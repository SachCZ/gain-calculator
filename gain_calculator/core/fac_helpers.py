"""
Module containing useful abstractions to encapsulate FAC
"""
import itertools
import re
import os
from pfac import fac
from pfac import crm
import classes


class Parser(object):
    """
    Context manager parsing FAC files and making them available in reasonable format

    :ivar transitions: FAC transitions file represented as list [(lower_index, upper_index, strength), ...]
    :ivar levels: FAC levels file represented as dict {name: index, ...}
    """

    def __init__(self, atom):  # type: (classes.Atom) -> None
        """
        Initialize the Parser using an atom with given maximal principal quantum number (n)
        :param atom: Atom instance
        """

        self.__initialize_fac()
        self.__atom = atom

        if self.__is_serialized():
            self.__init_from_cache()
        else:
            self.__init()

    __fac_temp_folder = os.path.join(os.path.abspath(os.path.dirname(__file__)), "fac_temp")

    def get_all_populations(self, temperature, density, population_total):  # type: (float, float, float) -> dict
        """
        Generate all level populations in given atom with max_n using temperature and density
        as params
        :param population_total: The sum of all populations over all levels
        :param temperature: temperature of plasma in eV
        :param density: density of plasma in cm^-3
        :return: returns a dict of populations indexed by energy level index, eg. {0: 0.25e-3, 1: 0.33e-4, ...}
        """

        self.__generate_populations(temperature=temperature, density=density, population_total=population_total)
        populations = self.__parse_population_file()
        self.__clean_population_files()
        return populations

    def get_weighted_oscillator_strength(self, lower, upper):
        # type: (classes.EnergyLevel, classes.EnergyLevel) -> float
        """
        Calculate the weighted oscillator strength (gf) between two energy levels
        :param lower: lower EnergyLevel instance
        :param upper: upper EnergyLevel instance
        :return: weighted oscillator strength value
        """
        return self.__parse_oscillator_strength(lower, upper)

    def __get_dir_name(self):
        return "-".join([str(self.__atom), str(self.__atom.config_groups.get_max_n())])

    def __create_dir(self):
        dir_name = self.__get_dir_name()
        dir_path = os.path.join(self.__fac_temp_folder, dir_name)
        try:
            os.mkdir(dir_path)
        except OSError as e:
            raise Exception("Failed to create directory to hold FAC files: {}".format(e.strerror))

        return dir_path

    def __init_filenames(self):
        self.__binary_filename = os.path.join(self.__dir_name, "fac_binary_temp")

        self.__excitation_binary_filename = self.__binary_filename + ".ce"
        self.__hamiltonian_binary_filename = self.__binary_filename + ".ham"
        self.__levels_binary_filename = self.__binary_filename + ".en"
        self.__transitions_binary_filename = self.__binary_filename + ".tr"
        self.__spec_binary_filename = self.__binary_filename + ".sp"

        self.__levels_filename = os.path.join(self.__dir_name, "levels.txt")
        self.__excitation_filename = os.path.join(self.__dir_name, "excitation.txt")
        self.__transitions_filename = os.path.join(self.__dir_name, "transitions.txt")
        self.__spec_filename = os.path.join(self.__dir_name, "populations.txt")

    def __clean_population_files(self):
        os.remove(self.__spec_binary_filename)
        os.remove(self.__spec_filename)

    def __generate_populations(self, temperature, density, population_total):  # type: (float, float, float) -> None

        # Keeping this in one method as I find it easier to manage FAC in one place
        electron_count = self.__atom.electron_count
        crm.ReinitCRM()
        crm.NormalizeMode(1)
        crm.AddIon(electron_count, 0, self.__binary_filename)
        crm.SetBlocks(-1)

        crm.SetAbund(electron_count, population_total)
        crm.SetEleDensity(density * 1e-10)
        crm.SetEleDist(0, temperature, -1, -1)
        crm.SetTRRates(0)
        crm.SetCERates(1)

        crm.InitBlocks()
        crm.SetIteration(1e-4, 0.5, 1500)

        crm.LevelPopulation()

        crm.SpecTable(self.__spec_binary_filename, -1)
        crm.PrintTable(self.__spec_binary_filename, self.__spec_filename)

    def __parse_levels_file(self):
        def __parse_line(line):  # type: (str) -> (int, str)
            match = re.search(r"\s+(?P<index>\d+)\s+\S+\s+\S+\s+\d+\s+\d+\s+\d+\s+\S+\s+\S+\s+(?P<name>\S+)", line)
            if match is None:
                return None
            return str(match.group('name')), int(match.group('index'))

        with open(self.__levels_filename, 'r') as f:
            levels = filter(None, map(__parse_line, f.readlines()))
            assert levels, "Fatal error, no levels parsed"
            return {level[0]: level[1] for level in levels}

    def __parse_transitions_file(self):
        def __parse_line(line):
            match = re.search(
                r"\s*(?P<upper_index>\d+)\s+\d+\s+(?P<lower_index>\d+)\s+\d+\s+\S+\s+(?P<strength>\S+)",
                line)
            if match is None:
                return None
            return int(match.group('lower_index')), int(match.group('upper_index')), float(match.group('strength'))

        with open(self.__transitions_filename, 'r') as f:
            transitions = filter(None, map(__parse_line, f.readlines()))
            assert transitions, "Fatal error, no transitions parsed"
            return transitions

    def __parse_population_file(self):
        def __parse_line(line):
            match = re.search(
                r"\s+(?P<index>\d+)\s+\d+\s+\S+\s+(?P<population>\S+)",
                line)
            if match is None:
                return None
            return int(match.group('index')), float(match.group('population'))

        with open(self.__spec_filename, 'r') as f:
            populations = filter(None, map(__parse_line, f.readlines()))
            assert populations, "Fatal error, no populations parsed"
            return {index: population for index, population in populations}

    def __generate_structure(self, groups):
        fac.Structure(self.__levels_binary_filename, self.__hamiltonian_binary_filename, groups)
        fac.MemENTable(self.__levels_binary_filename)
        fac.PrintTable(self.__levels_binary_filename, self.__levels_filename, 1)

    @staticmethod
    def __generate_group_combinations(groups):
        # TODO rewrite this - create config_group __gt__ etc

        group_combinations = list(itertools.combinations_with_replacement(groups, 2))

        def __fix_invalid(combination):
            if int(combination[0][-1]) > int(combination[1][-1]):
                return combination[1], combination[0]
            else:
                return combination

        return map(__fix_invalid, group_combinations)

    def __generate_transitions(self, groups):
        for group_combination in self.__generate_group_combinations(groups):
            fac.TransitionTable(self.__transitions_binary_filename, group_combination[0], group_combination[1])
        fac.PrintTable(self.__transitions_binary_filename, self.__transitions_filename, 1)

    def __generate_excitation(self, groups):
        for group_combination in self.__generate_group_combinations(groups):
            fac.CETable(self.__excitation_binary_filename, group_combination[0], group_combination[1])
        fac.PrintTable(self.__excitation_binary_filename, self.__excitation_filename, 1)

    def __generate_files(self):  # type: () -> None
        fac.SetAtom(self.__atom.symbol)
        config_groups = self.__atom.config_groups
        self.__configure_ion(config_groups)

        groups = config_groups.get_names()

        self.__generate_structure(groups)
        self.__generate_transitions(groups)
        self.__generate_excitation(groups)

    def __init(self):
        self.__dir_name = self.__create_dir()
        self.__init_filenames()
        self.__generate_files()
        self.levels = self.__parse_levels_file()
        self.transitions = self.__parse_transitions_file()

    def __init_from_cache(self):
        self.__dir_name = os.path.join(self.__fac_temp_folder, self.__get_dir_name())
        self.__init_filenames()
        self.levels = self.__parse_levels_file()
        self.transitions = self.__parse_transitions_file()

    def __is_serialized(self):
        return self.__get_dir_name() in os.listdir(self.__fac_temp_folder)

    @staticmethod
    def __configure_ion(config_groups):  # type: (classes.ConfigGroups) -> None
        for config_group in config_groups.all_groups:
            fac.Config(config_group.get_name(), config_group.config)

        fac.ConfigEnergy(0)
        fac.OptimizeRadial(config_groups.base_group.get_name())
        fac.ConfigEnergy(1)

    @staticmethod
    def __initialize_fac():
        fac.Reinit(0)

    def __get_level_index(self, energy_level):  # type: (classes.EnergyLevel) -> int
        return self.levels[energy_level.get_fac_repr()]

    def __indexes_match(self, lower_index, upper_index, lower_level, upper_level):
        #  type: (int, int, classes.EnergyLevel, classes.EnergyLevel) -> bool
        lower = self.__get_level_index(lower_level)
        upper = self.__get_level_index(upper_level)
        return lower == lower_index and upper == upper_index

    def __parse_oscillator_strength(self, lower_level, upper_level):
        # type: (classes.EnergyLevel, classes.EnergyLevel) -> float
        for lower_index, upper_index, strength in self.transitions:
            if self.__indexes_match(lower_index, upper_index, lower_level, upper_level):
                return strength

        raise Exception("Failed to find transition!")

    def get_population(self, energy_level, temperature, density, population_total):
        #  type: (classes.EnergyLevel, float, float, float) -> float
        """
        Calculate population of a single energy level at given temperature density with population total stating
        the sum of all energy levels.

        :param energy_level:
        :param temperature: temperature in eV
        :param density: electron density in cm^-3
        :param population_total:
        :return:
        """
        populations = self.get_all_populations(
            temperature=temperature,
            density=density,
            population_total=population_total)
        return populations[self.__get_level_index(energy_level)]
