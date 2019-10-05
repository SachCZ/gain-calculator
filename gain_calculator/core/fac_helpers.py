"""
Module containing useful abstractions to encapsulate FAC
"""
import itertools
import os, sys
import re
from pfac import fac
from pfac import crm
import classes


class Parser(object):
    """
    Context manager parsing FAC files and making them available in reasonable format

    :ivar transitions: FAC transitions file represented as list [(lower_index, upper_index, strength), ...]
    :ivar levels: FAC levels file represented as dict {name: index, ...}
    """

    def __init__(self, atom, principal_number_threshold):  # type: (classes.Atom, int) -> None
        """
        Initialize the Parser using an atom with given maximal principal quantum number (n)
        :param atom: Atom instance
        :param principal_number_threshold: maximum n to which configurations will be generated when calculating
            populations
        """

        self.__initialize_fac()
        self.__atom = atom
        self.__principal_number_threshold = principal_number_threshold

        if self.__is_serialized():
            self.__init_from_cache()
        else:
            self.__init()

    __fac_temp_folder = os.path.join(os.path.abspath(os.path.dirname(__file__)), "fac_temp")

    def get_all_populations(self, temperature, density):  # type: (float, float) -> dict
        """
        Generate all level populations in given atom with principal_number_threshold using temperature and density
        as params
        :param temperature: temperature of plasma in eV
        :param density: density of plasma in cm^-3
        :return: returns a dict of populations indexed by energy level index, eg. {0: 0.25e-3, 1: 0.33e-4, ...}
        """

        self.__generate_populations(temperature=temperature, density=density)
        populations = self.__parse_population_file()
        self.__clean_population_files()
        return populations

    def __get_dir_name(self):
        return "-".join([str(self.__atom), str(self.__principal_number_threshold)])

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

    def __generate_populations(self, temperature, density):  # type: (float, float) -> None
        electron_count = self.__atom.get_electron_count()
        dens = density * 1e-10

        crm.ReinitCRM()

        crm.AddIon(electron_count, 0.0, self.__binary_filename)
        crm.SetBlocks(-1)

        crm.SetEleDist(0, temperature, -1, -1)
        crm.SetTRRates(0)
        crm.SetCERates(1)
        crm.SetAbund(electron_count, 1.0)
        crm.SetEleDensity(dens)

        crm.InitBlocks()
        crm.SetIteration(1e-6, 0.5)

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

    def __generate_transitions(self, group_combinations):
        for group_combination in group_combinations:
            fac.TransitionTable(self.__transitions_binary_filename, group_combination[0], group_combination[1])
        fac.PrintTable(self.__transitions_binary_filename, self.__transitions_filename, 1)

    def __generate_excitation(self, groups):
        for group in groups:
            fac.CETable(self.__excitation_binary_filename, "base_group", group)
        fac.PrintTable(self.__excitation_binary_filename, self.__excitation_filename, 1)

    def __generate_files(self):  # type: () -> None
        fac.SetAtom(self.__atom.symbol)
        possible_configurations = self.__atom.get_possible_fac_configurations(self.__principal_number_threshold)
        self.__configure_ion(possible_configurations)

        groups = possible_configurations.keys()
        group_combinations = list(itertools.combinations_with_replacement(groups, 2))

        self.__generate_structure(groups)
        self.__generate_transitions(group_combinations)
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
    def __configure_ion(possible_configurations):
        for group in possible_configurations:
            fac.Config(group, possible_configurations[group])

        fac.ConfigEnergy(0)
        fac.OptimizeRadial(['base_group'])
        fac.ConfigEnergy(1)

    @staticmethod
    def __initialize_fac():
        fac.Reinit(0)
