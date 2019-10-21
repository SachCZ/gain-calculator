"""
Module containing useful abstractions to encapsulate FAC
"""
import re
import os
import uuid
import ray
from pfac import crm
import generator


@ray.remote(num_cpus=1)
class Parser(object):
    def __init__(self, files, electron_count):  # type: (generator.FacFiles, int) -> None
        self.__files = files
        self.__electron_count = electron_count
        self.levels = self.__parse_levels_file()
        self.transitions = self.__parse_transitions_file()

    def get_all_populations(self, temperature, density, population_total):  # type: (float, float, float) -> dict
        """
        Generate all level populations in given atom with max_n using temperature and density
        as params
        :param population_total: The sum of all populations over all levels
        :param temperature: temperature of plasma in eV
        :param density: density of plasma in cm^-3
        :return: returns a dict of populations indexed by energy level index, eg. {0: 0.25e-3, 1: 0.33e-4, ...}
        """

        self.__choose_population_filenames()
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

    def __parse_levels_file(self):
        def __parse_line(line):  # type: (str) -> (int, str)
            match = re.search(r"\s+(?P<index>\d+)\s+\S+\s+\S+\s+\d+\s+\d+\s+\d+\s+\S+\s+\S+\s+(?P<name>\S+)", line)
            if match is None:
                return None
            return str(match.group('name')), int(match.group('index'))

        with open(self.__files.levels_filename, 'r') as f:
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

        with open(self.__files.transitions_filename, 'r') as f:
            transitions = filter(None, map(__parse_line, f.readlines()))
            assert transitions, "Fatal error, no transitions parsed"
            return transitions

    def __choose_population_filenames(self):
        name = uuid.uuid4().hex[:6].upper()
        self.__spec_binary_filename = name + ".sp"
        self.__spec_filename = name + ".txt"

    def __generate_populations(self, temperature, density, population_total):  # type: (float, float, float) -> None

        # Keeping this in one method as I find it easier to manage FAC in one place
        electron_count = self.__electron_count
        crm.ReinitCRM()
        crm.NormalizeMode(1)
        crm.AddIon(electron_count, 0, self.__files.binary_filename)
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

    def __clean_population_files(self):
        os.remove(self.__spec_binary_filename)
        os.remove(self.__spec_filename)
