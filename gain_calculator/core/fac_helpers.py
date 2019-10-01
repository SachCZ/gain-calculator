"""
Module holding containing useful abstractions to encapsulate FAC
"""

import os
import re
from pfac import fac

import classes


class Parser(object):
    """
    Context manager parsing FAC files and making them available in reasonable format

    :ivar structure: FAC structure file represented as list [(lower_index, upper_index, strength), ...]
    :ivar levels: FAC levels file represented as dict {name: index, ...}
    """

    def __init__(self, atom, lower, upper):  # type: (str, classes.EnergyLevel, classes.EnergyLevel) -> None
        self.__generate_files(atom, lower, upper)
        self.levels = self.__parse_levels_file()
        self.structure = self.__parse_structure_file()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.__clean_files()

    __structure_filename = "structure.txt"
    __levels_filename = "levels.txt"
    __structure_binary_filename = "structure.b"
    __levels_binary_filename = "levels.b"

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

    def __parse_structure_file(self):
        def __parse_line(line):
            match = re.search(
                r"\s*(?P<upper_index>\d)\s+\d\s+(?P<lower_index>\d)\s+\d\s+[\d.+-E]*\s+(?P<strength>[\d.+-E]*)",
                line)
            if match is None:
                return None
            return int(match.group('lower_index')), int(match.group('upper_index')), float(match.group('strength'))

        with open(self.__structure_filename, 'r') as f:
            structure = filter(None, map(__parse_line, f.readlines()))
            assert structure, "Fatal error, no structure parsed"
            return structure

    def __generate_files(self, atom, lower, upper):
        fac.SetAtom(atom)
        fac.Config(lower.get_fac_style_configuration(), group='lower')
        fac.Config(upper.get_fac_style_configuration(), group='upper')
        fac.ConfigEnergy(0)
        fac.OptimizeRadial(['lower', 'upper'])
        fac.ConfigEnergy(1)
        fac.Structure(self.__levels_binary_filename, ['lower', 'upper'])
        fac.MemENTable(self.__levels_binary_filename)
        fac.PrintTable(self.__levels_binary_filename, self.__levels_filename, 1)
        fac.TransitionTable(self.__structure_binary_filename, ['lower'], ['upper'])
        fac.PrintTable(self.__structure_binary_filename, self.__structure_filename, 1)

    def __clean_files(self):
        try:
            os.remove(self.__structure_filename)
            os.remove(self.__levels_filename)
            os.remove(self.__structure_binary_filename)
            os.remove(self.__levels_binary_filename)
        except OSError as e:
            raise Exception("Failed to clean FAC file {} with following error: {}".format(e.filename, e.strerror))
