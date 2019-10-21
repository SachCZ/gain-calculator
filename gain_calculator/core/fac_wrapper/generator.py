import importlib

import os
import itertools
from multiprocessing import Process, Queue


def __generate_files(queue, atom, fac_temp_folder):
    queue.put(Generator(atom, fac_temp_folder).get_files())


def generate_files(atom, fac_temp_folder):
    queue = Queue()
    process = Process(target=__generate_files, args=(queue, atom, fac_temp_folder))
    process.start()
    process.join()
    return queue.get()


class FacFiles:
    def __init__(
            self,
            binary_filename,
            excitation_binary_filename,
            hamiltonian_binary_filename,
            levels_binary_filename,
            transitions_binary_filename,
            levels_filename,
            excitation_filename,
            transitions_filename,
            dir_name
    ):
        self.binary_filename = binary_filename
        self.excitation_binary_filename = excitation_binary_filename
        self.hamiltonian_binary_filename = hamiltonian_binary_filename
        self.levels_binary_filename = levels_binary_filename
        self.transitions_binary_filename = transitions_binary_filename
        self.levels_filename = levels_filename
        self.excitation_filename = excitation_filename
        self.transitions_filename = transitions_filename
        self.dir_name = dir_name


class Generator:
    """
    Class for generating FAC files and making them available in reasonable format

    :ivar Atom atom: atom instance to generate files for
    """

    def __init__(self, atom, fac_temp_folder):
        """
        Initialize the Parser using an atom with given maximal principal quantum number (n)
        :param atom: Atom instance
        """

        self.__atom = atom
        self.__fac_temp_folder = fac_temp_folder

        if self.__is_serialized():
            self.__init_from_cache()
        else:
            self.__init()

    def get_files(self):
        # TODO check everything is allright and ready to go
        return self.__files

    def __generate_structure(self, groups):
        self.fac.Structure(self.__files.levels_binary_filename, self.__files.hamiltonian_binary_filename, groups)
        self.fac.MemENTable(self.__files.levels_binary_filename)
        self.fac.PrintTable(self.__files.levels_binary_filename, self.__files.levels_filename, 1)

    def __generate_transitions(self, group_combinations):
        for lower, upper in group_combinations:
            self.fac.TransitionTable(self.__files.transitions_binary_filename, lower.get_name(), upper.get_name())
        self.fac.PrintTable(self.__files.transitions_binary_filename, self.__files.transitions_filename, 1)

    def __generate_excitation(self, group_combinations):
        for lower, upper in group_combinations:
            self.fac.CETable(self.__files.excitation_binary_filename, lower.get_name(), upper.get_name())
        self.fac.PrintTable(self.__files.excitation_binary_filename, self.__files.excitation_filename, 1)

    def __generate_files(self):
        self.fac.SetAtom(self.__atom.symbol)
        config_groups = self.__atom.config_groups
        self.__configure_ion(config_groups)

        group_names = config_groups.get_names()
        group_combinations = self.__generate_group_combinations(config_groups)

        self.__generate_structure(group_names)
        self.__generate_transitions(group_combinations)
        self.__generate_excitation(group_combinations)

    def __configure_ion(self, config_groups):
        for config_group in config_groups.all_groups:
            self.fac.Config(config_group.get_name(), config_group.config)

        self.fac.ConfigEnergy(0)
        self.fac.OptimizeRadial(config_groups.base_group.get_name())
        self.fac.ConfigEnergy(1)

    def __initialize_fac(self):
        self.fac.Reinit(0)

    def __init(self):
        # This monstrosity is here because FAC allocates 1.6 GB of memory when imported
        self.fac = importlib.import_module("pfac.fac")

        self.__initialize_fac()

        self.__dir_name = self.__create_dir()
        self.__init_filenames()
        self.__generate_files()

    def __create_dir(self):
        dir_name = self.__get_dir_name()
        dir_path = os.path.join(self.__fac_temp_folder, dir_name)
        try:
            os.mkdir(dir_path)
        except OSError as e:
            raise Exception("Failed to create directory to hold FAC files: {}".format(e.strerror))

        return dir_path

    def __init_filenames(self):
        binary_filename = os.path.join(self.__dir_name, "fac_binary_temp")
        self.__files = FacFiles(
            binary_filename=binary_filename,
            excitation_binary_filename=binary_filename + ".ce",
            hamiltonian_binary_filename=binary_filename + ".ham",
            levels_binary_filename=binary_filename + ".en",
            transitions_binary_filename=binary_filename + ".tr",
            levels_filename=os.path.join(self.__dir_name, "levels.txt"),
            excitation_filename=os.path.join(self.__dir_name, "excitation.txt"),
            transitions_filename=os.path.join(self.__dir_name, "transitions.txt"),
            dir_name=self.__dir_name
        )

    def __get_dir_name(self):
        return " up to n=".join([str(self.__atom), str(self.__atom.config_groups.get_max_n())])

    @staticmethod
    def __generate_group_combinations(config_groups):
        group_combinations = list(itertools.combinations_with_replacement(config_groups.all_groups, 2))

        def __fix_invalid(combination):
            first, second = combination
            if first > second:
                return second, first
            else:
                return combination

        return map(__fix_invalid, group_combinations)

    def __init_from_cache(self):
        self.__dir_name = os.path.join(self.__fac_temp_folder, self.__get_dir_name())
        self.__init_filenames()

    def __is_serialized(self):
        try:
            return self.__get_dir_name() in os.listdir(self.__fac_temp_folder)
        except OSError:
            return False
