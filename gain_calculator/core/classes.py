"""
Module containing the core classes of the GainCalculator project. End user should not be concerned about this.
A convention is used that principal quantum number is denoted n and orbital quantum number is denoted l, keep
this in mind.
"""

import re
import fac_helpers


class ConfigGroup:
    def __init__(self, index, config):  # type: (int, str) -> None
        self.index = index
        self.config = config

    def __repr__(self):
        return "[{}] ".format(self.index) + self.config

    def __eq__(self, other):
        return self.index == other.index and self.config == other.config

    def __lt__(self, other):
        return self.index < other.index

    def __le__(self, other):
        return self.index <= other.index

    def __gt__(self, other):
        return self.index > other.index

    def __ge__(self, other):
        return self.index >= other.index

    def __hash__(self):
        return hash(str(self))

    def get_name(self):
        return "group" + str(self.index)

    def get_electron_count(self):
        return sum([int(electrons) for electrons in re.findall(r'\*(\d+?)', self.config)])


class ConfigGroups:
    """
    This class represents a set of all  possible shell configurations up to max_n concerning only the
    number of electrons in each shell given by principal quantum number n. There is no end user use of this
    class aside from providing its instance as a param to other classes. Initialize it like this::

        ConfigGroups(
            base="1*2 2*8",
            max_n=4
        )

    :param str base: The base configuration from which all config groups are generated. Notation "1*2 2*8" means there
        are 2 electron in shell with :math:`n=1` and 8 electrons in shell with :math:`n=2`. Notation with asterisks
        and spaces must be used at all times.
    :param int max_n: Maximal principal quantum number :math:`n` up to which
        ConfigGroups are represented::

            ConfigGroups(base=1*2 2*8, max_n=4)

        represents::

            ["1*2 2*8", "1*2 2*7 3*1", "1*2 2*7 4*1"]
    """
    def __init__(self, base, max_n):  # type: (str, int) -> ConfigGroups
        self.base_group = ConfigGroup(0, base)
        begin_n = int(base[-3])
        base = base[:-1] + str(int(base[-1]) - 1)  # Subtract one from the last number in string
        self.all_groups = {ConfigGroup(n, base + " {}*1".format(n)) for n in range(begin_n + 1, max_n + 1)}
        self.all_groups.add(self.base_group)

    def get_names(self):
        return [group.get_name() for group in self.all_groups]

    def get_max_n(self):
        return max([group.index for group in self.all_groups])


class Atom:
    """
    This class represents a specific atom with given configuration. Initialize it like this::

        Atom(
            symbol="Fe",
            config_groups=ConfigGroups(
                base="1*2 2*8",
                max_n=4
            )
        )

    :param str symbol: A symbol of the element, eg. "Ge" for germanium
    :param ConfigGroups config_groups: Instance of class ConfigGroups representing possible shell configurations.
        For more info see :class:`ConfigGroups`.
    """
    def __init__(self, symbol, config_groups):  # type: (str, ConfigGroups) -> None
        self.symbol = symbol
        self.config_groups = config_groups
        self.electron_count = config_groups.base_group.get_electron_count()

    def __repr__(self):
        return " ".join([self.symbol, str(self.config_groups.base_group)])

    def get_population(self, energy_level, temperature, electron_density, population_total=1.0):
        """
        Get electron population on single energy level at given temperature and electron density.
        The calculation is done for all energy levels in given atom. Population
        total is the sum of populations of all energy levels. By default this is 1 and the function returns relative
        populations. Eg. 0.1 means 10% of all electrons are at energy_level. Example usage::

            atom = Atom(
                symbol="Fe",
                config_groups=ConfigGroups(
                    base="1*2 2*8",
                    max_n=4
                )
            )
            population = atom.get_population(
                energy_level=EnergyLevel(config="1s+2(0)0 2s+2(0)0 2p-2(0)0 2p+3(3)3 3s+1(1)4"),
                temperature=900.0  # eV
                electron_density=1e20  # cm^-3
            )

        :param EnergyLevel energy_level: the energy level instance
        :param float temperature: temperature in eV
        :param float electron_density: electron density in cm^-3
        :param float population_total: 1 by default
        :return: energy level electron population normalized by population_total
        """
        fac_parser = fac_helpers.Parser(self)
        return fac_parser.get_population(energy_level, temperature, electron_density, population_total)


class LevelTerm:
    """
    Simple comparable data class representing a term of form 2p+(1)1
    where the last number is 2*J (J = total angular momentum over all shells).
    The rest of the shell notation is documented elsewhere.

    :ivar shell: shell object instance
    :ivar j2: double the value of total angular momentum
    """

    def __init__(self, shell, j2):  # type: (Shell, int) -> None
        """
        Init using a shell and total angular momentum
        :param shell: shell object instance
        :param j2: double the value of total angular momentum
        """
        self.shell = shell
        self.j2 = j2
        assert j2 >= 0, "Total angular momentum (J2) in {} is negative".format(self)

    def __repr__(self):  # type: () -> str
        return str(self.shell) + str(self.j2)

    def __eq__(self, other):
        return str(self) == str(other)

    def __ne__(self, other):  # type: (LevelTerm, LevelTerm) -> bool
        return not self == other


class Shell:
    """
    Class representing a single shell of form 3d+4(2), where 3 is the principal quantum number
    d is the orbital quantum number (2) + represents spin direction (up),
    4 is electron count and 2 is double the value of angular momentum of the shell.

    :ivar n: Starting from 1
    :ivar l: Given by value, use Shell.orbital_numbers for convenience
    :ivar spin_direction: expects "+" or "-" representing spin up and down
    :ivar double_angular_momentum: 2*J value (of shell only) e.g. two electrons up give 2*J=2
    :ivar electron_count: number of electrons in the shell

    :cvar orbital_numbers: mapping of orbital letters to numbers
    :cvar orbital_letters: inverse mapping to orbital numbers
    :cvar spin_direction_sign: mapping from "-" to -1 and "+" to 1
    """

    def __init__(self, n, l, spin_direction,
                 double_angular_momentum, electron_count):  # type: (int, int, str, int, int) -> None
        """
        Initialize it based on quantum numbers spin and electron count. Consider using
        a factory method create from string for convenience.
        :param n: Starting from 1
        :param l: Given by value, use Shell.orbital_numbers for convenience
        :param spin_direction: expects "+" or "-" representing spin up and down
        :param double_angular_momentum: 2*J value (of shell only) e.g. two electrons up give 2*J=2
        :param electron_count: number of electrons in the shell
        """
        self.n = n
        self.l = l
        self.spin_direction = spin_direction
        self.double_angular_momentum = double_angular_momentum
        self.electron_count = electron_count

        self.__assert_state_is_fine(), "Invalid shell configuration {}".format(str(self))

    def __repr__(self):  # type: () -> str
        return str(self.n) + \
               self.orbital_letters[self.l] + \
               self.spin_direction + \
               str(self.electron_count) + \
               "(" + str(self.double_angular_momentum) + ")"

    def __eq__(self, other):  # type: (Shell, Shell) -> bool
        return str(self) == str(other)

    def __ne__(self, other):  # type: (Shell, Shell) -> bool
        return not self == other

    def get_latex_repr(self):  # type: () -> str
        """
        Method to generate shell representation in latex format
        :return: Representation of Shell analogous to "$[2p+]^4_0$"
        """
        return "$[{}{}{}]^{}_{}$".format(
            self.n,
            self.orbital_letters[self.l],
            self.spin_direction,
            self.electron_count,
            self.double_angular_momentum
        )

    def is_full(self):
        """
        Returns true if there is maximum number of electrons in the shell
        :return:
        """
        return self.electron_count == self.__get_max_electron_count()

    def __get_max_electron_count(self):
        j2 = 2 * self.l + Shell.spin_direction_sign[self.spin_direction]
        return j2 + 1

    def __get_max_double_angular_momentum(self):
        return 2 * self.l + self.__get_max_electron_count()

    @staticmethod
    def create_from_string(shell_repr):  # type: (str) -> Shell
        """
        Factory method to deconstruct a string of form 3d+4(2) and create a Shell instance based on it
        :param shell_repr: string to deconstruct
        :return: instance of Shell
        """
        orbitals = reduce(lambda string_so_far, orbital: string_so_far + orbital, Shell.orbital_letters.itervalues())
        pattern = re.compile(r"^\d[" + orbitals + r"][+-]\d\(\d\)$")
        assert pattern.match(shell_repr) is not None, "Invalid string to generate shell {}".format(shell_repr)

        return Shell(
            n=int(shell_repr[0]),
            l=Shell.orbital_numbers[shell_repr[1]],
            spin_direction=shell_repr[2],
            double_angular_momentum=int(shell_repr[5]),
            electron_count=int(shell_repr[3])
        )

    orbital_numbers = {"s": 0, "p": 1, "d": 2, "f": 3, "g": 4, "h": 5, "i": 6}
    orbital_letters = {number: letter for letter, number in orbital_numbers.iteritems()}
    spin_direction_sign = {"+": 1, "-": -1}

    def __assert_state_is_fine(self):
        assert self.n > 0, \
            "Principal quantum number (n) must be positive in {}".format(self)
        assert abs(self.l) < self.n, \
            "Absolute value of orbital quantum number must be smaller than principal quantum number in {}".format(self)
        assert self.spin_direction == "+" or self.spin_direction == "-", \
            "Spin direction must be either + or - in {}".format(self)
        assert self.electron_count <= self.__get_max_electron_count(), \
            "Electron count must be less or equal to {} in shell {}".format(self.__get_max_electron_count(), self)
        assert self.double_angular_momentum <= self.__get_max_double_angular_momentum(), \
            "Double angular momentum must be less or equal to {} in shell {}".format(
                self.__get_max_double_angular_momentum(), self)


class EnergyLevel:
    """
    This class represents a single energy level.
    Initialize it like this::

        EnergyLevel(config="1s+2(0)0 2s+2(0)0 2p-2(0)0 2p+3(3)3 3s+1(1)4")

    :param str config: A sequence of terms separated by space each representing single shell in jj coupling
        notation. A shell 5p+3(1)2 has principal quantum number :math:`n=5`, orbital quantum number :math:`l = p`.
        There are three electrons in this shell. The number in brackets is 2 times the total
        angular momentum :math:`2J`. Here :math:`2J` is 2 hence :math:`J` is :math:`\\frac{1}{2}`.
        The last number is 2 times total angular momentum when taking into account all previous
        shells here it is 2 meaning total angular momentum of the whole configuration up to this shell is 1.

    :ivar int degeneracy: A number representing the degeneracy of the level. It is usually denoted :math:`g`.
    """

    def __init__(self, config):  # type: (str) -> None
        self.configuration = map(lambda term_string: LevelTerm(
            shell=Shell.create_from_string(term_string[:-1]),
            j2=int(term_string[-1])
        ), config.split(" "))

        self.degeneracy = self.configuration[-1].j2 + 1

    def __repr__(self):  # type: () -> str
        return " ".join(map(lambda term: str(term.shell) + str(term.j2), self.configuration))

    def __eq__(self, other):  # type: (EnergyLevel, EnergyLevel) -> bool
        return str(self) == str(other)

    def __ne__(self, other):  # type: (EnergyLevel, EnergyLevel) -> bool
        return not self == other

    def get_fac_repr(self):
        """
        Returns a string of terms separated by dot. The terms are string representations of LevelTerm. Only terms
        with shells that are not full are listed
        :return:
        """
        return ".".join([str(term) for term in self.configuration if not term.shell.is_full()])


class Transition:
    """
    This class represents a transition between upper and lower energy levels.
    Initialize it like this::

        atom = Atom(
            symbol="Fe",
            config_groups=ConfigGroups(
                base="1*2 2*8",
                max_n=4
            )
        )
        Transition(
            atom=atom,
            lower=core.EnergyLevel("1s+2(0)0 2s+2(0)0 2p-1(1)1 2p+4(1)1 3s+1(1)2"),
            upper=core.EnergyLevel("1s+2(0)0 2s+2(0)0 2p-1(1)1 2p+4(6)1 3p+1(3)4")
        )

    :param Atom atom: The atom instance in which te transition occurs
    :param EnergyLevel lower: The lower energy level
    :param EnergyLevel upper: The upper energy level

    :ivar float weighted_oscillator_strength: The weighted oscillator strength of the transition :math:`gf`.
    """

    def __init__(self, atom, lower, upper):
        # type: (Atom, EnergyLevel, EnergyLevel) -> None
        self.lower = lower
        self.upper = upper
        self.atom = atom
        self.__fac_parser = fac_helpers.Parser(self.atom)
        self.weighted_oscillator_strength = self.__fac_parser.get_weighted_oscillator_strength(lower, upper)

    def get_populations(self, temperature, electron_density, population_total=1.0):
        # type: (float, float, float) -> {"lower": float, "upper": float}
        """
        Simple convenience wrapper around
        :func:`EnergyLevel.get_population`
        to get both lower and upper energy
        levels populations.

        :param float temperature: temperature in eV
        :param float electron_density: electron density in cm^-3
        :param float population_total: 1 by default
        :return: a dict with two keys: {upper: ..., lower: ...} - the values are the populations normalized by
            population_total
        """
        return {"lower": self.__fac_parser.get_population(self.lower, temperature, electron_density, population_total),
                "upper": self.__fac_parser.get_population(self.upper, temperature, electron_density, population_total)}
