import unittest
import gain_calculator.core as core
import copy


class TestShell(unittest.TestCase):
    def setUp(self):
        self.shell = core.Shell(
            n=2,
            l=1,
            spin_direction="+",
            double_angular_momentum=0,
            electron_count=4
        )

    def tearDown(self):
        del self.shell

    def test_string_init(self):
        another_shell = core.Shell.create_from_string("2p+4(0)")
        self.assertEqual(another_shell, self.shell)

    def test_repr(self):
        self.assertEqual("2p+4(0)", str(self.shell))

    def test_latex_repr(self):
        self.assertEqual("$[2p+]^4_0$", self.shell.get_latex_repr())

    def test_equals(self):
        another_shell = copy.deepcopy(self.shell)
        self.assertEqual(self.shell, another_shell)

    def test_is_full(self):
        self.assertEqual(True, self.shell.is_full())


class TestConfigGroups(unittest.TestCase):
    def setUp(self):
        self.config_groups = core.ConfigGroups(base="1*2 2*8", max_n=5)

    def tearDown(self):
        del self.config_groups

    def test_all_groups(self):
        expected = {core.ConfigGroup(5, '1*2 2*7 5*1'),
                    core.ConfigGroup(4, '1*2 2*7 4*1'),
                    core.ConfigGroup(3, '1*2 2*7 3*1'),
                    core.ConfigGroup(0, '1*2 2*8')}

        self.assertSetEqual(expected, self.config_groups.all_groups)

    def test_base_group(self):
        self.assertEqual(core.ConfigGroup(0, '1*2 2*8'), self.config_groups.base_group)

    def test_base_group_electron_count(self):
        self.assertEqual(10, self.config_groups.base_group.get_electron_count())


class TestAtom(unittest.TestCase):
    def setUp(self):
        self.atom = core.Atom(
            symbol="Ge",
            config_groups=core.ConfigGroups(base="1*2 2*8", max_n=3),
        )

    def tearDown(self):
        del self.atom

    def test_get_population(self):
        self.assertAlmostEqual(0.0071, self.atom.get_population(
            energy_level=core.EnergyLevel("1s+2(0)0 2s+2(0)0 2p-2(0)0 2p+3(3)3 3s+1(1)4"),
            temperature=900,
            electron_density=1e20,
        ), places=4)


class TestTransition(unittest.TestCase):
    def setUp(self):
        atom = core.Atom(
            symbol="Ge",
            config_groups=core.ConfigGroups(base="1*2 2*8", max_n=3),
        )
        self.transition = core.Transition(
            atom=atom,
            lower=core.EnergyLevel("1s+2(0)0 2s+2(0)0 2p-1(1)1 2p+4(1)1 3s+1(1)2"),
            upper=core.EnergyLevel("1s+2(0)0 2s+2(0)0 2p-1(1)1 2p+4(6)1 3p+1(3)4"),
        )

    def tearDown(self):
        del self.transition

    def test_get_weighted_oscillator_strength(self):
        self.assertAlmostEqual(0.52, self.transition.weighted_oscillator_strength, places=2)

    def test_get_population(self):
        expected = {"lower": 0.00056, "upper": 0.00295}
        populations = self.transition.get_populations(
            temperature=900,
            electron_density=1e20
        )
        self.assertAlmostEqual(expected['lower'], populations['lower'], places=4)
        self.assertAlmostEqual(expected['upper'], populations['upper'], places=4)


class TestEnergyLevel(unittest.TestCase):
    def setUp(self):
        self.energy_level = core.EnergyLevel("1s+2(0)0 2s+2(0)0 2p-2(0)0 2p+3(3)3 3s+1(1)4")

    def tearDown(self):
        del self.energy_level

    def test_equals(self):
        another_level = copy.deepcopy(self.energy_level)
        self.assertEqual(self.energy_level, another_level)

    def test_string_init(self):
        self.assertEqual("1s+2(0)0 2s+2(0)0 2p-2(0)0 2p+3(3)3 3s+1(1)4", str(self.energy_level))

    def test_configuration(self):
        self.assertSequenceEqual([
            core.LevelTerm(shell=core.Shell.create_from_string("1s+2(0)"), j2=0),
            core.LevelTerm(shell=core.Shell.create_from_string("2s+2(0)"), j2=0),
            core.LevelTerm(shell=core.Shell.create_from_string("2p-2(0)"), j2=0),
            core.LevelTerm(shell=core.Shell.create_from_string("2p+3(3)"), j2=3),
            core.LevelTerm(shell=core.Shell.create_from_string("3s+1(1)"), j2=4)
        ], self.energy_level.configuration)

    def test_get_degeneracy(self):
        self.assertEqual(5, self.energy_level.degeneracy)

    def test_get_fac_repr(self):
        self.assertEqual("2p+3(3)3.3s+1(1)4", self.energy_level.get_fac_repr())


if __name__ == '__main__':
    unittest.main()
