import unittest
import gain_calculator.core as core
import copy


class TestCoreClasses(unittest.TestCase):
    def setUp(self):
        self.shell = core.Shell(
            n=2,
            l=1,
            spin_direction="+",
            double_angular_momentum=0,
            electron_count=4
        )
        self.energy_level = core.EnergyLevel.create_from_string("1s+2(0)0 2s+2(0)0 2p-2(0)0 2p+3(3)3 3s+1(1)4")
        self.germanium = core.Atom(
            symbol="Ge",
            base_level=core.EnergyLevel.create_from_string("1s+2(0)0 2s+2(0)0 2p-2(0)0 2p+4(0)0"),
        )
        self.transition = core.Transition(
            atom=self.germanium,
            lower=core.EnergyLevel.create_from_string("1s+2(0)0 2s+2(0)0 2p-1(1)1 2p+4(1)1 3s+1(1)2"),
            upper=core.EnergyLevel.create_from_string("1s+2(0)0 2s+2(0)0 2p-1(1)1 2p+4(6)1 3p+1(3)4"),
            principal_number_threshold=8
        )

    def tearDown(self):
        del self.shell
        del self.energy_level
        del self.germanium
        del self.transition

    def test_atom_get_possible_fac_configuration(self):
        expected = {'group5': '1*2 2*7 5*1', 'group4': '1*2 2*7 4*1', 'group3': '1*2 2*7 3*1', 'base_group': '1*2 2*8'}
        self.assertDictEqual(expected, self.germanium.get_possible_fac_configurations(5))

    def test_shell_string_init(self):
        another_shell = core.Shell.create_from_string("2p+4(0)")
        self.assertEqual(another_shell, self.shell)

    def test_shell_repr(self):
        self.assertEqual("2p+4(0)", str(self.shell))

    def test_shell_latex_repr(self):
        self.assertEqual("$[2p+]^4_0$", self.shell.get_latex_repr())

    def test_shell_equals(self):
        another_shell = copy.deepcopy(self.shell)
        self.assertEqual(self.shell, another_shell)

    def test_energy_level_equals(self):
        another_level = copy.deepcopy(self.energy_level)
        self.assertEqual(self.energy_level, another_level)

    def test_energy_level_string_init(self):
        self.assertEqual("1s+2(0)0 2s+2(0)0 2p-2(0)0 2p+3(3)3 3s+1(1)4", str(self.energy_level))

    def test_transition_get_weighted_oscillator_strength(self):
        self.assertAlmostEqual(0.5, self.transition.weighted_oscillator_strength, places=1)

    def test_energy_level_configuration(self):
        self.assertSequenceEqual([
            core.LevelTerm(shell=core.Shell.create_from_string("1s+2(0)"), j2=0),
            core.LevelTerm(shell=core.Shell.create_from_string("2s+2(0)"), j2=0),
            core.LevelTerm(shell=core.Shell.create_from_string("2p-2(0)"), j2=0),
            core.LevelTerm(shell=core.Shell.create_from_string("2p+3(3)"), j2=3),
            core.LevelTerm(shell=core.Shell.create_from_string("3s+1(1)"), j2=4)

        ], self.energy_level.configuration)


if __name__ == '__main__':
    unittest.main()
