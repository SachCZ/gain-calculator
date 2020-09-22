#!/usr/bin/env python2.7

import numpy as np
import fac_parser as fp
from scipy.interpolate import CubicSpline


def calc_effective_coll_strength(collision_strength_data, temperature):
    transition_energy = collision_strength_data.transition_energy
    u = collision_strength_data['energy'] / temperature
    integrand = collision_strength_data['collision_strength'] * np.exp(-u)
    cubic_spline = CubicSpline(u, integrand)

    return cubic_spline.integrate(transition_energy / temperature, max(u))


def calc_excitation_coeff_c(collision_strength_data, temperature, lower_degeneration):
    effective = calc_effective_coll_strength(collision_strength_data, temperature)
    return 8.010e-8 / lower_degeneration / np.sqrt(temperature) * effective


def calc_deexcitation_coeff_c(collision_strength_data,
                              temperature,
                              upper_degeneration):
    excitation_coeff_c = calc_excitation_coeff_c(collision_strength_data, temperature, 1)
    return excitation_coeff_c / upper_degeneration * np.exp(collision_strength_data.transition_energy / temperature)


class Coefficients:
    def __init__(self, levels_filename, structure_filename, collision_strengths_filename, base_level, lower_level,
                 upper_level):
        self.__levels = fp.LevelsData(levels_filename)
        self.__structure_base_lower = fp.StructureData(structure_filename, self.__levels, base_level, lower_level)
        self.__structure_lower_upper = fp.StructureData(structure_filename, self.__levels, lower_level, upper_level)
        self.__col_str_base_lower = fp.CollisionStrengthData(collision_strengths_filename,
                                                             self.__levels, base_level, lower_level)
        self.__col_str_base_upper = fp.CollisionStrengthData(collision_strengths_filename,
                                                             self.__levels, base_level, upper_level)
        self.__col_str_lower_upper = fp.CollisionStrengthData(collision_strengths_filename,
                                                              self.__levels, lower_level, upper_level)

        self.A32 = self.__structure_lower_upper['transition_rate']
        self.A21 = self.__structure_base_lower['transition_rate']
        self.C31 = 0
        self.C32 = 0
        self.C13 = 0
        self.C12 = 0
        self.gf = self.__structure_lower_upper['gf']
        self.transition_energy = self.__col_str_lower_upper.transition_energy

    def recalculate_c(self, temperature, base_degeneration, upper_degeneration):
        self.C31 = calc_deexcitation_coeff_c(self.__col_str_base_upper, temperature, upper_degeneration)
        self.C32 = calc_deexcitation_coeff_c(self.__col_str_lower_upper, temperature, upper_degeneration)
        self.C13 = calc_excitation_coeff_c(self.__col_str_base_upper, temperature, base_degeneration)
        self.C12 = calc_excitation_coeff_c(self.__col_str_base_lower, temperature, base_degeneration)

    def __str__(self):
        return ''.join((
            "A32 = {:.3e}\n".format(self.A32),
            "A21 = {:.3e}\n".format(self.A21),
            "C31 = {:.3e}\n".format(self.C31),
            "C32 = {:.3e}\n".format(self.C32),
            "C13 = {:.3e}\n".format(self.C13),
            "C12 = {:.3e}\n".format(self.C12),
            "gf = {:.3}".format(self.gf)
        ))


if __name__ == '__main__':
    coeffs = Coefficients('ne.lev', 'ne.tr', 'ne.ce', '2p+4(0)0', '2p-1(1)1.3s+1(1)2', '2p-1(1)1.3p+1(3)4')
    coeffs.recalculate_c(900, 1, 5)

    print(coeffs)
