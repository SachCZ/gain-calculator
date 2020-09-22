#!/usr/bin/env python

import gain_predictor as gp
import ncutils as ncu
import numpy as np
from matplotlib import pyplot as plt
import sys
import os


def get_previous_gain():
    m_u = 1.6605e-24

    filename = "./run/data/variables.nc"
    positions = ncu.get_variable(filename, 'cell_position')
    densities = ncu.get_variable(filename, 'cell_density')
    ionizations = ncu.get_variable(filename, 'cell_ionization')
    temperatures = ncu.get_variable(filename, 'cell_temperature')
    ionTemperatures = ncu.get_variable(filename, 'cell_ionTemperature')
    nucleon_number = ncu.get_variable(filename, 'cell_nucleonNumber')

    electronDensities = np.divide(ionizations * densities, nucleon_number * m_u)

    gains = gp.calculate_gain(
        electron_density=electronDensities,
        temperature=temperatures,
        ion_temperature=ionTemperatures,
        ionization=ionizations,
        electron_number=10,  # Neon-like
        proton_number=26,  # Iron
        levels_filename=os.path.join(os.path.abspath(os.path.dirname(__file__)), 'neon_like_iron_levels'),
        structure_filename=os.path.join(os.path.abspath(os.path.dirname(__file__)), 'neon_like_iron_transitions'),
        excitation_filename=os.path.join(os.path.abspath(os.path.dirname(__file__)), 'neon_like_iron_excitation'),
        base_level='2p+4(0)0',
        lower_level='2p-1(1)1.3s+1(1)2',
        upper_level='2p-1(1)1.3p-1(1)0',
        base_degeneracy=1.0,
        lower_degeneracy=3.0,
        upper_degeneracy=1.0
    )
    return gains


def get_previous_inversion():
    m_u = 1.6605e-24

    filename = "./run/data/variables.nc"
    positions = ncu.get_variable(filename, 'cell_position')
    densities = ncu.get_variable(filename, 'cell_density')
    ionizations = ncu.get_variable(filename, 'cell_ionization')
    temperatures = ncu.get_variable(filename, 'cell_temperature')
    ionTemperatures = ncu.get_variable(filename, 'cell_ionTemperature')
    nucleon_number = ncu.get_variable(filename, 'cell_nucleonNumber')

    electronDensities = np.divide(ionizations * densities, nucleon_number * m_u)

    inversion = gp.calculate_inversion(
        electron_density=electronDensities,
        temperature=temperatures,
        ion_temperature=ionTemperatures,
        ionization=ionizations,
        electron_number=10,  # Neon-like
        proton_number=26,  # Iron
        levels_filename=os.path.join(os.path.abspath(os.path.dirname(__file__)), 'neon_like_iron_levels'),
        structure_filename=os.path.join(os.path.abspath(os.path.dirname(__file__)), 'neon_like_iron_transitions'),
        excitation_filename=os.path.join(os.path.abspath(os.path.dirname(__file__)), 'neon_like_iron_excitation'),
        base_level='2p+4(0)0',
        lower_level='2p-1(1)1.3s+1(1)2',
        upper_level='2p-1(1)1.3p-1(1)0',
        base_degeneracy=1.0,
        lower_degeneracy=3.0,
        upper_degeneracy=1.0
    )
    return inversion
