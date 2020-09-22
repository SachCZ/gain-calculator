import logging

import gain_calculator as gc
import os
import numpy as np

def create_data(filename, densities, temperatures):
    atom = gc.Atom(
        symbol="Fe",
        config_groups=gc.ConfigGroups(base="1*2 2*8", max_n=6),
        data_folder=os.path.join(os.path.abspath(os.path.dirname(__file__)), "atomic_data")
    )

    lower = gc.EnergyLevel("1s+2(0)0 2s+2(0)0 2p-1(1)1 2p+4(1)1 3s+1(1)2")
    upper = gc.EnergyLevel("1s+2(0)0 2s+2(0)0 2p-1(1)1 2p+4(6)1 3p-1(1)0")
    transition = gc.Transition(atom, lower, upper)

    calculator = gc.GainCalculator()
    calculator.init_by_calculation(
        transition,
        temperatures,
        densities,
        combine=zip,
        log=lambda current, total: gc.print_progress(current, total, "Generating populations:")
    )
    calculator.serialize(filename)


def main():
    gc.init(logging_handler=logging.StreamHandler())
    _, densities = np.genfromtxt("data/cell_density.csv", delimiter=",", unpack=True)
    _, ionizations = np.genfromtxt("data/cell_ionization.csv", delimiter=",", unpack=True)
    mass_unit = 1.6605e-24
    A = 55.845
    electron_densities = ionizations * densities / A / mass_unit
    _, temperatures = np.genfromtxt("data/cell_temperature.csv", delimiter=",", unpack=True)
    create_data("data/1D_Fe_data.npy", electron_densities, temperatures)


if __name__ == '__main__':
    main()
