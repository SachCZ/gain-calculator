"""
Script to generate relative population graphs showing the state variables dependence
"""
import multiprocessing

import gain_calculator as gc
from matplotlib import pyplot as plt
import numpy as np

if __name__ == "__main__":
    energy_level = gc.EnergyLevel("1s+2(0)0 2s+2(0)0 2p-1(1)1 2p+4(0)1 3s+1(1)2")
    atom = gc.Atom(
        symbol="Ge",
        config_groups=gc.ConfigGroups("1*2 2*8", 6)
    )
    temperatures = [650, 850, 1050, 1250, 1450, 1650, 1850]
    population_values = atom.get_populations(
        energy_level=energy_level,
        temperatures=temperatures,
        electron_densities=np.logspace(20, 23)
    )

    for temperature in temperatures:
        relevant = [element for element in population_values if element["temperature"] == temperature]
        densities = [element["electron_density"] for element in relevant]
        populations = np.asarray([element["population"] for element in relevant]) / energy_level.degeneracy
        plt.semilogx(densities, populations, label=u"$T_e = {}$ eV".format(temperature))

    plt.xlabel(u"$n_e$ [cm$^{-3}$]")
    plt.ylabel(u"$\\frac{N_i}{N}$ [-]")
    plt.grid()
    plt.legend()

    plt.savefig("populations.png")
