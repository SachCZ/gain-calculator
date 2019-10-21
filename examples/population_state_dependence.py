"""
Script to generate relative population graphs showing the state variables dependence
"""
import os
import gain_calculator as gc
from matplotlib import pyplot as plt
import numpy as np

if __name__ == "__main__":

    # Setup
    gc.init()
    energy_level = gc.EnergyLevel("1s+2(0)0 2s+2(0)0 2p-1(1)1 2p+4(0)1 3p-1(1)0")
    atom = gc.Atom(
        symbol="Ge",
        config_groups=gc.ConfigGroups(base="1*2 2*8", max_n=6),
        data_folder=os.path.join(os.path.abspath(os.path.dirname(__file__)), "atomic_data")
    )

    # The main calculation
    temperatures = [650.0, 850.0, 1050.0, 1250.0, 1450.0, 1650.0, 1850.0]
    population_values = atom.get_combined_populations(
        energy_level=energy_level,
        temperatures=temperatures,
        electron_densities=np.logspace(20, 23)
    )

    # Plotting
    for temperature in temperatures:
        relevant = population_values[population_values["temperature"] == temperature]
        densities = relevant["electron_density"]
        populations = relevant["population"] / energy_level.degeneracy
        plt.semilogx(densities, populations, label=u"$T_e = {}$ eV".format(temperature))

    plt.xlabel(u"$n_e$ [cm$^{-3}$]")
    plt.ylabel(u"$\\frac{N_i}{N}$ [-]")
    plt.grid()
    plt.legend()
    plt.tight_layout()
    plt.savefig("images/populations_Ge_2p-3p-.png")
