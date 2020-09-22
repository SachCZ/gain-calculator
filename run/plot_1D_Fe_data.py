# -*- coding: utf-8 -*-
import gain_calculator as gc
from matplotlib import pyplot as plt
import numpy as np
from external.oliva_gain import get_previous_gain
from external.oliva_gain import get_previous_inversion


def plot_data(filename):
    calculator = gc.GainCalculator(filename)
    # The main calculation
    densities = calculator.get_densities()
    temperatures = calculator.get_temperatures()
    _, ion_temperatures = np.genfromtxt("run/data/cell_ion_temperature.csv", delimiter=",", unpack=True)
    x, ionizations = np.genfromtxt("run/data/cell_ionization.csv", delimiter=",", unpack=True)

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
    ax1.plot(x, densities)
    ax1.set_title('Electron density')
    ax1.set_yscale("log")
    ax1.grid()
    ax1.set_xlabel(u"$x$ [μm]")
    ax1.set_ylabel(u"$n_e$ [cm$^{-3}$]")
    ax2.plot(x, temperatures)
    ax2.set_title('Electron temperature')
    ax2.grid()
    ax2.set_xlabel(u"$x$ [μm]")
    ax2.set_ylabel(u"$T_e$ [eV]")
    ax3.plot(x, ion_temperatures)
    ax3.set_title('Ion temperature')
    ax3.grid()
    ax3.set_xlabel(u"$x$ [μm]")
    ax3.set_ylabel(u"$T_i$ [eV]")
    ax4.plot(x, ionizations)
    ax4.set_title('Ionization')
    ax4.grid()
    ax4.set_xlabel(u"$x$ [μm]")
    ax4.set_ylabel(u"$Z$ [-]")
    plt.tight_layout()
    plt.savefig("run/images/hydro_result.eps")
    plt.clf()

    mask = tuple([x > 5])
    gains = calculator.get_gain(densities, temperatures, ion_temperatures, ionizations, 10, 26)
    #a = np.asarray(zip(x, gains))
    #np.savetxt("foo.csv", a, delimiter=",", )

    previous_gains = get_previous_gain()[mask]

    fig, (ax1, ax2) = plt.subplots(1, 2)
    ax1.plot(x, gains / max(gains), label="M-level")
    ax1.plot(x[mask], previous_gains / max(previous_gains), "--", label="3-level")
    ax1.legend()
    ax1.set_ylabel(u"$g_{r}$[-]")
    ax1.set_xlabel(u"$x$ [μm]")
    ax1.grid()

    ax2.plot(x, gains, label="M-level")
    ax2.plot(x[mask], previous_gains, "--", label="3-level")
    ax2.legend()
    ax2.set_ylabel(u"$g$ [cm$^{-1}$]")
    ax2.set_xlabel(u"$x$ [μm]")
    ax2.grid()
    plt.tight_layout()
    plt.savefig("run/images/absolute_gain.eps")
    plt.clf()

def plot_populations(filename):
    calculator = gc.GainCalculator(filename)
    # The main calculation
    densities = calculator.get_densities()
    temperatures = calculator.get_temperatures()
    _, ion_temperatures = np.genfromtxt("run/data/cell_ion_temperature.csv", delimiter=",", unpack=True)
    x, ionizations = np.genfromtxt("run/data/cell_ionization.csv", delimiter=",", unpack=True)

    mask = tuple([x > 5])
    inversion = calculator.get_population_at("upper", densities, temperatures) - calculator.get_population_at("lower", densities, temperatures)
    previous_inversion = get_previous_inversion()[mask]
    fig, (ax1, ax2) = plt.subplots(1, 2)
    ax1.plot(x, inversion / max(inversion), label="M-level")
    ax1.plot(x[mask], previous_inversion / max(previous_inversion), "--", label="3-level")
    ax1.set_xlabel(u"$x$ [μm]")
    ax1.set_ylabel(u"${R}/{R_0}$ [-]")
    ax1.grid()
    ax1.legend()
    ax2.plot(x, inversion, label="M-level")
    ax2.plot(x[mask], previous_inversion, "--", label="3-level")
    ax2.set_xlabel(u"$x$ [μm]")
    ax2.set_ylabel(u"$R$ [-]")
    ax2.grid()
    ax2.legend()
    plt.tight_layout()
    plt.savefig("run/images/population_compare.eps")

def main():
    my_file = "run/data/1D_Fe_data.npy"
    plot_data(my_file)
    plot_populations(my_file)


if __name__ == "__main__":
    main()
