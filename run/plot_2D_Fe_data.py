import gain_calculator as gc
from matplotlib import pyplot as plt
import numpy as np


def plot_data(filename):
    calculator = gc.GainCalculator(filename)
    # The main calculation
    original_densities = calculator.get_densities()
    original_temperatures = calculator.get_temperatures()

    densities, temperatures = np.meshgrid(original_densities, original_temperatures)

    gains = calculator.get_gain(densities, temperatures, temperatures, 25, 10, 26)
    positive_gains = gains.clip(0)

    i, j = np.unravel_index(positive_gains.argmax(), positive_gains.shape)
    print "temperature", original_temperatures[i]
    print "density", original_densities[j]
    print "max gain", np.max(gains)

    plt.contourf(temperatures, densities, gains, cmap=plt.get_cmap('plasma'))
    plt.colorbar()
    plt.xlabel(u"$T$ [eV]")
    plt.ylabel(u"$n_e$ [cm$^{-3}$]")
    plt.yscale("log")

    plt.tight_layout()

    plt.savefig("run/images/gain2d.eps")
    plt.clf()


def plot_populations(filename):
    calculator = gc.GainCalculator(filename)
    # The main calculation
    original_densities = calculator.get_densities()
    original_temperatures = calculator.get_temperatures()

    densities, temperatures = np.meshgrid(original_densities, original_temperatures)

    gains = calculator.get_gain(densities, temperatures, temperatures, 25, 10, 26)
    i, j = np.unravel_index(gains.argmax(), gains.shape)
    temperature = original_temperatures[i]

    upper_population = [calculator.get_population_at("upper", density, temperature) for density in original_densities]
    lower_population = [calculator.get_population_at("lower", density, temperature) for density in original_densities]

    plt.plot(original_densities, np.asarray(lower_population) * 100, label=u"3s+1(1)2 (lower)")
    plt.plot(original_densities, np.asarray(upper_population) * 100, label=u"3p-1(1)0 (upper)")
    plt.grid()
    plt.xlabel(u"$T$ [eV]")
    plt.ylabel(u"$N$ [%]")
    plt.xscale("log")
    plt.legend()

    plt.savefig("run/images/population.eps")


def main():
    my_file = "run/data/2D_Fe_data.npy"
    plot_data(my_file)
    plot_populations(my_file)


if __name__ == "__main__":
    main()
