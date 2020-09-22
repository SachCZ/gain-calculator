import itertools
import pickle
import numpy as np
from pfac.crm import *


class GainCalculator:
    def __init__(self, filename=None):
        self.get_gain = np.vectorize(self.get_gain)
        self.get_population_at = np.vectorize(self.get_population_at)
        if filename:
            f = open(filename, 'rb')
            self.data = pickle.load(f)
        else:
            self.data = None

    def init_by_calculation(self, transition, temperatures, densities, log=None, combine=itertools.product):
        generated_populations = transition.get_populations(
            temperatures,
            densities,
            combine=combine,
            log=log
        )
        self.data = {
            "upper": {
                "temperature": generated_populations["upper"]["temperature"],
                "electron_density": generated_populations["upper"]["electron_density"],
                "population": generated_populations["upper"]["population"],
            },
            "lower": {
                "temperature": generated_populations["lower"]["temperature"],
                "electron_density": generated_populations["lower"]["electron_density"],
                "population": generated_populations["lower"]["population"],
            },
            "oscillator_strength": transition.weighted_oscillator_strength,
            "transition_energy": transition.energy,
            "temperatures": temperatures,
            "densities": densities
        }

    def serialize(self, filename):
        f = open(filename, 'wb')
        pickle.dump(self.data, f, protocol=2)

    def get_population_at(self, level, density, temperature):
        mask = (self.data[level]["temperature"] < temperature + 1e-6) * (
                self.data[level]["temperature"] > temperature - 1e-6) * (
                       self.data[level]["electron_density"] <= density + 1e-6) * (
                       self.data[level]["electron_density"] >= density - 1e-6)
        return self.data[level]["population"][mask][0]

    def get_temperatures(self):
        return self.data["temperatures"]

    def get_densities(self):
        return self.data["densities"]

    def get_gain(
            self,
            electron_density,
            temperature,
            ion_temperature,
            ionizations,
            electron_number,
            proton_number
    ):
        relative_inversion = self.get_population_at("upper", electron_density, temperature) - self.get_population_at(
            "lower", electron_density, temperature)
        fractional_abundance = self.__get_abundance(temperature, electron_number, proton_number)
        ion_density = electron_density / ionizations * fractional_abundance
        doppler_fwhm = 0.6 * self.__get_doppler_width(ion_temperature, self.data["transition_energy"])
        lorenz_fwhm = self.__get_lorenz_width(temperature, electron_density)

        return self.__calculate_gain_coeff(
            relative_inversion,
            ion_density,
            doppler_fwhm,
            lorenz_fwhm,
            self.data["oscillator_strength"]
        )

    __constants = {
        'c': 2.99792458e10,
        'h': 6.626068760e-27,
        'e': 4.8032068e-10,
        'm_e': 9.1093897e-28,
        'k_b': 1.602115e-12,
        'm_u': 1.6605e-24
    }

    def __get_abundance(self, t, e, z):
        return FracAbund(z, t)[e]

    def __get_doppler_width(self, ion_temperature, transition_energy):
        nu = transition_energy * 1.6021773e-12 / self.__constants['h']
        return nu * np.sqrt(
            2 * self.__constants['k_b'] * ion_temperature / (
                    self.__constants['m_u'] * np.power(self.__constants['c'], 2)))

    def __get_lorenz_width(self, temperature, electron_density):
        return 4.6746988e-29 * np.divide(np.power(electron_density, 2), np.sqrt(temperature))  # TODO look into this

    def __calculate_gain_coeff(
            self,
            relative_inversion,
            ion_density,
            doppler_fwhm,
            lorenz_fwhm,
            weighted_oscillator_strength
    ):
        eps0 = 1.0 / (4.0 * np.pi)
        delta_n = relative_inversion * ion_density

        profile_function = 1 / (doppler_fwhm + lorenz_fwhm)

        return np.power(self.__constants['e'], 2) / (
                4.0 * self.__constants['c'] * eps0 * self.__constants[
            'm_e']) * weighted_oscillator_strength * delta_n * profile_function
