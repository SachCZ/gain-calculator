import numpy as np
import coefficients as cf
import populations as pop
import frac_abundance as fa
from matplotlib import pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from mpl_toolkits.mplot3d import Axes3D

constants = {
    'c': 2.99792458e10,
    'h': 6.626068760e-27,
    'e': 4.8032068e-10,
    'm_e': 9.1093897e-28,
    'k_b': 1.602115e-12,
    'm_u': 1.6605e-24
}


def doppler_width(ion_temperature, coeffs):
    nu = coeffs.transition_energy * 1.6021773e-12 / constants['h']
    return nu * np.sqrt(2 * constants['k_b'] * ion_temperature / (constants['m_u'] * np.power(constants['c'], 2)))


def lorenz_width(temperature, electron_density):
    return 4.6746988e-29 * np.divide(np.power(electron_density, 2), np.sqrt(temperature))  # TODO look into this


def calculate_gain(electron_density, temperature, ion_temperature, ionization, electron_number, proton_number,
                   levels_filename, structure_filename, excitation_filename, base_level, lower_level,
                   upper_level, base_degeneracy, lower_degeneracy, upper_degeneracy):
    coeffs = cf.Coefficients(levels_filename, structure_filename, excitation_filename, base_level, lower_level,
                             upper_level)
    coeffs.recalculate_c(temperature, base_degeneracy, upper_degeneracy)
    relative_inversion = pop.get_relative_inversion(
        coeffs=coeffs,
        electron_density=electron_density,
        temperature=temperature,
        base_degeneracy=base_degeneracy,
        lower_degeneracy=lower_degeneracy,
        upper_degeneracy=upper_degeneracy
    )

    fractional_abundance = fa.get_abundance(temperature, electron_number, proton_number)
    ion_density = electron_density / ionization * fractional_abundance
    doppler_fwhm = 0.6 * doppler_width(ion_temperature, coeffs)
    lorenz_fwhm = lorenz_width(temperature, electron_density)

    return __calculate_gain_coeff(relative_inversion, ion_density, doppler_fwhm, lorenz_fwhm, coeffs)


calculate_gain = np.vectorize(calculate_gain)


def calculate_inversion(electron_density, temperature, ion_temperature, ionization, electron_number, proton_number,
                        levels_filename, structure_filename, excitation_filename, base_level, lower_level,
                        upper_level, base_degeneracy, lower_degeneracy, upper_degeneracy):
    coeffs = cf.Coefficients(levels_filename, structure_filename, excitation_filename, base_level, lower_level,
                             upper_level)
    coeffs.recalculate_c(temperature, base_degeneracy, upper_degeneracy)
    return pop.get_relative_inversion(
        coeffs=coeffs,
        electron_density=electron_density,
        temperature=temperature,
        base_degeneracy=base_degeneracy,
        lower_degeneracy=lower_degeneracy,
        upper_degeneracy=upper_degeneracy
    )


calculate_inversion = np.vectorize(calculate_inversion)


def __calculate_gain_coeff(relative_inversion, ion_density, doppler_fwhm, lorenz_fwhm, coeffs):
    eps0 = 1.0 / (4.0 * np.pi)
    delta_n = relative_inversion * ion_density

    profile_function = 1.0 / (doppler_fwhm + lorenz_fwhm)

    return np.power(constants['e'], 2) / (
            4.0 * constants['c'] * eps0 * constants['m_e']) * coeffs.gf * delta_n * profile_function


if __name__ == '__main__':
    my_coeffs = cf.Coefficients('ne.lev', 'ne.tr', 'ne.ce', '2p+4(0)0', '2p-1(1)1.3s+1(1)2', '2p-1(1)1.3p+1(3)4')
    my_coeffs.recalculate_c(900, 1, 5)

    my_coeffs.gf = 0.521

    print "Holden:", __calculate_gain_coeff(9.7e-5, 5e20 / 22.0, 0.6 * 1.82e12, 3.58e11, my_coeffs)

    print "LW Holden:", lorenz_width(900, 5e20)

    densities = np.linspace(1e18, 0.3e21, 50)
    temperatures = np.linspace(300, 850, 50)
    temperatures, densities = np.meshgrid(temperatures, densities)

    gains = calculate_gain(
        electron_density=densities,
        temperature=temperatures,
        ion_temperature=temperatures,
        ionization=16,
        electron_number=10,  # neon-like
        proton_number=30,  # almost germanium, higher is impossible
        levels_filename='ne.lev',
        structure_filename='ne.tr',
        excitation_filename='ne.ce',
        base_level='2p+4(0)0',
        lower_level='2p-1(1)1.3s+1(1)2',
        upper_level='2p-1(1)1.3p+1(3)4',
        base_degeneracy=1,
        lower_degeneracy=3,
        upper_degeneracy=5
    )
    positive_gains = gains.clip(0)

    fig = plt.figure()
    ax = fig.gca(projection='3d')

    # Plot the surface.
    surf = ax.plot_surface(temperatures, densities, positive_gains, cmap=cm.coolwarm,
                           linewidth=0, antialiased=False)

    # Customize the z axis.
    ax.set_zlim(0, np.max(positive_gains))
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

    # Add a color bar which maps values to colors.
    fig.colorbar(surf, shrink=0.5, aspect=5)

    plt.show()
