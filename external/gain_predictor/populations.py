import coefficients as cf
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from mpl_toolkits.mplot3d import Axes3D


def get_relative_populations(coeffs, ne):
    n1 = 1.0  # Assumption
    n3 = coeffs.C13 * ne / (coeffs.A32 + coeffs.C31 * ne + coeffs.C32 * ne)
    n2 = (coeffs.C12 * ne + n3 * (coeffs.A32 + coeffs.C32 * ne)) / coeffs.A21

    total = n1 + n2 + n3
    n1 = n1 / total
    n2 = n2 / total
    n3 = n3 / total

    return Populations(n1, n2, n3)


class Populations:
    def __init__(self, n1, n2, n3):
        self.relative_N1 = n1
        self.relative_N2 = n2
        self.relative_N3 = n3


def get_relative_inversion(coeffs, electron_density, temperature, base_degeneracy, lower_degeneracy,
                           upper_degeneracy):

    coeffs.recalculate_c(temperature, base_degeneracy, upper_degeneracy)
    populations = get_relative_populations(coeffs, electron_density)
    inversion = populations.relative_N3 / upper_degeneracy - populations.relative_N2 / lower_degeneracy
    return inversion


get_relative_inversion = np.vectorize(get_relative_inversion)

if __name__ == '__main__':
    my_coeffs = cf.Coefficients('ne.lev', 'ne.tr', 'ne.ce', '2p+4(0)0', '2p-1(1)1.3s+1(1)2', '2p-1(1)1.3p+1(3)4')

    densities = np.linspace(1e18, 0.4e21, 200)
    temperatures = np.linspace(200, 1400, 15)
    temperatures, densities = np.meshgrid(temperatures, densities)

    inversions = get_relative_inversion(
        coeffs=my_coeffs,
        electron_density=densities,
        temperature=temperatures,
        base_degeneracy=1,
        lower_degeneracy=3,
        upper_degeneracy=5
    )
    positive_inversions = inversions.clip(0)

    fig = plt.figure()
    ax = fig.gca(projection='3d')

    # Plot the surface.
    surf = ax.plot_surface(temperatures, densities, positive_inversions, cmap=cm.coolwarm,
                           linewidth=0, antialiased=False)

    # Customize the z axis.
    ax.set_zlim(0, np.max(inversions))
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

    # Add a color bar which maps values to colors.
    fig.colorbar(surf, shrink=0.5, aspect=5)

    plt.show()
