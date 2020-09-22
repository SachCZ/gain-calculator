from pfac.crm import *
from matplotlib import pyplot as plt
import numpy as np

get_abundance = np.vectorize(lambda t, e, z: FracAbund(z, t)[e])

if __name__ == '__main__':
    Z = 29
    temperatures = np.logspace(0, 5, num=5000)

    plt.yscale('log')
    plt.xscale('log')

    electrons = [0, 5, 10, 16, 22, 26]
    for i in electrons:
        plt.plot(temperatures, get_abundance(temperatures, i, Z), label="Ge" + str(Z - i) + "+")

    plt.ylim([1e-15, 1])
    plt.legend()
    plt.show()
