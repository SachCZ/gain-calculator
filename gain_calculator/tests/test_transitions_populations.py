import gain_calculator.core as core
import numpy as np
from matplotlib import pyplot as plt

if __name__ == "__main__":
    atom = core.Atom(
        symbol="Ge",
        base_level=core.EnergyLevel.create_from_string("1s+2(0)0 2s+2(0)0 2p-2(0)0 2p+4(0)0"),
    )
    energy_level = core.EnergyLevel.create_from_string("1s+2(0)0 2s+2(0)0 2p-1(1)1 2p+4(0)1 3s+1(1)2")

    def get_population(density, temperature):
        return atom.get_population(
            energy_level=energy_level,
            principal_number_threshold=6,
            temperature=temperature,
            density=density
        ) / energy_level.get_degeneracy()

    for temp in [650, 850, 1050, 1250, 1450, 1650, 1850]:
        x, y = list(zip(*[(den, get_population(den, temp)) for den in np.logspace(20, 23)]))
        plt.semilogx(x, y)

    plt.show()