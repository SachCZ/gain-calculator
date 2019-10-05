import gain_calculator.core as core
import numpy as np
from matplotlib import pyplot as plt

if __name__ == "__main__":
    krypton = core.Atom(
        symbol="Kr",
        base_level=core.EnergyLevel.create_from_string("1s+2(0)0 2s+2(0)0 2p-2(0)0 2p+4(0)0"),
    )

    energy_levels = [
        #core.EnergyLevel.create_from_string("1s+2(0)0 2s+2(0)0 2p-1(1)1 2p+4(0)1 3s+1(1)2"),
        core.EnergyLevel.create_from_string("1s+2(0)0 2s+2(0)0 2p-2(0)0 2p+3(3)3 3p+1(3)6"),
        core.EnergyLevel.create_from_string("1s+2(0)0 2s+2(0)0 2p-2(0)0 2p+3(3)3 3p+1(3)2"),
        core.EnergyLevel.create_from_string("1s+2(0)0 2s+2(0)0 2p-2(0)0 2p+3(3)3 3p+1(3)4"),
        core.EnergyLevel.create_from_string("1s+2(0)0 2s+2(0)0 2p-2(0)0 2p+3(3)3 3p+1(3)0"),
        core.EnergyLevel.create_from_string("1s+2(0)0 2s+2(0)0 2p-1(1)1 2p+4(0)1 3p-1(1)2"),
        core.EnergyLevel.create_from_string("1s+2(0)0 2s+2(0)0 2p-1(1)1 2p+4(0)1 3p+1(3)2"),
        core.EnergyLevel.create_from_string("1s+2(0)0 2s+2(0)0 2p-1(1)1 2p+4(0)1 3p+1(3)4"),
    ]
    for energy_level in energy_levels:
        def get_population(density):
            return krypton.get_population(
                energy_level=energy_level,
                principal_number_threshold=8,
                temperature=900,
                density=density
            ) / float(energy_level.get_degeneracy())

        x, y = list(zip(*[(den, get_population(den)) for den in np.logspace(16, 22)]))
        plt.loglog(x, y)

    plt.show()