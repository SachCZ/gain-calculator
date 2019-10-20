.. _getting-started:

Getting started
===============

To get started using GainCalculator library you have to properly install it. For the installation instructions,
please refer to readme of the `project on GitHub <https://github.com/SachCZ/gain-calculator>`_.

After you have GainCalculator successfully installed on your machine you can use it like this::

    import gain_calculator as gc
    atom = gc.Atom(
        symbol="Ge",
        config_groups=gc.ConfigGroups(base="1*2 2*8", max_n=3),
    )

    print atom.get_combined_populations(
        energy_level=gc.EnergyLevel("1s+2(0)0 2s+2(0)0 2p-2(0)0 2p+3(3)3 3s+1(1)4"),
        temperatures=900,
        electron_densities=1e20,
    )["population"][0]

This will print the relative electron population on the specified energy level in given atom. For more usage examples
see :ref:`usage-examples` but first it is recommended to read the following paragraphs about the used physics
models and code architecture.

Used conventions
----------------

In the example above there is specified an energy level. This specification is done using the jj-coupling naming
convention. This means that shells are given by principal quantum number :math:`n`, orbital quantum number :math:`l`,
the way the spin is coupled (+ or -), number of electrons in given shell and total angular momentum of the
shell itself and total angular momentum of the whole system up to this shell. To get better understanding of
specifically how is this represented in the code see :ref:`API`.

The RAY library
---------------

GainCalculator relies heavily on the `RAY library <https://github.com/ray-project/ray>`_. It takes care of the
parallel computation. This is the reason why ``init()`` must be called before using the GainCalculator as it
initializes the RAY library. The RAY library then spawns multiple Python processes during computation to harness
multiple computational cores.

The FAC library
---------------

GainCalculator is roughly speaking just a parallelized convenience wrapper around the
`Flexible atomic code <https://github.com/flexible-atomic-code/fac>`_. It wraps its functionality to provide the user
with simple interface intended only to be used for gain calculations. FAC can be used for many other types of
calculations and as such its interface is way more complicated than it could be for gain calculations.
GainCalculator hides all the other logic as it is unnecessary burden for the user.



