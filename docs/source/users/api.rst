.. _API:

API
===
This document is the full description of user API. It should serve as a reference for anyone using the GainCalculator.
Although there are some code examples, these are very limited and only serve to demonstrate concepts.
If you are new to GainCalculator you should consider reading :ref:`getting-started`. If you are looking for
more code examples please visit :ref:`usage-examples`.

The API consists of several classes listed below:

.. autofunction:: gain_calculator.init

.. autoclass:: gain_calculator.EnergyLevel

.. autoclass:: gain_calculator.Atom
    :members: get_populations, get_combined_populations

.. autoclass:: gain_calculator.ConfigGroups

.. autoclass:: gain_calculator.Transition
    :members: get_populations

.. autofunction:: gain_calculator.print_progress
