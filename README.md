# Gain Calculator

> Small package to make atomic physics lasing gain calculations easier

This Python 2 package is a wrapper around the Flexible Atomic Code - 
[FAC](https://github.com/flexible-atomic-code/fac). It aims to greatly simplify the
calculation of gain based state variables and atomic data provided by
[FAC](https://github.com/flexible-atomic-code/fac). 

## Dependencies
A dependency that is not installed during the installation of this package is the
[FAC](https://github.com/flexible-atomic-code/fac). Before using this package please
install the Flexible Atomic Code with Python interface pfac. 

You can check that pfac is installed successfully using those commands:
```
python -c "import pfac"
echo $?
# should print 0
```


## Install
To install the package simply run:
```bash
pip install gain_calculator
```

## Usage

The usage of this package is demonstrated with following example:

```python
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
```

To get started or to see more examples please refer to the [documentation]().
