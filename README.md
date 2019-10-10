# Gain Calculator

> Small package to make atomic physics lasing gain calculations easier

This Python 2 package is a wrapper around the Flexible Atomic Code - 
[FAC](https://github.com/flexible-atomic-code/fac). It aims to greatly simplify the
calculation of gain based on atomic data provided by
[FAC](https://github.com/flexible-atomic-code/fac) and given state variables. 

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

For more examples please refer to the examples section in documentation.
The whole documentation consisting of getting started, more usage examples and
API reference can be found [here](#)