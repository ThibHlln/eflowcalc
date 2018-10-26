[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![PyPI Version](https://badge.fury.io/py/eflowcalc.svg)](https://pypi.python.org/pypi/eflowcalc)

# EFlowCalc - An open-source calculator of ecological stream flow characteristics in Python

EFlowCalc is an open-source calculator of ecological stream flow characteristics in Python. It is licensed under GNU GPL-3.0 (see [licence file](https://github.com/ThibHlln/eflowcalc/blob/master/LICENCE.md) provided). EFlowCalc currently gives access to 17 of the 171 different ecologically relevant stream flow characteristics inventoried by [Olden and Poff (2003)](https://doi.org/10.1002/rra.700). More characteristics will gradually be added to EFlowCalc. A key strength of EFlowCalc is the vectorisation of all calculations (using [numpy](https://github.com/numpy/numpy), and therefore C code in the background) which makes for very efficient computation of the stream flow characteristics.

## How to Install

EFlowCalc is available on PyPI, so you can simply use pip and the name of the package:

    python -m pip install eflowcalc

You can also use pip and a link to the GitHub repository directly:

	python -m pip install git+https://github.com/ThibHlln/eflowcalc.git

Alternatively, you can download the source code (*i.e.* the GitHub repository) and, from the downloaded directory itself, run the command:

    python setup.py install

## How to Use

A tutorial in the form of a [Jupyter notebook](https://github.com/ThibHlln/eflowcalc/blob/master/examples/api_usage_example.ipynb) is available to get started with the usage of EFlowCalc's API. The input file required for the tutorial is provided in the `examples/` folder.

## Stream Flow Characteristics Available

The stream flow characteristics currently available in EFlowCalc are as follows:
* Magnitude of flow events
    * MA26 - Variability of March flow
    * MA41 - Mean annual daily flow
    * ML17 - Base flow ratio
    * ML20 - Base flow 3
    * MH10 - Maximum October flow
* Frequency of flow events
    * FL2 - Variability in low pulse count
    * FH6 - Frequency of moderate floods
    * FH7 - Frequency of large floods
    * FH9 - Flood frequency
* Duration of flow events
    * DL9 - Variability in annual minimum of 30-day average flow
    * DH4 - Annual maximum of 30-day average flow
    * DH13 - Annual maximum of 30-day average flow normalised by median flow
    * DH16 - Variability in high-flow pulse count
* Timing of flow events
    * TA1 - Constancy by [Colwell (1974)](https://doi.org/10.2307/1940366)
    * TL1 - Timing of annual minimum flow
* Rate of change in flow events
    * RA2 - Variability in rise rate
    * RA7 - Rate of flow recession
    
These stream flow characteristics are amongst the 171 hydrological indices inventoried by [Olden and Poff (2003)](https://doi.org/10.1002/rra.700). The computations implemented in EFlowCalc are partially inspired by the work of [Henriksen et al. (2006)](https://doi.org/10.3133/ofr20061093), however EFlowCalc is neither endorsed by these authors nor by the U.S. Geological Survey.

## Dependencies

EFlowCalc requires the popular Python packages `numpy` and `pandas` to be installed on the Python implementation where `eflowcalc` is installed.

## Version History

* 0.0.1 [26 Oct 2018]: First version of EFlowCalc

## Acknowledgment

This tool was developed with the financial support of Ireland's Environmental Protection Agency (Grant Number 2014-W-LS-5).