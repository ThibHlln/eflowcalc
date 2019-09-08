[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![PyPI Version](https://badge.fury.io/py/eflowcalc.svg)](https://pypi.python.org/pypi/eflowcalc)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.2566757.svg)](https://doi.org/10.5281/zenodo.2566757)

# EFlowCalc - An efficient calculator of ecological streamflow characteristics in Python

EFlowCalc is an open-source calculator of ecological streamflow characteristics in Python. It is licensed under GNU GPL-3.0 (see [licence file](https://github.com/ThibHlln/eflowcalc/blob/master/LICENCE.md) provided). EFlowCalc currently gives access to 159 ecologically relevant streamflow characteristics inventoried by [Olden and Poff (2003)](https://doi.org/10.1002/rra.700). More characteristics are gradually added to EFlowCalc. A key strength of EFlowCalc is the vectorisation of all calculations (using [numpy](https://github.com/numpy/numpy), and therefore C code in the background) which makes for very efficient computation of the streamflow characteristics.

## How to Install

EFlowCalc is available on PyPI, so you can simply use pip and the name of the package:

    python -m pip install eflowcalc

You can also use pip and a link to the GitHub repository directly:

	python -m pip install git+https://github.com/ThibHlln/eflowcalc.git

Alternatively, you can download the source code (*i.e.* the GitHub repository) and, from the downloaded directory itself, run the command:

    python setup.py install

## How to Use

A tutorial in the form of a [Jupyter notebook](https://github.com/ThibHlln/eflowcalc/blob/master/examples/api_usage_example.ipynb) is available to get started with the usage of EFlowCalc's API. The input file required for the tutorial is provided in the `examples/` folder.

## How to Cite

If you are using EFlowCalc, please consider citing the software as follows (click on the link to get the DOI of a specific version):
* Hallouin, T. (XXXX). EFlowCalc: Ecological Streamflow Characteristics Calculator (Version X.X.X). Zenodo. https://doi.org/10.5281/zenodo.2566757

## Streamflow Characteristics Available

The streamflow characteristics currently available in EFlowCalc are as follows:
* Magnitude of flow events
    * Average flow events: MA1 to MA45
    * Low flow events: ML1 to ML22
    * High flow events: MH1 to MH23
* Frequency of flow events
    * Low flow events: FL1 to FL3
    * High flow events: FH1 to FH10
* Duration of flow events
    * Low flow events: DL1 to DL20
    * High flow events: DH1 to DH21
* Timing of flow events
    * Average flow events: TA1 and TA2
    * Low flow events: TL1 and TL2
    * High flow events: TH1 and TH2
* Rate of change in flow events
    * Average flow events: RA1 to RA9
    
These streamflow characteristics are amongst the 171 hydrological indices inventoried by [Olden and Poff (2003)](https://doi.org/10.1002/rra.700). The computations implemented in EFlowCalc are partially inspired by the work of [Henriksen et al. (2006)](https://doi.org/10.3133/ofr20061093), however EFlowCalc is neither endorsed by these authors nor by the U.S. Geological Survey.

## Dependencies

EFlowCalc requires the Python packages `numpy` and `pandas` to be installed on the Python interpreter where `eflowcalc` is installed.

## Version History

* 0.0.3 [08 Sep 2019]: [General enhancements](https://github.com/ThibHlln/eflowcalc/releases/tag/v0.0.3)
* 0.0.2 [16 Feb 2019]: [Version with 142 additional SFCs (159 SFCs)](https://github.com/ThibHlln/eflowcalc/releases/tag/v0.0.2)
* 0.0.1 [26 Oct 2018]: [First version of EFlowCalc (17 SFCs)](https://github.com/ThibHlln/eflowcalc/releases/tag/v0.0.1)

## Acknowledgment

This tool was developed with the financial support of Ireland's Environmental Protection Agency (Grant Number 2014-W-LS-5).