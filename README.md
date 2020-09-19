[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![PyPI Version](https://badge.fury.io/py/eflowcalc.svg)](https://pypi.python.org/pypi/eflowcalc)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.2566757.svg)](https://doi.org/10.5281/zenodo.2566757)

# EFlowCalc - An efficient calculator of ecological streamflow characteristics in Python

`eflowcalc` is an open-source `calculator` of streamflow characteristics
in Python. It is licensed under GNU GPL-3.0 (see licence file provided).
The package currently gives the Python scientific community access
to 159 ecologically relevant streamflow characteristics inventoried by 
[Olden and Poff (2003)](https://doi.org/10.1002/rra.700). A key strength 
of `eflowcalc` is the vectorisation of all calculations (using 
[numpy](https://github.com/numpy/numpy), and therefore C code in the
background) which makes for very efficient computation of the streamflow
characteristics. 

Please refer to the [online documentation](https://thibhlln.github.io/eflowcalc) 
for further details and learn how to use the package. 

If you are using `eflowcalc`, please consider citing the software as 
follows (click on the link to get the DOI of a specific version):
* Hallouin, T. (XXXX). EFlowCalc: Ecological Streamflow Characteristics 
  Calculator (Version X.X.X). Zenodo. https://doi.org/10.5281/zenodo.2566757

## Version History

* 0.0.3 [08 Sep 2019]: [General enhancements](https://github.com/ThibHlln/eflowcalc/releases/tag/v0.0.3)
* 0.0.2 [16 Feb 2019]: [Version with 142 additional SFCs (159 SFCs)](https://github.com/ThibHlln/eflowcalc/releases/tag/v0.0.2)
* 0.0.1 [26 Oct 2018]: [First version of EFlowCalc (17 SFCs)](https://github.com/ThibHlln/eflowcalc/releases/tag/v0.0.1)

## Acknowledgment

Early versions of this tool were developed with the financial support of 
Ireland’s Environmental Protection Agency (Grant Number 2014-W-LS-5).
