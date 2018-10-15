[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

# EFlowCalc - An open-source calculator of ecological stream flow characteristics in Python

EFlowCalc is an open-source calculator of ecological stream flow characteristics in Python. It is licensed under GNU GPL-3.0 (see [licence file](https://github.com/ThibHlln/eflowcalc/blob/master/LICENCE.md) provided). EFlowCalc currently gives access to 17 of the 171 different ecologically relevant stream flow characteristics provided by [EflowStats](https://github.com/USGS-R/EflowStats) developed by USGS. More characteristics will gradually be added to EFlowCalc. The main advantage of EFlowCalc is the vectorisation of all computations which makes for very efficient calculation of the stream flow characteristics using [numpy](https://github.com/numpy/numpy) (and hence C code under the hood).

## How to Install

The simplest way to install EFlowCalc is to use pip and a link to the GitHub repository:

	python -m pip install git+https://github.com/ThibHlln/eflowcalc.git

Alternatively, you can download the source code (*i.e.* the GitHub repository) and, from the downloaded directory itself, run the command:

    python setup.py install

## Dependencies

EFlowCalc requires the popular Python packages `numpy` and `pandas` to be installed on the Python implementation where `eflowcalc` is installed.

## Version History

* 0.0.1 [14 Oct 2018]: First version of eFlowCalc

## Acknowledgment

This tool was developed with the financial support of Ireland's Environmental Protection Agency (Grant Number 2014-W-LS-5).