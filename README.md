[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

# eFlowCalc - An open-source calculator of ecological stream flow characteristics in Python

eFlowCalc is an open-source calculator of ecological stream flow characteristics in Python. It is licensed under GNU GPL-3.0 (see [licence file](LICENCE.md) provided). eFlowCalc currently gives access to eighteen of the 171 different ecologically relevant stream flow characteristics provided by EflowStats developed by USGS. More characteristics will be gradually added to eFlowCalc. The advantage of eFlowCalc is the vectorisation of all computations which makes for very efficient calculation of the stream flow characteristics using numpy (and hence C code under the hood).

## How to Install

The simplest way to install eFlowCalc is to use pip and a link to the GitHub repository:

	python -m pip install git+https://github.com/ThibHlln/smartpy.git

Alternatively, you can download the source code (*i.e.* the GitHub repository) and, from the downloaded directory itself, run the command:

    python setup.py install

## Dependencies

eFlowCalc requires the popular Python packages `numpy` and `scipy` to be installed on the Python implementation where `eflowcalc` is installed.

## Version History

* 0.1.0 [14 Oct 2018]: First version of eFlowCalc

## Acknowledgment

This tool was developed with the financial support of Ireland's Environmental Protection Agency (Grant Number 2014-W-LS-5).