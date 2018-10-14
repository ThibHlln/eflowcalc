# -*- coding: utf-8 -*-

# This file is part of HydroEval: A Calculator for Ecological Stream Flow Characteristics
# Copyright (C) 2018  Thibault Hallouin (1)
#
# (1) Dooge Centre for Water Resources Research, University College Dublin, Ireland
#
# EFlowCalc is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# EFlowCalc is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with EFlowCalc. If not, see <http://www.gnu.org/licenses/>.

import numpy as np
from .tools import count_events


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# AVERAGE FLOWS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOW FLOWS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# FL2 - Variability in low pulse count
def fl2(flows, datetimes, hydro_years, drainage_area):
    quantile25 = np.reshape(np.quantile(flows[:, :], .25, axis=-1), (flows[:, :].shape[0], 1))
    # calculations per hydrological year
    info = np.zeros((flows.shape[0], hydro_years.shape[0],), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[:, hy] = count_events(flows[:, mask], threshold=quantile25, typ='low')
    # calculations for entire time series
    sfc = np.std(info[:, :], axis=-1) * 100 / np.mean(info[:, :], axis=-1)

    return sfc


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# HIGH FLOWS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# FH6 - Frequency of moderate floods
def fh6(flows, datetimes, hydro_years, drainage_area):
    median = np.median(flows[:, :], axis=-1)
    median_rs = np.reshape(median, (median.size, 1))
    # calculations per hydrological year
    info = np.zeros((flows.shape[0], hydro_years.shape[0],), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[:, hy] = count_events(flows[:, mask], median_rs * 3, typ='high')
    # calculations for entire time series
    sfc = np.mean(info[:, :], axis=-1)

    return sfc


# FH7 - Frequency of large floods
def fh7(flows, datetimes, hydro_years, drainage_area):
    median = np.median(flows[:, :], axis=-1)
    median_rs = np.reshape(median, (median.size, 1))
    # calculations per hydrological year
    info = np.zeros((flows.shape[0], hydro_years.shape[0],), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[:, hy] = count_events(flows[:, mask], median_rs * 7, typ='high')
    # calculations for entire time series
    sfc = np.mean(info[:, :], axis=-1)

    return sfc


# FH9 - Flood frequency
def fh9(flows, datetimes, hydro_years, drainage_area):
    quantile25 = np.reshape(np.quantile(flows[:, :], .25, axis=-1), (flows[:, :].shape[0], 1))
    # calculations per hydrological year
    info = np.zeros((flows.shape[0], hydro_years.shape[0],), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[:, hy] = count_events(flows[:, mask], threshold=quantile25, typ='high')
    # calculations for entire time series
    sfc = np.mean(info[:, :], axis=-1)

    return sfc


