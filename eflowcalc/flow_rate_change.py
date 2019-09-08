# -*- coding: utf-8 -*-

# This file is part of EFlowCalc: A Calculator of Ecological Streamflow Characteristics
# Copyright (C) 2019  Thibault Hallouin (1)
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
from .tools import count_reversals


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ALL FLOWS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# RA1 - Average rise rate
def ra1(flows, datetimes, hydro_years, drainage_area):
    # calculations for entire time series
    diffs = np.diff(flows, axis=0)
    rises = np.copy(diffs)
    rises[rises <= 0] = np.nan
    sfc = np.nanmean(rises, axis=0)

    return sfc


# RA2 - Variability in rise rate
def ra2(flows, datetimes, hydro_years, drainage_area):
    # calculations for entire time series
    diffs = np.diff(flows, axis=0)
    rises = np.copy(diffs)
    rises[rises <= 0] = np.nan
    sfc = np.nanstd(rises, ddof=1, axis=0) * 100 / np.nanmean(rises, axis=0)

    return sfc


# RA3 - Average fall rate
def ra3(flows, datetimes, hydro_years, drainage_area):
    # calculations for entire time series
    diffs = np.diff(flows, axis=0)
    falls = np.copy(diffs)
    falls[falls >= 0] = np.nan
    falls = np.abs(falls)
    sfc = np.nanmean(falls, axis=0)

    return sfc


# RA4 - Variability in fall rate
def ra4(flows, datetimes, hydro_years, drainage_area):
    # calculations for entire time series
    diffs = np.diff(flows, axis=0)
    falls = np.copy(diffs)
    falls[falls >= 0] = np.nan
    falls = np.abs(falls)
    sfc = np.nanstd(falls, ddof=1, axis=0) * 100 / np.nanmean(falls, axis=0)

    return sfc


# RA5 - Ratio of days with flow rise
def ra5(flows, datetimes, hydro_years, drainage_area):
    # calculations for entire time series
    diffs = np.diff(flows, axis=0)
    rises = np.copy(diffs)
    sfc = np.true_divide(np.sum(rises > 0, axis=0), flows.shape[0])

    return sfc


# RA6 - Rate of flow rise
def ra6(flows, datetimes, hydro_years, drainage_area):
    # calculations for entire time series
    cp_flows = np.copy(flows)
    cp_flows[cp_flows == 0.0] = 0.01  # replace 0 by 0.01 if necessary (to avoid log(0))
    diffs_log = np.diff(np.log(cp_flows), axis=0)
    diffs_log[diffs_log <= 0] = np.nan  # take rises only
    sfc = np.nanmedian(np.abs(diffs_log), axis=0)

    return sfc


# RA7 - Rate of flow recession
def ra7(flows, datetimes, hydro_years, drainage_area):
    # calculations for entire time series
    cp_flows = np.copy(flows)
    cp_flows[cp_flows == 0.0] = 0.01  # replace 0 by 0.01 if necessary (to avoid log(0))
    diffs_log = np.diff(np.log(cp_flows), axis=0)
    diffs_log[diffs_log >= 0] = np.nan  # take falls only
    sfc = np.nanmedian(np.abs(diffs_log), axis=0)

    return sfc


# RA8 - Average annual number of reversals
def ra8(flows, datetimes, hydro_years, drainage_area):
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[hy, :] = count_reversals(flows[mask, :])
    # calculations for entire time series
    sfc = np.mean(info[:, :], axis=0)

    return sfc


# RA9 - Variability in number of annual reversals
def ra9(flows, datetimes, hydro_years, drainage_area):
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[hy, :] = count_reversals(flows[mask, :])
    # calculations for entire time series
    sfc = np.std(info, ddof=1, axis=0) * 100 / np.mean(info, axis=0)

    return sfc
