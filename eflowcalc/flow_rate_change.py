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


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# AVERAGE FLOWS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# RA2 - Variability in rise rate
def ra2(flows, datetimes, hydro_years, drainage_area):
    # calculations per hydrological year
    # # NOTHING
    # calculations for entire time series
    diffs = np.diff(flows[:, :], axis=-1)
    rises = np.copy(diffs)
    rises[rises < 0] = np.nan
    sfc = np.nanstd(rises, axis=-1) * 100 / np.nanmean(rises, axis=-1)

    return sfc


# RA7 - Rate of flow recession
def ra7(flows, datetimes, hydro_years, drainage_area):
    # calculations per hydrological year
    # # NOTHING
    # calculations for entire time series
    cp_flows = np.copy(flows[:, :])
    cp_flows[cp_flows == 0.0] = 0.01  # replace 0 by 0.01 if necessary (to avoid log(0))
    diffs_log = np.diff(np.log(cp_flows), axis=-1)
    diffs_log[diffs_log > 0] = np.nan  # take falls only
    sfc = np.nanmedian(np.abs(diffs_log), axis=-1)

    return sfc

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOW FLOWS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# HIGH FLOWS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

