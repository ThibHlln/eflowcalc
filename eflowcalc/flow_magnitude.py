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
import pandas as pd
from .tools import calc_bfi, rolling_window


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# AVERAGE FLOWS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# MA26 - Variability of March flow
def ma26(flows, datetimes, hydro_years, drainage_area):
    # calculations per hydrological year
    info = np.zeros((flows.shape[0], hydro_years.shape[0],), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        month_march = flows[:, mask][:, 152:183] if sum(mask) == 366 else flows[:, mask][:, 151:182]
        info[:, hy] = np.std(month_march, axis=-1) / np.mean(month_march, axis=-1)
    # calculations for entire time series
    sfc = np.mean(info[:, :], axis=-1) * 100

    return sfc


# MA41 - Mean annual daily flow
def ma41(flows, datetimes, hydro_years, drainage_area):
    # calculations per hydrological year
    info = np.zeros((flows.shape[0], hydro_years.shape[0],), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[:, hy] = np.mean(flows[:, mask], axis=-1) / drainage_area
    # calculations for entire time series
    sfc = np.mean(info[:, :], axis=-1)

    return sfc


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOW FLOWS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ML17 - Base flow ratio
def ml17(flows, datetimes, hydro_years, drainage_area):
    # calculations per hydrological year
    info = np.zeros((flows.shape[0], hydro_years.shape[0],), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[:, hy] = calc_bfi(flows[:, mask])
    # calculations for entire time series
    sfc = np.mean(info[:, :], axis=-1)

    return sfc


# ML20 - Base flow 3
def ml20(flows, datetimes, hydro_years, drainage_area):
    # calculations per hydrological year
    # # NOTHING
    # calculations for entire time series
    nb_blocks = flows[:, :].shape[1] // 5
    flows_5day_blocks = flows[:, :][:, 0:5 * nb_blocks].reshape(flows.shape[0], nb_blocks, 5)
    blocks_min = np.amin(flows_5day_blocks, axis=-1)
    rolling_blocks_min = np.amin(rolling_window(blocks_min, 3), axis=-1)

    base_flow = np.zeros(blocks_min.shape, dtype=np.float64)
    base_flow[:] = np.nan

    check = blocks_min[:, 1:-1] * 0.90 < rolling_blocks_min
    base_flow[:, 1:-1][check] = blocks_min[:, 1:-1][check]
    base_flow[:, 0] = blocks_min[:, 0]
    base_flow[:, -1] = blocks_min[:, -1]

    base_flow_ = pd.DataFrame(base_flow).interpolate(method='linear', axis=1).values

    sfc = np.sum(base_flow_ * 5, axis=-1) / np.sum(flows[:, :], axis=-1)

    return sfc


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# HIGH FLOWS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# MH10 - Maximum October flow
def mh10(flows, datetimes, hydro_years, drainage_area):
    # calculations per hydrological year
    info = np.zeros((flows.shape[0], hydro_years.shape[0],), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[:, hy] = np.amax(flows[:, mask][:, 0:31], axis=-1)
    # calculations for entire time series
    sfc = np.mean(info[:, :], axis=-1)

    return sfc
