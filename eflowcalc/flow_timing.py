# -*- coding: utf-8 -*-

# This file is part of EFlowCalc: A Calculator for Ecological Stream Flow Characteristics
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
import math


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# AVERAGE FLOWS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# TA1 - Constancy by Colwell (1974)
def ta1(flows, datetimes, hydro_years, drainage_area):
    mean = np.mean(flows[:, :], axis=-1)
    log_mean = np.reshape(np.log10(mean), (flows[:, :].shape[0], 1))
    # calculations per hydrological year
    colwell = np.zeros((flows.shape[0], 365, 11), dtype=int)
    log_f = np.copy(flows[:, :])
    log_f[log_f == 0.0] = 0.01  # replace log10(0) by log10(0.01) if necessary
    log_f = np.log10(log_f, dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        mask_no_lpy = np.copy(mask)
        if np.sum(mask) == 366:
            replacement = np.ones((366,), dtype=bool)
            replacement[151] = False  # remove 29th of February
            mask_no_lpy[mask] = replacement
        colwell[:, :, 0] += (log_f[:, mask_no_lpy] < (0.10 * log_mean)) * 1
        colwell[:, :, 1] += \
            ((log_f[:, mask_no_lpy] >= (0.10 * log_mean)) & (log_f[:, mask_no_lpy] < (0.25 * log_mean))) * 1
        colwell[:, :, 2] += \
            ((log_f[:, mask_no_lpy] >= (0.25 * log_mean)) & (log_f[:, mask_no_lpy] < (0.50 * log_mean))) * 1
        colwell[:, :, 3] += \
            ((log_f[:, mask_no_lpy] >= (0.50 * log_mean)) & (log_f[:, mask_no_lpy] < (0.75 * log_mean))) * 1
        colwell[:, :, 4] += \
            ((log_f[:, mask_no_lpy] >= (0.75 * log_mean)) & (log_f[:, mask_no_lpy] < (1.00 * log_mean))) * 1
        colwell[:, :, 5] += \
            ((log_f[:, mask_no_lpy] >= (1.00 * log_mean)) & (log_f[:, mask_no_lpy] < (1.25 * log_mean))) * 1
        colwell[:, :, 6] += \
            ((log_f[:, mask_no_lpy] >= (1.25 * log_mean)) & (log_f[:, mask_no_lpy] < (1.50 * log_mean))) * 1
        colwell[:, :, 7] += \
            ((log_f[:, mask_no_lpy] >= (1.50 * log_mean)) & (log_f[:, mask_no_lpy] < (1.75 * log_mean))) * 1
        colwell[:, :, 8] += \
            ((log_f[:, mask_no_lpy] >= (1.75 * log_mean)) & (log_f[:, mask_no_lpy] < (2.00 * log_mean))) * 1
        colwell[:, :, 9] += \
            ((log_f[:, mask_no_lpy] >= (2.00 * log_mean)) & (log_f[:, mask_no_lpy] < (2.25 * log_mean))) * 1
        colwell[:, :, 10] += (log_f[:, mask_no_lpy] >= (2.25 * log_mean)) * 1
    # calculations for entire time series
    colwell_y = np.sum(colwell, axis=-2)  # sum up values in each column
    colwell_z = np.reshape(np.sum(colwell, axis=(-1, -2)), (flows.shape[0], 1))  # sum up val in each matrix
    colwell_yz = np.divide(colwell_y, colwell_z, dtype=np.float64)
    colwell_log_yz = np.log10(colwell_yz, where=(colwell_yz != 0))
    colwell_log_yz[colwell_yz == 0] = np.nan
    colwell_hy = \
        - np.nansum(np.multiply(colwell_yz, colwell_log_yz, dtype=np.float64), axis=-1)
    sfc = 1 - (colwell_hy / np.log10(11))

    return sfc


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOW FLOWS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# TL1 - Timing of annual minimum flow
def tl1(flows, datetimes, hydro_years, drainage_area):
    # calculations per hydrological year
    info = np.zeros((flows.shape[0], hydro_years.shape[0], 2), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        min_julian_day = pd.DatetimeIndex(datetimes[mask][np.argmin(flows[:, mask], axis=-1)]).dayofyear
        min_julian_day = min_julian_day * 2.0 * math.pi / 365.25
        info[:, hy, 0] = np.cos(min_julian_day)
        info[:, hy, 1] = np.sin(min_julian_day)
    # calculations for entire time series
    x = np.mean(info[:, :, 0], axis=-1)
    y = np.mean(info[:, :, 1], axis=-1)
    tl1_ = np.zeros((flows.shape[0],), dtype=np.float64)
    tl1_[:] = np.nan
    tl1_ = np.arctan(np.divide(y, x, where=(x != 0), dtype=np.float32)) * 180 / math.pi
    tl1_[x < 0] = tl1_[x < 0] + 180.0
    tl1_[(x == 0) & (y > 0)] = 90
    tl1_[(x == 0) & (y < 0)] = 270
    tl1_[tl1_ < 0] = tl1_[tl1_ < 0] + 360
    tl1_ = tl1_ * 365.25 / 360.0
    tl1_[tl1_ == 0] = 365.25
    sfc = np.rint(tl1_)

    return sfc

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# HIGH FLOWS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


