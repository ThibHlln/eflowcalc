# -*- coding: utf-8 -*-

# This file is part of EFlowCalc: A Calculator of Ecological Stream Flow Characteristics
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
    mean = np.mean(flows, axis=0)
    log_mean = np.log10(mean)
    # calculations per hydrological year
    colwell = np.zeros((365, 11, flows.shape[1]), dtype=int)
    log_f = np.copy(flows)
    log_f[log_f == 0.0] = 0.01  # replace log10(0) by log10(0.01) if necessary
    log_f = np.log10(log_f, dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        mask_no_lpy = np.copy(mask)
        if np.sum(mask) == 366:
            replacement = np.ones((366,), dtype=bool)
            replacement[151] = False  # remove 29th of February
            mask_no_lpy[mask] = replacement
        colwell[:, 0, :] += (log_f[mask_no_lpy, :] < (0.10 * log_mean)) * 1
        colwell[:, 1, :] += \
            ((log_f[mask_no_lpy, :] >= (0.10 * log_mean)) & (log_f[mask_no_lpy, :] < (0.25 * log_mean))) * 1
        colwell[:, 2, :] += \
            ((log_f[mask_no_lpy, :] >= (0.25 * log_mean)) & (log_f[mask_no_lpy, :] < (0.50 * log_mean))) * 1
        colwell[:, 3, :] += \
            ((log_f[mask_no_lpy, :] >= (0.50 * log_mean)) & (log_f[mask_no_lpy, :] < (0.75 * log_mean))) * 1
        colwell[:, 4, :] += \
            ((log_f[mask_no_lpy, :] >= (0.75 * log_mean)) & (log_f[mask_no_lpy, :] < (1.00 * log_mean))) * 1
        colwell[:, 5, :] += \
            ((log_f[mask_no_lpy, :] >= (1.00 * log_mean)) & (log_f[mask_no_lpy, :] < (1.25 * log_mean))) * 1
        colwell[:, 6, :] += \
            ((log_f[mask_no_lpy, :] >= (1.25 * log_mean)) & (log_f[mask_no_lpy, :] < (1.50 * log_mean))) * 1
        colwell[:, 7, :] += \
            ((log_f[mask_no_lpy, :] >= (1.50 * log_mean)) & (log_f[mask_no_lpy, :] < (1.75 * log_mean))) * 1
        colwell[:, 8, :] += \
            ((log_f[mask_no_lpy, :] >= (1.75 * log_mean)) & (log_f[mask_no_lpy, :] < (2.00 * log_mean))) * 1
        colwell[:, 9, :] += \
            ((log_f[mask_no_lpy, :] >= (2.00 * log_mean)) & (log_f[mask_no_lpy, :] < (2.25 * log_mean))) * 1
        colwell[:, 10, :] += (log_f[mask_no_lpy, :] >= (2.25 * log_mean)) * 1
    # calculations for entire time series
    colwell_y = np.sum(colwell, axis=0)  # sum up values in each column (i.e. state)
    colwell_z = np.sum(colwell, axis=(0, 1))  # sum up val in each matrix
    colwell_yz = np.divide(colwell_y, colwell_z, dtype=np.float64)
    colwell_log_yz = np.log10(colwell_yz, where=(colwell_yz != 0))
    colwell_log_yz[colwell_yz == 0] = np.nan
    colwell_hy = \
        - np.nansum(np.multiply(colwell_yz, colwell_log_yz, dtype=np.float64), axis=0)
    sfc = 1 - (colwell_hy / np.log10(11))

    return sfc


# TA2 - Predictability by Colwell (1974)
def ta2(flows, datetimes, hydro_years, drainage_area):
    mean = np.mean(flows, axis=0)
    log_mean = np.log10(mean)
    # calculations per hydrological year
    colwell = np.zeros((365, 11, flows.shape[1]), dtype=int)
    log_f = np.copy(flows)
    log_f[log_f == 0.0] = 0.01  # replace log10(0) by log10(0.01) if necessary
    log_f = np.log10(log_f, dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        mask_no_lpy = np.copy(mask)
        if np.sum(mask) == 366:
            replacement = np.ones((366,), dtype=bool)
            replacement[151] = False  # remove 29th of February
            mask_no_lpy[mask] = replacement
        colwell[:, 0, :] += (log_f[mask_no_lpy, :] < (0.10 * log_mean)) * 1
        colwell[:, 1, :] += \
            ((log_f[mask_no_lpy, :] >= (0.10 * log_mean)) & (log_f[mask_no_lpy, :] < (0.25 * log_mean))) * 1
        colwell[:, 2, :] += \
            ((log_f[mask_no_lpy, :] >= (0.25 * log_mean)) & (log_f[mask_no_lpy, :] < (0.50 * log_mean))) * 1
        colwell[:, 3, :] += \
            ((log_f[mask_no_lpy, :] >= (0.50 * log_mean)) & (log_f[mask_no_lpy, :] < (0.75 * log_mean))) * 1
        colwell[:, 4, :] += \
            ((log_f[mask_no_lpy, :] >= (0.75 * log_mean)) & (log_f[mask_no_lpy, :] < (1.00 * log_mean))) * 1
        colwell[:, 5, :] += \
            ((log_f[mask_no_lpy, :] >= (1.00 * log_mean)) & (log_f[mask_no_lpy, :] < (1.25 * log_mean))) * 1
        colwell[:, 6, :] += \
            ((log_f[mask_no_lpy, :] >= (1.25 * log_mean)) & (log_f[mask_no_lpy, :] < (1.50 * log_mean))) * 1
        colwell[:, 7, :] += \
            ((log_f[mask_no_lpy, :] >= (1.50 * log_mean)) & (log_f[mask_no_lpy, :] < (1.75 * log_mean))) * 1
        colwell[:, 8, :] += \
            ((log_f[mask_no_lpy, :] >= (1.75 * log_mean)) & (log_f[mask_no_lpy, :] < (2.00 * log_mean))) * 1
        colwell[:, 9, :] += \
            ((log_f[mask_no_lpy, :] >= (2.00 * log_mean)) & (log_f[mask_no_lpy, :] < (2.25 * log_mean))) * 1
        colwell[:, 10, :] += (log_f[mask_no_lpy, :] >= (2.25 * log_mean)) * 1
    # calculations for entire time series
    colwell_x = np.sum(colwell, axis=1)  # sum up values in each row (i.e. time)
    colwell_z = np.sum(colwell, axis=(0, 1))  # sum up val in each matrix
    colwell_xz = np.divide(colwell_x, colwell_z, dtype=np.float64)
    colwell_log_xz = np.log10(colwell_xz, where=(colwell_xz != 0))
    colwell_log_xz[colwell_xz == 0] = np.nan
    colwell_hx = \
        - np.nansum(np.multiply(colwell_xz, colwell_log_xz, dtype=np.float64), axis=0)
    colwell_mz = np.divide(colwell, colwell_z, dtype=np.float64)
    colwell_log_mz = np.log10(colwell_mz, where=(colwell != 0))
    colwell_log_mz[colwell == 0] = np.nan
    colwell_hxy = \
        - np.nansum(np.multiply(colwell_mz, colwell_log_mz, dtype=np.float64), axis=(0, 1))
    sfc = (1 - ((colwell_hxy - colwell_hx) / np.log10(11))) * 100

    return sfc


# TA3 not available


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOW FLOWS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# TL1 - Timing of annual minimum flow
def tl1(flows, datetimes, hydro_years, drainage_area):
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1], 2), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        min_julian_day = pd.DatetimeIndex(datetimes[mask][np.argmin(flows[mask, :], axis=0)]).dayofyear
        min_julian_day = min_julian_day * 2.0 * math.pi / 365.25
        info[hy, :, 0] = np.cos(min_julian_day)
        info[hy, :, 1] = np.sin(min_julian_day)
    # calculations for entire time series
    x = np.mean(info[:, :, 0], axis=0)
    y = np.mean(info[:, :, 1], axis=0)
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


# TL2 - Variability in timing of annual minimum flow
def tl2(flows, datetimes, hydro_years, drainage_area):
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1], 2), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        min_julian_day = pd.DatetimeIndex(datetimes[mask][np.argmin(flows[mask, :], axis=0)]).dayofyear
        min_julian_day = min_julian_day * 2.0 * math.pi / 365.25
        info[hy, :, 0] = np.cos(min_julian_day)
        info[hy, :, 1] = np.sin(min_julian_day)
    # calculations for entire time series
    x = np.mean(info[:, :, 0], axis=0)
    y = np.mean(info[:, :, 1], axis=0)
    tl2_ = np.sqrt(2 * (1 - np.sqrt((x * x) + (y * y))))
    sfc = tl2_ * 180 / math.pi * 365.25 / 360.0

    return sfc


# TL3 to TL4 not available


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# HIGH FLOWS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# TH1 - Timing of annual maximum flow
def th1(flows, datetimes, hydro_years, drainage_area):
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1], 2), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        max_julian_day = pd.DatetimeIndex(datetimes[mask][np.argmax(flows[mask, :], axis=0)]).dayofyear
        max_julian_day = max_julian_day * 2.0 * math.pi / 365.25
        info[hy, :, 0] = np.cos(max_julian_day)
        info[hy, :, 1] = np.sin(max_julian_day)
    # calculations for entire time series
    x = np.mean(info[:, :, 0], axis=0)
    y = np.mean(info[:, :, 1], axis=0)
    th1_ = np.zeros((flows.shape[0],), dtype=np.float64)
    th1_[:] = np.nan
    th1_ = np.arctan(np.divide(y, x, where=(x != 0), dtype=np.float32)) * 180 / math.pi
    th1_[x < 0] = th1_[x < 0] + 180.0
    th1_[(x == 0) & (y > 0)] = 90
    th1_[(x == 0) & (y < 0)] = 270
    th1_[th1_ < 0] = th1_[th1_ < 0] + 360
    th1_ = th1_ * 365.25 / 360.0
    th1_[th1_ == 0] = 365.25
    sfc = np.rint(th1_)

    return sfc


# TH2 - Variability in timing of annual maximum flow
def th2(flows, datetimes, hydro_years, drainage_area):
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1], 2), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        max_julian_day = pd.DatetimeIndex(datetimes[mask][np.argmax(flows[mask, :], axis=0)]).dayofyear
        max_julian_day = max_julian_day * 2.0 * math.pi / 365.25
        info[hy, :, 0] = np.cos(max_julian_day)
        info[hy, :, 1] = np.sin(max_julian_day)
    # calculations for entire time series
    x = np.mean(info[:, :, 0], axis=0)
    y = np.mean(info[:, :, 1], axis=0)
    th2_ = np.sqrt(2 * (1 - np.sqrt((x * x) + (y * y))))
    sfc = th2_ * 180 / math.pi * 365.25 / 360.0

    return sfc


# TH3 not available
