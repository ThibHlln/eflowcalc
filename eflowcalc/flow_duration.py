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
from .tools import rolling_window, calc_events_avg_duration


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# AVERAGE FLOWS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOW FLOWS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# DL9 - Variability in annual minimum of 30-day average flow
def dl9(flows, datetimes, hydro_years, drainage_area):
    roll_30 = np.mean(rolling_window(flows[:, :], 30), axis=-1)
    # calculations per hydrological year
    info = np.zeros((flows.shape[0], hydro_years.shape[0],), dtype=np.float64)
    i = 0
    for hy, mask in enumerate(hydro_years):
        if hy == 0:  # i.e. first year in period
            info[:, hy] = np.amin(roll_30[:, 0:(np.sum(mask) - 14)], axis=-1)
            i += np.sum(mask) - 14
        elif hy == (np.sum(mask) - 1):  # i.e. last year in period
            info[:, hy] = np.amin(roll_30[:, i:(i + np.sum(mask) - 15)], axis=-1)
            i += np.sum(mask)
        else:
            info[:, hy] = np.amin(roll_30[:, i:(i + np.sum(mask))], axis=-1)
            i += np.sum(mask)
    # calculations for entire time series
    sfc = np.std(info[:, :], axis=-1) * 100 / np.mean(info[:, :], axis=-1)

    return sfc


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# HIGH FLOWS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# DH4 - Annual maximum of 30-day average flow
def dh4(flows, datetimes, hydro_years, drainage_area):
    roll_30 = np.mean(rolling_window(flows[:, :], 30), axis=-1)
    # calculations per hydrological year
    info = np.zeros((flows.shape[0], hydro_years.shape[0],), dtype=np.float64)
    i = 0
    for hy, mask in enumerate(hydro_years):
        if hy == 0:  # i.e. first year in period
            info[:, hy] = np.amax(roll_30[:, 0:(np.sum(mask) - 14)], axis=-1)
            i += np.sum(mask) - 14
        elif hy == (np.sum(mask) - 1):  # i.e. last year in period
            info[:, hy] = np.amax(roll_30[:, i:(i + np.sum(mask) - 15)], axis=-1)
            i += np.sum(mask)
        else:
            info[:, hy] = np.amax(roll_30[:, i:(i + np.sum(mask))], axis=-1)
            i += np.sum(mask)
    # calculations for entire time series
    sfc = np.mean(info[:, :], axis=-1)

    return sfc


# DH13 - Annual maximum of 30-day average flow normalised by median flow (of the whole record)
def dh13(flows, datetimes, hydro_years, drainage_area):
    median = np.median(flows[:, :], axis=-1)
    roll_30 = np.mean(rolling_window(flows[:, :], 30), axis=-1)
    # calculations per hydrological year
    info = np.zeros((flows.shape[0], hydro_years.shape[0],), dtype=np.float64)
    i = 0
    for hy, mask in enumerate(hydro_years):
        if hy == 0:  # i.e. first year in period
            info[:, hy] = np.amax(roll_30[:, 0:(np.sum(mask) - 14)], axis=-1)
            i += np.sum(mask) - 14
        elif hy == (np.sum(mask) - 1):  # i.e. last year in period
            info[:, hy] = np.amax(roll_30[:, i:(i + np.sum(mask) - 15)], axis=-1)
            i += np.sum(mask)
        else:
            info[:, hy] = np.amax(roll_30[:, i:(i + np.sum(mask))], axis=-1)
            i += np.sum(mask)
    # calculations for entire time series
    sfc = np.mean(info[:, :], axis=-1) / median

    return sfc


# DH16 - Variability in high-flow pulse count
def dh16(flows, datetimes, hydro_years, drainage_area):
    quantile75 = np.reshape(np.quantile(flows[:, :], .75, axis=-1), (flows[:, :].shape[0], 1))
    # calculations per hydrological year
    info = np.zeros((flows.shape[0], hydro_years.shape[0],), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[:, hy] = calc_events_avg_duration(flows[:, mask], quantile75, typ='high')
    # calculations for entire time series
    sfc = np.std(info[:, :], axis=-1) * 100 / np.mean(info[:, :], axis=-1)

    return sfc
