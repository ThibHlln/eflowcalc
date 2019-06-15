# -*- coding: utf-8 -*-

# This file is part of EFlowCalc: A Calculator of Ecological Streamflow Characteristics
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
from .tools import rolling_window, calc_events_avg_duration


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOW FLOWS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# DL1 - Annual minimum average daily flow
def dl1(flows, datetimes, hydro_years, drainage_area):
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[hy, :] = np.amin(flows[mask, :], axis=0)
    # calculations for entire time series
    sfc = np.mean(info, axis=0)

    return sfc


# DL2 - Mean annual minimum of 3-day average flow
def dl2(flows, datetimes, hydro_years, drainage_area):
    roll_3 = np.mean(rolling_window(flows, 3), axis=1)
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    i = 0
    for hy, mask in enumerate(hydro_years):
        if hy == 0:  # i.e. first year in period
            info[hy, :] = np.amin(roll_3[0:(np.sum(mask) - 1), :], axis=0)
            i += np.sum(mask) - 1
        elif hy == (np.sum(mask) - 1):  # i.e. last year in period
            info[hy, :] = np.amin(roll_3[i:(i + np.sum(mask) - 2), :], axis=0)
            i += np.sum(mask)
        else:
            info[hy, :] = np.amin(roll_3[i:(i + np.sum(mask)), :], axis=0)
            i += np.sum(mask)
    # calculations for entire time series
    sfc = np.mean(info, axis=0)

    return sfc


# DL3 - Mean annual minimum of 7-day average flow
def dl3(flows, datetimes, hydro_years, drainage_area):
    roll_7 = np.mean(rolling_window(flows, 7), axis=1)
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    i = 0
    for hy, mask in enumerate(hydro_years):
        if hy == 0:  # i.e. first year in period
            info[hy, :] = np.amin(roll_7[0:(np.sum(mask) - 3), :], axis=0)
            i += np.sum(mask) - 3
        elif hy == (np.sum(mask) - 1):  # i.e. last year in period
            info[hy, :] = np.amin(roll_7[i:(i + np.sum(mask) - 4), :], axis=0)
            i += np.sum(mask)
        else:
            info[hy, :] = np.amin(roll_7[i:(i + np.sum(mask)), :], axis=0)
            i += np.sum(mask)
    # calculations for entire time series
    sfc = np.mean(info, axis=0)

    return sfc


# DL4 - Mean annual minimum of 30-day average flow
def dl4(flows, datetimes, hydro_years, drainage_area):
    roll_30 = np.mean(rolling_window(flows, 30), axis=1)
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    i = 0
    for hy, mask in enumerate(hydro_years):
        if hy == 0:  # i.e. first year in period
            info[hy, :] = np.amin(roll_30[0:(np.sum(mask) - 14), :], axis=0)
            i += np.sum(mask) - 14
        elif hy == (np.sum(mask) - 1):  # i.e. last year in period
            info[hy, :] = np.amin(roll_30[i:(i + np.sum(mask) - 15), :], axis=0)
            i += np.sum(mask)
        else:
            info[hy, :] = np.amin(roll_30[i:(i + np.sum(mask)), :], axis=0)
            i += np.sum(mask)
    # calculations for entire time series
    sfc = np.mean(info, axis=0)

    return sfc


# DL5 - Mean annual minimum of 90-day average flow
def dl5(flows, datetimes, hydro_years, drainage_area):
    roll_90 = np.mean(rolling_window(flows, 90), axis=1)
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    i = 0
    for hy, mask in enumerate(hydro_years):
        if hy == 0:  # i.e. first year in period
            info[hy, :] = np.amin(roll_90[0:(np.sum(mask) - 44), :], axis=0)
            i += np.sum(mask) - 44
        elif hy == (np.sum(mask) - 1):  # i.e. last year in period
            info[hy, :] = np.amin(roll_90[i:(i + np.sum(mask) - 45), :], axis=0)
            i += np.sum(mask)
        else:
            info[hy, :] = np.amin(roll_90[i:(i + np.sum(mask)), :], axis=0)
            i += np.sum(mask)
    # calculations for entire time series
    sfc = np.mean(info, axis=0)

    return sfc


# DL6 - Variability in annual minimum average daily flow
def dl6(flows, datetimes, hydro_years, drainage_area):
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[hy, :] = np.amin(flows[mask, :], axis=0)
    # calculations for entire time series
    sfc = np.std(info, ddof=1, axis=0) * 100 / np.mean(info, axis=0)

    return sfc


# DL7 - Variability in annual minimum of 3-day average flow
def dl7(flows, datetimes, hydro_years, drainage_area):
    roll_3 = np.mean(rolling_window(flows, 3), axis=1)
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    i = 0
    for hy, mask in enumerate(hydro_years):
        if hy == 0:  # i.e. first year in period
            info[hy, :] = np.amin(roll_3[0:(np.sum(mask) - 1), :], axis=0)
            i += np.sum(mask) - 1
        elif hy == (np.sum(mask) - 1):  # i.e. last year in period
            info[hy, :] = np.amin(roll_3[i:(i + np.sum(mask) - 2), :], axis=0)
            i += np.sum(mask)
        else:
            info[hy, :] = np.amin(roll_3[i:(i + np.sum(mask)), :], axis=0)
            i += np.sum(mask)
    # calculations for entire time series
    sfc = np.std(info, ddof=1, axis=0) * 100 / np.mean(info, axis=0)

    return sfc


# DL8 - Variability in annual minimum of 7-day average flow
def dl8(flows, datetimes, hydro_years, drainage_area):
    roll_7 = np.mean(rolling_window(flows, 7), axis=1)
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    i = 0
    for hy, mask in enumerate(hydro_years):
        if hy == 0:  # i.e. first year in period
            info[hy, :] = np.amin(roll_7[0:(np.sum(mask) - 3), :], axis=0)
            i += np.sum(mask) - 3
        elif hy == (np.sum(mask) - 1):  # i.e. last year in period
            info[hy, :] = np.amin(roll_7[i:(i + np.sum(mask) - 4), :], axis=0)
            i += np.sum(mask)
        else:
            info[hy, :] = np.amin(roll_7[i:(i + np.sum(mask)), :], axis=0)
            i += np.sum(mask)
    # calculations for entire time series
    sfc = np.std(info, ddof=1, axis=0) * 100 / np.mean(info, axis=0)

    return sfc


# DL9 - Variability in annual minimum of 30-day average flow
def dl9(flows, datetimes, hydro_years, drainage_area):
    roll_30 = np.mean(rolling_window(flows, 30), axis=1)
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    i = 0
    for hy, mask in enumerate(hydro_years):
        if hy == 0:  # i.e. first year in period
            info[hy, :] = np.amin(roll_30[0:(np.sum(mask) - 14), :], axis=0)
            i += np.sum(mask) - 14
        elif hy == (np.sum(mask) - 1):  # i.e. last year in period
            info[hy, :] = np.amin(roll_30[i:(i + np.sum(mask) - 15), :], axis=0)
            i += np.sum(mask)
        else:
            info[hy, :] = np.amin(roll_30[i:(i + np.sum(mask)), :], axis=0)
            i += np.sum(mask)
    # calculations for entire time series
    sfc = np.std(info, ddof=1, axis=0) * 100 / np.mean(info, axis=0)

    return sfc


# DL10 - Mean annual minimum of 90-day average flow
def dl10(flows, datetimes, hydro_years, drainage_area):
    roll_90 = np.mean(rolling_window(flows, 90), axis=1)
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    i = 0
    for hy, mask in enumerate(hydro_years):
        if hy == 0:  # i.e. first year in period
            info[hy, :] = np.amin(roll_90[0:(np.sum(mask) - 44), :], axis=0)
            i += np.sum(mask) - 14
        elif hy == (np.sum(mask) - 1):  # i.e. last year in period
            info[hy, :] = np.amin(roll_90[i:(i + np.sum(mask) - 45), :], axis=0)
            i += np.sum(mask)
        else:
            info[hy, :] = np.amin(roll_90[i:(i + np.sum(mask)), :], axis=0)
            i += np.sum(mask)
    # calculations for entire time series
    sfc = np.std(info, ddof=1, axis=0) * 100 / np.mean(info, axis=0)

    return sfc


# DL11 - Mean annual minimum average flow normalised by median flow (of the whole record)
def dl11(flows, datetimes, hydro_years, drainage_area):
    median = np.median(flows, axis=0)
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[hy, :] = np.amin(flows[mask, :], axis=0)
    # calculations for entire time series
    sfc = np.mean(info, axis=0) / median

    return sfc


# DL12 - Mean annual minimum of 7-day average flow normalised by median flow (of the whole record)
def dl12(flows, datetimes, hydro_years, drainage_area):
    median = np.median(flows, axis=0)
    roll_7 = np.mean(rolling_window(flows, 7), axis=1)
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    i = 0
    for hy, mask in enumerate(hydro_years):
        if hy == 0:  # i.e. first year in period
            info[hy, :] = np.amin(roll_7[0:(np.sum(mask) - 3), :], axis=0)
            i += np.sum(mask) - 3
        elif hy == (np.sum(mask) - 1):  # i.e. last year in period
            info[hy, :] = np.amin(roll_7[i:(i + np.sum(mask) - 4), :], axis=0)
            i += np.sum(mask)
        else:
            info[hy, :] = np.amin(roll_7[i:(i + np.sum(mask)), :], axis=0)
            i += np.sum(mask)
    # calculations for entire time series
    sfc = np.mean(info, axis=0) / median

    return sfc


# DL13 - Mean annual minimum of 30-day average flow normalised by median flow (of the whole record)
def dl13(flows, datetimes, hydro_years, drainage_area):
    median = np.median(flows, axis=0)
    roll_30 = np.mean(rolling_window(flows, 30), axis=1)
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    i = 0
    for hy, mask in enumerate(hydro_years):
        if hy == 0:  # i.e. first year in period
            info[hy, :] = np.amin(roll_30[0:(np.sum(mask) - 14), :], axis=0)
            i += np.sum(mask) - 14
        elif hy == (np.sum(mask) - 1):  # i.e. last year in period
            info[hy, :] = np.amin(roll_30[i:(i + np.sum(mask) - 15), :], axis=0)
            i += np.sum(mask)
        else:
            info[hy, :] = np.amin(roll_30[i:(i + np.sum(mask)), :], axis=0)
            i += np.sum(mask)
    # calculations for entire time series
    sfc = np.mean(info, axis=0) / median

    return sfc


# DL14 - Q75 exceedance value normalised by median flow (of the whole record)
def dl14(flows, datetimes, hydro_years, drainage_area):
    median = np.median(flows, axis=0)
    perc25 = np.percentile(flows, 25, axis=0)

    sfc = perc25 / median

    return sfc


# DL15 - Q90 value normalised by median flow (of the whole record)
def dl15(flows, datetimes, hydro_years, drainage_area):
    median = np.median(flows, axis=0)
    perc10 = np.percentile(flows, 10, axis=0)

    sfc = perc10 / median

    return sfc


# DL16 - Median low-flow pulse duration
def dl16(flows, datetimes, hydro_years, drainage_area):
    perc25 = np.percentile(flows, 25, axis=0)
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[hy, :] = calc_events_avg_duration(flows[mask, :], perc25, typ='low')
    # calculations for entire time series
    sfc = np.median(info, axis=0)

    return sfc


# DL17 - Variability in low-flow pulse duration
def dl17(flows, datetimes, hydro_years, drainage_area):
    perc25 = np.percentile(flows, 25, axis=0)
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[hy, :] = calc_events_avg_duration(flows[mask, :], perc25, typ='low')
    # calculations for entire time series
    sfc = np.std(info, ddof=1, axis=0) * 100 / np.mean(info, axis=0)

    return sfc


# DL18 - Average annual number of zero flow days
def dl18(flows, datetimes, hydro_years, drainage_area):
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[hy, :] = np.sum(flows[mask, :] == 0, axis=0)
    # calculations for entire time series
    sfc = np.mean(info, axis=0)

    return sfc


# DL19 - Variability in annual number of zero flow days
def dl19(flows, datetimes, hydro_years, drainage_area):
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[hy, :] = np.sum(flows[mask, :] == 0, axis=0)
    # calculations for entire time series
    mean_ = np.mean(info, axis=0)
    sfc = np.true_divide(np.std(info, axis=0) * 100, mean_, where=(mean_ != 0))
    sfc[mean_ == 0] = 0.0

    return sfc


# DL20 - Annual number of zero flow months
def dl20(flows, datetimes, hydro_years, drainage_area):
    # calculations per month for each year
    zero_flow = np.array(pd.DataFrame(flows != 0, index=datetimes).groupby(lambda x: (x.year, x.month)).sum())
    # calculations for entire time series
    sfc = np.sum(zero_flow == 0, axis=0)

    return sfc


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# HIGH FLOWS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# DH1 - Annual maximum average daily flow
def dh1(flows, datetimes, hydro_years, drainage_area):
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[hy, :] = np.amax(flows[mask, :], axis=0)
    # calculations for entire time series
    sfc = np.mean(info, axis=0)

    return sfc


# DH2 - Mean annual maximum of 3-day average flow
def dh2(flows, datetimes, hydro_years, drainage_area):
    roll_3 = np.mean(rolling_window(flows, 3), axis=1)
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    i = 0
    for hy, mask in enumerate(hydro_years):
        if hy == 0:  # i.e. first year in period
            info[hy, :] = np.amax(roll_3[0:(np.sum(mask) - 1), :], axis=0)
            i += np.sum(mask) - 1
        elif hy == (np.sum(mask) - 1):  # i.e. last year in period
            info[hy, :] = np.amax(roll_3[i:(i + np.sum(mask) - 2), :], axis=0)
            i += np.sum(mask)
        else:
            info[hy, :] = np.amax(roll_3[i:(i + np.sum(mask)), :], axis=0)
            i += np.sum(mask)
    # calculations for entire time series
    sfc = np.mean(info, axis=0)

    return sfc


# DH3 - Mean annual maximum of 7-day average flow
def dh3(flows, datetimes, hydro_years, drainage_area):
    roll_7 = np.mean(rolling_window(flows, 7), axis=1)
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    i = 0
    for hy, mask in enumerate(hydro_years):
        if hy == 0:  # i.e. first year in period
            info[hy, :] = np.amax(roll_7[0:(np.sum(mask) - 3), :], axis=0)
            i += np.sum(mask) - 3
        elif hy == (np.sum(mask) - 1):  # i.e. last year in period
            info[hy, :] = np.amax(roll_7[i:(i + np.sum(mask) - 4), :], axis=0)
            i += np.sum(mask)
        else:
            info[hy, :] = np.amax(roll_7[i:(i + np.sum(mask)), :], axis=0)
            i += np.sum(mask)
    # calculations for entire time series
    sfc = np.mean(info, axis=0)

    return sfc


# DH4 - Annual maximum of 30-day average flow
def dh4(flows, datetimes, hydro_years, drainage_area):
    roll_30 = np.mean(rolling_window(flows, 30), axis=1)
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    i = 0
    for hy, mask in enumerate(hydro_years):
        if hy == 0:  # i.e. first year in period
            info[hy, :] = np.amax(roll_30[0:(np.sum(mask) - 14), :], axis=0)
            i += np.sum(mask) - 14
        elif hy == (np.sum(mask) - 1):  # i.e. last year in period
            info[hy, :] = np.amax(roll_30[i:(i + np.sum(mask) - 15), :], axis=0)
            i += np.sum(mask)
        else:
            info[hy, :] = np.amax(roll_30[i:(i + np.sum(mask)), :], axis=0)
            i += np.sum(mask)
    # calculations for entire time series
    sfc = np.mean(info, axis=0)

    return sfc


# DH5 - Mean annual maximum of 90-day average flow
def dh5(flows, datetimes, hydro_years, drainage_area):
    roll_90 = np.mean(rolling_window(flows, 90), axis=1)
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    i = 0
    for hy, mask in enumerate(hydro_years):
        if hy == 0:  # i.e. first year in period
            info[hy, :] = np.amax(roll_90[0:(np.sum(mask) - 44), :], axis=0)
            i += np.sum(mask) - 44
        elif hy == (np.sum(mask) - 1):  # i.e. last year in period
            info[hy, :] = np.amax(roll_90[i:(i + np.sum(mask) - 45), :], axis=0)
            i += np.sum(mask)
        else:
            info[hy, :] = np.amax(roll_90[i:(i + np.sum(mask)), :], axis=0)
            i += np.sum(mask)
    # calculations for entire time series
    sfc = np.mean(info, axis=0)

    return sfc


# DH6 - Variability in annual maximum average daily flow
def dh6(flows, datetimes, hydro_years, drainage_area):
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[hy, :] = np.amax(flows[mask, :], axis=0)
    # calculations for entire time series
    sfc = np.std(info, ddof=1, axis=0) * 100 / np.mean(info, axis=0)

    return sfc


# DH7 - Variability in annual maximum of 3-day average flow
def dh7(flows, datetimes, hydro_years, drainage_area):
    roll_3 = np.mean(rolling_window(flows, 3), axis=1)
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    i = 0
    for hy, mask in enumerate(hydro_years):
        if hy == 0:  # i.e. first year in period
            info[hy, :] = np.amax(roll_3[0:(np.sum(mask) - 1), :], axis=0)
            i += np.sum(mask) - 1
        elif hy == (np.sum(mask) - 1):  # i.e. last year in period
            info[hy, :] = np.amax(roll_3[i:(i + np.sum(mask) - 2), :], axis=0)
            i += np.sum(mask)
        else:
            info[hy, :] = np.amax(roll_3[i:(i + np.sum(mask)), :], axis=0)
            i += np.sum(mask)
    # calculations for entire time series
    sfc = np.std(info, ddof=1, axis=0) * 100 / np.mean(info, axis=0)

    return sfc


# DH8 - Variability in annual maximum of 7-day average flow
def dh8(flows, datetimes, hydro_years, drainage_area):
    roll_7 = np.mean(rolling_window(flows, 7), axis=1)
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    i = 0
    for hy, mask in enumerate(hydro_years):
        if hy == 0:  # i.e. first year in period
            info[hy, :] = np.amax(roll_7[0:(np.sum(mask) - 3), :], axis=0)
            i += np.sum(mask) - 3
        elif hy == (np.sum(mask) - 1):  # i.e. last year in period
            info[hy, :] = np.amax(roll_7[i:(i + np.sum(mask) - 4), :], axis=0)
            i += np.sum(mask)
        else:
            info[hy, :] = np.amax(roll_7[i:(i + np.sum(mask)), :], axis=0)
            i += np.sum(mask)
    # calculations for entire time series
    sfc = np.std(info, ddof=1, axis=0) * 100 / np.mean(info, axis=0)

    return sfc


# DH9 - Variability in annual maximum of 30-day average flow
def dh9(flows, datetimes, hydro_years, drainage_area):
    roll_30 = np.mean(rolling_window(flows, 30), axis=1)
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    i = 0
    for hy, mask in enumerate(hydro_years):
        if hy == 0:  # i.e. first year in period
            info[hy, :] = np.amax(roll_30[0:(np.sum(mask) - 14), :], axis=0)
            i += np.sum(mask) - 14
        elif hy == (np.sum(mask) - 1):  # i.e. last year in period
            info[hy, :] = np.amax(roll_30[i:(i + np.sum(mask) - 15), :], axis=0)
            i += np.sum(mask)
        else:
            info[hy, :] = np.amax(roll_30[i:(i + np.sum(mask)), :], axis=0)
            i += np.sum(mask)
    # calculations for entire time series
    sfc = np.std(info, ddof=1, axis=0) * 100 / np.mean(info, axis=0)

    return sfc


# DH10 - Mean annual maximum of 90-day average flow
def dh10(flows, datetimes, hydro_years, drainage_area):
    roll_90 = np.mean(rolling_window(flows, 90), axis=1)
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    i = 0
    for hy, mask in enumerate(hydro_years):
        if hy == 0:  # i.e. first year in period
            info[hy, :] = np.amax(roll_90[0:(np.sum(mask) - 44), :], axis=0)
            i += np.sum(mask) - 44
        elif hy == (np.sum(mask) - 1):  # i.e. last year in period
            info[hy, :] = np.amax(roll_90[i:(i + np.sum(mask) - 45), :], axis=0)
            i += np.sum(mask)
        else:
            info[hy, :] = np.amax(roll_90[i:(i + np.sum(mask)), :], axis=0)
            i += np.sum(mask)
    # calculations for entire time series
    sfc = np.std(info, ddof=1, axis=0) * 100 / np.mean(info, axis=0)

    return sfc


# DH11 - Mean annual minimum average flow normalised by median flow (of the whole record)
def dh11(flows, datetimes, hydro_years, drainage_area):
    median = np.median(flows, axis=0)
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[hy, :] = np.amax(flows[mask, :], axis=0)
    # calculations for entire time series
    sfc = np.mean(info, axis=0) / median

    return sfc


# DH12 - Mean annual minimum of 7-day average flow normalised by median flow (of the whole record)
def dh12(flows, datetimes, hydro_years, drainage_area):
    median = np.median(flows, axis=0)
    roll_7 = np.mean(rolling_window(flows, 7), axis=1)
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    i = 0
    for hy, mask in enumerate(hydro_years):
        if hy == 0:  # i.e. first year in period
            info[hy, :] = np.amax(roll_7[0:(np.sum(mask) - 3), :], axis=0)
            i += np.sum(mask) - 3
        elif hy == (np.sum(mask) - 1):  # i.e. last year in period
            info[hy, :] = np.amax(roll_7[i:(i + np.sum(mask) - 4), :], axis=0)
            i += np.sum(mask)
        else:
            info[hy, :] = np.amax(roll_7[i:(i + np.sum(mask)), :], axis=0)
            i += np.sum(mask)
    # calculations for entire time series
    sfc = np.mean(info, axis=0) / median

    return sfc


# DH13 - Annual maximum of 30-day average flow normalised by median flow (of the whole record)
def dh13(flows, datetimes, hydro_years, drainage_area):
    median = np.median(flows, axis=0)
    roll_30 = np.mean(rolling_window(flows, 30), axis=1)
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    i = 0
    for hy, mask in enumerate(hydro_years):
        if hy == 0:  # i.e. first year in period
            info[hy, :] = np.amax(roll_30[0:(np.sum(mask) - 14), :], axis=0)
            i += np.sum(mask) - 14
        elif hy == (np.sum(mask) - 1):  # i.e. last year in period
            info[hy, :] = np.amax(roll_30[i:(i + np.sum(mask) - 15), :], axis=0)
            i += np.sum(mask)
        else:
            info[hy, :] = np.amax(roll_30[i:(i + np.sum(mask)), :], axis=0)
            i += np.sum(mask)
    # calculations for entire time series
    sfc = np.mean(info, axis=0) / median

    return sfc


# DH14 - Annual maximum of 30-day average flow normalised by median flow (of the whole record)
def dh14(flows, datetimes, hydro_years, drainage_area):
    # calculations per month for each year
    mean = np.array(pd.DataFrame(flows, index=datetimes).groupby(lambda x: (x.year, x.month)).mean())
    # calculations for entire time series
    perc95 = np.percentile(mean, 95, axis=0)
    sfc = perc95 / np.mean(mean, axis=0)

    return sfc


# DH15 - Median high-flow pulse annual average duration
def dh15(flows, datetimes, hydro_years, drainage_area):
    perc75 = np.percentile(flows, 75, axis=0)
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[hy, :] = calc_events_avg_duration(flows[mask, :], perc75, typ='high')
    # calculations for entire time series
    sfc = np.median(info, axis=0)

    return sfc


# DH16 - Variability in high-flow pulse annual average duration
def dh16(flows, datetimes, hydro_years, drainage_area):
    perc75 = np.percentile(flows, 75, axis=0)
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[hy, :] = calc_events_avg_duration(flows[mask, :], perc75, typ='high')
    # calculations for entire time series
    sfc = np.std(info, ddof=1, axis=0) * 100 / np.mean(info, axis=0)

    return sfc


# DH17 - Average duration of high flow events for entire record above median flow
def dh17(flows, datetimes, hydro_years, drainage_area):
    median = np.median(flows, axis=0)
    # calculations for entire time series
    sfc = calc_events_avg_duration(flows, median, typ='high')

    return sfc


# DH18 - Average duration of high flow events for entire record above three times the median flow
def dh18(flows, datetimes, hydro_years, drainage_area):
    median = np.median(flows, axis=0)
    # calculations for entire time series
    sfc = calc_events_avg_duration(flows, 3 * median, typ='high')

    return sfc


# DH19 - Average duration of high flow events for entire record above seven times the median flow
def dh19(flows, datetimes, hydro_years, drainage_area):
    median = np.median(flows, axis=0)
    # calculations for entire time series
    sfc = calc_events_avg_duration(flows, 7 * median, typ='high')

    return sfc


# DH20 - Average duration of high flow events for entire record above Q25
def dh20(flows, datetimes, hydro_years, drainage_area):
    perc75 = np.percentile(flows, 75, axis=0)
    # calculations for entire time series
    sfc = calc_events_avg_duration(flows, perc75, typ='high')

    return sfc


# DH21 - Average duration of high flow events for entire record above Q75
def dh21(flows, datetimes, hydro_years, drainage_area):
    perc25 = np.percentile(flows, 25, axis=0)
    # calculations for entire time series
    sfc = calc_events_avg_duration(flows, perc25, typ='high')

    return sfc


# DH22 to DH24 not available
