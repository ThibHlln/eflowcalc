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
import warnings
from .tools import count_events, count_days


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOW FLOWS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# FL1 - Mean annual low flow pulse count (below 25th percentile)
def fl1(flows, datetimes, hydro_years, drainage_area):
    perc25 = np.percentile(flows, 25, axis=0)
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[hy, :] = count_events(flows[mask, :], threshold=perc25, typ='low')
    # calculations for entire time series
    sfc = np.mean(info, axis=0)

    return sfc


# FL2 - Variability in annual low flow pulse count
def fl2(flows, datetimes, hydro_years, drainage_area):
    perc25 = np.percentile(flows, 25, axis=0)
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[hy, :] = count_events(flows[mask, :], threshold=perc25, typ='low')
    # calculations for entire time series
    sfc = np.std(info, ddof=1, axis=0) * 100 / np.mean(info, axis=0)

    return sfc


# FL3 - Mean annual low flow pulse count (below 5% of the mean flow)
def fl3(flows, datetimes, hydro_years, drainage_area):
    mean_ = np.mean(flows, axis=0)
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[hy, :] = count_events(flows[mask, :], threshold=mean_ * 0.05, typ='low')
    # calculations for entire time series
    info[info <= 0] = np.nan
    with warnings.catch_warnings():  # info could be an empty slice that would rather an unnecessary warning
        warnings.simplefilter("ignore", category=RuntimeWarning)
        sfc = np.nanmean(info, axis=0)
    sfc[np.isnan(sfc)] = 0.0

    return sfc


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# HIGH FLOWS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# FH1 - Mean annual high flow pulse count (above 75th percentile)
def fh1(flows, datetimes, hydro_years, drainage_area):
    perc75 = np.percentile(flows, 75, axis=0)
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[hy, :] = count_events(flows[mask, :], threshold=perc75, typ='high')
    # calculations for entire time series
    sfc = np.mean(info, axis=0)

    return sfc


# FH2 - Variability in annual high flow pulse count
def fh2(flows, datetimes, hydro_years, drainage_area):
    perc75 = np.percentile(flows, 75, axis=0)
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[hy, :] = count_events(flows[mask, :], threshold=perc75, typ='high')
    # calculations for entire time series
    sfc = np.std(info, ddof=1, axis=0) * 100 / np.mean(info, axis=0)

    return sfc


# FH3 - Mean number of days per year with moderate floods
def fh3(flows, datetimes, hydro_years, drainage_area):
    median_ = np.median(flows, axis=0)
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[hy, :] = count_days(flows[mask, :], threshold=median_ * 3, typ='high')
    # calculations for entire time series
    sfc = np.mean(info, axis=0)

    return sfc


# FH4 - Mean number of days per year with large floods
def fh4(flows, datetimes, hydro_years, drainage_area):
    median_ = np.median(flows, axis=0)
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[hy, :] = count_days(flows[mask, :], threshold=median_ * 7, typ='high')
    # calculations for entire time series
    sfc = np.mean(info, axis=0)

    return sfc


# FH5 - Frequency of events above median flow for entire record
def fh5(flows, datetimes, hydro_years, drainage_area):
    median_ = np.median(flows, axis=0)
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[hy, :] = count_events(flows[mask, :], median_, typ='high')
    # calculations for entire time series
    info[info == 0] = np.nan
    sfc = np.nanmean(info, axis=0)

    return sfc


# FH6 - Frequency of moderate floods
def fh6(flows, datetimes, hydro_years, drainage_area):
    median_ = np.median(flows, axis=0)
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[hy, :] = count_events(flows[mask, :], median_ * 3, typ='high')
    # calculations for entire time series
    sfc = np.mean(info, axis=0)

    return sfc


# FH7 - Frequency of large floods
def fh7(flows, datetimes, hydro_years, drainage_area):
    median_ = np.median(flows, axis=0)
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[hy, :] = count_events(flows[mask, :], median_ * 7, typ='high')
    # calculations for entire time series
    sfc = np.mean(info, axis=0)

    return sfc


# FH8 - Frequency of events beyond Q25
def fh8(flows, datetimes, hydro_years, drainage_area):
    perc75 = np.percentile(flows, 75, axis=0)
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[hy, :] = count_events(flows[mask, :], threshold=perc75, typ='high')
    # calculations for entire time series
    sfc = np.mean(info, axis=0)

    return sfc


# FH9 - Frequency of events beyond Q75
def fh9(flows, datetimes, hydro_years, drainage_area):
    perc25 = np.percentile(flows, 25, axis=0)
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[hy, :] = count_events(flows[mask, :], threshold=perc25, typ='high')
    # calculations for entire time series
    sfc = np.mean(info, axis=0)

    return sfc


# FH10 - Frequency of events beyond Q75
def fh10(flows, datetimes, hydro_years, drainage_area):
    # calculations per hydrological year
    min_ = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        min_[hy, :] = np.amin(flows[mask, :], axis=0)
    median_ = np.median(min_, axis=0)
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[hy, :] = count_events(flows[mask, :], threshold=median_, typ='high')
    # calculations for entire time series
    info[info <= 0] = np.nan
    sfc = np.nanmean(info, axis=0)

    return sfc


# FH11 not available
