# This file is part of EFlowCalc:
# A Calculator of Ecological Streamflow Characteristics
# Copyright (C) 2020  Thibault Hallouin (1)
#
# (1) Dooge Centre for Water Resources Research,
#     University College Dublin, Ireland
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


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ALL FLOWS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def ra1(flows, datetimes, hydro_years, drainage_area):
    """Average rise rate.

    :Calculation Details:
        Identify the flow rises (i.e. the increases in daily flows from
        one day to the next). Calculate the mean of these rises.

    """
    # calculations for entire time series
    diffs = np.diff(flows, axis=0)
    rises = np.copy(diffs)
    rises[rises <= 0] = np.nan
    sfc = np.nanmean(rises, axis=0)

    return sfc


def ra2(flows, datetimes, hydro_years, drainage_area):
    """Variability in rise rate.

    :Calculation Details:
        Identify the flow rises (i.e. the increases in daily flows from
        one day to the next). Calculate the standard deviation and the
        mean of these rises, multiply the former by 100, and divide the
        result by the latter.

    """
    # calculations for entire time series
    diffs = np.diff(flows, axis=0)
    rises = np.copy(diffs)
    rises[rises <= 0] = np.nan
    sfc = np.nanstd(rises, ddof=1, axis=0) * 100 / np.nanmean(rises, axis=0)

    return sfc


def ra3(flows, datetimes, hydro_years, drainage_area):
    """Average fall rate.

    :Calculation Details:
        Identify the flow falls (i.e. the decreases in daily flows from
        one day to the next). Calculate the mean of these falls.

    """
    # calculations for entire time series
    diffs = np.diff(flows, axis=0)
    falls = np.copy(diffs)
    falls[falls >= 0] = np.nan
    falls = np.abs(falls)
    sfc = np.nanmean(falls, axis=0)

    return sfc


def ra4(flows, datetimes, hydro_years, drainage_area):
    """Variability in fall rate.

    :Calculation Details:
        Identify the flow falls (i.e. the decreases in daily flows from
        one day to the next). Calculate the standard deviation and the
        mean of these falls, multiply the former by 100, and divide the
        result by the latter.

    """
    # calculations for entire time series
    diffs = np.diff(flows, axis=0)
    falls = np.copy(diffs)
    falls[falls >= 0] = np.nan
    falls = np.abs(falls)
    sfc = np.nanstd(falls, ddof=1, axis=0) * 100 / np.nanmean(falls, axis=0)

    return sfc


def ra5(flows, datetimes, hydro_years, drainage_area):
    """Ratio of days with flow rise.

    :Calculation Details:
        Calculate the number of days with flow rises (i.e. the increases
        in daily flows from one day to the next). Divide this number by
        the total number of days in the daily flow record.

    """
    # calculations for entire time series
    diffs = np.diff(flows, axis=0)
    rises = np.copy(diffs)
    sfc = np.true_divide(np.sum(rises > 0, axis=0), flows.shape[0])

    return sfc


def ra6(flows, datetimes, hydro_years, drainage_area):
    """Change of flow rises.

    :Calculation Details:
        Compute the natural logarithm of the daily flow record. Identify
        the rises in log-transformed flows (i.e. the increases in daily
        log flows from one day to the next). Calculate the median of
        these rises.

    """
    # calculations for entire time series
    cp_flows = np.copy(flows)
    # replace 0 by 0.01 if necessary (to avoid log(0))
    cp_flows[cp_flows == 0.0] = 0.01
    diffs_log = np.diff(np.log(cp_flows), axis=0)
    diffs_log[diffs_log <= 0] = np.nan  # take rises only
    sfc = np.nanmedian(np.abs(diffs_log), axis=0)

    return sfc


def ra7(flows, datetimes, hydro_years, drainage_area):
    """Change of flow falls.

    :Calculation Details:
        Compute the natural logarithm of the daily flow record. Identify
        the falls in log-transformed flows (i.e. the decreases in daily
        log flows from one day to the next). Calculate the median of
        these falls.

    """
    # calculations for entire time series
    cp_flows = np.copy(flows)
    # replace 0 by 0.01 if necessary (to avoid log(0))
    cp_flows[cp_flows == 0.0] = 0.01
    diffs_log = np.diff(np.log(cp_flows), axis=0)
    diffs_log[diffs_log >= 0] = np.nan  # take falls only
    sfc = np.nanmedian(np.abs(diffs_log), axis=0)

    return sfc


def ra8(flows, datetimes, hydro_years, drainage_area):
    """Average annual number of reversals.

    :Calculation Details:
        Count the number of flow reversals (change from a period of
        increasing daily flows from one day to the next to a period of
        decreasing daily flows from one day to the next) for each
        hydrological year in the daily flow record. Calculate the mean
        of these numbers.

    """
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[hy, :] = count_reversals(flows[mask, :])
    # calculations for entire time series
    sfc = np.mean(info[:, :], axis=0)

    return sfc


def ra9(flows, datetimes, hydro_years, drainage_area):
    """Variability in annual number of reversals.

    :Calculation Details:
        Count the number of flow reversals (change from a period of
        increasing daily flows from one day to the next to a period of
        decreasing daily flows from one day to the next) for each
        hydrological year in the daily flow record. Calculate the
        standard deviation and the mean of these numbers. Multiply the
        former by 100, and divide the result by the latter.

    """
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[hy, :] = count_reversals(flows[mask, :])
    # calculations for entire time series
    sfc = np.std(info, ddof=1, axis=0) * 100 / np.mean(info, axis=0)

    return sfc
