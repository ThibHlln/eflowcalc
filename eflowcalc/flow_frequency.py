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
import warnings
from .tools import count_events, count_days


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOW FLOWS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def fl1(flows, datetimes, hydro_years, drainage_area):
    """Mean annual low-flow pulse count (below 25th percentile).

    :Calculation Details:
        Identify for each hydrological year the flow events below the
        25th percentile of the whole daily flow record. Count the
        number of these flow events for each hydrological year.
        Calculate the mean of these numbers.

    """
    perc25 = np.percentile(flows, 25, axis=0)
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[hy, :] = count_events(flows[mask, :], threshold=perc25, typ='low')
    # calculations for entire time series
    sfc = np.mean(info, axis=0)

    return sfc


def fl2(flows, datetimes, hydro_years, drainage_area):
    """Variability in annual low-flow pulse count (below 25th percentile).

    :Calculation Details:
        Identify for each hydrological year the flow events below the
        25th percentile of the whole daily flow record. Count the
        number of these flow events for each hydrological year.
        Calculate the standard deviation and the mean of these number of
        events. Multiply the former by 100 and divide the result by the
        latter.

    """
    perc25 = np.percentile(flows, 25, axis=0)
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[hy, :] = count_events(flows[mask, :], threshold=perc25, typ='low')
    # calculations for entire time series
    sfc = np.std(info, ddof=1, axis=0) * 100 / np.mean(info, axis=0)

    return sfc


def fl3(flows, datetimes, hydro_years, drainage_area):
    """Mean annual low-flow pulse count (below 5% of overall mean flow).

    :Calculation Details:
        Identify for each hydrological year the flow events below 5% of
        the mean flow of the whole daily flow record. Count the number
        of these flow events for each hydrological year. Calculate the
        mean of these number of events.

    """
    mean_ = np.mean(flows, axis=0)
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[hy, :] = count_events(flows[mask, :], threshold=mean_ * 0.05,
                                   typ='low')
    # calculations for entire time series
    info[info <= 0] = np.nan
    with warnings.catch_warnings():
        # info could be an empty slice that would rather an
        # unnecessary warning
        warnings.simplefilter("ignore", category=RuntimeWarning)
        sfc = np.nanmean(info, axis=0)
    sfc[np.isnan(sfc)] = 0.0

    return sfc


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# HIGH FLOWS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def fh1(flows, datetimes, hydro_years, drainage_area):
    """Mean annual high-flow pulse count (above 75th percentile).

    *Note:* this streamflow characteristic is equivalent to `fh8`.

    :Calculation Details:
        Identify for each hydrological year the flow events above the
        75th percentile of the whole daily flow record. Count the
        number of these flow events for each hydrological year.
        Calculate the mean of these number of events.

    """
    perc75 = np.percentile(flows, 75, axis=0)
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]),
                    dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[hy, :] = count_events(flows[mask, :], threshold=perc75,
                                   typ='high')
    # calculations for entire time series
    sfc = np.mean(info, axis=0)

    return sfc


def fh2(flows, datetimes, hydro_years, drainage_area):
    """Variability in annual low-flow pulse count (above 75th percentile).

    :Calculation Details:
        Identify for each hydrological year the flow events above the
        75th percentile of the whole daily flow record. Count the
        number of these flow events for each hydrological year.
        Calculate the standard deviation and the mean of these number of
        events. Multiply the former by 100 and divide the result by the
        latter.

    """
    perc75 = np.percentile(flows, 75, axis=0)
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[hy, :] = count_events(flows[mask, :], threshold=perc75,
                                   typ='high')
    # calculations for entire time series
    sfc = np.std(info, ddof=1, axis=0) * 100 / np.mean(info, axis=0)

    return sfc


def fh3(flows, datetimes, hydro_years, drainage_area):
    """Mean number of days per year with moderate floods.

    :Calculation Details:
        Count for each hydrological year the number of days where flow
        above 3 times the median flow value of the whole record.
        Calculate the mean of these number of days.

    """
    median_ = np.median(flows, axis=0)
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[hy, :] = count_days(flows[mask, :], threshold=median_ * 3,
                                 typ='high')
    # calculations for entire time series
    sfc = np.mean(info, axis=0)

    return sfc


def fh4(flows, datetimes, hydro_years, drainage_area):
    """Mean number of days per year with large floods.

    :Calculation Details:
        Count for each hydrological year the number of days where flow
        above 7 times the median flow value of the whole record.
        Calculate the mean of these number of days.

    """
    median_ = np.median(flows, axis=0)
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[hy, :] = count_days(flows[mask, :], threshold=median_ * 7,
                                 typ='high')
    # calculations for entire time series
    sfc = np.mean(info, axis=0)

    return sfc


def fh5(flows, datetimes, hydro_years, drainage_area):
    """Frequency of events above overall median daily flow.

    :Calculation Details:
        Count for each hydrological year the number of flow events
        above the median flow value of the whole record. Calculate the
        mean of these number of days.

    """
    median_ = np.median(flows, axis=0)
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[hy, :] = count_events(flows[mask, :], median_, typ='high')
    # calculations for entire time series
    info[info == 0] = np.nan
    sfc = np.nanmean(info, axis=0)

    return sfc


def fh6(flows, datetimes, hydro_years, drainage_area):
    """Frequency of moderate floods.

    :Calculation Details:
        Count for each hydrological year the number of flow events
        above 3 times the median flow value of the whole record.
        Calculate the mean of these number of events.

    """
    median_ = np.median(flows, axis=0)
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[hy, :] = count_events(flows[mask, :], median_ * 3, typ='high')
    # calculations for entire time series
    sfc = np.mean(info, axis=0)

    return sfc


def fh7(flows, datetimes, hydro_years, drainage_area):
    """Frequency of large floods.

    :Calculation Details:
        Count for each hydrological year the number of flow events
        above 7 times the median flow value of the whole record.
        Calculate the mean of these number of events.

    """
    median_ = np.median(flows, axis=0)
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[hy, :] = count_events(flows[mask, :], median_ * 7, typ='high')
    # calculations for entire time series
    sfc = np.mean(info, axis=0)

    return sfc


def fh8(flows, datetimes, hydro_years, drainage_area):
    """Frequency of events above Q25.

    *Note:* this streamflow characteristic is equivalent to `fh1`.

    :Calculation Details:
        Identify for each hydrological year the flow events above the
        75th percentile of the whole daily flow record. Count the
        number of these flow events for each hydrological year.
        Calculate the mean of these number of events.

    """
    perc75 = np.percentile(flows, 75, axis=0)
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[hy, :] = count_events(flows[mask, :], threshold=perc75,
                                   typ='high')
    # calculations for entire time series
    sfc = np.mean(info, axis=0)

    return sfc


def fh9(flows, datetimes, hydro_years, drainage_area):
    """Frequency of events above Q75.

    :Calculation Details:
        Identify for each hydrological year the flow events above the
        25th percentile of the whole daily flow record. Count the
        number of these flow events for each hydrological year.
        Calculate the mean of these number of events.

    """
    perc25 = np.percentile(flows, 25, axis=0)
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[hy, :] = count_events(flows[mask, :], threshold=perc25,
                                   typ='high')
    # calculations for entire time series
    sfc = np.mean(info, axis=0)

    return sfc


def fh10(flows, datetimes, hydro_years, drainage_area):
    """Frequency of events beyond median of annual minima.

    :Calculation Details:
        Take the minimum daily flow for each hydrological year. Compute
        the median of these minima. Identify for each hydrological year
        the flow events above this median. Count the number of these
        flow events for each hydrological year. Calculate the mean of
        these number of events.

    """
    # calculations per hydrological year
    min_ = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        min_[hy, :] = np.amin(flows[mask, :], axis=0)
    median_ = np.median(min_, axis=0)
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[hy, :] = count_events(flows[mask, :], threshold=median_,
                                   typ='high')
    # calculations for entire time series
    info[info <= 0] = np.nan
    sfc = np.nanmean(info, axis=0)

    return sfc


# FH11 not available
