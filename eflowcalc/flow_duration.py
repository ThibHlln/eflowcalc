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
import pandas as pd
from .tools import rolling_window, calc_events_avg_duration


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOW FLOWS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def dl1(flows, datetimes, hydro_years, drainage_area):
    """Mean annual minimum daily flow.

    :Calculation Details:
        Take the minimum flow for each hydrological year. Calculate
        the mean of these minimum values.

    """
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[hy, :] = np.amin(flows[mask, :], axis=0)
    # calculations for entire time series
    sfc = np.mean(info, axis=0)

    return sfc


def dl2(flows, datetimes, hydro_years, drainage_area):
    """Mean annual minimum of 3-day daily flow.

    :Calculation Details:
        Compute the 3-day rolling mean flows. Take the minimum of these
        rolling means for each hydrological year. Calculate the mean of
        these minimum values.

    """
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


def dl3(flows, datetimes, hydro_years, drainage_area):
    """Mean annual minimum of 7-day daily flow.

    :Calculation Details:
        Compute the 7-day rolling mean flows. Take the minimum of these
        rolling means for each hydrological year. Calculate the mean of
        these minimum values.

    """
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


def dl4(flows, datetimes, hydro_years, drainage_area):
    """Mean annual minimum of 30-day daily flow.

    :Calculation Details:
        Compute the 30-day rolling mean flows. Take the minimum of these
        rolling means for each hydrological year. Calculate the mean of
        these minimum values.

    """
    roll_30 = np.mean(rolling_window(flows, 30), axis=1)
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    i = 0
    for hy, mask in enumerate(hydro_years):
        if hy == 0:  # i.e. first year in period
            info[hy, :] = np.amin(roll_30[0:(np.sum(mask) - 14), :], axis=0)
            i += np.sum(mask) - 14
        elif hy == (np.sum(mask) - 1):  # i.e. last year in period
            info[hy, :] = np.amin(roll_30[i:(i + np.sum(mask) - 15), :],
                                  axis=0)
            i += np.sum(mask)
        else:
            info[hy, :] = np.amin(roll_30[i:(i + np.sum(mask)), :], axis=0)
            i += np.sum(mask)
    # calculations for entire time series
    sfc = np.mean(info, axis=0)

    return sfc


def dl5(flows, datetimes, hydro_years, drainage_area):
    """Mean annual minimum of 90-day daily flow.

    :Calculation Details:
        Compute the 90-day rolling mean flows. Take the minimum of these
        rolling means for each hydrological year. Calculate the mean of
        these minimum values.

    """
    roll_90 = np.mean(rolling_window(flows, 90), axis=1)
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    i = 0
    for hy, mask in enumerate(hydro_years):
        if hy == 0:  # i.e. first year in period
            info[hy, :] = np.amin(roll_90[0:(np.sum(mask) - 44), :], axis=0)
            i += np.sum(mask) - 44
        elif hy == (np.sum(mask) - 1):  # i.e. last year in period
            info[hy, :] = np.amin(roll_90[i:(i + np.sum(mask) - 45), :],
                                  axis=0)
            i += np.sum(mask)
        else:
            info[hy, :] = np.amin(roll_90[i:(i + np.sum(mask)), :], axis=0)
            i += np.sum(mask)
    # calculations for entire time series
    sfc = np.mean(info, axis=0)

    return sfc


def dl6(flows, datetimes, hydro_years, drainage_area):
    """Variability in annual minimum daily flow.

    *Note:* this streamflow characteristic is equivalent to `ml21`.

    :Calculation Details:
        Take the minimum flow for each hydrological year. Calculate
        the standard deviation and mean of these minimum values.
        Multiply the former by 100 and divide the result by the latter.

    """
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[hy, :] = np.amin(flows[mask, :], axis=0)
    # calculations for entire time series
    sfc = np.std(info, ddof=1, axis=0) * 100 / np.mean(info, axis=0)

    return sfc


def dl7(flows, datetimes, hydro_years, drainage_area):
    """Variability in annual minimum of 3-day daily flow.

    :Calculation Details:
        Compute the 3-day rolling mean flows. Take the minimum of these
        rolling means for each hydrological year. Calculate the standard
        deviation and mean of these minimum values. Multiply the former
        by 100 and divide the result by the latter.

    """
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


def dl8(flows, datetimes, hydro_years, drainage_area):
    """Variability in annual minimum of 7-day daily flow.

    :Calculation Details:
        Compute the 7-day rolling mean flows. Take the minimum of these
        rolling means for each hydrological year. Calculate the standard
        deviation and the mean of these minimum values. Multiply the
        former by 100 and divide the result by the latter.

    """
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


def dl9(flows, datetimes, hydro_years, drainage_area):
    """Variability in annual minimum of 30-day daily flow.

    :Calculation Details:
        Compute the 30-day rolling mean flows. Take the minimum of these
        rolling means for each hydrological year. Calculate the standard
        deviation and the mean of these minimum values. Multiply the
        former by 100 and divide the result by the latter.

    """
    roll_30 = np.mean(rolling_window(flows, 30), axis=1)
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    i = 0
    for hy, mask in enumerate(hydro_years):
        if hy == 0:  # i.e. first year in period
            info[hy, :] = np.amin(roll_30[0:(np.sum(mask) - 14), :], axis=0)
            i += np.sum(mask) - 14
        elif hy == (np.sum(mask) - 1):  # i.e. last year in period
            info[hy, :] = np.amin(roll_30[i:(i + np.sum(mask) - 15), :],
                                  axis=0)
            i += np.sum(mask)
        else:
            info[hy, :] = np.amin(roll_30[i:(i + np.sum(mask)), :], axis=0)
            i += np.sum(mask)
    # calculations for entire time series
    sfc = np.std(info, ddof=1, axis=0) * 100 / np.mean(info, axis=0)

    return sfc


def dl10(flows, datetimes, hydro_years, drainage_area):
    """Variability in annual minimum of 90-day daily flow.

    :Calculation Details:
        Compute the 90-day rolling mean flows. Take the minimum of these
        rolling means for each hydrological year. Calculate the standard
        deviation and the mean of these minimum values. Multiply the
        former by 100 and divide the result by the latter.

    """
    roll_90 = np.mean(rolling_window(flows, 90), axis=1)
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    i = 0
    for hy, mask in enumerate(hydro_years):
        if hy == 0:  # i.e. first year in period
            info[hy, :] = np.amin(roll_90[0:(np.sum(mask) - 44), :], axis=0)
            i += np.sum(mask) - 14
        elif hy == (np.sum(mask) - 1):  # i.e. last year in period
            info[hy, :] = np.amin(roll_90[i:(i + np.sum(mask) - 45), :],
                                  axis=0)
            i += np.sum(mask)
        else:
            info[hy, :] = np.amin(roll_90[i:(i + np.sum(mask)), :], axis=0)
            i += np.sum(mask)
    # calculations for entire time series
    sfc = np.std(info, ddof=1, axis=0) * 100 / np.mean(info, axis=0)

    return sfc


def dl11(flows, datetimes, hydro_years, drainage_area):
    """Mean annual minimum daily flow normalised by overall median
    daily flow.

    :Calculation Details:
        Take the minimum flow for each hydrological year. Calculate
        the mean of these minimum values. Divide this mean value by the
        median of the whole daily flow record.

    """
    median = np.median(flows, axis=0)
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[hy, :] = np.amin(flows[mask, :], axis=0)
    # calculations for entire time series
    sfc = np.mean(info, axis=0) / median

    return sfc


def dl12(flows, datetimes, hydro_years, drainage_area):
    """Mean annual minimum of 7-day daily flow normalised by overall
    median daily flow.

    :Calculation Details:
        Compute the 7-day rolling mean flows. Take the minimum of these
        rolling means for each hydrological year. Calculate the mean of
        these minimum values. Divide this mean value by the median of
        the whole daily flow record.

    """
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


def dl13(flows, datetimes, hydro_years, drainage_area):
    """Mean annual minimum of 30-day daily flow normalised by overall
    median daily flow.

    :Calculation Details:
        Compute the 30-day rolling mean flows. Take the minimum of these
        rolling means for each hydrological year. Calculate the mean of
        these minimum values. Divide this mean value by the median of
        the whole daily flow record.

    """
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
            info[hy, :] = np.amin(roll_30[i:(i + np.sum(mask) - 15), :],
                                  axis=0)
            i += np.sum(mask)
        else:
            info[hy, :] = np.amin(roll_30[i:(i + np.sum(mask)), :], axis=0)
            i += np.sum(mask)
    # calculations for entire time series
    sfc = np.mean(info, axis=0) / median

    return sfc


def dl14(flows, datetimes, hydro_years, drainage_area):
    """Q75 exceedance value normalised by overall median daily flow.

    :Calculation Details:
        Calculate the 25th percentile of the whole daily flow record.
        Divide this percentile by the median of the whole daily flow
        record.

    """
    median = np.median(flows, axis=0)
    perc25 = np.percentile(flows, 25, axis=0)

    sfc = perc25 / median

    return sfc


def dl15(flows, datetimes, hydro_years, drainage_area):
    """Q90 exceedance value normalised by overall median daily flow.

    :Calculation Details:
        Calculate the 10th percentile of the whole daily flow record.
        Divide this percentile by the median of the whole daily flow
        record.

    """

    median = np.median(flows, axis=0)
    perc10 = np.percentile(flows, 10, axis=0)

    sfc = perc10 / median

    return sfc


def dl16(flows, datetimes, hydro_years, drainage_area):
    """Median low-flow pulse duration.

    :Calculation Details:
        Identify for each hydrological year the flow events below the
        25th percentile of the whole daily flow record. Calculate the
        mean duration of these flow events for hydrological year. Take
        the median of these mean values.

    """
    perc25 = np.percentile(flows, 25, axis=0)
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[hy, :] = calc_events_avg_duration(flows[mask, :], perc25,
                                               typ='low')
    # calculations for entire time series
    sfc = np.median(info, axis=0)

    return sfc


def dl17(flows, datetimes, hydro_years, drainage_area):
    """Variability in low-flow pulse duration.

    :Calculation Details:
        Identify for each hydrological year the flow events below the
        25th percentile of the whole daily flow record. Calculate the
        mean duration of these flow events for each hydrological year.
        Calculate the standard deviation and the mean of these mean
        values. Multiply the former by 100 and divide the result by the
        latter.

    """
    perc25 = np.percentile(flows, 25, axis=0)
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[hy, :] = calc_events_avg_duration(flows[mask, :], perc25,
                                               typ='low')
    # calculations for entire time series
    sfc = np.std(info, ddof=1, axis=0) * 100 / np.mean(info, axis=0)

    return sfc


def dl18(flows, datetimes, hydro_years, drainage_area):
    """Mean annual number of zero flow days.

    :Calculation Details:
        Calculate the number of days with mean daily flow of zero for
        each hydrological year. Calculate the mean of these numbers.

    """
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[hy, :] = np.sum(flows[mask, :] == 0, axis=0)
    # calculations for entire time series
    sfc = np.mean(info, axis=0)

    return sfc


def dl19(flows, datetimes, hydro_years, drainage_area):
    """Variability in annual number of zero flow days.

    :Calculation Details:
        Calculate the number of days with mean daily flow of zero for
        each hydrological year. Calculate the standard deviation and the
        mean of these numbers. Multiply the former by 100 and divide the
        result by the latter.

    """
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[hy, :] = np.sum(flows[mask, :] == 0, axis=0)
    # calculations for entire time series
    mean_ = np.mean(info, axis=0)
    sfc = np.true_divide(np.std(info, axis=0) * 100, mean_, where=(mean_ != 0))
    sfc[mean_ == 0] = 0.0

    return sfc


def dl20(flows, datetimes, hydro_years, drainage_area):
    """Mean annual number of zero flow days.

    :Calculation Details:
        Calculate the monthly flows for the whole flow record. Compute
        the number of months with monthly flow of zero.

    """
    # calculations per month for each year
    zero_flow = np.array(pd.DataFrame(flows != 0, index=datetimes).groupby(
        lambda x: (x.year, x.month)).sum())
    # calculations for entire time series
    sfc = np.sum(zero_flow == 0, axis=0)

    return sfc


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# HIGH FLOWS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def dh1(flows, datetimes, hydro_years, drainage_area):
    """Mean annual maximum daily flow.

    :Calculation Details:
        Take the maximum flow for each hydrological year. Calculate
        the mean of these maximum values.

    """
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[hy, :] = np.amax(flows[mask, :], axis=0)
    # calculations for entire time series
    sfc = np.mean(info, axis=0)

    return sfc


def dh2(flows, datetimes, hydro_years, drainage_area):
    """Mean annual maximum of 3-day daily flow.

    :Calculation Details:
        Compute the 3-day rolling mean flows. Take the maximum of these
        rolling means for each hydrological year. Calculate the mean of
        these maximum values.

    """
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


def dh3(flows, datetimes, hydro_years, drainage_area):
    """Mean annual maximum of 7-day daily flow.

    :Calculation Details:
        Compute the 7-day rolling mean flows. Take the maximum of these
        rolling means for each hydrological year. Calculate the mean of
        these maximum values.

    """
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


def dh4(flows, datetimes, hydro_years, drainage_area):
    """Mean annual maximum of 30-day daily flow.

    :Calculation Details:
        Compute the 30-day rolling mean flows. Take the maximum of these
        rolling means for each hydrological year. Calculate the mean of
        these maximum values.

    """
    roll_30 = np.mean(rolling_window(flows, 30), axis=1)
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    i = 0
    for hy, mask in enumerate(hydro_years):
        if hy == 0:  # i.e. first year in period
            info[hy, :] = np.amax(roll_30[0:(np.sum(mask) - 14), :], axis=0)
            i += np.sum(mask) - 14
        elif hy == (np.sum(mask) - 1):  # i.e. last year in period
            info[hy, :] = np.amax(roll_30[i:(i + np.sum(mask) - 15), :],
                                  axis=0)
            i += np.sum(mask)
        else:
            info[hy, :] = np.amax(roll_30[i:(i + np.sum(mask)), :], axis=0)
            i += np.sum(mask)
    # calculations for entire time series
    sfc = np.mean(info, axis=0)

    return sfc


def dh5(flows, datetimes, hydro_years, drainage_area):
    """Mean annual maximum of 90-day daily flow.

    :Calculation Details:
        Compute the 90-day rolling mean flows. Take the maximum of these
        rolling means for each hydrological year. Calculate the mean of
        these maximum values.

    """
    roll_90 = np.mean(rolling_window(flows, 90), axis=1)
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    i = 0
    for hy, mask in enumerate(hydro_years):
        if hy == 0:  # i.e. first year in period
            info[hy, :] = np.amax(roll_90[0:(np.sum(mask) - 44), :], axis=0)
            i += np.sum(mask) - 44
        elif hy == (np.sum(mask) - 1):  # i.e. last year in period
            info[hy, :] = np.amax(roll_90[i:(i + np.sum(mask) - 45), :],
                                  axis=0)
            i += np.sum(mask)
        else:
            info[hy, :] = np.amax(roll_90[i:(i + np.sum(mask)), :], axis=0)
            i += np.sum(mask)
    # calculations for entire time series
    sfc = np.mean(info, axis=0)

    return sfc


def dh6(flows, datetimes, hydro_years, drainage_area):
    """Variability in annual maximum daily flow.

    :Calculation Details:
        Take the maximum flow for each hydrological year. Calculate
        the standard deviation and the mean of these maximum values.
        Multiply the former by 100 and divide the result by the latter.

    """
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[hy, :] = np.amax(flows[mask, :], axis=0)
    # calculations for entire time series
    sfc = np.std(info, ddof=1, axis=0) * 100 / np.mean(info, axis=0)

    return sfc


def dh7(flows, datetimes, hydro_years, drainage_area):
    """Variability in annual maximum of 3-day daily flow.

    :Calculation Details:
        Compute the 3-day rolling mean flows. Take the maximum of these
        rolling means for each hydrological year. Calculate the standard
        deviation and the mean of these maximum values. Multiply the
        former by 100 and divide the result by the latter.

    """
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


def dh8(flows, datetimes, hydro_years, drainage_area):
    """Variability in annual maximum of 7-day daily flow.

    :Calculation Details:
        Compute the 7-day rolling mean flows. Take the maximum of these
        rolling means for each hydrological year. Calculate the standard
        deviation and the mean of these maximum values. Multiply the
        former by 100 and divide the result by the latter.

    """
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


def dh9(flows, datetimes, hydro_years, drainage_area):
    """Variability in annual maximum of 30-day daily flow.

    :Calculation Details:
        Compute the 30-day rolling mean flows. Take the maximum of these
        rolling means for each hydrological year. Calculate the standard
        deviation and the mean of these maximum values. Multiply the
        former by 100 and divide the result by the latter.

    """
    roll_30 = np.mean(rolling_window(flows, 30), axis=1)
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    i = 0
    for hy, mask in enumerate(hydro_years):
        if hy == 0:  # i.e. first year in period
            info[hy, :] = np.amax(roll_30[0:(np.sum(mask) - 14), :], axis=0)
            i += np.sum(mask) - 14
        elif hy == (np.sum(mask) - 1):  # i.e. last year in period
            info[hy, :] = np.amax(roll_30[i:(i + np.sum(mask) - 15), :],
                                  axis=0)
            i += np.sum(mask)
        else:
            info[hy, :] = np.amax(roll_30[i:(i + np.sum(mask)), :], axis=0)
            i += np.sum(mask)
    # calculations for entire time series
    sfc = np.std(info, ddof=1, axis=0) * 100 / np.mean(info, axis=0)

    return sfc


def dh10(flows, datetimes, hydro_years, drainage_area):
    """Variability in annual maximum of 90-day daily flow.

    :Calculation Details:
        Compute the 90-day rolling mean flows. Take the maximum of these
        rolling means for each hydrological year. Calculate the standard
        deviation and the mean of these maximum values. Multiply the
        former by 100 and divide the result by the latter.

    """
    roll_90 = np.mean(rolling_window(flows, 90), axis=1)
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    i = 0
    for hy, mask in enumerate(hydro_years):
        if hy == 0:  # i.e. first year in period
            info[hy, :] = np.amax(roll_90[0:(np.sum(mask) - 44), :], axis=0)
            i += np.sum(mask) - 44
        elif hy == (np.sum(mask) - 1):  # i.e. last year in period
            info[hy, :] = np.amax(roll_90[i:(i + np.sum(mask) - 45), :],
                                  axis=0)
            i += np.sum(mask)
        else:
            info[hy, :] = np.amax(roll_90[i:(i + np.sum(mask)), :], axis=0)
            i += np.sum(mask)
    # calculations for entire time series
    sfc = np.std(info, ddof=1, axis=0) * 100 / np.mean(info, axis=0)

    return sfc


def dh11(flows, datetimes, hydro_years, drainage_area):
    """Mean annual maximum daily flow normalised by overall median
    daily flow.

    :Calculation Details:
        Take the maximum flow for each hydrological year. Calculate
        the mean of these maximum values. Divide this mean value by the
        median of the whole daily flow record.

    """
    median = np.median(flows, axis=0)
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[hy, :] = np.amax(flows[mask, :], axis=0)
    # calculations for entire time series
    sfc = np.mean(info, axis=0) / median

    return sfc


def dh12(flows, datetimes, hydro_years, drainage_area):
    """Mean annual maximum of 7-day daily flow normalised by overall
    median daily flow.

    :Calculation Details:
        Compute the 7-day rolling mean flows. Take the maximum of these
        rolling means for each hydrological year. Calculate the mean of
        these maximum values. Divide this mean value by the median of
        the whole daily flow record.

    """
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


def dh13(flows, datetimes, hydro_years, drainage_area):
    """Mean annual maximum of 30-day daily flow normalised by overall
    median daily flow.

    :Calculation Details:
        Compute the 30-day rolling mean flows. Take the maximum of these
        rolling means for each hydrological year. Calculate the mean of
        these maximum values. Divide this mean value by the median of
        the whole daily flow record.

    """
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
            info[hy, :] = np.amax(roll_30[i:(i + np.sum(mask) - 15), :],
                                  axis=0)
            i += np.sum(mask)
        else:
            info[hy, :] = np.amax(roll_30[i:(i + np.sum(mask)), :], axis=0)
            i += np.sum(mask)
    # calculations for entire time series
    sfc = np.mean(info, axis=0) / median

    return sfc


def dh14(flows, datetimes, hydro_years, drainage_area):
    """Q5 exceedance monthly mean flow value normalised by overall
    mean monthly mean flow.

    :Calculation Details:
        Compute the monthly mean flow for the whole daily flow record.
        Calculate the 95th percentile of these mean values. Divide this
        percentile by the mean of these monthly mean values.

    """
    # calculations per month for each year
    mean = np.array(pd.DataFrame(flows, index=datetimes).groupby(
        lambda x: (x.year, x.month)).mean())
    # calculations for entire time series
    perc95 = np.percentile(mean, 95, axis=0)
    sfc = perc95 / np.mean(mean, axis=0)

    return sfc


def dh15(flows, datetimes, hydro_years, drainage_area):
    """Median high-flow pulse duration.

    :Calculation Details:
        Identify for each hydrological year the flow events above the
        75th percentile of the whole daily flow record. Calculate the
        mean duration of these flow events for hydrological year. Take
        the median of these mean values.

    """
    perc75 = np.percentile(flows, 75, axis=0)
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[hy, :] = calc_events_avg_duration(flows[mask, :], perc75,
                                               typ='high')
    # calculations for entire time series
    sfc = np.median(info, axis=0)

    return sfc


def dh16(flows, datetimes, hydro_years, drainage_area):
    """Variability in high-flow pulse duration.

    :Calculation Details:
        Identify for each hydrological year the flow events above the
        75th percentile of the whole daily flow record. Calculate the
        mean duration of these flow events for hydrological year.
        Calculate the standard deviation and the mean of these mean
        values. Multiply the former by 100 and divide the result by the
        latter.

    """
    perc75 = np.percentile(flows, 75, axis=0)
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[hy, :] = calc_events_avg_duration(flows[mask, :], perc75,
                                               typ='high')
    # calculations for entire time series
    sfc = np.std(info, ddof=1, axis=0) * 100 / np.mean(info, axis=0)

    return sfc


def dh17(flows, datetimes, hydro_years, drainage_area):
    """Mean duration of flow events above the median flow for the whole
    record.

    :Calculation Details:
        Identify the flow events across the whole daily flow record
        which are above the median of the whole daily flow record.
        Calculate the mean duration of these flow events.

    """
    median = np.median(flows, axis=0)
    # calculations for entire time series
    sfc = calc_events_avg_duration(flows, median, typ='high')

    return sfc


def dh18(flows, datetimes, hydro_years, drainage_area):
    """Mean duration of flow events above three times the median flow
    for the whole record.

    :Calculation Details:
        Identify the flow events across the whole daily flow record
        which are three times above the median of the whole daily flow
        record. Calculate the mean duration of these flow events.

    """
    median = np.median(flows, axis=0)
    # calculations for entire time series
    sfc = calc_events_avg_duration(flows, 3 * median, typ='high')

    return sfc


def dh19(flows, datetimes, hydro_years, drainage_area):
    """Mean duration of flow events above seven times the median flow
    for the whole record.

    :Calculation Details:
        Identify the flow events across the whole daily flow record
        which are seven times above the median of the whole daily flow
        record. Calculate the mean duration of these flow events.

    """
    median = np.median(flows, axis=0)
    # calculations for entire time series
    sfc = calc_events_avg_duration(flows, 7 * median, typ='high')

    return sfc


def dh20(flows, datetimes, hydro_years, drainage_area):
    """Mean duration of flow events above Q25 for the whole record.

    :Calculation Details:
        Identify the flow events across the whole daily flow record
        which are above the 75th percentile of the whole daily flow
        record. Calculate the mean duration of these flow events.

    """
    perc75 = np.percentile(flows, 75, axis=0)
    # calculations for entire time series
    sfc = calc_events_avg_duration(flows, perc75, typ='high')

    return sfc


def dh21(flows, datetimes, hydro_years, drainage_area):
    """Mean duration of flow events above Q75 for the whole record.

    :Calculation Details:
        Identify the flow events across the whole daily flow record
        which are above the 25th percentile of the whole daily flow
        record. Calculate the mean duration of these flow events.

    """
    perc25 = np.percentile(flows, 25, axis=0)
    # calculations for entire time series
    sfc = calc_events_avg_duration(flows, perc25, typ='high')

    return sfc


# DH22 to DH24 not available
