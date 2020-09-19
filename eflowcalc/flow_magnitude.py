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
from .tools import calc_bfi, rolling_window, calc_events_avg_volume_above


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# AVERAGE FLOWS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def ma1(flows, datetimes, hydro_years, drainage_area):
    """Mean daily flow for the entire daily flow record.

    :Calculation Details:
        Calculate the mean for the whole daily flow record.

    """
    # calculations for entire time series
    sfc = np.mean(flows, axis=0)

    return sfc


def ma2(flows, datetimes, hydro_years, drainage_area):
    """Median daily flow for the entire daily flow record.

    :Calculation Details:
        Calculate the median for the whole daily flow record.

    """
    # calculations for entire time series
    sfc = np.median(flows, axis=0)

    return sfc


def ma3(flows, datetimes, hydro_years, drainage_area):
    """Mean annual coefficient of variation for daily flows.

    :Calculation Details:
        Calculate the standard deviation and the mean of the daily flows
        for each hydrological year. Divide the values of the former by
        the values of the latter to get the coefficients of variation
        for each hydrological year. Calculate the mean of these
        coefficients.

    """
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[hy, :] = (np.std(flows[mask, :], ddof=1, axis=0)
                       / np.mean(flows[mask, :], axis=0))
    # calculations for entire time series
    sfc = np.mean(info, axis=0) * 100

    return sfc


def ma4(flows, datetimes, hydro_years, drainage_area):
    """Variability in the percentiles on log10 for the entire daily
    flow record.

    :Calculation Details:
        Compute the decimal logarithm of the daily flow values.
        Calculate the 5th, 10th, 15th, ..., 95th percentiles of these
        log-transformed values. Calculate the standard deviation and the
        mean of the 19 percentile values. Multiply the former by 100 and
        divide the result by the latter.

    """
    log_f = np.copy(flows)
    # replace log10(0) by log10(0.01) if necessary
    log_f[log_f == 0.0] = 0.01
    log_f = np.log10(log_f, dtype=np.float64)
    # calculations for entire time series
    perc = np.percentile(log_f, (5, 10, 15, 20, 25, 30, 35, 40, 45, 50,
                                 55, 60, 65, 70, 75, 80, 85, 90, 95), axis=0)
    sfc = np.std(perc, ddof=1, axis=0) * 100 / np.mean(perc, axis=0)

    return sfc


def ma5(flows, datetimes, hydro_years, drainage_area):
    """Skewness of the daily flow values.

    :Calculation Details:
        Calculate the mean flow for the whole daily flow record. Divide
        by the median flow for the whole daily flow record.

    """
    # calculations for entire time series
    sfc = np.mean(flows, axis=0) / np.median(flows, axis=0)

    return sfc


def ma6(flows, datetimes, hydro_years, drainage_area):
    """Ratio between 90th and 10th percentiles for the entire record.

    :Calculation Details:
        Calculate the 90th percentile and the 10th percentile of the
        daily flow record. Divide the former value by the latter value.

    """
    # calculations for entire time series
    sfc = np.percentile(flows, 90, axis=0) / np.percentile(flows, 10, axis=0)

    return sfc


def ma7(flows, datetimes, hydro_years, drainage_area):
    """Ratio between 80th and 20th percentiles for the entire record.

    :Calculation Details:
        Calculate the 80th percentile and the 20th percentile of the
        daily flow record. Divide the former value by the latter value.

    """
    # calculations for entire time series
    sfc = np.percentile(flows, 80, axis=0) / np.percentile(flows, 20, axis=0)

    return sfc


def ma8(flows, datetimes, hydro_years, drainage_area):
    """Ratio between 75th and 25th percentiles for the entire record.

    :Calculation Details:
        Calculate the 75th percentile and the 25th percentile of the
        daily flow record. Divide the former value by the latter value.

    """
    # calculations for entire time series
    sfc = np.percentile(flows, 75, axis=0) / np.percentile(flows, 25, axis=0)

    return sfc


def ma9(flows, datetimes, hydro_years, drainage_area):
    """Spread in 90th-10th percentile range on decimal logarithm
    transformed daily flows.

    :Calculation Details:
        Compute the decimal logarithm of the daily flow values.
        Calculate the 90th percentile and the 10th percentile of these
        log-transformed flow values. Subtract the latter to the former.
        Divide this value by the median of the log-transformed flow
        values.

    """
    log_f = np.copy(flows)
    # replace log10(0) by log10(0.01) if necessary
    log_f[log_f == 0.0] = 0.01
    log_f = np.log10(log_f, dtype=np.float64)
    med_f = np.log10(np.median(flows, axis=0), dtype=np.float64)
    # calculations for entire time series
    sfc = (np.percentile(log_f, 90, axis=0)
           - np.percentile(log_f, 10, axis=0)) / med_f

    return sfc


def ma10(flows, datetimes, hydro_years, drainage_area):
    """Spread in 80th-20th percentile range on decimal logarithm
    transformed daily flows.

    :Calculation Details:
        Compute the decimal logarithm of the daily flow values.
        Calculate the 80th percentile and the 20th percentile of these
        log-transformed flow values. Subtract the latter to the former.
        Divide this value by the median of the log-transformed flow
        values.

    """
    log_f = np.copy(flows)
    # replace log10(0) by log10(0.01) if necessary
    log_f[log_f == 0.0] = 0.01
    log_f = np.log10(log_f, dtype=np.float64)
    med_f = np.log10(np.median(flows, axis=0), dtype=np.float64)
    # calculations for entire time series
    sfc = (np.percentile(log_f, 80, axis=0)
           - np.percentile(log_f, 20, axis=0)) / med_f

    return sfc


def ma11(flows, datetimes, hydro_years, drainage_area):
    """Spread in 75th-25th percentile range on decimal logarithm
    transformed daily flows.

    :Calculation Details:
        Compute the decimal logarithm of the daily flow values.
        Calculate the 75th percentile and the 25th percentile of these
        log-transformed flow values. Subtract the latter to the former.
        Divide this value by the median of the log-transformed flow
        values.

    """
    log_f = np.copy(flows)
    # replace log10(0) by log10(0.01) if necessary
    log_f[log_f == 0.0] = 0.01
    log_f = np.log10(log_f, dtype=np.float64)
    med_f = np.log10(np.median(flows, axis=0), dtype=np.float64)
    # calculations for entire time series
    sfc = (np.percentile(log_f, 75, axis=0)
           - np.percentile(log_f, 25, axis=0)) / med_f

    return sfc


def ma12(flows, datetimes, hydro_years, drainage_area):
    """Mean of January flows.

    :Calculation Details:
        Compute the mean daily flow values for each month in the record.
        Take the mean value for January.

    """
    # calculations per month
    info = np.array(pd.DataFrame(flows, index=datetimes).groupby(
        lambda x: x.month).mean())
    # calculations for entire time series
    sfc = info[0, :]

    return sfc


def ma13(flows, datetimes, hydro_years, drainage_area):
    """Mean of February flows.

    :Calculation Details:
        Compute the mean daily flow values for each month in the record.
        Take the mean value for February.

    """
    # calculations per month
    info = np.array(pd.DataFrame(flows, index=datetimes).groupby(
        lambda x: x.month).mean())
    # calculations for entire time series
    sfc = info[1, :]

    return sfc


def ma14(flows, datetimes, hydro_years, drainage_area):
    """Mean of March flows.

    :Calculation Details:
        Compute the mean daily flow values for each month in the record.
        Take the mean value for March.

    """
    # calculations per month
    info = np.array(pd.DataFrame(flows, index=datetimes).groupby(
        lambda x: x.month).mean())
    # calculations for entire time series
    sfc = info[2, :]

    return sfc


def ma15(flows, datetimes, hydro_years, drainage_area):
    """Mean of April flows.

    :Calculation Details:
        Compute the mean daily flow values for each month in the record.
        Take the mean value for April.

    """
    # calculations per month
    info = np.array(pd.DataFrame(flows, index=datetimes).groupby(
        lambda x: x.month).mean())
    # calculations for entire time series
    sfc = info[3, :]

    return sfc


def ma16(flows, datetimes, hydro_years, drainage_area):
    """Mean of May flows.

    :Calculation Details:
        Compute the mean daily flow values for each month in the record.
        Take the mean value for May.

    """
    # calculations per month
    info = np.array(pd.DataFrame(flows, index=datetimes).groupby(
        lambda x: x.month).mean())
    # calculations for entire time series
    sfc = info[4, :]

    return sfc


def ma17(flows, datetimes, hydro_years, drainage_area):
    """Mean of June flows.

    :Calculation Details:
        Compute the mean daily flow values for each month in the record.
        Take the mean value for June.

    """
    # calculations per month
    info = np.array(pd.DataFrame(flows, index=datetimes).groupby(
        lambda x: x.month).mean())
    # calculations for entire time series
    sfc = info[5, :]

    return sfc


def ma18(flows, datetimes, hydro_years, drainage_area):
    """Mean of July flows.

    :Calculation Details:
        Compute the mean daily flow values for each month in the record.
        Take the mean value for July.

    """
    # calculations per month
    info = np.array(pd.DataFrame(flows, index=datetimes).groupby(
        lambda x: x.month).mean())
    # calculations for entire time series
    sfc = info[6, :]

    return sfc


def ma19(flows, datetimes, hydro_years, drainage_area):
    """Mean of August flows.

    :Calculation Details:
        Compute the mean daily flow values for each month in the record.
        Take the mean value for August.

    """
    # calculations per month
    info = np.array(pd.DataFrame(flows, index=datetimes).groupby(
        lambda x: x.month).mean())
    # calculations for entire time series
    sfc = info[7, :]

    return sfc


def ma20(flows, datetimes, hydro_years, drainage_area):
    """Mean of September flows.

    :Calculation Details:
        Compute the mean daily flow values for each month in the record.
        Take the mean value for September.

    """
    # calculations per month
    info = np.array(pd.DataFrame(flows, index=datetimes).groupby(
        lambda x: x.month).mean())
    # calculations for entire time series
    sfc = info[8, :]

    return sfc


def ma21(flows, datetimes, hydro_years, drainage_area):
    """Mean of October flows.

    :Calculation Details:
        Compute the mean daily flow values for each month in the record.
        Take the mean value for October.

    """
    # calculations per month
    info = np.array(pd.DataFrame(flows, index=datetimes).groupby(
        lambda x: x.month).mean())
    # calculations for entire time series
    sfc = info[9, :]

    return sfc


def ma22(flows, datetimes, hydro_years, drainage_area):
    """Mean of November flows.

    :Calculation Details:
        Compute the mean daily flow values for each month in the record.
        Take the mean value for November.

    """
    # calculations per month
    info = np.array(pd.DataFrame(flows, index=datetimes).groupby(
        lambda x: x.month).mean())
    # calculations for entire time series
    sfc = info[10, :]

    return sfc


def ma23(flows, datetimes, hydro_years, drainage_area):
    """Mean of December flows.

    :Calculation Details:
        Compute the mean daily flow values for each month in the record.
        Take the mean value for December.

    """
    # calculations per month
    info = np.array(pd.DataFrame(flows, index=datetimes).groupby(
        lambda x: x.month).mean())
    # calculations for entire time series
    sfc = info[11, :]

    return sfc


def ma24(flows, datetimes, hydro_years, drainage_area):
    """Variability in January flows.

    :Calculation Details:
        Compute the standard deviation and the mean of the daily flow
        values for each month in the record. Multiply the former values
        by 100 and divide the results by the latter values to obtain
        the coefficients of variations for each month. Take the mean
        value for January.

    """
    # calculations per month for each year
    std_ = pd.DataFrame(flows, index=datetimes).groupby(
        lambda x: (x.year, x.month)).std()
    mean_ = pd.DataFrame(flows, index=datetimes).groupby(
        lambda x: (x.year, x.month)).mean()
    # calculations per month for the entire series
    cv_ = (std_ / mean_).groupby(lambda x: x[1]).mean()
    info = np.array(cv_)
    # calculations for entire time series
    sfc = info[0, :] * 100

    return sfc


def ma25(flows, datetimes, hydro_years, drainage_area):
    """Variability in February flows.

    :Calculation Details:
        Compute the standard deviation and the mean of the daily flow
        values for each month in the record. Multiply the former values
        by 100 and divide the results by the latter values to obtain
        the coefficients of variations for each month. Take the mean
        value for February.

    """
    # calculations per month for each year
    std_ = pd.DataFrame(flows, index=datetimes).groupby(
        lambda x: (x.year, x.month)).std()
    mean_ = pd.DataFrame(flows, index=datetimes).groupby(
        lambda x: (x.year, x.month)).mean()
    # calculations per month for the entire series
    cv_ = (std_ / mean_).groupby(lambda x: x[1]).mean()
    info = np.array(cv_)
    # calculations for entire time series
    sfc = info[1, :] * 100

    return sfc


def ma26(flows, datetimes, hydro_years, drainage_area):
    """Variability in March flows.

    :Calculation Details:
        Compute the standard deviation and the mean of the daily flow
        values for each month in the record. Multiply the former values
        by 100 and divide the results by the latter values to obtain
        the coefficients of variations for each month. Take the mean
        value for March.

    """
    # calculations per month for each year
    std_ = pd.DataFrame(flows, index=datetimes).groupby(
        lambda x: (x.year, x.month)).std()
    mean_ = pd.DataFrame(flows, index=datetimes).groupby(
        lambda x: (x.year, x.month)).mean()
    # calculations per month for the entire series
    cv_ = (std_ / mean_).groupby(lambda x: x[1]).mean()
    info = np.array(cv_)
    # calculations for entire time series
    sfc = info[2, :] * 100

    return sfc


def ma27(flows, datetimes, hydro_years, drainage_area):
    """Variability in April flows.

    :Calculation Details:
        Compute the standard deviation and the mean of the daily flow
        values for each month in the record. Multiply the former values
        by 100 and divide the results by the latter values to obtain
        the coefficients of variations for each month. Take the mean
        value for April.

    """
    # calculations per month for each year
    std_ = pd.DataFrame(flows, index=datetimes).groupby(
        lambda x: (x.year, x.month)).std()
    mean_ = pd.DataFrame(flows, index=datetimes).groupby(
        lambda x: (x.year, x.month)).mean()
    # calculations per month for the entire series
    cv_ = (std_ / mean_).groupby(lambda x: x[1]).mean()
    info = np.array(cv_)
    # calculations for entire time series
    sfc = info[3, :] * 100

    return sfc


def ma28(flows, datetimes, hydro_years, drainage_area):
    """Variability in May flows.

    :Calculation Details:
        Compute the standard deviation and the mean of the daily flow
        values for each month in the record. Multiply the former values
        by 100 and divide the results by the latter values to obtain
        the coefficients of variations for each month. Take the mean
        value for May.

    """
    # calculations per month for each year
    std_ = pd.DataFrame(flows, index=datetimes).groupby(
        lambda x: (x.year, x.month)).std()
    mean_ = pd.DataFrame(flows, index=datetimes).groupby(
        lambda x: (x.year, x.month)).mean()
    # calculations per month for the entire series
    cv_ = (std_ / mean_).groupby(lambda x: x[1]).mean()
    info = np.array(cv_)
    # calculations for entire time series
    sfc = info[4, :] * 100

    return sfc


def ma29(flows, datetimes, hydro_years, drainage_area):
    """Variability in June flows.

    :Calculation Details:
        Compute the standard deviation and the mean of the daily flow
        values for each month in the record. Multiply the former values
        by 100 and divide the results by the latter values to obtain
        the coefficients of variations for each month. Take the mean
        value for June.

    """
    # calculations per month for each year
    std_ = pd.DataFrame(flows, index=datetimes).groupby(
        lambda x: (x.year, x.month)).std()
    mean_ = pd.DataFrame(flows, index=datetimes).groupby(
        lambda x: (x.year, x.month)).mean()
    # calculations per month for the entire series
    cv_ = (std_ / mean_).groupby(lambda x: x[1]).mean()
    info = np.array(cv_)
    # calculations for entire time series
    sfc = info[5, :] * 100

    return sfc


def ma30(flows, datetimes, hydro_years, drainage_area):
    """Variability in July flows.

    :Calculation Details:
        Compute the standard deviation and the mean of the daily flow
        values for each month in the record. Multiply the former values
        by 100 and divide the results by the latter values to obtain
        the coefficients of variations for each month. Take the mean
        value for July.

    """
    # calculations per month for each year
    std_ = pd.DataFrame(flows, index=datetimes).groupby(
        lambda x: (x.year, x.month)).std()
    mean_ = pd.DataFrame(flows, index=datetimes).groupby(
        lambda x: (x.year, x.month)).mean()
    # calculations per month for the entire series
    cv_ = (std_ / mean_).groupby(lambda x: x[1]).mean()
    info = np.array(cv_)
    # calculations for entire time series
    sfc = info[6, :] * 100

    return sfc


def ma31(flows, datetimes, hydro_years, drainage_area):
    """Variability in August flows.

    :Calculation Details:
        Compute the standard deviation and the mean of the daily flow
        values for each month in the record. Multiply the former values
        by 100 and divide the results by the latter values to obtain
        the coefficients of variations for each month. Take the mean
        value for August.

    """
    # calculations per month for each year
    std_ = pd.DataFrame(flows, index=datetimes).groupby(
        lambda x: (x.year, x.month)).std()
    mean_ = pd.DataFrame(flows, index=datetimes).groupby(
        lambda x: (x.year, x.month)).mean()
    # calculations per month for the entire series
    cv_ = (std_ / mean_).groupby(lambda x: x[1]).mean()
    info = np.array(cv_)
    # calculations for entire time series
    sfc = info[7, :] * 100

    return sfc


def ma32(flows, datetimes, hydro_years, drainage_area):
    """Variability in September flows.

    :Calculation Details:
        Compute the standard deviation and the mean of the daily flow
        values for each month in the record. Multiply the former values
        by 100 and divide the results by the latter values to obtain
        the coefficients of variations for each month. Take the mean
        value for September.

    """
    # calculations per month for each year
    std_ = pd.DataFrame(flows, index=datetimes).groupby(
        lambda x: (x.year, x.month)).std()
    mean_ = pd.DataFrame(flows, index=datetimes).groupby(
        lambda x: (x.year, x.month)).mean()
    # calculations per month for the entire series
    cv_ = (std_ / mean_).groupby(lambda x: x[1]).mean()
    info = np.array(cv_)
    # calculations for entire time series
    sfc = info[8, :] * 100

    return sfc


def ma33(flows, datetimes, hydro_years, drainage_area):
    """Variability in October flows.

    :Calculation Details:
        Compute the standard deviation and the mean of the daily flow
        values for each month in the record. Multiply the former values
        by 100 and divide the results by the latter values to obtain
        the coefficients of variations for each month. Take the mean
        value for October.

    """
    # calculations per month for each year
    std_ = pd.DataFrame(flows, index=datetimes).groupby(
        lambda x: (x.year, x.month)).std()
    mean_ = pd.DataFrame(flows, index=datetimes).groupby(
        lambda x: (x.year, x.month)).mean()
    # calculations per month for the entire series
    cv_ = (std_ / mean_).groupby(lambda x: x[1]).mean()
    info = np.array(cv_)
    # calculations for entire time series
    sfc = info[9, :] * 100

    return sfc


def ma34(flows, datetimes, hydro_years, drainage_area):
    """Variability in November flows.

    :Calculation Details:
        Compute the standard deviation and the mean of the daily flow
        values for each month in the record. Multiply the former values
        by 100 and divide the results by the latter values to obtain
        the coefficients of variations for each month. Take the mean
        value for November.

    """
    # calculations per month for each year
    std_ = pd.DataFrame(flows, index=datetimes).groupby(
        lambda x: (x.year, x.month)).std()
    mean_ = pd.DataFrame(flows, index=datetimes).groupby(
        lambda x: (x.year, x.month)).mean()
    # calculations per month for the entire series
    cv_ = (std_ / mean_).groupby(lambda x: x[1]).mean()
    info = np.array(cv_)
    # calculations for entire time series
    sfc = info[10, :] * 100

    return sfc


def ma35(flows, datetimes, hydro_years, drainage_area):
    """Variability in December flows.

    :Calculation Details:
        Compute the standard deviation and the mean of the daily flow
        values for each month in the record. Multiply the former values
        by 100 and divide the results by the latter values to obtain
        the coefficients of variations for each month. Take the mean
        value for December.

    """
    # calculations per month for each year
    std_ = pd.DataFrame(flows, index=datetimes).groupby(
        lambda x: (x.year, x.month)).std()
    mean_ = pd.DataFrame(flows, index=datetimes).groupby(
        lambda x: (x.year, x.month)).mean()
    # calculations per month for the entire series
    cv_ = (std_ / mean_).groupby(lambda x: x[1]).mean()
    info = np.array(cv_)
    # calculations for entire time series
    sfc = info[11, :] * 100

    return sfc


def ma36(flows, datetimes, hydro_years, drainage_area):
    """Skewness in mean monthly flow using min and max.

    :Calculation Details:
        Compute the mean of the daily flow values for each month in the
        record. Calculate the minimum, maximum, and median of these mean
        values. Subtract this maximum to this minimum, and divide the
        result by this median.

    """
    # calculations per month for each year
    mean_ = np.array(pd.DataFrame(flows, index=datetimes).groupby(
        lambda x: (x.year, x.month)).mean())
    # calculations for entire time series
    sfc = (np.amax(mean_, axis=0)
           - np.amin(mean_, axis=0)) / np.median(mean_, axis=0)

    return sfc


def ma37(flows, datetimes, hydro_years, drainage_area):
    """Skewness in mean monthly flow using 25th and 75th percentiles.

    :Calculation Details:
        Compute the mean of the daily flow values for each month in the
        record. Calculate the 75th percentile, 25th percentile, and
        median of these mean values. Subtract this 75th percentile to
        this 25th percentile, and divide the result by this median.

    """
    # calculations per month for each year
    mean_ = np.array(pd.DataFrame(flows, index=datetimes).groupby(
        lambda x: (x.year, x.month)).mean())
    # calculations for entire time series
    sfc = (np.percentile(mean_, 75, axis=0)
           - np.percentile(mean_, 25, axis=0)) / np.median(mean_, axis=0)

    return sfc


def ma38(flows, datetimes, hydro_years, drainage_area):
    """Skewness in mean monthly flow using 10th and 90th percentiles.

    :Calculation Details:
        Compute the mean of the daily flow values for each month in the
        record. Calculate the 90th percentile, 10th percentile, and
        median of these mean values. Subtract this 90th percentile to
        this 10th percentile, and divide the result by this median.

    """
    # calculations per month for each year
    mean_ = np.array(pd.DataFrame(flows, index=datetimes).groupby(
        lambda x: (x.year, x.month)).mean())
    # calculations for entire time series
    sfc = (np.percentile(mean_, 90, axis=0)
           - np.percentile(mean_, 10, axis=0)) / np.median(mean_, axis=0)

    return sfc


def ma39(flows, datetimes, hydro_years, drainage_area):
    """Variability in mean monthly flow.

    :Calculation Details:
        Compute the mean of the daily flow values for each month in the
        record. Compute the standard deviation and the mean of these
        mean values. Multiply the former by 100, and divide the result
        by the latter.

    """
    # calculations per month for each year
    mean_ = np.array(pd.DataFrame(flows, index=datetimes).groupby(
        lambda x: (x.year, x.month)).mean())
    # calculations for entire time series
    sfc = np.std(mean_, ddof=1, axis=0) * 100 / np.mean(mean_, axis=0)

    return sfc


def ma40(flows, datetimes, hydro_years, drainage_area):
    """Skewness in mean monthly flow.

    :Calculation Details:
        Compute the mean of the daily flow values for each month in the
        record. Compute the mean and the median of these mean values.
        Subtract this median to this mean, and divide the result by this
        median.

    """
    # calculations per month for each year
    mean_ = np.array(pd.DataFrame(flows, index=datetimes).groupby(
        lambda x: (x.year, x.month)).mean())
    # calculations for entire time series
    sfc = (np.mean(mean_, axis=0)
           - np.median(mean_, axis=0)) / np.median(mean_, axis=0)

    return sfc


def ma41(flows, datetimes, hydro_years, drainage_area):
    """Mean annual daily flow normalised by drainage area.

    :Calculation Details:
        Calculate the mean daily flow value for each hydrological year.
        Divide these mean values by the catchment drainage area.
        Calculate the mean of these normalised values.

    """
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[hy, :] = np.mean(flows[mask, :], axis=0) / drainage_area
    # calculations for entire time series
    sfc = np.mean(info, axis=0)

    return sfc


def ma42(flows, datetimes, hydro_years, drainage_area):
    """Skewness in mean annual flow using min and max.

    :Calculation Details:
        Compute the mean of the daily flow values for each hydrological
        year. Calculate the minimum, maximum, and median of these mean
        values. Subtract this maximum to this minimum, and divide the
        result by this median.

    """
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[hy, :] = np.mean(flows[mask, :], axis=0)
    # calculations for entire time series
    sfc = (np.amax(info, axis=0)
           - np.amin(info, axis=0)) / np.median(info, axis=0)

    return sfc


def ma43(flows, datetimes, hydro_years, drainage_area):
    """Skewness in mean annual flow using 25th and 75th percentiles.

    :Calculation Details:
        Compute the mean of the daily flow values for each hydrological
        year. Calculate the 75th percentile, 25th percentile, and median
        of these mean values. Subtract this 75th percentile to this 25th
        percentile, and divide the result by this median.

    """
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[hy, :] = np.mean(flows[mask, :], axis=0)
    # calculations for entire time series
    sfc = (np.percentile(info, 75, axis=0)
           - np.percentile(info, 25, axis=0)) / np.median(info, axis=0)

    return sfc


def ma44(flows, datetimes, hydro_years, drainage_area):
    """Skewness in mean annual flow using 10th and 90th percentiles.

    :Calculation Details:
        Compute the mean of the daily flow values for each hydrological
        year. Calculate the 90th percentile, 10th percentile, and median
        of these mean values. Subtract this 90th percentile to this 10th
        percentile, and divide the result by this median.

    """
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[hy, :] = np.mean(flows[mask, :], axis=0)
    # calculations for entire time series
    sfc = (np.percentile(info, 90, axis=0)
           - np.percentile(info, 10, axis=0)) / np.median(info, axis=0)

    return sfc


# MA45 - Skewness in mean annual flow
def ma45(flows, datetimes, hydro_years, drainage_area):
    """Skewness in mean annual flow.

    :Calculation Details:
        Compute the mean of the daily flow values for each hydrological
        year. Calculate the mean, and median of these mean values.
        Subtract this mean to this median, and divide the result by this
        median.

    """
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[hy, :] = np.mean(flows[mask, :], axis=0)
    # calculations for entire time series
    sfc = (np.mean(info, axis=0)
           - np.median(info, axis=0)) / np.median(info, axis=0)

    return sfc


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOW FLOWS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def ml1(flows, datetimes, hydro_years, drainage_area):
    """Mean minimum flow for January.

    :Calculation Details:
        Compute the minimum flow values for each month in the daily
        record. Take the mean of these minima for January.

    """
    # calculations per month for each year
    min_ = pd.DataFrame(flows, index=datetimes).groupby(
        lambda x: (x.year, x.month)).min()
    info = np.array(min_.groupby(lambda x: x[1]).mean())
    # calculations for entire time series
    sfc = info[0, :]

    return sfc


def ml2(flows, datetimes, hydro_years, drainage_area):
    """Mean minimum flow for February.

    :Calculation Details:
        Compute the minimum flow values for each month in the daily
        record. Take the mean of these minima for February.

    """
    # calculations per month for each year
    min_ = pd.DataFrame(flows, index=datetimes).groupby(
        lambda x: (x.year, x.month)).min()
    info = np.array(min_.groupby(lambda x: x[1]).mean())
    # calculations for entire time series
    sfc = info[1, :]

    return sfc


def ml3(flows, datetimes, hydro_years, drainage_area):
    """Mean minimum flow for March.

    :Calculation Details:
        Compute the minimum flow values for each month in the daily
        record. Take the mean of these minima for March.

    """
    # calculations per month for each year
    min_ = pd.DataFrame(flows, index=datetimes).groupby(
        lambda x: (x.year, x.month)).min()
    info = np.array(min_.groupby(lambda x: x[1]).mean())
    # calculations for entire time series
    sfc = info[2, :]

    return sfc


def ml4(flows, datetimes, hydro_years, drainage_area):
    """Mean minimum flow for April.

    :Calculation Details:
        Compute the minimum flow values for each month in the daily
        record. Take the mean of these minima for April.

    """
    # calculations per month for each year
    min_ = pd.DataFrame(flows, index=datetimes).groupby(
        lambda x: (x.year, x.month)).min()
    info = np.array(min_.groupby(lambda x: x[1]).mean())
    # calculations for entire time series
    sfc = info[3, :]

    return sfc


def ml5(flows, datetimes, hydro_years, drainage_area):
    """Mean minimum flow for May.

    :Calculation Details:
        Compute the minimum flow values for each month in the daily
        record. Take the mean of these minima for May.

    """
    # calculations per month for each year
    min_ = pd.DataFrame(flows, index=datetimes).groupby(
        lambda x: (x.year, x.month)).min()
    info = np.array(min_.groupby(lambda x: x[1]).mean())
    # calculations for entire time series
    sfc = info[4, :]

    return sfc


def ml6(flows, datetimes, hydro_years, drainage_area):
    """Mean minimum flow for June.

    :Calculation Details:
        Compute the minimum flow values for each month in the daily
        record. Take the mean of these minima for June.

    """
    # calculations per month for each year
    min_ = pd.DataFrame(flows, index=datetimes).groupby(
        lambda x: (x.year, x.month)).min()
    info = np.array(min_.groupby(lambda x: x[1]).mean())
    # calculations for entire time series
    sfc = info[5, :]

    return sfc


def ml7(flows, datetimes, hydro_years, drainage_area):
    """Mean minimum flow for July.

    :Calculation Details:
        Compute the minimum flow values for each month in the daily
        record. Take the mean of these minima for July.

    """
    # calculations per month for each year
    min_ = pd.DataFrame(flows, index=datetimes).groupby(
        lambda x: (x.year, x.month)).min()
    info = np.array(min_.groupby(lambda x: x[1]).mean())
    # calculations for entire time series
    sfc = info[6, :]

    return sfc


def ml8(flows, datetimes, hydro_years, drainage_area):
    """Mean minimum flow for August.

    :Calculation Details:
        Compute the minimum flow values for each month in the daily
        record. Take the mean of these minima for August.

    """
    # calculations per month for each year
    min_ = pd.DataFrame(flows, index=datetimes).groupby(
        lambda x: (x.year, x.month)).min()
    info = np.array(min_.groupby(lambda x: x[1]).mean())
    # calculations for entire time series
    sfc = info[7, :]

    return sfc


def ml9(flows, datetimes, hydro_years, drainage_area):
    """Mean minimum flow for September.

    :Calculation Details:
        Compute the minimum flow values for each month in the daily
        record. Take the mean of these minima for September.

    """
    # calculations per month for each year
    min_ = pd.DataFrame(flows, index=datetimes).groupby(
        lambda x: (x.year, x.month)).min()
    info = np.array(min_.groupby(lambda x: x[1]).mean())
    # calculations for entire time series
    sfc = info[8, :]

    return sfc


def ml10(flows, datetimes, hydro_years, drainage_area):
    """Mean minimum flow for October.

    :Calculation Details:
        Compute the minimum flow values for each month in the daily
        record. Take the mean of these minima for October.

    """
    # calculations per month for each year
    min_ = pd.DataFrame(flows, index=datetimes).groupby(
        lambda x: (x.year, x.month)).min()
    info = np.array(min_.groupby(lambda x: x[1]).mean())
    # calculations for entire time series
    sfc = info[9, :]

    return sfc


def ml11(flows, datetimes, hydro_years, drainage_area):
    """Mean minimum flow for November.

    :Calculation Details:
        Compute the minimum flow values for each month in the daily
        record. Take the mean of these minima for November.

    """
    # calculations per month for each year
    min_ = pd.DataFrame(flows, index=datetimes).groupby(
        lambda x: (x.year, x.month)).min()
    info = np.array(min_.groupby(lambda x: x[1]).mean())
    # calculations for entire time series
    sfc = info[10, :]

    return sfc


def ml12(flows, datetimes, hydro_years, drainage_area):
    """Mean minimum flow for December.

    :Calculation Details:
        Compute the minimum flow values for each month in the daily
        record. Take the mean of these minima for December.

    """
    # calculations per month for each year
    min_ = pd.DataFrame(flows, index=datetimes).groupby(
        lambda x: (x.year, x.month)).min()
    info = np.array(min_.groupby(lambda x: x[1]).mean())
    # calculations for entire time series
    sfc = info[11, :]

    return sfc


def ml13(flows, datetimes, hydro_years, drainage_area):
    """Variability in minimum monthly flow.

    :Calculation Details:
        Compute the minimum flow values for each month in the daily
        record. Calculate the standard deviation and the mean of these
        minimum values. Multiply the former by 100 and divide the result
        by the latter.

    """
    # calculations per month for each year
    min_ = np.array(pd.DataFrame(flows, index=datetimes).groupby(
        lambda x: (x.year, x.month)).min())
    # calculations for entire time series
    sfc = np.std(min_, ddof=1, axis=0) * 100 / np.mean(min_, axis=0)

    return sfc


def ml14(flows, datetimes, hydro_years, drainage_area):
    """Mean of minimum annual flow normalised by their annual median flow.

    :Calculation Details:
        Compute the minimum and median flow values for each hydrological
        year in the record, and divide the former by the latter.
        Calculate the mean of these ratios.

    """
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[hy, :] = (np.amin(flows[mask, :], axis=0)
                       / np.median(flows[mask, :], axis=0))
    # calculations for entire time series
    sfc = np.mean(info, axis=0)

    return sfc


def ml15(flows, datetimes, hydro_years, drainage_area):
    """Mean of minimum annual flow normalised by their annual mean flow.

    :Calculation Details:
        Compute the minimum and mean flow values for each hydrological
        year in the record, and divide the former by the latter.
        Calculate the mean of these ratios.

    """
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[hy, :] = (np.amin(flows[mask, :], axis=0)
                       / np.mean(flows[mask, :], axis=0))
    # calculations for entire time series
    sfc = np.mean(info, axis=0)

    return sfc


def ml16(flows, datetimes, hydro_years, drainage_area):
    """Median of minimum annual flow normalised by their annual median flow.

    :Calculation Details:
        Compute the minimum and median flow values for each hydrological
        year in the record, and divide the former by the latter.
        Calculate the median of these ratios.

    """
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[hy, :] = (np.amin(flows[mask, :], axis=0)
                       / np.median(flows[mask, :], axis=0))
    # calculations for entire time series
    sfc = np.median(info, axis=0)

    return sfc


def ml17(flows, datetimes, hydro_years, drainage_area):
    """Base flow 1.

    :Calculation Details:
        Compute the 7-day rolling mean flows for each hydrological year,
        take the minimum rolling mean value, and divide this minimum by
        the mean flow for the whole daily record to obtain the base flow
        index for the hydrological year. Compute the mean of the base
        flow indices.

    """
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[hy, :] = calc_bfi(flows[mask, :])
    # calculations for entire time series
    sfc = np.mean(info, axis=0)

    return sfc


def ml18(flows, datetimes, hydro_years, drainage_area):
    """Variability in base flow 1.

    :Calculation Details:
        Compute the 7-day rolling mean flows for each hydrological year,
        take the minimum rolling mean value, and divide this minimum by
        the mean flow for the whole daily record to obtain the base flow
        index for the hydrological year. Compute the standard deviation
        and the mean of the base flow indices. Multiply the former by
        100 and divide the result by the latter.

    """
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[hy, :] = calc_bfi(flows[mask, :])
    # calculations for entire time series
    sfc = np.std(info, ddof=1, axis=0) * 100 / np.mean(info, axis=0)

    return sfc


def ml19(flows, datetimes, hydro_years, drainage_area):
    """Base flow 2.

    :Calculation Details:
        Calculate the minimum and mean flow values for each hydrological
        year, and divide the former by the latter. Compute the mean of
        these ratios and multiply by 100.

    """
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[hy, :] = (np.amin(flows[mask, :], axis=0)
                       / np.mean(flows[mask, :], axis=0))
    # calculations for entire time series
    sfc = np.mean(info, axis=0) * 100

    return sfc


def ml20(flows, datetimes, hydro_years, drainage_area):
    """Base flow 3.

    :Calculation Details:
        Divide the whole record into 5-day blocks. Calculate the minimum
        daily flow values in each block. For each block, base flow is
        taken as this minimum value if 90% of this minimum is less than
        the minimum for the blocks on either side of the block.
        Otherwise, base flow is linearly interpolated from the blocks
        that do meet the condition previously mentioned. Calculate the
        sum of five times these base flow values, and the sum of all
        daily flows. Divide the former by the latter.

    """
    # calculations for entire time series

    # divide the flow record into 5-day blocks
    nb_blocks = flows.shape[0] // 5
    flows_5day_blocks = flows[0:5 * nb_blocks, :].reshape(nb_blocks, 5,
                                                          flows.shape[1])

    # compute the minimum in each block
    blocks_min = np.amin(flows_5day_blocks, axis=1)
    # compute the 3-block rolling minimum as a mean to compare blocks
    # with their preceding and succeeding neighbouring blocks
    rolling_blocks_min = np.amin(rolling_window(blocks_min, 3), axis=1)

    # compute base flow values
    base_flow = np.zeros(blocks_min.shape, dtype=np.float64)
    base_flow[:] = np.nan

    # take base flow value as block minimum if 90% of it is less than
    # the minimum of the blocks on either side of the block (using the
    # 3-block rolling minimum to do so)
    check = blocks_min[1:-1, :] * 0.90 < rolling_blocks_min
    base_flow[1:-1, :][check] = blocks_min[1:-1, :][check]
    base_flow[0, :] = blocks_min[0, :]
    base_flow[-1, :] = blocks_min[-1, :]

    # for those blocks not meeting the condition above, infer base flow
    # value through linear interpolation
    base_flow_ = pd.DataFrame(base_flow).interpolate(
        method='linear', axis=0).values

    sfc = np.sum(base_flow_ * 5, axis=0) / np.sum(flows, axis=0)

    return sfc


def ml21(flows, datetimes, hydro_years, drainage_area):
    """Variability in annual minimum daily flow.

    *Note:* this streamflow characteristic is equivalent to `dl6`.

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


# ML22 - Mean in annual minimum daily flow
def ml22(flows, datetimes, hydro_years, drainage_area):
    """Mean annual minimum daily flow.

    :Calculation Details:
        Take the minimum flow for each hydrological year. Calculate
        the mean of these minimum values.

    """
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[hy, :] = np.amin(flows[mask, :], axis=0) / drainage_area
    # calculations for entire time series
    sfc = np.mean(info, axis=0)

    return sfc


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# HIGH FLOWS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def mh1(flows, datetimes, hydro_years, drainage_area):
    """Mean maximum flow for January.

    :Calculation Details:
        Compute the maximum flow values for each month in the daily
        record. Take the mean of these maxima for January.

    """
    # calculations per month for each year
    max_ = pd.DataFrame(flows, index=datetimes).groupby(
        lambda x: (x.year, x.month)).max()
    info = np.array(max_.groupby(lambda x: x[1]).mean())
    # calculations for entire time series
    sfc = info[0, :]

    return sfc


def mh2(flows, datetimes, hydro_years, drainage_area):
    """Mean maximum flow for February.

    :Calculation Details:
        Compute the maximum flow values for each month in the daily
        record. Take the mean of these maxima for February.

    """
    # calculations per month for each year
    max_ = pd.DataFrame(flows, index=datetimes).groupby(
        lambda x: (x.year, x.month)).max()
    info = np.array(max_.groupby(lambda x: x[1]).mean())
    # calculations for entire time series
    sfc = info[1, :]

    return sfc


def mh3(flows, datetimes, hydro_years, drainage_area):
    """Mean maximum flow for March.

    :Calculation Details:
        Compute the maximum flow values for each month in the daily
        record. Take the mean of these maxima for March.

    """
    # calculations per month for each year
    max_ = pd.DataFrame(flows, index=datetimes).groupby(
        lambda x: (x.year, x.month)).max()
    info = np.array(max_.groupby(lambda x: x[1]).mean())
    # calculations for entire time series
    sfc = info[2, :]

    return sfc


def mh4(flows, datetimes, hydro_years, drainage_area):
    """Mean maximum flow for April.

    :Calculation Details:
        Compute the maximum flow values for each month in the daily
        record. Take the mean of these maxima for April.

    """
    # calculations per month for each year
    max_ = pd.DataFrame(flows, index=datetimes).groupby(
        lambda x: (x.year, x.month)).max()
    info = np.array(max_.groupby(lambda x: x[1]).mean())
    # calculations for entire time series
    sfc = info[3, :]

    return sfc


def mh5(flows, datetimes, hydro_years, drainage_area):
    """Mean maximum flow for May.

    :Calculation Details:
        Compute the maximum flow values for each month in the daily
        record. Take the mean of these maxima for May.

    """
    # calculations per month for each year
    max_ = pd.DataFrame(flows, index=datetimes).groupby(
        lambda x: (x.year, x.month)).max()
    info = np.array(max_.groupby(lambda x: x[1]).mean())
    # calculations for entire time series
    sfc = info[4, :]

    return sfc


def mh6(flows, datetimes, hydro_years, drainage_area):
    """Mean maximum flow for June.

    :Calculation Details:
        Compute the maximum flow values for each month in the daily
        record. Take the mean of these maxima for June.

    """
    # calculations per month for each year
    max_ = pd.DataFrame(flows, index=datetimes).groupby(
        lambda x: (x.year, x.month)).max()
    info = np.array(max_.groupby(lambda x: x[1]).mean())
    # calculations for entire time series
    sfc = info[5, :]

    return sfc


def mh7(flows, datetimes, hydro_years, drainage_area):
    """Mean maximum flow for July.

    :Calculation Details:
        Compute the maximum flow values for each month in the daily
        record. Take the mean of these maxima for July.

    """
    # calculations per month for each year
    max_ = pd.DataFrame(flows, index=datetimes).groupby(
        lambda x: (x.year, x.month)).max()
    info = np.array(max_.groupby(lambda x: x[1]).mean())
    # calculations for entire time series
    sfc = info[6, :]

    return sfc


def mh8(flows, datetimes, hydro_years, drainage_area):
    """Mean maximum flow for August.

    :Calculation Details:
        Compute the maximum flow values for each month in the daily
        record. Take the mean of these maxima for August.

    """
    # calculations per month for each year
    max_ = pd.DataFrame(flows, index=datetimes).groupby(
        lambda x: (x.year, x.month)).max()
    info = np.array(max_.groupby(lambda x: x[1]).mean())
    # calculations for entire time series
    sfc = info[7, :]

    return sfc


def mh9(flows, datetimes, hydro_years, drainage_area):
    """Mean maximum flow for September.

    :Calculation Details:
        Compute the maximum flow values for each month in the daily
        record. Take the mean of these maxima for September.

    """
    # calculations per month for each year
    max_ = pd.DataFrame(flows, index=datetimes).groupby(
        lambda x: (x.year, x.month)).max()
    info = np.array(max_.groupby(lambda x: x[1]).mean())
    # calculations for entire time series
    sfc = info[8, :]

    return sfc


def mh10(flows, datetimes, hydro_years, drainage_area):
    """Mean maximum flow for October.

    :Calculation Details:
        Compute the maximum flow values for each month in the daily
        record. Take the mean of these maxima for October.

    """
    # calculations per month for each year
    max_ = pd.DataFrame(flows, index=datetimes).groupby(
        lambda x: (x.year, x.month)).max()
    info = np.array(max_.groupby(lambda x: x[1]).mean())
    # calculations for entire time series
    sfc = info[9, :]

    return sfc


def mh11(flows, datetimes, hydro_years, drainage_area):
    """Mean maximum flow for November.

    :Calculation Details:
        Compute the maximum flow values for each month in the daily
        record. Take the mean of these maxima for November.

    """
    # calculations per month for each year
    max_ = pd.DataFrame(flows, index=datetimes).groupby(
        lambda x: (x.year, x.month)).max()
    info = np.array(max_.groupby(lambda x: x[1]).mean())
    # calculations for entire time series
    sfc = info[10, :]

    return sfc


def mh12(flows, datetimes, hydro_years, drainage_area):
    """Mean maximum flow for December.

    :Calculation Details:
        Compute the maximum flow values for each month in the daily
        record. Take the mean of these maxima for December.

    """
    # calculations per month for each year
    max_ = pd.DataFrame(flows, index=datetimes).groupby(
        lambda x: (x.year, x.month)).max()
    info = np.array(max_.groupby(lambda x: x[1]).mean())
    # calculations for entire time series
    sfc = info[11, :]

    return sfc


def mh13(flows, datetimes, hydro_years, drainage_area):
    """Variability in maximum monthly flow.

    :Calculation Details:
        Compute the maximum flow values for each month in the daily
        record. Calculate the standard deviation and the mean of these
        maximum values. Multiply the former by 100 and divide the result
        by the latter.

    """
    # calculations per month for each year
    max_ = np.array(pd.DataFrame(flows, index=datetimes).groupby(
        lambda x: (x.year, x.month)).max())
    # calculations for entire time series
    sfc = np.std(max_, ddof=1, axis=0) * 100 / np.mean(max_, axis=0)

    return sfc


def mh14(flows, datetimes, hydro_years, drainage_area):
    """Median of maximum annual flow normalised by their annual median flow.

    :Calculation Details:
        Compute the maximum and median flow values for each hydrological
        year in the record, and divide the former by the latter.
        Calculate the median of these ratios.

    """
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[hy, :] = (np.amax(flows[mask, :], axis=0)
                       / np.median(flows[mask, :], axis=0))
    # calculations for entire time series
    sfc = np.median(info, axis=0)

    return sfc


def mh15(flows, datetimes, hydro_years, drainage_area):
    """Q1 exceedance value normalised by overall median daily flow.

    :Calculation Details:
        Calculate the 99th percentile of the whole daily flow record.
        Divide this percentile by the median of the whole daily flow
        record.

    """
    # calculations for entire time series
    sfc = np.percentile(flows, 99, axis=0) / np.median(flows, axis=0)

    return sfc


def mh16(flows, datetimes, hydro_years, drainage_area):
    """Q10 exceedance value normalised by overall median daily flow.

    :Calculation Details:
        Calculate the 90th percentile of the whole daily flow record.
        Divide this percentile by the median of the whole daily flow
        record.

    """
    # calculations for entire time series
    sfc = np.percentile(flows, 90, axis=0) / np.median(flows, axis=0)

    return sfc


def mh17(flows, datetimes, hydro_years, drainage_area):
    """Q25 exceedance value normalised by overall median daily flow.

    :Calculation Details:
        Calculate the 75th percentile of the whole daily flow record.
        Divide this percentile by the median of the whole daily flow
        record.

    """
    # calculations for entire time series
    sfc = np.percentile(flows, 75, axis=0) / np.median(flows, axis=0)

    return sfc


def mh18(flows, datetimes, hydro_years, drainage_area):
    """Variability in log-transformed annual maximum daily flow.

    :Calculation Details:
        Compute the maximum daily flow for each hydrological year in
        the record. Compute the decimal logarithm of the maxima.
        Calculate the standard deviation and the mean of these
        log-transformed maxima. Multiply the former by 100, and divide
        the result by the latter.

    """
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[hy, :] = np.amax(flows[mask, :], axis=0)
    # calculations for entire time series
    log_f = np.copy(info)
    # replace log10(0) by log10(0.01) if necessary
    log_f[log_f == 0.0] = 0.01
    log_f = np.log10(log_f, dtype=np.float64)
    sfc = np.std(log_f, ddof=1, axis=0) * 100 / np.mean(log_f, axis=0)

    return sfc


def mh19(flows, datetimes, hydro_years, drainage_area):
    """Skewness in annual maximum daily flow.

    :Calculation Details:
        Compute the maximum daily flow for each hydrological year in
        the record. Calculate the skewness using the following
        equation: ::

            N**2 * Σ(log(Q)) - 3 * N * Σ(log(Q)) * Σ(log(Q)**2) + 2 * Σ(log(Q)**3)
            ----------------------------------------------------------------------
                             N * (N-1) * (N-2) * (σ(log(Q)))**3

        where *N* is the number of years in the record, *Q* is the
        vector of annual maxima, *log* in the decimal logarithm, *Σ* is
        the arithmetic sum, and *σ* is the standard deviation.

    """
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[hy, :] = np.amax(flows[mask, :], axis=0)
    # calculations for entire time series
    nb_years = hydro_years.shape[0]
    log_f = np.copy(info)
    # replace log10(0) by log10(0.01) if necessary
    log_f[log_f == 0.0] = 0.01
    log_f = np.log10(log_f, dtype=np.float64)
    sum_log_f = np.sum(log_f, axis=0)
    sum_log_f_2 = np.sum(log_f ** 2, axis=0)
    sum_log_f_3 = np.sum(log_f ** 3, axis=0)
    std_log_f = np.std(log_f, ddof=1, axis=0)
    sfc = (((nb_years ** 2) * sum_log_f_3 - 3 * nb_years * sum_log_f_2
            * sum_log_f + 2 * (sum_log_f ** 3))
           / (nb_years * (nb_years - 1) * (nb_years - 2) * (std_log_f ** 3)))

    return sfc


def mh20(flows, datetimes, hydro_years, drainage_area):
    """Mean annual maximum daily flow normalised by drainage area.

    :Calculation Details:
        Take the maximum for each hydrological year in the daily flow
        record, and divide them by the drainage area. Compute the mean
        of these ratios.

    """
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[hy, :] = np.amax(flows[mask, :], axis=0) / drainage_area
    # calculations for entire time series
    sfc = np.mean(info, axis=0)

    return sfc


def mh21(flows, datetimes, hydro_years, drainage_area):
    """Mean flood volume (above overall median) normalised by overall
    median.

    :Calculation Details:
        Identify the events above the median of the whole daily flow
        record. Calculate the mean volume of these events. Divide this
        mean by the median of the whole daily flow record.

    """
    median = np.median(flows, axis=0)
    # calculations for entire time series
    sfc = calc_events_avg_volume_above(flows, median) / median

    return sfc


def mh22(flows, datetimes, hydro_years, drainage_area):
    """Mean flood volume (above twice overall median) normalised by
    overall median.

    :Calculation Details:
        Identify the events above twice the median of the whole daily
        flow record. Calculate the mean volume of these events. Divide
        this mean by the median of the whole daily flow record.

    """
    median = np.median(flows, axis=0)
    # calculations for entire time series
    sfc = calc_events_avg_volume_above(flows, 3 * median) / median

    return sfc


def mh23(flows, datetimes, hydro_years, drainage_area):
    """Mean flood volume (above 3 times overall median) normalised by
    overall median.

    :Calculation Details:
        Identify the events above third times the median of the whole
        daily flow record. Calculate the mean volume of these events.
        Divide this mean by the median of the whole daily flow record.

    """
    median = np.median(flows, axis=0)
    # calculations for entire time series
    sfc = calc_events_avg_volume_above(flows, 7 * median) / median

    return sfc


# MH24 to MH27 not available
