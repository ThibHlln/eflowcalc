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
from .tools import calc_bfi, rolling_window, calc_events_avg_volume_above


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# AVERAGE FLOWS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# MA1 - Mean of the entire record
def ma1(flows, datetimes, hydro_years, drainage_area):
    # calculations for entire time series
    sfc = np.mean(flows, axis=0)

    return sfc


# MA2 - Median of the entire record
def ma2(flows, datetimes, hydro_years, drainage_area):
    # calculations for entire time series
    sfc = np.median(flows, axis=0)

    return sfc


# MA3 - Mean of the annual coefficients of variation
def ma3(flows, datetimes, hydro_years, drainage_area):
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[hy, :] = np.std(flows[mask, :], ddof=1, axis=0) / np.mean(flows[mask, :], axis=0)
    # calculations for entire time series
    sfc = np.mean(info, axis=0) * 100

    return sfc


# MA4 - Variability in the percentiles on log10 for the entire record
def ma4(flows, datetimes, hydro_years, drainage_area):
    log_f = np.copy(flows)
    log_f[log_f == 0.0] = 0.01  # replace log10(0) by log10(0.01) if necessary
    log_f = np.log10(log_f, dtype=np.float64)
    # calculations for entire time series
    perc = np.percentile(log_f, (5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95), axis=0)
    sfc = np.std(perc, ddof=1, axis=0) * 100 / np.mean(perc, axis=0)

    return sfc


# MA5 - Skewness of the entire record
def ma5(flows, datetimes, hydro_years, drainage_area):
    # calculations for entire time series
    sfc = np.mean(flows, axis=0) / np.median(flows, axis=0)

    return sfc


# MA6 - Ratio for range between 90th and 10th percentile for the entire record
def ma6(flows, datetimes, hydro_years, drainage_area):
    # calculations for entire time series
    sfc = np.percentile(flows, 90, axis=0) / np.percentile(flows, 10, axis=0)

    return sfc


# MA7 - Ratio for range between 80th and 20th percentile for the entire record
def ma7(flows, datetimes, hydro_years, drainage_area):
    # calculations for entire time series
    sfc = np.percentile(flows, 80, axis=0) / np.percentile(flows, 20, axis=0)

    return sfc


# MA8 - Ratio for range between 75th and 25th percentile for the entire record
def ma8(flows, datetimes, hydro_years, drainage_area):
    # calculations for entire time series
    sfc = np.percentile(flows, 75, axis=0) / np.percentile(flows, 25, axis=0)

    return sfc


# MA9 - Spread in range between 90th and 10th percentile on log10 for the entire record
def ma9(flows, datetimes, hydro_years, drainage_area):
    log_f = np.copy(flows)
    log_f[log_f == 0.0] = 0.01  # replace log10(0) by log10(0.01) if necessary
    log_f = np.log10(log_f, dtype=np.float64)
    med_f = np.log10(np.median(flows, axis=0), dtype=np.float64)
    # calculations for entire time series
    sfc = (np.percentile(log_f, 90, axis=0) - np.percentile(log_f, 10, axis=0)) / med_f

    return sfc


# MA10 - Spread in range between 80th and 20th percentile on log10 for the entire record
def ma10(flows, datetimes, hydro_years, drainage_area):
    log_f = np.copy(flows)
    log_f[log_f == 0.0] = 0.01  # replace log10(0) by log10(0.01) if necessary
    log_f = np.log10(log_f, dtype=np.float64)
    med_f = np.log10(np.median(flows, axis=0), dtype=np.float64)
    # calculations for entire time series
    sfc = (np.percentile(log_f, 80, axis=0) - np.percentile(log_f, 20, axis=0)) / med_f

    return sfc


# MA11 - Spread in range between 75th and 25th percentile on log10 for the entire record
def ma11(flows, datetimes, hydro_years, drainage_area):
    log_f = np.copy(flows)
    log_f[log_f == 0.0] = 0.01  # replace log10(0) by log10(0.01) if necessary
    log_f = np.log10(log_f, dtype=np.float64)
    med_f = np.log10(np.median(flows, axis=0), dtype=np.float64)
    # calculations for entire time series
    sfc = (np.percentile(log_f, 75, axis=0) - np.percentile(log_f, 25, axis=0)) / med_f

    return sfc


# MA12 - Mean of January flows
def ma12(flows, datetimes, hydro_years, drainage_area):
    # calculations per month
    info = np.array(pd.DataFrame(flows, index=datetimes).groupby(lambda x: x.month).mean())
    # calculations for entire time series
    sfc = info[0, :]

    return sfc


# MA13 - Mean of February flows
def ma13(flows, datetimes, hydro_years, drainage_area):
    # calculations per month
    info = np.array(pd.DataFrame(flows, index=datetimes).groupby(lambda x: x.month).mean())
    # calculations for entire time series
    sfc = info[1, :]

    return sfc


# MA14 - Mean of March flows
def ma14(flows, datetimes, hydro_years, drainage_area):
    # calculations per month
    info = np.array(pd.DataFrame(flows, index=datetimes).groupby(lambda x: x.month).mean())
    # calculations for entire time series
    sfc = info[2, :]

    return sfc


# MA15 - Mean of April flows
def ma15(flows, datetimes, hydro_years, drainage_area):
    # calculations per month
    info = np.array(pd.DataFrame(flows, index=datetimes).groupby(lambda x: x.month).mean())
    # calculations for entire time series
    sfc = info[3, :]

    return sfc


# MA16 - Mean of May flows
def ma16(flows, datetimes, hydro_years, drainage_area):
    # calculations per month
    info = np.array(pd.DataFrame(flows, index=datetimes).groupby(lambda x: x.month).mean())
    # calculations for entire time series
    sfc = info[4, :]

    return sfc


# MA17 - Mean of June flows
def ma17(flows, datetimes, hydro_years, drainage_area):
    # calculations per month
    info = np.array(pd.DataFrame(flows, index=datetimes).groupby(lambda x: x.month).mean())
    # calculations for entire time series
    sfc = info[5, :]

    return sfc


# MA18 - Mean of July flows
def ma18(flows, datetimes, hydro_years, drainage_area):
    # calculations per month
    info = np.array(pd.DataFrame(flows, index=datetimes).groupby(lambda x: x.month).mean())
    # calculations for entire time series
    sfc = info[6, :]

    return sfc


# MA19 - Mean of August flows
def ma19(flows, datetimes, hydro_years, drainage_area):
    # calculations per month
    info = np.array(pd.DataFrame(flows, index=datetimes).groupby(lambda x: x.month).mean())
    # calculations for entire time series
    sfc = info[7, :]

    return sfc


# MA20 - Mean of September flows
def ma20(flows, datetimes, hydro_years, drainage_area):
    # calculations per month
    info = np.array(pd.DataFrame(flows, index=datetimes).groupby(lambda x: x.month).mean())
    # calculations for entire time series
    sfc = info[8, :]

    return sfc


# MA21 - Mean of October flows
def ma21(flows, datetimes, hydro_years, drainage_area):
    # calculations per month
    info = np.array(pd.DataFrame(flows, index=datetimes).groupby(lambda x: x.month).mean())
    # calculations for entire time series
    sfc = info[9, :]

    return sfc


# MA22 - Mean of November flows
def ma22(flows, datetimes, hydro_years, drainage_area):
    # calculations per month
    info = np.array(pd.DataFrame(flows, index=datetimes).groupby(lambda x: x.month).mean())
    # calculations for entire time series
    sfc = info[10, :]

    return sfc


# MA23 - Mean of December flows
def ma23(flows, datetimes, hydro_years, drainage_area):
    # calculations per month
    info = np.array(pd.DataFrame(flows, index=datetimes).groupby(lambda x: x.month).mean())
    # calculations for entire time series
    sfc = info[11, :]

    return sfc


# MA24 - Variability of January flow
def ma24(flows, datetimes, hydro_years, drainage_area):
    # calculations per month for each year
    std_ = pd.DataFrame(flows, index=datetimes).groupby(lambda x: (x.year, x.month)).std()
    mean_ = pd.DataFrame(flows, index=datetimes).groupby(lambda x: (x.year, x.month)).mean()
    # calculations per month for the entire series
    cv_ = (std_ / mean_).groupby(lambda x: x[1]).mean()
    info = np.array(cv_)
    # calculations for entire time series
    sfc = info[0, :] * 100

    return sfc


# MA25 - Variability of February flow
def ma25(flows, datetimes, hydro_years, drainage_area):
    # calculations per month for each year
    std_ = pd.DataFrame(flows, index=datetimes).groupby(lambda x: (x.year, x.month)).std()
    mean_ = pd.DataFrame(flows, index=datetimes).groupby(lambda x: (x.year, x.month)).mean()
    # calculations per month for the entire series
    cv_ = (std_ / mean_).groupby(lambda x: x[1]).mean()
    info = np.array(cv_)
    # calculations for entire time series
    sfc = info[1, :] * 100

    return sfc


# MA26 - Variability of March flow
def ma26(flows, datetimes, hydro_years, drainage_area):
    # calculations per month for each year
    std_ = pd.DataFrame(flows, index=datetimes).groupby(lambda x: (x.year, x.month)).std()
    mean_ = pd.DataFrame(flows, index=datetimes).groupby(lambda x: (x.year, x.month)).mean()
    # calculations per month for the entire series
    cv_ = (std_ / mean_).groupby(lambda x: x[1]).mean()
    info = np.array(cv_)
    # calculations for entire time series
    sfc = info[2, :] * 100

    return sfc


# MA27 - Variability of April flow
def ma27(flows, datetimes, hydro_years, drainage_area):
    # calculations per month for each year
    std_ = pd.DataFrame(flows, index=datetimes).groupby(lambda x: (x.year, x.month)).std()
    mean_ = pd.DataFrame(flows, index=datetimes).groupby(lambda x: (x.year, x.month)).mean()
    # calculations per month for the entire series
    cv_ = (std_ / mean_).groupby(lambda x: x[1]).mean()
    info = np.array(cv_)
    # calculations for entire time series
    sfc = info[3, :] * 100

    return sfc


# MA28 - Variability of May flow
def ma28(flows, datetimes, hydro_years, drainage_area):
    # calculations per month for each year
    std_ = pd.DataFrame(flows, index=datetimes).groupby(lambda x: (x.year, x.month)).std()
    mean_ = pd.DataFrame(flows, index=datetimes).groupby(lambda x: (x.year, x.month)).mean()
    # calculations per month for the entire series
    cv_ = (std_ / mean_).groupby(lambda x: x[1]).mean()
    info = np.array(cv_)
    # calculations for entire time series
    sfc = info[4, :] * 100

    return sfc


# MA29 - Variability of June flow
def ma29(flows, datetimes, hydro_years, drainage_area):
    # calculations per month for each year
    std_ = pd.DataFrame(flows, index=datetimes).groupby(lambda x: (x.year, x.month)).std()
    mean_ = pd.DataFrame(flows, index=datetimes).groupby(lambda x: (x.year, x.month)).mean()
    # calculations per month for the entire series
    cv_ = (std_ / mean_).groupby(lambda x: x[1]).mean()
    info = np.array(cv_)
    # calculations for entire time series
    sfc = info[5, :] * 100

    return sfc


# MA30 - Variability of July flow
def ma30(flows, datetimes, hydro_years, drainage_area):
    # calculations per month for each year
    std_ = pd.DataFrame(flows, index=datetimes).groupby(lambda x: (x.year, x.month)).std()
    mean_ = pd.DataFrame(flows, index=datetimes).groupby(lambda x: (x.year, x.month)).mean()
    # calculations per month for the entire series
    cv_ = (std_ / mean_).groupby(lambda x: x[1]).mean()
    info = np.array(cv_)
    # calculations for entire time series
    sfc = info[6, :] * 100

    return sfc


# MA31 - Variability of August flow
def ma31(flows, datetimes, hydro_years, drainage_area):
    # calculations per month for each year
    std_ = pd.DataFrame(flows, index=datetimes).groupby(lambda x: (x.year, x.month)).std()
    mean_ = pd.DataFrame(flows, index=datetimes).groupby(lambda x: (x.year, x.month)).mean()
    # calculations per month for the entire series
    cv_ = (std_ / mean_).groupby(lambda x: x[1]).mean()
    info = np.array(cv_)
    # calculations for entire time series
    sfc = info[7, :] * 100

    return sfc


# MA32 - Variability of September flow
def ma32(flows, datetimes, hydro_years, drainage_area):
    # calculations per month for each year
    std_ = pd.DataFrame(flows, index=datetimes).groupby(lambda x: (x.year, x.month)).std()
    mean_ = pd.DataFrame(flows, index=datetimes).groupby(lambda x: (x.year, x.month)).mean()
    # calculations per month for the entire series
    cv_ = (std_ / mean_).groupby(lambda x: x[1]).mean()
    info = np.array(cv_)
    # calculations for entire time series
    sfc = info[8, :] * 100

    return sfc


# MA33 - Variability of October flow
def ma33(flows, datetimes, hydro_years, drainage_area):
    # calculations per month for each year
    std_ = pd.DataFrame(flows, index=datetimes).groupby(lambda x: (x.year, x.month)).std()
    mean_ = pd.DataFrame(flows, index=datetimes).groupby(lambda x: (x.year, x.month)).mean()
    # calculations per month for the entire series
    cv_ = (std_ / mean_).groupby(lambda x: x[1]).mean()
    info = np.array(cv_)
    # calculations for entire time series
    sfc = info[9, :] * 100

    return sfc


# MA34 - Variability of November flow
def ma34(flows, datetimes, hydro_years, drainage_area):
    # calculations per month for each year
    std_ = pd.DataFrame(flows, index=datetimes).groupby(lambda x: (x.year, x.month)).std()
    mean_ = pd.DataFrame(flows, index=datetimes).groupby(lambda x: (x.year, x.month)).mean()
    # calculations per month for the entire series
    cv_ = (std_ / mean_).groupby(lambda x: x[1]).mean()
    info = np.array(cv_)
    # calculations for entire time series
    sfc = info[10, :] * 100

    return sfc


# MA35 - Variability of December flow
def ma35(flows, datetimes, hydro_years, drainage_area):
    # calculations per month for each year
    std_ = pd.DataFrame(flows, index=datetimes).groupby(lambda x: (x.year, x.month)).std()
    mean_ = pd.DataFrame(flows, index=datetimes).groupby(lambda x: (x.year, x.month)).mean()
    # calculations per month for the entire series
    cv_ = (std_ / mean_).groupby(lambda x: x[1]).mean()
    info = np.array(cv_)
    # calculations for entire time series
    sfc = info[11, :] * 100

    return sfc


# MA36 - Skewness of mean monthly flow using min and max
def ma36(flows, datetimes, hydro_years, drainage_area):
    # calculations per month for each year
    mean_ = np.array(pd.DataFrame(flows, index=datetimes).groupby(lambda x: (x.year, x.month)).mean())
    # calculations for entire time series
    sfc = (np.amax(mean_, axis=0) - np.amin(mean_, axis=0)) / np.median(mean_, axis=0)

    return sfc


# MA37 - Skewness of mean monthly flow using 25th and 75th percentiles
def ma37(flows, datetimes, hydro_years, drainage_area):
    # calculations per month for each year
    mean_ = np.array(pd.DataFrame(flows, index=datetimes).groupby(lambda x: (x.year, x.month)).mean())
    # calculations for entire time series
    sfc = (np.percentile(mean_, 75, axis=0) - np.percentile(mean_, 25, axis=0)) / np.median(mean_, axis=0)

    return sfc


# MA38 - Skewness of mean monthly flow using 10th and 90th percentiles
def ma38(flows, datetimes, hydro_years, drainage_area):
    # calculations per month for each year
    mean_ = np.array(pd.DataFrame(flows, index=datetimes).groupby(lambda x: (x.year, x.month)).mean())
    # calculations for entire time series
    sfc = (np.percentile(mean_, 90, axis=0) - np.percentile(mean_, 10, axis=0)) / np.median(mean_, axis=0)

    return sfc


# MA39 - Variability in mean monthly flow
def ma39(flows, datetimes, hydro_years, drainage_area):
    # calculations per month for each year
    mean_ = np.array(pd.DataFrame(flows, index=datetimes).groupby(lambda x: (x.year, x.month)).mean())
    # calculations for entire time series
    sfc = np.std(mean_, ddof=1, axis=0) * 100 / np.mean(mean_, axis=0)

    return sfc


# MA40 - Skewness in mean monthly flow
def ma40(flows, datetimes, hydro_years, drainage_area):
    # calculations per month for each year
    mean_ = np.array(pd.DataFrame(flows, index=datetimes).groupby(lambda x: (x.year, x.month)).mean())
    # calculations for entire time series
    sfc = (np.mean(mean_, axis=0) - np.median(mean_, axis=0)) / np.median(mean_, axis=0)

    return sfc


# MA41 - Mean annual daily flow normalised with drainage area
def ma41(flows, datetimes, hydro_years, drainage_area):
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[hy, :] = np.mean(flows[mask, :], axis=0) / drainage_area
    # calculations for entire time series
    sfc = np.mean(info, axis=0)

    return sfc


# MA42 - Skewness of mean annual flow using min and max
def ma42(flows, datetimes, hydro_years, drainage_area):
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[hy, :] = np.mean(flows[mask, :], axis=0)
    # calculations for entire time series
    sfc = (np.amax(info, axis=0) - np.amin(info, axis=0)) / np.median(info, axis=0)

    return sfc


# MA43 - Skewness of mean annual flow using 25th and 75th percentiles
def ma43(flows, datetimes, hydro_years, drainage_area):
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[hy, :] = np.mean(flows[mask, :], axis=0)
    # calculations for entire time series
    sfc = (np.percentile(info, 75, axis=0) - np.percentile(info, 25, axis=0)) / np.median(info, axis=0)

    return sfc


# MA44 - Skewness of mean annual flow using 10th and 90th percentiles
def ma44(flows, datetimes, hydro_years, drainage_area):
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[hy, :] = np.mean(flows[mask, :], axis=0)
    # calculations for entire time series
    sfc = (np.percentile(info, 90, axis=0) - np.percentile(info, 10, axis=0)) / np.median(info, axis=0)

    return sfc


# MA45 - Skewness in mean annual flow
def ma45(flows, datetimes, hydro_years, drainage_area):
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[hy, :] = np.mean(flows[mask, :], axis=0)
    # calculations for entire time series
    sfc = (np.mean(info, axis=0) - np.median(info, axis=0)) / np.median(info, axis=0)

    return sfc


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOW FLOWS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ML1 - Minimum January flow
def ml1(flows, datetimes, hydro_years, drainage_area):
    # calculations per month for each year
    min_ = pd.DataFrame(flows, index=datetimes).groupby(lambda x: (x.year, x.month)).min()
    info = np.array(min_.groupby(lambda x: x[1]).mean())
    # calculations for entire time series
    sfc = info[0, :]

    return sfc


# ML2 - Minimum February flow
def ml2(flows, datetimes, hydro_years, drainage_area):
    # calculations per month for each year
    min_ = pd.DataFrame(flows, index=datetimes).groupby(lambda x: (x.year, x.month)).min()
    info = np.array(min_.groupby(lambda x: x[1]).mean())
    # calculations for entire time series
    sfc = info[1, :]

    return sfc


# ML3 - Minimum March flow
def ml3(flows, datetimes, hydro_years, drainage_area):
    # calculations per month for each year
    min_ = pd.DataFrame(flows, index=datetimes).groupby(lambda x: (x.year, x.month)).min()
    info = np.array(min_.groupby(lambda x: x[1]).mean())
    # calculations for entire time series
    sfc = info[2, :]

    return sfc


# ML4 - Minimum April flow
def ml4(flows, datetimes, hydro_years, drainage_area):
    # calculations per month for each year
    min_ = pd.DataFrame(flows, index=datetimes).groupby(lambda x: (x.year, x.month)).min()
    info = np.array(min_.groupby(lambda x: x[1]).mean())
    # calculations for entire time series
    sfc = info[3, :]

    return sfc


# ML5 - Minimum May flow
def ml5(flows, datetimes, hydro_years, drainage_area):
    # calculations per month for each year
    min_ = pd.DataFrame(flows, index=datetimes).groupby(lambda x: (x.year, x.month)).min()
    info = np.array(min_.groupby(lambda x: x[1]).mean())
    # calculations for entire time series
    sfc = info[4, :]

    return sfc


# ML6 - Minimum June flow
def ml6(flows, datetimes, hydro_years, drainage_area):
    # calculations per month for each year
    min_ = pd.DataFrame(flows, index=datetimes).groupby(lambda x: (x.year, x.month)).min()
    info = np.array(min_.groupby(lambda x: x[1]).mean())
    # calculations for entire time series
    sfc = info[5, :]

    return sfc


# ML7 - Minimum July flow
def ml7(flows, datetimes, hydro_years, drainage_area):
    # calculations per month for each year
    min_ = pd.DataFrame(flows, index=datetimes).groupby(lambda x: (x.year, x.month)).min()
    info = np.array(min_.groupby(lambda x: x[1]).mean())
    # calculations for entire time series
    sfc = info[6, :]

    return sfc


# ML8 - Minimum August flow
def ml8(flows, datetimes, hydro_years, drainage_area):
    # calculations per month for each year
    min_ = pd.DataFrame(flows, index=datetimes).groupby(lambda x: (x.year, x.month)).min()
    info = np.array(min_.groupby(lambda x: x[1]).mean())
    # calculations for entire time series
    sfc = info[7, :]

    return sfc


# ML9 - Minimum September flow
def ml9(flows, datetimes, hydro_years, drainage_area):
    # calculations per month for each year
    min_ = pd.DataFrame(flows, index=datetimes).groupby(lambda x: (x.year, x.month)).min()
    info = np.array(min_.groupby(lambda x: x[1]).mean())
    # calculations for entire time series
    sfc = info[8, :]

    return sfc


# ML10 - Minimum October flow
def ml10(flows, datetimes, hydro_years, drainage_area):
    # calculations per month for each year
    min_ = pd.DataFrame(flows, index=datetimes).groupby(lambda x: (x.year, x.month)).min()
    info = np.array(min_.groupby(lambda x: x[1]).mean())
    # calculations for entire time series
    sfc = info[9, :]

    return sfc


# ML11 - Minimum November flow
def ml11(flows, datetimes, hydro_years, drainage_area):
    # calculations per month for each year
    min_ = pd.DataFrame(flows, index=datetimes).groupby(lambda x: (x.year, x.month)).min()
    info = np.array(min_.groupby(lambda x: x[1]).mean())
    # calculations for entire time series
    sfc = info[10, :]

    return sfc


# ML12 - Minimum December flow
def ml12(flows, datetimes, hydro_years, drainage_area):
    # calculations per month for each year
    min_ = pd.DataFrame(flows, index=datetimes).groupby(lambda x: (x.year, x.month)).min()
    info = np.array(min_.groupby(lambda x: x[1]).mean())
    # calculations for entire time series
    sfc = info[11, :]

    return sfc


# ML13 - Variability in minimum monthly flow
def ml13(flows, datetimes, hydro_years, drainage_area):
    # calculations per month for each year
    min_ = np.array(pd.DataFrame(flows, index=datetimes).groupby(lambda x: (x.year, x.month)).min())
    # calculations for entire time series
    sfc = np.std(min_, ddof=1, axis=0) * 100 / np.mean(min_, axis=0)

    return sfc


# ML14 - Mean of minimum annual flow normalised by annual median flow
def ml14(flows, datetimes, hydro_years, drainage_area):
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[hy, :] = np.amin(flows[mask, :], axis=0) / np.median(flows[mask, :], axis=0)
    # calculations for entire time series
    sfc = np.mean(info, axis=0)

    return sfc


# ML15 - Mean of minimum annual flow normalised by annual mean flow
def ml15(flows, datetimes, hydro_years, drainage_area):
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[hy, :] = np.amin(flows[mask, :], axis=0) / np.mean(flows[mask, :], axis=0)
    # calculations for entire time series
    sfc = np.mean(info, axis=0)

    return sfc


# ML16 - Median of minimum annual flow normalised by annual median flow
def ml16(flows, datetimes, hydro_years, drainage_area):
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[hy, :] = np.amin(flows[mask, :], axis=0) / np.median(flows[mask, :], axis=0)
    # calculations for entire time series
    sfc = np.median(info, axis=0)

    return sfc


# ML17 - Base flow 1
def ml17(flows, datetimes, hydro_years, drainage_area):
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[hy, :] = calc_bfi(flows[mask, :])
    # calculations for entire time series
    sfc = np.mean(info, axis=0)

    return sfc


# ML18 - Variability in base flow 1
def ml18(flows, datetimes, hydro_years, drainage_area):
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[hy, :] = calc_bfi(flows[mask, :])
    # calculations for entire time series
    sfc = np.std(info, ddof=1, axis=0) * 100 / np.mean(info, axis=0)

    return sfc


# ML19 - Base flow 2
def ml19(flows, datetimes, hydro_years, drainage_area):
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[hy, :] = np.amin(flows[mask, :], axis=0) / np.mean(flows[mask, :], axis=0)
    # calculations for entire time series
    sfc = np.mean(info, axis=0) * 100

    return sfc


# ML20 - Base flow 3
def ml20(flows, datetimes, hydro_years, drainage_area):
    # calculations for entire time series
    nb_blocks = flows.shape[0] // 5
    flows_5day_blocks = flows[0:5 * nb_blocks, :].reshape(nb_blocks, 5, flows.shape[1])
    blocks_min = np.amin(flows_5day_blocks, axis=1)
    rolling_blocks_min = np.amin(rolling_window(blocks_min, 3), axis=1)

    base_flow = np.zeros(blocks_min.shape, dtype=np.float64)
    base_flow[:] = np.nan

    check = blocks_min[1:-1, :] * 0.90 < rolling_blocks_min
    base_flow[1:-1, :][check] = blocks_min[1:-1, :][check]
    base_flow[0, :] = blocks_min[0, :]
    base_flow[-1, :] = blocks_min[-1, :]

    base_flow_ = pd.DataFrame(base_flow).interpolate(method='linear', axis=0).values

    sfc = np.sum(base_flow_ * 5, axis=0) / np.sum(flows, axis=0)

    return sfc


# ML21 - Variability in annual minimum daily flow
def ml21(flows, datetimes, hydro_years, drainage_area):
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[hy, :] = np.amin(flows[mask, :], axis=0)
    # calculations for entire time series
    sfc = np.std(info, ddof=1, axis=0) * 100 / np.mean(info, axis=0)

    return sfc


# ML22 - Mean in annual minimum daily flow
def ml22(flows, datetimes, hydro_years, drainage_area):
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[hy, :] = np.amin(flows[mask, :], axis=0) / drainage_area
    # calculations for entire time series
    sfc = np.mean(info, axis=0)

    return sfc


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# HIGH FLOWS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# MH1 - Maximum January flow
def mh1(flows, datetimes, hydro_years, drainage_area):
    # calculations per month for each year
    max_ = pd.DataFrame(flows, index=datetimes).groupby(lambda x: (x.year, x.month)).max()
    info = np.array(max_.groupby(lambda x: x[1]).mean())
    # calculations for entire time series
    sfc = info[0, :]

    return sfc


# MH2 - Maximum February flow
def mh2(flows, datetimes, hydro_years, drainage_area):
    # calculations per month for each year
    max_ = pd.DataFrame(flows, index=datetimes).groupby(lambda x: (x.year, x.month)).max()
    info = np.array(max_.groupby(lambda x: x[1]).mean())
    # calculations for entire time series
    sfc = info[1, :]

    return sfc


# MH3 - Maximum March flow
def mh3(flows, datetimes, hydro_years, drainage_area):
    # calculations per month for each year
    max_ = pd.DataFrame(flows, index=datetimes).groupby(lambda x: (x.year, x.month)).max()
    info = np.array(max_.groupby(lambda x: x[1]).mean())
    # calculations for entire time series
    sfc = info[2, :]

    return sfc


# MH4 - Maximum April flow
def mh4(flows, datetimes, hydro_years, drainage_area):
    # calculations per month for each year
    max_ = pd.DataFrame(flows, index=datetimes).groupby(lambda x: (x.year, x.month)).max()
    info = np.array(max_.groupby(lambda x: x[1]).mean())
    # calculations for entire time series
    sfc = info[3, :]

    return sfc


# MH5 - Maximum May flow
def mh5(flows, datetimes, hydro_years, drainage_area):
    # calculations per month for each year
    max_ = pd.DataFrame(flows, index=datetimes).groupby(lambda x: (x.year, x.month)).max()
    info = np.array(max_.groupby(lambda x: x[1]).mean())
    # calculations for entire time series
    sfc = info[4, :]

    return sfc


# MH6 - Maximum June flow
def mh6(flows, datetimes, hydro_years, drainage_area):
    # calculations per month for each year
    max_ = pd.DataFrame(flows, index=datetimes).groupby(lambda x: (x.year, x.month)).max()
    info = np.array(max_.groupby(lambda x: x[1]).mean())
    # calculations for entire time series
    sfc = info[5, :]

    return sfc


# MH7 - Maximum July flow
def mh7(flows, datetimes, hydro_years, drainage_area):
    # calculations per month for each year
    max_ = pd.DataFrame(flows, index=datetimes).groupby(lambda x: (x.year, x.month)).max()
    info = np.array(max_.groupby(lambda x: x[1]).mean())
    # calculations for entire time series
    sfc = info[6, :]

    return sfc


# MH8 - Maximum August flow
def mh8(flows, datetimes, hydro_years, drainage_area):
    # calculations per month for each year
    max_ = pd.DataFrame(flows, index=datetimes).groupby(lambda x: (x.year, x.month)).max()
    info = np.array(max_.groupby(lambda x: x[1]).mean())
    # calculations for entire time series
    sfc = info[7, :]

    return sfc


# MH9 - Maximum September flow
def mh9(flows, datetimes, hydro_years, drainage_area):
    # calculations per month for each year
    max_ = pd.DataFrame(flows, index=datetimes).groupby(lambda x: (x.year, x.month)).max()
    info = np.array(max_.groupby(lambda x: x[1]).mean())
    # calculations for entire time series
    sfc = info[8, :]

    return sfc


# MH10 - Maximum October flow
def mh10(flows, datetimes, hydro_years, drainage_area):
    # calculations per month for each year
    max_ = pd.DataFrame(flows, index=datetimes).groupby(lambda x: (x.year, x.month)).max()
    info = np.array(max_.groupby(lambda x: x[1]).mean())
    # calculations for entire time series
    sfc = info[9, :]

    return sfc


# MH11 - Maximum November flow
def mh11(flows, datetimes, hydro_years, drainage_area):
    # calculations per month for each year
    max_ = pd.DataFrame(flows, index=datetimes).groupby(lambda x: (x.year, x.month)).max()
    info = np.array(max_.groupby(lambda x: x[1]).mean())
    # calculations for entire time series
    sfc = info[10, :]

    return sfc


# MH12 - Maximum December flow
def mh12(flows, datetimes, hydro_years, drainage_area):
    # calculations per month for each year
    max_ = pd.DataFrame(flows, index=datetimes).groupby(lambda x: (x.year, x.month)).max()
    info = np.array(max_.groupby(lambda x: x[1]).mean())
    # calculations for entire time series
    sfc = info[11, :]

    return sfc


# MH13 - Variability in maximum monthly flow
def mh13(flows, datetimes, hydro_years, drainage_area):
    # calculations per month for each year
    max_ = np.array(pd.DataFrame(flows, index=datetimes).groupby(lambda x: (x.year, x.month)).max())
    # calculations for entire time series
    sfc = np.std(max_, ddof=1, axis=0) * 100 / np.mean(max_, axis=0)

    return sfc


# MH14 - Mean of maximum annual flow normalised by annual median flow
def mh14(flows, datetimes, hydro_years, drainage_area):
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[hy, :] = np.amax(flows[mask, :], axis=0) / np.median(flows[mask, :], axis=0)
    # calculations for entire time series
    sfc = np.median(info, axis=0)

    return sfc


# MH15 - 1% exceedance value normalised by median flow for entire record
def mh15(flows, datetimes, hydro_years, drainage_area):
    # calculations for entire time series
    sfc = np.percentile(flows, 99, axis=0) / np.median(flows, axis=0)

    return sfc


# MH16 - 10% exceedance value normalised by median flow for entire record
def mh16(flows, datetimes, hydro_years, drainage_area):
    # calculations for entire time series
    sfc = np.percentile(flows, 90, axis=0) / np.median(flows, axis=0)

    return sfc


# MH17 - 25% exceedance value normalised by median flow for entire record
def mh17(flows, datetimes, hydro_years, drainage_area):
    # calculations for entire time series
    sfc = np.percentile(flows, 75, axis=0) / np.median(flows, axis=0)

    return sfc


# MH18 - Variability in annual maximum daily flow on log10
def mh18(flows, datetimes, hydro_years, drainage_area):
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[hy, :] = np.amax(flows[mask, :], axis=0)
    # calculations for entire time series
    log_f = np.copy(info)
    log_f[log_f == 0.0] = 0.01  # replace log10(0) by log10(0.01) if necessary
    log_f = np.log10(log_f, dtype=np.float64)
    sfc = np.std(log_f, ddof=1, axis=0) * 100 / np.mean(log_f, axis=0)

    return sfc


# MH19 - Skewness in annual maximum daily flow
def mh19(flows, datetimes, hydro_years, drainage_area):
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[hy, :] = np.amax(flows[mask, :], axis=0)
    # calculations for entire time series
    nb_years = hydro_years.shape[0]
    log_f = np.copy(info)
    log_f[log_f == 0.0] = 0.01  # replace log10(0) by log10(0.01) if necessary
    log_f = np.log10(log_f, dtype=np.float64)
    sum_log_f = np.sum(log_f, axis=0)
    sum_log_f_2 = np.sum(log_f ** 2, axis=0)
    sum_log_f_3 = np.sum(log_f ** 3, axis=0)
    std_log_f = np.std(log_f, ddof=1, axis=0)
    sfc = ((nb_years ** 2) * sum_log_f_3 - 3 * nb_years * sum_log_f_2 * sum_log_f + 2 * (sum_log_f ** 3)) / \
          (nb_years * (nb_years - 1) * (nb_years - 2) * (std_log_f ** 3))

    return sfc


# MH20 - Mean in annual maximum daily flow
def mh20(flows, datetimes, hydro_years, drainage_area):
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1]), dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        info[hy, :] = np.amax(flows[mask, :], axis=0) / drainage_area
    # calculations for entire time series
    sfc = np.mean(info, axis=0)

    return sfc


# MH21 - Average volume of flood (above median flow) normalised by median flow over entire period
def mh21(flows, datetimes, hydro_years, drainage_area):
    median = np.median(flows, axis=0)
    # calculations for entire time series
    sfc = calc_events_avg_volume_above(flows, median) / median

    return sfc


# MH22 - Average volume of flood (above twice median flow) normalised by median flow over entire period
def mh22(flows, datetimes, hydro_years, drainage_area):
    median = np.median(flows, axis=0)
    # calculations for entire time series
    sfc = calc_events_avg_volume_above(flows, 3 * median) / median

    return sfc


# MH23 - Average volume of flood (above three times median flow) normalised by median flow over entire period
def mh23(flows, datetimes, hydro_years, drainage_area):
    median = np.median(flows, axis=0)
    # calculations for entire time series
    sfc = calc_events_avg_volume_above(flows, 7 * median) / median

    return sfc


# MH24 to MH27 not available
