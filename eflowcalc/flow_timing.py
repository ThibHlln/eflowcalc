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
import math


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# AVERAGE FLOWS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def ta1(flows, datetimes, hydro_years, drainage_area):
    """Constancy by `Colwell (1974) <https://doi.org/10.2307/1940366>`_
    applied to flows.

    :Calculation Details:
        Compute the decimal logarithm of the daily flow values.
        Calculate the decimal logarithm of the overall mean daily flow
        for the entire record. Compute the Colwell matrix featuring 365
        rows for 365 d in a year (ignoring last day of February for leap
        years) and 11 columns for 11 flow states (break points are 0.10,
        0.25, 0.50, 0.75, 1.00, 1.25, 1.50, 1.75, 2.00, and 2.25 times
        the log mean daily flow calculated previously) for each
        hydrological year, incrementally adding to the tally in each
        cell from year to year.

        Calculate Y, the sum of each column (vector), and Z, the sum of
        the whole matrix (scalar). Divide the elements of vector Y by
        scalar Z. Multiply the elements of the new vector by their
        respective decimal log-transformed value; sum the elements of
        the vector to obtain a scalar; and multiply by minus one to
        obtain the uncertainty with respect to state H(Y). Divide H(Y)
        by the decimal log of the number of states (11), and subtract
        this ratio from one.

    """
    mean = np.mean(flows, axis=0)
    log_mean = np.log10(mean)
    # calculations per hydrological year
    colwell = np.zeros((365, 11, flows.shape[1]), dtype=int)
    log_f = np.copy(flows)
    # replace log10(0) by log10(0.01) if necessary
    log_f[log_f == 0.0] = 0.01
    log_f = np.log10(log_f, dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        mask_no_lpy = np.copy(mask)
        if np.sum(mask) == 366:
            replacement = np.ones((366,), dtype=bool)
            replacement[151] = False  # remove 29th of February
            mask_no_lpy[mask] = replacement
        colwell[:, 0, :] += (log_f[mask_no_lpy, :] < (0.10 * log_mean)) * 1
        colwell[:, 1, :] += \
            ((log_f[mask_no_lpy, :] >= (0.10 * log_mean))
             & (log_f[mask_no_lpy, :] < (0.25 * log_mean))) * 1
        colwell[:, 2, :] += \
            ((log_f[mask_no_lpy, :] >= (0.25 * log_mean))
             & (log_f[mask_no_lpy, :] < (0.50 * log_mean))) * 1
        colwell[:, 3, :] += \
            ((log_f[mask_no_lpy, :] >= (0.50 * log_mean))
             & (log_f[mask_no_lpy, :] < (0.75 * log_mean))) * 1
        colwell[:, 4, :] += \
            ((log_f[mask_no_lpy, :] >= (0.75 * log_mean))
             & (log_f[mask_no_lpy, :] < (1.00 * log_mean))) * 1
        colwell[:, 5, :] += \
            ((log_f[mask_no_lpy, :] >= (1.00 * log_mean))
             & (log_f[mask_no_lpy, :] < (1.25 * log_mean))) * 1
        colwell[:, 6, :] += \
            ((log_f[mask_no_lpy, :] >= (1.25 * log_mean))
             & (log_f[mask_no_lpy, :] < (1.50 * log_mean))) * 1
        colwell[:, 7, :] += \
            ((log_f[mask_no_lpy, :] >= (1.50 * log_mean))
             & (log_f[mask_no_lpy, :] < (1.75 * log_mean))) * 1
        colwell[:, 8, :] += \
            ((log_f[mask_no_lpy, :] >= (1.75 * log_mean))
             & (log_f[mask_no_lpy, :] < (2.00 * log_mean))) * 1
        colwell[:, 9, :] += \
            ((log_f[mask_no_lpy, :] >= (2.00 * log_mean))
             & (log_f[mask_no_lpy, :] < (2.25 * log_mean))) * 1
        colwell[:, 10, :] += (log_f[mask_no_lpy, :] >= (2.25 * log_mean)) * 1
    # calculations for entire time series
    # sum up values in each column (i.e. state)
    colwell_y = np.sum(colwell, axis=0)
    # sum up val in each matrix
    colwell_z = np.sum(colwell, axis=(0, 1))
    colwell_yz = np.divide(colwell_y, colwell_z, dtype=np.float64)
    colwell_log_yz = np.log10(colwell_yz, where=(colwell_yz != 0))
    colwell_log_yz[colwell_yz == 0] = np.nan
    colwell_hy = - np.nansum(np.multiply(colwell_yz, colwell_log_yz,
                                         dtype=np.float64),
                             axis=0)
    sfc = 1 - (colwell_hy / np.log10(11))

    return sfc


def ta2(flows, datetimes, hydro_years, drainage_area):
    """Predictability by `Colwell (1974) <https://doi.org/10.2307/1940366>`_
    applied to flows.

    :Calculation Details:
        Compute the decimal logarithm of the daily flow values.
        Calculate the decimal logarithm of the overall mean daily flow
        for the entire record. Compute the Colwell matrix featuring 365
        rows for 365 d in a year (ignoring last day of February for leap
        years) and 11 columns for 11 flow states (break points are 0.10,
        0.25, 0.50, 0.75, 1.00, 1.25, 1.50, 1.75, 2.00, and 2.25 times
        the log mean daily flow calculated previously) for each
        hydrological year, incrementally adding to the tally in each
        cell from year to year.

        Calculate X, the sum of each row (vector), and Z, the sum of
        the whole matrix (scalar). Divide the elements of vector X by
        scalar Z. Multiply the elements of the new vector by their
        respective decimal log-transformed value; sum the elements of
        the vector to obtain a scalar; and multiply by minus one to
        obtain the uncertainty with respect to time H(X).

        Take the Colwell matrix N. Divide the elements of matrix N by
        scalar Z. Multiply the elements of the new matrix by their
        respective decimal log-transformed value; sum the elements of
        the matrix to obtain a scalar; and multiply by minus one to
        obtain the uncertainty with respect to the interaction of time
        and state H(XY).

        Subtract H(X) to H(XY), and divide the result by the decimal log
        of the number of states (11), and subtract this ratio from one.

    """
    mean = np.mean(flows, axis=0)
    log_mean = np.log10(mean)
    # calculations per hydrological year
    colwell = np.zeros((365, 11, flows.shape[1]), dtype=int)
    log_f = np.copy(flows)
    # replace log10(0) by log10(0.01) if necessary
    log_f[log_f == 0.0] = 0.01
    log_f = np.log10(log_f, dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        mask_no_lpy = np.copy(mask)
        if np.sum(mask) == 366:
            replacement = np.ones((366,), dtype=bool)
            replacement[151] = False  # remove 29th of February
            mask_no_lpy[mask] = replacement
        colwell[:, 0, :] += (log_f[mask_no_lpy, :] < (0.10 * log_mean)) * 1
        colwell[:, 1, :] += \
            ((log_f[mask_no_lpy, :] >= (0.10 * log_mean))
             & (log_f[mask_no_lpy, :] < (0.25 * log_mean))) * 1
        colwell[:, 2, :] += \
            ((log_f[mask_no_lpy, :] >= (0.25 * log_mean))
             & (log_f[mask_no_lpy, :] < (0.50 * log_mean))) * 1
        colwell[:, 3, :] += \
            ((log_f[mask_no_lpy, :] >= (0.50 * log_mean))
             & (log_f[mask_no_lpy, :] < (0.75 * log_mean))) * 1
        colwell[:, 4, :] += \
            ((log_f[mask_no_lpy, :] >= (0.75 * log_mean))
             & (log_f[mask_no_lpy, :] < (1.00 * log_mean))) * 1
        colwell[:, 5, :] += \
            ((log_f[mask_no_lpy, :] >= (1.00 * log_mean))
             & (log_f[mask_no_lpy, :] < (1.25 * log_mean))) * 1
        colwell[:, 6, :] += \
            ((log_f[mask_no_lpy, :] >= (1.25 * log_mean))
             & (log_f[mask_no_lpy, :] < (1.50 * log_mean))) * 1
        colwell[:, 7, :] += \
            ((log_f[mask_no_lpy, :] >= (1.50 * log_mean))
             & (log_f[mask_no_lpy, :] < (1.75 * log_mean))) * 1
        colwell[:, 8, :] += \
            ((log_f[mask_no_lpy, :] >= (1.75 * log_mean))
             & (log_f[mask_no_lpy, :] < (2.00 * log_mean))) * 1
        colwell[:, 9, :] += \
            ((log_f[mask_no_lpy, :] >= (2.00 * log_mean))
             & (log_f[mask_no_lpy, :] < (2.25 * log_mean))) * 1
        colwell[:, 10, :] += (log_f[mask_no_lpy, :] >= (2.25 * log_mean)) * 1
    # calculations for entire time series
    # sum up values in each row (i.e. time)
    colwell_x = np.sum(colwell, axis=1)
    # sum up val in each matrix
    colwell_z = np.sum(colwell, axis=(0, 1))
    colwell_xz = np.divide(colwell_x, colwell_z, dtype=np.float64)
    colwell_log_xz = np.log10(colwell_xz, where=(colwell_xz != 0))
    colwell_log_xz[colwell_xz == 0] = np.nan
    colwell_hx = - np.nansum(np.multiply(colwell_xz, colwell_log_xz,
                                         dtype=np.float64),
                             axis=0)
    colwell_nz = np.divide(colwell, colwell_z, dtype=np.float64)
    colwell_log_nz = np.log10(colwell_nz, where=(colwell != 0))
    colwell_log_nz[colwell == 0] = np.nan
    colwell_hxy = - np.nansum(np.multiply(colwell_nz, colwell_log_nz,
                                          dtype=np.float64),
                              axis=(0, 1))
    sfc = (1 - ((colwell_hxy - colwell_hx) / np.log10(11))) * 100

    return sfc


# TA3 not available


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOW FLOWS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def tl1(flows, datetimes, hydro_years, drainage_area):
    """Timing of annual minimum flow.

    :Calculation Details:
        Determine the day of the year (i.e. number) in the Julian
        calendar where the flow is minimal in each hydrological year in
        the daily flow record. Map these dates onto a circular scale as
        angles in radians. Determine the x and y components of these
        angles (i.e. cosine, and sine, respectively). Compute the mean
        of the x and y components separately. Use these means to compute
        the corresponding mean angle, and convert this mean angle back
        to a day of the year in the Julian calendar (using the
        arc-tangent).

    """
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1], 2),
                    dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        min_julian_day = pd.DatetimeIndex(
            datetimes[mask][np.argmin(flows[mask, :], axis=0)]).dayofyear
        min_julian_day = min_julian_day * 2.0 * math.pi / 365.25
        info[hy, :, 0] = np.cos(min_julian_day)
        info[hy, :, 1] = np.sin(min_julian_day)
    # calculations for entire time series
    x = np.mean(info[:, :, 0], axis=0)
    y = np.mean(info[:, :, 1], axis=0)
    tl1_ = np.zeros((flows.shape[0],), dtype=np.float64)
    tl1_[:] = np.nan
    tl1_ = np.arctan(np.divide(y, x, where=(x != 0),
                               dtype=np.float32)) * 180 / math.pi
    tl1_[x < 0] = tl1_[x < 0] + 180.0
    tl1_[(x == 0) & (y > 0)] = 90
    tl1_[(x == 0) & (y < 0)] = 270
    tl1_[tl1_ < 0] = tl1_[tl1_ < 0] + 360
    tl1_ = tl1_ * 365.25 / 360.0
    tl1_[tl1_ == 0] = 365.25
    sfc = np.rint(tl1_)

    return sfc


def tl2(flows, datetimes, hydro_years, drainage_area):
    """Variability in timing of annual minimum flow.

    :Calculation Details:
        Determine the day of the year (i.e. number) in the Julian
        calendar where the flow is minimal in each hydrological year in
        the daily flow record. Map these dates onto a circular scale as
        angles in radians. Determine the x and y components of these
        angles (i.e. cosine, and sine, respectively). Compute the mean
        of the x and y components separately. Compute the variability
        using the following equation: ::

            √(2 * (1 - √(X*X + Y*Y))

        where √ is the square root, and X and Y are the mean of the
        x and y components of the angles, respectively.

        Convert this angular variability back to a day of the year in
        the Julian calendar.

    """
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1], 2),
                    dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        min_julian_day = pd.DatetimeIndex(
            datetimes[mask][np.argmin(flows[mask, :], axis=0)]).dayofyear
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


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# HIGH FLOWS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def th1(flows, datetimes, hydro_years, drainage_area):
    """Timing of annual maximum flow.

    :Calculation Details:
        Determine the day of the year (i.e. number) in the Julian
        calendar where the flow is maximal in each hydrological year in
        the daily flow record. Map these dates onto a circular scale as
        angles in radians. Determine the x and y components of these
        angles (i.e. cosine, and sine, respectively). Compute the mean
        of the x and y components separately. Use these means to compute
        the corresponding mean angle, and convert this mean angle back
        to a day of the year in the Julian calendar (using the
        arc-tangent).

    """
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1], 2),
                    dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        max_julian_day = pd.DatetimeIndex(
            datetimes[mask][np.argmax(flows[mask, :], axis=0)]).dayofyear
        max_julian_day = max_julian_day * 2.0 * math.pi / 365.25
        info[hy, :, 0] = np.cos(max_julian_day)
        info[hy, :, 1] = np.sin(max_julian_day)
    # calculations for entire time series
    x = np.mean(info[:, :, 0], axis=0)
    y = np.mean(info[:, :, 1], axis=0)
    th1_ = np.zeros((flows.shape[0],), dtype=np.float64)
    th1_[:] = np.nan
    th1_ = np.arctan(np.divide(y, x, where=(x != 0),
                               dtype=np.float32)) * 180 / math.pi
    th1_[x < 0] = th1_[x < 0] + 180.0
    th1_[(x == 0) & (y > 0)] = 90
    th1_[(x == 0) & (y < 0)] = 270
    th1_[th1_ < 0] = th1_[th1_ < 0] + 360
    th1_ = th1_ * 365.25 / 360.0
    th1_[th1_ == 0] = 365.25
    sfc = np.rint(th1_)

    return sfc


def th2(flows, datetimes, hydro_years, drainage_area):
    """Variability in timing of annual maximum flow.

    :Calculation Details:
        Determine the day of the year (i.e. number) in the Julian
        calendar where the flow is maximal in each hydrological year in
        the daily flow record. Map these dates onto a circular scale as
        angles in radians. Determine the x and y components of these
        angles (i.e. cosine, and sine, respectively). Compute the mean
        of the x and y components separately. Compute the variability
        using the following equation: ::

            √(2 * (1 - √(X*X + Y*Y))

        where √ is the square root, and X and Y are the mean of the
        x and y components of the angles, respectively.

        Convert this angular variability back to a day of the year in
        the Julian calendar.

    """
    # calculations per hydrological year
    info = np.zeros((hydro_years.shape[0], flows.shape[1], 2),
                    dtype=np.float64)
    for hy, mask in enumerate(hydro_years):
        max_julian_day = pd.DatetimeIndex(
            datetimes[mask][np.argmax(flows[mask, :], axis=0)]).dayofyear
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
