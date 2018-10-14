# -*- coding: utf-8 -*-

# This file is part of HydroEval: A Calculator for Ecological Stream Flow Characteristics
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

from datetime import datetime, timedelta
import numpy as np


def calculator(sfc_s, datetimes, simulation_s, drainage_area, axis=1, hydro_year='01/10'):
    # check the format of the different arguments given
    if not isinstance(datetimes, np.ndarray):
        raise Exception('The DateTimes given are not in a NumPy array.')
    if not datetimes.ndim == 1:
        raise Exception('The DateTimes array is not unidimensional.')
    if not np.issubdtype(datetimes.dtype, np.dtype(datetime)):
        raise Exception('The DateTimes array does not contain DateTime objects.')
    if not isinstance(simulation_s, np.ndarray):
        raise Exception('The Simulations given are not in a NumPy array.')
    if not ((axis == 0) or (axis == 1)):
        raise Exception('The index for axis is not 0 or 1.')

    # check the dimensions of the simulation data provided
    if simulation_s.ndim == 1:
        my_simu = np.reshape(simulation_s, (1, simulation_s.size))
    elif simulation_s.ndim == 2:
        if axis == 0:
            my_simu = simulation_s.T
        else:
            my_simu = simulation_s
    else:
        raise Exception('The simulation array contains more than 2 dimensions.')

    # trim heads and tails of time series to only include full hydrological years
    start = datetime.strptime('{}/{} 00:00:00'.format(hydro_year, datetimes[0].year), '%d/%m/%Y %H:%M:%S')
    head = start if datetimes[0] <= start else start.replace(year=datetimes[0].year + 1)

    end = datetime.strptime('{}/{} 00:00:00'.format(hydro_year, datetimes[-1].year), '%d/%m/%Y %H:%M:%S') - \
        timedelta(days=1)
    tail = end if datetimes[-1] >= end else end.replace(year=datetimes[-1].year - 1)

    my_time = datetimes[(datetimes >= head) & (datetimes <= tail)]
    my_simu = my_simu[:, (datetimes >= head) & (datetimes <= tail)]

    # check that there is no missing data
    if np.isnan(my_simu).any():
        raise Exception('The simulation(s) time series contain(s) invalid values (NaN).')
    if not my_simu.shape[1] == (my_time[-1] - my_time[0]).days + 1:
        raise Exception('The simulation(s) time series is (are) not complete (missing days)')

    # determine mask for each hydrological year
    my_masks = np.zeros(((my_time[-1].year - my_time[0].year), my_simu.shape[1]), dtype=bool)
    for hy, hydro_year in enumerate(range(my_time[0].year, my_time[-1].year, 1)):
        start_hydro_year = datetime.strptime("{}-10-01 00:00:00".format(hydro_year), "%Y-%m-%d %H:%M:%S")
        end_hydro_year = datetime.strptime("{}-09-30 00:00:00".format(hydro_year + 1), "%Y-%m-%d %H:%M:%S")
        my_masks[hy, :] = (my_time >= start_hydro_year) & (my_time <= end_hydro_year)

    # calculate the SFC(s)
    if hasattr(sfc_s, '__iter__'):
        calc_sfc = np.zeros((my_simu.shape[0], len(sfc_s)), dtype=np.float32)
        calc_sfc[:] = np.nan
        for i, sfc in enumerate(sfc_s):
            calc_sfc[:, i] = sfc(my_simu, my_time, my_masks, drainage_area)
    else:
        calc_sfc = np.zeros((my_simu.shape[0], 1), dtype=np.float32)
        calc_sfc[:] = np.nan
        calc_sfc[:, 0] = sfc_s(my_simu, my_time, my_masks, drainage_area)

    return calc_sfc
