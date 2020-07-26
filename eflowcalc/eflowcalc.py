# -*- coding: utf-8 -*-

# This file is part of EFlowCalc:
# A Calculator of Ecological Streamflow Characteristics
# Copyright (C) 2019  Thibault Hallouin (1)
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

from datetime import datetime, timedelta
import numpy as np


def calculator(sfcs, datetimes, streamflows, drainage_area, axis=0,
               hydro_year='01/10', years=None):
    # check the format of the different arguments given,
    # if not compliant, abort
    if not isinstance(datetimes, np.ndarray):
        raise Exception('The DateTimes given are not in a NumPy array.')
    if not datetimes.ndim == 1:
        raise Exception('The DateTimes array is not uni-dimensional.')
    if not np.issubdtype(datetimes.dtype, np.dtype(datetime)):
        raise Exception('The DateTimes array does not contain '
                        'DateTime objects.')
    if not isinstance(streamflows, np.ndarray):
        raise Exception('The streamflow data given is not in a NumPy array.')
    if not ((axis == 0) or (axis == 1)):
        raise Exception('The index for axis must be 0 or 1.')
    try:
        datetime.strptime(hydro_year, '%d/%m')
    except ValueError:
        raise Exception('The \'hydro_year\' argument does not match the '
                        'format \'%d/%m\' or it is semantically incorrect.')

    # check the dimensions of the streamflow data provided
    if streamflows.ndim == 1:
        # if a flat array, reshape to a 2D array
        my_streamflow = np.reshape(streamflows, (streamflows.size, 1))
    elif streamflows.ndim == 2:
        # if a 2D array, transpose if necessary
        if axis == 0:
            my_streamflow = streamflows
        else:
            my_streamflow = streamflows.T
    else:
        # if the array is neither 1D nor 2D, abort
        raise Exception('The streamflow array contains more '
                        'than 2 dimensions.')

    # subset only full hydrological years and determine mask
    # for each hydrological year
    if years and hasattr(years, '__iter__'):
        # i.e. user provided specific hydrological years

        # subset time series to only include hydrological years requested
        my_subset = np.zeros((my_streamflow.shape[0],), dtype=bool)
        for year_ in years:
            start_hydro_year = datetime.strptime(
                '{}/{} 00:00:00'.format(hydro_year, year_),
                "%d/%m/%Y %H:%M:%S"
            )
            end_hydro_year = datetime.strptime(
                '{}/{} 00:00:00'.format(hydro_year, year_ + 1),
                "%d/%m/%Y %H:%M:%S"
            ) - timedelta(days=1)
            my_subset += ((datetimes >= start_hydro_year)
                          & (datetimes <= end_hydro_year))

        my_time = datetimes[my_subset]
        my_streamflow = my_streamflow[my_subset, :]

        # determine mask for each hydrological year requested
        my_masks_hy = np.zeros((len(years), my_streamflow.shape[0]),
                               dtype=bool)
        for y, year_ in enumerate(years):
            start_hydro_year = datetime.strptime(
                '{}/{} 00:00:00'.format(hydro_year, year_),
                "%d/%m/%Y %H:%M:%S"
            )
            end_hydro_year = datetime.strptime(
                '{}/{} 00:00:00'.format(hydro_year, year_ + 1),
                "%d/%m/%Y %H:%M:%S"
            ) - timedelta(days=1)
            my_masks_hy[y, :] = ((my_time >= start_hydro_year)
                                 & (my_time <= end_hydro_year))

            # check that there is no invalid or missing data
            if np.isnan(my_streamflow[my_masks_hy[y, :], :]).any():
                raise Exception('The hydrological year {} contain(s) '
                                'invalid values (NaN).'.format(hydro_year))
            if not (my_streamflow[my_masks_hy[y, :], :].shape[0]
                    == (end_hydro_year - start_hydro_year).days + 1):
                raise Exception('The hydrological year {} is not complete '
                                '(missing days).'.format(hydro_year))

    else:
        # i.e. user did not provide specific hydrological years,
        # so work on the whole time series

        # trim head and tail of time series to only include
        # full hydrological years
        start = datetime.strptime(
            '{}/{} 00:00:00'.format(hydro_year, datetimes[0].year),
            '%d/%m/%Y %H:%M:%S'
        )
        head = start if datetimes[0] <= start else start.replace(
            year=datetimes[0].year + 1)

        end = datetime.strptime(
            '{}/{} 00:00:00'.format(hydro_year, datetimes[-1].year),
            '%d/%m/%Y %H:%M:%S'
        ) - timedelta(days=1)
        tail = end if datetimes[-1] >= end else end.replace(
            year=datetimes[-1].year - 1)

        my_time = datetimes[(datetimes >= head) & (datetimes <= tail)]
        my_streamflow = my_streamflow[(datetimes >= head)
                                      & (datetimes <= tail), :]

        # check that there is no invalid or missing data
        if np.isnan(my_streamflow).any():
            raise Exception('The simulation(s) time series contain(s) '
                            'invalid values (NaN).')
        if not my_streamflow.shape[0] == (my_time[-1] - my_time[0]).days + 1:
            raise Exception('The simulation(s) time series is (are) '
                            'not complete (missing days).')

        # determine mask for each hydrological year in the whole series
        # (from start to end)
        my_masks_hy = np.zeros(
            ((my_time[-1].year - my_time[0].year), my_streamflow.shape[0]),
            dtype=bool
        )
        for y, year_ in enumerate(range(my_time[0].year, my_time[-1].year, 1)):
            start_hydro_year = datetime.strptime(
                '{}/{} 00:00:00'.format(hydro_year, year_),
                "%d/%m/%Y %H:%M:%S"
            )
            end_hydro_year = datetime.strptime(
                '{}/{} 00:00:00'.format(hydro_year, year_ + 1),
                "%d/%m/%Y %H:%M:%S"
            ) - timedelta(days=1)
            my_masks_hy[y, :] = ((my_time >= start_hydro_year)
                                 & (my_time <= end_hydro_year))

    # calculate the requested streamflow characteristic(s)
    if hasattr(sfcs, '__iter__'):
        calc_sfc = np.zeros((len(sfcs), my_streamflow.shape[1]),
                            dtype=np.float32)
        calc_sfc[:] = np.nan
        for i, sfc in enumerate(sfcs):
            calc_sfc[i, :] = sfc(my_streamflow, my_time, my_masks_hy,
                                 drainage_area)
    else:
        calc_sfc = np.zeros((1, my_streamflow.shape[1]), dtype=np.float32)
        calc_sfc[:] = np.nan
        calc_sfc[0, :] = sfcs(my_streamflow, my_time, my_masks_hy,
                              drainage_area)

    # return array in its original orientation
    if axis == 0:
        return calc_sfc
    else:
        return calc_sfc.T
