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

from datetime import datetime, timedelta
import numpy as np


def calculator(sfcs, datetimes, streamflows, drainage_area,
               hydro_year='01/10', years=None, axis=0):
    """Calculate streamflow characteristics for one time series stored
    in a 1D array (or several time series of equal length stored  in a
    2D array). Typically used for streamflow time series for a same
    period and for a same catchment.

    :Parameters:

        sfcs: (sequence of) `eflowcalc` SFC functions
            The (sequence of) streamflow characteristic(s) to be
            calculated for the given *streamflows* series.

            *Parameter example:* ::

                sfcs=ma41

            *Parameter example:* ::

                sfcs=(ma41, dh4, ra7)

            *Parameter example:* ::

                sfcs=everything

        datetimes: array-like object
            The array of datetimes corresponding to the date and time
            of the *streamflows* values. Must be uni-dimensional.

        streamflows: array-like object
            The array of daily streamflow values in cubic metres per
            second on which to calculate the given *sfcs*. Note, the
            array can be one- or two-dimensional. If it is 2D, the time
            dimension must be the one specified through *axis*.

        drainage_area: `int` or `float`
            The drainage area of the catchment in square kilometres
            for which *streamflows* are provided.

        hydro_year: `str`, optional
            The day and month of the beginning of the hydrological
            (or water) year. Typically '01/10' (i.e. 1st of October)
            for the Northern Hemisphere, and '01/07' (i.e. 1st of July)
            for the Southern Hemisphere. If not provided, set to default
            value '01/10'.

        years: sequence of `int`, optional
            A sequence of years to use to subset the *streamflows*
            series. If not provided, no subset is carried out and the
            whole series is considered by the calculator.

        axis: `int`, optional
            The axis along which the *streamflows* time dimension is
            if *streamflows* is a 2D array. If not provided, set to
            default value 0 (i.e. time along first axis).

    **Examples**

    >>> from datetime import datetime, timedelta
    >>> times = [
    ...     datetime(2010, 1, 1) + timedelta(days=d)
    ...     for d in range(3652)
    ... ]
    >>> import numpy
    >>> numpy.random.seed(7)
    >>> flows = numpy.random.uniform(3, 50, 3652)
    >>> import eflowcalc as efc
    >>> print(
    ...     efc.calculator(
    ...         efc.ma1,
    ...         times, flows, drainage_area=147.
    ...     )
    ... )
    [26.29431]
    >>> print(
    ...     efc.calculator(
    ...         [efc.ma1, efc.dh7],
    ...         times, flows, drainage_area=147.
    ...     )
    ... )
    [[26.29431  ]
     [ 2.8495195]]

    Computations on multiple streamflow series at once are possible.

    >>> flows = numpy.random.uniform(3, 50, (3652, 3))
    >>> print(
    ...     efc.calculator(
    ...         [efc.ma1, efc.dh7],
    ...         times, flows, drainage_area=147.
    ...     )
    ... )
    [[26.548698  26.465912  26.271872 ]
     [ 2.4092593  2.8745487  2.2962854]]
    >>> flows = numpy.random.uniform(3, 50, (3, 3652))
    >>> print(
    ...     efc.calculator(
    ...         [efc.ma1, efc.dh7],
    ...         times, flows, drainage_area=147.,
    ...         axis=1
    ...     )
    ... )
    [[26.409723   4.8691554]
     [26.398853   3.6545587]
     [26.251245   3.4656363]]

    """
    # check the format of the different arguments given,
    # if not compliant, abort
    datetimes = np.asarray(datetimes)
    if not datetimes.shape:
        raise TypeError('datetimes must be an array')
    if not np.issubdtype(datetimes.dtype, np.dtype(datetime)):
        raise TypeError('datetimes array must contain datetimes')
    if not datetimes.ndim == 1:
        raise ValueError('datetimes array must be uni-dimensional')
    streamflows = np.asarray(streamflows)
    if not datetimes.shape:
        raise TypeError('streamflows must be an array')
    if not np.issubdtype(streamflows.dtype, np.number):
        raise TypeError('streamflows array must contain numerical values')
    if axis not in (0, 1):
        raise IndexError('index for axis must be 0 or 1')
    try:
        datetime.strptime(hydro_year, '%d/%m')
    except ValueError:
        raise ValueError('\'hydro_year\' argument does not match format '
                         '\'%d/%m\' or is semantically incorrect')

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
        raise ValueError('streamflows array contains more than 2 dimensions')

    # check streamflows' time dimension compatible with datetimes
    if my_streamflow.shape[0] != datetimes.shape[0]:
        raise RuntimeError(
            'streamflows time length ({}) not equal to datetimes '
            'length ({})'.format(my_streamflow.shape[0], datetimes.shape[0])
        )

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
                raise ValueError('hydrological year {} with invalid '
                                 'values (NaN)'.format(hydro_year))
            if not (my_streamflow[my_masks_hy[y, :], :].shape[0]
                    == (end_hydro_year - start_hydro_year).days + 1):
                raise ValueError('hydrological year {} not complete '
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
            raise ValueError('streamflow series with invalid values (NaN)')
        if not my_streamflow.shape[0] == (my_time[-1] - my_time[0]).days + 1:
            raise Exception('streamflow series not complete (missing days)')

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
        if streamflows.ndim == 1:
            calc_sfc = calc_sfc[0]

    # return array in its original orientation
    if axis == 0:
        return calc_sfc
    else:
        return calc_sfc.T
