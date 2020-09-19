.. currentmodule:: eflowcalc
.. default-role:: obj

Tutorial
========

This notebook contains a simple example of the usage of the API of
`eflowcalc` to calculate stream flow characteristics from streamflow
time series.

.. code-block:: python
   :caption: Importing the package and checking its version.

   >>> import eflowcalc
   >>> print(eflowcalc.__version__)
   0.0.3


.. rubric:: Load streamflow time series

An example file is provided in the folder *examples/* in order for anyone
to reproduce this tutorial. Because this is a NetCDF file, we are going
to use the Python package `netCDF4`, but `eflowcalc` is independent of the
file format you are working with because it only requires `numpy` arrays
as inputs for streamflow time series and datetime series. The datetime
series must be made of datetime objects from the datetime package.

.. code-block:: python
   :caption: Reading in the streamflow dataset.

   >>> from netCDF4 import Dataset
   >>> import numpy as np
   >>> from datetime import datetime, timedelta
   >>> with Dataset('examples/catchment.sim.flow.nc', 'r', format='NETCDF4') as f:
   ...     streamflow = f.variables['flow'][:]  # streamflow time series
   ...     timestamps = f.variables['time'][:]  # timestamp series for the period
   ...     time_units = f.variables['time'].units
   >>> print(streamflow.shape, timestamps.shape)
   (20, 4383) (4383,)

In NetCDF files, the time dimension is stored as numerics representing
times elapsed since a reference time (timestamps), but `eflowcalc`
requires datetime objects, so we need to convert the timestamps to
datetimes first.

.. code-block:: python
   :caption: Converting timestamps to datetimes.

   >>> import cftime
   >>> print(timestamps[0], timestamps[-1])
   1096588800.0 1475193600.0
   >>> print(time_units)
   seconds since 1970-01-01 00:00:00.0
   >>> datetimes = cftime.num2date(timestamps, time_units)
   >>> print(datetimes[0], datetimes[-1])
   2004-10-01 00:00:00 2016-09-30 00:00:00


.. rubric:: Calculate one or more streamflow characteristics

Now that the dataset is loaded in memory, it is time to use `eflowcalc`
to calculate the streamflow characteristics from the hydrograph(s). To
do so, import `eflowcalc`, which will give you access to its `calculator`
Python function as well as all streamflow characteristics implemented
in `eflowcalc` (as Python functions as well).

By default, `eflowcalc` expects the time dimension to be on `axis=0`.
In this example, this is not the case, so we need to specify explicitly
that it is on `axis=1`.

.. code-block:: python
   :caption: Calculating only one streamflow characteristic (e.g. `ma41` here).

   >>> from eflowcalc import calculator, ma41
   >>> my_sfc = calculator(ma41, my_dt, my_flow, 1246, axis=1)
   >>> print(my_sfc[0])
   [0.00307121]


.. code-block:: python
   :caption: Calculating multiple streamflow characteristics at once.

   >>> from eflowcalc import calculator, ma41, dh4, ra7
   >>> my_sfcs = calculator((ma41, dh4, ra7), my_dt, my_flow, 1246, axis=1)
   >>> print(my_sfcs[0, :])
   [3.0712057e-03 8.4707642e+00 3.3340059e-02]

It is important to be aware that `eflowcalc` requires strictly continuous
time series of daily streamflow. Moreover, all streamflow characteristics
are only computed on full hydrological years (it will automatically trim
the head and tail of the time series to guarantee so). By default, a
hydrological year starts on the 1st of October, which is commonly the
case in the Northern hemisphere. However, it can be changed to suit any
location using the keyword argument *hydro_year*. For example, if working
on a catchment in the Southern hemisphere, it is likely that the
hydrological year starts on the 1st of July, see example below.

.. code-block:: python
   :caption: Changing the definition for the hydrological year.

   >>> from eflowcalc import calculator, ma41
   >>> my_sfc = calculator(ma41, my_dt, my_flow, 1246, hydro_year='01/07', axis=1)
   >>> print(my_sfc[0])
   [0.00312234]
