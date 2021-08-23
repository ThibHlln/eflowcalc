.. default-role:: obj

v0.1.0
------

Released on 2021-04-26.

.. rubric:: General

* drop support for Python 2.7.x

.. rubric:: API changes

* `calculator` now supports any array-like object, e.g. `numpy.ndarray`, `list`
* `calculator` now returns an array with the same number of dimensions as
  streamflow input array

.. rubric:: Tests

* add `unittest` test suite to check that new commits do not change results
* add `doctest` tests to check that the API as presented in the docs is working
* add GitHub workflow to run tests

.. rubric:: Documentation

* add a documentation website generated with `sphinx`


v0.0.3
------

Released on 2019-09-08.

.. rubric:: Algorithm

* check for the format of the keyword parameter *hydro_year* in `calculator`
  to make sure it is a valid 2-digit day/2-digit month format
* add a conditional statement on the return of `calculator` to send back
  the output array in the same orientation as the input streamflow array
* for `dl19`: deal with zero divide issues and returns zeros in place of NaNs
  when the denominator is null
* for `fl3`: deal with empty slice issues and replaces NaN values by zeros
  when the slice is empty and avoids raising unnecessary warnings
* for `count_reversal`: deal with special case where the null difference
  happens at the very beginning of the period, where NaN values were forward
  replaced (i.e. null difference) as long as the difference remained null.
  Now, after forward replacing, it backward replaces these NaN values
  for this special case at the beginning of the period if necessary.

.. rubric:: Bug fixes

* remove hard coded month and year of start of hydrological year that prevented
  `calculator` from honouring the use of its keyword parameter *hydro_year*
  (if a value different from '01/10' was required)
* for `count_reversal`: correct error on the array index used to replace null
  difference by the most recent non-null value (low likelihood of having
  produced numerically erroneous outputs without running into an exception
  first with previous versions of `eflowcalc` though)

.. rubric:: Documentation

* update the tutorial notebook to correct typos

v0.0.2
------

Released on 2019-02-16.

.. rubric:: General

* add single-sourcing package version

.. rubric:: Scope

* add 142 streamflow characteristics of the 171 inventoried by
  `Olden and Poff (2003) <https://doi.org/10.1002/rra.700)>`_

.. rubric:: Functionality

* change default array orientation to *axis=0*
* add possibility to select specific years where to calculate the SFCs using
  the keyword parameter *years*
* add gathering functions for each flow category (`magnitude`, `frequency`,
  `timing`, `duration`, `rate_change`) and for all SFCs (`everything`)

.. rubric:: Algorithm

* change the degrees of freedom to n-1 in calculation of std
  **[break backwards compatibility]**


v0.0.1
------

Released on 2018-10-26.

* first release