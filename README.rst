A streamflow characteristics calculator in Python
-------------------------------------------------

.. image:: https://img.shields.io/pypi/v/eflowcalc?color=blue
   :target: https://pypi.python.org/pypi/eflowcalc
   :alt: PyPI Version
.. image:: https://zenodo.org/badge/153001813.svg
   :target: https://zenodo.org/badge/latestdoi/153001813
   :alt: DOI
.. image:: https://img.shields.io/badge/License-GPL%20v3-blue.svg
   :target: https://www.gnu.org/licenses/gpl-3.0
   :alt: License: GPL v3
.. image:: https://img.shields.io/github/workflow/status/ThibHlln/eflowcalc/Tests?label=tests
   :target: https://github.com/ThibHlln/eflowcalc/actions/workflows/tests.yml
   :alt: GitHub Actions Test Workflow Status

`eflowcalc` is an open-source `calculator` of ecological streamflow
characteristics in Python. It is licensed under GNU GPL-3.0.
The package currently gives the Python scientific community access
to 159 ecologically relevant streamflow characteristics inventoried by
`Olden and Poff (2003) <https://doi.org/10.1002/rra.700>`_. A key strength
of `eflowcalc` is the vectorisation of all calculations (using
`numpy <https://github.com/numpy/numpy>`_, and therefore C code in the
background) which makes for very efficient computation of the streamflow
characteristics.

If you are using `eflowcalc`, please consider citing the software as
follows (click on the link to get the DOI of a specific version):

.. pull-quote::

   *Hallouin, T. (XXXX). EFlowCalc: Ecological Streamflow Characteristics Calculator (Version X.X.X). Zenodo.* `<https://doi.org/10.5281/zenodo.2566757>`_

.. rubric:: Brief overview of the API

.. code-block:: python

   from datetime import datetime, timedelta
   import numpy as np
   import eflowcalc as efc

   datetimes = [datetime(2010, 1, 1) + timedelta(days=d) for d in range(3652)]
   streamflows = np.random.uniform(3, 50, 3652)
   drainage_area = 120.7

   ma41 = efc.calculator(efc.ma41, datetimes, streamflows, drainage_area)

   ma41, dh4, ra7 = efc.calculator((efc.ma41, efc.dh4, efc.ra7),
                                   datetimes, streamflows, drainage_area)

.. rubric:: Streamflow characteristics available

The streamflow characteristics currently available in `eflowcalc` are
as follows:

* Magnitude of flow events
   * Average flow events: `ma1`, `ma2`, `ma3`, `ma4`, `ma5`, `ma6`, `ma7`,
     `ma8`, `ma9`, `ma10`, `ma11`, `ma12`, `ma13`, `ma14`, `ma15`, `ma16`,
     `ma17`, `ma18`, `ma19`, `ma20`, `ma21`, `ma22`, `ma23`, `ma24`, `ma25`,
     `ma26`, `ma27`, `ma28`, `ma29`, `ma30`, `ma31`, `ma32`, `ma33`, `ma34`,
     `ma35`, `ma36`, `ma37`, `ma38`, `ma39`, `ma40`, `ma41`, `ma42`, `ma43`,
     `ma44`, `ma45`
   * Low flow events: `ml1`, `ml2`, `ml3`, `ml4`, `ml5`, `ml6`, `ml7`, `ml8`,
     `ml9`, `ml10`, `ml11`, `ml12`, `ml13`, `ml14`, `ml15`, `ml16`, `ml17`,
     `ml18`, `ml19`, `ml20`, `ml21`, `ml22`
   * High flow events: `mh1`, `mh2`, `mh3`, `mh4`, `mh5`, `mh6`, `mh7`, `mh8`,
     `mh9`, `mh10`, `mh11`, `mh12`, `mh13`, `mh14`, `mh15`, `mh16`, `mh17`,
     `mh18`, `mh19`, `mh20`, `mh21`, `mh22`, `mh23`
* Frequency of flow events
   * Low flow events: `fl1`, `fl2`, `fl3`
   * High flow events: `fh1`, `fh2`, `fh3`, `fh4`, `fh5`, `fh6`, `fh7`, `fh8`,
     `fh9`, `fh10`
* Duration of flow events
   * Low flow events: `dl1`, `dl2`, `dl3`, `dl4`, `dl5`, `dl6`, `dl7`, `dl8`,
     `dl9`, `dl10`, `dl11`, `dl12`, `dl13`, `dl14`, `dl15`, `dl16`, `dl17`,
     `dl18`, `dl19`, `dl20`
   * High flow events: `dh1`, `dh2`, `dh3`, `dh4`, `dh5`, `dh6`, `dh7`, `dh8`,
     `dh9`, `dh10`, `dh11`, `dh12`, `dh13`, `dh14`, `dh15`, `dh16`, `dh17`,
     `dh18`, `dh19`, `dh20`, `dh21`
* Timing of flow events
   * Average flow events: `ta1`, `ta2`
   * Low flow events: `tl1`, `tl2`
   * High flow events: `th1`, `th2`
* Rate of change in flow events
   * Average flow events: `ra1`, `ra2`, `ra3`, `ra4`, `ra5`, `ra6`, `ra7`,
     `ra8`, `ra9`

.. note::
   The computations implemented in `eflowcalc` are partially inspired
   by the work of `Henriksen et al. (2006)
   <https://doi.org/10.3133/ofr20061093>`_, however `eflowcalc` is
   neither endorsed by these authors nor by the U.S. Geological Survey.

.. rubric:: Acknowledgement

Early versions of this tool were developed with the financial support of
Ireland's Environmental Protection Agency (Grant Number 2014-W-LS-5).
