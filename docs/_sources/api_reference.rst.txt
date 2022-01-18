.. currentmodule:: eflowcalc
.. default-role:: obj

API Reference
=============

`eflowcalc` comes with a `calculator` function that pre-processes the input
arrays and calculates the desired streamflow characteristic(s) (SFC(s)).

.. toctree::
   :maxdepth: 1

   functions/eflowcalc.calculator.rst

.. rubric:: Streamflow Characteristics

`eflowcalc` features 156 SFCs categorised into five categories focussing
on flow duration, flow magnitude, flow frequency, flow timing, and flow
rate of change, respectively.  Each category can be subdivided into low
flow, average flow, and or high flow events.

.. toctree::
   :maxdepth: 2

   functions/eflowcalc.flow_duration.rst
   functions/eflowcalc.flow_magnitude.rst
   functions/eflowcalc.flow_frequency.rst
   functions/eflowcalc.flow_rate_change.rst
   functions/eflowcalc.flow_timing.rst

.. toctree::
   :maxdepth: 1

   bundles/eflowcalc.bundles.rst
