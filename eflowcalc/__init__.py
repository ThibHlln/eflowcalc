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

from .eflowcalc import calculator

from .flow_magnitude import ma26, ma41, ml17, ml20, mh10
from .flow_frequency import fl2, fh6, fh7, fh9
from .flow_duration import dl9, dh4, dh13, dh16
from .flow_timing import ta1, tl1
from .flow_rate_change import ra2, ra7
