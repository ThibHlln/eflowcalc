# -*- coding: utf-8 -*-

# This file is part of EFlowCalc: A Calculator of Ecological Stream Flow Characteristics
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

from .flow_magnitude import \
    ma1, ma2, ma3, ma4, ma5, ma6, ma7, ma8, ma9, ma10, ma11, ma12, ma13, ma14, ma15, ma16, ma17, ma18, ma19, ma20, \
    ma21, ma22, ma23, ma24, ma25, ma26, ma27, ma28, ma29, ma30, ma31, ma32, ma33, ma34, ma35, ma36, ma37, ma38, ma39, \
    ma40, ma41, ma42, ma43, ma44, ma45, \
    ml1, ml2, ml3, ml4, ml5, ml6, ml7, ml8, ml9, ml10, ml11, ml12, ml13, ml14, ml15, ml16, ml17, ml18, ml19, ml20, \
    ml21, ml22, \
    mh1, mh2, mh3, mh4, mh5, mh6, mh7, mh8, mh9, mh10, mh11, mh12, mh13, mh14, mh15, mh16, mh17, mh18, mh19, mh20, \
    mh21, mh22, mh23
from .flow_frequency import \
    fl1, fl2, fl3, \
    fh1, fh2, fh3, fh4, fh5, fh6, fh7, fh8, fh9, fh10
from .flow_duration import \
    dl1, dl2, dl3, dl4, dl5, dl6, dl7, dl8, dl9, dl10, dl11, dl12, dl13, dl14, dl15, dl16, dl17, dl18, dl19, dl20, \
    dh1, dh2, dh3, dh4, dh5, dh6, dh7, dh8, dh9, dh10, dh11, dh12, dh13, dh14, dh15, dh16, dh17, dh18, dh19, dh20, dh21
from .flow_timing import \
    ta1, ta2, \
    tl1, tl2, \
    th1, th2
from .flow_rate_change import \
    ra1, ra2, ra3, ra4, ra5, ra6, ra7, ra8, ra9

from .version import __version__

magnitude = (
    ma1, ma2, ma3, ma4, ma5, ma6, ma7, ma8, ma9, ma10, ma11, ma12, ma13, ma14, ma15, ma16, ma17, ma18, ma19, ma20, ma21,
    ma22, ma23, ma24, ma25, ma26, ma27, ma28, ma29, ma30, ma31, ma32, ma33, ma34, ma35, ma36, ma37, ma38, ma39, ma40,
    ma41, ma42, ma43, ma44, ma45, ml1, ml2, ml3, ml4, ml5, ml6, ml7, ml8, ml9, ml10, ml11, ml12, ml13, ml14, ml15, ml16,
    ml17, ml18, ml19, ml20, ml21, ml22, mh1, mh2, mh3, mh4, mh5, mh6, mh7, mh8, mh9, mh10, mh11, mh12, mh13, mh14, mh15,
    mh16, mh17, mh18, mh19, mh20, mh21, mh22, mh23
)

frequency = (
    fl1, fl2, fl3, fh1, fh2, fh3, fh4, fh5, fh6, fh7, fh8, fh9, fh10
)

duration = (
    dl1, dl2, dl3, dl4, dl5, dl6, dl7, dl8, dl9, dl10, dl11, dl12, dl13, dl14, dl15, dl16, dl17, dl18, dl19, dl20,
    dh1, dh2, dh3, dh4, dh5, dh6, dh7, dh8, dh9, dh10, dh11, dh12, dh13, dh14, dh15, dh16, dh17, dh18, dh19, dh20, dh21
)

timing = (
    ta1, ta2, tl1, tl2, th1, th2
)

rate_change = (
    ra1, ra2, ra3, ra4, ra5, ra6, ra7, ra8, ra9
)

everything = magnitude + frequency + duration + timing + rate_change
