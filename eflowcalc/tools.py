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

import numpy as np


def rolling_window(a, window):
    # By Erik Rigtorp (http://www.rigtorp.se/2011/01/01/rolling-statistics-numpy.html)
    shape = a.shape[:-1] + (a.shape[-1] - window + 1, window)
    strides = a.strides + (a.strides[-1],)
    return np.lib.stride_tricks.as_strided(a, shape=shape, strides=strides)


def count_events(arr, threshold, typ='high'):
    if typ == 'high':
        m = arr > threshold
    else:
        m = arr < threshold
    return np.sum((np.diff(m * 1, axis=-1) > 0), axis=-1) + m[:, 0]


def calc_events_avg_duration(arr, threshold, typ='high'):
    if typ == 'high':
        m = arr > threshold
    else:
        m = arr < threshold
    count = np.sum((np.diff(m * 1, axis=-1) > 0), axis=-1) + m[:, 0]
    avg_duration = np.true_divide(np.sum(m * 1, axis=-1), count, where=(count != 0))
    avg_duration[count == 0] = 0.0
    return avg_duration


def calc_bfi(arr):
    return np.amin(np.mean(rolling_window(arr, 7), axis=-1), axis=-1) / np.mean(arr, axis=-1)
