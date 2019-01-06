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

import numpy as np


def rolling_window(a, window):
    # From Erik Rigtorp (http://www.rigtorp.se/2011/01/01/rolling-statistics-numpy.html)
    shape = a.shape[:-1] + (a.shape[-1] - window + 1, window)
    strides = a.strides + (a.strides[-1],)
    return np.lib.stride_tricks.as_strided(a, shape=shape, strides=strides)


def count_events(arr, threshold, typ='high'):
    if typ == 'high':
        m = arr > threshold
    else:
        m = arr < threshold
    return np.sum((np.diff(m * 1, axis=-1) > 0), axis=-1) + m[:, 0]


def count_days(arr, threshold, typ='high'):
    if typ == 'high':
        m = arr > threshold
    else:
        m = arr < threshold
    return np.sum(m, axis=-1)


def count_reversals(arr):
    diff = np.diff(arr, axis=-1)
    diff[diff == 0] = np.nan

    # https://stackoverflow.com/a/41191127
    mask = np.isnan(diff)
    idx = np.where(~mask, np.arange(mask.shape[-1]), 0)
    np.maximum.accumulate(idx, axis=-1, out=idx)
    diff[mask] = diff[np.nonzero(mask)[0], idx[mask]]

    diff_pos = (diff > 0)
    diff_neg = (diff < 0)
    events_pos = np.sum((np.diff(diff_pos * 1, axis=-1) > 0), axis=-1) + diff_pos[:, 0]
    events_neg = np.sum((np.diff(diff_neg * 1, axis=-1) > 0), axis=-1) + diff_neg[:, 0]

    return events_pos + events_neg - 1


def calc_events_avg_duration(arr, threshold, typ='high'):
    if typ == 'high':
        m = arr > threshold
    else:
        m = arr < threshold
    count = np.sum((np.diff(m * 1, axis=-1) > 0), axis=-1) + m[:, 0]
    avg_duration = np.true_divide(np.sum(m * 1, axis=-1), count, where=(count != 0))
    avg_duration[count == 0] = 0.0
    return avg_duration


def calc_events_avg_volume_above(arr, threshold):
    m = arr > threshold
    count = np.sum((np.diff(m * 1, axis=-1) > 0), axis=-1) + m[:, 0]
    above = arr - threshold
    above[~m] = 0.0
    avg_volume = np.true_divide(np.sum(above, axis=-1), count, where=(count != 0))
    avg_volume[count == 0] = 0.0
    return avg_volume


def calc_bfi(arr):
    return np.amin(np.mean(rolling_window(arr, 7), axis=-1), axis=-1) / np.mean(arr, axis=-1)
