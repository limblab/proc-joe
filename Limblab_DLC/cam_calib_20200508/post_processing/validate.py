#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 12:43:48 2020

@author: minyoungpark
"""

import numpy as np
import pandas as pd
from numpy import array as arr


def read_single_labeled_2d_data(data_path, bp_interested, offset):
    data = pd.read_csv(data_path, header=[1,2], index_col=0)
    length = len(data.index)
    index = arr(data.index)
    
    coords = np.zeros((length, len(bp_interested), 2))
        
    for bp_idx, bp in enumerate(bp_interested):
        bp_coords = arr(data[bp])
        coords[index, bp_idx, :] = bp_coords[:, :2] + [offset[0], offset[1]]

    return {'length': length,
            'coords': coords
            }

def read_labeled_2d_data(data_paths, bp_interested, offsets):
    data = []
    for data_path, offset in zip(data_paths, offsets):
        data.append(read_single_labeled_2d_data(data_path, bp_interested, offset))
    
    return data
    
    