#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 11:20:11 2019

@author: minyoungpark
"""

import os
import numpy as np
import pandas as pd
from tqdm import trange
from numpy import array as arr


def read_single_2d_data(data_path, offset, bp_interested):
    data = pd.read_csv(data_path, header=[1,2], index_col=0)
    length = len(data.index)
    index = arr(data.index)
    
    coords = np.zeros((length, len(bp_interested), 2))
    scores = np.zeros((length, len(bp_interested)))
        
    for bp_idx, bp in enumerate(bp_interested):
        bp_coords = arr(data[bp])
        coords[index, bp_idx, :] = bp_coords[:, :2] + [offset[0], offset[1]]
        scores[index, bp_idx] = bp_coords[:, 2]

    return {'length': length,
            'coords': coords,
            'scores': scores}

    
def load_offsets_dict(config, vid_indices):
    offsets_dict = dict()
    for vid_idx in vid_indices:
        # if record_dict is None:
        if 'cameras' not in config or vid_idx not in config['cameras']:
            offsets_dict[vid_idx] = [0, 0]
        else:
            offsets_dict[vid_idx] = config['cameras'][vid_idx]['offset']
        # else:
        #     offsets_dict[cname] = record_dict['cameras'][cname]['video']['ROIPosition']

    return offsets_dict


def load_2d_data(config, vid_indices, bp_interested):
    paths_to_2d_data = config['paths_to_2d_data']
    
    offsets_dict = load_offsets_dict(config, vid_indices)
    
    # TODO: If there is any frame dropping, do interpolation. Now, just assume
    #       that there isn't any.
    all_points_raw = []
    all_scores = []
    
    # all_points_raw = np.zeros((length, len(cam_names), len(bodyparts), 2))
    # all_scores = np.zeros((length, len(cam_names), len(bodyparts)))
    
    for ix_cam, (vid_idx, data_path) in \
            enumerate(zip(vid_indices, paths_to_2d_data)):
        out = read_single_2d_data(data_path, offsets_dict[vid_idx], bp_interested)
        all_points_raw.append(out['coords'])
        all_scores.append(out['scores'])
        
    all_points_raw = np.stack(all_points_raw, axis=1)
    all_scores = np.stack(all_scores, axis=1)
    
    return {'points': all_points_raw,
            'scores': all_scores}


def read_single_labeled_2d_data(data_path, bp_interested, offset):
    data = pd.read_csv(data_path, header=[1,2], index_col=0)
    length = len(data.index)
    indices = arr(data.index)
    
    for i, index in enumerate(indices):
        indices[i] = index.split('/')[-1]
    
    coords = np.zeros((length, len(bp_interested), 2))
        
    for bp_idx, bp in enumerate(bp_interested):
        bp_coords = arr(data[bp])
        coords[:, bp_idx, :] = bp_coords[:, :] + [offset[0], offset[1]]

    return {'length': length,
            'coords': coords,
            'indices': indices
            }

def load_labeled_2d_data(config, vid_indices, bp_interested):
    paths_to_2d_data = config['paths_to_labeled_2d_data']
    
    offsets_dict = load_offsets_dict(config, vid_indices)
    
    # TODO: If there is any frame dropping, do interpolation. Now, just assume
    #       that there isn't any.
    all_points_raw = []
    all_indices = []
    all_lengths = []
    
    # all_points_raw = np.zeros((length, len(cam_names), len(bodyparts), 2))
    # all_scores = np.zeros((length, len(cam_names), len(bodyparts)))
    
    for ix_cam, (vid_idx, data_path) in \
            enumerate(zip(vid_indices, paths_to_2d_data)):
        out = read_single_labeled_2d_data(data_path, bp_interested, offsets_dict[vid_idx])
        all_points_raw.append(out['coords'])
        all_indices.append(out['indices'])
        all_lengths.append(out['length'])

    min_len = min(all_lengths)
    # amin_len = amin(all_lengths)
    
    for i in range(len(all_lengths)):
        all_points_raw[i] = all_points_raw[i][:min_len]
        
    # for j in :
    #     if 

    all_points_raw = np.stack(all_points_raw, axis=1)
    
    return {'points': all_points_raw,
            'indices': all_indices}
    
    
def add_static_points(config, labels, static, snapshots):
    data_paths = config['paths_to_2d_data']
    path_to_save = config['path_to_save_static_data']
    if not os.path.exists(path_to_save):
        print('Path to save does not exist.')
        folder_input = input('Do you want to create this path (folder)? (y/n) ')
        if folder_input is 'y':
            os.mkdir(path_to_save)
        elif folder_input is 'n':
            return
        else:
            print('Wrong input.')
            return
        
    for i, (snapshot, data_path) in enumerate(zip(snapshots, data_paths)):
        data = pd.read_csv(data_path, header=[0,1,2], index_col=0)
        for label in labels:
            if np.isnan(static[label][i][0]):
                x = np.zeros(len(data))
                y = np.zeros(len(data))
                likelihood = np.zeros(len(data))
            else:
                x = np.ones(len(data)) * static[label][i][0]
                y = np.ones(len(data)) * static[label][i][1]
                likelihood = np.ones(len(data))
            data = data.join(pd.DataFrame(x,
                                          columns=pd.MultiIndex.from_product([[snapshot],[label],['x']]),
                                          index=data.index))
            data = data.join(pd.DataFrame(y,
                                          columns=pd.MultiIndex.from_product([[snapshot],[label],['y']]),
                                          index=data.index))
            data = data.join(pd.DataFrame(likelihood,
                                          columns=pd.MultiIndex.from_product([[snapshot],[label],['likelihood']]),
                                          index=data.index))
            data.to_csv(os.path.join(path_to_save, 'cam_' + str(i) +'.csv'), mode='w')