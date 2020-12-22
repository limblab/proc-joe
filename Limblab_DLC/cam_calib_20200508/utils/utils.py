#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 23:07:13 2019

@author: minyoungpark
"""

import os
import os.path
import toml
import numpy as np

DEFAULT_CONFIG = {
    'filter': {
        'enabled': False,
        'medfilt': 13,
        'offset_threshold': 25,
        'score_threshold': 0.8,
        'spline': True
    },
    'triangulation': {
        'optim': False
    }
}

def full_path(path):
    path_user = os.path.expanduser(path)
    path_full = os.path.abspath(path_user)
    path_norm = os.path.normpath(path_full)
    return path_norm


def load_config(config_filename):
    if config_filename is None:
        config_filename = 'config.toml'

    if os.path.exists(config_filename):
        config = toml.load(config_filename)
    else:
        config = dict()
    
    for k, v in DEFAULT_CONFIG.items():
        if k not in config:
            config[k] = v
        elif isinstance(v, dict): # handle nested defaults
            for k2, v2 in v.items():
                if k2 not in config[k]:
                    config[k][k2] = v2

    return config



def align_frames(ts_list):
    # this function outputs frame numbers from a list of timestamps. The timestamps 
    # are relative to only to itself, not to the other lists
    # (ts_list[0][0] likely happens at the same time as ts_list[1][0], 
    # but the ts_list entries will not be the same)
    all_cams_frame_list = []
    is_good_frame = []
    good_frame_nums = []
    
    dur_list = np.zeros(len(ts_list))
    n_frames = np.zeros(len(ts_list))
    
    for vid_idx in range(len(ts_list)):
        n_frames[vid_idx] = len(ts_list[vid_idx])
        dur_list[vid_idx] = ts_list[vid_idx][-1]-ts_list[vid_idx][0]
        dt_list = ts_list[vid_idx][1:]-ts_list[vid_idx][0:-1]
        dt = np.median(ts_list[vid_idx][1:]-ts_list[vid_idx][0:-1])

    max_dur_idx = np.argmax(dur_list)
    # condition 1: the four lists have the same timestamp duration. In this case,
    # subtract the first timestamp, and then match based on timestamps (accounting for a bit of noise)
    if(np.all(abs(dur_list-dur_list[0]) < dt*0.5)):
        # use longest length to make ground truth array[0,1,2,3,4....], then
        # subtract first timestamp from each list, and determine the idx for each frame
        master_framenums = np.arange(0,np.ceil(dur_list[max_dur_idx]/dt)+1,1).astype(int)
        master_frame_ts = master_framenums*dt
        
        for vid_idx in range(len(ts_list)):
            ts = ts_list[vid_idx] - ts_list[vid_idx][0]
            frame_list = np.zeros_like(master_framenums) - 1
            # match lists up!
            for t_idx in range(len(ts)):
                # find master_frame_ts that is closest to ts[t_idx]
                master_idx = np.argmin(abs(master_frame_ts-ts[t_idx]))
                dist = abs(master_frame_ts[master_idx] - ts[t_idx])
                # if the distance to that frame is sufficiently small, append frame num to frame_list
                if(dist < dt*0.5):
                    frame_list[master_framenums[master_idx]] = t_idx
                else:
                    print("hmmm")
                
            all_cams_frame_list.append(frame_list)
            # get list of good frames (where there were no dropouts)
            if(vid_idx==0):
                is_good_frame = frame_list >= 0 
            else:
                is_good_frame = np.all([is_good_frame==1, frame_list>=0],axis=0)
    # condition 2: the four lists do not have the same timestamp duration. hmm
    else:
        print("timestamps do not have same duration, something is not right")
    
    good_frame_nums = np.argwhere(is_good_frame)
    good_frame_nums = np.reshape(good_frame_nums,(len(good_frame_nums),))
    return all_cams_frame_list, good_frame_nums
    
    

def get_framenums(vid_indices, videos):

    # info files have the same filename as the real videos, just with .xiinfo appended at end
    # check to see if info files exist first
    append_str= '.xiinfo'
    ts_list = []
    n_frames = np.zeros(len(vid_indices))
    counter = 0
    for vid_idx, vid in zip(vid_indices, videos):  
        # sanity check -- see if filenames has two .avis....if yes, remove one of them
        if(".avi" in vid[-4:] and ".avi" in vid[-8:-4]):
            vid = vid[:-4]
        
        if(not ".avi" in vid[-4:]):
            vid = vid + ".avi"
        
        if(not os.path.exists(vid+append_str)):
            print('xiinfo file not found. This is necessary to align frames')
            continue
        temp_list = np.array([])
        # load in xiinfo store timestamps in ts_list
        info_file = open(vid+append_str,'r')
        lines = info_file.readlines()
        
        for line in lines:
            if("timestamp" in line): # find number and store
                temp_list=np.append(temp_list,(int(line[line.index('\"')+1:-4])))
        
        # pop first ts_list entry if it's nonsense
        if(temp_list[1] - temp_list[0] > 40000*10):
            temp_list = temp_list[1:]
            
        ts_list.append(temp_list)
        
        n_frames[counter] = len(temp_list)
        counter = counter + 1
        
                
    all_cams_frame_list,good_frame_nums = align_frames(ts_list)   
    
    return all_cams_frame_list, good_frame_nums