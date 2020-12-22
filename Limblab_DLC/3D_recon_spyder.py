# -*- coding: utf-8 -*-
"""
Created on Fri Dec 18 11:41:30 2020

@author: Joseph Sombeck
"""
# imports and other stuff
import sys

   
import numpy as np
import toml
import glob
import os
import yaml
import pandas as pd

# add folder with calibration functions to path
calib_folder = r'D:\Lab\GIT\proc-joe\Limblab_DLC\cam_calib_20200508'
sys.path.append(calib_folder)
from utils.utils import load_config
from calibration.intrinsic import calibrate_intrinsic
from calibration.extrinsic import calibrate_extrinsic
from triangulation.triangulate import reconstruct_3d
from utils.vis_utils import generate_three_dim_video
from utils.vis_utils import generate_three_dim_pictures
from utils.triangulation_utils import add_static_points




# set project folder 
project_folder = r'D:\Lab\Data\DLC_videos\Han_20201203_rwFreeReach'

# determine if we are using filtered data or not
use_filtered_data = True
remove_triangulation = True

# using reference frame or use an arbitrary frame?
use_reference_frame = False

# set number of cameras
n_cams = 4



# setup config file 
# load basic toml folder and then fill out relevant entries
parsed_toml = toml.load(calib_folder + r'\config_master.toml')

# upate calib video path and prefix and extension
parsed_toml['calibration']['calib_video_path'] = project_folder + r'\videos\calib'
parsed_toml['calibration']['calib_video_prefix'] = 'Calib_20201203_0000'

# update paths to 2d data while removing (or keeping) videos with filtered
if(use_filtered_data):
    vid_list = glob.glob(project_folder + r'\videos\*filtered.csv')
else:
    vid_list_temp = glob.glob(project_folder + r'\videos\*.csv')
    vid_list = []
    for vid_name in vid_list_temp:
        if vid_name.find('filtered') == -1:
            vid_list.append(vid_name)
    
parsed_toml['paths_to_2d_data'] = vid_list
    
num_vids = len(vid_list)
num_camera_sets = int(num_vids/n_cams)

# update path to save static data
parsed_toml['path_to_save_static_data'] = project_folder + r'\videos'

# update output video path
parsed_toml['output_video_path'] = project_folder + r'\reconstructed-3d-data'

# update triangulation data (or remove if desired)
if(remove_triangulation):
    parsed_toml['triangulation'].pop('axes', None)
    parsed_toml['triangulation'].pop('reference_point', None)
# else, likely leave alone. Could fill this in if we change how we do reference points

# update reconstruction output path and threshold
parsed_toml['triangulation']['reconstruction_threshold'] = 0.7 # 0.7 default
parsed_toml['triangulation']['reconstruction_output_path'] = project_folder + r'\reconstructed-3d-data'

# update labeling scheme and bodyparts interested
# only need to update if using something other than base arm points
parsed_toml['labeling']['scheme'] = []
parsed_toml['labeling']['bodyparts_interested'] = ['shoulder','elbow1','elbow2','wrist1','wrist2','hand1','hand2','hand3','pointX','pointY','pointZ']

recon_config_file = project_folder + r'\recon_config.toml'

with open(recon_config_file,'w+') as file:
    toml.dump(parsed_toml,file)




config = load_config(recon_config_file)

#%% If you already ran calibration you don't need to run these.
calibrate_intrinsic(config)
calibrate_extrinsic(config)

# add static points if reference frame is provided 
if(use_reference_frame):
    labels = ['pointX', 'pointY', 'pointZ']
    
    snapshots = vid_list

    # initialize static
    static = {label : [] for label in labels}
        
    # get labeled reference point for each camera and store in a new file, also copy over other tracking data
    # make sure to overwrite pointX, pointY, and pointZ if already exist
    # labeled reference point is in each csv file in each folder in project_folder + 'labeled-data'
    labeled_folder = project_folder + r'\labeled-data'
    
    video_files = config['paths_to_2d_data']
    pointX_data = []; pointY_data = []; pointZ_data = [];
    
    for vid_idx in range(len(video_files)):
        video_name = os.path.split(video_files[vid_idx])[-1]
        DLC_idx = video_name.rindex("DLC")
    
        video_prefix = video_name[:DLC_idx]
        
        path_to_folder = labeled_folder + '\\' + video_prefix
        csv_file = glob.glob(path_to_folder+ '\\*.csv')[0]
        data = pd.read_csv(csv_file, header=[1,2], index_col=0)
    
        # get (x,y) for each label Store in Static properly
        
        for label in labels:
            x = data[(label, 'x')].values[~np.isnan(data[(label, 'x')].values)]
            y = data[(label, 'y')].values[~np.isnan(data[(label, 'y')].values)]
            x = x[0] # only get first labels
            y = y[0]
            static[label].append([x,y])

    
    add_static_points(config, labels, static, snapshots)
    
    
#%% edit config to only use specific videos for 3D reconstruction?
    
# make 3D recon folder if needed
if(not os.path.isdir(parsed_toml['triangulation']['reconstruction_output_path'])):
    os.mkdir(parsed_toml['triangulation']['reconstruction_output_path'])
    
vid_list_all = config['paths_to_2d_data']

if(num_vids % n_cams != 0):
    print("wrong number of videos or cameras")
    
for i_set in range(num_camera_sets):
    config["paths_to_2d_data"] = vid_list_all[(i_set)*n_cams:(i_set+1)*n_cams]
    #recovery = reconstruct_3d(config,i_set)

config["paths_to_2d_data"] = vid_list_all


#%% Try 3D stick figure video

vid_list_all = config['paths_to_2d_data']
for i_set in range(num_camera_sets):
    config["paths_to_2d_data"] = vid_list_all[(i_set)*n_cams:(i_set+1)*n_cams]
    generate_three_dim_video(config,i_set)
    
config["paths_to_2d_data"] = vid_list_all

#%% Testing if the generated 3D results makes sense, by calculating the distance between wrist and hand

