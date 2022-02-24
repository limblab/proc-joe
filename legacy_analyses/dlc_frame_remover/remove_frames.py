# -*- coding: utf-8 -*-
"""
Created on Thu Feb 17 11:28:40 2022

@author: Joseph Sombeck
"""

#%% import useful things



import numpy as np
import matplotlib.pyplot as plt

import remove_frames_utils as utils

import pandas as pd


#%% subsample labeled images, store in new folder.
path = 'D:\\Lab\\Data\\DLCgrasp_training_frames\\20210921_cam_0_b2'
#path = 'D:\\Lab\\Data\\DLCgrasp_training_frames\\20211105_cam_3_b'

out = utils.run_cluster_removal(path, ratio_keep=0.35, des_width=300, des_height=300, n_pca_comp=12)

frames_keep = out[0]
SSW = out[1]
pix_vals_pca = out[2]
all_im = out[3]
n_clusters = out[4]
class_pred = out[5]
clust_centers = out[6]

new_df, new_image_fnames = utils.remove_frames(path, frames_keep)


new_path = path + '_subsampled'
utils.write_subsampled_data(path,new_path, new_df, new_image_fnames)





#%% the following methods make plots related to the subsampling
#%% plot some frames from the same cluster
    
for c in range(15):

    in_clust_idx = np.argwhere(class_pred==c)[0:4]
    f, ax = plt.subplots(nrows=2,ncols=2)
    ax=ax.reshape((-1,))
    counter = 0
    for idx in in_clust_idx:
        ax[counter].imshow(all_im[idx[0]])
        counter=counter+1
        
        
#%% plot 2D pca projection with clusters colored

plt.figure()
for i in range(n_clusters):
    in_cluster = class_pred == i
    plt.plot(pix_vals_pca[in_cluster==1,0],pix_vals_pca[in_cluster==1,1],'.')
    
    



#%% this runs for various numbers of clusters.

path = 'D:\\Lab\\Data\\DLCgrasp_training_frames\\20210921_cam_0_b2'
des_width = 250
des_height = 250

pca_n_components = 15

ratio_keeps = [0.2,0.3,0.4, 0.5, 0.6, 0.8]

SSW_all = []
for i in range(len(ratio_keeps)):
    out = utils.run_cluster_removal(path, ratio_keep=ratio_keeps[i], des_width=250, des_height=250, n_pca_comp=15)
    SSW_all.append(out[1])
    

    
# plot mean SSW for each run
    
plt.figure()
SSW_mean = []
for i in range(len(ratio_keeps)):
    SSW_mean.append(np.mean(SSW_all[i]))
    
plt.plot(ratio_keeps, SSW_mean)