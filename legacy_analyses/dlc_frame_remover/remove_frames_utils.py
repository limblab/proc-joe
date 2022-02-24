# -*- coding: utf-8 -*-
"""
Created on Thu Feb 17 11:40:40 2022

@author: Joseph Sombeck
"""

from PIL import Image

import numpy as np
import glob
import matplotlib.pyplot as plt

import remove_frames_utils as utils

from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
import pandas as pd
import os
import shutil


def get_pixel_distance(fr1, fr2): 
    # expects two matrices (pixel idx, rgb), returns a matrix of the distance at each pixel. 
    
    dist = np.sqrt(np.sum((fr1-fr2)**2,axis=-1))
    return dist




def get_diff_bounds(frame_diff, axis, des):
    
        
    flattened_diff = np.mean(np.sum(frame_diff,axis=axis),axis=1)
     
    max_idx = 0
    max_val = np.sum(flattened_diff[max_idx:max_idx+des-1])
    
    for i in range(1,len(flattened_diff)):
        curr_val = np.sum(flattened_diff[i:i+des-1])
        if(curr_val > max_val):
            max_idx = i
            max_val = curr_val
    
    
    return [max_idx, max_idx+des]





def check_bounds(bounds, max_val):
    
    if(bounds[0] < 0):
        bounds = bounds + -1*bounds[0]
    if(bounds[1] > max_val):
        bounds = bounds - (bounds[1]-(max_val+1))

    return bounds


def get_SSE(data):
    # samples are along first axis, features along second
    mean_val = np.mean(data, axis=0).reshape((1,-1))
    
    return np.mean(np.square(data-mean_val))


def load_h5_data(fpath):
    label_fname = glob.glob(fpath+'\\*.h5')[0]
    
    # get labeled pixels for each image. Use this to crop images
    df = pd.read_hdf(label_fname)

    return df


def remove_frames(fpath, frames_keep):
    df = load_h5_data(fpath)
    fnames = glob.glob(fpath+'\\*.png')
        
    # get keep_mask
    keep_mask = np.zeros((len(fnames),),dtype=bool)
    keep_mask[frames_keep] = True
    
    # get fnames (index of dataframe)
    fnames = df.index.to_numpy()
    
    # drop frames
    fnames_drop = fnames[keep_mask==False]
    new_df = df.drop(fnames_drop,axis=0)
    new_fnames = fnames[keep_mask==True]
    
    
    return new_df, new_fnames
    
def write_subsampled_data(path,new_path, new_df, new_image_fnames):
    

    label_fname = os.path.basename(glob.glob(path+'\\*.h5')[0])
    
    # check if subsampled folder exists, if not: make it
    if(os.path.isdir(new_path) == False):
        os.mkdir(new_path)
      
    # write labels as h5
    new_df.to_hdf(new_path+os.path.sep+label_fname,key='new_df',mode='w')

    # copy images from path to new path
    for f in new_image_fnames:
        # get basename
        n = os.path.basename(f)
    
        #copy file 
        src = path + os.path.sep + n
        dest = new_path + os.path.sep + n
        shutil.copyfile(src, dest)
    
    
    return 0

def run_cluster_removal(fpath, ratio_keep, des_width, des_height, n_pca_comp):
    fnames = glob.glob(fpath+'\\*.png')
    
    df = load_h5_data(fpath)
    scorer = df.columns.get_level_values(0)
    bodyparts=df.columns.get_level_values(1)
    
    # use bodyparts[0], which is the wrist
    y_vals = df[scorer[0]][bodyparts[0]]['x'].to_numpy()
    x_vals = df[scorer[0]][bodyparts[0]]['y'].to_numpy()
    
    x_vals = np.round(x_vals).astype(int)
    y_vals = np.round(y_vals).astype(int)
    
    
    # get pixels for all frames, only keep those within x and y bounds
    
    pix_vals = np.zeros((des_width*des_height*3,len(fnames)))
    all_im = []
    for i in range(len(fnames)):
        im = Image.open(fnames[i])
        temp_vals = np.array(im)
        
        # get slice vals
        x_bounds = np.array([x_vals[i]-des_width/2, x_vals[i]+des_width/2]).astype(int)
        y_bounds = np.array([y_vals[i]-des_height/2, y_vals[i]+des_height/2]).astype(int)
        
        
        # check for out of bounds
        x_bounds = utils.check_bounds(x_bounds,temp_vals.shape[0])
        y_bounds = utils.check_bounds(y_bounds,temp_vals.shape[1])
        
        
        temp_vals = temp_vals[x_bounds[0]:x_bounds[1], y_bounds[0]:y_bounds[1],:]
        all_im.append(temp_vals)
        pix_vals[:,i] = temp_vals.reshape((-1,))
            
    pix_vals = np.transpose(pix_vals)
    
    # we have way too many features, do pca to reduce dimensionality
    
    pca = PCA(n_components = n_pca_comp)
    pix_vals_pca = pca.fit_transform(pix_vals)
    
    
    # now do clustering to remove a percentage of frames
    
    n_clusters = np.ceil(pix_vals_pca.shape[0]*ratio_keep).astype(int)
    
    kmeans = KMeans(n_clusters=n_clusters,random_state=0).fit(pix_vals_pca)
    
    class_pred = kmeans.predict(pix_vals_pca)
    clust_centers = kmeans.cluster_centers_
    
    
    # for each cluster center, get closest frame and store that
    frames_keep = []
    for i in range(n_clusters):
        in_clust = np.argwhere(class_pred == i)
        
        dist_to_cent = utils.get_pixel_distance(pix_vals_pca[in_clust,:], clust_centers[i,:])
        
        frames_keep.append(in_clust[np.argmin(dist_to_cent)][0])
        
    
    # for each cluster, get SSE within a cluster
        
    SSW = np.zeros((n_clusters,))
    
    for c in range(n_clusters):
        in_clust = class_pred == c
        SSW[c] = utils.get_SSE(pix_vals_pca[in_clust==1,:])
        
        

    return (frames_keep, SSW, pix_vals_pca, all_im, n_clusters, class_pred, clust_centers)