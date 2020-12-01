#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 23:07:13 2019

@author: minyoungpark
"""

import os
import cv2
import numpy as np
import pandas as pd

from tqdm import trange
from numpy import array as arr
from scipy import optimize

from utils.triangulation_utils import load_2d_data, load_labeled_2d_data

from utils.calibration_utils import get_video_path, load_intrinsics, load_extrinsics


def expand_matrix(mtx):
    z = np.zeros((4,4))
    z[0:3,0:3] = mtx[0:3,0:3]
    z[3,3] = 1
    return z


def reproject_points(p3d, points2d, camera_mats):
    proj = np.dot(camera_mats, p3d)
    proj = proj[:, :2] / proj[:, 2, None]
    return proj


def reprojection_error(p3d, points2d, camera_mats):
    proj = np.dot(camera_mats, p3d)
    proj = proj[:, :2] / proj[:, 2, None]
    errors = np.linalg.norm(proj - points2d, axis=1)
    return np.mean(errors)


def distort_points_cams(points, camera_mats):
    out = []
    for i in range(len(points)):
        point = np.append(points[i], 1)
        mat = camera_mats[i]
        new = mat.dot(point)[:2]
        out.append(new)
    return np.array(out)


def reprojection_error_und(p3d, points2d, camera_mats, camera_mats_dist):
    proj = np.dot(camera_mats, p3d)
    proj = proj[:, :2] / proj[:, 2, None]
    proj_d = distort_points_cams(proj, camera_mats_dist)
    points2d_d = distort_points_cams(points2d, camera_mats_dist)
    errors = np.linalg.norm(proj_d - points2d_d, axis=1)
    return np.mean(errors)


def triangulate_simple(points, camera_mats):
    num_cams = len(camera_mats)
    A = np.zeros((num_cams*2, 4))
    for i in range(num_cams):
        x, y = points[i]
        mat = camera_mats[i]
        A[(i*2):(i*2+1)] = x*mat[2]-mat[0]
        A[(i*2+1):(i*2+2)] = y*mat[2]-mat[1]
    u, s, vh = np.linalg.svd(A, full_matrices=True)
    p3d = vh[-1]
    p3d = p3d / p3d[3]
    return p3d


def triangulate_points(the_points, cam_mats):
    p3ds = []
    errors = []
    for ptnum in range(the_points.shape[0]):
        points = the_points[ptnum]
        good = ~np.isnan(points[:, 0])
        p3d = triangulate_simple(points[good], cam_mats[good])
        err = reprojection_error(p3d, points[good], cam_mats[good])
        p3ds.append(p3d)
        errors.append(err)
    p3ds = np.array(p3ds)
    errors = np.array(errors)
    return p3ds, errors


def optim_error_fun(points, camera_mats):
    def fun(x):
        p3d = np.array([x[0], x[1], x[2], 1])
        proj = np.dot(camera_mats, p3d)
        resid = points - proj[:, :2] / proj[:, 2, None]
        return resid.flatten()
        # return np.linalg.norm(resid, axis=1)
    return fun


def triangulate_optim(points, camera_mats, max_error=20):
    try:
        p3d = triangulate_simple(points, camera_mats)
        # error = reprojection_error(p3d, points, camera_mats)
    except np.linalg.linalg.LinAlgError:
        return np.array([0,0,0,0])

    fun = optim_error_fun(points, camera_mats)
    try:
        res = optimize.least_squares(fun, p3d[:3], loss='huber', f_scale=1e-3)
        x = res.x
        p3d = np.array([x[0], x[1], x[2], 1])
    except ValueError:
        pass

    return p3d


def proj(u, v):
    """Project u onto v"""
    return u * np.dot(v,u) / np.dot(u,u)


def ortho(u, v):
    """Orthagonalize u with respect to v"""
    return u - proj(v, u)


def get_median(all_points_3d, ix):
    pts = all_points_3d[:, ix]
    pts = pts[~np.isnan(pts[:, 0])]
    return np.median(pts, axis=0)


def correct_coordinate_frame(config, all_points_3d, bodyparts):
    """Given a config and a set of points and bodypart names, this function will rotate the coordinate frame to match the one in config"""
    bp_interested = config['labeling']['bodyparts_interested']
    bp_index = dict(zip(bp_interested, range(len(bp_interested))))
    axes_mapping = dict(zip('xyz', range(3)))
    ref_point = config['triangulation']['reference_point']
    axes_spec = config['triangulation']['axes']
    a_dirx, a_l, a_r = axes_spec[0]
    b_dirx, b_l, b_r = axes_spec[1]
    a_dir = axes_mapping[a_dirx]
    b_dir = axes_mapping[b_dirx]
    ## find the missing direction
    done = np.zeros(3, dtype='bool')
    done[a_dir] = True
    done[b_dir] = True
    c_dir = np.where(~done)[0][0]
    a_lv = get_median(all_points_3d, bp_index[a_l])
    a_rv = get_median(all_points_3d, bp_index[a_r])
    b_lv = get_median(all_points_3d, bp_index[b_l])
    b_rv = get_median(all_points_3d, bp_index[b_r])
    a_diff = a_rv - a_lv
    b_diff = ortho(b_rv - b_lv, a_diff)
    M = np.zeros((3,3))
    M[a_dir] = a_diff
    M[b_dir] = b_diff
    if (a_dir==0 and b_dir==1) or (a_dir==1 and b_dir==2) or (a_dir==2 and b_dir==0):
        M[c_dir] = np.cross(a_diff, b_diff)
    else:
        M[c_dir] = np.cross(b_diff, a_diff)
    M /= np.linalg.norm(M, axis=1)[:,None]
    center = get_median(all_points_3d, bp_index[ref_point])
    # all_points_3d_adj = np.dot(all_points_3d - center, M.T)
    all_points_3d_adj = (all_points_3d - center).dot(M.T)
    center_new = get_median(all_points_3d_adj, bp_index[ref_point])
    all_points_3d_adj = all_points_3d_adj - center_new
    recovery = {'center_new': center_new,
                'registration_mat': M,
                'center': center}
    return all_points_3d_adj, recovery


def undistort_points(all_points_raw, cam_names, intrinsics):
    all_points_und = np.zeros(all_points_raw.shape)

    for ix_cam, cam_name in enumerate(cam_names):
        calib = intrinsics[cam_name]
        points = all_points_raw[:, ix_cam].reshape(-1, 1, 2)
        points_new = cv2.undistortPoints(
            points, arr(calib['camera_mat']), arr(calib['dist_coeff']))
        all_points_und[:, ix_cam] = points_new.reshape(
            all_points_raw[:, ix_cam].shape)

    return all_points_und

def reconstruct_3d(config, **kwargs):
    path, videos, vid_indices = get_video_path(config)
    bp_interested = config['labeling']['bodyparts_interested']
    reconstruction_threshold = config['triangulation']['reconstruction_threshold']
    if config['triangulation'].get('reconstruction_output_path') is None:
        output_path = kwargs.get('output_path', '')
    else:
        output_path = config['triangulation']['reconstruction_output_path']
    try:
        intrinsics = load_intrinsics(path, vid_indices)
    except:
        print("Intrinsic calibration output does not exist.")
        return
    try:
        extrinsics = load_extrinsics(path)
    except:
        print("Extrinsic calibration output does not exist.")
        return
    # intrinsics, extrinsics = load_calib_new(config)
    cam_mats = []
    cam_mats_dist = []
    for vid_idxs in vid_indices:
        mat = arr(extrinsics[vid_idxs])
        left = arr(intrinsics[vid_idxs]['camera_mat'])
        cam_mats.append(mat)
        cam_mats_dist.append(left)
    cam_mats = arr(cam_mats)
    cam_mats_dist = arr(cam_mats_dist)
    out = load_2d_data(config, vid_indices, bp_interested)
    all_points_raw = out['points']
    all_scores = out['scores']
    all_points_und = undistort_points(all_points_raw, vid_indices, intrinsics)
    length = all_points_raw.shape[0]
    shape = all_points_raw.shape
    all_points_3d = np.zeros((shape[0], shape[2], 3))
    all_points_3d.fill(np.nan)
    errors = np.zeros((shape[0], shape[2]))
    errors.fill(np.nan)
    scores_3d = np.zeros((shape[0], shape[2]))
    scores_3d.fill(np.nan)
    num_cams = np.zeros((shape[0], shape[2]))
    num_cams.fill(np.nan)
    all_points_und[all_scores < reconstruction_threshold] = np.nan
    for i in trange(all_points_und.shape[0], ncols=70):
        for j in range(all_points_und.shape[2]):
            pts = all_points_und[i, :, j, :]
            good = ~np.isnan(pts[:, 0])
            if np.sum(good) >= 2:
                # TODO: make triangulation type configurable
                # p3d = triangulate_optim(pts[good], cam_mats[good])
                p3d = triangulate_simple(pts[good], cam_mats[good])
                all_points_3d[i, j] = p3d[:3]
                errors[i,j] = reprojection_error_und(p3d, pts[good], cam_mats[good], cam_mats_dist[good])
                num_cams[i,j] = np.sum(good)
                scores_3d[i,j] = np.min(all_scores[i, :, j][good])
    if 'reference_point' in config['triangulation'] and 'axes' in config['triangulation']:
        all_points_3d_adj, recovery = correct_coordinate_frame(config, all_points_3d, bp_interested)
    else:
        all_points_3d_adj = all_points_3d
    dout = pd.DataFrame()
    for bp_num, bp in enumerate(bp_interested):
        for ax_num, axis in enumerate(['x','y','z']):
            dout[bp + '_' + axis] = all_points_3d_adj[:, bp_num, ax_num]
        dout[bp + '_error'] = errors[:, bp_num]
        dout[bp + '_ncams'] = num_cams[:, bp_num]
        dout[bp + '_score'] = scores_3d[:, bp_num]
    dout['fnum'] = np.arange(length)
    dout.to_csv(os.path.join(output_path, 'output_3d_data.csv'), index=False)
    if 'reference_point' in config['triangulation'] and 'axes' in config['triangulation']:
        return recovery
    else:
        return None
    
# def reconstruct_3d(config, **kwargs):
#     path, videos, vid_indices = get_video_path(config)
#     bp_interested = config['labeling']['bodyparts_interested']
#     reconstruction_threshold = config['triangulation']['reconstruction_threshold']
#     if config['triangulation'].get('reconstruction_output_path') is None:
#         output_path = kwargs.get('output_path', '')
#     else:
#         output_path = config['triangulation']['reconstruction_output_path']

#     try:
#         intrinsics = load_intrinsics(path, vid_indices)
#     except:
#         print("Intrinsic calibration output does not exist.")
#         return

#     try:
#         extrinsics = load_extrinsics(path)
#     except:
#         print("Extrinsic calibration output does not exist.")
#         return

#     cam_mats = []
#     cam_mats_dist = []

#     for vid_idxs in vid_indices:
#         mat = arr(extrinsics[vid_idxs])
#         left = arr(intrinsics[vid_idxs]['camera_mat'])
#         cam_mats.append(mat)
#         cam_mats_dist.append(left)

#     cam_mats = arr(cam_mats)
#     cam_mats_dist = arr(cam_mats_dist)

#     out = load_2d_data(config, vid_indices, bp_interested)
    
#     all_points_raw = out['points']
#     all_scores = out['scores']

#     all_points_und = undistort_points(all_points_raw, vid_indices, intrinsics)

#     length = all_points_raw.shape[0]
#     shape = all_points_raw.shape

#     all_points_3d = np.zeros((shape[0], shape[2], 3))
#     all_points_3d.fill(np.nan)

#     errors = np.zeros((shape[0], shape[2]))
#     errors.fill(np.nan)

#     scores_3d = np.zeros((shape[0], shape[2]))
#     scores_3d.fill(np.nan)

#     num_cams = np.zeros((shape[0], shape[2]))
#     num_cams.fill(np.nan)

#     all_points_und[all_scores < reconstruction_threshold] = np.nan

#     for i in trange(all_points_und.shape[0], ncols=70):
#         for j in range(all_points_und.shape[2]):
#             pts = all_points_und[i, :, j, :]
#             good = ~np.isnan(pts[:, 0])
#             if np.sum(good) >= 2:
#                 # TODO: make triangulation type configurable
#                 # p3d = triangulate_optim(pts[good], cam_mats[good])
#                 p3d = triangulate_simple(pts[good], cam_mats[good])
#                 all_points_3d[i, j] = p3d[:3]
#                 errors[i,j] = reprojection_error_und(p3d, pts[good], cam_mats[good], cam_mats_dist[good])
#                 num_cams[i,j] = np.sum(good)
#                 scores_3d[i,j] = np.min(all_scores[i, :, j][good])

#     if 'reference_point' in config['triangulation'] and 'axes' in config['triangulation']:
#         #all_points_3d_adj, recovery = correct_coordinate_frame(config, all_points_3d, bp_interested)
#         all_points_3d_adj = correct_coordinate_frame(config, all_points_3d, bp_interested)
#     else:
#         all_points_3d_adj = all_points_3d
#     dout = pd.DataFrame()
#     for bp_num, bp in enumerate(bp_interested):
#         for ax_num, axis in enumerate(['x','y','z']):
#             dout[bp + '_' + axis] = all_points_3d_adj[:, bp_num, ax_num]
#         dout[bp + '_error'] = errors[:, bp_num]
#         dout[bp + '_ncams'] = num_cams[:, bp_num]
#         dout[bp + '_score'] = scores_3d[:, bp_num]
#     dout['fnum'] = np.arange(length)
#     dout.to_csv(os.path.join(output_path, 'output_3d_data.csv'), index=False)
#     if 'reference_point' in config['triangulation'] and 'axes' in config['triangulation']:
#         return recovery
#         # return None
#     else:
#         return None


def validate_3d(config, **kwargs):
    path, videos, vid_indices = get_video_path(config)
    bp_interested = config['labeling']['bodyparts_interested']
    reconstruction_threshold = config['triangulation']['reconstruction_threshold']
    if config['triangulation'].get('reconstruction_output_path') is None:
        output_path = kwargs.get('output_path', '')
    else:
        output_path = config['triangulation']['reconstruction_output_path']

    try:
        intrinsics = load_intrinsics(path, vid_indices)
    except:
        print("Intrinsic calibration output does not exist.")
        return

    try:
        extrinsics = load_extrinsics(path)
    except:
        print("Extrinsic calibration output does not exist.")
        return

    cam_mats = []
    cam_mats_dist = []

    for vid_idxs in vid_indices:
        mat = arr(extrinsics[vid_idxs])
        left = arr(intrinsics[vid_idxs]['camera_mat'])
        cam_mats.append(mat)
        cam_mats_dist.append(left)

    cam_mats = arr(cam_mats)
    cam_mats_dist = arr(cam_mats_dist)

    out = load_labeled_2d_data(config, vid_indices, bp_interested)
    
    all_points_raw = out['points']

    all_points_und = undistort_points(all_points_raw, vid_indices, intrinsics)

    length = all_points_raw.shape[0]
    shape = all_points_raw.shape

    all_points_3d = np.zeros((shape[0], shape[2], 3))
    all_points_3d.fill(np.nan)

    errors = np.zeros((shape[0], shape[2]))
    errors.fill(np.nan)

    scores_3d = np.zeros((shape[0], shape[2]))
    scores_3d.fill(np.nan)

    num_cams = np.zeros((shape[0], shape[2]))
    num_cams.fill(np.nan)

    for i in trange(all_points_und.shape[0], ncols=70):
        for j in range(all_points_und.shape[2]):
            pts = all_points_und[i, :, j, :]
            good = ~np.isnan(pts[:, 0])
            if np.sum(good) >= 2:
                # TODO: make triangulation type configurable
                # p3d = triangulate_optim(pts[good], cam_mats[good])
                p3d = triangulate_simple(pts[good], cam_mats[good])
                all_points_3d[i, j] = p3d[:3]
                errors[i,j] = reprojection_error_und(p3d, pts[good], cam_mats[good], cam_mats_dist[good])
                num_cams[i,j] = np.sum(good)

    if 'reference_point' in config['triangulation'] and 'axes' in config['triangulation']:
        all_points_3d_adj = correct_coordinate_frame(config, all_points_3d, bp_interested)
    else:
        all_points_3d_adj = all_points_3d

    dout = pd.DataFrame()
    for bp_num, bp in enumerate(bp_interested):
        for ax_num, axis in enumerate(['x','y','z']):
            dout[bp + '_' + axis] = all_points_3d_adj[:, bp_num, ax_num]
        dout[bp + '_error'] = errors[:, bp_num]
        dout[bp + '_ncams'] = num_cams[:, bp_num]

    dout['fnum'] = np.arange(length)

    dout.to_csv(os.path.join(output_path, 'validate_3d_data.csv'), index=False)
