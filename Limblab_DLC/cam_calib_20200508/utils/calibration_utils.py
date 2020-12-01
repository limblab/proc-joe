#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 27 21:24:28 2019

@author: minyoungpark
"""

import os
import re
import cv2
import toml
import fnmatch
import numpy as np
import pandas as pd

from cv2 import aruco
from tqdm import tqdm
from glob import glob
from numpy import array as arr

def get_video_params(vid):
    cap = cv2.VideoCapture(vid)
    
    params = dict()
    params['width'] = int(cap.get(cv2.CAP_PROP_FRAME_WIDTH))
    params['height'] = int(cap.get(cv2.CAP_PROP_FRAME_HEIGHT))
    params['nframes'] = int(cap.get(cv2.CAP_PROP_FRAME_COUNT))
    params['fps'] = cap.get(cv2.CAP_PROP_FPS)
    
    cap.release()
    return params


def get_video_path(config):
    path = config['calibration']['calib_video_path']
    video_extension = config['calibration']['video_extension']
    video_prefix = config['calibration']['calib_video_prefix']
    
    if not os.path.exists(path):
        print('Path does not exist.')

    videos = glob(os.path.join(path, '*' + video_prefix + '*.' +  video_extension))
    videos = sorted(videos)
    
    if not videos:
        print('Calibration videos are not found. Please check if "calib_video_path",' 
              ' "calib_video_prefix", and "video_extension" are correct in the config file.')

    # Trim after the prefix, find numbers and store them into a string inside an array,
    # And then choose the first element.
    vid_indices = [re.findall(r'\d+', vid.split(video_prefix)[1])[0]
                   for vid in videos]
    
    return path, videos, vid_indices


def load_intrinsics(path, vid_indices):
    intrinsics = {}
    for vid_idx in vid_indices:
        intrinsic_path = os.path.join(path, 'intrinsics_{}.toml'.format(vid_idx))
        intrinsics[vid_idx] = toml.load(intrinsic_path)
    return intrinsics


def load_extrinsics(path):
    extrinsics = toml.load(os.path.join(path, 'extrinsics.toml'))
    return extrinsics


ARUCO_DICTS = {
    (4, 50): aruco.DICT_4X4_50,
    (5, 50): aruco.DICT_5X5_50,
    (6, 50): aruco.DICT_6X6_50,
    (7, 50): aruco.DICT_7X7_50,

    (4, 100): aruco.DICT_4X4_100,
    (5, 100): aruco.DICT_5X5_100,
    (6, 100): aruco.DICT_6X6_100,
    (7, 100): aruco.DICT_7X7_100,

    (4, 250): aruco.DICT_4X4_250,
    (5, 250): aruco.DICT_5X5_250,
    (6, 250): aruco.DICT_6X6_250,
    (7, 250): aruco.DICT_7X7_250,

    (4, 1000): aruco.DICT_4X4_1000,
    (5, 1000): aruco.DICT_5X5_1000,
    (6, 1000): aruco.DICT_6X6_1000,
    (7, 1000): aruco.DICT_7X7_1000
}

class Checkerboard:
    def __init__(self, squaresX, squaresY, squareLength):
        self.squaresX = squaresX
        self.squaresY = squaresY
        self.squareLength = squareLength

        objp = np.zeros((squaresX * squaresY, 3), np.float32)
        objp[:, :2] = np.mgrid[0:squaresY, 0:squaresX].T.reshape(-1, 2)
        objp *= squareLength
        self.chessboardCorners = objp
        self.objPoints = objp

    def getChessboardSize(self):
        size = (self.squaresX, self.squaresY)
        return size

    def getGridSize(self):
        return self.getChessboardSize()

    def getSquareLength(self):
        return self.squareLength


def get_calibration_board(config):
    board_size = config['calibration']['board_size']
    board_type = config['calibration']['board_type'].lower()

    if board_type in ['aruco', 'charuco']:
        dkey = (config['calibration']['board_marker_bits'],
                config['calibration']['board_marker_dict_number'])
        dictionary = aruco.getPredefinedDictionary(ARUCO_DICTS[dkey])
        if board_type == 'aruco':
            board = aruco.GridBoard_create(
                board_size[0], board_size[1],
                config['calibration']['board_marker_length'],
                config['calibration']['board_marker_separation_length'],
                dictionary)
        elif board_type == 'charuco':
            board = aruco.CharucoBoard_create(
                board_size[0], board_size[1],
                config['calibration']['board_square_side_length'],
                config['calibration']['board_marker_length'],
                dictionary)
    elif board_type == 'checkerboard':
        board = Checkerboard(board_size[0], board_size[1],
                             config['calibration']['board_square_side_length'])
    else:
        raise ValueError("board_type should be one of "
                         "'aruco', 'charuco', or 'checkerboard' not '{}'".format(
                             board_type))

    return board


def get_board_type(board):
    if isinstance(board, cv2.aruco_GridBoard):
        return 'aruco'
    elif isinstance(board, cv2.aruco_CharucoBoard):
        return 'charuco'
    elif isinstance(board, Checkerboard):
        return 'checkerboard'
    else:
        return None


def get_board_size(board):
    board_type = get_board_type(board)
    if board_type == 'charuco':
        return board.getChessboardSize()
    else:
        return board.getGridSize()


def get_expected_corners(board):
    board_size = get_board_size(board)
    board_type = get_board_type(board)
    if board_type == 'charuco':
        return (board_size[0]-1)*(board_size[1]-1)
    else:
        return board_size[0]*board_size[1]


def get_calibration_board_image(config):
    board = get_calibration_board(config)
    numx, numy = get_board_size(board)
    size = numx*200, numy*200
    img = board.draw(size)
    return img


def undistort_images(config, images_path, cam_num):
    intrinsics = load_intrinsics(config['calibration']['calib_video_path'], cam_num)
    mtx = np.array(intrinsics[cam_num]['camera_mat'])
    dist = np.array(intrinsics[cam_num]['dist_coeff'])
    resolution = tuple(config['video']['resolution'])
    newcameramtx, _ = cv2.getOptimalNewCameraMatrix(mtx, dist, resolution, 1, resolution)
    
    images = []
    for file in os.listdir(images_path):
        if fnmatch.fnmatch(file, '*.png'):
            images.append(file)
    
    path_to_save_undistorted_images = os.path.join(images_path, 'undistorted')
    
    if not os.path.exists(path_to_save_undistorted_images):
        os.mkdir(path_to_save_undistorted_images)
    
    for image in tqdm(images):
        img = cv2.imread(os.path.join(images_path, image))
        dst = cv2.undistort(img, mtx, dist, None, newcameramtx)
        cv2.imwrite(os.path.join(path_to_save_undistorted_images, image), dst)