#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 16 15:42:04 2019

@author: minyoungpark
"""

from mayavi import mlab
mlab.options.offscreen = True



import os
import cv2
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import animation
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Line3DCollection, Poly3DCollection
from tqdm import tqdm, trange
from matplotlib.pyplot import get_cmap
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas

from utils.calibration_utils import get_video_path, get_calibration_board, \
    load_intrinsics, load_extrinsics
from utils.triangulation_utils import load_2d_data
from calibration.extrinsic import detect_aruco
from triangulation.triangulate import triangulate_points, undistort_points


def draw_dropout_histogram(config):
    data_paths = config['paths_to_2d_data']
    bp_interested = config['labeling']['bodyparts_interested']
    
    path, videos, vid_indices = get_video_path(config)
    
    from matplotlib import pyplot as plt
    data = load_2d_data(config, vid_indices, bp_interested)
    l = len(data['scores'])
    fig = plt.figure()
    # print(data['scores'].shape)
    for i in range(4):
        ax = fig.add_subplot(2,2,i+1)
        ax.hist(data['scores'][:,i,:])
        ax.set_title('Cam '+str(i))
        ax.set_xlabel('Likelihood')
        ax.set_ylabel('Percentage')
        ax.set_ylim([0, l])
        ax.set_yticks([0, l/5, 2*l/5, 3*l/5, 4*l/5, l])
        ax.set_yticklabels(['0', '20', '40', '60', '80', '100'])
        ax.legend(config['labeling']['bodyparts_interested'],loc='upper center',prop={'size': 6})
    plt.show()


# Function to extract frames 
def extract_frames(vidpath, times, every_n_frames, path_to_save): 
    if not os.path.exists(vidpath):
        print('Video does not exist.')
        return
    
    if not os.path.exists(path_to_save):
        print(path_to_save + ' does not exsit.')
        folder_input = input('Do you want to create this path (folder)? (y/n) ')
        if folder_input is 'y':
            os.mkdir(path_to_save)
        elif folder_input is 'n':
            return
        else:
            print('Wrong input.')
            return
            
    cap = cv2.VideoCapture(vidpath) 
    fps = cap.get(cv2.CAP_PROP_FPS)
    frame_num = int(cap.get(cv2.CAP_PROP_FRAME_COUNT))
    duration = frame_num/fps
    count = 0
    
    for start_time, end_time in times:
        frame_times = np.arange(start_time, end_time, every_n_frames/fps)
        count += len(frame_times)
        
    with tqdm(total=count) as pbar:
        for start_time, end_time in times:
            frame_times = np.arange(start_time, end_time, every_n_frames/fps)
            frame_counts = list(map(int, frame_times * fps))
            for frame_count in frame_counts:
                cap.set(cv2.CAP_PROP_POS_FRAMES, frame_count)
                ret, frame = cap.read()
                cv2.imwrite(os.path.join(path_to_save, 'img' + str(frame_count).zfill(6) + '.png'), frame)
                pbar.update(1)
            
    print('\n{} frames were extracted.'.format(count))
            

from mpl_toolkits.mplot3d.proj3d import proj_transform
from matplotlib.text import Annotation

class Annotation3D(Annotation):
    '''Annotate the point xyz with text s'''

    def __init__(self, s, xyz, *args, **kwargs):
        Annotation.__init__(self,s, xy=(0,0), *args, **kwargs)
        self._verts3d = xyz        

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.xy=(xs,ys)
        Annotation.draw(self, renderer)
        
        
def annotate3D(ax, s, *args, **kwargs):
    '''add anotation text s to to Axes3d ax'''

    tag = Annotation3D(s, *args, **kwargs)
    ax.add_artist(tag)
    

def transform_corners(points_3d, w, h):
    def proj(u, v):
        """Project u onto v"""
        return u * np.dot(v,u) / np.dot(u,u)
    
    def ortho(u, v):
        """Orthagonalize u with respect to v"""
        return u - proj(v, u)
    
    id0 = points_3d[0]
    idw = points_3d[w-2]
    idh = points_3d[(w-1)*(h-2)]
    
    xAxis = idw - id0
    yAxis = idh - id0
    yAxis_ortho = ortho(yAxis, xAxis)
    
    M = np.zeros((3, 3))
    M[0] = xAxis
    M[1] = yAxis_ortho
    M[2] = np.cross(xAxis, yAxis_ortho)
    M /= np.linalg.norm(M, axis=1)[:,None]
    
    center = np.mean(points_3d, axis=0)
    points_3d_transform = (points_3d - center).dot(M.T)
    
    return points_3d_transform


def plot_poses(config, scale_factor=1):
    '''Creates a plot showing the location and orientation of all cameras.
    Creates a plot showing the location and orientation of all cameras given based on translation
    and rotation vectors. If your cameras are very close together or far apart you can change the
    scaling factor as necessary.
    Arguments:
        pose_estimation_config {dict} -- see help(ncams.camera_tools). Should have following keys:
            serials {list of numbers} -- list of camera serials.
            world_locations {list of np.arrays} -- world locations of each camera.
            world_orientations {list of np.arrays} -- world orientation of each camera.
    '''
    from utils.calibration_utils import load_intrinsics, load_extrinsics
    path = config['calibration']['calib_video_path']
    # intrinscis_dict = load_intrinsics(path, ['1', '2', '3', '4'])
    extrinsics_dict = load_extrinsics(path)
    
    num_cameras = 4
    world_orientations = []
    world_locations = []
    
    for i in range(num_cameras):
        extrinsics = np.array(extrinsics_dict[str(i+1)])
        world_orientations.append(extrinsics[:3, :3])
        world_locations.append(extrinsics[:3, -1])
        
    labels = ['cam_'+str(i) for i in range(num_cameras)]
    # world_locations = pose_estimation_config['world_locations']
    # world_orientations = pose_estimation_config['world_orientations']
    # serials = pose_estimation_config['serials']
    # num_cameras = len(serials)

    # Only accepts list format so check if this is true only when a single camera is present
    # if num_cameras == 1:  # AS: Not sure if needed anymore
    #     if isinstance(world_locations, np.ndarray):
    #         world_locations = [world_locations]
    #     if isinstance(world_orientations, np.ndarray):
    #         world_orientations = [world_orientations]

    # Create a figure with axes
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    plt.setp(ax.get_xticklabels(), visible=False)
    plt.setp(ax.get_yticklabels(), visible=False)
    plt.setp(ax.get_zticklabels(), visible=False)
    
    # Keep the verts for setting the axes later
    cam_verts = [[] for _ in range(num_cameras)]
    for icam in range(num_cameras):
        # Get the vertices to plot appropriate to the translation and rotation
        cam_verts[icam], cam_center = create_camera(
            scale_factor=scale_factor,
            rotation_vector=world_orientations[icam],
            translation_vector=world_locations[icam])

        # Plot it and change the color according to it's number
        ax.add_collection3d(Poly3DCollection(
            cam_verts[icam], facecolors='C'+str(icam), linewidths=1, edgecolors='k', alpha=1))

        # Give each camera a label
        ax.text(np.asscalar(cam_center[0]), np.asscalar(cam_center[1]), np.asscalar(cam_center[2]),
                'cam_' + str(icam))

    # mpl is weird about maintaining aspect ratios so this has to be done
    ax_min = np.min(np.hstack(cam_verts))
    ax_max = np.max(np.hstack(cam_verts))

    # Set the axes and viewing angle
    # Note that this is reversed so that the cameras are looking towards us
    ax.set_xlim([ax_max, ax_min])
    ax.set_ylim([ax_min, ax_max])
    ax.set_zlim([ax_min, ax_max])
    ax.view_init(elev=105, azim=-90)

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')


#################### Camera plotting helper functions
def create_camera(scale_factor=1, rotation_vector=None, translation_vector=None):
        '''Create a typical camera shape.
        [description]
        Keyword Arguments:
            scale_factor {number} -- [description] (default: {1})
            rotation_vector {[type]} -- [description] (default: {None})
            translation_vector {[type]} -- [description] (default: {None})
        Output:
            camera_vertices {np.array} -- [description]
            cam_center {np.array} -- [description]
        '''
        # Lines:
        # Back of camera body
        #  Front of camera body/back of lens
        # Back of camera body
        cam_points = np.array([
            [0, 0, 0],       [1, 0, 0],       [1, 1, 0],       [0, 1, 0],
            [0.2, 0.2, 0.5], [0.8, 0.2, 0.5], [0.8, 0.8, 0.5], [0.2, 0.8, 0.5],
            [0.2, 0.2, 1],   [0.8, 0.2, 1],   [0.8, 0.8, 1],   [0.2, 0.8, 1]])
    
        # Set the origin as the back of the lens
        centering_vector = [0.5, 0.5, 0.5]
        cam_points = cam_points - centering_vector
    
        # Scale the points
        cam_points = cam_points * scale_factor
    
        # Move the camera
        cam_points = move_camera(cam_points, rotation_vector, translation_vector)
    
        # Get the vertices & center
        camera_vertices = get_camera_vertices(cam_points)
        cam_center = np.mean(cam_points[4:8, :], 0)
        cam_center[1] = cam_center[1] + scale_factor
    
        return camera_vertices, cam_center


def move_camera(cam_points, rotation_vector=None, translation_vector=None):
        '''Applies the appropriate rotation and translation to the camera points.
        [description]
        Arguments:
            cam_points {[type]} -- [description]
        Keyword Arguments:
            rotation_vector {np.array} -- [description] (default: {None})
            translation_vector {np.array} -- [description] (default: {None})
        '''
        # Check rotation vector format
        if rotation_vector is None:
            rotation_vector = np.identity(3) # Assume it's not rotating
        elif rotation_vector.shape == (3, 1) or rotation_vector.shape == (1, 3):
            # Make matrix if necessary
            rotation_vector = cv2.Rodrigues(rotation_vector)[0] # Convert to matrix
    
        if translation_vector is None:
            translation_vector = np.zeros((3, 1)) # Assume there is no translation
        elif translation_vector.shape == (1, 3):
            translation_vector = np.transpose(translation_vector) # Format
    
        # Create the translation vector
        translation_vector = np.matmul(-np.transpose(rotation_vector), translation_vector)
    
        # Rotate and then translate
        cam_points = np.transpose(np.matmul(np.transpose(rotation_vector), np.transpose(cam_points)))
        cam_points = cam_points - np.transpose(translation_vector)
    
        return cam_points


def get_camera_vertices(cam_points):
    '''Manual mapping of the camera points from in create_camera.
    [description]
    Arguments:
        cam_points {list} -- 12-element array.
    Output:
        cam_verts {list 9x4} -- [description]
    '''
    cam_verts = [
        [cam_points[0], cam_points[4], cam_points[5], cam_points[1]],
        [cam_points[1], cam_points[5], cam_points[6], cam_points[2]],
        [cam_points[2], cam_points[6], cam_points[7], cam_points[3]],
        [cam_points[3], cam_points[7], cam_points[4], cam_points[0]], # Sides of lenses
        [cam_points[4], cam_points[8], cam_points[9], cam_points[5]],
        [cam_points[5], cam_points[9], cam_points[10], cam_points[6]],
        [cam_points[6], cam_points[10], cam_points[11], cam_points[7]],
        [cam_points[7], cam_points[11], cam_points[8], cam_points[4]],  # Sides of body
        [cam_points[8], cam_points[9], cam_points[10], cam_points[11]]]  # Back of body

    return cam_verts



def reconstruction_validation(config, img_paths):
    images = [cv2.imread(img_path) for img_path in img_paths]
    
    board = get_calibration_board(config)
    (w,h) = board.getChessboardSize()
    num_corners = (w-1)*(h-1)
    
    path, videos, vid_indices = get_video_path(config)
    num_cams = len(vid_indices)
    intrinsics = load_intrinsics(path, vid_indices)
    
    all_detectedCorners = np.zeros((num_corners, num_cams, 2))
    all_detectedCorners.fill(np.nan)
    all_detectedIds = []
    images_with_corners = []
    
    for i, (image, vid_idx) in enumerate(zip(images, vid_indices)):
        gray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
        intrinsic = intrinsics[vid_idx]
        detectedCorners, detectedIds = detect_aruco(gray, intrinsic, board)
        images_with_corners.append(cv2.aruco.drawDetectedCornersCharuco(image, 
                                                     detectedCorners,
                                                     detectedIds))
        all_detectedIds.append(detectedIds)
        
        for coord, j in zip(detectedCorners, detectedIds):
            all_detectedCorners[j[0]][i] = coord
    
    extrinsics = load_extrinsics(path)
    
    cam_mats = np.array([extrinsics[vid_idx] for vid_idx in vid_indices])
    
    undistorted_corners = undistort_points(all_detectedCorners, vid_indices, intrinsics)
    
    points_3d, _ = triangulate_points(undistorted_corners, cam_mats)
    points_3d = points_3d[:, :3]
    
    points_3d_transform = transform_corners(points_3d, w, h)
    
    square_length = config['calibration']['board_square_side_length']
    
    x = np.linspace((-w/2+1)*square_length, (w/2-1)*square_length, w-1)
    y = np.linspace((-h/2+1)*square_length, (h/2-1)*square_length, h-1)
    xv, yv = np.meshgrid(x, y)
    
    pseudo_real_corners = np.zeros((num_corners, 3))
    pseudo_real_corners[:, 0] = np.reshape(xv, (-1))
    pseudo_real_corners[:, 1] = np.reshape(yv, (-1))

    subplot_w = int(num_cams // np.floor(np.sqrt(num_cams)))
    subplot_h = num_cams / subplot_w
    
    fig = plt.figure()
    
    for i in range(num_cams):
        ax = fig.add_subplot(subplot_w, subplot_h, i+1)
        ax.imshow(images_with_corners[i])
        ax.axis('off')
    
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    
    labels = ['id_'+str(i) for i in range(num_corners)]
    
    for i, (pred, real) in enumerate(zip(points_3d_transform, pseudo_real_corners)):
        ax.scatter(pred[0], pred[1], pred[2], c='r', alpha=0.5)
        ax.scatter(real[0], real[1], real[2], c='b', alpha=0.5)
        annotate3D(ax, s=labels[i], xyz=pred, fontsize=8, xytext=(-3,3),
               textcoords='offset points', ha='center',va='bottom')
    
    ax_min, ax_max = -np.floor(max(w, h)*square_length)//2, np.floor(max(w, h)*square_length//2)
    ax.set_xlim([ax_min, ax_max])
    ax.set_ylim([ax_min, ax_max])
    ax.set_zlim([ax_min, ax_max])
    ax.view_init(elev=50, azim=-70)
    
    ax.set_xlabel('x (mm)')
    ax.set_ylabel('y (mm)')
    ax.set_zlabel('z (mm)')    
        
    return points_3d, points_3d_transform


def generate_three_dim_video(config, **kwargs):
#def generate_three_dim_video(config, azimuth, elevation **kwargs):
   
    path, videos, vid_indices = get_video_path(config)
    if config.get('output_video_path') is None:
        output_path = kwargs.get('output_path', '')
    else:
        output_path = config['output_video_path']
    
    try:
        connections = config['labeling']['scheme']
    except KeyError:
        connections = []
        
    bp_interested = config['labeling']['bodyparts_interested']
    bp_dict = dict(zip(bp_interested, range(len(bp_interested))))
    
    if config['triangulation'].get('reconstruction_output_path') is None:
        reconstruction_output_path = kwargs.get('output_path', '')
    else:
        reconstruction_output_path = config['triangulation']['reconstruction_output_path']
    
    reconstruction_filename = 'output_3d_data.csv'
    data = pd.read_csv(os.path.join(reconstruction_output_path, reconstruction_filename))
    
    all_points = np.array([np.array(data.loc[:, (bp+'_x', bp+'_y', bp+'_z')])
                            for bp in bp_interested])
    
    all_errors = np.array([np.array(data.loc[:, bp+'_error'])
                            for bp in bp_interested])
    
    
    all_errors[np.isnan(all_errors)] = 10000
    
    good = (all_errors < 150)
    all_points[~good] = np.nan
    
    all_points_flat = all_points.reshape(-1, 3)
    check = ~np.isnan(all_points_flat[:, 0])
    low, high = np.percentile(all_points_flat[check], [5, 95], axis=0)
    
    fps = config['video']['fps']
    resolution = config['video']['resolution']
    resolution = (640, 480)
    
    fourcc = cv2.VideoWriter_fourcc(*'mp4v')
    three_dim_video_filename = 'output_3d_video.mp4'
    writer = cv2.VideoWriter(os.path.join(output_path, three_dim_video_filename), 
                              fourcc, fps, resolution)
    
    # cmap = plt.get_cmap('tab10')
    
    all_x_points = all_points[:,:,0]
    all_y_points = all_points[:,:,1]
    all_z_points = all_points[:,:,2]
    x_low, x_high = np.percentile(all_x_points[np.isfinite(all_x_points)], [3, 97])
    y_low, y_high = np.percentile(all_y_points[np.isfinite(all_y_points)], [3, 97])
    z_low, z_high = np.percentile(all_z_points[np.isfinite(all_z_points)], [3, 97])
    
    fig = plt.figure()
    ax = Axes3D(fig)
    
    all_points = np.transpose(all_points, (1,0,2))
    ax.view_init(azim=-48, elev=28)
    
    
    def draw_lines(points):
        ax.set_xlim([x_low-100, x_high+100])
        ax.set_ylim([y_low-100, y_high+100])
        ax.set_zlim([z_low-100, z_high+100])
        ax.set_xlabel('x (mm)')
        ax.set_ylabel('y (mm)')
        ax.set_zlabel('z (mm)')    
        ax.scatter(points[:, 0], points[:, 1], points[:, 2], 'bo', alpha=0.5)
        for connection in connections:
            xs = [points[bp_dict[connection[0]]][0], points[bp_dict[connection[1]]][0]]
            ys = [points[bp_dict[connection[0]]][1], points[bp_dict[connection[1]]][1]]
            zs = [points[bp_dict[connection[0]]][2], points[bp_dict[connection[1]]][2]]
            ax.plot(xs, ys, zs, 'g-')
    
    # points = all_points[0]
    # draw_lines(points)
    # print('ax.azim {}'.format(ax.azim))
    # print('ax.elev {}'.format(ax.elev))
    # # For testing
    # all_points = all_points[:600]
    canvas = FigureCanvas(fig)
    for i in trange(len(all_points), ncols=70):
        ax.cla()
        points = all_points[i]
        draw_lines(points)
        canvas.draw()
        img = np.array(canvas.renderer._renderer)
        img = cv2.cvtColor(img, cv2.COLOR_RGB2BGR)
        writer.write(img)
    plt.close()
    writer.release()


def generate_three_dim_pictures(config, **kwargs):
   
    path, videos, vid_indices = get_video_path(config)
    if config.get('output_video_path') is None:
        output_path = kwargs.get('output_path', '')
    else:
        output_path = config['output_video_path']
    
    try:
        connections = config['labeling']['scheme']
    except KeyError:
        connections = []
        
    bp_interested = config['labeling']['bodyparts_interested']
    bp_dict = dict(zip(bp_interested, range(len(bp_interested))))
    
    if config['triangulation'].get('reconstruction_output_path') is None:
        reconstruction_output_path = kwargs.get('output_path', '')
    else:
        reconstruction_output_path = config['triangulation']['reconstruction_output_path']
    
    reconstruction_filename = 'output_3d_data.csv'
    data = pd.read_csv(os.path.join(reconstruction_output_path, reconstruction_filename))
    
    all_points = np.array([np.array(data.loc[:, (bp+'_x', bp+'_y', bp+'_z')])
                            for bp in bp_interested])
    
    all_errors = np.array([np.array(data.loc[:, bp+'_error'])
                            for bp in bp_interested])
    
    
    all_errors[np.isnan(all_errors)] = 10000
    
    good = (all_errors < 150)
    all_points[~good] = np.nan
    
    all_points_flat = all_points.reshape(-1, 3)
    check = ~np.isnan(all_points_flat[:, 0])
    low, high = np.percentile(all_points_flat[check], [5, 95], axis=0)
    
    fps = config['video']['fps']
    resolution = config['video']['resolution']
    resolution = (640, 480)
    
    fourcc = cv2.VideoWriter_fourcc(*'mp4v')
    three_dim_video_filename = 'output_3d_video.mp4'
    writer = cv2.VideoWriter(os.path.join(output_path, three_dim_video_filename), 
                              fourcc, fps, resolution)
    
    # cmap = plt.get_cmap('tab10')
    
    all_x_points = all_points[:,:,0]
    all_y_points = all_points[:,:,1]
    all_z_points = all_points[:,:,2]
    x_low, x_high = np.percentile(all_x_points[np.isfinite(all_x_points)], [3, 97])
    y_low, y_high = np.percentile(all_y_points[np.isfinite(all_y_points)], [3, 97])
    z_low, z_high = np.percentile(all_z_points[np.isfinite(all_z_points)], [3, 97])
    
    fig = plt.figure()
    ax = Axes3D(fig)
    
    all_points = np.transpose(all_points, (1,0,2))
    ax.view_init(azim=-48, elev=28)
    
    
    def draw_lines(points):
        ax.set_xlim([x_low-100, x_high+100])
        ax.set_ylim([y_low-100, y_high+100])
        ax.set_zlim([z_low-100, z_high+100])
        ax.set_xlabel('x (mm)')
        ax.set_ylabel('y (mm)')
        ax.set_zlabel('z (mm)')    
        ax.scatter(points[:, 0], points[:, 1], points[:, 2], 'bo', alpha=0.5)
        for connection in connections:
            xs = [points[bp_dict[connection[0]]][0], points[bp_dict[connection[1]]][0]]
            ys = [points[bp_dict[connection[0]]][1], points[bp_dict[connection[1]]][1]]
            zs = [points[bp_dict[connection[0]]][2], points[bp_dict[connection[1]]][2]]
            ax.plot(xs, ys, zs, 'g-')
    
    points = all_points[1000]
    draw_lines(points)
            
            
    
# def generate_valid_three_dim_video(config, **kwargs):
#     path, videos, vid_indices = get_video_path(config)
#     if config.get('output_video_path') is None:
#         output_path = kwargs.get('output_path', '')
#     else:
#         output_path = config['output_video_path']

#     try:
#         scheme = config['labeling']['scheme']
#     except KeyError:
#         scheme = []
        
#     if config['triangulation'].get('reconstruction_output_path') is None:
#         reconstruction_output_path = kwargs.get('output_path', '')
#     else:
#         reconstruction_output_path = config['triangulation']['reconstruction_output_path']
    
#     reconstruction_filename = 'validate_3d_data.csv'
#     data = pd.read_csv(os.path.join(reconstruction_output_path, reconstruction_filename))
#     cols = [x for x in data.columns if '_error' in x]

#     if len(scheme) == 0:
#         bodyparts = [c.replace('_error', '') for c in cols]
#     else:
#         bodyparts = sorted(set([x for dx in scheme for x in dx]))

#     bp_dict = dict(zip(bodyparts, range(len(bodyparts))))

#     all_points = np.array([np.array(data.loc[:, (bp+'_x', bp+'_y', bp+'_z')])
#                             for bp in bodyparts])

#     all_errors = np.array([np.array(data.loc[:, bp+'_error'])
#                             for bp in bodyparts])

#     # all_scores = np.array([np.array(data.loc[:, bp+'_score'])
#     #                         for bp in bodyparts])

#     # if config['triangulation']['optim']:
#     #     all_errors[np.isnan(all_errors)] = 0
#     # else:
#     all_errors[np.isnan(all_errors)] = 10000
    
# #    good = (all_errors < 100)
#     good = (all_errors < 150)
#     all_points[~good] = np.nan

#     all_points_flat = all_points.reshape(-1, 3)
#     check = ~np.isnan(all_points_flat[:, 0])
#     low, high = np.percentile(all_points_flat[check], [5, 95], axis=0)

#     nparts = len(bodyparts)
#     framedict = dict(zip(data['fnum'], data.index))

# #    writer = skvideo.io.FFmpegWriter(outname, inputdict={
# #        # '-hwaccel': 'auto',
# #        '-framerate': str(fps),
# #    }, outputdict={
# #        '-vcodec': 'h264', '-qp': '30'
# #    })
    
#     fps = config['video']['fps']
#     resolution = config['video']['resolution']
#     resolution = (resolution[0], resolution[1])
    
#     # OpenCV vid writer
#     fourcc = cv2.VideoWriter_fourcc(*'XVID')
#     three_dim_video_filename = 'output_3d_valid_video.mp4'
#     writer = cv2.VideoWriter(os.path.join(output_path, three_dim_video_filename), 
#                               fourcc, fps, resolution)
    
#     cmap = get_cmap('tab10')


#     points = all_points[:, 20]
#     points[0] = low
#     points[1] = high

#     s = np.arange(points.shape[0])
#     good = ~np.isnan(points[:, 0])

#     fig = mlab.figure(bgcolor=(1,1,1), size=resolution)
#     fig.scene.anti_aliasing_frames = 2

#     low, high = np.percentile(points[good, 0], [10,90])
#     scale_factor = (high - low) / 10.0

#     mlab.clf()
#     mlab.axes(ranges=(
#             np.min(points[:, 0])-100,
#             np.max(points[:, 0])+100,
            
#             np.min(points[:, 1])-100,
#             np.max(points[:, 1])+100,
            
#             np.min(points[:, 2])-100,
#             np.max(points[:, 2])+100,
#             ))
# #    pts = mlab.points3d(points[:, 0], -points[:, 1], points[:, 2], s,
# #                        scale_mode='none', scale_factor=scale_factor)
#     pts = mlab.points3d(points[:, 0], -points[:, 1], points[:, 2], s,
#                         scale_mode='none', scale_factor=1)
#     lines = connect_all(points, scheme, bp_dict, cmap)
#     mlab.orientation_axes()

#     view = list(mlab.view())

#     mlab.view(focalpoint='auto', distance='auto')

# #    for framenum in trange(30, ncols=70):
#     for framenum in trange(data.shape[0], ncols=70):
#         fig.scene.disable_render = True

#         if framenum in framedict:
#             points = all_points[:, framenum]
#         else:
#             points = np.ones((nparts, 3))*np.nan

#         s = np.arange(points.shape[0])
#         good = ~np.isnan(points[:, 0])

#         new = np.vstack([points[:, 0], -points[:, 1], points[:, 2]]).T
#         pts.mlab_source.points = new
#         update_all_lines(lines, points, scheme, bp_dict)

#         fig.scene.disable_render = False


#         # mlab.view(*view, reset_roll=False)
#         mlab.view(120, 0)
# #        writer.writeFrame(img)
#         img = mlab.screenshot()
        
#         #OpenCV
#         writer.write(img)

#     mlab.close(all=True)
# #    writer.close()
