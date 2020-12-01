#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 07:02:53 2020

@author: minyoungpark
"""

#%% Filtering
import pandas as pd
import os
import numpy as np
from scipy import interpolate

def interpolate_low_prob(x, prob):
    x_val = np.array([])
    idx_high_prob = np.array([])
    idx_to_interp = np.array([])
    for i in range(len(x)):
        if prob[i] > 0.7:
            x_val = np.append(x_val, x[i])
            idx_high_prob = np.append(idx_high_prob, i)
        else:
            idx_to_interp = np.append(idx_to_interp, i)

    x_interp = interpolate.griddata(np.expand_dims(idx_high_prob, axis=1), x_val,
                                np.expand_dims(idx_to_interp, axis=1), method='nearest')
        
    return idx_to_interp, x_interp

#x[ii.astype(int)] = np.squeeze(xx)

from scipy import signal

video_names = ['cam_0', 'cam_1', 'cam_2', 'cam_3']
joints = ['Wrist', 'CMC_thumb', 'MCP_thumb', 'MCP1', 'MCP2', 'MCP3', 'MCP4',
          'IP_thumb', 'PIP1', 'PIP2', 'PIP3', 'PIP4', 'Dip1', 'Dip2', 'Dip3', 'Dip4',
          'Tip_thumb', 'Tip1', 'Tip2', 'Tip3', 'Tip4']
video_folder_path = '/media/minyoungpark/Min/pop_0317_2/dlc/3d_filt'

snapshots = ['DLC_resnet50_Pop_freeReach_0317_cam_0Mar23shuffle1_1030000',
             'DLC_resnet50_Pop_freeReach_0317_cam_1Mar19shuffle1_1030000',
             'DLC_resnet50_Pop_freeReach_0317_cam_2Mar23shuffle1_950000',
             'DLC_resnet50_Pop_freeReach_0317_cam_3Mar23shuffle1_1030000']

data = []
for i in range(4):
    data.append(pd.read_hdf(os.path.join(video_folder_path , video_names[i] +'.h5'), 'df_with_missing'))
# data = pd.read_hdf(os.path.join(video_folder_path , video_names[i] +'.h5'), 'df_with_missing')
# for video_name, snapshot in zip(video_names, snapshots):
#     data = pd.read_hdf(os.path.join(video_folder_path, video_name+snapshot+'.h5'), 
#                        'df_with_missing')
#     prefilt = pd.read_csv(video_folder_path + video_name + '_prefiltered.csv')
# #    coord = pd.read_hdf(project_name+'videos/'+video_name+snapshot+'.h5')
#     for joint in joints:
#     #    ii, xx = interpolate_low_prob(coord[snapshot, joint, 'x'], 
#     #                                  coord[snapshot, joint, 'likelihood'])
#     #    coord[snapshot, joint, 'x'][ii.astype(int)]
# #        xmed = signal.medfilt(coord[snapshot, joint, 'x'], 13)
# #        ymed = signal.medfilt(coord[snapshot, joint, 'y'], 13)
# #        for i in range(len(coord)):
# #            if coord[snapshot, joint, 'likelihood'][i] < 0.25:
# #                coord[snapshot, joint, 'x'][i] = xmed[i]
# #                coord[snapshot, joint, 'y'][i] = ymed[i]
#         coord[snapshot, joint, 'x'] = signal.savgol_filter(prefilt[joint + '_x'], 13, 2)
#         coord[snapshot, joint, 'y'] = signal.savgol_filter(prefilt[joint + '_y'], 13, 2)
        
#     coord.to_hdf(video_folder_path+video_name+snapshot+'.h5', key='df_with_missing', mode='w')

#%%
import numpy as np
import matplotlib.pyplot as plt
from numpy import array as arr

def kalman_xy(x, P, measurement, R, dt, segmaux, segmauy):
    """
    Parameters:    
    x: initial state 4-tuple of location and velocity: (x0, x1, x0_dot, x1_dot)
    P: initial uncertainty convariance matrix
    measurement: observed position
    R: measurement noise 
    motion: external motion added to state vector x
    Q: motion noise (same shape as P)
    """
    F = np.matrix([[1, 0, dt, 0],
                   [0, 1, 0, dt],
                   [0, 0, 1, 0],
                   [0, 0, 0, 1]])
    
    H = np.matrix([[1, 0, 0, 0],
                   [0, 1, 0, 0]])
    motion = np.matrix([[0.5*dt**2, 0.5*dt**2, dt, dt]]).T # G
    # Q = np.matrix(np.eye(4))
    # Q = np.matrix([[0.25*dt**4, 0.25*dt**4 , 0.5*dt**3, 0.5*dt**3 ],
    #                 [0.25*dt**4 , 0.25*dt**4.  , 0.5*dt**3 , 0.5*dt**3  ],
    #                 [0.5*dt**3, 0.5*dt**3 , 1*dt**2, 1*dt**2 ],
    #                 [0.5*dt**3 , 0.5*dt**3.  , 1*dt**2 , 1*dt**2  ]])
    Q = np.matrix([[segmaux**2, 0         ],
                   [0,          segmauy**2]])
    
    B = np.matrix([[motion.item(0), 0        ], # define the input matrix
                  [0,         motion.item(1)],
                  [motion.item(2),        0        ],
                  [0,         motion.item(3)       ]])
    
    return kalman(x, P, measurement, R, motion, Q, F, H, B)

def kalman(x, P, measurement, R, motion, Q, F, H, B):
    '''
    Parameters:
    x: initial state
    P: initial uncertainty convariance matrix
    measurement: observed position (same shape as H*x)
    R: measurement noise (same shape as H)
    motion: external motion added to state vector x
    Q: motion noise (same shape as P)
    F: next state function: x_prime = F*x
    H: measurement function: position = H*x

    Return: the updated and predicted new values for (x, P)

    See also http://en.wikipedia.org/wiki/Kalman_filter

    This version of kalman can be applied to many different situations by
    appropriately defining F and H 
    '''
    # UPDATE x, P based on measurement m    
    # distance between measured and current position-belief
    y = np.matrix(measurement).T - H * x
    S = H * P * H.T + R  # residual convariance
    K = P * H.T * S.I    # Kalman gain
    x = x + K*y
    I = np.matrix(np.eye(F.shape[0])) # identity matrix
    P = (I - K*H)*P

    # PREDICT x, P based on motion
    x = F*x + motion
    # P = F*P*F.T + Q
    P = F*P*F.T + B*Q*B.T

    return x, P


# for idx, joint in enumerate(joints):
coord = arr(data[0][snapshots[0]]['Wrist'])[:200, :]
x = np.expand_dims(np.array([coord[0, 0], coord[0, 1], 0, 0]), axis=-1)
P = np.matrix(np.eye(4))*0.01 # initial uncertainty
# N = len(coord)
observed_x = coord[:, 0]
observed_y = coord[:, 1]
result = []
R = 0.0001**2
dt = 1/30

segmaux = 7 # standard deviation ax
segmauy = 4 # standard deviation ay


for meas in zip(observed_x, observed_y):
    x, P = kalman_xy(x, P, meas, R, dt, segmaux, segmauy)
    result.append((x[:2]).tolist())
kalman_x, kalman_y = zip(*result)

from scipy.signal import butter,filtfilt, freqz, freqs

n = len(coord)  # total number of samples
fs = 30.0       # sample rate, Hz
T = n/fs         # Sample Period
cutoff = 10      # desired cutoff frequency of the filter, Hz, slightly higher than actual 1.2 Hz
nyq = 0.5 * fs  # Nyquist Frequency
order = 5       # sin wave can be approx represented as quadratic

def butter_lowpass_filter(data, cutoff, fs, order):
    normal_cutoff = cutoff / nyq
    # Get the filter coefficients 
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    y = filtfilt(b, a, data)
    output = {'y': y,
              'b': b,
              'a': a}
    return output

data_3d_path = '/media/minyoungpark/Min/pop_0317_2/dlc/3d_filt/output_3d_data.csv'
data_3d = pd.read_csv(data_3d_path)
coord_3d = arr([data_3d['Dip1_x'], data_3d['Dip1_y'], data_3d['Dip1_z']])

def lpf_3d(data, cutoff, fs, order):
    normal_cutoff = cutoff / nyq
    # Get the filter coefficients 
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    y = filtfilt(b, a, data)
    output = {'y': y,
              'b': b,
              'a': a}
    return output

lpf_3d = lpf_3d(coord_3d, cutoff, fs, order)

lpf_x = butter_lowpass_filter(observed_x, cutoff, fs, order)
lpf_y = butter_lowpass_filter(observed_y, cutoff, fs, order)


#%%
x = np.expand_dims(np.array([lpf_x['y'][0], lpf_y['y'][0], 0, 0]), axis=-1)
P = np.matrix(np.eye(4))*0.01 # initial uncertainty
result = []
R = 0.1**2
for meas in zip(lpf_x['y'], lpf_y['y']):
    x, P = kalman_xy(x, P, meas, R, dt, segmaux, segmauy)
    result.append((x[:2]).tolist())
lpf_kalman_x, lpf_kalman_y = zip(*result)

plt.plot(observed_x, observed_y, 'ko')
plt.plot(lpf_x['y'], lpf_y['y'], 'g-')
plt.plot(kalman_x, kalman_y, 'b-')
# plt.plot(lpf_kalman_x, lpf_kalman_y, 'm-')
plt.legend(['Observed', 'LP filtered', 'Kalman filtered', 'LPF + Kalman'])
plt.title('Trajectory of the first 200 frames')
plt.show()

#%%
w_x, h_x = freqz(lpf_x['b'], lpf_x['a'])
w_y, h_y = freqz(lpf_y['b'], lpf_y['a'])
plt.semilogx((fs * 0.5 / np.pi)*w_x, 20 * np.log10(abs(h_x)))
plt.title('Butterworth filter frequency response')
plt.xlabel('Frequency [Hz]')
plt.ylabel('Amplitude [dB]')
plt.margins(0, 0.1)
plt.grid(which='both', axis='both')
plt.axvline(10, color='green') # cutoff frequency
plt.show()

#%%
from scipy import signal
joints = ['Wrist', 'CMC_thumb', 'MCP_thumb', 'MCP1', 'MCP2', 'MCP3', 'MCP4',
          'IP_thumb', 'PIP1', 'PIP2', 'PIP3', 'PIP4', 'Dip1', 'Dip2', 'Dip3', 'Dip4',
          'Tip_thumb', 'Tip1', 'Tip2', 'Tip3', 'Tip4']

fs = 30
fig = plt.figure()

for idx, joint in enumerate(joints):
    ax = fig.add_subplot(3, 7, idx+1)
    
    coord = arr(data[snapshots[0]][joint])[:,:]
    observed_x = coord[:, 0]
    observed_y = coord[:, 1]
    
    f_x, Pxx_spec_x = signal.welch(observed_x, fs, 'flattop', 1024, scaling='spectrum')
    f_y, Pxx_spec_y = signal.welch(observed_y, fs, 'flattop', 1024, scaling='spectrum')

    ax.semilogy(f_x, np.sqrt(Pxx_spec_x))
    ax.semilogy(f_y, np.sqrt(Pxx_spec_y))
    ax.legend(['x', 'y'])
    ax.set_xlabel('frequency [Hz]')
    ax.set_ylabel('Linear spectrum [V RMS]')
    ax.set_title(joint)
    
fig.subplots_adjust(hspace=.5, wspace=.5)
plt.show()

#%%

plt.show()




#%%
import numpy as np
from numpy.random import normal as normrnd
from control.matlab import ss, lsim

np.seterr(divide='ignore', invalid='ignore')

def my_kalman_filter(data, Ts, stds, accels):

    # Ts = 0.1; # define the sample time
    
    # A = np.matrix([[1, 0, 0, Ts, 0,  0], # define the state matrix
    #                 [0, 1, 0, 0,  Ts, 0],
    #                 [0, 0, 1, 0,  0, Ts], 
    #                 [0, 0, 0, 1,  0,  0],
    #                 [0, 0, 0, 0,  1,  0],
    #                 [0, 0, 0, 0,  0,  1]]) 
    A = np.matrix([[1, 0, dt, 0],
                   [0, 1, 0,  dt],
                   [0, 0, 1,  0],
                   [0, 0, 0,  1]])
    # A = np.matrix([[1, 0, 0, 0, Ts, 0,  0], # define the state matrix
    #                [0, 1, 0, 0,  Ts, 0],
    #                [0, 0, 1, 0,  0, Ts], 
    #                [0, 0, 0, 1,  0,  0],
    #                [0, 0, 0, 0,  1,  0],
    #                [0, 0, 0, 0,  0,  1]]) 
    
    # C = np.matrix([[1, 0, 0, 0, 0, 0], # define the output matrix
    #                [0, 1, 0, 0, 0, 0],
    #                [0, 0, 1, 0, 0, 0]])
    C = np.matrix([[1, 0, 0, 0],
                   [0, 1, 0, 0]])
    
    # B = np.matrix([[0.5*Ts**2, 0,         0], # define the input matrix
    #                [0,         0.5*Ts**2, 0],
    #                [0,         0,         0.5*Ts**2],
    #                [Ts,        0,         0],
    #                [0,         Ts,        0],
    #                [0,         0,         Ts]])
    B = np.matrix([[0.5*Ts**2, 0        ], # define the input matrix
                   [0,         0.5*Ts**2],
                   [Ts,        0        ],
                   [0,         Ts       ]])
    
    # x0 = np.matrix([0, 0, 0, 0, 0, 0]).T # define the initial conditions
    x0 = np.matrix([0, 0, 0, 0]).T # define the initial conditions
    # sys = ss(A, B, np.eye(6), np.zeros(B.shape), Ts) # define a system to generate true data
    sys = ss(A, B, np.eye(4), np.zeros(B.shape), Ts) # define a system to generate true data
    # t = np.arange(0, 40 ,Ts) # define the time interval
    t = Ts*np.arange(len(data))
    
    # assuming that the uncertanities in the accelerations are equal, we define
    # them as follow:
    segmaux = 5 # standard deviation ax
    segmauy = 5 # standard deviation ay
    segmaualpha = 5 # standard deviation angular acceleration
    
    # In practice, these values are determined experimentally.
    # define the input(accelerations):
    ux = np.concatenate([np.zeros((1, 30)), 
                        25*np.ones((1, 20)),
                        -20*np.ones((1, 20)), 
                        15*np.ones((1, len(t)-70))], axis=1) + normrnd(0, segmaux, (1, len(t)))
    
    uy = np.concatenate([np.zeros((1, 10)),
                         60*np.ones((1,60)),
                         -20*np.ones((1,len(t)-70))], axis=1) + normrnd(0, segmauy, (1, len(t)))
    ualpha = np.concatenate([np.zeros((1, 30)),
                              25*np.ones((1, 20)),
                              -20*np.ones((1, 20)),
                              15*np.ones((1,len(t)-70))], axis=1) + normrnd(0, segmaualpha, (1, len(t)))
    u = np.concatenate([ux, uy, ualpha], axis=0)
    # generating the true data:
    Xtrue = lsim(sys, u, t, x0)
    xtrue = Xtrue[0][:, 0]
    ytrue = Xtrue[0][:, 1]
    thtrue = Xtrue[0][:, 2]
    vxtrue = Xtrue[0][:, 3]
    vytrue = Xtrue[0][:, 4]
    wtrue = Xtrue[0][:, 5]
    # defining V:
    # measurmentsV = np.matrix([[200**2, 0,      0], 
    #                           [0,      200**2, 0],
    #                           [0,      0,      300**2]])
    measurmentsV = np.matrix([[stds[0]**2, 0     ], 
                              [0,      stds[1]**2]])
    # generating measurment data by adding noise to the true data:
    # xm = xtrue + normrnd(0, 200, (len(xtrue),))
    # ym = ytrue + normrnd(0, 200, (len(ytrue),))
    # thm = thtrue+normrnd(0, 300, (len(ytrue),))
    
    xm = data[:, 0]
    ym = data[:, 1]
    # thm = data[:, 2]
    
    # initializing the matricies for the for loop (this will make the matlab run
    # the for loop faster.
    # Xest = np.zeros((6, len(t)))
    Xest = np.zeros((4, len(t)))
    Xest[:, 0] = x0.T
    # defining R and Q
    R = measurmentsV * C * C.T
    Q = np.matrix([[segmaux**2, 0,          0],
                   [0,          segmauy**2, 0], 
                   [0,          0,          segmaualpha**2]])
    # Initializing P 
    P = B * Q * B.T
    
    for i in range(1, len(t)):
        P = A * P * A.T + B * Q * B.T # predicting P
        Xest[:, i] = np.squeeze(A * np.expand_dims(Xest[:, i-1], axis=-1) + \
        B * np.expand_dims(u[:, i-1], axis=-1)) # Predicitng the state
        # num = P * C.T
        # den = np.tile((C * P * C.T + R), (2,1))
        # K = num/den # calculating the Kalman gains
        K = P * C.T * (C * P * C.T + R).I
        K = np.nan_to_num(K)
        # Xest[:, i] = Xest[:, i] + np.squeeze(K * (np.matrix([xm[i], ym[i], thm[i]]).T - 
        #                                C * np.expand_dims(Xest[:, i], axis=-1))) # Correcting: estimating the state
        
        Xest[:, i] = Xest[:, i] + np.squeeze(K * (np.matrix([xm[i], ym[i]]).T - 
                                       C * np.expand_dims(Xest[:, i], axis=-1))) # Correcting: estimating the state
        # P = (np.eye(6) - K * C) * P # Correcting: estimating P
        P = (np.eye(4) - K * C) * P # Correcting: estimating P



from numpy import array as arr
coords = []
for i in range(4):
    coords.append(arr(data[i][snapshots[i]]['MCP1'])[:200,:2])
# coord = arr(data[snapshots[0]]['MCP1'])[:200,:]
# x = np.expand_dims(np.array([coord[0, 0], coord[0, 1], 0, 0]), axis=-1)
#%%
from matplotlib import pyplot as plt
fig = plt.figure()

ax1 = fig.add_subplot(311)
ax1.plot(t,Xest[0,:],'r',t,xm,'g',t,xtrue,'b')
ax1.set_xlabel('time [sec]');
ax1.set_ylabel('displacementx [m/s]');
ax1.set_title('displacementx');
# ax1.legend('estimated displacementx','measured displacementx','true displacementx');
ax2 = fig.add_subplot(312)
ax2.plot(t,Xest[1,:],'r',t,ym,'g',t,ytrue,'b')
ax2.set_xlabel('time [sec]');
ax2.set_ylabel('displacementy [m/s]');
ax2.set_title('displacementy');
# legend('estimated displacementy','measured displacementy','true displacementy');
# t=0:0.1:40;
ax3 = fig.add_subplot(313)
ax3.plot(t,Xest[2,:],'r',t,thm,'g',t,thtrue,'b')
ax3.set_xlabel('time [sec]');
ax3.set_ylabel('angle');
ax3.set_title('angle theta');
# legend('estimated angle theta','measured angle theta','true angle theta');
# t=0:0.1:40;
# figure
# hold on 
# %simple animation:
# for i=1:1:length(t)
# axis([min(xtrue)-500 max(xtrue)+500 min(ytrue)-500 max(ytrue)+500]);
# %viscircles([xtrue(i) ytrue(i)],20,'color','b')
# %viscircles([Xest(1,i) Xest(2,i)],20,'color','r')
# plot(xtrue(i),ytrue(i),'bo');
# plot(Xest(1,i),Xest(2,i),'rx');
# pause(0.1)
# end
