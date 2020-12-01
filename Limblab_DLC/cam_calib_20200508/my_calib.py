#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 25 14:31:01 2019

@author: minyoungpark
"""


from utils.utils import load_config
from calibration.intrinsic import calibrate_intrinsic
from calibration.extrinsic import calibrate_extrinsic
from triangulation.triangulate import reconstruct_3d
from utils.vis_utils import generate_three_dim_video
from utils.vis_utils import generate_three_dim_pictures

#%%
#config = load_config('config_20200804_FR.toml' )
config = load_config('config_Test_20201123_15in.toml' )

#%% If you already ran calibration you don't need to run these.
#%%
calibrate_intrinsic(config)
#%%
calibrate_extrinsic(config)

#%%
from utils.triangulation_utils import add_static_points
    
import numpy as np

labels = ['pointX', 'pointY', 'pointZ']

snapshots = ['Han_20201123_test_15in00001DLC_resnet50_TestNov23shuffle1_1030000filtered.csv',
             'Han_20201123_test_15in00001DLC_resnet50_TestNov23shuffle1_1030000filtered.csv']

#2020.7.2
#      POINT X-X          X-Y         Y-X        Y-Y         Z-X            Z-Y   
#cam1: 450.4066209	426.7973505	482.023147	339.8519037	401.0057988	120.5122537
#cam2: 438.5002129	383.2673147	639.7415945	287.6776584	393.220902	27.74087388
#cam3: 507.2576849	280.9696124	666.5737787	463.7638673	534.0898691	37.80294296
#cam4: 554.2140073	410.0994989	589.4312491	564.3845581	544.1519382	225.6282324
#static = {'pointX': [[450.4066209, 426.7973505],[438.5002129, 383.2673147],[507.2576849, 280.9696124],[554.2140073, 410.0994989]],
#          'pointY': [[482.023147,  339.8519037],[639.7415945, 287.6776584],[666.5737787, 463.7638673],[589.4312491, 564.3845581]],
#          'pointZ': [[401.0057988, 120.5122537],[393.220902,  27.74087388],[534.0898691, 37.80294296],[544.1519382, 225.6282324]]}


#2020.8.4_FR (Probably the wrong one?)
#      POINT X-X          X-Y         Y-X        Y-Y         Z-X            Z-Y   
#cam1: 508.8817229	267.1377897	667.7777467	462.039809	526.4760327	74.07798101
#cam2: 456.7275814	368.0002233	527.5915553	537.4893153	426.3346535	227.4486637
#cam3: 327.6524372	322.4532157	539.1979602	255.4095576	315.2851605	17.82766254
#cam4: 552.4286331	330.8305919	580.9893248	243.8201125	498.6282604	73.1201644
#                        CAM1                        CAM2                     CAM3                       CAM4
#static = {'pointX': [[508.8817229, 267.1377897],[456.7275814, 368.0002233],[327.6524372, 322.4532157],[552.4286331, 330.8305919]],
#          'pointY': [[667.7777467, 462.039809],[527.5915553, 537.4893153],[539.1979602, 255.4095576],[580.9893248, 243.8201125]],
#          'pointZ': [[526.4760327, 74.07798101],[426.3346535, 227.4486637],[315.2851605, 17.82766254],[498.6282604, 73.1201644]]}


#2020.8.4_FR (Tried again on 2020.8.27)
#      POINT X-X          X-Y         Y-X        Y-Y         Z-X            Z-Y   
#cam1: 506.8593188	270.2912942	669.1746763	466.3880205	524.985907	72.54669628
#cam2: 457.7902575	372.1308333	529.4067594	539.0107952	423.3332613	218.0877915
#cam3: 325.8480212	321.1805457	540.3628991	257.5143654	316.0973449	21.77742742
#cam4: 551.487489	332.0127751	572.5356485	268.9786674	496.6288531	73.691097
#                        CAM1                        CAM2                     CAM3                       CAM4
#static = {'pointX': [[506.8593188, 270.2912942],[457.7902575, 372.1308333],[325.8480212, 321.1805457],[551.487489,  332.0127751]],
#          'pointY': [[669.1746763, 466.3880205],[529.4067594, 539.0107952],[540.3628991, 257.5143654],[572.5356485, 268.9786674]],
#          'pointZ': [[524.985907,  72.54669628],[423.3332613, 218.0877915],[316.0973449, 21.77742742],[496.6288531, 73.691097]]}


#2020-08-07-RT
#	X	X	Y	Y	Z	Z
#	x	y	x	y	x	y
#cam1	557.0009009	348.1218377	585.5283056	255.8629969	496.9112612	70.13138334
#cam2	289.3930809	331.2071218	524.1873264	260.4859635	276.3803878	11.54748636
#cam3	479.8047276	220.4174975	650.8757764	436.3106552	498.4806063	15.729867
#cam4	385.0349836	181.3921891	526.2136407	518.348984	369.435132	99.49296814

#static = {'pointX': [[557.0009009, 348.1218377],[289.3930809, 331.2071218],[479.8047276, 220.4174975],[385.0349836,  181.3921891]],
#          'pointY': [[585.5283056, 255.8629969],[524.1873264, 260.4859635],[650.8757764, 436.3106552],[526.2136407, 518.348984]],
#          'pointZ': [[496.9112612,  70.13138334],[276.3803878, 11.54748636],[498.4806063, 15.729867],[369.435132, 99.49296814]]}


#2020-08-07-RT3D
#	arm1	arm1	arm2	arm2	shoulder1	shoulder1
#	x	y	x	y	x	y
#cam1	556.1550072	318.5502591	590.5872853	246.4064383	504.2333179	69.87269495
#cam2	349.213698	308.9553891	561.6655504	248.6295545	340.2959659	18.8668104
#cam3	520.4926231	284.3971305	680.2702582	470.2770525	536.0246524	86.21005152
#cam4	463.3357761	383.7667022	540.0112439	548.5373884	430.7079174	241.835517

#static = {'pointX': [[556.1550072, 318.5502591],[349.213698, 308.9553891],[520.4926231, 284.3971305],[463.3357761,  383.7667022]],
#          'pointY': [[590.5872853, 246.4064383],[561.6655504, 248.6295545],[680.2702582, 470.2770525],[540.0112439, 548.5373884]],
#          'pointZ': [[504.2333179,  69.87269495],[340.2959659, 18.8668104],[536.0246524, 86.21005152],[430.7079174, 241.835517]]}

#2020-09-22-RT2D
#	X	X	Y	Y	Z	Z
#	x	y	x	y	x	y
#cam1	517.1872033	263.5312303	588.2434289	153.8049501	476.5836459	82.26534891
#cam2	452.2257636	202.5520421	563.3692135	393.9248696	426.4640368	103.9214309
#cam3	435.7077227	371.1354913	613.3987771	322.7378273	426.0676736	240.6734928
#cam4	567.2712472	337.4798826	708.4608901	480.0537377	571.4238838	251.6587271

#static = {'pointX': [[517.1872033, 263.5312303],[452.2257636, 202.5520421],[435.7077227, 371.1354913],[567.2712472, 337.4798826]],
#          'pointY': [[588.2434289, 153.8049501],[563.3692135, 393.9248696],[613.3987771, 322.7378273],[708.4608901, 480.0537377]],
#          'pointZ': [[476.5836459,  82.26534891],[426.4640368, 103.9214309],[426.0676736, 240.6734928],[571.4238838, 251.6587271]]}


#2020-09-22-RT3D
#	arm1	arm1	arm2	arm2	shoulder1	shoulder1
#	x	y	x	y	x	y
#9	546.9585833	264.9640822	597.4650171	162.3512794	509.6068259	83.47848455
#10	462.7466311	196.8663011	570.7345505	408.9389621	439.3275642	99.93738544
#11	455.1261485	370.0948878	624.0439784	320.7136674	446.2715849	236.9358728
#12	574.3066612	326.7770386	705.9304092	474.7806308	577.2316333	238.4428788

#2020-07-22-RT2D Rocket-Chris
#	X	    X	    Y	    Y	    Z	        Z
#	x	    y	    x	    y	    x	        y
#1 406.4153767	216.5187243	668.1026678	442.7691948	424.5881053	69.31962309
#2 507.2779366	346.8194973	597.4961629	566.0413557	499.6894876	243.110695
#3 379.9852334	446.6277599	443.0163239	208.5927001	324.3695652	253.8267769
#4 176.8703511	373.8421794	523.2740916	164.8417463	133.1750465	161.6830495

#static = {'pointX': [[406.4153767, 216.5187243],[507.2779366, 346.8194973],[379.9852334, 446.6277599],[176.8703511, 373.8421794]],
#          'pointY': [[668.1026678, 442.7691948],[597.4961629, 566.0413557],[443.0163239, 208.5927001],[523.2740916, 164.8417463]],
#          'pointZ': [[424.5881053, 69.31962309],[499.6894876, 243.110695],[324.3695652, 253.8267769],[133.1750465, 161.6830495]]}



#2020-11-25 Test
#	X	X	Y	Y	Z	Z
#	x	y	x	y	x	y
#5in						
#cam1	443.1085196	467.9040349	495.6460439	422.4729002	437.9277924	180.9704713
#cam2	220.9928146	405.0866184	384.0414749	330.2677526	208.1804985	104.4070038
#static = {'pointX': [[443.1085196, 467.9040349],[220.9928146, 405.0866184]],
#          'pointY': [[495.6460439, 422.4729002],[384.0414749, 330.2677526]],
#          'pointZ': [[437.9277924, 180.9704713],[208.1804985, 104.4070038]]}

#2020-11-25 Test
#	X	X	Y	Y	Z	Z
#	x	y	x	y	x	y
#6in						
#cam1	542.4749267	475.5029479	599.2452249	416.1128572	544.2887054	186.179833
#cam2	341.2335451	363.1431765	499.5970366	311.6208215	329.8525244	90.88841407

#static = {'pointX': [[542.4749267, 475.5029479],[341.2335451, 363.1431765]],
#         'pointY': [[599.2452249, 416.1128572],[499.5970366, 311.6208215]],
#         'pointZ': [[544.2887054, 186.179833],[329.8525244, 90.88841407]]}


#2020-11-25 Test
#	X	X	Y	Y	Z	Z
#	x	y	x	y	x	y
#7in						
#cam1	525.7048116	477.1799594	582.723203	421.8385795	524.0278001	175.317887
#cam2	292.6002113	368.174211	483.7795238	311.1558196	282.5381422	83.08225381

#static = {'pointX': [[525.7048116, 477.1799594],[292.6002113, 368.174211]],
#          'pointY': [[582.723203, 421.8385795],[483.7795238, 311.1558196]],
#          'pointZ': [[524.0278001, 175.317887],[282.5381422, 83.08225381]]}


#2020-11-25 Test
#	X	X	Y	Y	Z	Z
#	x	y	x	y	x	y
#9.5in						
#cam1	485.4565353	460.4098443	540.7979152	416.8075449	482.1025122	175.317887
#cam2	369.7427409	361.466165	589.4312491	314.5098426	358.0036603	103.206392

#static = {'pointX': [[485.4565353,460.4098443],[369.7427409,361.466165]],
#          'pointY': [[540.7979152,416.8075449],[589.4312491,314.5098426]],
#          'pointZ': [[482.1025122,175.317887],[358.0036603,103.206392]]}

#2020-11-25 Test
#	X	X	Y	Y	Z	Z
#	x	y	x	y	x	y
#15in						
#cam1	497.1956159	477.1799594	602.8473412	420.1615679	493.8415928	182.0259331
#cam2	286.2408452	338.6374849	606.4223067	301.0989688	274.6480682	95.18916682

static = {'pointX': [[497.1956159,477.1799594],[286.2408452,338.6374849]],
          'pointY': [[602.8473412,420.1615679],[606.4223067,301.0989688]],
          'pointZ': [[493.8415928,182.0259331],[274.6480682,95.18916682]]}










### EXAMPLE
#                        CAM1                        CAM2                     CAM3                       CAM4
#static = {'pointX': [[450.4066209, 426.7973505],[438.5002129, 383.2673147],[507.2576849, 280.9696124],[554.2140073, 410.0994989]],
#          'pointY': [[482.023147,  339.8519037],[639.7415945, 287.6776584],[666.5737787, 463.7638673],[589.4312491, 564.3845581]],
#         'pointZ': [[401.0057988, 120.5122537],[393.220902,  27.74087388],[534.0898691, 37.80294296],[544.1519382, 225.6282324]]}

add_static_points(config, labels, static, snapshots)
#%%
"""
MAKE A COPY OF THE CONFIG FILE, AND CHANGE THE 2D DLC TRACKING PATHS TO
THE FILES WITH STATIC POINTS
"""
#%%
#config = load_config('config_20200804_FR_static.toml' )
config = load_config('config_Test_20201123_15in_static.toml' )

#%%
# You can set 3D reconstruction output path through argument,
# reconstruct_3d(config, output_path='/media/minyoungpark/Min/1101')
# or skip it and set it in config file (reconstruction_output_folder_path).
recovery = reconstruct_3d(config)

# =============================================================================
# #%% Save 3d recovery json file
# import numpy as np
# from json import JSONEncoder
# import json
# class NumpyArrayEncoder(JSONEncoder):
#     def default(self, obj):
#         if isinstance(obj, np.ndarray):
#             return obj.tolist()
#         return JSONEncoder.default(self, obj)
# with open("pop_0610_anipose.json", "w") as write_file:
#     json.dump(recovery, write_file, cls=NumpyArrayEncoder)
#     
# #%% Load 3d recovery json file
# import numpy as np
# from json import JSONEncoder
# import json
# with open("pop_0317_3.json", "r") as read_file:
#     print("Converting JSON encoded data into Numpy array")
#     recovery = json.load(read_file)
# recovery['registration_mat'] = np.array(recovery['registration_mat'])
# recovery['center'] = np.array(recovery['center'])
# =============================================================================


#%% Save 3d recovery json file
import numpy as np
from json import JSONEncoder
import json
class NumpyArrayEncoder(JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return JSONEncoder.default(self, obj)
with open("pop_0610_anipose.json", "w") as write_file:
    json.dump(recovery, write_file, cls=NumpyArrayEncoder)
#%% Load 3d recovery json file
import numpy as np
from json import JSONEncoder
import json
with open("pop_0317_3.json", "r") as read_file:
    print("Converting JSON encoded data into Numpy array")
    recovery = json.load(read_file)
recovery['registration_mat'] = np.array(recovery['registration_mat'])
recovery['center'] = np.array(recovery['center'])

#%% generate 3D picture to find the optimal azimuth and elevation value
generate_three_dim_pictures(config)

#%% Try 3D stick figure video

generate_three_dim_video(config)

#%% Testing if the generated 3D results makes sense, by calculating the distance between wrist and hand
"""
import pandas as pd
import numpy as np
from numpy import array as arr

data_path = 'C:/Users/dongq/DeepLabCut/Han-Qiwei-2020-02-21/3D-data/output_3d_data_rotate4.csv'

df = pd.read_csv(data_path)
wrist = np.empty((len(df), 3))
hand = np.empty((len(df), 3))

wrist[:,0] = arr(df['wrist1_x'])
wrist[:,1] = arr(df['wrist1_y'])
wrist[:,2] = arr(df['wrist1_z'])


hand[:,0] = arr(df['hand2_x'])
hand[:,1] = arr(df['hand2_y'])
hand[:,2] = arr(df['hand2_z'])

dist = wrist - hand
dist_finite = dist[np.isfinite(dist[:,0]), :]
dist_3d = np.linalg.norm(dist_finite, axis=1)
print(np.median(dist_3d))
"""

#%% TEMP test, Matrix product the world coordinates

