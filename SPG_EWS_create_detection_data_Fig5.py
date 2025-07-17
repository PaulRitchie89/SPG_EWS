# -*- coding: utf-8 -*-
"""
Created on Thu Jul 17 09:14:27 2025

@author: Paul

Script to calculate detection data (max/min of detection time series, year of
abrupt shift) that is used for Figure 5 in EWS SPG paper.
Require regimeshifts Python library https://github.com/BeatrizArellano/regimeshifts
Last accessed 17th July 2025
"""

import numpy as np
from scipy.io import loadmat
from regimeshifts import regime_shifts as rs
from scipy.io import savemat
from warnings import filterwarnings
filterwarnings(action='ignore', category=DeprecationWarning, message='`np.bool` is a deprecated alias')




########### CMIP6 models ####################


################# Specify variable of interest ######################
var = 'msftbarot'

################ Specify region of interest #######################
region = 'Atlantic'

################# Specify model details here ################################## 
model = 'MRI-ESM2-0' 
variant_id = 'r1i1p1f1'
date_ends = '2100'

# Specify experiments
experiment_hist = 'historical'
experiment_ssp = 'ssp245'

# Initialise empty arrays
x_ts = []
x_ts_control = []

# File name of the data
fname = var+'_'+model+'_'+experiment_hist+'_'+experiment_ssp+'_'+variant_id+'_1850-'+date_ends+'_'+region+'.mat'
   
LONS = loadmat('data/'+fname)['LONS'][0]
LATS = loadmat('data/'+fname)['LATS'][0]
ts = loadmat('data/'+fname)['region_data']
    

# Obtain dimension sizes
nt = int(len(ts[:,0]))
nx = int(len(ts[0,:]))

t = np.linspace(1850,2100,nt)



detection_index = np.zeros(nx)
tip_time = np.zeros(nx)

for i in range(nx):
    if ts[0,i] == 0:
        tip_time[i], detection_index[i] = np.nan, np.nan
    else:
        ts2 = rs.Regime_shift(ts[:,i].data)
        detection = ts2.as_detect()
        detection_index[i] = detection[np.argmax(np.abs(detection))]
        
        if detection_index[i] < 0:
            bef_rs = ts2.before_drs()
        else:
            bef_rs = ts2.before_urs()
        
        tip_time[i] = t[len(bef_rs)]
        


mdic = {"LONS": LONS, "LATS": LATS, "detection_index": detection_index, "tip_time": tip_time}

outname = 'data/'++var+'_'+model+'_'+experiment_ssp+'_detection_data_'+region+'.mat'

savemat(outname, mdic)