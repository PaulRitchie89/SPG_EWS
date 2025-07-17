# -*- coding: utf-8 -*-
"""
Created on Fri Feb  9 11:54:31 2024

@author: Paul

Script to plot one column (change model to plot other column) of Figure 6 in EWS
SPG paper showing each cluster's mean time series for barotropic streamfunction
along with detection time series and EWS (AR(1) and variance) 
Require regimeshifts Python library https://github.com/BeatrizArellano/regimeshifts
Last accessed 17th July 2025
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat

from regimeshifts import regime_shifts as rs
from regimeshifts import ews

from warnings import filterwarnings
filterwarnings(action='ignore', category=DeprecationWarning, message='`np.bool` is a deprecated alias')

from scipy.stats import norm

from numpy import genfromtxt
import numpy.ma as ma

from matplotlib import rc

fontsize = 14
rc('font', **{'size' : fontsize})
rc('text', usetex=True)

def sig_test(ts,sig_ts,test='positive',**kwargs):
    ts = ts.loc[ts.first_valid_index():ts.last_valid_index()]
    high_tail = sig_ts.quantile(0.975)
    low_tail = sig_ts.quantile(0.025)    
    time_ind = np.nan
    std_err = norm.ppf(0.95)*np.sqrt(sig_ts.loc[:,'Time series'].var(ddof=0))
    if test == 'positive':
        for j in range(len(ts)):
            if all(i >= float(std_err) for i in ts.loc[ts.first_valid_index()+j:,'Time series']-ts.loc[ts.first_valid_index(),'Time series']):
                time_ind = len(ts) - j
                break
    elif test == 'negative':
        for j in range(len(ts)):
            if all(i <= float(low_tail) for i in ts.loc[j:,'Time series']):
                time_ind = j
                break
    return (low_tail,high_tail,std_err,time_ind)



########### CMIP6 models ####################


################# Specify variable of interest ######################
var = 'msftbarot'

############## Specify time points of interest ###########################
years_start = 1
years = 1


################ Specify region of interest #######################
region = 'Atlantic'

################# Specify model details here ################################## 
model = 'CESM2-WACCM' 
variant_id = 'r1i1p1f1'
date_ends = '2100'

# Specify experiments
experiment_hist = 'historical'
experiment_ssp = 'ssp126'

# Initialise empty arrays
x_ts = []
x_ts_control = []

# File name of the data
fname = var+'_'+model+'_'+experiment_hist+'_'+experiment_ssp+'_'+variant_id+'_1850-'+date_ends+'_'+region+'.mat'
   
LONS = loadmat('data/'+fname)['LONS'][0]
LATS = loadmat('data/'+fname)['LATS'][0]
ts = loadmat('data/'+fname)['region_data']

date_ends2 = '0499'
fname2 = var+'_'+model+'_piControl_'+variant_id+'_0001-'+date_ends2+'_'+region+'.mat'

ts_control = loadmat('data/'+fname2)['region_data']

filename = var+'_'+model+'_'+experiment_ssp+'_clusters_scaled.csv'

my_data = genfromtxt('data/'+filename, delimiter=',')

LONS2 = my_data[1:,1]
LATS2 = my_data[1:,2]
CLUSTER = my_data[1:,7]#7]

nc = int(np.max(CLUSTER))

t = np.linspace(1850,int(date_ends),int(date_ends)-1850+1)


detection_index2 = np.zeros(nc)
tip_time2 = np.zeros(nc)

ar1_tau = np.zeros(nc)
var_tau = np.zeros(nc)

ar1_time = np.zeros(nc)
var_time = np.zeros(nc)

wL = 100 ## Window length specified in number of points in the series
margin = 40
bW = 60

lw = 1.5

n2 = np.min([500,int(date_ends2)])

colours = ['tab:blue','tab:orange','tab:green','tab:red','tab:purple','tab:brown']

fig,ax = plt.subplots(4,1,sharex=True,figsize=(6.4,8.5))
ax[0].set_ylabel('Normalised '+var)
ax[1].set_ylabel('Detection metric')
ax[2].set_ylabel('Change in AR(1)')
ax[3].set_ylabel('Variance')
ax[3].set_xlabel('Year')
ax[0].set_title(model+', '+experiment_ssp, fontsize=fontsize)
fig.tight_layout()
ax[0].set_xlim(1850,2100)
ax[1].set_ylim(-0.75,0.75)


for j in range(nc):
    LONS3 = LONS2[CLUSTER==j+1]
    LATS3 = LATS2[CLUSTER==j+1]
    
    x_ts = np.zeros((int(date_ends)-1850+1,len(LONS3)))
    x_ts_control = np.zeros((n2,len(LONS3)))
    
    for k in range(len(LONS3)):
        idx = np.nanargmin(np.abs(LONS-LONS3[k])+np.abs(LATS-LATS3[k]))
        x_ts[:,k] = ts[:,idx]
        x_ts_control[:,k] = ts_control[:n2,idx]
    
    weights = np.cos(np.deg2rad(LATS3))
    

    x_TS = ma.average(x_ts,axis=1,weights=weights)
    

    x_TS_CONTROL = ma.average(x_ts_control,axis=1,weights=weights)


    
    ts2 = rs.Regime_shift(x_TS)
    detection = ts2.as_detect()
    detection_index2[j] = detection[np.argmax(np.abs(detection))]
    
    if detection_index2[j] < 0:
        bef_rs = ts2.before_drs()
    else:
        bef_rs = ts2.before_urs()
    
    tip_time2[j] = t[len(bef_rs)]
    

    ax[0].plot(t,(x_TS-np.min(x_TS))/(np.max(x_TS)-np.min(x_TS)))
    ax[1].plot(t,detection)
    
    if (len(bef_rs) > wL + margin):
        
        
        series = ews.Ews(bef_rs)
        series = series.rename(columns={0:'Time series'}) ## The Ews class returns an extended Dataframe object, if we provided a series, it sets 0 for the column name. 
                        
        ar1 = series.ar1(detrend=True,bW=bW,wL=wL) ### Computing lag-1 autocorrelation using the ar1() method
        var = series.var(detrend=True,bW=bW,wL=wL) ## Computing variance
        
        a = series.gaussian_det(bW)
        
        ar1_tau[j] = ar1.kendall
        var_tau[j] = var.kendall
        

        ts3 = rs.Regime_shift(x_TS_CONTROL[:n2])
        series2 = ews.Ews(ts3.iloc[:499])
        series2 = series2.rename(columns={0:'Time series'})
        
        
        sig_ar1 = series2.ar1(detrend=True,bW=bW,wL=wL) ### Computing lag-1 autocorrelation using the ar1() method
        sig_var = series2.var(detrend=True,bW=bW,wL=wL) ## Computing variance
        
        ar1_low, ar1_high, ar1_err, ar1_time[j] = sig_test(ar1,sig_ar1,test='positive')
        var_low, var_high, var_err, var_time[j] = sig_test(var,sig_var,test='positive')
        n=len(bef_rs)
        ax[2].plot(t[:n],ar1-ar1.iloc[int(wL-1)],c=colours[j],lw=lw,label='$\\tau$ = '+"{:.2f}".format(ar1_tau[j]))
        ax[3].plot(t[:n],var,c=colours[j],lw=lw,label='$\\tau$ = '+"{:.2f}".format(var_tau[j]))
        
        ax[2].fill_between([1850,2100],-ar1_err,ar1_err,alpha=0.3,facecolor=colours[j],edgecolor='none')
        ax[3].fill_between([1850,2100],var.iloc[int(wL-1)]-var_err,var.iloc[int(wL-1)]+var_err,alpha=0.3,facecolor=colours[j],edgecolor='none')
        
ax[2].legend(frameon=False)
ax[3].legend(frameon=False)

fig.subplots_adjust(left=0.175, right=0.95, bottom=0.075, top=0.95, hspace=0.1)

import seaborn as sns
sns.despine()

plt.text(-0.21,0.95,'$\\textbf{(i)}$',transform=ax[0].transAxes)
plt.text(-0.21,0.95,'$\\textbf{(j)}$',transform=ax[1].transAxes)
plt.text(-0.21,0.95,'$\\textbf{(k)}$',transform=ax[2].transAxes)
plt.text(-0.21,0.95,'$\\textbf{(l)}$',transform=ax[3].transAxes)
