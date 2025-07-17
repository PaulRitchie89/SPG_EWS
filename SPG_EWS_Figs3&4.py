# -*- coding: utf-8 -*-
"""
Created on Tue Feb 13 11:54:25 2024

@author: Paul

Script to plot Figures 3 & 4 in EWS SPG paper showing box mean time series for
variable of interest along with detection time series and EWS (AR(1) and variance)
Require regimeshifts Python library https://github.com/BeatrizArellano/regimeshifts
Last accessed 17th July 2025
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from scipy.io import loadmat
import seaborn as sns
import numpy.ma as ma

from regimeshifts import regime_shifts as rs
from regimeshifts import ews

from warnings import filterwarnings
filterwarnings(action='ignore', category=DeprecationWarning, message='`np.bool` is a deprecated alias')

from scipy.stats import norm

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


def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap


########### CMIP6 models ####################


################# Specify variable of interest ######################
var = 'msftbarot'

############## Specify time points of interest ###########################
years_start = 1
years = 1


################ Specify region of interest #######################
region = 'Atlantic'

################# Specify model details here ################################## 
models = ['NorESM2-LM','CESM2-WACCM','MRI-ESM2-0'] 
variant_ids = ['r1i1p1f1']
date_ends = ['2100']
control_ends = ['0501','0499','0701']

colours = ['tab:blue','tab:orange','tab:green']

# Specify experiments
experiment_hist = 'historical'
experiment_ssps = ['ssp126','ssp126','ssp245']

t = np.linspace(1850,2100,251)
fig, ax = plt.subplots(4,1,sharex=True,figsize=(6.4,8.5))

# ax[0].set_ylabel('Change in barotropic\nstreamfunction (Sv)',fontsize=fontsize)
ax[0].set_ylabel('Change in air\ntemperature ($^o$C)',fontsize=fontsize)
ax[1].set_ylabel('Detection index',fontsize=fontsize)
ax[2].set_ylabel('AR(1)',fontsize=fontsize)
ax[3].set_ylabel('Variance',fontsize=fontsize)
ax[3].set_xlabel('Year',fontsize=fontsize)

ax[0].set_xlim(t[0],t[-1])


    
for l in range(len(models)):

    
    # File name of the data
    fname = var+'_'+models[l]+'_'+experiment_hist+'_'+experiment_ssps[l]+'_'+variant_ids[0]+'_1850-'+date_ends[0]+'_'+region+'_corrected.mat'#'_monthly.nc'
    fname2 = var+'_'+models[l]+'_piControl_'+variant_ids[0]+'_0001-'+control_ends[l]+'_'+region+'_corrected.mat'
    
    LONS = loadmat('data/'+fname)['LONS'][0]
    LATS = loadmat('data/'+fname)['LATS'][0]
    x_tas = loadmat('data/'+fname)['region_data']
    
    x_tas_control = loadmat('data/'+fname2)['region_data']



    
    weights = np.cos(np.deg2rad(LATS))
    bW = 60
    wL = 100
    margin = 40
    
    lag=1
    
    
    x_tas_mean = np.zeros(251)
    x_tas_control_mean = np.zeros(499)

    for i in range(251):
        x_tas_mean[i] = ma.average(x_tas[i,:].data,weights=weights)#


    for i in range(499):
        x_tas_control_mean[i] = ma.average(x_tas_control[i,:].data,weights=weights)#



    ts = rs.Regime_shift(x_tas_mean)
    detection = ts.as_detect()
    detection_index = detection[np.argmax(np.abs(detection))]
        
    if detection_index < 0:
        bef_rs = ts.before_drs()
    else:
        bef_rs = ts.before_urs()
    
    if l==0:
        zo=5
    else:
        zo=1
    
    if var == 'tas':
        ax[0].plot(t,ts-ts[0],zorder=zo)
    elif var == 'msftbarot':
        ax[0].plot(t,ts-ts[0],zorder=zo)
    else:
        ax[0].plot(t,ts)
    
    ax[1].plot(t,detection,label=models[l]+', '+experiment_ssps[l],zorder=zo)
    

    series = ews.Ews(bef_rs)
    series = series.rename(columns={0:'Time series'}) ## The Ews class returns an extended Dataframe object, if we provided a series, it sets 0 for the column name. 
    
    if (len(bef_rs) > wL + margin):
              
        ar1 = series.ar1(detrend=True,bW=bW,wL=wL) ### Computing lag-1 autocorrelation using the ar1() method
        vari = series.var(detrend=True,bW=bW,wL=wL) ## Computing variance
        
        ar1_tau = ar1.kendall
        var_tau = vari.kendall
        
        ts2 = rs.Regime_shift(x_tas_control_mean)
        series2 = ews.Ews(ts2)
        series2 = series2.rename(columns={0:'Time series'})
        
        sig_ar1 = series2.ar1(detrend=True,bW=bW,wL=wL) ### Computing lag-1 autocorrelation using the ar1() method
        sig_var = series2.var(detrend=True,bW=bW,wL=wL) ## Computing variance
        
        ar1_low, ar1_high, ar1_err, ar1_time = sig_test(ar1,sig_ar1,test='positive')
        var_low, var_high, var_err, var_time = sig_test(vari,sig_var,test='positive')
   
        
        ax[2].plot(t[:len(ar1)],ar1,label='$\\tau$ = '+"{:.2f}".format(ar1_tau),color=colours[l],zorder=zo)
        ax[3].plot(t[:len(vari)],vari,label='$\\tau$ = '+"{:.2f}".format(var_tau),color=colours[l],zorder=zo)
        ax[2].fill_between([1850,2100],ar1.loc[wL-1,'Time series']-ar1_err,ar1.loc[wL-1,'Time series']+ar1_err,alpha=0.3,color=colours[l])#,zorder=zo)
        ax[3].fill_between([1850,2100],vari.loc[wL-1,'Time series']-var_err,vari.loc[wL-1,'Time series']+var_err,alpha=0.3,color=colours[l])#,zorder=zo)



ax[1].legend(frameon=False)
ax[2].legend(frameon=False)
ax[3].legend(frameon=False)


fig.subplots_adjust(left=0.175, right=0.95, bottom=0.075, top=0.975, hspace=0.1)

sns.despine()

plt.text(-0.21,0.95,'$\\textbf{(e)}$',transform=ax[0].transAxes)
plt.text(-0.21,0.95,'$\\textbf{(f)}$',transform=ax[1].transAxes)
plt.text(-0.21,0.95,'$\\textbf{(g)}$',transform=ax[2].transAxes)
plt.text(-0.21,0.95,'$\\textbf{(h)}$',transform=ax[3].transAxes)