# -*- coding: utf-8 -*-
"""
Created on Tue May 14 15:02:39 2024

@author: Paul

Script to plot Figure 2 in EWS SPG paper showing ensemble mean changes in key
variables for the SPG
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from scipy.io import loadmat
import numpy.ma as ma
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.feature as cfeature



def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap

plt.rcParams["font.family"] = "serif"



# Define colour intervals, and colour maps of 30 year mean plots
vr_tas = np.concatenate((np.arange(-2, 0, 0.2),np.arange(0.2, 2.01, 0.2)))
vr_sos = np.concatenate((np.arange(-1, 0, 0.1),np.arange(0.1, 1.01, 0.1)))
vr_mld = np.concatenate((np.arange(-200, 0, 20),np.arange(20, 201, 20)))
vr_obs = np.concatenate((np.arange(-5, 0, 0.5),np.arange(0.5, 5.01, 0.5)))


cmap_name_tas = 'RdBu_r'
cmap_name_sos = 'RdBu_r'
cmap_name_mld = 'BrBG'
cmap_name_obs = 'RdBu'

cmap_tas = plt.cm.get_cmap(cmap_name_tas,len(vr_tas)+1)
cmap_sos = plt.cm.get_cmap(cmap_name_sos,len(vr_sos)+1)
cmap_mld = plt.cm.get_cmap(cmap_name_mld,len(vr_mld)+1)
cmap_obs = plt.cm.get_cmap(cmap_name_obs,len(vr_obs)+1)

cmap_tas = truncate_colormap(cmap_tas, minval=0.1, maxval=0.9)
cmap_sos = truncate_colormap(cmap_sos, minval=0.1, maxval=0.9)
cmap_mld = truncate_colormap(cmap_mld, minval=0.1, maxval=0.9)
cmap_obs = truncate_colormap(cmap_obs, minval=0.1, maxval=0.9)

cmaps = [cmap_tas,cmap_sos,cmap_mld,cmap_obs]

norm_tas = colors.BoundaryNorm(vr_tas, cmap_tas.N)
norm_sos = colors.BoundaryNorm(vr_sos, cmap_sos.N)
norm_mld = colors.BoundaryNorm(vr_mld, cmap_mld.N)
norm_obs = colors.BoundaryNorm(vr_obs, cmap_obs.N)

norms = [norm_tas,norm_sos,norm_mld,norm_obs]

# Define projection used for maps
proj = 'merc'

########### CMIP6 models ####################


################# Specify variable of interest ######################
var = ['tas','sos','mlotst','msftbarot']

############## Specify time points of interest ###########################
years = 20
start_year = 1850
crit_year = [2040,2035,2040]

################ Specify region of interest #######################
region = 'Atlantic'

################# Specify model details here ################################## 
models = ['CESM2-WACCM','NorESM2-LM','MRI-ESM2-0'] 
variant_ids = ['r1i1p1f1']
date_ends = ['2100']

# Specify experiments
experiment_hist = 'historical'
experiment_ssp = ['ssp126','ssp126','ssp245']

# The projection keyword determines how the plot will look
fig, ax = plt.subplots(2, 2, subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(13, 5.1))

lon_low = -80
lon_high = 0
lat_low = 30
lat_high = 65
project = ccrs.PlateCarree()

ax_pos_tas = [0.05,0.54, 0.37, 0.45]
ax_pos_sos = [0.53,0.54, 0.37, 0.45]
ax_pos_mld = [0.05,0.05, 0.37, 0.45]
ax_pos_obs = [0.53,0.05, 0.37, 0.45]
ax_pos = [ax_pos_tas,ax_pos_sos,ax_pos_mld,ax_pos_obs]

cbar_pos_tas = [0.425, 0.558, 0.015, 0.415]
cbar_pos_sos = [0.905, 0.558, 0.015, 0.415]
cbar_pos_mld = [0.425, 0.067, 0.015, 0.415]
cbar_pos_obs = [0.905, 0.067, 0.015, 0.415]
cbar_pos = [cbar_pos_tas,cbar_pos_sos,cbar_pos_mld,cbar_pos_obs]

cbar_lab_tas = 'Temperature change ($^o$C)'
cbar_lab_sos = 'Sea surface salinity\nchange (psu)'
cbar_lab_mld = 'Mixed layer depth (m)'
cbar_lab_obs = 'Ocean barotropic\nstreamfunction change (Sv)'
cbar_lab = [cbar_lab_tas,cbar_lab_sos,cbar_lab_mld,cbar_lab_obs]

cbar_ticks_tas = np.arange(-2,2.01,step=1)
cbar_ticks_sos = np.arange(-1,1.01,step=0.5)
cbar_ticks_mld = np.arange(-200,200.1,step=100)
cbar_ticks_obs = np.arange(-5,5.01,step=2.5)
cbar_ticks = [cbar_ticks_tas,cbar_ticks_sos,cbar_ticks_mld,cbar_ticks_obs]


for j in range(len(var)):

    # Initialise empty arrays
    x_ts_start = []
    x_ts_end = []

        
    for i in range(len(models)):
    
        # File name of the data
        fname = var[j]+'_'+models[i]+'_'+experiment_hist+'_'+experiment_ssp[i]+'_'+variant_ids[0]+'_1850-'+date_ends[0]+'_'+region+'_uniform_grid.mat'#'_monthly.nc'

        LONS = loadmat('data/'+fname)['LONS'][0]
        LATS = loadmat('data/'+fname)['LATS'][0]
        
        ts_start = loadmat('data/'+fname)['region_data'][crit_year[i]-years-start_year:crit_year[i]-start_year,:]
        ts_end = loadmat('data/'+fname)['region_data'][crit_year[i]-start_year:crit_year[i]+years-start_year,:]    
        
        x_ts_start.append(ma.average(ts_start,axis=0))
        x_ts_end.append(ma.average(ts_end,axis=0))
        
    
    # Store T data
    x_tas_start = np.stack(x_ts_start)
    x_tas_end = np.stack(x_ts_end)
    
    x_tas_change = x_tas_end - x_tas_start
    
    res = np.zeros(len(LONS))
    for k in range(len(LONS)):
        res[k] = np.all(x_tas_change[:,k]>0) | np.all(x_tas_change[:,k]<0)


    x_tas_mean_change = ma.average(x_tas_change,axis=0)

    
    LONS = LONS.reshape(36,81)
    LATS = LATS.reshape(36,81)
    DATA = (x_tas_mean_change).reshape(36,81)
    RES = (res).reshape(36,81)


    im = ax[int(j/2),j%2].pcolor(LONS,LATS,DATA,cmap=cmaps[j],norm=norms[j],transform=project)
    im2 = ax[int(j/2),j%2].scatter(LONS[RES==0],LATS[RES==0],c='k',transform=project,s=0.2)
    
    if j>0:
        ax[int(j/2),j%2].add_feature(cfeature.LAND, color='silver', zorder=2)
        ax[int(j/2),j%2].add_feature(cfeature.COASTLINE, zorder=3)
    


    
    ax[int(j/2),j%2].coastlines(resolution='50m',zorder=5)
    
    

    ax[int(j/2),j%2].set_yticks(np.arange(30,90,step=5), crs=project)
    ax[int(j/2),j%2].yaxis.set_major_formatter(LATITUDE_FORMATTER)
    ax[int(j/2),j%2].set_xticks(np.arange(-80,10,step=10), crs=project)
    ax[int(j/2),j%2].xaxis.set_major_formatter(LONGITUDE_FORMATTER)
    
    if j%2 == 1:
        ax[int(j/2),j%2].set_yticklabels([])
    if int(j/2) == 0:
        ax[int(j/2),j%2].set_xticklabels([])

    
    
    ax[int(j/2),j%2].grid(color='black', linestyle='dotted')
    
    ax[int(j/2),j%2].set_xlim(lon_low,lon_high)
    ax[int(j/2),j%2].set_ylim(lat_low,lat_high)
    
    ax[int(j/2),j%2].plot([-50,-10],[45,45],'k',transform=ccrs.PlateCarree())
    ax[int(j/2),j%2].plot([-50,-10],[60,60],'k',transform=ccrs.PlateCarree())
    ax[int(j/2),j%2].plot([-50,-50],[45,60],'k',transform=ccrs.PlateCarree())
    ax[int(j/2),j%2].plot([-10,-10],[45,60],'k',transform=ccrs.PlateCarree())
    


    ax[int(j/2),j%2].set_position(ax_pos[j])
    cbar_ax = fig.add_axes(cbar_pos[j])
    cbar = plt.colorbar(im, cax=cbar_ax,spacing="proportional")
    cbar.set_label(cbar_lab[j])
    cbar.set_ticks(cbar_ticks[j])