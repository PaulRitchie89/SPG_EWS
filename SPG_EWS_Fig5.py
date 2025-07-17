# -*- coding: utf-8 -*-
"""
Created on Fri Feb  2 13:30:25 2024

@author: Paul

Script to plot one column (change model name for other columns) of Figure 5 in
EWS SPG paper showing colour plots of the detection index (max/min of detection
time series), year of abrupt shift and the clusters based on these two statistics 
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors

from scipy.io import loadmat

from warnings import filterwarnings
filterwarnings(action='ignore', category=DeprecationWarning, message='`np.bool` is a deprecated alias')

from numpy import genfromtxt


def create_polar_map(fig,ax,LONS,LATS,stat,cmap,norm,project,map_project,model):
    ax.set_title(model)
    ax.gridlines()
    ax.set_extent([-180, 180, 60, 90], project)
    ax.add_feature(cfeature.LAND, color='silver', zorder=2)
    ax.add_feature(cfeature.COASTLINE, zorder=3)
    theta  = np.linspace(0, 2*np.pi, 100)
    center = [0.5, 0.5]
    radius =  0.5
    verts  = np.vstack([np.sin(theta), np.cos(theta)]).T
    circle = mpath.Path(verts * radius + center)
    ax.set_boundary(circle, transform=ax.transAxes)
    im = ax.scatter(LONS,LATS,c=stat,cmap=cmap,norm=norm,transform=project,s=5,zorder=1)
    
    return im

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap


v2  = [-1, -0.9, -0.7, -0.5, -0.2, 0.2, 0.5, 0.7, 0.9, 1]
cmap2_name = 'RdBu'
cmap2 = plt.cm.get_cmap(cmap2_name,len(v2)+1)
colors2 = list(cmap2(np.arange(len(v2)+1)))
norm2 = colors.BoundaryNorm(v2, cmap2.N)


v3  = [1850, 1900, 1950, 2000, 2020, 2040, 2060, 2080, 2100]
cmap3_name = 'hot'
cmap3 = plt.cm.get_cmap(cmap3_name,100)
cmap3 = truncate_colormap(cmap3, minval=0.0, maxval=0.9, n=100)
colors3 = list(cmap3(np.arange(len(v3)+1)))
norm3 = colors.BoundaryNorm(v3, cmap3.N)



import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.feature as cfeature
import matplotlib.path as mpath

project = ccrs.PlateCarree()
map_project = ccrs.NorthPolarStereo()

########### CMIP6 models ####################


################# Specify variable of interest ######################
var = 'msftbarot'

################ Specify region of interest #######################
region = 'Atlantic'

################# Specify model details here ################################## 
models = ['NorESM2-LM']

# Specify experiments
experiment_ssp = 'ssp126'

# File name of the data
filename = var+'_'+models[0]+'_'+experiment_ssp+'_clusters_scaled.csv'

my_data = genfromtxt('data/'+filename, delimiter=',')

LONS2 = my_data[1:,1]
LATS2 = my_data[1:,2]
CLUSTER = my_data[1:,7]

v4  = np.linspace(0.5,np.max(CLUSTER)+0.5,int(np.max(CLUSTER)+1))
cmap4_name = 'tab10'
cmap4 = plt.cm.get_cmap(cmap4_name)
cmap4 = truncate_colormap(cmap4, minval=0.0, maxval=np.max(CLUSTER)/10-0.01, n=100)
colors4 = list(cmap4(np.arange(len(v4)+1)))
norm4 = colors.BoundaryNorm(v4, cmap4.N)
    
fname = var+'_'+models[0]+'_'+experiment_ssp+'_detection_data_'+region+'.mat'
         
LONS = loadmat('data/'+fname)['LONS'][0]
LATS = loadmat('data/'+fname)['LATS'][0]
detection_index = loadmat('data/'+fname)['detection_index']
tip_time = loadmat('data/'+fname)['tip_time']

fig, ax = plt.subplots(3,1, figsize=(6.5, 7.25),subplot_kw={'projection':project})

    
for j in range(3):
    ax[j].add_feature(cfeature.LAND, color='silver', zorder=2)
    ax[j].add_feature(cfeature.COASTLINE, zorder=3)
    
    ax[j].set_yticks(np.arange(30,90,step=5), crs=project)
    ax[j].set_xticks(np.arange(-80,10,step=10), crs=project)
    ax[j].xaxis.set_major_formatter(LONGITUDE_FORMATTER)
    ax[j].yaxis.set_major_formatter(LATITUDE_FORMATTER)
    ax[j].grid(color='black', linestyle='dotted')
    
    ax[j].set_xlim(-80,0)
    ax[j].set_ylim(30,65)

im = ax[0].scatter(LONS,LATS,c=detection_index[0,:],cmap=cmap2,norm=norm2,transform=project,s=5,zorder=1)
im2 = ax[1].scatter(LONS,LATS,c=tip_time[0,:],cmap=cmap3,norm=norm3,transform=project,s=5,zorder=1)
im3 = ax[2].scatter(LONS2,LATS2,c=CLUSTER,cmap=cmap4,norm=norm4,transform=project,s=5,zorder=1)


fig.subplots_adjust(hspace=0.2, wspace=0, bottom=0.06, top=0.94, left=0.02, right=0.9)

cbar_ax = fig.add_axes([0.85, 0.675, 0.04, 0.265])
cbar = plt.colorbar(im, cax=cbar_ax,spacing="proportional")
cbar.set_label('Detection index')
cbar.set_ticks(np.linspace(-1,1,5))

cbar_ax2 = fig.add_axes([0.85, 0.365, 0.04, 0.265])
cbar2 = plt.colorbar(im2, cax=cbar_ax2,spacing="proportional")
cbar2.set_label('Year of abrupt shift')
cbar2.set_ticks(np.linspace(1850,2100,6))

cbar_ax3 = fig.add_axes([0.85, 0.055, 0.04, 0.265])
cbar3 = plt.colorbar(im3, cax=cbar_ax3,spacing="proportional")
cbar3.set_label('Cluster')
cbar3.set_ticks(np.linspace(1,int(np.max(CLUSTER)),int(np.max(CLUSTER))))

fig.suptitle(models[0]+', '+experiment_ssp)