
import os
import numpy as np
import gdal
import matplotlib.pyplot as plt
import geopandas as gpd
import pandas as pd
from shapely.geometry import Polygon,Point

# Import custom modules
from CustomFunctions import gis as gis
from CustomFunctions import basic_functions as bf

# Set figure properties
params={'font.sans-serif':'Arial',
        'font.size':7,
        'axes.labelsize':7,
        'axes.titlesize':14,
        'axes.linewidth':0.5,        
        'xtick.labelsize':7,
        'xtick.major.width':0.5,
        'xtick.major.size':5,
        'xtick.direction':'in',
        'ytick.labelsize':7,
        'ytick.major.width':0.5,
        'ytick.major.size':5,
        'ytick.direction':'in',
        'legend.fontsize':10,
        'savefig.dpi':150}
plt.rcParams.update(params)

'''============================================================================
IMPORT DATA
============================================================================'''

# BC land basemap
gdf_bm=gpd.read_file(r'Z:\!Workgrp\Forest Carbon\Data\Basemaps\bcbound.shp')

# Import TSA maps
zTSA=gis.OpenGdal(r'Z:\!Workgrp\Forest Carbon\Data\BC1ha\Admin\tsa.tif')
tsa_key=pd.read_excel(r'Z:\!Workgrp\Forest Carbon\Data\BC1ha\Admin\tsa.tif.vat.xlsx')
tsa_boundaries=gpd.read_file(r'C:\Users\rhember\Documents\Data\TSA\tsa_boundaries.shp')


# Roads
#gdf_rd=gpd.read_file(r'Z:\!Workgrp\Forest Carbon\Data\Basemaps\roads.shp')

# Districts
#gdf_d=gpd.read_file(r'Z:\!Workgrp\Forest Carbon\Data\Districts\district.shp')



'''============================================================================
ISOLATE ROI
============================================================================'''

# Pick the TSAs to include
#tsaList=['Kamloops TSA','100 Mile House TSA','Williams Lake TSA']
tsaList=list(tsa_key['Name'])
ind_tsa=tsa_key.VALUE[np.isin(tsa_key.Name,tsaList)].values

tsa_mask1=np.zeros((zTSA['m'],zTSA['n']))
ind=(np.isin(zTSA['Data'],ind_tsa))
tsa_mask1[ind]=1

# Define extent based on mask
xlim=[np.min(zTSA['X'][ind])-5000,np.max(zTSA['X'][ind])+5000]
ylim=[np.min(zTSA['Y'][ind])-5000,np.max(zTSA['Y'][ind])+5000]

tsa_boundaries1=tsa_boundaries[np.isin(tsa_boundaries.Name,tsaList)]
tsa_boundaries1=tsa_boundaries1.reset_index(drop=True)

# Isolate lakes inside the ROI
#lakes=gpd.overlay(gdf_bm[(gdf_bm['TAG']=='lake') | (gdf_bm['TAG']=='river')],tsa_boundaries[np.isin(tsa_boundaries.Name,tsaList)],how='intersection')




'''============================================================================
PLOT A TSA
============================================================================'''

# Create and adjust the colormap
cm=plt.cm.get_cmap('viridis',2)
cm.colors[0,:]=[1,1,1,1]

plt.close('all')
fig,ax=plt.subplots(1)
mngr=plt.get_current_fig_manager()
mngr.window.setGeometry(100,100,700,600)
gdf_bm.plot(ax=ax,facecolor=[0.8,0.8,0.8],edgecolor=[0,0,0],label='Districts',linewidth=0.25,alpha=1)
tsa_boundaries1.plot(ax=ax,color="none",edgecolor=[0,0,0],linewidth=0.5,facecolor="none")
ax.set(position=[0.04,0.03,0.92,0.95],xlim=xlim,ylim=ylim)
#plt.imshow(tsa_mask1[0::2,0::2],clim=(0,1),extent=extent,cmap=cm)
#lakes.plot(ax=ax,color=[0.8,0.85,1],label='Water')
#roads.plot(ax=ax,color=[0,0,0],label='Road')





'''============================================================================
PLOT BGC Zones within the TSA mask
============================================================================'''

# BGC classification
zBGC=gis.OpenGdal(r'Z:\!Workgrp\Forest Carbon\Data\BC1ha\VRI\becz.tif')
zBGC['key']=pd.read_excel(r'Z:\!Workgrp\Forest Carbon\Data\BC1ha\VRI\becz.tif.vat.xlsx')

u=np.unique(zBGC['Data'][tsa_mask1==1])

lab=[]
for i in range(len(u)):
    try:
        lab.append(zBGC['key'].ZONE[zBGC['key'].VALUE==u[i]].values[0])
    except:
        lab.append('')

z1=np.ones(zBGC['Data'].shape)
for i in range(u.size):
    z1[(tsa_mask1==1) & (zBGC['Data']==u[i])]=i
L=i+1
z1[(tsa_mask1!=1)]=L

cm=plt.cm.get_cmap('viridis',L)
cm.colors[-1,:]=[1,1,1,1]


plt.close('all')
fig,ax=plt.subplots(1,2)
mngr=plt.get_current_fig_manager()
mngr.window.setGeometry(100,100,750,750)
im=ax[0].imshow(z1[0::2,0::2],clim=(0,L+1),extent=zTSA['Extent'],cmap=cm)
tsa_boundaries1.plot(ax=ax[0],color=None,edgecolor=[0,0,0],facecolor='none')
ax[0].set(position=[0.04,0.02,0.92,0.96],xlim=xlim,ylim=ylim,aspect='auto')

cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,L+1,1),ticks=np.arange(0.5,L+0.5,1))
cb.ax.set(yticklabels=lab)
cb.ax.tick_params(labelsize=11,length=0)
for i in range(0,L):
    ax[1].plot([0,100],[i/L,i/L],'k-',linewidth=0.5)
ax[1].set(position=[0.06,0.04,0.025,0.4])












