
#%% Import modules

import os
import numpy as np
import gc
import gdal
import matplotlib.pyplot as plt
import matplotlib.colors
import geopandas as gpd
import pandas as pd
from shapely.geometry import Polygon,Point

# Import custom modules
import fcgadgets.utilities.utilities_general as gu
import fcgadgets.utilities.utilities_gis as gis

#%% Set figure properties
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

#%% Import data

# BC land basemap
gdf_bc_boundary=gpd.read_file(r'C:\Users\rhember\Documents\Data\Basemaps\Basemaps.gdb',layer='NRC_POLITICAL_BOUNDARIES_1M_SP')
gdf_bm=gpd.read_file(r'C:\Users\rhember\Documents\Data\Basemaps\bc_land.shp')

# Import TSA maps
zTSA=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Admin\tsa.tif')
tsa_key=pd.read_excel(r'C:\Users\rhember\Documents\Data\BC1ha\Admin\lut_tsa.xlsx')
gdf_tsa=gpd.read_file(r'C:\Users\rhember\Documents\Data\TSA\tsa_boundaries.shp')

# Roads
#gdf_road=gpd.read_file(r'Z:\!Workgrp\Forest Carbon\Data\Basemaps\roads.shp')

# Districts
#gdf_d=gpd.read_file(r'Z:\!Workgrp\Forest Carbon\Data\Districts\district.shp')


#%% Area of interest

# Pick the TSAs to include
#tsaList=['Soo TSA']
#tsaList=['Kamloops TSA','100 Mile House TSA','Williams Lake TSA']
#tsaList=['Williams Lake TSA']
tsaList=list(tsa_key['Name'])

# Index to TSAs in query
iROI=tsa_key.VALUE[np.isin(tsa_key.Name,tsaList)].values

tsa_mask1=zTSA.copy()
tsa_mask1['Data']=np.zeros(zTSA['Data'].shape)
ind=(np.isin(zTSA['Data'],iROI))
tsa_mask1['Data'][ind]=1

# Define extent based on mask
xlim=[np.min(zTSA['X'][ind])-5000,np.max(zTSA['X'][ind])+5000]
ylim=[np.min(zTSA['Y'][ind])-5000,np.max(zTSA['Y'][ind])+5000]

# Clip mask
tsa_mask1=gis.ClipRaster(tsa_mask1,xlim,ylim)
gc.collect()

gdf_tsa_roi=gdf_tsa[np.isin(gdf_tsa.Name,tsaList)]
gdf_tsa_roi=gdf_tsa_roi.reset_index(drop=True)

#%% Isolate lakes inside the ROI

gdf_lakes_roi=gpd.overlay(gdf_bm[(gdf_bm['TAG']=='lake')],gdf_tsa[np.isin(gdf_tsa.Name,tsaList)],how='intersection')
#gdf_rivers_roi=gpd.overlay(gdf_bm[(gdf_bm['TAG']=='river')],gdf_tsa[np.isin(gdf_tsa.Name,tsaList)],how='intersection')

# Roads
try:
    gdf_road_roi=gdf_road.cx[tsa_mask1['minx']:tsa_mask1['maxx'],tsa_mask1['miny']:tsa_mask1['maxy']]
    gdf_road_roi=gdf_road_roi.reset_index(drop=True)
    gdf_road_roi=gpd.sjoin(gdf_road_roi,gdf_tsa_roi,how='left')
    gdf_road_roi=gdf_road_roi.groupby('index_right')
except:
    pass

#%% Import land cover scheme level 2 (Treed=4)

zLC2_tmp=gis.OpenGeoTiff(r'Z:\!Workgrp\Forest Carbon\Data\BC1ha\VRI\lc2.tif')
zLC2=zTSA.copy()
zLC2['Data']=zLC2_tmp['Data']
del zLC2_tmp
zLC2=gis.ClipRaster(zLC2,xlim,ylim)
gc.collect()

#%% Import Base Thematic Map

zBTM=gis.OpenGeoTiff(r'Z:\!Workgrp\Forest Carbon\Data\BC1ha\LandUseLandCover\landuse.btm.tif')
zBTM=gis.ClipRaster(zBTM,xlim,ylim)
dBTM=pd.read_csv(r'Z:\!Workgrp\Forest Carbon\Data\BC1ha\LandUseLandCover\landuse_btm_category_metadata.csv')
cl=np.column_stack( (dBTM['C1'].values,dBTM['C2'].values,dBTM['C3'].values) )
zBTM.Data1,zBTM.lab1,zBTM.cl1=gis.CompressCats(zBTM.Data,dBTM['Raster Value'].values,dBTM['PLU Label'].values,cl)


#%% PLOT TSA

# Grid and labels
lab=[]
z1=np.ones(zLC2['Data'].shape)
z1[(tsa_mask1['Data']==1) & (zLC2['Data']==4)]=0; lab.append('Treed')
z1[(tsa_mask1['Data']==1) & (zLC2['Data']!=4)]=1; lab.append('Non-treed')
z1[(tsa_mask1['Data']!=1) & (zLC2['Data']!=1)]=2; lab.append('hidden')
z1[(tsa_mask1['Data']!=1) & (zLC2['Data']==1)]=3; lab.append('hidden')

# Number of colours and number of colours excluded from colorbar
N_color=4
N_hidden=2

# Colormap
cm=np.vstack( ((0.55,0.55,0.55,1),(0.8,0.8,0.8,1),(0.93,0.93,0.93,1),(1,1,1,1)) )
cm=matplotlib.colors.ListedColormap(cm)

plt.close('all')
fig,ax=plt.subplots(1,2)
mngr=plt.get_current_fig_manager()
mngr.window.setGeometry(100,100,700,600)
im=ax[0].matshow(z1[0::1,0::1],clim=(0,N_color),extent=zLC2['Extent'],cmap=cm)
gdf_bc_boundary.plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
gdf_lakes_roi.plot(ax=ax[0],facecolor=[0.82,0.88,1],edgecolor=[0.7*0.82,0.7*0.88,0.7*1],linewidth=0.25,label='Water')
gdf_tsa_roi.plot(ax=ax[0],color=None,edgecolor=[0,0,0],facecolor='none')
gdf_road_roi.plot(ax=ax[0],facecolor='none',edgecolor=[0,0,0],label='Roads',linewidth=1.5,alpha=1,zorder=1)
ax[0].set(position=[0.04,0.02,0.92,0.96],xlim=xlim,ylim=ylim,aspect='auto')
ax[0].grid(False)

cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_color-1,1),ticks=np.arange(0.5,N_color+1.5,1))
cb.ax.set(yticklabels=lab)
cb.ax.tick_params(labelsize=6,length=0)
for i in range(0,N_color):
    ax[1].plot([0,100],[i/(N_color-N_hidden),i/(N_color-N_hidden)],'k-',linewidth=0.5)
ax[1].set(position=[0.8,0.04,0.03,0.06])

#gu.PrintFig(r'G:\My Drive\Figures\CutYear','png',500)




#%% PLOT BGC Zones within the TSA mask

# BGC classification
zBGC_tmp=gis.OpenGeoTiff(r'Z:\!Workgrp\Forest Carbon\Data\BC1ha\VRI\becz.tif')
zBGC=zTSA.copy()
zBGC['Data']=zBGC_tmp['Data']
del zBGC_tmp
zBGC['key']=pd.read_excel(r'Z:\!Workgrp\Forest Carbon\Data\BC1ha\VRI\becz.tif.vat.xlsx')
zBGC=gis.ClipRaster(zBGC,xlim,ylim)
gc.collect()

# Grid
u=np.unique(zBGC['Data'])
z1=np.ones(zBGC['Data'].shape)
for i in range(u.size):
    z1[(zBGC['Data']==u[i])]=i
L=i+1
z1[(tsa_mask1['Data']!=1)]=L

# Labels
lab=[]
for i in range(len(u)):
    try:
        lab.append(zBGC['key'].ZONE[zBGC['key'].VALUE==u[i]].values[0])
    except:
        lab.append('')

# Colormap
cm=plt.cm.get_cmap('viridis',L)
#cm.colors=np.concatenate((cm.colors,np.array([1,1,1,1]).reshape(1,4)),axis=0)

# Plot
plt.close('all')
fig,ax=plt.subplots(1,2)
mngr=plt.get_current_fig_manager()
mngr.window.setGeometry(100,100,750,750)
im=ax[0].matshow(z1[0::1,0::1],clim=(0,L+1),extent=zBGC['Extent'],cmap=cm)
gdf_tsa_roi.plot(ax=ax[0],color=None,edgecolor=[0,0,0],facecolor='none')
ax[0].set(position=[0.04,0.02,0.92,0.96],xlim=xlim,ylim=ylim,aspect='auto')
#ax[0].set(position=[0.04,0.02,0.92,0.96],xlim=[zTSA.minx,zTSA.maxx],ylim=[zTSA.miny,zTSA.maxy],aspect='auto')

cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,L+1,1),ticks=np.arange(0.5,L+0.5,1))
cb.ax.set(yticklabels=lab)
cb.ax.tick_params(labelsize=11,length=0)
for i in range(0,L):
    ax[1].plot([0,100],[i/L,i/L],'k-',linewidth=0.5)
ax[1].set(position=[0.06,0.04,0.025,0.4])



#%% Plot Burn Severity within the TSA mask

zBSH=zTSA.copy()
zBSH=gis.ClipRaster(zBSH,xlim,ylim)
zBSH['Data']=0*zBSH['Data']
for i in range(1,5):
    zBSH_tmp=gis.OpenGeoTiff(r'Z:\!Workgrp\Forest Carbon\Data\BC1ha\Disturbances\VEG_BURN_SEVERITY_SP_High_year' + str(i) + '.tif')    
    zBSH_tmp=gis.ClipRaster(zBSH_tmp,xlim,ylim)    
    zBSH['Data']=np.maximum(zBSH['Data'],zBSH_tmp['Data'])
    del zBSH_tmp
    gc.collect()

zBSM=zTSA.copy()
zBSM=gis.ClipRaster(zBSM,xlim,ylim)
zBSM['Data']=0*zBSM['Data']
for i in range(1,5):
    zBSM_tmp=gis.OpenGeoTiff(r'Z:\!Workgrp\Forest Carbon\Data\BC1ha\Disturbances\VEG_BURN_SEVERITY_SP_Medium_year' + str(i) + '.tif')
    zBSM_tmp=gis.ClipRaster(zBSM_tmp,xlim,ylim)    
    zBSM['Data']=np.maximum(zBSM['Data'],zBSM_tmp['Data'])
    del zBSM_tmp
    gc.collect()

zBSL=zTSA.copy()
zBSL=gis.ClipRaster(zBSL,xlim,ylim)
zBSL['Data']=0*zBSL['Data']
for i in range(1,5):
    zBSL_tmp=gis.OpenGeoTiff(r'Z:\!Workgrp\Forest Carbon\Data\BC1ha\Disturbances\VEG_BURN_SEVERITY_SP_Low_year' + str(i) + '.tif')
    zBSL_tmp=gis.ClipRaster(zBSL_tmp,xlim,ylim)    
    zBSL['Data']=np.maximum(zBSL['Data'],zBSL_tmp['Data'])
    del zBSL_tmp
    gc.collect()

# Grid and labels
lab=[]
z1=8*np.ones(zLC2['Data'].shape)
z1[(tsa_mask1['Data']==1) & (zLC2['Data']==4) & (zBSH['Data']==2017)]=0; lab.append('Treed, High severity')
z1[(tsa_mask1['Data']==1) & (zLC2['Data']==4) & (zBSM['Data']==2017)]=1; lab.append('Treed, Medium severity')
z1[(tsa_mask1['Data']==1) & (zLC2['Data']==4) & (zBSL['Data']==2017)]=2; lab.append('Treed, Low severity')
z1[(tsa_mask1['Data']==1) & (zLC2['Data']!=4) & (zBSH['Data']==2017)]=3; lab.append('Non-treed, High severity')
z1[(tsa_mask1['Data']==1) & (zLC2['Data']!=4) & (zBSM['Data']==2017)]=4; lab.append('Non-treed, Medium severity')
z1[(tsa_mask1['Data']==1) & (zLC2['Data']!=4) & (zBSL['Data']==2017)]=5; lab.append('Non-treed, Low severity')
z1[(tsa_mask1['Data']==1) & (zLC2['Data']==4) & (zBSL['Data']!=2017) & (zBSM['Data']!=2017) & (zBSH['Data']!=2017)]=6; lab.append('Treed, unburned')
z1[(tsa_mask1['Data']==1) & (zLC2['Data']!=4) & (zBSL['Data']!=2017) & (zBSM['Data']!=2017) & (zBSH['Data']!=2017)]=7; lab.append('Non-treed, unburned')
#z1[(tsa_mask1['Data']!=1) & (zLC2['Data']!=1)]=7; lab.append('hidden')
z1[(tsa_mask1['Data']!=1) & (zLC2['Data']==1)]=8; lab.append('hidden')

# Number of colours and number of colours excluded from colorbar
N_color=9
N_hidden=1

# Colormap
cm=np.vstack( ( (0.25,0,0,1),(0.625,0,0,1),(1,0,0,1), \
                (1,0.5,0,1),(1,0.75,0,1),(1,1,0,1), \
                (0.93,0.93,0.93,1),(0.8,0.8,0.8,1),(1,1,1,1) ) )
cm=matplotlib.colors.ListedColormap(cm)


plt.close('all')
fig,ax=plt.subplots(1,2)
mngr=plt.get_current_fig_manager()
mngr.window.setGeometry(100,100,700,600)
im=ax[0].matshow(z1,clim=(0,N_color),extent=zLC2['Extent'],cmap=cm)
gdf_bc_boundary.plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
gdf_lakes_roi.plot(ax=ax[0],facecolor=[0.82,0.88,1],edgecolor=[0.7*0.82,0.7*0.88,0.7*1],linewidth=0.25,label='Water')
#gdf_rivers_roi.plot(ax=ax[0],linecolor=[0,0,0.7],label='Water',linewidth=0.25)
gdf_tsa_roi.plot(ax=ax[0],color=None,edgecolor=[0,0,0],facecolor='none')
rc1=gdf_road_roi.plot(ax=ax[0],facecolor='none',edgecolor=[0,0,0],label='Roads',linewidth=1.5,alpha=1,zorder=1)
rc2=gdf_road_roi.plot(ax=ax[0],facecolor='none',edgecolor=[1,1,1],label='Roads',linewidth=1,alpha=1,zorder=2)
ax[0].set(position=[0.04,0.02,0.92,0.96],xlim=xlim,ylim=ylim,aspect='auto')
ax[0].grid(False)

cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_color,1),ticks=np.arange(0.5,N_color+1.5,1))
cb.ax.set(yticklabels=lab)
cb.ax.tick_params(labelsize=6,length=0)
for i in range(0,N_color):
    ax[1].plot([0,100],[i/(N_color-N_hidden),i/(N_color-N_hidden)],'k-',linewidth=0.5)
ax[1].set(position=[0.78,0.04,0.03,0.2])

#gu.PrintFig(r'G:\My Drive\Figures\CaribooPlantOverBurns','png',300)


# Add invenotry data

# Inventory data
sgrd=gu.ipickle(r'C:\Users\rhember\Documents\Data\FCI_Projects\FCI_RollupByTile1\TileBC_0909\Inputs\Geospatial\Coordinates.pkl')
atu=gu.ipickle(r'C:\Users\rhember\Documents\Data\FCI_Projects\FCI_RollupByTile1\TileBC_0909\Inputs\Geospatial\RSLT_ACTIVITY_TREATMENT_SVW.pkl')
pest=gu.ipickle(r'C:\Users\rhember\Documents\Data\FCI_Projects\FCI_RollupByTile1\TileBC_0909\Inputs\Geospatial\PEST_INFESTATION_POLY.pkl')
vri=gu.ipickle(r'C:\Users\rhember\Documents\Data\FCI_Projects\FCI_RollupByTile1\TileBC_0909\Inputs\Geospatial\VEG_COMP_LYR_R1_POLY.pkl')
lut_at=gu.ipickle(r'C:\Users\rhember\Documents\Data\ForestInventory\Results\20200430\LUTs_RSLT_ACTIVITY_TREATMENT_SVW.pkl')
lut_pest=gu.ipickle(r'C:\Users\rhember\Documents\Data\ForestInventory\Disturbances\20200430\LUTs_PEST_INFESTATION_POLY.pkl')
xg=sgrd['X'].flatten()
yg=sgrd['Y'].flatten()


ind=np.where( (atu['SILV_BASE_CODE']==lut_at['SILV_BASE_CODE']['PL']) & (atu['Year']>=2017) )[0]

ax[0].plot(xg[atu['IdxToGrid'][ind]],yg[atu['IdxToGrid'][ind]],'o',markeredgecolor=[0,1,0],markerfacecolor=[0.7,1,0.7],markersize=3,linewidth=0.25)




#%% PLOT stand age within the TSA mask

zAge_tmp=gis.OpenGeoTiff(r'Z:\!Workgrp\Forest Carbon\Data\BC1ha\VRI\age1.tif')
zAge=zTSA.copy()
zAge['Data']=zAge_tmp['Data']
zAge=gis.ClipRaster(zAge,xlim,ylim)
del zAge_tmp
gc.collect()

# Grid
bin=np.arange(0,210,10); bw=10;
z1=(bin.size)*np.ones(zAge['Data'].shape)
for i in range(bin.size):
    ind=np.where(np.abs(zAge['Data']-bin[i])<=bw/2)
    z1[ind]=i
z1[(tsa_mask1['Data']==1) & (zLC2['Data']!=4)]=i+1
z1[(tsa_mask1['Data']!=1)]=i+2
L=i+2

ind=np.where( (tsa_mask1['Data']==1) & (zLC2['Data']==4) )
np.mean(zAge['Data'][ind])

# Labels
lab=bin.astype(str)

# Colormap
#cm=plt.cm.get_cmap('viridis',i)
cm=plt.cm.get_cmap('plasma',i)
cm=np.vstack( (cm.colors,(0.9,0.9,0.9,1),(1,1,1,1)) )
cm=matplotlib.colors.ListedColormap(cm)
                              
# Plot
plt.close('all')
fig,ax=plt.subplots(1,2)
mngr=plt.get_current_fig_manager()
mngr.window.setGeometry(100,100,950,750)
im=ax[0].matshow(z1[0::1,0::1],clim=(0,L+1),extent=zBGC['Extent'],cmap=cm)
gdf_tsa_roi.plot(ax=ax[0],color=None,edgecolor=[0,0,0],facecolor='none')
ax[0].set(position=[0.04,0.02,0.92,0.96],xlim=xlim,ylim=ylim,aspect='auto')
ax[0].grid(False)
#ax[0].set(position=[0.04,0.02,0.92,0.96],xlim=[zTSA.minx,zTSA.maxx],ylim=[zTSA.miny,zTSA.maxy],aspect='auto')

cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,L,1),ticks=np.arange(0.5,L+1.5,1))
cb.ax.set(yticklabels=lab)
cb.ax.tick_params(labelsize=6,length=0)
for i in range(0,L):
    ax[1].plot([0,100],[i/(L-1),i/(L-1)],'k-',linewidth=0.5)
ax[1].set(position=[0.8,0.04,0.025,0.5])



#%% PLOT harvested year within the TSA mask

zH=zTSA.copy()
zH=gis.ClipRaster(zH,xlim,ylim)
zH['Data']=0*zH['Data']
for i in range(1,5):
    zH_tmp=gis.OpenGeoTiff(r'Z:\!Workgrp\Forest Carbon\Data\BC1ha\Disturbances\VEG_CONSOLIDATED_CUT_BLOCKS_SP_year' + str(i) + '.tif')    
    zH_tmp=gis.ClipRaster(zH_tmp,xlim,ylim)    
    zH['Data']=np.maximum(zH['Data'],zH_tmp['Data'])
    del zH_tmp
    gc.collect()

zLC2_tmp=gis.OpenGeoTiff(r'Z:\!Workgrp\Forest Carbon\Data\BC1ha\VRI\lc2.tif')
zLC2=zTSA.copy()
zLC2['Data']=zLC2_tmp['Data']
del zLC2_tmp
zLC2=gis.ClipRaster(zLC2,xlim,ylim)
gc.collect()

# Grid
bw=5; bin=np.arange(1960,2025,bw); 
z1=(bin.size)*np.ones(zH['Data'].shape)
for i in range(bin.size):
    ind=np.where(np.abs(zH['Data']-bin[i])<=bw/2)
    z1[ind]=i
z1[(tsa_mask1['Data']==1) & (zLC2['Data']==4) & (zH['Data']==0)]=i+1
z1[(tsa_mask1['Data']==1) & (zLC2['Data']!=4)]=i+2
z1[(tsa_mask1['Data']!=1)]=i+3
L=i+3

# Labels
lab=bin.astype(str)

# Colormap
#cm=plt.cm.get_cmap('viridis',i)
cm=plt.cm.get_cmap('plasma',i)
cm=np.vstack( (cm.colors,(0.8,0.8,0.8,1),(0.9,0.9,0.9,1),(1,1,1,1)) )
cm=matplotlib.colors.ListedColormap(cm)
                              
# Plot
plt.close('all')
fig,ax=plt.subplots(1,2)
mngr=plt.get_current_fig_manager()
mngr.window.setGeometry(100,100,950,750)
im=ax[0].matshow(z1[0::1,0::1],clim=(0,L+1),extent=zBGC['Extent'],cmap=cm)
gdf_tsa_roi.plot(ax=ax[0],color=None,edgecolor=[0,0,0],facecolor='none')
ax[0].set(position=[0.04,0.02,0.92,0.96],xlim=xlim,ylim=ylim,aspect='auto')
ax[0].grid(False)
#ax[0].set(position=[0.04,0.02,0.92,0.96],xlim=[zTSA.minx,zTSA.maxx],ylim=[zTSA.miny,zTSA.maxy],aspect='auto')

cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,L-1,1),ticks=np.arange(0.5,L+1.5,1))
cb.ax.set(yticklabels=lab)
cb.ax.tick_params(labelsize=6,length=0)
for i in range(0,L):
    ax[1].plot([0,100],[i/(L-2),i/(L-2)],'k-',linewidth=0.5)
ax[1].set(position=[0.8,0.04,0.025,0.35])
gu.PrintFig(r'G:\My Drive\Figures\CutYear','png',500)


#%% Plot AOS year within the TSA mask

dSC=gu.ReadExcel(r'Z:\!Workgrp\Forest Carbon\Data\BC1ha\Disturbances\PEST_INFESTATION_POLY_PEST_SEVERITY_CODE.xlsx')
sco=['L','M','S','V']
tv=np.arange(1990,2021,1)
A=np.zeros((tv.size,len(sco)))
for iT in range(tv.size):
    print(tv[iT])
    z=gis.OpenGeoTiff(r'Z:\!Workgrp\Forest Carbon\Data\BC1ha\Disturbances\PEST_INFESTATION_POLY_IDL_SeverityClass_' + str(tv[iT]) + '.tif')
    for iSC in range(len(sco)):
        ind=np.where(dSC['PEST_SEVERITY_CODE']==sco[iSC])[0]
        ind=np.where(z['Data']==dSC['ID'][ind])[0]
        A[iT,iSC]=ind.size

plt.close('all')
plt.plot(tv,A[:,1],'-s')        
plt.plot(tv,A[:,2],'-o')
plt.plot(tv,A[:,3],'-d')


zLC2_tmp=gis.OpenGeoTiff(r'Z:\!Workgrp\Forest Carbon\Data\BC1ha\VRI\lc2.tif')
zLC2=zTSA.copy()
zLC2['Data']=zLC2_tmp['Data']
del zLC2_tmp
zLC2=gis.ClipRaster(zLC2,xlim,ylim)
gc.collect()

# Grid
bw=5; bin=np.arange(1960,2025,bw); 
z1=(bin.size)*np.ones(zH['Data'].shape)
for i in range(bin.size):
    ind=np.where(np.abs(zH['Data']-bin[i])<=bw/2)
    z1[ind]=i
z1[(tsa_mask1['Data']==1) & (zLC2['Data']==4) & (zH['Data']==0)]=i+1
z1[(tsa_mask1['Data']==1) & (zLC2['Data']!=4)]=i+2
z1[(tsa_mask1['Data']!=1)]=i+3
L=i+3

# Labels
lab=bin.astype(str)

# Colormap
#cm=plt.cm.get_cmap('viridis',i)
cm=plt.cm.get_cmap('plasma',i)
cm=np.vstack( (cm.colors,(0.8,0.8,0.8,1),(0.9,0.9,0.9,1),(1,1,1,1)) )
cm=matplotlib.colors.ListedColormap(cm)
                              
# Plot
plt.close('all')
fig,ax=plt.subplots(1,2)
mngr=plt.get_current_fig_manager()
mngr.window.setGeometry(100,100,950,750)
im=ax[0].matshow(z1[0::1,0::1],clim=(0,L+1),extent=zBGC['Extent'],cmap=cm)
gdf_tsa_roi.plot(ax=ax[0],color=None,edgecolor=[0,0,0],facecolor='none')
ax[0].set(position=[0.04,0.02,0.92,0.96],xlim=xlim,ylim=ylim,aspect='auto')
ax[0].grid(False)
#ax[0].set(position=[0.04,0.02,0.92,0.96],xlim=[zTSA.minx,zTSA.maxx],ylim=[zTSA.miny,zTSA.maxy],aspect='auto')

cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,L-1,1),ticks=np.arange(0.5,L+1.5,1))
cb.ax.set(yticklabels=lab)
cb.ax.tick_params(labelsize=6,length=0)
for i in range(0,L):
    ax[1].plot([0,100],[i/(L-2),i/(L-2)],'k-',linewidth=0.5)
ax[1].set(position=[0.8,0.04,0.025,0.35])
gu.PrintFig(r'G:\My Drive\Figures\CutYear','png',500)




#%% PLOT wildfire within the TSA mask

zF1=zTSA.copy()
zF1=gis.ClipRaster(zF1,xlim,ylim)
zF1['Data']=0*zF1['Data']
for i in range(1,5):
    zF1_tmp=gis.OpenGeoTiff(r'Z:\!Workgrp\Forest Carbon\Data\BC1ha\Disturbances\PROT_HISTORICAL_FIRE_POLYS_SP_year' + str(i) + '.tif')    
    zF1_tmp=gis.ClipRaster(zF1_tmp,xlim,ylim)
    ind=np.where(zF1_tmp['Data']==2017)
    zF1['Data'][ind]=zF1['Data'][ind]+zF1_tmp['Data'][ind]
    del zF1_tmp    
    gc.collect()

# Grid
z1=np.ones(zF1['Data'].shape)
z1[(zF1['Data']==2017)]=2
L=3
z1[(tsa_mask1['Data']!=1)]=L

# Labels
lab=[]
for i in range(len(u)):
    try:
        lab.append(zF1['key'].ZONE[zF1['key'].VALUE==u[i]].values[0])
    except:
        lab.append('')

# Colormap
cm=plt.cm.get_cmap('viridis',L)
#cm.colors=np.concatenate((cm.colors,np.array([1,1,1,1]).reshape(1,4)),axis=0)
cm.colors[0,:]=[0.8,0.8,0.8,1]
cm.colors[1,:]=[0.5,0,0,1]
cm.colors[2,:]=[1,1,1,1]

# Plot
plt.close('all')
fig,ax=plt.subplots(1,2)
mngr=plt.get_current_fig_manager()
mngr.window.setGeometry(100,100,750,750)
im=ax[0].matshow(z1[0::1,0::1],clim=(0,L+1),extent=zF1['Extent'],cmap=cm)
gdf_tsa_roi.plot(ax=ax[0],color=None,edgecolor=[0,0,0],facecolor='none')
ax[0].set(position=[0.04,0.02,0.92,0.96],xlim=xlim,ylim=ylim,aspect='auto')
#ax[0].set(position=[0.04,0.02,0.92,0.96],xlim=[zTSA.minx,zTSA.maxx],ylim=[zTSA.miny,zTSA.maxy],aspect='auto')

cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,L+1,1),ticks=np.arange(0.5,L+0.5,1))
cb.ax.set(yticklabels=lab)
cb.ax.tick_params(labelsize=11,length=0)
for i in range(0,L):
    ax[1].plot([0,100],[i/L,i/L],'k-',linewidth=0.5)
ax[1].set(position=[0.06,0.04,0.025,0.4])





#%% Map of N deposition from 2017 aerial fertilization

pth=r'C:\Users\rhember\Documents\Data\FCI_Projects\FertilizationSummary\Inputs\Nitrogen Deposition\ndep_2017.tif'
z0=gis.OpenGeoTiff(pth)

Mask=np.zeros(z0['Data'].shape)
id=np.where(zBTM.lab1=='Old Forest (FO)')[0]; Mask[(zBTM['Data1']==id)]=1
id=np.where(zBTM.lab1=='Range Lands (RANG)')[0]; Mask[(zBTM['Data1']==id)]=1
id=np.where(zBTM.lab1=='Recently Burned (BURN)')[0]; Mask[(zBTM['Data1']==id)]=1
id=np.where(zBTM.lab1=='Recently Logged (LOG)')[0]; Mask[(zBTM['Data1']==id)]=1
id=np.where(zBTM.lab1=='Selectively Logged (LOGS)')[0]; Mask[(zBTM['Data1']==id)]=1
id=np.where(zBTM.lab1=='Shurbs (SHRB)')[0]; Mask[(zBTM['Data1']==id)]=1
id=np.where(zBTM.lab1=='Young Forest (FY)')[0]; Mask[(zBTM['Data1']==id)]=1

N0=np.sum(z0['Data'].flatten())/25000
N1=np.sum(z0['Data'].flatten()*Mask.flatten())/25000
N1/N0

N1=np.sum(z0['Data'].flatten()*Mask.flatten())
3.667*35*0.8*N1/1000/2500

#Mask=np.ones(z0['Data'].shape)

# Grid
L=8
z1=L*np.ones(zTSA['Data'].shape,dtype='int8')
ind=np.where((Mask==1) & (z0['Data']>0.001) & (z0['Data']<=0.01)); z1[ind]=1
ind=np.where((Mask==1) & (z0['Data']>0.01) & (z0['Data']<=0.05)); z1[ind]=2
ind=np.where((Mask==1) & (z0['Data']>0.05) & (z0['Data']<=0.1)); z1[ind]=3
ind=np.where((Mask==1) & (z0['Data']>0.1) & (z0['Data']<=0.2)); z1[ind]=4
ind=np.where((Mask==1) & (z0['Data']>0.2) & (z0['Data']<=0.5)); z1[ind]=5
ind=np.where((Mask==1) & (z0['Data']>0.5) & (z0['Data']<=1)); z1[ind]=6
ind=np.where((Mask==1) & (z0['Data']>1) ); z1[ind]=7
z1[(Mask==0)]=8
z1[(z0['Data']==0)]=8
#id=np.where(zBTM.lab1=='Unknown')[0]; z1[(zBTM['Data1']==id)]=8
#id=np.where(zBTM.lab1=='Outside B.C. (NONE)')[0]; z1[(zBTM['Data1']==id)]=8

pth=r'C:\Users\rhember\Documents\Data\FCI_Projects\FertilizationSummary\Inputs\Nitrogen Deposition\ndep_2017_cat.tif'
zo=z0.copy()
zo['Data']=z1.astype('int16')
gis.SaveGeoTiff(zo,pth)

# Labels
lab=['0.005','0.025','0.075','0.150','0.350','0.750','1.000']

# Colormap
cm=plt.cm.get_cmap('plasma',7)
cm=np.vstack( (cm.colors,(1,1,1,1)) ) # ,(1,1,1,1)
cm=matplotlib.colors.ListedColormap(cm)

from rasterio.plot import show
show(z1)
                              
# Plot
plt.close('all')
fig,ax=plt.subplots(1,2)
mngr=plt.get_current_fig_manager()
mngr.window.setGeometry(100,100,600,600)
im=ax[0].matshow(z1[0::50,0::50])
im=ax[0].matshow(z1[0::50,0::50],clim=(0,L),extent=zTSA['Extent'],cmap=cm)
#gdf_tsa_roi.plot(ax=ax[0],color=None,edgecolor=[0,0,0],facecolor='none')
gdf_bc_boundary.plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.5)
ax[0].set(position=[0,0,1,1],xlim=xlim,ylim=ylim,aspect='auto',xticks=xlim,yticks=ylim)
ax[0].grid(False)
ax[0].tick_params(left=False, right=False)

cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,L-1,1),ticks=np.arange(0.5,L+1.5,1))
cb.ax.set(yticklabels=lab)
cb.ax.tick_params(labelsize=13,length=0)
for i in range(0,L):
    ax[1].plot([0,100],[i/(L-2),i/(L-2)],'k-',linewidth=0.5)
ax[1].set(position=[0.8,0.54,0.03,0.3])
gu.PrintFig(r'G:\My Drive\Figures\Fertilization\ndep_from_aerial_fert_2017','png',500)

