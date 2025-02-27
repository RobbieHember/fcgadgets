#%% Import modules
import os
import numpy as np
import time
import matplotlib.pyplot as plt
import geopandas as gpd
import pandas as pd
import pyproj
from pyproj import Proj, transform
import scipy.io as spio
import fiona
import rasterio
from rasterio import features
from shapely import geometry
from scipy.interpolate import griddata
import cv2
import matplotlib.patches as patches
import matplotlib.colors
from pyproj import Proj, transform
#from rasterio.transform import from_origin
#from shapely.geometry import Point,Polygon
import fcgadgets.macgyver.util_general as gu
import fcgadgets.macgyver.util_gis as gis
import fcgadgets.macgyver.util_query_gdb as qgdb
import fcgadgets.bc1ha.bc1ha_utils as u1ha
import fcgadgets.cbrunner.cbrun_util as cbu
gp=gu.SetGraphics('Manuscript')

#%% Import paths and look-up-tables
meta=u1ha.Init()
zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
gdf=u1ha.Import_GDBs_ProvinceWide(meta)
srs=gis.ImportSRSs()

#%% Select ROIs
cm=plt.cm.get_cmap('binary',2)
plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(16,12))
im=ax.matshow(zRef['Data'],extent=zRef['Extent'],clim=(0,2));
gdf['bc_bound']['gdf'].plot(ax=ax,color='None',edgecolor=[1,1,1],lw=0.25,fc='None')
gdf['popp']['gdf'].plot(ax=ax,marker='s',markersize=5,facecolor=[1,1,0],edgecolor=[0,0,0],linewidth=0.25,label='Opening',alpha=1)
#ax[0].set(position=[0.075,0.08,0.8,0.88],xlim=zTSA['xlim'],ylim=zTSA['ylim'],aspect='auto');
#ax[0].xaxis.tick_bottom()
#ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].tick_params(length=1.5)
#cb=plt.colorbar(im,cax=ax[1])
#ax[1].set(position=[0.89,0.2,0.03,0.6])

def on_click(event):
	if event.inaxes is not None:
		print(event.xdata,event.ydata)
	else:
		print('Clicked ouside axes bounds but inside plot window')

flg=1
if flg==1:
	fig.canvas.callbacks.connect('button_press_event', on_click)
	plt.show()

# project coordinates
flg=0
if flg==1:
	lat=49.562844; lon=-115.763263; x,y=transform(srs['Proj']['Geographic'],srs['Proj']['BC1ha'],lon,lat); # Cranbrook
	lat=59.139484; lon=-122.419605; x,y=transform(srs['Proj']['Geographic'],srs['Proj']['BC1ha'],lon,lat); # NE BC
	lat=56.194609; lon=-120.912890; x,y=transform(srs['Proj']['Geographic'],srs['Proj']['BC1ha'],lon,lat); # Site C Dam
	lat=50.033426; lon=-125.294933; x,y=transform(srs['Proj']['Geographic'],srs['Proj']['BC1ha'],lon,lat); # Campbell River
	lat=54.790580; lon=-127.181148; x,y=transform(srs['Proj']['Geographic'],srs['Proj']['BC1ha'],lon,lat); # Smithers	
	#lon,lat=transform(srs['Proj']['BC1ha'],srs['Proj']['Geographic'],x,y)

#%% Select ROIs
cm=plt.cm.get_cmap('binary',2)
plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(16,12))
im=ax.matshow(zRef['Data'],extent=zRef['Extent'],clim=(0,2));
gdf['bc_bound']['gdf'].plot(ax=ax,color='None',edgecolor=[1,1,1],lw=0.25,fc='None')
gdf['popp']['gdf'].plot(ax=ax,marker='s',markersize=5,facecolor=[1,1,0],edgecolor=[0,0,0],linewidth=0.25,label='Opening',alpha=1)

Sites=[]
Sites.append({'x':1738380,'y':558653,'Name':'Cranbrook'})
Sites.append({'x':1314813,'y':1265000,'Name':'Fort St John'})
Sites.append({'x':1180000,'y':397000,'Name':'Victoria'})
Sites.append({'x':1265338,'y':798420,'Name':'Williams Lake'})
Sites.append({'x':1050512,'y':557753,'Name':'Campbell River'})
Sites.append({'x':924231,'y':1088485,'Name':'Smithers'})
w=30000
for i in range(len(Sites)):
	xy=(Sites[i]['x']-w/2 ,Sites[i]['y']-w/2) 
	rect=patches.Rectangle(xy,w,w,fc='None',ec=[0.55,1,0],lw=1)
	ax.add_patch(rect)
	ax.text(Sites[i]['x']+1.25*w,Sites[i]['y']+1.25*w,str(i+1),color=[1,1,1],fontweight='bold',fontsize=10,ha='center',va='center')

#%% Plot LC and LU for ROIs 

z=u1ha.Import_Raster(meta,[],['lc_comp1_2019','lu_comp1_2019'])

#%% Plot LC/LU ROIs Example 1

Loci=[0,0,1,1,2,2]
plt.close('all'); fig,ax=plt.subplots(8,figsize=gu.cm2inch(18,17.5))
for i in range(6):
	xlim=[Sites[Loci[i]]['x']-w/2,Sites[Loci[i]]['x']+w/2]
	ylim=[Sites[Loci[i]]['y']-w/2,Sites[Loci[i]]['y']+w/2]
	zRef0=gis.ClipRasterByXYLimits(zRef,xlim,ylim)
	
	if np.isin(i,[0,2,4,6])==True:
		lab=list(meta['LUT']['Derived']['lc_comp1'].keys())
		z0=gis.ClipRasterByXYLimits(z['lc_comp1_2019'],xlim,ylim)
		z1=len(lab)*np.ones(z0['Data'].shape,dtype='int8')
		for k in meta['LUT']['Derived']['lc_comp1'].keys():
			ind=np.where(z0['Data']==meta['LUT']['Derived']['lc_comp1'][k])
			if ind[0].size>0:
				z1[ind]=meta['LUT']['Derived']['lc_comp1'][k]
			else:
				z1[0,meta['LUT']['Derived']['lc_comp1'][k]]=meta['LUT']['Derived']['lc_comp1'][k]		
		ind=np.where(zRef0['Data']==0)
		if ind[0].size>0:
			z1[ind]=meta['LUT']['Derived']['lc_comp1']['Water']
		#else:
		z1[-1,-1]=meta['LUT']['Derived']['lc_comp1'][k]+1

		N_vis=len(lab)
		N_hidden=1
		N_tot=N_vis+N_hidden
		cm=np.vstack( ((0.3,0.6,0,1),(0.65,1,0,1),(1,1,0.5,1),(0.8,0.7,1,1),(0.75,0.5,0.5,1),(0.75,0.75,0.75,1),(0.8,0,0,1),(0.95,0.95,0.95,1),(0.75,0.85,0.95,1),(1,1,1,1)) )
		cm=matplotlib.colors.ListedColormap(cm)		
		im1=ax[i].matshow(z1,extent=z0['Extent'],cmap=cm)
		
	else:
		z0=gis.ClipRasterByXYLimits(z['lu_comp1_2019'],xlim,ylim)
		lab2=list(meta['LUT']['Derived']['lu_comp1'].keys())
		z2=len(lab2)*np.ones(z0['Data'].shape,dtype='int8')
		for k in meta['LUT']['Derived']['lu_comp1'].keys():
			ind=np.where(z0['Data']==meta['LUT']['Derived']['lu_comp1'][k])
			z2[ind]=meta['LUT']['Derived']['lu_comp1'][k]
		ind=np.where(zRef0==0)
		if ind[0].size>0:
			z2[ind]=meta['LUT']['Derived']['lu_comp1'][k]+1
		else:
			z2[-1,-1]=meta['LUT']['Derived']['lu_comp1'][k]+1
		for k in meta['LUT']['Derived']['lu_comp1'].keys():
			z2[0,meta['LUT']['Derived']['lu_comp1'][k]]=meta['LUT']['Derived']['lu_comp1'][k]
		
		N_vis2=len(lab2)
		N_hidden2=1
		N_tot2=N_vis2+N_hidden2
		cm2=np.vstack( ((0.75,0,0,1),(0.9,0.8,0.4,1),(0.5,0.8,0.1,1),(0,0.2,0.5,1),(0.425,0,0,1),(0.2,0.5,0,1),(0.35,0.87,0.94,1),(0.85,0.95,0.7,1),(0.8,0.6,1,1),(1,0.5,0,1),(1,1,0,1),(0.4,0.7,0.79,1),(0.6,0.6,0.6,1),(1,1,1,1)) )
		cm2=matplotlib.colors.ListedColormap(cm2)
		im2=ax[i].matshow(z2,extent=z0['Extent'],cmap=cm2)			
		
	ax[i].set(xlim=z0['xlim'],ylim=z0['ylim'],xticks=xlim,yticks=ylim,xticklabels=['',''],yticklabels=['',''])
	#ax[i].yaxis.set_ticks_position('both'); ax[i].xaxis.set_ticks_position('both'); 
	ax[i].grid(meta['Graphics']['Map']['Map Grid Vis']); #ax[i].axis(meta['Graphics']['Map']['Map Axis Vis'])
	ax[i].text(Sites[Loci[i]]['x'],Sites[Loci[i]]['y'],Sites[Loci[i]]['Name'],color=[0,0,0],fontweight='bold',fontsize=10,ha='center',va='center')
ax[0].set(position=[-0.025,0.668,0.4,0.325])
ax[1].set(position=[0.43,0.668,0.4,0.325])
ax[2].set(position=[-0.025,0.337,0.4,0.325])
ax[3].set(position=[0.43,0.337,0.4,0.325])
ax[4].set(position=[-0.025,0.005,0.4,0.325])
ax[5].set(position=[0.43,0.005,0.4,0.325])

zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
cb=plt.colorbar(im1,cax=ax[6],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
cb.ax.set(yticklabels=lab)
cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
cb.outline.set_edgecolor('w')
for j in range(cb_bnd.size):
	ax[6].plot([0,100],[cb_bnd[j],cb_bnd[j]],'w-',linewidth=2)
ax[6].set(position=[0.34,0.98-(0.018*len(lab)),0.02,0.018*len(lab)])

zmn=np.min(z2); zmx=np.max(z2); cb_ivl=(zmx-zmn)/N_tot2; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden2*cb_ivl,cb_ivl)
cb_ticks=np.arange(zmn+cb_ivl/2,N_tot2-1,cb_ivl)
cb=plt.colorbar(im2,cax=ax[7],cmap=cm2,boundaries=cb_bnd,ticks=cb_ticks)
cb.ax.set(yticklabels=lab2)
cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
cb.outline.set_edgecolor('w')
for j in range(cb_bnd.size):
	ax[7].plot([0,100],[cb_bnd[j],cb_bnd[j]],'w-',linewidth=2)
ax[7].set(position=[0.795,0.98-(0.018*len(lab2)),0.02,0.018*len(lab2)])

gu.axletters(ax,plt,-0.04,0.96,FontColor=meta['Graphics']['gp']['cla'],FontSize=8,LetterStyle='NoPar',FontWeight='Bold',Skip=[6,7])
gu.axletters(ax,plt,-0.04,0.96,FontColor=meta['Graphics']['gp']['cla'],FontSize=8,LetterStyle='NoPar',FontWeight='Bold',Skip=[6,7])
gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\LandCover\LCLU_Examples1','png',900)

#%% Plot LC / LU for ROIs Example 2
Loci=[3,3,4,4,5,5]
plt.close('all'); fig,ax=plt.subplots(8,figsize=gu.cm2inch(18,17.5))
for i in range(6):
	xlim=[Sites[Loci[i]]['x']-w/2,Sites[Loci[i]]['x']+w/2]
	ylim=[Sites[Loci[i]]['y']-w/2,Sites[Loci[i]]['y']+w/2]
	zRef0=gis.ClipRasterByXYLimits(zRef,xlim,ylim)
	
	if np.isin(i,[0,2,4,6])==True:
		lab=list(meta['LUT']['Derived']['lc_comp1'].keys())
		z0=gis.ClipRasterByXYLimits(z['lc_comp1_2019'],xlim,ylim)
		z1=len(lab)*np.ones(z0['Data'].shape,dtype='int8')
		for k in meta['LUT']['Derived']['lc_comp1'].keys():
			ind=np.where(z0['Data']==meta['LUT']['Derived']['lc_comp1'][k])
			if ind[0].size>0:
				z1[ind]=meta['LUT']['Derived']['lc_comp1'][k]
			else:
				z1[0,meta['LUT']['Derived']['lc_comp1'][k]]=meta['LUT']['Derived']['lc_comp1'][k]		
		ind=np.where(zRef0['Data']==0)
		if ind[0].size>0:
			z1[ind]=meta['LUT']['Derived']['lc_comp1']['Water']
		#else:
		z1[-1,-1]=meta['LUT']['Derived']['lc_comp1'][k]+1

		N_vis=len(lab)
		N_hidden=1
		N_tot=N_vis+N_hidden
		cm=np.vstack( ((0.3,0.6,0,1),(0.65,1,0,1),(1,1,0.5,1),(0.8,0.7,1,1),(0.75,0.5,0.5,1),(0.75,0.75,0.75,1),(0.8,0,0,1),(0.95,0.95,0.95,1),(0.75,0.85,0.95,1),(1,1,1,1)) )
		cm=matplotlib.colors.ListedColormap(cm)		
		im1=ax[i].matshow(z1,extent=z0['Extent'],cmap=cm)			
		
	else:
		z0=gis.ClipRasterByXYLimits(z['lu_comp1_2019'],xlim,ylim)
		lab2=list(meta['LUT']['Derived']['lu_comp1'].keys())
		z2=len(lab2)*np.ones(z0['Data'].shape,dtype='int8')
		for k in meta['LUT']['Derived']['lu_comp1'].keys():
			ind=np.where(z0['Data']==meta['LUT']['Derived']['lu_comp1'][k])
			z2[ind]=meta['LUT']['Derived']['lu_comp1'][k]
		ind=np.where(zRef0==0)
		if ind[0].size>0:
			z2[ind]=meta['LUT']['Derived']['lu_comp1'][k]+1
		else:
			z2[-1,-1]=meta['LUT']['Derived']['lu_comp1'][k]+1
		for k in meta['LUT']['Derived']['lu_comp1'].keys():
			z2[0,meta['LUT']['Derived']['lu_comp1'][k]]=meta['LUT']['Derived']['lu_comp1'][k]
		
		N_vis2=len(lab2)
		N_hidden2=1
		N_tot2=N_vis2+N_hidden2
		cm2=np.vstack( ((0.75,0,0,1),(0.9,0.8,0.4,1),(0.5,0.8,0.1,1),(0,0.2,0.5,1),(0.425,0,0,1),(0.2,0.5,0,1),(0.35,0.87,0.94,1),(0.85,0.95,0.7,1),(0.8,0.6,1,1),(1,0.5,0,1),(1,1,0,1),(0.4,0.7,0.79,1),(0.6,0.6,0.6,1),(1,1,1,1)) )
		cm2=matplotlib.colors.ListedColormap(cm2)
		im2=ax[i].matshow(z2,extent=z0['Extent'],cmap=cm2)			
		
	ax[i].set(xlim=z0['xlim'],ylim=z0['ylim'],xticks=xlim,yticks=ylim,xticklabels=['',''],yticklabels=['',''])
	#ax[i].yaxis.set_ticks_position('both'); ax[i].xaxis.set_ticks_position('both'); 
	ax[i].grid(meta['Graphics']['Map']['Map Grid Vis']); #ax[i].axis(meta['Graphics']['Map']['Map Axis Vis'])
	ax[i].text(Sites[Loci[i]]['x'],Sites[Loci[i]]['y'],Sites[Loci[i]]['Name'],color=[0,0,0],fontweight='bold',fontsize=10,ha='center',va='center')
ax[0].set(position=[-0.025,0.668,0.4,0.325])
ax[1].set(position=[0.43,0.668,0.4,0.325])
ax[2].set(position=[-0.025,0.337,0.4,0.325])
ax[3].set(position=[0.43,0.337,0.4,0.325])
ax[4].set(position=[-0.025,0.005,0.4,0.325])
ax[5].set(position=[0.43,0.005,0.4,0.325])

zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
cb=plt.colorbar(im1,cax=ax[6],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
cb.ax.set(yticklabels=lab)
cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
cb.outline.set_edgecolor('w')
for j in range(cb_bnd.size):
	ax[6].plot([0,100],[cb_bnd[j],cb_bnd[j]],'w-',linewidth=2)
ax[6].set(position=[0.34,0.98-(0.018*len(lab)),0.02,0.018*len(lab)])

zmn=np.min(z2); zmx=np.max(z2); cb_ivl=(zmx-zmn)/N_tot2; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden2*cb_ivl,cb_ivl)
cb_ticks=np.arange(zmn+cb_ivl/2,N_tot2-1,cb_ivl)
cb=plt.colorbar(im2,cax=ax[7],cmap=cm2,boundaries=cb_bnd,ticks=cb_ticks)
cb.ax.set(yticklabels=lab2)
cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
cb.outline.set_edgecolor('w')
for j in range(cb_bnd.size):
	ax[7].plot([0,100],[cb_bnd[j],cb_bnd[j]],'w-',linewidth=2)
ax[7].set(position=[0.795,0.98-(0.018*len(lab2)),0.02,0.018*len(lab2)])

gu.axletters(ax,plt,-0.04,0.96,FontColor=meta['Graphics']['gp']['cla'],FontSize=8,LetterStyle='NoPar',FontWeight='Bold',Skip=[6,7],StartLetterIndex=6)
gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\LandCover\LCLU_Examples2','png',900)

#%% Plot deforestation for ROIs

z={}
z['Pre10']=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Disturbances\\LandUseChange_Comp1_Early_Type.tif')
z['Pos10']=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Disturbances\\LandUseChange_Comp1_1019_Type.tif')

lut_chng=meta['LUT']['Derived']['lclu_chng_comp1']

#%% Plot LUC Examples 1
lab=['FL-CL','FL-PA','FL-RC','FL-TR','FL-EM','FL-WM','FL-ND']
lab2=['FL-CL','FL-PA','FL-RC','FL-TR','FL-EM','FL-WM','FL-ND','']
Loci=[0,0,1,1,2,2]
plt.close('all'); fig,ax=plt.subplots(7,figsize=gu.cm2inch(18,17.5))
for i in range(6):
	xlim=[Sites[Loci[i]]['x']-w/2,Sites[Loci[i]]['x']+w/2]
	ylim=[Sites[Loci[i]]['y']-w/2,Sites[Loci[i]]['y']+w/2]
	zRef0=gis.ClipRasterByXYLimits(zRef,xlim,ylim)
	
	if np.isin(i,[0,2,4,6])==True:
		z0=gis.ClipRasterByXYLimits(z['Pre10'],xlim,ylim)
		z1=len(lab)*np.ones(z0['Data'].shape,dtype='int8')
		for k in lut_chng.keys():
			ind=np.where(z0['Data']==lut_chng[k])
			if ind[0].size>0:
				z1[ind]=lut_chng[k]
		ind=np.where(z0['Data']==0)
		if ind[0].size>0:
			z1[ind]=lut_chng['FL-ND']+1
		else:
			z1[-1,-1]=lut_chng['FL-ND']+1
		ind=np.where(zRef0['Data']==0)
		if ind[0].size>0:
			z1[ind]=lut_chng['FL-ND']+2
		else:
			z1[1,3]=lut_chng['FL-ND']+2
		z1[0,0:7]=np.arange(1,8)
		
		N_vis=len(lab)
		N_hidden=2
		N_tot=N_vis+N_hidden
		cm=np.vstack( ((0.95,0.9,0.7,1),(0.75,0.9,0.55,1),(0.85,0.4,0.4,1),(0,0.3,0.7,1),(0.65,0.1,0,1),(0.4,0.6,0.8,1),(0.87,0.87,0.87,1),(0.15,0.15,0.15,1),(1,1,1,1)) )
		cm=matplotlib.colors.ListedColormap(cm)		
		im1=ax[i].matshow(z1,extent=z0['Extent'],cmap=cm)
		
	else:
		z0=gis.ClipRasterByXYLimits(z['Pos10'],xlim,ylim)
		z1=len(lab)*np.ones(z0['Data'].shape,dtype='int8')
		for k in lut_chng.keys():
			ind=np.where(z0['Data']==lut_chng[k])
			if ind[0].size>0:
				z1[ind]=lut_chng[k]
		ind=np.where(z0['Data']==0)
		if ind[0].size>0:
			z1[ind]=lut_chng['FL-ND']+1
		else:
			z1[-1,-1]=lut_chng['FL-ND']+1
		ind=np.where(zRef0['Data']==0)
		if ind[0].size>0:
			z1[ind]=lut_chng['FL-ND']+2
		else:
			z1[1,3]=lut_chng['FL-ND']+2
		z1[0,0:7]=np.arange(1,8)
		
		N_vis=len(lab)
		N_hidden=2
		N_tot=N_vis+N_hidden
		cm=np.vstack( ((0.95,0.9,0.7,1),(0.75,0.9,0.55,1),(0.85,0.4,0.4,1),(0,0.3,0.7,1),(0.65,0.1,0,1),(0.4,0.6,0.8,1),(0.87,0.87,0.87,1),(0.15,0.15,0.15,1),(1,1,1,1)) )
		cm=matplotlib.colors.ListedColormap(cm)		
		im1=ax[i].matshow(z1,extent=z0['Extent'],cmap=cm)		   
		
	ax[i].set(xlim=z0['xlim'],ylim=z0['ylim'],xticks=xlim,yticks=ylim,xticklabels=['',''],yticklabels=['',''])
	#ax[i].yaxis.set_ticks_position('both'); ax[i].xaxis.set_ticks_position('both'); 
	ax[i].grid(meta['Graphics']['Map']['Map Grid Vis']); #ax[i].axis(meta['Graphics']['Map']['Map Axis Vis'])
	ax[i].text(Sites[Loci[i]]['x'],Sites[Loci[i]]['y'],Sites[Loci[i]]['Name'],color=[0.8,0.6,1],fontweight='bold',fontsize=10,ha='center',va='center')
ax[0].set(position=[-0.025,0.668,0.4,0.325])
ax[1].set(position=[0.43,0.668,0.4,0.325])
ax[2].set(position=[-0.025,0.337,0.4,0.325])
ax[3].set(position=[0.43,0.337,0.4,0.325])
ax[4].set(position=[-0.025,0.005,0.4,0.325])
ax[5].set(position=[0.43,0.005,0.4,0.325])

zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
cb=plt.colorbar(im1,cax=ax[6],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
cb.ax.set(yticklabels=lab2)
cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
cb.outline.set_edgecolor('w')
for j in range(cb_bnd.size):
	ax[6].plot([0,100],[cb_bnd[j],cb_bnd[j]],'w-',linewidth=2)
ax[6].set(position=[0.34,0.98-(0.018*len(lab2)),0.02,0.018*len(lab2)])

gu.axletters(ax,plt,-0.04,0.96,FontColor=meta['Graphics']['gp']['cla'],FontSize=8,LetterStyle='NoPar',FontWeight='Bold',Skip=[6])
gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\LandCover\LUC_Examples1','png',900)

#%% Plot LUC Examples 2 
lab=['FL-CL','FL-PA','FL-RC','FL-TR','FL-EM','FL-WM','FL-ND']
lab2=['FL-CL','FL-PA','FL-RC','FL-TR','FL-EM','FL-WM','FL-ND','']
Loci=[3,3,4,4,5,5]
plt.close('all'); fig,ax=plt.subplots(7,figsize=gu.cm2inch(18,17.5))
for i in range(6):
	xlim=[Sites[Loci[i]]['x']-w/2,Sites[Loci[i]]['x']+w/2]
	ylim=[Sites[Loci[i]]['y']-w/2,Sites[Loci[i]]['y']+w/2]
	zRef0=gis.ClipRasterByXYLimits(zRef,xlim,ylim)
	
	if np.isin(i,[0,2,4,6])==True:
		z0=gis.ClipRasterByXYLimits(z['Pre10'],xlim,ylim)
		z1=len(lab)*np.ones(z0['Data'].shape,dtype='int8')
		for k in lut_chng.keys():
			ind=np.where(z0['Data']==lut_chng[k])
			if ind[0].size>0:
				z1[ind]=lut_chng[k]
		ind=np.where(z0['Data']==0)
		if ind[0].size>0:
			z1[ind]=lut_chng['FL-ND']+1
		else:
			z1[-1,-1]=lut_chng['FL-ND']+1
		ind=np.where(zRef0['Data']==0)
		if ind[0].size>0:
			z1[ind]=lut_chng['FL-ND']+2
		else:
			z1[1,3]=lut_chng['FL-ND']+2
		z1[0,0:7]=np.arange(1,8)
		
		N_vis=len(lab)
		N_hidden=2
		N_tot=N_vis+N_hidden
		cm=np.vstack( ((0.95,0.9,0.7,1),(0.75,0.9,0.55,1),(0.85,0.4,0.4,1),(0,0.3,0.7,1),(0.65,0.1,0,1),(0.4,0.6,0.8,1),(0.87,0.87,0.87,1),(0.15,0.15,0.15,1),(1,1,1,1)) )
		cm=matplotlib.colors.ListedColormap(cm)		
		im1=ax[i].matshow(z1,extent=z0['Extent'],cmap=cm)
		
	else:
		z0=gis.ClipRasterByXYLimits(z['Pos10'],xlim,ylim)
		z1=len(lab)*np.ones(z0['Data'].shape,dtype='int8')
		for k in lut_chng.keys():
			ind=np.where(z0['Data']==lut_chng[k])
			if ind[0].size>0:
				z1[ind]=lut_chng[k]
		ind=np.where(z0['Data']==0)
		if ind[0].size>0:
			z1[ind]=lut_chng['FL-ND']+1
		else:
			z1[-1,-1]=lut_chng['FL-ND']+1
		ind=np.where(zRef0['Data']==0)
		if ind[0].size>0:
			z1[ind]=lut_chng['FL-ND']+2
		else:
			z1[1,3]=lut_chng['FL-ND']+2
		z1[0,0:7]=np.arange(1,8)
		
		N_vis=len(lab)
		N_hidden=2
		N_tot=N_vis+N_hidden
		cm=np.vstack( ((0.95,0.9,0.7,1),(0.75,0.9,0.55,1),(0.85,0.4,0.4,1),(0,0.3,0.7,1),(0.65,0.1,0,1),(0.4,0.6,0.8,1),(0.87,0.87,0.87,1),(0.15,0.15,0.15,1),(1,1,1,1)) )
		cm=matplotlib.colors.ListedColormap(cm)		
		im1=ax[i].matshow(z1,extent=z0['Extent'],cmap=cm)		   
		
	ax[i].set(xlim=z0['xlim'],ylim=z0['ylim'],xticks=xlim,yticks=ylim,xticklabels=['',''],yticklabels=['',''])
	#ax[i].yaxis.set_ticks_position('both'); ax[i].xaxis.set_ticks_position('both'); 
	ax[i].grid(meta['Graphics']['Map']['Map Grid Vis']); #ax[i].axis(meta['Graphics']['Map']['Map Axis Vis'])
	ax[i].text(Sites[Loci[i]]['x'],Sites[Loci[i]]['y'],Sites[Loci[i]]['Name'],color=[0.8,0.6,1],fontweight='bold',fontsize=10,ha='center',va='center')
ax[0].set(position=[-0.025,0.668,0.4,0.325])
ax[1].set(position=[0.43,0.668,0.4,0.325])
ax[2].set(position=[-0.025,0.337,0.4,0.325])
ax[3].set(position=[0.43,0.337,0.4,0.325])
ax[4].set(position=[-0.025,0.005,0.4,0.325])
ax[5].set(position=[0.43,0.005,0.4,0.325])

lab.append('')
zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
cb=plt.colorbar(im1,cax=ax[6],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
cb.ax.set(yticklabels=lab)
cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
cb.outline.set_edgecolor('w')
for j in range(cb_bnd.size):
	ax[6].plot([0,100],[cb_bnd[j],cb_bnd[j]],'w-',linewidth=2)
ax[6].set(position=[0.34,0.98-(0.018*len(lab)),0.02,0.018*len(lab)])

gu.axletters(ax,plt,-0.04,0.96,FontColor=meta['Graphics']['gp']['cla'],FontSize=8,LetterStyle='NoPar',FontWeight='Bold',Skip=[6],StartLetterIndex=6)
gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\LandCover\LUC_Examples2','png',900)

