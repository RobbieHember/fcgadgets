#%% Import modules
import os
import numpy as np
import gc
import matplotlib.pyplot as plt
import matplotlib.colors
from matplotlib.colors import LightSource
import geopandas as gpd
import pandas as pd
from shapely.geometry import Polygon,Point,box,shape
from rasterio import features
import copy
import time
import cv2
import fcgadgets.macgyver.util_general as gu
import fcgadgets.macgyver.util_gis as gis
import fcgadgets.bc1ha.bc1ha_utils as u1ha
import fcgadgets.macgyver.util_query_gdb as qgdb
import fcexplore.field_plots.fp_util as ufp
import fcgadgets.cbrunner.cbrun_util as cbu

#%%
def PlotVectorBaseMaps(meta,roi,ax):
	#ec_bound=[0,0,0]
	ec_bound=[0.5,0.5,0.5]

	ec_road=[0.4,0,0]

	#ec_cities=[0.75,0.3,0] # Orange
	#fc_cities=[1,0.6,0]
	ec_cities=[0,0,0] # Black
	fc_cities=[0,0,0] # Black

	#ms_cities=10
	ms_cities=5
	fw_cities='bold'

	#symb_cities='o'
	symb_cities='s'

	if meta['Graphics']['Map']['Show Bound Within']=='On':
		roi['gdf']['bound ROI'].plot(ax=ax,edgecolor=ec_bound,facecolor='none',linewidth=0.25)
	if meta['Graphics']['Map']['Show Bound Land Mask']=='On':
		roi['gdf']['bc_bound'].plot(ax=ax,edgecolor=ec_bound,facecolor='none',linewidth=0.25)
	if meta['Graphics']['Map']['Show Lakes']=='On':
		roi['gdf']['lakes'].plot(ax=ax,facecolor=[0.5,0.8,1],label='Lakes',linewidth=0.25)
	if meta['Graphics']['Map']['Show Rivers']=='On':
		roi['gdf']['rivers'].plot(ax=ax,color=[0.5,0.8,1],label='Rivers',linewidth=0.4)
		roi['gdf']['riversecond'].plot(ax=ax,facecolor=[0.5,0.8,1],edgecolor=[0.5,0.8,1],label='Rivers',linewidth=0.6)
		roi['gdf']['rivermajor'].plot(ax=ax,facecolor=[0.5,0.8,1],edgecolor=[0.5,0.8,1],label='Rivers',linewidth=1)
	if meta['Graphics']['Map']['Show Roads']=='On':
		roi['gdf']['road'].plot(ax=ax,edgecolor=ec_road,linewidth=0.25,label='Road',alpha=1,zorder=1)
	if meta['Graphics']['Map']['Show Rail']=='On':
		roi['gdf']['rail'].plot(ax=ax,edgecolor=[0,0,0],linewidth=0.5,label='Rail',alpha=1,zorder=1)
		roi['gdf']['rail'].plot(ax=ax,edgecolor=[1,1,1],linewidth=0.25,label='Rail',alpha=1,zorder=1)
	if meta['Graphics']['Map']['Show TPFs']=='On':
		roi['gdf']['tpf'].plot(ax=ax,marker='s',edgecolor=[0.4,0,1],linewidth=0.25,label='TPF',alpha=1,zorder=1,markersize=15)
		if meta['Graphics']['Map']['Show Symbol Labels']=='On':
			for x,y,label in zip(roi['gdf']['tpf'].geometry.x,roi['gdf']['tpf'].geometry.y,roi['gdf']['tpf']['COMPANY_NAME']):
				ax.annotate(label,xy=(x,y),xytext=(4,3),textcoords="offset points",color=[0,0,0],fontsize=6)
	if meta['Graphics']['Map']['Show Cities']=='On':
		ind=np.where( (roi['gdf']['cities']['Level']==1) & (roi['gdf']['cities']['Territory']=='BC') )[0]
		roi['gdf']['cities'].iloc[ind].plot(ax=ax,marker=symb_cities,edgecolor=ec_cities,facecolor=fc_cities,lw=1,markersize=ms_cities,alpha=1,zorder=2)
		if meta['Graphics']['Map']['Show Symbol Labels']=='On':
			for x,y,label in zip(roi['gdf']['cities'].iloc[ind].geometry.x,roi['gdf']['cities'].iloc[ind].geometry.y,roi['gdf']['cities'].iloc[ind]['City Name']):
				ax.annotate(label,xy=(x,y),xytext=(3,2),textcoords="offset points",color=[0,0,0],fontsize=5,fontweight=fw_cities)
	return roi

#%%
def Plot_LandCoverComp1(meta,roi,nam):
	z0=roi['grd'][nam]['Data']
	lab=list(meta['LUT']['Derived']['lc_comp1'].keys())

	z1=len(lab)*np.ones(z0.shape,dtype='int8')
	for k in meta['LUT']['Derived']['lc_comp1'].keys():
		ind=np.where(z0==meta['LUT']['Derived']['lc_comp1'][k])
		if ind[0].size>0:
			z1[ind]=meta['LUT']['Derived']['lc_comp1'][k]
		else:
			z1[0,meta['LUT']['Derived']['lc_comp1'][k]]=meta['LUT']['Derived']['lc_comp1'][k]
	ind=np.where(roi['grd']['Data']==0)
	if ind[0].size>0:
		z1[ind]=meta['LUT']['Derived']['lc_comp1'][k]+1
	else:
		z1[-1,-1]=meta['LUT']['Derived']['lc_comp1'][k]+1

	N_vis=len(lab)
	N_hidden=1
	N_tot=N_vis+N_hidden

	cm=np.vstack( ((0.3,0.6,0,1),(0.65,1,0,1),(1,1,0.5,1),(0.8,0.7,1,1),(0.75,0.5,0.5,1),(0.75,0.75,0.75,1),(0.8,0,0,1),(0.95,0.95,0.95,1),(0.75,0.85,0.95,1),(1,1,1,1)) )
	cm=matplotlib.colors.ListedColormap(cm)

	plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
	im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
	roi=PlotVectorBaseMaps(meta,roi,ax[0])
	ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

	zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
	cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
	cb=plt.colorbar(im,cax=ax[1],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
	#ax[1].set(position=[0.71,0.6,0.05,0.14])
	cb.ax.set(yticklabels=lab)
	cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
	cb.outline.set_edgecolor('w')
	for i in range(cb_bnd.size):
		ax[1].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)
	pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
	ax[1].set(position=pos2)
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_' + nam,'png',900)
	
	return fig,ax

#%%
def Plot_Peatland(meta,roi,nam):
	z0=roi['grd'][nam]['Data']
	lab=list(meta['LUT']['Raw']['Peat']['Name'])

	z1=len(lab)*np.ones(z0.shape,dtype='int8')
	for i in range(meta['LUT']['Raw']['Peat']['ID'].size):
		ind=np.where(z0==meta['LUT']['Raw']['Peat']['ID'][i])
		if ind[0].size>0:
			z1[ind]=i
		else:
			z1[0,i]=i
	ind=np.where( (roi['grd']['lc_comp1_2019']['Data']==1) & (z0==0) ); z1[ind]=i+1
	ind=np.where( (roi['grd']['lc_comp1_2019']['Data']!=1) & (z0==0) ); z1[ind]=i+2
	ind=np.where(roi['grd']['Data']==0); z1[ind]=i+3

	lab.append('Upland forest')
	lab.append('Non-forest')
	N_vis=len(lab)
	N_hidden=1
	N_tot=N_vis+N_hidden

	cm=np.vstack( (np.column_stack((meta['LUT']['Raw']['Peat']['R'],meta['LUT']['Raw']['Peat']['G'],meta['LUT']['Raw']['Peat']['B'],np.ones(meta['LUT']['Raw']['Peat']['B'].size))),(0.8,0.8,0.8,1),(0.95,0.95,0.95,1),(1,1,1,1)) )
	cm=matplotlib.colors.ListedColormap(cm)

	plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
	im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
	roi=PlotVectorBaseMaps(meta,roi,ax[0])
	ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

	zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
	cb_ticks=np.arange(zmn+cb_ivl/2,N_vis-1,cb_ivl)
	cb=plt.colorbar(im,cax=ax[1],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
	#ax[1].set(position=[0.71,0.6,0.05,0.14])
	cb.ax.set(yticklabels=lab)
	cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
	cb.outline.set_edgecolor('w')
	for i in range(cb_bnd.size):
		ax[1].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)
	pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
	ax[1].set(position=pos2)
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_' + nam,'png',900)
	
	return fig,ax

#%%
def Plot_gromo_GN(meta,roi,nam):
	z0=roi['grd'][nam]['Data']

	bw=0.5; bin=np.arange(-1,6+bw,bw)
	N_vis=bin.size
	N_hidden=2
	N_tot=N_vis+N_hidden
	z1=N_vis*np.ones(z0.shape)
	for i in range(bin.size):
		ind=np.where(np.abs(z0-bin[i])<=bw/2)
		z1[ind]=i
	z1[(z0<=bin[0])]=1
	z1[(z0>=bin[i])]=i
	z1[(roi['grd']['lc_comp1_2019']['Data']!=meta['LUT']['Derived']['lc_comp1']['Forest'])]=i+1
	z1[(roi['grd']['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Water']) | (roi['grd']['lc_comp1_2019']['Data']==0)]=i+2
	z1[(roi['grd']['Data']==0)]=i+2
	for i in range(bin.size):
		z1[0,i]=i
	lab=np.round(bin,decimals=2).astype(str)
	lab=np.append(lab,'')
	lab=np.append(lab,'')

	cm=plt.cm.get_cmap('RdYlGn',N_vis)
	cm.colors=cm(np.arange(0,cm.N))
	cm=np.vstack( ((50/255,30/255,0/255,1),(125/255,85/255,0/255,1),(175/255,155/255,105/255,1),cm.colors[6:,:],(0/255,60/255,0/255,1),(5/255,40/255,5/255,1),(10/255,20/255,10/255,1),(0.85,0.85,0.85,1),(1,1,1,1)) )
	cm=matplotlib.colors.ListedColormap(cm)

	plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
	im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
	roi=PlotVectorBaseMaps(meta,roi,ax[0])
	ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

	zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
	cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
	cb=plt.colorbar(im,cax=ax[1],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
	#ax[1].set(position=[0.71,0.6,0.05,0.14])
	cb.ax.set(yticklabels=lab)
	cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
	cb.outline.set_edgecolor('w')
	for i in range(cb_bnd.size):
		ax[1].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)
	pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
	ax[1].set(position=pos2)
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_' + nam,'png',900)
	return fig,ax

#%%
def Plot_gromo_G(meta,roi,nam):
	z0=roi['grd'][nam]['Data']

	bw=0.5; bin=np.arange(0,7+bw,bw)
	N_vis=bin.size
	N_hidden=2
	N_tot=N_vis+N_hidden
	z1=N_vis*np.ones(z0.shape)
	for i in range(bin.size):
		ind=np.where(np.abs(z0-bin[i])<=bw/2)
		z1[ind]=i
	z1[(z0<=bin[0])]=1
	z1[(z0>=bin[i])]=i
	z1[(roi['grd']['lc_comp1_2019']['Data']!=meta['LUT']['Derived']['lc_comp1']['Forest'])]=i+1
	z1[(roi['grd']['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Water']) | (roi['grd']['lc_comp1_2019']['Data']==0)]=i+2
	z1[(roi['grd']['Data']==0)]=i+2
	for i in range(bin.size):
		z1[0,i]=i
	lab=np.round(bin,decimals=2).astype(str)
	lab=np.append(lab,'')
	lab=np.append(lab,'')

	cm=plt.cm.get_cmap('YlGn',N_vis)
	cm.colors=cm(np.arange(0,cm.N))
	cm1=np.vstack([cm.colors[0,:],cm.colors[2,:],cm.colors[4,:],cm.colors[6:,:]])
	cm=np.vstack( ((100/255,100/255,25/255,1),(170/255,165/255,140/255,1),(225/255,225/255,175/255,1),cm1,(0.85,0.85,0.85,1),(1,1,1,1)) )
	#cm=np.vstack( (cm.colors,(0.85,0.85,0.85,1),(1,1,1,1)) )
	cm=matplotlib.colors.ListedColormap(cm)

	plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
	im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
	roi=PlotVectorBaseMaps(meta,roi,ax[0])
	ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

	zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
	cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
	cb=plt.colorbar(im,cax=ax[1],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
	cb.ax.set(yticklabels=lab)
	cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
	cb.outline.set_edgecolor('w')
	for i in range(cb_bnd.size):
		ax[1].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)
	pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
	ax[1].set(position=pos2)
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_' + nam,'png',900)
	return fig,ax

#%%
def Plot_gromo_M(meta,roi,nam):
	z0=roi['grd'][nam]['Data']

	bw=0.5; bin=np.arange(0,7+bw,bw)
	N_vis=bin.size
	N_hidden=2
	N_tot=N_vis+N_hidden
	z1=N_vis*np.ones(z0.shape)
	for i in range(bin.size):
		ind=np.where(np.abs(z0-bin[i])<=bw/2)
		z1[ind]=i
	z1[(z0<=bin[0])]=1
	z1[(z0>=bin[i])]=i
	z1[(roi['grd']['lc_comp1_2019']['Data']!=meta['LUT']['Derived']['lc_comp1']['Forest'])]=i+1
	z1[(roi['grd']['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Water']) | (roi['grd']['lc_comp1_2019']['Data']==0)]=i+2
	z1[(roi['grd']['Data']==0)]=i+2
	for i in range(bin.size):
		z1[0,i]=i
	lab=np.round(bin,decimals=2).astype(str)
	lab=np.append(lab,'')
	lab=np.append(lab,'')

	cm=plt.cm.get_cmap('YlOrBr',N_vis)
	cm.colors=cm(np.arange(0,cm.N))
	cm=np.vstack( (cm.colors,(0.85,0.85,0.85,1),(1,1,1,1)) )
	cm=matplotlib.colors.ListedColormap(cm)

	plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
	im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
	roi=PlotVectorBaseMaps(meta,roi,ax[0])
	ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

	zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
	cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
	cb=plt.colorbar(im,cax=ax[1],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
	cb.ax.set(yticklabels=lab)
	cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
	cb.outline.set_edgecolor('w')
	for i in range(cb_bnd.size):
		ax[1].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)
	pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
	ax[1].set(position=pos2)
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_' + nam,'png',900)
	return fig,ax

#%%
def Plot_VectorVariable(meta,roi):
	gdf=gpd.read_file(meta['Paths']['BCFCS_NMC']['Data'] + '\\geos.geojson')
	fig,ax=plt.subplots(1,figsize=gu.cm2inch(12,12))
	#roi=PlotVectorBaseMaps(meta,roi,ax)	
	roi['gdf']['bc_bound'].plot(ax=ax,facecolor=[0.8,0.8,0.8],edgecolor=[0,0,0],label='Political Boundary',linewidth=0.25,alpha=1)
	gdf.plot(ax=ax,markersize=2,facecolor=[0.75,0,0],edgecolor=None,linewidth=0.75,alpha=1)
	#ax.grid(color='k',linestyle='-',linewidth=0.25)	
	#ax.set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
	#ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])
	# ax.grid(color='k',linestyle='-',linewidth=0.25)
	ax.set(position=[0.01,0.01,0.98,0.98],xticks=[],yticks=[])
	ax.axis(meta['Graphics']['Map']['Map Axis Vis'])
	#if meta['Graphics']['Print Figures']=='On':
	#	gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_' + nam,'png',900)
	
	return fig,ax

#%%
def Plot_ClimateNormalsPanels(meta,roi):

	plt.close('all'); fig,ax=plt.subplots(16,figsize=gu.cm2inch(24,12))
	vL=[   'tmean_ann_n','prcp_ann_n','etp_ann_n','runoff_ann_n','melt_ann_n','wsp_ann_n','ws_ann_n','cwd_ann_n']
	dim0=[ 1.0,           30.5,        30.5,       30.5,          30.5,        1.0,        1.0,       30.5]
	bw0=[  1.0,            0.5,         0.2,        0.5,           0.2,       50.0,       20.0,        0.2]
	ymin=[-5.0,            0.0,         0.0,        0.0,           0.0,        0,          0,          0.0]
	ymax=[10.0,            6.0,         3.0,        6.0,           3.0,      400,        200,          3.0]
	for i in range(len(vL)):
		v=vL[i]
		bw=bw0[i]
		bin=np.arange(ymin[i],ymax[i]+bw,bw)
		z1=bin.size*np.ones( roi['grd'][v]['Data'].shape)
		for j in range(bin.size):
			ind=np.where(np.abs(roi['grd'][v]['Data']/dim0[i]-bin[j])<=bw/2)
			if ind[0].size>0:
				z1[ind]=j
			else:
				z1[0,i]=j

		#plt.close('all'); plt.matshow(z1)
		z1[(roi['grd'][v]['Data']/dim0[i]<ymin[i])]=0
		z1[(roi['grd'][v]['Data']/dim0[i]>ymax[i])]=j
		z1[(roi['grd']['lc_comp1_2019']==0) | (roi['grd']['lc_comp1_2019']==meta['LUT']['Derived']['lc_comp1']['Water']) | (roi['grd']['Data']!=1)]=j+1
		lab=np.round(bin,decimals=1).astype(str)
		lab=np.append(lab,'Hidden')
		N_vis=bin.size
		N_hidden=1
		N_tot=N_vis+N_hidden
		cm=plt.cm.get_cmap('plasma',bin.size)
		cm=np.vstack( (cm.colors,(1,1,1,1)) )
		cm=matplotlib.colors.ListedColormap(cm)
		im=ax[i].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
		#roi=PlotVectorBaseMaps(meta,roi,ax[0])
		ax[i].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
		ax[i].yaxis.set_ticks_position('both'); ax[i].xaxis.set_ticks_position('both'); ax[i].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[i].axis(meta['Graphics']['Map']['Map Axis Vis'])

		zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
		cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
		cb=plt.colorbar(im,cax=ax[i+8],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
		ax[i+8].set(position=[0.71,0.6,0.05,0.14])
		cb.ax.set(yticklabels=lab)
		cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
		cb.outline.set_edgecolor('w')
		for j in range(cb_bnd.size):
			ax[i+8].plot([0,100],[cb_bnd[j],cb_bnd[j]],'w-',linewidth=2)

	ax[0].set(position=[0,0.5,0.25,0.5])
	ax[1].set(position=[0.25,0.5,0.25,0.5])
	ax[2].set(position=[0.5,0.5,0.25,0.5])
	ax[3].set(position=[0.75,0.5,0.25,0.5])
	ax[4].set(position=[0,0,0.25,0.5])
	ax[5].set(position=[0.25,0,0.25,0.5])
	ax[6].set(position=[0.5,0,0.25,0.5])
	ax[7].set(position=[0.75,0,0.25,0.5])

	for i in range(8):
		p=[ax[i].get_position().xmin+0.03,ax[i].get_position().ymin-0.03,0.22,0.5]
		ax[i].set(position=p)

	ah=0.38
	ax[8].set(position=[0.0,0.53,0.012,ah])
	ax[9].set(position=[0.25,0.53,0.012,ah])
	ax[10].set(position=[0.5,0.53,0.012,ah])
	ax[11].set(position=[0.75,0.53,0.012,ah])
	ax[12].set(position=[0.0,0.03,0.012,ah])
	ax[13].set(position=[0.25,0.03,0.012,ah])
	ax[14].set(position=[0.5,0.03,0.012,ah])
	ax[15].set(position=[0.75,0.03,0.012,ah])

	vL2=['Temperature (degC)','Precipiration (mm/d)','Potential evapotranspiration (mm/d)','Runnoff (mm/d)','Snowmelt (mm/d)','Snowpack water (mm)','Soil water (mm)','Climatic water deficit (mm/d)']
	gu.axletters(ax,plt,-0.025,1.01,FontColor=meta['Graphics']['gp']['cla'],LetterStyle='NoPar',FontWeight='Bold',Skip=[8,9,10,11,12,13,14,15],Labels=vL2)

	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_ClimateNormalPanels','png',900)
	return fig,ax

#%%
def Plot_LandUseComp1(meta,roi,nam):
	
	z0=roi['grd'][nam]['Data']
	lab=list(meta['LUT']['Derived']['lu_comp1'].keys())
	z1=len(lab)*np.ones(z0.shape,dtype='int8')
	for k in meta['LUT']['Derived']['lu_comp1'].keys():
		ind=np.where(z0==meta['LUT']['Derived']['lu_comp1'][k])
		z1[ind]=meta['LUT']['Derived']['lu_comp1'][k]
	ind=np.where(roi['grd']['Data']==0)
	if ind[0].size>0:
		z1[ind]=meta['LUT']['Derived']['lu_comp1'][k]+1
	else:
		z1[-1,-1]=meta['LUT']['Derived']['lu_comp1'][k]+1
	for k in meta['LUT']['Derived']['lu_comp1'].keys():
		z1[0,meta['LUT']['Derived']['lu_comp1'][k]]=meta['LUT']['Derived']['lu_comp1'][k]
	
	N_vis=len(lab)
	N_hidden=1
	N_tot=N_vis+N_hidden
	cm=np.vstack( ((0.65,0,0,1),(1,0.95,0.5,1),(0.5,0.8,0.1,1),(0,0.2,0.5,1),(0.25,0,0,1),(0.2,0.5,0,1),
				   (0.35,0.87,0.94,1),(0.85,0.95,0.7,1),(0.8,0.6,1,1),(1,0.75,0.2,1),(1,0.4,0.15,1),(0.4,0.7,0.79,1),(0.6,0.6,0.6,1),(1,1,1,1)) )
	cm=matplotlib.colors.ListedColormap(cm)

	plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
	im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
	roi=PlotVectorBaseMaps(meta,roi,ax[0])
	ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

	zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
	cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
	cb=plt.colorbar(im,cax=ax[1],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
	ax[1].set(position=[0.71,0.6,0.05,0.14])
	cb.ax.set(yticklabels=lab)
	cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
	cb.outline.set_edgecolor('w')
	for i in range(cb_bnd.size):
		ax[1].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)
	pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
	ax[1].set(position=pos2)
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_' + nam,'png',900)
	return fig,ax

#%%
def Plot_LandUseComp1Panels(meta,roi):
	namL=['lu_comp1_2019','lu_comp1_2049s1','lu_comp1_2049s2','lu_comp1_2049s3']
	cnt=0
	plt.close('all'); 
	fig,ax=plt.subplots(5,figsize=gu.cm2inch(2*meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*(2*meta['Graphics']['Map']['Fig Width'])*roi['grd']['yxrat']))
	for nam in namL:
		z0=roi['grd'][nam]['Data']
		lab=list(meta['LUT']['Derived']['lu_comp1'].keys())
		z1=len(lab)*np.ones(z0.shape,dtype='int8')
		for k in meta['LUT']['Derived']['lu_comp1'].keys():
			ind=np.where(z0==meta['LUT']['Derived']['lu_comp1'][k])
			z1[ind]=meta['LUT']['Derived']['lu_comp1'][k]
		ind=np.where(roi['grd']['Data']==0)
		if ind[0].size>0:
			z1[ind]=meta['LUT']['Derived']['lu_comp1'][k]+1
		else:
			z1[-1,-1]=meta['LUT']['Derived']['lu_comp1'][k]+1
		for k in meta['LUT']['Derived']['lu_comp1'].keys():
			z1[0,meta['LUT']['Derived']['lu_comp1'][k]]=meta['LUT']['Derived']['lu_comp1'][k]
		
		N_vis=len(lab)
		N_hidden=1
		N_tot=N_vis+N_hidden
		cm=np.vstack( ((0.65,0,0,1),(0.9,0.8,0.4,1),(0.5,0.8,0.1,1),(0,0.2,0.5,1),(0.25,0,0,1),(0.2,0.5,0,1),(0.35,0.87,0.94,1),(0.85,0.95,0.7,1),(0.8,0.6,1,1),(1,1,0,1),(1,0,0,1),(0.4,0.7,0.79,1),(0.6,0.6,0.6,1),(1,1,1,1)) )
		cm=matplotlib.colors.ListedColormap(cm)	
		im=ax[cnt].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
		roi=PlotVectorBaseMaps(meta,roi,ax[cnt])
		ax[cnt].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
		ax[cnt].yaxis.set_ticks_position('both'); ax[cnt].xaxis.set_ticks_position('both'); ax[cnt].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[cnt].axis(meta['Graphics']['Map']['Map Axis Vis'])
		cnt=cnt+1
	gu.axletters(ax,plt,-0.025,0.9,FontColor=meta['Graphics']['gp']['cla'],LetterStyle='NoPar',FontWeight='Bold',Skip=[4])

	ax[0].set(position=[0,0.5,0.5,0.5])
	ax[1].set(position=[0.5,0.5,0.5,0.5])
	ax[2].set(position=[0,0,0.5,0.5])
	ax[3].set(position=[0.5,0,0.5,0.5])

	zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
	cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
	cb=plt.colorbar(im,cax=ax[4],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
	cb.ax.set(yticklabels=lab)
	cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
	cb.outline.set_edgecolor('w')
	for i in range(cb_bnd.size):
		ax[4].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)
	#pos2=[0.78,0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
	pos2=[0.86,0.74,0.02,0.25]
	ax[4].set(position=pos2)
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_LandUseComp1Panels','png',900)
	return fig,ax

#%%
def Plot_LandCoverLandUseChange(meta,roi,nam):
	#nam='luc1_early_type'
	z0=roi['grd'][nam]['Data']
	lab=list(meta['LUT']['Derived']['lclu_chng_comp1'].keys())

	z1=(len(lab)+1)*np.ones(z0.shape,dtype='int8')
	for k in meta['LUT']['Derived']['lc_comp1'].keys():
		ind=np.where(z0==meta['LUT']['Derived']['lc_comp1'][k])
		if ind[0].size>0:
			z1[ind]=meta['LUT']['Derived']['lc_comp1'][k]
		else:
			z1[0,meta['LUT']['Derived']['lc_comp1'][k]]=meta['LUT']['Derived']['lc_comp1'][k]
	ind=np.where(roi['grd']['Data']==0)
	if ind[0].size>0:
		z1[ind]=meta['LUT']['Derived']['lc_comp1'][k]+1
	else:
		z1[-1,-1]=meta['LUT']['Derived']['lc_comp1'][k]+1
	#z0[0,0:7]=np.arange(1,8)
	z1[(roi['grd']['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Water'])]=(len(lab)+2)
	
	N_vis=len(lab)
	N_hidden=2
	N_tot=N_vis+N_hidden

	cm=np.vstack( ((0.95,0.9,0.7,1),(0.75,0.9,0.55,1),(0.85,0.4,0.4,1),(0.75,0.65,1,1),(0.65,0.1,0,1),(0.45,0.75,1,1),(0.87,0.87,0.87,1),(0.15,0.15,0.15,1),(1,1,1,1)) )
	cm=matplotlib.colors.ListedColormap(cm)

	plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
	im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
	roi=PlotVectorBaseMaps(meta,roi,ax[0])
	ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

	zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
	cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
	cb=plt.colorbar(im,cax=ax[1],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
	ax[1].set(position=[0.71,0.6,0.05,0.14])
	cb.ax.set(yticklabels=lab)
	cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
	cb.outline.set_edgecolor('w')
	for i in range(cb_bnd.size):
		ax[1].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)
	pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
	ax[1].set(position=pos2)
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_' + nam,'png',900)
	
	return fig,ax

#%%
def Plot_Geomorphons(meta,roi,vnam):

	lab=list(meta['LUT']['Derived']['geomorph'].keys())
	z1=roi['grd'][vnam]['Data']
	ind=np.where(z1<=0); z1[ind]=len(lab)+1
	ind=np.where(roi['grd']['Data']==0)
	if ind[0].size>0:
		z1[ind]=len(lab)+1
	ind=np.where(roi['grd']['lc_comp1_2019']['Data']==0); z1[ind]=11
	ind=np.where(roi['grd']['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Water']); z1[ind]=11
	z1[0:10,0]=np.arange(1,11)

	N_vis=len(lab)
	N_hidden=1
	N_tot=N_vis+N_hidden
	cm=np.vstack( ((0.9,0.9,0.9,1),(0.25,0,0,1),
				   (0.75,0,0,1),(1,0.5,0.5,1),
				   (1,1,0.5,1),
				   (0.65,0.87,0.24,1),(0.85,0.95,0.7,1),(0.65,0.25,1,1),
				   (0.8,0.9,1,1),(0,0.25,0.5,1),(1,1,1,1)) )
	cm=matplotlib.colors.ListedColormap(cm)

	plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
	im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
	roi=PlotVectorBaseMaps(meta,roi,ax[0])
	ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

	zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
	cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
	cb=plt.colorbar(im,cax=ax[1],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
	ax[1].set(position=[0.71,0.6,0.05,0.14])
	cb.ax.set(yticklabels=lab)
	cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
	cb.outline.set_edgecolor('w')
	for i in range(cb_bnd.size):
		ax[1].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)
	pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
	ax[1].set(position=pos2)
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_' + vnam,'png',900)
	return fig,ax

#%%
def Plot_UplandWetlandForest(meta,roi,nam):

	z0=roi['grd'][nam]['Data']
	lab=list(meta['LUT']['Derived']['upwetf'].keys())
	z1=len(lab)*np.ones(z0.shape,dtype='int8')
	for k in meta['LUT']['Derived']['upwetf'].keys():
		ind=np.where(z0==meta['LUT']['Derived']['upwetf'][k])
		if ind[0].size>0:
			z1[ind]=meta['LUT']['Derived']['upwetf'][k]
		else:
			z1[0,meta['LUT']['Derived']['upwetf'][k]]=meta['LUT']['Derived']['upwetf'][k]
	ind=np.where(roi['grd']['Data']==0)
	if ind[0].size>0:
		z1[ind]=meta['LUT']['Derived']['upwetf'][k]+1
	else:
		z1[-1,-1]=meta['LUT']['Derived']['upwetf'][k]+1

	N_vis=len(lab)
	N_hidden=1
	N_tot=N_vis+N_hidden

	cm=np.vstack( ((0.3,0.6,0,1),(0.75,0.95,1,1),(0.2,0.65,1,1),(0.9,0.9,0.9,1),(1,1,1,1)) )
	cm=matplotlib.colors.ListedColormap(cm)

	plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
	im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
	roi=PlotVectorBaseMaps(meta,roi,ax[0])
	ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

	zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
	cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
	cb=plt.colorbar(im,cax=ax[1],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
	ax[1].set(position=[0.71,0.6,0.05,0.14])
	cb.ax.set(yticklabels=lab)
	cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
	cb.outline.set_edgecolor('w')
	for i in range(cb_bnd.size):
		ax[1].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)
	pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
	ax[1].set(position=pos2)
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_' + nam,'png',900)
	plt.show()
	return fig,ax

#%%
def Plot_AccessZones(meta,roi,nam):

	meta['Graphics']['Map']['Show Rail']='On'
	meta['Graphics']['Map']['Show Cities']='On'
	meta['Graphics']['Map']['Show Symbol Labels']='On'

	z1=roi['grd'][nam]['Data']
	lab=list(meta['LUT']['Derived'][nam].keys())
# 	z1=len(lab)*np.ones(z0.shape,dtype='int8')
# 	for k in meta['LUT']['Derived']['access'].keys():
# 		ind=np.where(z0==meta['LUT']['Derived']['access'][k])
# 		if ind[0].size>0:
# 			z1[ind]=meta['LUT']['Derived']['access'][k]
# 		else:
# 			z1[0,meta['LUT']['Derived']['access'][k]]=meta['LUT']['Derived']['access'][k]

	ind=np.where( (roi['grd']['Data']==0) | (z1==0) )
	z1[ind]=meta['LUT']['Derived']['access']['Non-forest']
	z1[0,0:5]=np.arange(1,6,1)

	N_vis=len(lab)-1
	N_hidden=1
	N_tot=N_vis+N_hidden
	lab=lab[0:-1]
	cm=np.vstack( ((0.7,0.75,1,1),(0.75,0.6,1,1),(0.4,0.8,0,1),(0.95,0.95,0.95,1),(1,1,1,1)) )
	cm=matplotlib.colors.ListedColormap(cm)

	plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
	im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
	roi=PlotVectorBaseMaps(meta,roi,ax[0])
	ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

	zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
	cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
	cb=plt.colorbar(im,cax=ax[1],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
	ax[1].set(position=[0.71,0.6,0.05,0.14])
	cb.ax.set(yticklabels=lab)
	cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
	cb.outline.set_edgecolor('w')
	for i in range(cb_bnd.size):
		ax[1].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)
	pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
	ax[1].set(position=pos2)
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_' + nam,'png',900)
	return fig,ax

#%%
def Plot_BurnSeverity(meta,roi,vnam):
	z0=roi['grd'][vnam]['Data']	
	lab=list(meta['LUT']['Derived']['burnsev_comp1'].keys())
	z1=z0
	
	ind=np.where(roi['grd']['Data']==0)
	z1[ind]=6
	z1[0,0:6]=np.arange(1,7)

	N_vis=4
	N_hidden=2
	N_tot=N_vis+N_hidden
	
	lab.append('')
	
	cm=np.vstack( ((0.85,1,0.55,1),(1,0.9,0.65,1),(1,0.5,0.4,1),(0.5,0,0,1),(0.92,0.92,0.92,1),(1,1,1,1)) )
	cm=matplotlib.colors.ListedColormap(cm)

	plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
	im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
	roi=PlotVectorBaseMaps(meta,roi,ax[0])
	ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

	zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
	cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
	cb=plt.colorbar(im,cax=ax[1],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
	#ax[1].set(position=[0.71,0.6,0.05,0.14])
	cb.ax.set(yticklabels=lab)
	cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
	cb.outline.set_edgecolor('w')
	for i in range(cb_bnd.size):
		ax[1].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)
	pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
	ax[1].set(position=pos2)
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_' + vnam,'png',900)
	
	return fig,ax

#%%
def Plot_HarvestYear(meta,roi,vnam):
	z0=roi['grd'][vnam]['Data']
	bw=10; bin=np.arange(1855,2025+bw,bw);
	
	N_vis=bin.size
	N_hidden=3
	N_tot=N_vis+N_hidden
	
	z1=N_vis*np.ones(z0.shape)
	for i in range(bin.size):
		ind=np.where(np.abs(z0-bin[i])<=bw/2)
		z1[ind]=i
	z1[(roi['grd']['Data']==1) & (z0==0)]=i+2
	z1[(roi['grd']['Data']==1) & (roi['grd']['lc_comp1_2019']['Data']!=meta['LUT']['Derived']['lc_comp1']['Forest'])]=i+1	
	z1[(roi['grd']['Data']==0) | (roi['grd']['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Water'])]=i+3
	z1[0,0:N_tot]=np.arange(0,N_tot,1)

	lab=bin.astype(str)

	cm=plt.cm.get_cmap('plasma',i)
	cm=np.vstack( (cm.colors,(0.25,0.25,0.25,1),(0.75,0.75,0.75,1),(1,1,1,1)) )
	cm=matplotlib.colors.ListedColormap(cm)

	plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
	im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm,clim=(0,N_tot)) #
	roi=PlotVectorBaseMaps(meta,roi,ax[0])
	ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])
	
	cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_vis+1,1),ticks=np.arange(0.5,N_vis,1))
	cb.ax.set(yticklabels=lab)
	cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
	cb.outline.set_edgecolor('w')
	for i in range(0,N_vis):
		ax[1].plot([0,100],[i,i],'w-',linewidth=1.5)
	pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
	ax[1].set(position=pos2)
		
	gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_' + vnam,'png',900)
	return fig,ax

#%%
def Plot_ReserveComp1(meta,roi,vnam):

	lab=list(meta['LUT']['Derived']['reserve_comp1'].keys())
	z1=roi['grd'][vnam]['Data']
	z1[roi['grd']['Data']==0]=9
	z1[0,0:9]=np.arange(1,10)

	N_vis=8
	N_hidden=1
	N_tot=N_vis+N_hidden

	cm=np.vstack( ((0.6,1,0,1),(1,0.5,0,1),(0.5,0.8,1,1),(0,0.5,0,1),(0.75,0,1,1),(0.4,0.4,0.4,1),(0.6,0.6,0.6,1),(0.83,0.83,0.83,1),(1,1,1,1)) )
	cm=matplotlib.colors.ListedColormap(cm)

	plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
	im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
	roi=PlotVectorBaseMaps(meta,roi,ax[0])
	ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'],xticklabels='',yticklabels='')
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

	zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
	cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
	cb=plt.colorbar(im,cax=ax[1],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
	ax[1].set(position=[0.71,0.6,0.05,0.14])
	cb.ax.set(yticklabels=lab)
	cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
	cb.outline.set_edgecolor('w')
	for i in range(cb_bnd.size):
		ax[1].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)
	pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
	ax[1].set(position=pos2)

	gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_' + vnam,'png',900)

	return fig,ax

#%%
def Plot_REARs(meta,roi,vnam):
	if 'rears' in roi['grd'].keys():
		z1=roi['grd'][vnam]['Data']
	else:
		z1=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\REARs_CompPlusProp.tif')['Data']
		z1=gis.UpdateGridCellsize(z1,meta['Graphics']['Map']['RGSF'])['Data']
	z1[(roi['grd']['Data']==1) & (z1==0) & (roi['grd']['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Forest'])]=3
	z1[(roi['grd']['Data']==1) & (roi['grd']['lc_comp1_2019']['Data']!=meta['LUT']['Derived']['lc_comp1']['Forest'])]=4
	z1[(roi['grd']['Data']==0) | (roi['grd']['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Water'])]=5

	lab=['Protected','Protected (proposed)','Unprotected forest','Non-forest land']

	N_vis=4
	N_hidden=1
	N_tot=N_vis+N_hidden

	cm=np.vstack( ((0.7,0.6,1,1),(0.25,0.85,0.9,1),(0,0.4,0,1),(0.83,0.86,0.9,1),(1,1,1,1)) )
	cm=matplotlib.colors.ListedColormap(cm)

	plt.close('all')
	fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
	im=ax[0].matshow(z1[0::1,0::1],extent=roi['grd']['Extent'],cmap=cm)
	roi=PlotVectorBaseMaps(meta,roi,ax[0])
	ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

	zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
	cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
	cb=plt.colorbar(im,cax=ax[1],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
	ax[1].set(position=[0.71,0.6,0.05,0.14])
	cb.ax.set(yticklabels=lab)
	cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
	cb.outline.set_edgecolor('w')
	for i in range(cb_bnd.size):
		ax[1].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)

	pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
	ax[1].set(position=pos2)

	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_' + vnam,'png',900)

	return fig,ax

#%%
def Plot_Harvest_Po(meta,roi):
	z=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Disturbances\\HarvestProbability.tif')
	z=gis.UpdateGridCellsize(z,meta['Graphics']['Map']['RGSF'])
	sf=1000
	# Apply scale factor to convert to (%/yr)
	#z['Data']=z['Data'].astype('float')/1000
	bw=0.1; bin=np.arange(0,0.5+bw,bw)
	z1=np.ones( z['Data'].shape)
	for i in range(bin.size):
		ind=np.where(np.abs( z['Data']-bin[i]*sf)<=bw*sf/2)
		if ind[0].size==0:
			z1[ind[0][i],ind[1][i]]=i
		else:
			z1[ind]=i
	ind=np.where( z['Data']>=bin[i]*sf ); z1[ind]=i
	ind=np.where( roi['grd']['Data']!=1 ); z1[ind]=i+1

	lab=["%.2f" % x for x in bin]
	lab=np.append(lab,np.array(['Water']))

	N_vis=bin.size
	N_hidden=1
	N_tot=N_vis+N_hidden

	cm=plt.cm.get_cmap('viridis',bin.size)
	cm=np.vstack( (cm.colors,(1,1,1,1)) )
	cm=matplotlib.colors.ListedColormap(cm)

	plt.close('all'); fig,ax=plt.subplots(1,3,figsize=gu.cm2inch(18,18*0.5*roi['grd']['yxrat']))
	im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
	if meta['Graphics']['Map']['Show Bound Land Mask']=='On':
		roi['gdf']['bc_bound'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
	ax[0].set(position=[0,0,0.5,1],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'],aspect='auto')
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

	zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
	cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
	cb=plt.colorbar(im,cax=ax[1],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
	ax[1].set(position=[0.71,0.6,0.05,0.14])
	cb.ax.set(yticklabels=lab)
	cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
	cb.outline.set_edgecolor('w')
	for i in range(cb_bnd.size):
		ax[1].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)
	pos2=[0.39,0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],0.02,N_vis*meta['Graphics']['Map']['Legend Text Space']]
	ax[1].set(position=pos2)

	# Model
	V_Merch=np.arange(1,1200)
	beta=[0.005,-0.04,400]
	Po=beta[0]*(1/(1+np.exp(beta[1]*(V_Merch-beta[2]))))
	ax[2].plot(V_Merch,Po*100,'k-',linewidth=1,label='Harvest on-the-fly model 1')
	beta=[0.005,-0.04,500]
	Po=beta[0]*(1/(1+np.exp(beta[1]*(V_Merch-beta[2]))))
	ax[2].plot(V_Merch,Po*100,'g--',linewidth=1,label='Harvest on-the-fly model 2')
	ax[2].set(position=[0.55,0.12,0.4,0.86],xticks=np.arange(0,1300,100),xlabel='Merchantable volume (m$^3$ ha$^-$$^1$)',ylabel='Annual probability of harvest (%)',xlim=[0,800],ylim=[0,0.7])
	ax[2].legend(loc='upper left',bbox_to_anchor=(0.06,0.92),frameon=False,facecolor='w')
	ax[2].yaxis.set_ticks_position('both'); ax[2].xaxis.set_ticks_position('both'); ax[2].tick_params(length=meta['Graphics']['gp']['tickl'])

	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\taz_ann_prob_harvest','png',900)
	plt.show()
	return

#%%
def AddScalebar(ax,roi):

	dx=roi['grd']['xmax']-roi['grd']['xmin']
	dy=roi['grd']['ymax']-roi['grd']['ymin']

	if dy>1000000:
		x0=0.9 # Starting relative x
		y0=0.04 # Starting relative y (from bottom)
		sb_width=100*1000 # Distance of scalebar (meters)
		y_text_offset=15000
		y_offset=1*y_text_offset
	elif (dy>200000) & (dy<200000):
		x0=0.9 # Starting relative x
		y0=0.07 # Starting relative y (from bottom)
		sb_width=50*1000 # Distance of scalebar (meters)
		y_text_offset=5000
		y_offset=1*y_text_offset
	else:
		x0=0.9 # Starting relative x
		y0=0.07 # Starting relative y (from bottom)
		sb_width=50*1000 # Distance of scalebar (meters)
		y_text_offset=300
		y_offset=5*y_text_offset

	sb_x=[roi['grd']['xmin']+x0*dx-sb_width,roi['grd']['xmin']+x0*dx]
	sb_y=[roi['grd']['ymin']+y0*dy,roi['grd']['ymin']+y0*dy]
	ax.plot(sb_x,sb_y,'k-')
	ax.plot([sb_x[0],sb_x[0]],[sb_y[0]-y_offset,sb_y[0]],'k-')
	ax.plot([sb_x[1],sb_x[1]],[sb_y[0]-y_offset,sb_y[0]],'k-')
	ax.text(np.mean(sb_x),sb_y[0]+y_text_offset,str(int(sb_width/1000)) + ' km',ha='center')
	return

#%%
def Plot_BGC_Zone(meta,roi,vnam):

	z0=roi['grd'][vnam]['Data']
	z0[(roi['grd']['Data']==0) | (roi['grd']['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Water'])]=0

	lab0=np.array(list(meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE'].keys()))
	cl0=np.column_stack([meta['LUT']['Raw']['bgc_zone']['R'],meta['LUT']['Raw']['bgc_zone']['G'],meta['LUT']['Raw']['bgc_zone']['B']])
	id0=np.array(list(meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE'].values()))

	uid=np.unique(z0[z0!=0])
	id1=np.zeros(uid.size)
	lab1=np.array(['' for _ in range(uid.size)],dtype=object)
	cl1=np.zeros((uid.size,3))
	z1=(uid.size+1)*np.ones(z0.shape)
	for i in range(uid.size):
		ind=np.where(z0==uid[i])
		z1[ind]=i+1
		ind=np.where(id0==uid[i])[0][0]
		id1[i]=id0[ind]
		lab1[i]=lab0[ind]
		cl1[i,:]=cl0[ind,:]

	d=gu.CountByCategories(z1.flatten(),'Percent')

	N_vis=uid.size
	N_hidden=1
	N_tot=N_vis+N_hidden
	#lab1=np.append(lab1,[''])

	# Colormap
	cm=plt.cm.get_cmap('viridis',N_vis);
	for i in range(N_vis):
		cm.colors[i,0:3]=cl1[i,:]
	cm=np.vstack( (cm.colors,(1,1,1,1)) )
	cm=matplotlib.colors.ListedColormap(cm)

	if (meta['Graphics']['Map']['Show Inset Map']=='On'):
		N_panel=3
	else:
		N_panel=2

	plt.close('all'); fig,ax=plt.subplots(1,N_panel,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
	im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
	roi=PlotVectorBaseMaps(meta,roi,ax[0])
	ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

	# Legend
	zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
	cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
	cb=plt.colorbar(im,cax=ax[1],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
	ax[1].set(position=[0.71,0.6,0.05,0.14])
	cb.ax.set(yticklabels=lab1)
	cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
	cb.outline.set_edgecolor('w')
	for i in range(cb_bnd.size):
		ax[1].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)
	pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
	ax[1].set(position=pos2)

	# Scalebar
	if meta['Graphics']['Map']['Show Scalebar']=='On':
		AddScalebar(ax[0],roi)

	# Inset
	if meta['Graphics']['Map']['Show Inset Map']=='On':
		meta['gdf']['bc_bound']['gdf'].plot(ax=ax[2],edgecolor=None,facecolor=[0.8,0.8,0.8],linewidth=0.25)
		minx,miny,maxx,maxy=meta['gdf']['bc_bound']['gdf'].geometry.total_bounds
		roi['gdf']['bound within'].plot(ax=ax[2],edgecolor=None,facecolor=[0,0,0],linewidth=0.25)
		# Lower right
		#pos=[meta['Graphics']['Map']['Legend X'],0.1,1-meta['Graphics']['Map']['Legend X']-0.01,1-meta['Graphics']['Map']['Legend X']+0.05]
		# Lower left
		pos=[0,0,1-meta['Graphics']['Map']['Legend X']-0.01,1-meta['Graphics']['Map']['Legend X']+0.05]
		ax[2].set(position=pos,xlim=[minx,maxx],ylim=[miny,maxy])
		ax[2].yaxis.set_ticks_position('both'); ax[2].xaxis.set_ticks_position('both'); ax[2].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[2].axis(meta['Graphics']['Map']['Map Axis Vis'])

	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_' + vnam,'png',900)

	return fig,ax

#%%
def Plot_Ownership(meta,roi,vnam):

	z1=roi['grd'][vnam]['Data']
	lab1=list(meta['LUT']['F_OWN']['OWNERSHIP_DESCRIPTION'].keys())	
	cl1=np.array([[0.95,0.9,1],[0.9,0.85,1],[0.85,0.8,1],[0.2,0.5,0],[0.75,0.7,1],[0.7,0.65,1],[0.3,0.6,0],[0.6,0.55,1],[0.55,0.5,1],
				  [1,1,0.8],[0.95,0.95,0.7],[0.95,0.95,0.6],[0.95,0.95,0.5],[0.95,0.95,0.4],[0.9,0.9,0.3],[0.9,0.9,0.2],[0.85,0.85,0.1],				  
				  [0,0.75,0],[0,0.5,0],[0,0.25,0],[0.85,0.65,0.75],[0.7,0.5,0.6],[0.6,0.8,0.8],[0.7,0.9,0.8],[0.8,1,0.9],[0.9,1,0.65],[0.9,1,0.6],
				  [0.85,0.95,1],[0.8,0.9,1],[0.75,0.85,1],[0.7,0.8,1],[0.65,0.75,1],
				  [1,0.9,0.8],[1,0.8,0.7],[1,0.7,0.6],[0.9,0.9,0.9]])
		
	N_vis=len(lab1)
	N_hidden=1
	N_tot=N_vis+N_hidden
	#ind=np.where(z1==0); z1[ind]=N_vis+1
	#lab1=lab1[1:]
	#cl1=cl1[1:,:]
	
	z1[(roi['grd']['Data']==0) | (roi['grd']['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Water'])]=N_vis+1
	for i in range(len(lab1)):
		z1[i,0]=i+1
	#lab1=np.append(lab1,[''])
	cm=plt.cm.get_cmap('viridis',N_vis);
	for i in range(N_vis):
		cm.colors[i,0:3]=cl1[i,:]
	cm=np.vstack( (cm.colors,(1,1,1,1)) )
	cm=matplotlib.colors.ListedColormap(cm)

	plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
	im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
	#roi=PlotVectorBaseMaps(meta,roi,ax[0])
	ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

	zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
	cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
	cb=plt.colorbar(im,cax=ax[1],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
	ax[1].set(position=[0.71,0.01,0.03,0.98])
	cb.ax.set(yticklabels=lab1)
	cb.ax.tick_params(labelsize=4,length=0)
	cb.outline.set_edgecolor('w')
	for i in range(cb_bnd.size):
		ax[1].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)
	pos2=[0.66,0.45,0.01,0.54]
	ax[1].set(position=pos2)
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_' + vnam,'png',900)
	return fig,ax

#%%
def Plot_ASET(meta,roi,vnam):

	N_vis=meta['LUT']['Raw']['ASET']['ID'].size
	N_hidden=2
	N_tot=N_vis+N_hidden

	lab=list(meta['LUT']['Raw']['ASET']['Name'])
	z1=roi['grd'][vnam]['Data'].copy()
	z1[(z1==0) | (z1>N_vis)]=N_vis+1
	z1[roi['grd']['Data']==0]=N_vis+2
	z1[0,0:N_vis]=meta['LUT']['Raw']['ASET']['ID']
	np.unique(z1)

	# Colormap
	cl=np.array([meta['LUT']['Raw']['ASET']['c1'],meta['LUT']['Raw']['ASET']['c2'],meta['LUT']['Raw']['ASET']['c3']]).T
	cl=np.column_stack((cl,np.ones(cl.shape[0])))
	cm=np.vstack( (cl,(0.25,0.25,0.25,1),(1,1,1,1)) )
	cm=matplotlib.colors.ListedColormap(cm)

	plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
	im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
	roi=PlotVectorBaseMaps(meta,roi,ax[0])
	ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'],xticklabels='',yticklabels='')
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

	zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
	cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
	cb=plt.colorbar(im,cax=ax[1],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)

	try:
		cb.ax.set(yticklabels=lab)
	except:
		lab.append('')
		cb.ax.set(yticklabels=lab)
	cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
	cb.outline.set_edgecolor('w')
	for i in range(cb_bnd.size):
		ax[1].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)

	pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
	pos2[0]=0.69
	ax[1].set(position=pos2)

	gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_' + vnam,'png',900)

	return fig,ax

#%%
def Plot_ASET_Wildfire(meta,roi,vnam):

	vnam2='aset_post17_wf'
	N_vis=meta['LUT']['Raw']['ASET']['ID'].size
	N_hidden=3
	N_tot=N_vis+N_hidden

	lab=list(meta['LUT']['Raw']['ASET']['Name'])
	z1=roi['grd'][vnam]['Data'].copy()
	z1[(z1==0) | (z1>N_vis) | (roi['grd']['fire_yl']['Data']<2017)]=N_vis+2
	ind=np.where( (roi['grd'][vnam]['Data']==0) & (roi['grd']['fire_yl']['Data']>=2017) )
	z1[ind]=N_vis+1
	z1[roi['grd']['Data']==0]=N_vis+3
	z1[0,0:N_vis]=meta['LUT']['Raw']['ASET']['ID']
	#np.unique(z1)
	lab.append('')
	lab.append('')

	# Colormap
	cl=np.array([meta['LUT']['Raw']['ASET']['c1'],meta['LUT']['Raw']['ASET']['c2'],meta['LUT']['Raw']['ASET']['c3']]).T
	cl=np.column_stack((cl,np.ones(cl.shape[0])))
	cm=np.vstack( (cl,(0.25,0.25,0.25,1),(0.95,0.95,0.95,1),(1,1,1,1)) )
	cm=matplotlib.colors.ListedColormap(cm)

	plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
	im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
	roi=PlotVectorBaseMaps(meta,roi,ax[0])
	ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'],xticklabels='',yticklabels='')
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

	zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
	cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
	cb=plt.colorbar(im,cax=ax[1],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)

	try:
		cb.ax.set(yticklabels=lab)
	except:
		lab.append('')
		cb.ax.set(yticklabels=lab)
	cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
	cb.outline.set_edgecolor('w')
	for i in range(cb_bnd.size):
		ax[1].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)

	pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
	pos2[0]=0.69
	ax[1].set(position=pos2)

	gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_' + vnam2,'png',900)

	return fig,ax

#%%
def Plot_SalvageLogging(meta,roi,vnam):

	lab=['Salvage (>50% dead)','Salvage (>10% dead)','Harvest (<10% dead)','No harvesting','Non-forest land','Outside']
	lab=lab[0:-1]
	z1=roi['grd'][vnam]['Data']
	z1[roi['grd']['Data']==0]=6

	# Number of colours and number of colours excluded from colorbar
	N_vis=4
	N_hidden=2
	N_tot=N_vis+N_hidden

	# Colormap
	cm=np.vstack( ((0.9,0,0,1),(1,0.85,0,1),(0.5,0.5,0.5,1),(0.75,0.75,0.75,1), (0.93,0.93,0.93,1),(1,1,1,1)) )
	cm=matplotlib.colors.ListedColormap(cm)

	plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
	im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
	roi=PlotVectorBaseMaps(meta,roi,ax[0])
	ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'],xticklabels='',yticklabels='')
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

	zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
	cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
	cb=plt.colorbar(im,cax=ax[1],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)

	try:
		cb.ax.set(yticklabels=lab)
	except:
		lab.append('')
		cb.ax.set(yticklabels=lab)
	cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
	cb.outline.set_edgecolor('w')
	for i in range(cb_bnd.size):
		ax[1].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)

	pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
	ax[1].set(position=pos2)

	gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_' + vnam,'png',900)

	return fig,ax

#%%
def Plot_formask(meta,roi,vnam):

	z0=roi['grd'][vnam]['Data']

	bw=10; bin=np.arange(0,100+bw,bw)

	N_vis=bin.size
	N_hidden=2
	N_tot=N_vis+N_hidden

	z1=N_vis*np.ones(z0.shape)
	for i in range(N_vis):
		ind=np.where(np.abs(z0-bin[i])<=bw/2)
		if ind[0].size>0:
			z1[ind]=i+1
	ind=np.where(z0>=bin[i]); z1[ind]=i+1
	z1[1,1]=i+2
	z1[(roi['grd']['Data']==0) | (roi['grd']['lc_comp1_2019']['Data']==0) | (roi['grd']['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Water'])]=i+3

	for i in range(N_vis):
		z1[0,i]=i+1

	lab=['']*(N_tot-1)
	lab[0:N_vis]=bin.astype(str)

	cm=plt.cm.get_cmap('viridis',N_vis)
	cm=np.vstack( (cm.colors,(0.9,0.9,0.9,1),(1,1,1,1)) )
	cm=matplotlib.colors.ListedColormap(cm)

	plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
	im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
	roi=PlotVectorBaseMaps(meta,roi,ax[0])
	ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

	zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
	cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
	cb=plt.colorbar(im,cax=ax[1],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
	ax[1].set(position=[0.71,0.6,0.05,0.14])
	cb.ax.set(yticklabels=lab)
	cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
	cb.outline.set_edgecolor('w')
	for i in range(cb_bnd.size):
		ax[1].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)

	pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
	ax[1].set(position=pos2)

	gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_' + vnam,'png',900)

	return fig,ax

#%%
def Plot_NDEP(meta,roi,vnam):

	z0=roi['grd'][vnam]['Data']

	bw=0.5; bin=np.arange(0,10+bw,bw)
	N_vis=bin.size
	N_hidden=2
	N_tot=N_vis+N_hidden

	z1=N_vis*np.ones(z0.shape)
	for i in range(N_vis):
		ind=np.where(np.abs(z0-bin[i])<=bw/2)
		if ind[0].size>0:
			z1[ind]=i+1
	ind=np.where(z0>=bin[i]); z1[ind]=i+1
	z1[1,1]=i+2
	z1[(roi['grd']['Data']==0) | (roi['grd']['lc_comp1_2019']['Data']==0) | (roi['grd']['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Water'])]=i+3

	for i in range(N_vis):
		z1[0,i]=i+1

	lab=['']*(N_tot-1)
	lab[0:N_vis]=bin.astype(str)

	#cmn='cividis'
	cmn='viridis'
	cm=plt.cm.get_cmap(cmn,N_vis)
	cm=np.vstack( (cm.colors,(0.9,0.9,0.9,1),(1,1,1,1)) )
	cm=matplotlib.colors.ListedColormap(cm)

	plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
	im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
	roi=PlotVectorBaseMaps(meta,roi,ax[0])
	ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

	zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
	cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
	cb=plt.colorbar(im,cax=ax[1],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
	ax[1].set(position=[0.71,0.6,0.05,0.14])
	cb.ax.set(yticklabels=lab)
	cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
	cb.outline.set_edgecolor('w')
	for i in range(cb_bnd.size):
		ax[1].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)

	pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
	ax[1].set(position=pos2)

	gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_' + vnam,'png',900)

	return fig,ax

#%%
def Plot_MAP(meta,roi,vnam):

	z0=12*roi['grd'][vnam]['Data']

	bw=200; bin=np.arange(0,2400+bw,bw)

	N_vis=bin.size
	N_hidden=2
	N_tot=N_vis+N_hidden

	z1=N_vis*np.ones(z0.shape)
	for i in range(N_vis):
		ind=np.where(np.abs(z0-bin[i])<=bw/2)
		if ind[0].size>0:
			z1[ind]=i+1
	ind=np.where(z0>=bin[i]); z1[ind]=i+1
	z1[1,1]=i+2
	z1[(roi['grd']['Data']==0) | (roi['grd']['lc_comp1_2019']['Data']==0) | (roi['grd']['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Water'])]=i+3

	for i in range(N_vis):
		z1[0,i]=i+1

	lab=['']*(N_tot-1)
	lab[0:N_vis]=bin.astype(str)

	cm=plt.cm.get_cmap('viridis',N_vis)
	cm=np.vstack( (cm.colors,(0.9,0.9,0.9,1),(1,1,1,1)) )
	cm=matplotlib.colors.ListedColormap(cm)

	plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
	im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
	roi=PlotVectorBaseMaps(meta,roi,ax[0])
	ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

	zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
	cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
	cb=plt.colorbar(im,cax=ax[1],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
	ax[1].set(position=[0.71,0.6,0.05,0.14])
	cb.ax.set(yticklabels=lab)
	cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
	cb.outline.set_edgecolor('w')
	for i in range(cb_bnd.size):
		ax[1].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)

	pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
	ax[1].set(position=pos2)

	gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_' + vnam,'png',900)

	return fig,ax

#%%
def Plot_MAT(meta,roi,vnam):

	z0=roi['grd'][vnam]['Data']
	bw=1; bin=np.arange(-5,12+bw,bw)
	N_vis=bin.size
	N_hidden=2
	N_tot=N_vis+N_hidden

	z1=N_vis*np.ones(z0.shape)
	for i in range(N_vis):
		ind=np.where(np.abs(z0-bin[i])<=bw/2)
		if ind[0].size>0:
			z1[ind]=i+1
	ind=np.where(z0<bin[0]); z1[ind]=0
	ind=np.where(z0>=bin[i]); z1[ind]=i+1
	z1[1,1]=i+2
	z1[(roi['grd']['Data']==0) | (roi['grd']['lc_comp1_2019']['Data']==0) | (roi['grd']['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Water'])]=i+3

	for i in range(N_vis):
		z1[0,i]=i+1

	lab=['']*(N_tot-1)
	lab[0:N_vis]=np.array(bin).astype(str)

	cm=plt.cm.get_cmap('viridis',N_vis)
	cm=np.vstack( (cm.colors,(0.9,0.9,0.9,1),(1,1,1,1)) )
	cm=matplotlib.colors.ListedColormap(cm)

	plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
	im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
	roi=PlotVectorBaseMaps(meta,roi,ax[0])
	ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

	zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
	cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
	cb=plt.colorbar(im,cax=ax[1],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
	ax[1].set(position=[0.71,0.6,0.05,0.14])
	cb.ax.set(yticklabels=lab)
	cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
	cb.outline.set_edgecolor('w')
	for i in range(cb_bnd.size):
		ax[1].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)

	pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
	ax[1].set(position=pos2)

	gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_' + vnam,'png',900)
	plt.show()
	return fig,ax

#%%
def Plot_Fire2023(meta,roi):

	z1=4*np.ones(roi['grd']['Data'].shape,dtype='int8')

	ind=np.where( (roi['grd']['lc_comp1_2019']['Data']==1) ); z1[ind]=1
	ind=np.where( (roi['grd']['lc_comp1_2019']['Data']==1) & (roi['grd']['fire_2023']['Data']>0) ); z1[ind]=2
	ind=np.where( (roi['grd']['lc_comp1_2019']['Data']!=1) ); z1[ind]=3
	ind=np.where( (roi['grd']['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Water']) ); z1[ind]=4

	lab=['Forest land','Forest land (affected)','Non-forest land']

	N_vis=len(lab)
	N_hidden=1
	N_tot=N_vis+N_hidden

	cm=np.vstack( ((0.8,0.8,0.8,1),(0.85,0,0,1),(0.9,0.9,0.9,1),(1,1,1,1)) )
	cm=matplotlib.colors.ListedColormap(cm)

	plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
	im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
	roi=PlotVectorBaseMaps(meta,roi,ax[0])
	gdf['cities']['gdf'][gdf['cities']['gdf']['Territory']=='BC'].plot(ax=ax[0],marker='s',edgecolor=[0,0,0],facecolor=[1,1,0],lw=0.25,markersize=9,alpha=1,zorder=2)
	for x,y,label in zip(gdf['cities']['gdf'].geometry.x,gdf['cities']['gdf'].geometry.y,gdf['cities']['gdf']['City Name']):
		ax[0].annotate(label,xy=(x,y),xytext=(3,2),textcoords="offset points",color=[0,0,0],fontsize=4)

	ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

	zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
	cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
	cb=plt.colorbar(im,cax=ax[1],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
	ax[1].set(position=[0.71,0.6,0.05,0.14])
	cb.ax.set(yticklabels=lab)
	cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
	cb.outline.set_edgecolor('w')
	for i in range(cb_bnd.size):
		ax[1].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)
	pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
	ax[1].set(position=pos2)
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_fire2023','png',900)

	return fig,ax

#%% Plot FECA year
def Plot_FECA_Year(meta,roi,vnam):

	z0=roi['grd'][vnam]['Data']
	bw=5; bin=np.arange(1975,2025+bw,bw);
	z1=(bin.size)*np.ones(z0.shape)
	for i in range(bin.size):
		ind=np.where(np.abs(z0-bin[i])<=bw/2)
		z1[ind]=i
	z1[(z0==0)]=i+1
	z1[(roi['grd']['Data']==0) | (roi['grd']['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Water'])]=i+2
	for i in range(bin.size):
		z1[0,i]=i

	N_vis=bin.size
	N_hidden=2
	N_tot=N_vis+N_hidden
	lab=bin.astype(str)
	lab=np.append(lab,'Hidden')
	lab=np.append(lab,'Hidden')
	
	#cm=plt.cm.get_cmap('viridis',N_vis)
	cm=plt.cm.get_cmap('plasma',N_vis)
	cm=np.vstack( (cm.colors,(0.8,0.8,0.8,1),(1,1,1,1)) )
	cm=matplotlib.colors.ListedColormap(cm)

	plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
	im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
	roi=PlotVectorBaseMaps(meta,roi,ax[0])
	ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

	zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
	cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
	cb=plt.colorbar(im,cax=ax[1],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
	ax[1].set(position=[0.71,0.6,0.05,0.14])
	cb.ax.set(yticklabels=lab)
	cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
	cb.outline.set_edgecolor('w')
	for i in range(cb_bnd.size):
		ax[1].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)
		
	pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
	ax[1].set(position=pos2)
	gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_' + vnam,'png',900)
	return fig,ax

#%%
def Plot_Spc1_NTEMS(meta,roi,vnam):

	# Compress categories
	z0=roi['grd'][vnam]['Data']
	lab0=np.array(list(meta['LUT']['Derived']['spc1_ntems'].keys()))
	id0=np.array(list(meta['LUT']['Derived']['spc1_ntems'].values()))
	cl0=np.column_stack([meta['LUT']['Raw']['spc1_ntems']['R'],meta['LUT']['Raw']['spc1_ntems']['G'],meta['LUT']['Raw']['spc1_ntems']['B']])
	z1,lab1,cl1=gis.CompressCats(z0,id0,lab0,cl0)
	z1=z1+1

	n=int(np.max(z1))
	z1[0,0]=n+1
	z1[(roi['grd']['Data']==0) | (roi['grd']['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Water'])]=n+2

	lab1=np.append(lab1,[''])

	N_vis=n
	N_hidden=2
	N_tot=N_vis+N_hidden

	# Colormap
	cm=plt.cm.get_cmap('viridis',N_vis);
	for i in range(N_vis):
		cm.colors[i,0:3]=cl1[i,:]
	cm=np.vstack( (cm.colors,(0.9,0.9,0.9,1),(1,1,1,1)) )
	cm=matplotlib.colors.ListedColormap(cm)

	plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
	im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
	roi=PlotVectorBaseMaps(meta,roi,ax[0])
	ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

	zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
	cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
	cb=plt.colorbar(im,cax=ax[1],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
	ax[1].set(position=[0.71,0.6,0.05,0.14])
	cb.ax.set(yticklabels=lab1)
	cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
	cb.outline.set_edgecolor('w')
	for i in range(cb_bnd.size):
		ax[1].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)

	pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
	ax[1].set(position=pos2)

	gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_' + vnam,'png',900)
	plt.show()
	return fig,ax

#%%
def Plot_SI(meta,roi):

	#z0=roi['grd']['SITE_INDEX']['Data']
	#z0=roi['grd']['si_spl_fd']['Data']
	#z0=roi['grd']['si_spl_hw']['Data']
	z0=roi['grd']['si_spl_ntems']['Data']

	bw=2
	bin=np.arange(10,40+bw,bw)
	id=np.arange(1,bin.size+1)

	N_vis=bin.size
	N_hidden=2
	N_tot=N_vis+N_hidden

	z1=N_vis*np.ones(z0.shape)
	for i in range(N_vis):
		ind=np.where(np.abs(z0-bin[i])<=bw/2)
		if ind[0].size>0:
			z1[ind]=id[i]
		else:
			z1[0,i]=id[i]
	ind=np.where(z0>bin[i]); z1[ind]=id[-1]
	z1[(z0==0) | (roi['grd']['lc_comp1_2019']['Data']!=meta['LUT']['Derived']['lc_comp1']['Forest'])]=id[-1]+1
	#z1[1,1]=id[-1]+1
	z1[(roi['grd']['Data']==0) | (roi['grd']['lc_comp1_2019']['Data']==0) | (roi['grd']['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Water'])]=id[-1]+2

	lab=['']*(N_vis+1)
	lab[0:N_vis]=bin.astype(str)

	cm=plt.cm.get_cmap('viridis',N_vis)
	cm=np.vstack( (cm.colors,(0.9,0.9,0.9,1),(1,1,1,1)) )
	cm=matplotlib.colors.ListedColormap(cm)

	plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
	im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
	roi=PlotVectorBaseMaps(meta,roi,ax[0])
	ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

	zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
	cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
	cb=plt.colorbar(im,cax=ax[1],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
	ax[1].set(position=[0.71,0.6,0.05,0.14])
	cb.ax.set(yticklabels=lab)
	cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
	cb.outline.set_edgecolor('w')
	for i in range(cb_bnd.size):
		ax[1].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)

	pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
	ax[1].set(position=pos2)

	#gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_SITE_INDEX','png',900)
	#gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_si_spl_hw','png',900)
	gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_si_spl_ntems','png',900)

	return fig,ax

#%%
def Plot_Age(meta,roi,vnam):
	z0=roi['grd'][vnam]['Data']
	bw=20
	bin=np.arange(0,220,bw)
	id=np.arange(1,bin.size+1)

	N_vis=bin.size
	N_hidden=2
	N_tot=N_vis+N_hidden

	z1=N_vis*np.ones(z0.shape)
	for i in range(N_vis):
		ind=np.where(np.abs(z0-bin[i])<=bw/2)
		if ind[0].size>0:
			z1[ind]=id[i]
		else:
			z1[0,i]=id[i]
	ind=np.where(z0>bin[i]); z1[ind]=id[-1]
	z1[(z0==0) | (roi['grd']['lc_comp1_2019']['Data']!=meta['LUT']['Derived']['lc_comp1']['Forest'])]=id[-1]+1
	z1[(roi['grd']['Data']==0) | (roi['grd']['lc_comp1_2019']['Data']==0) | (roi['grd']['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Water'])]=id[-1]+2

	lab=['']*(N_tot-1)
	lab[0:N_vis]=bin.astype(str)

	cm=plt.cm.get_cmap('viridis',N_vis)
	cm=np.vstack( (cm.colors,(0.9,0.9,0.9,1),(1,1,1,1)) )
	cm=matplotlib.colors.ListedColormap(cm)

	plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
	im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
	roi=PlotVectorBaseMaps(meta,roi,ax[0])

	ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

	zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
	cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
	cb=plt.colorbar(im,cax=ax[1],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
	ax[1].set(position=[0.71,0.6,0.05,0.14])
	cb.ax.set(yticklabels=lab)
	cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
	cb.outline.set_edgecolor('w')
	for i in range(cb_bnd.size):
		ax[1].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)
	pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
	ax[1].set(position=pos2)
	gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_' + vnam,'png',900)
	return fig,ax

#%%
def Plot_Height(meta,roi,vnam):
	z0=roi['grd'][vnam]['Data']
	bw=5
	bin=np.arange(0,40+bw,bw)
	id=np.arange(1,bin.size+1)

	N_vis=bin.size
	N_hidden=2
	N_tot=N_vis+N_hidden

	z1=N_vis*np.ones(z0.shape)
	for i in range(N_vis):
		ind=np.where(np.abs(z0-bin[i])<=bw/2)
		if ind[0].size>0:
			z1[ind]=id[i]
		else:
			z1[0,i]=id[i]
	ind=np.where(z0>bin[i]); z1[ind]=id[-1]
	z1[(z0==0) | (roi['grd']['lc_comp1_2019']['Data']!=meta['LUT']['Derived']['lc_comp1']['Forest'])]=id[-1]+1
	z1[(roi['grd']['Data']==0) | (roi['grd']['lc_comp1_2019']['Data']==0) | (roi['grd']['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Water'])]=id[-1]+2

	lab=['']*(N_tot-1)
	lab[0:N_vis]=bin.astype(str)

	cm=plt.cm.get_cmap('viridis',N_vis)
	cm=np.vstack( (cm.colors,(0.9,0.9,0.9,1),(1,1,1,1)) )
	cm=matplotlib.colors.ListedColormap(cm)

	plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
	im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
	roi=PlotVectorBaseMaps(meta,roi,ax[0])

	ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

	zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
	cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
	cb=plt.colorbar(im,cax=ax[1],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
	ax[1].set(position=[0.71,0.6,0.05,0.14])
	cb.ax.set(yticklabels=lab)
	cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
	cb.outline.set_edgecolor('w')
	for i in range(cb_bnd.size):
		ax[1].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)
	pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
	ax[1].set(position=pos2)
	gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_' + vnam,'png',900)
	return fig,ax

#%%
def Plot_BroadleafDeciduousFrac(meta,roi,vnam):
	z0=roi['grd'][vnam]['Data']
	bw=5
	bin=np.arange(0,100+bw,bw)
	id=np.arange(1,bin.size+1)

	N_vis=bin.size
	N_hidden=2
	N_tot=N_vis+N_hidden

	z1=(N_vis+1)*np.ones(z0.shape)
	for i in range(N_vis):
		ind=np.where(np.abs(z0-bin[i])<=bw/2)
		if ind[0].size>0:
			z1[ind]=id[i]
		else:
			z1[0,i]=id[i]
	ind=np.where(z0>bin[i]); z1[ind]=id[-1]
	z1[(roi['grd']['lc_comp1_2019']['Data']!=meta['LUT']['Derived']['lc_comp1']['Forest'])]=N_vis+1
	z1[(roi['grd']['Data']==0)]=N_vis+2

	lab=['']*(N_tot-1)
	lab[0:N_vis]=bin.astype(str)

	cm0=plt.cm.get_cmap('turbo',N_vis)
	cm=np.zeros((N_vis,4))
	for i in range(N_vis):
		cm[i,:]=cm0(i)
	cm[0,:]=[0,0,0,1]
	cm=np.vstack( (cm,(0.0,0.0,0.0,1),(1,1,1,1)) )
	cm=matplotlib.colors.ListedColormap(cm)

	plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
	im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
	roi=PlotVectorBaseMaps(meta,roi,ax[0])

	ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

	zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
	cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
	cb=plt.colorbar(im,cax=ax[1],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
	ax[1].set(position=[0.71,0.6,0.05,0.14])
	cb.ax.set(yticklabels=lab)
	cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
	cb.outline.set_edgecolor('w')
	for i in range(cb_bnd.size):
		ax[1].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)
	pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
	ax[1].set(position=pos2)
	gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_' + vnam,'png',900)
	return fig,ax

#%%
def Plot_PlantedMask(meta,roi,vnam):

	lab=[]
	z1=2*np.ones(roi['grd'][vnam]['Data'].shape)
	z1[(roi['grd']['Data']==1) & (roi['grd'][vnam]['Data']>0)]=1; lab.append('Planted')
	z1[(roi['grd']['Data']==1) & (roi['grd'][vnam]['Data']==0)]=2; lab.append('Not planted')
	z1[(roi['grd']['Data']!=1)]=3;

	# Number of colours and number of colours excluded from colorbar
	N_vis=2
	N_hidden=1
	N_tot=N_vis+N_hidden

	# Colormap
	cm=np.vstack( ((0,0.6,0,1),(0.83,0.83,0.83,1),(1,1,1,1)) )
	cm=matplotlib.colors.ListedColormap(cm)

	plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
	im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
	roi=PlotVectorBaseMaps(meta,roi,ax[0])
	ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'],xticklabels='',yticklabels='')
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

	zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
	cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
	cb=plt.colorbar(im,cax=ax[1],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
	ax[1].set(position=[0.71,0.6,0.05,0.14])
	cb.ax.set(yticklabels=lab)
	cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
	cb.outline.set_edgecolor('w')
	for i in range(cb_bnd.size):
		ax[1].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)

	pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
	ax[1].set(position=pos2)

	gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_' + vnam,'png',900)

	return fig,ax

#%%
def Plot_Elev(meta,roi,vnam):
	z0=roi['grd'][vnam]['Data']
	bw=200; bin=np.arange(0,3200+bw,bw);
	
	N_vis=bin.size
	N_hidden=1
	N_tot=N_vis+N_hidden
	
	z1=N_vis*np.ones(z0.shape)
	for i in range(bin.size):
		ind=np.where(np.abs(z0-bin[i])<=bw/2)
		z1[ind]=i
	z1[(z1>bin[i])]=i
	#z1[(roi['grd']['Data']==1) & (z0==0)]=i+2
	#z1[(roi['grd']['Data']==1) & (roi['grd']['lc_comp1_2019']['Data']!=meta['LUT']['Derived']['lc_comp1']['Forest'])]=i+1	
	z1[(roi['grd']['Data']==0)]=i+1
	z1[0,0:N_tot]=np.arange(0,N_tot,1)

	lab=bin.astype(str)

	cm=plt.cm.get_cmap('Greys',i)
	cmc=np.zeros((N_vis,4))
	cnt=0
	for ivl in np.linspace(0,1,N_vis):
		cmc[cnt,:]=np.array(cm(ivl))
		cnt=cnt+1
	cm=np.vstack( (cmc,(1,1,1,1)) )
	cm=matplotlib.colors.ListedColormap(cm)

	plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
	im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm,clim=(0,N_tot)) #
	roi=PlotVectorBaseMaps(meta,roi,ax[0])
	ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])
	
	cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_vis+1,1),ticks=np.arange(0.5,N_vis,1))
	cb.ax.set(yticklabels=lab)
	cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
	cb.outline.set_edgecolor('w')
	for i in range(0,N_vis):
		ax[1].plot([0,100],[i,i],'w-',linewidth=1.5)
	pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
	ax[1].set(position=pos2)
		
	gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_' + vnam,'png',900)
	return fig,ax

#%%
def Plot_Slope(meta,roi,vnam):
	z0=roi['grd'][vnam]['Data']
	bw=2; bin=np.arange(0,20+bw,bw);
	
	N_vis=bin.size
	N_hidden=2
	N_tot=N_vis+N_hidden
	
	z1=N_vis*np.ones(z0.shape)
	for i in range(bin.size):
		ind=np.where(np.abs(z0-bin[i])<=bw/2)
		z1[ind]=i
	z1[(z1>bin[i])]=i
	#z1[(roi['grd']['Data']==1) & (z0==0)]=i+2
	#z1[(roi['grd']['Data']==1) & (roi['grd']['lc_comp1_2019']['Data']!=meta['LUT']['Derived']['lc_comp1']['Forest'])]=i+1	
	z1[z0>bin[-1]]=i+1
	z1[(roi['grd']['Data']==0)]=i+2
	z1[0,0:N_tot]=np.arange(0,N_tot,1)

	lab=bin.astype(str)

	cm=plt.cm.get_cmap('jet',i)
	cm.colors=cm(np.arange(0,cm.N))
	cm=np.vstack( (cm.colors,(0,0,0,1),(1,1,1,1)) )
	cm=matplotlib.colors.ListedColormap(cm)

	plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
	im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm,clim=(0,N_tot)) #
	roi=PlotVectorBaseMaps(meta,roi,ax[0])
	ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])
	
	cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_vis+1,1),ticks=np.arange(0.5,N_vis,1))
	cb.ax.set(yticklabels=lab)
	cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
	cb.outline.set_edgecolor('w')
	for i in range(0,N_vis):
		ax[1].plot([0,100],[i,i],'w-',linewidth=1.5)
	pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
	ax[1].set(position=pos2)
		
	gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_' + vnam,'png',900)
	return fig,ax

#%%
def Plot_Infastructure1(meta,roi,vnam):
	z0=roi['grd'][vnam]['Data']
	#bw=500; bin=np.arange(0,5000+bw,bw);
	bw=1; bin=np.arange(120,160+bw,bw);
	
	N_vis=bin.size
	N_hidden=1
	N_tot=N_vis+N_hidden
	
	z1=N_vis*np.ones(z0.shape)
	for i in range(bin.size):
		ind=np.where(np.abs(z0-bin[i])<=bw/2)
		z1[ind]=i
	z1[(z1>bin[i])]=i
	z1[(roi['grd']['Data']==0)]=i+1
	z1[0,0:N_tot]=np.arange(0,N_tot,1)

	lab=bin.astype(str)

	cm=plt.cm.get_cmap('Greys',i)
	cmc=np.zeros((N_vis,4))
	cnt=0
	for ivl in np.linspace(0,1,N_vis):
		cmc[cnt,:]=np.array(cm(ivl))
		cnt=cnt+1
	cm=np.vstack( (cmc,(1,1,1,1)) )
	cm=matplotlib.colors.ListedColormap(cm)

	plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
	im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm,clim=(0,N_tot)) #

	a=gis.OpenGeoTiff(r'D:\Data\dem_CV_p.tif')
	ax[0].matshow(a['Data'],extent=a['Extent'],cmap=cm,clim=(100,170))
	
	#roi['gdf']['bound within'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
	roi['gdf']['bc_bound'].plot(ax=ax[0],edgecolor=[0.6,0.6,0.6],facecolor='none',linewidth=0.25)
	roi['gdf']['lakes'].plot(ax=ax[0],facecolor=[0.82,0.88,1],label='Lakes',linewidth=0.25)
	#roi['gdf']['rivers'].plot(ax=ax[0],color=[0.6,0.8,1],label='Rivers',linewidth=0.25)

	roi['gdf']['road'].plot(ax=ax[0],edgecolor=[0.9,0.45,0],linewidth=0.25,label='Road',alpha=1,zorder=1)

	roi['gdf']['hydrol'].plot(ax=ax[0],linestyle='--',edgecolor=[0.27,0.45,0.8],linewidth=0.5,alpha=1,zorder=1,label='Road')

	roi['gdf']['rail'].plot(ax=ax[0],edgecolor=[0,0,0],linewidth=0.5,alpha=1,zorder=1,label='Rail')
	roi['gdf']['rail'].plot(ax=ax[0],edgecolor=[1,1,1],linewidth=0.25,alpha=1,zorder=1,label='Rail')

	roi['gdf']['tpf'].plot(ax=ax[0],marker='^',edgecolor=[0,0.75,0.75],facecolor=[0,1,1],linewidth=0.25,label='TPF',alpha=1,zorder=1,markersize=5)

	#for x,y,label in zip(roi['gdf']['tpf'].geometry.x,roi['gdf']['tpf'].geometry.y,roi['gdf']['tpf']['COMPANY_NAME']):
	#	ax[0].annotate(label,xy=(x,y),xytext=(4,3),textcoords="offset points",color=[0,0,0],fontsize=4)
	ind=np.where( (roi['gdf']['cities']['Territory']=='BC') & (roi['gdf']['cities']['Level']==1) )[0]
	roi['gdf']['cities'].iloc[ind].plot(ax=ax[0],marker='s',edgecolor=[0.5,0,0],facecolor=[1,0,0],lw=0.5,markersize=7,alpha=1,zorder=2)
	for x,y,label in zip(roi['gdf']['cities'].iloc[ind].geometry.x,roi['gdf']['cities'].iloc[ind].geometry.y,roi['gdf']['cities'].iloc[ind]['City Name']):
		ax[0].annotate(label,xy=(x,y),xytext=(2,1.5),textcoords="offset points",color=[0,0,0],fontsize=4)
	#roi['gdf']['popp'].plot(ax=ax[0],marker='o',edgecolor=[0.75,0.3,0],facecolor=[1,0.6,0],lw=0.25,markersize=6,alpha=1,zorder=2)
	#for x,y,label in zip(roi['gdf']['popp'].geometry.x,roi['gdf']['popp'].geometry.y,roi['gdf']['popp']['NAME']):
	#	ax[0].annotate(label,xy=(x,y),xytext=(2,1.5),textcoords="offset points",color=[0,0,0],fontsize=4)
	
	ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])
	
	cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_vis+1,1),ticks=np.arange(0.5,N_vis,1))
	cb.ax.set(yticklabels=lab)
	cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
	cb.outline.set_edgecolor('w')
	for i in range(0,N_vis):
		ax[1].plot([0,100],[i,i],'w-',linewidth=1.5)
	pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
	ax[1].set(position=pos2)
		
	gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_infastructure1','png',900)
	return fig,ax

#%%
def Plot_CHM(meta,roi,vnam):
	#z0=roi['grd'][vnam]['Data']
	chm=gis.OpenGeoTiff(r'C:\Data\LiDAR\chm_p.tif')
	z0=chm['Data'].copy()

	bw=5; bin=np.arange(0,40+bw,bw)
	N_vis=bin.size
	N_hidden=1
	N_tot=N_vis+N_hidden
	
	z1=N_vis*np.ones(z0.shape)
	for i in range(bin.size):
		ind=np.where(np.abs(z0-bin[i])<=bw/2)
		z1[ind]=i
	z1[(z0>bin[i])]=i
	#z1[(roi['grd']['Data']==0)]=i+1
	z1[0,0:N_tot]=np.arange(0,N_tot,1)

	lab=bin.astype(str)

	cm=plt.cm.get_cmap('viridis',i)
	cmc=np.zeros((N_vis,4))
	cnt=0
	for ivl in np.linspace(0,1,N_vis):
		cmc[cnt,:]=np.array(cm(ivl))
		cnt=cnt+1
	cm=np.vstack( (cmc,(1,1,1,1)) )
	cm=matplotlib.colors.ListedColormap(cm)

	plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
	#im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm,clim=(0,N_tot)) #
	#ax[0].matshow(chm['Data'],extent=chm['Extent'],cmap=cm,clim=(0,N_tot))
	im=ax[0].matshow(z1,extent=chm['Extent'],cmap=cm,clim=(0,N_tot))

	prop=gis.ConvertGeographicToGeojson(r'G:\My Drive\Property\Rainville\Property Boundary.xlsx',r'G:\My Drive\Property\Rainville\Property Boundary.geojson','Polygon')
	prop.plot(ax=ax[0],facecolor='none',edgecolor=[0,1,1],label='Property',linewidth=1,linestyle='-')

	#prop=gis.ConvertGeographicToGeojson(r'G:\My Drive\Property\Cowichan Lake Road Lot 3\Coordinates Geographic.xlsx',r'G:\My Drive\Property\Cowichan Lake Road Lot 3\Coordinates Geographic.geojson','Polygon')
	#prop.plot(ax=ax[0],facecolor='none',edgecolor=[1,0,0],label='Property',linewidth=1.25,linestyle='-')

	#row=gis.ConvertGeographicToGeojson(r'G:\My Drive\Property\Cowichan Lake Road Lot 3\ROW.xlsx',r'G:\My Drive\Property\Cowichan Lake Road Lot 3\ROW.geojson','Line')
	#row.plot(ax=ax[0],facecolor='none',edgecolor=[0.7,0.2,1],label='Property',linewidth=1,linestyle='-')

	#phot=gis.ConvertGeographicToGeojson(r'G:\My Drive\Property\Cowichan Lake Road Lot 3\Photos.xlsx',r'G:\My Drive\Property\Cowichan Lake Road Lot 3\Photos.geojson','Point')
	#phot.plot(ax=ax[0],marker='o',markersize=7,facecolor='none',edgecolor=[1,0,0],label='Property',linewidth=1,linestyle='-')

	#roi['gdf']['bound within'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
	roi['gdf']['bc_bound'].plot(ax=ax[0],edgecolor=[0.6,0.6,0.6],facecolor='none',linewidth=0.25)
	roi['gdf']['lakes'].plot(ax=ax[0],facecolor=[0.82,0.88,1],label='Lakes',linewidth=0.25)
	#roi['gdf']['rivers'].plot(ax=ax[0],color=[0.6,0.8,1],label='Rivers',linewidth=0.25)
	roi['gdf']['road'].plot(ax=ax[0],edgecolor=[0.9,0.45,0],linewidth=0.25,label='Road',alpha=1,zorder=1)
	roi['gdf']['hydrol'].plot(ax=ax[0],linestyle='--',edgecolor=[0.27,0.45,0.8],linewidth=0.5,alpha=1,zorder=1,label='Road')
	roi['gdf']['rail'].plot(ax=ax[0],edgecolor=[0,0,0],linewidth=0.5,alpha=1,zorder=1,label='Rail')
	roi['gdf']['rail'].plot(ax=ax[0],edgecolor=[1,1,1],linewidth=0.25,alpha=1,zorder=1,label='Rail')
	roi['gdf']['tpf'].plot(ax=ax[0],marker='^',edgecolor=[0,0.75,0.75],facecolor=[0,1,1],linewidth=0.25,label='TPF',alpha=1,zorder=1,markersize=5)
	#for x,y,label in zip(roi['gdf']['tpf'].geometry.x,roi['gdf']['tpf'].geometry.y,roi['gdf']['tpf']['COMPANY_NAME']):
	#	ax[0].annotate(label,xy=(x,y),xytext=(4,3),textcoords="offset points",color=[0,0,0],fontsize=4)
	ind=np.where( (roi['gdf']['cities']['Territory']=='BC') & (roi['gdf']['cities']['Level']==1) )[0]
	roi['gdf']['cities'].iloc[ind].plot(ax=ax[0],marker='s',edgecolor=[0.5,0,0],facecolor=[1,0,0],lw=0.5,markersize=7,alpha=1,zorder=2)
	for x,y,label in zip(roi['gdf']['cities'].iloc[ind].geometry.x,roi['gdf']['cities'].iloc[ind].geometry.y,roi['gdf']['cities'].iloc[ind]['City Name']):
		ax[0].annotate(label,xy=(x,y),xytext=(2,1.5),textcoords="offset points",color=[0,0,0],fontsize=4)
	#roi['gdf']['popp'].plot(ax=ax[0],marker='o',edgecolor=[0.75,0.3,0],facecolor=[1,0.6,0],lw=0.25,markersize=6,alpha=1,zorder=2)
	#for x,y,label in zip(roi['gdf']['popp'].geometry.x,roi['gdf']['popp'].geometry.y,roi['gdf']['popp']['NAME']):
	#	ax[0].annotate(label,xy=(x,y),xytext=(2,1.5),textcoords="offset points",color=[0,0,0],fontsize=4)

	#ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=chm['xlim'],ylim=chm['ylim'])
	ax[0].set(position=[0.1,0.1,0.55,0.8],xlim=chm['xlim'],ylim=chm['ylim'])
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both');
	#ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])
	#ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

	cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_vis+1,1),ticks=np.arange(0.5,N_vis,1))
	cb.ax.set(yticklabels=lab)
	cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
	cb.outline.set_edgecolor('w')
	for i in range(0,N_vis):
		ax[1].plot([0,100],[i,i],'w-',linewidth=1.5)

	pos2=[0.75,0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
	ax[1].set(position=pos2)
		
	#gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_chm_a','png',900)
	#gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_chm_b','png',900)
	#gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_chm_c','png',900)
	gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_chm_Rainville','png',900)
	return fig,ax

#%%
def Plot_hcl(meta,roi):

	chm=gis.OpenGeoTiff(r'C:\Data\LiDAR\chm_p.tif')
	#chm=gis.OpenGeoTiff(r'D:\Data\chm_p2.tif')

	N_vis=5
	N_hidden=1
	N_tot=N_vis+N_hidden

	chm1=N_vis*np.ones(chm['Data'].shape,dtype='int8')
	ind=np.where( (chm['Data']<0.5) ); chm1[ind]=1
	ind=np.where( (chm['Data']>=0.5) & (chm['Data']<2) ); chm1[ind]=2
	ind=np.where( (chm['Data']>=2.0) & (chm['Data']<10) ); chm1[ind]=3
	ind=np.where( (chm['Data']>=10) & (chm['Data']<20) ); chm1[ind]=4
	ind=np.where( (chm['Data']>=20) ); chm1[ind]=5
	chm1[:,0]=6

	lab=['< 0.5m','> 0.5m and < 2m','>2m and < 10m','>10m and < 20m','>20m','','']
	cm=np.vstack( ((0.5,0.5,0.5,1),(0.85,1,0.85,1),(0.25,0.9,0.25,1),(0,0.5,0,1),(0,0.25,0,1),(1,1,1,1)) )
	cm=matplotlib.colors.ListedColormap(cm)

	plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
	#im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm,clim=(0,N_tot)) #
	#ax[0].matshow(chm['Data'],extent=chm['Extent'],cmap=cm,clim=(0,N_tot))
	im=ax[0].matshow(chm1,extent=chm['Extent'],cmap=cm)

	prop=gis.ConvertGeographicToGeojson(r'G:\My Drive\Property\Rainville\Property Boundary.xlsx',r'G:\My Drive\Property\Rainville\Property Boundary.geojson','Polygon')
	prop.plot(ax=ax[0],facecolor='none',edgecolor=[1,1,0],label='Property',linewidth=1,linestyle='-')

	#prop=gis.ConvertGeographicToGeojson(r'G:\My Drive\Property\Cowichan Lake Road Lot 3\Coordinates Geographic.xlsx',r'G:\My Drive\Property\Cowichan Lake Road Lot 3\Coordinates Geographic.geojson','Polygon')
	#prop.plot(ax=ax[0],facecolor='none',edgecolor=[1,0,0],label='Property',linewidth=1.25,linestyle='-')

	#row=gis.ConvertGeographicToGeojson(r'G:\My Drive\Property\Cowichan Lake Road Lot 3\ROW.xlsx',r'G:\My Drive\Property\Cowichan Lake Road Lot 3\ROW.geojson','Line')
	#row.plot(ax=ax[0],facecolor='none',edgecolor=[0.7,0.2,1],label='Property',linewidth=1,linestyle='-')

	#phot=gis.ConvertGeographicToGeojson(r'G:\My Drive\Property\Cowichan Lake Road Lot 3\Photos.xlsx',r'G:\My Drive\Property\Cowichan Lake Road Lot 3\Photos.geojson','Point')
	#phot.plot(ax=ax[0],marker='o',markersize=7,facecolor='none',edgecolor=[1,1,0],label='Property',linewidth=1,linestyle='-')

	#roi['gdf']['bound within'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
	roi['gdf']['bc_bound'].plot(ax=ax[0],edgecolor=[0.6,0.6,0.6],facecolor='none',linewidth=0.25)
	roi['gdf']['lakes'].plot(ax=ax[0],facecolor=[0.82,0.88,1],label='Lakes',linewidth=0.25)
	#roi['gdf']['rivers'].plot(ax=ax[0],color=[0.6,0.8,1],label='Rivers',linewidth=0.25)
	roi['gdf']['road'].plot(ax=ax[0],edgecolor=[0.9,0.45,0],linewidth=0.25,label='Road',alpha=1,zorder=1)
	roi['gdf']['hydrol'].plot(ax=ax[0],linestyle='--',edgecolor=[0.27,0.45,0.8],linewidth=0.5,alpha=1,zorder=1,label='Road')
	roi['gdf']['rail'].plot(ax=ax[0],edgecolor=[0,0,0],linewidth=0.5,alpha=1,zorder=1,label='Rail')
	roi['gdf']['rail'].plot(ax=ax[0],edgecolor=[1,1,1],linewidth=0.25,alpha=1,zorder=1,label='Rail')
	roi['gdf']['tpf'].plot(ax=ax[0],marker='^',edgecolor=[0,0.75,0.75],facecolor=[0,1,1],linewidth=0.25,label='TPF',alpha=1,zorder=1,markersize=5)
	#for x,y,label in zip(roi['gdf']['tpf'].geometry.x,roi['gdf']['tpf'].geometry.y,roi['gdf']['tpf']['COMPANY_NAME']):
	#	ax[0].annotate(label,xy=(x,y),xytext=(4,3),textcoords="offset points",color=[0,0,0],fontsize=4)
	ind=np.where( (roi['gdf']['cities']['Territory']=='BC') & (roi['gdf']['cities']['Level']==1) )[0]
	roi['gdf']['cities'].iloc[ind].plot(ax=ax[0],marker='s',edgecolor=[0.5,0,0],facecolor=[1,0,0],lw=0.5,markersize=7,alpha=1,zorder=2)
	for x,y,label in zip(roi['gdf']['cities'].iloc[ind].geometry.x,roi['gdf']['cities'].iloc[ind].geometry.y,roi['gdf']['cities'].iloc[ind]['City Name']):
		ax[0].annotate(label,xy=(x,y),xytext=(2,1.5),textcoords="offset points",color=[0,0,0],fontsize=4)
	#roi['gdf']['popp'].plot(ax=ax[0],marker='o',edgecolor=[0.75,0.3,0],facecolor=[1,0.6,0],lw=0.25,markersize=6,alpha=1,zorder=2)
	#for x,y,label in zip(roi['gdf']['popp'].geometry.x,roi['gdf']['popp'].geometry.y,roi['gdf']['popp']['NAME']):
	#	ax[0].annotate(label,xy=(x,y),xytext=(2,1.5),textcoords="offset points",color=[0,0,0],fontsize=4)
	
	#ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
	ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=chm['xlim'],ylim=chm['ylim'])
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])
	cb=plt.colorbar(im,cax=ax[1])

	w=N_vis/N_tot
	ticks=np.arange(1+w/2,N_tot+1,w)
	cb=plt.colorbar(im,cax=ax[1],ticks=ticks)
	cb.ax.set(yticklabels=lab)
	cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
	cb.outline.set_edgecolor('w')
	for i,ti in enumerate(ticks):
		ax[1].plot([0,100],[ti+w/2,ti+w/2],'w-',linewidth=1.5)

	pos2=[0.75,0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
	ax[1].set(position=pos2)
		
	#gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_HeightClass_a','png',900)
	#gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_HeightClass_b','png',900)
	gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_HeightClass_Rainville','png',900)
	return fig,ax

#%%
def Plot_zmin(meta,roi,vnam):
	#z0=roi['grd'][vnam]['Data']
	zmin=gis.OpenGeoTiff(r'D:\Data\zmin_p.tif')
	z0=zmin['Data'].copy()

	bw=0.25; bin=np.arange(124,138+bw,bw)
	N_vis=bin.size
	N_hidden=1
	N_tot=N_vis+N_hidden
	
	z1=N_vis*np.ones(z0.shape)
	for i in range(bin.size):
		ind=np.where(np.abs(z0-bin[i])<=bw/2)
		z1[ind]=i
	z1[(z0>bin[i])]=i
	z1[0,0:N_tot]=np.arange(0,N_tot,1)

	lab=bin.astype(str)

	cm=plt.cm.get_cmap('viridis',i)
	cmc=np.zeros((N_vis,4))
	cnt=0
	for ivl in np.linspace(0,1,N_vis):
		cmc[cnt,:]=np.array(cm(ivl))
		cnt=cnt+1
	cm=np.vstack( (cmc,(1,1,1,1)) )
	cm=matplotlib.colors.ListedColormap(cm)

	plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
	#im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm,clim=(0,N_tot)) #
	#ax[0].matshow(chm['Data'],extent=chm['Extent'],cmap=cm,clim=(0,N_tot))
	im=ax[0].matshow(z1,extent=zmin['Extent'],cmap=cm,clim=(0,N_tot))

	prop=gis.ConvertGeographicToGeojson(r'G:\My Drive\Property\Cowichan Lake Road Lot 3\Coordinates Geographic.xlsx',r'G:\My Drive\Property\Cowichan Lake Road Lot 3\Coordinates Geographic.geojson','Polygon')
	prop.plot(ax=ax[0],facecolor='none',edgecolor=[1,0,0],label='Property',linewidth=1.25,linestyle='-')

	row=gis.ConvertGeographicToGeojson(r'G:\My Drive\Property\Cowichan Lake Road Lot 3\ROW.xlsx',r'G:\My Drive\Property\Cowichan Lake Road Lot 3\ROW.geojson','Line')
	row.plot(ax=ax[0],facecolor='none',edgecolor=[0.7,0.2,1],label='Property',linewidth=1,linestyle='-')

	phot=gis.ConvertGeographicToGeojson(r'G:\My Drive\Property\Cowichan Lake Road Lot 3\Photos.xlsx',r'G:\My Drive\Property\Cowichan Lake Road Lot 3\Photos.geojson','Point')
	phot.plot(ax=ax[0],marker='o',markersize=7,facecolor='none',edgecolor=[1,0,0],label='Property',linewidth=1,linestyle='-')

	#roi['gdf']['bound within'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
	roi['gdf']['bc_bound'].plot(ax=ax[0],edgecolor=[0.6,0.6,0.6],facecolor='none',linewidth=0.25)
	roi['gdf']['lakes'].plot(ax=ax[0],facecolor=[0.82,0.88,1],label='Lakes',linewidth=0.25)
	#roi['gdf']['rivers'].plot(ax=ax[0],color=[0.6,0.8,1],label='Rivers',linewidth=0.25)
	roi['gdf']['road'].plot(ax=ax[0],edgecolor=[0.9,0.45,0],linewidth=0.25,label='Road',alpha=1,zorder=1)
	roi['gdf']['hydrol'].plot(ax=ax[0],linestyle='--',edgecolor=[0.27,0.45,0.8],linewidth=0.5,alpha=1,zorder=1,label='Road')
	roi['gdf']['rail'].plot(ax=ax[0],edgecolor=[0,0,0],linewidth=0.5,alpha=1,zorder=1,label='Rail')
	roi['gdf']['rail'].plot(ax=ax[0],edgecolor=[1,1,1],linewidth=0.25,alpha=1,zorder=1,label='Rail')
	roi['gdf']['tpf'].plot(ax=ax[0],marker='^',edgecolor=[0,0.75,0.75],facecolor=[0,1,1],linewidth=0.25,label='TPF',alpha=1,zorder=1,markersize=5)
	#for x,y,label in zip(roi['gdf']['tpf'].geometry.x,roi['gdf']['tpf'].geometry.y,roi['gdf']['tpf']['COMPANY_NAME']):
	#	ax[0].annotate(label,xy=(x,y),xytext=(4,3),textcoords="offset points",color=[0,0,0],fontsize=4)
	ind=np.where( (roi['gdf']['cities']['Territory']=='BC') & (roi['gdf']['cities']['Level']==1) )[0]
	roi['gdf']['cities'].iloc[ind].plot(ax=ax[0],marker='s',edgecolor=[0.5,0,0],facecolor=[1,0,0],lw=0.5,markersize=7,alpha=1,zorder=2)
	for x,y,label in zip(roi['gdf']['cities'].iloc[ind].geometry.x,roi['gdf']['cities'].iloc[ind].geometry.y,roi['gdf']['cities'].iloc[ind]['City Name']):
		ax[0].annotate(label,xy=(x,y),xytext=(2,1.5),textcoords="offset points",color=[0,0,0],fontsize=4)
	#roi['gdf']['popp'].plot(ax=ax[0],marker='o',edgecolor=[0.75,0.3,0],facecolor=[1,0.6,0],lw=0.25,markersize=6,alpha=1,zorder=2)
	#for x,y,label in zip(roi['gdf']['popp'].geometry.x,roi['gdf']['popp'].geometry.y,roi['gdf']['popp']['NAME']):
	#	ax[0].annotate(label,xy=(x,y),xytext=(2,1.5),textcoords="offset points",color=[0,0,0],fontsize=4)

	ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=zmin['xlim'],ylim=zmin['ylim'])
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

	cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_vis+1,1),ticks=np.arange(0.5,N_vis,1))
	cb.ax.set(yticklabels=lab)
	cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
	cb.outline.set_edgecolor('w')
	for i in range(0,N_vis):
		ax[1].plot([0,100],[i,i],'w-',linewidth=1.5)

	pos2=[0.75,0.1,0.02,0.8]
	ax[1].set(position=pos2)
		
	#gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_zmin_b','png',900)
	return fig,ax

#%%
def Plot_InfastructureLumber(meta,roi,vnam):
	z0=roi['grd'][vnam]['Data']
	bw=500; bin=np.arange(0,5000+bw,bw);
	
	N_vis=bin.size
	N_hidden=1
	N_tot=N_vis+N_hidden
	
	z1=N_vis*np.ones(z0.shape)
	for i in range(bin.size):
		ind=np.where(np.abs(z0-bin[i])<=bw/2)
		z1[ind]=i
	z1[(z1>bin[i])]=i
	z1[(roi['grd']['Data']==0)]=i+1
	z1[0,0:N_tot]=np.arange(0,N_tot,1)

	lab=bin.astype(str)

	cm=plt.cm.get_cmap('Greys',i)
	cmc=np.zeros((N_vis,4))
	cnt=0
	for ivl in np.linspace(0,1,N_vis):
		cmc[cnt,:]=np.array(cm(ivl))
		cnt=cnt+1
	cm=np.vstack( (cmc,(1,1,1,1)) )
	cm=matplotlib.colors.ListedColormap(cm)

	plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
	im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm,clim=(0,N_tot)) #
	
	#roi['gdf']['bound within'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
	roi['gdf']['bc_bound'].plot(ax=ax[0],edgecolor=[0.6,0.6,0.6],facecolor='none',linewidth=0.25)
	roi['gdf']['lakes'].plot(ax=ax[0],facecolor=[0.82,0.88,1],label='Lakes',linewidth=0.25)
	#roi['gdf']['rivers'].plot(ax=ax[0],color=[0.6,0.8,1],label='Rivers',linewidth=0.25)

	roi['gdf']['road'].plot(ax=ax[0],edgecolor=[0.9,0.45,0],linewidth=0.25,label='Road',alpha=1,zorder=1)

	roi['gdf']['hydrol'].plot(ax=ax[0],linestyle='--',edgecolor=[0.27,0.45,0.8],linewidth=0.5,alpha=1,zorder=1,label='Road')

	roi['gdf']['rail'].plot(ax=ax[0],edgecolor=[0,0,0],linewidth=0.5,alpha=1,zorder=1,label='Rail')
	roi['gdf']['rail'].plot(ax=ax[0],edgecolor=[1,1,1],linewidth=0.25,alpha=1,zorder=1,label='Rail')

	#mtypeL=['LBR','PLP','PLT']
	ind=np.where( (roi['gdf']['tpf']['PRODUCT_CODE']=='LBR') )[0]
	y=roi['gdf']['tpf'].iloc[ind]['EST_AN_CAP_MLN_BOARD_FT']/453
	ms=50+500*y
	roi['gdf']['tpf'].iloc[ind].plot(ax=ax[0],marker='o',edgecolor='g',facecolor='g',lw=0.25,markersize=ms,alpha=0.35,zorder=2)
	for x,y,label in zip(roi['gdf']['tpf'].iloc[ind].geometry.x,roi['gdf']['tpf'].iloc[ind].geometry.y,roi['gdf']['tpf'].iloc[ind]['COMPANY_NAME']):
		ax[0].annotate(label,xy=(x,y),xytext=(5,5),textcoords="offset points",color='g',fontsize=3)

	ind=np.where( (roi['gdf']['cities']['Territory']=='BC') & (roi['gdf']['cities']['Level']==1) )[0]
	roi['gdf']['cities'].iloc[ind].plot(ax=ax[0],marker='s',edgecolor=[0.5,0,0],facecolor=[1,0,0],lw=0.5,markersize=7,alpha=1,zorder=2)
	for x,y,label in zip(roi['gdf']['cities'].iloc[ind].geometry.x,roi['gdf']['cities'].iloc[ind].geometry.y,roi['gdf']['cities'].iloc[ind]['City Name']):
		ax[0].annotate(label,xy=(x,y),xytext=(2,1.5),textcoords="offset points",color=[0,0,0],fontsize=4)

	#roi['gdf']['popp'].plot(ax=ax[0],marker='o',edgecolor=[0.75,0.3,0],facecolor=[1,0.6,0],lw=0.25,markersize=6,alpha=1,zorder=2)
	#for x,y,label in zip(roi['gdf']['popp'].geometry.x,roi['gdf']['popp'].geometry.y,roi['gdf']['popp']['NAME']):
	#	ax[0].annotate(label,xy=(x,y),xytext=(2,1.5),textcoords="offset points",color=[0,0,0],fontsize=4)
	
	ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])
	
	cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_vis+1,1),ticks=np.arange(0.5,N_vis,1))
	cb.ax.set(yticklabels=lab)
	cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
	cb.outline.set_edgecolor('w')
	for i in range(0,N_vis):
		ax[1].plot([0,100],[i,i],'w-',linewidth=1.5)
	pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
	ax[1].set(position=pos2)
		
	gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_infastructureLumber','png',900)
	return fig,ax

#%%
def Plot_InfastructurePulp(meta,roi,vnam):
	z0=roi['grd'][vnam]['Data']
	bw=500; bin=np.arange(0,5000+bw,bw);
	
	N_vis=bin.size
	N_hidden=1
	N_tot=N_vis+N_hidden
	
	z1=N_vis*np.ones(z0.shape)
	for i in range(bin.size):
		ind=np.where(np.abs(z0-bin[i])<=bw/2)
		z1[ind]=i
	z1[(z1>bin[i])]=i
	z1[(roi['grd']['Data']==0)]=i+1
	z1[0,0:N_tot]=np.arange(0,N_tot,1)

	lab=bin.astype(str)

	cm=plt.cm.get_cmap('Greys',i)
	cmc=np.zeros((N_vis,4))
	cnt=0
	for ivl in np.linspace(0,1,N_vis):
		cmc[cnt,:]=np.array(cm(ivl))
		cnt=cnt+1
	cm=np.vstack( (cmc,(1,1,1,1)) )
	cm=matplotlib.colors.ListedColormap(cm)

	plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
	im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm,clim=(0,N_tot)) #
	
	#roi['gdf']['bound within'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
	roi['gdf']['bc_bound'].plot(ax=ax[0],edgecolor=[0.6,0.6,0.6],facecolor='none',linewidth=0.25)
	roi['gdf']['lakes'].plot(ax=ax[0],facecolor=[0.82,0.88,1],label='Lakes',linewidth=0.25)
	#roi['gdf']['rivers'].plot(ax=ax[0],color=[0.6,0.8,1],label='Rivers',linewidth=0.25)

	roi['gdf']['road'].plot(ax=ax[0],edgecolor=[0.9,0.45,0],linewidth=0.25,label='Road',alpha=1,zorder=1)

	roi['gdf']['hydrol'].plot(ax=ax[0],linestyle='--',edgecolor=[0.27,0.45,0.8],linewidth=0.5,alpha=1,zorder=1,label='Road')

	roi['gdf']['rail'].plot(ax=ax[0],edgecolor=[0,0,0],linewidth=0.5,alpha=1,zorder=1,label='Rail')
	roi['gdf']['rail'].plot(ax=ax[0],edgecolor=[1,1,1],linewidth=0.25,alpha=1,zorder=1,label='Rail')

	#mtypeL=['LBR','PLP','PLT']
	ind=np.where( (roi['gdf']['tpf']['PRODUCT_CODE']=='PLP') )[0]
	y=roi['gdf']['tpf'].iloc[ind]['EST_AN_CAP_000_TONNES']
	ms=0.05+0.75*2*y
	roi['gdf']['tpf'].iloc[ind].plot(ax=ax[0],marker='o',edgecolor='g',facecolor='g',lw=0.25,markersize=ms,alpha=0.35,zorder=2)
	for x,y,label in zip(roi['gdf']['tpf'].iloc[ind].geometry.x,roi['gdf']['tpf'].iloc[ind].geometry.y,roi['gdf']['tpf'].iloc[ind]['COMPANY_NAME']):
		ax[0].annotate(label,xy=(x,y),xytext=(5,5),textcoords="offset points",color='g',fontsize=3)

	ind=np.where( (roi['gdf']['cities']['Territory']=='BC') & (roi['gdf']['cities']['Level']==1) )[0]
	roi['gdf']['cities'].iloc[ind].plot(ax=ax[0],marker='s',edgecolor=[0.5,0,0],facecolor=[1,0,0],lw=0.5,markersize=7,alpha=1,zorder=2)
	for x,y,label in zip(roi['gdf']['cities'].iloc[ind].geometry.x,roi['gdf']['cities'].iloc[ind].geometry.y,roi['gdf']['cities'].iloc[ind]['City Name']):
		ax[0].annotate(label,xy=(x,y),xytext=(2,1.5),textcoords="offset points",color=[0,0,0],fontsize=4)

	#roi['gdf']['popp'].plot(ax=ax[0],marker='o',edgecolor=[0.75,0.3,0],facecolor=[1,0.6,0],lw=0.25,markersize=6,alpha=1,zorder=2)
	#for x,y,label in zip(roi['gdf']['popp'].geometry.x,roi['gdf']['popp'].geometry.y,roi['gdf']['popp']['NAME']):
	#	ax[0].annotate(label,xy=(x,y),xytext=(2,1.5),textcoords="offset points",color=[0,0,0],fontsize=4)
	
	ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])
	
	cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_vis+1,1),ticks=np.arange(0.5,N_vis,1))
	cb.ax.set(yticklabels=lab)
	cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
	cb.outline.set_edgecolor('w')
	for i in range(0,N_vis):
		ax[1].plot([0,100],[i,i],'w-',linewidth=1.5)
	pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
	ax[1].set(position=pos2)
		
	gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_infastructurePulp','png',900)
	return fig,ax

#%%
def Plot_InfastructurePanel(meta,roi,vnam):
	z0=roi['grd'][vnam]['Data']
	bw=500; bin=np.arange(0,5000+bw,bw);
	
	N_vis=bin.size
	N_hidden=1
	N_tot=N_vis+N_hidden
	
	z1=N_vis*np.ones(z0.shape)
	for i in range(bin.size):
		ind=np.where(np.abs(z0-bin[i])<=bw/2)
		z1[ind]=i
	z1[(z1>bin[i])]=i
	z1[(roi['grd']['Data']==0)]=i+1
	z1[0,0:N_tot]=np.arange(0,N_tot,1)

	lab=bin.astype(str)

	cm=plt.cm.get_cmap('Greys',i)
	cmc=np.zeros((N_vis,4))
	cnt=0
	for ivl in np.linspace(0,1,N_vis):
		cmc[cnt,:]=np.array(cm(ivl))
		cnt=cnt+1
	cm=np.vstack( (cmc,(1,1,1,1)) )
	cm=matplotlib.colors.ListedColormap(cm)

	plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
	im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm,clim=(0,N_tot)) #
	
	#roi['gdf']['bound within'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
	roi['gdf']['bc_bound'].plot(ax=ax[0],edgecolor=[0.6,0.6,0.6],facecolor='none',linewidth=0.25)
	roi['gdf']['lakes'].plot(ax=ax[0],facecolor=[0.82,0.88,1],label='Lakes',linewidth=0.25)
	#roi['gdf']['rivers'].plot(ax=ax[0],color=[0.6,0.8,1],label='Rivers',linewidth=0.25)

	roi['gdf']['road'].plot(ax=ax[0],edgecolor=[0.9,0.45,0],linewidth=0.25,label='Road',alpha=1,zorder=1)

	roi['gdf']['hydrol'].plot(ax=ax[0],linestyle='--',edgecolor=[0.27,0.45,0.8],linewidth=0.5,alpha=1,zorder=1,label='Road')

	roi['gdf']['rail'].plot(ax=ax[0],edgecolor=[0,0,0],linewidth=0.5,alpha=1,zorder=1,label='Rail')
	roi['gdf']['rail'].plot(ax=ax[0],edgecolor=[1,1,1],linewidth=0.25,alpha=1,zorder=1,label='Rail')

	ind=np.where( (np.isin(roi['gdf']['tpf']['PRODUCT_CODE'],['PLY','VNR','OSB','PNL'])==True) )[0]
	y=roi['gdf']['tpf'].iloc[ind]['EST_AN_CAP_MLN_SQ_FT']
	ms=10+1000*(y/885)
	roi['gdf']['tpf'].iloc[ind].plot(ax=ax[0],marker='o',edgecolor='g',facecolor='g',lw=0.25,markersize=ms,alpha=0.35,zorder=2)
	for x,y,label in zip(roi['gdf']['tpf'].iloc[ind].geometry.x,roi['gdf']['tpf'].iloc[ind].geometry.y,roi['gdf']['tpf'].iloc[ind]['COMPANY_NAME']):
		ax[0].annotate(label,xy=(x,y),xytext=(5,5),textcoords="offset points",color='g',fontsize=3)

	ind=np.where( (roi['gdf']['cities']['Territory']=='BC') & (roi['gdf']['cities']['Level']==1) )[0]
	roi['gdf']['cities'].iloc[ind].plot(ax=ax[0],marker='s',edgecolor=[0.5,0,0],facecolor=[1,0,0],lw=0.5,markersize=7,alpha=1,zorder=2)
	for x,y,label in zip(roi['gdf']['cities'].iloc[ind].geometry.x,roi['gdf']['cities'].iloc[ind].geometry.y,roi['gdf']['cities'].iloc[ind]['City Name']):
		ax[0].annotate(label,xy=(x,y),xytext=(2,1.5),textcoords="offset points",color=[0,0,0],fontsize=4)

	#roi['gdf']['popp'].plot(ax=ax[0],marker='o',edgecolor=[0.75,0.3,0],facecolor=[1,0.6,0],lw=0.25,markersize=6,alpha=1,zorder=2)
	#for x,y,label in zip(roi['gdf']['popp'].geometry.x,roi['gdf']['popp'].geometry.y,roi['gdf']['popp']['NAME']):
	#	ax[0].annotate(label,xy=(x,y),xytext=(2,1.5),textcoords="offset points",color=[0,0,0],fontsize=4)
	
	ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])
	
	cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_vis+1,1),ticks=np.arange(0.5,N_vis,1))
	cb.ax.set(yticklabels=lab)
	cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
	cb.outline.set_edgecolor('w')
	for i in range(0,N_vis):
		ax[1].plot([0,100],[i,i],'w-',linewidth=1.5)
	pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
	ax[1].set(position=pos2)
		
	gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_infastructurePanels','png',900)
	return fig,ax

#%%
def Plot_InfastructurePellet(meta,roi,vnam):
	z0=roi['grd'][vnam]['Data']
	bw=500; bin=np.arange(0,5000+bw,bw);
	
	N_vis=bin.size
	N_hidden=1
	N_tot=N_vis+N_hidden
	
	z1=N_vis*np.ones(z0.shape)
	for i in range(bin.size):
		ind=np.where(np.abs(z0-bin[i])<=bw/2)
		z1[ind]=i
	z1[(z1>bin[i])]=i
	z1[(roi['grd']['Data']==0)]=i+1
	z1[0,0:N_tot]=np.arange(0,N_tot,1)

	lab=bin.astype(str)

	cm=plt.cm.get_cmap('Greys',i)
	cmc=np.zeros((N_vis,4))
	cnt=0
	for ivl in np.linspace(0,1,N_vis):
		cmc[cnt,:]=np.array(cm(ivl))
		cnt=cnt+1
	cm=np.vstack( (cmc,(1,1,1,1)) )
	cm=matplotlib.colors.ListedColormap(cm)

	plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
	im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm,clim=(0,N_tot)) #
	
	#roi['gdf']['bound within'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
	roi['gdf']['bc_bound'].plot(ax=ax[0],edgecolor=[0.6,0.6,0.6],facecolor='none',linewidth=0.25)
	roi['gdf']['lakes'].plot(ax=ax[0],facecolor=[0.82,0.88,1],label='Lakes',linewidth=0.25)
	#roi['gdf']['rivers'].plot(ax=ax[0],color=[0.6,0.8,1],label='Rivers',linewidth=0.25)

	roi['gdf']['road'].plot(ax=ax[0],edgecolor=[0.9,0.45,0],linewidth=0.25,label='Road',alpha=1,zorder=1)

	roi['gdf']['hydrol'].plot(ax=ax[0],linestyle='--',edgecolor=[0.27,0.45,0.8],linewidth=0.5,alpha=1,zorder=1,label='Road')

	roi['gdf']['rail'].plot(ax=ax[0],edgecolor=[0,0,0],linewidth=0.5,alpha=1,zorder=1,label='Rail')
	roi['gdf']['rail'].plot(ax=ax[0],edgecolor=[1,1,1],linewidth=0.25,alpha=1,zorder=1,label='Rail')

	ind=np.where( (np.isin(roi['gdf']['tpf']['PRODUCT_CODE'],['PLT'])==True) )[0]
	y=roi['gdf']['tpf'].iloc[ind]['EST_AN_CAP_000_TONNES']
	ms=1+2*y
	roi['gdf']['tpf'].iloc[ind].plot(ax=ax[0],marker='o',edgecolor='g',facecolor='g',lw=0.25,markersize=ms,alpha=0.35,zorder=2)
	for x,y,label in zip(roi['gdf']['tpf'].iloc[ind].geometry.x,roi['gdf']['tpf'].iloc[ind].geometry.y,roi['gdf']['tpf'].iloc[ind]['COMPANY_NAME']):
		ax[0].annotate(label,xy=(x,y),xytext=(5,5),textcoords="offset points",color='g',fontsize=3)

	ind=np.where( (roi['gdf']['cities']['Territory']=='BC') & (roi['gdf']['cities']['Level']==1) )[0]
	roi['gdf']['cities'].iloc[ind].plot(ax=ax[0],marker='s',edgecolor=[0.5,0,0],facecolor=[1,0,0],lw=0.5,markersize=7,alpha=1,zorder=2)
	for x,y,label in zip(roi['gdf']['cities'].iloc[ind].geometry.x,roi['gdf']['cities'].iloc[ind].geometry.y,roi['gdf']['cities'].iloc[ind]['City Name']):
		ax[0].annotate(label,xy=(x,y),xytext=(2,1.5),textcoords="offset points",color=[0,0,0],fontsize=4)

	#roi['gdf']['popp'].plot(ax=ax[0],marker='o',edgecolor=[0.75,0.3,0],facecolor=[1,0.6,0],lw=0.25,markersize=6,alpha=1,zorder=2)
	#for x,y,label in zip(roi['gdf']['popp'].geometry.x,roi['gdf']['popp'].geometry.y,roi['gdf']['popp']['NAME']):
	#	ax[0].annotate(label,xy=(x,y),xytext=(2,1.5),textcoords="offset points",color=[0,0,0],fontsize=4)
	
	ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])
	
	cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_vis+1,1),ticks=np.arange(0.5,N_vis,1))
	cb.ax.set(yticklabels=lab)
	cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
	cb.outline.set_edgecolor('w')
	for i in range(0,N_vis):
		ax[1].plot([0,100],[i,i],'w-',linewidth=1.5)
	pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
	ax[1].set(position=pos2)
		
	gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_infastructurePellet','png',900)
	return fig,ax

#%%
def Plot_InfastructureChipper(meta,roi,vnam):
	z0=roi['grd'][vnam]['Data']
	bw=500; bin=np.arange(0,5000+bw,bw);
	
	N_vis=bin.size
	N_hidden=1
	N_tot=N_vis+N_hidden
	
	z1=N_vis*np.ones(z0.shape)
	for i in range(bin.size):
		ind=np.where(np.abs(z0-bin[i])<=bw/2)
		z1[ind]=i
	z1[(z1>bin[i])]=i
	z1[(roi['grd']['Data']==0)]=i+1
	z1[0,0:N_tot]=np.arange(0,N_tot,1)

	lab=bin.astype(str)

	cm=plt.cm.get_cmap('Greys',i)
	cmc=np.zeros((N_vis,4))
	cnt=0
	for ivl in np.linspace(0,1,N_vis):
		cmc[cnt,:]=np.array(cm(ivl))
		cnt=cnt+1
	cm=np.vstack( (cmc,(1,1,1,1)) )
	cm=matplotlib.colors.ListedColormap(cm)

	plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
	im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm,clim=(0,N_tot)) #
	
	#roi['gdf']['bound within'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
	roi['gdf']['bc_bound'].plot(ax=ax[0],edgecolor=[0.6,0.6,0.6],facecolor='none',linewidth=0.25)
	roi['gdf']['lakes'].plot(ax=ax[0],facecolor=[0.82,0.88,1],label='Lakes',linewidth=0.25)
	#roi['gdf']['rivers'].plot(ax=ax[0],color=[0.6,0.8,1],label='Rivers',linewidth=0.25)

	roi['gdf']['road'].plot(ax=ax[0],edgecolor=[0.9,0.45,0],linewidth=0.25,label='Road',alpha=1,zorder=1)

	roi['gdf']['hydrol'].plot(ax=ax[0],linestyle='--',edgecolor=[0.27,0.45,0.8],linewidth=0.5,alpha=1,zorder=1,label='Road')

	roi['gdf']['rail'].plot(ax=ax[0],edgecolor=[0,0,0],linewidth=0.5,alpha=1,zorder=1,label='Rail')
	roi['gdf']['rail'].plot(ax=ax[0],edgecolor=[1,1,1],linewidth=0.25,alpha=1,zorder=1,label='Rail')

	ind=np.where( (np.isin(roi['gdf']['tpf']['PRODUCT_CODE'],['CHP'])==True) )[0]
	y=roi['gdf']['tpf'].iloc[ind]['EST_AN_CAP_000_BDUS']
	ms=1+0.7*2*y
	roi['gdf']['tpf'].iloc[ind].plot(ax=ax[0],marker='o',edgecolor='g',facecolor='g',lw=0.25,markersize=ms,alpha=0.35,zorder=2)
	for x,y,label in zip(roi['gdf']['tpf'].iloc[ind].geometry.x,roi['gdf']['tpf'].iloc[ind].geometry.y,roi['gdf']['tpf'].iloc[ind]['COMPANY_NAME']):
		ax[0].annotate(label,xy=(x,y),xytext=(5,5),textcoords="offset points",color='g',fontsize=3)

	ind=np.where( (roi['gdf']['cities']['Territory']=='BC') & (roi['gdf']['cities']['Level']==1) )[0]
	roi['gdf']['cities'].iloc[ind].plot(ax=ax[0],marker='s',edgecolor=[0.5,0,0],facecolor=[1,0,0],lw=0.5,markersize=7,alpha=1,zorder=2)
	for x,y,label in zip(roi['gdf']['cities'].iloc[ind].geometry.x,roi['gdf']['cities'].iloc[ind].geometry.y,roi['gdf']['cities'].iloc[ind]['City Name']):
		ax[0].annotate(label,xy=(x,y),xytext=(2,1.5),textcoords="offset points",color=[0,0,0],fontsize=4)

	#roi['gdf']['popp'].plot(ax=ax[0],marker='o',edgecolor=[0.75,0.3,0],facecolor=[1,0.6,0],lw=0.25,markersize=6,alpha=1,zorder=2)
	#for x,y,label in zip(roi['gdf']['popp'].geometry.x,roi['gdf']['popp'].geometry.y,roi['gdf']['popp']['NAME']):
	#	ax[0].annotate(label,xy=(x,y),xytext=(2,1.5),textcoords="offset points",color=[0,0,0],fontsize=4)
	
	ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])
	
	cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_vis+1,1),ticks=np.arange(0.5,N_vis,1))
	cb.ax.set(yticklabels=lab)
	cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
	cb.outline.set_edgecolor('w')
	for i in range(0,N_vis):
		ax[1].plot([0,100],[i,i],'w-',linewidth=1.5)
	pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
	ax[1].set(position=pos2)
		
	gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_infastructureChipper','png',900)
	return fig,ax

#%%
def Plot_InfastructureBioenergy(meta,roi,vnam):
	z0=roi['grd'][vnam]['Data']
	bw=500; bin=np.arange(0,5000+bw,bw);
	
	N_vis=bin.size
	N_hidden=1
	N_tot=N_vis+N_hidden
	
	z1=N_vis*np.ones(z0.shape)
	for i in range(bin.size):
		ind=np.where(np.abs(z0-bin[i])<=bw/2)
		z1[ind]=i
	z1[(z1>bin[i])]=i
	z1[(roi['grd']['Data']==0)]=i+1
	z1[0,0:N_tot]=np.arange(0,N_tot,1)

	lab=bin.astype(str)

	cm=plt.cm.get_cmap('Greys',i)
	cmc=np.zeros((N_vis,4))
	cnt=0
	for ivl in np.linspace(0,1,N_vis):
		cmc[cnt,:]=np.array(cm(ivl))
		cnt=cnt+1
	cm=np.vstack( (cmc,(1,1,1,1)) )
	cm=matplotlib.colors.ListedColormap(cm)

	plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
	im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm,clim=(0,N_tot)) #
	
	#roi['gdf']['bound within'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
	roi['gdf']['bc_bound'].plot(ax=ax[0],edgecolor=[0.6,0.6,0.6],facecolor='none',linewidth=0.25)
	roi['gdf']['lakes'].plot(ax=ax[0],facecolor=[0.82,0.88,1],label='Lakes',linewidth=0.25)
	#roi['gdf']['rivers'].plot(ax=ax[0],color=[0.6,0.8,1],label='Rivers',linewidth=0.25)

	roi['gdf']['road'].plot(ax=ax[0],edgecolor=[0.9,0.45,0],linewidth=0.25,label='Road',alpha=1,zorder=1)

	roi['gdf']['hydrol'].plot(ax=ax[0],linestyle='--',edgecolor=[0.27,0.45,0.8],linewidth=0.5,alpha=1,zorder=1,label='Road')

	roi['gdf']['rail'].plot(ax=ax[0],edgecolor=[0,0,0],linewidth=0.5,alpha=1,zorder=1,label='Rail')
	roi['gdf']['rail'].plot(ax=ax[0],edgecolor=[1,1,1],linewidth=0.25,alpha=1,zorder=1,label='Rail')

	#ind=np.where( (np.isin(roi['gdf']['bioe']['PRODUCT_CODE'],['CHP'])==True) )[0]
	y=roi['gdf']['bioe']['Capacity Electrical (MW)']
	ms=5+1.5*y
	roi['gdf']['bioe'].plot(ax=ax[0],marker='o',edgecolor='g',facecolor='g',lw=0.25,markersize=ms,alpha=0.35,zorder=2)
	for x,y,label in zip(roi['gdf']['bioe'].geometry.x,roi['gdf']['bioe'].geometry.y,roi['gdf']['bioe']['Name']):
		ax[0].annotate(label,xy=(x,y),xytext=(5,5),textcoords="offset points",color='g',fontsize=3)

	ind=np.where( (roi['gdf']['cities']['Territory']=='BC') & (roi['gdf']['cities']['Level']==1) )[0]
	roi['gdf']['cities'].iloc[ind].plot(ax=ax[0],marker='s',edgecolor=[0.5,0,0],facecolor=[1,0,0],lw=0.5,markersize=7,alpha=1,zorder=2)
	for x,y,label in zip(roi['gdf']['cities'].iloc[ind].geometry.x,roi['gdf']['cities'].iloc[ind].geometry.y,roi['gdf']['cities'].iloc[ind]['City Name']):
		ax[0].annotate(label,xy=(x,y),xytext=(2,1.5),textcoords="offset points",color=[0,0,0],fontsize=4)

	ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])
	
	cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_vis+1,1),ticks=np.arange(0.5,N_vis,1))
	cb.ax.set(yticklabels=lab)
	cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
	cb.outline.set_edgecolor('w')
	for i in range(0,N_vis):
		ax[1].plot([0,100],[i,i],'w-',linewidth=1.5)
	pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
	ax[1].set(position=pos2)
		
	gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_infastructureBioenergy','png',900)
	return fig,ax

#%%
def Plot_WaterManagement(meta,roi,vnam):
	z0=roi['grd'][vnam]['Data']
	bw=200; bin=np.arange(600,3000+bw,bw);
	
	N_vis=bin.size
	N_hidden=1
	N_tot=N_vis+N_hidden
	
	z1=N_vis*np.ones(z0.shape)
	for i in range(bin.size):
		ind=np.where(np.abs(z0-bin[i])<=bw/2)
		z1[ind]=i
	z1[(z1>bin[i])]=i
	z1[(roi['grd']['Data']==0)]=i+1
	z1[0,0:N_tot]=np.arange(0,N_tot,1)

	lab=bin.astype(str)

	cm=plt.cm.get_cmap('Greys',i)
	cmc=np.zeros((N_vis,4))
	cnt=0
	for ivl in np.linspace(0,1,N_vis):
		cmc[cnt,:]=np.array(cm(ivl))
		cnt=cnt+1
	cm=np.vstack( (cmc,(1,1,1,1)) )
	cm=matplotlib.colors.ListedColormap(cm)

	if meta['Graphics']['Map']['Show Inset Map']=='On':
		N_panel=3
	else:
		N_panel=2

	plt.close('all');
	fig,ax=plt.subplots(1,N_panel,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
	im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm,clim=(0,N_tot)) #

	#roi['gdf']['bound within'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
	roi['gdf']['bc_bound'].plot(ax=ax[0],edgecolor=[0.6,0.6,0.6],facecolor='none',linewidth=0.25)
	roi['gdf']['lakes'].plot(ax=ax[0],facecolor=[0.5,0.8,1],label='Lakes',linewidth=0.25)
	roi['gdf']['rivers'].plot(ax=ax[0],color=[0.5,0.8,1],label='Rivers',linewidth=0.25)
	roi['gdf']['riversecond'].plot(ax=ax[0],facecolor=[0.5,0.8,1],edgecolor=[0.5,0.8,1],label='Rivers',linewidth=0.5)
	roi['gdf']['rivermajor'].plot(ax=ax[0],facecolor=[0.5,0.8,1],edgecolor=[0.5,0.8,1],label='Rivers',linewidth=0.75)

	#roi['gdf']['road'].plot(ax=ax[0],edgecolor=[0.9,0.45,0],linewidth=0.25,label='Road',alpha=1,zorder=1)
	#roi['gdf']['rail'].plot(ax=ax[0],edgecolor=[0.5,0.2,0],linewidth=0.5,alpha=1,zorder=1,label='Rail')
	#roi['gdf']['hydrol'].plot(ax=ax[0],linestyle='--',edgecolor=[0,0,0.8],linewidth=0.5,alpha=1,zorder=1,label='Road')
	#roi['gdf']['tpf'].plot(ax=ax[0],marker='^',edgecolor=[0,0.75,0.75],facecolor=[0,1,1],linewidth=0.25,label='TPF',alpha=1,zorder=1,markersize=5)

	# Add cities
	ind=np.where( (roi['gdf']['cities']['Territory']=='BC') & (roi['gdf']['cities']['Level']==1) )[0]
	roi['gdf']['cities'].iloc[ind].plot(ax=ax[0],marker='s',edgecolor=[0.5,0,0],facecolor=[1,0,0],lw=0.5,markersize=7,alpha=1,zorder=2)
	for x,y,label in zip(roi['gdf']['cities'].iloc[ind].geometry.x,roi['gdf']['cities'].iloc[ind].geometry.y,roi['gdf']['cities'].iloc[ind]['City Name']):
		ax[0].annotate(label,xy=(x,y),xytext=(2,1.5),textcoords="offset points",color=[0,0,0],fontsize=4)

	# Add wetlands
	if meta['Graphics']['Vector Import']['Wetland']=='On':
		roi['gdf']['wetland'].plot(ax=ax[0],facecolor=[0.7,0.6,1],edgecolor=None,linewidth=0.25,label='Wetlands')

	# Add water management features
	if meta['Graphics']['Vector Import']['Water Management']=='On':
		roi['gdf']['wlic'].plot(ax=ax[0],facecolor=[1,0.65,0.25],edgecolor=[1,0.65,0.25],label='Anthropogenic water bodies',linewidth=0.25,alpha=0.5)
		roi['gdf']['mmwb'].plot(ax=ax[0],facecolor=[0.5,0.9,0],edgecolor=[0.25,0.45,0],label='Anthropogenic water bodies',linewidth=0.25,alpha=0.5)
		roi['gdf']['flooda'].plot(ax=ax[0],facecolor=[0.25,1,1],edgecolor=[0,0.85,.85],linewidth=0.25,label='Flooded areas',alpha=0.25)
		roi['gdf']['dams'].plot(ax=ax[0],facecolor=[0,0,1],edgecolor=[0,0,1],linewidth=0.5,label='Dams')

		# Add streamflow observation sites
		flg=1
		if flg==1:
			roi['gdf']['hydat'].plot(ax=ax[0],marker='^',edgecolor=[0.05,0.25,0.5],facecolor=[0.05,0.25,0.5],linewidth=0.25,markersize=3)
			for x,y,label in zip(roi['gdf']['hydat'].geometry.x,roi['gdf']['hydat'].geometry.y,roi['gdf']['hydat']['Station Name']):
				ax[0].annotate(label,xy=(x,y),xytext=(2,1.5),textcoords="offset points",color=[0.05,0.25,0.5],fontsize=4)

	# Add watersheds
	flg=1
	if flg==1:
		A_th=1e6
		wL=[5,6,7]
		for w in wL:
			nam='wshed' + str(w)
			roi=u1ha.Import_GDB_Over_ROI(meta,roi,[nam])
			roi['gdf'][nam]['gdf']=roi['gdf'][nam]['gdf'].explode()
			iBig=np.where(roi['gdf'][nam]['gdf'].geometry.area>A_th)[0]
			roi['gdf'][nam]['gdf'].iloc[iBig].plot(ax=ax[0],edgecolor=[0.2,0.2,0.2],facecolor='none',linewidth=0.25)
	
			# Add labels
			flg2=1
			if flg2==1:
				roi['gdf'][nam]['gdf']['coords']=roi['gdf'][nam]['gdf']['geometry'].apply(lambda x: x.representative_point().coords[:])
				roi['gdf'][nam]['gdf']['coords']=[coords[0] for coords in roi['gdf'][nam]['gdf']['coords']]
				for idx, row in roi['gdf'][nam]['gdf'].iterrows():
					if row.geometry.area>A_th:
						#print(row['GNIS_NAME'] + ' ' + str(row.geometry.area))
						ax[0].annotate(row['GNIS_NAME'],xy=row['coords'],horizontalalignment='center',verticalalignment='center',color=[0,0,0],fontsize=3)
						#ax[0].annotate(row['GNIS_ID'],xy=row['coords'],horizontalalignment='center',verticalalignment='center',color=[0,0,1])

	ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])
	
	cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_vis+1,1),ticks=np.arange(0.5,N_vis,1))
	cb.ax.set(yticklabels=lab)
	cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
	cb.outline.set_edgecolor('w')
	for i in range(0,N_vis):
		ax[1].plot([0,100],[i,i],'w-',linewidth=1.5)
	pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
	ax[1].set(position=pos2)

	if meta['Graphics']['Map']['Show Inset Map']=='On':
		meta['gdf']['bc_bound']['gdf'].plot(ax=ax[2],edgecolor='k',facecolor=[0.8,0.8,0.8],linewidth=0.25)
		minx, miny, maxx, maxy = meta['gdf']['bc_bound']['gdf'].geometry.total_bounds
		roi['gdf']['bound within'].plot(ax=ax[2],edgecolor=[0,0,0],facecolor=[0,0,0],linewidth=0.25)
		ax[2].set(position=[meta['Graphics']['Map']['Legend X'],0.1,1-meta['Graphics']['Map']['Legend X']-0.01,1-meta['Graphics']['Map']['Legend X']+0.05],xlim=[minx,maxx],ylim=[miny,maxy])
		ax[2].yaxis.set_ticks_position('both'); ax[2].xaxis.set_ticks_position('both'); ax[2].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[2].axis(meta['Graphics']['Map']['Map Axis Vis'])

	gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_watermanagement','png',900)
	return fig,ax

#%%
def PlotWatershedBoundaries(meta,roi,wL,A_th):
	#A_th=1e6
	#wL=[5,6,7]
	#roi=u1ha.Import_Raster(meta,roi,['elev'])

	fig,ax=Plot_Elev(meta,roi,'elev')

	for w in wL:
		nam='wshed' + str(w)
		roi=u1ha.Import_GDB_Over_ROI(meta,roi,[nam])
		roi['gdf'][nam]['gdf']=roi['gdf'][nam]['gdf'].explode()
		iBig=np.where(roi['gdf'][nam]['gdf'].geometry.area>A_th)[0]
		roi['gdf'][nam]['gdf'].iloc[iBig].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)

		# Add labels
		flg=1
		if flg==1:
			roi['gdf'][nam]['gdf']['coords']=roi['gdf'][nam]['gdf']['geometry'].apply(lambda x: x.representative_point().coords[:])
			roi['gdf'][nam]['gdf']['coords']=[coords[0] for coords in roi['gdf'][nam]['gdf']['coords']]
			for idx, row in roi['gdf'][nam]['gdf'].iterrows():
				if row.geometry.area>A_th:
					#print(row['GNIS_NAME'] + ' ' + str(row.geometry.area))
					ax[0].annotate(row['GNIS_NAME'],xy=row['coords'],horizontalalignment='center',verticalalignment='center',color=[0,0,0],fontsize=3)
					#ax[0].annotate(row['GNIS_ID'],xy=row['coords'],horizontalalignment='center',verticalalignment='center',color=[0,0,1])

	# Add streamflow observation sites
	roi['gdf']['streamflow_grdc'].plot(ax=ax[0],marker='^',edgecolor=[0,0,0.5],facecolor=[0,0,1],linewidth=0.25)

	# Save
	gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_watersheds','png',900)

	# Select watersheds
	flg=0
	if flg==1:
		ind=np.where(roi['gdf']['wshed']['gdf']['GNIS_NAME']=='Spius Creek')[0]
		ind=np.where(roi['gdf']['wshed']['gdf']['GNIS_NAME']=='Salmon River')[0]
		roi['gdf']['wshed']['gdf'].iloc[ind].plot(ax=ax[0],edgecolor=[0,0,1],facecolor='none',linewidth=1)
		roi['gdf']['rivers'].plot(ax=ax[0],color=[0.4,0.8,1],label='Rivers',linewidth=1)

	return fig,ax

#%%
def Plot_Settlements(meta,roi,vnam):
	z0=roi['grd'][vnam]['Data']
	bw=200; bin=np.arange(600,3000+bw,bw);
	
	N_vis=bin.size
	N_hidden=1
	N_tot=N_vis+N_hidden
	
	z1=N_vis*np.ones(z0.shape)
	for i in range(bin.size):
		ind=np.where(np.abs(z0-bin[i])<=bw/2)
		z1[ind]=i
	z1[(z1>bin[i])]=i
	z1[(roi['grd']['Data']==0)]=i+1
	z1[0,0:N_tot]=np.arange(0,N_tot,1)

	lab=bin.astype(str)

	cm=plt.cm.get_cmap('Greys',i)
	cmc=np.zeros((N_vis,4))
	cnt=0
	for ivl in np.linspace(0,1,N_vis):
		cmc[cnt,:]=np.array(cm(ivl))
		cnt=cnt+1
	cm=np.vstack( (cmc,(1,1,1,1)) )
	cm=matplotlib.colors.ListedColormap(cm)

	if meta['Graphics']['Map']['Show Inset Map']=='On':
		N_panel=3
	else:
		N_panel=2

	plt.close('all'); fig,ax=plt.subplots(1,N_panel,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
	im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm,clim=(0,N_tot)) #
	
	#roi['gdf']['bound within'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
	roi['gdf']['bc_bound'].plot(ax=ax[0],edgecolor=[0.6,0.6,0.6],facecolor='none',linewidth=0.25)
	roi['gdf']['lakes'].plot(ax=ax[0],facecolor=[0.5,0.8,1],label='Lakes',linewidth=0.25)
	roi['gdf']['rivers'].plot(ax=ax[0],color=[0.5,0.8,1],label='Rivers',linewidth=0.25)
	roi['gdf']['rivermajor'].plot(ax=ax[0],facecolor=[0.5,0.8,1],edgecolor=[0.5,0.8,1],label='Rivers',linewidth=0.5)

	roi['gdf']['fnc'].plot(ax=ax[0],marker='^',edgecolor=[0.5,0,0],facecolor=[1,0,0],linewidth=0.25,label='TPF',alpha=1,zorder=1,markersize=5)
	for x,y,label in zip(roi['gdf']['fnc'].geometry.x,roi['gdf']['fnc'].geometry.y,roi['gdf']['fnc']['FIRST_NATION_BC_NAME' ]):
		ax[0].annotate(label,xy=(x,y),xytext=(2,1.5),textcoords="offset points",color=[0,0,0],fontsize=4)

	ind=np.where( (roi['gdf']['cities']['Territory']=='BC') & (roi['gdf']['cities']['Level']==1) )[0]
	roi['gdf']['cities'].iloc[ind].plot(ax=ax[0],marker='s',edgecolor=[0.5,0,0],facecolor=[1,0,0],lw=0.5,markersize=9,alpha=1,zorder=2)
	for x,y,label in zip(roi['gdf']['cities'].iloc[ind].geometry.x,roi['gdf']['cities'].iloc[ind].geometry.y,roi['gdf']['cities'].iloc[ind]['City Name']):
		ax[0].annotate(label,xy=(x,y),xytext=(2,1.5),textcoords="offset points",color=[0,0,0],fontsize=4)

	ind=np.where( (roi['gdf']['cities']['Territory']=='BC') & (roi['gdf']['cities']['Level']==2) )[0]
	roi['gdf']['cities'].iloc[ind].plot(ax=ax[0],marker='o',edgecolor=[0.5,0,0],facecolor=[1,0,0],lw=0.25,markersize=3,alpha=1,zorder=2)
	for x,y,label in zip(roi['gdf']['cities'].iloc[ind].geometry.x,roi['gdf']['cities'].iloc[ind].geometry.y,roi['gdf']['cities'].iloc[ind]['City Name']):
		ax[0].annotate(label,xy=(x,y),xytext=(2,1.5),textcoords="offset points",color=[0,0,0],fontsize=4)

	#roi['gdf']['popp'].plot(ax=ax[0],marker='o',edgecolor=[0.75,0.3,0],facecolor=[1,0.6,0],lw=0.25,markersize=6,alpha=1,zorder=2)
	#for x,y,label in zip(roi['gdf']['popp'].geometry.x,roi['gdf']['popp'].geometry.y,roi['gdf']['popp']['NAME']):
	#	ax[0].annotate(label,xy=(x,y),xytext=(2,1.5),textcoords="offset points",color=[0,0,0],fontsize=4)

	ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])
	
	cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_vis+1,1),ticks=np.arange(0.5,N_vis,1))
	cb.ax.set(yticklabels=lab)
	cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
	cb.outline.set_edgecolor('w')
	for i in range(0,N_vis):
		ax[1].plot([0,100],[i,i],'w-',linewidth=1.5)
	pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
	ax[1].set(position=pos2)

	if meta['Graphics']['Map']['Show Inset Map']=='On':
		meta['gdf']['bc_bound']['gdf'].plot(ax=ax[2],edgecolor='k',facecolor=[0.8,0.8,0.8],linewidth=0.25)
		minx, miny, maxx, maxy = meta['gdf']['bc_bound']['gdf'].geometry.total_bounds
		roi['gdf']['bound within'].plot(ax=ax[2],edgecolor=[0,0,0],facecolor=[0,0,0],linewidth=0.25)
		ax[2].set(position=[meta['Graphics']['Map']['Legend X'],0.1,1-meta['Graphics']['Map']['Legend X']-0.01,1-meta['Graphics']['Map']['Legend X']+0.05],xlim=[minx,maxx],ylim=[miny,maxy])
		ax[2].yaxis.set_ticks_position('both'); ax[2].xaxis.set_ticks_position('both'); ax[2].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[2].axis(meta['Graphics']['Map']['Map Axis Vis'])

	gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_settlements','png',900)
	return fig,ax
 
#%%
def Plot_LC20_CEC(meta,roi,vnam):

	dCEC=gu.ReadExcel(meta['Paths']['Model']['Parameters'] + '\\LUT_lc_cec_Compressed.xlsx')

	N_vis=dCEC['Name'].size
	N_hidden=1
	N_tot=N_vis+N_hidden

	z1=roi['grd'][vnam]['Data']
	z1[(roi['grd']['Data']==0)]=np.max(z1)
	z1[0,0:N_vis]=np.arange(1,N_vis+1,1)
	
	# Labels
	lab=list(meta['LUT']['Derived']['lc_cec_c'].keys())
	#lab.append('')
	#lab.append('')

	# Colormap
	cm=plt.cm.get_cmap('viridis',N_vis);
	for i in range(N_vis):
		cm.colors[i,0:3]=[ dCEC['R'][i]/255,dCEC['G'][i]/255,dCEC['B'][i]/255 ]
	cm=np.vstack( (cm.colors,(1,1,1,1)) )
	cm=matplotlib.colors.ListedColormap(cm)

	plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
	im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
	roi=PlotVectorBaseMaps(meta,roi,ax[0])
	ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'],xticklabels='',yticklabels='')
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

	zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
	cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
	cb=plt.colorbar(im,cax=ax[1],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
	ax[1].set(position=[0.71,0.6,0.05,0.14])
	cb.ax.set(yticklabels=lab)
	cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
	cb.outline.set_edgecolor('w')
	for i in range(cb_bnd.size):
		ax[1].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)

	pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
	ax[1].set(position=pos2)
	gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_' + vnam,'png',900)

	return fig,ax

#%%
def Plot_TreeDensityClass(meta,roi,vnam):

	lab=['Treed sparse','Treed open','Treed dense','Shrubland','Grassland','Non-treed']
	z1=roi['grd'][vnam]['Data']
	z1[roi['grd']['Data']==0]=np.max(z1)
	z1[0,0:6]=np.arange(1,7,1)

	N_vis=6
	N_hidden=1
	N_tot=N_vis+N_hidden

	# Colormap
	cm=np.vstack( ((0.75,0.95,0.6,1),(0,0.7,0,1),(0,0.25,0,1),(0.8,0.6,0.6,1),(1,1,0.5,1),(0.85,0.85,0.85,1),(1,1,1,1)) )
	cm=matplotlib.colors.ListedColormap(cm)

	plt.close('all')
	fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
	im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
	roi=PlotVectorBaseMaps(meta,roi,ax[0])
	ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

	zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
	cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
	cb=plt.colorbar(im,cax=ax[1],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
	ax[1].set(position=[0.71,0.6,0.05,0.14])
	try:
		cb.ax.set(yticklabels=lab)
		cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
	except:
		lab.append('')
		cb.ax.set(yticklabels=lab)
		cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
	cb.outline.set_edgecolor('w')
	for i in range(cb_bnd.size):
		ax[1].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)

	pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
	ax[1].set(position=pos2)

	gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_' + vnam,'png',900)

	return fig,ax

#%%
def Plot_TreeDensityGT30(meta,roi,vnam):
	lab=['> 3000','< 3000','Non-treed']
	z1=roi['grd'][vnam]['Data']
	z1[(roi['grd'][vnam]['Data']==0) & (roi['grd']['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Forest'])]=2
	z1[(roi['grd']['lc_comp1_2019']['Data']!=meta['LUT']['Derived']['lc_comp1']['Forest'])]=3
	z1[(roi['grd']['Data']==0)]=4
	z1[0,0:4]=np.arange(1,5,1)

	N_vis=3
	N_hidden=1
	N_tot=N_vis+N_hidden

	# Colormap
	cm=np.vstack( ((0.2,0.2,0.2,1),(0.8,0.8,0.8,1),(0.93,0.93,0.93,1),(1,1,1,1)) )
	cm=matplotlib.colors.ListedColormap(cm)

	plt.close('all')
	fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
	im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
	roi=PlotVectorBaseMaps(meta,roi,ax[0])
	ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

	zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
	cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
	cb=plt.colorbar(im,cax=ax[1],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
	ax[1].set(position=[0.71,0.6,0.05,0.14])
	try:
		cb.ax.set(yticklabels=lab)
		cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
	except:
		lab.append('')
		cb.ax.set(yticklabels=lab)
		cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
	cb.outline.set_edgecolor('w')
	for i in range(cb_bnd.size):
		ax[1].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)

	pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
	ax[1].set(position=pos2)

	gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_' + vnam,'png',900)

	return fig,ax

#%%
def Plot_SI(meta,roi,vnam):
	z=roi['grd'][vnam]['Data']
	if roi['Type']=='ByTSA':
		z[roi['grd']['Data']==0]=0

	plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
	im=ax[0].matshow(roi['grd']['SITE_INDEX']['Data'],extent=roi['grd']['Extent'],cmap='magma',clim=[5,22])

	#roi['gdf']['tsa'].plot(ax=ax[0],color=None,edgecolor=[0,0,0],facecolor='none')
	ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'],xticklabels='',yticklabels='')
	ax[0].grid(False)

	cb=plt.colorbar(im,cax=ax[1])#,boundaries=np.arange(0,N_vis-(N_hidden-1),1),ticks=np.arange(0.5,N_vis+1.5,1))
	pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
	ax[1].set(position=pos2)
	gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + vnam,'png',900)
	
	return fig,ax

#%%
def Plot_PFI(meta,roi):

	# Grid
	bw=20; bin=np.arange(0,140+bw,bw);
	z1=(bin.size)*np.ones( roi['grd']['pfi_c']['Data'].shape)
	for i in range(bin.size):
		ind=np.where(np.abs( roi['grd']['pfi_c']['Data']-bin[i])<=bw/2)
		if ind[0].size>0:
			z1[ind]=i
		else:
			z1[0,i]=i
	z1[(roi['grd']['Data']==1) & ( roi['grd']['pfi_c']['Data']==0)]=i+1
	z1[(roi['grd']['Data']!=1)]=i+2
	L=i+2

	lab=bin.astype(str)

	# Colormap
	cm=plt.cm.get_cmap('cividis',i)
	cm=np.vstack( (cm.colors,(0.9,0.9,0.9,1),(1,1,1,1)) )
	cm=matplotlib.colors.ListedColormap(cm)

	N_vis=bin.size+3
	N_hidden=3

	# Plot
	plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
	im=ax[0].matshow(z1,clim=(0,L+1),extent=roi['grd']['Extent'],cmap=cm)
	roi=PlotVectorBaseMaps(meta,roi,ax[0])
	ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

	cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_vis-(N_hidden-1),1),ticks=np.arange(0.5,N_vis-N_hidden,1))
	cb.ax.set(yticklabels=lab)
	cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
	cb.outline.set_edgecolor('w')
	for i in range(0,N_vis):
		ax[1].plot([0,100],[i,i],'w-',linewidth=1.5)
	pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
	ax[1].set(position=pos2)

	return fig,ax

#%%

def Plot_PFI(meta,roi):

	zB=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\ForestInventory\LiDAR\Boundary TSA\From GQ 20230302\BoundaryTSA_PFI_LnTotalWSV.tif')
	zB['Data']=0.5*0.5*zB['Data']
	zB['Data']=zB['Data'].astype('int16')

	# Grid
	bw=10; bin=np.arange(0,160+bw,bw);
	z1=(bin.size)*np.ones( zB['Data'].shape)
	for i in range(bin.size):
		ind=np.where(np.abs( zB['Data']-bin[i])<=bw/2)
		if ind[0].size>0:
			z1[ind]=i
		else:
			z1[0,i]=i
	ind=np.where(zB['Data']<0); z1[ind]=0
	ind=np.where(zB['Data']>bin[i]); z1[ind]=i
	z1[0,0]=i+1
	z1[1,1]=i+2
	L=i+2

	lab=bin.astype(str)

	# Colormap
	cm=plt.cm.get_cmap('cividis',i)
	cm=np.vstack( (cm.colors,(0.9,0.9,0.9,1),(1,1,1,1)) )
	cm=matplotlib.colors.ListedColormap(cm)

	N_vis=bin.size+3
	N_hidden=3

	# Plot
	plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
	im=ax[0].matshow(z1,clim=(0,L+1),extent=roi['grd']['Extent'],cmap=cm)
	roi=PlotVectorBaseMaps(meta,roi,ax[0])
	ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

	cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_vis-(N_hidden-1),1),ticks=np.arange(0.5,N_vis-N_hidden,1))
	cb.ax.set(yticklabels=lab)
	cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
	cb.outline.set_edgecolor('w')
	for i in range(0,N_vis):
		ax[1].plot([0,100],[i,i],'w-',linewidth=1.5)
	ax[1].set(position=meta['Graphics']['ax2 pos long']);
	#ind=np.where(roi['gdf']['cc']['gdf']['HARVEST_YEAR']==1995)
	roi['gdf']['cc']['gdf'].plot(ax=ax[0],facecolor='None',edgecolor=[0.25,0.75,1],linewidth=1,label='Opening')
	roi['gdf']['fcres']['gdf'].plot(ax=ax[0],facecolor='None',edgecolor=[0.65,1,0.1],linewidth=1,linestyle='--',label='Opening')
	#for idx, row in roi['gdf']['cc']['gdf'].iterrows():
		#print(row.keys())
		#break
	#	plt.text(row['geometry'].centroid.x,row['geometry'].centroid.y,str(row['HARVEST_YEAR']),ha='center',fontsize=5,color=[0.25,0.75,1])
	return fig,ax

#%%

def Plot_SoilOrganicCarbon_GSOC(meta,roi,vnam):

	# Grid
	bw=10; bin=np.arange(0,100+bw,bw);
	z1=(bin.size)*np.ones( roi['grd'][vnam]['Data'].shape)
	for i in range(bin.size):
		ind=np.where(np.abs( roi['grd'][vnam]['Data']-bin[i])<=bw/2)
		z1[ind]=i
	ind=np.where(roi['grd'][vnam]['Data']>=bin[i]); z1[ind]=i
	z1[(roi['grd']['Data']==1) & ( roi['grd'][vnam]['Data']==0)]=i+1
	z1[(roi['grd']['Data']!=1)]=i+2
	L=i+2

	lab=bin.astype(str)

	# Colormap
	cm=plt.cm.get_cmap('viridis',i)
	#cm=plt.cm.get_cmap('plasma',i)
	cm=np.vstack( (cm.colors,(0.9,0.9,0.9,1),(1,1,1,1)) )
	cm=matplotlib.colors.ListedColormap(cm)

	N_vis=bin.size+3
	N_hidden=3

	# Plot
	plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
	im=ax[0].matshow(z1,clim=(0,L+1),extent=roi['grd']['Extent'],cmap=cm)
	roi=PlotVectorBaseMaps(meta,roi,ax[0])
	ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

	cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_vis-(N_hidden-1),1),ticks=np.arange(0.5,N_vis-N_hidden,1))
	cb.ax.set(yticklabels=lab)
	cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
	cb.outline.set_edgecolor('w')
	for i in range(0,N_vis):
		ax[1].plot([0,100],[i,i],'w-',linewidth=1.5)

	pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
	ax[1].set(position=pos2)
	gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_' + vnam,'png',900)

	return fig,ax

#%%
def Plot_CrownCover(meta,roi,vnam):

	bw=10; bin=np.arange(0,100+bw,bw);
	z1=(bin.size)*np.ones( roi['grd'][vnam]['Data'].shape)
	for i in range(bin.size):
		ind=np.where(np.abs( roi['grd'][vnam]['Data']-bin[i])<=bw/2)
		z1[ind]=i
	ind=np.where(roi['grd'][vnam]['Data']>bin[i])
	z1[ind]=0
	z1[(roi['grd']['Data']==1) & ( roi['grd'][vnam]['Data']==0)]=i+1
	z1[(roi['grd']['Data']!=1)]=i+2
	L=i+2

	lab=bin.astype(str)

	# Colormap
	cm=plt.cm.get_cmap('viridis',i)
	#cm=plt.cm.get_cmap('plasma',i)
	cm=np.vstack( (cm.colors,(0.9,0.9,0.9,1),(1,1,1,1)) )
	cm=matplotlib.colors.ListedColormap(cm)

	N_vis=bin.size+3
	N_hidden=3

	# Plot
	plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
	im=ax[0].matshow(z1,clim=(0,L+1),extent=roi['grd']['Extent'],cmap=cm)
	ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'],xticklabels='',yticklabels='')
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

	cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_vis-(N_hidden-1),1),ticks=np.arange(0.5,N_vis-N_hidden,1))
	cb.ax.set(yticklabels=lab)
	cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
	cb.outline.set_edgecolor('w')
	for i in range(0,N_vis):
		ax[1].plot([0,100],[i,i],'w-',linewidth=1.5)
	pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
	ax[1].set(position=pos2)

	gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_' + vnam,'png',900)

	return fig,ax

#%%
def Plot_SPH(meta,roi,vnam):

	# Grid
	bw=400; bin=np.arange(0,2400,bw);
	z1=(bin.size)*np.ones( roi['grd'][vnam]['Data'].shape)
	for i in range(bin.size):
		ind=np.where(np.abs( roi['grd'][vnam]['Data']-bin[i])<=bw/2)
		z1[ind]=i
	ind=np.where(roi['grd'][vnam]['Data']>bin[i])
	z1[ind]=i
	z1[(roi['grd']['Data']==1) & ( roi['grd'][vnam]['Data']==0)]=i+1
	z1[(roi['grd']['Data']!=1)]=i+2
	L=i+2

	lab=bin.astype(str)

	# Colormap
	#cm=plt.cm.get_cmap('viridis',i)
	cm=plt.cm.get_cmap('plasma',i)
	cm=np.vstack( (cm.colors,(0.9,0.9,0.9,1),(1,1,1,1)) )
	cm=matplotlib.colors.ListedColormap(cm)

	N_vis=bin.size+3
	N_hidden=3

	# Plot
	plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
	im=ax[0].matshow(z1,clim=(0,L+1),extent=roi['grd']['Extent'],cmap=cm)
	#roi['gdf']['tsa'].plot(ax=ax[0],color=None,edgecolor=[0,0,0],facecolor='none')
	ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
	ax[0].grid(False)

	cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_vis-(N_hidden-1),1),ticks=np.arange(0.5,N_vis-N_hidden,1))
	cb.ax.set(yticklabels=lab)
	cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
	for i in range(0,N_vis):
		ax[1].plot([0,100],[i/(N_vis-N_hidden),i/(N_vis-N_hidden)],'k-',linewidth=0.5)
	ax[1].set(position=meta['Graphics']['ax2 pos long']);

	# cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_vis-(N_hidden-1),1),ticks=np.arange(0.5,N_vis+1.5,1))
	# cb.ax.set(yticklabels=lab)
	# cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
	# for i in range(0,N_vis):
	#	 ax[1].plot([0,100],[i/(N_vis-N_hidden),i/(N_vis-N_hidden)],'k-',linewidth=0.5)
	# ax[1].set(position=meta['Graphics']['ax2 pos long']);
	# pos2=copy.copy(meta['Graphics']['pos2'])
	# pos2[1]=0.6
	# pos2[3]=0.24
	# ax[1].set(position=pos2)
	gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_' + vnam,'png',900)

	return fig,ax

#%%

def Plot_LUC_Year(meta,roi,vnam):

	z0=roi['grd'][vnam]['Data']

	bw=10; bin=np.arange(1850,2020+bw,bw);
	z1=(bin.size)*np.ones( z0.shape)
	for i in range(bin.size):
		ind=np.where(np.abs( z0-bin[i])<=bw/2)
		z1[ind]=i
	z1[(roi['grd']['Data']==0) | (roi['grd']['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Water'])]=i+2
	#z1[(roi['grd']['Data']!=1)]=i+2
	L=i+2

	lab=bin.astype(str)

	# Colormap
	#cm=plt.cm.get_cmap('viridis',i)
	cm=plt.cm.get_cmap('plasma',i)
	cm=np.vstack( (cm.colors,(0.85,0.85,0.85,1),(1,1,1,1)) )
	cm=matplotlib.colors.ListedColormap(cm)

	N_vis=bin.size+3
	N_hidden=3

	# Plot
	plt.close('all')
	fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
	im=ax[0].matshow(z1,clim=(0,L+1),extent=roi['grd']['Extent'],cmap=cm)
	roi=PlotVectorBaseMaps(meta,roi,ax[0])
	ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

	cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_vis-(N_hidden-1),1),ticks=np.arange(0.5,N_vis-N_hidden,1))
	cb.ax.set(yticklabels=lab)
	cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
	cb.outline.set_edgecolor('w')
	for i in range(0,N_vis):
		ax[1].plot([0,100],[i,i],'w-',linewidth=1.5)
	pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
	ax[1].set(position=pos2)
	gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_' + vnam,'png',900)
	return fig,ax

#%%
def Plot_WildfireYear(meta,roi,vnam):

	bw=10; bin=np.arange(1910,2020+bw,bw);
	z1=(bin.size)*np.ones( roi['grd'][vnam]['Data'].shape)
	for i in range(bin.size):
		ind=np.where(np.abs( roi['grd'][vnam]['Data']-bin[i])<=bw/2)
		z1[ind]=i
	z1[(roi['grd']['Data']==0) | (roi['grd']['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Water'])]=i+2
	#z1[(roi['grd']['Data']!=1)]=i+2
	L=i+2
	lab=bin.astype(str)

	# Colormap
	#cm=plt.cm.get_cmap('viridis',i)
	cm=plt.cm.get_cmap('plasma',i)
	cm=np.vstack( (cm.colors,(0.85,0.85,0.85,1),(1,1,1,1)) )
	cm=matplotlib.colors.ListedColormap(cm)

	N_vis=bin.size+3
	N_hidden=3

	# Plot
	plt.close('all')
	fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
	im=ax[0].matshow(z1,clim=(0,L+1),extent=roi['grd']['Extent'],cmap=cm)
	roi=PlotVectorBaseMaps(meta,roi,ax[0])
	ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

	cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_vis-(N_hidden-1),1),ticks=np.arange(0.5,N_vis-N_hidden,1))
	cb.ax.set(yticklabels=lab)

	cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
	cb.outline.set_edgecolor('w')
	for i in range(0,N_vis):
		ax[1].plot([0,100],[i,i],'w-',linewidth=1.5)

	pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
	ax[1].set(position=pos2)
	gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_' + vnam,'png',900)

	return fig,ax

#%%
def Plot_GFC_LossYear(meta,roi,vnam):

	z0=roi['grd'][vnam]['Data']
	#z0=roi['grd']['gfcly_filt']['Data']

	bw=1; bin=np.arange(2001,2022,bw);
	z1=(bin.size)*np.ones( z0.shape)
	for i in range(bin.size):
		ind=np.where(z0==bin[i])
		z1[ind]=i
	z1[(z0==0)]=i+1
	#z1[(roi['grd']['Data']==0)]=i+2
	z1[(roi['grd']['Data']==0) | (roi['grd']['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Water'])]=i+2
	L=i+2

	lab=np.array(bin).astype(str)

	# Colormap
	#cm=plt.cm.get_cmap('viridis',i)
	cm=plt.cm.get_cmap('plasma',i)
	cm=np.vstack( (cm.colors,(0.85,0.85,0.85,1),(1,1,1,1)) )
	cm=matplotlib.colors.ListedColormap(cm)

	N_vis=bin.size+3
	N_hidden=3

	plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
	im=ax[0].matshow(z1,clim=(0,L+1),extent=roi['grd']['Extent'],cmap=cm)
	roi=PlotVectorBaseMaps(meta,roi,ax[0])
	ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'],xticklabels='',yticklabels='')
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

	cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_vis-(N_hidden-1),1),ticks=np.arange(0.5,N_vis-N_hidden,1))
	cb.ax.set(yticklabels=lab)

	cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
	cb.outline.set_edgecolor('w')
	for i in range(0,N_vis):
		ax[1].plot([0,100],[i,i],'w-',linewidth=1.5)

	pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
	ax[1].set(position=pos2)

	gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_' + vnam,'png',900)
	#gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_gfcly_filt','png',900)

	return fig,ax

#%%
def Plot_SoilWaterContent(meta,roi,vnam):
	bw=20; bin=np.arange(0,200+bw,bw);
	z1=bin.size*np.ones( roi['grd'][vnam]['Data'].shape)
	for i in range(bin.size):
		ind=np.where(np.abs( roi['grd'][vnam]['Data']-bin[i])<=bw/2)
		z1[ind]=i
	z1[(roi['grd'][vnam]['Data']>200)]=i
	z1[(roi['grd']['lc_comp1_2019']['Data']==0) | (roi['grd']['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Water']) | (roi['grd']['Data']!=1)]=i+1
	#z1[(roi['grd']['lc_comp1_2019']==0) | (roi['grd']['Data']!=1)]=i+1
	for i in range(bin.size):
		z1[0,i]=i

	lab=bin.astype(str)
	lab=np.append(lab,'Hidden')

	N_vis=bin.size
	N_hidden=1
	N_tot=N_vis+N_hidden
	cm0=gu.ReadExcel(r'C:\Data\Colormaps\colormap_ws.xlsx')
	cm=np.column_stack((cm0['cl1'],cm0['cl2'],cm0['cl3'],np.ones(cm0['bin'].size)))
	cm=np.vstack( (cm,(1,1,1,1)) ) # (0.83137,0.81569,0.78431,1)
	cm=matplotlib.colors.ListedColormap(cm)

	plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
	im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
	roi=PlotVectorBaseMaps(meta,roi,ax[0])
	ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

	zmn=np.min(z1);
	zmx=np.max(z1);
	cb_ivl=(zmx-zmn)/N_tot;
	cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
	cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
	cb=plt.colorbar(im,cax=ax[1],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
	ax[1].set(position=[0.71,0.6,0.05,0.14])
	try:
		cb.ax.set(yticklabels=lab)
	except:
		lab=np.append(lab,'Hidden')
		cb.ax.set(yticklabels=lab)
	cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
	cb.outline.set_edgecolor('w')
	for i in range(cb_bnd.size):
		ax[1].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)

	pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
	ax[1].set(position=pos2)

	gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_' + vnam,'png',300)
	return fig,ax



#%% Empty forest mask
# *** For use with custom vector layers ***
def Plot_ForestMask(meta,roi):

	z1=3*np.ones(roi['grd']['Data'].shape,dtype='int8')
	ind=np.where( (roi['grd']['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (roi['grd']['Data']>0) )
	z1[ind]=1
	ind=np.where( (roi['grd']['lc_comp1_2019']['Data']!=meta['LUT']['Derived']['lc_comp1']['Forest']) & (roi['grd']['Data']>0) )
	z1[ind]=2

	if (meta['Graphics']['Map']['Show Inset Map']=='On'):
		N_panel=3
	else:
		N_panel=2

	# Number of colours and number of colours excluded from colorbar
	N_vis=2
	N_hidden=1
	N_tot=N_vis+N_hidden
	lab=np.array(['Forest','Non-forest'])

	cm=np.vstack( ((0.65,0.75,0.5,1),(0.9,0.9,0.9,1),(1,1,1,1)) )
	cm=matplotlib.colors.ListedColormap(cm)

	plt.close('all'); fig,ax=plt.subplots(1,N_panel,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
	im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
	roi=PlotVectorBaseMaps(meta,roi,ax[0])
	ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

	zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
	cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
	cb=plt.colorbar(im,cax=ax[1],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
	ax[1].set(position=[0.65,0.75,0.03,0.07])
	cb.ax.set(yticklabels=lab)
	cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
	cb.outline.set_edgecolor('w')
	for i in range(cb_bnd.size):
		ax[1].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)
	#ax[1].axis(meta['Graphics']['Map']['Map Axis Vis']) # Exclude colorbar

	# Scalebar
	if meta['Graphics']['Map']['Show Scalebar']=='On':
		AddScalebar(ax[0],roi)

	if meta['Graphics']['Map']['Show Inset Map']=='On':
		meta['gdf']['bc_bound']['gdf'].plot(ax=ax[2],edgecolor=None,facecolor=[0.8,0.8,0.8],linewidth=0.25)
		minx, miny, maxx, maxy = meta['gdf']['bc_bound']['gdf'].geometry.total_bounds
		roi['gdf']['bound within'].plot(ax=ax[2],edgecolor=None,facecolor=[0,0,0],linewidth=0.25)
		ax[2].set(position=[meta['Graphics']['Map']['Legend X'],0.1,1-meta['Graphics']['Map']['Legend X']-0.01,1-meta['Graphics']['Map']['Legend X']+0.05],xlim=[minx,maxx],ylim=[miny,maxy])
		ax[2].yaxis.set_ticks_position('both'); ax[2].xaxis.set_ticks_position('both'); ax[2].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[2].axis(meta['Graphics']['Map']['Map Axis Vis'])

	# *** Add Vector layers ***
	
	# points=[]
	# for k in range(x.size):
	#	 points.append(Point(x[k],y[k]))
	# gdf1=gpd.GeoDataFrame({'geometry':points,'ID':np.ones(x.size)})
	# gdf1.crs=roi['crs']
	# gdf1=gpd.overlay(gdf1,roi['gdf']['bound within'],how='intersection')
	# gdf1.plot(ax=ax[0],markersize=0.7,facecolor=[0.75,0,0],edgecolor=None,linewidth=0.75,alpha=1,label='Proposed sample')

	gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_forest_mask','png',900)
	return fig,ax

#%%
def Plot_FieldPlots(meta,roi):

	z1=roi['grd']['Data'].copy()
	#ind=np.where(roi['grd']['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Forest'])
	#z1[ind]=2

	ind=np.where( (roi['grd']['bgcz']['Data']==meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE']['CWH']) )
	z1[ind]=2

	ind=np.where( (roi['grd']['bgcz']['Data']==meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE']['ESSF']) | (roi['grd']['bgcz']['Data']==meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE']['SWB']) )
	z1[ind]=3

	ind=np.where( (roi['grd']['bgcz']['Data']==meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE']['SBS']) | \
			(roi['grd']['bgcz']['Data']==meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE']['MS']) | \
			(roi['grd']['bgcz']['Data']==meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE']['SBPS']) | \
			(roi['grd']['bgcz']['Data']==meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE']['IDF']) )
	z1[ind]=4

	#ind=np.where( (roi['grd']['bgcz']['Data']==meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE']['BWBS']) | (roi['grd']['bgcz']['Data']==meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE']['SWB']) )
	#z1[ind]=5

	cm=np.vstack( ((1,1,1,1),
				(0.95,0.95,0.95,1),
				(0.8,0.95,0.65,1),
				(0.8,0.9,1,1),
				(1,1,0.8,1)) )
	cm=matplotlib.colors.ListedColormap(cm)

	meta['Graphics']['Map']['Fig Width']=14

	plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
	ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
	#roi=PlotVectorBaseMaps(meta,roi,ax[0])
	roi['gdf']['bc_bound'].plot(ax=ax[0],edgecolor=[0.5,0.5,0.5],facecolor='none',linewidth=0.25)
	ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

	flg=0
	if flg==1:
		# Grid
		ivl=32
		ind=np.where(z1[0::ivl,0::ivl]==2)
		x=roi['grd']['X'][0::ivl,0::ivl][ind]
		y=roi['grd']['Y'][0::ivl,0::ivl][ind]
	
		points=[]
		for k in range(x.size):
			points.append(Point(x[k],y[k]))
		gdf1=gpd.GeoDataFrame({'geometry':points,'ID':np.ones(x.size)})
		gdf1.crs=roi['crs']
		gdf1=gpd.overlay(gdf1,roi['gdf']['bound within'],how='intersection')
		gdf1.plot(ax=ax[0],markersize=0.7,facecolor=[0.75,0,0],edgecolor=None,linewidth=0.75,alpha=1,label='Proposed sample')

	# CN plots
	flg=0
	if flg==1:
		ind=np.where( (fp['PTF CNY']==1) & (np.isnan(fp['Year t1'])==True) )[0]
		x=fp['X BC'][ind]
		y=fp['Y BC'][ind]
		points=[]
		for k in range(x.size):
			points.append(Point(x[k],y[k]))
		gdf1=gpd.GeoDataFrame({'geometry':points,'ID':np.ones(x.size)})
		gdf1.crs=roi['crs']
		gdf1=gpd.overlay(gdf1,roi['gdf']['bound within'],how='intersection')
		gdf1.plot(ax=ax[0],markersize=0.7,facecolor=[0,0,0],edgecolor=None,linewidth=0.75,alpha=1,label='Single measurement')

	ind=np.where( (fp['PTF CN']==1) )[0]
	#ind=np.where( (fp['PTF CN']==1) & (np.isnan(fp['Year t1'])==False) )[0]
	#ind=np.where( (np.isnan(fp['Year t1'])==False) )[0]
	x=fp['X BC'][ind]
	y=fp['Y BC'][ind]
	points=[]
	for k in range(x.size):
		points.append(Point(x[k],y[k]))
	gdf1=gpd.GeoDataFrame({'geometry':points,'ID':np.ones(x.size)})
	gdf1.crs=roi['crs']
	gdf1=gpd.overlay(gdf1,roi['gdf']['bound within'],how='intersection')
	gdf1.plot(ax=ax[0],marker='o',markersize=3,facecolor=[1,1,1],edgecolor=[0,0,0],linewidth=0.5,alpha=1,label='With remeasurement')

	# Number of colours and number of colours excluded from colorbar
	N_vis=4
	N_hidden=1
	N_tot=N_vis+N_hidden

	# Colormap
	cm=np.vstack( ((0.85,1,0.5,1),(0,0.5,0,1),(0.6,0.6,0.6,1),(0.83,0.83,0.83,1),(1,1,1,1)) )
	cm=matplotlib.colors.ListedColormap(cm)

	plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
	im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
	roi=PlotVectorBaseMaps(meta,roi,ax[0])
	ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'],xticklabels='',yticklabels='')
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

	zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
	cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
	cb=plt.colorbar(im,cax=ax[1],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
	ax[1].set(position=[0.71,0.6,0.05,0.14])
	cb.ax.set(yticklabels=lab)
	cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
	cb.outline.set_edgecolor('w')
	for i in range(cb_bnd.size):
		ax[1].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)

	pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
	ax[1].set(position=pos2)

	#ax[0].legend(loc='lower left',facecolor=[1,1,1],frameon=False)
	ax[1].axis(meta['Graphics']['Map']['Map Axis Vis'])
	#gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_FieldPlots','png',900)
	gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Field Plots\Simple Biomass Model vs Observations\Fig1_Map','png',900)
	return
#%%
def Plot_FieldPlotsSBM(meta,roi,fp):

	z1=5*np.ones(roi['grd']['Data'].shape,dtype='int8')
	ind=np.where(roi['grd']['Data']==1)
	z1[ind]=4

	ind=np.where( (roi['grd']['bgcz']['Data']==meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE']['SBS']) | \
			(roi['grd']['bgcz']['Data']==meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE']['MS']) | \
			(roi['grd']['bgcz']['Data']==meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE']['SBPS']) | \
			(roi['grd']['bgcz']['Data']==meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE']['IDF']) | \
			(roi['grd']['bgcz']['Data']==meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE']['BWBS']) )
	z1[ind]=1

	ind=np.where( (roi['grd']['bgcz']['Data']==meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE']['ESSF']) | (roi['grd']['bgcz']['Data']==meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE']['SWB']) )
	z1[ind]=2

	ind=np.where( (roi['grd']['bgcz']['Data']==meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE']['CWH']) )
	z1[ind]=3

	cm=np.vstack( ((1,1,0.8,1),(0.8,0.9,1,1),(0.8,0.95,0.65,1),(0.95,0.95,0.95,1),(1,1,1,1)) )
	cm=matplotlib.colors.ListedColormap(cm)

	meta['Graphics']['Map']['Fig Width']=14

	plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
	im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
	roi['gdf']['bc_bound'].plot(ax=ax[0],edgecolor=[0.5,0.5,0.5],facecolor='none',linewidth=0.25)
	ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

	ind=np.where( (fp['PTF CN']==1) )[0]
	#ind=np.where( (fp['PTF CN']==1) & (np.isnan(fp['Year t1'])==False) )[0]
	#ind=np.where( (np.isnan(fp['Year t1'])==False) )[0]
	x=fp['X BC'][ind]
	y=fp['Y BC'][ind]
	points=[]
	for k in range(x.size):
		points.append(Point(x[k],y[k]))
	gdf1=gpd.GeoDataFrame({'geometry':points,'ID':np.ones(x.size)})
	gdf1.crs=roi['crs']
	gdf1=gpd.overlay(gdf1,roi['gdf']['bound within'],how='intersection')
	gdf1.plot(ax=ax[0],marker='o',markersize=3,facecolor=[1,1,1],edgecolor=[0,0,0],linewidth=0.5,alpha=1,label='With remeasurement')

	# Number of colours and number of colours excluded from colorbar
	N_vis=3
	N_hidden=2
	N_tot=N_vis+N_hidden
	lab=['Subhumid continental','Humid continental','Humid maritime','']

	zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
	cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
	cb=plt.colorbar(im,cax=ax[1],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
	ax[1].set(position=[0.71,0.6,0.05,0.14])
	cb.ax.set(yticklabels=lab)
	cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
	cb.outline.set_edgecolor('w')
	for i in range(cb_bnd.size):
		ax[1].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)

	pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
	ax[1].set(position=pos2)

	#ax[1].axis(meta['Graphics']['Map']['Map Axis Vis'])

	#gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Field Plots\Simple Biomass Model vs Observations\Fig1_Map','png',900)
	return

#%%
def Plot_FieldPlotsGROMO(meta,roi,fp):

	z1=3*np.ones(roi['grd']['Data'].shape)
	ind=np.where(roi['grd']['Data']==1); z1[ind]=2
	ind=np.where(roi['grd']['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Forest']); z1[ind]=1
	cm=np.vstack( ((0.6,0.8,0.5,1),(0.95,0.95,0.95,1),(1,1,1,1)) )
	cm=matplotlib.colors.ListedColormap(cm)

	meta['Graphics']['Map']['Fig Width']=15
	plt.close('all'); fig,ax=plt.subplots(1,5,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
	im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
	roi['gdf']['bc_bound'].plot(ax=ax[0],edgecolor=[0.5,0.5,0.5],facecolor='none',linewidth=0.25)
	ax[0].set(position=[0,0.5,0.5,0.5],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

	ikp=np.where((fp['Age Mean t0']>0) & (fp['Age Mean t0']<1000) & \
					(fp['Csw G Tot']>0) & (fp['Csw G Tot']<25) & \
					(np.isnan(fp['tmean_wyr_n'])==False) & \
					(fp['Occ_Harv']==0) & \
					(fp['Occ_Wildfire']>=0) & \
					(fp['Occ_IBM']>=0) )[0]
	#ind=np.where( (fp['PTF CN']==1) )[0]
	#ind=np.where( (fp['PTF CN']==1) & (np.isnan(fp['Year t1'])==False) )[0]
	#ind=np.where( (np.isnan(fp['Year t1'])==False) )[0]
	x=fp['X BC'][ikp]
	y=fp['Y BC'][ikp]
	points=[]
	for k in range(x.size):
		points.append(Point(x[k],y[k]))
	gdf1=gpd.GeoDataFrame({'geometry':points,'ID':np.ones(x.size)})
	gdf1.crs=roi['crs']
	gdf1=gpd.overlay(gdf1,roi['gdf']['bound within'],how='intersection')
	gdf1.plot(ax=ax[0],marker='.',markersize=3,facecolor=[1,1,1],edgecolor=[0,0,0],linewidth=0.5,alpha=1,label='With remeasurement')

	# Number of colours and number of colours excluded from colorbar
	N_vis=2
	N_hidden=1
	N_tot=N_vis+N_hidden
	lab=['Forest','Non-forest']

	zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
	cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
	cb=plt.colorbar(im,cax=ax[1],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
	ax[1].set(position=[0.71,0.6,0.05,0.14])
	cb.ax.set(yticklabels=lab)
	cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
	cb.outline.set_edgecolor('w')
	for i in range(cb_bnd.size):
		ax[1].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)

	pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
	ax[1].set(position=[0.01,0.57,0.03,0.05])

	# Add climate space
	aw=0.4; ms=2
	ax[2].plot(fp['Year t0'][ikp],fp['Age Mean t0'][ikp],'k.',ms=ms)
	ax[2].set(position=[0.55,0.55,aw,aw],xlabel='Time, years',ylabel='Stand age, years',ylim=[0,300])
	ax[2].yaxis.set_ticks_position('both'); ax[2].xaxis.set_ticks_position('both'); ax[2].tick_params(length=meta['Graphics']['gp']['tickl'])

	ax[3].plot(fp['tmean_wyr_n'][ikp],fp['cwd_mjjas_n'][ikp]/30.6,'k.',ms=ms)
	ax[3].set(position=[0.075,0.08,aw,aw],xlabel='Air temperature ($\circ$C)',ylabel='Climatic water deficit (mm d$^{-1}$)',xlim=[-3,11],ylim=[0,4.5])
	ax[3].yaxis.set_ticks_position('both'); ax[3].xaxis.set_ticks_position('both'); ax[3].tick_params(length=meta['Graphics']['gp']['tickl'])

	ax[4].plot(fp['CO2'][ikp],fp['ndep'][ikp],'k.',ms=ms)
	ax[4].set(position=[0.55,0.08,aw,aw],xlabel='Carbon dioxide (ppm)',ylabel='Nitrogen deposition (kg N ha$^{-1}$ yr$^{-1}$)',xlim=[300,415],ylim=[0,5])
	ax[4].yaxis.set_ticks_position('both'); ax[4].xaxis.set_ticks_position('both'); ax[4].tick_params(length=meta['Graphics']['gp']['tickl'])

	gu.axletters(ax,plt,0.03,0.92,FontColor=meta['Graphics']['gp']['cla'],LetterStyle='Default',FontWeight=meta['Graphics']['Modelling']['AxesFontWeight'],Skip=1)
	#ax[1].axis(meta['Graphics']['Map']['Map Axis Vis'])
	gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Field Plots\gromo\Fig1_Map','png',900)

	return

#%%
def Plot_FieldPlotsPanels(meta,roi):
	meta2,fp,soc=ufp.ImportFieldPlotData(meta,type='Stand')

	plt.close('all'); fig,ax=plt.subplots(1,4,figsize=gu.cm2inch(17.5,15))

	# VRI plots
	z1=roi['grd']['Data'].copy()
	ind=np.where(roi['grd']['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Forest'])
	z1[ind]=2
	cm=np.vstack( ((1,1,1,1),(0.95,0.95,0.95,1),(0.8,0.8,0.8,1)) )
	cm=matplotlib.colors.ListedColormap(cm)
	ms=2
	im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
	#roi=PlotVectorBaseMaps(meta,roi,ax[0])
	ax[0].set(position=[0,0.5,0.5,0.5],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])
	
	ind=np.where( (fp['PTF CNV']==1) & (fp['PTF CN']==0) )[0]
	x=fp['X BC'][ind]; y=fp['Y BC'][ind]
	points=[]
	for k in range(x.size):
		points.append(Point(x[k],y[k]))
	gdf1=gpd.GeoDataFrame({'geometry':points,'ID':np.ones(x.size)})
	gdf1.crs=roi['crs']
	gdf1=gpd.overlay(gdf1,roi['gdf']['bound within'],how='intersection')
	gdf1.plot(ax=ax[0],markersize=0.7,facecolor=[0.75,0,0],edgecolor=None,linewidth=0.75,alpha=1,label='Proposed sample')
	
	# CMI plots
	z1=roi['grd']['Data'].copy()
	ind=np.where(roi['grd']['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Forest'])
	z1[ind]=2
	cm=np.vstack( ((1,1,1,1),(0.95,0.95,0.95,1),(0.8,0.8,0.8,1)) )
	cm=matplotlib.colors.ListedColormap(cm)
	ms=2
	im=ax[1].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
	#roi=PlotVectorBaseMaps(meta,roi,ax[1])
	ax[1].set(position=[0.5,0.5,0.5,0.5],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
	ax[1].yaxis.set_ticks_position('both'); ax[1].xaxis.set_ticks_position('both'); ax[1].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[1].axis(meta['Graphics']['Map']['Map Axis Vis'])
	
	ind=np.where( (fp['PTF CN']==1) & (np.isnan(fp['Year t1'])==True) )[0]
	x=fp['X BC'][ind]; y=fp['Y BC'][ind]
	points=[]
	for k in range(x.size):
		points.append(Point(x[k],y[k]))
	gdf1=gpd.GeoDataFrame({'geometry':points,'ID':np.ones(x.size)})
	gdf1.crs=roi['crs']
	gdf1=gpd.overlay(gdf1,roi['gdf']['bound within'],how='intersection')
	gdf1.plot(ax=ax[1],markersize=0.7,facecolor=[0.75,0,0],edgecolor=None,linewidth=0.75,alpha=1,label='Single measurement')
	
	ind=np.where( (fp['PTF CN']==1) & (np.isnan(fp['Year t1'])==False) )[0]
	x=fp['X BC'][ind]; y=fp['Y BC'][ind]
	points=[]
	for k in range(x.size):
		points.append(Point(x[k],y[k]))
	gdf1=gpd.GeoDataFrame({'geometry':points,'ID':np.ones(x.size)})
	gdf1.crs=roi['crs']
	gdf1=gpd.overlay(gdf1,roi['gdf']['bound within'],how='intersection')
	gdf1.plot(ax=ax[1],marker='s',markersize=0.7,facecolor=[0.25,0.5,0.75],edgecolor=None,linewidth=0.75,alpha=1,label='With remeasurement')
	ax[1].legend(loc='lower left',facecolor=[1,1,1],frameon=False)
		  
	# YSM plots
	z1=roi['grd']['Data'].copy()
	ind=np.where(roi['grd']['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Forest'])
	z1[ind]=2
	cm=np.vstack( ((1,1,1,1),(0.95,0.95,0.95,1),(0.8,0.8,0.8,1)) )
	cm=matplotlib.colors.ListedColormap(cm)
	ms=2
	im=ax[2].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
	#roi=PlotVectorBaseMaps(meta,roi,ax[2])
	ax[2].set(position=[0,0,0.5,0.5],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
	ax[2].yaxis.set_ticks_position('both'); ax[2].xaxis.set_ticks_position('both'); ax[2].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[2].axis(meta['Graphics']['Map']['Map Axis Vis'])
	
	ind=np.where( (fp['PTF YSM']==1) & (np.isnan(fp['Year t1'])==True) )[0]
	x=fp['X BC'][ind]; y=fp['Y BC'][ind]
	points=[]
	for k in range(x.size):
		points.append(Point(x[k],y[k]))
	gdf1=gpd.GeoDataFrame({'geometry':points,'ID':np.ones(x.size)})
	gdf1.crs=roi['crs']
	gdf1=gpd.overlay(gdf1,roi['gdf']['bound within'],how='intersection')
	gdf1.plot(ax=ax[2],markersize=0.7,facecolor=[0.75,0,0],edgecolor=None,linewidth=0.75,alpha=1,label='Proposed sample')
	
	ind=np.where( (fp['PTF YSM']==1) & (np.isnan(fp['Year t1'])==False) )[0]
	x=fp['X BC'][ind]; y=fp['Y BC'][ind]
	points=[]
	for k in range(x.size):
		points.append(Point(x[k],y[k]))
	gdf1=gpd.GeoDataFrame({'geometry':points,'ID':np.ones(x.size)})
	gdf1.crs=roi['crs']
	gdf1=gpd.overlay(gdf1,roi['gdf']['bound within'],how='intersection')
	gdf1.plot(ax=ax[2],marker='s',markersize=0.7,facecolor=[0.25,0.5,0.75],edgecolor=None,linewidth=0.75,alpha=1,label='Proposed sample')
	
	flg=0
	if flg==1:
		# SOC
		z1=roi['grd']['Data'].copy()
		ind=np.where(roi['grd']['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Forest'])
		z1[ind]=2
		cm=np.vstack( ((1,1,1,1),(0.95,0.95,0.95,1),(0.8,0.8,0.8,1)) )
		cm=matplotlib.colors.ListedColormap(cm)
		ms=2
		im=ax[3].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
		#roi=PlotVectorBaseMaps(meta,roi,ax[3])
		ax[3].set(position=[0.5,0,0.5,0.5],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
		ax[3].yaxis.set_ticks_position('both'); ax[3].xaxis.set_ticks_position('both'); ax[3].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[3].axis(meta['Graphics']['Map']['Map Axis Vis'])
		
		soils=gu.ipickle(r'C:\Data\Soils\Shaw et al 2018 Database\SITES.pkl')
		x=soils['x']; y=soils['y']
		points=[]
		for k in range(x.size):
			points.append(Point(x[k],y[k]))
		gdf1=gpd.GeoDataFrame({'geometry':points,'ID':np.ones(x.size)})
		gdf1.crs=roi['crs']
		gdf1=gpd.overlay(gdf1,roi['gdf']['bound within'],how='intersection')
		gdf1.plot(ax=ax[3],markersize=0.7,facecolor=[0.75,0,0],edgecolor=None,linewidth=0.75,alpha=1,label='Proposed sample')
	else:
		# PSPs
		z1=roi['grd']['Data'].copy()
		ind=np.where(roi['grd']['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Forest'])
		z1[ind]=2
		cm=np.vstack( ((1,1,1,1),(0.95,0.95,0.95,1),(0.8,0.8,0.8,1)) )
		cm=matplotlib.colors.ListedColormap(cm)
		ms=2
		im=ax[3].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
		#roi=PlotVectorBaseMaps(meta,roi,ax[2])
		ax[3].set(position=[0.5,0,0.5,0.5],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
		ax[3].yaxis.set_ticks_position('both'); ax[3].xaxis.set_ticks_position('both');
		ax[3].grid(meta['Graphics']['Map']['Map Grid Vis'])
		ax[3].axis(meta['Graphics']['Map']['Map Axis Vis'])
		
		ind=np.where( (fp['PTF PSP']==1) & (np.isnan(fp['Year t1'])==True) )[0]
		x=fp['X BC'][ind]; y=fp['Y BC'][ind]
		points=[]
		for k in range(x.size):
			points.append(Point(x[k],y[k]))
		gdf1=gpd.GeoDataFrame({'geometry':points,'ID':np.ones(x.size)})
		gdf1.crs=roi['crs']
		gdf1=gpd.overlay(gdf1,roi['gdf']['bound within'],how='intersection')
		gdf1.plot(ax=ax[3],markersize=0.7,facecolor=[0.75,0,0],edgecolor=None,linewidth=0.75,alpha=1,label='Proposed sample')
		
		ind=np.where( (fp['PTF PSP']==1) & (np.isnan(fp['Year t1'])==False) )[0]
		x=fp['X BC'][ind]; y=fp['Y BC'][ind]
		points=[]
		for k in range(x.size):
			points.append(Point(x[k],y[k]))
		gdf1=gpd.GeoDataFrame({'geometry':points,'ID':np.ones(x.size)})
		gdf1.crs=roi['crs']
		gdf1=gpd.overlay(gdf1,roi['gdf']['bound within'],how='intersection')
		gdf1.plot(ax=ax[3],marker='s',markersize=0.7,facecolor=[0.25,0.5,0.75],edgecolor=None,linewidth=0.75,alpha=1,label='Proposed sample')
	
	gu.axletters(ax,plt,0.0,0.9,FontColor=meta['Graphics']['gp']['cla'],LetterStyle='NoPar',FontWeight='Bold')
	
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_FieldPlotsFull','png',900)

	return

#%%
def Plot_RangeTenure(meta,roi,vnam):

	lab=['Forest with grazing tenure','Forest with haycutting tenure','Forest with no range tenure','Non-forest land','Non land']
	lab=lab[0:-1]
	z1=roi['grd'][vnam]['Data']
	z1[roi['grd']['Data']==0]=5

	# Number of colours and number of colours excluded from colorbar
	N_vis=4
	N_hidden=1
	N_tot=N_vis+N_hidden

	# Colormap
	cm=np.vstack( ((0.85,1,0.5,1),(0,0.5,0,1),(0.6,0.6,0.6,1),(0.83,0.83,0.83,1),(1,1,1,1)) )
	cm=matplotlib.colors.ListedColormap(cm)

	plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
	im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
	roi=PlotVectorBaseMaps(meta,roi,ax[0])
	ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'],xticklabels='',yticklabels='')
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

	zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
	cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
	cb=plt.colorbar(im,cax=ax[1],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
	ax[1].set(position=[0.71,0.6,0.05,0.14])
	cb.ax.set(yticklabels=lab)
	cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
	cb.outline.set_edgecolor('w')
	for i in range(cb_bnd.size):
		ax[1].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)

	pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
	ax[1].set(position=pos2)

	gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_' + vnam,'png',900)

	return fig,ax

#%%

def Plot_biomass_glob(meta,roi,vnam):

	# Grid
	bw=20; bin=np.arange(0,200+bw,bw);
	z1=(bin.size)*np.ones( roi['grd'][vnam]['Data'].shape)
	for i in range(bin.size):
		ind=np.where(np.abs( roi['grd'][vnam]['Data']-bin[i])<=bw/2)
		if ind[0].size>0:
			z1[ind]=i
		else:
			z1[0,i]=i
	ind=np.where(roi['grd'][vnam]['Data']>bin[i]); z1[ind]=i
	z1[(roi['grd']['Data']==1) & ( roi['grd'][vnam]['Data']==0)]=i+1
	z1[(roi['grd']['Data']!=1) | (roi['grd']['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Water'])]=i+2
	L=i+2

	lab=bin.astype(str)

	# Colormap
	cm=plt.cm.get_cmap('viridis',i)
	cm=np.vstack( (cm.colors,(0.9,0.9,0.9,1),(1,1,1,1)) )
	cm=matplotlib.colors.ListedColormap(cm)

	N_vis=bin.size+3
	N_hidden=3

	# Plot
	plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
	im=ax[0].matshow(z1,clim=(0,L+1),extent=roi['grd']['Extent'],cmap=cm)
	roi=PlotVectorBaseMaps(meta,roi,ax[0])
	ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

	cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_vis-(N_hidden-1),1),ticks=np.arange(0.5,N_vis-N_hidden,1))
	cb.ax.set(yticklabels=lab)
	cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
	cb.outline.set_edgecolor('w')
	for i in range(0,N_vis):
		ax[1].plot([0,100],[i,i],'w-',linewidth=1.5)

	pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
	ax[1].set(position=pos2)

	gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_' + vnam,'png',900)

	return fig,ax

#%%
def Plot_DistanceFrom(meta,roi,vnam):

	#z=roi['grd']['d2fac']['Data']
	z=roi['grd'][vnam]['Data']
	if roi['Type']=='ByTSA':
		z[roi['grd']['Data']==0]=0

	N_vis=10
	# Plot
	plt.close('all')
	fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
	im=ax[0].matshow(z,extent=roi['grd']['Extent'],cmap='viridis')

	#roi['gdf']['tsa'].plot(ax=ax[0],color=None,edgecolor=[0,0,0],facecolor='none')
	ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'],xticklabels='',yticklabels='')
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

	# Add relief shading
	#z=roi['grd']['elev']['Data']; dx,dy=roi['grd']['Cellsize'],roi['grd']['Cellsize']
	#ls=LightSource(azdeg=90,altdeg=45)
	#ve=0.1
	#hs=ls.hillshade(z,vert_exag=ve,dx=dx,dy=dy)
	#ax[0].matshow(hs,extent=roi['grd']['elev']['Extent'],cmap='Greys',alpha=0.4,clim=(np.min(hs),np.max(hs)))

	cb=plt.colorbar(im,cax=ax[1])#,boundaries=np.arange(0,N_vis-(N_hidden-1),1),ticks=np.arange(0.5,N_vis+1.5,1))
	pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
	ax[1].set(position=pos2)
	gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_' + vnam,'png',900)

	return fig,ax

#%%
def Plot_LU_Change_FromCEC(meta,roi,type):

	#type='def'
	if type=='def':
		lab=['Forest to Settlement','Forest to Cropland','Forest to Barren Ground','Land']# Deforestation
	elif type=='aff':
		lab=['Settlement to Forest','Cropland to Forest','Barren Ground to Forest','Land'] # Afforestation
	else:
		pass

	z0=roi['grd']['luc_' + type + '_cec']['Data']
	z1=z0
	z1[(z0==0)]=4
	z1[(roi['grd']['Data']!=1) | (roi['grd']['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Water'])]=5

	# Number of colours and number of colours excluded from colorbar
	N_vis=3
	N_hidden=2
	N_tot=N_vis+N_hidden

	# Colormap
	cm=np.vstack( ((0.5,0.0,0.0,1),(1.0,0,0,1),(1.0,0.5,0,1),(0.83,0.83,0.83,1),(1,1,1,1)) )
	cm=matplotlib.colors.ListedColormap(cm)

	plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
	im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
	if meta['Graphics']['Map']['Show Bound Within']=='On':
		roi['gdf']['bound within'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
	if meta['Graphics']['Map']['Show Bound Land Mask']=='On':
		roi['gdf']['bc_bound'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
	if meta['Graphics']['Map']['Show Lakes']=='On':
		roi['gdf']['lakes'].plot(ax=ax[0],facecolor=[0.82,0.88,1],label='Lakes',linewidth=2.25)
	if meta['Graphics']['Map']['Show Rivers']=='On':
		roi['gdf']['rivers'].plot(ax=ax[0],color=[0.6,0.8,1],label='Rivers',linewidth=0.25)
	if meta['Graphics']['Map']['Show Roads']=='On':
		roi['gdf']['road'].plot(ax=ax[0],edgecolor='k',linewidth=0.25,label='Road',alpha=1,zorder=1)
	ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'],xticklabels='',yticklabels='')
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

	zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
	cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
	cb=plt.colorbar(im,cax=ax[1],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
	ax[1].set(position=[0.71,0.6,0.05,0.14])
	cb.ax.set(yticklabels=lab)
	cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
	cb.outline.set_edgecolor('w')
	for i in range(cb_bnd.size):
		ax[1].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)

	pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
	ax[1].set(position=pos2)
	gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_lu_change_' + type + '_cec','png',900)

	return fig,ax

#%% Forest health outbreak maps
def Plot_OutbreakMap(meta,roi):
	pest_cd='IBS'
	tv=np.arange(1999,2022,1)
	zG=np.zeros(zRef['Data'].shape,dtype='int8')
	for iT in range(tv.size):
	
		z=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Disturbances\PEST_INFESTATION_POLY_' + pest_cd + '_SeverityClass_' + str(tv[iT]) + '.tif')
	
		zG[(z['Data']>1)]=1
	
		# Grid and labels
		z1=6*np.ones(zRef['Data'].shape)
		z1[(z['Data']==1)]=0
		z1[(z['Data']==2)]=1
		z1[(z['Data']==3)]=2
		z1[(z['Data']==4)]=3
		z1[(z['Data']==5)]=4
		z1[(z['Data']==0) & (zRef['Data']==1)]=6
		z1[(z['Data']==0) & (zG==1)]=5
		z1[(zRef['Data']==0)]=7
	
		ind=np.where(z1==7)
		for i in range(7):
			z1[ind[0][i],ind[1][i]]=i
		#z1[ind[0][1],ind[1][1]]=4
	
		lab=['Trace','Light','Moderate','Severe','Very severe','Previously affected','Unaffected','N/A']
	
		# Number of colours and number of colours excluded from colorbar
		N_vis=7
		N_hidden=1
		N_tot=N_vis+N_hidden
	
		# Colormap
		cm=np.vstack( ((1,1,0,1),(1,0.5,0,1),(1,0,0,1),(0.7,0,0,1),(0.35,0,0,1),(0.8,0.8,0.8,1),(0.93,0.93,0.93,1),(1,1,1,1)) )
		cm=matplotlib.colors.ListedColormap(cm)
	
		plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['figwidth'],(1-meta['Graphics']['sidespace'])*meta['Graphics']['figwidth']*zRef['yxrat']))
		im=ax[0].matshow(z1,extent=zRef['Extent'],cmap=cm)
		gdf['bc_bound']['gdf'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
		ax[0].set(position=meta['Graphics']['ax1 pos'],xlim=zRef['xlim'],ylim=zRef['ylim'],aspect='auto')
		ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both')
		ax[0].grid(meta['Graphics']['ax1 gridvis'])
		ax[0].axis(meta['Graphics']['ax1 vis'])
	
		cb_ivl=(N_tot-1)/N_tot
		cbivls=np.arange( cb_ivl , N_tot , cb_ivl)
		cbivls_low=np.arange( 0 , N_tot , cb_ivl)
	
		cb=plt.colorbar(im,cax=ax[1],cmap=cm,
						boundaries=np.arange(0,cbivls_low[N_vis],cb_ivl),
						ticks=np.arange(cb_ivl/2,N_tot-1,cb_ivl) )
		ax[1].set(position=[0.71,0.65,0.045,0.2])
	
		cb.ax.set(yticklabels=lab)
		cb.ax.tick_params(labelsize=11,length=0)
		cb.outline.set_edgecolor('w')
		for i in range(cbivls.size):
			ax[1].plot([0,100],[cbivls[i],cbivls[i]],'w-',linewidth=2)
	
		gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Beetles\Maps\Map_' + pest_cd + '_Severity_' + str(tv[iT]),'png',300)
	return

#%%
def Plot_FromModel(meta,roi,zRef,geos,md,vNam):
	
	z0=copy.deepcopy(zRef)
	z0['Data']=np.zeros(zRef['Data'].shape,dtype='int16')
	z0['Data'][geos['iMask']]=md[vNam]
	z0=gis.ClipToRaster(z0,roi['grd'])

	if vNam=='C_Biomass':
		bw=50; bin=np.arange(0,500+bw,bw); cm_nam='RdYlGn_r'
	elif (vNam=='C_Litter'):
		bw=20; bin=np.arange(0,200+bw,bw); cm_nam='RdYlGn_r'
	elif (vNam=='C_DeadWood'):
		bw=5; bin=np.arange(0,50+bw,bw); cm_nam='RdYlGn_r'
	elif (vNam=='C_Forest'):
		bw=50; bin=np.arange(0,700+bw,bw); cm_nam='viridis'
	elif (vNam=='C_ToMillTotal'):
		bw=5; bin=np.arange(0,50+bw,bw); cm_nam='RdYlGn_r'
	elif (vNam=='C_G_Gross'):
		bw=0.2; bin=np.arange(0,3+bw,bw); cm_nam='RdYlGn_r'
	elif (vNam=='E_Domestic_ForestSector_NEE'):
		z0['Data']=z0['Data']/3.667
		bw=0.5; bin=np.arange(-6,6+bw,bw); cm_nam='RdYlGn_r'
	elif (vNam=='E_NAB'):
		bw=2; bin=np.arange(-16,30+bw,bw); cm_nam='RdYlGn_r'

	z1=(bin.size)*np.ones(z0['Data'].shape)
	for i in range(bin.size):
		ind=np.where(np.abs(z0['Data']-bin[i])<=bw/2)
		z1[ind]=i
	ind=np.where(z0['Data']>bin[i]); z1[ind]=i
	z1[(z0['Data']==0) | (roi['grd']['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Water'])]=i+2
	#z1[(roi['grd']['Data']!=1)]=i+2
	L=i+2
	lab=np.round(bin,decimals=1).astype(str)

	# Colormap
	#cm=plt.cm.get_cmap('viridis',i)
	#cm=plt.cm.get_cmap('plasma',i)
	cm=plt.cm.get_cmap(cm_nam,i)
	cm.colors=cm(np.arange(0,cm.N))
	cm=np.vstack( (cm.colors,(0.85,0.85,0.85,1),(1,1,1,1)) )
	cm=matplotlib.colors.ListedColormap(cm)

	N_vis=bin.size+3
	N_hidden=3

	# Plot
	plt.close('all')
	fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
	im=ax[0].matshow(z1,clim=(0,L+1),extent=roi['grd']['Extent'],cmap=cm)
	roi=PlotVectorBaseMaps(meta,roi,ax[0])
	ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

	cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_vis-(N_hidden-1),1),ticks=np.arange(0.5,N_vis-N_hidden,1))
	cb.ax.set(yticklabels=lab)

	cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
	cb.outline.set_edgecolor('w')
	for i in range(0,N_vis):
		ax[1].plot([0,100],[i,i],'w-',linewidth=1.5)

	pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
	ax[1].set(position=pos2)
	gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_FromModel_' + vNam,'png',900)

	return fig,ax

#%%
def PlotWildfireBoundaries(meta,ax,roi,t0,t1):
	df=roi['gdf']['wf']['gdf']

	# Import names
	dNam=gu.ReadExcel(r'C:\Data\Wildfire\Wildfire Names\WildfireNames.xlsx')

	# Plot elevation
	fig,ax=Plot_Elev(meta,roi,'elev'); #del roi['grd'][vList[0]]

	# Add wildfirefire
	#YearStart=2017;ikp=np.where(df['FIRE_YEAR']==Year)[0]
	#YearStart=0;ikp=np.where(df['FIRE_YEAR']>=Year)[0]
	#YearStart=2017;
	ikp=np.where( (df['FIRE_YEAR']>=t0) & (df['FIRE_YEAR']<=t1) & (df['SHAPE_Area']/10000>50000) )[0]
	df0=df.iloc[ikp]
	df0.plot(ax=ax[0],facecolor=[1,0.4,0],edgecolor=[0.5,0.1,0],alpha=0.35,linewidth=0.5)

	# Add labels
	flg=1
	if flg==1:
		df0['coords']=df0['geometry'].apply(lambda x: x.representative_point().coords[:])
		df0['coords']=[coords[0] for coords in df0['coords']]
		for idx, row in df0.iterrows():
			ind=np.where(dNam['Code']==row['FIRE_NUMBER'])[0]
			try:
				ax[0].annotate(dNam['Name'][ind[0]],xy=row['coords'],horizontalalignment='center',color=[0.5,0.1,0],fontsize=6)
			except:
				ax[0].annotate(row['FIRE_NUMBER'],xy=row['coords'],horizontalalignment='center',color=[0.5,0.1,0],fontsize=6)

			#ax[0].annotate(row['FIRE_LABEL'],xy=row['coords'],horizontalalignment='center',color=[0,0,0],fontsize=6)

	gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_wildfire_bounds','png',900)
	return

#%%
def Plot_RoadsFromForestTenures(meta,roi):
	# Import data
	vList=['elev']; roi=u1ha.Import_Raster(meta,roi,vList)
	roi=u1ha.Import_GDB_Over_ROI(meta,roi,['road_ften'])
	# Plot
	meta['Graphics']['Map']['Show Roads']='On'
	fig,ax=p1ha.Plot_Elev(meta,roi,vList[0]);
	roi['gdf']['road_ften']['gdf'].plot(ax=ax[0],edgecolor=[1,1,0],facecolor='none',linewidth=0.25)
	gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_road_ften','png',900)
	return

#%%
def Plot3DTerrain(meta,roi):
	# Use subsampling of at least every other hectare to avoid being really slow and glitchy
	# vList=['elev']; roi=u1ha.Import_Raster(meta,roi,vList);
	ivl=3
	cm0=np.flip(np.vstack( ((0.32,0.19,0.19,1),(0.41,0.28,0.27,1),(0.49,0.38,0.36,1),(0.58,0.47,0.44,1),(0.66,0.57,0.53,1),(0.75,0.66,0.61,1),(0.83,0.76,0.7,1),(0.92,0.85,0.78,1),(1,0.95,0.87,1),(0.83,0.87,0.85,1),(0.67,0.78,0.82,1),(0.5,0.7,0.8,1),(0.33,0.61,0.78,1),(0.17,0.53,0.76,1),(0,0.45,0.74,1)) ),axis=0)
	cm0=matplotlib.colors.LinearSegmentedColormap.from_list('twi',cm0,N=30)
	#cmap=cm0(twi['Data']/np.amax(twi['Data']))
	z=roi['grd']['Data'][0::ivl,0::ivl]*roi['grd']['elev']['Data'][0::ivl,0::ivl]
	zmin=300#np.amin(roi['grd']['elev']['Data'][0::ivl,0::ivl])
	zmax=np.amax(z)
	cmap=cm0( (z-zmin)/(zmax-zmin) )
	
	e=60; a=360-2*45
	plt.close('all'); fig=plt.figure(figsize=gu.cm2inch(12,12))
	ax=fig.add_subplot(111, projection='3d')
	#ax=plt.axes(projection='3d')
	ax.plot_surface(roi['grd']['X'][0::ivl,0::ivl],roi['grd']['Y'][0::ivl,0::ivl],z,rstride=2,cstride=2,facecolors=cmap,linewidth=0,antialiased=False,shade=False) # ,
	ax.set(position=[0,0,1,1],xlim=[roi['grd']['xmin'],roi['grd']['xmax']],ylim=[roi['grd']['ymin'],roi['grd']['ymax']])
	ax.set_zlim(300,35000)
	ax.set_axis_on()
	ax.view_init(elev=e,azim=a)
	plt.tight_layout()
	ax.set(position=[-0.5,-0.2,2,2])
	return

#%%
def Plot_SiteVisits(meta,roi,vList):
	#vList=['lu_comp1_2019']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_LandUseComp1(meta,roi,vList[0]);
	vList=['elev']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_Elev(meta,roi,vList[0]);
	sites=gu.ReadExcel(r'C:\Data\Site Field Visits\Site Field Visits.xlsx')
	srs=gis.ImportSRSs()
	x,y=srs['Proj']['BC1ha'](sites['Lon'],sites['Lat'])
	points=[]
	for k in range(x.size):
		points.append(Point(x[k],y[k]))
	gdf_sites=gpd.GeoDataFrame({'geometry':points,'ID':np.ones(x.size)})
	gdf_sites.crs=roi['crs']
	gdf_sites['Name']=sites['Name']
	gdf_sites.plot(ax=ax[0],markersize=30,fc='r',ec='k',lw=0.5)
	for i in range(x.size):
		ax[0].annotate(sites['Name'][i],xy=[x[i]+1000,y[i]+1000],horizontalalignment='center',color=[0,0,0],fontsize=6)
	
	df.iloc[ikp].plot(ax=ax[0],facecolor='none',edgecolor=[0.6,0,0],linewidth=1)
	df['coords']=df['geometry'].apply(lambda x: x.representative_point().coords[:])
	df['coords']=[coords[0] for coords in df['coords']]
	
	roi['gdf']['road'].plot(ax=ax[0],edgecolor=[0.9,0.45,0],linewidth=0.5,label='Road',alpha=1,zorder=1)
	return

#%% Add watershed to water management basemap
def WaterManagementWithSelectWatersheds(meta,roi):
	vList=['elev']; roi=u1ha.Import_Raster(meta,roi,vList);
	fig,ax=p1ha.Plot_WaterManagement(meta,roi,vList[0]);
	lw=0.5

	roi=u1ha.Import_GDB_Over_ROI(meta,roi,['wshed8'])
	roi=u1ha.Import_GDB_Over_ROI(meta,roi,['wshed7'])

	# Select watersheds
	#orderL=[6]; roi=u1ha.Import_GDB_Over_ROI(meta,roi,['wshed'],OrderList=orderL)
	ind=np.where(roi['gdf']['wshed8']['gdf']['GNIS_NAME']=='Clearwater River')[0]
	roi['gdf']['wshed8']['gdf'].iloc[ind].plot(ax=ax[0],edgecolor=[0.4,0.6,0.2],facecolor=[0.8,1,0.4],alpha=0.35,linewidth=lw)
	ind=np.where(roi['gdf']['wshed7']['gdf']['GNIS_NAME']=='Bonaparte River')[0]
	roi['gdf']['wshed7']['gdf'].iloc[ind].plot(ax=ax[0],edgecolor=[0.4,0.6,0.2],facecolor=[0.8,1,0.4],alpha=0.35,linewidth=lw)

	ind=np.where(roi['gdf']['wshed6']['gdf']['GNIS_NAME']=='Deadman River')[0]
	roi['gdf']['wshed6']['gdf'].iloc[ind].plot(ax=ax[0],edgecolor=[0.4,0.6,0.2],facecolor=[0.8,1,0.4],alpha=0.35,linewidth=lw)

	ind=np.where(roi['gdf']['wshed6']['gdf']['GNIS_NAME']=='San Jose River')[0]
	roi['gdf']['wshed6']['gdf'].iloc[ind].plot(ax=ax[0],edgecolor=[0.4,0.6,0.2],facecolor=[0.8,1,0.4],alpha=0.35,linewidth=lw)

	ind=np.where(roi['gdf']['wshed6']['gdf']['GNIS_NAME']=='McKinley Creek')[0]
	roi['gdf']['wshed6']['gdf'].iloc[ind].plot(ax=ax[0],edgecolor=[0.4,0.6,0.2],facecolor=[0.8,1,0.4],alpha=0.35,linewidth=lw)

	ind=np.where(roi['gdf']['wshed5']['gdf']['GNIS_NAME']=='Dog Creek')[0]
	roi['gdf']['wshed5']['gdf'].iloc[ind].plot(ax=ax[0],edgecolor=[0.4,0.6,0.2],facecolor=[0.8,1,0.4],alpha=0.5,linewidth=1)
	ind=np.where(roi['gdf']['wshed5']['gdf']['GNIS_NAME']=='Indian Meadows Creek')[0]
	roi['gdf']['wshed5']['gdf'].iloc[ind].plot(ax=ax[0],edgecolor=[0.4,0.6,0.2],facecolor=[0.8,1,0.4],alpha=0.5,linewidth=1)
	ind=np.where(roi['gdf']['wshed5']['gdf']['GNIS_NAME']=='Canoe Creek')[0]
	roi['gdf']['wshed5']['gdf'].iloc[ind].plot(ax=ax[0],edgecolor=[0.4,0.6,0.2],facecolor=[0.8,1,0.4],alpha=0.5,linewidth=1)
	ind=np.where(roi['gdf']['wshed5']['gdf']['GNIS_NAME']=='Valenzuela Creek')[0]
	roi['gdf']['wshed5']['gdf'].iloc[ind].plot(ax=ax[0],edgecolor=[0.4,0.6,0.2],facecolor=[0.8,1,0.4],alpha=0.5,linewidth=1)
	ind=np.where(roi['gdf']['wshed5']['gdf']['GNIS_NAME']=='Big Bar Creek')[0]
	roi['gdf']['wshed5']['gdf'].iloc[ind].plot(ax=ax[0],edgecolor=[0.4,0.6,0.2],facecolor=[0.8,1,0.4],alpha=0.5,linewidth=1)
	ind=np.where(roi['gdf']['wshed5']['gdf']['GNIS_NAME']=='China Gulch')[0]
	roi['gdf']['wshed5']['gdf'].iloc[ind].plot(ax=ax[0],edgecolor=[0.4,0.6,0.2],facecolor=[0.8,1,0.4],alpha=0.5,linewidth=1)

	#ind=np.where(roi['gdf']['wshed']['gdf']['GNIS_NAME']=='Spius Creek')[0]
	#roi['gdf']['wshed']['gdf'].iloc[ind].plot(ax=ax[0],edgecolor=[0.4,0.6,0.2],facecolor=[0.8,1,0.4],alpha=0.5,linewidth=1)
	#ind=np.where(roi['gdf']['wshed']['gdf']['GNIS_NAME']=='Salmon River')[0]
	#roi['gdf']['wshed']['gdf'].iloc[ind].plot(ax=ax[0],edgecolor=[0.4,0.6,0.2],facecolor=[0.8,1,0.4],alpha=0.5,linewidth=1)
	#ind=np.where(roi['gdf']['wshed']['gdf']['GNIS_NAME']=='Canimred Creek')[0]
	#ind=np.where(roi['gdf']['wshed']['gdf']['GNIS_NAME']=='Deadman River')[0]
	#roi['gdf']['wshed']['gdf'].iloc[ind].plot(ax=ax[0],edgecolor=[0.4,0.6,0.2],facecolor=[0.8,1,0.4],alpha=0.5,linewidth=1)
	#ind=np.where(roi['gdf']['wshed']['gdf']['GNIS_NAME']=='Adams River')[0]
	#roi['gdf']['wshed']['gdf'].iloc[ind].plot(ax=ax[0],edgecolor=[0.4,0.6,0.2],facecolor=[0.8,1,0.4],alpha=0.5,linewidth=1)
	#ind=np.where(roi['gdf']['wshed']['gdf']['GNIS_NAME']=='Otter Creek')[0]
	#roi['gdf']['wshed']['gdf'].iloc[ind].plot(ax=ax[0],edgecolor=[0.4,0.6,0.2],facecolor=[0.8,1,0.4],alpha=0.5,linewidth=1)
	#ind=np.where(roi['gdf']['wshed']['gdf']['GNIS_NAME']=='Fishtrap Creek')[0]
	#roi['gdf']['wshed']['gdf'].iloc[ind].plot(ax=ax[0],edgecolor=[0.4,0.6,0.2],facecolor=[0.8,1,0.4],alpha=0.5,linewidth=1)

	#orderL=[4]; roi=u1ha.Import_GDB_Over_ROI(meta,roi,['wshed'],OrderList=orderL)
	#ind=np.where(roi['gdf']['wshed']['gdf']['GNIS_NAME']=='Pennask Creek')[0]
	#roi['gdf']['wshed']['gdf'].iloc[ind].plot(ax=ax[0],edgecolor=[0.4,0.6,0.2],facecolor=[0.8,1,0.4],alpha=0.5,linewidth=1)

	gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_ClearwaterRiver_' + nam,'png',900)
	return

#%%
def Plot_Map_FNM_ForTimeSpan(meta,roi,t0,t1):

	vnam='feca_yr'
	roi=u1ha.Import_Raster(meta,roi,[vnam])
	gdf=gpd.read_file(meta['Paths']['GDB']['GDB'] + '\\atu_FE_CA.geojson')

	z1=4*np.ones(roi['grd']['Data'].shape,dtype='int8')
	ind=np.where( (roi['grd']['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (roi['grd']['Data']>0) )
	z1[ind]=1
	ind=np.where( (roi['grd']['lc_comp1_2019']['Data']!=meta['LUT']['Derived']['lc_comp1']['Forest']) & (roi['grd']['Data']>0) )
	z1[ind]=2
	ind=np.where( (roi['grd'][vnam]['Data']>=t0) & (roi['grd'][vnam]['Data']<=t1) )
	z1[ind]=3

	if (meta['Graphics']['Map']['Show Inset Map']=='On'):
		N_panel=3
	else:
		N_panel=2

	# Number of colours and number of colours excluded from colorbar
	N_vis=3
	N_hidden=1
	N_tot=N_vis+N_hidden
	lab=np.array(['Forest','Non-forest','Treatment area'])

	cm=np.vstack( ((0.65,0.75,0.5,1),(0.9,0.9,0.9,1),(1,0.15,0.0,1),(1,1,1,1)) )
	cm=matplotlib.colors.ListedColormap(cm)

	plt.close('all'); fig,ax=plt.subplots(1,N_panel,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
	im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
	roi=PlotVectorBaseMaps(meta,roi,ax[0])
	ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])
	gdf[(gdf['Year']>=t0) & (gdf['Year']<=t1)].plot(ax=ax[0],facecolor='r',edgecolor='r',linewidth=1)

	zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
	cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
	cb=plt.colorbar(im,cax=ax[1],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
	ax[1].set(position=[meta['Graphics']['Map']['Legend X'],meta['Graphics']['Map']['Legend Y'],0.03,0.1])

	cb.ax.set(yticklabels=lab)
	cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
	cb.outline.set_edgecolor('w')
	for i in range(cb_bnd.size):
		ax[1].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)
	#ax[1].axis(meta['Graphics']['Map']['Map Axis Vis']) # Exclude colorbar

	# Scalebar
	if meta['Graphics']['Map']['Show Scalebar']=='On':
		AddScalebar(ax[0],roi)

	if meta['Graphics']['Map']['Show Inset Map']=='On':
		meta['gdf']['bc_bound']['gdf'].plot(ax=ax[2],edgecolor=None,facecolor=[0.8,0.8,0.8],linewidth=0.25)
		minx, miny, maxx, maxy = meta['gdf']['bc_bound']['gdf'].geometry.total_bounds
		roi['gdf']['bound within'].plot(ax=ax[2],edgecolor=None,facecolor=[0,0,0],linewidth=0.25)
		ax[2].set(position=[meta['Graphics']['Map']['Legend X'],0.1,1-meta['Graphics']['Map']['Legend X']-0.01,1-meta['Graphics']['Map']['Legend X']+0.05],xlim=[minx,maxx],ylim=[miny,maxy])
		ax[2].yaxis.set_ticks_position('both'); ax[2].xaxis.set_ticks_position('both'); ax[2].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[2].axis(meta['Graphics']['Map']['Map Axis Vis'])

	# *** Add Vector layers ***
	
	# points=[]
	# for k in range(x.size):
	#	 points.append(Point(x[k],y[k]))
	# gdf1=gpd.GeoDataFrame({'geometry':points,'ID':np.ones(x.size)})
	# gdf1.crs=roi['crs']
	# gdf1=gpd.overlay(gdf1,roi['gdf']['bound within'],how='intersection')
	# gdf1.plot(ax=ax[0],markersize=0.7,facecolor=[0.75,0,0],edgecolor=None,linewidth=0.75,alpha=1,label='Proposed sample')

	gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\Map_FNM_TreatmentAreas_' + str(t0) + '_to_' + str(t1),'png',900)
	return fig,ax

#%%
def Plot_Map_NOSE_ForTimeSpan(meta,roi,t0,t1):

	meta['Graphics']['Map']['Show Rail']='Off'
	meta['Graphics']['Map']['Show Roads']='Off'
	meta['Graphics']['Map']['Show Cities']='On'
	meta['Graphics']['Map']['Show Inset Map']='Off'

	# Import vector
	gdf=gpd.read_file(meta['Paths']['GDB']['GDB'] + '\\atu_PL_PL.geojson')

	# Import raster
	vnam='aset_no_all'
	roi=u1ha.Import_Raster(meta,roi,[vnam])

	z1=4*np.ones(roi['grd']['Data'].shape,dtype='int8')
	ind=np.where( (roi['grd']['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (roi['grd']['Data']>0) )
	z1[ind]=1
	ind=np.where( (roi['grd']['lc_comp1_2019']['Data']!=meta['LUT']['Derived']['lc_comp1']['Forest']) & (roi['grd']['Data']>0) )
	z1[ind]=2
	ind=np.where( (roi['grd'][vnam]['Data']>0) )
	z1[ind]=3

	if (meta['Graphics']['Map']['Show Inset Map']=='On'):
		N_panel=3
	else:
		N_panel=2

	# Number of colours and number of colours excluded from colorbar
	N_vis=3
	N_hidden=1
	N_tot=N_vis+N_hidden
	lab=np.array(['Forest','Non-forest','Treatment area'])

	cm=np.vstack( ((0.8,0.9,0.65,1),(0.9,0.9,0.9,1),(1,0,0.0,1),(1,1,1,1)) )
	cm=matplotlib.colors.ListedColormap(cm)

	plt.close('all'); fig,ax=plt.subplots(1,N_panel,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
	im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
	roi=PlotVectorBaseMaps(meta,roi,ax[0])
	ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

	gdf[(gdf['Year']>=t0) & (gdf['Year']<=t1)].plot(ax=ax[0],facecolor='r',edgecolor='r',linewidth=1)

	zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
	cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
	cb=plt.colorbar(im,cax=ax[1],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
	ax[1].set(position=[meta['Graphics']['Map']['Legend X'],meta['Graphics']['Map']['Legend Y'],0.03,0.1])

	cb.ax.set(yticklabels=lab)
	cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
	cb.outline.set_edgecolor('w')
	for i in range(cb_bnd.size):
		ax[1].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)
	#ax[1].axis(meta['Graphics']['Map']['Map Axis Vis']) # Exclude colorbar

	# Scalebar
	if meta['Graphics']['Map']['Show Scalebar']=='On':
		AddScalebar(ax[0],roi)

	if meta['Graphics']['Map']['Show Inset Map']=='On':
		meta['gdf']['bc_bound']['gdf'].plot(ax=ax[2],edgecolor=None,facecolor=[0.8,0.8,0.8],linewidth=0.25)
		minx,miny,maxx,maxy=meta['gdf']['bc_bound']['gdf'].geometry.total_bounds
		roi['gdf']['bound within'].plot(ax=ax[2],edgecolor=None,facecolor=[0,0,0],linewidth=0.25)
		ax[2].set(position=[meta['Graphics']['Map']['Legend X'],0.1,1-meta['Graphics']['Map']['Legend X']-0.01,1-meta['Graphics']['Map']['Legend X']+0.05],xlim=[minx,maxx],ylim=[miny,maxy])
		ax[2].yaxis.set_ticks_position('both'); ax[2].xaxis.set_ticks_position('both'); ax[2].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[2].axis(meta['Graphics']['Map']['Map Axis Vis'])

	gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_nose_' + str(t0) + 'to' + str(t1),'png',900)
	return fig,ax

#%%
def Plot_HoldridgeClimateZones(meta,roi,vnam):

	lut=gu.ReadExcel(r'G:\My Drive\Code_Python\fcgadgets\cbrunner\Parameters\LUT_HoldgridgeLifeZone.xlsx',sheet_name='Sheet1',skiprows=0)

	z0=roi['grd'][vnam]['Data']
	z0[(roi['grd']['Data']==0) | (roi['grd']['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Water'])]=0

	lab0=lut['Name']
	cl0=np.column_stack([lut['c1'],lut['c2'],lut['c3']])
	id0=lut['ID']

# 	uid=np.unique(z0[z0!=0])
# 	id1=np.zeros(uid.size)
# 	lab1=np.array(['' for _ in range(uid.size)],dtype=object)
# 	cl1=np.zeros((uid.size,3))
# 	z1=(uid.size+1)*np.ones(z0.shape)
# 	for i in range(uid.size):
# 		ind=np.where(z0==uid[i])
# 		z1[ind]=i+1
# 		ind=np.where(id0==uid[i])[0][0]
# 		id1[i]=id0[ind]
# 		lab1[i]=lab0[ind]
# 		cl1[i,:]=cl0[ind,:]

	uid=np.unique(z0[z0!=0])
	id1=np.zeros(uid.size)
	lab1=np.array(['' for _ in range(uid.size)],dtype=object)
	cl1=np.zeros((uid.size,3))
	z1=(uid.size+1)*np.ones(z0.shape)
	for i in range(uid.size):
		ind=np.where(z0==uid[i])
		z1[ind]=i+1
		ind=np.where(id0==uid[i])[0][0]
		id1[i]=id0[ind]
		lab1[i]=lab0[ind]
		cl1[i,:]=cl0[ind,:]

	#d=gu.CountByCategories(z1.flatten(),'Percent')

	N_vis=uid.size
	N_hidden=1
	N_tot=N_vis+N_hidden
	#lab1=np.append(lab1,[''])

	# Colormap
	cm=plt.cm.get_cmap('viridis',N_vis);
	for i in range(N_vis):
		cm.colors[i,0:3]=cl1[i,:]
	cm=np.vstack( (cm.colors,(1,1,1,1)) )
	cm=matplotlib.colors.ListedColormap(cm)

	if (meta['Graphics']['Map']['Show Inset Map']=='On'):
		N_panel=3
	else:
		N_panel=2

	plt.close('all'); fig,ax=plt.subplots(1,N_panel,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
	im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
	roi=PlotVectorBaseMaps(meta,roi,ax[0])
	ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

	# Legend
	zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
	cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
	cb=plt.colorbar(im,cax=ax[1],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)

# 	zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
# 	cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
# 	cb=plt.colorbar(im,cax=ax[1],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
	ax[1].set(position=[0.71,0.6,0.05,0.14])
	cb.ax.set(yticklabels=lab1)
	cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
	cb.outline.set_edgecolor('w')
	for i in range(cb_bnd.size):
		ax[1].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)
	pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
	ax[1].set(position=pos2)

	# Scalebar
	if meta['Graphics']['Map']['Show Scalebar']=='On':
		AddScalebar(ax[0],roi)

	# Inset
	if meta['Graphics']['Map']['Show Inset Map']=='On':
		meta['gdf']['bc_bound']['gdf'].plot(ax=ax[2],edgecolor=None,facecolor=[0.8,0.8,0.8],linewidth=0.25)
		minx,miny,maxx,maxy=meta['gdf']['bc_bound']['gdf'].geometry.total_bounds
		roi['gdf']['bound within'].plot(ax=ax[2],edgecolor=None,facecolor=[0,0,0],linewidth=0.25)
		# Lower right
		#pos=[meta['Graphics']['Map']['Legend X'],0.1,1-meta['Graphics']['Map']['Legend X']-0.01,1-meta['Graphics']['Map']['Legend X']+0.05]
		# Lower left
		pos=[0,0,1-meta['Graphics']['Map']['Legend X']-0.01,1-meta['Graphics']['Map']['Legend X']+0.05]
		ax[2].set(position=pos,xlim=[minx,maxx],ylim=[miny,maxy])
		ax[2].yaxis.set_ticks_position('both'); ax[2].xaxis.set_ticks_position('both'); ax[2].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[2].axis(meta['Graphics']['Map']['Map Axis Vis'])

	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_' + vnam,'png',900)

	return fig,ax

#%%