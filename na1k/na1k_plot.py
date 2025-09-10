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
import na1k.na1k_util as u1k
import fcgadgets.macgyver.util_general as gu
import fcgadgets.macgyver.util_gis as gis
import fcgadgets.macgyver.util_query_gdb as qgdb
import fcexplore.field_plots.Processing.fp_util as ufp
import fcgadgets.cbrunner.cbrun_util as cbu

#%%
def Plot_tmean_ann_Normal(meta,zRef):
	z0=gis.OpenGeoTiff(meta['Paths']['na1k'] + '\\Summaries\\Normals\\na1k_tmean_norm_1971to2000_ann.tif')['Data'].astype('float')*meta['SF']['tmean']
	#plt.close('all'); plt.matshow(z0,clim=[-5,15])
	bw=2; bin=np.arange(-16,20+bw,bw)
	N_vis=bin.size
	N_hidden=1
	N_tot=N_vis+N_hidden
	z1=N_vis*np.ones(z0.shape)
	for i in range(bin.size):
		ind=np.where(np.abs(z0-bin[i])<=bw/2)
		z1[ind]=i
	z1[(z0<=bin[0])]=1
	z1[(z0>=bin[i])]=i
	z1[(zRef['Data']==0)]=i+1
	for i in range(bin.size):
		z1[0,i]=i
	lab=bin.astype(str)

	#cm0=gu.ReadExcel(r'C:\Data\Colormaps\colormap_ws_zscore.xlsx')
	#cm=np.column_stack((cm0['cl1'],cm0['cl2'],cm0['cl3'],np.ones(cm0['bin'].size)))
	#cm=np.vstack( (cm,(1,1,1,1)) ) # (0.83137,0.81569,0.78431,1)
	#cm=matplotlib.colors.ListedColormap(cm)

	cm=plt.cm.get_cmap('jet',N_vis)
	cm.colors=cm(np.arange(0,cm.N))
	cm=np.vstack( (cm.colors,(1,1,1,1)) )
	cm=matplotlib.colors.ListedColormap(cm)

	plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*zRef['yxrat']))
	im=ax[0].matshow(z1,extent=zRef['Extent'],cmap=cm,clim=(0,N_tot)) #
	meta['pb'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
	#ax[0].text(-2.8e6,4.6e6,str(yr),fontweight='bold',fontsize=14)
	ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=zRef['xlim'],ylim=zRef['ylim'])
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])
	cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_vis+1,1),ticks=np.arange(0.5,N_vis,1))
	cb.ax.set(yticklabels=lab)
	cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
	cb.outline.set_edgecolor('w')
	for i in range(0,N_vis):
		ax[1].plot([0,100],[i,i],'w-',linewidth=1.5)
	ax[1].set(position=[0.01,0.01,0.03,0.65])
	gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\Summaries\\Normals\\na1k_tmean_norm_1971to2000_ann','png',900)
	return

#%%
def Plot_rswd_mjjas_Normal(meta,zRef):
	zRef=gis.OpenGeoTiff(meta['Paths']['na1k Ref Grid'])
	z0=gis.OpenGeoTiff(meta['Paths']['na1k'] + '\\Summaries\\Normals\\na1k_rswd_norm_1971to2000_mjjas.tif')['Data'].astype('float')*meta['SF']['rswd']
	#plt.close('all'); plt.matshow(z0,clim=[0,15])
	bw=1;
	bin=np.arange(12,27+bw,bw)
	N_vis=bin.size
	N_hidden=1
	N_tot=N_vis+N_hidden
	z1=N_vis*np.ones(z0.shape)
	for i in range(bin.size):
		ind=np.where(np.abs(z0-bin[i])<=bw/2)
		z1[ind]=i
	z1[(z0<=bin[0])]=1
	z1[(z0>=bin[i])]=i
	z1[(zRef['Data']==0)]=i+1
	for i in range(bin.size):
		z1[0,i]=i
	lab=bin.astype(str)

	cm=plt.cm.get_cmap('jet',N_vis)
	cm.colors=cm(np.arange(0,cm.N))
	cm=np.vstack( (cm.colors,(1,1,1,1)) )
	cm=matplotlib.colors.ListedColormap(cm)

	plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*zRef['yxrat']))
	im=ax[0].matshow(z1,extent=zRef['Extent'],cmap=cm,clim=(0,N_tot)) #
	meta['pb'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
	#ax[0].text(-2.8e6,4.6e6,str(yr),fontweight='bold',fontsize=14)
	ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=zRef['xlim'],ylim=zRef['ylim'])
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])
	cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_vis+1,1),ticks=np.arange(0.5,N_vis,1))
	cb.ax.set(yticklabels=lab)
	cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
	cb.outline.set_edgecolor('w')
	for i in range(0,N_vis):
		ax[1].plot([0,100],[i,i],'w-',linewidth=1.5)
	ax[1].set(position=[0.01,0.01,0.03,0.4])
	gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\Summaries\\Normals\\na1k_rswd_norm_1971to2000_mjjas','png',900)
	return

#%%
def Plot_prcp_ann_Normal(meta,zRef):
	zRef=gis.OpenGeoTiff(meta['Paths']['na1k Ref Grid'])
	DIM=30.5
	z0=gis.OpenGeoTiff(meta['Paths']['na1k'] + '\\Summaries\\Normals\\na1k_prcp_norm_1971to2000_ann.tif')['Data'].astype('float')*meta['SF']['prcp']/DIM
	#plt.close('all'); plt.matshow(z0,clim=[0,15])
	bw=0.5;
	bin=np.arange(0,6+bw,bw)
	N_vis=bin.size
	N_hidden=1
	N_tot=N_vis+N_hidden
	z1=N_vis*np.ones(z0.shape)
	for i in range(bin.size):
		ind=np.where(np.abs(z0-bin[i])<=bw/2)
		z1[ind]=i
	z1[(z0<=bin[0])]=1
	z1[(z0>=bin[i])]=i
	z1[(zRef['Data']==0)]=i+1
	for i in range(bin.size):
		z1[0,i]=i
	lab=np.round(bin,decimals=1).astype(str)

	cm=plt.cm.get_cmap('viridis',N_vis)
	cm.colors=cm(np.arange(0,cm.N))
	cm=np.vstack( (cm.colors,(1,1,1,1)) )
	cm=matplotlib.colors.ListedColormap(cm)

	plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*zRef['yxrat']))
	im=ax[0].matshow(z1,extent=zRef['Extent'],cmap=cm,clim=(0,N_tot)) #
	meta['pb'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
	#ax[0].text(-2.8e6,4.6e6,str(yr),fontweight='bold',fontsize=14)
	ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=zRef['xlim'],ylim=zRef['ylim'])
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])
	cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_vis+1,1),ticks=np.arange(0.5,N_vis,1))
	cb.ax.set(yticklabels=lab)
	cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
	cb.outline.set_edgecolor('w')
	for i in range(0,N_vis):
		ax[1].plot([0,100],[i,i],'w-',linewidth=1.5)
	ax[1].set(position=[0.01,0.01,0.03,0.5])
	gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\Summaries\\Normals\\na1k_prcp_norm_1971to2000_ann','png',900)
	return

#%%
def Plot_prcp_mjjas_Normal(meta,zRef):
	zRef=gis.OpenGeoTiff(meta['Paths']['na1k Ref Grid'])
	DIM=30.5
	z0=gis.OpenGeoTiff(meta['Paths']['na1k'] + '\\Summaries\\Normals\\na1k_prcp_norm_1971to2000_mjjas.tif')['Data'].astype('float')*meta['SF']['prcp']/DIM
	#plt.close('all'); plt.matshow(z0,clim=[0,30])
	bw=0.5;
	bin=np.arange(0,6+bw,bw)
	N_vis=bin.size
	N_hidden=1
	N_tot=N_vis+N_hidden
	z1=N_vis*np.ones(z0.shape)
	for i in range(bin.size):
		ind=np.where(np.abs(z0-bin[i])<=bw/2)
		z1[ind]=i
	z1[(z0<=bin[0])]=1
	z1[(z0>=bin[i])]=i
	z1[(zRef['Data']==0)]=i+1
	for i in range(bin.size):
		z1[0,i]=i
	lab=np.round(bin,decimals=1).astype(str)

	cm=plt.cm.get_cmap('viridis',N_vis)
	cm.colors=cm(np.arange(0,cm.N))
	cm=np.vstack( (cm.colors,(1,1,1,1)) )
	cm=matplotlib.colors.ListedColormap(cm)

	plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*zRef['yxrat']))
	im=ax[0].matshow(z1,extent=zRef['Extent'],cmap=cm,clim=(0,N_tot)) #
	meta['pb'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
	#ax[0].text(-2.8e6,4.6e6,str(yr),fontweight='bold',fontsize=14)
	ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=zRef['xlim'],ylim=zRef['ylim'])
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])
	cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_vis+1,1),ticks=np.arange(0.5,N_vis,1))
	cb.ax.set(yticklabels=lab)
	cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
	cb.outline.set_edgecolor('w')
	for i in range(0,N_vis):
		ax[1].plot([0,100],[i,i],'w-',linewidth=1.5)
	ax[1].set(position=[0.01,0.01,0.03,0.5])
	gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\Summaries\\Normals\\na1k_prcp_norm_1971to2000_mjjas','png',900)
	return

#%%
def Plot_ea_mjjas_Normal(meta,zRef):
	zRef=gis.OpenGeoTiff(meta['Paths']['na1k Ref Grid'])

	z0=gis.OpenGeoTiff(meta['Paths']['na1k'] + '\\Summaries\\Normals\\na1k_ea_biasadj_norm_1971to2000_mjjas.tif')['Data'].astype('float')*meta['SF']['ea']
	#plt.close('all'); plt.matshow(z0,clim=[0,15])
	bw=1;
	bin=np.arange(0,26+bw,bw)
	N_vis=bin.size
	N_hidden=1
	N_tot=N_vis+N_hidden
	z1=N_vis*np.ones(z0.shape)
	for i in range(bin.size):
		ind=np.where(np.abs(z0-bin[i])<=bw/2)
		z1[ind]=i
	z1[(z0<=bin[0])]=1
	z1[(z0>=bin[i])]=i
	z1[(zRef['Data']==0)]=i+1
	for i in range(bin.size):
		z1[0,i]=i
	lab=bin.astype(str)

	cm0=gu.ReadExcel(r'C:\Data\Colormaps\colormap_ea_norm_mjjas.xlsx')
	cm=np.column_stack((cm0['cl1'],cm0['cl2'],cm0['cl3'],np.ones(cm0['bin'].size)))
	cm=np.vstack( (cm,(1,1,1,1)) ) # (0.83137,0.81569,0.78431,1)
	cm=matplotlib.colors.ListedColormap(cm)

	plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*zRef['yxrat']))
	im=ax[0].matshow(z1,extent=zRef['Extent'],cmap=cm,clim=(0,N_tot)) #
	meta['pb'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
	#ax[0].text(-2.8e6,4.6e6,str(yr),fontweight='bold',fontsize=14)
	ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=zRef['xlim'],ylim=zRef['ylim'])
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])
	cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_vis+1,1),ticks=np.arange(0.5,N_vis,1))
	cb.ax.set(yticklabels=lab)
	cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
	cb.outline.set_edgecolor('w')
	for i in range(0,N_vis):
		ax[1].plot([0,100],[i,i],'w-',linewidth=1.5)
	ax[1].set(position=[0.01,0.01,0.03,0.7])
	gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\Summaries\\Normals\\na1k_ea_biasadj_norm_1971to2000_mjjas','png',900)
	return

#%%
def Plot_vpd_mjjas_Normal(meta,zRef):
	zRef=gis.OpenGeoTiff(meta['Paths']['na1k Ref Grid'])
	z0=gis.OpenGeoTiff(meta['Paths']['na1k'] + '\\Summaries\\Normals\\na1k_vpd_norm_1971to2000_mjjas.tif')['Data'].astype('float')*meta['SF']['vpd']
	#plt.close('all'); plt.matshow(z0,clim=[0,15])
	bw=1;
	bin=np.arange(0,26+bw,bw)
	N_vis=bin.size
	N_hidden=1
	N_tot=N_vis+N_hidden
	z1=N_vis*np.ones(z0.shape)
	for i in range(bin.size):
		ind=np.where(np.abs(z0-bin[i])<=bw/2)
		z1[ind]=i
	z1[(z0<=bin[0])]=1
	z1[(z0>=bin[i])]=i
	z1[(zRef['Data']==0)]=i+1
	for i in range(bin.size):
		z1[0,i]=i
	lab=bin.astype(str)

	cm=plt.cm.get_cmap('jet',N_vis)
	cm.colors=cm(np.arange(0,cm.N))
	cm=np.vstack( (cm.colors,(1,1,1,1)) )
	cm=matplotlib.colors.ListedColormap(cm)

	plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*zRef['yxrat']))
	im=ax[0].matshow(z1,extent=zRef['Extent'],cmap=cm,clim=(0,N_tot)) #
	meta['pb'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
	#ax[0].text(-2.8e6,4.6e6,str(yr),fontweight='bold',fontsize=14)
	ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=zRef['xlim'],ylim=zRef['ylim'])
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])
	cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_vis+1,1),ticks=np.arange(0.5,N_vis,1))
	cb.ax.set(yticklabels=lab)
	cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
	cb.outline.set_edgecolor('w')
	for i in range(0,N_vis):
		ax[1].plot([0,100],[i,i],'w-',linewidth=1.5)
	ax[1].set(position=[0.01,0.01,0.03,0.65])
	gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\Summaries\\Normals\\na1k_vpd_norm_1971to2000_mjjas','png',900)
	return

#%%
def Plot_etp_mjjas_Normal(meta,zRef):
	zRef=gis.OpenGeoTiff(meta['Paths']['na1k Ref Grid'])
	DIM=30.5
	z0=gis.OpenGeoTiff(meta['Paths']['na1k'] + '\\Summaries\\Normals\\na1k_etp_norm_1971to2000_mjjas.tif')['Data'].astype('float')*meta['SF']['etp']/DIM
	#plt.close('all'); plt.matshow(z0,clim=[0,15])
	bw=0.5;
	bin=np.arange(0,13+bw,bw)
	N_vis=bin.size
	N_hidden=1
	N_tot=N_vis+N_hidden
	z1=N_vis*np.ones(z0.shape)
	for i in range(bin.size):
		ind=np.where(np.abs(z0-bin[i])<=bw/2)
		z1[ind]=i
	z1[(z0<=bin[0])]=1
	z1[(z0>=bin[i])]=i
	z1[(zRef['Data']==0)]=i+1
	for i in range(bin.size):
		z1[0,i]=i
	lab=bin.astype(str)

	#cm0=gu.ReadExcel(r'C:\Data\Colormaps\colormap_etp_norm_mjjas.xlsx')
	cm0=gu.ReadExcel(r'C:\Data\Colormaps\colormap_ea_norm_mjjas.xlsx')
	cm=np.column_stack((cm0['cl1'],cm0['cl2'],cm0['cl3'],np.ones(cm0['bin'].size)))
	cm=np.vstack( (cm,(1,1,1,1)) ) # (0.83137,0.81569,0.78431,1)
	cm=matplotlib.colors.ListedColormap(cm)

	plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*zRef['yxrat']))
	im=ax[0].matshow(z1,extent=zRef['Extent'],cmap=cm,clim=(0,N_tot)) #
	meta['pb'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
	#ax[0].text(-2.8e6,4.6e6,str(yr),fontweight='bold',fontsize=14)
	ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=zRef['xlim'],ylim=zRef['ylim'])
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])
	cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_vis+1,1),ticks=np.arange(0.5,N_vis,1))
	cb.ax.set(yticklabels=lab)
	cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
	cb.outline.set_edgecolor('w')
	for i in range(0,N_vis):
		ax[1].plot([0,100],[i,i],'w-',linewidth=1.5)
	ax[1].set(position=[0.01,0.01,0.03,0.7])
	gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\Summaries\\Normals\\na1k_etp_norm_1971to2000_mjjas','png',900)
	return

#%%
def Plot_runoff_ann_Normal(meta,zRef):
	zRef=gis.OpenGeoTiff(meta['Paths']['na1k Ref Grid'])
	DIM=30.5
	z0=gis.OpenGeoTiff(meta['Paths']['na1k'] + '\\Summaries\\Normals\\na1k_runoff_norm_1971to2000_ann.tif')['Data'].astype('float')*meta['SF']['runoff']/DIM
	#plt.close('all'); plt.matshow(z0,clim=[0,6])
	bw=1;
	bin=np.arange(0,15+bw,bw)
	N_vis=bin.size
	N_hidden=1
	N_tot=N_vis+N_hidden
	z1=N_vis*np.ones(z0.shape)
	for i in range(bin.size):
		ind=np.where(np.abs(z0-bin[i])<=bw/2)
		z1[ind]=i
	z1[(z0<=bin[0])]=1
	z1[(z0>=bin[i])]=i
	z1[(zRef['Data']==0)]=i+1
	for i in range(bin.size):
		z1[0,i]=i
	lab=bin.astype(str)

	cm=plt.cm.get_cmap('Blues',N_vis)
	cm.colors=cm(np.arange(0,cm.N))
	cm=np.vstack( (cm.colors,(1,1,1,1)) )
	cm=matplotlib.colors.ListedColormap(cm)

	plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*zRef['yxrat']))
	im=ax[0].matshow(z1,extent=zRef['Extent'],cmap=cm,clim=(0,N_tot)) #
	meta['pb'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
	#ax[0].text(-2.8e6,4.6e6,str(yr),fontweight='bold',fontsize=14)
	ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=zRef['xlim'],ylim=zRef['ylim'])
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])
	cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_vis+1,1),ticks=np.arange(0.5,N_vis,1))
	cb.ax.set(yticklabels=lab)
	cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
	cb.outline.set_edgecolor('w')
	for i in range(0,N_vis):
		ax[1].plot([0,100],[i,i],'w-',linewidth=1.5)
	ax[1].set(position=[0.01,0.01,0.03,0.44])
	gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\Summaries\\Normals\\na1k_runoff_norm_1971to2000_ann','png',900)
	return

#%%
def Plot_melt_jfmam_Normal(meta,zRef):
	zRef=gis.OpenGeoTiff(meta['Paths']['na1k Ref Grid'])
	DIM=30.5
	z0=gis.OpenGeoTiff(meta['Paths']['na1k'] + '\\Summaries\\Normals\\na1k_melt_norm_1971to2000_jfmam.tif')['Data'].astype('float')*meta['SF']['melt']/DIM
	#plt.close('all'); plt.matshow(z0,clim=[0,6])
	bw=1;
	bin=np.arange(0,15+bw,bw)
	N_vis=bin.size
	N_hidden=1
	N_tot=N_vis+N_hidden
	z1=N_vis*np.ones(z0.shape)
	for i in range(bin.size):
		ind=np.where(np.abs(z0-bin[i])<=bw/2)
		z1[ind]=i
	z1[(z0<=bin[0])]=1
	z1[(z0>=bin[i])]=i
	z1[(zRef['Data']==0)]=i+1
	for i in range(bin.size):
		z1[0,i]=i
	lab=np.round(bin,decimals=1).astype(str)

	cm=plt.cm.get_cmap('Blues',N_vis)
	cm.colors=cm(np.arange(0,cm.N))
	cm=np.vstack( (cm.colors,(1,1,1,1)) )
	cm=matplotlib.colors.ListedColormap(cm)

	plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*zRef['yxrat']))
	im=ax[0].matshow(z1,extent=zRef['Extent'],cmap=cm,clim=(0,N_tot)) #
	meta['pb'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
	#ax[0].text(-2.8e6,4.6e6,str(yr),fontweight='bold',fontsize=14)
	ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=zRef['xlim'],ylim=zRef['ylim'])
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])
	cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_vis+1,1),ticks=np.arange(0.5,N_vis,1))
	cb.ax.set(yticklabels=lab)
	cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
	cb.outline.set_edgecolor('w')
	for i in range(0,N_vis):
		ax[1].plot([0,100],[i,i],'w-',linewidth=1.5)
	ax[1].set(position=[0.01,0.01,0.03,0.44])
	gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\Summaries\\Normals\\na1k_melt_norm_1971to2000_jfmam','png',900)
	return

#%%
def Plot_cwd_mjjas_Normal(meta,zRef):
	zRef=gis.OpenGeoTiff(meta['Paths']['na1k Ref Grid'])
	DIM=30.5
	z0=gis.OpenGeoTiff(meta['Paths']['na1k'] + '\\Summaries\\Normals\\na1k_cwd_norm_1971to2000_mjjas.tif')['Data'].astype('float')*meta['SF']['cwd']/DIM
	#plt.close('all'); plt.matshow(z0,clim=[0,12])
	bw=0.5;
	bin=np.arange(0,10+bw,bw)
	N_vis=bin.size
	N_hidden=1
	N_tot=N_vis+N_hidden
	z1=N_vis*np.ones(z0.shape)
	for i in range(bin.size):
		ind=np.where(np.abs(z0-bin[i])<=bw/2)
		z1[ind]=i
	z1[(z0<=bin[0])]=1
	z1[(z0>=bin[i])]=i
	z1[(zRef['Data']==0)]=i+1
	for i in range(bin.size):
		z1[0,i]=i
	lab=bin.astype(str)

	cm0=gu.ReadExcel(r'C:\Data\Colormaps\colormap_cwd_norm_mjjas.xlsx')
	cm=np.column_stack((cm0['cl1'],cm0['cl2'],cm0['cl3'],np.ones(cm0['bin'].size)))
	cm=np.vstack( (cm,(1,1,1,1)) ) # (0.83137,0.81569,0.78431,1)
	cm=matplotlib.colors.ListedColormap(cm)

	plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*zRef['yxrat']))
	im=ax[0].matshow(z1,extent=zRef['Extent'],cmap=cm,clim=(0,N_tot)) #
	meta['pb'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
	#ax[0].text(-2.8e6,4.6e6,str(yr),fontweight='bold',fontsize=14)
	ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=zRef['xlim'],ylim=zRef['ylim'])
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])
	cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_vis+1,1),ticks=np.arange(0.5,N_vis,1))
	cb.ax.set(yticklabels=lab)
	cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
	cb.outline.set_edgecolor('w')
	for i in range(0,N_vis):
		ax[1].plot([0,100],[i,i],'w-',linewidth=1.5)
	ax[1].set(position=[0.01,0.01,0.03,0.65])
	gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\Summaries\\Normals\\na1k_cwd_norm_1971to2000_mjjas','png',900)
	return

#%%
def Plot_cmi_mjjas_Normal(meta,zRef):
	zRef=gis.OpenGeoTiff(meta['Paths']['na1k Ref Grid'])
	DIM=30.5
	z0=gis.OpenGeoTiff(meta['Paths']['na1k'] + '\\Summaries\\Normals\\na1k_cmi_norm_1971to2000_mjjas.tif')['Data'].astype('float')*meta['SF']['cmi']/DIM
	#plt.close('all'); plt.matshow(z0,clim=[0,12])
	bw=0.5;
	bin=np.arange(-6,6+bw,bw)
	N_vis=bin.size
	N_hidden=1
	N_tot=N_vis+N_hidden
	z1=N_vis*np.ones(z0.shape)
	for i in range(bin.size):
		ind=np.where(np.abs(z0-bin[i])<=bw/2)
		z1[ind]=i
	z1[(z0<=bin[0])]=1
	z1[(z0>=bin[i])]=i
	z1[(zRef['Data']==0)]=i+1
	for i in range(bin.size):
		z1[0,i]=i
	lab=bin.astype(str)

# 	cm0=gu.ReadExcel(r'C:\Data\Colormaps\colormap_cwd_norm_mjjas.xlsx')
# 	cm=np.column_stack((cm0['cl1'],cm0['cl2'],cm0['cl3'],np.ones(cm0['bin'].size)))
# 	cm=np.vstack( (cm,(1,1,1,1)) ) # (0.83137,0.81569,0.78431,1)
# 	cm=matplotlib.colors.ListedColormap(cm)
	cm=plt.cm.get_cmap('RdBu',N_vis)
	cm.colors=cm(np.arange(0,cm.N))
	cm=np.vstack( (cm.colors,(1,1,1,1)) )
	cm=matplotlib.colors.ListedColormap(cm)

	plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*zRef['yxrat']))
	im=ax[0].matshow(z1,extent=zRef['Extent'],cmap=cm,clim=(0,N_tot)) #
	meta['pb'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
	#ax[0].text(-2.8e6,4.6e6,str(yr),fontweight='bold',fontsize=14)
	ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=zRef['xlim'],ylim=zRef['ylim'])
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])
	cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_vis+1,1),ticks=np.arange(0.5,N_vis,1))
	cb.ax.set(yticklabels=lab)
	cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
	cb.outline.set_edgecolor('w')
	for i in range(0,N_vis):
		ax[1].plot([0,100],[i,i],'w-',linewidth=1.5)
	ax[1].set(position=[0.01,0.01,0.03,0.65])
	gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\Summaries\\Normals\\na1k_cmi_norm_1971to2000_mjjas','png',900)
	return

#%%
def Plot_wsp_djf_Normal(meta,zRef):
	zRef=gis.OpenGeoTiff(meta['Paths']['na1k Ref Grid'])
	z0=gis.OpenGeoTiff(meta['Paths']['na1k'] + '\\Summaries\\Normals\\na1k_wsp_norm_1971to2000_djf.tif')['Data'].astype('float')*meta['SF']['wsp']
	#plt.close('all'); plt.matshow(z0,clim=[0,300])
	bw=20;
	bin=np.arange(0,340+bw,bw)
	N_vis=bin.size
	N_hidden=1
	N_tot=N_vis+N_hidden
	z1=N_vis*np.ones(z0.shape)
	for i in range(bin.size):
		ind=np.where(np.abs(z0-bin[i])<=bw/2)
		z1[ind]=i
	z1[(z0<=bin[0])]=0
	z1[(z0>=bin[i])]=i
	z1[(zRef['Data']==0)]=i+1
	for i in range(bin.size):
		z1[0,i]=i
	lab=bin.astype(str)

	#cm0=gu.ReadExcel(r'C:\Data\Colormaps\colormap_etp_norm_mjjas.xlsx')
	cm0=gu.ReadExcel(r'C:\Data\Colormaps\colormap_wsp.xlsx')
	cm=np.column_stack((cm0['cl1'],cm0['cl2'],cm0['cl3'],np.ones(cm0['bin'].size)))
	cm=np.vstack( (cm,(1,1,1,1)) ) # (0.83137,0.81569,0.78431,1)
	cm=matplotlib.colors.ListedColormap(cm)

	plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*zRef['yxrat']))
	im=ax[0].matshow(z1,extent=zRef['Extent'],cmap=cm,clim=(0,N_tot)) #
	meta['pb'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
	#ax[0].text(-2.8e6,4.6e6,str(yr),fontweight='bold',fontsize=14)
	ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=zRef['xlim'],ylim=zRef['ylim'])
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])
	cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_vis+1,1),ticks=np.arange(0.5,N_vis,1))
	cb.ax.set(yticklabels=lab)
	cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
	cb.outline.set_edgecolor('w')
	for i in range(0,N_vis):
		ax[1].plot([0,100],[i,i],'w-',linewidth=1.5)
	ax[1].set(position=[0.01,0.01,0.03,0.44])
	gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\Summaries\\Normals\\na1k_wsp_norm_1971to2000_djf','png',900)
	return

#%%
def Plot_ws_mjjas_Normal(meta,zRef):
	zRef=gis.OpenGeoTiff(meta['Paths']['na1k Ref Grid'])
	z0=gis.OpenGeoTiff(meta['Paths']['na1k'] + '\\Summaries\\Normals\\na1k_ws_norm_1971to2000_mjjas.tif')['Data'].astype('float')*meta['SF']['ws']
	bw=10;
	bin=np.arange(0,200+bw,bw)
	N_vis=bin.size
	N_hidden=1
	N_tot=N_vis+N_hidden
	z1=N_vis*np.ones(z0.shape)
	for i in range(bin.size):
		ind=np.where(np.abs(z0-bin[i])<=bw/2)
		z1[ind]=i
	z1[(z0<=bin[0])]=1
	z1[(z0>=bin[i])]=i
	z1[(zRef['Data']==0)]=i+1
	for i in range(bin.size):
		z1[0,i]=i
	lab=bin.astype(str)

	cm0=gu.ReadExcel(r'C:\Data\Colormaps\colormap_ws.xlsx')
	cm=np.column_stack((cm0['cl1'],cm0['cl2'],cm0['cl3'],np.ones(cm0['bin'].size)))
	cm=np.vstack( (cm,(1,1,1,1)) ) # (0.83137,0.81569,0.78431,1)
	cm=matplotlib.colors.ListedColormap(cm)

	plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*zRef['yxrat']))
	im=ax[0].matshow(z1,extent=zRef['Extent'],cmap=cm,clim=(0,N_tot)) #
	meta['pb'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
	#ax[0].text(-2.8e6,4.6e6,str(yr),fontweight='bold',fontsize=14)
	ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=zRef['xlim'],ylim=zRef['ylim'])
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])
	cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_vis+1,1),ticks=np.arange(0.5,N_vis,1))
	cb.ax.set(yticklabels=lab)
	cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
	cb.outline.set_edgecolor('w')
	for i in range(0,N_vis):
		ax[1].plot([0,100],[i,i],'w-',linewidth=1.5)
	ax[1].set(position=[0.01,0.01,0.03,0.65])
	gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\Summaries\\Normals\\na1k_ws_norm_1971to2000_mjjas','png',900)
	return

#%%
def Plot_tmean_anom(meta,zRef):
	for yr in range(1998,2024):
		z0=gis.OpenGeoTiff(meta['Paths']['na1k'] + '\\Summaries\\Time Series\\Anomalies\\na1k_tmean_anom_mjjas_' + str(yr) + '.tif')['Data'].astype('float')*meta['SF']['tmean']
		#plt.close('all'); plt.matshow(z0,clim=[-5,5])
		bw=0.5; bin=np.arange(-4,4+bw,bw)
		N_vis=bin.size
		N_hidden=1
		N_tot=N_vis+N_hidden
		z1=N_vis*np.ones(z0.shape)
		for i in range(bin.size):
			ind=np.where(np.abs(z0-bin[i])<=bw/2)
			z1[ind]=i
		z1[(z0<=bin[0])]=1
		z1[(z0>=bin[i])]=i
		z1[(zRef['Data']==0)]=i+1
		for i in range(bin.size):
			z1[0,i]=i
		lab=np.round(bin,decimals=2).astype(str)

		cm=plt.cm.get_cmap('RdYlBu_r',N_vis)
		cm.colors=cm(np.arange(0,cm.N))
		cm=np.vstack( (cm.colors,(1,1,1,1)) )
		cm=matplotlib.colors.ListedColormap(cm)
		#cm0=gu.ReadExcel(r'C:\Data\Colormaps\colormap_ws_zscore.xlsx')
		#cm=np.column_stack((cm0['cl1'],cm0['cl2'],cm0['cl3'],np.ones(cm0['bin'].size)))
		#cm=np.vstack( (cm,(1,1,1,1)) ) # (0.83137,0.81569,0.78431,1)
		#cm=matplotlib.colors.ListedColormap(cm)

		plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*zRef['yxrat']))
		im=ax[0].matshow(z1,extent=zRef['Extent'],cmap=cm,clim=(0,N_tot)) #
		meta['pb'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
		ax[0].text(-2.8e6,4.6e6,str(yr),fontweight='bold',fontsize=14)
		ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=zRef['xlim'],ylim=zRef['ylim'])
		ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])
		cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_vis+1,1),ticks=np.arange(0.5,N_vis,1))
		cb.ax.set(yticklabels=lab)
		cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
		cb.outline.set_edgecolor('w')
		for i in range(0,N_vis):
			ax[1].plot([0,100],[i,i],'w-',linewidth=1.5)
		ax[1].set(position=[0.01,0.01,0.03,0.5])
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\Summaries\\Time Series\\Anomalies\\na1k_tmean_anom_mjjas_' + str(yr),'png',500)
	return

#%%
def Plot_prcp_anom(meta,zRef):
	DIM=30.5
	for yr in range(1998,2024):
		z0=gis.OpenGeoTiff(meta['Paths']['na1k'] + '\\Summaries\\Time Series\\Anomalies\\na1k_prcp_anom_mjjas_' + str(yr) + '.tif')['Data'].astype('float')*meta['SF']['prcp']/DIM
		#plt.close('all'); plt.matshow(z0,clim=[-2,2])
		bw=0.25; bin=np.arange(-2.5,2.5+bw,bw)
		N_vis=bin.size
		N_hidden=1
		N_tot=N_vis+N_hidden
		z1=N_vis*np.ones(z0.shape)
		for i in range(bin.size):
			ind=np.where(np.abs(z0-bin[i])<=bw/2)
			z1[ind]=i
		z1[(z0<=bin[0])]=1
		z1[(z0>=bin[i])]=i
		z1[(zRef['Data']==0)]=i+1
		for i in range(bin.size):
			z1[0,i]=i
		lab=np.round(bin,decimals=2).astype(str)

		cm=plt.cm.get_cmap('RdYlGn',N_vis)
		cm.colors=cm(np.arange(0,cm.N))
		cm=np.vstack( (cm.colors,(1,1,1,1)) )
		cm=matplotlib.colors.ListedColormap(cm)
		#cm0=gu.ReadExcel(r'C:\Data\Colormaps\colormap_ws_zscore.xlsx')
		#cm=np.column_stack((cm0['cl1'],cm0['cl2'],cm0['cl3'],np.ones(cm0['bin'].size)))
		#cm=np.vstack( (cm,(1,1,1,1)) ) # (0.83137,0.81569,0.78431,1)
		#cm=matplotlib.colors.ListedColormap(cm)

		plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*zRef['yxrat']))
		im=ax[0].matshow(z1,extent=zRef['Extent'],cmap=cm,clim=(0,N_tot)) #
		meta['pb'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
		ax[0].text(-2.8e6,4.6e6,str(yr),fontweight='bold',fontsize=14)
		ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=zRef['xlim'],ylim=zRef['ylim'])
		ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])
		cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_vis+1,1),ticks=np.arange(0.5,N_vis,1))
		cb.ax.set(yticklabels=lab)
		cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
		cb.outline.set_edgecolor('w')
		for i in range(0,N_vis):
			ax[1].plot([0,100],[i,i],'w-',linewidth=1.5)
		ax[1].set(position=[0.01,0.01,0.03,0.65])
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\Summaries\\Time Series\\Anomalies\\na1k_prcp_anom_mjjas_' + str(yr),'png',500)
	return

#%%
def Plot_etp_anom(meta,zRef):
	DIM=30.5
	for yr in range(1998,2024):
		z0=gis.OpenGeoTiff(meta['Paths']['na1k'] + '\\Summaries\\Time Series\\Anomalies\\na1k_etp_anom_mjjas_' + str(yr) + '.tif')['Data'].astype('float')*meta['SF']['etp']/DIM
		#plt.close('all'); plt.matshow(z0,clim=[-2,2])
		bw=0.5; bin=np.arange(-4,4+bw,bw)
		N_vis=bin.size
		N_hidden=1
		N_tot=N_vis+N_hidden
		z1=N_vis*np.ones(z0.shape)
		for i in range(bin.size):
			ind=np.where(np.abs(z0-bin[i])<=bw/2)
			z1[ind]=i
		z1[(z0<=bin[0])]=1
		z1[(z0>=bin[i])]=i
		z1[(zRef['Data']==0)]=i+1
		for i in range(bin.size):
			z1[0,i]=i
		lab=np.round(bin,decimals=2).astype(str)

		cm=plt.cm.get_cmap('RdYlBu_r',N_vis)
		cm.colors=cm(np.arange(0,cm.N))
		cm=np.vstack( (cm.colors,(1,1,1,1)) )
		cm=matplotlib.colors.ListedColormap(cm)
		#cm0=gu.ReadExcel(r'C:\Data\Colormaps\colormap_ws_zscore.xlsx')
		#cm=np.column_stack((cm0['cl1'],cm0['cl2'],cm0['cl3'],np.ones(cm0['bin'].size)))
		#cm=np.vstack( (cm,(1,1,1,1)) ) # (0.83137,0.81569,0.78431,1)
		#cm=matplotlib.colors.ListedColormap(cm)

		plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*zRef['yxrat']))
		im=ax[0].matshow(z1,extent=zRef['Extent'],cmap=cm,clim=(0,N_tot)) #
		meta['pb'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
		ax[0].text(-2.8e6,4.6e6,str(yr),fontweight='bold',fontsize=14)
		ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=zRef['xlim'],ylim=zRef['ylim'])
		ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])
		cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_vis+1,1),ticks=np.arange(0.5,N_vis,1))
		cb.ax.set(yticklabels=lab)
		cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
		cb.outline.set_edgecolor('w')
		for i in range(0,N_vis):
			ax[1].plot([0,100],[i,i],'w-',linewidth=1.5)
		ax[1].set(position=[0.01,0.01,0.03,0.65])
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\Summaries\\Time Series\\Anomalies\\na1k_etp_anom_mjjas_' + str(yr),'png',900)
	return

#%%
def Plot_cmi_anom(meta,zRef):
	DIM=30.5
	for yr in range(1998,2024):
		z0=gis.OpenGeoTiff(meta['Paths']['na1k'] + '\\Summaries\\Time Series\\Anomalies\\na1k_cmi_anom_mjjas_' + str(yr) + '.tif')['Data'].astype('float')*meta['SF']['cmi']/DIM
		#plt.close('all'); plt.matshow(z0,clim=[-2,2])
		bw=0.5; bin=np.arange(-4,4+bw,bw)
		N_vis=bin.size
		N_hidden=1
		N_tot=N_vis+N_hidden
		z1=N_vis*np.ones(z0.shape)
		for i in range(bin.size):
			ind=np.where(np.abs(z0-bin[i])<=bw/2)
			z1[ind]=i
		z1[(z0<=bin[0])]=1
		z1[(z0>=bin[i])]=i
		z1[(zRef['Data']==0)]=i+1
		for i in range(bin.size):
			z1[0,i]=i
		lab=np.round(bin,decimals=2).astype(str)

		cm=plt.cm.get_cmap('RdYlGn',N_vis)
		cm.colors=cm(np.arange(0,cm.N))
		cm=np.vstack( (cm.colors,(1,1,1,1)) )
		cm=matplotlib.colors.ListedColormap(cm)
		#cm0=gu.ReadExcel(r'C:\Data\Colormaps\colormap_ws_zscore.xlsx')
		#cm=np.column_stack((cm0['cl1'],cm0['cl2'],cm0['cl3'],np.ones(cm0['bin'].size)))
		#cm=np.vstack( (cm,(1,1,1,1)) ) # (0.83137,0.81569,0.78431,1)
		#cm=matplotlib.colors.ListedColormap(cm)

		plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*zRef['yxrat']))
		im=ax[0].matshow(z1,extent=zRef['Extent'],cmap=cm,clim=(0,N_tot)) #
		meta['pb'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
		ax[0].text(-2.8e6,4.6e6,str(yr),fontweight='bold',fontsize=14)
		ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=zRef['xlim'],ylim=zRef['ylim'])
		ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])
		cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_vis+1,1),ticks=np.arange(0.5,N_vis,1))
		cb.ax.set(yticklabels=lab)
		cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
		cb.outline.set_edgecolor('w')
		for i in range(0,N_vis):
			ax[1].plot([0,100],[i,i],'w-',linewidth=1.5)
		ax[1].set(position=[0.01,0.01,0.03,0.65])
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\Summaries\\Time Series\\Anomalies\\na1k_cmi_anom_mjjas_' + str(yr),'png',500)
	return

#%%
def Plot_ws_zscore(meta,zRef):
	for yr in range(1998,2024):
		z0=gis.OpenGeoTiff(meta['Paths']['na1k'] + '\\Summaries\\Time Series\\Zscores\\na1k_ws_zscore_mjjas_' + str(yr) + '.tif')['Data'].astype('float')*meta['SF']['ws']
		#plt.close('all'); plt.matshow(z0,clim=[-5,5])
		bw=0.5; bin=np.arange(-5,5+bw,bw)
		N_vis=bin.size
		N_hidden=1
		N_tot=N_vis+N_hidden
		z1=N_vis*np.ones(z0.shape)
		for i in range(bin.size):
			ind=np.where(np.abs(z0-bin[i])<=bw/2)
			z1[ind]=i
		z1[(z0<=bin[0])]=1
		z1[(z0>=bin[i])]=i
		z1[(zRef['Data']==0)]=i+1
		for i in range(bin.size):
			z1[0,i]=i
		lab=bin.astype(str)
		cm0=gu.ReadExcel(r'C:\Data\Colormaps\colormap_ws_zscore.xlsx')
		cm=np.column_stack((cm0['cl1'],cm0['cl2'],cm0['cl3'],np.ones(cm0['bin'].size)))
		cm=np.vstack( (cm,(1,1,1,1)) ) # (0.83137,0.81569,0.78431,1)
		cm=matplotlib.colors.ListedColormap(cm)
		plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*zRef['yxrat']))
		im=ax[0].matshow(z1,extent=zRef['Extent'],cmap=cm,clim=(0,N_tot)) #
		meta['pb'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
		ax[0].text(-2.8e6,4.6e6,str(yr),fontweight='bold',fontsize=14)
		ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=zRef['xlim'],ylim=zRef['ylim'])
		ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])
		cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_vis+1,1),ticks=np.arange(0.5,N_vis,1))
		cb.ax.set(yticklabels=lab)
		cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
		cb.outline.set_edgecolor('w')
		for i in range(0,N_vis):
			ax[1].plot([0,100],[i,i],'w-',linewidth=1.5)
		ax[1].set(position=[0.01,0.01,0.03,0.65])
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\Summaries\\Time Series\\Zscores\\na1k_ws_zscore_mjjas_' + str(yr),'png',500)
	return


#%%
def Plot_HoldridgeClimateZones(meta,roi,vnam):

	lut=gu.ReadExcel(r'G:\My Drive\Code_Python\fcgadgets\cbrunner\Parameters\LUT_HoldgridgeLifeZone.xlsx',sheet_name='Sheet1',skiprows=0)

	z0=gis.OpenGeoTiff(meta['Paths']['na1k'] + '\\Summaries\\Normals\\na1k_HoldridgeLifeZones_1971to2000.tif')['Data']
	lab0=lut['Name']
	cl0=np.column_stack([lut['c1'],lut['c2'],lut['c3']])
	id0=lut['ID']

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
	#lab1=np.append([''],lab1)
	#np.unique(z1)

	# Colormap
	cm=plt.cm.get_cmap('viridis',N_vis);
	for i in range(N_vis):
		cm.colors[i,0:3]=cl1[i,:]
	cm=np.vstack( (cm.colors,(1,1,1,1)) )
	cm=matplotlib.colors.ListedColormap(cm)

	plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*zRef['yxrat']))
	im=ax[0].matshow(z1,extent=zRef['Extent'],cmap=cm)
	meta['pb'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
	ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=zRef['xlim'],ylim=zRef['ylim'])
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])
	zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
	cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
	cb=plt.colorbar(im,cax=ax[1],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
	ax[1].set(position=[0.71,0.6,0.05,0.14])
	cb.ax.set(yticklabels=lab1)
	cb.ax.tick_params(labelsize=5,length=0)
	cb.outline.set_edgecolor('w')
	for i in range(cb_bnd.size):
		ax[1].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)
	ax[1].set(position=[0.001,0.001,0.03,0.3])

	gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\Summaries\\na1k_HoldridgeLifeZones','png',500)

	return