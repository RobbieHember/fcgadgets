#%% Import modules
import os
import numpy as np
import time
import matplotlib.pyplot as plt
import pandas as pd
import fiona
import fcgadgets.macgyver.util_general as gu
import fcgadgets.macgyver.util_gis as gis
import fcgadgets.glob50k.glob50k_util as ug50
gp=gu.SetGraphics('Manuscript')

#%% Import paths and look-up-tables
meta,gdf=ug50.Init()

#%% Land cover class
#z=gis.OpenGeoTiff(meta['Paths']['glob50k'] + '\\Downloads\\lcc_copern.tif')
z=gis.OpenGeoTiff(meta['Paths']['glob50k'] + '\\lcc1.tif')
#plt.close('all'); plt.matshow(z['Data'],clim=[0,100])

#%% Forest Land Mask
z1=copy.deepcopy(z)
z1['Data']=3*np.ones(z['Data'].shape,dtype='int8')
ind=np.where( (z['Data']>0) ); z1['Data'][ind]=1
ind=np.where( (z['Data']==80) | (z['Data']==200) ); z1['Data'][ind]=3
ind=np.where( (z['Data']>=111) & (z['Data']<=116) ); z1['Data'][ind]=2
ind=np.where( (z['Data']>=121) & (z['Data']<=126) ); z1['Data'][ind]=2
#plt.close('all'); plt.matshow(z1['Data'],clim=[0,2])
gis.SaveGeoTiff(z1,meta['Paths']['glob50k'] + '\\ForestLandMask.tif')
print(np.where(z1['Data']==2)[0].size/1500)

#%% Plot

N_vis=2
N_hidden=1
N_tot=N_vis+N_hidden
lab=['Non-forest','Forest']

# Colormap
cm=np.vstack( ((0.8,0.8,0.8,1),(0.4,0.8,0,1),(1,1,1,1)) )
cm=matplotlib.colors.ListedColormap(cm)

plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*meta['Graphics']['Map']['yxrat']))
im=ax[0].matshow(z1['Data'],extent=z1['Extent'],cmap=cm)
gdf.plot(ax=ax[0],edgecolor=[0,0,0],facecolor='None',lw=0.25,alpha=1)
ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=z1['xlim'],ylim=z1['ylim'],xticklabels='',yticklabels='')
ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

zmn=np.min(z1['Data'])
zmx=np.max(z1['Data'])
cb_ivl=(zmx-zmn)/N_tot
cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
cb=plt.colorbar(im,cax=ax[1],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
#ax[1].set(position=[0.71,0.6,0.05,0.14])
cb.ax.set(yticklabels=lab)
cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
cb.outline.set_edgecolor('w')
for i in range(cb_bnd.size):
	ax[1].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)
pos2=[meta['Graphics']['Map']['Legend X'],0.45-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
ax[1].set(position=pos2)
#if meta['Graphics']['Print Figures']=='On':
#	gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_' + nam,'png',900)

#%%
