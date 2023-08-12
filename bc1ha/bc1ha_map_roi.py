
#%% Import modules

import numpy as np
import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib.colors
from matplotlib.colors import LightSource
from shapely.geometry import Polygon,Point
import copy
import time
import fcgadgets.macgyver.util_general as gu
import fcgadgets.macgyver.util_gis as gis
import fcgadgets.bc1ha.bc1ha_util as u1ha
import fcgadgets.macgyver.util_query_gdb as qgdb

#%% Import parameters

meta=u1ha.Init()
meta['Graphics']['Map']['RGSF']=1
meta['Graphics']['Map']['Fig Width']=15.5
meta['Graphics']['Map']['Side Space']=0.25
meta['Graphics']['Map']['Map Position']=[0,0,1-meta['Graphics']['Map']['Side Space']-0.01,1]
meta['Graphics']['Map']['Map Axis Vis']='off'
meta['Graphics']['Map']['Map Grid Vis']=False
meta['Graphics']['Map']['Legend X']=1-meta['Graphics']['Map']['Side Space']+meta['Graphics']['Map']['Map Position'][0]#+0.01,0.6,0.03,0.35]
meta['Graphics']['Map']['Legend Width']=0.0275
meta['Graphics']['Map']['Legend Font Size']=7
meta['Graphics']['Map']['Legend Text Space']=0.035
meta['Graphics']['Map']['Show Bound Land Mask']='On'
meta['Graphics']['Map']['Show Bound Within']='Off'
meta['Graphics']['Map']['Show Lakes']='On'
meta['Graphics']['Map']['Show Rivers']='On'
meta['Graphics']['Map']['Show Roads']='On'

meta['Graphics']['Plot Style']='Manuscript'
meta['Graphics']['gp']=gu.SetGraphics(meta['Graphics']['Plot Style'])
meta['Graphics']['Print Figures']='On'

#%% Import base maps

gdf=u1ha.Import_GDBs_ProvinceWide(meta)

#%% Define region of interest

# Identify how ROI is defined
roi={}

roi['Type']='Prov'
#roi['Type']='ByRegDis'
#roi['Type']='ByTSA'
#roi['Type']='ByLatLon'

#roi['Name']='Arrowsmith TSA'
#roi['Name']='Fort St John TSA'
#roi['Name']='Boundary TSA'
#roi['Name']='Cassiar TSA'
#roi['Name']='Merritt TSA'
#roi['Name']='Fort Nelson TSA'
#roi['Name']='Prince George TSA'
#roi['Name']='Williams Lake TSA'
#roi['Name']='North Island TSA'
#roi['Name']='Quesnel TSA'
#roi['Name']='100 Mile House TSA'

roi['Name']='Prov'
#roi['Name']='CAPITAL'
#roi['Name']='STRATHCONA'
#roi['Name']='COWICHAN VALLEY'

if roi['Type']=='ByTSA':
    meta['Graphics']['Print Figure Path']=r'C:\Users\rhember\OneDrive - Government of BC\Figures' + '\\TSAs\\' + roi['Name']
elif roi['Type']=='ByRegDis':
    meta['Graphics']['Print Figure Path']=r'C:\Users\rhember\OneDrive - Government of BC\Figures' + '\\Regional Districts\\' + roi['Name']
elif roi['Type']=='Prov':
    meta['Graphics']['Print Figure Path']=r'C:\Users\rhember\OneDrive - Government of BC\Figures' + '\\Province\\' + roi['Name']

t0=time.time()
if roi['Type']=='ByTSA':
    # Search: gdf['tsa']['key']
    #roi['List']=['Arrowsmith TSA']
    #roi['List']=['Fort St. John TSA']
    #roi['List']=['Boundary TSA']
    #roi['List']=['Cassiar TSA']
    #roi['List']=['Fort Nelson TSA']
    #roi['List']=['Merritt TSA']
    #roi['List']=['Prince George TSA']
    #roi['List']=['Kamloops TSA','100 Mile House TSA','Williams Lake TSA']
    #roi['List']=['Williams Lake TSA']
    #roi['List']=['North Island TSA']
    #roi['List']=['Quesnel TSA']
    #roi['List']=['100 Mile House TSA']
    #roi['List']=['Kootenay Lake TSA']
    roi['List']=['100 Mile House TSA']
    #roi['List']=['Okanagan TSA']
    #roi['List']=['Merritt TSA','Kamloops TSA','100 Mile House TSA','Okanagan TSA','Williams Lake TSA','Lillooet TSA','Boundary TSA'] # ,'Arrow TSA','Revelstoke TSA'
    #roi['List']=list(gdf['tsa']['key']['Name'])

    # Western spruce budworm study (do not change!)
    #roi['List']=['Merritt TSA','Kamloops TSA','100 Mile House TSA','Okanagan TSA','Williams Lake TSA','Lillooet TSA','Boundary TSA'] # ,'Arrow TSA','Revelstoke TSA'

elif roi['Type']=='ByRegDis':
    # Search: gdf['regdis']['gdf']['REGIONAL_DISTRICT_NAME'].unique()
    roi['List']=['CAPITAL']
    #roi['List']=['COWICHAN VALLEY']
    #roi['List']=['STRATHCONA']

elif roi['Type']=='ByLatLon':
    # Hanceville fire
    flg=0
    if flg==1:
        roi['Centre']=[-122.92,51.92]
        roi['Radius']=40*1000
        meta['Graphics']['Print Figure Path']=r'C:\Users\rhember\OneDrive - Government of BC\Figures\Reforestation\Reforestation Hanceville ROI'

    # Elephant Hill fire
    flg=0
    if flg==1:
        roi['Centre']=[-121.15,51.15]
        roi['Radius']=45*1000 # metres
        meta['Graphics']['Print Figure Path']=r'C:\Users\rhember\OneDrive - Government of BC\Figures\Reforestation\Reforestation Elephant Hill ROI'

    # Yung Lake, Elephant Hill fire
    flg=0
    if flg==1:
        roi['Centre']=[-121.07,51.308]
        roi['Radius']=10*1000 # metres

    # Yahk
    flg=0
    if flg==1:
        #roi['Centre']=[-116.086296,49.106646]
        #roi['Radius']=100*1000 # metres
        roi['Lower Left']=[-117.5231,48.9926]
        roi['Upper Right']=[-115.5786,49.6713]

    # IWB
    flg=0
    if flg==1:
        roi['Centre']=[-117.5,50.8]
        roi['Radius']=100*1000 # metres

    # LICS
    flg=1
    if flg==1:
        roi['Name']='ROI_EX'
        lics=gu.ipickle(r'D:\Data\FCI_Projects\LICS Site List\LICS Sites.pkl')
        id=5; roi['Centre']=[lics[id]['Lon'],lics[id]['Lat']]; roi['Name']=lics[id]['Loc']
        roi['Radius']=15*1000 # metres

elif roi['Type']=='Prov':
    pass

# Prepare region of interest
# Roads is empty because it takes a long time
roi=u1ha.DefineROI(meta,roi,gdf)
t1=time.time()
print((t1-t0)/60)

#%% Simplify

gdf['road']['gdf'].to_file(r'C:\Users\rhember\Documents\Data\Geodatabases\test1.geojson', driver="GeoJSON")  
a=gdf['road']['gdf'].simplify(100)
a.to_file(r'C:\Users\rhember\Documents\Data\Geodatabases\test2.geojson', driver="GeoJSON")  

#%% Import rasters over ROI

vList=u1ha.GetRasterListFromSpreadsheet(r'C:\Users\rhember\Documents\Data\BC1ha\RasterInclusion.xlsx')
roi=u1ha.Import_Raster(meta,roi,vList)

#%% Import required vector geodatabases

#vList=['cc','fcres']
#vList=['wf'] # 'op','cc','fcres','ogsr'
#roi=u1ha.Import_GDB_Over_ROI(meta,roi,vList)

#%% Resource roads

#gdf['rd']={}
#gdf['rd']['gdf']=gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='FTEN_ROAD_SEGMENT_LINES_SVW')
#gdf['rd']['gdf']=gpd.overlay(gdf['rd']['gdf'],roi['gdf']['bound within'],how='intersection')

#%% Plot FECA year

def Plot_FECA_Year(meta,roi):

    z0=roi['grd']['feca_yr']['Data']
    bw=5; bin=np.arange(1975,2025+bw,bw);
    z1=(bin.size)*np.ones( z0.shape)
    for i in range(bin.size):
        ind=np.where(np.abs(z0-bin[i])<=bw/2)
        z1[ind]=i
    z1[(z0==0)]=i+1
    z1[(roi['grd']['lcc1_c']['Data']==meta['LUT']['Derived']['lcc1']['Water'])]=i+2

    N_vis=bin.size+3
    N_hidden=3
    N_tot=N_vis+N_hidden

    lab=bin.astype(str)

    #cm=plt.cm.get_cmap('viridis',N_vis)
    cm=plt.cm.get_cmap('plasma',N_vis)
    cm=np.vstack( (cm.colors,(0,0,0,1),(1,1,1,1)) )
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
    gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_feca_yr','png',900)
    return fig,ax

fig,ax=Plot_FECA_Year(meta,roi)

#%% Land Cover Class 1

fig,ax=u1ha.Plot_LandCoverClass1(meta,roi)

#%%

u1ha.Plot_REARs(meta,roi)

#%% Plot BGC Zones

fig,ax=u1ha.Plot_BGC_Zone(meta,roi)

#%% Species leading NTEMS

def Plot_Spc1(meta,roi):

    # Compress categories
    z0=roi['grd']['spc1_ntems']['Data']
    lab0=np.array(list(meta['LUT']['Derived']['spc1_ntems'].keys()))
    id0=np.array(list(meta['LUT']['Derived']['spc1_ntems'].values()))
    cl0=np.column_stack([meta['LUT']['Raw']['spc1_ntems']['R'],meta['LUT']['Raw']['spc1_ntems']['G'],meta['LUT']['Raw']['spc1_ntems']['B']])
    z1,lab1,cl1=gis.CompressCats(z0,id0,lab0,cl0)
    z1=z1+1

    n=int(np.max(z1))
    z1[0,0]=n+1
    z1[(roi['grd']['Data']==0) | (roi['grd']['lcc1_c']['Data']==meta['LUT']['Derived']['lcc1']['Water'])]=n+2

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

    gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_spc1_ntems','png',900)
    plt.show()
    return fig,ax

fig,ax=Plot_Spc1(meta,roi)

#%% Plot mean annual temp

fig,ax=u1ha.Plot_MAT(meta,roi)

#%% Plot mean annual precip

fig,ax=u1ha.Plot_MAP(meta,roi)

#%% PLOT age from VRI

def Plot_PROJ_AGE_1(meta,roi):

    z0=roi['grd']['PROJ_AGE_1']['Data']

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
    ind=np.where(roi['grd']['PROJ_AGE_1']['Data']>bin[i]); z1[ind]=id[-1]
    z1[(roi['grd']['PROJ_AGE_1']['Data']==0) | (roi['grd']['lcc1_c']['Data']!=meta['LUT']['Derived']['lcc1']['Forest'])]=id[-1]+1
    z1[(roi['grd']['Data']==0) | (roi['grd']['lcc1_c']['Data']==0) | (roi['grd']['lcc1_c']['Data']==meta['LUT']['Derived']['lcc1']['Water'])]=id[-1]+2

    lab=['']*(N_tot-1)
    lab[0:N_vis]=bin.astype(str)

    cm=plt.cm.get_cmap('viridis',N_vis)
    cm=np.vstack( (cm.colors,(0.9,0.9,0.9,1),(1,1,1,1)) )
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

    gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_PROJ_AGE_1','png',900)

    return fig,ax

fig,ax=Plot_PROJ_AGE_1(meta,roi)

#%% Plot Site Index from VRI

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
    z1[(z0==0) | (roi['grd']['lcc1_c']['Data']!=meta['LUT']['Derived']['lcc1']['Forest'])]=id[-1]+1
    #z1[1,1]=id[-1]+1
    z1[(roi['grd']['Data']==0) | (roi['grd']['lcc1_c']['Data']==0) | (roi['grd']['lcc1_c']['Data']==meta['LUT']['Derived']['lcc1']['Water'])]=id[-1]+2

    lab=['']*(N_vis+1)
    lab[0:N_vis]=bin.astype(str)

    cm=plt.cm.get_cmap('viridis',N_vis)
    cm=np.vstack( (cm.colors,(0.9,0.9,0.9,1),(1,1,1,1)) )
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

fig,ax=Plot_SI(meta,roi)



#%% PLOT age from NTEMS

def Plot_Age_NTEMS(meta,roi):

    z0=roi['grd']['age_ntem']['Data']

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
    z1[(z0==0) | (roi['grd']['lcc1_c']['Data']!=meta['LUT']['Derived']['lcc1']['Forest'])]=id[-1]+1
    z1[(roi['grd']['Data']==0) | (roi['grd']['lcc1_c']['Data']==0) | (roi['grd']['lcc1_c']['Data']==meta['LUT']['Derived']['lcc1']['Water'])]=id[-1]+2

    lab=['']*(N_tot-1)
    lab[0:N_vis]=bin.astype(str)

    cm=plt.cm.get_cmap('viridis',N_vis)
    cm=np.vstack( (cm.colors,(0.9,0.9,0.9,1),(1,1,1,1)) )
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

    gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_Age_NTEM','png',900)

    return fig,ax

fig,ax=Plot_Age_NTEMS(meta,roi)

#%%

def Plot_RegenType(meta,roi,gdf):

    lab=['harvested and planted','Harvested (pre-1988), not planted','Harvested (1988-2018), not planted','Harvested (post-2018), not planted','Straight planting','No harvesting / no planting']
    z1=7*np.ones(roi['grd']['Data'].shape,dtype='int8')
    for i in range(1,7):
        z1[(roi['grd']['regentype']['Data']==i)]=i
    ind=np.where( (roi['grd']['Data']==0) ); z1[ind]=7

    N_vis=len(lab)
    N_hidden=1
    N_tot=N_vis+N_hidden
    cm=np.vstack( ((0.85,0.7,1,1),(1,1,0,1),(1,0.25,0,1),(0.7,0.7,0.7,1),(0.7,1,0.2,1),(0.9,0.9,0.9,1),(1,1,1,1)) )
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
    #gdf['cities'].plot(ax=ax[0],marker='s',edgecolor=[0.75,0.3,0],facecolor=[1,0.6,0],lw=1,markersize=20,alpha=1,zorder=2)
    #for x,y,label in zip(gdf['cities'].geometry.x,gdf['cities'].geometry.y,gdf['cities'].Name):
    #    ax[0].annotate(label,xy=(x,y),xytext=(4,3),textcoords="offset points",color=[0.75,0.3,0],fontsize=4)
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
    gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_regentype','png',900)

    return fig,ax

Plot_RegenType(meta,roi,gdf)

#%% Plot planted area

def Plot_PlantedMask(meta,roi):

    lab=[]
    z1=2*np.ones(roi['grd']['plam']['Data'].shape)
    z1[(roi['grd']['Data']==1) & (roi['grd']['plam']['Data']>0)]=1; lab.append('Planted')
    z1[(roi['grd']['Data']==1) & (roi['grd']['plam']['Data']==0)]=2; lab.append('Not planted')
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

    gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_planted','png',900)

    return fig,ax

fig,ax=Plot_PlantedMask(meta,roi)

#%% Get Fuel Treatments

fre={}
fre['Path']=r'C:\Users\rhember\Documents\Data\ForestInventory\Results\20230430\Results.gdb'
fre['Layer']='RSLT_ACTIVITY_TREATMENT_SVW'; # fiona.listlayers(fre['Path'])
fre['crs']=gdf['bc_bound']['gdf'].crs
fre['Keep Geom']='On'
fre['Select Openings']=np.array([])
fre['SBC']=np.array([])
fre['FSC']=np.array([])
fre['SOC1']=np.array(['FRE'])
fre['ROI']=[]
fre['gdf']=qgdb.Query_Openings(fre,roi)

fre={}
fre['Path']=r'C:\Users\rhember\Documents\Data\ForestInventory\Results\20230430\Results.gdb'
fre['Layer']='RSLT_ACTIVITY_TREATMENT_SVW'; # fiona.listlayers(fre['Path'])
fre['crs']=gdf['bc_bound']['gdf'].crs
fre['Keep Geom']='On'
fre['Select Openings']=np.array([])
fre['SBC']=np.array([])
fre['FSC']=np.array([])
fre['SOC1']=np.array(['FRE'])
fre['ROI']=[]
fre['gdf']=qgdb.Query_Openings(fre,roi)

#%% Plot elevation (and fire)

def Plot_Elev(meta,roi):

    z=roi['grd']['elev']['Data']
    if (roi['Type']=='ByTSA') | (roi['Type']=='ByRegDis'):
        z[roi['grd']['Data']==0]=0

    plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
    im=ax[0].matshow(roi['grd']['elev']['Data'],extent=roi['grd']['Extent'],cmap='Greys')
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

    # Plot fuel treatments
    #fre['gdf'].plot(ax=ax[0],facecolor='None',edgecolor=[1,0.5,0],linewidth=1.5,label='Opening',alpha=1)

    # Plot lumber mills, pulp mills and chipper mills
    #mtypeL=['LBR','PLP','PLT']
    #mtypeL=list(gdf['tpf']['gdf']['PRODUCT_CODE'].unique())
    #ax=u1ha.PlotMills(ax,gdf['tpf']['gdf'],mtypeL,labels='On')

    #gdf['cities'].plot(ax=ax[0],marker='s',edgecolor=[0.75,0.3,0],facecolor=[1,0.6,0],lw=1,markersize=20,alpha=1,zorder=2)
    #for x,y,label in zip(gdf['cities'].geometry.x,gdf['cities'].geometry.y,gdf['cities'].Name):
    #    ax[0].annotate(label,xy=(x,y),xytext=(4,3),textcoords="offset points",color=[0.75,0.3,0],fontsize=4)

    roi['gdf']['popp'].plot(ax=ax[0],marker='s',edgecolor=None,facecolor=[0.9,0.45,0],lw=0.5,markersize=3,alpha=1,zorder=2)
    for x,y,label in zip(roi['gdf']['popp'].geometry.x,roi['gdf']['popp'].geometry.y,roi['gdf']['popp'].NAME):
        ax[0].annotate(label,xy=(x,y),xytext=(2,1.5),textcoords="offset points",color=[0.9,0.45,0],fontsize=4)

    ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'],xticklabels='',yticklabels='')
    ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

    # Add relief shading
    #z=roi['grd']['elev']['Data']; dx,dy=roi['grd']['Cellsize'],roi['grd']['Cellsize']
    #ls=LightSource(azdeg=90,altdeg=45)
    #ve=0.1
    #hs=ls.hillshade(z,vert_exag=ve,dx=dx,dy=dy)
    #ax[0].matshow(hs,extent=roi['grd']['elev']['Extent'],cmap='Greys',alpha=0.4,clim=(np.min(hs),np.max(hs)))

    plt.colorbar(im,cax=ax[1])#,boundaries=np.arange(0,N_vis-(N_hidden-1),1),ticks=np.arange(0.5,N_vis+1.5,1))
    pos=[meta['Graphics']['Map']['Legend X'],0.5,meta['Graphics']['Map']['Legend Width'],0.4]
    ax[1].set(position=pos);

    gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_elev','png',900)
    return fig,ax

plt.close('all')
fig,ax=Plot_Elev(meta,roi)

#roi['gdf']['ogsr']['gdf'].plot(ax=ax[0],facecolor='None',edgecolor=[1,0,0],linewidth=0.5,label='Opening',alpha=1)
#roi['gdf']['op']['gdf'].plot(ax=ax[0],facecolor='None',edgecolor=[0,0,0],linewidth=0.5,label='Opening',alpha=1)

# gp=gpd.read_file(r'C:\Users\rhember\Documents\Data\GroundPlots\PSP-NADB\ground_plots.geojson')
# ind=np.where( (sl['Age_t0']>125) & (sl['Csw_L_t0']<50) )[0]
# gp.iloc[ind].plot(ax=ax[0],facecolor='None',marker='s',edgecolor=[1,0,0],linewidth=1,markersize=12,label='Opening',alpha=1)

# ind=np.where( (sl['Age_t0']>125) & (sl['Csw_L_t0']>250) )[0]
# gp.iloc[ind].plot(ax=ax[0],facecolor='None',marker='s',edgecolor=[0,1,0],linewidth=1,markersize=12,label='Opening',alpha=1)

# a=gpd.read_file(r'C:\Users\rhember\Documents\Data\GroundPlots\DellaSala et al 2022 IWB\data\v10\outputs.gdb')
# a=a.to_crs(roi['crs'])
# a.plot(ax=ax[0],facecolor='None',marker='^',edgecolor=[1,1,0],linewidth=1.25,markersize=14,label='Opening',alpha=1)


#ind=np.where(roi['gdf']['op']['gdf']['OPENING_ID']==1760606)[0]
#roi['gdf']['op']['gdf'].iloc[ind].plot(ax=ax[0],facecolor='None',edgecolor=[1,0,1],linewidth=2.5,label='Opening',alpha=1)

#%% Plot BTM

def Plot_ROI_BTM(meta,roi):

    # Grid
    bin=np.unique(roi['grd']['btm']['Compressed']['Data'])

    N_bin=bin.size
    N_hidden=1
    N_vis=N_bin+N_hidden

    z1=np.ones(roi['grd']['btm']['Compressed']['Data'].shape,dtype=int)
    for i in range(bin.size):
        z1[(roi['grd']['btm']['Compressed']['Data']==bin[i])]=i

    z1[(roi['grd']['Data']!=1)]=N_bin

    # Labels
    lab=[]
    for i in range(N_bin):
        try:
            lab.append(roi['grd']['btm']['Compressed']['lab'][i])
        except:
            lab.append('')

    # Colormap
    cm=plt.cm.get_cmap('viridis',N_bin);
    for i in range(N_bin):
        cm.colors[i,0:3]=roi['grd']['btm']['Compressed']['cl1'][i,:]
    cm=np.vstack( (cm.colors,(1,1,1,1)) )
    cm=matplotlib.colors.ListedColormap(cm)

    plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
    im=ax[0].matshow(z1,clim=(0,N_vis),extent=roi['grd']['Extent'],cmap=cm)
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

    # Add relief shading
    # z=roi['grd']['elev']['Data']; dx,dy=roi['grd']['elev']['Cellsize'],roi['grd']['elev']['Cellsize']
    # ls=LightSource(azdeg=90,altdeg=45)
    # ve=0.1
    # hs=ls.hillshade(z,vert_exag=ve,dx=dx,dy=dy)
    # ax[0].matshow(hs,extent=roi['grd']['elev']['Extent'],cmap='Greys',alpha=0.4,clim=(np.min(hs),np.max(hs)))

    #cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_vis,1),ticks=np.arange(0.5,N_vis+1.5,1))
    #cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_vis-1,1),ticks=np.arange(0.5,N_vis+1.5,1))
    cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_vis-(N_hidden-1),1),ticks=np.arange(0.5,N_vis-N_hidden,1))
    cb.ax.set(yticklabels=lab)
    cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
    cb.outline.set_edgecolor('w')
    for i in range(0,N_vis):
        ax[1].plot([0,100],[i,i],'w-',linewidth=1.5)
    pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
    ax[1].set(position=pos2)

    return fig,ax

plt.close('all')
fig,ax=Plot_ROI_BTM(meta,roi)
#roi['gdf']['tpf'].plot(ax=ax[0],marker='^',edgecolor='c',facecolor=[0.5,1,1],markersize=75,zorder=2)
#roi['gdf']['ogsr']['gdf'].plot(ax=ax[0],facecolor='None',edgecolor=[1,1,0],linewidth=0.5,label='Opening',alpha=1)
#op['gdf'].plot(ax=ax[0],facecolor='None',ls='-',edgecolor=[0.8,0.8,1],linewidth=2,label='Wildfire',alpha=1)
#fcinv['gdf'].plot(ax=ax[0],facecolor='None',ls='-',edgecolor=[0.4,0.4,1],linewidth=1,label='Wildfire',alpha=1)
#fcres['gdf'].plot(ax=ax[0],facecolor='None',ls='--',edgecolor=[1,0.5,0.25],linewidth=1,label='Wildfire',alpha=1)
#atu['gdf'].plot(ax=ax[0],facecolor='None',ls='--',edgecolor=[1,0,0],linewidth=2,label='Planting',alpha=1)
#wf['gdf'].plot(ax=ax[0],facecolor='None',ls='-',edgecolor=[1,0.5,0.5],linewidth=2,label='Planting',alpha=1)
#wfp.plot(ax=ax[0],facecolor='None',edgecolor=[0,0,0],linewidth=1,label='Wildfire',alpha=1)
#atup['gdf overlay'].plot(ax=ax[0],facecolor=[0,0,0],edgecolor=[0,0.5,0],linewidth=1.25,label='Planting',alpha=0.25)
gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_btm','png',900)

#%% PLOT Land Cover 2020 from CEC

def Plot_CEC_LC20(meta,roi):

    dCEC=gu.ReadExcel(r'C:\Users\rhember\Documents\Data\BC1ha\LUTs\LUT_lcc_cec_Compressed.xlsx')

    bin=np.arange(1,dCEC['Name'].size+1,1)
    N_vis=bin.size+1
    N_hidden=1
    N_tot=N_vis+N_hidden

    z0=roi['grd']['lcc_cec_2020']['Data']
    z1=bin.size*np.ones(z0.shape,dtype=int)
    for i in range(bin.size):
        z1[(z0==bin[i])]=i
    z1[(roi['grd']['Data']!=1)]=i+1

    for i in range(bin.size):
        z1[0,i]=i

    # Labels
    lab=list(meta['LUT']['Derived']['lcc_cec_c'].keys())
    #lab.append('')
    #lab.append('')

    # Colormap
    cm=plt.cm.get_cmap('viridis',bin.size);
    for i in range(bin.size):
        cm.colors[i,0:3]=[ dCEC['R'][i]/255,dCEC['G'][i]/255,dCEC['B'][i]/255 ]
    cm=np.vstack( (cm.colors,(1,1,1,1)) )
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

    cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_vis-(N_hidden-1),1),ticks=np.arange(0.5,N_vis-N_hidden,1))
    cb.ax.set(yticklabels=lab)
    cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
    cb.outline.set_edgecolor('w')
    for i in range(0,N_vis):
        ax[1].plot([0,100],[i,i],'w-',linewidth=1.5)
    pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
    ax[1].set(position=pos2)
    gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_lc20_cec','png',900)

    return fig,ax

fig,ax=Plot_CEC_LC20(meta,roi)

#%% Plot LUC (Afforestation or Deforestation) from CEC

def Plot_LUC(meta,roi):

    type='def'
    if type=='def':
        lab=['Forest to Settlement','Forest to Cropland','Forest to Barren Ground','Land']# Deforestation
    elif type=='aff':
        lab=['Settlement to Forest','Cropland to Forest','Barren Ground to Forest','Land'] # Afforestation
    else:
        pass

    z0=roi['grd']['luc_' + type + '_cec']['Data']
    z1=z0
    z1[(z0==0)]=4
    z1[(roi['grd']['Data']!=1) | (roi['grd']['lcc1_c']['Data']==meta['LUT']['Derived']['lcc1']['Water'])]=5

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
    gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_luc_' + type + '_cec','png',900)

    return fig,ax

fig,ax=Plot_LUC(meta,roi)

#%% Density Class

def Plot_DensityClass(meta,roi,gdf):

    lab=['Treed sparse','Treed open','Treed dense','Shrubland','Grassland','Non-treed']
    z1=roi['grd']['denseclass']['Data']
    z1[roi['grd']['Data']==0]=np.max(z1)

    N_vis=6
    N_hidden=1
    N_tot=N_vis+N_hidden

    # Colormap
    cm=np.vstack( ((0.75,0.95,0.6,1),(0,0.7,0,1),(0,0.25,0,1),(0.8,0.6,0.6,1),(1,1,0.5,1),(0.85,0.85,0.85,1),(1,1,1,1)) )
    cm=matplotlib.colors.ListedColormap(cm)

    plt.close('all')
    fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
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

    gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_denseclass','png',900)

    return fig,ax

Plot_DensityClass(meta,roi,gdf)

#%%

Plot_REARs(meta,roi,gdf)

#%% Plot distance from road

def Plot_DistanceFromRoad(meta,roi):

    z=roi['grd']['d2fac']['Data']
    if roi['Type']=='ByTSA':
        z[roi['grd']['Data']==0]=0

    # Plot
    plt.close('all')
    fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
    im=ax[0].matshow(z,extent=roi['grd']['Extent'],cmap='Greys')

    #roi['gdf']['tsa'].plot(ax=ax[0],color=None,edgecolor=[0,0,0],facecolor='none')
    ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'],xticklabels='',yticklabels='')
    ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

    # Add relief shading
    z=roi['grd']['elev']['Data']; dx,dy=roi['grd']['Cellsize'],roi['grd']['Cellsize']
    ls=LightSource(azdeg=90,altdeg=45)
    ve=0.1
    hs=ls.hillshade(z,vert_exag=ve,dx=dx,dy=dy)
    ax[0].matshow(hs,extent=roi['grd']['elev']['Extent'],cmap='Greys',alpha=0.4,clim=(np.min(hs),np.max(hs)))

    cb=plt.colorbar(im,cax=ax[1])#,boundaries=np.arange(0,N_vis-(N_hidden-1),1),ticks=np.arange(0.5,N_vis+1.5,1))
    pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
    ax[1].set(position=pos2)

    return fig,ax

plt.close('all')
fig,ax=Plot_DistanceFromRoad(meta,roi)

#roi['gdf']['road'].plot(ax=ax[0],edgecolor='y',linewidth=2,label='Road',alpha=1,zorder=1)
roi['gdf']['road'].plot(ax=ax[0],edgecolor=[0.7,0.9,1],linewidth=1,label='Road',alpha=1,zorder=1)

# Plot lumber mills, pulp mills and chipper mills
mtypeL=['LBR','PLP','PLT']
ax=u1ha.PlotMills(ax,gdf['tpf']['gdf'],mtypeL,labels='On')

gdf['cities'].plot(ax=ax[0],marker='s',edgecolor=[0.75,0.3,0],facecolor=[1,0.6,0],lw=1,markersize=20,alpha=1,zorder=2)
for x,y,label in zip(gdf['cities'].geometry.x,gdf['cities'].geometry.y,gdf['cities'].Name):
    ax[0].annotate(label,xy=(x,y),xytext=(4,3),textcoords="offset points",color=[0.75,0.3,0])

#%% Plot SI

def Plot_SI(meta,roi):

    z=roi['grd']['SITE_INDEX']['Data']
    if roi['Type']=='ByTSA':
        z[roi['grd']['Data']==0]=0

    # Plot
    plt.close('all')
    fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
    im=ax[0].matshow(roi['grd']['SITE_INDEX']['Data'],extent=roi['grd']['Extent'],cmap='magma',clim=[5,22])

    #roi['gdf']['tsa'].plot(ax=ax[0],color=None,edgecolor=[0,0,0],facecolor='none')
    ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'],xticklabels='',yticklabels='')
    ax[0].grid(False)

    # # Add relief shading
    # z=roi['grd']['elev']['Data']; dx,dy=roi['grd']['Cellsize'],roi['grd']['Cellsize']
    # ls=LightSource(azdeg=90,altdeg=45)
    # ve=0.1
    # hs=ls.hillshade(z,vert_exag=ve,dx=dx,dy=dy)
    # ax[0].matshow(hs,extent=roi['grd']['elev']['Extent'],cmap='Greys',alpha=0.4,clim=(np.min(hs),np.max(hs)))

    cb=plt.colorbar(im,cax=ax[1])#,boundaries=np.arange(0,N_vis-(N_hidden-1),1),ticks=np.arange(0.5,N_vis+1.5,1))
    pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
    ax[1].set(position=pos2)

    return fig,ax

plt.close('all')
fig,ax=Plot_SI(meta,roi)
#roi['gdf']['road'].plot(ax=ax[0],edgecolor='y',linewidth=1.25,label='Road',alpha=1)
#roi['gdf']['road'].plot(ax=ax[0],edgecolor='k',linewidth=0.5,label='Road',alpha=1)
#roi['gdf']['ogsr']['gdf'].plot(ax=ax[0],facecolor='None',edgecolor=[1,1,0],linewidth=1,label='Opening',alpha=1)
#roi['gdf']['op']['gdf'].plot(ax=ax[0],facecolor='None',edgecolor=[0,1,1],linewidth=1,label='Opening',alpha=1)

#%% PLOT PFI stemwood carbon

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

plt.close('all')
fig,ax=Plot_PFI(meta,roi)

#gdf['cities'].plot(ax=ax[0],marker='s',edgecolor=[0.75,0.3,0],facecolor=[1,0.6,0],lw=1,markersize=20,alpha=1,zorder=2)
#for x,y,label in zip(gdf['cities'].geometry.x,gdf['cities'].geometry.y,gdf['cities'].Name):
#    ax[0].annotate(label,xy=(x,y),xytext=(4,3),textcoords="offset points",color=[0.75,0.3,0])
#roi['gdf']['cc']['gdf'].plot(ax=ax[0],facecolor=None,edgecolor=[0.5,1,1],linewidth=0.25,label='Opening',alpha=0.3)

gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_pfi_c','png',900)

#%% PLOT PFI stemwood carbon (20 m resolution(

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
    #    plt.text(row['geometry'].centroid.x,row['geometry'].centroid.y,str(row['HARVEST_YEAR']),ha='center',fontsize=5,color=[0.25,0.75,1])
    return fig,ax

#%% PLOT GLOB Biomass

def Plot_biomass_glob(meta,roi):

    # Grid
    bw=20; bin=np.arange(0,200+bw,bw);
    z1=(bin.size)*np.ones( roi['grd']['biomass_glob']['Data'].shape)
    for i in range(bin.size):
        ind=np.where(np.abs( roi['grd']['biomass_glob']['Data']-bin[i])<=bw/2)
        if ind[0].size>0:
            z1[ind]=i
        else:
            z1[0,i]=i
    ind=np.where(roi['grd']['biomass_glob']['Data']>bin[i]); z1[ind]=i
    z1[(roi['grd']['Data']==1) & ( roi['grd']['biomass_glob']['Data']==0)]=i+1
    z1[(roi['grd']['Data']!=1) | (roi['grd']['lcc1_c']['Data']==meta['LUT']['Derived']['lcc1']['Water'])]=i+2
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

    gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_biomass_glob','png',900)

    return fig,ax

fig,ax=Plot_biomass_glob(meta,roi)

#%% PLOT soc

def Plot_SOC_WithinROI(meta,roi):

    # Grid
    bw=10; bin=np.arange(0,100+bw,bw);
    z1=(bin.size)*np.ones( roi['grd']['gsoc']['Data'].shape)
    for i in range(bin.size):
        ind=np.where(np.abs( roi['grd']['gsoc']['Data']-bin[i])<=bw/2)
        z1[ind]=i
    ind=np.where(roi['grd']['gsoc']['Data']>=bin[i]); z1[ind]=i
    z1[(roi['grd']['Data']==1) & ( roi['grd']['gsoc']['Data']==0)]=i+1
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
    gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_gsoc','png',900)

    return fig,ax

fig,ax=Plot_SOC_WithinROI(meta,roi)

#%% Plot Crown Cover Percent from VRI

def Plot_CrownCover(meta,roi):

    # Grid
    bw=10; bin=np.arange(0,100+bw,bw);
    z1=(bin.size)*np.ones( roi['grd']['crownc']['Data'].shape)
    for i in range(bin.size):
        ind=np.where(np.abs( roi['grd']['crownc']['Data']-bin[i])<=bw/2)
        z1[ind]=i
    ind=np.where(roi['grd']['crownc']['Data']>bin[i])
    z1[ind]=0
    z1[(roi['grd']['Data']==1) & ( roi['grd']['crownc']['Data']==0)]=i+1
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

    gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_crownc','png',900)

    return fig,ax

fig,ax=Plot_CrownCover(meta,roi)

#%% Plot Live SPH from VRI

def Plot_SPH_Live(meta,roi):

    # Grid
    bw=400; bin=np.arange(0,2400,bw);
    z1=(bin.size)*np.ones( roi['grd']['sphlive']['Data'].shape)
    for i in range(bin.size):
        ind=np.where(np.abs( roi['grd']['sphlive']['Data']-bin[i])<=bw/2)
        z1[ind]=i
    ind=np.where(roi['grd']['sphlive']['Data']>bin[i])
    z1[ind]=i
    z1[(roi['grd']['Data']==1) & ( roi['grd']['sphlive']['Data']==0)]=i+1
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
    #     ax[1].plot([0,100],[i/(N_vis-N_hidden),i/(N_vis-N_hidden)],'k-',linewidth=0.5)
    # ax[1].set(position=meta['Graphics']['ax2 pos long']);
    # pos2=copy.copy(meta['Graphics']['pos2'])
    # pos2[1]=0.6
    # pos2[3]=0.24
    # ax[1].set(position=pos2)

    return fig,ax

plt.close('all')
fig,ax=Plot_SPH_Live(meta,roi)
#wfp.plot(ax=ax[0],facecolor='None',edgecolor=[0,0,0],linewidth=1,label='Wildfire',alpha=1)
#atup['gdf overlay'].plot(ax=ax[0],facecolor='None',edgecolor=[0,0,0],linewidth=1.25,label='Planting',alpha=1)
#gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\sphlive','png',300)

#%% Plot harvested year

def Plot_HarvestYear(meta,roi):

    nam='ntem'
    z0=roi['grd']['harv_yr_' + nam]['Data']

    bw=10; bin=np.arange(1855,2025+bw,bw);
    z1=(bin.size)*np.ones( z0.shape)
    for i in range(bin.size):
        ind=np.where(np.abs( z0-bin[i])<=bw/2)
        z1[ind]=i
    z1[(roi['grd']['Data']==0) | (roi['grd']['lcc1_c']['Data']==meta['LUT']['Derived']['lcc1']['Water'])]=i+2
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
    gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_harv_yr_' + nam,'png',900)
    return fig,ax

fig,ax=Plot_HarvestYear(meta,roi)

#%% Plot Land Use Change year

def Plot_LUC_Year(meta,roi):

    z0=roi['grd']['luc1_yr']['Data']

    bw=10; bin=np.arange(1850,2020+bw,bw);
    z1=(bin.size)*np.ones( z0.shape)
    for i in range(bin.size):
        ind=np.where(np.abs( z0-bin[i])<=bw/2)
        z1[ind]=i
    z1[(roi['grd']['Data']==0) | (roi['grd']['lcc1_c']['Data']==meta['LUT']['Derived']['lcc1']['Water'])]=i+2
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
    gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_luc1_yr','png',900)
    return fig,ax

fig,ax=Plot_LUC_Year(meta,roi)

#%% Plot wildfire year

def Plot_FireYear(meta,roi):

    # Grid
    bw=10; bin=np.arange(1910,2020+bw,bw);
    z1=(bin.size)*np.ones( roi['grd']['fire_yr']['Data'].shape)
    for i in range(bin.size):
        ind=np.where(np.abs( roi['grd']['fire_yr']['Data']-bin[i])<=bw/2)
        z1[ind]=i
    z1[(roi['grd']['Data']==0) | (roi['grd']['lcc1_c']['Data']==meta['LUT']['Derived']['lcc1']['Water'])]=i+2
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
    gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_fire_yr','png',900)

    return fig,ax

fig,ax=Plot_FireYear(meta,roi)

#%% Plot Global Forest Change Loss Year

def Plot_GFC_LossYear(meta,roi):

    z0=roi['grd']['gfcly']['Data']
    #z0=roi['grd']['gfcly_filt']['Data']

    bw=1; bin=np.arange(2001,2022,bw);
    z1=(bin.size)*np.ones( z0.shape)
    for i in range(bin.size):
        ind=np.where(z0==bin[i])
        z1[ind]=i
    z1[(z0==0)]=i+1
    #z1[(roi['grd']['Data']==0)]=i+2
    z1[(roi['grd']['Data']==0) | (roi['grd']['lcc1_c']['Data']==meta['LUT']['Derived']['lcc1']['Water'])]=i+2
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

    cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_vis-(N_hidden-1),1),ticks=np.arange(0.5,N_vis-N_hidden,1))
    cb.ax.set(yticklabels=lab)

    cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
    cb.outline.set_edgecolor('w')
    for i in range(0,N_vis):
        ax[1].plot([0,100],[i,i],'w-',linewidth=1.5)

    pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
    ax[1].set(position=pos2)

    gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_gfcly','png',900)
    #gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_gfcly_filt','png',900)

    return fig,ax

fig,ax=Plot_GFC_LossYear(meta,roi)

#%% Plot Burn Severity within the TSA mask

def Plot_ROI_BSR(bsr):

    z1=np.zeros(bsr['grd']['Data'].shape)
    ic=np.where(bsr['key']['Code']=='High')[0]; z1[np.where(bsr['grd']['Data']==bsr['key']['ID'][ic])]=1
    ic=np.where(bsr['key']['Code']=='Medium')[0]; z1[np.where(bsr['grd']['Data']==bsr['key']['ID'][ic])]=2
    ic=np.where(bsr['key']['Code']=='Low')[0]; z1[np.where(bsr['grd']['Data']==bsr['key']['ID'][ic])]=3
    ind=np.where(z1==0);
    z1[ind]=3
    ind=np.where(roi['Mask']['Data']==0)
    z1[ind]=4

    lab=['High','Medium','Low','Unburned','']

    # Number of colours and number of colours excluded from colorbar
    N_vis=5
    N_hidden=0

    # Colormap
    cm=np.vstack( ( (0.5,0,0,1),(1,0.25,0.25,1),(1,0.75,0.75,1),(0.96,0.96,0.96,1),(1,1,1,1) ) )
    cm=matplotlib.colors.ListedColormap(cm)

    fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
    im=ax[0].matshow(z1,clim=(0,N_vis),extent=lc2['grd']['Extent'],cmap=cm)
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
    ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
    ax[0].grid(False)

    cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_vis,1),ticks=np.arange(0.5,N_vis+1.5,1))
    cb.ax.set(yticklabels=lab)
    cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
    for i in range(0,N_vis):
        ax[1].plot([0,100],[i/(N_vis-N_hidden-1),i/(N_vis-N_hidden-1)],'k-',linewidth=0.5)
    pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
    ax[1].set(position=pos2)

    return fig,ax

plt.close('all')
fig,ax=Plot_ROI_BSR(bsr)
#wfp.plot(ax=ax[0],facecolor='None',edgecolor=[0,0,0],linewidth=1.5,label='Wildfire',alpha=1)
#atup['gdf overlay'].plot(ax=ax[0],facecolor=[0,0,0],edgecolor=[0,0.5,0],linewidth=1.25,label='Planting',alpha=0.25)
#gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\Planted areas and BSR','png',300)


#%% Plot soil water content

def Plot_SoilWaterContent(meta,roi):

    bw=20; bin=np.arange(0,200+bw,bw);
    z1=bin.size*np.ones( roi['grd']['ws_gs_n']['Data'].shape)
    for i in range(bin.size):
        ind=np.where(np.abs( roi['grd']['ws_gs_n']['Data']-bin[i])<=bw/2)
        if ind[0].size>0:
            z1[ind]=i
        else:
            z1[0,i]=i
    z1[(roi['grd']['ws_gs_n']['Data']>200)]=i
    z1[(roi['grd']['lcc1']==0) | (roi['grd']['lcc1']==meta['LUT']['Derived']['lcc1']['Water']) | (roi['grd']['Data']!=1)]=i+1

    lab=bin.astype(str)
    lab=np.append(lab,'Hidden')

    N_vis=bin.size
    N_hidden=1
    N_tot=N_vis+N_hidden

    cm=np.column_stack((roi['grd']['ws_gs_n']['cm']['cl1'],roi['grd']['ws_gs_n']['cm']['cl2'],roi['grd']['ws_gs_n']['cm']['cl3'],np.ones(roi['grd']['ws_gs_n']['cm']['bin'].size)))
    cm=np.vstack( (cm,(1,1,1,1)) ) # (0.83137,0.81569,0.78431,1)
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
    ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
    ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

    zmn=np.min(z1);
    zmx=np.max(z1);
    cb_ivl=(zmx-zmn)/N_tot;
    cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
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

    gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_SoilWaterContent_Normal','png',300)

    return fig,ax

fig,ax=Plot_SoilWaterContent(meta,roi)

#%% Plot ground plot sample grid

tmp=roi['grd']['Data'].copy()
ind=np.where(roi['grd']['lcc1_c']['Data']==meta['LUT']['lcc1']['Forest'])
tmp[ind]=2

cm=np.vstack( ((1,1,1,1),(0.95,0.95,0.95,1),(0.8,0.8,0.8,1)) )
cm=matplotlib.colors.ListedColormap(cm)

ms=2
plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
mp=ax.matshow(tmp,extent=roi['grd']['Extent'],cmap=cm,label='Forest mask')
roi['gdf']['bound'].plot(ax=ax,facecolor='None',edgecolor=[0,0,0],label='Political Boundary',linewidth=0.25,alpha=1)
roi['gdf']['bound within'].plot(ax=ax,facecolor='None',edgecolor=[0,0,0],label='Political Boundary',linewidth=0.25,alpha=1)

ind=np.where(tmp==2)
A=ind[0].size/1e3

ivl=32
ind=np.where(tmp[0::ivl,0::ivl]==2)
x=roi['grd']['X'][0::ivl,0::ivl][ind]
y=roi['grd']['Y'][0::ivl,0::ivl][ind]
rho1=x.size/A
print(A)
print(rho1)

points=[]
for k in range(x.size):
    points.append(Point(x[k],y[k]))
gdf1=gpd.GeoDataFrame({'geometry':points,'ID':np.ones(x.size)})
gdf1.crs=roi['crs']
gdf1=gpd.overlay(gdf1,roi['gdf']['bound within'],how='intersection')
gdf1.plot(ax=ax,markersize=0.7,facecolor=[0.75,0,0],edgecolor=None,linewidth=0.75,alpha=1,label='Proposed sample')

ax.set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'],aspect='auto')
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both')
ax.grid(meta['Graphics']['Map']['Map Grid Vis'])
ax.axis(meta['Graphics']['Map']['Map Axis Vis'])
#gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Biomass\SamplingPowerMap2','png',900)

#%% Salvage mask from timber cruise

fig,ax=u1ha.Plot_SalvageLogging(meta,roi)

#%% Plot consolidate forest range

def Plot_RangeTenure(meta,roi):

    lab=['Forest with grazing tenure','Forest with haycutting tenure','Forest with no range tenure','Non-forest land','Non land']
    lab=lab[0:-1]
    z1=roi['grd']['rangecon']['Data']
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

    gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_range_consol','png',900)

    return fig,ax

fig,ax=Plot_RangeTenure(meta,roi)

#%% Plot reserves

def Plot_ReservesConsolidated(meta,roi):

    lab=list(meta['LUT']['rescon'].keys())
    lab=lab[0:-1]
    z1=roi['grd']['rescon']['Data']
    z1[roi['grd']['Data']==0]=8

    # Number of colours and number of colours excluded from colorbar
    N_vis=6
    N_hidden=2
    N_tot=N_vis+N_hidden

    # Colormap
    cm=np.vstack( ((0.6,1,0,1),(1,0.5,0,1),(0.5,0.8,1,1),(0,0.5,0,1),(0.75,0,1,1), (0.4,0.4,0.4,1),(0.83,0.83,0.83,1),(1,1,1,1)) )
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

    gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_reserves_consol','png',900)

    return fig,ax

fig,ax=Plot_ReservesConsolidated(meta,roi)