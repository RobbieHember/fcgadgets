
#%% Import modules

import numpy as np
import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib.colors
from matplotlib.colors import LightSource
from shapely.geometry import Polygon,Point
import copy
import time
import fcgadgets.macgyver.utilities_general as gu
import fcgadgets.macgyver.utilities_gis as gis
import fcgadgets.bc1ha.bc1ha_utilities as u1ha
from fcgadgets.macgyver import utilities_query_gdb as qgdb

#%% Import parameters

meta=u1ha.Init()
meta=u1ha.ImportLUTs(meta)
meta['Graphics']={}
meta['Graphics']['figwidth']=20
meta['Graphics']['sidespace']=0.45
meta['Graphics']['ax1 pos']=[0.01,0.01,1-meta['Graphics']['sidespace']-0.01,0.98]
meta['Graphics']['ax1 vis']='off'
meta['Graphics']['ax1 gridvis']=False
#meta['Graphics']['leg x']=meta['Graphics']['ax1 pos'][0]+0.005 # Upper left
meta['Graphics']['leg x']=meta['Graphics']['ax1 pos'][0]+meta['Graphics']['ax1 pos'][2]+0.005 # Upper right side
meta['Graphics']['leg text space']=0.035
meta['Graphics']['leg w']=0.03
gp=gu.SetGraphics('Manuscript')

#%% Import base maps (full province)

gdf=u1ha.Import_GDBs_ProvinceWide(meta)

#%% Define region of interest

# Identify how ROI is defined
flg_roi='ByRegDis'
#flg_roi='ByTSA'
#flg_roi='ByLatLon'

roi={}
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

roi['Name']='CAPITAL'
#roi['Name']='STRATHCONA'
#roi['Name']='ROI_CowichanValley'

if flg_roi=='ByTSA':
    meta['Paths']['Figures']=r'C:\Users\rhember\OneDrive - Government of BC\Figures' + '\\TSAs\\' + roi['Name']
elif flg_roi=='ByRegDis':
    meta['Paths']['Figures']=r'C:\Users\rhember\OneDrive - Government of BC\Figures' + '\\Regional Districts\\' + roi['Name']

t0=time.time()
if flg_roi=='ByTSA':
    roi['Type']='ByTSA'
    # Search: gdf['tsa']['key']
    #roi['List']=['Arrowsmith TSA']
    #roi['List']=['Fort St. John TSA']
    #roi['List']=['Boundary TSA']
    #roi['List']=['Cassiar TSA']
    #roi['List']=['Fort Nelson TSA']
    #roi['List']=['Merritt TSA']
    #roi['List']=['Prince George TSA']
    #roi['List']=['Kamloops TSA','100 Mile House TSA','Williams Lake TSA']
    roi['List']=['Williams Lake TSA']
    #roi['List']=['North Island TSA']
    #roi['List']=['Quesnel TSA']
    #roi['List']=['100 Mile House TSA']
    #roi['List']=['Kootenay Lake TSA']
    #roi['List']=['100 Mile House TSA']
    #roi['List']=['Okanagan TSA']
    #roi['List']=['Merritt TSA','Kamloops TSA','100 Mile House TSA','Okanagan TSA','Williams Lake TSA','Lillooet TSA','Boundary TSA'] # ,'Arrow TSA','Revelstoke TSA'
    #roi['List']=list(gdf['tsa']['key']['Name'])

    # Western spruce budworm study (do not change!)
    #roi['List']=['Merritt TSA','Kamloops TSA','100 Mile House TSA','Okanagan TSA','Williams Lake TSA','Lillooet TSA','Boundary TSA'] # ,'Arrow TSA','Revelstoke TSA'

elif flg_roi=='ByRegDis':

    roi['Type']='ByRegDis'
    # Search: gdf['regdis']['gdf']['REGIONAL_DISTRICT_NAME'].unique()
    roi['List']=['CAPITAL']
    #roi['List']=['COWICHAN VALLEY']
    #roi['List']=['STRATHCONA']

elif flg_roi=='ByLatLon':

    roi['Type']='ByLatLon'

    # Hanceville fire
    flg=0
    if flg==1:
        roi['Centre']=[-122.92,51.92]
        roi['Radius']=40*1000
        meta['Paths']['Figures']=r'C:\Users\rhember\OneDrive - Government of BC\Figures\Reforestation\Reforestation Hanceville ROI'

    # Elephant Hill fire
    flg=0
    if flg==1:
        roi['Centre']=[-121.15,51.15]
        roi['Radius']=45*1000 # metres
        meta['Paths']['Figures']=r'C:\Users\rhember\OneDrive - Government of BC\Figures\Reforestation\Reforestation Elephant Hill ROI'

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

# Prepare region of interest
# Roads is empty because it takes a long time
roi=u1ha.DefineROI(meta,roi,gdf)
t1=time.time()
print((t1-t0)/60)

#%% Import rasters over ROI

vList=['lcc1','cut_yr','rears']#,'spc1_ntems','si_spl_ntems','SITE_INDEX','si_spl_fd','si_spl_hw'] #,'age_ntem','globbio','PROJ_AGE_1','cut_yr','fire_yr','gfcly','ws_gs_n','elev','becz'] 'tmean_ann_n'
#vList=['lcc1','elev','becz','plam','regentype','salvage','cut_yr_con1']
#vList=['lc2','lcc1','btm','elev','becz','denseclass','age1','cut_yr','globbio','gsoc','rears','ws_gs_n','gfcly','pfi_c','lc20_cec','plam']
#vList=['lc20_cec','gfcly','btm','cut_yr','plam']
#vList=['ws_gs_n','ta_ann_n']
#vList=['lc2','btm','elev','wsb_mask','wsb_treat','becz','age1','cut_yr','d2road','d2fac','idw_mask','protected','si']
#vList=['temp_norm','ws_gs_n','lc2','btm','elev','bgcz','cut_yr']
#vList=['lc2','btm','elev','soc','age1','si','temp_norm','ws_gs_n','cut_yr']
#vList=['lc2','btm','elev','bgcz','cut_yr','bsr','wf','age1','sphlive','sphdead']
#vList=['cut_yr']
#vList=['si']
roi=u1ha.Import_Raster_Over_ROI(meta,roi,vList)

#%% Import required vector geodatabases

#vList=['cc','fcres']
#vList=['wf'] # 'op','cc','fcres','ogsr'
#roi=u1ha.Import_GDB_Over_ROI(meta,roi,vList)

#%% Resource roads

#gdf['rd']={}
#gdf['rd']['gdf']=gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='FTEN_ROAD_SEGMENT_LINES_SVW')
#gdf['rd']['gdf']=gpd.overlay(gdf['rd']['gdf'],roi['gdf']['bound within'],how='intersection')

#%% Plot BGC Zones

def Plot_BGC_Zone(roi):

    z0=roi['grd']['bgc_zone']['Data']
    lab0=list(meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE'].keys())
    cl0=np.column_stack([meta['LUT']['Raw']['bgc_zone']['R'],meta['LUT']['Raw']['bgc_zone']['G'],meta['LUT']['Raw']['bgc_zone']['B']])
    id0=list(meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE'].values())
    z1,lab1,cl1=gis.CompressCats(z0,id0,lab0,cl0)

    N_vis=int(np.max(z1))
    N_hidden=2
    N_tot=N_vis+N_hidden
    ind=np.where(z1==0); z1[ind]=N_vis+1
    lab1=lab1[1:]
    cl1=cl1[1:,:]

    z1[0,0]=N_vis+1
    z1[(roi['grd']['Data']==0) | (roi['grd']['lcc1']['Data']==meta['LUT']['Derived']['lcc1']['Water'])]=N_vis+2

    lab1=np.append(lab1,[''])

    # Colormap
    cm=plt.cm.get_cmap('viridis',N_vis);
    for i in range(N_vis):
        cm.colors[i,0:3]=cl1[i,:]
    cm=np.vstack( (cm.colors,(0.9,0.9,0.9,1),(1,1,1,1)) )
    cm=matplotlib.colors.ListedColormap(cm)

    plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['figwidth'],(1-meta['Graphics']['sidespace'])*meta['Graphics']['figwidth']*roi['grd']['yxrat']))
    im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
    #roi['gdf']['bound within'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
    roi['gdf']['bc_bound'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
    roi['gdf']['lakes'].plot(ax=ax[0],facecolor=[0.82,0.88,1],label='Lakes',linewidth=2.25)
    roi['gdf']['rivers'].plot(ax=ax[0],color=[0.6,0.8,1],label='Rivers',linewidth=0.25)
    roi['gdf']['road'].plot(ax=ax[0],edgecolor='k',linewidth=0.25,label='Road',alpha=1,zorder=1)

    ax[0].set(position=meta['Graphics']['ax1 pos'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
    ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['ax1 gridvis']); ax[0].axis(meta['Graphics']['ax1 vis'])

    zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
    cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
    cb=plt.colorbar(im,cax=ax[1],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
    ax[1].set(position=[0.71,0.6,0.05,0.14])
    cb.ax.set(yticklabels=lab1)
    cb.ax.tick_params(labelsize=11,length=0)
    cb.outline.set_edgecolor('w')
    for i in range(cb_bnd.size):
        ax[1].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)

    pos2=[meta['Graphics']['leg x'],0.99-N_vis*meta['Graphics']['leg text space'],meta['Graphics']['leg w'],N_vis*meta['Graphics']['leg text space']]
    ax[1].set(position=pos2)

    gu.PrintFig(meta['Paths']['Figures'] + '\\' + roi['Name'] + '_bgc_zone','png',900)

    return fig,ax

fig,ax=Plot_BGC_Zone(roi)

#%% Species leading NTEMS

def Plot_Spc1(roi):

    # Compress categories
    z0=roi['grd']['spc1_ntems']['Data']
    lab0=np.array(list(meta['LUT']['Derived']['spc1_ntems'].keys()))
    id0=np.array(list(meta['LUT']['Derived']['spc1_ntems'].values()))
    cl0=np.column_stack([meta['LUT']['Raw']['spc1_ntems']['R'],meta['LUT']['Raw']['spc1_ntems']['G'],meta['LUT']['Raw']['spc1_ntems']['B']])
    z1,lab1,cl1=gis.CompressCats(z0,id0,lab0,cl0)
    z1=z1+1

    n=int(np.max(z1))
    z1[0,0]=n+1
    z1[(roi['grd']['Data']==0) | (roi['grd']['lcc1']['Data']==meta['LUT']['Derived']['lcc1']['Water'])]=n+2

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

    plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['figwidth'],(1-meta['Graphics']['sidespace'])*meta['Graphics']['figwidth']*roi['grd']['yxrat']))
    im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
    #roi['gdf']['bound within'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
    roi['gdf']['bc_bound'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
    roi['gdf']['lakes'].plot(ax=ax[0],facecolor=[0.82,0.88,1],label='Lakes',linewidth=2.25)
    roi['gdf']['rivers'].plot(ax=ax[0],color=[0.6,0.8,1],label='Rivers',linewidth=0.25)
    roi['gdf']['road'].plot(ax=ax[0],edgecolor='k',linewidth=0.25,label='Road',alpha=1,zorder=1)

    ax[0].set(position=meta['Graphics']['ax1 pos'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
    ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['ax1 gridvis']); ax[0].axis(meta['Graphics']['ax1 vis'])

    zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
    cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
    cb=plt.colorbar(im,cax=ax[1],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
    ax[1].set(position=[0.71,0.6,0.05,0.14])
    cb.ax.set(yticklabels=lab1)
    cb.ax.tick_params(labelsize=11,length=0)
    cb.outline.set_edgecolor('w')
    for i in range(cb_bnd.size):
        ax[1].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)

    pos2=[meta['Graphics']['leg x'],0.99-N_vis*meta['Graphics']['leg text space'],meta['Graphics']['leg w'],N_vis*meta['Graphics']['leg text space']]
    ax[1].set(position=pos2)

    gu.PrintFig(meta['Paths']['Figures'] + '\\' + roi['Name'] + '_spc1_ntems','png',900)

    return fig,ax

fig,ax=Plot_Spc1(roi)

#%% Land Cover Class 1

def Plot_LandCoverClass1(meta,roi,gdf):

    lab=list(meta['LUT']['Derived']['lcc1'].keys())
    z1=len(lab)*np.ones(roi['grd']['Data'].shape,dtype='int8')
    for k in meta['LUT']['Derived']['lcc1'].keys():
        ind=np.where(roi['grd']['lcc1']['Data']==meta['LUT']['Derived']['lcc1'][k])
        if ind[0].size>0:
            z1[ind]=meta['LUT']['Derived']['lcc1'][k]
        else:
            z1[0,meta['LUT']['Derived']['lcc1'][k]]=meta['LUT']['Derived']['lcc1'][k]
    z1[roi['grd']['Data']==0]=meta['LUT']['Derived']['lcc1'][k]+1

    N_vis=len(lab)
    N_hidden=1
    N_tot=N_vis+N_hidden

    # Colormap
    cm=np.vstack( ((0.3,0.6,0,1),(0.65,1,0,1),(1,1,0.5,1),(0.75,0.4,1,1),(0.75,0.5,0.5,1),(0.75,0.75,0.75,1),(0.8,0,0,1),(0.95,0.95,0.95,1),(0.75,0.85,0.95,1),(1,1,1,1)) )
    cm=matplotlib.colors.ListedColormap(cm)

    plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['figwidth'],(1-meta['Graphics']['sidespace'])*meta['Graphics']['figwidth']*roi['grd']['yxrat']))
    im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
    roi['gdf']['bound within'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
    roi['gdf']['bc_bound'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
    roi['gdf']['lakes'].plot(ax=ax[0],facecolor=[0.82,0.88,1],label='Lakes',linewidth=2.25)
    roi['gdf']['rivers'].plot(ax=ax[0],color=[0.6,0.8,1],label='Rivers',linewidth=0.25)
    roi['gdf']['road'].plot(ax=ax[0],edgecolor='k',linewidth=0.25,label='Road',alpha=1,zorder=1)

    #gdf['cities'].plot(ax=ax[0],marker='s',edgecolor=[0.75,0.3,0],facecolor=[1,0.6,0],lw=1,markersize=20,alpha=1,zorder=2)
    #for x,y,label in zip(gdf['cities'].geometry.x,gdf['cities'].geometry.y,gdf['cities'].Name):
    #    ax[0].annotate(label,xy=(x,y),xytext=(4,3),textcoords="offset points",color=[0.75,0.3,0],fontsize=4)

    ax[0].set(position=meta['Graphics']['ax1 pos'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
    ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['ax1 gridvis']); ax[0].axis(meta['Graphics']['ax1 vis'])

    zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
    cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
    cb=plt.colorbar(im,cax=ax[1],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
    ax[1].set(position=[0.71,0.6,0.05,0.14])
    cb.ax.set(yticklabels=lab)
    cb.ax.tick_params(labelsize=11,length=0)
    cb.outline.set_edgecolor('w')
    for i in range(cb_bnd.size):
        ax[1].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)

    pos2=[meta['Graphics']['leg x'],0.99-N_vis*meta['Graphics']['leg text space'],meta['Graphics']['leg w'],N_vis*meta['Graphics']['leg text space']]
    ax[1].set(position=pos2)

    gu.PrintFig(meta['Paths']['Figures'] + '\\' + roi['Name'] + '_lcc1','png',900)

    return fig,ax

fig,ax=Plot_LandCoverClass1(meta,roi,gdf)

#%% Plot mean annual temp

def Plot_MAT(roi):

    z0=roi['grd']['tmean_ann_n']['Data']

    bw=5; bin=np.arange(20,110+bw,bw)

    N_shown=bin.size
    N_hidden=2
    N_tot=N_shown+N_hidden

    z1=N_shown*np.ones(z0.shape)
    for i in range(N_shown):
        ind=np.where(np.abs(z0-bin[i])<=bw/2)
        if ind[0].size>0:
            z1[ind]=i+1
    ind=np.where(z0>=bin[i]); z1[ind]=i+1
    z1[1,1]=i+2
    z1[(roi['grd']['Data']==0) | (roi['grd']['lcc1']['Data']==0) | (roi['grd']['lcc1']['Data']==meta['LUT']['Derived']['lcc1']['Water'])]=i+3

    for i in range(N_shown):
        z1[0,i]=i+1

    lab=['']*(N_tot-1)
    lab[0:N_shown]=np.array(bin/10).astype(str)

    cm=plt.cm.get_cmap('viridis',N_shown)
    cm=np.vstack( (cm.colors,(0.9,0.9,0.9,1),(1,1,1,1)) )
    cm=matplotlib.colors.ListedColormap(cm)

    plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['figwidth'],(1-meta['Graphics']['sidespace'])*meta['Graphics']['figwidth']*roi['grd']['yxrat']))
    im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)

    #roi['gdf']['bound within'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
    roi['gdf']['bc_bound'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
    roi['gdf']['lakes'].plot(ax=ax[0],facecolor=[0.82,0.88,1],label='Lakes',linewidth=2.25)
    roi['gdf']['rivers'].plot(ax=ax[0],color=[0.6,0.8,1],label='Rivers',linewidth=0.25)
    roi['gdf']['road'].plot(ax=ax[0],edgecolor='k',linewidth=0.25,label='Road',alpha=1,zorder=1)

    ax[0].set(position=meta['Graphics']['ax1 pos'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
    ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['ax1 gridvis']); ax[0].axis(meta['Graphics']['ax1 vis'])

    zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
    cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
    cb=plt.colorbar(im,cax=ax[1],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
    ax[1].set(position=[0.71,0.6,0.05,0.14])
    cb.ax.set(yticklabels=lab)
    cb.ax.tick_params(labelsize=11,length=0)
    cb.outline.set_edgecolor('w')
    for i in range(cb_bnd.size):
        ax[1].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)

    ax[1].set(position=[0.56,0.98-0.8,0.03,0.8])

    gu.PrintFig(meta['Paths']['Figures'] + '\\' + roi['Name'] + '_tmean_ann_n','png',900)

    return fig,ax

fig,ax=Plot_MAT(roi)

#%% Plot mean annual precip

def Plot_MAP(roi):

    z0=roi['grd']['prcp_ann_n']['Data']

    bw=250; bin=np.arange(0,4000+bw,bw)

    N_shown=bin.size
    N_hidden=2
    N_tot=N_shown+N_hidden

    z1=N_shown*np.ones(z0.shape)
    for i in range(N_shown):
        ind=np.where(np.abs(z0-bin[i])<=bw/2)
        if ind[0].size>0:
            z1[ind]=i+1
    ind=np.where(z0>=bin[i]); z1[ind]=i+1
    z1[1,1]=i+2
    z1[(roi['grd']['Data']==0) | (roi['grd']['lcc1']['Data']==0) | (roi['grd']['lcc1']['Data']==meta['LUT']['Derived']['lcc1']['Water'])]=i+3

    for i in range(N_shown):
        z1[0,i]=i+1

    lab=['']*(N_tot-1)
    lab[0:N_shown]=bin.astype(str)

    cm=plt.cm.get_cmap('viridis',N_shown)
    cm=np.vstack( (cm.colors,(0.9,0.9,0.9,1),(1,1,1,1)) )
    cm=matplotlib.colors.ListedColormap(cm)

    plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['figwidth'],(1-meta['Graphics']['sidespace'])*meta['Graphics']['figwidth']*roi['grd']['yxrat']))
    im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)

    roi['gdf']['bound within'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
    roi['gdf']['bc_bound'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
    roi['gdf']['lakes'].plot(ax=ax[0],facecolor=[0.82,0.88,1],label='Lakes',linewidth=2.25)
    roi['gdf']['rivers'].plot(ax=ax[0],color=[0.6,0.8,1],label='Rivers',linewidth=0.25)
    roi['gdf']['road'].plot(ax=ax[0],edgecolor='k',linewidth=0.25,label='Road',alpha=1,zorder=1)

    ax[0].set(position=meta['Graphics']['ax1 pos'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
    ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['ax1 gridvis']); ax[0].axis(meta['Graphics']['ax1 vis'])

    zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
    cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
    cb=plt.colorbar(im,cax=ax[1],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
    ax[1].set(position=[0.71,0.6,0.05,0.14])
    cb.ax.set(yticklabels=lab)
    cb.ax.tick_params(labelsize=11,length=0)
    cb.outline.set_edgecolor('w')
    for i in range(cb_bnd.size):
        ax[1].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)

    ax[1].set(position=[0.56,0.05,0.03,0.9])

    gu.PrintFig(meta['Paths']['Figures'] + '\\' + roi['Name'] + '_prcp_ann_n','png',900)

    return fig,ax

fig,ax=Plot_MAP(roi)

#%% PLOT age from VRI

def Plot_PROJ_AGE_1(roi):

    z0=roi['grd']['PROJ_AGE_1']['Data']

    bw=20
    bin=np.arange(0,220,bw)
    id=np.arange(1,bin.size+1)

    N_shown=bin.size
    N_hidden=2
    N_tot=N_shown+N_hidden

    z1=N_shown*np.ones(z0.shape)
    for i in range(N_shown):
        ind=np.where(np.abs(z0-bin[i])<=bw/2)
        if ind[0].size>0:
            z1[ind]=id[i]
        else:
            z1[0,i]=id[i]
    ind=np.where(roi['grd']['PROJ_AGE_1']['Data']>bin[i]); z1[ind]=id[-1]
    z1[(roi['grd']['PROJ_AGE_1']['Data']==0) | (roi['grd']['lcc1']['Data']!=meta['LUT']['Derived']['lcc1']['Forest'])]=id[-1]+1
    z1[(roi['grd']['Data']==0) | (roi['grd']['lcc1']['Data']==0) | (roi['grd']['lcc1']['Data']==meta['LUT']['Derived']['lcc1']['Water'])]=id[-1]+2

    lab=['']*(N_tot-1)
    lab[0:N_shown]=bin.astype(str)

    cm=plt.cm.get_cmap('viridis',N_shown)
    cm=np.vstack( (cm.colors,(0.9,0.9,0.9,1),(1,1,1,1)) )
    cm=matplotlib.colors.ListedColormap(cm)

    plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['figwidth'],(1-meta['Graphics']['sidespace'])*meta['Graphics']['figwidth']*roi['grd']['yxrat']))
    im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)

    roi['gdf']['bound within'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
    roi['gdf']['bc_bound'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
    roi['gdf']['lakes'].plot(ax=ax[0],facecolor=[0.82,0.88,1],label='Lakes',linewidth=2.25)
    roi['gdf']['rivers'].plot(ax=ax[0],color=[0.6,0.8,1],label='Rivers',linewidth=0.25)
    roi['gdf']['road'].plot(ax=ax[0],edgecolor='k',linewidth=0.25,label='Road',alpha=1,zorder=1)

    ax[0].set(position=meta['Graphics']['ax1 pos'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
    ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['ax1 gridvis']); ax[0].axis(meta['Graphics']['ax1 vis'])

    zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
    cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
    cb=plt.colorbar(im,cax=ax[1],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
    ax[1].set(position=[0.71,0.6,0.05,0.14])
    cb.ax.set(yticklabels=lab)
    cb.ax.tick_params(labelsize=11,length=0)
    cb.outline.set_edgecolor('w')
    for i in range(cb_bnd.size):
        ax[1].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)

    pos2=meta['Graphics']['ax2 pos']
    pos2[1]=0.5
    pos2[2]=0.03
    pos2[3]=0.4
    ax[1].set(position=pos2)

    gu.PrintFig(meta['Paths']['Figures'] + '\\' + roi['Name'] + '_PROJ_AGE_1','png',900)

    return fig,ax

fig,ax=Plot_PROJ_AGE_1(roi)

#%% Plot Site Index from VRI

def Plot_SI(roi):

    #z0=roi['grd']['SITE_INDEX']['Data']
    #z0=roi['grd']['si_spl_fd']['Data']
    z0=roi['grd']['si_spl_hw']['Data']
    #z0=roi['grd']['si_spl_ntems']['Data']

    bw=2
    bin=np.arange(10,40+bw,bw)
    id=np.arange(1,bin.size+1)

    N_shown=bin.size
    N_hidden=2
    N_tot=N_shown+N_hidden

    z1=N_shown*np.ones(z0.shape)
    for i in range(N_shown):
        ind=np.where(np.abs(z0-bin[i])<=bw/2)
        if ind[0].size>0:
            z1[ind]=id[i]
        else:
            z1[0,i]=id[i]
    ind=np.where(z0>bin[i]); z1[ind]=id[-1]
    z1[(z0==0) | (roi['grd']['lcc1']['Data']!=meta['LUT']['Derived']['lcc1']['Forest'])]=id[-1]+1
    #z1[1,1]=id[-1]+1
    z1[(roi['grd']['Data']==0) | (roi['grd']['lcc1']['Data']==0) | (roi['grd']['lcc1']['Data']==meta['LUT']['Derived']['lcc1']['Water'])]=id[-1]+2

    lab=['']*(N_shown+1)
    lab[0:N_shown]=bin.astype(str)

    cm=plt.cm.get_cmap('viridis',N_shown)
    cm=np.vstack( (cm.colors,(0.9,0.9,0.9,1),(1,1,1,1)) )
    cm=matplotlib.colors.ListedColormap(cm)

    plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['figwidth'],(1-meta['Graphics']['sidespace'])*meta['Graphics']['figwidth']*roi['grd']['yxrat']))
    im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)

    #roi['gdf']['bound within'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
    roi['gdf']['bc_bound'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
    roi['gdf']['lakes'].plot(ax=ax[0],facecolor=[0.82,0.88,1],label='Lakes',linewidth=2.25)
    roi['gdf']['rivers'].plot(ax=ax[0],color=[0.6,0.8,1],label='Rivers',linewidth=0.25)
    roi['gdf']['road'].plot(ax=ax[0],edgecolor='k',linewidth=0.25,label='Road',alpha=1,zorder=1)

    ax[0].set(position=meta['Graphics']['ax1 pos'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
    ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['ax1 gridvis']); ax[0].axis(meta['Graphics']['ax1 vis'])

    zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
    cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
    cb=plt.colorbar(im,cax=ax[1],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
    ax[1].set(position=[0.71,0.6,0.05,0.14])
    cb.ax.set(yticklabels=lab)
    cb.ax.tick_params(labelsize=11,length=0)
    cb.outline.set_edgecolor('w')
    for i in range(cb_bnd.size):
        ax[1].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)

    pos2=meta['Graphics']['ax2 pos']
    pos2[1]=0.1
    pos2[2]=0.03
    pos2[3]=0.8
    ax[1].set(position=pos2)

    #gu.PrintFig(meta['Paths']['Figures'] + '\\' + roi['Name'] + '_SITE_INDEX','png',900)
    gu.PrintFig(meta['Paths']['Figures'] + '\\' + roi['Name'] + '_si_spl_hw','png',900)
    #gu.PrintFig(meta['Paths']['Figures'] + '\\' + roi['Name'] + '_si_spl_ntems','png',900)

    return fig,ax

fig,ax=Plot_SI(roi)



#%% PLOT age from NTEMS

def Plot_Age_NTEMS(roi):

    z0=roi['grd']['age_ntem']['Data']

    bw=20
    bin=np.arange(0,220,bw)
    id=np.arange(1,bin.size+1)

    N_shown=bin.size
    N_hidden=2
    N_tot=N_shown+N_hidden

    z1=N_shown*np.ones(z0.shape)
    for i in range(N_shown):
        ind=np.where(np.abs(z0-bin[i])<=bw/2)
        if ind[0].size>0:
            z1[ind]=id[i]
        else:
            z1[0,i]=id[i]
    ind=np.where(z0>bin[i]); z1[ind]=id[-1]
    z1[(z0==0) | (roi['grd']['lcc1']['Data']!=meta['LUT']['Derived']['lcc1']['Forest'])]=id[-1]+1
    z1[(roi['grd']['Data']==0) | (roi['grd']['lcc1']['Data']==0) | (roi['grd']['lcc1']['Data']==meta['LUT']['Derived']['lcc1']['Water'])]=id[-1]+2

    lab=['']*(N_tot-1)
    lab[0:N_shown]=bin.astype(str)

    cm=plt.cm.get_cmap('viridis',N_shown)
    cm=np.vstack( (cm.colors,(0.9,0.9,0.9,1),(1,1,1,1)) )
    cm=matplotlib.colors.ListedColormap(cm)

    plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['figwidth'],(1-meta['Graphics']['sidespace'])*meta['Graphics']['figwidth']*roi['grd']['yxrat']))
    im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)

    roi['gdf']['bound within'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
    roi['gdf']['bc_bound'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
    roi['gdf']['lakes'].plot(ax=ax[0],facecolor=[0.82,0.88,1],label='Lakes',linewidth=2.25)
    roi['gdf']['rivers'].plot(ax=ax[0],color=[0.6,0.8,1],label='Rivers',linewidth=0.25)
    roi['gdf']['road'].plot(ax=ax[0],edgecolor='k',linewidth=0.25,label='Road',alpha=1,zorder=1)

    ax[0].set(position=meta['Graphics']['ax1 pos'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
    ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['ax1 gridvis']); ax[0].axis(meta['Graphics']['ax1 vis'])

    zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
    cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
    cb=plt.colorbar(im,cax=ax[1],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
    ax[1].set(position=[0.71,0.6,0.05,0.14])
    cb.ax.set(yticklabels=lab)
    cb.ax.tick_params(labelsize=11,length=0)
    cb.outline.set_edgecolor('w')
    for i in range(cb_bnd.size):
        ax[1].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)

    pos2=meta['Graphics']['ax2 pos']
    pos2[1]=0.5
    pos2[2]=0.03
    pos2[3]=0.4
    ax[1].set(position=pos2)

    gu.PrintFig(meta['Paths']['Figures'] + '\\' + roi['Name'] + '_Age_NTEM','png',900)

    return fig,ax

fig,ax=Plot_Age_NTEMS(roi)

#%%

def Plot_RegenType(meta,roi,gdf):

    lab=['Harvested and planted','Harvested (pre-1988), not planted','Harvested (1988-2018), not planted','Harvested (post-2018), not planted','Straight planting','No harvesting / no planting']
    z1=7*np.ones(roi['grd']['Data'].shape,dtype='int8')
    for i in range(1,7):
        z1[(roi['grd']['regentype']['Data']==i)]=i
    ind=np.where( (roi['grd']['Data']==0) ); z1[ind]=7

    N_vis=len(lab)
    N_hidden=1
    N_tot=N_vis+N_hidden
    cm=np.vstack( ((0.85,0.7,1,1),(1,1,0,1),(1,0.25,0,1),(0.7,0.7,0.7,1),(0.7,1,0.2,1),(0.9,0.9,0.9,1),(1,1,1,1)) )
    cm=matplotlib.colors.ListedColormap(cm)

    plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['figwidth'],(1-meta['Graphics']['sidespace'])*meta['Graphics']['figwidth']*roi['grd']['yxrat']))
    im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
    roi['gdf']['bound within'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
    roi['gdf']['lakes'].plot(ax=ax[0],facecolor=[0.82,0.88,1],label='Lakes',linewidth=2.25)
    roi['gdf']['rivers'].plot(ax=ax[0],color=[0.6,0.8,1],label='Rivers',linewidth=0.25)
    roi['gdf']['road'].plot(ax=ax[0],edgecolor='k',linewidth=0.25,label='Road',alpha=1,zorder=1)
    #gdf['cities'].plot(ax=ax[0],marker='s',edgecolor=[0.75,0.3,0],facecolor=[1,0.6,0],lw=1,markersize=20,alpha=1,zorder=2)
    #for x,y,label in zip(gdf['cities'].geometry.x,gdf['cities'].geometry.y,gdf['cities'].Name):
    #    ax[0].annotate(label,xy=(x,y),xytext=(4,3),textcoords="offset points",color=[0.75,0.3,0],fontsize=4)
    ax[0].set(position=meta['Graphics']['ax1 pos'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
    ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['ax1 gridvis']); ax[0].axis(meta['Graphics']['ax1 vis'])

    zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
    cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
    cb=plt.colorbar(im,cax=ax[1],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
    ax[1].set(position=[0.71,0.6,0.05,0.14])
    cb.ax.set(yticklabels=lab)
    cb.ax.tick_params(labelsize=11,length=0)
    cb.outline.set_edgecolor('w')
    for i in range(cb_bnd.size):
        ax[1].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)
    pos2=meta['Graphics']['ax2 pos']
    pos2[1]=0.6
    pos2[3]=0.26
    ax[1].set(position=pos2)
    gu.PrintFig(meta['Paths']['Figures'] + '\\' + roi['Name'] + '_regentype','png',900)

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

    plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['figwidth'],(1-meta['Graphics']['sidespace'])*meta['Graphics']['figwidth']*roi['grd']['yxrat']))
    im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
    roi['gdf']['rivers'].plot(ax=ax[0],color=[0.6,0.8,1],label='Rivers',linewidth=0.25)
    roi['gdf']['lakes'].plot(ax=ax[0],color=[0.82,0.88,1],label='Lakes',linewidth=2.25)
    roi['gdf']['bound within'].plot(ax=ax[0],color=None,edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
    roi['gdf']['road'].plot(ax=ax[0],edgecolor='k',linewidth=0.5,label='Road',alpha=1,zorder=1)
    ax[0].set(position=meta['Graphics']['ax1 pos'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'],xticklabels='',yticklabels='')
    ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['ax1 gridvis']); ax[0].axis(meta['Graphics']['ax1 vis'])

    zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
    cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
    cb=plt.colorbar(im,cax=ax[1],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
    ax[1].set(position=[0.71,0.6,0.05,0.14])
    cb.ax.set(yticklabels=lab)
    cb.ax.tick_params(labelsize=11,length=0)
    cb.outline.set_edgecolor('w')
    for i in range(cb_bnd.size):
        ax[1].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)

    pos2=meta['Graphics']['ax2 pos']
    pos2[1]=0.7
    pos2[2]=0.027
    pos2[3]=0.09
    ax[1].set(position=pos2)

    gu.PrintFig(meta['Paths']['Figures'] + '\\' + roi['Name'] + '_planted','png',900)

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

    plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['figwidth'],(1-meta['Graphics']['sidespace'])*meta['Graphics']['figwidth']*roi['grd']['yxrat']))
    im=ax[0].matshow(roi['grd']['elev']['Data'],extent=roi['grd']['Extent'],cmap='Greys')
    roi['gdf']['bc_bound'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
    #roi['gdf']['bound within'].plot(ax=ax[0],color=None,edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
    roi['gdf']['rivers'].plot(ax=ax[0],color=[0.6,0.8,1],label='Rivers',linewidth=0.25)
    roi['gdf']['lakes'].plot(ax=ax[0],color=[0.82,0.88,1],label='Lakes',linewidth=2.25)
    #roi['gdf']['road'].plot(ax=ax[0],edgecolor='y',linewidth=2,label='Road',alpha=1,zorder=1)
    #roi['gdf']['road'].plot(ax=ax[0],edgecolor=[0,0,0],linewidth=0.5,label='Road',alpha=1,zorder=1)
    #roi['gdf']['wf']['gdf'].plot(ax=ax[0],facecolor=[1,.5,0],edgecolor=[1,.5,0],linewidth=0.5,label='Opening',alpha=0.2)

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

    ax[0].set(position=meta['Graphics']['ax1 pos'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'],xticklabels='',yticklabels='')
    ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['ax1 gridvis']); ax[0].axis(meta['Graphics']['ax1 vis'])

    # Add relief shading
    #z=roi['grd']['elev']['Data']; dx,dy=roi['grd']['Cellsize'],roi['grd']['Cellsize']
    #ls=LightSource(azdeg=90,altdeg=45)
    #ve=0.1
    #hs=ls.hillshade(z,vert_exag=ve,dx=dx,dy=dy)
    #ax[0].matshow(hs,extent=roi['grd']['elev']['Extent'],cmap='Greys',alpha=0.4,clim=(np.min(hs),np.max(hs)))

    plt.colorbar(im,cax=ax[1])#,boundaries=np.arange(0,N_vis-(N_hidden-1),1),ticks=np.arange(0.5,N_vis+1.5,1))
    pos=[meta['Graphics']['leg x'],0.5,meta['Graphics']['leg w'],0.4]
    ax[1].set(position=pos);

    gu.PrintFig(meta['Paths']['Figures'] + '\\' + roi['Name'] + '_elev','png',900)
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

def Plot_ROI_BTM(roi):

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

    plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['figwidth'],(1-meta['Graphics']['sidespace'])*meta['Graphics']['figwidth']*roi['grd']['yxrat']))
    im=ax[0].matshow(z1,clim=(0,N_vis),extent=roi['grd']['Extent'],cmap=cm)
    roi['gdf']['bc_bound'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
    roi['gdf']['rivers'].plot(ax=ax[0],color=[0.6,0.8,1],label='Rivers',linewidth=0.25)
    #roi['gdf']['lakes'].plot(ax=ax[0],color=[0.82,0.88,1],label='Lakes',linewidth=2.25)
    roi['gdf']['bound within'].plot(ax=ax[0],color=None,edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
    #roi['gdf']['road'].plot(ax=ax[0],edgecolor='y',linewidth=2.25,label='Road',alpha=1,zorder=1)
    #roi['gdf']['road'].plot(ax=ax[0],edgecolor='k',linewidth=0.2,label='Road',alpha=1,zorder=1)
    ax[0].set(position=meta['Graphics']['ax1 pos'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'],xticklabels='',yticklabels='')
    ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['ax1 gridvis']); ax[0].axis(meta['Graphics']['ax1 vis'])

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
    cb.ax.tick_params(labelsize=6,length=0)
    cb.outline.set_edgecolor('w')
    for i in range(0,N_vis):
        ax[1].plot([0,100],[i,i],'w-',linewidth=1.5)
    pos2=meta['Graphics']['ax2 pos']
    pos2[1]=0.35
    pos2[3]=0.55
    ax[1].set(position=pos2)

    return fig,ax

plt.close('all')
fig,ax=Plot_ROI_BTM(roi)
#roi['gdf']['tpf'].plot(ax=ax[0],marker='^',edgecolor='c',facecolor=[0.5,1,1],markersize=75,zorder=2)
#roi['gdf']['ogsr']['gdf'].plot(ax=ax[0],facecolor='None',edgecolor=[1,1,0],linewidth=0.5,label='Opening',alpha=1)
#op['gdf'].plot(ax=ax[0],facecolor='None',ls='-',edgecolor=[0.8,0.8,1],linewidth=2,label='Wildfire',alpha=1)
#fcinv['gdf'].plot(ax=ax[0],facecolor='None',ls='-',edgecolor=[0.4,0.4,1],linewidth=1,label='Wildfire',alpha=1)
#fcres['gdf'].plot(ax=ax[0],facecolor='None',ls='--',edgecolor=[1,0.5,0.25],linewidth=1,label='Wildfire',alpha=1)
#atu['gdf'].plot(ax=ax[0],facecolor='None',ls='--',edgecolor=[1,0,0],linewidth=2,label='Planting',alpha=1)
#wf['gdf'].plot(ax=ax[0],facecolor='None',ls='-',edgecolor=[1,0.5,0.5],linewidth=2,label='Planting',alpha=1)
#wfp.plot(ax=ax[0],facecolor='None',edgecolor=[0,0,0],linewidth=1,label='Wildfire',alpha=1)
#atup['gdf overlay'].plot(ax=ax[0],facecolor=[0,0,0],edgecolor=[0,0.5,0],linewidth=1.25,label='Planting',alpha=0.25)
gu.PrintFig(meta['Paths']['Figures'] + '\\' + roi['Name'] + '_btm','png',900)

#%% PLOT Land Cover 2020 from CEC

def Plot_CEC_LC20(roi):

    dCEC=gu.ReadExcel(r'C:\Users\rhember\Documents\Data\BC1ha\LUTs\LUT_lcc_cec_Compressed.xlsx')

    bin=np.arange(1,dCEC['Name'].size+1,1)
    N_vis=bin.size+1
    N_hidden=1
    N_tot=N_vis+N_hidden

    z1=bin.size*np.ones(roi['grd']['lc20_cec']['Data'].shape,dtype=int)
    for i in range(bin.size):
        z1[(roi['grd']['lc20_cec']['Data']==bin[i])]=i
    z1[(roi['grd']['Data']!=1)]=i+1

    for i in range(bin.size):
        z1[0,i]=i

    # Labels
    lab=list(meta['LUT']['Derived']['lcc_cec_c'].keys())
    lab.append('')
    lab.append('')

    # Colormap
    cm=plt.cm.get_cmap('viridis',bin.size);
    for i in range(bin.size):
        cm.colors[i,0:3]=[ dCEC['R'][i]/255,dCEC['G'][i]/255,dCEC['B'][i]/255 ]
    cm=np.vstack( (cm.colors,(1,1,1,1)) )
    cm=matplotlib.colors.ListedColormap(cm)

    plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['figwidth'],(1-meta['Graphics']['sidespace'])*meta['Graphics']['figwidth']*roi['grd']['yxrat']))
    im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
    roi['gdf']['bc_bound'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
    roi['gdf']['rivers'].plot(ax=ax[0],color=[0.6,0.8,1],label='Rivers',linewidth=0.25)
    roi['gdf']['lakes'].plot(ax=ax[0],color=[0.82,0.88,1],label='Lakes',linewidth=2.25)
    roi['gdf']['bound within'].plot(ax=ax[0],color=None,edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
    roi['gdf']['road'].plot(ax=ax[0],edgecolor=[0,0,0],linewidth=0.5,label='Road',alpha=1,zorder=1)
    ax[0].set(position=meta['Graphics']['ax1 pos'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'],xticklabels='',yticklabels='')
    ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['ax1 gridvis']); ax[0].axis(meta['Graphics']['ax1 vis'])

    cb_ivl=(N_tot-1)/N_tot
    cbivls=np.arange( cb_ivl , N_tot , cb_ivl)
    cbivls_low=np.arange( 0 , N_tot , cb_ivl)
    cb=plt.colorbar(im,cax=ax[1],cmap=cm,
                    boundaries=np.arange(0,cbivls_low[N_vis],cb_ivl),
                    ticks=np.arange(cb_ivl/2,N_tot-1,cb_ivl) )
    ax[1].set(position=[0.56,0.4,0.025,0.5])

    cb.ax.set(yticklabels=lab)
    cb.ax.tick_params(labelsize=10,length=0)
    cb.outline.set_edgecolor('w')
    for i in range(cbivls.size):
        ax[1].plot([0,100],[cbivls[i],cbivls[i]],'w-',linewidth=2)
    gu.PrintFig(meta['Paths']['Figures'] + '\\' + roi['Name'] + '_lc20_cec','png',900)

    return fig,ax

fig,ax=Plot_CEC_LC20(roi)

#%% Plot deforestation from CEC

def Plot_Afforestation(meta,roi):

    lab=['Settlement to Forest','Cropland to Forest','Barren Ground to Forest','Land']
    z1=roi['grd']['aff_cec']['Data']
    z1[(roi['grd']['aff_cec']['Data']==0)]=4
    z1[(roi['grd']['Data']!=1)]=5

    # Number of colours and number of colours excluded from colorbar
    N_vis=3
    N_hidden=2
    N_tot=N_vis+N_hidden

    # Colormap
    cm=np.vstack( ((0.5,0.0,0.0,1),(1.0,0,0,1),(1.0,0.5,0,1),(0.83,0.83,0.83,1),(1,1,1,1)) )
    cm=matplotlib.colors.ListedColormap(cm)

    #plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['figwidth'],(1-meta['Graphics']['sidespace'])*meta['Graphics']['figwidth']*roi['grd']['yxrat']))
    plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['figwidth'],(1-meta['Graphics']['sidespace'])*meta['Graphics']['figwidth']*roi['grd']['yxrat']))
    im=ax[0].matshow(z1[0::1,0::1],extent=roi['grd']['Extent'],cmap=cm)
    roi['gdf']['bc_bound'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
    roi['gdf']['rivers'].plot(ax=ax[0],color=[0.6,0.8,1],label='Rivers',linewidth=0.25)
    roi['gdf']['lakes'].plot(ax=ax[0],color=[0.82,0.88,1],label='Lakes',linewidth=2.25)
    roi['gdf']['bound within'].plot(ax=ax[0],color=None,edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
    roi['gdf']['road'].plot(ax=ax[0],edgecolor='k',linewidth=0.25,label='Road',alpha=1,zorder=1)
    ax[0].set(position=meta['Graphics']['ax1 pos'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'],xticklabels='',yticklabels='')
    ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['ax1 gridvis']); ax[0].axis(meta['Graphics']['ax1 vis'])

    zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
    cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
    cb=plt.colorbar(im,cax=ax[1],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
    ax[1].set(position=[0.71,0.6,0.05,0.14])
    cb.ax.set(yticklabels=lab)
    cb.ax.tick_params(labelsize=11,length=0)
    cb.outline.set_edgecolor('w')
    for i in range(cb_bnd.size):
        ax[1].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)

    pos2=meta['Graphics']['ax2 pos']
    pos2[1]=0.7
    pos2[2]=0.027
    pos2[3]=0.14
    ax[1].set(position=pos2)
    gu.PrintFig(meta['Paths']['Figures'] + '\\' + roi['Name'] + '_afforestation_cec','png',900)

    return fig,ax

fig,ax=Plot_Afforestation(meta,roi)

#%% Plot deforestation from CEC

def Plot_Deforestation(meta,roi):

    lab=['Forest to Settlement','Forest to Cropland','Forest to Barren Ground','Land']
    z1=roi['grd']['def_cec']['Data']
    z1[(roi['grd']['def_cec']['Data']==0)]=4
    z1[(roi['grd']['Data']!=1)]=5

    # Number of colours and number of colours excluded from colorbar
    N_vis=3
    N_hidden=2
    N_tot=N_vis+N_hidden

    # Colormap
    cm=np.vstack( ((0.5,0.0,0.0,1),(1.0,0,0,1),(1.0,0.5,0,1),(0.83,0.83,0.83,1),(1,1,1,1)) )
    cm=matplotlib.colors.ListedColormap(cm)

    plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['figwidth'],(1-meta['Graphics']['sidespace'])*meta['Graphics']['figwidth']*roi['grd']['yxrat']))
    im=ax[0].matshow(z1[0::1,0::1],extent=roi['grd']['Extent'],cmap=cm)
    roi['gdf']['bc_bound'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
    roi['gdf']['rivers'].plot(ax=ax[0],color=[0.6,0.8,1],label='Rivers',linewidth=0.25)
    roi['gdf']['lakes'].plot(ax=ax[0],color=[0.82,0.88,1],label='Lakes',linewidth=2.25)
    roi['gdf']['bound within'].plot(ax=ax[0],color=None,edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
    roi['gdf']['road'].plot(ax=ax[0],edgecolor='k',linewidth=0.25,label='Road',alpha=1,zorder=1)
    ax[0].set(position=meta['Graphics']['ax1 pos'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'],xticklabels='',yticklabels='')
    ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['ax1 gridvis']); ax[0].axis(meta['Graphics']['ax1 vis'])

    zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
    cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
    cb=plt.colorbar(im,cax=ax[1],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
    ax[1].set(position=[0.71,0.6,0.05,0.14])
    cb.ax.set(yticklabels=lab)
    cb.ax.tick_params(labelsize=11,length=0)
    cb.outline.set_edgecolor('w')
    for i in range(cb_bnd.size):
        ax[1].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)

    pos2=meta['Graphics']['ax2 pos']
    pos2[1]=0.7
    pos2[2]=0.027
    pos2[3]=0.14
    ax[1].set(position=pos2)
    gu.PrintFig(meta['Paths']['Figures'] + '\\' + roi['Name'] + '_deforestation_cec','png',900)

    return fig,ax

fig,ax=Plot_Deforestation(meta,roi)

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
    fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['figwidth'],(1-meta['Graphics']['sidespace'])*meta['Graphics']['figwidth']*roi['grd']['yxrat']))
    im=ax[0].matshow(z1[0::1,0::1],extent=roi['grd']['Extent'],cmap=cm)
    roi['gdf']['bc_bound'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
    #roi['gdf']['bound within'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
    roi['gdf']['lakes'].plot(ax=ax[0],facecolor=[0.82,0.88,1],label='Lakes',linewidth=2.25)
    roi['gdf']['rivers'].plot(ax=ax[0],color=[0.6,0.8,1],label='Rivers',linewidth=0.25)
    roi['gdf']['road'].plot(ax=ax[0],edgecolor='k',linewidth=0.25,label='Road',alpha=1,zorder=1)
    ax[0].set(position=meta['Graphics']['ax1 pos'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
    ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['ax1 gridvis']); ax[0].axis(meta['Graphics']['ax1 vis'])

    zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
    cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
    cb=plt.colorbar(im,cax=ax[1],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
    ax[1].set(position=[0.71,0.6,0.05,0.14])
    cb.ax.set(yticklabels=lab)
    cb.ax.tick_params(labelsize=11,length=0)
    cb.outline.set_edgecolor('w')
    for i in range(cb_bnd.size):
        ax[1].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)

    pos2=[meta['Graphics']['leg x'],0.99-N_vis*meta['Graphics']['leg text space'],meta['Graphics']['leg w'],N_vis*meta['Graphics']['leg text space']]
    ax[1].set(position=pos2)

    gu.PrintFig(meta['Paths']['Figures'] + '\\' + roi['Name'] + '_denseclass','png',900)

    return fig,ax

Plot_DensityClass(meta,roi,gdf)

#%% Resource extraction area restrictions

def Plot_REARs(meta,roi,gdf):

    z1=roi['grd']['rears']['Data'].copy()
    z1[(roi['grd']['Data']==1) & (z1==0) & (roi['grd']['lcc1']['Data']==meta['LUT']['Derived']['lcc1']['Forest'])]=3
    z1[(roi['grd']['Data']==1) & (roi['grd']['lcc1']['Data']!=meta['LUT']['Derived']['lcc1']['Forest'])]=4
    z1[(roi['grd']['Data']==0)]=5

    lab=['Protected (completed)','Protected (proposed)','Unprotected forest','Non-forest land']

    # Number of colours and number of colours excluded from colorbar
    N_vis=4
    N_hidden=1
    N_tot=N_vis+N_hidden

    cm=np.vstack( ((0.7,0.6,1,1),(0.25,0.85,0.9,1),(0,0.4,0,1),(0.83,0.86,0.9,1),(1,1,1,1)) )
    cm=matplotlib.colors.ListedColormap(cm)

    plt.close('all')
    fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['figwidth'],(1-meta['Graphics']['sidespace'])*meta['Graphics']['figwidth']*roi['grd']['yxrat']))
    im=ax[0].matshow(z1[0::1,0::1],extent=roi['grd']['Extent'],cmap=cm)
    roi['gdf']['bc_bound'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
    #roi['gdf']['bound within'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
    roi['gdf']['lakes'].plot(ax=ax[0],facecolor=[0.82,0.88,1],label='Lakes',linewidth=2.25)
    roi['gdf']['rivers'].plot(ax=ax[0],color=[0.6,0.8,1],label='Rivers',linewidth=0.25)
    roi['gdf']['road'].plot(ax=ax[0],edgecolor='k',linewidth=0.25,label='Road',alpha=1,zorder=1)
    ax[0].set(position=meta['Graphics']['ax1 pos'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
    ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['ax1 gridvis']); ax[0].axis(meta['Graphics']['ax1 vis'])

    zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
    cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
    cb=plt.colorbar(im,cax=ax[1],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
    ax[1].set(position=[0.71,0.6,0.05,0.14])
    cb.ax.set(yticklabels=lab)
    cb.ax.tick_params(labelsize=11,length=0)
    cb.outline.set_edgecolor('w')
    for i in range(cb_bnd.size):
        ax[1].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)

    pos2=[meta['Graphics']['leg x'],0.99-N_vis*meta['Graphics']['leg text space'],meta['Graphics']['leg w'],N_vis*meta['Graphics']['leg text space']]
    ax[1].set(position=pos2)
    gu.PrintFig(meta['Paths']['Figures'] + '\\' + roi['Name'] + '_rears','png',900)

    return fig,ax

Plot_REARs(meta,roi,gdf)

#%% Plot distance from road

def Plot_DistanceFromRoad(meta,roi):

    z=roi['grd']['d2fac']['Data']
    if roi['Type']=='ByTSA':
        z[roi['grd']['Data']==0]=0

    # Plot
    plt.close('all')
    fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['figwidth'],(1-meta['Graphics']['sidespace'])*meta['Graphics']['figwidth']*roi['grd']['yxrat']))
    im=ax[0].matshow(z,extent=roi['grd']['Extent'],cmap='Greys')

    #roi['gdf']['tsa'].plot(ax=ax[0],color=None,edgecolor=[0,0,0],facecolor='none')
    ax[0].set(position=meta['Graphics']['ax1 pos'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'],xticklabels='',yticklabels='')
    ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['ax1 gridvis']); ax[0].axis(meta['Graphics']['ax1 vis'])

    # Add relief shading
    z=roi['grd']['elev']['Data']; dx,dy=roi['grd']['Cellsize'],roi['grd']['Cellsize']
    ls=LightSource(azdeg=90,altdeg=45)
    ve=0.1
    hs=ls.hillshade(z,vert_exag=ve,dx=dx,dy=dy)
    ax[0].matshow(hs,extent=roi['grd']['elev']['Extent'],cmap='Greys',alpha=0.4,clim=(np.min(hs),np.max(hs)))

    cb=plt.colorbar(im,cax=ax[1])#,boundaries=np.arange(0,N_vis-(N_hidden-1),1),ticks=np.arange(0.5,N_vis+1.5,1))
    #cb.ax.set(yticklabels=lab)
    #cb.ax.tick_params(labelsize=6,length=0)
    #for i in range(0,N_vis):
    #    ax[1].plot([0,100],[i/(N_vis-N_hidden),i/(N_vis-N_hidden)],'k-',linewidth=0.5)
    ax[1].set(position=meta['Graphics']['ax2 pos long']);

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
    fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['figwidth'],(1-meta['Graphics']['sidespace'])*meta['Graphics']['figwidth']*roi['grd']['yxrat']))
    im=ax[0].matshow(roi['grd']['SITE_INDEX']['Data'],extent=roi['grd']['Extent'],cmap='magma',clim=[5,22])

    #roi['gdf']['tsa'].plot(ax=ax[0],color=None,edgecolor=[0,0,0],facecolor='none')
    ax[0].set(position=meta['Graphics']['ax1 pos'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'],xticklabels='',yticklabels='')
    ax[0].grid(False)

    # # Add relief shading
    # z=roi['grd']['elev']['Data']; dx,dy=roi['grd']['Cellsize'],roi['grd']['Cellsize']
    # ls=LightSource(azdeg=90,altdeg=45)
    # ve=0.1
    # hs=ls.hillshade(z,vert_exag=ve,dx=dx,dy=dy)
    # ax[0].matshow(hs,extent=roi['grd']['elev']['Extent'],cmap='Greys',alpha=0.4,clim=(np.min(hs),np.max(hs)))

    cb=plt.colorbar(im,cax=ax[1])#,boundaries=np.arange(0,N_vis-(N_hidden-1),1),ticks=np.arange(0.5,N_vis+1.5,1))
    #cb.ax.set(yticklabels=lab)
    #cb.ax.tick_params(labelsize=6,length=0)
    #for i in range(0,N_vis):
    #    ax[1].plot([0,100],[i/(N_vis-N_hidden),i/(N_vis-N_hidden)],'k-',linewidth=0.5)
    ax[1].set(position=meta['Graphics']['ax2 pos long']);

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
    plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['figwidth'],(1-meta['Graphics']['sidespace'])*meta['Graphics']['figwidth']*roi['grd']['yxrat']))
    im=ax[0].matshow(z1[0::1,0::1],clim=(0,L+1),extent=roi['grd']['Extent'],cmap=cm)
    roi['gdf']['bound within'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
    roi['gdf']['lakes'].plot(ax=ax[0],facecolor=[0.82,0.88,1],label='Lakes',linewidth=2.25)
    roi['gdf']['rivers'].plot(ax=ax[0],color=[0.6,0.8,1],label='Rivers',linewidth=0.25)
    roi['gdf']['road'].plot(ax=ax[0],edgecolor=[1,0.6,0],linewidth=0.25,label='Road',alpha=1,zorder=1)
    ax[0].set(position=meta['Graphics']['ax1 pos'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
    ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['ax1 gridvis']); ax[0].axis(meta['Graphics']['ax1 vis'])

    cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_vis-(N_hidden-1),1),ticks=np.arange(0.5,N_vis-N_hidden,1))
    cb.ax.set(yticklabels=lab)
    cb.ax.tick_params(labelsize=10,length=0)
    cb.outline.set_edgecolor('w')
    for i in range(0,N_vis):
        ax[1].plot([0,100],[i,i],'w-',linewidth=1.5)
    ax[1].set(position=meta['Graphics']['ax2 pos long']);

    return fig,ax

plt.close('all')
fig,ax=Plot_PFI(meta,roi)

#gdf['cities'].plot(ax=ax[0],marker='s',edgecolor=[0.75,0.3,0],facecolor=[1,0.6,0],lw=1,markersize=20,alpha=1,zorder=2)
#for x,y,label in zip(gdf['cities'].geometry.x,gdf['cities'].geometry.y,gdf['cities'].Name):
#    ax[0].annotate(label,xy=(x,y),xytext=(4,3),textcoords="offset points",color=[0.75,0.3,0])
#roi['gdf']['cc']['gdf'].plot(ax=ax[0],facecolor=None,edgecolor=[0.5,1,1],linewidth=0.25,label='Opening',alpha=0.3)

gu.PrintFig(meta['Paths']['Figures'] + '\\' + roi['Name'] + '_pfi_c','png',900)

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
    plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['figwidth'],(1-meta['Graphics']['sidespace'])*meta['Graphics']['figwidth']*roi['grd']['yxrat']))
    im=ax[0].matshow(z1,clim=(0,L+1),extent=roi['grd']['Extent'],cmap=cm)
    roi['gdf']['bound within'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
    roi['gdf']['lakes'].plot(ax=ax[0],facecolor=[0.82,0.88,1],label='Lakes',linewidth=2.25)
    roi['gdf']['rivers'].plot(ax=ax[0],color=[0.6,0.8,1],label='Rivers',linewidth=0.25)
    roi['gdf']['road'].plot(ax=ax[0],edgecolor=[1,0.6,0],linewidth=1,label='Road',alpha=1,zorder=1)
    ax[0].set(position=meta['Graphics']['ax1 pos'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
    ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['ax1 gridvis']); ax[0].axis(meta['Graphics']['ax1 vis'])

    cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_vis-(N_hidden-1),1),ticks=np.arange(0.5,N_vis-N_hidden,1))
    cb.ax.set(yticklabels=lab)
    cb.ax.tick_params(labelsize=10,length=0)
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

def Plot_GlobBio(meta,roi):

    # Grid
    bw=20; bin=np.arange(0,200+bw,bw);
    z1=(bin.size)*np.ones( roi['grd']['globbio']['Data'].shape)
    for i in range(bin.size):
        ind=np.where(np.abs( roi['grd']['globbio']['Data']-bin[i])<=bw/2)
        if ind[0].size>0:
            z1[ind]=i
        else:
            z1[0,i]=i
    ind=np.where(roi['grd']['globbio']['Data']>bin[i]); z1[ind]=i
    z1[(roi['grd']['Data']==1) & ( roi['grd']['globbio']['Data']==0)]=i+1
    z1[(roi['grd']['Data']!=1)]=i+2
    L=i+2

    lab=bin.astype(str)

    # Colormap
    cm=plt.cm.get_cmap('viridis',i)
    cm=np.vstack( (cm.colors,(0.9,0.9,0.9,1),(1,1,1,1)) )
    cm=matplotlib.colors.ListedColormap(cm)

    N_vis=bin.size+3
    N_hidden=3

    # Plot
    plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['figwidth'],(1-meta['Graphics']['sidespace'])*meta['Graphics']['figwidth']*roi['grd']['yxrat']))
    im=ax[0].matshow(z1[0::1,0::1],clim=(0,L+1),extent=roi['grd']['Extent'],cmap=cm)
    roi['gdf']['bound within'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
    roi['gdf']['lakes'].plot(ax=ax[0],facecolor=[0.82,0.88,1],label='Lakes',linewidth=2.25)
    roi['gdf']['rivers'].plot(ax=ax[0],color=[0.6,0.8,1],label='Rivers',linewidth=0.25)
    roi['gdf']['road'].plot(ax=ax[0],edgecolor='k',linewidth=0.25,label='Road',alpha=1,zorder=1)

    ax[0].set(position=meta['Graphics']['ax1 pos'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
    ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['ax1 gridvis']); ax[0].axis(meta['Graphics']['ax1 vis'])

    cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_vis-(N_hidden-1),1),ticks=np.arange(0.5,N_vis-N_hidden,1))
    cb.ax.set(yticklabels=lab)
    cb.ax.tick_params(labelsize=10,length=0)
    cb.outline.set_edgecolor('w')
    for i in range(0,N_vis):
        ax[1].plot([0,100],[i,i],'w-',linewidth=1.5)
    ax[1].set(position=meta['Graphics']['ax2 pos long']);

    gu.PrintFig(meta['Paths']['Figures'] + '\\' + roi['Name'] + '_globbio','png',900)

    return fig,ax

fig,ax=Plot_GlobBio(meta,roi)

#%% PLOT soc

def Plot_SOC_WithinROI(roi):

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
    plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['figwidth'],(1-meta['Graphics']['sidespace'])*meta['Graphics']['figwidth']*roi['grd']['yxrat']))
    im=ax[0].matshow(z1[0::1,0::1],clim=(0,L+1),extent=roi['grd']['Extent'],cmap=cm)
    roi['gdf']['bound within'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
    roi['gdf']['lakes'].plot(ax=ax[0],facecolor=[0.82,0.88,1],label='Lakes',linewidth=2.25)
    roi['gdf']['rivers'].plot(ax=ax[0],color=[0.6,0.8,1],label='Rivers',linewidth=0.25)
    roi['gdf']['road'].plot(ax=ax[0],edgecolor='k',linewidth=0.25,label='Road',alpha=1,zorder=1)
    ax[0].set(position=meta['Graphics']['ax1 pos'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
    ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['ax1 gridvis']); ax[0].axis(meta['Graphics']['ax1 vis'])

    cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_vis-(N_hidden-1),1),ticks=np.arange(0.5,N_vis-N_hidden,1))
    cb.ax.set(yticklabels=lab)
    cb.ax.tick_params(labelsize=10,length=0)
    cb.outline.set_edgecolor('w')
    for i in range(0,N_vis):
        ax[1].plot([0,100],[i,i],'w-',linewidth=1.5)

    #ax[1].set(position=meta['Graphics']['ax2 pos long']);
    pos2=copy.copy(meta['Graphics']['ax2 pos'])
    pos2[1]=0.55
    pos2[2]=0.024
    pos2[3]=0.4
    ax[1].set(position=pos2)
    gu.PrintFig(meta['Paths']['Figures'] + '\\' + roi['Name'] + '_gsoc','png',900)

    return fig,ax

fig,ax=Plot_SOC_WithinROI(roi)

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
    plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['figwidth'],(1-meta['Graphics']['sidespace'])*meta['Graphics']['figwidth']*roi['grd']['yxrat']))
    im=ax[0].matshow(z1[0::1,0::1],clim=(0,L+1),extent=roi['grd']['Extent'],cmap=cm)
    ax[0].set(position=meta['Graphics']['ax1 pos'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'],xticklabels='',yticklabels='')
    ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['ax1 gridvis']); ax[0].axis(meta['Graphics']['ax1 vis'])

    cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_vis-(N_hidden-1),1),ticks=np.arange(0.5,N_vis-N_hidden,1))
    cb.ax.set(yticklabels=lab)
    cb.ax.tick_params(labelsize=10,length=0)
    cb.outline.set_edgecolor('w')
    for i in range(0,N_vis):
        ax[1].plot([0,100],[i,i],'w-',linewidth=1.5)
    pos2=meta['Graphics']['ax2 pos']
    pos2[1]=0.25
    pos2[2]=0.027
    pos2[3]=0.5
    ax[1].set(position=pos2)

    gu.PrintFig(meta['Paths']['Figures'] + '\\' + roi['Name'] + '_crownc','png',900)

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
    plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['figwidth'],(1-meta['Graphics']['sidespace'])*meta['Graphics']['figwidth']*roi['grd']['yxrat']))
    im=ax[0].matshow(z1[0::1,0::1],clim=(0,L+1),extent=roi['grd']['Extent'],cmap=cm)
    #roi['gdf']['tsa'].plot(ax=ax[0],color=None,edgecolor=[0,0,0],facecolor='none')
    ax[0].set(position=meta['Graphics']['ax1 pos'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
    ax[0].grid(False)

    cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_vis-(N_hidden-1),1),ticks=np.arange(0.5,N_vis-N_hidden,1))
    cb.ax.set(yticklabels=lab)
    cb.ax.tick_params(labelsize=6,length=0)
    for i in range(0,N_vis):
        ax[1].plot([0,100],[i/(N_vis-N_hidden),i/(N_vis-N_hidden)],'k-',linewidth=0.5)
    ax[1].set(position=meta['Graphics']['ax2 pos long']);

    # cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_vis-(N_hidden-1),1),ticks=np.arange(0.5,N_vis+1.5,1))
    # cb.ax.set(yticklabels=lab)
    # cb.ax.tick_params(labelsize=6,length=0)
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
#gu.PrintFig(meta['Paths']['Figures'] + '\\sphlive','png',300)

#%% Plot harvested year within the TSA mask

def Plot_HarvestYear(roi):

    z0=roi['grd']['cut_yr_con1']['Data']
    #z0=roi['grd']['cut_yr']['Data']
    #z0=roi['grd']['cut_yr_ntems']['Data']

    bw=10; bin=np.arange(1855,2025+bw,bw);
    z1=(bin.size)*np.ones( z0.shape)
    for i in range(bin.size):
        ind=np.where(np.abs( z0-bin[i])<=bw/2)
        z1[ind]=i
    z1[(roi['grd']['Data']==0)]=i+2
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
    fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['figwidth'],(1-meta['Graphics']['sidespace'])*meta['Graphics']['figwidth']*roi['grd']['yxrat']))
    im=ax[0].matshow(z1[0::1,0::1],clim=(0,L+1),extent=roi['grd']['Extent'],cmap=cm)
    #roi['gdf']['bound within'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
    roi['gdf']['bc_bound'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
    roi['gdf']['lakes'].plot(ax=ax[0],facecolor=[0.82,0.88,1],label='Lakes',linewidth=2.25)
    roi['gdf']['rivers'].plot(ax=ax[0],color=[0.6,0.8,1],label='Rivers',linewidth=0.25)
    roi['gdf']['road'].plot(ax=ax[0],edgecolor='k',linewidth=0.25,label='Road',alpha=1,zorder=1)
    ax[0].set(position=meta['Graphics']['ax1 pos'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
    ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['ax1 gridvis']); ax[0].axis(meta['Graphics']['ax1 vis'])

    cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_vis-(N_hidden-1),1),ticks=np.arange(0.5,N_vis-N_hidden,1))
    cb.ax.set(yticklabels=lab)

    cb.ax.tick_params(labelsize=10,length=0)
    cb.outline.set_edgecolor('w')
    for i in range(0,N_vis):
        ax[1].plot([0,100],[i,i],'w-',linewidth=1.5)

    pos2=[meta['Graphics']['leg x'],0.99-N_vis*meta['Graphics']['leg text space'],meta['Graphics']['leg w'],N_vis*meta['Graphics']['leg text space']]
    ax[1].set(position=pos2)

    gu.PrintFig(meta['Paths']['Figures'] + '\\' + roi['Name'] + '_cut_yr_con1','png',900)
    #gu.PrintFig(meta['Paths']['Figures'] + '\\' + roi['Name'] + '_cut_yr','png',900)
    #gu.PrintFig(meta['Paths']['Figures'] + '\\' + roi['Name'] + '_cut_yr_ntems','png',900)

    return fig,ax

fig,ax=Plot_HarvestYear(roi)

#%% Plot wildfire year

def Plot_FireYear(roi):

    # Grid
    bw=5; bin=np.arange(1915,2025,bw);
    z1=(bin.size)*np.ones( roi['grd']['fire_yr']['Data'].shape)
    for i in range(bin.size):
        ind=np.where(np.abs( roi['grd']['fire_yr']['Data']-bin[i])<=bw/2)
        z1[ind]=i
    z1[(roi['grd']['Data']==0)]=i+2
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
    fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['figwidth'],(1-meta['Graphics']['sidespace'])*meta['Graphics']['figwidth']*roi['grd']['yxrat']))
    im=ax[0].matshow(z1[0::1,0::1],clim=(0,L+1),extent=roi['grd']['Extent'],cmap=cm)
    roi['gdf']['bound within'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
    roi['gdf']['bc_bound'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
    roi['gdf']['lakes'].plot(ax=ax[0],facecolor=[0.82,0.88,1],label='Lakes',linewidth=2.25)
    roi['gdf']['rivers'].plot(ax=ax[0],color=[0.6,0.8,1],label='Rivers',linewidth=0.25)
    roi['gdf']['road'].plot(ax=ax[0],edgecolor='k',linewidth=0.25,label='Road',alpha=1,zorder=1)
    ax[0].set(position=meta['Graphics']['ax1 pos'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
    ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['ax1 gridvis']); ax[0].axis(meta['Graphics']['ax1 vis'])

    cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_vis-(N_hidden-1),1),ticks=np.arange(0.5,N_vis-N_hidden,1))
    cb.ax.set(yticklabels=lab)

    cb.ax.tick_params(labelsize=10,length=0)
    cb.outline.set_edgecolor('w')
    for i in range(0,N_vis):
        ax[1].plot([0,100],[i,i],'w-',linewidth=1.5)

    pos2=copy.copy(meta['Graphics']['ax2 pos'])
    pos2[1]=0.1
    pos2[2]=0.03
    pos2[3]=0.8
    ax[1].set(position=pos2)
    gu.PrintFig(meta['Paths']['Figures'] + '\\' + roi['Name'] + '_fire_yr','png',900)

    return fig,ax

fig,ax=Plot_FireYear(roi)

#%% Plot Global Forest Change Loss Year

def Plot_LossYear(roi):

    # Grid
    bw=1; bin=np.arange(1,22,bw);
    z1=(bin.size)*np.ones( roi['grd']['gfcly']['Data'].shape)
    for i in range(bin.size):
        ind=np.where(roi['grd']['gfcly']['Data']==bin[i])
        z1[ind]=i
    z1[(roi['grd']['gfcly']['Data']==0)]=i+1
    z1[(roi['grd']['Data']==0)]=i+2
    #z1[(roi['grd']['Data']!=1)]=i+2
    L=i+2

    lab=np.array(2000+bin).astype(str)

    # Colormap
    #cm=plt.cm.get_cmap('viridis',i)
    cm=plt.cm.get_cmap('plasma',i)
    cm=np.vstack( (cm.colors,(0.85,0.85,0.85,1),(1,1,1,1)) )
    cm=matplotlib.colors.ListedColormap(cm)

    N_vis=bin.size+3
    N_hidden=3

    plt.close('all');
    fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['figwidth'],(1-meta['Graphics']['sidespace'])*meta['Graphics']['figwidth']*roi['grd']['yxrat']))
    im=ax[0].matshow(z1[0::1,0::1],clim=(0,L+1),extent=roi['grd']['Extent'],cmap=cm)
    roi['gdf']['bound within'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
    roi['gdf']['bc_bound'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
    roi['gdf']['lakes'].plot(ax=ax[0],facecolor=[0.82,0.88,1],label='Lakes',linewidth=2.25)
    roi['gdf']['rivers'].plot(ax=ax[0],color=[0.6,0.8,1],label='Rivers',linewidth=0.25)
    #roi['gdf']['road'].plot(ax=ax[0],edgecolor='k',linewidth=0.25,label='Road',alpha=1,zorder=1)
    ax[0].set(position=meta['Graphics']['ax1 pos'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'],xticklabels='',yticklabels='')
    ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['ax1 gridvis']); ax[0].axis(meta['Graphics']['ax1 vis'])

    cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_vis-(N_hidden-1),1),ticks=np.arange(0.5,N_vis-N_hidden,1))
    cb.ax.set(yticklabels=lab)

    cb.ax.tick_params(labelsize=10,length=0)
    cb.outline.set_edgecolor('w')
    for i in range(0,N_vis):
        ax[1].plot([0,100],[i,i],'w-',linewidth=1.5)

    pos2=copy.copy(meta['Graphics']['ax2 pos'])
    pos2[1]=0.05
    pos2[2]=0.03
    pos2[3]=0.9
    ax[1].set(position=pos2)

    return fig,ax

plt.close('all')
fig,ax=Plot_LossYear(roi)
#roi['gdf']['wf']['gdf'].plot(ax=ax[0],facecolor=[1,0.6,0],edgecolor=[0.75,0,0],linewidth=0.25,label='Opening',alpha=0.3)
#roi['gdf']['cc']['gdf'].plot(ax=ax[0],facecolor=[0.5,0.75,1],edgecolor=[0.5,0.75,1],linewidth=0.25,label='Opening',alpha=0.3)

#gdf['cities'].plot(ax=ax[0],marker='s',edgecolor=[0.75,0.3,0],facecolor=[1,0.6,0],lw=1,markersize=20,alpha=1,zorder=2)
#for x,y,label in zip(gdf['cities'].geometry.x,gdf['cities'].geometry.y,gdf['cities'].Name):
#    ax[0].annotate(label,xy=(x,y),xytext=(4,3),textcoords="offset points",color=[0.75,0.3,0])

gu.PrintFig(meta['Paths']['Figures'] + '\\' + roi['Name'] + '_gfcly','png',900)

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

    fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['figwidth'],(1-meta['Graphics']['sidespace'])*meta['Graphics']['figwidth']*roi['grd']['yxrat']))
    im=ax[0].matshow(z1,clim=(0,N_vis),extent=lc2['grd']['Extent'],cmap=cm)
    #bm['gdf_bc_bound'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
    roi['gdf_lakes'].plot(ax=ax[0],facecolor=[0.82,0.88,1],edgecolor=[0.7*0.82,0.7*0.88,0.7*1],linewidth=0.25,label='Water')
    roi['gdf_rivers'].plot(ax=ax[0],linecolor=[0,0,0.7],label='Water',linewidth=0.25)
    roi['gdf_bound'].plot(ax=ax[0],color=None,edgecolor=[0,0,0],facecolor='none')
    roi['gdf_roads'].plot(ax=ax[0],facecolor='none',edgecolor=[0,0,0],label='Roads',linewidth=0.75,alpha=1,zorder=1)
    ax[0].set(position=meta['Graphics']['ax1 pos'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
    ax[0].grid(False)

    cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_vis,1),ticks=np.arange(0.5,N_vis+1.5,1))
    cb.ax.set(yticklabels=lab)
    cb.ax.tick_params(labelsize=6,length=0)
    for i in range(0,N_vis):
        ax[1].plot([0,100],[i/(N_vis-N_hidden-1),i/(N_vis-N_hidden-1)],'k-',linewidth=0.5)
    pos2=meta['Graphics']['ax2 pos']
    pos2[1]=0.8
    pos2[3]=0.14
    ax[1].set(position=pos2)

    return fig,ax

plt.close('all')
fig,ax=Plot_ROI_BSR(bsr)
#wfp.plot(ax=ax[0],facecolor='None',edgecolor=[0,0,0],linewidth=1.5,label='Wildfire',alpha=1)
#atup['gdf overlay'].plot(ax=ax[0],facecolor=[0,0,0],edgecolor=[0,0.5,0],linewidth=1.25,label='Planting',alpha=0.25)
#gu.PrintFig(meta['Paths']['Figures'] + '\\Planted areas and BSR','png',300)


#%% Plot soil water content

def Plot_SoilWaterContent(roi):

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

    plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['figwidth'],(1-meta['Graphics']['sidespace'])*meta['Graphics']['figwidth']*roi['grd']['yxrat']))
    im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
    roi['gdf']['bound within'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
    roi['gdf']['lakes'].plot(ax=ax[0],facecolor=[0.82,0.88,1],label='Lakes',linewidth=2.25)
    roi['gdf']['rivers'].plot(ax=ax[0],color=[0.6,0.8,1],label='Rivers',linewidth=0.25)
    roi['gdf']['road'].plot(ax=ax[0],edgecolor='k',linewidth=0.25,label='Road',alpha=1,zorder=1)
    ax[0].set(position=meta['Graphics']['ax1 pos'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
    ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['ax1 gridvis']); ax[0].axis(meta['Graphics']['ax1 vis'])

    zmn=np.min(z1);
    zmx=np.max(z1);
    cb_ivl=(zmx-zmn)/N_tot;
    cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
    cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
    cb=plt.colorbar(im,cax=ax[1],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
    ax[1].set(position=[0.71,0.6,0.05,0.14])
    cb.ax.set(yticklabels=lab)
    cb.ax.tick_params(labelsize=11,length=0)
    cb.outline.set_edgecolor('w')
    for i in range(cb_bnd.size):
        ax[1].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)

    pos2=meta['Graphics']['ax2 pos']
    pos2[1]=0.45
    pos2[2]=0.027
    pos2[3]=0.45
    ax[1].set(position=pos2)

    gu.PrintFig(meta['Paths']['Figures'] + '\\' + roi['Name'] + '_SoilWaterContent_Normal','png',300)

    return fig,ax

fig,ax=Plot_SoilWaterContent(roi)

#%%

#%% Plot ground plot sample grid

tmp=roi['grd']['Data'].copy()
ind=np.where(roi['grd']['lcc1']['Data']==meta['LUT']['lcc1']['Forest'])
tmp[ind]=2

cm=np.vstack( ((1,1,1,1),(0.95,0.95,0.95,1),(0.8,0.8,0.8,1)) )
cm=matplotlib.colors.ListedColormap(cm)

ms=2
plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(meta['Graphics']['figwidth'],(1-meta['Graphics']['sidespace'])*meta['Graphics']['figwidth']*roi['grd']['yxrat']))
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

ax.set(position=meta['Graphics']['ax1 pos'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'],aspect='auto')
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both')
ax.grid(meta['Graphics']['ax1 gridvis'])
ax.axis(meta['Graphics']['ax1 vis'])
#gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Biomass\SamplingPowerMap2','png',900)

#%% Salvage mask from timber cruise

def Plot_SalvageMask(meta,roi):

    lab=['Salvage (>50% dead)','Salvage (>10% dead)','Harvest (<10% dead)','No harvesting','Non-forest land','Outside']
    lab=lab[0:-1]
    z1=roi['grd']['salvage']['Data']
    z1[roi['grd']['Data']==0]=6

    # Number of colours and number of colours excluded from colorbar
    N_vis=4
    N_hidden=2
    N_tot=N_vis+N_hidden

    # Colormap
    cm=np.vstack( ((0.9,0,0,1),(1,0.85,0,1),(0.85,1,0.65,1),(0.6,0.6,0.6,1), (0.93,0.93,0.93,1),(1,1,1,1)) )
    cm=matplotlib.colors.ListedColormap(cm)

    plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['figwidth'],(1-meta['Graphics']['sidespace'])*meta['Graphics']['figwidth']*roi['grd']['yxrat']))
    im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
    roi['gdf']['rivers'].plot(ax=ax[0],color=[0.6,0.8,1],label='Rivers',linewidth=0.25)
    roi['gdf']['lakes'].plot(ax=ax[0],color=[0.82,0.88,1],label='Lakes',linewidth=2.25)
    roi['gdf']['bound within'].plot(ax=ax[0],color=None,edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
    roi['gdf']['road'].plot(ax=ax[0],edgecolor='k',linewidth=0.5,label='Road',alpha=1,zorder=1)
    ax[0].set(position=meta['Graphics']['ax1 pos'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'],xticklabels='',yticklabels='')
    ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['ax1 gridvis']); ax[0].axis(meta['Graphics']['ax1 vis'])

    zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
    cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
    cb=plt.colorbar(im,cax=ax[1],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
    ax[1].set(position=[0.71,0.6,0.05,0.14])
    cb.ax.set(yticklabels=lab)
    cb.ax.tick_params(labelsize=11,length=0)
    cb.outline.set_edgecolor('w')
    for i in range(cb_bnd.size):
        ax[1].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)

    pos2=meta['Graphics']['ax2 pos']
    pos2[1]=0.7
    pos2[2]=0.027
    pos2[3]=0.22
    ax[1].set(position=pos2)

    gu.PrintFig(meta['Paths']['Figures'] + '\\' + roi['Name'] + '_salvage_mask','png',900)

    return fig,ax

fig,ax=Plot_SalvageMask(meta,roi)

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

    plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['figwidth'],(1-meta['Graphics']['sidespace'])*meta['Graphics']['figwidth']*roi['grd']['yxrat']))
    im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
    roi['gdf']['rivers'].plot(ax=ax[0],color=[0.6,0.8,1],label='Rivers',linewidth=0.25)
    roi['gdf']['lakes'].plot(ax=ax[0],color=[0.82,0.88,1],label='Lakes',linewidth=2.25)
    roi['gdf']['bound within'].plot(ax=ax[0],color=None,edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
    roi['gdf']['road'].plot(ax=ax[0],edgecolor='k',linewidth=0.5,label='Road',alpha=1,zorder=1)
    ax[0].set(position=meta['Graphics']['ax1 pos'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'],xticklabels='',yticklabels='')
    ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['ax1 gridvis']); ax[0].axis(meta['Graphics']['ax1 vis'])

    zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
    cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
    cb=plt.colorbar(im,cax=ax[1],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
    ax[1].set(position=[0.71,0.6,0.05,0.14])
    cb.ax.set(yticklabels=lab)
    cb.ax.tick_params(labelsize=11,length=0)
    cb.outline.set_edgecolor('w')
    for i in range(cb_bnd.size):
        ax[1].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)

    pos2=meta['Graphics']['ax2 pos']
    pos2[1]=0.7
    pos2[2]=0.03
    pos2[3]=0.2
    ax[1].set(position=pos2)

    gu.PrintFig(meta['Paths']['Figures'] + '\\' + roi['Name'] + '_range_consol','png',900)

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

    plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['figwidth'],(1-meta['Graphics']['sidespace'])*meta['Graphics']['figwidth']*roi['grd']['yxrat']))
    im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
    roi['gdf']['rivers'].plot(ax=ax[0],color=[0.6,0.8,1],label='Rivers',linewidth=0.25)
    roi['gdf']['lakes'].plot(ax=ax[0],color=[0.82,0.88,1],label='Lakes',linewidth=2.25)
    roi['gdf']['bound within'].plot(ax=ax[0],color=None,edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
    roi['gdf']['road'].plot(ax=ax[0],edgecolor='k',linewidth=0.5,label='Road',alpha=1,zorder=1)
    ax[0].set(position=meta['Graphics']['ax1 pos'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'],xticklabels='',yticklabels='')
    ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['ax1 gridvis']); ax[0].axis(meta['Graphics']['ax1 vis'])

    zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
    cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
    cb=plt.colorbar(im,cax=ax[1],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
    ax[1].set(position=[0.71,0.6,0.05,0.14])
    cb.ax.set(yticklabels=lab)
    cb.ax.tick_params(labelsize=11,length=0)
    cb.outline.set_edgecolor('w')
    for i in range(cb_bnd.size):
        ax[1].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)

    pos2=meta['Graphics']['ax2 pos']
    pos2[1]=0.7
    pos2[2]=0.027
    pos2[3]=0.24
    ax[1].set(position=pos2)

    gu.PrintFig(meta['Paths']['Figures'] + '\\' + roi['Name'] + '_reserves_consol','png',900)

    return fig,ax

fig,ax=Plot_ReservesConsolidated(meta,roi)