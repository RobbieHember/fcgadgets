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
import fcgadgets.bc1ha.bc1ha_plot as p1ha
import fcgadgets.macgyver.util_query_gdb as qgdb

#%% Import parameters
meta=u1ha.Init()
meta['Graphics']['Map']['RGSF']=1
meta['Graphics']['Map']['Fig Width']=15.5
meta['Graphics']['Map']['Side Space']=0.25
meta['Graphics']['Map']['Map Position']=[0,0,1-meta['Graphics']['Map']['Side Space']-0.01,1]
meta['Graphics']['Map']['Map Axis Vis']='off'
meta['Graphics']['Map']['Map Grid Vis']=False
#meta['Graphics']['Map']['Legend X']=1-meta['Graphics']['Map']['Side Space']+meta['Graphics']['Map']['Map Position'][0]#+0.01,0.6,0.03,0.35]
#meta['Graphics']['Map']['Legend X']=0.55 # Good for province
meta['Graphics']['Map']['Legend Width']=0.0275
meta['Graphics']['Map']['Legend Font Size']=7
meta['Graphics']['Map']['Legend Text Space']=0.035
meta['Graphics']['Map']['Show Bound Land Mask']='On'
meta['Graphics']['Map']['Show Bound Within']='Off'
meta['Graphics']['Map']['Show Lakes']='Off'
meta['Graphics']['Map']['Show Rivers']='On'
meta['Graphics']['Map']['Show Rail']='On'
meta['Graphics']['Map']['Show Roads']='On'
meta['Graphics']['Map']['Show TPF and Cities']='On'
meta['Graphics']['Map']['Show Symbol Labels']='Off'

meta['Graphics']['Plot Style']='Manuscript'
meta['Graphics']['gp']=gu.SetGraphics(meta['Graphics']['Plot Style'])
meta['Graphics']['Print Figures']='On'

#%% Define region of interest
roi={}
#roi['Type']='Prov'; roi['Name']='Prov'

#roi['Type']='ByRegDis'
#roi['Name']='CAPITAL'
#roi['Name']='STRATHCONA'
#roi['Name']='COMOX VALLEY'
#roi['Name']='COWICHAN VALLEY'

#roi['Type']='ByTSA'; roi['Name']='Arrowsmith TSA'
#roi['Type']='ByTSA'; roi['Name']='Fort St John TSA'
#roi['Type']='ByTSA'; roi['Name']='Boundary TSA'
#roi['Type']='ByTSA'; roi['Name']='Cassiar TSA'
#roi['Type']='ByTSA'; roi['Name']='Merritt TSA'
#roi['Type']='ByTSA'; roi['Name']='Fort Nelson TSA'
#roi['Type']='ByTSA'; roi['Name']='Prince George TSA'
roi['Type']='ByTSA'; roi['Name']='Williams Lake TSA'
#roi['Type']='ByTSA'; roi['Name']='North Island TSA'
#roi['Type']='ByTSA'; roi['Name']='Quesnel TSA'
#roi['Type']='ByTSA'; roi['Name']='100 Mile House TSA'

#roi['Type']='LICS'
if roi['Type']=='ByTSA':
    meta['Graphics']['Print Figure Path']=r'C:\Users\rhember\OneDrive - Government of BC\Figures' + '\\TSAs\\' + roi['Name']
elif roi['Type']=='ByRegDis':
    meta['Graphics']['Print Figure Path']=r'C:\Users\rhember\OneDrive - Government of BC\Figures' + '\\Regional Districts\\' + roi['Name']
elif roi['Type']=='Prov':
    meta['Graphics']['Print Figure Path']=r'C:\Users\rhember\OneDrive - Government of BC\Figures' + '\\Province\\' + roi['Name']
elif roi['Type']=='LICS':
    meta['Graphics']['Print Figure Path']=r'C:\Users\rhember\OneDrive - Government of BC\Figures\\LICS'

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

elif roi['Type']=='ByRegDis':
    roi['List']=roi['Name']
    # Search: gdf['regdis']['gdf']['REGIONAL_DISTRICT_NAME'].unique()
    #roi['List']=['CAPITAL']
    #roi['List']=['COWICHAN VALLEY']
    #roi['List']=['STRATHCONA']
elif roi['Type']=='LICS':
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
        lics=gu.ipickle(r'D:\Data\FCI_Projects\LICS Site List\LICS Sites.pkl')
        id=5; 
        roi['Name']=lics[id]['Loc']
        roi['Centre']=[lics[id]['Lon'],lics[id]['Lat']]; roi['Name']=lics[id]['Loc']
        roi['Radius']=15*1000 # metres

elif roi['Type']=='Prov':
    pass

# Import base maps
gdf=u1ha.Import_GDBs_ProvinceWide(meta)

# Prepare region of interest
roi=u1ha.DefineROI(meta,roi,gdf)
t1=time.time()
print((t1-t0)/60)

#%% Import rasters over ROI
vList=u1ha.GetRasterListFromSpreadsheet(r'C:\Users\rhember\Documents\Data\BC1ha\RasterInclusion.xlsx')
roi=u1ha.Import_Raster(meta,roi,vList)

#%% Forest mask for costum vector layers
fig,ax=p1ha.Plot_ForestMask(meta,roi)

#%% Import required vector geodatabases
#vList=['cc','fcres']
#vList=['wf'] # 'op','cc','fcres','ogsr'
#roi=u1ha.Import_GDB_Over_ROI(meta,roi,vList)

#%% Land Cover Compilation 1
fig,ax=p1ha.Plot_LandCoverComp1(meta,roi,'lc_comp1_1800')
fig,ax=p1ha.Plot_LandCoverComp1(meta,roi,'lc_comp1_2019')
fig,ax=p1ha.Plot_LandCoverComp1(meta,roi,'lc_comp1_2049s1')
fig,ax=p1ha.Plot_LandCoverComp1(meta,roi,'lc_comp1_2049s2')
fig,ax=p1ha.Plot_LandCoverComp1(meta,roi,'lc_comp1_2049s3')
fig,ax=p1ha.Plot_LandCoverComp1(meta,roi,'lc_comp1_2049s4')

#%% Land Use Compilation 1
fig,ax=p1ha.Plot_LandUseComp1(meta,roi,'lu_comp1_2019')
fig,ax=p1ha.Plot_LandUseComp1(meta,roi,'lu_comp1_2049s1')
fig,ax=p1ha.Plot_LandUseComp1(meta,roi,'lu_comp1_2049s2')
fig,ax=p1ha.Plot_LandUseComp1(meta,roi,'lu_comp1_2049s3')
fig,ax=p1ha.Plot_LandUseComp1(meta,roi,'lu_comp1_2049s4')
p1ha.Plot_LandUseComp1Panels(meta,roi) # Plot panels

#%% Other Land Cover estimates
fig,ax=p1ha.Plot_LandCoverComp1(meta,roi,'lc_ntems_2019_recl')
fig,ax=p1ha.Plot_LandCoverComp1(meta,roi,'lc_vri_recl')
fig,ax=p1ha.Plot_LC20_CEC(meta,roi)

#%% Plot LUC (Afforestation or Deforestation) from CEC
fig,ax=p1ha.Plot_LandCoverLandUseChange(meta,roi,'luc1_early_type')
fig,ax=p1ha.Plot_LandCoverLandUseChange(meta,roi,'luc1_1019_type')
#fig,ax=p1ha.Plot_LU_Change_FromCEC(meta,roi,'aff')
#fig,ax=p1ha.Plot_LU_Change_FromCEC(meta,roi,'def')

#%% Plot Land Use Change year
fig,ax=p1ha.Plot_LUC_Year(meta,roi)

#%% Plot elevation
fig,ax=p1ha.Plot_Elev(meta,roi)

#%% Plot BGC Zones
fig,ax=p1ha.Plot_BGC_Zone(meta,roi)

#%% Upland-wetland forest mask
fig,ax=p1ha.Plot_UplandWetlandForest(meta,roi,'upwetf_ntems')
fig,ax=p1ha.Plot_UplandWetlandForest(meta,roi,'upwetf_vri')

#%% Tree density class
p1ha.Plot_TreeDensityClass(meta,roi)

#%% Species leading NTEMS
fig,ax=p1ha.Plot_Spc1(meta,roi)

#%% Plot mean annual temp
fig,ax=p1ha.Plot_MAT(meta,roi)

#%% Plot mean annual precip
fig,ax=p1ha.Plot_MAP(meta,roi)

#%% Plot age from VRI
fig,ax=p1ha.Plot_PROJ_AGE_1(meta,roi)

#%% Plot age from NTEMS
fig,ax=p1ha.Plot_Age_NTEMS(meta,roi)

#%% Plot site index from VRI
fig,ax=p1ha.Plot_SI(meta,roi)

#%% Plot Crown Cover Percent from VRI
fig,ax=p1ha.Plot_CrownCover(meta,roi)

#%% Plot Live SPH from VRI
fig,ax=p1ha.Plot_SPH_Live(meta,roi)

#%% Plot wildfire year
fig,ax=p1ha.Plot_WildfireYear(meta,roi)

#%% Plot harvest year
fig,ax=p1ha.Plot_HarvestYear(meta,roi,'comp1')
fig,ax=p1ha.Plot_HarvestYear(meta,roi,'comp2')

#%% Plot planted area mask
fig,ax=p1ha.Plot_PlantedMask(meta,roi)

#%% Plot fertilization year
fig,ax=p1ha.Plot_FECA_Year(meta,roi)

#%% Plot PFI stemwood carbon (20 m)
fig,ax=p1ha.Plot_PFI(meta,roi)

#%% Plot GLOB Biomass
fig,ax=p1ha.Plot_biomass_glob(meta,roi)

#%% Plot soil organic carbon
fig,ax=p1ha.Plot_SoilOrganicCarbon_GSOC(meta,roi)

#%% Plot Global Forest Change Loss Year
fig,ax=p1ha.Plot_GFC_LossYear(meta,roi)

#%% Plot Burn Severity
fig,ax=p1ha.Plot_BurnSeverity(meta,roi)

#%% Plot soil water content
fig,ax=p1ha.Plot_SoilWaterContent(meta,roi)

#%% Plot ground plot sample grid
p1ha.Plot_GroundPlots(meta,roi)

#%% Plot ground plots 2 (panels for paper)
p1ha.Plot_GroundPlotsPanels(meta,roi)

#%% Salvage mask from timber cruise
fig,ax=p1ha.Plot_SalvageLogging(meta,roi)

#%% Ownership
fig,ax=p1ha.Plot_Ownership(meta,roi)

#%% Plot forest range
fig,ax=p1ha.Plot_RangeTenure(meta,roi)

#%% Plot harvest retention (reserves) compilation 1
fig,ax=p1ha.Plot_HarvestRetentionComp1(meta,roi)

#%% Plot distance from road
fig,ax=p1ha.Plot_DistanceFromRoad(meta,roi)

#%% Obsolete
#fig,ax=p1ha.Plot_REARs(meta,roi)

#%% Geomorphons
fig,ax=p1ha.Plot_Geomorphons(meta,roi)

#%% Add vector layers

roi['gdf']['road'].plot(ax=ax[0],facecolor='None',edgecolor=[0,0,0],linewidth=0.25,label='Opening',alpha=1)
roi['gdf']['rail'].plot(ax=ax[0],facecolor='None',edgecolor=[0.25,0.25,0.25],linewidth=0.25,label='Opening',alpha=1)
roi['gdf']['hydrol'].plot(ax=ax[0],facecolor='None',edgecolor=[0.5,0.5,0.5],linewidth=0.25,label='Opening',alpha=1)
roi['gdf']['ogp'].plot(ax=ax[0],facecolor='None',edgecolor=[0.75,0.75,0.75],linewidth=0.25,label='Opening',alpha=1)

roi['gdf']['tpf'].plot(ax=ax[0],marker='^',markersize=4,facecolor=[1,1,0],edgecolor=[0.85,0,0],linewidth=0.25,label='Opening',alpha=1)
roi['gdf']['popp'].plot(ax=ax[0],marker='s',markersize=3,facecolor=[0.75,0.75,0.75],edgecolor=[0,0,0],linewidth=0.25,label='Opening',alpha=1)

gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_ownership','png',900)

# gdf=gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='TA_PARK_ECORES_PA_SVW')
# gdf.plot(ax=ax[0],facecolor='None',edgecolor=[0,0,1],linewidth=1,label='Opening',alpha=1)

# gdf=gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='TA_PROTECTED_LANDS_SV')
# gdf.plot(ax=ax[0],facecolor='None',edgecolor=[0,1,0],linewidth=1,label='Opening',alpha=1)

# gdf=gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='WCL_CONSERVATION_AREAS_NGO_SP')
# gdf.plot(ax=ax[0],facecolor='None',edgecolor=[1,0.7,0],linewidth=1,label='Opening',alpha=1)

# gdf=gpd.read_file(meta['Paths']['GDB']['GDB'] + '\\FTEN_ROAD_SEGMENT_LINES_SVW.geojson')
# gdf.plot(ax=ax[0],facecolor='None',edgecolor=[1,0.7,0],linewidth=1,label='Opening',alpha=1)

# gdf.plot(facecolor='None',edgecolor=[1,0,0],linewidth=0.5,label='Opening',alpha=1)

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

#%%
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

# #roi['gdf']['road'].plot(ax=ax[0],edgecolor='y',linewidth=2,label='Road',alpha=1,zorder=1)
# roi['gdf']['road'].plot(ax=ax[0],edgecolor=[0.7,0.9,1],linewidth=1,label='Road',alpha=1,zorder=1)

# # Plot lumber mills, pulp mills and chipper mills
# mtypeL=['LBR','PLP','PLT']
# ax=u1ha.PlotMills(ax,gdf['tpf']['gdf'],mtypeL,labels='On')

# gdf['cities'].plot(ax=ax[0],marker='s',edgecolor=[0.75,0.3,0],facecolor=[1,0.6,0],lw=1,markersize=20,alpha=1,zorder=2)
# for x,y,label in zip(gdf['cities'].geometry.x,gdf['cities'].geometry.y,gdf['cities'].Name):
#     ax[0].annotate(label,xy=(x,y),xytext=(4,3),textcoords="offset points",color=[0.75,0.3,0])



#%%

# def Plot_RegenType(meta,roi,gdf):

#     lab=['harvested and planted','Harvested (pre-1988), not planted','Harvested (1988-2018), not planted','Harvested (post-2018), not planted','Straight planting','No harvesting / no planting']
#     z1=7*np.ones(roi['grd']['Data'].shape,dtype='int8')
#     for i in range(1,7):
#         z1[(roi['grd']['regentype']['Data']==i)]=i
#     ind=np.where( (roi['grd']['Data']==0) ); z1[ind]=7

#     N_vis=len(lab)
#     N_hidden=1
#     N_tot=N_vis+N_hidden
#     cm=np.vstack( ((0.85,0.7,1,1),(1,1,0,1),(1,0.25,0,1),(0.7,0.7,0.7,1),(0.7,1,0.2,1),(0.9,0.9,0.9,1),(1,1,1,1)) )
#     cm=matplotlib.colors.ListedColormap(cm)

#     plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
#     im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
#     if meta['Graphics']['Map']['Show Bound Within']=='On':
#         roi['gdf']['bound within'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
#     if meta['Graphics']['Map']['Show Bound Land Mask']=='On':
#         roi['gdf']['bc_bound'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
#     if meta['Graphics']['Map']['Show Lakes']=='On':
#         roi['gdf']['lakes'].plot(ax=ax[0],facecolor=[0.82,0.88,1],label='Lakes',linewidth=2.25)
#     if meta['Graphics']['Map']['Show Rivers']=='On':
#         roi['gdf']['rivers'].plot(ax=ax[0],color=[0.6,0.8,1],label='Rivers',linewidth=0.25)
#     if meta['Graphics']['Map']['Show Roads']=='On':
#         roi['gdf']['road'].plot(ax=ax[0],edgecolor='k',linewidth=0.25,label='Road',alpha=1,zorder=1)
#     #gdf['cities'].plot(ax=ax[0],marker='s',edgecolor=[0.75,0.3,0],facecolor=[1,0.6,0],lw=1,markersize=20,alpha=1,zorder=2)
#     #for x,y,label in zip(gdf['cities'].geometry.x,gdf['cities'].geometry.y,gdf['cities'].Name):
#     #    ax[0].annotate(label,xy=(x,y),xytext=(4,3),textcoords="offset points",color=[0.75,0.3,0],fontsize=4)
#     ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
#     ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

#     zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
#     cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
#     cb=plt.colorbar(im,cax=ax[1],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
#     ax[1].set(position=[0.71,0.6,0.05,0.14])
#     cb.ax.set(yticklabels=lab)
#     cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
#     cb.outline.set_edgecolor('w')
#     for i in range(cb_bnd.size):
#         ax[1].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)
#     pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
#     ax[1].set(position=pos2)
#     gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_regentype','png',900)

#     return fig,ax

# Plot_RegenType(meta,roi,gdf)