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
meta['Graphics']['Map']['Show Rivers']='Off'
meta['Graphics']['Map']['Show Rail']='Off'
meta['Graphics']['Map']['Show Roads']='Off'
meta['Graphics']['Map']['Show Cities']='Off'
meta['Graphics']['Map']['Show TPFs']='Off'
meta['Graphics']['Map']['Show Symbol Labels']='Off'

meta['Graphics']['Plot Style']='Manuscript'
meta['Graphics']['gp']=gu.SetGraphics(meta['Graphics']['Plot Style'])
meta['Graphics']['Print Figures']='On'

# Define region of interest
roi={}
#roi['Type']='Prov'; roi['Name']='Prov'

#roi['Type']='LICS'

#roi['Type']='ByRegDis'; roi['Name']='CAPITAL'
#roi['Type']='ByRegDis'; roi['Name']='STRATHCONA'
#roi['Type']='ByRegDis'; roi['Name']='COMOX VALLEY'
#roi['Type']='ByRegDis'; roi['Name']='COWICHAN VALLEY'

#roi['Type']='ByWatershed'; roi['Name']='Chilko-Taseko River'
#roi['Type']='ByWatershed'; roi['Name']='Big Creek'
#roi['Type']='ByWatershed'; roi['Name']='West Arm'
#roi['Type']='ByWatershed'; roi['Name']='Cowichan River'
#roi['Type']='ByWatershed'; roi['Name']='Campbell River'
#roi['Type']='ByWatershed'; roi['Name']='Elephant Hill Fire'
#roi['Type']='ByWatershed'; roi['Name']='Spius Creek'

roi['Type']='ByTSA'; roi['Name']='South'
#roi['Type']='ByTSA'; roi['Name']='Arrowsmith TSA'
#roi['Type']='ByTSA'; roi['Name']='Fort St John TSA'
#roi['Type']='ByTSA'; roi['Name']='Dawson Creek TSA'
#roi['Type']='ByTSA'; roi['Name']='Boundary TSA'
#roi['Type']='ByTSA'; roi['Name']='Cassiar TSA'
#roi['Type']='ByTSA'; roi['Name']='Fort Nelson TSA'
#roi['Type']='ByTSA'; roi['Name']='Kamloops TSA'
#roi['Type']='ByTSA'; roi['Name']='Kalum TSA'
#roi['Type']='ByTSA'; roi['Name']='Kootenay Lake TSA'
#roi['Type']='ByTSA'; roi['Name']='Merritt TSA'
#roi['Type']='ByTSA'; roi['Name']='North Island TSA'
#roi['Type']='ByTSA'; roi['Name']='Prince George TSA'
#roi['Type']='ByTSA'; roi['Name']='Quesnel TSA'
#roi['Type']='ByTSA'; roi['Name']='Williams Lake TSA'
#roi['Type']='ByTSA'; roi['Name']='100 Mile House TSA'

#roi['Type']='LICS'
if roi['Type']=='ByTSA':
	meta['Graphics']['Print Figure Path']=r'C:\Users\rhember\OneDrive - Government of BC\Figures' + '\\TSAs\\' + roi['Name']
elif roi['Type']=='ByRegDis':
	meta['Graphics']['Print Figure Path']=r'C:\Users\rhember\OneDrive - Government of BC\Figures' + '\\Regional Districts\\' + roi['Name']
elif roi['Type']=='Prov':
	meta['Graphics']['Print Figure Path']=r'C:\Users\rhember\OneDrive - Government of BC\Figures' + '\\Province\\' + roi['Name']
elif roi['Type']=='ByWatershed':
	meta['Graphics']['Print Figure Path']=r'C:\Users\rhember\OneDrive - Government of BC\Figures\\Watershed\\' + roi['Name']

t0=time.time()
if roi['Type']=='ByTSA':
	# Search: list(gdf['tsa']['gdf']['Name'].unique())
	roi['List']=roi['Name']
	#roi['List']=['Arrowsmith TSA']
	#roi['List']=['Boundary TSA']
	#roi['List']=['Cassiar TSA']
	#roi['List']=['Dawson Creek TSA']
	#roi['List']=['Fort Nelson TSA']
	#roi['List']=['Fort St. John TSA']
	#roi['List']=['Kalum TSA']
	#roi['List']=['Kamloops TSA']
	#roi['List']=['Kamloops TSA','100 Mile House TSA','Williams Lake TSA']
	roi['List']=['Kamloops TSA','Lillooet TSA','Merritt TSA','Okanagan TSA','Boundary TSA','Fraser TSA','Soo TSA','Arrow TSA','100 Mile House TSA'] # South
	#roi['List']=['Kootenay Lake TSA']
	#roi['List']=['Merritt TSA']
	#roi['List']=['North Island TSA']
	#roi['List']=['Okanagan TSA']
	#roi['List']=['Prince George TSA']
	#roi['List']=['Quesnel TSA']
	#roi['List']=['Williams Lake TSA']
	#roi['List']=['100 Mile House TSA']
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
elif roi['Type']=='ByWatershed':
	#roi['List']=[33206] # Big Creek,33159
	#roi['List']=[12227] # Cowichan River
	#roi['List']=[39532] # Campbell River
	#roi['List']=[13741,23935] # Chilko River / Taseko River
	#roi['List']=[2366,8354,15037] # West Arm
	roi['List']=[23016] # Spius Creek
	#roi['List']=[34069,24249,16771] # Elephant Hill Fire
elif roi['Type']=='LICS':
	#roi['Name']='Hanceville Fire'; roi['Centre']=[-122.92,51.92]; roi['Radius']=40*1000 # Hanceville fire
	#meta['Graphics']['Print Figure Path']=r'C:\Users\rhember\OneDrive - Government of BC\Figures\Reforestation\Reforestation Hanceville ROI'
	#roi['Name']='Plateau Fire'; roi['Centre']=[-123.6,52.75]; roi['Radius']=65*1000 # Plateau fire
	#meta['Graphics']['Print Figure Path']=r'C:\Users\rhember\OneDrive - Government of BC\Figures\LICS\Plateau Fire'
	roi['Name']='Elephant Hill Fire'; roi['Centre']=[-121.15,51.12]; roi['Radius']=42*1000 # Elephant Hill fire
	meta['Graphics']['Print Figure Path']=r'C:\Users\rhember\OneDrive - Government of BC\Figures\LICS\Elephant Hill Fire'
elif roi['Type']=='Prov':
	pass

# Import base maps
gdf=u1ha.Import_GDBs_ProvinceWide(meta)

# Prepare region of interest
roi=u1ha.DefineROI(meta,roi,gdf)
t1=time.time()
print((t1-t0)/60)

# Import rasters over ROI
vList=u1ha.GetRasterListFromSpreadsheet(r'C:\Data\BC1ha\RasterInclusion.xlsx')
roi=u1ha.Import_Raster(meta,roi,vList)

#%% Plot everything
vList=['access']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_AccessZones(meta,roi,vList[0]); del roi['grd'][vList[0]]
vList=['age_ntems']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_Age(meta,roi,vList[0]); del roi['grd'][vList[0]]
vList=['aset']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_ASET(meta,roi,vList[0]); del roi['grd'][vList[0]]
vList=['aset_post18']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_ASET(meta,roi,vList[0]); del roi['grd'][vList[0]]
vList=['age_vri23']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_Age(meta,roi,vList[0]); del roi['grd'][vList[0]]
vList=['bdfrac']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_BroadleafDeciduousFrac(meta,roi,vList[0]); del roi['grd'][vList[0]]
vList=['bdfrac_2049s4']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_BroadleafDeciduousFrac(meta,roi,vList[0]); del roi['grd'][vList[0]]
vList=['bgcz']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_BGC_Zone(meta,roi,vList[0]); del roi['grd'][vList[0]]
vList=['biomass_glob']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_biomass_glob(meta,roi,vList[0]); del roi['grd'][vList[0]]
vList=['bsr_sc']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_BurnSeverity(meta,roi,vList[0]); del roi['grd'][vList[0]]
#vList=['crownc']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_CrownCover(meta,roi,vList[0]); del roi['grd'][vList[0]]
vList=['d2road']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_DistanceFrom(meta,roi,vList[0]); del roi['grd'][vList[0]]
vList=['d2fac']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_DistanceFrom(meta,roi,vList[0]); del roi['grd'][vList[0]]
vList=['elev']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_Elev(meta,roi,vList[0]); del roi['grd'][vList[0]]
vList=['elev']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_Infastructure1(meta,roi,vList[0]); del roi['grd'][vList[0]]
vList=['feca_yr']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_FECA_Year(meta,roi,vList[0]); del roi['grd'][vList[0]]
vList=['fire_yr']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_WildfireYear(meta,roi,vList[0]); del roi['grd'][vList[0]]
vList=['fire_2023']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_WildfireYear(meta,roi,vList[0]); del roi['grd'][vList[0]]
vList=['gfcly']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_GFC_LossYear(meta,roi,vList[0]); del roi['grd'][vList[0]]
vList=['geomorph']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_Geomorphons(meta,roi,vList[0]); del roi['grd'][vList[0]]
vList=['gsoc']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_SoilOrganicCarbon_GSOC(meta,roi,vList[0]); del roi['grd'][vList[0]]
vList=['harv_yr_comp1']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_HarvestYear(meta,roi,vList[0]); del roi['grd'][vList[0]]
vList=['harv_yr_comp2']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_HarvestYear(meta,roi,vList[0]); del roi['grd'][vList[0]]
#vList=['harv_salv']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_SalvageLogging(meta,roi,vList[0]); del roi['grd'][vList[0]]
vList=['harvret1']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_HarvestRetentionComp1(meta,roi,vList[0]); del roi['grd'][vList[0]]
vList=['lc_comp1_1800']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_LandCoverComp1(meta,roi,vList[0]); del roi['grd'][vList[0]]
vList=['lc_comp1_2019']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_LandCoverComp1(meta,roi,vList[0]); 
vList=['lc_comp1_2049s1']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_LandCoverComp1(meta,roi,vList[0]); del roi['grd'][vList[0]]
vList=['lc_comp1_2049s2']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_LandCoverComp1(meta,roi,vList[0]); del roi['grd'][vList[0]]
#vList=['lc_comp1_2049s3']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_LandCoverComp1(meta,roi,vList[0]); del roi['grd'][vList[0]]
#vList=['lc_ntems_2019_recl']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_LandCoverComp1(meta,roi,vList[0]); del roi['grd'][vList[0]]
vList=['lc_vri_recl']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_LandCoverComp1(meta,roi,vList[0]); del roi['grd'][vList[0]]
vList=['lc_cec_2020']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_LC20_CEC(meta,roi,vList[0]); del roi['grd'][vList[0]]
vList=['lu_comp1_2019']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_LandUseComp1(meta,roi,vList[0]); del roi['grd'][vList[0]]
vList=['lu_comp1_2049s1']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_LandUseComp1(meta,roi,vList[0]); del roi['grd'][vList[0]]
vList=['lu_comp1_2049s2']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_LandUseComp1(meta,roi,vList[0]); del roi['grd'][vList[0]]
vList=['lu_comp1_2049s3']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_LandUseComp1(meta,roi,vList[0]); del roi['grd'][vList[0]]
vList=['lu_comp1_2049s4']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_LandUseComp1(meta,roi,vList[0]); del roi['grd'][vList[0]]
vList=['luc1_hist_type']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_LandCoverLandUseChange(meta,roi,vList[0]); del roi['grd'][vList[0]]
vList=['luc1_1019_type']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_LandCoverLandUseChange(meta,roi,vList[0]); del roi['grd'][vList[0]]
vList=['luc1_fut_s1_type']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_LandCoverLandUseChange(meta,roi,vList[0]); del roi['grd'][vList[0]]
vList=['luc1_hist_yr']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_LUC_Year(meta,roi,vList[0]); del roi['grd'][vList[0]]
vList=['own']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_Ownership(meta,roi,vList[0]); del roi['grd'][vList[0]]
vList=['plam']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_PlantedMask(meta,roi,vList[0]); del roi['grd'][vList[0]]
vList=['prcp_ann_n']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_MAP(meta,roi,vList[0]); del roi['grd'][vList[0]]
vList=['rangecon']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_RangeTenure(meta,roi,vList[0]); del roi['grd'][vList[0]]
vList=['rears']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_REARs(meta,roi,vList[0]); del roi['grd'][vList[0]]
vList=['spc1_ntems']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_Spc1_NTEMS(meta,roi,vList[0]); del roi['grd'][vList[0]]
vList=['tdc_wsg']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_TreeDensityClass(meta,roi,vList[0]); del roi['grd'][vList[0]]
vList=['tmean_ann_n']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_MAT(meta,roi,vList[0]); del roi['grd'][vList[0]]
vList=['upwetf_ntems']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_UplandWetlandForest(meta,roi,vList[0]); del roi['grd'][vList[0]]
vList=['upwetf_vri']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_UplandWetlandForest(meta,roi,vList[0]); del roi['grd'][vList[0]]
#vList=['si_vri23']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_SI(meta,roi); del roi['grd'][vList[0]]
#vList=['sphl_vri23']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_SPH(meta,roi,vList[0]); del roi['grd'][vList[0]]
#vList=['sphd_vri23']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_SPH(meta,roi,vList[0]); del roi['grd'][vList[0]]
vList=['elev']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_WaterManagement(meta,roi,vList[0]); del roi['grd'][vList[0]]
vList=['ws_mjjas_n']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_SoilWaterContent(meta,roi,vList[0]); del roi['grd'][vList[0]]
#vList=['']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_PFI(meta,roi,vList[0]); del roi['grd'][vList[0]]
#vList=['']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.(meta,roi,vList[0]); del roi['grd'][vList[0]]

#%% Plot forest tenure roads
def PlotRoadsForest(meta):
	# Import data
	vList=['elev']; roi=u1ha.Import_Raster(meta,roi,vList)
	roi=u1ha.Import_GDB_Over_ROI(meta,roi,['road_ften'])
	# Plot
	meta['Graphics']['Map']['Show Roads']='On'
	fig,ax=p1ha.Plot_Elev(meta,roi,vList[0]);
	roi['gdf']['road_ften']['gdf'].plot(ax=ax[0],edgecolor=[1,1,0],facecolor='none',linewidth=0.25)
	gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_road_ften','png',900)
return

#%% Mills
vList=['elev']; roi=u1ha.Import_Raster(meta,roi,vList);
fig,ax=p1ha.Plot_InfastructureLumber(meta,roi,vList[0]);
fig,ax=p1ha.Plot_InfastructurePulp(meta,roi,vList[0]);
fig,ax=p1ha.Plot_InfastructurePanel(meta,roi,vList[0]);
fig,ax=p1ha.Plot_InfastructurePellet(meta,roi,vList[0]);
fig,ax=p1ha.Plot_InfastructureChipper(meta,roi,vList[0]);

# # Plot lumber mills, pulp mills and chipper mills
# mtypeL=['LBR','PLP','PLT']
# ax=u1ha.PlotMills(ax,gdf['tpf']['gdf'],mtypeL,labels='On')

# #%% Plot major timber producing facilities
# def PlotMills(ax,gdf,mtypeL,labels):
# 	lw=0.75
# 	for mtype in mtypeL:
# 		#mtype='PLP'
# 		ind=np.where( (gdf['PRODUCT_CODE']==mtype) )[0]
# 		if mtype=='LBR':
# 			y=gdf.iloc[ind]['EST_AN_CAP_MLN_BOARD_FT']/453
# 			ms=500*y
# 			gdf.iloc[ind].plot(ax=ax,marker='o',edgecolor='g',facecolor='g',lw=lw,markersize=ms,alpha=0.35,zorder=2)
# 		elif (mtype=='PLY') | (mtype=='VNR') | (mtype=='OSB') | (mtype=='PNL'):
# 			# One PNL (panel) mill WestPine MDF in Quesnell - it is MDF
# 			y=gdf.iloc[ind]['EST_AN_CAP_MLN_SQ_FT']/885
# 			ms=300*y
# 			gdf.iloc[ind].plot(ax=ax,marker='o',edgecolor='c',facecolor='c',lw=lw,markersize=ms,alpha=0.35,zorder=2)
# 		elif mtype=='PLT':
# 			y=gdf.iloc[ind]['EST_AN_CAP_000_TONNES']*2
# 			ms=1.0*y
# 			gdf.iloc[ind].plot(ax=ax,marker='o',edgecolor='r',facecolor='r',lw=lw,markersize=ms,alpha=0.35,zorder=2)
# 		elif mtype=='CHP':
# 			y=gdf.iloc[ind]['EST_AN_CAP_000_BDUS']*2
# 			ms=0.75*y
# 			gdf.iloc[ind].plot(ax=ax,marker='o',edgecolor='y',facecolor='y',lw=lw,markersize=ms,alpha=0.35,zorder=2)
# 		elif mtype=='PLP':
# 			y=gdf.iloc[ind]['EST_AN_CAP_000_TONNES']*2
# 			ms=0.75*y
# 			gdf.iloc[ind].plot(ax=ax,marker='o',edgecolor='k',facecolor='k',lw=lw,markersize=ms,alpha=0.35,zorder=2)
# 		elif mtype=='LVL':
# 			# Laminated veneer lumber
# 			y=gdf.iloc[ind]['EST_AN_CAP_MLN_CUBIC_FT']#/885
# 			ms=200*y
# 			gdf.iloc[ind].plot(ax=ax,marker='o',edgecolor='k',facecolor='g',lw=lw,markersize=ms,alpha=0.35,zorder=2)

# 		if labels=='On':
# 			for x,y,label in zip(gdf.iloc[ind].geometry.x,gdf.iloc[ind].geometry.y,gdf.iloc[ind]['COMPANY_NAME']):
# 				ax.annotate(label,xy=(x,y),xytext=(5,5),textcoords="offset points")
# 	return ax

#%% Plot vector fire
def PlotWildfireBoundaries(meta):
	roi=u1ha.Import_GDB_Over_ROI(meta,roi,['wf'])
	df=roi['gdf']['wf']['gdf']
	Year=2017
	ikp=np.where(df['FIRE_YEAR']==Year)[0]
	vList=['elev']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_Elev(meta,roi,vList[0])
	df.iloc[ikp].plot(ax=ax[0],facecolor='none',edgecolor=[0.6,0,0],linewidth=1)
	df['coords']=df['geometry'].apply(lambda x: x.representative_point().coords[:])
	df['coords']=[coords[0] for coords in df['coords']]
	for idx, row in df.iloc[ikp].iterrows():
		ax[0].annotate(row['FIRE_NUMBER'],xy=row['coords'],horizontalalignment='center',color=[0,0,0],fontsize=12)
	return

#%% Plot watersheds
def PlotWatershedBoundaries(meta):
	roi=u1ha.Import_GDB_Over_ROI(meta,roi,['wshed'],OrderList=[6])
	vList=['elev']; roi=u1ha.Import_Raster(meta,roi,vList);
	fig,ax=p1ha.Plot_Elev(meta,roi,vList[0])

	roi['gdf']['wshed']['gdf'].plot(ax=ax[0],edgecolor=[0,0,1],facecolor='none',linewidth=1)
	roi['gdf']['rivers'].plot(ax=ax[0],color=[0.4,0.8,1],label='Rivers',linewidth=1)

	roi['gdf']['wshed']['gdf']['coords']=roi['gdf']['wshed']['gdf']['geometry'].apply(lambda x: x.representative_point().coords[:])
	roi['gdf']['wshed']['gdf']['coords']=[coords[0] for coords in roi['gdf']['wshed']['gdf']['coords']]
	for idx, row in roi['gdf']['wshed']['gdf'].iterrows():
		ax[0].annotate(row['GNIS_NAME'],xy=row['coords'],horizontalalignment='center',verticalalignment='center',color=[0,0,1])
		#ax[0].annotate(row['GNIS_ID'],xy=row['coords'],horizontalalignment='center',verticalalignment='center',color=[0,0,1])

	# Select watersheds
	ind=np.where(roi['gdf']['wshed']['gdf']['GNIS_NAME']=='Spius Creek')[0]
	ind=np.where(roi['gdf']['wshed']['gdf']['GNIS_NAME']=='Salmon River')[0]
	roi['gdf']['wshed']['gdf'].iloc[ind].plot(ax=ax[0],edgecolor=[0,0,1],facecolor='none',linewidth=1)
	roi['gdf']['rivers'].plot(ax=ax[0],color=[0.4,0.8,1],label='Rivers',linewidth=1)
	return

#%% Add watershed
def WaterManagementWithSelectWatersheds():
	vList=['elev']; roi=u1ha.Import_Raster(meta,roi,vList);
	fig,ax=p1ha.Plot_WaterManagement(meta,roi,vList[0]);

	# Select watersheds
	#roi=u1ha.Import_GDB_Over_ROI(meta,roi,['wshed'],OrderList=[6])
	nam='Spius Creek'
	ind=np.where(roi['gdf']['wshed']['gdf']['GNIS_NAME']==nam)[0]
	#ind=np.where(roi['gdf']['wshed']['gdf']['GNIS_NAME']=='Salmon River')[0]
	roi['gdf']['wshed']['gdf'].iloc[ind].plot(ax=ax[0],edgecolor=[0.2,0.4,0.6],facecolor=[0.6,0.8,1],alpha=0.5,linewidth=1)
	gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_watermanagement_' + nam,'png',900)
	return

#%% Plot Streamflow sites
def PlotStreamflowSites(meta):
	d=gu.ReadExcel(r'C:\Data\Streamflow\GRDC_Stations.xlsx')
	srs=gis.ImportSRSs()
	x,y=srs['Proj']['BC1ha'](d['long'],d['lat'])
	
	#vList=['elev']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_Elev(meta,roi,vList[0]);
	ax[0].plot(x,y,'co')
	for i in range(x.size):
		ax[0].annotate(d['grdc_no'][i],xy=[x[i],y[i]],horizontalalignment='center',color=[0,0,0])
	return

#%%
vList=['elev']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_Elev(meta,roi,vList[0]);
ax[0].contour(np.flip(roi['grd']['elev']['Data'],axis=0),extent=roi['grd']['Extent'],cmap="viridis",levels=[400,600,800,1000,1200,1400],linewidth=0.25)

#%% Streamflow simulations from PCIC
a=pd.read_excel(r'C:\Data\Streamflow\PCIC\SPIUS.csv.ascii.xlsx')
y=a[' CanESM2_rcp45_r1i1p1'].values
tv=gu.tvec('d',1945,2099)

yH=np.zeros(12)
yF=np.zeros(12)
for mo in range(12):
	ind=np.where( (tv[:,0]>=1971) & (tv[:,0]<=2000) & (tv[:,1]==mo+1) )[0]
	yH[mo]=np.sum(y[ind])/30
	ind=np.where( (tv[:,0]>=2071) & (tv[:,0]<=2099) & (tv[:,1]==mo+1) )[0]
	yF[mo]=np.sum(y[ind])/30
fig,ax=plt.subplots(1,figsize=gu.cm2inch(10,10))
plt.plot(yH)
plt.plot(yF,'r--')


#%% Streamflow
# m3/s
# Time is provided in days since 1700-01-01
nc=gu.ReadNC(r'C:\Data\Streamflow\GRDC-Monthly.nc')
tvd=gu.tvec('d',1700,2025)
tvd[nc['time'][0],:]
#indS=np.where(nc['id']==4207300)[0]
indS=np.where(nc['id']==4207305)[0]
Area=nc['area'][indS]
sf=nc['runoff_mean'][:,indS[0]]
sf[sf<0]=np.nan
mu=gu.IvlMean(sf,12)
#plt.close('all'); plt.bar(np.arange(12),mu)

# Gap-fill
mut=np.tile(mu,(int(sf.size/12)))
sf_gf=sf.copy()
sf_gf[np.isnan(sf)]=mut[np.isnan(sf)]

sfA=gu.BlockMean(sf,12)
sfgfA=gu.BlockMean(sf_gf,12)
MissingA=gu.BlockMissing(sf,12)
sfgfA[MissingA>3]=np.nan

tStart=1911; tv=np.arange(tStart,tStart+sfA.size,1)
#plt.bar(tv,Missing)
plt.close('all');
plt.plot(tv,sfgfA,'-rs')
plt.plot(tv,sfA,'-bo')


#%% 3D terrrain
# vList=['elev']; roi=u1ha.Import_Raster(meta,roi,vList);
ivl=1
cm0=np.flip(np.vstack( ((0.32,0.19,0.19,1),(0.41,0.28,0.27,1),(0.49,0.38,0.36,1),(0.58,0.47,0.44,1),(0.66,0.57,0.53,1),(0.75,0.66,0.61,1),(0.83,0.76,0.7,1),(0.92,0.85,0.78,1),(1,0.95,0.87,1),(0.83,0.87,0.85,1),(0.67,0.78,0.82,1),(0.5,0.7,0.8,1),(0.33,0.61,0.78,1),(0.17,0.53,0.76,1),(0,0.45,0.74,1)) ),axis=0)
cm0=matplotlib.colors.LinearSegmentedColormap.from_list('twi',cm0,N=30)
#cmap=cm0(twi['Data']/np.amax(twi['Data']))
z=roi['grd']['Data'][0::ivl,0::ivl]*roi['grd']['elev']['Data'][0::ivl,0::ivl]
zmin=300#np.amin(roi['grd']['elev']['Data'][0::ivl,0::ivl])
zmax=np.amax(z)
cmap=cm0( (z-zmin)/(zmax-zmin) )

e=20; a=360-2*45
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

#%% Plot model results
flg=0
if flg==1:
	pth=r'D:\FCI_Projects\TSA_DawsonCreek'
	md=gu.ipickle(pth + '\\Outputs\MapData_Scn1.pkl')
	geos=gu.ipickle(r'D:\FCI_Projects\TSA_DawsonCreek\geos.pkl')
	zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
	vnam='E_CO2e_AGHGB_WSub'
	#
	fig,ax=p1ha.Plot_FromModel(meta,roi,zRef,geos,md,vnam)

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
#fig,ax=p1ha.Plot_LandCoverComp1(meta,roi,'lc_comp1_2049s4')

#%% Land Use Compilation 1
fig,ax=p1ha.Plot_LandUseComp1(meta,roi,'lu_comp1_2019')
fig,ax=p1ha.Plot_LandUseComp1(meta,roi,'lu_comp1_2049s1')
fig,ax=p1ha.Plot_LandUseComp1(meta,roi,'lu_comp1_2049s2')
fig,ax=p1ha.Plot_LandUseComp1(meta,roi,'lu_comp1_2049s3')
#fig,ax=p1ha.Plot_LandUseComp1(meta,roi,'lu_comp1_2049s4')
p1ha.Plot_LandUseComp1Panels(meta,roi) # Plot panels

#%% Other Land Cover estimates
fig,ax=p1ha.Plot_LandCoverComp1(meta,roi,'lc_ntems_2019_recl')
fig,ax=p1ha.Plot_LandCoverComp1(meta,roi,'lc_vri_recl')
fig,ax=p1ha.Plot_LC20_CEC(meta,roi)

#%% Plot LUC (Afforestation or Deforestation) from CEC
fig,ax=p1ha.Plot_LandCoverLandUseChange(meta,roi,'luc1_hist_type')
fig,ax=p1ha.Plot_LandCoverLandUseChange(meta,roi,'luc1_1019_type')
fig,ax=p1ha.Plot_LandCoverLandUseChange(meta,roi,'luc1_fut_s1_type')
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

#%% Salvage mask from timber cruise
fig,ax=p1ha.Plot_SalvageLogging(meta,roi,'harv_salv')

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


