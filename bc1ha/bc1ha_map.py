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
import fcgadgets.bc1ha.bc1ha_utils as u1ha
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
meta['Graphics']['Map']['Show Lakes']='On'
meta['Graphics']['Map']['Show Rivers']='On'
meta['Graphics']['Map']['Show Rail']='Off'
meta['Graphics']['Map']['Show Roads']='Off'
meta['Graphics']['Map']['Show Cities']='On'
meta['Graphics']['Map']['Show TPFs']='Off'
meta['Graphics']['Map']['Show Symbol Labels']='On'
meta['Graphics']['Map']['Show Inset Map']='On'
meta['Graphics']['Map']['Show Scalebar']='On'
meta['Graphics']['Vector Import']={'Water Management':'On','Wetland':'Off'}

meta['Graphics']['Plot Style']='Manuscript'
meta['Graphics']['gp']=gu.SetGraphics(meta['Graphics']['Plot Style'])
meta['Graphics']['Print Figures']='On'

# Define region of interest
roi={}
#roi['Type']='Prov'
#roi['Type']='FromMask'
#roi['Type']='ByRegDis'
roi['Type']='ByWatershed'
#roi['Type']='ByTSA'
#roi['Type']='LICS'

t0=time.time()
if roi['Type']=='ByTSA':
	# Search: list(gdf['tsa']['gdf']['Name'].unique())
	#roi['List']=roi['Name']

	# North:
	#roi['Name']='North'; roi['List']=['Kalum TSA','Kispiox TSA','Bulkley TSA','Morice TSA','MacKenzie TSA','Prince George TSA','Lakes TSA','Fort St. John TSA','Dawson Creek TSA'] # North
	# Central:
	#roi['Type']='ByTSA'; roi['Name']='Central'; roi['List']=['100 Mile House TSA','Quesnel TSA','Williams Lake TSA','GBR North TSA','GBR South TSA','Mid Coast TSA'] # Central
	# South:
	roi['Name']='South'; roi['List']=['Kamloops TSA','Lillooet TSA','Merritt TSA','Okanagan TSA','Boundary TSA','Fraser TSA','Soo TSA','Arrow TSA','100 Mile House TSA'] # South

	# Individual:
	#roi['Name']='Arrowsmith TSA'; roi['List']=['Arrowsmith TSA']
	#roi['Name']='Boundary TSA'; roi['List']=['Boundary TSA']
	#roi['Name']='Cassiar TSA'; roi['List']=['Cassiar TSA']
	#roi['Name']='Dawson Creek TSA'; roi['List']=['Dawson Creek TSA']
	#roi['List']=['Fort Nelson TSA']
	#roi['Name']='Fort St John TSA'; roi['List']=['Fort St. John TSA']
	#roi['List']=['Kalum TSA']
	#roi['List']=['Kamloops TSA']
	#roi['List']=['Kamloops TSA','100 Mile House TSA','Williams Lake TSA']
	#roi['List']=['Kootenay Lake TSA']
	#roi['Name']='Merritt TSA' roi['List']=['Merritt TSA']
	#roi['Name']='North Island TSA' roi['List']=['North Island TSA']
	#roi['List']=['Okanagan TSA']
	#roi['Name']='Prince George TSA'; roi['List']=['Prince George TSA']
	#roi['Name']='Quesnel TSA'; roi['List']=['Quesnel TSA']
	#roi['Name']='Williams Lake TSA'; roi['List']=['Williams Lake TSA']
	#roi['Name']='100 Mile House TSA'; roi['List']=['100 Mile House TSA']
	#roi['List']=['Merritt TSA','Kamloops TSA','100 Mile House TSA','Okanagan TSA','Williams Lake TSA','Lillooet TSA','Boundary TSA'] # ,'Arrow TSA','Revelstoke TSA'
	#roi['List']=list(gdf['tsa']['key']['Name'])
	# Western spruce budworm study (do not change!)
	#roi['List']=['Merritt TSA','Kamloops TSA','100 Mile House TSA','Okanagan TSA','Williams Lake TSA','Lillooet TSA','Boundary TSA'] # ,'Arrow TSA','Revelstoke TSA'

elif roi['Type']=='ByRegDis':
	roi['List']=roi['Name']
	# Search: gdf['regdis']['gdf']['REGIONAL_DISTRICT_NAME'].unique()
	# roi['Name']='CAPITAL'; roi['List']=['CAPITAL']
	# roi['Name']='STRATHCONA'; roi['List']=['STRATHCONA']
	# roi['Name']='COMOX VALLEY'
	# roi['Name']='COWICHAN VALLEY'; roi['List']=['COWICHAN VALLEY']

elif roi['Type']=='ByWatershed':
	#roi['Name']='Chilko-Taseko River'; roi['List']=[13741,23935] # Chilko River / Taseko River
	#roi['Name']='Big Creek'; roi['List']=[33206] # Big Creek,33159
	#roi['Name']='West Arm'; roi['List']=[2366,8354,15037] # West Arm
	#roi['Name']='Cowichan River'; roi['List']=[12227] # Cowichan River
	#roi['Name']='Campbell River'; roi['List']=[39532] # Campbell River
	#roi['Name']='Elephant Hill Fire'; roi['List']=[34069,24249,16771] # Elephant Hill Fire
	#roi['Name']='Spius Creek'; roi['List']=[23016] # Spius Creek
	#roi['Name']='Adams River'; roi['List']=[39257] # Adams River
	roi['Name']='Nicola River'; roi['List']=[39400]

elif roi['Type']=='FromMask':
	#roi['Name']='CDF'
	#roi['Name']='IDF'
	roi['Name']='Lower Thompson'

elif roi['Type']=='LICS':
	#roi['Name']='Hanceville Fire'; roi['Centre']=[-122.92,51.92]; roi['Radius']=40*1000 # Hanceville fire
	#meta['Graphics']['Print Figure Path']=r'C:\Users\rhember\OneDrive - Government of BC\Figures\Reforestation\Reforestation Hanceville ROI'
	#roi['Name']='Plateau Fire'; roi['Centre']=[-123.6,52.75]; roi['Radius']=65*1000 # Plateau fire
	#meta['Graphics']['Print Figure Path']=r'C:\Users\rhember\OneDrive - Government of BC\Figures\LICS\Plateau Fire'
	roi['Name']='Elephant Hill Fire'; roi['Centre']=[-121.1,51.1]; roi['Radius']=70*1000 # Elephant Hill fire
	meta['Graphics']['Print Figure Path']=r'C:\Users\rhember\OneDrive - Government of BC\Figures\LICS\Elephant Hill Fire'

elif roi['Type']=='Prov':
	roi['Name']='Prov'
else:
	pass

# Import base maps
meta['gdf']=u1ha.Import_GDBs_ProvinceWide(meta)

# Prepare region of interest
meta,roi=u1ha.DefineROI(meta,roi,meta['gdf'])
t1=time.time()
print((t1-t0)/60)

# Import rasters over ROI
vList=u1ha.GetRasterListFromSpreadsheet(r'C:\Data\BC1ha\RasterInclusion.xlsx')
roi=u1ha.Import_Raster(meta,roi,vList)

#%% Plot everything
vList=['access']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_AccessZones(meta,roi,vList[0]); del roi['grd'][vList[0]]
vList=['age_ntems']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_Age(meta,roi,vList[0]); del roi['grd'][vList[0]]
vList=['aset']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_ASET(meta,roi,vList[0]); del roi['grd'][vList[0]]
vList=['aset_post17']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_ASET(meta,roi,vList[0]); del roi['grd'][vList[0]]
vList=['aset_post17','fire_yl']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_ASET_Wildfire(meta,roi,vList[0]); del roi['grd'][vList[0]]
vList=['age_vri23']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_Age(meta,roi,vList[0]); del roi['grd'][vList[0]]
vList=['bdfrac']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_BroadleafDeciduousFrac(meta,roi,vList[0]); del roi['grd'][vList[0]]
vList=['bdfrac_2049s4']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_BroadleafDeciduousFrac(meta,roi,vList[0]); del roi['grd'][vList[0]]
vList=['bgcz']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_BGC_Zone(meta,roi,vList[0]); del roi['grd'][vList[0]]
vList=['biomass_glob']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_biomass_glob(meta,roi,vList[0]); del roi['grd'][vList[0]]
vList=['bsr_sc']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_BurnSeverity(meta,roi,vList[0]); del roi['grd'][vList[0]]
#vList=['crownc']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_CrownCover(meta,roi,vList[0]); del roi['grd'][vList[0]]
vList=['d2road']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_DistanceFrom(meta,roi,vList[0]); del roi['grd'][vList[0]]
vList=['d2fac']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_DistanceFrom(meta,roi,vList[0]); del roi['grd'][vList[0]]
vList=['elev']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_Elev(meta,roi,vList[0])
vList=['elev']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_Infastructure1(meta,roi,vList[0])
vList=['elev']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_Settlements(meta,roi,vList[0])
vList=['elev']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_WaterManagement(meta,roi,vList[0])
vList=['elev']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.PlotWatershedBoundaries(meta,roi,[5,6,7],1e6)
vList=['feca_yr']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_FECA_Year(meta,roi,vList[0]); del roi['grd'][vList[0]]
vList=['fire_yl']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_WildfireYear(meta,roi,vList[0]); del roi['grd'][vList[0]]
vList=['fire_2023']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_WildfireYear(meta,roi,vList[0]); del roi['grd'][vList[0]]
vList=['forfrac5m_Tom23']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_formask(meta,roi,vList[0]); del roi['grd'][vList[0]]
vList=['gfcly']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_GFC_LossYear(meta,roi,vList[0]); del roi['grd'][vList[0]]
vList=['geomorph']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_Geomorphons(meta,roi,vList[0]); del roi['grd'][vList[0]]
vList=['gromo_net50']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_gromo_GN(meta,roi,vList[0]); del roi['grd'][vList[0]]
vList=['gromo_net']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_gromo_GN(meta,roi,vList[0]); del roi['grd'][vList[0]]
vList=['gromo_g']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_gromo_G(meta,roi,vList[0]); del roi['grd'][vList[0]]
vList=['gromo_m']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_gromo_M(meta,roi,vList[0]); del roi['grd'][vList[0]]
vList=['gsoc']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_SoilOrganicCarbon_GSOC(meta,roi,vList[0]); del roi['grd'][vList[0]]
vList=['h_vri23'];roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_Height(meta,roi,vList[0]); del roi['grd'][vList[0]]
vList=['h_Tolan23'];roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_Height(meta,roi,vList[0]); del roi['grd'][vList[0]]
vList=['h_Lang22'];roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_Height(meta,roi,vList[0]); del roi['grd'][vList[0]]
vList=['harv_yr_comp1']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_HarvestYear(meta,roi,vList[0]); del roi['grd'][vList[0]]
vList=['harv_yr_comp2']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_HarvestYear(meta,roi,vList[0]); del roi['grd'][vList[0]]
#vList=['harv_salv']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_SalvageLogging(meta,roi,vList[0]); del roi['grd'][vList[0]]
vList=['reserve_comp1']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_ReserveComp1(meta,roi,vList[0]); del roi['grd'][vList[0]]
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
vList=['peat']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_Peatland(meta,roi,vList[0]); del roi['grd'][vList[0]]
vList=['plam']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_PlantedMask(meta,roi,vList[0]); del roi['grd'][vList[0]]
vList=['pl_yl']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_PlantedYearLast(meta,roi,vList[0]); del roi['grd'][vList[0]]
vList=['prcp_ann_n']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_MAP(meta,roi,vList[0]); del roi['grd'][vList[0]]
vList=['rangecon']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_RangeTenure(meta,roi,vList[0]); del roi['grd'][vList[0]]
vList=['rears']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_REARs(meta,roi,vList[0]); del roi['grd'][vList[0]]
vList=['slope']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_Slope(meta,roi,vList[0]); del roi['grd'][vList[0]]
vList=['spc1_ntems']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_Spc1_NTEMS(meta,roi,vList[0]); del roi['grd'][vList[0]]
vList=['tdc_wsg']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_TreeDensityClass(meta,roi,vList[0]); del roi['grd'][vList[0]]
vList=['tdgt30_vri23']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_TreeDensityGT30(meta,roi,vList[0]); del roi['grd'][vList[0]]
vList=['tmean_ann_n']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_MAT(meta,roi,vList[0]); del roi['grd'][vList[0]]
vList=['upwetf_ntems']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_UplandWetlandForest(meta,roi,vList[0]); del roi['grd'][vList[0]]
vList=['upwetf_vri']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_UplandWetlandForest(meta,roi,vList[0]); del roi['grd'][vList[0]]
#vList=['si_vri23']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_SI(meta,roi); del roi['grd'][vList[0]]
#vList=['sphl_vri23']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_SPH(meta,roi,vList[0]); del roi['grd'][vList[0]]
#vList=['sphd_vri23']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_SPH(meta,roi,vList[0]); del roi['grd'][vList[0]]
vList=['ws_mjjas_n']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_SoilWaterContent(meta,roi,vList[0]); del roi['grd'][vList[0]]
#vList=['']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_PFI(meta,roi,vList[0]); del roi['grd'][vList[0]]
#vList=['']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.(meta,roi,vList[0]); del roi['grd'][vList[0]]

#%% Forest mask for costum vector layers
fig,ax=p1ha.Plot_ForestMask(meta,roi)

#%% Import wildfire perimiters
roi=u1ha.Import_GDB_Over_ROI(meta,roi,['wf']); t0=2017; t1=2100;
vList=['elev']; roi=u1ha.Import_Raster(meta,roi,vList); p1ha.PlotWildfireBoundaries(meta,ax,roi,t0,t1)

#%% Climate panel
vL=['tmean_ann_n','prcp_ann_n','etp_ann_n','runoff_ann_n','melt_ann_n','wsp_ann_n','ws_ann_n','cwd_ann_n']
roi=u1ha.Import_Raster(meta,roi,vL)
p1ha.Plot_ClimateNormalsPanels(meta,roi)

#%% Plot forest tenure roads
p1ha.Plot_RoadsFromForestTenures(meta,roi)

#%% Plot mills
vList=['elev']; roi=u1ha.Import_Raster(meta,roi,vList);
fig,ax=p1ha.Plot_InfastructureLumber(meta,roi,vList[0]);
fig,ax=p1ha.Plot_InfastructurePulp(meta,roi,vList[0]);
fig,ax=p1ha.Plot_InfastructurePanel(meta,roi,vList[0]);
fig,ax=p1ha.Plot_InfastructurePellet(meta,roi,vList[0]);
fig,ax=p1ha.Plot_InfastructureChipper(meta,roi,vList[0]);
fig,ax=p1ha.Plot_InfastructureBioenergy(meta,roi,vList[0]);

#%% Site visits
p1ha.Plot_SiteVisits(meta,roi,vList)

#%% Plot elevation contours
vList=['elev']; roi=u1ha.Import_Raster(meta,roi,vList); fig,ax=p1ha.Plot_Elev(meta,roi,vList[0]);
ax[0].contour(np.flip(roi['grd']['elev']['Data'],axis=0),extent=roi['grd']['Extent'],cmap="viridis",levels=[400,600,800,1000,1200,1400],linewidth=0.25)

#%% 3D terrrain
p1ha.Plot3DTerrain(meta,roi)

#%% Plot model results
flg=0
if flg==1:
	md=gu.ipickle(pth + '\\Outputs\MapData_Scn1.pkl')
	geos=gu.ipickle(pth + '\\geos.pkl')
	zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
	vNam='C_Biomass'; fig,ax=p1ha.Plot_FromModel(meta,roi,zRef,geos,md,vNam)
	vNam='C_Forest'; fig,ax=p1ha.Plot_FromModel(meta,roi,zRef,geos,md,vNam)
	vNam='C_G_Gross'; fig,ax=p1ha.Plot_FromModel(meta,roi,zRef,geos,md,vNam)
	vNam='E_Domestic_ForestSector_NEE'; fig,ax=p1ha.Plot_FromModel(meta,roi,zRef,geos,md,vNam)
	#vNam='C_ToMillTotal'; fig,ax=p1ha.Plot_FromModel(meta,roi,zRef,geos,md,vNam)
	vNam='E_NAB'; fig,ax=p1ha.Plot_FromModel(meta,roi,zRef,geos,md,vNam)

#%% Import required vector geodatabases
#vList=['cc','fcres']
#vList=['wf'] # 'op','cc','fcres','ogsr'
#roi=u1ha.Import_GDB_Over_ROI(meta,roi,vList)

#%%
plt.close('all')
plt.matshow(roi['grd']['elev']['Data'],clim=[1090,3000])

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
