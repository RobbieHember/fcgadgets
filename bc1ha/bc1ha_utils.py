#%% Import modules
import os
import numpy as np
import gc
import matplotlib.pyplot as plt
import matplotlib.colors
import geopandas as gpd
import pandas as pd
from scipy.interpolate import griddata
import copy
from shapely.geometry import Polygon,Point,box,shape
import pyproj
import rasterio
from rasterio import features
from scipy.interpolate import NearestNDInterpolator
import fiona
import time
import cv2
import gzip
from datetime import datetime
import netCDF4 as nc
from scipy.interpolate import griddata
import scipy.ndimage
from perlin_noise import PerlinNoise
from numpy.random import default_rng
import fcgadgets.macgyver.util_general as gu
import fcgadgets.macgyver.util_gis as gis
import fcgadgets.macgyver.util_query_gdb as qgdb
import fcgadgets.bc1ha.bc1ha_utils as u1ha
import fcgadgets.gaia.gaia_util as gaia
import fcexplore.field_plots.Processing.fp_util as ufp
import fcgadgets.cbrunner.cbrun_util as cbu

#%%
# srs=gis.ImportSRSs()
# pts=[Point(srs['Proj']['BC1ha'](-123.370254,48.417787))]
# gdf=gpd.GeoDataFrame({'geometry':pts,'Name':'545 Superior Street'})
# #gdf_xy.crs=gdf_bc_boundary.crs
# gdf.to_file(r'C:\Users\rhember\Government of BC\FCCS Team - General\temp\offi_alb.shp')

# import os, fiona
# fiona.supported_drivers['KML'] = 'rw' #Enable kml driver
# gdf.to_file(r'C:\Users\rhember\Government of BC\FCCS Team - General\temp\offi_alb.kml',driver='KML')

#%% Create ID for categorical variable
def CreateIdForCategoricalVariable(meta,lNam,vNam,df):
	v=list(df.keys())[0]
	df['ID_' + vNam]=np.zeros(df[v].size)
	for k in meta['LUT'][lNam][vNam].keys():
		ind=np.where(df[vNam]==k)[0]
		df['ID_' + vNam][ind]=meta['LUT'][lNam][vNam][k]
	return df

#%% Initialize project
def Init(*argv):
	meta={}
	# Set paths
	meta['Paths']={}
	meta['Paths']['bc1ha']=r'C:\Data\BC1ha'
	meta['Paths']['bc1ha Ref Grid']=meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandMask.tif'
	meta['Paths']['bc1ha Ref Grid Buf']=meta['Paths']['bc1ha'] + '\\LandCoverUse\LandMaskBuf.tif'
	meta['Paths']['bc1ha Ref Grid Buf Rev']=meta['Paths']['bc1ha'] + '\\LandCoverUse\LandMaskBufRev.tif'

	meta['Paths']['bc5k']=r'C:\Data\BC5k'
	meta['Paths']['bc5k Ref Grid']=r'C:\Data\BC5k' + '\\LandMask.tif'
	meta['Paths']['GDB']={}
	meta['Paths']['GDB']['GDB']=r'C:\Data\Geodatabases'
	meta['Paths']['GDB']['LandCover']=meta['Paths']['GDB']['GDB'] + '\\LandCover\\20230607\\LandCover.gdb'
	meta['Paths']['GDB']['LandUse']=meta['Paths']['GDB']['GDB'] + '\\LandUse\\20240529\\LandUse.gdb'
	meta['Paths']['GDB']['Disturbance']=meta['Paths']['GDB']['GDB'] + '\\Disturbances\\20250405\\Disturbances.gdb'
	#meta['Paths']['GDB']['Results']=meta['Paths']['GDB']['GDB'] + '\\Results\\20240528\\Results.gdb'
	meta['Paths']['GDB']['Results']=meta['Paths']['GDB']['GDB'] + '\\Results\\20250405\\Results.gdb'
	meta['Paths']['GDB']['VRI 2002']=meta['Paths']['GDB']['GDB'] + '\\VRI\\2002\\VRI2002_VEG_COMP_LYR_R1_POLY_FINAL_DELIVERYV4.gdb'
	meta['Paths']['GDB']['VRI 2015']=meta['Paths']['GDB']['GDB'] + '\\VRI\\2015\\VEG_COMP_LYR_R1_POLY.gdb'
	meta['Paths']['GDB']['VRI 2023']=meta['Paths']['GDB']['GDB'] + '\\VRI\\20230401\\VRI.gdb'
	#meta['Paths']['GDB']['VRI 2024']=meta['Paths']['GDB']['GDB'] + '\\VRI\\20240529\\VRI.gdb'
	#meta['Paths']['GDB']['VRI 2025']=meta['Paths']['GDB']['GDB'] + '\\VRI\\20250404\\VRI.gdb'

	meta['Paths']['DB']={}
	meta['Paths']['DB']['Taz']=r'C:\Users\rhember\Government of BC\External NRS Data Science and Modelling - General\Data\Events DB'
	meta['Paths']['DB']['na1k']=r'D:\Data\na1k'
	meta['Paths']['DB']['Beetles']=meta['Paths']['DB']['Taz'] + '\\Beetles'
	meta['Paths']['DB']['Harvest']=meta['Paths']['DB']['Taz'] + '\\Harvest'
	meta['Paths']['DB']['Waste']=r'C:\Data\Waste Wood\FromWasteSystem'
	meta['Paths']['DB']['HBS']=r'C:\Data\Harvest\HBS'
	meta['Paths']['DB']['Surface Climate']=r'C:\Users\rhember\Government of BC\External NRS Data Science and Modelling - General\Data\Surface Climate DB'
	meta['Paths']['DB']['CO2']=r'C:\Data\Carbon Dioxide Concentration\Carbon Dioxide Concentration Mid Year.xlsx'
	meta['Paths']['DB']['NDEP']=r'C:\Data\Nitrogen Deposition'
	meta['Paths']['DB']['CMIP6']=r'C:\Data\Climate\CMIP6\Monthly'
	meta['Paths']['DB']['20thCR']=r'C:\Data\Climate\Reanlysis\20th Century'
	meta['Paths']['DB']['CRU']=r'C:\Data\Climate\CRU'
	meta['Paths']['DB']['Field Plot Observations']=r'C:\Users\rhember\Government of BC\External NRS Data Science and Modelling - General\Data\FIeld Plot Observations DB'
	meta['Paths']['DB']['FAIR']=r'C:\Users\rhember\Government of BC\External NRS Data Science and Modelling - General\Data\FAIR DB'
	meta['Paths']['DB']['Nutrient Applications']=r'C:\Users\rhember\Government of BC\External NRS Data Science and Modelling - General\Data\Forest Nutrient Addition Experiments DB'

	meta['Paths']['Model']={}
	meta['Paths']['Model']['Code']=r'G:\My Drive\Code_Python\fcgadgets\cbrunner'
	meta['Paths']['Model']['Parameters']=meta['Paths']['Model']['Code'] + '\\Parameters'
	meta['Paths']['Model']['gromo']=r'D:\Data\na1k\gromo'

	meta['Paths']['Projects']={}
	meta['Paths']['Projects']['Demos']=r'D:\Modelling Projects'
	meta['Paths']['Projects']['BC-FCS']={}
	meta['Paths']['Projects']['BC-FCS']['Documentation']=r'C:\Users\rhember\Government of BC\External Forest Carbon Estimation - Documents\General\Projects\BC-FCS\Documentation\R2025'

	meta['Paths']['Field Plots']={}
	meta['Paths']['Field Plots']['DB']=r'C:\Data\Field Plots\PSP-NADB2'
	meta['Paths']['Field Plots']['Raw Data']={}
	#meta['Paths']['Field Plots']['Raw Data']['BC']=meta['Paths']['Field Plots']['DB'] + '\\Given\BC\Received 2023-03-02'
	#meta['Paths']['Field Plots']['Raw Data']['BC']=meta['Paths']['Field Plots']['DB'] + '\\Given\BC\\Recieved 2024-03-19'
	meta['Paths']['Field Plots']['Raw Data']['BC']=meta['Paths']['Field Plots']['DB'] + '\\Given\BC\\Recieved 2024-12-04'
	
	meta['Graphics']={'Plot Style':{},'Map':{},'Flowchart':{}}
	meta['Graphics']['Plot Style']='Web' # Manuscript
	meta['Graphics']['gp']=gu.SetGraphics(meta['Graphics']['Plot Style'])
	meta['Graphics']['gp']['AxesLetterStyle']='Caps'
	meta['Graphics']['gp']['AxesLetterFontWeight']='Bold'
	meta['Graphics']['Print Figures']='Off'
	meta['Graphics']['Print Figure Path']=r'C:\Users\rhember\OneDrive - Government of BC\Figures\BCFCS'

	meta['Graphics']['Colours']={'rgb':{},'hex':{}}
	meta['Graphics']['Colours']['rgb']['Blue Light']=matplotlib.colors.hex2color('#c5dcfa') # OLD: #e6f5ff
	meta['Graphics']['Colours']['rgb']['Blue Dark']=matplotlib.colors.hex2color('#004f99')
	meta['Graphics']['Colours']['rgb']['Green Light']=matplotlib.colors.hex2color('#d9ffcc')
	meta['Graphics']['Colours']['rgb']['Green Dark']=matplotlib.colors.hex2color('#134d00')
	meta['Graphics']['Colours']['rgb']['Green Neon']=matplotlib.colors.hex2color('#b8d911')
	meta['Graphics']['Colours']['rgb']['Red Light']=matplotlib.colors.hex2color('#ffc2b3')
	meta['Graphics']['Colours']['rgb']['Red Dark']=matplotlib.colors.hex2color('#801a00')
	meta['Graphics']['Colours']['rgb']['Brown Light']=matplotlib.colors.hex2color('#ecd9c5')
	meta['Graphics']['Colours']['rgb']['Brown Dark']=matplotlib.colors.hex2color('#4d3319')
	meta['Graphics']['Colours']['rgb']['Purple Light']=matplotlib.colors.hex2color('#e0ccff')
	meta['Graphics']['Colours']['rgb']['Purple Dark']=matplotlib.colors.hex2color('#5900b3')
	meta['Graphics']['Colours']['rgb']['Orange']=matplotlib.colors.hex2color('#d99b1e')
	meta['Graphics']['Colours']['rgb']['Hover']=matplotlib.colors.hex2color('#fffbb5')
	for k in meta['Graphics']['Colours']['rgb'].keys():
		meta['Graphics']['Colours']['hex'][k]=matplotlib.colors.to_hex(meta['Graphics']['Colours']['rgb'][k])

	# Defaults assume province-wide map
	meta['Graphics']['Map']['RGSF']=1
	meta['Graphics']['Map']['Fig Width']=9.75
	meta['Graphics']['Map']['Side Space']=0
	meta['Graphics']['Map']['Map Position']=[0,0,1-meta['Graphics']['Map']['Side Space']-0.01,1]
	meta['Graphics']['Map']['Map Axis Vis']='off'
	meta['Graphics']['Map']['Map Grid Vis']=False
	meta['Graphics']['Map']['Legend X']=0.72
	meta['Graphics']['Map']['Legend Y']=0.75
	meta['Graphics']['Map']['Legend Width']=0.03
	meta['Graphics']['Map']['Legend Font Size']=7
	meta['Graphics']['Map']['Legend Text Space']=0.035
	meta['Graphics']['Map']['Show Bound Land Mask']='On'
	meta['Graphics']['Map']['Show Bound Within']='Off'
	meta['Graphics']['Map']['Show Lakes']='Off'
	meta['Graphics']['Map']['Show Rivers']='Off'
	meta['Graphics']['Map']['Show Roads']='Off'
	meta['Graphics']['Map']['Show Rail']='Off'
	meta['Graphics']['Map']['Show Cities']='Off'
	meta['Graphics']['Map']['Show TPFs']='Off'
	meta['Graphics']['Map']['Show Symbol Labels']='Off'
	meta['Graphics']['Map']['Show Inset Map']='Off'
	meta['Graphics']['Map']['Show Scalebar']='Off'
	meta['Graphics']['Vector Import']={'Water Management':'Off','Wetland':'Off'}
	
	meta['Graphics']['Flowchart']={}
	meta['Graphics']['Flowchart']['Font Name']='Arial'
	meta['Graphics']['Flowchart']['Font Color']='#6b7d2a'
	meta['Graphics']['Flowchart']['Font Size']='10'
	meta['Graphics']['Flowchart']['Penwidth']='0.5'
	meta['Graphics']['Flowchart']['Node Background Color']='#f0fca2'
	meta['Graphics']['Flowchart']['Cluster Background Color']='#f4f5f2'
	meta['Graphics']['Flowchart']['Cluster Background Color 2']='#d1d1d1'
	
	meta['Graphics']['Field Plot Comparison']={}
	meta['Graphics']['Field Plot Comparison']['bl']=np.array([0.7,0.8,0.95])
	meta['Graphics']['Field Plot Comparison']['bd']=np.array([0.27,0.49,0.77])
	meta['Graphics']['Field Plot Comparison']['gl']=np.array([0.7,0.95,0]) # 0.85,1,0.65
	meta['Graphics']['Field Plot Comparison']['gd']=np.array([0.4,0.75,0])
	meta['Graphics']['Field Plot Comparison']['rl']=np.array([0.8,0.6,0.4])
	meta['Graphics']['Field Plot Comparison']['rd']=np.array([0.4,0.3,0.2])
	
	meta['Graphics']['Modelling']={}
	meta['Graphics']['Modelling']['AxesLetterStyle']='Default' #'NoPar'
	meta['Graphics']['Modelling']['AxesFontWeight']='Bold'
	
	# Tracking parameters
	meta['Graphics']['Fig Count']=1
	meta['Graphics']['Tab Count']=1
	
	# Initiate geospatial info
	meta['Geos']={}
	
	# Import variable info
	meta['Geos']['Variable Info']=gu.ReadExcel(meta['Paths']['Model']['Code'] + '\\Parameters\\Table_BCFCS_BC1haRasterVariableList.xlsx',sheet_name='Sheet1',skiprows=0)
	
	# Import coordinate reference system
	meta['Geos']['crs']=gpd.read_file(meta['Paths']['GDB']['GDB'] + '//lakes.geojson').crs
	
	# Import LUTs
	if 'LUT' not in meta:
		meta=ImportLUTs(meta)

	# Import model parameters
	meta=cbu.ImportParameters(meta)

	# Import field plot parameters
	if 'Include Field Plots' in argv:
		meta['Field Plots'],data,soc=ufp.ImportGroundPlotData(meta,type='Just Parameters',include_soil='True')
		meta['Field Plots']['Data']=data
		meta['Field Plots']['SOC']=soc
		del meta['Data']

	# Climate
	meta['Climate']={'SF':{}}
	meta['Climate']['SF']['cwd']=0.1
	meta['Climate']['SF']['cmi']=0.1
	meta['Climate']['SF']['ea']=0.01
	meta['Climate']['SF']['es']=0.01
	meta['Climate']['SF']['eta']=0.1
	meta['Climate']['SF']['etp']=0.1
	meta['Climate']['SF']['prcp']=0.1
	meta['Climate']['SF']['rswd']=0.01
	meta['Climate']['SF']['tmean']=0.01
	meta['Climate']['SF']['tmax']=0.01
	meta['Climate']['SF']['tmin']=0.01
	meta['Climate']['SF']['vpd']=0.01
	meta['Climate']['SF']['ws']=0.1
	meta['Climate']['SF']['wsp']=0.1
	meta['Climate']['SF']['runoff']=0.1
	meta['Climate']['SF']['melt']=0.1
	meta['Climate']['SF']['ndep']=0.01
	meta['Climate']['SF']['ndwi']=0.001
	meta['Climate']['SF']['lst']=0.01
	meta['Climate']['SF']['ca']=0.1
	meta['Climate']['SF']['albedo']=0.001
	meta['Climate']['SF']['AbsorptionRSW']=0.1
	meta['Climate']['SF']['ndvi']=0.001
	meta['Climate']['Missing Number']=-99

	# Gromo (growth and mortality fitted against field plots)
	if 'Modules' not in meta:
		meta['Modules']={}
		meta['Modules']['gromo']={}
		meta['Modules']['gromo']['Species CD']=['PL','SW','BL','FDC','FDI','HWC','HWI','CW','SE','AT','SB','SS','AE','O']
		meta['Modules']['gromo']['SF']=0.01

	# Ancillary tables
	meta['Tables']={}
	fileL=os.listdir(meta['Paths']['Model']['Parameters'])
	for file in fileL:
		if file[0:6]=='Table_':
			meta['Tables'][file[6:-5]]=gu.ReadExcel(meta['Paths']['Model']['Parameters'] + '\\' + file,sheet_name='Sheet1',skiprows=0)

	return meta

#%% Import variables for ROI
def Import_Raster(meta,roi,vList,*argv):
	d={}
	for v in vList:
		if roi!=[]:
			if 'grd' in roi.keys():
				if v in roi['grd'].keys():
					continue
		if v=='access':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Management\\AccessZones.tif')
		elif v=='age_ntems':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Age\\Age_NTEM_2019.tif')
		elif v=='age_vri02':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2002\\PROJ_AGE_1.tif')
		elif v=='age_vri15':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2015\\PROJ_AGE_1.tif')
		elif v=='age_vri23':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\PROJ_AGE_1.tif')
		elif v=='alr':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\OATS_ALR_POLYS\\ALR_POLY_ID.tif')
		elif v=='aspect':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Terrain\\aspect.tif')
		elif v=='aset':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_ASET_All.tif')
		elif v=='aset_post17':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_ASET_Post2017.tif')
		elif v=='aset_no_all':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_NonOb_ASET_All.tif')
		elif v=='aset_no_2024':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_NonOb_ASET_2024.tif')
		elif v=='BEC_ZONE_CODE':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\BEC_ZONE_CODE.tif')
		elif v=='bdfrac':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\BroadleafDeciduousFraction.tif')
		elif v=='bdfrac_2049s4':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\BDF_Comp1_2049_Scn4.tif')
		elif v=='bgcz':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\BEC_BIOGEOCLIMATIC_POLY\\ZONE_GapFilled.tif')
		elif v=='biomass_glob':
			d[v]=gis.OpenGeoTiff(r'C:\Data\Biomass\GlobBiomass\N80W140_agb pc8.tif')
			d[v]['Data']=np.squeeze(d[v]['Data'])
			d[v]['Data']=0.5*d[v]['Data']
		elif v=='biomass_sw_vri23':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\WHOLE_STEM_BIOMASS_PER_HA.tif')
		elif v=='biomass_bk_vri23':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\BARK_BIOMASS_PER_HA.tif')
		elif v=='biomass_br_vri23':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\BRANCH_BIOMASS_PER_HA.tif')
		elif v=='biomass_f_vri23':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\FOLIAGE_BIOMASS_PER_HA.tif')
		elif v=='bsr_yr':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\PROT_HISTORICAL_FIRE_POLYS_SP\\BurnSevComp_YearLast.tif')
		elif v=='bsr_sc':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\PROT_HISTORICAL_FIRE_POLYS_SP\\BurnSevComp_SevClassLast.tif')
		elif v=='citym':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\BC_MAJOR_CITIES_POINTS_500M\\NAME.tif')
		elif v=='cmi_mjjas_n':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\Climate\Summaries\\Normals\\bc1ha_cmi_mjjas_norm_1971to2000.tif')
			d[v]['Data']=d[v]['Data'].astype('float')*meta['Climate']['SF']['cmi']
		elif v=='crownc':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\CROWN_CLOSURE.tif')
		elif v=='cwd_ann_n':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\Climate\Summaries\\Normals\\bc1ha_cwd_ann_norm_1971to2000.tif')
			d[v]['Data']=d[v]['Data'].astype('float')*meta['Climate']['SF']['cwd']
		elif v=='cwd_mjjas_n':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\Climate\Summaries\\Normals\\bc1ha_cwd_mjjas_norm_1971to2000.tif')
			d[v]['Data']=d[v]['Data'].astype('float')*meta['Climate']['SF']['cwd']
		elif v=='cmi_ann_n':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\Climate\Summaries\\Normals\\bc1ha_cmi_ann_norm_1971to2000.tif')
			d[v]['Data']=d[v]['Data'].astype('float')*meta['Climate']['SF']['cmi']
		elif v=='cwd_mjjas_n':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\Climate\Summaries\\Normals\\bc1ha_cmi_mjjas_norm_1971to2000.tif')
			d[v]['Data']=d[v]['Data'].astype('float')*meta['Climate']['SF']['cmi']
		elif v=='d2road':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Management\\DistanceFromRoads.tif')
		elif v=='d2fac':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Management\\DistanceFromForestryFacility.tif')
		elif v=='d2set':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Management\\DistanceFromSettlement.tif')
		elif v=='elev':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Terrain\\elevation.tif')
		elif v=='etp_ann_n':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Climate\\Summaries\\Normals\\bc1ha_etp_ann_norm_1971to2000.tif')
			d[v]['Data']=d[v]['Data'].astype('float')*meta['Climate']['SF']['etp']
		elif v=='etp_mjjas_n':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Climate\\Summaries\\Normals\\bc1ha_etp_mjjas_norm_1971to2000.tif')
			d[v]['Data']=d[v]['Data'].astype('float')*meta['Climate']['SF']['etp']
		elif v=='ezcan':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Ecozones_Canada\\Ecozones_Canada.tif')
		elif v=='feca_yr':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\FE-CA_YearLast.tif')
		elif v=='feca_mask':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\FE-CA_MaskAll.tif')
		elif v=='fire_yf':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\PROT_HISTORICAL_FIRE_POLYS_SP\\WildfireComp1_YearFirst.tif')
		elif v=='fire_yl':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\PROT_HISTORICAL_FIRE_POLYS_SP\\WildfireComp1_YearLast.tif')
		elif v=='fire_2023':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\PROT_HISTORICAL_FIRE_POLYS_SP\\WildfireComp1_YEAR_2023.tif')
		elif v=='forfrac5m_Tom23':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\CanopyHeight\\ForestMask_2020_ThreshHeight5m.tif')
		elif v=='forfrac2m_Tom23':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\CanopyHeight\\ForestMask_2020_ThreshHeight2m.tif')
		elif v=='fc_yr':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_FOREST_COVER_INV_SVW\\REFERENCE_YEAR.tif')
		elif v=='fc_ssc':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_FOREST_COVER_INV_SVW\\ID_STOCKING_STATUS_CODE.tif')
		elif v=='fc_stc':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_FOREST_COVER_INV_SVW\\ID_STOCKING_TYPE_CODE.tif')
		elif v=='gfc_tc2000':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Global Forest Change\\bc1ha_gfc_treecover2000.tif')
		elif v=='gfc_ly':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Global Forest Change\\bc1ha_gfc_lossyear.tif')
		elif v=='gfc_loss':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Global Forest Change\\bc1ha_gfc_loss.tif')
		elif v=='gfc_gain':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Global Forest Change\\bc1ha_gfc_gain.tif')
		elif v=='gromo_g_n':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\gromo\\Normals\\GrowthGross_ActualAge_1971to2000mean.tif')
			d[v]['Data']=d[v]['Data'].astype('float')*0.01
		elif v=='gromo_net':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\gromo\\GrowthNet_Actual_2023.tif')
			d[v]['Data']=d[v]['Data'].astype('float')*0.01
		elif v=='gromo_net50':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\gromo\\GrowthNet_A50_2023.tif')
			d[v]['Data']=d[v]['Data'].astype('float')*0.01
		elif v=='gromo_m':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\gromo\\Mortality_Actual_2023.tif')
			d[v]['Data']=d[v]['Data'].astype('float')*0.01
		elif v=='gsoc':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Soil\\gsoc2010_bc1ha.tif')
		elif v=='geomorph':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Terrain\\FromBurtWilliam\\geomorphons_search50.tif')
		elif v=='h_vri02':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2002\\PROJ_HEIGHT_1.tif')
		elif v=='h_vri15':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2015\\PROJ_HEIGHT_1.tif')
		elif v=='h_vri23':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\PROJ_HEIGHT_1.tif')
		elif v=='h_Tolan23':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\CanopyHeight\\h_Tolan23.tif')
		elif v=='h_Lang22':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\CanopyHeight\\h_Lang22.tif')
		elif v=='harv_yr_cc':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VEG_CONSOLIDATED_CUT_BLOCKS_SP\\HARVEST_YEAR_YearLast.tif')
		elif v=='harv_yr_ntems':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Disturbances\\Harvest_NTEM_Year.tif')
		elif v=='harv_yr_early':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Disturbances\\Harvest_EarlyReconstruction_Year.tif')
		elif v=='harv_yr_comp1':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Disturbances\\Harvest_Comp1_Year.tif')
		elif v=='harv_yr_comp2':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Disturbances\\Harvest_Comp2_Year.tif')
		elif v=='harv_salv':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Disturbances\\HarvestSalvageMask_FromCruise.tif')
		elif v=='harv_prob':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Disturbances\\HarvestProbability.tif')
		elif v=='hlz':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Climate\\Summaries\\Normals\\bc1ha_HoldridgeLifeZones_1971to2000.tif')
		elif v=='kd_yl':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\SP-KD_YearLast.tif')
		elif v=='ibm_yr':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\PEST_INFESTATION_POLY\\PEST_SEVERITY_CODE_IBM_YearLast.tif')
		elif v=='ibd_yr':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\PEST_INFESTATION_POLY\\PEST_SEVERITY_CODE_IBD_YearLast.tif')
		elif v=='ibb_yr':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\PEST_INFESTATION_POLY\\PEST_SEVERITY_CODE_IBB_YearLast.tif')
		elif v=='ibs_yr':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\PEST_INFESTATION_POLY\\PEST_SEVERITY_CODE_IBS_YearLast.tif')
		elif v=='idw_yr':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\PEST_INFESTATION_POLY\\PEST_SEVERITY_CODE_IDW_YearLast.tif')
		elif v=='idl_yr':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\PEST_INFESTATION_POLY\\PEST_SEVERITY_CODE_IDL_YearLast.tif')
		elif v=='lc_vri_l2':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\BCLCS_LEVEL_2.tif')
		elif v=='lc_vri_l3':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\BCLCS_LEVEL_3.tif')
		elif v=='lc_vri_l4':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\BCLCS_LEVEL_4.tif')
		elif v=='lc_vri_l5':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\BCLCS_LEVEL_5.tif')
		elif v=='lc_vri_l5_02':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2002\\BCLCS_LEVEL_5.tif')
		elif v=='lc_comp1_1800':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCover_Comp1_1800.tif')
		elif v=='lc_comp1_2019':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCover_Comp1_2019.tif')
		elif v=='lc_comp1_vri':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCover_Comp1_2022_FromVRI2022.tif')
		elif v=='lc_comp1_2049s1':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCover_Comp1_2049_Scn1.tif')
		elif v=='lc_comp1_2049s2':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCover_Comp1_2049_Scn2.tif')
		elif v=='lc_comp1_2049s3':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCover_Comp1_2049_Scn3.tif')
		elif v=='lc_comp1_2049s4':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCover_Comp1_2049_Scn4.tif')
		elif v=='lc_cec_2020':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCover_CEC_2020_Compressed.tif')
		elif v=='lc_cec_2010':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCover_CEC_2010_Compressed.tif')
		elif v=='lc_ntems_2019':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCover_NTEMS_2019.tif')
		elif v=='lc_ntems_2019_recl':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCover_Comp1_2019_FromNTEMS.tif')
		elif v=='lc_vri_recl':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCover_Comp1_2022_FromVRI2022.tif')
		elif v=='lu_comp1_2019':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandUse_Comp1_2019.tif')
		elif v=='lu_comp1_2049s1':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandUse_Comp1_2049_Scn1.tif')
		elif v=='lu_comp1_2049s2':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandUse_Comp1_2049_Scn2.tif')
		elif v=='lu_comp1_2049s3':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandUse_Comp1_2049_Scn3.tif')
		elif v=='lu_comp1_2049s4':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandUse_Comp1_2049_Scn4.tif')
		elif v=='luc1_hist_yr':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Disturbances\\LandUseChange_Comp1_1800to2019_Year.tif')
		elif v=='luc1_hist_type':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Disturbances\\LandUseChange_Comp1_1800to2019_Type.tif')
		elif v=='luc1_1019_type':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Disturbances\\LandUseChange_Comp1_2010to2019_Type.tif')
		elif v=='luc1_1019_yr':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Disturbances\\LandUseChange_Comp1_2010to2019_Year.tif')
		elif v=='luc1_fut_s1_yr':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Disturbances\\LandUseChange_Comp1_2020to2049_Scn1_Year.tif')
		elif v=='luc1_fut_s1_type':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Disturbances\\LandUseChange_Comp1_2020to2049_Scn1_Type.tif')
		elif v=='luc_aff_cec':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Disturbances\\LandUseChange_Afforestation_10to20_CEC.tif')
		elif v=='melt_ann_n':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Climate\\Summaries\\Normals\\bc1ha_melt_ann_norm_1971to2000.tif')
			d[v]['Data']=d[v]['Data'].astype('float')*meta['Climate']['SF']['melt']
		elif v=='mines':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\HSP_MJR_MINES_PERMTTD_AREAS_SP\\STATUS_TYPE.tif')
		elif v=='munic':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\TA_MUNICIPALITIES_SVW\MUNICIPALITY_NAME.tif')
		elif v=='ogma':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RMP_OGMA_LEGAL_ALL_SVW\\LEGAL_OGMA_PROVID.tif')
		elif v=='ogdef':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\OGSR_TAP_PRIORITY_DEF_AREA_SP\\PRIORITY_DEFERRAL_ID.tif')
		elif v=='own':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\F_OWN\\OWNERSHIP_DESCRIPTION.tif')
		elif v=='pfi_c':
			d[v]=gis.OpenGeoTiff(r'C:\Data\ForestInventory\LiDAR\Boundary TSA\From GQ 20230302\BoundaryTSA_PFI_LnTotalWSV_bc1ha.tif')
			d[v]['Data']=np.squeeze(d[v]['Data'])
			d[v]['Data']=0.5*0.5*d[v]['Data']
		elif v=='pdead_cruise':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Disturbances\\HarvestPercentDead_FromCruise.tif')
		elif v=='pdead_vri23':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\STAND_PERCENTAGE_DEAD.tif')
		elif v=='plam':
			d[v]=gis.OpenGeoTiff(r'C:\Data\BC1ha\Results\Planting_FromRESULTS_MaskCount.tif')
		elif v=='pl_yl':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_YearLast.tif')
		elif v=='popp':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\NRC_POPULATED_PLACES_1M_SP\\NAME.tif')
		elif v=='prcp_ann_n':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\Climate\\Summaries\\Normals\\bc1ha_prcp_ann_norm_1971to2000.tif')
			d[v]['Data']=d[v]['Data'].astype('float')*meta['Climate']['SF']['prcp']
		elif v=='prcp_mjjas_n':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\Climate\\Summaries\\Normals\\bc1ha_prcp_mjjas_norm_1971to2000.tif')
			d[v]['Data']=d[v]['Data'].astype('float')*meta['Climate']['SF']['prcp']
		elif v=='protected':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\PROTECTED_LANDS_DESIGNATION.tif')
		elif v=='park':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\TA_PARK_ECORES_PA_SVW\\ADMIN_AREA_SID.tif')
		elif v=='parknat':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\CLAB_NATIONAL_PARKS\\ENGLISH_NAME.tif')
		elif v=='peat':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\peat.tif')
		elif v=='rangecon':
			d[v]=gis.OpenGeoTiff(r'C:\Data\BC1ha\LandCoverUse\rangeten_consol.tif')
		elif v=='rears':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\REARs_CompPlusProp.tif')
		elif v=='refg':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
		elif v=='refyear':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\REFERENCE_YEAR.tif')
		elif v=='regentype':
			d[v]=gis.OpenGeoTiff(r'C:\Data\BC1ha\Results\Harvest_Regen_Type.tif')
		elif v=='reserve_comp1':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Management\\ReserveComp1.tif')
		elif v=='rd':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\TA_REGIONAL_DISTRICTS_SVW\\REGIONAL_DISTRICT_NAME.tif')
		elif v=='road_atl':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\DRA_DGTL_ROAD_ATLAS_MPAR_SP\\DIGITAL_ROAD_ATLAS_LINE_ID.tif')
		elif v=='road_ften_s':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\FTEN_ROAD_SEGMENT_LINES_SVW\\LIFE_CYCLE_STATUS_CODE.tif')
		elif v=='roadrehab_yl':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\RD-DE_YearLast.tif')
		elif v=='rswd_ann_n':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\Climate\\Summaries\\bc1ha_rswd_ann_norm_1971to2000.tif')
			d[v]['Data']=d[v]['Data'].astype('float')*meta['Climate']['SF']['rswd']
		elif v=='rswd_mjjas_n':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\Climate\\Summaries\\bc1ha_rswd_mjjas_norm_1971to2000.tif')
			d[v]['Data']=d[v]['Data'].astype('float')*meta['Climate']['SF']['rswd']
		elif v=='runoff_ann_n':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Climate\\Summaries\\Normals\\bc1ha_runoff_ann_norm_1971to2000.tif')
			d[v]['Data']=d[v]['Data'].astype('float')*meta['Climate']['SF']['runoff']
		elif v=='si_vri23':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\SITE_INDEX.tif')
		elif v=='slope':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Terrain\\Slope.tif')
		elif v=='slope':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Terrain\\slope.tif')
		elif v=='sp_rip':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\SP-Rip_YearLast.tif')
		elif v=='sp_disc':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\SP-DISC_YearLast.tif')
		elif v=='sp_broad':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\SP-BU-BROAD_YearLast.tif')
		elif v=='sp_wing':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\SP-WING_YearLast.tif')
		elif v=='sp_mound':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\SP-MOUND_YearLast.tif')
		elif v=='sphl_vri23':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI\\sphlive.tif')
		elif v=='sphd_vri23':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI\\sphdead.tif')
		elif v=='soc':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Soil\\soc_tot_forest_Shawetal2018.tif')
		elif v=='spc1_ntems':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Species\\NTEMS_Spc1_Filled_RD_CAPITAL.tif')
		elif v=='spc1_vri02':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2002\\SPECIES_CD_1.tif')
		elif v=='spc2_vri02':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2002\\SPECIES_CD_2.tif')
		elif v=='spc3_vri02':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2002\\SPECIES_CD_3.tif')
		elif v=='spc1_pct_vri02':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2002\\SPECIES_PCT_1.tif')
		elif v=='spc2_pct_vri02':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2002\\SPECIES_PCT_2.tif')
		elif v=='spc3_pct_vri02':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2002\\SPECIES_PCT_3.tif')
		elif v=='spc1_vri15':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2015\\SPECIES_CD_1.tif')
		elif v=='spc1_vri23':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\SPECIES_CD_1.tif')
		elif v=='spc2_vri23':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\SPECIES_CD_2.tif')
		elif v=='spc3_vri23':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\SPECIES_CD_3.tif')
		elif v=='spc4_vri23':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\SPECIES_CD_4.tif')
		elif v=='spc5_vri23':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\SPECIES_CD_5.tif')
		elif v=='spc6_vri23':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\SPECIES_CD_6.tif')
		elif v=='spc1_pct_vri23':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\SPECIES_PCT_1.tif')
		elif v=='spc2_pct_vri23':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\SPECIES_PCT_2.tif')
		elif v=='spc3_pct_vri23':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\SPECIES_PCT_3.tif')
		elif v=='spc4_pct_vri23':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\SPECIES_PCT_4.tif')
		elif v=='spc5_pct_vri23':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\SPECIES_PCT_5.tif')
		elif v=='spc6_pct_vri23':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\SPECIES_PCT_6.tif')
		elif v=='si_spl_ntems':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\SPL\\SI_From_NTEMS_Spc1_Filled_RD_CAPITAL_GF.tif')
		elif v=='si_spl_fd':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\SPL\\Site_Prod_Fd.tif')
		elif v=='si_spl_hw':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\SPL\\Site_Prod_Hw.tif')
		elif v=='ssc':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_OPENING_SVW\\SILV_SYSTEM_CODE.tif')
		elif v=='tdc_vri02':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\TreeDensityClass_2002.tif')
		elif v=='tdc_vri15':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\TreeDensityClass_VRI2015.tif')
		elif v=='tdc_vri23':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\TreeDensityClass_Current.tif')
		elif v=='tdgt30_vri23':
			tmp=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\VRI_LIVE_STEMS_PER_HA.tif')
			d[v]=copy.deepcopy(tmp)
			d[v]['Data']=0*d[v]['Data']
			d[v]['Data'][(tmp['Data']>3000)]=1
		elif v=='tdc_wsg':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\TreeDensityClass_WithShrubsGrasses.tif')
		elif v=='tm':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\FTEN_CUT_BLOCK_POLY_SVW\\TIMBER_MARK.tif')
		elif v=='tmean_ann_n':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\Climate\Summaries\\Normals\\BC1ha_tmean_ann_norm_1971to2000.tif')
			d[v]['Data']=d[v]['Data'].astype('float')*meta['Climate']['SF']['tmean']
		elif v=='tmean_mjjas_n':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\Climate\Summaries\\Normals\\BC1ha_tmean_mjjas_norm_1971to2000.tif')
			d[v]['Data']=d[v]['Data'].astype('float')*meta['Climate']['SF']['tmean']
		elif v=='transl':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\GBA_TRANSMISSION_LINES_SP\\TRANSMISSION_LINE_ID.tif')
		elif v=='twi':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Terrain\\twi.tif')
		elif v=='upwetf_ntems':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoveruse\\UplandWetlandForest_NTEMS.tif')
		elif v=='upwetf_vri':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoveruse\\UplandWetlandForest_VRI.tif')
		elif v=='v_vri02':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2002\\VRI_LIVE_STEMS_PER_HA.tif')
		elif v=='wbt':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\FWA_WETLANDS_POLY\\WATERBODY_TYPE.tif')
		elif v=='wshed':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Terrain\\WatershedsNamed.tif')
		elif v=='wbmm':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\FWA_MANMADE_WATERBODIES_POLY\\WATERBODY_TYPE.tif')
		elif v=='ws_ann_n':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Climate\\Summaries\\Normals\\bc1ha_ws_ann_norm_1971to2000.tif')
			d[v]['Data']=d[v]['Data'].astype('float')*meta['Climate']['SF']['ws']
		elif v=='ws_mjjas_n':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Climate\\Summaries\\Normals\\bc1ha_ws_mjjas_norm_1971to2000.tif')
			d[v]['Data']=d[v]['Data'].astype('float')*meta['Climate']['SF']['ws']
		elif v=='wsp_ann_n':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Climate\\Summaries\\Normals\\bc1ha_wsp_ann_norm_1971to2000.tif')
			d[v]['Data']=d[v]['Data'].astype('float')*meta['Climate']['SF']['wsp']
		elif v=='wsp_mjjas_n':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Climate\\Summaries\\Normals\\bc1ha_wsp_mjjas_norm_1971to2000.tif')
			d[v]['Data']=d[v]['Data'].astype('float')*meta['Climate']['SF']['wsp']
		else:
			pass

		if meta['Graphics']['Map']['RGSF']!=1:
			d[v]=gis.UpdateGridCellsize(d[v],meta['Graphics']['Map']['RGSF'])

		if roi!=[]:
			if 'points' in roi.keys():
				if v==vList[0]:
					iPoints=gis.GetGridIndexToPoints(d[vList[0]],roi['points']['x'],roi['points']['y'])
				d[v]=d[v]['Data'][iPoints]

	if 'Extract Grid' in argv:
		for v in vList:
			d[v]=d[v]['Data']

	if roi!=[]:
		if 'grd' in roi.keys():
			for v in vList:
				if v in roi['grd'].keys():
					continue
				# Add to ROI structure and clip to ROI
				roi['grd'][v]=d[v]
				try:
					roi['grd'][v]=gis.ClipToRaster(roi['grd'][v],roi['grd'])
				except:
					print('Failed, Try using RGSF = 1')
			return roi
		elif 'points' in roi.keys():
			return d

	else:
		return d

#%% Look up tables
def ImportLUTs(meta):

	if 'LUT' not in meta.keys():
		meta['LUT']={}

	# Source databases
	uL=np.unique(meta['Geos']['Variable Info']['Layer Name'])
	for iL in range(uL.size):
		indL=np.where( (meta['Geos']['Variable Info']['Layer Name']==uL[iL]) & (meta['Geos']['Variable Info']['LUT Required']=='Yes') )[0]
		if indL.size==0:
			# Not categorical variables
			continue
		try:
			meta['LUT'][uL[iL]]=gu.ipickle(meta['Paths']['bc1ha'] + '\\LUTs\\LUTs_' + uL[iL] + '.pkl')
		except:
			print(uL[iL])

	# Override pest severity so that it is in order
	meta['LUT']['PEST_INFESTATION_POLY']['PEST_SEVERITY_CODE']={'T':1,'L':2,'M':3,'S':4,'V':5,'G':6}

	# Raw and Derived layer
	lNam='Derived'
	meta['LUT'][lNam]={}
	meta['LUT']['Raw']={}

	# Peatland map
	meta['LUT']['Raw']['Peat']=gu.ReadExcel(meta['Paths']['Model']['Parameters'] + '\\LUT_Peatlands.xlsx',sheet_name='Sheet1',skiprows=0)

	# Geomorphons
	d=gu.ReadExcel(meta['Paths']['Model']['Parameters'] + '\\LUT_Geomorphons.xlsx',sheet_name='Sheet1',skiprows=0)
	vNam='geomorph'
	meta['LUT'][lNam][vNam]={}
	for i in range(d['ID'].size):
		meta['LUT'][lNam][vNam][ d['Name'][i] ]=d['ID'][i]

	# Land cover - VRI level 2
	vNam='lc_vri_l2'
	meta['LUT'][lNam][vNam]={}
	meta['LUT'][lNam][vNam]['Water']=1
	meta['LUT'][lNam][vNam]['Land']=2
	meta['LUT'][lNam][vNam]['Non-treed']=3
	meta['LUT'][lNam][vNam]['Treed']=4
		
	# Land cover - VRI level 4
	vNam='lc_vri_l4'
	meta['LUT'][lNam][vNam]={}
	meta['LUT'][lNam][vNam]['Water']=1
	meta['LUT'][lNam][vNam]['Land']=2
	meta['LUT'][lNam][vNam]['Non-treed']=3
	meta['LUT'][lNam][vNam]['Treed']=4

	# Land Cover Compilation 1
	d=gu.ReadExcel(meta['Paths']['Model']['Parameters'] + '\\LUT_lc_comp1.xlsx',sheet_name='Sheet1',skiprows=0)
	vNam='lc_comp1'
	meta['LUT']['Raw'][vNam]=d
	meta['LUT'][lNam][vNam]={}
	for i in range(d['Value'].size):
		meta['LUT'][lNam][vNam][ d['Name'][i] ]=d['Value'][i]
	
	# Upland-wetland forest mask
	d=gu.ReadExcel(meta['Paths']['Model']['Parameters'] + '\\LUT_uplandwetland.xlsx',sheet_name='Sheet1',skiprows=0)
	vNam='upwetf'
	meta['LUT']['Raw'][vNam]=d
	meta['LUT'][lNam][vNam]={}
	for i in range(d['Value'].size):
		meta['LUT'][lNam][vNam][ d['Name'][i] ]=d['Value'][i]

	# Land Cover Class - NTEMS
	d=gu.ReadExcel(meta['Paths']['Model']['Parameters'] + '\\LUT_lc_NTEMS.xlsx',sheet_name='Sheet1',skiprows=0)
	vNam='lc_ntems_2019'
	meta['LUT']['Raw'][vNam]=d
	meta['LUT'][lNam][vNam]={}
	for i in range(d['Value'].size):
		meta['LUT'][lNam][vNam][ d['Name'][i] ]=d['Value'][i]

	# Land cover 2020 (CEC Given)
	d=gu.ReadExcel(meta['Paths']['Model']['Parameters'] + '\\LUT_lc_cec.xlsx',sheet_name='Sheet1',skiprows=0)
	vNam='lc_cec'
	meta['LUT']['Raw'][vNam]=d
	meta['LUT'][lNam][vNam]={}
	for i in range(d['Value'].size):
		meta['LUT'][lNam][vNam][ d['Name'][i] ]=d['Value'][i]

	# Land cover 2020 (CEC Compressed)
	d=gu.ReadExcel(meta['Paths']['Model']['Parameters'] + '\\LUT_lc_cec_Compressed.xlsx',sheet_name='Sheet1',skiprows=0)
	vNam='lc_cec_c'
	meta['LUT']['Raw'][vNam]=d
	meta['LUT'][lNam][vNam]={}
	for i in range(d['Value'].size):
		meta['LUT'][lNam][vNam][ d['Name'][i] ]=d['Value'][i]

	# Land Use Compilation 1
	d=gu.ReadExcel(meta['Paths']['Model']['Parameters'] + '\\LUT_lu_comp1.xlsx',sheet_name='Sheet1',skiprows=0)
	vNam='lu_comp1'
	meta['LUT']['Raw'][vNam]=d
	meta['LUT'][lNam][vNam]={}
	for i in range(d['Value'].size):
		meta['LUT'][lNam][vNam][ d['Name'][i] ]=d['Value'][i]

	# Land Cover Land Use Compilation 1
	d=gu.ReadExcel(meta['Paths']['Model']['Parameters'] + '\\LUT_lclu_change_comp1.xlsx',sheet_name='Sheet1',skiprows=0)
	vNam='lclu_chng_comp1'
	meta['LUT']['Raw'][vNam]=d
	meta['LUT'][lNam][vNam]={}
	for i in range(d['Value'].size):
		meta['LUT'][lNam][vNam][ d['Name'][i] ]=d['Value'][i]

	# Species leading - NTEMS
	d=gu.ReadExcel(meta['Paths']['Model']['Parameters'] + '\\LUT_spc1_NTEMS.xlsx',sheet_name='Sheet1',skiprows=0)
	vNam='spc1_ntems'
	meta['LUT']['Raw'][vNam]=d
	meta['LUT'][lNam][vNam]={}
	for i in range(d['Value'].size):
		meta['LUT'][lNam][vNam][ d['Name'][i] ]=d['Value'][i]

	# BGC Zone / NDT Zone Combo
	d=gu.ReadExcel(meta['Paths']['Model']['Parameters'] + '\\LUT_bgc_zone.xlsx',sheet_name='Sheet1',skiprows=0)
	vNam='bgc_zone'
	meta['LUT']['Raw'][vNam]=d

	# BGC Zone / NDT Zone Combo
	d=gu.ReadExcel(meta['Paths']['Model']['Parameters'] + '\\LUT_bgcz_ndt_combo.xlsx',sheet_name='Sheet1',skiprows=0)
	vNam='bgc-ndt'
	meta['LUT']['Raw'][vNam]=d
	meta['LUT'][lNam][vNam]={}
	for i in range(d['ID'].size):
		meta['LUT'][lNam][vNam][ d['BGC-NDT'][i] ]=d['ID'][i]

	# Tree Density Class
	vNam='tdc'
	meta['LUT'][lNam][vNam]={'Sparse':1,'Open':2,'Dense':3}

	# Access zones
	d=gu.ReadExcel(meta['Paths']['Model']['Parameters'] + '\\LUT_access.xlsx',sheet_name='Sheet1',skiprows=0)
	vNam='access'
	meta['LUT']['Raw'][vNam]=d
	meta['LUT'][lNam][vNam]={}
	for i in range(d['ID'].size):
		meta['LUT'][lNam][vNam][ d['Name'][i] ]=d['ID'][i]

	# Forest cover reserve compilation
	d=gu.ReadExcel(meta['Paths']['Model']['Parameters'] + '\\LUT_reserve_comp1.xlsx',sheet_name='Sheet1',skiprows=0)
	vNam='reserve_comp1'
	meta['LUT']['Raw'][vNam]=d
	meta['LUT'][lNam][vNam]={}
	for i in range(d['ID'].size):
		meta['LUT'][lNam][vNam][ d['Name'][i] ]=d['ID'][i]
	
	# Harvest regeneration type
	# Change name to ASET
# 	d=gu.ReadExcel(meta['Paths']['Model']['Parameters'] + '\\LUT_HarvestRegenType.xlsx')
# 	vNam='HarvestRegenType'
# 	meta['LUT']['Raw'][vNam]=d
# 	meta['LUT'][lNam][vNam]={}
# 	for i in range(d['ID'].size):
# 		meta['LUT'][lNam][vNam][ d['Name'][i] ]=d['ID'][i]

	# Planting regeneration type (used in NOSE modelling)
	d=gu.ReadExcel(meta['Paths']['Model']['Parameters'] + '\\LUT_ASET.xlsx',sheet_name='Sheet1',skiprows=0)
	vNam='ASET'
	meta['LUT']['Raw'][vNam]=d
	meta['LUT'][lNam][vNam]={}
	for i in range(d['ID'].size):
		meta['LUT'][lNam][vNam][ d['Name'][i] ]=d['ID'][i]

	# Forest nutrient management
	d=gu.ReadExcel(meta['Paths']['Model']['Parameters'] + '\\LUT_FNM.xlsx',sheet_name='Sheet1',skiprows=0)
	vNam='FNM'
	meta['LUT']['Raw'][vNam]=d
	meta['LUT'][lNam][vNam]={}
	for i in range(d['ID'].size):
		meta['LUT'][lNam][vNam][ d['Name'][i] ]=d['ID'][i]

	d=gu.ReadExcel(meta['Paths']['Model']['Parameters'] + '\\LUT_EcozoneCanada.xlsx',sheet_name='Sheet1',skiprows=0)
	vNam='ezcan'
	meta['LUT']['Raw'][vNam]=d
	meta['LUT'][lNam][vNam]={}
	for i in range(d['ID'].size):
		meta['LUT'][lNam][vNam][ d['Name'][i] ]=d['ID'][i]

	# Burn severity class compelation 1
	d=gu.ReadExcel(meta['Paths']['Model']['Parameters'] + '\\LUT_burnsev_comp1.xlsx',sheet_name='Sheet1',skiprows=0)
	vNam='burnsev_comp1'
	meta['LUT']['Raw'][vNam]=d
	meta['LUT'][lNam][vNam]={}
	for i in range(d['Value'].size):
		meta['LUT'][lNam][vNam][ d['Name'][i] ]=d['Value'][i]

	# Log grades
	d=gu.ReadExcel(meta['Paths']['Model']['Parameters'] + '\\Parameters_LogGrades.xlsx',sheet_name='Sheet1',skiprows=0)
	vNam='LogGrade'
	meta['LUT'][lNam][vNam]={}
	for i in range(d['ID'].size):
		meta['LUT'][lNam][vNam][ d['Code'][i] ]=d['ID'][i]

	# Forest Cover Stocking Type
	#d=gu.ReadExcel(r'C:\Data\BC1ha\Results\stocktype.tif.vat.xlsx')
	#d=gu.ReadExcel(r'C:\Code_Python\fcgadgets\cbrunner\Parameters\Parameters_ForestCover_StockingType.xlsx')
	# nam='fcst'
	# meta['LUT'][nam]={}
	# for i in range(d['VALUE'].size):
	#	 meta['LUT'][nam][ d['STOCKING_T'][i] ]=d['VALUE'][i]

	return meta

#%% Convert LUT number to string
def lut_n2s(dc,numb):
	if numb!=-999:
		vals=np.fromiter(dc.values(),dtype=float)
		keys=np.fromiter(dc.keys(),dtype='<U70')
		ind=np.where(vals==numb)[0]
		s=keys[ind]
	else:
		s=np.array(['Unidentified'],ndmin=1)
	return s

#%% Simplify geodatabases used for mapping
def SimplifyProvincialGDBs(meta,gdf):
	vL=['road','rivers','lakes','ogp','ogf']
	for v in vL:
		a=gdf[v]['gdf'].simplify(100)
		a.to_file(meta['Paths']['GDB']['GDB'] + '\\' + v + '.geojson',driver='GeoJSON')
	return

#%% Simplify road and line geodatabases
def SimplifyRoadGDBs(meta):
	vL=['FTEN_ROAD_SEGMENT_LINES_SVW','OG_ROAD_SEGMENT_PERMIT_SP','DRP_OIL_GAS_PIPELINES_BC_SP','GBA_TRANSMISSION_LINES_SP']	
	for v in vL:
		df=gpd.read_file(meta['Paths']['GDB']['LandUse'],layer=v)
		df=df.simplify(100)
		df.to_file(meta['Paths']['GDB']['GDB'] + '\\' + v + '.geojson',driver='GeoJSON')
	return

#%% Major rivers
def GetRiversMajor(meta):
	r=gpd.read_file(meta['Paths']['GDB']['LandCover'],layer='FWA_RIVERS_POLY')
	#ord=np.flip(np.argsort(r['AREA_HA']))
	#L=list(r['GNIS_NAME_1'][ord])
	#L[0:300]
	rmL=['Fraser River','Skeena River','Columbia River','Thompson River','Peace River','Liard River','North Thompson River','South Thompson River', \
	  'Lillooet River']
	rm=r[(np.isin(r['GNIS_NAME_1'],rmL)==True)]
	rm=rm.simplify(100)
	rm.to_file(meta['Paths']['GDB']['GDB'] + '\\rivermajor.geojson',driver='GeoJSON')

	rmL=['Nechako River','Nicola River','Nass River','Quesnel River','Chilcotin River','Liard River','Similkameen River','Fort Nelson River','Eagle River', \
	  'Granby River','Chilko River','Kootenay River']
	rm=r[(np.isin(r['GNIS_NAME_1'],rmL)==True)]
	rm=rm.simplify(100)
	rm.to_file(meta['Paths']['GDB']['GDB'] + '\\riversecond.geojson',driver='GeoJSON')
	return

#%% Rasterize from source
# Only used when the GDB can be opened in python - big files (e.g. VRI) done differently
def RasterizeFromSource(meta,zRef,lNam,vNam):
	if os.path.exists(meta['Paths']['bc1ha'] + '\\' + lNam)==False:
		os.mkdir(meta['Paths']['bc1ha'] + '\\' + lNam)

	ind=np.where( (meta['Geos']['Variable Info']['Layer Name']==lNam) & (meta['Geos']['Variable Info']['Variable Name']==vNam) )[0]

	pthin=meta['Paths']['GDB'][ meta['Geos']['Variable Info']['Geodatabase Name'][ind[0]] ] # fiona.listlayers(pthin)
	df=gpd.read_file(pthin,layer=lNam)

	df=df[df.geometry!=None]
	df=df.reset_index()

	if meta['Geos']['Variable Info']['LUT Required'][ind]=='Yes':
		df=CreateIdForCategoricalVariable(meta,lNam,vNam,df)
		shapes=((geom,value) for geom, value in zip(df['geometry'],df['ID_' + vNam]))
	elif meta['Geos']['Variable Info']['Date Conversion Required'][ind]=='Yes':
		df[vNam + '_Year']=df[vNam].yearnp.zeros(len(df),dtype='int16')
		for i in range(len(df)):
			try:
				df[vNam + '_Year'][i]=df[vNam][i].year
			except:
				continue
		shapes=((geom,value) for geom, value in zip(df['geometry'],df[vNam + '_Year']))
		vNam=vNam + '_Year'
	else:
		shapes=((geom,value) for geom, value in zip(df['geometry'],df[vNam]))

	z0=np.zeros(zRef['Data'].shape,dtype=float)
	burned=features.rasterize(shapes=shapes,fill=0,out=z0,transform=zRef['Transform'])

	z1=copy.deepcopy(zRef)
	z1['Data']=z0.astype(meta['Geos']['Variable Info']['Precision'][ind[0]])
	gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\' + lNam + '\\' + vNam + '.tif')

	return

#%% Get list of variable labels from GDB
def GetVariablesFromGDB(meta,lNam):
	ind=np.where( (meta['Geos']['Variable Info']['Layer Name']==lNam) )[0]
	if ind.size>1:
		ind=ind[0]
	try:
		pthin=meta['Paths']['GDB'][ meta['Geos']['Variable Info']['Geodatabase Name'][ind] ] # fiona.listlayers(pthin)
	except:
		pthin=meta['Paths']['GDB'][ meta['Geos']['Variable Info']['Geodatabase Name'][ind[0]] ] # fiona.listlayers(pthin)
	with fiona.open(pthin,layer=lNam) as source:
		for feat in source:
			prp=dict(feat['properties'].items())
			break
	L=list(prp.keys())
	print(L)
	return prp

#%% Region of interest
def DefineROI(meta,roi,gdf):

	if roi['Type']=='Prov':
		meta['Graphics']['Map']['Legend X']=0.7
		meta['Graphics']['Map']['Show Lakes']='Off'
		meta['Graphics']['Map']['Show Rivers']='Off'
	elif roi['Type']=='LICS':
		meta['Graphics']['Map']['Legend X']=0.75

	roi['crs']=gdf['bc_bound']['gdf'].crs

	#--------------------------------------------------------------------------
	# Raster data
	#--------------------------------------------------------------------------
	zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
	if roi['Type']=='ByTSA':
		tsa=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\FADM_TSA\\TSA_NUMBER_DESCRIPTION.tif')
		tsa=gis.UpdateGridCellsize(tsa,meta['Graphics']['Map']['RGSF'])
		List=[]
		for k in roi['List']:
			List.append(meta['LUT']['FADM_TSA']['TSA_NUMBER_DESCRIPTION'][k])
		roi['grd']=tsa.copy()
		roi['grd']['Data']=np.zeros(tsa['Data'].shape)
		ind=(np.isin(tsa['Data'],List))
		roi['grd']['Data'][ind]=1
		# Define extent based on mask
		xlim=[np.min(tsa['X'][ind])-5000,np.max(tsa['X'][ind])+5000]
		ylim=[np.min(tsa['Y'][ind])-5000,np.max(tsa['Y'][ind])+5000]
		roi['grd']=gis.ClipRasterByXYLimits(roi['grd'],xlim,ylim)
		gc.collect()

	elif roi['Type']=='ByRegDis':
		df=gdf['regdis']['gdf'][np.isin(gdf['regdis']['gdf']['REGIONAL_DISTRICT_NAME'],roi['List'])==True]
		df['ID']=1.0
		df=df[df.geometry!=None]
		df=df.reset_index()
		shapes=((geom,value) for geom, value in zip(df.geometry,df.ID))
		z0=np.zeros(zRef['Data'].shape,dtype=float)
		burned=features.rasterize(shapes=shapes,fill=0,out=z0,transform=zRef['Transform'])
		roi['grd']=copy.deepcopy(zRef)
		roi['grd']['Data']=burned.astype('int8')
		roi['grd']['Data'][(zRef['Data']==0)]=0
		ind=np.where(burned>0)
		xlim=[np.min(zRef['X'][ind])-5000,np.max(zRef['X'][ind])+5000]
		ylim=[np.min(zRef['Y'][ind])-5000,np.max(zRef['Y'][ind])+5000]
		# Clip mask
		roi['grd']=gis.ClipRasterByXYLimits(roi['grd'],xlim,ylim)
		gc.collect()

	elif roi['Type']=='ByWatershed':
		pthin=meta['Paths']['GDB']['GDB'] + '\\WatershedsNamed.gdb'
		df=gpd.read_file(pthin,layer='FWA_NAMED_WATERSHEDS_POLY')
		df=df[np.isin(df['GNIS_ID'],roi['List'])==True]
		df=df.reset_index()
		#roi['gdf']['wshed']=df
		df['ID']=1.0
		shapes=((geom,value) for geom, value in zip(df.geometry,df.ID))
		z0=np.zeros(zRef['Data'].shape,dtype=float)
		burned=features.rasterize(shapes=shapes,fill=0,out=z0,transform=zRef['Transform'])
		roi['grd']=copy.deepcopy(zRef)
		roi['grd']['Data']=burned.astype('int8')
		roi['grd']['Data'][(zRef['Data']==0)]=0

		# Save
		gis.SaveGeoTiff(roi['grd'],r'C:\Data\BC1ha\Masks' + '\\' + roi['Name'] + '.tif')

		ind=np.where(burned>0)
		xlim=[np.min(zRef['X'][ind])-5000,np.max(zRef['X'][ind])+5000]
		ylim=[np.min(zRef['Y'][ind])-5000,np.max(zRef['Y'][ind])+5000]
		# Clip mask
		roi['grd']=gis.ClipRasterByXYLimits(roi['grd'],xlim,ylim)
		gc.collect()

	elif roi['Type']=='LICS':
		srs=gis.ImportSRSs()
		if 'Centre' in roi.keys():
			xc,yc=srs['Proj']['BC1ha'](roi['Centre'][0],roi['Centre'][1])
			xc=np.round(xc/100)*100
			yc=np.round(yc/100)*100
			xlim=[xc-roi['Radius'],xc+roi['Radius']]
			ylim=[yc-roi['Radius'],yc+roi['Radius']]
		else:
			ll=srs['Proj']['BC1ha'](roi['Lower Left'][0],roi['Lower Left'][1])
			ur=srs['Proj']['BC1ha'](roi['Upper Right'][0],roi['Upper Right'][1])
			xlim=[ll[0],ur[0]]
			ylim=[ll[1],ur[1]]
		roi['grd']=gis.ClipRasterByXYLimits(zRef,xlim,ylim)
		roi['grd']['Data']=0*roi['grd']['Data']+1
		roi['grd']['yxrat']=np.diff(ylim)[0]/np.diff(xlim)[0]

	elif roi['Type']=='FromMask':
		pthin=r'C:\Data\BC1ha\Masks' + '\\' + roi['Name']
		zMask=gis.OpenGeoTiff(pthin + '\\Mask.tif')
		roi['grd']=copy.deepcopy(zRef)
		roi['grd']['Data']=zMask['Data'].astype('int8')
		roi['grd']['Data'][(zRef['Data']==0)]=0
		ind=np.where(zMask['Data']>0)
		xlim=[np.min(zRef['X'][ind])-5000,np.max(zRef['X'][ind])+5000]
		ylim=[np.min(zRef['Y'][ind])-5000,np.max(zRef['Y'][ind])+5000]
		# Clip mask
		roi['grd']=gis.ClipRasterByXYLimits(roi['grd'],xlim,ylim)
		gc.collect()

	elif roi['Type']=='Prov':
		roi['grd']=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
		roi['grd']=gis.UpdateGridCellsize(roi['grd'],meta['Graphics']['Map']['RGSF'])

	#--------------------------------------------------------------------------
	# Vector layers
	#--------------------------------------------------------------------------

	roi['gdf']={}
	geom=box(roi['grd']['Extent'][0],roi['grd']['Extent'][2],roi['grd']['Extent'][1],roi['grd']['Extent'][3]) #Takes: box(W, S, E, N)

	roi['gdf']['bound']=gpd.GeoDataFrame({"id":1,"geometry":[geom]})

	roi['gdf']['bound ROI']=gis.DigitizeBinaryMask(roi['grd'])

	if roi['Type']=='ByTSA':
		roi['gdf']['bound within']=gdf['tsa']['gdf'].iloc[np.isin(gdf['tsa']['gdf'].Name,roi['List'])]
	elif roi['Type']=='ByRegDis':
		roi['gdf']['bound within']=gdf['regdis']['gdf'].iloc[ np.isin(gdf['regdis']['gdf']['REGIONAL_DISTRICT_NAME'],roi['List']) ]
	elif roi['Type']=='ByWatershed':
		roi['gdf']['wshed']=df
		roi['gdf']['bound within']=df
	elif roi['Type']=='FromMask':
		roi['gdf']['bound within']=roi['gdf']['bound ROI']
	else:
		roi['gdf']['bound within']=roi['gdf']['bound']

	# Adjust spatial extent and overlay with boundary within bounding box
	for k in gdf.keys():
		#print('Clipping ' + k)

		if k=='road':
			# This takes a long time, only proceed if showing roads
			if meta['Graphics']['Map']['Show Roads']=='Off':
				continue

		# Confine extent to bounding box
		roi['gdf'][k]=gdf[k]['gdf'].cx[roi['grd']['xmin']:roi['grd']['xmax'],roi['grd']['ymin']:roi['grd']['ymax']]
		roi['gdf'][k]=roi['gdf'][k].reset_index(drop=True)

		if (roi['Type']=='ByTSA') | (roi['Type']=='ByRegDis') | (roi['Type']=='ByWatershed') | (roi['Type']=='FromMask'):
			# Confine extent to irregular boundary of ROI
			try:
				# Crashing on TSA
				# This can be really slow for complex irregular shapse (e.g. water licenses in IDF mask)
				roi['gdf'][k]=gpd.overlay(roi['gdf'][k],roi['gdf']['bound ROI'],how='intersection')
			except:
				print(k + ' not working.')

	# Set path for figures
	if roi['Type']=='ByTSA':
		meta['Graphics']['Print Figure Path']=r'C:\Users\rhember\OneDrive - Government of BC\Figures' + '\\TSAs\\' + roi['Name']
	elif roi['Type']=='ByRegDis':
		meta['Graphics']['Print Figure Path']=r'C:\Users\rhember\OneDrive - Government of BC\Figures' + '\\Regional Districts\\' + roi['Name']
	elif roi['Type']=='Prov':
		meta['Graphics']['Print Figure Path']=r'C:\Users\rhember\OneDrive - Government of BC\Figures' + '\\Province\\' + roi['Name']
	elif roi['Type']=='ByWatershed':
		meta['Graphics']['Print Figure Path']=r'C:\Users\rhember\OneDrive - Government of BC\Figures\\Watershed\\' + roi['Name']
	elif roi['Type']=='FromMask':
		meta['Graphics']['Print Figure Path']=r'C:\Users\rhember\OneDrive - Government of BC\Figures\\LICS\\' + roi['Name']

	# Get bounds in geographic coordinates
	srs=gis.ImportSRSs()
	lon0,lat0=pyproj.transform(srs['Proj']['BC1ha'],srs['Proj']['Geographic'],roi['grd']['xmin'],roi['grd']['ymin'])
	lon1,lat1=pyproj.transform(srs['Proj']['BC1ha'],srs['Proj']['Geographic'],roi['grd']['xmax'],roi['grd']['ymax'])
	roi['grd']['Geographic Bounds']=[lon0,lat0,lon1,lat1]

	return meta,roi

#%% Import province-wide basemaps
def Import_GDBs_ProvinceWide(meta):
	flg=0
	if flg==1:
		fiona.listlayers(meta['Paths']['GDB']['LandCover'])
		fiona.listlayers(meta['Paths']['GDB']['LandUse'])
	gdf={}
	gdf['bc_bound']={'gdf':gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='NRC_POLITICAL_BOUNDARIES_1M_SP')}
	gdf['rivers']={'gdf':gpd.read_file(meta['Paths']['GDB']['GDB'] + '\\rivers.geojson')}
	gdf['rivermajor']={'gdf':gpd.read_file(meta['Paths']['GDB']['GDB'] + '\\rivermajor.geojson')}
	gdf['riversecond']={'gdf':gpd.read_file(meta['Paths']['GDB']['GDB'] + '\\riversecond.geojson')}
	gdf['lakes']={'gdf':gpd.read_file(meta['Paths']['GDB']['GDB'] + '\\lakes.geojson')} 
	gdf['tsa']={'gdf':gpd.read_file(meta['Paths']['GDB']['GDB'] + '\\tsa.geojson')}
	gdf['regdis']={'gdf':gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='TA_REGIONAL_DISTRICTS_SVW')}
	gdf['road']={'gdf':gpd.read_file(meta['Paths']['GDB']['GDB'] + '\\road.geojson')}
	gdf['rail']={'gdf':gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='DRP_RAILWAYS_1M_SP')}
	gdf['tpf']={'gdf':gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='GSR_TMBR_PRCSSING_FAC_SV')}
	gdf['fnc']={'gdf':gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='FN_COMMUNITY_LOCATIONS_SP')}
	gdf['hydrol']={'gdf':gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='GBA_TRANSMISSION_LINES_SP')}
	gdf['skih']={'gdf':gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='REC_TENURE_ALPINE_SKI_AREAS_SP')}
	gdf['ogf']={'gdf':gpd.read_file(meta['Paths']['GDB']['GDB'] + '\\ogf.geojson')}
	gdf['ogp']={'gdf':gpd.read_file(meta['Paths']['GDB']['GDB'] + '\\ogp.geojson')}	
	gdf['popp']={'gdf':gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='NRC_POPULATED_PLACES_1M_SP')}
	gdf['cities']={'gdf':gis.ImportCities(r'C:\Data\Cities\Cities.xlsx','GDF')}
	gdf['cities']['gdf'].crs=meta['Geos']['crs']
	gdf['bioe']={'gdf':gis.ImportBioenergy(meta['Paths']['Model']['Parameters'] + '\\Parameters_Facilities.xlsx','GDF')}
	gdf['bioe']['gdf'].crs=meta['Geos']['crs']

	if meta['Graphics']['Vector Import']['Wetland']=='On':
		gdf['wetland']={'gdf':gpd.read_file(meta['Paths']['GDB']['GDB'] + '\\wetlands_digitized1.geojson')}

	if meta['Graphics']['Vector Import']['Water Management']=='On':
		pth=meta['Paths']['GDB']['GDB'] + '\\WaterManagement.gdb'
		gdf['wlic']={'gdf':gpd.read_file(pth,layer='WLS_LICENCE_WITH_PARCELS_ISP')}
		gdf['dams']={'gdf':gpd.read_file(pth,layer='WRIS_DAMS_PUBLIC_SVW')}
		gdf['mmwb']={'gdf':gpd.read_file(pth,layer='FWA_MANMADE_WATERBODIES_POLY')}
		gdf['flooda']={'gdf':gpd.read_file(pth,layer='WLS_NON_TRIM_FLDAREA_LINES_SP')}
		gdf['hydat']={'gdf':gpd.read_file(meta['Paths']['GDB']['GDB'] + '\\streamflow_stations_hydat.geojson')}

		# Streamflow from GRDC
		d=gu.ReadExcel(r'C:\Data\Streamflow\GRDC\GRDC_Stations.xlsx')
		srs=gis.ImportSRSs()
		d['X']=np.zeros(d['long'].size)
		d['Y']=np.zeros(d['lat'].size)
		for i in range(d['lat'].size):
			d['X'][i],d['Y'][i]=srs['Proj']['BC1ha'](d['long'][i],d['lat'][i])
		points=[]
		for i in range(d['X'].size):
			points.append(Point(d['X'][i],d['Y'][i]))
		gdf['streamflow_grdc']={'gdf':gpd.GeoDataFrame({'geometry':points,'grdc_no':d['grdc_no'],'station':d['station'],'river':d['river']})}

	return gdf

#%% Clip geodataframe to ROI

def ClipGDF_ByROI(gdf_in,roi):

	gdf_out=gdf_in.cx[roi['grd']['xmin']:roi['grd']['xmax'],roi['grd']['ymin']:roi['grd']['ymax']]
	gdf_out=gdf_out.reset_index(drop=True)
	gdf_out=gpd.sjoin(gdf_out,roi['gdf']['bound'],how='left')

	# This grouped by may not be necessary - it prevents the file from working in overlays
	#gdf_out=gdf_out.groupby('index_right')

	return gdf_out

#%% Import geodatabases for ROI
def Import_GDB_Over_ROI(meta,roi,vList,**kwargs):
	for nam in vList:
		if nam=='op':
			roi['gdf'][nam]={}
			roi['gdf'][nam]['Path']=meta['Paths']['GDB']['Results']
			roi['gdf'][nam]['Layer']='RSLT_OPENING_SVW'; # fiona.listlayers(op['Path'])
			roi['gdf'][nam]['crs']=roi['crs']
			roi['gdf'][nam]['Keep Geom']='On'
			roi['gdf'][nam]['Select Openings']=np.array([])
			roi['gdf'][nam]['SBC']=np.array([])
			roi['gdf'][nam]['FSC']=np.array([])
			roi['gdf'][nam]['gdf']=qgdb.Query_Openings(roi['gdf'][nam],roi)
		elif nam=='fcinv':
			roi['gdf'][nam]={}
			roi['gdf'][nam]['Path']=meta['Paths']['GDB']['Results']
			roi['gdf'][nam]['Layer']='RSLT_FOREST_COVER_INV_SVW'; # fiona.listlayers(fcinv['Path'])
			roi['gdf'][nam]['crs']=roi['crs']
			roi['gdf'][nam]['Keep Geom']='On'
			roi['gdf'][nam]['Select Openings']=np.array([])
			roi['gdf'][nam]['SBC']=np.array([])
			roi['gdf'][nam]['FSC']=np.array([])
			roi['gdf'][nam]['gdf']=qgdb.Query_Openings(roi['gdf'][nam],roi)
		elif nam=='fcres':
			roi['gdf'][nam]={}
			roi['gdf'][nam]['Path']=meta['Paths']['GDB']['Results']
			roi['gdf'][nam]['Layer']='RSLT_FOREST_COVER_RESERVE_SVW'; # fiona.listlayers(fcres['Path'])
			roi['gdf'][nam]['crs']=roi['crs']
			roi['gdf'][nam]['Keep Geom']='On'
			roi['gdf'][nam]['Select Openings']=np.array([])
			roi['gdf'][nam]['SBC']=np.array([])
			roi['gdf'][nam]['STC']=np.array([])
			roi['gdf'][nam]['SMC']=np.array([])
			roi['gdf'][nam]['SOC1']=np.array([])
			roi['gdf'][nam]['FSC']=np.array([])
			roi['gdf'][nam]['gdf']=qgdb.Query_Openings(roi['gdf'][nam],roi)
		# Import atu within ROI
		elif nam=='pl':
			roi['gdf'][nam]={}
			roi['gdf'][nam]['Path']=meta['Paths']['GDB']['Results']
			roi['gdf'][nam]['Layer']='RSLT_ACTIVITY_TREATMENT_SVW'; # fiona.listlayers(atu['Path'])
			roi['gdf'][nam]['crs']=roi['crs']
			roi['gdf'][nam]['Keep Geom']='On'
			roi['gdf'][nam]['Select Openings']=np.array([])
			roi['gdf'][nam]['SBC']=np.array(['PL'])
			roi['gdf'][nam]['FSC']=np.array([])
			roi['gdf'][nam]['gdf']=qgdb.Query_Openings(roi['gdf'][nam],roi)
		elif nam=='wf':
			roi['gdf'][nam]={}
			roi['gdf'][nam]['Path']=meta['Paths']['GDB']['Disturbance']
			roi['gdf'][nam]['Layer']='PROT_HISTORICAL_FIRE_POLYS_SP'; # fiona.listlayers(wf['Path'])
			roi['gdf'][nam]['crs']=roi['crs']
			roi['gdf'][nam]['Keep Geom']='On'
			roi['gdf'][nam]['Select Openings']=np.array([])
			roi['gdf'][nam]['SBC']=np.array([])
			roi['gdf'][nam]['FSC']=np.array([])
			roi['gdf'][nam]['Year Start']=1800
			roi['gdf'][nam]['Year End']=2100
			roi['gdf'][nam]['gdf']=qgdb.Query_Wildfire(roi['gdf'][nam],roi)
		elif nam=='cc':
			roi['gdf'][nam]={}
			roi['gdf'][nam]['Path']=meta['Paths']['GDB']['Disturbance']
			roi['gdf'][nam]['Layer']='VEG_CONSOLIDATED_CUT_BLOCKS_SP'; # fiona.listlayers(cc['Path'])
			roi['gdf'][nam]['crs']=roi['crs']
			roi['gdf'][nam]['Keep Geom']='On'
			roi['gdf'][nam]['Select Openings']=np.array([])
			roi['gdf'][nam]['SBC']=np.array([])
			roi['gdf'][nam]['FSC']=np.array([])
			roi['gdf'][nam]['Year Start']=1900
			roi['gdf'][nam]['Year End']=2025
			roi['gdf'][nam]['gdf']=qgdb.Query_ConsolidatedCutblocks(roi['gdf'][nam],roi)
		elif nam=='vri':
			roi['gdf'][nam]={}
			roi['gdf'][nam]['Path']=meta['Paths']['GDB']['VRI']
			roi['gdf'][nam]['Layer']='VEG_COMP_LYR_R1_POLY';
			roi['gdf'][nam]['crs']=roi['crs']
			roi['gdf'][nam]['Keep Geom']='On'
			roi['gdf'][nam]['Select Openings']=np.array([])
			roi['gdf'][nam]['SBC']=np.array([])
			roi['gdf'][nam]['FSC']=np.array([])
			roi['gdf'][nam]['Year Start']=1900
			roi['gdf'][nam]['Year End']=2025
			roi['gdf'][nam]['gdf']=qgdb.Query_VRI(roi['gdf'][nam],roi)
		elif nam=='ogsr':
			roi['gdf'][nam]={}
			roi['gdf'][nam]['Path']=meta['Paths']['GDB']['LandUse']
			roi['gdf'][nam]['Layer']='OGSR_TAP_PRIORITY_DEF_AREA_SP';# fiona.listlayers(roi['gdf'][nam]['Path'])
			roi['gdf'][nam]['crs']=roi['crs']
			roi['gdf'][nam]['Keep Geom']='On'
			roi['gdf'][nam]['Select Openings']=np.array([])
			roi['gdf'][nam]['SBC']=np.array([])
			roi['gdf'][nam]['FSC']=np.array([])
			roi['gdf'][nam]['Year Start']=1900
			roi['gdf'][nam]['Year End']=2025
			roi['gdf'][nam]['gdf']=qgdb.Query_Openings(roi['gdf'][nam],roi)
		elif nam=='wshed':
			pthin=meta['Paths']['GDB']['GDB'] + '\\WatershedsNamed.gdb'
			df=gpd.read_file(pthin,layer='FWA_NAMED_WATERSHEDS_POLY')
			df=df[np.isin(df['STREAM_ORDER'],kwargs['OrderList'])==True]
			df=df.reset_index()
			roi['gdf'][nam]={}
			roi['gdf'][nam]['gdf']=df.cx[roi['grd']['xmin']:roi['grd']['xmax'],roi['grd']['ymin']:roi['grd']['ymax']]
		elif (nam=='wshed4') | (nam=='wshed5') | (nam=='wshed6') | (nam=='wshed7') | (nam=='wshed8') | (nam=='wshed9'):
			if nam not in roi['gdf'].keys():
				pthin=meta['Paths']['GDB']['GDB'] + '\\WatershedsNamed\\WatershedsNamed_ord' + str(nam[-1]) + '.geojson'
				df=gpd.read_file(pthin)
				#df=df[np.isin(df['STREAM_ORDER'],kwargs['OrderList'])==True]
				#df=df.reset_index()
				roi['gdf'][nam]={}
				roi['gdf'][nam]['gdf']=df.cx[roi['grd']['xmin']:roi['grd']['xmax'],roi['grd']['ymin']:roi['grd']['ymax']]
		elif nam=='road_ften':
			pthin=meta['Paths']['GDB']['GDB'] + '\\FTEN_ROAD_SEGMENT_LINES_SVW.geojson'
			df=gpd.read_file(pthin)
			roi['gdf'][nam]={}
			roi['gdf'][nam]['gdf']=df.cx[roi['grd']['xmin']:roi['grd']['xmax'],roi['grd']['ymin']:roi['grd']['ymax']]
		else:
			pass

		if (roi['Type']=='ByTSA') | (roi['Type']=='ByWatershed') | (roi['Type']=='FromMask'):
			#if nam not in roi['gdf'].keys():
			try:
				roi['gdf'][nam]['gdf']=gpd.overlay(roi['gdf'][nam]['gdf'],roi['gdf']['bound ROI'],how='intersection')
			except:
				# This will crash if it has already bee
				pass

	return roi

#%% Planting (Non-obligation by project type)
# Important layer for NOSE modelling project
def DerivePlantingMaskNonOb(meta,zRef):
	zMask=copy.deepcopy(zRef)
	zMask['Data']=np.zeros(zRef['Data'].shape,dtype='int16')
	for iEY in range(6):
		zFSC=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(iEY+1) + '_SILV_FUND_SOURCE_CODE.tif')
		ind=np.where( (np.isin(zFSC['Data'],meta['Param']['Raw']['FSC']['NO List ID'])==True) )
		zMask['Data'][ind]=1
	gis.SaveGeoTiff(zMask,meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_NonOb_MaskAll.tif')
	return

#%%
def RasterizeSilviculture(meta,sbc,stc,smc,soc1,vNam,YearLast):

	flg=0
	if flg==1:
		YearLast=2023
		sbc=np.array(['SU']);stc=np.array(['RA','RE','RG','SR','SU']); smc=np.array(['PLOT']);soc1=np.array([]);vNam='SU-PLOT'
		sbc=np.array(['DN']);stc=np.array([]); smc=np.array([]);soc1=np.array([]);vNam='DN'
		sbc=np.array(['FE']);stc=np.array(['CA']); smc=np.array([]);soc1=np.array([]);vNam='FE-CA'
		sbc=np.array(['SP']);stc=np.array([]); smc=np.array(['CABLE','GUARD','HARV','MDOWN','PUSH']);soc1=np.array([]);vNam='SP-KD'
		sbc=meta,np.array(['SP']),np.array(['BU']),np.array(['BROAD']),np.array([]),np.array(['SP-BU-BROAD'])
		sbc=np.array(['SP']);stc=np.array(['ME']);smc=np.array(['DISC','EXCAV','MOUND','LRIP','RRIP','WING']);soc1=np.array([]);vNam='SP-Rip'
		sbc=np.array(['SP']);stc=np.array([]);smc=np.array(['BROAD']);soc1=np.array([]);vNam=np.array(['SP-BROAD'])

	tv=np.arange(1960,YearLast+1,1)
	zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])

	#--------------------------------------------------------------------------
	# Start with areas that reported spatial from RESULTS (takes 15 min)
	#--------------------------------------------------------------------------
	flg=1
	if flg==1:
		t0=time.time()
		ats={}
		ats['Path']=meta['Paths']['GDB']['Results']
		ats['Layer']='RSLT_ACTIVITY_TREATMENT_SVW'; # fiona.listlayers(ats['Path'])
		ats['crs']=meta['Geos']['crs']
		ats['Keep Geom']='On'
		ats['Select Openings']=np.array([])
		ats['SBC']=sbc
		ats['STC']=stc
		ats['SMC']=smc
		ats['FSC']=np.array([])
		ats['SOC1']=soc1
		ats['ROI']=[]
		ats['gdf']=qgdb.Query_Openings(ats,[])
		ats['gdf']=ats['gdf'][ats['gdf'].geometry!=None]
		ats['gdf']=ats['gdf'].reset_index()
		ats['gdf']['Year']=np.zeros(len(ats['gdf']))
		for i in range(ats['gdf']['Year'].size):
			ats['gdf']['Year'][i]=int(ats['gdf']['ATU_COMPLETION_DATE'][i][0:4])
		ats['gdf']['Month']=np.zeros(len(ats['gdf']))
		for i in range(ats['gdf']['Month'].size):
			ats['gdf']['Month'][i]=int(ats['gdf']['ATU_COMPLETION_DATE'][i][5:7])
		ats['gdf']['Day']=np.zeros(len(ats['gdf']))
		for i in range(ats['gdf']['Day'].size):
			ats['gdf']['Day'][i]=int(ats['gdf']['ATU_COMPLETION_DATE'][i][8:10])
		ats['gdf']['DOY']=np.zeros(len(ats['gdf']))
		for i in range(ats['gdf']['Day'].size):
			ats['gdf']['DOY'][i]=datetime(int(ats['gdf']['Year'][i]),int(ats['gdf']['Month'][i]),int(ats['gdf']['Day'][i])).timetuple().tm_yday

		AreaPlanted=ats['gdf']['ACTUAL_TREATMENT_AREA']
		ats['gdf']=u1ha.CreateIdForCategoricalVariable(meta,'RSLT_ACTIVITY_TREATMENT_SVW','SILV_FUND_SOURCE_CODE',ats['gdf'])
		ats['gdf']=u1ha.CreateIdForCategoricalVariable(meta,'RSLT_ACTIVITY_TREATMENT_SVW','SILV_BASE_CODE',ats['gdf'])
		ats['gdf']=u1ha.CreateIdForCategoricalVariable(meta,'RSLT_ACTIVITY_TREATMENT_SVW','SILV_TECHNIQUE_CODE',ats['gdf'])
		ats['gdf']=u1ha.CreateIdForCategoricalVariable(meta,'RSLT_ACTIVITY_TREATMENT_SVW','SILV_METHOD_CODE',ats['gdf'])
		ats['gdf']=ats['gdf'][ (ats['gdf']['RESULTS_IND']=='Y') & (ats['gdf']['SILV_METHOD_CODE']!='LAYOT') ].copy()
		ats['gdf']=ats['gdf'][ats['gdf'].geometry!=None]
		ats['gdf']=ats['gdf'].reset_index()

	#--------------------------------------------------------------------------
	# Add areas where FC is artificial
	#--------------------------------------------------------------------------
	flg=1
	if flg==1:
		at={}
		at['Path']=meta['Paths']['GDB']['Results']
		at['Layer']='RSLT_ACTIVITY_TREATMENT_SVW'; # fiona.listlayers(at['Path'])
		at['crs']=meta['Geos']['crs']
		at['Keep Geom']='Off'
		at['Select Openings']=np.array([])
		at['SBC']=sbc
		at['STC']=stc
		at['SMC']=smc
		at['FSC']=np.array([])
		at['SOC1']=soc1
		at['ROI']=[]
		at['gdf']=qgdb.Query_Openings(at,[])
		at['gdf']=u1ha.CreateIdForCategoricalVariable(meta,'RSLT_ACTIVITY_TREATMENT_SVW','SILV_FUND_SOURCE_CODE',at['gdf'])
		at['gdf']=u1ha.CreateIdForCategoricalVariable(meta,'RSLT_ACTIVITY_TREATMENT_SVW','SILV_BASE_CODE',at['gdf'])
		at['gdf']=u1ha.CreateIdForCategoricalVariable(meta,'RSLT_ACTIVITY_TREATMENT_SVW','SILV_TECHNIQUE_CODE',at['gdf'])
		at['gdf']=u1ha.CreateIdForCategoricalVariable(meta,'RSLT_ACTIVITY_TREATMENT_SVW','SILV_METHOD_CODE',at['gdf'])
		# Make sure to remove entries that we know did not occur (planned or layout)
		ikp=np.where(  (at['gdf']['RESULTS_IND']=='Y') & (at['gdf']['SILV_METHOD_CODE']!='LAYOT') )[0]
		for k in at['gdf'].keys():
			at['gdf'][k]=at['gdf'][k][ikp]
		at['gdf']['Year']=np.zeros(at['gdf']['ACTIVITY_TREATMENT_UNIT_ID'].size)
		for i in range(at['gdf']['Year'].size):
			at['gdf']['Year'][i]=int(at['gdf']['ATU_COMPLETION_DATE'][i][0:4])
		at['gdf']['Month']=np.zeros(at['gdf']['ACTIVITY_TREATMENT_UNIT_ID'].size)
		for i in range(at['gdf']['Month'].size):
			at['gdf']['Month'][i]=int(at['gdf']['ATU_COMPLETION_DATE'][i][5:7])
		at['gdf']['Day']=np.zeros(at['gdf']['ACTIVITY_TREATMENT_UNIT_ID'].size)
		for i in range(at['gdf']['Day'].size):
			at['gdf']['Day'][i]=int(at['gdf']['ATU_COMPLETION_DATE'][i][8:10])
		at['gdf']['DOY']=np.zeros(at['gdf']['ACTIVITY_TREATMENT_UNIT_ID'].size)
		for i in range(at['gdf']['Day'].size):
			try:
				at['gdf']['DOY'][i]=datetime(int(at['gdf']['Year'][i]),int(at['gdf']['Month'][i]),int(at['gdf']['Day'][i])).timetuple().tm_yday
			except:
				pass
	#np.sum(ats['gdf']['ACTUAL_TREATMENT_AREA'])
	#ind=np.where(at['gdf']['SILV_METHOD_CODE']=='BROAD')[0]
	#np.sum(at['gdf']['ACTUAL_TREATMENT_AREA'])
	#ind.size

	#--------------------------------------------------------------------------
	# Import opening ID with spatial
	#--------------------------------------------------------------------------

	zOP1=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_OPENING_SVW\\OPENING_ID.tif')['Data']
	zOP2=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_OPENING_SVW\\OPENING_ID_2.tif')['Data']

	zFC_OID=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_FOREST_COVER_INV_SVW\\OPENING_ID.tif')['Data']
	zFC_STC=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_FOREST_COVER_INV_SVW\\ID_STOCKING_TYPE_CODE.tif')['Data']

	zVRI_OID=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\OPENING_ID.tif')['Data']

	# Reduce the size of rasters
	indOP1=np.where( (zOP1!=0) )
	zOP1s=zOP1[indOP1]
	indOP2=np.where( (zOP2!=0) )
	zOP2s=zOP2[indOP2]

	indFC=np.where( (zFC_OID!=0) )
	zFC_OIDs=zFC_OID[indFC]
	zFC_STCs=zFC_STC[indFC]

	indVRI=np.where( (zVRI_OID!=0) )
	zVRI_OIDs=zVRI_OID[indVRI]

	# Unique indices to Opening ID
	uOP1s=gu.IndicesFromUniqueArrayValues(zOP1s)
	uOP2s=gu.IndicesFromUniqueArrayValues(zOP2s)
	uFCs=gu.IndicesFromUniqueArrayValues(zFC_OIDs)
	uVRIs=gu.IndicesFromUniqueArrayValues(zVRI_OIDs)

	# Index to activity and year
	dP={}
	for iT in range(tv.size):
		dP[tv[iT]]={}
		for iS in range(4):
			dP[tv[iT]][iS]={}
			dP[tv[iT]][iS]={'IndexToGrid':np.array([],dtype=int),
					   'DOY':np.array([]),
					   'Source FC':np.array([]),
					   'ID_SILV_FUND_SOURCE_CODE':np.array([]),
					   'ACTIVITY_TREATMENT_UNIT_ID':np.array([]),
					   'ID_SILV_BASE_CODE':np.array([]),
					   'ID_SILV_TECHNIQUE_CODE':np.array([]),
					   'ID_SILV_METHOD_CODE':np.array([])}

	N_MissingArea=0 # 2 entries with missing area
	for iAT in range(at['gdf']['Year'].size):
		print(iAT)
		Year=at['gdf']['Year'][iAT].astype(int)
		if (Year<tv[0]) | (Year>tv[-1]):
			continue
		DOY=at['gdf']['DOY'][iAT].astype(int)
		ID=at['gdf']['OPENING_ID'][iAT]
		FSC=at['gdf']['ID_SILV_FUND_SOURCE_CODE'][iAT]
		ATUID=at['gdf']['ACTIVITY_TREATMENT_UNIT_ID'][iAT]
		SBC=at['gdf']['ID_SILV_BASE_CODE'][iAT]
		STC=at['gdf']['ID_SILV_TECHNIQUE_CODE'][iAT]
		SMC=at['gdf']['ID_SILV_METHOD_CODE'][iAT]
		A_Planted=at['gdf']['ACTUAL_TREATMENT_AREA'][iAT]
		NumTreesPlanted=at['gdf']['ACTUAL_PLANTED_NUMBER'][iAT]
		SPH_Planted=NumTreesPlanted/A_Planted
		if np.isnan(A_Planted)==True:
			N_MissingArea=N_MissingArea+1
			continue

		iS=0
		flg=1
		try:
			indArt=np.where(zFC_STCs[uFCs[ID]]==meta['LUT']['RSLT_FOREST_COVER_INV_SVW']['STOCKING_TYPE_CODE']['ART'])[0]
		except:
			flg=0
		if flg==1:
			A_Art=indArt.size
			if A_Art>0:
				fA_fc=np.sum(A_Planted)/A_Art
				if (np.abs(fA_fc-1.0)<0.02):
					ind=uFCs[ID][indArt]
					dP[Year][iS]['IndexToGrid']=np.append(dP[Year][iS]['IndexToGrid'],ind)
					dP[Year][iS]['DOY']=np.append(dP[Year][iS]['DOY'],DOY*np.ones(ind.size))
					dP[Year][iS]['Source FC']=np.append(dP[Year][iS]['Source FC'],1*np.ones(ind.size))
					dP[Year][iS]['ID_SILV_FUND_SOURCE_CODE']=np.append(dP[Year][iS]['ID_SILV_FUND_SOURCE_CODE'],FSC*np.ones(ind.size))
					dP[Year][iS]['ACTIVITY_TREATMENT_UNIT_ID']=np.append(dP[Year][iS]['ACTIVITY_TREATMENT_UNIT_ID'],ATUID*np.ones(ind.size))
					dP[Year][iS]['ID_SILV_BASE_CODE']=np.append(dP[Year][iS]['ID_SILV_BASE_CODE'],SBC*np.ones(ind.size))
					dP[Year][iS]['ID_SILV_TECHNIQUE_CODE']=np.append(dP[Year][iS]['ID_SILV_TECHNIQUE_CODE'],STC*np.ones(ind.size))
					dP[Year][iS]['ID_SILV_METHOD_CODE']=np.append(dP[Year][iS]['ID_SILV_METHOD_CODE'],STC*np.ones(ind.size))
					#print('1')
					continue

		iS=1
		flg=1
		try:
			ind=uOP1s[ID]
		except:
			flg=0
		if flg==1:
			ind2=ind[0:int(np.minimum(ind.size,np.round(A_Planted)))]
			dP[Year][iS]['IndexToGrid']=np.append(dP[Year][iS]['IndexToGrid'],ind2)
			dP[Year][iS]['DOY']=np.append(dP[Year][iS]['DOY'],DOY*np.ones(ind2.size))
			dP[Year][iS]['Source FC']=np.append(dP[Year][iS]['Source FC'],2*np.ones(ind2.size))
			dP[Year][iS]['ID_SILV_FUND_SOURCE_CODE']=np.append(dP[Year][iS]['ID_SILV_FUND_SOURCE_CODE'],FSC*np.ones(ind2.size))
			dP[Year][iS]['ACTIVITY_TREATMENT_UNIT_ID']=np.append(dP[Year][iS]['ACTIVITY_TREATMENT_UNIT_ID'],ATUID*np.ones(ind2.size))
			dP[Year][iS]['ID_SILV_BASE_CODE']=np.append(dP[Year][iS]['ID_SILV_BASE_CODE'],SBC*np.ones(ind2.size))
			dP[Year][iS]['ID_SILV_TECHNIQUE_CODE']=np.append(dP[Year][iS]['ID_SILV_TECHNIQUE_CODE'],STC*np.ones(ind2.size))
			dP[Year][iS]['ID_SILV_METHOD_CODE']=np.append(dP[Year][iS]['ID_SILV_METHOD_CODE'],STC*np.ones(ind2.size))
			#print('2')
			continue

		iS=2
		flg=1
		try:
			ind=uOP2s[ID]
		except:
			flg=0
		if flg==1:
			ind2=ind[0:int(np.minimum(ind.size,np.round(A_Planted)))]
			dP[Year][iS]['IndexToGrid']=np.append(dP[Year][iS]['IndexToGrid'],ind2)
			dP[Year][iS]['DOY']=np.append(dP[Year][iS]['DOY'],DOY*np.ones(ind2.size))
			dP[Year][iS]['Source FC']=np.append(dP[Year][iS]['Source FC'],3*np.ones(ind2.size))
			dP[Year][iS]['ID_SILV_FUND_SOURCE_CODE']=np.append(dP[Year][iS]['ID_SILV_FUND_SOURCE_CODE'],FSC*np.ones(ind2.size))
			dP[Year][iS]['ACTIVITY_TREATMENT_UNIT_ID']=np.append(dP[Year][iS]['ACTIVITY_TREATMENT_UNIT_ID'],ATUID*np.ones(ind2.size))
			dP[Year][iS]['ID_SILV_BASE_CODE']=np.append(dP[Year][iS]['ID_SILV_BASE_CODE'],SBC*np.ones(ind2.size))
			dP[Year][iS]['ID_SILV_TECHNIQUE_CODE']=np.append(dP[Year][iS]['ID_SILV_TECHNIQUE_CODE'],STC*np.ones(ind2.size))
			dP[Year][iS]['ID_SILV_METHOD_CODE']=np.append(dP[Year][iS]['ID_SILV_METHOD_CODE'],STC*np.ones(ind2.size))
			#print('3')
			continue
		iS=3
		flg=1
		try:
			ind=uVRIs[ID]
			if ind.size==1:
				continue
		except:
			flg=0
		if flg==1:
			ind2=ind[0:int(np.minimum(ind.size,np.round(A_Planted)))]
			dP[Year][iS]['IndexToGrid']=np.append(dP[Year][iS]['IndexToGrid'],ind2)
			dP[Year][iS]['DOY']=np.append(dP[Year][iS]['DOY'],DOY*np.ones(ind2.size))
			dP[Year][iS]['Source FC']=np.append(dP[Year][iS]['Source FC'],4*np.ones(ind2.size))
			dP[Year][iS]['ID_SILV_FUND_SOURCE_CODE']=np.append(dP[Year][iS]['ID_SILV_FUND_SOURCE_CODE'],FSC*np.ones(ind2.size))
			dP[Year][iS]['ACTIVITY_TREATMENT_UNIT_ID']=np.append(dP[Year][iS]['ACTIVITY_TREATMENT_UNIT_ID'],ATUID*np.ones(ind2.size))
			dP[Year][iS]['ID_SILV_BASE_CODE']=np.append(dP[Year][iS]['ID_SILV_BASE_CODE'],SBC*np.ones(ind2.size))
			dP[Year][iS]['ID_SILV_TECHNIQUE_CODE']=np.append(dP[Year][iS]['ID_SILV_TECHNIQUE_CODE'],STC*np.ones(ind2.size))
			dP[Year][iS]['ID_SILV_METHOD_CODE']=np.append(dP[Year][iS]['ID_SILV_METHOD_CODE'],STC*np.ones(ind2.size))
			print('From VRI')
	 
		#print('Missing')
	#gu.opickle(r'C:\Data\BC1ha\RSLT_ACTIVITY_TREATMENT_SVW\dP.pkl',dP)
	#dP=gu.ipickle(r'C:\Data\BC1ha\RSLT_ACTIVITY_TREATMENT_SVW\dP.pkl')

	#--------------------------------------------------------------------------
	# Pack
	#--------------------------------------------------------------------------

	# Initialize rasters
	N_Year=6
	vL=['ID_SILV_FUND_SOURCE_CODE','ID_SILV_BASE_CODE','ID_SILV_TECHNIQUE_CODE','ID_SILV_METHOD_CODE','DOY'] # 'ACTIVITY_TREATMENT_UNIT_ID',
	zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
	zPac={'Year':{},'DOY':{},'ID_SILV_FUND_SOURCE_CODE':{},'ID_SILV_BASE_CODE':{},'ID_SILV_TECHNIQUE_CODE':{},'ID_SILV_METHOD_CODE':{}} # ,'ACTIVITY_TREATMENT_UNIT_ID':{}
	for iY in range(N_Year):
		for k in zPac.keys():
			if k=='ACTIVITY_TREATMENT_UNIT_ID':
				zPac[k][iY+1]=np.zeros(zRef['Data'].shape,dtype='int32')
			else:
				zPac[k][iY+1]=np.zeros(zRef['Data'].shape,dtype='int16')

	for iT in range(tv.size):
		print(tv[iT])
		
		zYr={}
		for k in zPac.keys():
			zYr[k]=np.zeros(zRef['Data'].shape,dtype=float)

		# Add activities without spatial
		iS=0
		iA=indFC[0][dP[tv[iT]][iS]['IndexToGrid']]
		iB=indFC[1][dP[tv[iT]][iS]['IndexToGrid']]
		zYr['DOY'][ iA,iB  ]=dP[tv[iT]][iS]['DOY']
		#zYr['ACTIVITY_TREATMENT_UNIT_ID'][ iA,iB  ]=dP[tv[iT]][iS]['ACTIVITY_TREATMENT_UNIT_ID']
		zYr['ID_SILV_FUND_SOURCE_CODE'][ iA,iB  ]=dP[tv[iT]][iS]['ID_SILV_FUND_SOURCE_CODE']
		zYr['ID_SILV_BASE_CODE'][ iA,iB  ]=dP[tv[iT]][iS]['ID_SILV_BASE_CODE']
		zYr['ID_SILV_TECHNIQUE_CODE'][ iA,iB  ]=dP[tv[iT]][iS]['ID_SILV_TECHNIQUE_CODE']
		zYr['ID_SILV_METHOD_CODE'][ iA,iB  ]=dP[tv[iT]][iS]['ID_SILV_METHOD_CODE']
		
		iS=1
		iA=indOP1[0][dP[tv[iT]][iS]['IndexToGrid']]
		iB=indOP1[1][dP[tv[iT]][iS]['IndexToGrid']]
		zYr['DOY'][ iA,iB  ]=dP[tv[iT]][iS]['DOY']
		#zYr['ACTIVITY_TREATMENT_UNIT_ID'][ iA,iB  ]=dP[tv[iT]][iS]['ACTIVITY_TREATMENT_UNIT_ID']
		zYr['ID_SILV_FUND_SOURCE_CODE'][ iA,iB  ]=dP[tv[iT]][iS]['ID_SILV_FUND_SOURCE_CODE']
		zYr['ID_SILV_BASE_CODE'][ iA,iB  ]=dP[tv[iT]][iS]['ID_SILV_BASE_CODE']
		zYr['ID_SILV_TECHNIQUE_CODE'][ iA,iB  ]=dP[tv[iT]][iS]['ID_SILV_TECHNIQUE_CODE']
		zYr['ID_SILV_METHOD_CODE'][ iA,iB  ]=dP[tv[iT]][iS]['ID_SILV_METHOD_CODE']

		iS=2
		iA=indOP2[0][dP[tv[iT]][iS]['IndexToGrid']]
		iB=indOP2[1][dP[tv[iT]][iS]['IndexToGrid']]
		zYr['DOY'][ iA,iB  ]=dP[tv[iT]][iS]['DOY']
		#zYr['ACTIVITY_TREATMENT_UNIT_ID'][ iA,iB  ]=dP[tv[iT]][iS]['ACTIVITY_TREATMENT_UNIT_ID']
		zYr['ID_SILV_FUND_SOURCE_CODE'][ iA,iB  ]=dP[tv[iT]][iS]['ID_SILV_FUND_SOURCE_CODE']
		zYr['ID_SILV_BASE_CODE'][ iA,iB  ]=dP[tv[iT]][iS]['ID_SILV_BASE_CODE']
		zYr['ID_SILV_TECHNIQUE_CODE'][ iA,iB  ]=dP[tv[iT]][iS]['ID_SILV_TECHNIQUE_CODE']
		zYr['ID_SILV_METHOD_CODE'][ iA,iB  ]=dP[tv[iT]][iS]['ID_SILV_METHOD_CODE']
		
		# Add activity layer with spatial
		ats0=ats['gdf'][ (ats['gdf']['Year']==tv[iT]) ].copy()
		ats0=ats0[ats0.geometry!=None]; #ats0=ats0.reset_index()
		if len(ats0)>0:
			for v in vL:
				shapes=((geom,value) for geom, value in zip(ats0['geometry'],ats0[v]))
				burned=features.rasterize(shapes=shapes,fill=0,out=zYr[v],transform=zRef['Transform'])

		# Populate packed grids
		ind=np.where( (zPac['Year'][1]==0) & (zYr['ID_SILV_FUND_SOURCE_CODE']!=0) ) # (zCounter<=1) & 
		zPac['Year'][1][ind]=tv[iT]
		for k in zYr.keys():
			if k=='Year':
				continue
			zPac[k][1][ind]=zYr[k][ind]

		ind=np.where( (zPac['Year'][1]!=0) & (zPac['Year'][1]!=tv[iT]) & (zPac['Year'][2]==0) & (zYr['ID_SILV_FUND_SOURCE_CODE']!=0) )
		zPac['Year'][2][ind]=tv[iT]
		for k in zYr.keys():
			if k=='Year':
				continue
			zPac[k][2][ind]=zYr[k][ind]

		ind=np.where( (zPac['Year'][2]!=0) & (zPac['Year'][2]!=tv[iT]) & (zPac['Year'][3]==0) & (zYr['ID_SILV_FUND_SOURCE_CODE']!=0) )
		zPac['Year'][3][ind]=tv[iT]
		for k in zYr.keys():
			if k=='Year':
				continue
			zPac[k][3][ind]=zYr[k][ind]

		ind=np.where( (zPac['Year'][3]!=0) & (zPac['Year'][3]!=tv[iT]) & (zPac['Year'][4]==0) & (zYr['ID_SILV_FUND_SOURCE_CODE']!=0) )
		zPac['Year'][4][ind]=tv[iT]
		for k in zYr.keys():
			if k=='Year':
				continue
			zPac[k][4][ind]=zYr[k][ind]

		ind=np.where( (zPac['Year'][4]!=0) & (zPac['Year'][4]!=tv[iT]) & (zPac['Year'][5]==0) & (zYr['ID_SILV_FUND_SOURCE_CODE']!=0) )
		zPac['Year'][5][ind]=tv[iT]
		for k in zYr.keys():
			if k=='Year':
				continue
			zPac[k][5][ind]=zYr[k][ind]

		ind=np.where( (zPac['Year'][5]!=0) & (zPac['Year'][5]!=tv[iT]) & (zYr['ID_SILV_FUND_SOURCE_CODE']!=0) )
		zPac['Year'][6][ind]=tv[iT]
		for k in zYr.keys():
			if k=='Year':
				continue
			zPac[k][6][ind]=zYr[k][ind]

	# Save to file
	for iY in range(N_Year):
		z1=copy.deepcopy(zRef)
		z1['Data']=zPac['Year'][iY+1].astype('int16')
		gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + vNam + '_' + str(iY+1) + '_Year.tif')
		z1=copy.deepcopy(zRef)
		z1['Data']=zPac['DOY'][iY+1].astype('int8')
		gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + vNam + '_' + str(iY+1) + '_DOY.tif')
		#z1=copy.deepcopy(zRef)
		#z1['Data']=zPac['ACTIVITY_TREATMENT_UNIT_ID'][iY+1].astype('int32')
		#gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + vNam + '_' + str(iY+1) + '_ACTIVITY_TREATMENT_UNIT_ID.tif')
		z1=copy.deepcopy(zRef)
		z1['Data']=zPac['ID_SILV_FUND_SOURCE_CODE'][iY+1].astype('int16')
		gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + vNam + '_' + str(iY+1) + '_SILV_FUND_SOURCE_CODE.tif')
		z1=copy.deepcopy(zRef)
		z1['Data']=zPac['ID_SILV_BASE_CODE'][iY+1].astype('int16')
		gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + vNam + '_' + str(iY+1) + '_SILV_BASE_CODE.tif')
		z1=copy.deepcopy(zRef)
		z1['Data']=zPac['ID_SILV_TECHNIQUE_CODE'][iY+1].astype('int16')
		gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + vNam + '_' + str(iY+1) + '_SILV_TECHNIQUE_CODE.tif')
		z1=copy.deepcopy(zRef)
		z1['Data']=zPac['ID_SILV_METHOD_CODE'][iY+1].astype('int16')
		gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + vNam + '_' + str(iY+1) + '_SILV_METHOD_CODE.tif')

	# Mask of occurrence
	z1=copy.deepcopy(zRef)
	z1['Data']=np.zeros(zRef['Data'].shape,dtype='int8')
	for iY in range(N_Year):
		z=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + vNam + '_' + str(iY+1) + '_Year.tif')
		ind=np.where(z['Data']>0)
		z1['Data'][ind]=1
	gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + vNam + '_MaskAll.tif')

	# Count
	z1=copy.deepcopy(zRef)
	z1['Data']=np.zeros(zRef['Data'].shape,dtype='int8')
	for iY in range(N_Year):
		z=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + vNam + '_' + str(iY+1) + '_Year.tif')
		ind=np.where(z['Data']>0)
		z1['Data'][ind]=z1['Data'][ind]+1
	gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + vNam + '_Count.tif')

	# Last year
	z1=copy.deepcopy(zRef)
	z1['Data']=np.zeros(zRef['Data'].shape,dtype='int16')
	for iY in range(N_Year):
		z=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + vNam + '_' + str(iY+1) + '_Year.tif')
		ind=np.where( (z['Data']>z1['Data']) )
		z1['Data'][ind]=z['Data'][ind]
	gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + vNam + '_YearLast.tif')

	return

#lNam='RSLT_ACTIVITY_TREATMENT_SVW'
#vNam='SP-ME'
#tv,N=u1ha.TimeSeriesOccurrenceFromPackedEvents(meta,lNam,vNam,6)

flg=0
if flg==1:
	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(15.5,6)); barw=0.35;
	ax.bar(tv,N/1000,0.8)
	ax.set(xticks=np.arange(0,3000,5),xlabel='Time, years',ylabel='Implementation (hectares x 1000 per year)',xlim=[1960,2023],ylim=[0,40])
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=gp['tickl'])
	plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
	plt.tight_layout()
	gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Site Prep\SE-ME_TimeSeries','png',900)
	
# z=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\FRE_1_SILV_TECHNIQUE_CODE.tif')
# u,N=gu.CountByCategories(z['Data'],'Percent')
# plt.bar(u[1:],N[1:])
# meta['LUT']['RSLT_ACTIVITY_TREATMENT_SVW']['SILV_TECHNIQUE_CODE']

#%% Create BGC Zone / NDT Combination
def DeriveBGCZoneNDTCombo(meta):
	zBGC=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\BEC_ZONE_CODE.tif')
	zNDT=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\ndt.tif')

	n=np.zeros(5)
	for i in range(1,5):
		ind=np.where(zNDT['Data']==i)
		n[i]=ind[0].size/1e6

	uBGC=np.unique(zBGC['Data'][(zBGC['Data']>0) & (zBGC['Data']<255)])
	uNDT=np.unique(zNDT['Data'][zNDT['Data']>0])
	cnt=1
	z=np.zeros(zBGC['Data'].shape,dtype='int16')
	d={'ID':np.zeros(100),'BGC':np.array(['empty' for _ in range(100)],dtype=object),'NDT':np.array(['empty' for _ in range(100)],dtype=object),'BGC-NDT':np.array(['empty' for _ in range(100)],dtype=object),'Area':np.zeros(100)}
	for i in range(uBGC.size):
		for j in range(uNDT.size):
			print(cnt)
			ind=np.where( (zBGC['Data']==uBGC[i]) & (zNDT['Data']==uNDT[j]) )
			z[ind]=cnt
			d['ID'][cnt-1]=cnt
			d['BGC'][cnt-1]=uBGC[i]
			d['NDT'][cnt-1]=uNDT[j]
			d['BGC-NDT'][cnt-1]=u1ha.lut_n2s(lut['bgcz'],uBGC[i])[0] + str(uNDT[j])
			d['Area'][cnt-1]=ind[0].size
			cnt=cnt+1

	ind=np.where(d['BGC']!='empty')[0]
	for k in d.keys():
		d[k]=d[k][ind]

	df=pd.DataFrame(d)
	df.to_excel(meta['Paths']['bc1ha'] + '\\VRI 2023\\lut_bgcz_ndt_combo.xlsx')

	z1=copy.deepcopy(zRef)
	z1['Data']=z
	gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\VRI 2023\\bgcz_ndt_combo.tif')
	return

#%%
def DeriveDistanceFromRoads(meta):
	df=gpd.read_file(meta['Paths']['GDB']['GDB'] + '\\road.geojson')
	df=df[df.geometry!=None]
	df['ID']=np.ones(len(df))
	bwD=1; binD=np.arange(bwD,51,bwD)
	hiD=(binD[-1]+bwD)
	z=hiD*np.ones(zRef['Data'].shape,dtype='int16')
	for iD in range(binD.size):
		print(binD[iD])
		df0=df.copy()
		df0['geometry']=df0.geometry.buffer(1000*binD[iD])
		z0=np.zeros(zRef['Data'].shape,dtype=float)
		shapes=((geom,value) for geom, value in zip(df0['geometry'],df0['ID']))
		burned=features.rasterize(shapes=shapes,fill=0,out=z0,transform=zRef['Transform'])
		z[(burned>0) & (z==hiD)]=binD[iD]
	plt.matshow(z)
	z1=copy.deepcopy(zRef)
	z1['Data']=z
	gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\Management\DistanceFromRoads.tif')	
	return

#%%
def DeriveDistanceFromFacility(meta):	
	df=gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='GSR_TMBR_PRCSSING_FAC_SV')	
	df=df[df.geometry!=None]
	df['ID']=np.ones(len(df))	
	bwD=1; binD=np.arange(bwD,51,bwD)
	hiD=(binD[-1]+bwD)
	z=hiD*np.ones(zRef['Data'].shape,dtype='int16')
	for iD in range(binD.size):
		print(binD[iD])
		df0=df.copy()
		df0['geometry']=df0.geometry.buffer(1000*binD[iD])
		z0=np.zeros(zRef['Data'].shape,dtype=float)
		shapes=((geom,value) for geom, value in zip(df0['geometry'],df0['ID']))
		burned=features.rasterize(shapes=shapes,fill=0,out=z0,transform=zRef['Transform'])
		z[(burned>0) & (z==hiD)]=binD[iD]
	plt.matshow(z)
	z1=copy.deepcopy(zRef)
	z1['Data']=z
	gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\Management\\DistanceFromForestryFacility.tif')
	return

#%%
def DeriveDistanceFromSettlement(meta):
	df=gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='BC_MAJOR_CITIES_POINTS_500M')  
	df=df[df.geometry!=None]	
	df['ID']=np.ones(len(df))
	bwD=1; binD=np.arange(bwD,151,bwD)
	hiD=(binD[-1]+bwD)
	z=hiD*np.ones(zRef['Data'].shape,dtype='int16')
	for iD in range(binD.size):
		print(binD[iD])
		df0=df.copy()
		df0['geometry']=df0.geometry.buffer(1000*binD[iD])
		z0=np.zeros(zRef['Data'].shape,dtype=float)
		shapes=((geom,value) for geom, value in zip(df0['geometry'],df0['ID']))
		burned=features.rasterize(shapes=shapes,fill=0,out=z0,transform=zRef['Transform'])
		z[(burned>0) & (z==hiD)]=binD[iD]
	plt.matshow(z)
	z1=copy.deepcopy(zRef)
	z1['Data']=z
	gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\Management\DistanceFromSettlement.tif')
	return

#%% Derive Upland-wetland Forest Mask for NTEMS and VRI
def DeriveUplandWetlandForstMask(meta):
	vList=['refg','lc_vri_l4','lc_vri_l3','lc_ntems_2019']
	z0=Import_Raster(meta,[],vList)
	
	# NTEMS
	z1=copy.deepcopy(z0['refg'])
	z1['Data']=np.zeros(z0['refg']['Data'].shape,dtype='int8')
	ind=np.where( (np.isin(z0['lc_ntems_2019']['Data'],[meta['LUT']['Derived']['lc_ntems_2019']['Coniferous'],meta['LUT']['Derived']['lc_ntems_2019']['Broadleaf'],meta['LUT']['Derived']['lc_ntems_2019']['Mixedwood']])==True) )
	z1['Data'][ind]=meta['LUT']['Derived']['upwetf']['Upland Forest']
	ind=np.where( (np.isin(z0['lc_ntems_2019']['Data'],[meta['LUT']['Derived']['lc_ntems_2019']['Wetland-treed']])==True) )
	z1['Data'][ind]=meta['LUT']['Derived']['upwetf']['Wetland Forest']
	ind=np.where( (np.isin(z0['lc_ntems_2019']['Data'],[meta['LUT']['Derived']['lc_ntems_2019']['Wetland']])==True) )
	z1['Data'][ind]=meta['LUT']['Derived']['upwetf']['Wetland']
	ind=np.where( (z0['refg']['Data']==1) & (z1['Data']==0) )
	z1['Data'][ind]=meta['LUT']['Derived']['upwetf']['Other Land']
	ind=np.where( (z1['Data']==0) )
	z1['Data'][ind]=meta['LUT']['Derived']['upwetf']['Other Land']+1
	gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\LandCoverUse\\UplandWetlandForest_NTEMS.tif')

	# Digitize wetlands (this takes 1.5 hours)
	z1=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\UplandWetlandForest_NTEMS.tif')
	z2=copy.deepcopy(z1)
	z2['Data']=np.zeros(z1['Data'].shape)
	z2['Data'][z1['Data']==meta['LUT']['Derived']['upwetf']['Wetland']]=1
	z2['Data'][z1['Data']==meta['LUT']['Derived']['upwetf']['Wetland Forest']]=1
	gdf=gis.DigitizeBinaryMask(z2)
	gdf.to_file(meta['Paths']['GDB']['GDB'] + '\\wetlands_digitized.geojson',driver='GeoJSON')
	#bdf_bk=gdf.copy()
	#gdf.geometry=gdf.geometry.simplify(200)
	#gdf.to_file(meta['Paths']['GDB']['GDB'] + '\\wetlands_digitized1.geojson',driver='GeoJSON')
	#gdf.geometry=gdf.geometry.simplify(500)
	#gdf.to_file(meta['Paths']['GDB']['GDB'] + '\\wetlands_digitized2.geojson',driver='GeoJSON')

	# VRI
	z1=copy.deepcopy(z0['refg'])
	z1['Data']=np.zeros(z0['refg']['Data'].shape,dtype='int8')
	ind=np.where( (np.isin(z0['lc_vri_l4']['Data'],[meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_4']['TC'],meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_4']['TM'],meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_4']['TB']])==True) )
	z1['Data'][ind]=meta['LUT']['Derived']['upwetf']['Upland Forest']
	ind=np.where( (z0['lc_vri_l3']['Data']==meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_3']['W']) & (np.isin(z0['lc_vri_l4']['Data'],[meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_4']['TC'],meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_4']['TM'],meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_4']['TB']])==True) )
	z1['Data'][ind]=meta['LUT']['Derived']['upwetf']['Wetland Forest']
	ind=np.where( (z0['lc_vri_l3']['Data']==meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_3']['W']) & (np.isin(z0['lc_vri_l4']['Data'],[meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_4']['BL'],meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_4']['BM'],meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_4']['BY'],meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_4']['HG'],meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_4']['HF'],meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_4']['HE'],meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_4']['SI'],meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_4']['SL'],meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_4']['ST']])==True) )
	z1['Data'][ind]=meta['LUT']['Derived']['upwetf']['Wetland']
	ind=np.where( (z0['refg']['Data']==1) & (z1['Data']==0) )
	z1['Data'][ind]=meta['LUT']['Derived']['upwetf']['Other Land']
	ind=np.where( (z1['Data']==0) )
	z1['Data'][ind]=meta['LUT']['Derived']['upwetf']['Other Land']+1
	gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\LandCoverUse\\UplandWetlandForest_VRI.tif')
	return

#%% Clip rasters to standard grid
# *** This also compresses files that come out of Arc crazy big ***
def ClipToBC1ha(meta):

	# Forest Cover ID
	fin=r'C:\Data\BC1ha\RSLT_FOREST_COVER_INV_SVW\fcid_orig.tif'
	fout=r'C:\Data\BC1ha\RSLT_FOREST_COVER_INV_SVW\fcid.tif'
	gis.ClipToRaster_ByFile(fin,fout,meta['Paths']['bc1ha Ref Grid'])

	# VRI Feature ID
	fin=r'C:\Data\BC1ha\VRI 2023\vri_feaid.tif'
	fout=r'C:\Data\BC1ha\VRI 2023\vri_feaid.tif'
	gis.ClipToRaster_ByFile(fin,fout,meta['Paths']['bc1ha Ref Grid'])

	# VRI Feature ID
	fin=r'C:\Data\BC1ha\VRI 2015\featid.tif'
	fout=r'C:\Data\BC1ha\VRI 2015\vri_feaid.tif'
	gis.ClipToRaster_ByFile(fin,fout,meta['Paths']['bc1ha Ref Grid'])

	# DEM
	fin=r'C:\Data\BC1ha\Terrain\elevation.tif'
	fout=r'C:\Data\BC1ha\Terrain\elevation.tif'
	gis.ClipToRaster_ByFile(fin,fout,meta['Paths']['bc1ha Ref Grid'])

	# Aspect
	fin=r'C:\Data\BC1ha\Terrain\aspect.tif'
	fout=r'C:\Data\BC1ha\Terrain\aspect.tif'
	gis.ClipToRaster_ByFile(fin,fout,meta['Paths']['bc1ha Ref Grid'])

	# Slope
	fin=r'C:\Data\BC1ha\Terrain\slope.tif'
	fout=r'C:\Data\BC1ha\Terrain\slope.tif'
	gis.ClipToRaster_ByFile(fin,fout,meta['Paths']['bc1ha Ref Grid'])

	sL=['At','Ba','Bl','Cw','Ep','Fd','Hm','Hw','Lt','Lw','Pa','Pl','Pw','Py','Sb','Se','Ss','Sx','Sw','Yc']
	for s in sL:
		fin=r'C:\Data\SiteProductivityLayer\Site_Prod_BC_Geotiffs\Site_Prod_' + s + '.tif'
		fout=r'C:\Data\BC1ha\SPL\Site_Prod_' + s + '.tif'
		#gis.ClipToRaster_ByFile(fin,fout,meta['Paths']['bc1ha Ref Grid'])
		z=gis.OpenGeoTiff(fout)
		z['Data']=z['Data'].astype('int16')
		gis.SaveGeoTiff(z,fout)
	return

#%% Import packed event data and generate time series of occurrence
def TimeSeriesOccurrenceFromPackedEvents(meta,lNam,vNam,nPack):
	tv=np.arange(1900,2030,1)
	N=np.zeros(tv.size)
	for i in range(nPack):
		z=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\' + lNam + '\\' + vNam + '_' + str(i+1) + '_Year.tif')['Data']
		u,c=np.unique(z,return_counts=True)
		for j in range(u.size):
			ind=np.where(tv==u[j])[0]
			N[ind]=N[ind]+c[j]
	return tv,N

#%% Get list of rasters
def GetRasterListFromSpreadsheet(path):
	d=gu.ReadExcel(path,sheet_name='Sheet1',skiprows=0)
	vList=[]
	for i in range(d['Name'].size):
		if d['Included'][i]==1:
			vList.append(d['Name'][i])
	return vList

#%%
def DeriveHarvest_SILV_SYSTEM_CODE(meta,zRef):
	# Import opening data (takes 26 min)
	t0=time.time()
	ops={}
	ops['Path']=meta['Paths']['GDB']['Results']
	ops['Layer']='RSLT_OPENING_SVW'; 
	ops['crs']=meta['Geos']['crs']
	ops['Keep Geom']='Off'
	ops['Select Openings']=np.array([])
	ops['SBC']=np.array([])
	ops['STC']=np.array([])
	ops['SMC']=np.array([])
	ops['FSC']=np.array([])
	ops['SOC1']=np.array([])
	ops['ROI']=[]
	ops['gdf']=qgdb.Query_Openings(ops,[])
	print((time.time()-t0)/60)
	
	zOP1=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_OPENING_SVW\\OPENING_ID.tif')['Data']
	zOP2=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_OPENING_SVW\\OPENING_ID_2.tif')['Data']

	z1=copy.deepcopy(zRef)
	z1['Data']=np.zeros(zRef['Data'].shape,dtype='int8')
	for k in meta['LUT']['RSLT_OPENING_SVW']['DENUDATION_1_SILV_SYSTEM_CODE'].keys():
		indOP=np.where(ops['gdf']['DENUDATION_1_SILV_SYSTEM_CODE']==k)[0]
		ind=np.where(np.isin(zOP1,ops['gdf']['OPENING_ID'][indOP])==True)
		z1['Data'][ind]=meta['LUT']['RSLT_OPENING_SVW']['DENUDATION_1_SILV_SYSTEM_CODE'][k]
		ind=np.where(np.isin(zOP2,ops['gdf']['OPENING_ID'][indOP])==True)
		z1['Data'][ind]=meta['LUT']['RSLT_OPENING_SVW']['DENUDATION_1_SILV_SYSTEM_CODE'][k]
		indOP=np.where(ops['gdf']['DENUDATION_2_SILV_SYSTEM_CODE']==k)[0]
		ind=np.where(np.isin(zOP1,ops['gdf']['OPENING_ID'][indOP])==True)
		z1['Data'][ind]=meta['LUT']['RSLT_OPENING_SVW']['DENUDATION_1_SILV_SYSTEM_CODE'][k]
		ind=np.where(np.isin(zOP2,ops['gdf']['OPENING_ID'][indOP])==True)
		z1['Data'][ind]=meta['LUT']['RSLT_OPENING_SVW']['DENUDATION_1_SILV_SYSTEM_CODE'][k]

	#plt.close('all'); plt.matshow(z1['Data'])
	gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\RSLT_OPENING_SVW\\SILV_SYSTEM_CODE.tif')
	return

#%%
def RasterizeHarvest_CC(meta,zRef):
	lNam='VEG_CONSOLIDATED_CUT_BLOCKS_SP'
	vNam='HARVEST_YEAR'

	if os.path.exists(meta['Paths']['bc1ha'] + '\\' + lNam)==False:
		os.mkdir(meta['Paths']['bc1ha'] + '\\' + lNam)

	ind=np.where( (meta['Geos']['Variable Info']['Layer Name']==lNam) & (meta['Geos']['Variable Info']['Variable Name']==vNam) )[0]

	pthin=meta['Paths']['GDB'][ meta['Geos']['Variable Info']['Geodatabase Name'][ind[0]] ]

	df=gpd.read_file(pthin,layer=lNam)
	df=df[df.geometry!=None]
	df=df.reset_index()

	zYearFirst=copy.deepcopy(zRef)
	zYearFirst['Data']=np.zeros(zRef['Data'].shape,dtype='int16')
	
	zYearLast=copy.deepcopy(zRef)
	zYearLast['Data']=np.zeros(zRef['Data'].shape,dtype='int16')

	uYear=df[vNam].unique()
	tv=np.arange(np.min(uYear),np.max(uYear)+1,1)

	for iT in range(tv.size):

		df0=df[df[vNam]==tv[iT]].copy()
		shapes=((geom,value) for geom, value in zip(df0.geometry,df0[vNam]))

		z0=np.zeros(zRef['Data'].shape,dtype=float)
		if len(df0)>0:
			burned=features.rasterize(shapes=shapes,fill=0,out=z0,transform=zRef['Transform'])

		z1=copy.deepcopy(zRef)
		z1['Data']=z0.astype(meta['Geos']['Variable Info']['Precision'][ind[0]])
		gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\' + lNam + '\\' + vNam + '_' + str(tv[iT]) + '.tif')

		# Update by year grid
		zYearLast['Data'][burned>0]=tv[iT]

	# Year of first occurrence
	z1=copy.deepcopy(zRef)
	z1['Data']=zYearLast['Data'].astype('int16')
	gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\' + lNam + '\\' + vNam + '_YearFirst.tif')
	
	# Year of last occurrence
	z1=copy.deepcopy(zRef)
	z1['Data']=zYearLast['Data'].astype('int16')
	gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\' + lNam + '\\' + vNam + '_YearLast.tif')

	# Mask of occurrence
	z1=copy.deepcopy(zRef)
	z1['Data']=np.zeros(zRef['Data'].shape,dtype='int8')
	ind0=np.where(zYearLast['Data']>0)
	z1['Data'][ind0]=1
	gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\' + lNam + '\\' + vNam + '_MaskAll.tif')

	# Pack into smaller number of layers

	# Initialize rasters
	N_Year=6
	z={'Year':{}}
	for iY in range(N_Year):
		z['Year'][iY+1]=np.zeros(zRef['Data'].shape,dtype='int16')

	for iT in range(tv.size):
		print(tv[iT])
		z0=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\' + lNam + '\\' + vNam + '_' + str(tv[iT]) + '.tif')['Data']
		ind=np.where((z['Year'][1]==0) & (z0!=0))
		z['Year'][1][ind]=tv[iT];
		ind=np.where((z['Year'][1]!=0) & (z['Year'][1]!=tv[iT]) & (z['Year'][2]==0) & (z0!=0))
		z['Year'][2][ind]=tv[iT];
		ind=np.where((z['Year'][2]!=0) & (z['Year'][2]!=tv[iT]) & (z['Year'][3]==0) & (z0!=0))
		z['Year'][3][ind]=tv[iT];
		ind=np.where((z['Year'][3]!=0) & (z['Year'][3]!=tv[iT]) & (z['Year'][4]==0) & (z0!=0))
		z['Year'][4][ind]=tv[iT];
		ind=np.where((z['Year'][4]!=0) & (z['Year'][4]!=tv[iT]) & (z['Year'][5]==0) & (z0!=0))
		z['Year'][5][ind]=tv[iT];
		ind=np.where((z['Year'][5]!=0) & (z['Year'][5]!=tv[iT]) & (z0!=0))
		z['Year'][6][ind]=tv[iT];

	for iY in range(N_Year):
		z1=copy.deepcopy(zRef)
		z1['Data']=z['Year'][iY+1].astype('int16')
		gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\' + lNam + '\\' + vNam + '_' + str(iY+1) + '_Year.tif')
	return

#%%
def RasterizeCruisePercentDead(meta):
	dC=gu.ipickle(r'C:\Data\ECAS\Received 2023-04-04\rh_COMP_DATA_CleanCombined.pkl')

	vTM=list(meta['LUT']['FTEN_CUT_BLOCK_POLY_SVW']['TIMBER_MARK'].values())
	kTM=list(meta['LUT']['FTEN_CUT_BLOCK_POLY_SVW']['TIMBER_MARK'].keys())

	zTM=gis.OpenGeoTiff(meta['Paths'] + '\\FTEN_CUT_BLOCK_POLY_SVW\\TIMBER_MARK.tif')

	u=np.unique(dC['PRIMARY_MARK'])
	iZ=gu.IndicesFromUniqueArrayValues(zTM['Data'].flatten())

	z1=copy.deepcopy(zRef)
	z1['Data']=np.zeros(zRef['Data'].shape,dtype='int16').flatten()
	for k in iZ.keys():
		ind1=np.where(vTM==k)[0]
		if ind1.size==0:
			continue
		ind2=np.where( (dC['PRIMARY_MARK']==kTM[ind1[0]]) )[0]
		if ind2.size>0:
			z1['Data'][iZ[k]]=np.nanmean(dC['Pct Dead Net'][ind2])
			#print(np.nanmean(d['Pct Dead Net'][ind2]))

	z1['Data']=np.reshape(z1['Data'],zRef['Data'].shape)
	#plt.matshow(z1['Data'],clim=[0,100])
	#plt.hist(z1['Data'][0::20,0::20].flatten(),np.arange(0,105,5))
	gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\Disturbances\\HarvestPercentDead_FromCruise.tif')
	return

#%%
def DeriveRangeTenureMask(meta):
	zLCC1=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCoverClass1.tif')
	zR=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\rangeten.tif')

	d={1:'Forest with grazing tenure',2:'Forest with haycutting tenure',3:'Forest with no range tenure',4:'Non-forest land',5:'Non land'}

	z1=copy.deepcopy(zRef)
	z1['Data']=np.zeros(zRef['Data'].shape,dtype='int8')
	ind=np.where( (zLCC1['Data']==1) & (zR['Data']==1) | (zLCC1['Data']==1) & (zR['Data']==2) | (zLCC1['Data']==1) & (zR['Data']==3) ); z1['Data'][ind]=1
	ind=np.where( (zLCC1['Data']==1) & (zR['Data']==4) | (zLCC1['Data']==1) & (zR['Data']==5) | (zLCC1['Data']==1) & (zR['Data']==6) ); z1['Data'][ind]=2
	ind=np.where( (zLCC1['Data']==1) & (zR['Data']==0) ); z1['Data'][ind]=3
	ind=np.where( (zLCC1['Data']!=1) ); z1['Data'][ind]=4
	ind=np.where( (zRef['Data']!=1) ); z1['Data'][ind]=5
	gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\LandCoverUse\\rangeten_consol.tif')
	return

#%%
def DeriveCrownLandMask(meta):
	z0=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\f_own.tif')
	z1=copy.deepcopy(zRef)
	z1['Data']=np.zeros(zRef['Data'].shape,dtype='int8')
	ind=np.where( (zRef['Data']==1) & (z0['Data']>=9) )
	z1['Data'][ind]=1
	# plt.close('all'); plt.matshow(z1['Data'])# Confirm that it worked
	gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\LandCoverUse\\CrownForestMask.tif')
	return

#%%
def RasterizeEcozonesOfCanada(meta):
	zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])

	srs=gis.ImportSRSs()
	crs=pyproj.CRS(srs['String']['Geographic'])

	pthin=r'C:\Data\Ecozones\nef_ca_ter_ecozone_v2_2.geojson'
	df=gpd.read_file(pthin)
	df=df[df.geometry!=None]
	df=df.reset_index()
	df.crs=pyproj.CRS(srs['String']['Geographic'])
	df=df.to_crs({'init':'epsg:3005'})
	# Used to create LUT: df.drop(columns='geometry').to_excel(r'C:\Data\Ecozones\table.xlsx')

	shapes=((geom,value) for geom, value in zip(df['geometry'],df['ECOZONE_ID']))
	z0=np.zeros(zRef['Data'].shape,dtype=float)
	burned=features.rasterize(shapes=shapes,fill=0,out=z0,transform=zRef['Transform'])
	z1=copy.deepcopy(zRef)
	z1['Data']=np.zeros(zRef['Data'].shape,dtype='int16')
	ind=np.where( (burned>0) &(zRef['Data']==1) ); z1['Data'][ind]=burned[ind]
	plt.close('all'); plt.matshow(z1['Data'],clim=[0,15])
	gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\Ecozones_Canada\\Ecozones_Canada.tif')
	return


#%%
def DeriveBurnSeverityCompilation(meta):
	zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])

	# Import provincial data and re-classify according to comp1
	zYearH=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VEG_BURN_SEVERITY_SP\\FIRE_YEAR.tif')
	zBSH0=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VEG_BURN_SEVERITY_SP\\BURN_SEVERITY_RATING.tif')
	zBSH=copy.deepcopy(zBSH0)
	for k in meta['LUT']['Derived']['burnsev_comp1'].keys():
		ind=np.where( (zBSH0['Data']==meta['LUT']['VEG_BURN_SEVERITY_SP']['BURN_SEVERITY_RATING'][k]) )
		zBSH['Data'][ind]=meta['LUT']['Derived']['burnsev_comp1'][k]
	#d=gu.CountByCategories(zBSH['Data'][zYearH['Data']>0],'Percent')
	#plt.bar(n,p)

	# Import BC current year data
	zBSC=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VEG_BURN_SEVERITY_SP\\BURN_SEVERITY_RATING_Current_2023.tif')
	#d=gu.CountByCategories(zBSC['Data'][zBSC['Data']>0],'Percent')
	#plt.bar(n,p)

# 	# Import national data, scale and classify according to (https://catalogue.data.gov.bc.ca/dataset/fire-burn-severity-same-year)
# 	zYn=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Disturbances\\BurnSev_Can_Year.tif')
# 	zBn=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Disturbances\\BurnSev_Can_dNBR.tif')
# 	#plt.hist(zBn['Data'][zBn['Data']>0].flatten()[0::50])
# 	zBn2=copy.deepcopy(zBn)
# 	zBn2['Data']=(zBn2['Data'].astype('float')+275)/5
# 	zBn2['Data']=zBn2['Data'].astype('int16')
# 	#plt.hist(zBn2['Data'][zBn2['Data']>0].flatten()[0::50])
# 	zBn3=copy.deepcopy(zRef)
# 	zBn3['Data']=np.zeros(zRef['Data'].shape,'int8')
# 	ind=np.where( (zBn2['Data']<76) ); zBn3['Data'][ind]=meta['LUT']['Derived']['burnsev_comp1']['Unburned']
# 	ind=np.where( (zBn2['Data']>=76) & (zBn2['Data']<110) ); zBn3['Data'][ind]=meta['LUT']['Derived']['burnsev_comp1']['Low']
# 	ind=np.where( (zBn2['Data']>=110) & (zBn2['Data']<187) ); zBn3['Data'][ind]=meta['LUT']['Derived']['burnsev_comp1']['Medium']
# 	ind=np.where( (zBn2['Data']>=187) ); zBn3['Data'][ind]=meta['LUT']['Derived']['burnsev_comp1']['High']
	
	#n,p=gu.CountByCategories(zBn3['Data'][zYn['Data']>0],'Percent')
	#plt.bar(n,p)
	
	# # Use distribution of BC classes to derive classes from dNBR for Canada (exclude unburned)
	# p=gu.CountByCategories(zBSH['Data'][zBSH['Data']>0],'Percent')
	# p_ord=np.array([23.83,19.68,38.49,15.23])
	# p_ord=p_ord/np.sum(p_ord)*100 # normalize without unknown class	
	# # Get percentile values for national dNBR
	# cpn=np.flip(np.percentile(zBn['Data'][zBn['Data']>0],np.cumsum(p_ord)))
	
	zLastY=copy.deepcopy(zRef)
	zLastY['Data']=np.zeros(zRef['Data'].shape,dtype='int16')
	zLastS=copy.deepcopy(zRef)
	zLastS['Data']=5*np.ones(zRef['Data'].shape,dtype='int8')
	
	N_Year=6
	for iY in range(N_Year):
		zY=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\PROT_HISTORICAL_FIRE_POLYS_SP\\WildfireComp1_' + str(iY+1) + '_Year.tif')

		# Initialize
		zBS=copy.deepcopy(zRef)
		zBS['Data']=5*np.ones(zRef['Data'].shape,dtype='int8')
		
# 		# Add national data (converted from dNBR to severity class)
# 		uY=np.unique(zYn['Data'][zYn['Data']>0])
# 		for iU in range(uY.size):	
# 			ind=np.where( (zY['Data']==uY[iU]) & (zYn['Data']==uY[iU]) )
# 			zBS['Data'][ind]=zBn3['Data'][ind]
# 			zLastY['Data'][ind]=uY[iU]
# 			zLastS['Data'][ind]=zBn3['Data'][ind]
			
			# for iC in range(cpn.size):
			#	 ind=np.where( (zY['Data']==uY[iU]) & (zYn['Data']==uY[iU]) & (zBn['Data']<cpn[iC]) )
			#	 zBS['Data'][ind]=list(meta['LUT']['Derived']['burnsev_comp1'].values())[iC]
			#	 zLastY['Data'][ind]=uY[iU]
			#	 zLastS['Data'][ind]=list(meta['LUT']['Derived']['burnsev_comp1'].values())[iC]
		
		# Add BC historical data
		uY=np.unique(zYearH['Data'][zYearH['Data']>0])
		for iU in range(uY.size):
			ind=np.where( (zY['Data']==uY[iU]) & (zYearH['Data']==uY[iU]) )
			zBS['Data'][ind]=zBSH['Data'][ind]
			zLastY['Data'][ind]=uY[iU]
			zLastS['Data'][ind]=zBSH['Data'][ind]
		#plt.matshow(zBS['Data'],clim=[0,4])

		# Add BC current year data
		yr=2023
		ind=np.where( (zBSC['Data']>0) )
		zBS['Data'][ind]=zBSC['Data'][ind]
		zLastY['Data'][ind]=yr
		zLastS['Data'][ind]=zBSC['Data'][ind]
		#plt.matshow(zBS['Data'],clim=[0,4])

		# Save
		gis.SaveGeoTiff(zBS,meta['Paths']['bc1ha'] + '\\PROT_HISTORICAL_FIRE_POLYS_SP\\WildfireComp1_' + str(iY+1) + '_SevClass.tif')

	gis.SaveGeoTiff(zLastY,meta['Paths']['bc1ha'] + '\\PROT_HISTORICAL_FIRE_POLYS_SP\\BurnSevComp_YearLast.tif')
	gis.SaveGeoTiff(zLastS,meta['Paths']['bc1ha'] + '\\PROT_HISTORICAL_FIRE_POLYS_SP\\BurnSevComp_SevClassLast.tif')

	flg=0
	if flg==1:
		# Check that it is working for a specific year
		yr=2023
		lNam='PROT_HISTORICAL_FIRE_POLYS_SP'
		vNam='FIRE_YEAR'
		N_Year=6
		N=np.zeros(5)
		for iY in range(N_Year):
			zY=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\' + lNam + '\\' + vNam + '_' + str(iY+1) + '_Year.tif')
			zS=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\' + lNam + '\\' + vNam + '_' + str(iY+1) + '_SevClass.tif')
			for j in range(5):
				ind=np.where( (zY['Data']==yr) & (zS['Data']==j+1) )[0]
				N[j]=N[j]+ind.size
		plt.bar(np.arange(5),N)
	return
	
#%% 
def DeriveForestCoverReserveComp1(meta):
	zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
	lut1=meta['LUT']['Derived']['reserve_comp1']
	lut0=meta['LUT']['RSLT_FOREST_COVER_RESERVE_SVW']['SILV_RESERVE_CODE']
	#d=gu.ReadExcel(meta['Paths']['bc1ha'] + '\\Results\\SILV_RESERVE_CODE.xlsx')
	zLCC1=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCover_Comp1_2019.tif')
	zH=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VEG_CONSOLIDATED_CUT_BLOCKS_SP\\HARVEST_YEAR_MaskAll.tif')
	zFCR=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_FOREST_COVER_RESERVE_SVW\\SILV_RESERVE_CODE.tif')
	zIBR=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_FOREST_COVER_INV_SVW\\ID_SILV_RESERVE_CODE.tif')
		
	zC=np.maximum(zFCR['Data'],zIBR['Data'])
	u,N=gu.CountByCategories(zC,'Percent')
	
	flg=0
	if flg==1:
		# Look at overlap
		ind1=np.where(zFCR['Data']>0)
		ind2=np.where(zIBR['Data']>0)
		ind3=np.where( (zFCR['Data']>0) & (zIBR['Data']>0) )
		ind4=np.where( (zFCR['Data']>0) & (zIBR['Data']==0) )
		ind5=np.where( (zFCR['Data']==0) & (zIBR['Data']>0) )
	
		print(ind1[0].size)
		print(ind2[0].size)
		print(ind3[0].size)
		print(ind4[0].size)
		print(ind5[0].size)
	
	z1=copy.deepcopy(zRef)
	z1['Data']=lut1['Unharvested forest']*np.ones(zRef['Data'].shape,dtype='int8')
	ind=np.where( (zC==lut0['D']) ); z1['Data'][ind]=lut1['Dispersed']
	ind=np.where( (zC==lut0['G']) ); z1['Data'][ind]=lut1['Group']
	ind=np.where( (zC==lut0['R']) ); z1['Data'][ind]=lut1['Riparian']
	ind=np.where( (zC==lut0['W']) ); z1['Data'][ind]=lut1['Wildlife trees']
	ind=np.where( (zC==lut0['O']) ); z1['Data'][ind]=lut1['Other']
	ind=np.where( (z1['Data']==lut1['Unharvested forest']) & (zH['Data']>0) ); z1['Data'][ind]=lut1['Harvested with no reserves']
	ind=np.where( (zRef['Data']==1) & (zH['Data']==0) & (zLCC1['Data']==meta['LUT']['Derived']['lc_comp1']['Forest']) ); z1['Data'][ind]=lut1['Unharvested forest']
	ind=np.where( (zRef['Data']==1) & (zLCC1['Data']!=meta['LUT']['Derived']['lc_comp1']['Forest']) ); z1['Data'][ind]=lut1['Non-forest']
	gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\Management\\ForestCoverReserveComp1.tif')
	return

#%% Consolidate Tree Density Class
# *** Ensure harvested areas are listed as dense ***
def DeriveTreeDensityClass(meta):

	# 2002 VRI
	flg=0
	if flg==1:
		z=Import_Raster(meta,[],['refg','lc_comp1_2019','lc_vri_l5_02','harv_yr_comp1'])
		z1=z['refg'].copy()
		z1['Data']=np.zeros(z['lc_vri_l5_02']['Data'].shape,dtype='int8')
		ind=np.where( (z['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (z['lc_vri_l5_02']['Data']==meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_5']['SP']) ); z1['Data'][ind]=meta['LUT']['Derived']['tdc']['Sparse']
		ind=np.where( (z['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (z['lc_vri_l5_02']['Data']==meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_5']['OP']) ); z1['Data'][ind]=meta['LUT']['Derived']['tdc']['Open']
		ind=np.where( (z['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (z['lc_vri_l5_02']['Data']==meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_5']['DE']) ); z1['Data'][ind]=meta['LUT']['Derived']['tdc']['Dense']
		#ind=np.where( (z['harv_yr_comp1']['Data']==1) ); z1['Data'][ind]=meta['LUT']['Derived']['tdc']['Dense']
		ind=np.where( (z['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (z1['Data']==0) ); z1['Data'][ind]=meta['LUT']['Derived']['tdc']['Open']
		# plt.matshow(z1,clim=[0,3])
		gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\LandCoverUse\\TreeDensityClass_2002.tif')

	# 2015 VRI
	flg=1
	if flg==1:
		z=Import_Raster(meta,[],['refg','lc_comp1_2019','harv_yr_comp1'])
		zL5=gis.OpenGeoTiff(r'C:\Data\BC1ha\VRI 2015\BCLCS_LEVEL_5.tif')
		z1=z['refg'].copy()
		z1['Data']=np.zeros(z['refg']['Data'].shape,dtype='int8')
		ind=np.where( (z['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (zL5['Data']==meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_5']['SP']) ); z1['Data'][ind]=meta['LUT']['Derived']['tdc']['Sparse']
		ind=np.where( (z['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (zL5['Data']==meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_5']['OP']) ); z1['Data'][ind]=meta['LUT']['Derived']['tdc']['Open']
		ind=np.where( (z['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (zL5['Data']==meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_5']['DE']) ); z1['Data'][ind]=meta['LUT']['Derived']['tdc']['Dense']
		#ind=np.where( (z['harv_yr_comp1']['Data']==1) ); z1['Data'][ind]=meta['LUT']['Derived']['tdc']['Dense']
		ind=np.where( (z['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (z1['Data']==0) ); z1['Data'][ind]=meta['LUT']['Derived']['tdc']['Open']
		# plt.matshow(z1,clim=[0,3])
		gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\LandCoverUse\\TreeDensityClass_VRI2015.tif')

	# # 2023 VRI
	# z=Import_Raster(meta,[],['refg','lc_comp1_2019','lc_vri_l5','harv_yr_comp1'])
		
	# z1=z['refg'].copy()
	# z1['Data']=np.zeros(z['lc_vri_l5']['Data'].shape,dtype='int8')
	# ind=np.where( (z['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (z['lc_vri_l5']['Data']==meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_5']['SP']) ); z1['Data'][ind]=meta['LUT']['Derived']['tdc']['Sparse']
	# ind=np.where( (z['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (z['lc_vri_l5']['Data']==meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_5']['OP']) ); z1['Data'][ind]=meta['LUT']['Derived']['tdc']['Open']
	# ind=np.where( (z['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (z['lc_vri_l5']['Data']==meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_5']['DE']) ); z1['Data'][ind]=meta['LUT']['Derived']['tdc']['Dense']
	# ind=np.where( (z['harv_yr_comp1']['Data']==1) ); z1['Data'][ind]=meta['LUT']['Derived']['tdc']['Dense']
	# ind=np.where( (z['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (z1['Data']==0) ); z1['Data'][ind]=meta['LUT']['Derived']['tdc']['Open']
	# # plt.matshow(z1,clim=[0,3])
	# gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\LandCoverUse\\TreeDensityClass_Current.tif')
	
	# # With shrubs and grasses
	# z1=z['refg'].copy()
	# z1['Data']=np.zeros(z['lc_vri_l5']['Data'].shape,dtype='int8')
	# ind=np.where( (z['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (z['lc_vri_l5']['Data']==meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_5']['SP']) ); z1['Data'][ind]=1
	# ind=np.where( (z['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (z['lc_vri_l5']['Data']==meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_5']['OP']) ); z1['Data'][ind]=2
	# ind=np.where( (z['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (z['lc_vri_l5']['Data']==meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_5']['DE']) ); z1['Data'][ind]=3
	# ind=np.where( (z['harv_yr_comp1']['Data']==1) ); z1['Data'][ind]=3
	# ind=np.where( (z['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (z1['Data']==0) ); z1['Data'][ind]=2
	# ind=np.where( (z['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Shrub']) ); z1['Data'][ind]=4
	# ind=np.where( (z['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Herb']) ); z1['Data'][ind]=5
	# zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
	# ind=np.where( (zRef['Data']==1) & (z1['Data']==0) ); z1['Data'][ind]=6
	# ind=np.where( (zRef['Data']==0) ); z1['Data'][ind]=7
	# gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\LandCoverUse\\TreeDensityClass_WithShrubsGrasses.tif')
	
	# return

#%% Land mask for BC
def GenerateLandMaskBC(meta):
	zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
	df=gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='NRC_POLITICAL_BOUNDARIES_1M_SP')
	df=df[df.geometry!=None]
	df['ID']=1
	shapes=((geom,value) for geom, value in zip(df.geometry,df.ID))
	z=np.zeros(zRef['Data'].shape,dtype=float)
	burned=features.rasterize(shapes=shapes,fill=0,out=z,transform=zRef['Transform'])
	# plt.close('all'); plt.matshow(burned) # Confirm that it worked
	zOut=copy.deepcopy(zRef)
	zOut['Data']=np.zeros(zRef['Data'].shape,dtype='int8')
	zOut['Data'][burned>0]=1
	zOut['Data']=zOut['Data'].astype('int8')
	gis.SaveGeoTiff(zOut,meta['Paths']['bc1ha'] + '\\LandCoverUse\LandMask.tif')

	# Buffer
	zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
	kernel=np.ones((5,5))
	mask=cv2.dilate(zRef['Data'].astype(np.uint8),kernel,iterations=1)
	z1=copy.deepcopy(zRef)
	z1['Data']=mask.astype('int8')
	#plt.matshow(zRef['Data'])
	#plt.matshow(z1['Data'])
	gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\LandCoverUse\LandMaskBuf.tif')

	# Reverse buffer
	kernel=np.ones((3,3))
	mask=cv2.erode(zRef['Data'].astype(np.uint8),kernel,iterations=1)
	z1=copy.deepcopy(zRef)
	z1['Data']=mask.astype('int8')
	plt.matshow(z1['Data'])
	gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\LandCoverUse\LandMaskBufRev.tif')

	# Create 1km version
	f1k=meta['Paths']['bc1ha'] + '\\LandCoverUse\LandMask_1k.tif'
	zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
	gis.SaveGeoTiff(zRef,f1k)
	gis.ResampleRaster(f1k,0.1)
	return

#%% Harvest early reconstruction
def DeriveHarvestEarlyYears(meta):
	
	zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])	
	lut_lc=meta['LUT']['Derived']['lc_comp1']
	lut_lu=meta['LUT']['Derived']['lu_comp1']
	lut_ntems=meta['LUT']['Derived']['lc_ntems_2019']
	lut_lc_vri_l3=meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_3']

	vList=['bgcz','fire_yr','lc_comp1_2019','lu_comp1_2019','harv_yr_comp1','age_vri23','d2road','d2fac']  # ,'d2fac','d2set'
	z0=Import_Raster(meta,[],vList,'Extract Grid')

	# Regional mask
	MaskR=np.zeros(zRef['Data'].shape,dtype='int8')
	ind=np.where( (z0['bgcz']==meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE']['CWH']) | (z0['bgcz']==meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE']['CDF']) )
	MaskR[ind]=1
	ind=np.where( (np.isin(z0['bgcz'],[meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE']['BG'],
									   meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE']['ICH'],
									   meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE']['IDF'],
									   meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE']['MS'],
									   meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE']['PP'],
									   meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE']['SBPS'],
									   meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE']['SBS']])==True) )
	MaskR[ind]=2
	ind=np.where( (np.isin(z0['bgcz'],[meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE']['BWBS']])==True) )
	MaskR[ind]=3

	# Mask area where early harvesting will occur
	dTH_Road_Coast=30
	#dTH_Fac_Coast=150
	dTH_Road_Int=18
	dTH_Fac_Int=18
	dTH_Road_Bor=6
	dTH_Fac_Bor=6
	aTH_Coast=150
	aTH_Interior=95
	MaskEH=np.zeros(zRef['Data'].shape,dtype='int8')	
	ind=np.where( (z0['lc_comp1_2019']==lut_lc['Forest']) & (z0['harv_yr_comp1']==0) & (MaskR==1) & (z0['age_vri']<aTH_Coast) & (z0['fire_yr']==0) & (z0['d2road']<=dTH_Road_Coast) ); MaskEH[ind]=1
	ind=np.where( (z0['lc_comp1_2019']==lut_lc['Forest']) & (z0['harv_yr_comp1']==0) & (MaskR==2) & (z0['age_vri']<aTH_Interior) & (z0['fire_yr']==0) & (z0['d2road']<=dTH_Road_Int) & (z0['d2fac']<dTH_Fac_Int) ); MaskEH[ind]=1
	ind=np.where( (z0['lc_comp1_2019']==lut_lc['Forest']) & (z0['harv_yr_comp1']==0) & (MaskR==3) & (z0['age_vri']<aTH_Interior) & (z0['fire_yr']==0) & (z0['d2road']<=dTH_Road_Bor) & (z0['d2fac']<dTH_Fac_Bor) ); MaskEH[ind]=1
	
	# Calculate year of harvest
	zDTR=np.zeros(zRef['Data'].shape,dtype='int8')
	ind=np.where(MaskEH==1)
	zDTR[ind]=z0['d2road'][ind]
	#plt.matshow(zDTR,clim=[0,5])
		
	zH=copy.deepcopy(zRef)
	zH['Data']=np.zeros(zRef['Data'].shape,dtype='int16')
	zH['Data'][ind]=1860+(1975-1860)*(1-np.exp(-0.3*(zDTR[ind].astype(float))))+np.random.randint(low=-10,high=10,size=ind[0].size)
	#zH['Data'][0:4000,:]=0
	#plt.close('all'); plt.plot(zDTR[ind][0::20],zH['Data'][ind][0::20],'k.')	
	plt.close('all'); plt.matshow(zH['Data'],clim=[1865,1975])
	gis.SaveGeoTiff(zH,meta['Paths']['bc1ha'] + '\\Disturbances\\Harvest_EarlyReconstruction_Year.tif')
	return

#%%
def DeriveHarvestCompilation(meta):
	zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
	vList=['lc_comp1_2019','harv_yr_cc','harv_yr_ntems','harv_yr_early'] 
	z0=Import_Raster(meta,[],vList,'Extract Grid')
		
	a=gu.IndicesFromUniqueArrayValues(z0['harv_yr_cc'].flatten())
	b=gu.IndicesFromUniqueArrayValues(z0['harv_yr_early'].flatten())
	
	tv=np.arange(1800,2022,1)
	A1=np.zeros(tv.size)
	A2=np.zeros(tv.size)
	for i in tv:
		it=np.where(tv==i)[0]
		try:
			A1[it]=a[i].size
		except:
			pass
	for i in tv:
		it=np.where(tv==i)[0]
		try:
			A2[it]=b[i].size
		except:
			pass
	plt.close('all')
	plt.plot(tv,A1,'-k')
	plt.plot(tv,A1+A2,'-r')
	
	# Harvest consolidated 1 (w/o early reconstruction)
	z=copy.deepcopy(zRef)
	z['Data']=np.zeros(zRef['Data'].shape,dtype='int16')
	ind=np.where( (z0['lc_comp1_2019']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (z0['harv_yr_cc']>0) ); z['Data'][ind]=z0['harv_yr_cc'][ind]
	ind=np.where( (z0['lc_comp1_2019']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (z0['harv_yr_ntems']>0) ); z['Data'][ind]=z0['harv_yr_ntems'][ind]
	z['Data']=z['Data'].astype('int16')
	gis.SaveGeoTiff(z,meta['Paths']['bc1ha'] + '\\Disturbances\\Harvest_Comp1_Year.tif')

	# Harvest consolidated 2 (with early reconstruction)
	z=copy.deepcopy(zRef)
	z['Data']=np.zeros(zRef['Data'].shape,dtype='int16')
	ind=np.where( (z0['lc_comp1_2019']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (z0['harv_yr_early']>0) ); z['Data'][ind]=z0['harv_yr_early'][ind]
	ind=np.where( (z0['lc_comp1_2019']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (z0['harv_yr_cc']>0) ); z['Data'][ind]=z0['harv_yr_cc'][ind]
	ind=np.where( (z0['lc_comp1_2019']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (z0['harv_yr_ntems']>0) ); z['Data'][ind]=z0['harv_yr_ntems'][ind]
	z['Data']=z['Data'].astype('int16')
	gis.SaveGeoTiff(z,meta['Paths']['bc1ha'] + '\\Disturbances\\Harvest_Comp2_Year.tif')
	return

#%% Salvage harvest mask
def HarvestSalvageMaskFromCruise(meta):
	zH=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VEG_CONSOLIDATED_CUT_BLOCKS_SP\\HARVEST_YEAR_YearLast.tif')
	zLCC1=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCover_Comp1_2019.tif')
	zPD=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Disturbances\\HarvestPercentDead_FromCruise.tif')
	
	z=copy.deepcopy(zRef)
	z['Data']=6*np.ones(zRef['Data'].shape,dtype='int8')
	ind=np.where( (zPD['Data']>=50) ); z['Data'][ind]=1 # Salvage high
	ind=np.where( (zPD['Data']>=10) & (zPD['Data']<50) ); z['Data'][ind]=2 # Salvage low
	ind=np.where( (zPD['Data']<10) & (zH['Data']>=2015) ); z['Data'][ind]=3 # Non-salvage harvested forest
	ind=np.where( (zH['Data']<2015) & (zLCC1['Data']==1) ); z['Data'][ind]=4 # Unharvested forest
	ind=np.where( (zRef['Data']==1) & (zLCC1['Data']!=1) ); z['Data'][ind]=5 # Non-forest land
	ind=np.where( (zRef['Data']==0) ); z['Data'][ind]=6 # Non land
	gis.SaveGeoTiff(z,meta['Paths'] + '\\Disturbances\\HarvestSalvageMask_FromCruise.tif')
	return

#%% Generate sparse inputs
# This should speed import by an order of magnitude
def GenerateSparseInputs(meta,rgsf,mask):

	flg=0
	if flg==1:
		rgsf=20
		mask='BCFCS_NOSE'

	z=Import_Raster(meta,[],['refg','lc_comp1_2019'])
	
	if mask=='Province':
		zMask=z['refg']
	elif mask=='BCFCS_NOSE':
		# Non-obligation stand establishment
		zMask=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_NonOb_MaskAll.tif')
	elif mask=='BCFCS_NMC':
		# Nutrient management
		zMask=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\FE-CA_MaskAll.tif')
	elif mask=='BCFCS_LUC':
		# Land use change
		zMask=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Disturbances\\LandCover_Comp1_DeforestationMaskAll.tif')
	elif mask=='BCFCS_EvalAtPlots':
		# Evaluate at plots
		zMask=gu.ipickle(r'D:\Data\FCI_Projects\BCFCS_EvalAtPlots\geos.pkl')['Grid']
	elif mask=='BCFCS_EvalAtCN':
		# Evaluate at plots
		zMask=gu.ipickle(r'D:\Data\FCI_Projects\BCFCS_EvalAtCN\geos.pkl')['Grid']
	elif mask=='BCFCS_Eval':
		# Evaluate
		zMask=gu.ipickle(r'D:\Data\FCI_Projects\BCFCS_Eval\geos.pkl')['Grid']
	elif mask=='BCFCS_EvalCoast':
		# Evaluate Coast
		zMask=gu.ipickle(r'D:\Data\FCI_Projects\BCFCS_EvalCoast\geos.pkl')['Grid']
	elif mask=='BCFCS_EvalInterior':
		# Evaluate Interior
		zMask=gu.ipickle(r'D:\Data\FCI_Projects\BCFCS_EvalInterior\geos.pkl')['Grid']
	elif mask=='TSA_DawsonCreek':
		# Mask from timber supply area
		zTSA=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\FADM_TSA\\TSA_NUMBER_DESCRIPTION.tif')
		zMask=copy.deepcopy(zTSA)
		zMask['Data']=np.zeros(zTSA['Data'].shape,'int8')
		ind=np.where( (zTSA['Data']==meta['LUT']['FADM_TSA']['TSA_NUMBER_DESCRIPTION']['Dawson Creek TSA']) )
		zMask['Data'][ind]=1
	elif mask=='Landscape_NicolaRiverWatershed':
		# Mask straight from BC1ha map
		zMask=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Masks\\Nicola River.tif')
	elif mask=='BCFCS_CWH':
		# Mask from timber supply area
		zBGC=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\BEC_BIOGEOCLIMATIC_POLY\\ZONE.tif')
		zMask=copy.deepcopy(zBGC)
		zMask['Data']=np.zeros(zBGC['Data'].shape,'int8')
		ind=np.where( (zBGC['Data']==meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE']['CWH']) )
		zMask['Data'][ind]=1
	elif mask=='BCFCS_SBS':
		# Mask from timber supply area
		zBGC=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\BEC_BIOGEOCLIMATIC_POLY\\ZONE.tif')
		zMask=copy.deepcopy(zBGC)
		zMask['Data']=np.zeros(zBGC['Data'].shape,'int8')
		ind=np.where( (zBGC['Data']==meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE']['SBS']) )
		zMask['Data'][ind]=1
	else:
		pass

	# Update grid subsampling
	if rgsf>1:
		z['lc_comp1_2019']=gis.UpdateGridCellsize(z['lc_comp1_2019'],rgsf)
		zMask=gis.UpdateGridCellsize(zMask,rgsf)
	
	# Index to ROI subsample
	if mask=='BCFCS_LUC':
		iMask=np.where( (zMask['Data']==1) )
	else:
		iMask=np.where( (zMask['Data']==1) & (z['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Forest']) )
	
	z0=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Disturbances\\Harvest_Comp2_Year.tif')
	z0=gis.UpdateGridCellsize(z0,rgsf)
	gu.opickle(meta['Paths']['bc1ha'] + '\\Sparse\\RGSF' + str(rgsf) + '_Mask' + mask + '_Harvest_Comp2_Year.pkl',z0['Data'][iMask])
	
	# Land cover / land use
	vL=['LandCover_Comp1_1800','LandCover_Comp1_2019','LandUse_Comp1_2019','LandUse_Comp1_2049_Scn1','LandUse_Comp1_2049_Scn2']
	for v in vL:
		z0=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\' + v + '.tif')
		z0=gis.UpdateGridCellsize(z0,rgsf)
		gu.opickle(meta['Paths']['bc1ha'] + '\\Sparse\\RGSF' + str(rgsf) + '_Mask' + mask + '_' + v + '.pkl',z0['Data'][iMask])
	
	# Land use change
	#if mask=='BCFCS_LUC':
	vL=['LandUseChange_Comp1_1800to2019_Type','LandUseChange_Comp1_1800to2019_Year','LandUseChange_Comp1_1800to2019_Type','LandUseChange_Comp1_2020to2049_Scn1_Type','LandUseChange_Comp1_2020to2049_Scn1_Year',
		'LandUseChange_Comp1_2020to2049_Scn2_Type','LandUseChange_Comp1_2020to2049_Scn2_Year']
	for v in vL:
		z0=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Disturbances\\' + v + '.tif')
		z0=gis.UpdateGridCellsize(z0,rgsf)
		gu.opickle(meta['Paths']['bc1ha'] + '\\Sparse\\RGSF' + str(rgsf) + '_Mask' + mask + '_' + v + '.pkl',z0['Data'][iMask])
	
	# Wildfire
	lNam='PROT_HISTORICAL_FIRE_POLYS_SP'
	vNam='WildfireComp1'
	N_Year=6
	d=[None]*N_Year
	for iY in range(N_Year):
		z0=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\' + lNam + '\\' + vNam + '_' + str(iY+1) + '_Year.tif')
		z0=gis.UpdateGridCellsize(z0,rgsf)
		d[iY]=z0['Data'][iMask]
	gu.opickle(meta['Paths']['bc1ha'] + '\\Sparse\\RGSF' + str(rgsf) + '_Mask' + mask + '_' + lNam + '_Year.pkl',d)
	d=[None]*N_Year
	for iY in range(N_Year):
		z0=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\' + lNam + '\\' + vNam + '_' + str(iY+1) + '_DOY.tif')
		z0=gis.UpdateGridCellsize(z0,rgsf)
		d[iY]=z0['Data'][iMask]
	gu.opickle(meta['Paths']['bc1ha'] + '\\Sparse\\RGSF' + str(rgsf) + '_Mask' + mask + '_' + lNam + '_DOY.pkl',d)
	d=[None]*N_Year
	for iY in range(N_Year):
		z0=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\' + lNam + '\\' + vNam + '_' + str(iY+1) + '_SevClass.tif')
		z0=gis.UpdateGridCellsize(z0,rgsf)
		d[iY]=z0['Data'][iMask]
	gu.opickle(meta['Paths']['bc1ha'] + '\\Sparse\\RGSF' + str(rgsf) + '_Mask' + mask + '_' + lNam + '_SevClass.pkl',d)
	   
	# Insect compilation 1
	N_Year=10
	d=[None]*N_Year
	for iY in range(N_Year):
		z0=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\PEST_INFESTATION_POLY\\InsectComp1_' + str(iY+1) + '_Year.tif')
		z0=gis.UpdateGridCellsize(z0,rgsf)
		d[iY]=z0['Data'][iMask]
	gu.opickle(meta['Paths']['bc1ha'] + '\\Sparse\\RGSF' + str(rgsf) + '_Mask' + mask + '_InsectComp1_Year.pkl',d)
	d=[None]*N_Year
	for iY in range(N_Year):
		z0=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\PEST_INFESTATION_POLY\\InsectComp1_' + str(iY+1) + '_Type.tif')
		z0=gis.UpdateGridCellsize(z0,rgsf)
		d[iY]=z0['Data'][iMask]
	gu.opickle(meta['Paths']['bc1ha'] + '\\Sparse\\RGSF' + str(rgsf) + '_Mask' + mask + '_InsectComp1_Type.pkl',d)

	# Disturbance from RESULTS ATU layer
	vL=['C','D','E','F','I','L','P','R','S','W','B']
	for v in vL:
		N_Year=1
		d=[None]*N_Year
		z0=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\DisturbanceFromATU_' + v + '_Year.tif')
		z0=gis.UpdateGridCellsize(z0,rgsf)
		d[0]=z0['Data'][iMask]
		gu.opickle(meta['Paths']['bc1ha'] + '\\Sparse\\RGSF' + str(rgsf) + '_Mask' + mask + '_DisturbanceFromATU_' + v + '_Year.pkl',d)
	
	# Harvest
	N_Year=3
	d=[None]*N_Year
	for iY in range(N_Year):
		z0=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VEG_CONSOLIDATED_CUT_BLOCKS_SP\\HARVEST_YEAR_' + str(iY+1) + '_Year.tif')
		z0=gis.UpdateGridCellsize(z0,rgsf)
		d[iY]=z0['Data'][iMask]
	gu.opickle(meta['Paths']['bc1ha'] + '\\Sparse\\RGSF' + str(rgsf) + '_Mask' + mask + '_VEG_CONSOLIDATED_CUT_BLOCKS_SP_Year.pkl',d)
	
	# Planting year
	N_Year=6
	d=[None]*N_Year
	for iY in range(N_Year):
		z0=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(iY+1) + '_Year.tif')
		z0=gis.UpdateGridCellsize(z0,rgsf)
		d[iY]=z0['Data'][iMask]
	gu.opickle(meta['Paths']['bc1ha'] + '\\Sparse\\RGSF' + str(rgsf) + '_Mask' + mask + '_PL_All_Year.pkl',d)

	d1=[None]*N_Year
	for iY in range(N_Year):
		z0=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(iY+1) + '_Year.tif')
		z0=gis.UpdateGridCellsize(z0,rgsf)
		d1[iY]=z0['Data'][iMask]

	# Planting DOY
	d=[None]*N_Year
	for iY in range(N_Year):
		z0=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(iY+1) + '_DOY.tif')
		z0=gis.UpdateGridCellsize(z0,rgsf)
		d[iY]=z0['Data'][iMask]
	gu.opickle(meta['Paths']['bc1ha'] + '\\Sparse\\RGSF' + str(rgsf) + '_Mask' + mask + '_PL_All_DOY.pkl',d)

	# Planting FSC
	d=[None]*N_Year
	for iY in range(N_Year):
		z0=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(iY+1) + '_SILV_FUND_SOURCE_CODE.tif')
		z0=gis.UpdateGridCellsize(z0,rgsf)
		d[iY]=z0['Data'][iMask]
	gu.opickle(meta['Paths']['bc1ha'] + '\\Sparse\\RGSF' + str(rgsf) + '_Mask' + mask + '_PL_All_SILV_FUND_SOURCE_CODE.pkl',d)
	d=[None]*N_Year
	for iY in range(N_Year):
		z0=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(iY+1) + '_SPH_Planted.tif')
		z0=gis.UpdateGridCellsize(z0,rgsf)
		d[iY]=z0['Data'][iMask]
	gu.opickle(meta['Paths']['bc1ha'] + '\\Sparse\\RGSF' + str(rgsf) + '_Mask' + mask + '_PL_All_SPH_Planted.pkl',d)
	d=[None]*N_Year
	for iY in range(N_Year):
		z0=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(iY+1) + '_ASET.tif')
		z0=gis.UpdateGridCellsize(z0,rgsf)
		d[iY]=z0['Data'][iMask]
	gu.opickle(meta['Paths']['bc1ha'] + '\\Sparse\\RGSF' + str(rgsf) + '_Mask' + mask + '_PL_All_ASET.pkl',d)
	for iSpc in range(6):
		d=[None]*N_Year
		for iY in range(N_Year):
			z0=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(iY+1) + '_PL_SPECIES_CD' + str(iSpc+1) + '.tif')
			z0=gis.UpdateGridCellsize(z0,rgsf)
			d[iY]=z0['Data'][iMask]
		gu.opickle(meta['Paths']['bc1ha'] + '\\Sparse\\RGSF' + str(rgsf) + '_Mask' + mask + '_PL_SPECIES_CD' + str(iSpc+1) + '.pkl',d)
		d=[None]*N_Year
		for iY in range(N_Year):
			z0=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(iY+1) + '_PL_SPECIES_PCT' + str(iSpc+1) + '.tif')
			z0=gis.UpdateGridCellsize(z0,rgsf)
			d[iY]=z0['Data'][iMask]
		gu.opickle(meta['Paths']['bc1ha'] + '\\Sparse\\RGSF' + str(rgsf) + '_Mask' + mask + '_PL_SPECIES_PCT' + str(iSpc+1) + '.pkl',d)
		d=[None]*N_Year
		for iY in range(N_Year):
			z0=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(iY+1) + '_PL_SPECIES_GW' + str(iSpc+1) + '.tif')
			z0=gis.UpdateGridCellsize(z0,rgsf)
			d[iY]=z0['Data'][iMask]
		gu.opickle(meta['Paths']['bc1ha'] + '\\Sparse\\RGSF' + str(rgsf) + '_Mask' + mask + '_PL_SPECIES_GW' + str(iSpc+1) + '.pkl',d)
	
	# Nutrient management
	N_Year=3
	d=[None]*N_Year
	dDOY=[None]*N_Year
	dFSC=[None]*N_Year
	for iY in range(N_Year):
		z0=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\FE-CA_' + str(iY+1) + '_Year.tif')
		z0=gis.UpdateGridCellsize(z0,rgsf)
		d[iY]=z0['Data'][iMask]
		z0=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\FE-CA_' + str(iY+1) + '_DOY.tif')
		z0=gis.UpdateGridCellsize(z0,rgsf)
		dDOY[iY]=z0['Data'][iMask]
		z0=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\FE-CA_' + str(iY+1) + '_SILV_FUND_SOURCE_CODE.tif')
		z0=gis.UpdateGridCellsize(z0,rgsf)
		dFSC[iY]=z0['Data'][iMask]
	gu.opickle(meta['Paths']['bc1ha'] + '\\Sparse\\RGSF' + str(rgsf) + '_Mask' + mask + '_FE-CA_Year.pkl',d)
	gu.opickle(meta['Paths']['bc1ha'] + '\\Sparse\\RGSF' + str(rgsf) + '_Mask' + mask + '_FE-CA_DOY.pkl',dDOY)
	gu.opickle(meta['Paths']['bc1ha'] + '\\Sparse\\RGSF' + str(rgsf) + '_Mask' + mask + '_FE-CA_SILV_FUND_SOURCE_CODE.pkl',dFSC)
	
	# Knockdown
	N_Year=3
	d=[None]*N_Year
	dDOY=[None]*N_Year
	dFSC=[None]*N_Year
	for iY in range(N_Year):
		z0=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\SP-KD_' + str(iY+1) + '_Year.tif')
		z0=gis.UpdateGridCellsize(z0,rgsf)
		d[iY]=z0['Data'][iMask]
		z0=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\SP-KD_' + str(iY+1) + '_DOY.tif')
		z0=gis.UpdateGridCellsize(z0,rgsf)
		dDOY[iY]=z0['Data'][iMask]
		z0=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\SP-KD_' + str(iY+1) + '_SILV_FUND_SOURCE_CODE.tif')
		z0=gis.UpdateGridCellsize(z0,rgsf)
		dFSC[iY]=z0['Data'][iMask]
	gu.opickle(meta['Paths']['bc1ha'] + '\\Sparse\\RGSF' + str(rgsf) + '_Mask' + mask + '_SP-KD_Year.pkl',d)
	gu.opickle(meta['Paths']['bc1ha'] + '\\Sparse\\RGSF' + str(rgsf) + '_Mask' + mask + '_SP-KD_DOY.pkl',dDOY)
	gu.opickle(meta['Paths']['bc1ha'] + '\\Sparse\\RGSF' + str(rgsf) + '_Mask' + mask + '_SP-KD_SILV_FUND_SOURCE_CODE.pkl',dFSC)

	# Ripping
	N_Year=3
	d=[None]*N_Year
	dDOY=[None]*N_Year
	dFSC=[None]*N_Year
	for iY in range(N_Year):
		z0=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\SP-Rip_' + str(iY+1) + '_Year.tif')
		z0=gis.UpdateGridCellsize(z0,rgsf)
		d[iY]=z0['Data'][iMask]
		z0=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\SP-Rip_' + str(iY+1) + '_DOY.tif')
		z0=gis.UpdateGridCellsize(z0,rgsf)
		dDOY[iY]=z0['Data'][iMask]
		z0=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\SP-Rip_' + str(iY+1) + '_SILV_FUND_SOURCE_CODE.tif')
		z0=gis.UpdateGridCellsize(z0,rgsf)
		dFSC[iY]=z0['Data'][iMask]
	gu.opickle(meta['Paths']['bc1ha'] + '\\Sparse\\RGSF' + str(rgsf) + '_Mask' + mask + '_SP-Rip_Year.pkl',d)
	gu.opickle(meta['Paths']['bc1ha'] + '\\Sparse\\RGSF' + str(rgsf) + '_Mask' + mask + '_SP-Rip_DOY.pkl',dDOY)
	gu.opickle(meta['Paths']['bc1ha'] + '\\Sparse\\RGSF' + str(rgsf) + '_Mask' + mask + '_SP-Rip_SILV_FUND_SOURCE_CODE.pkl',dFSC)

	# Pile burn
	N_Year=3
	d=[None]*N_Year
	dDOY=[None]*N_Year
	dFSC=[None]*N_Year
	for iY in range(N_Year):
		z0=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\SP-PBURN_' + str(iY+1) + '_Year.tif')
		z0=gis.UpdateGridCellsize(z0,rgsf)
		d[iY]=z0['Data'][iMask]
		z0=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\SP-PBURN_' + str(iY+1) + '_DOY.tif')
		z0=gis.UpdateGridCellsize(z0,rgsf)
		dDOY[iY]=z0['Data'][iMask]
		z0=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\SP-PBURN_' + str(iY+1) + '_SILV_FUND_SOURCE_CODE.tif')
		z0=gis.UpdateGridCellsize(z0,rgsf)
		dFSC[iY]=z0['Data'][iMask]
	gu.opickle(meta['Paths']['bc1ha'] + '\\Sparse\\RGSF' + str(rgsf) + '_Mask' + mask + '_SP-PBURN_Year.pkl',d)
	gu.opickle(meta['Paths']['bc1ha'] + '\\Sparse\\RGSF' + str(rgsf) + '_Mask' + mask + '_SP-PBURN_DOY.pkl',dDOY)
	gu.opickle(meta['Paths']['bc1ha'] + '\\Sparse\\RGSF' + str(rgsf) + '_Mask' + mask + '_SP-PBURN_SILV_FUND_SOURCE_CODE.pkl',dFSC)

	# Burn at landings
	N_Year=3
	d=[None]*N_Year
	dDOY=[None]*N_Year
	dFSC=[None]*N_Year
	for iY in range(N_Year):
		z0=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\SP-BU-LAND_' + str(iY+1) + '_Year.tif')
		z0=gis.UpdateGridCellsize(z0,rgsf)
		d[iY]=z0['Data'][iMask]
		z0=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\SP-BU-LAND_' + str(iY+1) + '_DOY.tif')
		z0=gis.UpdateGridCellsize(z0,rgsf)
		dDOY[iY]=z0['Data'][iMask]
		z0=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\SP-BU-LAND_' + str(iY+1) + '_SILV_FUND_SOURCE_CODE.tif')
		z0=gis.UpdateGridCellsize(z0,rgsf)
		dFSC[iY]=z0['Data'][iMask]
	gu.opickle(meta['Paths']['bc1ha'] + '\\Sparse\\RGSF' + str(rgsf) + '_Mask' + mask + '_SP-BU-LAND_Year.pkl',d)
	gu.opickle(meta['Paths']['bc1ha'] + '\\Sparse\\RGSF' + str(rgsf) + '_Mask' + mask + '_SP-BU-LAND_DOY.pkl',dDOY)
	gu.opickle(meta['Paths']['bc1ha'] + '\\Sparse\\RGSF' + str(rgsf) + '_Mask' + mask + '_SP-BU-LAND_SILV_FUND_SOURCE_CODE.pkl',dFSC)

	# Piling
	N_Year=3
	d=[None]*N_Year
	dDOY=[None]*N_Year
	dFSC=[None]*N_Year
	for iY in range(N_Year):
		z0=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\SP-PILE_' + str(iY+1) + '_Year.tif')
		z0=gis.UpdateGridCellsize(z0,rgsf)
		d[iY]=z0['Data'][iMask]
		z0=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\SP-PILE_' + str(iY+1) + '_DOY.tif')
		z0=gis.UpdateGridCellsize(z0,rgsf)
		dDOY[iY]=z0['Data'][iMask]
		z0=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\SP-PILE_' + str(iY+1) + '_SILV_FUND_SOURCE_CODE.tif')
		z0=gis.UpdateGridCellsize(z0,rgsf)
		dFSC[iY]=z0['Data'][iMask]
	gu.opickle(meta['Paths']['bc1ha'] + '\\Sparse\\RGSF' + str(rgsf) + '_Mask' + mask + '_SP-PILE_Year.pkl',d)
	gu.opickle(meta['Paths']['bc1ha'] + '\\Sparse\\RGSF' + str(rgsf) + '_Mask' + mask + '_SP-PILE_DOY.pkl',dDOY)
	gu.opickle(meta['Paths']['bc1ha'] + '\\Sparse\\RGSF' + str(rgsf) + '_Mask' + mask + '_SP-PILE_SILV_FUND_SOURCE_CODE.pkl',dFSC)

	# Prescribed burning
	N_Year=3
	d=[None]*N_Year
	dDOY=[None]*N_Year
	dFSC=[None]*N_Year
	for iY in range(N_Year):
		z0=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\SP-BU-BROAD_' + str(iY+1) + '_Year.tif')
		z0=gis.UpdateGridCellsize(z0,rgsf)
		d[iY]=z0['Data'][iMask]
		z0=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\SP-BU-BROAD_' + str(iY+1) + '_DOY.tif')
		z0=gis.UpdateGridCellsize(z0,rgsf)
		dDOY[iY]=z0['Data'][iMask]
	gu.opickle(meta['Paths']['bc1ha'] + '\\Sparse\\RGSF' + str(rgsf) + '_Mask' + mask + '_SP-BU-BROAD_Year.pkl',d)
	gu.opickle(meta['Paths']['bc1ha'] + '\\Sparse\\RGSF' + str(rgsf) + '_Mask' + mask + '_SP-BU-BROAD_DOY.pkl',dDOY)
	gu.opickle(meta['Paths']['bc1ha'] + '\\Sparse\\RGSF' + str(rgsf) + '_Mask' + mask + '_SP_BU-BROAD_SILV_FUND_SOURCE_CODE.pkl',dFSC)
	
	return

#%% Derive Artificial Stand Establishment Type (used for non-ob stand establishment)
def DeriveASETComp1(meta,YearLast):
	#YearLast=2023
	# Notes:
	# Must run RasterizePlanting first, but not MaskPlantingNonOb not necessary

	meta=cbu.Load_LUTs_Modelling(meta)
	zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
	zD=Import_Raster(meta,[],['harv_yr_comp1','roadrehab_yl','pdead_cruise','spc1_vri02'])

	ptNam=np.array(list(meta['LUT']['Derived']['ASET'].keys()))
	ptID=np.array(list(meta['LUT']['Derived']['ASET'].values()))
	cd=meta['LUT']['RSLT_ACTIVITY_TREATMENT_SVW']['SILV_TECHNIQUE_CODE']
	stcInclude=[cd['PL'],cd['RP'],cd['FP'],cd['RR']]
	psevInclude=[meta['LUT']['PEST_INFESTATION_POLY']['PEST_SEVERITY_CODE']['L'],meta['LUT']['PEST_INFESTATION_POLY']['PEST_SEVERITY_CODE']['M'],meta['LUT']['PEST_INFESTATION_POLY']['PEST_SEVERITY_CODE']['S'],meta['LUT']['PEST_INFESTATION_POLY']['PEST_SEVERITY_CODE']['V']]

	# Track annual summary of area treated
	d={}
	d['tv']=np.arange(1960,YearLast+1,1)
	d['A Tot']=np.zeros(d['tv'].size)
	d['A All']=np.zeros((d['tv'].size,ptNam.size))
	d['A NOSE']=np.zeros((d['tv'].size,ptNam.size))

	# Import planting events
	zP={}
	for iP in range(6):
		zP[iP]={}
		zP[iP]['Type']=np.zeros(zRef['Data'].shape,dtype='int8')
		zP[iP]['Year']=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(iP+1) + '_Year.tif')['Data']
		zP[iP]['DOY']=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(iP+1) + '_DOY.tif')['Data']
		zP[iP]['STC']=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(iP+1) + '_SILV_TECHNIQUE_CODE.tif')['Data']
		zP[iP]['SMC']=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(iP+1) + '_SILV_METHOD_CODE.tif')['Data']
		zP[iP]['FSC']=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(iP+1) + '_SILV_FUND_SOURCE_CODE.tif')['Data']
		zP[iP]['OBJ']=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(iP+1) + '_SILV_OBJECTIVE_CODE_1.tif')['Data']
		zP[iP]['SP1']=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(iP+1) + '_PL_SPECIES_CD1.tif')['Data']

	# Import knockdown events
	zKD={}
	for iKD in range(3):
		zKD[iKD]={}
		zKD[iKD]['Year']=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\SP-KD_' + str(iKD+1) + '_Year.tif')['Data']
		zKD[iKD]['DOY']=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\SP-KD_' + str(iKD+1) + '_DOY.tif')['Data']

	# Import direct seeding events
	zDS={}
	for iDS in range(3):
		zDS[iDS]={}
		zDS[iDS]['Year']=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\DS_All_' + str(iDS+1) + '_Year.tif')['Data']

	# Initialize last disturbance type
	LastDisturbance=np.zeros(zRef['Data'].shape,dtype='int8')

	# History of beetles
	BH=np.zeros(zRef['Data'].shape,dtype='int8')
	BH_counter=np.zeros(zRef['Data'].shape,dtype='int8')

	# History of harvesting
	HH=np.zeros(zRef['Data'].shape,dtype='int8')
	HH_counter=np.zeros(zRef['Data'].shape,dtype='int8')

	for iT in range(d['tv'].size):
		print(d['tv'][iT])

		# Import disturbance events (confined to the previous year)
		# *** This may miss events that occur in the same year as planting ***
		fire_oc=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\PROT_HISTORICAL_FIRE_POLYS_SP\\PROT_HISTORICAL_FIRE_POLYS_SP_' + str(d['tv'][iT]) + '_Year.tif')['Data']
		fire_doy=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\PROT_HISTORICAL_FIRE_POLYS_SP\\PROT_HISTORICAL_FIRE_POLYS_SP_' + str(d['tv'][iT]) + '_DOY.tif')['Data']

		ibm_oc=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\PEST_INFESTATION_POLY\\PEST_SEVERITY_CODE_IBM_' + str(d['tv'][iT]) + '.tif')['Data']
		ibs_oc=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\PEST_INFESTATION_POLY\\PEST_SEVERITY_CODE_IBS_' + str(d['tv'][iT]) + '.tif')['Data']
		ibd_oc=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\PEST_INFESTATION_POLY\\PEST_SEVERITY_CODE_IBD_' + str(d['tv'][iT]) + '.tif')['Data']

		try:
			# If you only check pervious year, you can miss instances where KD and PL
			# occur during the same year.
			#zH=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VEG_CONSOLIDATED_CUT_BLOCKS_SP\\HARVEST_YEAR_' + str(d['tv'][iT]-1) + '.tif')['Data']
			harv_oc=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VEG_CONSOLIDATED_CUT_BLOCKS_SP\\HARVEST_YEAR_' + str(d['tv'][iT]) + '.tif')['Data']
		except:
			# 2023 not yet available
			print('Harvest year missing!')
			harv_oc=np.zeros(zRef['Data'].shape)

		# Extract knockdown for current time step
		kd_oc=np.zeros(zRef['Data'].shape,dtype='int8')
		kd_doy=np.zeros(zRef['Data'].shape,dtype='int16')
		for iP in range(3):
			ind=np.where( (zKD[iP]['Year']==d['tv'][iT]) )
			kd_oc[ind]=1
			kd_doy[ind]=zKD[iP]['DOY'][ind]

		# Extract planting for current time step
		pl_oc=np.zeros(zRef['Data'].shape,dtype='int8')
		pl_doy=np.zeros(zRef['Data'].shape,dtype='int16')
		pl_stc=np.zeros(zRef['Data'].shape,dtype='int16')
		pl_smc=np.zeros(zRef['Data'].shape,dtype='int16')
		pl_fsc=np.zeros(zRef['Data'].shape,dtype='int16')
		pl_obj=np.zeros(zRef['Data'].shape,dtype='int16')
		pl_spc=np.zeros(zRef['Data'].shape,dtype='int16')
		for iP in range(6):
			ind=np.where( (zP[iP]['Year']==d['tv'][iT]) & (np.isin(zP[iP]['STC'],stcInclude)==True) )
			pl_oc[ind]=1
			pl_doy[ind]=zP[iP]['DOY'][ind]
			pl_stc[ind]=zP[iP]['STC'][ind]
			pl_smc[ind]=zP[iP]['STC'][ind]
			pl_fsc[ind]=zP[iP]['FSC'][ind]
			pl_obj[ind]=zP[iP]['OBJ'][ind]
			pl_spc[ind]=zP[iP]['SP1'][ind]

		# Populate type of last disturbance

		# Insects have no DOY, run first
		ind=np.where(np.isin(ibm_oc,psevInclude)==True)
		LastDisturbance[ind]=meta['LUT']['Event']['Mountain Pine Beetle'] # Severity greater than light
		BH[ind]=1
		BH_counter[ind]=BH_counter[ind]+1
		ind_stop=np.where(BH_counter>20)
		BH[ind_stop]=0
		BH_counter[ind_stop]=0

		ind=np.where(np.isin(ibs_oc,psevInclude)==True)
		LastDisturbance[ind]=meta['LUT']['Event']['Spruce Beetle']
		BH[ind]=1
		BH_counter[ind]=BH_counter[ind]+1
		ind_stop=np.where(BH_counter>20)
		BH[ind_stop]=0
		BH_counter[ind_stop]=0

		ind=np.where(np.isin(ibd_oc,psevInclude)==True)
		LastDisturbance[ind]=meta['LUT']['Event']['Douglas-fir Beetle']
		BH[ind]=1
		BH_counter[ind]=BH_counter[ind]+1
		ind_stop=np.where(BH_counter>20)
		BH[ind_stop]=0
		BH_counter[ind_stop]=0

		# Harvest has no DOY, run second
		ind=np.where(harv_oc>0)
		LastDisturbance[ind]=meta['LUT']['Event']['Harvest']
		HH[ind]=1
		HH_counter[ind]=HH_counter[ind]+1
		ind_stop=np.where(HH_counter>10)
		HH[ind_stop]=0
		HH_counter[ind_stop]=0

		# Wilfire
		ind=np.where( (fire_oc>0) )
		LastDisturbance[ind]=meta['LUT']['Event']['Wildfire']

		# If you only check pervious year, you can miss instances where KD and PL
		# occur during the same year.
		ind=np.where( (kd_oc>0) )
		LastDisturbance[ind]=meta['LUT']['Event']['Knockdown']

		ind=np.where(zD['roadrehab_yl']['Data']==d['tv'][iT])
		LastDisturbance[ind]=meta['LUT']['Event']['Road Deactivation']

		ind=np.where( (pl_oc==1) )
		d['A Tot'][iT]=ind[0].size

		# Initialize a grid where planting will be tracked. In many cases, the
		# classification of each type will only be recorded if the tracking grid
		# has not already been classified. But there are exceptions (e.g. UNDER).
		zType0=np.zeros(zRef['Data'].shape,dtype='int8')

		#----------------------------------------------------------------------
		# Replanting
		#----------------------------------------------------------------------
		ind=np.where( (pl_oc==1) & (pl_stc==meta['LUT']['RSLT_ACTIVITY_TREATMENT_SVW']['SILV_TECHNIQUE_CODE']['RP']) & (np.isin(pl_fsc,meta['Param']['Raw']['FSC']['NO List ID'])==True) )
		d['A NOSE'][iT,meta['LUT']['Derived']['ASET']['Replanting']-1]=ind[0].size
		ind=np.where( (pl_oc==1) & (pl_stc==meta['LUT']['RSLT_ACTIVITY_TREATMENT_SVW']['SILV_TECHNIQUE_CODE']['RP']) )
		d['A All'][iT,meta['LUT']['Derived']['ASET']['Replanting']-1]=ind[0].size
		zType0[ind]=meta['LUT']['Derived']['ASET']['Replanting']

		#----------------------------------------------------------------------
		# Fill planting
		#----------------------------------------------------------------------
		ind=np.where( (pl_oc==1) & (pl_stc==meta['LUT']['RSLT_ACTIVITY_TREATMENT_SVW']['SILV_TECHNIQUE_CODE']['FP']) & (np.isin(pl_fsc,meta['Param']['Raw']['FSC']['NO List ID'])==True) )
		d['A NOSE'][iT,meta['LUT']['Derived']['ASET']['Fill Planting']-1]=ind[0].size
		ind=np.where( (pl_oc==1) & (pl_stc==meta['LUT']['RSLT_ACTIVITY_TREATMENT_SVW']['SILV_TECHNIQUE_CODE']['FP']) )
		d['A All'][iT,meta['LUT']['Derived']['ASET']['Fill Planting']-1]=ind[0].size
		zType0[ind]=meta['LUT']['Derived']['ASET']['Fill Planting']

		#----------------------------------------------------------------------
		# Ecosystem restoration
		#----------------------------------------------------------------------
		ind=np.where( (pl_oc==1) & (np.isin(pl_fsc,meta['Param']['Raw']['FSC']['NO List ID'])==True) & (pl_obj==meta['LUT']['RSLT_ACTIVITY_TREATMENT_SVW']['SILV_OBJECTIVE_CODE_1']['ER']) |
			(pl_oc==1) & (np.isin(pl_fsc,meta['Param']['Raw']['FSC']['NO List ID'])==True) & (pl_obj==meta['LUT']['RSLT_ACTIVITY_TREATMENT_SVW']['SILV_OBJECTIVE_CODE_1']['HER']) |
			(pl_oc==1) & (np.isin(pl_fsc,meta['Param']['Raw']['FSC']['NO List ID'])==True) &
			(np.isin(pl_spc,[meta['LUT']['VEG_COMP_LYR_R1_POLY']['SPECIES_CD_1']['AT'],meta['LUT']['VEG_COMP_LYR_R1_POLY']['SPECIES_CD_1']['ACT'],meta['LUT']['VEG_COMP_LYR_R1_POLY']['SPECIES_CD_1']['PA']])==True) )
		d['A NOSE'][iT,meta['LUT']['Derived']['ASET']['Ecosystem Restoration']-1]=ind[0].size
		d['A All'][iT,meta['LUT']['Derived']['ASET']['Ecosystem Restoration']-1]=ind[0].size
		zType0[ind]=meta['LUT']['Derived']['ASET']['Ecosystem Restoration']

		#----------------------------------------------------------------------
		# Road rehab
		#----------------------------------------------------------------------
		# Missing plenty, being retroactively fixed in RESULTS ***
		ind=np.where( (pl_oc==1) & (pl_stc==meta['LUT']['RSLT_ACTIVITY_TREATMENT_SVW']['SILV_TECHNIQUE_CODE']['RR']) & (np.isin(pl_fsc,meta['Param']['Raw']['FSC']['NO List ID'])==True) |
			   (pl_oc==1) & (LastDisturbance==meta['LUT']['Event']['Road Deactivation']) & (np.isin(pl_fsc,meta['Param']['Raw']['FSC']['NO List ID'])==True))
		d['A NOSE'][iT,meta['LUT']['Derived']['ASET']['Road Rehabilitation']-1]=ind[0].size
		ind=np.where( (pl_oc==1) & (pl_stc==meta['LUT']['RSLT_ACTIVITY_TREATMENT_SVW']['SILV_TECHNIQUE_CODE']['RR']) |
					  (pl_oc==1) & (LastDisturbance==meta['LUT']['Event']['Road Deactivation']) )
		d['A All'][iT,meta['LUT']['Derived']['ASET']['Road Rehabilitation']-1]=ind[0].size
		zType0[ind]=meta['LUT']['Derived']['ASET']['Road Rehabilitation']

		#----------------------------------------------------------------------
		# Back to back planting
		#----------------------------------------------------------------------
		ind=np.where( (pl_oc==1) & (zType0==0) & (LastDisturbance==meta['LUT']['Event']['Planting']) & (np.isin(pl_fsc,meta['Param']['Raw']['FSC']['NO List ID'])==True) )
		d['A NOSE'][iT,meta['LUT']['Derived']['ASET']['Back-to-back Planting']-1]=ind[0].size
		ind=np.where( (pl_oc==1) & (zType0==0) & (LastDisturbance==meta['LUT']['Event']['Planting']) )
		d['A All'][iT,meta['LUT']['Derived']['ASET']['Back-to-back Planting']-1]=ind[0].size
		zType0[ind]=meta['LUT']['Derived']['ASET']['Back-to-back Planting']

		#----------------------------------------------------------------------
		# Salvage with beetles
		#----------------------------------------------------------------------
		ind=np.where( (pl_oc==1) &
			   (np.isin(pl_fsc,meta['Param']['Raw']['FSC']['NO List ID'])==True) &
			   (d['tv'][iT]>2000) &
			   (zD['harv_yr_comp1']['Data']>1987) &
			   (zType0==0) &
			   (LastDisturbance==meta['LUT']['Event']['Harvest']) &
			   (BH>0) )
		d['A NOSE'][iT,meta['LUT']['Derived']['ASET']['Salvage and Planting Post Beetle']-1]=ind[0].size
		d['A All'][iT,meta['LUT']['Derived']['ASET']['Salvage and Planting Post Beetle']-1]=ind[0].size
		zType0[ind]=meta['LUT']['Derived']['ASET']['Salvage and Planting Post Beetle']

		# Licensees - based on percent dead after 2014
		ind=np.where( (pl_oc==1) &
			   (np.isin(pl_fsc,meta['Param']['Raw']['FSC']['NO List ID'])==False) &
			   (d['tv'][iT]>=2015) &
			   (zD['harv_yr_comp1']['Data']>1987) &
			   (zType0==0) &
			   (LastDisturbance==meta['LUT']['Event']['Harvest']) &
			   (zD['pdead_cruise']['Data']>25) &
			   (BH>0) )
		d['A All'][iT,meta['LUT']['Derived']['ASET']['Salvage and Planting Post Beetle']-1]=d['A All'][iT,meta['LUT']['Derived']['ASET']['Salvage and Planting Post Beetle']-1]+ind[0].size
		zType0[ind]=meta['LUT']['Derived']['ASET']['Salvage and Planting Post Beetle']

		# Licensees - based on beetle history and leading species=PL, before 2015
		ind=np.where( (pl_oc==1) &
			   (np.isin(pl_fsc,meta['Param']['Raw']['FSC']['NO List ID'])==False) &
			   (d['tv'][iT]<2015) &
			   (zD['harv_yr_comp1']['Data']>1987) &
			   (np.isin(zD['spc1_vri02']['Data'],[meta['LUT']['VEG_COMP_LYR_R1_POLY']['SPECIES_CD_1']['PL'],meta['LUT']['VEG_COMP_LYR_R1_POLY']['SPECIES_CD_1']['PLI']])==True) &
			   (zType0==0) &
			   (LastDisturbance==meta['LUT']['Event']['Harvest']) &
			   (BH>0) )
		d['A All'][iT,meta['LUT']['Derived']['ASET']['Salvage and Planting Post Beetle']-1]=d['A All'][iT,meta['LUT']['Derived']['ASET']['Salvage and Planting Post Beetle']-1]+ind[0].size
		zType0[ind]=meta['LUT']['Derived']['ASET']['Salvage and Planting Post Beetle']

		#----------------------------------------------------------------------
		# Salvage with no beetle history
		#----------------------------------------------------------------------
		ind=np.where( (pl_oc==1) &
			   (np.isin(pl_fsc,meta['Param']['Raw']['FSC']['NO List ID'])==True) &
			   (d['tv'][iT]>2000) &
			   (zD['harv_yr_comp1']['Data']>1987) &
			   (zType0==0) &
			   (LastDisturbance==meta['LUT']['Event']['Harvest']) )
		d['A NOSE'][iT,meta['LUT']['Derived']['ASET']['Salvage and Planting Post Other']-1]=ind[0].size
		d['A All'][iT,meta['LUT']['Derived']['ASET']['Salvage and Planting Post Other']-1]=ind[0].size
		zType0[ind]=meta['LUT']['Derived']['ASET']['Salvage and Planting Post Other']
		ind=np.where( (pl_oc==1) &
			   (np.isin(pl_fsc,meta['Param']['Raw']['FSC']['NO List ID'])==False) &
			   (d['tv'][iT]>2000) &
			   (zD['harv_yr_comp1']['Data']>1987) &
			   (zType0==0) &
			   (LastDisturbance==meta['LUT']['Event']['Harvest']) &
			   (zD['pdead_cruise']['Data']>20) )
		d['A All'][iT,meta['LUT']['Derived']['ASET']['Salvage and Planting Post Other']-1]=d['A All'][iT,meta['LUT']['Derived']['ASET']['Salvage and Planting Post Other']-1]+ind[0].size
		zType0[ind]=meta['LUT']['Derived']['ASET']['Salvage and Planting Post Other']

		#----------------------------------------------------------------------
		# Knockdown
		#----------------------------------------------------------------------
		ind=np.where( (pl_oc==1) & (zType0==0) & (LastDisturbance==meta['LUT']['Event']['Knockdown']) & (np.isin(pl_fsc,meta['Param']['Raw']['FSC']['NO List ID'])==True) )
		d['A NOSE'][iT,meta['LUT']['Derived']['ASET']['Knockdown and Planting']-1]=ind[0].size
		ind=np.where( (pl_oc==1) & (zType0==0) & (LastDisturbance==meta['LUT']['Event']['Knockdown']) )
		d['A All'][iT,meta['LUT']['Derived']['ASET']['Knockdown and Planting']-1]=ind[0].size
		zType0[ind]=meta['LUT']['Derived']['ASET']['Knockdown and Planting']

		#----------------------------------------------------------------------
		# Straight fire (prior to use of SMC='UNDER')
		#----------------------------------------------------------------------
		ind=np.where( (pl_oc==1) & (zType0==0) & (LastDisturbance==meta['LUT']['Event']['Wildfire']) & (np.isin(pl_fsc,meta['Param']['Raw']['FSC']['NO List ID'])==True) )
		d['A NOSE'][iT,meta['LUT']['Derived']['ASET']['Underplanting']-1]=ind[0].size
		ind=np.where( (pl_oc==1) & (zType0==0) & (LastDisturbance==meta['LUT']['Event']['Wildfire']) )
		d['A All'][iT,meta['LUT']['Derived']['ASET']['Underplanting']-1]=ind[0].size
		zType0[ind]=meta['LUT']['Derived']['ASET']['Underplanting']

		#----------------------------------------------------------------------
		# Straight fire (use of SMC='UNDER')
		#----------------------------------------------------------------------
		ind=np.where( (pl_oc==1) & (pl_smc==meta['LUT']['RSLT_ACTIVITY_TREATMENT_SVW']['SILV_METHOD_CODE']['UNDER']) & (LastDisturbance==meta['LUT']['Event']['Wildfire']) & (np.isin(pl_fsc,meta['Param']['Raw']['FSC']['NO List ID'])==True) )
		d['A NOSE'][iT,meta['LUT']['Derived']['ASET']['Underplanting']-1]=d['A NOSE'][iT,meta['LUT']['Derived']['ASET']['Underplanting']-1]+ind[0].size
		ind=np.where( (pl_oc==1) & (pl_smc==meta['LUT']['RSLT_ACTIVITY_TREATMENT_SVW']['SILV_METHOD_CODE']['UNDER']) & (LastDisturbance==meta['LUT']['Event']['Wildfire']) )
		d['A All'][iT,meta['LUT']['Derived']['ASET']['Underplanting']-1]=d['A All'][iT,meta['LUT']['Derived']['ASET']['Underplanting']-1]+ind[0].size
		zType0[ind]=meta['LUT']['Derived']['ASET']['Underplanting']

		#----------------------------------------------------------------------
		# Straight insect
		#----------------------------------------------------------------------

		# Some straight insect is very likely salvage harvesting
		ind=np.where( (pl_oc==1) &
			   (zType0==0) &
			   (np.isin(LastDisturbance,[meta['LUT']['Event']['Douglas-fir Beetle'],meta['LUT']['Event']['Spruce Beetle'],meta['LUT']['Event']['Mountain Pine Beetle']])==True) &
			   (HH>0) &
			   (np.isin(pl_fsc,meta['Param']['Raw']['FSC']['NO List ID'])==True) )
		d['A NOSE'][iT,meta['LUT']['Derived']['ASET']['Salvage and Planting Post Beetle']-1]=d['A NOSE'][iT,meta['LUT']['Derived']['ASET']['Salvage and Planting Post Beetle']-1]+ind[0].size
		ind=np.where( (pl_oc==1) &
			   (zType0==0) &
			   (np.isin(LastDisturbance,[meta['LUT']['Event']['Douglas-fir Beetle'],meta['LUT']['Event']['Spruce Beetle'],meta['LUT']['Event']['Mountain Pine Beetle']])==True) &
			   (HH>0) )
		d['A All'][iT,meta['LUT']['Derived']['ASET']['Salvage and Planting Post Beetle']-1]=d['A All'][iT,meta['LUT']['Derived']['ASET']['Salvage and Planting Post Beetle']-1]+ind[0].size
		zType0[ind]=meta['LUT']['Derived']['ASET']['Salvage and Planting Post Beetle']

		# Classify remaining non-obligation as straight to insects
		ind=np.where( (pl_oc==1) &
			   (zType0==0) &
			   (np.isin(LastDisturbance,[meta['LUT']['Event']['Douglas-fir Beetle'],meta['LUT']['Event']['Spruce Beetle'],meta['LUT']['Event']['Mountain Pine Beetle']])==True) &
			   (np.isin(pl_fsc,meta['Param']['Raw']['FSC']['NO List ID'])==True) )
		d['A NOSE'][iT,meta['LUT']['Derived']['ASET']['Straight-to-planting Post Beetles']-1]=d['A NOSE'][iT,meta['LUT']['Derived']['ASET']['Straight-to-planting Post Beetles']-1]+ind[0].size
		zType0[ind]=meta['LUT']['Derived']['ASET']['Straight-to-planting Post Beetles']

		# Classify remaining licensee planting as 'Harvest and Planting'
		ind=np.where( (pl_oc==1) &
			   (zType0==0) &
			   (np.isin(LastDisturbance,[meta['LUT']['Event']['Douglas-fir Beetle'],meta['LUT']['Event']['Spruce Beetle'],meta['LUT']['Event']['Mountain Pine Beetle']])==True) &
			   (np.isin(pl_fsc,meta['Param']['Raw']['FSC']['NO List ID'])==False) )
		d['A All'][iT,meta['LUT']['Derived']['ASET']['Harvest and Planting']-1]=d['A All'][iT,meta['LUT']['Derived']['ASET']['Harvest and Planting']-1]+ind[0].size
		zType0[ind]=meta['LUT']['Derived']['ASET']['Harvest and Planting']

		#----------------------------------------------------------------------
		# Straight other
		# Many should be classified as Ecosystem Restoration now - consider retiring
		#----------------------------------------------------------------------
		# n/a

		#----------------------------------------------------------------------
		# Agroforestry and Short Rotation Crops
		#----------------------------------------------------------------------
		# n/a

		#----------------------------------------------------------------------
		# Harvest and Planting NSR backlog
		#----------------------------------------------------------------------
		ind=np.where( (pl_oc==1) &
			   (zType0==0) &
			   (LastDisturbance==meta['LUT']['Event']['Harvest']) &
			   (zD['harv_yr_comp1']['Data']<=1987) &
			   (np.isin(pl_fsc,meta['Param']['Raw']['FSC']['NO List ID'])==True) )
		d['A NOSE'][iT,meta['LUT']['Derived']['ASET']['Harvest and Planting NSR Backlog']-1]=ind[0].size
		d['A All'][iT,meta['LUT']['Derived']['ASET']['Harvest and Planting NSR Backlog']-1]=ind[0].size
		zType0[ind]=meta['LUT']['Derived']['ASET']['Harvest and Planting NSR Backlog']

		#----------------------------------------------------------------------
		# Harvest and Planting
		#----------------------------------------------------------------------
		#ind=np.where( (pl_oc==1) &
		#	   (zType0==0) &
		#	   (LastDisturbance==meta['LUT']['Event']['Harvest']) &
		#	   (np.isin(pl_fsc,meta['Param']['Raw']['FSC']['NO List ID'])==True) )
		#d['A NOSE'][iT,meta['LUT']['Derived']['ASET']['Harvest and Planting']-1]=d['A All'][iT,meta['LUT']['Derived']['ASET']['Harvest and Planting']-1]+ind[0].size
		ind=np.where( (pl_oc==1) &
			   (zType0==0) &
			   (LastDisturbance==meta['LUT']['Event']['Harvest']) &
			   (np.isin(pl_fsc,meta['Param']['Raw']['FSC']['NO List ID'])==False) )
		d['A All'][iT,meta['LUT']['Derived']['ASET']['Harvest and Planting']-1]=d['A All'][iT,meta['LUT']['Derived']['ASET']['Harvest and Planting']-1]+ind[0].size
		zType0[ind]=meta['LUT']['Derived']['ASET']['Harvest and Planting']

		# Even if harvest is not the last disturbance, assume the rest is licensee
		# harvesting
		ind=np.where( (pl_oc==1) &
			   (zType0==0) &
			   (np.isin(pl_fsc,meta['Param']['Raw']['FSC']['NO List ID'])==False) )
		d['A All'][iT,meta['LUT']['Derived']['ASET']['Harvest and Planting']-1]=d['A All'][iT,meta['LUT']['Derived']['ASET']['Harvest and Planting']-1]+ind[0].size
		zType0[ind]=meta['LUT']['Derived']['ASET']['Harvest and Planting']

		#----------------------------------------------------------------------
		# Unknown
		#----------------------------------------------------------------------
		ind=np.where( (pl_oc==1) &
			   (zType0==0) &
			   (np.isin(pl_fsc,meta['Param']['Raw']['FSC']['NO List ID'])==True) )
		d['A NOSE'][iT,meta['LUT']['Derived']['ASET']['Unknown']-1]=ind[0].size
		ind=np.where( (pl_oc==1) &
			   (zType0==0) )
		d['A All'][iT,meta['LUT']['Derived']['ASET']['Unknown']-1]=ind[0].size
		zType0[ind]=meta['LUT']['Derived']['ASET']['Unknown']

		#----------------------------------------------------------------------
		# Direct seeding
		#----------------------------------------------------------------------
		for iDS in range(3):
			ind=np.where( (zDS[iDS]['Year']==d['tv'][iT]) & (np.isin(pl_fsc,meta['Param']['Raw']['FSC']['NO List ID'])==True) )
			d['A NOSE'][iT,meta['LUT']['Derived']['ASET']['Direct Seeding']-1]=ind[0].size
			ind=np.where( (zDS[iDS]['Year']==d['tv'][iT]) )
			d['A All'][iT,meta['LUT']['Derived']['ASET']['Direct Seeding']-1]=ind[0].size
			zType0[ind]=meta['LUT']['Derived']['ASET']['Direct Seeding']

		# Update last disturbance if planting occurs
		ind=np.where( (pl_oc==1) & (pl_doy>fire_doy) & (pl_doy>kd_doy) )
		LastDisturbance[ind]=meta['LUT']['Event']['Planting']

		# Package
		for iP in range(6):
			ind=np.where( (zP[iP]['Year']==d['tv'][iT]) )
			zP[iP]['Type'][ind]=zType0[ind]

	# Save summary
	gu.opickle(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\ASET_SummaryByTime.pkl',d)

	# Save summary to spreadsheet
	time=str(datetime.today().year) + '-' + str(datetime.today().month) + '-' + str(datetime.today().day)
	ptL=list(meta['LUT']['Derived']['ASET'].keys())
	d0={'Year':d['tv']}
	for i in range(d['A All'].shape[1]):
		d0[ptL[i]]=d['A All'][:,i]
	df=pd.DataFrame.from_dict(d0)
	df.to_excel(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\ASET_SummaryByTime_All_' + time + '.xlsx',index=False)
	d0={'Year':d['tv']}
	for i in range(d['A NOSE'].shape[1]):
		d0[ptL[i]]=d['A NOSE'][:,i]
	df=pd.DataFrame.from_dict(d0)
	df.to_excel(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\ASET_SummaryByTime_NOSE_' + time + '.xlsx',index=False)

	# Save packed rasters
	for i in range(6):
		z1=copy.deepcopy(zRef)
		z1['Data']=zP[i]['Type']
		gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(i+1) + '_ASET.tif')

	return

#%%
def RasterizeASET(meta,YearLast):
	# Import data
	zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
	zH=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VEG_CONSOLIDATED_CUT_BLOCKS_SP\\HARVEST_YEAR_YearLast.tif')
	zPD=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Disturbances\\HarvestPercentDead_FromCruise.tif')

	# All stand establishment
	z1=copy.deepcopy(zRef)
	z1['Data']=np.zeros(zRef['Data'].shape,dtype='int16')
	for i in range(6):
		zY=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(i+1) + '_Year.tif')
		zT=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(i+1) + '_ASET.tif')
		ind=np.where(zY['Data']>0)
		z1['Data'][ind]=zT['Data'][ind]
	# Add salvage
	ind=np.where( (zH['Data']>=2017) & (zPD['Data']>=15) & (z1['Data']!=meta['LUT']['Derived']['ASET']['Replanting']) )
	z1['Data'][ind]=meta['LUT']['Derived']['ASET']['Salvage and Planting Post Beetle']
	gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_ASET_All.tif')

	# All stand establishment (post 2017)
	z1=copy.deepcopy(zRef)
	z1['Data']=np.zeros(zRef['Data'].shape,dtype='int16')
	for i in range(6):
		zY=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(i+1) + '_Year.tif')
		zT=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(i+1) + '_ASET.tif')
		ind=np.where(zY['Data']>=2017)
		z1['Data'][ind]=zT['Data'][ind]
	# Add salvage
	ind=np.where( (zH['Data']>=2017) & (zPD['Data']>=15) & (z1['Data']!=meta['LUT']['Derived']['ASET']['Replanting']) )
	z1['Data'][ind]=meta['LUT']['Derived']['ASET']['Salvage and Planting Post Beetle']
	gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_ASET_Post2017.tif')

	# Non-obligation stand establishment
	z1=copy.deepcopy(zRef)
	z1['Data']=np.zeros(zRef['Data'].shape,dtype='int16')
	z2=copy.deepcopy(zRef)
	z2['Data']=np.zeros(zRef['Data'].shape,dtype='int16')
	for i in range(6):
		zY=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(i+1) + '_Year.tif')
		zT=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(i+1) + '_ASET.tif')
		zFSC=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(i+1) + '_SILV_FUND_SOURCE_CODE.tif')
		ind=np.where( (zY['Data']>0) & (np.isin(zFSC['Data'],meta['Param']['Raw']['FSC']['NO List ID'])==True) )
		z1['Data'][ind]=zT['Data'][ind]
		ind=np.where( (zY['Data']>0) & (np.isin(zFSC['Data'],meta['Param']['Raw']['FSC']['NO List ID'])==True) & (zY['Data']==YearLast) )
		z2['Data'][ind]=zT['Data'][ind]
	gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_NonOb_ASET_All.tif')
	gis.SaveGeoTiff(z2,meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_NonOb_ASET_' + str(YearLast) + '.tif')

	return

#%%
def DerivePlantingStatsByTime(meta):
	# ASET is not included here because a summary time series is created during production

	YearLast=2023

	zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
	iMask=np.where(zRef['Data']==1)
	zBGC=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\BEC_ZONE_CODE.tif')['Data'][iMask]

	ds={}
	ds['Year']=np.arange(1960,YearLast+1,1)
	ds['Global']={}
	ds['Global']['GW Mean']=np.zeros(ds['Year'].size)
	ds['Global']['SPH Mean']=np.zeros(ds['Year'].size)
	ds['ByFSC']={}
	ds['ByFSC']['FSC']=np.array(list(meta['LUT']['RSLT_ACTIVITY_TREATMENT_SVW']['SILV_FUND_SOURCE_CODE'].keys()))
	ds['ByFSC']['Area']=np.zeros((ds['Year'].size,ds['ByFSC']['FSC'].size))

	ds['ByBGC']={}
	ds['ByBGC']['ZONE']=np.array(list(meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE'].keys()))
	ds['ByBGC']['GW Mean']=np.zeros((ds['Year'].size,ds['ByBGC']['ZONE'].size))
	ds['ByBGC']['SPH Mean']=np.zeros((ds['Year'].size,ds['ByBGC']['ZONE'].size))

	# Need to track average over 6 compacted files and then calculate weighted average at the end
	trackG=[None]*6
	for i in range(6):
		trackG[i]={}
		trackG[i]['GW Mean']=np.zeros(ds['Year'].size)
		trackG[i]['SPH Mean']=np.zeros(ds['Year'].size)
		trackG[i]['N']=np.zeros(ds['Year'].size)
	trackBGC=[None]*6
	for i in range(6):
		trackBGC[i]={}
		trackBGC[i]['GW Mean']=np.zeros((ds['Year'].size,ds['ByBGC']['ZONE'].size))
		trackBGC[i]['SPH Mean']=np.zeros((ds['Year'].size,ds['ByBGC']['ZONE'].size))
		trackBGC[i]['N']=np.zeros((ds['Year'].size,ds['ByBGC']['ZONE'].size))

	for i in range(6):
		print(i)
		z={}
		z['yr']=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(i+1) + '_Year.tif')['Data'][iMask]
		ikp=np.where(z['yr']>0)
		z['yr']=z['yr'][ikp]
		z['BGC']=zBGC[ikp]
		z['FSC']=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(i+1) + '_SILV_FUND_SOURCE_CODE.tif')['Data'][iMask][ikp]
		z['SPH']=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(i+1) + '_SPH_Planted.tif')['Data'][iMask][ikp]
		z['S1']=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(i+1) + '_PL_SPECIES_CD1.tif')['Data'][iMask][ikp]
		z['P1']=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(i+1) + '_PL_SPECIES_PCT1.tif')['Data'][iMask][ikp]
		z['GW1']=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(i+1) + '_PL_SPECIES_GW1.tif')['Data'][iMask][ikp]
		z['S2']=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(i+1) + '_PL_SPECIES_CD2.tif')['Data'][iMask][ikp]
		z['P2']=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(i+1) + '_PL_SPECIES_PCT2.tif')['Data'][iMask][ikp]
		z['GW2']=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(i+1) + '_PL_SPECIES_GW2.tif')['Data'][iMask][ikp]
		z['S3']=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(i+1) + '_PL_SPECIES_CD3.tif')['Data'][iMask][ikp]
		z['P3']=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(i+1) + '_PL_SPECIES_PCT3.tif')['Data'][iMask][ikp]
		z['GW3']=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(i+1) + '_PL_SPECIES_GW3.tif')['Data'][iMask][ikp]
		z['S4']=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(i+1) + '_PL_SPECIES_CD4.tif')['Data'][iMask][ikp]
		z['P4']=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(i+1) + '_PL_SPECIES_PCT4.tif')['Data'][iMask][ikp]
		z['GW4']=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(i+1) + '_PL_SPECIES_GW4.tif')['Data'][iMask][ikp]

		for iT in range(ds['Year'].size):
			indT=np.where( (z['yr']==ds['Year'][iT]) )[0]

			# By Funding Source
			fsc0=z['FSC'][indT]
			uFSC=np.unique(fsc0)
			uFSC=uFSC[uFSC>0]
			for iFSC in range(uFSC.size):
				ind1=np.where( (fsc0==uFSC[iFSC]) )[0]
				ind2=np.where(ds['ByFSC']['FSC']==u1ha.lut_n2s(meta['LUT']['RSLT_ACTIVITY_TREATMENT_SVW']['SILV_FUND_SOURCE_CODE'],uFSC[iFSC])[0])[0]
				ds['ByFSC']['Area'][iT,ind2]=ds['ByFSC']['Area'][iT,ind2]+ind1.size

			# Global
			f1=z['P1'][indT].astype('float')/100
			f2=z['P2'][indT].astype('float')/100
			f3=z['P3'][indT].astype('float')/100
			f4=z['P4'][indT].astype('float')/100
			gw1=np.nan_to_num(z['GW1'][indT].astype('float'))
			gw2=np.nan_to_num(z['GW2'][indT].astype('float'))
			gw3=np.nan_to_num(z['GW3'][indT].astype('float'))
			gw4=np.nan_to_num(z['GW4'][indT].astype('float'))
			trackG[i]['GW Mean'][iT]=np.nanmean(( (f1*gw1)+(f2*gw2)+(f3*gw3)+(f4*gw4) )/(f1+f2+f3+f4))
			trackG[i]['SPH Mean'][iT]=np.nanmean(z['SPH'][indT].astype('float'))
			trackG[i]['N'][iT]=indT.size

			# By BGC zone
			bgc0=z['BGC'][indT]
			sph0=z['SPH'][indT]
			uBGC=np.unique(bgc0)
			uBGC=uBGC[uBGC>0]
			for iBGC in range(uBGC.size):
				ind1=np.where( (bgc0==uBGC[iBGC]) )[0]
				ind2=np.where(ds['ByBGC']['ZONE']==u1ha.lut_n2s(meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE'],uBGC[iBGC])[0])[0]
				trackBGC[i]['GW Mean'][iT,ind2]=np.mean(( (f1[ind1]*gw1[ind1])+(f2[ind1]*gw2[ind1])+(f3[ind1]*gw3[ind1])+(f4[ind1]*gw4[ind1]) )/(f1[ind1]+f2[ind1]+f3[ind1]+f4[ind1]))
				trackBGC[i]['SPH Mean'][iT,ind2]=np.mean(sph0[ind1].astype('float'))
				trackBGC[i]['N'][iT,ind2]=ind1.size

	# Calculate global mean genetic worth
	gw_sum=np.zeros(ds['Year'].size)
	sph_sum=np.zeros(ds['Year'].size)
	n=np.zeros(ds['Year'].size)
	for i in range(6):
		gw_sum=gw_sum+trackG[i]['N']*trackG[i]['GW Mean']
		sph_sum=sph_sum+trackG[i]['N']*trackG[i]['SPH Mean']
		n=n+trackG[i]['N']
	ds['Global']['GW Mean']=gw_sum/n
	ds['Global']['SPH Mean']=sph_sum/n

	# Calculate mean genetic worth by BGC zone
	gw_sum=np.zeros((ds['Year'].size,ds['ByBGC']['ZONE'].size))
	sph_sum=np.zeros((ds['Year'].size,ds['ByBGC']['ZONE'].size))
	n=np.zeros((ds['Year'].size,ds['ByBGC']['ZONE'].size))
	for i in range(6):
		gw_sum=gw_sum+trackBGC[i]['N']*trackBGC[i]['GW Mean']
		sph_sum=sph_sum+trackBGC[i]['N']*trackBGC[i]['SPH Mean']
		n=n+trackBGC[i]['N']
	ds['ByBGC']['GW Mean']=gw_sum/n
	ds['ByBGC']['SPH Mean']=sph_sum/n

	flg=0
	if flg==1:
		iT=np.where( (ds['Year']>=1985) )[0]
		plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(11,6.5)); ms=2; lw=0.5
		ax.plot(ds['Year'][iT],ds['Global']['GW Mean'][iT],'-ko',ms=ms,lw=lw,label='Province')
		ind=np.where(ds['ByBGC']['ZONE']=='CWH')[0]
		iGF=np.where(np.isnan(ds['ByBGC']['GW Mean'][iT,ind])==False)[0]
		y=np.interp(ds['Year'][iT],ds['Year'][iT[iGF]],ds['ByBGC']['GW Mean'][iT[iGF],ind])
		#ax.plot(ds['Year'][iT],y,'-gs',ms=ms,lw=lw,label='CWH')
		ind=np.where(ds['ByBGC']['ZONE']=='SBS')[0]
		#ax.plot(ds['Year'][iT],ds['ByBGC']['GW Mean'][iT,ind],'-c^',ms=ms,lw=lw,label='SBS')
		ax.set(xticks=np.arange(0,3000,5),yticks=np.arange(0,3000,2),ylabel='Average genetic worth (%)',xlabel='Time, years',xlim=[1984,2024],ylim=[0,18])
		ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
		ax.legend(loc='upper left',facecolor=[1,1,1],frameon=False);
		gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Reforestation\GeneticWorth_ts','png',900)

		iT=np.where( (ds['Year']>=1985) )[0]
		plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(11,6.5)); ms=2; lw=0.5
		ax.plot(ds['Year'][iT],ds['Global']['SPH Mean'][iT],'-ko',ms=ms,lw=lw,label='Province')
		ind=np.where(ds['ByBGC']['ZONE']=='CWH')[0]
		iGF=np.where(np.isnan(ds['ByBGC']['SPH Mean'][iT,ind])==False)[0]
		y=np.interp(ds['Year'][iT],ds['Year'][iT[iGF]],ds['ByBGC']['SPH Mean'][iT[iGF],ind])
		#np.mean(y[-10:])
		ax.plot(ds['Year'][iT],y,'-gs',ms=ms,lw=lw,label='CWH')
		ind=np.where(ds['ByBGC']['ZONE']=='SBS')[0]
		#np.mean(ds['ByBGC']['SPH Mean'][iT,ind][-10:])
		ax.plot(ds['Year'][iT],ds['ByBGC']['SPH Mean'][iT,ind],'-c^',ms=ms,lw=lw,label='SBS')
		ax.set(xticks=np.arange(0,3000,5),yticks=np.arange(0,3000,100),ylabel='Planting density (SPH)',xlabel='Time, years',xlim=[1984,2024],ylim=[0,1600])
		ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
		ax.legend(loc='lower right',facecolor=[1,1,1],frameon=False);
		gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Reforestation\PlantingDensity_ts','png',900)

	# Remove FSCs with no activity
	ind=np.where(np.sum(ds['ByFSC']['Area'],axis=0)>0)[0]
	ds['Area Unique']=ds['ByFSC']['Area'][:,ind]
	ds['FSC Unique']=ds['ByFSC']['FSC'][ind]

	# Save
	gu.opickle(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\Planting_SummaryByTime.pkl',ds)
	return

#%% Calculate time series summaries for insect comp 1
def CalcInsectComp1_TimeSeries(meta):
	zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
	tv=np.arange(1951,2030,1)
	id=meta['Param']['Raw']['InsectComp1']['ID']
	name=meta['Param']['Raw']['InsectComp1']['Insect Name']
	sev=meta['Param']['Raw']['InsectComp1']['Severity Code']
	uNam=np.unique(name)
	sevL=np.array(['L','M','S','V'])
	d={'tv':tv,'data':{}}
	for nam in uNam:
		d['data'][nam]={}
		for s in sevL:
			d['data'][nam][s]=np.zeros(tv.size)
	N_Year=10
	for iY in range(N_Year):
		yr0=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\PEST_INFESTATION_POLY\\InsectComp1_' + str(iY+1) + '_Year.tif')['Data']
		id0=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\PEST_INFESTATION_POLY\\InsectComp1_' + str(iY+1) + '_Type.tif')['Data']
		ind=np.where(yr0>0)
		yr0=yr0[ind]
		id0=id0[ind]
		uID=np.unique(id0)
		for u in uID:
			iID=np.where(id==u)[0]
			iSev=np.where( sevL==sev[iID[0]] )[0]
			ind1=np.where(id0==u)[0]
			yr1=yr0[ind1]
			idx=gu.IndicesFromUniqueArrayValues(yr1)
			for i in idx.keys():
				iT=np.where(tv==i)[0]
				d['data'][ name[iID[0]] ][ sevL[iSev[0]] ][iT]=d['data'][ name[iID[0]] ][ sevL[iSev[0]] ][iT]+idx[i].size
	gu.opickle(meta['Paths']['bc1ha'] + '\\PEST_INFESTATION_POLY\\InsectComp1_SummaryTS.pkl',d)

	flg=0
	if flg==1:
		tv=np.arange(1951,2024,1)
		d=gu.ipickle(meta['Paths']['bc1ha'] + '\\PEST_INFESTATION_POLY\\InsectComp1_SummaryTS.pkl')
		plt.plot(d['tv'],d['data']['Mountain Pine Beetle']['M'],'ko')

	return

#%% Calculate map of last year for each severe event from insect comp 1
def CalcInsectComp1_YearLast(meta):
	zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
	iL=['IBM','IBS','IBB','iBD','IDL','IDW']
	z1={}
	for i in iL:
		z1[i]=copy.deepcopy(zRef)
		z1[i]['Data']=0*z1[i]['Data']
		z1[i]['Data']=z1[i]['Data'].astype('int16')
	N_Year=10
	for iY in range(N_Year):
		print(iY)
		yr0=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\PEST_INFESTATION_POLY\\InsectComp1_' + str(iY+1) + '_Year.tif')['Data']
		id0=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\PEST_INFESTATION_POLY\\InsectComp1_' + str(iY+1) + '_Type.tif')['Data']
		for i in iL:
			ind0=np.where( (meta['Param']['Raw']['InsectComp1']['Severity Code']=='S') & (meta['Param']['Raw']['InsectComp1']['Insect Code']==i) )[0]
			if ind0.size==0:
				continue
			id=meta['Param']['Raw']['InsectComp1']['ID'][ind0[0]]
			ind1=np.where( (yr0>0) & (id0==id) )
			z1[i]['Data'][ind1]=yr0[ind1]
	#plt.matshow(z1['IDW']['Data'],clim=[1980,2022])
	for i in iL:
		gis.SaveGeoTiff(z1[i],meta['Paths']['bc1ha'] + '\\PEST_INFESTATIOn_POLY\\InsectComp1_' + i + '_YearLast_Severe.tif')

	return

#%% Generate a random sample of points from areas with disagreement in forest
def GenerateSampleForGroundTruthing(meta):
	import pyproj
	import geopandas as gpd
	from shapely.geometry import Polygon,Point
	
	ind=np.where( (z['lc2']['Data']!=meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_2']['T']) & \
				 (z['lcc_cec_2020']['Data']==meta['LUT']['Derived']['lcc_cec_c']['Forest']) & \
				 (np.isin(z['lcc_ntem_2019']['Data'],[meta['LUT']['Derived']['lcc_ntem']['Coniferous'],meta['LUT']['Derived']['lcc_ntem']['Broadleaf'],meta['LUT']['Derived']['lcc_ntem']['Mixedwood'],meta['LUT']['Derived']['lcc_ntem']['Wetland-treed']])==True) )
	
	srs=gis.ImportSRSs()
	
	lon,lat=pyproj.transform(srs['Proj']['BC1ha'],srs['Proj']['Geographic'],zRef['X'][ind],zRef['Y'][ind])
	lon_s=lon[0::30000]
	lat_s=lat[0::30000]
	print(lon_s.size)
	
	flg=1
	if flg==1:
		points=[]
		for k in range(lon_s.size):
			points.append(Point(lon_s[k],lat_s[k]))
		gdf_xy=gpd.GeoDataFrame({'geometry':points,'ID':np.arange(lon_s.size,dtype='int32')})
		#gdf_xy.crs=gdf_bc_boundary.crs
		gdf_xy.to_file(meta['Paths']['bc1ha'] + '\\forest.geojson',driver='GeoJSON')
	return

#%%
def CalculateStatsByTDC(meta):
	vList=['lc_comp1_2019','tdc','crownc']
	z0=Import_Raster(meta,[],vList)

	cdL=list(meta['LUT']['Derived']['tdc'].keys())
	vL=['Frequency (%)','Crown cover from VRI (%)','Tree density from plots (stems/ha)','Tree biomass from plots (MgC/ha)','Net growth from plots (MgC/ha/yr)']

	d={}
	for v in vL:
		d[v]={}
		for cd in cdL:
			d[v][cd]=0.0

	A_tot=np.where( (z0['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Forest'])  )[0].size
	for cd in cdL:
		ind=np.where( (z0['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (z0['tdc']['Data']==meta['LUT']['Derived']['tdc'][cd]) )
		d['Frequency (%)'][cd]=np.round(ind[0].size/A_tot*100,decimals=1)
		d['Crown cover from VRI (%)'][cd]=np.round(np.mean(z0['crownc']['Data'][ind]),decimals=1)

	# PSP stuff moved to PSP summary script
	meta,gpt,soc=ufp.ImportGroundPlotData(meta,type='Stand',include_soil='True')
	vList=['tdc']
	roi={'points':{'x':gpt['X'],'y':gpt['Y']}}
	z=Import_Raster(meta,roi,vList)
	for cd in cdL:
		ind=np.where( (gpt['PTF CN']==1) & (z['tdc']==meta['LUT']['Derived']['tdc'][cd]) )[0]
		d['Tree density from plots (stems/ha)'][cd]=int(np.nanmean(gpt['N L t0'][ind]))
		d['Tree biomass from plots (MgC/ha)'][cd]=np.round(np.nanmean(gpt['Ctot L t0'][ind]),decimals=0)
		d['Net growth from plots (MgC/ha/yr)'][cd]=np.round(np.nanmean(gpt['Ctot Net'][ind]),decimals=1)

	# OAF1 (systematic effect
	#vn='OAF1 (%)'
	#d[vn]={'Sparse':45,'Open':30,'Dense':15}
	
	# Add labels and convert to dataframe
	
	#lab=[' (61-100%)',' (26-60%)',' (10-25%)']
	#df=pd.DataFrame(d)
	df=pd.DataFrame.from_dict(d).T
	#df=df.reset_index()
	#df=df.rename(columns={'index':'Tree Density Class'})
	df.to_excel(meta['Paths']['Model']['Parameters'] + '\\Parameters_ByTreeDensityClass.xlsx')
	return df

#%% Species groups
def DeriveSpeciesGroups(meta):
	z0=u1ha.Import_Raster(meta,[],['refg','lc_comp1_2019'])
	z1=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\SPECIES_CD_1.tif')['Data']
	z2=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\SPECIES_CD_2.tif')['Data']
	z3=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\SPECIES_CD_3.tif')['Data']
	z4=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\SPECIES_CD_4.tif')['Data']
	ikp=np.where( (z0['lc_comp1_2019']['Data']==1) & (z1>0) )
	X=np.flip(np.array([z1[ikp],z2[ikp],z3[ikp],z4[ikp]],dtype='int8').T,axis=1)
	u,cnts=np.unique(X,axis=0,return_counts=True)

	#zF=0*z0['refg']['Data'].copy()
	#ind=np.where(z1==meta['LUT']['VEG_COMP_LYR_R1_POLY']['SPECIES_CD_1']['AT'])
	#zF[ind]=1
	#plt.matshow(zF)
	#ind[0].size

	ord=np.argsort(cnts)
	cntss=np.flip(cnts[ord])
	us=np.flip(u[ord,:])
	
	n=100
	cds=np.array(['empty' for _ in range(n)],dtype=object)
	for i in range(n):
		s=''
		for j in range(4):
			s0=u1ha.lut_n2s(meta['LUT']['VEG_COMP_LYR_R1_POLY']['SPECIES_CD_1'],us[i,j])
			if s0.size>0:
				if j==0:
					s=s+s0[0]
				else:
					s=s+'-'+s0[0]
		cds[i]=s
	d={}
	d['Species Group']=cds
	d['Frequency (%)']=cntss[0:n]/ikp[0].size*100
	df=pd.DataFrame.from_dict(d)
	df.to_excel(r'C:\Data\Species Composition\SpeciesGroups.xlsx')
	return

#%% Deciduous fraction
def DeriveBroadleafDeciduousFraction(meta):
	zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
	vList=['spc1_vri23','spc2_vri23','spc3_vri23','spc4_vri23','spc5_vri23','spc6_vri23','spc1_pct_vri23','spc2_pct_vri23','spc3_pct_vri23','spc4_pct_vri23','spc5_pct_vri23','spc6_pct_vri23']
	z0=u1ha.Import_Raster(meta,[],vList,'Extract Grid')
	sL=['A','AC','ACB','ACT','AD','AT','AX','D','DG','DM','DR','E','EA','EB','EE','EP','ES','EW','EXP','M','MB','MR','MV','Q','QG','RA','W','WA','WB','WD','WP','WS','WT']
	sidL=np.zeros(len(sL))
	for i in range(len(sL)):
		sidL[i]=meta['LUT']['VEG_COMP_LYR_R1_POLY']['SPECIES_CD_1'][sL[i]]
	z1=copy.deepcopy(zRef)
	ind=np.where( (np.isin(z0['spc1_vri23'],sidL)==True) ); z1['Data'][ind]=z1['Data'][ind]+z0['spc1_pct_vri23'][ind]
	ind=np.where( (np.isin(z0['spc2_vri23'],sidL)==True) ); z1['Data'][ind]=z1['Data'][ind]+z0['spc2_pct_vri23'][ind]
	ind=np.where( (np.isin(z0['spc3_vri23'],sidL)==True) ); z1['Data'][ind]=z1['Data'][ind]+z0['spc3_pct_vri23'][ind]
	ind=np.where( (np.isin(z0['spc4_vri23'],sidL)==True) ); z1['Data'][ind]=z1['Data'][ind]+z0['spc4_pct_vri23'][ind]
	ind=np.where( (np.isin(z0['spc5_vri23'],sidL)==True) ); z1['Data'][ind]=z1['Data'][ind]+z0['spc5_pct_vri23'][ind]
	ind=np.where( (np.isin(z0['spc6_vri23'],sidL)==True) ); z1['Data'][ind]=z1['Data'][ind]+z0['spc6_pct_vri23'][ind]
	#plt.close('all'); plt.matshow(z1['Data'],clim=[0,100])
	gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\VRI 2023\\BroadleafDeciduousFraction.tif')

	# Future scenario
	vList=['lu_comp1_2049s4']
	z0=u1ha.Import_Raster(meta,[],vList,'Extract Grid')
	ind=np.where(z0['lu_comp1_2049s4']==meta['LUT']['Derived']['lu_comp1']['Bioenergy Production'])
	z1['Data'][ind]=85
	gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\LandCoverUse\\BDF_Comp1_2049_Scn4.tif')

	return

#%% Access zones
def DeriveAccessZones(meta):
	zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
	z=Import_Raster(meta,[],['lc_comp1_2019','lu_comp1_2019','d2road','d2fac','slope'],'Extract Grid')
	
	slp_th=20
	tL=[meta['LUT']['Derived']['lu_comp1']['Timber Production Passive'],meta['LUT']['Derived']['lu_comp1']['Timber Production Intense']]
	z1=copy.deepcopy(zRef)
	z1['Data']=4*z1['Data']
	ind=np.where( (z['lc_comp1_2019']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (np.isin(z['lu_comp1_2019'],tL)==True) & (z['d2road']<60) & (z['slope']<slp_th) ); z1['Data'][ind]=1
	ind=np.where( (z['lc_comp1_2019']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (np.isin(z['lu_comp1_2019'],tL)==True) & (z['d2road']<30) & (z['d2fac']<50) & (z['slope']<slp_th) ); z1['Data'][ind]=2
	ind=np.where( (z['lc_comp1_2019']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (np.isin(z['lu_comp1_2019'],tL)==True) & (z['d2road']<10) & (z['d2fac']<25) & (z['slope']<slp_th) ); z1['Data'][ind]=3
	#plt.close('all'); plt.matshow(z1['Data'])
	gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\Management\\AccessZones.tif')
	return

#%%
def Calc_AIL_Silv_ByFSC(meta,cd,YearLast):
	# cd: FE-CA
	#YearLast=2023

	zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
	d={}
	d['Year']=np.arange(1970,YearLast+1,1)
	d['FSC']=np.unique(np.array(list(meta['LUT']['RSLT_ACTIVITY_TREATMENT_SVW']['SILV_FUND_SOURCE_CODE'].keys())))
	d['Area']=np.zeros( (d['Year'].size,d['FSC'].size) )
	for i in range(10):
		try:
			zYr=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + cd + '_' + str(i+1) + '_Year.tif')['Data']
			zFSC=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + cd + '_' + str(i+1) + '_SILV_FUND_SOURCE_CODE.tif')['Data']
		except:
			zYr=np.zeros(zRef['Data'].shape)

		ikp=np.where(zYr>0)
		if ikp[0].size==0:
			continue

		zYr=zYr[ikp]
		zFSC=zFSC[ikp]
		uY=gu.IndicesFromUniqueArrayValues(zYr)
		for iY in uY:
			zFSC_yr=zFSC[uY[iY]]
			uFSC=np.unique(zFSC_yr)
			uFSC=uFSC[uFSC>0]
			for iFSC in range(uFSC.size):
				ind1=np.where( (zFSC_yr==uFSC[iFSC]) )[0]
				ind2=np.where(d['FSC']==u1ha.lut_n2s(meta['LUT']['RSLT_ACTIVITY_TREATMENT_SVW']['SILV_FUND_SOURCE_CODE'],uFSC[iFSC])[0])[0]
				iT=np.where(d['Year']==iY)[0]
				d['Area'][iT,ind2]=d['Area'][iT,ind2]+ind1.size

	# Remove FSCs with no activity
	ind=np.where(np.sum(d['Area'],axis=0)>0)[0]
	d['Area']=d['Area'][:,ind]
	d['FSC']=d['FSC'][ind]

	# Summarize
	d['Area Total']=np.sum(d['Area'],axis=1)
	iNO=np.where( np.isin(d['FSC'],meta['Param']['Raw']['FSC']['NO List Name'])==True )[0]
	d['Area NO']=np.sum(d['Area'][:,iNO],axis=1)
	iL=np.where( np.isin(d['FSC'],meta['Param']['Raw']['FSC']['NO List Name'])==False )[0]
	d['Area Licensee']=np.sum(d['Area'][:,iL],axis=1)

	# Save
	gu.opickle(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + cd + '_SummaryByTimeAndFSC.pkl',d)

	flg=0
	if flg==1:
		d=gu.ipickle(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_SummaryByTimeAndFSC.pkl')
		plt.plot(d['Year'],d['Area NO']/1e3,'-k.')

	return

#%%
def Calc_AIL_Silv_ByBGC(meta,cd,YearLast):
	# cd: FE-CA
	#YearLast=2023

	zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
	zBGC=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\BEC_BIOGEOCLIMATIC_POLY\\ZONE_GapFilled.tif')['Data']

	d={}
	d['Year']=np.arange(1970,YearLast+1,1)
	d['BGC']=np.array(list(meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE'].keys()))
	d['Area']=np.zeros( (d['Year'].size,d['BGC'].size) )
	for i in range(10):
		try:
			zYr=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + cd + '_' + str(i+1) + '_Year.tif')['Data']
		except:
			zYr=np.zeros(zRef['Data'].shape)

		ikp=np.where(zYr>0)
		if ikp[0].size==0:
			continue

		zYr=zYr[ikp]
		zBGCy=zBGC[ikp]
		uYr=gu.IndicesFromUniqueArrayValues(zYr)
		for iY in uYr.keys():
			uBGC=np.unique(zBGCy[uYr[iY]])
			uBGC=uBGC[uBGC>0]
			for iBGC in range(uBGC.size):
				ind1=np.where( (zBGCy[uYr[iY]]==uBGC[iBGC]) )[0]
				ind2=np.where(d['BGC']==u1ha.lut_n2s(meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE'],uBGC[iBGC])[0])[0]
				iT=np.where(d['Year']==iY)[0]
				d['Area'][iT,ind2]=d['Area'][iT,ind2]+ind1.size

	d['A Total']=np.sum(d['Area'],axis=1)
	gu.opickle(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + cd + '_SummaryByTimeAndBGCZone.pkl',d)
	return

#%% Rasterize Watersheds (Named)
def RasterizeWatersheds(meta):
	# Import reference grid
	zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])

	pthin=meta['Paths']['GDB']['GDB'] + '\\WatershedsNamed.gdb' # fiona.listlayers(pthin)
	df=gpd.read_file(pthin,layer='FWA_NAMED_WATERSHEDS_POLY')

	# Save list without spatial data to save time
	flg=0
	if flg==1:
		df0=pd.DataFrame(df.drop(columns='geometry'))
		df0.to_excel(meta['Paths']['GDB']['GDB'] + '\\WatershedsNamed\\WatershedsNamed.xlsx')

	df=df[df.geometry!=None]
	df=df.reset_index()

	#z1=copy.deepcopy(zRef)
	#z1['Data']=np.zeros(zRef['Data'].shape,dtype='int32')
	ord=np.flip(np.arange(1,11))
	for iO in ord:
		print(iO)
		df0=df[df['STREAM_ORDER']==iO]
		df0=df0.reset_index()
		df0.to_file(meta['Paths']['GDB']['GDB'] + '\\WatershedsNamed\\WatershedsNamed_ord' + str(iO) + '.geojson',driver='GeoJSON')
		shapes=((geom,value) for geom, value in zip(df0['geometry'],df0['NAMED_WATERSHED_ID']))
		z0=np.zeros(zRef['Data'].shape,dtype=float)
		burned=features.rasterize(shapes=shapes,fill=0,out=z0,transform=zRef['Transform'])
		ind=np.where( (burned>0) )
		#z1['Data'][ind]=burned[ind]

		z2=copy.deepcopy(zRef)
		z2['Data']=np.zeros(zRef['Data'].shape,dtype='int16')
		z2['Data'][ind]=burned[ind]
		gis.SaveGeoTiff(z2,meta['Paths']['bc1ha'] + '\\Terrain\\WatershedsNamed_Order' + str(iO) + '.tif')
	#plt.close('all'); plt.matshow(z1['Data'])
	#gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\Terrain\\WatershedsNamed.tif')
	return

#%% Find Watersheds
def FindWatersheds(meta):
	# Import data
	pthin=meta['Paths']['GDB']['GDB'] + '\\WatershedsNamed.gdb'
	df=gpd.read_file(pthin,layer='FWA_NAMED_WATERSHEDS_POLY') # fiona.listlayers(pthin)
	df=df[df.geometry!=None]
	df=df.reset_index()

	# Global search
	Name='Nicola River'
	ind=np.where(df['GNIS_NAME']==Name)[0]
	df['AREA_HA'][ind]
	print(df['GNIS_ID'][ind])

	# TSA search
	tsa=gpd.read_file(meta['Paths']['GDB']['GDB'] + '\\tsa.geojson')
	tsaF=tsa[tsa['Name']=='Kootenay Lake TSA']
	tsaF=tsaF.reset_index()
	dfo=gpd.overlay(tsaF,df,how='intersection')

	# List watersheds within TSA at each stream order
	dO={}
	for Ord in range(1,11):
		ind=np.where(dfo['STREAM_ORDER']==Ord)[0]
		dO[Ord]=np.unique(dfo['GNIS_NAME'][ind])
	
	Name='Kokanee Creek'
	ind=np.where(dfo['GNIS_NAME']==Name)[0]
	dfo['AREA_HA'][ind]
	print(dfo['GNIS_ID'][ind])

	# Get IDs
	namL=['Big Creek','Gaspard Creek']
	print(dfo['GNIS_ID'][np.isin(dfo['GNIS_NAME'],namL)==True].values)

	# Get IDs
	namL=['Kokanee Creek','Redfish Creek','Laird Creek']
	print(dfo['GNIS_ID'][np.isin(dfo['GNIS_NAME'],namL)==True].values)

	return

#%% Export attribute water management to spreadsheet
def ExportWaterManagementToSpreadsheet(meta):
	d={}
	d['mmwb']=gpd.read_file(meta['Paths']['GDB']['GDB'] + '\\WaterManagement.gdb',layer='FWA_MANMADE_WATERBODIES_POLY')
	d['dams']=gpd.read_file(meta['Paths']['GDB']['GDB'] + '\\WaterManagement.gdb',layer='WRIS_DAMS_PUBLIC_SVW')
	d['flooda']=gpd.read_file(meta['Paths']['GDB']['GDB'] + '\\WaterManagement.gdb',layer='WLS_NON_TRIM_FLDAREA_LINES_SP')
	pd.DataFrame(d['dams'].drop(columns='geometry')).to_excel(meta['Paths']['GDB']['GDB'] + '\\WaterManagement_WRIS_DAMS_PUBLIC_SVW.xlsx',index=False)
	pd.DataFrame(d['mmwb'].drop(columns='geometry')).to_excel(meta['Paths']['GDB']['GDB'] + '\\WaterManagement_FWA_MANMADE_WATERBODIES_POLY.xlsx',index=False)
	pd.DataFrame(d['flooda'].drop(columns='geometry')).to_excel(meta['Paths']['GDB']['GDB'] + '\\WaterManagement_WLS_NON_TRIM_FLDAREA_LINES_SP.xlsx',index=False)
	return

#%%
def Calc_AIL_ByFSC_And_ASET(meta,YearLast):
	zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])

	d={}
	d['tv']=np.arange(2018,YearLast+1,1)
	d['Data']={}
	for fsc in meta['LUT']['RSLT_ACTIVITY_TREATMENT_SVW']['SILV_FUND_SOURCE_CODE'].keys():
		d['Data'][fsc]={}
		for aset in meta['LUT']['Derived']['ASET'].keys():
			d['Data'][fsc][aset]=np.zeros(d['tv'].size)
	
	for i in range(6):
		zY=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(i+1) + '_Year.tif')['Data']
		zT=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(i+1) + '_ASET.tif')['Data']
		zF=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(i+1) + '_SILV_FUND_SOURCE_CODE.tif')['Data']
		for iT,yr in enumerate(d['tv']):
			ind0=np.where( (zY==yr) )
			zT0=zT[ind0]
			zF0=zF[ind0]
			for fsc in meta['LUT']['RSLT_ACTIVITY_TREATMENT_SVW']['SILV_FUND_SOURCE_CODE'].keys():
				for aset in meta['LUT']['Derived']['ASET'].keys():
					ind1=np.where( (zT0==meta['LUT']['Derived']['ASET'][aset]) & \
						(zF0==meta['LUT']['RSLT_ACTIVITY_TREATMENT_SVW']['SILV_FUND_SOURCE_CODE'][fsc]) )
					d['Data'][fsc][aset][iT]=d['Data'][fsc][aset][iT]+ind1[0].size
	pthout=meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\ASET_SummaryByTimeAndFSC.pkl'
	gu.opickle(pthout,d)

	return d

#%%
# def DeriveConservationCompolation(meta):
#	 vList=['refg','lc_comp1_2019','protected','park','ogma','ogdef','parknat']
#	 z=Import_Raster(meta,[],vList)
	
#	 #zLCC1=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCoverClass1.tif')
#	 #zProt=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\TA_PROTECTED_LANDS_SV\\PROTECTED_LANDS_DESIGNATION.tif')
#	 #zPark=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\TA_PARK_ECORES_PA_SVW\\ADMIN_AREA_SID.tif')
#	 #zOGMA=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RMP_OGMA_LEGAL_ALL_SVW\\LEGAL_OGMA_PROVID.tif')
#	 #zOGD=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\OGSR_TAP_PRIORITY_DEF_AREA_SP\\PRIORITY_DEFERRAL_ID.tif')
#	 #zNP=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\CLAB_NATIONAL_PARKS\\ENGLISH_NAME.tif')
#	 #zFCR=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\RSLT_FOREST_COVER_RESERVE_SVW.tif')
	
#	 # ListCon=['No Timber Harvesting Areas','Landscape Corridors',
#	 #   'Critical Deer Winter Range','Sensitive Watershed','Water Management Units','High Value Wetlands for Moose','Telkwa Caribou Recovery Area',
#	 #   'Caribou Migration Corridor','High Biodiversity Emphasis Areas','Scenic Corridors']
	
#	 # Generate random new protected areas
#	 flg=0
#	 if flg==1:
#		 gdf=u1ha.Import_GDBs_ProvinceWide()
#		 N=2500 # Number within bounding box
#		 Dbuf=6200 # metres
	
#		 x=np.random.uniform(zRef['xmin'],zRef['xmax'],size=N)
#		 y=np.random.uniform(zRef['ymin'],zRef['ymax'],size=N)
#		 points=[]
#		 for k in range(x.size):
#			 points.append(Point(x[k],y[k]))
#		 gdf_xy=gpd.GeoDataFrame({'geometry':points,'ID':1})
#		 gdf_xy.crs=gdf['bc_bound']['gdf'].crs
	
#		 gdf_xy=gpd.sjoin(gdf_xy,gdf['bc_bound']['gdf'],op='within')
#		 #gdf_xy.plot(ax=ax[0],markersize=8)
	
#		 gdf_xyb=gdf_xy.geometry.buffer(Dbuf)
#		 gdf_xyb=gpd.GeoDataFrame({'geometry':gdf_xyb,'ID':np.arange(0,gdf_xyb.size,1)})
#		 gdf_xyb.crs=gdf['bc_bound']['gdf'].crs
#		 #gdf_xyb.plot(ax=ax[0],color='r')
	
#		 shapes=((geom,value) for geom, value in zip(gdf_xyb['geometry'],gdf_xyb['ID']))
#		 z=np.zeros(zRef['Data'].shape,dtype=float)
#		 burned=features.rasterize(shapes=shapes,fill=0,out=z,transform=zRef['Transform'])
	
#		 zRn=copy.deepcopy(zRef)
#		 zRn['Data']=z.astype('int32')
#		 gis.SaveGeoTiff(zRn,meta['Paths']['bc1ha'] + '\\LandCoverUse\\REARs_RandomAreas.tif')
	
#	 else:
#		 zRn=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\REARs_RandomAreas.tif')
	
#	 # # Protected areas with random areas
	
#	 # Area treed
#	 ind0=np.where( (zLCC1['Data']==meta['LUT']['Derived']['lc_comp1']['Forest']) );
#	 A_treed=ind0[0].size
	
#	 # Everything completed+proposed without random additions (Comp=1, Prop=2)
#	 zCP=copy.deepcopy(zRef)
#	 zCP['Data']=np.zeros(zRef['Data'].shape,dtype='int16')
#	 ind=np.where( (zLCC1['Data']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (zPark['Data']>0) |
#				   (zLCC1['Data']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (zOGMA['Data']>0) |
#				   (zLCC1['Data']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (zNP['Data']>0) |
#				   (zLCC1['Data']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (zProt['Data']>0))
#	 zCP['Data'][ind]=1
#	 ind=np.where( (zCP['Data']==0) & (zLCC1['Data']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (zOGD['Data']>0) )
#	 zCP['Data'][ind]=2
#	 ind=np.where( (zLCC1['Data']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (zCP['Data']==0) ); zCP['Data'][ind]=3
#	 ind=np.where( (zLCC1['Data']!=meta['LUT']['Derived']['lc_comp1']['Forest']) & (zCP['Data']==0) ); zCP['Data'][ind]=4
#	 ind=np.where( (zRef['Data']==0) ); zCP['Data'][ind]=5
	
#	 ind=np.where( (zLCC1['Data']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (zCP['Data']>=1) & (zCP['Data']<=2) ); A_CPsep=ind[0].size
#	 print('Area Protected from completed+proposed (Sep) (%): ' + str(A_CPsep/A_treed*100))
	
#	 zCP['Data']=zCP['Data'].astype('int16')
#	 gis.SaveGeoTiff(zCP,meta['Paths']['bc1ha'] + '\\LandCoverUse\\REARs_CompPlusProp.tif')
	
#	 # Everything completed+proposed with random additions (Comp=1, Prop=2)
#	 zCP=copy.deepcopy(zRef)
#	 zCP['Data']=np.zeros(zRef['Data'].shape,dtype='int16')
#	 ind=np.where( (zLCC1['Data']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (zPark['Data']>0) |
#				   (zLCC1['Data']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (zOGMA['Data']>0) |
#				   (zLCC1['Data']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (zNP['Data']>0) |
#				   (zLCC1['Data']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (zProt['Data']>0) )
#	 zCP['Data'][ind]=1
#	 ind=np.where( (zCP['Data']==0) & (zLCC1['Data']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (zOGD['Data']>0) | (zCP['Data']==0) & (zLCC1['Data']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (zRn['Data']>0) )
#	 zCP['Data'][ind]=2
#	 ind=np.where( (zLCC1['Data']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (zCP['Data']==0) ); zCP['Data'][ind]=3
#	 ind=np.where( (zLCC1['Data']!=meta['LUT']['Derived']['lc_comp1']['Forest']) & (zCP['Data']==0) ); zCP['Data'][ind]=4
#	 ind=np.where( (zRef['Data']==0) ); zCP['Data'][ind]=5
	
#	 ind=np.where( (zLCC1['Data']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (zCP['Data']>=1) & (zCP['Data']<=2) ); A_CPsep=ind[0].size
#	 print('Area Protected from completed+proposed (Sep) (%): ' + str(A_CPsep/A_treed*100))
	
#	 zCP['Data']=zCP['Data'].astype('int32')
#	 gis.SaveGeoTiff(zCP,meta['Paths']['bc1ha'] + '\\LandCoverUse\\REARs_CompPlusPropWithRandomAdditions.tif')
#	 return