#%% Import modules
import os
import numpy as np
import gc
import matplotlib.pyplot as plt
import geopandas as gpd
import pandas as pd
from scipy.interpolate import griddata
import copy
from shapely.geometry import Polygon,Point,box,shape
import rasterio
from rasterio import features
from scipy.interpolate import NearestNDInterpolator
import fiona
import time
import cv2
import gzip
import netCDF4 as nc
from scipy.interpolate import griddata
import scipy.ndimage
from perlin_noise import PerlinNoise
from numpy.random import default_rng
import fcgadgets.macgyver.util_general as gu
import fcgadgets.macgyver.util_gis as gis
import fcgadgets.macgyver.util_query_gdb as qgdb
import fcgadgets.gaia.gaia_util as gaia
import fcexplore.field_plots.Processing.fp_util as ufp
import fcgadgets.cbrunner.cbrun_util as cbu

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
	meta['Paths']['GDB']['LandCover']=r'C:\Data\Geodatabases\LandCover\20230607\LandCover.gdb'
	meta['Paths']['GDB']['LandUse']=r'C:\Data\Geodatabases\LandUse\20230501\LandUse.gdb'
	meta['Paths']['GDB']['Disturbance']=r'C:\Data\Geodatabases\Disturbances\20230501\Disturbances.gdb'
	meta['Paths']['GDB']['Results']=r'C:\Data\Geodatabases\Results\20230430\Results.gdb'
	meta['Paths']['GDB']['VRI 2002']=r'C:\Data\Geodatabases\VRI\2002\VRI2002_VEG_COMP_LYR_R1_POLY_FINAL_DELIVERYV4.gdb'
	meta['Paths']['GDB']['VRI 2015']=r'C:\Data\Geodatabases\VRI\2015\VEG_COMP_LYR_R1_POLY.gdb'
	meta['Paths']['GDB']['VRI 2023']=r'C:\Data\Geodatabases\VRI\20230401\VRI.gdb'
	meta['Paths']['DB']={}
	meta['Paths']['DB']['Harvest']=r'C:\Data\Harvest'
	meta['Paths']['DB']['Waste']=r'C:\Data\Waste Wood\FromWasteSystem'
	meta['Paths']['DB']['HBS']=r'C:\Data\Harvest\HBS'
	meta['Paths']['DB']['CO2']=r'C:\Data\Carbon Dioxide Concentration\Carbon Dioxide Concentration Mid Year.xlsx'
	meta['Paths']['DB']['NDEP']=r'C:\Data\Nitrogen Deposition'
	meta['Paths']['Model']={}
	#meta['Paths']['Model']['Code']=r'C:\Code_Python\fcgadgets\cbrunner'
	meta['Paths']['Model']['Code']=r'G:\My Drive\Code_Python\fcgadgets\cbrunner'
	meta['Paths']['Model']['Parameters']=meta['Paths']['Model']['Code'] + '\\Parameters'
	meta['Paths']['Model']['Taz Datasets']=r'C:\Data\Taz Datasets'
	meta['Paths']['GP']={}
	meta['Paths']['GP']['DB']=r'C:\Data\Field Plots\PSP-NADB2'
	meta['Paths']['GP']['Raw Data']={}
	#meta['Paths']['GP']['Raw Data']['BC']=meta['Paths']['GP']['DB'] + '\\Given\BC\Received 2023-03-02'
	meta['Paths']['GP']['Raw Data']['BC']=meta['Paths']['GP']['DB'] + '\\Given\BC\\Recieved 2024-03-19'
	
	meta['Graphics']={'Plot Style':{},'Map':{},'Flowchart':{}}
	meta['Graphics']['Plot Style']='Web' # Manuscript
	meta['Graphics']['gp']=gu.SetGraphics(meta['Graphics']['Plot Style'])
	meta['Graphics']['gp']['AxesLetterStyle']='Caps'
	meta['Graphics']['gp']['AxesLetterFontWeight']='Bold'
	meta['Graphics']['Print Figures']='Off'
	meta['Graphics']['Print Figure Path']=r'C:\Users\rhember\OneDrive - Government of BC\Figures\BCFCS'
	
	# Defaults assume province-wide map
	meta['Graphics']['Map']['RGSF']=1
	meta['Graphics']['Map']['Fig Width']=9.75
	meta['Graphics']['Map']['Side Space']=0
	meta['Graphics']['Map']['Map Position']=[0,0,1-meta['Graphics']['Map']['Side Space']-0.01,1]
	meta['Graphics']['Map']['Map Axis Vis']='off'
	meta['Graphics']['Map']['Map Grid Vis']=False
	meta['Graphics']['Map']['Legend X']=0.72
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
	
	meta['Graphics']['Flowchart']={}
	meta['Graphics']['Flowchart']['Font Name']='Arial'
	meta['Graphics']['Flowchart']['Font Color']='#6b7d2a'
	meta['Graphics']['Flowchart']['Font Size']='10'
	meta['Graphics']['Flowchart']['Penwidth']='0.5'
	meta['Graphics']['Flowchart']['Node Background Color']='#f0fca2'
	meta['Graphics']['Flowchart']['Cluster Background Color']='#f4f5f2'
	meta['Graphics']['Flowchart']['Cluster Background Color 2']='#d1d1d1'
	
	meta['Graphics']['GP Comp']={}
	meta['Graphics']['GP Comp']['bl']=np.array([0.55,0.75,1])
	meta['Graphics']['GP Comp']['bd']=np.array([0.27,0.49,0.77])
	meta['Graphics']['GP Comp']['gl']=np.array([0.7,0.95,0]) # 0.85,1,0.65
	meta['Graphics']['GP Comp']['gd']=np.array([0.4,0.75,0])
	meta['Graphics']['GP Comp']['rl']=np.array([0.8,0.6,0.4])
	meta['Graphics']['GP Comp']['rd']=np.array([0.4,0.3,0.2])
	
	meta['Graphics']['Modelling']={}
	meta['Graphics']['Modelling']['AxesLetterStyle']='NoPar'
	meta['Graphics']['Modelling']['AxesFontWeight']='Bold'
	
	# Tracking parameters
	meta['Graphics']['Fig Count']=1
	meta['Graphics']['Tab Count']=1
	
	# Initiate geospatial info
	meta['Geos']={}
	
	# Import variable info
	meta['Geos']['Variable Info']=gu.ReadExcel(meta['Paths']['Model']['Code'] + '\\Parameters\\Table_BCFCS_BC1haRasterVariableList.xlsx')
	
	# Import coordinate reference system
	meta['Geos']['crs']=gpd.read_file(meta['Paths']['GDB']['GDB'] + '//lakes.geojson').crs
	
	# Import LUTs
	if 'LUT' not in meta:
		meta=ImportLUTs(meta)

	# Import model parameters
	meta=cbu.ImportParameters(meta)

	# Import field plot parameters
	if 'Include GPs' in argv:
		meta['GP'],data,soc=ufp.ImportGroundPlotData(meta,type='Just Parameters',include_soil='True')
		meta['GP']['Data']=data
		meta['GP']['SOC']=soc
		del meta['Data']

	# Climate
	meta['Climate']={'SF':{}}
	meta['Climate']['SF']['cwd']=0.01
	meta['Climate']['SF']['cmi']=0.01
	meta['Climate']['SF']['ea']=0.01
	meta['Climate']['SF']['es']=0.01
	meta['Climate']['SF']['eta']=0.01
	meta['Climate']['SF']['etp']=0.01
	meta['Climate']['SF']['prcp']=0.01
	meta['Climate']['SF']['rswd']=0.01
	meta['Climate']['SF']['tmean']=0.01
	meta['Climate']['SF']['tmax']=0.01
	meta['Climate']['SF']['tmin']=0.01
	meta['Climate']['SF']['vpd']=0.01
	meta['Climate']['SF']['ws']=0.01
	meta['Climate']['SF']['wsp']=0.01
	meta['Climate']['SF']['runoff']=0.01
	meta['Climate']['SF']['melt']=0.01
	meta['Climate']['Missing Number']=-99

	# Ancillary tables
	meta['Tables']={}
	fileL=os.listdir(meta['Paths']['Model']['Parameters'])
	for file in fileL:
		if file[0:6]=='Table_':
			meta['Tables'][file[6:-5]]=gu.ReadExcel(meta['Paths']['Model']['Parameters'] + '\\' + file)

	return meta

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
			pass

	# Override pest severity so that it is in order
	meta['LUT']['PEST_INFESTATION_POLY']['PEST_SEVERITY_CODE']={'T':1,'L':2,'M':3,'S':4,'V':5,'G':6}

	# Raw and Derived layer
	lNam='Derived'
	meta['LUT'][lNam]={}
	meta['LUT']['Raw']={}

	# Geomorphons
	d=gu.ReadExcel(meta['Paths']['Model']['Parameters'] + '\\LUT_Geomorphons.xlsx')
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
	d=gu.ReadExcel(meta['Paths']['Model']['Parameters'] + '\\LUT_lc_comp1.xlsx')
	vNam='lc_comp1'
	meta['LUT']['Raw'][vNam]=d
	meta['LUT'][lNam][vNam]={}
	for i in range(d['Value'].size):
		meta['LUT'][lNam][vNam][ d['Name'][i] ]=d['Value'][i]
	
	# Upland-wetland forest mask
	d=gu.ReadExcel(meta['Paths']['Model']['Parameters'] + '\\LUT_uplandwetland.xlsx')
	vNam='upwetf'
	meta['LUT']['Raw'][vNam]=d
	meta['LUT'][lNam][vNam]={}
	for i in range(d['Value'].size):
		meta['LUT'][lNam][vNam][ d['Name'][i] ]=d['Value'][i]

	# Land Cover Class - NTEMS
	d=gu.ReadExcel(meta['Paths']['Model']['Parameters'] + '\\LUT_lc_NTEMS.xlsx')
	vNam='lc_ntems_2019'
	meta['LUT']['Raw'][vNam]=d
	meta['LUT'][lNam][vNam]={}
	for i in range(d['Value'].size):
		meta['LUT'][lNam][vNam][ d['Name'][i] ]=d['Value'][i]

	# Land cover 2020 (CEC Given)
	d=gu.ReadExcel(meta['Paths']['Model']['Parameters'] + '\\LUT_lc_cec.xlsx')
	vNam='lc_cec'
	meta['LUT']['Raw'][vNam]=d
	meta['LUT'][lNam][vNam]={}
	for i in range(d['Value'].size):
		meta['LUT'][lNam][vNam][ d['Name'][i] ]=d['Value'][i]

	# Land cover 2020 (CEC Compressed)
	d=gu.ReadExcel(meta['Paths']['Model']['Parameters'] + '\\LUT_lc_cec_Compressed.xlsx')
	vNam='lc_cec_c'
	meta['LUT']['Raw'][vNam]=d
	meta['LUT'][lNam][vNam]={}
	for i in range(d['Value'].size):
		meta['LUT'][lNam][vNam][ d['Name'][i] ]=d['Value'][i]

	# Land Use Compilation 1
	d=gu.ReadExcel(meta['Paths']['Model']['Parameters'] + '\\LUT_lu_comp1.xlsx')
	vNam='lu_comp1'
	meta['LUT']['Raw'][vNam]=d
	meta['LUT'][lNam][vNam]={}
	for i in range(d['Value'].size):
		meta['LUT'][lNam][vNam][ d['Name'][i] ]=d['Value'][i]

	# Land Cover Land Use Compilation 1
	d=gu.ReadExcel(meta['Paths']['Model']['Parameters'] + '\\LUT_lclu_change_comp1.xlsx')
	vNam='lclu_chng_comp1'
	meta['LUT']['Raw'][vNam]=d
	meta['LUT'][lNam][vNam]={}
	for i in range(d['Value'].size):
		meta['LUT'][lNam][vNam][ d['Name'][i] ]=d['Value'][i]

	# Species leading - NTEMS
	d=gu.ReadExcel(meta['Paths']['Model']['Parameters'] + '\\LUT_spc1_NTEMS.xlsx')
	vNam='spc1_ntems'
	meta['LUT']['Raw'][vNam]=d
	meta['LUT'][lNam][vNam]={}
	for i in range(d['Value'].size):
		meta['LUT'][lNam][vNam][ d['Name'][i] ]=d['Value'][i]

	# BGC Zone / NDT Zone Combo
	d=gu.ReadExcel(meta['Paths']['Model']['Parameters'] + '\\LUT_bgc_zone.xlsx')
	vNam='bgc_zone'
	meta['LUT']['Raw'][vNam]=d

	# BGC Zone / NDT Zone Combo
	d=gu.ReadExcel(meta['Paths']['Model']['Parameters'] + '\\LUT_bgcz_ndt_combo.xlsx')
	vNam='bgc-ndt'
	meta['LUT']['Raw'][vNam]=d
	meta['LUT'][lNam][vNam]={}
	for i in range(d['ID'].size):
		meta['LUT'][lNam][vNam][ d['BGC-NDT'][i] ]=d['ID'][i]

	# Tree Density Class
	vNam='tdc'
	meta['LUT'][lNam][vNam]={'Sparse':1,'Open':2,'Dense':3}

	# Access zones
	d=gu.ReadExcel(meta['Paths']['Model']['Parameters'] + '\\LUT_access.xlsx')
	vNam='access'
	meta['LUT']['Raw'][vNam]=d
	meta['LUT'][lNam][vNam]={}
	for i in range(d['ID'].size):
		meta['LUT'][lNam][vNam][ d['Name'][i] ]=d['ID'][i]

	# Forest cover reserve compilation
	d=gu.ReadExcel(meta['Paths']['Model']['Parameters'] + '\\LUT_reserve_comp1.xlsx')
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
	d=gu.ReadExcel(meta['Paths']['Model']['Parameters'] + '\\LUT_ASET.xlsx')
	vNam='ASET'
	meta['LUT']['Raw'][vNam]=d
	meta['LUT'][lNam][vNam]={}
	for i in range(d['ID'].size):
		meta['LUT'][lNam][vNam][ d['Name'][i] ]=d['ID'][i]

	d=gu.ReadExcel(meta['Paths']['Model']['Parameters'] + '\\LUT_EcozoneCanada.xlsx')
	vNam='ezcan'
	meta['LUT']['Raw'][vNam]=d
	meta['LUT'][lNam][vNam]={}
	for i in range(d['ID'].size):
		meta['LUT'][lNam][vNam][ d['Name'][i] ]=d['ID'][i]

	# Burn severity class compelation 1
	d=gu.ReadExcel(meta['Paths']['Model']['Parameters'] + '\\LUT_burnsev_comp1.xlsx')
	vNam='burnsev_comp1'
	meta['LUT']['Raw'][vNam]=d
	meta['LUT'][lNam][vNam]={}
	for i in range(d['Value'].size):
		meta['LUT'][lNam][vNam][ d['Name'][i] ]=d['Value'][i]

	# Log grades
	d=gu.ReadExcel(meta['Paths']['Model']['Parameters'] + '\\Parameters_LogGrades.xlsx')
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

#%% Build LUTs for variables with categorical data from source geodatabases
# Takes 7.5 hours
# Raw look-up-table spreadsheets stored with model parameters, while the processed
# pickle files are are stored with bc1ha data

#fiona.listlayers(meta['Paths']['GDB']['LandCover'])
#fiona.listlayers(meta['Paths']['GDB']['LandUse'])
#fiona.listlayers(meta['Paths']['GDB']['VRI'])
#fiona.listlayers(meta['Paths']['GDB']['Disturbance'])

def BuildLUTsFromSourceGDBs(meta):
	# Unique layers
	uL=np.unique(meta['Geos']['Variable Info']['Layer Name'])

	d={}
	for iL in range(uL.size):

		if uL[iL]=='VEG_BURN_SEVERITY_SP': # 'BEC_NATURAL_DISTURBANCE_SV':
			break

		t_start=time.time()

		# Index to all variables from layer
		indL=np.where( (meta['Geos']['Variable Info']['Layer Name']==uL[iL]) & (meta['Geos']['Variable Info']['LUT Required']=='Yes') )[0]

		if indL.size==0:
			# Not categorical variables
			continue

		# Path to file
		pthL=meta['Paths']['GDB'][ meta['Geos']['Variable Info']['Geodatabase Name'][indL][0] ]

		# Initilize dictionary to store variable codes
		d[uL[iL]]={}
		for iV in range(indL.size):
			vNam=meta['Geos']['Variable Info']['Variable Name'][indL[iV]]
			d[uL[iL]][vNam]=[]

		with fiona.open(pthL,layer=uL[iL]) as source:
			for feat in source:
				prop=dict(feat['properties'].items())  
				
				for iV in range(indL.size):
					vNam=meta['Geos']['Variable Info']['Variable Name'][indL[iV]]
					# Only continue if the variable is populated
					if prop[vNam]==None:
						continue
					# Add to code list
					d[uL[iL]][vNam].append(prop[vNam])

		t_ela=time.time()-t_start
		print(t_ela)

		# Create LUT for each variable with categorical data
		lut={}
		for iV in range(indL.size):
			vNam=meta['Geos']['Variable Info']['Variable Name'][indL[iV]]
			u=np.unique(d[uL[iL]][vNam])
			lut[vNam]={}
			for iU in range(u.size):
				lut[vNam][u[iU]]=iU+1#,dtype=meta['Geos']['Variable Info']['Precision'][indL[iV]])

		# Save
		gu.opickle(meta['Paths']['bc1ha'] + '\\LUTs\\LUTs_' + uL[iL] + '.pkl',lut)

	#--------------------------------------------------------------------------
	# Standardize species codes
	#--------------------------------------------------------------------------

	# Import data
	lut_vri=gu.ipickle(meta['Paths']['bc1ha'] + '\\LUTs\\LUTs_VEG_COMP_LYR_R1_POLY.pkl')
	lut_fci=gu.ipickle(meta['Paths']['bc1ha'] + '\\LUTs\\LUTs_RSLT_FOREST_COVER_INV_SVW.pkl')
	lut_pl=gu.ipickle(meta['Paths']['bc1ha'] + '\\LUTs\\LUTs_RSLT_PLANTING_SVW.pkl')

	# Add everything to a list
	cd=list(lut_vri['SPECIES_CD_1'].keys())
	cd=cd+list(lut_vri['SPECIES_CD_2'].keys())
	cd=cd+list(lut_vri['SPECIES_CD_3'].keys())
	cd=cd+list(lut_vri['SPECIES_CD_4'].keys())
	cd=cd+list(lut_vri['SPECIES_CD_5'].keys())
	cd=cd+list(lut_vri['SPECIES_CD_6'].keys())
	cd=cd+list(lut_fci['I_SPECIES_CODE_1'].keys())
	cd=cd+list(lut_fci['I_SPECIES_CODE_2'].keys())
	cd=cd+list(lut_fci['I_SPECIES_CODE_3'].keys())
	cd=cd+list(lut_fci['I_SPECIES_CODE_4'].keys())
	cd=cd+list(lut_fci['I_SPECIES_CODE_5'].keys())
	cd=cd+list(lut_pl['SILV_TREE_SPECIES_CODE'].keys())

	# Get unique list
	uCode=np.unique(cd)

	lut_vri['SPECIES_CD_1']={}
	lut_vri['SPECIES_CD_2']={}
	lut_vri['SPECIES_CD_3']={}
	lut_vri['SPECIES_CD_4']={}
	lut_vri['SPECIES_CD_5']={}
	lut_vri['SPECIES_CD_6']={}
	for i in range(len(uCode)):
		lut_vri['SPECIES_CD_1'][uCode[i]]=i+1
		lut_vri['SPECIES_CD_2'][uCode[i]]=i+1
		lut_vri['SPECIES_CD_3'][uCode[i]]=i+1
		lut_vri['SPECIES_CD_4'][uCode[i]]=i+1
		lut_vri['SPECIES_CD_5'][uCode[i]]=i+1
		lut_vri['SPECIES_CD_6'][uCode[i]]=i+1
	gu.opickle(meta['Paths']['bc1ha'] + '\\LUTs\\LUTs_VEG_COMP_LYR_R1_POLY.pkl',lut_vri)

	lut_fci['I_SPECIES_CODE_1']={}
	lut_fci['I_SPECIES_CODE_2']={}
	lut_fci['I_SPECIES_CODE_3']={}
	lut_fci['I_SPECIES_CODE_4']={}
	lut_fci['I_SPECIES_CODE_5']={}
	for i in range(len(uCode)):
		lut_fci['I_SPECIES_CODE_1'][uCode[i]]=i+1
		lut_fci['I_SPECIES_CODE_2'][uCode[i]]=i+1
		lut_fci['I_SPECIES_CODE_3'][uCode[i]]=i+1
		lut_fci['I_SPECIES_CODE_4'][uCode[i]]=i+1
		lut_fci['I_SPECIES_CODE_5'][uCode[i]]=i+1
	gu.opickle(meta['Paths']['bc1ha'] + '\\LUTs\\LUTs_RSLT_FOREST_COVER_INV_SVW.pkl',lut_fci)

	lut_pl['SILV_TREE_SPECIES_CODE']={}
	for i in range(len(uCode)):
		lut_pl['SILV_TREE_SPECIES_CODE'][uCode[i]]=i+1
	gu.opickle(meta['Paths']['bc1ha'] + '\\LUTs\\LUTs_RSLT_PLANTING_SVW.pkl',lut_pl)

	# Export excel versions of each LUT
	for k1 in meta['LUT'].keys():
		if (k1=='Derived') | (k1=='Raw'):
			continue
		for k2 in meta['LUT'][k1].keys():
			df=pd.DataFrame({'ID':meta['LUT'][k1][k2].values(),'Name':meta['LUT'][k1][k2].keys()})
			df.to_excel(meta['Paths']['bc1ha'] + '\\LUTs\\LUT_' + k1 + '_' + k2 + '.xlsx',index=False)

	return

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
		a.to_file(meta['Paths']['GDB']['GDB'] + '\\' + v + '.geojson',driver="GeoJSON")
	return

#%% Simplify road and line geodatabases
def SimplifyRoadGDBs(meta):
	vL=['FTEN_ROAD_SEGMENT_LINES_SVW','OG_ROAD_SEGMENT_PERMIT_SP','DRP_OIL_GAS_PIPELINES_BC_SP','GBA_TRANSMISSION_LINES_SP']	
	for v in vL:
		df=gpd.read_file(meta['Paths']['GDB']['LandUse'],layer=v)
		df=df.simplify(100)
		df.to_file(meta['Paths']['GDB']['GDB'] + '\\' + v + '.geojson',driver="GeoJSON")
	return

#%% Major rivers
def GetRiversMajor(meta):
	r=gpd.read_file(meta['Paths']['GDB']['LandCover'],layer='FWA_RIVERS_POLY')
	#ord=np.flip(np.argsort(r['AREA_HA']))
	#L=list(r['GNIS_NAME_1'][ord])
	#L[0:100]
	rmL=['Fraser River','Skeena River','Columbia River','Thompson River','Peace River','Liard River','North Thompson River','South Thompson River','Lillooet River']
	rm=r[(np.isin(r['GNIS_NAME_1'],rmL)==True)]
	rm=rm.simplify(100)
	rm.to_file(meta['Paths']['GDB']['GDB'] + '\\rivermajor.geojson',driver="GeoJSON")
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

#%% Rasterize VRI
# *** Takes 7.5 hours ***
def RasterizeVRI(meta,year):

	# Pick year
	year=2023
	#year=2015
	t0=time.time()
	lNam='VEG_COMP_LYR_R1_POLY'

	# Import reference grid
	zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])

	# Import feature ID
	zFID=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI ' + str(year) + '\\vri_feaid.tif')
	fid=zFID['Data'].flatten()
	iu=gu.IndicesFromUniqueArrayValues(fid)

	# Index to variables
	indL=np.where( (meta['Geos']['Variable Info']['Layer Name']==lNam) )[0]

	# Just specific variables
	vL=['WHOLE_STEM_BIOMASS_PER_HA','BRANCH_BIOMASS_PER_HA','FOLIAGE_BIOMASS_PER_HA','BARK_BIOMASS_PER_HA']
	indL=np.where( (np.isin(meta['Geos']['Variable Info']['Variable Name'],vL)==True) )[0]
	#indL=np.where( (meta['Geos']['Variable Info']['Variable Name']=='CROWN_CLOSURE') )[0]

	# Initialize variables
	d={}
	for i in indL:
		vNam=meta['Geos']['Variable Info']['Variable Name'][i]
		Prec=meta['Geos']['Variable Info']['Precision'][i]
		Date=meta['Geos']['Variable Info']['Date Conversion Required'][i]
		if Date=='Yes':
			d[vNam + '_Year']=np.zeros(fid.size,dtype=Prec)
			d[vNam + '_Month']=np.zeros(fid.size,dtype=Prec)
		d[vNam]=np.zeros(fid.size,dtype=Prec)

	# Keep track of instances where there is no crosswalk between GDB and rasterized feature ID
	cn=0
	with fiona.open(meta['Paths']['GDB']['VRI ' + str(year)],layer=lNam) as source:
		for feat in source:
			prp=dict(feat['properties'].items())
			# if prp['BCLCS_LEVEL_2']!='T':
			#	 continue
			try:
				ind=iu[ prp['FEATURE_ID'] ]
			except:
				cn=cn+1
				print(cn)

			for i in indL:
				vNam=meta['Geos']['Variable Info']['Variable Name'][i]
				Cat=meta['Geos']['Variable Info']['LUT Required'][i]
				Date=meta['Geos']['Variable Info']['Date Conversion Required'][i]

				if prp[vNam]==None:
					continue

				if Date=='Yes':
					d[vNam + '_Year'][ind]=int(prp[vNam][0:4])
					d[vNam + '_Month'][ind]=int(prp[vNam][5:7])
				else:
					if Cat=='Yes':
						# Categorical, requires LUT
						try:
							d[vNam][ind]=meta['LUT'][lNam][vNam][ prp[vNam] ]
						except:
							continue
					else:
						# Numerical
						d[vNam][ind]=prp[vNam]

	# Save
	for i in indL:
		vNam=meta['Geos']['Variable Info']['Variable Name'][i]
		z=copy.deepcopy(zRef)
		z['Data']=np.reshape(d[vNam],zRef['Data'].shape)
		gis.SaveGeoTiff(z,meta['Paths']['bc1ha'] + '\\VRI ' + str(year) + '\\' + vNam + '.tif')

	print((time.time()-t0)/60/60)

	return

#%% Rasterize VRI 2002
# *** Takes 7.5 hours ***
def RasterizeVRI_2002(meta):
	t0=time.time()
	lNam='VEG_COMP_LYR_R1_POLY'

	# Import reference grid
	zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])

	# Import feature ID
	zFID=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2002\\featid.tif')
	fid=zFID['Data'].flatten()
	iu=gu.IndicesFromUniqueArrayValues(fid)

	# Index to variables
	indL=np.where( (meta['Geos']['Variable Info']['Layer Name']==lNam) )[0]

	# Initialize variables
	d={}
	for i in indL:
		vNam=meta['Geos']['Variable Info']['Variable Name'][i]
		Prec=meta['Geos']['Variable Info']['Precision'][i]
		Date=meta['Geos']['Variable Info']['Date Conversion Required'][i]
		if Date=='Yes':
			d[vNam + '_Year']=np.zeros(fid.size,dtype=Prec)
			d[vNam + '_Month']=np.zeros(fid.size,dtype=Prec)
		d[vNam]=np.zeros(fid.size,dtype=Prec)

	# Keep track of instances where there is no crosswalk between GDB and rasterized feature ID
	cn=0

	with fiona.open(meta['Paths']['GDB']['VRI 2002'],layer='VEG_COMP_LYR_R1_POLY_FINALV4') as source:
		for feat in source:
			prp=dict(feat['properties'].items())
			
			# if prp['BCLCS_LEVEL_2']!='T':
			#	 continue

			try:
				ind=iu[ prp['FEATURE_ID'] ]
			except:
				cn=cn+1
				print(cn)

			for i in indL:
				vNam=meta['Geos']['Variable Info']['Variable Name'][i]
				Cat=meta['Geos']['Variable Info']['LUT Required'][i]
				Date=meta['Geos']['Variable Info']['Date Conversion Required'][i]

				if vNam not in prp.keys():
					continue

				if prp[vNam]==None:
					continue

				if Date=='Yes':
					d[vNam + '_Year'][ind]=int(prp[vNam][0:4])
					d[vNam + '_Month'][ind]=int(prp[vNam][5:7])
				else:
					if Cat=='Yes':
						# Categorical, requires LUT
						try:
							d[vNam][ind]=meta['LUT'][lNam][vNam][ prp[vNam] ]
						except:
							# Older 2002 elements missing from 2023 LUT elements
							pass
					else:
						# Numerical
						d[vNam][ind]=prp[vNam]

	# Save
	for i in indL:
		vNam=meta['Geos']['Variable Info']['Variable Name'][i]
		z=copy.deepcopy(zFID)
		z['Data']=np.reshape(d[vNam],zFID['Data'].shape)
		z=gis.ClipToRaster(z,zRef)
		gis.SaveGeoTiff(z,meta['Paths']['bc1ha'] + '\\VRI 2002\\' + vNam + '.tif')
		
	print((time.time()-t0)/60/60)

	return

#%% Region of interest
def DefineROI(meta,roi,gdf):

	if roi['Type']=='Prov':
		meta['Graphics']['Map']['Legend X']=0.7
		meta['Graphics']['Map']['Show Lakes']='Off'
		meta['Graphics']['Map']['Show Rivers']='Off'
	elif roi['Type']=='LICS':
		meta['Graphics']['Map']['Legend X']=0.75

	roi['crs']=gdf['bc_bound']['gdf'].crs

	# Raster data
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
	elif roi['Type']=='Prov':
		roi['grd']=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
		roi['grd']=gis.UpdateGridCellsize(roi['grd'],meta['Graphics']['Map']['RGSF'])

	# Vector layers
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
	else:
		roi['gdf']['bound within']=roi['gdf']['bound']

	for k in gdf.keys():
		roi['gdf'][k]=gdf[k]['gdf'].cx[roi['grd']['xmin']:roi['grd']['xmax'],roi['grd']['ymin']:roi['grd']['ymax']]
		roi['gdf'][k]=roi['gdf'][k].reset_index(drop=True)
		if (roi['Type']=='ByTSA') | (roi['Type']=='ByRegDis') | (roi['Type']=='ByWatershed'):
			roi['gdf'][k]=gpd.overlay(roi['gdf'][k],roi['gdf']['bound ROI'],how='intersection')

	return roi

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
		elif nam=='road_ften':
			pthin=meta['Paths']['GDB']['GDB'] + '\\FTEN_ROAD_SEGMENT_LINES_SVW.geojson'
			df=gpd.read_file(pthin)
			roi['gdf'][nam]={}
			roi['gdf'][nam]['gdf']=df.cx[roi['grd']['xmin']:roi['grd']['xmax'],roi['grd']['ymin']:roi['grd']['ymax']]
		else:
			pass

		if (roi['Type']=='ByTSA') | (roi['Type']=='ByWatershed'):
			roi['gdf'][nam]['gdf']=gpd.overlay(roi['gdf'][nam]['gdf'],roi['gdf']['bound ROI'],how='intersection')

	return roi

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
		elif v=='aset_post18':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_ASET.tif')
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
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Disturbances\\BurnSevComp_YearLast.tif')
		elif v=='bsr_sc':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Disturbances\\BurnSevComp_SevClassLast.tif')
		elif v=='citym':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\BC_MAJOR_CITIES_POINTS_500M\\NAME.tif')
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
		elif v=='fire_yr':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\PROT_HISTORICAL_FIRE_POLYS_SP\\FIRE_YEAR_YearLast.tif')
		elif v=='fire_2023':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\PROT_HISTORICAL_FIRE_POLYS_SP\\FIRE_YEAR_2023.tif')
		elif v=='fc_yr':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_FOREST_COVER_INV_SVW\\REFERENCE_YEAR.tif')
		elif v=='fc_ssc':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_FOREST_COVER_INV_SVW\\ID_STOCKING_STATUS_CODE.tif')
		elif v=='fc_stc':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_FOREST_COVER_INV_SVW\\ID_STOCKING_TYPE_CODE.tif')
		elif v=='gfcly':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Disturbances\\GlobalForestChange_LossYear_2021.tif')
		elif v=='gfcly_filt':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Disturbances\\GlobalForestChange_LossYear_2021_Filtered.tif')
		elif v=='gsoc':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Soil\\gsoc2010_bc1ha.tif')
		elif v=='geomorph':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Terrain\\FromBurtWilliam\\geomorphons_search50.tif')
		elif v=='h_vri23':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\PROJ_HEIGHT_1.tif')
		elif v=='h_vri15':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2015\\PROJ_HEIGHT_1.tif')
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
		elif v=='kd_yr':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\SP_KD_YearLast.tif')
		elif v=='ibm_yr':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\PEST_INFESTATION_POLY\\InsectComp1_IBM_YearLast_Severe.tif')
		elif v=='ibd_yr':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\PEST_INFESTATION_POLY\\InsectComp1_IBD_YearLast_Severe.tif')
		elif v=='ibb_yr':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\PEST_INFESTATION_POLY\\InsectComp1_IBB_YearLast_Severe.tif')
		elif v=='ibs_yr':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\PEST_INFESTATION_POLY\\InsectComp1_IBS_YearLast_Severe.tif')
		elif v=='idw_yr':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\PEST_INFESTATION_POLY\\InsectComp1_IDW_YearLast_Severe.tif')
		elif v=='idl_yr':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\PEST_INFESTATION_POLY\\InsectComp1_IDL_YearLast_Severe.tif')
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
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCover_NTEMS_2019_ReclassAsComp1.tif')
		elif v=='lc_vri_recl':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCover_VRI_ReclassAsComp1.tif')
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
		elif v=='rswd_ann_n':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\Climate\\Summaries\\bc1ha_rswd_ann_norm_1971to2000.tif')
			d[v]['Data']=d[v]['Data'].astype('float')*meta['Climate']['SF']['rswd']
		elif v=='rswd_mjjas_n':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\Climate\\Summaries\\bc1ha_rswd_mjjas_norm_1971to2000.tif')
			d[v]['Data']=d[v]['Data'].astype('float')*meta['Climate']['SF']['rswd']
		elif v=='runoff_ann_n':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Climate\\Summaries\\Normals\\bc1ha_runoff_ann_norm_1971to2000.tif')
			d[v]['Data']=d[v]['Data'].astype('float')*meta['Climate']['SF']['runoff']
		elif v=='slope':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Terrain\\Slope.tif')
		elif v=='si_vri23':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\SITE_INDEX.tif')
		elif v=='slope':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Terrain\\slope.tif')
		elif v=='sp_me':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\SP-ME_MaskAll.tif')
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
		elif v=='tdc_wsg':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\TreeDensityClass_WithShrubsGrasses.tif')
		elif v=='tm':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\FTEN_CUT_BLOCK_POLY_SVW\\TIMBER_MARK.tif')
		elif v=='tmean_ann_n':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\Climate\Summaries\\BC1ha_tmean_ann_norm_1971to2000_si_hist_v1.tif')
			d[v]['Data']=d[v]['Data'].astype('float')*meta['Climate']['SF']['tmean']
		elif v=='tmean_mjjas_n':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\Climate\Summaries\\BC1ha_tmean_mjjas_norm_1971to2000_si_hist_v1.tif')
			d[v]['Data']=d[v]['Data'].astype('float')*meta['Climate']['SF']['tmean']
		elif v=='transl':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\GBA_TRANSMISSION_LINES_SP\\TRANSMISSION_LINE_ID.tif')
		elif v=='twi':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Terrain\\twi.tif')
		elif v=='upwetf_ntems':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoveruse\\UplandWetlandForest_NTEMS.tif')
		elif v=='upwetf_vri':
			d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoveruse\\UplandWetlandForest_VRI.tif')
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
				roi['grd'][v]=gis.ClipToRaster(roi['grd'][v],roi['grd'])
			return roi
		elif 'points' in roi.keys():
			return d

	else:
		return d

#%%
def GapFillBGCZ(meta):
	zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
	z=Import_Raster(meta,[],['lc_comp1_2019','bgcz'])
	ivl=5
	iGap=np.where( (zRef['Data']==1) & (z['bgcz']['Data']==0) )
	iCal=np.where( (zRef['Data'][0::ivl,0::ivl]==1) & (z['bgcz']['Data'][0::ivl,0::ivl]>0) )
	xy=np.column_stack([zRef['X'][0::ivl,0::ivl][iCal],zRef['Y'][0::ivl,0::ivl][iCal]])
	vals=z['bgcz']['Data'][0::ivl,0::ivl][iCal]
	zFill=griddata(xy,vals,(zRef['X'][iGap],zRef['Y'][iGap]),method='nearest')
	z1=copy.deepcopy(zRef)
	z1['Data']=np.zeros(zRef['Data'].shape,dtype='int8')
	z1['Data']=z['bgcz']['Data']
	z1['Data'][iGap]=zFill.astype('int8')
	z1['Data'][zRef['Data']==0]=0
	#plt.close('all'); plt.matshow(z1['Data'])
	gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\BEC_BIOGEOCLIMATIC_POLY\\ZONE_GapFilled.tif')
	return

#%%
def RasterizeWildfire(meta,zRef):
	# Rasterize historical polygons
	lNam='PROT_HISTORICAL_FIRE_POLYS_SP'
	vNam='FIRE_YEAR'
	if os.path.exists(meta['Paths']['bc1ha'] + '\\' + lNam)==False:
		os.mkdir(meta['Paths']['bc1ha'] + '\\' + lNam)
	ind=np.where( (meta['Geos']['Variable Info']['Layer Name']==lNam) & (meta['Geos']['Variable Info']['Variable Name']==vNam) )[0]
	pthin=meta['Paths']['GDB'][ meta['Geos']['Variable Info']['Geodatabase Name'][ind[0]] ]
	df=gpd.read_file(pthin,layer=lNam)
	df=df[df.geometry!=None]
	df=df.reset_index()
	zYearLast=copy.deepcopy(zRef)
	zYearLast['Data']=np.zeros(zRef['Data'].shape,dtype='int16')
	uYear=df[vNam].unique()
	tv=np.arange(np.min(uYear),np.max(uYear),1)
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

	#zYearLast=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\' + lNam + '\\' + vNam + '_YearLast.tif')

	# Get current year
	zC=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\' + lNam + '\\' + vNam + '_2022.tif')
	ind=np.where( (zYearLast['Data']==0) & (zC['Data']>0) )
	zYearLast['Data'][ind]=2022
	zC=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\' + lNam + '\\' + vNam + '_2023.tif')
	ind=np.where( (zYearLast['Data']==0) & (zC['Data']>0) )
	zYearLast['Data'][ind]=2023

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
	tv=np.arange(1920,2025,1)

	# Initialize rasters
	N_Year=6
	z={'Year':{}}
	for iY in range(N_Year):
		z['Year'][iY+1]=np.zeros(zRef['Data'].shape,dtype='int16')

	for iT in range(tv.size):
		print(tv[iT])
		try:
			z0=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\' + lNam + '\\' + vNam + '_' + str(tv[iT]) + '.tif')['Data']
		except:
			continue
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

	# Plot time series to confirm it worked
	flg=0
	if flg==1:
		lNam='PROT_HISTORICAL_FIRE_POLYS_SP'
		vNam='FIRE_YEAR'
		tv,N=u1ha.TimeSeriesOccurrenceFromPackedEvents(meta,lNam,vNam,6)
		plt.plot(tv,N,'-bo')
	return

#%% Current-year wildfire
# Current year + preveous year (sometimes missing)
def RasterizeWildfireCurrentYear(meta,yr):
	zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
	# fiona.listlayers(pthin)
	# Current year
	yr=2023
	pthin=r'C:\Data\Wildfire\Wildfire.gdb'
	df=gpd.read_file(pthin,layer='PROT_CURRENT_FIRE_POLYS_SP')
	ind=np.where(df['FIRE_YEAR']==yr)[0]
	df=df.iloc[ind]
	df=df[df.geometry!=None]
	df=df.reset_index()
	shapes=((geom,value) for geom, value in zip(df['geometry'],df['FIRE_YEAR']))
	z0=np.zeros(zRef['Data'].shape,dtype=float)
	burned=features.rasterize(shapes=shapes,fill=0,out=z0,transform=zRef['Transform'])
	z1=copy.deepcopy(zRef)
	z1['Data']=np.zeros(zRef['Data'].shape,dtype='int16')
	ind=np.where(burned>0); z1['Data'][ind]=1
	plt.close(); plt.matshow(z1['Data'])
	print(np.sum(z1['Data']))
	#gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\prot_current_fire_polys\prot_current_fire_polys.tif')
	gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\PROT_HISTORICAL_FIRE_POLYS_SP\\FIRE_YEAR_' + str(yr) + '.tif')
	
# 	# Previous year when missing
# 	import geotable
# 	import pyproj
# 	srs=gis.ImportSRSs()
# 	crs=pyproj.CRS(srs['String']['Geographic'])
# 	
# 	a=geotable.load(r'C:\Data\Wildfire\BC Fire Perimeters 2020-2022.kmz')
# 	df=gpd.GeoDataFrame(data=a.Name,geometry=a.geometry_object)
# 	df['ID']=np.ones(len(df))
# 	df.crs=pyproj.CRS(srs['String']['Geographic'])
# 	df=df.to_crs({'init':'epsg:3005'})
# 	
# 	shapes=((geom,value) for geom, value in zip(df['geometry'],df['ID']))
# 	z0=np.zeros(zRef['Data'].shape,dtype=float)
# 	burned=features.rasterize(shapes=shapes,fill=0,out=z0,transform=zRef['Transform'])
# 	
# 	# Fix by removing other years
# 	tv=np.arange(2020,2021+1,1)
# 	for iT in range(tv.size):
# 		zF0=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\PROT_HISTORICAL_FIRE_POLYS_SP\\FIRE_YEAR_' + str(tv[iT]) + '.tif')
# 		ind=np.where(zF0['Data']>0)
# 		burned[ind]=0
# 	
# 	z1=copy.deepcopy(zRef)
# 	z1['Data']=np.zeros(zRef['Data'].shape,dtype='int16')
# 	ind=np.where(burned>0); z1['Data'][ind]=1
# 	print(np.sum(z1['Data']))
# 	gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\PROT_HISTORICAL_FIRE_POLYS_SP\\FIRE_YEAR_2022.tif')
	return

#%%
def RasterizeInsects(meta,zRef):
	lNam='PEST_INFESTATION_POLY'
	vNam='PEST_SEVERITY_CODE'

	if os.path.exists(meta['Paths']['bc1ha'] + '\\' + lNam)==False:
		os.mkdir(meta['Paths']['bc1ha'] + '\\' + lNam)

	indI=np.where( (meta['Geos']['Variable Info']['Layer Name']==lNam) & (meta['Geos']['Variable Info']['Variable Name']==vNam) )[0]

	pthin=meta['Paths']['GDB'][ meta['Geos']['Variable Info']['Geodatabase Name'][indI[0]] ]

	df=gpd.read_file(pthin,layer=lNam)
	df=df[df.geometry!=None]
	df=df.reset_index()

	df=CreateIdForCategoricalVariable(meta,lNam,vNam,df)

	pestL=['IBM']#,'IBS','IBB','IBD','IDW','IDL']

	tv=np.arange(1951,2023,1)

	for pest in pestL:

		zYearLast=copy.deepcopy(zRef)
		zYearLast['Data']=np.zeros(zRef['Data'].shape,dtype='int16')

		for iT in range(tv.size):

			df0=df[ (df['PEST_SPECIES_CODE']==pest) & (df['CAPTURE_YEAR']==tv[iT]) ].copy()
			df0=df0[df0.geometry!=None]
			df0=df0.reset_index()

			if len(df0)>0:
				shapes=((geom,value) for geom, value in zip(df0.geometry,df0['ID_PEST_SEVERITY_CODE']))
				z0=np.zeros(zRef['Data'].shape,dtype=float)
				burned=features.rasterize(shapes=shapes,fill=0,out=z0,transform=zRef['Transform'])
			else:
				z0=np.zeros(zRef['Data'].shape,dtype=float)

			z1=copy.deepcopy(zRef)
			z1['Data']=z0.astype(meta['Geos']['Variable Info']['Precision'][indI[0]])
			gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\' + lNam + '\\' + vNam + '_' + pest + '_' + str(tv[iT]) + '.tif')

			# Update by year grid
			zYearLast['Data'][z0>0]=tv[iT]

		# Year of last occurrence
		z1=copy.deepcopy(zRef)
		z1['Data']=zYearLast['Data'].astype('int16')
		gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\' + lNam + '\\' + vNam + '_' + pest + '_YearLast.tif')

		# Mask of occurrence
		z1=copy.deepcopy(zRef)
		z1['Data']=np.zeros(zRef['Data'].shape,dtype='int8')
		ind0=np.where(zYearLast['Data']>0)
		z1['Data'][ind0]=1
		gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\' + lNam + '\\' + vNam + '_' + pest + '_MaskAll.tif')

		# Pack into smaller number of layers

		# Initialize rasters
		N_Year=10
		z={'Year':{},'Severity':{}}
		for iY in range(N_Year):
			z['Year'][iY+1]=np.zeros(zRef['Data'].shape,dtype='int16')
			z['Severity'][iY+1]=np.zeros(zRef['Data'].shape,dtype='int16')

		for iT in range(tv.size):
			print(tv[iT])
			z0=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\' + lNam + '\\' + vNam + '_' + pest + '_' + str(tv[iT]) + '.tif')['Data']
			ind=np.where((z['Year'][1]==0) & (z0!=0))
			z['Year'][1][ind]=tv[iT]
			z['Severity'][1][ind]=z0[ind]
			ind=np.where((z['Year'][1]!=0) & (z['Year'][1]!=tv[iT]) & (z['Year'][2]==0) & (z0!=0))
			z['Year'][2][ind]=tv[iT]
			z['Severity'][2][ind]=z0[ind]
			ind=np.where((z['Year'][2]!=0) & (z['Year'][2]!=tv[iT]) & (z['Year'][3]==0) & (z0!=0))
			z['Year'][3][ind]=tv[iT]
			z['Severity'][3][ind]=z0[ind]
			ind=np.where((z['Year'][3]!=0) & (z['Year'][3]!=tv[iT]) & (z['Year'][4]==0) & (z0!=0))
			z['Year'][4][ind]=tv[iT]
			z['Severity'][4][ind]=z0[ind]
			ind=np.where((z['Year'][4]!=0) & (z['Year'][4]!=tv[iT]) & (z['Year'][5]==0) & (z0!=0))
			z['Year'][5][ind]=tv[iT]
			z['Severity'][5][ind]=z0[ind]
			ind=np.where((z['Year'][5]!=0) & (z['Year'][5]!=tv[iT]) & (z['Year'][6]==0) & (z0!=0))
			z['Year'][6][ind]=tv[iT]
			z['Severity'][6][ind]=z0[ind]
			ind=np.where((z['Year'][6]!=0) & (z['Year'][6]!=tv[iT]) & (z['Year'][7]==0) & (z0!=0))
			z['Year'][7][ind]=tv[iT]
			z['Severity'][7][ind]=z0[ind]
			ind=np.where((z['Year'][7]!=0) & (z['Year'][7]!=tv[iT]) & (z['Year'][8]==0) & (z0!=0))
			z['Year'][8][ind]=tv[iT]
			z['Severity'][8][ind]=z0[ind]
			ind=np.where((z['Year'][8]!=0) & (z['Year'][8]!=tv[iT]) & (z['Year'][9]==0) & (z0!=0))
			z['Year'][9][ind]=tv[iT]
			z['Severity'][9][ind]=z0[ind];
			ind=np.where((z['Year'][9]!=0) & (z['Year'][9]!=tv[iT]) & (z0!=0))
			z['Year'][10][ind]=tv[iT]
			z['Severity'][10][ind]=z0[ind]

		for iY in range(N_Year):
			z1=copy.deepcopy(zRef)
			z1['Data']=z['Year'][iY+1].astype('int16')
			gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\' + lNam + '\\' + vNam + '_' + pest + '_' + str(iY+1) + '_Year.tif')
			z1=copy.deepcopy(zRef)
			z1['Data']=z['Severity'][iY+1].astype('int16')
			gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\' + lNam + '\\' + vNam + '_' + pest + '_' + str(iY+1) + '_Severity.tif')
	return

#%% Rasterize planting
def RasterizePlanting(meta):
	
	tv=np.arange(1960,2024,1)
	zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
	
	# Start with planting with spatial from RESULTS (takes 15 min)
	flg=0
	if flg==1:
		t0=time.time()
		ats={}
		ats['Path']=meta['Paths']['GDB']['Results']
		ats['Layer']='RSLT_ACTIVITY_TREATMENT_SVW'; # fiona.listlayers(ats['Path'])
		ats['crs']=meta['Geos']['crs']
		ats['Keep Geom']='On'
		ats['Select Openings']=np.array([])
		ats['SBC']=np.array(['PL'])
		ats['STC']=np.array([])
		ats['SMC']=np.array([])
		ats['FSC']=np.array([])
		ats['SOC1']=np.array([])
		ats['ROI']=[]
		ats['gdf']=qgdb.Query_Openings(ats,[])
		ats['gdf']=ats['gdf'][ats['gdf'].geometry!=None]
		ats['gdf']=ats['gdf'].reset_index()	
		ats['gdf']['Year']=np.zeros(len(ats['gdf']))
		for i in range(ats['gdf']['Year'].size):
			ats['gdf']['Year'][i]=int(ats['gdf']['ATU_COMPLETION_DATE'][i][0:4])
		AreaPlanted=ats['gdf']['ACTUAL_TREATMENT_AREA']
		NumTreesPlanted=ats['gdf']['ACTUAL_PLANTED_NUMBER']
		ats['gdf']['SPH_Planted']=NumTreesPlanted/AreaPlanted
		ats['gdf']=CreateIdForCategoricalVariable(meta,'RSLT_ACTIVITY_TREATMENT_SVW','SILV_FUND_SOURCE_CODE',ats['gdf'])
		ats['gdf']=CreateIdForCategoricalVariable(meta,'RSLT_ACTIVITY_TREATMENT_SVW','SILV_TECHNIQUE_CODE',ats['gdf'])	
		ats['gdf']=ats['gdf'][ (ats['gdf']['RESULTS_IND']=='Y') & (ats['gdf']['SILV_METHOD_CODE']!='LAYOT') ].copy()
		ats['gdf']=ats['gdf'][ats['gdf'].geometry!=None]
		ats['gdf']=ats['gdf'].reset_index()
		#ats['gdf'].to_file(r'C:\Data\BC1ha\RSLT_ACTIVITY_TREATMENT_SVW\ats.geojson',driver="GeoJSON")
	else:
		ats={}
		ats['gdf']=gpd.read_file(r'C:\Data\BC1ha\RSLT_ACTIVITY_TREATMENT_SVW\ats.geojson')

	# Add areas where FC is artificial
	flg=0
	if flg==1:
		at={}
		at['Path']=meta['Paths']['GDB']['Results']
		at['Layer']='RSLT_ACTIVITY_TREATMENT_SVW'; # fiona.listlayers(at['Path'])
		at['crs']=meta['Geos']['crs']
		at['Keep Geom']='Off'
		at['Select Openings']=np.array([])
		at['SBC']=np.array(['PL'])
		at['STC']=np.array([])
		at['SMC']=np.array([])
		at['FSC']=np.array([])
		at['SOC1']=np.array([])
		at['ROI']=[]
		at['gdf']=qgdb.Query_Openings(at,[])
		at['gdf']=CreateIdForCategoricalVariable(meta,'RSLT_ACTIVITY_TREATMENT_SVW','SILV_FUND_SOURCE_CODE',at['gdf'])
		at['gdf']=CreateIdForCategoricalVariable(meta,'RSLT_ACTIVITY_TREATMENT_SVW','SILV_TECHNIQUE_CODE',at['gdf'])
		# Make sure to remove entries that we know did not occur (planned or layout)
		ikp=np.where(  (at['gdf']['RESULTS_IND']=='Y') & (at['gdf']['SILV_METHOD_CODE']!='LAYOT') )[0]
		for k in at['gdf'].keys():
			at['gdf'][k]=at['gdf'][k][ikp]
		at['gdf']['Year']=np.zeros(at['gdf']['ACTIVITY_TREATMENT_UNIT_ID'].size)
		for i in range(at['gdf']['Year'].size):
			at['gdf']['Year'][i]=int(at['gdf']['ATU_COMPLETION_DATE'][i][0:4])
		gu.opickle(r'C:\Data\BC1ha\RSLT_ACTIVITY_TREATMENT_SVW\at.pkl',at)
	else:
		at=gu.ipickle(r'C:\Data\BC1ha\RSLT_ACTIVITY_TREATMENT_SVW\at.pkl')

	#ind=np.where(at['gdf']['Year']==2010)[0]
	#at['gdf']['OPENING_ID'][ind]

	flg=0
	if flg==1:
		A=np.zeros(tv.size)
		N=np.zeros(tv.size)
		for iT in range(tv.size):
			ind=np.where(at['gdf']['Year']==tv[iT])[0]
			A[iT]=np.nansum(at['gdf']['ACTUAL_TREATMENT_AREA'][ind])
			N[iT]=np.nansum(at['gdf']['ACTUAL_PLANTED_NUMBER'][ind])
		plt.close('all'); plt.plot(tv,A,'-bo')
		#plt.close('all'); plt.plot(tv,N/A,'-bo')
		
	# Import opening ID with spatial   
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
	
	# Index to planting and year
	dP={}
	for iT in range(tv.size):
		dP[tv[iT]]={}
		for iS in range(4):
			dP[tv[iT]][iS]={}
			dP[tv[iT]][iS]={'IndexToGrid':np.array([],dtype=int),
					   'Source FC':np.array([]),
					   'SPH_Planted':np.array([]),
					   'ID_SILV_FUND_SOURCE_CODE':np.array([]),
					   'ACTIVITY_TREATMENT_UNIT_ID':np.array([]),
					   'ID_SILV_TECHNIQUE_CODE':np.array([])}
	
	N_MissingArea=0 # 2 entries with missing area
	for iAT in range(at['gdf']['Year'].size):
		print(iAT)
		Year=at['gdf']['Year'][iAT].astype(int)
		if (Year<tv[0]) | (Year>tv[-1]):
			continue
		ID=at['gdf']['OPENING_ID'][iAT]		
		FSC=at['gdf']['ID_SILV_FUND_SOURCE_CODE'][iAT]
		ATUID=at['gdf']['ACTIVITY_TREATMENT_UNIT_ID'][iAT]
		STC=at['gdf']['ID_SILV_TECHNIQUE_CODE'][iAT]
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
					dP[Year][iS]['Source FC']=np.append(dP[Year][iS]['Source FC'],1*np.ones(ind.size))
					dP[Year][iS]['SPH_Planted']=np.append(dP[Year][iS]['SPH_Planted'],SPH_Planted*np.ones(ind.size))
					dP[Year][iS]['ID_SILV_FUND_SOURCE_CODE']=np.append(dP[Year][iS]['ID_SILV_FUND_SOURCE_CODE'],FSC*np.ones(ind.size))
					dP[Year][iS]['ACTIVITY_TREATMENT_UNIT_ID']=np.append(dP[Year][iS]['ACTIVITY_TREATMENT_UNIT_ID'],ATUID*np.ones(ind.size))
					dP[Year][iS]['ID_SILV_TECHNIQUE_CODE']=np.append(dP[Year][iS]['ID_SILV_TECHNIQUE_CODE'],STC*np.ones(ind.size))
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
			dP[Year][iS]['Source FC']=np.append(dP[Year][iS]['Source FC'],2*np.ones(ind2.size))
			dP[Year][iS]['SPH_Planted']=np.append(dP[Year][iS]['SPH_Planted'],SPH_Planted*np.ones(ind2.size))
			dP[Year][iS]['ID_SILV_FUND_SOURCE_CODE']=np.append(dP[Year][iS]['ID_SILV_FUND_SOURCE_CODE'],FSC*np.ones(ind2.size))
			dP[Year][iS]['ACTIVITY_TREATMENT_UNIT_ID']=np.append(dP[Year][iS]['ACTIVITY_TREATMENT_UNIT_ID'],ATUID*np.ones(ind2.size))
			dP[Year][iS]['ID_SILV_TECHNIQUE_CODE']=np.append(dP[Year][iS]['ID_SILV_TECHNIQUE_CODE'],STC*np.ones(ind2.size))
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
			dP[Year][iS]['Source FC']=np.append(dP[Year][iS]['Source FC'],3*np.ones(ind2.size))
			dP[Year][iS]['SPH_Planted']=np.append(dP[Year][iS]['SPH_Planted'],SPH_Planted*np.ones(ind2.size))
			dP[Year][iS]['ID_SILV_FUND_SOURCE_CODE']=np.append(dP[Year][iS]['ID_SILV_FUND_SOURCE_CODE'],FSC*np.ones(ind2.size))
			dP[Year][iS]['ACTIVITY_TREATMENT_UNIT_ID']=np.append(dP[Year][iS]['ACTIVITY_TREATMENT_UNIT_ID'],ATUID*np.ones(ind2.size))
			dP[Year][iS]['ID_SILV_TECHNIQUE_CODE']=np.append(dP[Year][iS]['ID_SILV_TECHNIQUE_CODE'],STC*np.ones(ind2.size))
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
			dP[Year][iS]['Source FC']=np.append(dP[Year][iS]['Source FC'],4*np.ones(ind2.size))
			dP[Year][iS]['SPH_Planted']=np.append(dP[Year][iS]['SPH_Planted'],SPH_Planted*np.ones(ind2.size))
			dP[Year][iS]['ID_SILV_FUND_SOURCE_CODE']=np.append(dP[Year][iS]['ID_SILV_FUND_SOURCE_CODE'],FSC*np.ones(ind2.size))
			dP[Year][iS]['ACTIVITY_TREATMENT_UNIT_ID']=np.append(dP[Year][iS]['ACTIVITY_TREATMENT_UNIT_ID'],ATUID*np.ones(ind2.size))
			dP[Year][iS]['ID_SILV_TECHNIQUE_CODE']=np.append(dP[Year][iS]['ID_SILV_TECHNIQUE_CODE'],STC*np.ones(ind2.size))
			print('From VRI')
	 
		#print('Missing')
	gu.opickle(r'C:\Data\BC1ha\RSLT_ACTIVITY_TREATMENT_SVW\dP.pkl',dP)
	#dP=gu.ipickle(r'C:\Data\BC1ha\RSLT_ACTIVITY_TREATMENT_SVW\dP.pkl')

	# Pack

	# Initialize rasters
	N_Year=6
	vL=['ACTIVITY_TREATMENT_UNIT_ID','ID_SILV_FUND_SOURCE_CODE','SPH_Planted','ID_SILV_TECHNIQUE_CODE']
	zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])	
	zPac={'Year':{},'ACTIVITY_TREATMENT_UNIT_ID':{},'ID_SILV_FUND_SOURCE_CODE':{},'SPH_Planted':{},'ID_SILV_TECHNIQUE_CODE':{}}
	for iY in range(N_Year):
		for k in zPac.keys():
			if k=='ACTIVITY_TREATMENT_UNIT_ID':
				zPac[k][iY+1]=np.zeros(zRef['Data'].shape,dtype='int32')
			else:
				zPac[k][iY+1]=np.zeros(zRef['Data'].shape,dtype='int16')

	vNam='PL_All'
	for iT in range(tv.size):
		print(tv[iT])
		
		zYr={}
		for k in zPac.keys():
			zYr[k]=np.zeros(zRef['Data'].shape,dtype=float)

		# Add activities without spatial
		iS=0
		iA=indFC[0][dP[tv[iT]][iS]['IndexToGrid']]
		iB=indFC[1][dP[tv[iT]][iS]['IndexToGrid']]
		zYr['ACTIVITY_TREATMENT_UNIT_ID'][ iA,iB  ]=dP[tv[iT]][iS]['ACTIVITY_TREATMENT_UNIT_ID']
		zYr['ID_SILV_FUND_SOURCE_CODE'][ iA,iB  ]=dP[tv[iT]][iS]['ID_SILV_FUND_SOURCE_CODE']
		zYr['SPH_Planted'][ iA,iB  ]=dP[tv[iT]][iS]['SPH_Planted']
		zYr['ID_SILV_TECHNIQUE_CODE'][ iA,iB  ]=dP[tv[iT]][iS]['ID_SILV_TECHNIQUE_CODE']
		
		iS=1
		iA=indOP1[0][dP[tv[iT]][iS]['IndexToGrid']]
		iB=indOP1[1][dP[tv[iT]][iS]['IndexToGrid']]
		zYr['ACTIVITY_TREATMENT_UNIT_ID'][ iA,iB  ]=dP[tv[iT]][iS]['ACTIVITY_TREATMENT_UNIT_ID']
		zYr['ID_SILV_FUND_SOURCE_CODE'][ iA,iB  ]=dP[tv[iT]][iS]['ID_SILV_FUND_SOURCE_CODE']
		zYr['SPH_Planted'][ iA,iB  ]=dP[tv[iT]][iS]['SPH_Planted']
		zYr['ID_SILV_TECHNIQUE_CODE'][ iA,iB  ]=dP[tv[iT]][iS]['ID_SILV_TECHNIQUE_CODE']
		
		iS=2
		iA=indOP2[0][dP[tv[iT]][iS]['IndexToGrid']]
		iB=indOP2[1][dP[tv[iT]][iS]['IndexToGrid']]
		zYr['ACTIVITY_TREATMENT_UNIT_ID'][ iA,iB  ]=dP[tv[iT]][iS]['ACTIVITY_TREATMENT_UNIT_ID']
		zYr['ID_SILV_FUND_SOURCE_CODE'][ iA,iB  ]=dP[tv[iT]][iS]['ID_SILV_FUND_SOURCE_CODE']
		zYr['SPH_Planted'][ iA,iB  ]=dP[tv[iT]][iS]['SPH_Planted']
		zYr['ID_SILV_TECHNIQUE_CODE'][ iA,iB  ]=dP[tv[iT]][iS]['ID_SILV_TECHNIQUE_CODE']

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
		z1['Data']=zPac['ACTIVITY_TREATMENT_UNIT_ID'][iY+1].astype('int32')
		gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + vNam + '_' + str(iY+1) + '_ACTIVITY_TREATMENT_UNIT_ID.tif')
		z1=copy.deepcopy(zRef)
		z1['Data']=zPac['ID_SILV_FUND_SOURCE_CODE'][iY+1].astype('int16')
		gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + vNam + '_' + str(iY+1) + '_SILV_FUND_SOURCE_CODE.tif')
		z1=copy.deepcopy(zRef)
		z1['Data']=zPac['SPH_Planted'][iY+1].astype('int16')
		gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + vNam + '_' + str(iY+1) + '_SPH_Planted.tif')
		z1=copy.deepcopy(zRef)
		z1['Data']=zPac['ID_SILV_TECHNIQUE_CODE'][iY+1].astype('int16')
		gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + vNam + '_' + str(iY+1) + '_SILV_TECHNIQUE_CODE.tif')
	
	return

#%% Rasterize direct seeding

def RasterizeDirectSeeding(meta):
	
	tv=np.arange(1960,2024,1)
	zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
	
	# Start with planting with spatial from RESULTS (takes 15 min)
	flg=0
	if flg==1:
		t0=time.time()
		ats={}
		ats['Path']=meta['Paths']['GDB']['Results']
		ats['Layer']='RSLT_ACTIVITY_TREATMENT_SVW'; # fiona.listlayers(ats['Path'])
		ats['crs']=meta['Geos']['crs']
		ats['Keep Geom']='On'
		ats['Select Openings']=np.array([])
		ats['SBC']=np.array(['DS'])
		ats['STC']=np.array([])
		ats['SMC']=np.array([])
		ats['FSC']=np.array([])
		ats['SOC1']=np.array([])
		ats['ROI']=[]
		ats['gdf']=qgdb.Query_Openings(ats,[])
		ats['gdf']=ats['gdf'][ats['gdf'].geometry!=None]
		ats['gdf']=ats['gdf'].reset_index()	
		ats['gdf']['Year']=np.zeros(len(ats['gdf']))
		for i in range(ats['gdf']['Year'].size):
			ats['gdf']['Year'][i]=int(ats['gdf']['ATU_COMPLETION_DATE'][i][0:4])
		
		AreaPlanted=ats['gdf']['ACTUAL_TREATMENT_AREA']
		NumTreesPlanted=ats['gdf']['ACTUAL_PLANTED_NUMBER']
		ats['gdf']['SPH_Planted']=NumTreesPlanted/AreaPlanted
		ats['gdf']=CreateIdForCategoricalVariable(meta,'RSLT_ACTIVITY_TREATMENT_SVW','SILV_FUND_SOURCE_CODE',ats['gdf'])
		ats['gdf']=CreateIdForCategoricalVariable(meta,'RSLT_ACTIVITY_TREATMENT_SVW','SILV_TECHNIQUE_CODE',ats['gdf'])	
		ats['gdf']=ats['gdf'][ (ats['gdf']['RESULTS_IND']=='Y') & (ats['gdf']['SILV_METHOD_CODE']!='LAYOT') ].copy()
		ats['gdf']=ats['gdf'][ats['gdf'].geometry!=None]
		ats['gdf']=ats['gdf'].reset_index()
		#ats['gdf'].to_file(r'C:\Data\BC1ha\RSLT_ACTIVITY_TREATMENT_SVW\ats_ds.geojson',driver="GeoJSON")
	else:
		ats={}
		ats['gdf']=gpd.read_file(r'C:\Data\BC1ha\RSLT_ACTIVITY_TREATMENT_SVW\ats_ds.geojson')

	# Add areas where FC is artificial
	flg=0
	if flg==1:
		at={}
		at['Path']=meta['Paths']['GDB']['Results']
		at['Layer']='RSLT_ACTIVITY_TREATMENT_SVW'; # fiona.listlayers(at['Path'])
		at['crs']=meta['Geos']['crs']
		at['Keep Geom']='Off'
		at['Select Openings']=np.array([])
		at['SBC']=np.array(['DS'])
		at['STC']=np.array([])
		at['SMC']=np.array([])
		at['FSC']=np.array([])
		at['SOC1']=np.array([])
		at['ROI']=[]
		at['gdf']=qgdb.Query_Openings(at,[])
		at['gdf']=CreateIdForCategoricalVariable(meta,'RSLT_ACTIVITY_TREATMENT_SVW','SILV_FUND_SOURCE_CODE',at['gdf'])
		at['gdf']=CreateIdForCategoricalVariable(meta,'RSLT_ACTIVITY_TREATMENT_SVW','SILV_TECHNIQUE_CODE',at['gdf'])
		# Make sure to remove entries that we know did not occur (planned or layout)
		ikp=np.where(  (at['gdf']['RESULTS_IND']=='Y') & (at['gdf']['SILV_METHOD_CODE']!='LAYOT') )[0]
		for k in at['gdf'].keys():
			at['gdf'][k]=at['gdf'][k][ikp]
		at['gdf']['Year']=np.zeros(at['gdf']['ACTIVITY_TREATMENT_UNIT_ID'].size)
		for i in range(at['gdf']['Year'].size):
			at['gdf']['Year'][i]=int(at['gdf']['ATU_COMPLETION_DATE'][i][0:4])
		gu.opickle(r'C:\Data\BC1ha\RSLT_ACTIVITY_TREATMENT_SVW\at_ds.pkl',at)
	else:
		at=gu.ipickle(r'C:\Data\BC1ha\RSLT_ACTIVITY_TREATMENT_SVW\at_ds.pkl')

	# Import opening ID with spatial   
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
	
	# Index to planting and year
	dDS={}
	for iT in range(tv.size):
		dDS[tv[iT]]={}
		for iS in range(4):
			dDS[tv[iT]][iS]={}
			dDS[tv[iT]][iS]={'IndexToGrid':np.array([],dtype=int),
					   'Source FC':np.array([]),
					   'SPH_Planted':np.array([]),
					   'ID_SILV_FUND_SOURCE_CODE':np.array([]),
					   'ACTIVITY_TREATMENT_UNIT_ID':np.array([]),
					   'ID_SILV_TECHNIQUE_CODE':np.array([])}
	
	N_MissingArea=0 # 2 entries with missing area
	for iAT in range(at['gdf']['Year'].size):
		print(iAT)
		Year=at['gdf']['Year'][iAT].astype(int)
		if (Year<tv[0]) | (Year>tv[-1]):
			continue
		ID=at['gdf']['OPENING_ID'][iAT]		
		FSC=at['gdf']['ID_SILV_FUND_SOURCE_CODE'][iAT]
		ATUID=at['gdf']['ACTIVITY_TREATMENT_UNIT_ID'][iAT]
		STC=at['gdf']['ID_SILV_TECHNIQUE_CODE'][iAT]
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
					dDS[Year][iS]['IndexToGrid']=np.append(dDS[Year][iS]['IndexToGrid'],ind)
					dDS[Year][iS]['Source FC']=np.append(dDS[Year][iS]['Source FC'],1*np.ones(ind.size))
					dDS[Year][iS]['SPH_Planted']=np.append(dDS[Year][iS]['SPH_Planted'],SPH_Planted*np.ones(ind.size))
					dDS[Year][iS]['ID_SILV_FUND_SOURCE_CODE']=np.append(dDS[Year][iS]['ID_SILV_FUND_SOURCE_CODE'],FSC*np.ones(ind.size))
					dDS[Year][iS]['ACTIVITY_TREATMENT_UNIT_ID']=np.append(dDS[Year][iS]['ACTIVITY_TREATMENT_UNIT_ID'],ATUID*np.ones(ind.size))
					dDS[Year][iS]['ID_SILV_TECHNIQUE_CODE']=np.append(dDS[Year][iS]['ID_SILV_TECHNIQUE_CODE'],STC*np.ones(ind.size))
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
			dDS[Year][iS]['IndexToGrid']=np.append(dDS[Year][iS]['IndexToGrid'],ind2)
			dDS[Year][iS]['Source FC']=np.append(dDS[Year][iS]['Source FC'],2*np.ones(ind2.size))
			dDS[Year][iS]['SPH_Planted']=np.append(dDS[Year][iS]['SPH_Planted'],SPH_Planted*np.ones(ind2.size))
			dDS[Year][iS]['ID_SILV_FUND_SOURCE_CODE']=np.append(dDS[Year][iS]['ID_SILV_FUND_SOURCE_CODE'],FSC*np.ones(ind2.size))
			dDS[Year][iS]['ACTIVITY_TREATMENT_UNIT_ID']=np.append(dDS[Year][iS]['ACTIVITY_TREATMENT_UNIT_ID'],ATUID*np.ones(ind2.size))
			dDS[Year][iS]['ID_SILV_TECHNIQUE_CODE']=np.append(dDS[Year][iS]['ID_SILV_TECHNIQUE_CODE'],STC*np.ones(ind2.size))
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
			dDS[Year][iS]['IndexToGrid']=np.append(dDS[Year][iS]['IndexToGrid'],ind2)
			dDS[Year][iS]['Source FC']=np.append(dDS[Year][iS]['Source FC'],3*np.ones(ind2.size))
			dDS[Year][iS]['SPH_Planted']=np.append(dDS[Year][iS]['SPH_Planted'],SPH_Planted*np.ones(ind2.size))
			dDS[Year][iS]['ID_SILV_FUND_SOURCE_CODE']=np.append(dDS[Year][iS]['ID_SILV_FUND_SOURCE_CODE'],FSC*np.ones(ind2.size))
			dDS[Year][iS]['ACTIVITY_TREATMENT_UNIT_ID']=np.append(dDS[Year][iS]['ACTIVITY_TREATMENT_UNIT_ID'],ATUID*np.ones(ind2.size))
			dDS[Year][iS]['ID_SILV_TECHNIQUE_CODE']=np.append(dDS[Year][iS]['ID_SILV_TECHNIQUE_CODE'],STC*np.ones(ind2.size))
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
			dDS[Year][iS]['IndexToGrid']=np.append(dDS[Year][iS]['IndexToGrid'],ind2)
			dDS[Year][iS]['Source FC']=np.append(dDS[Year][iS]['Source FC'],4*np.ones(ind2.size))
			dDS[Year][iS]['SPH_Planted']=np.append(dDS[Year][iS]['SPH_Planted'],SPH_Planted*np.ones(ind2.size))
			dDS[Year][iS]['ID_SILV_FUND_SOURCE_CODE']=np.append(dDS[Year][iS]['ID_SILV_FUND_SOURCE_CODE'],FSC*np.ones(ind2.size))
			dDS[Year][iS]['ACTIVITY_TREATMENT_UNIT_ID']=np.append(dDS[Year][iS]['ACTIVITY_TREATMENT_UNIT_ID'],ATUID*np.ones(ind2.size))
			dDS[Year][iS]['ID_SILV_TECHNIQUE_CODE']=np.append(dDS[Year][iS]['ID_SILV_TECHNIQUE_CODE'],STC*np.ones(ind2.size))
			print('From VRI')
	 
		#print('Missing')
	gu.opickle(r'C:\Data\BC1ha\RSLT_ACTIVITY_TREATMENT_SVW\dDS.pkl',dDS)
	#dDS=gu.ipickle(r'C:\Data\BC1ha\RSLT_ACTIVITY_TREATMENT_SVW\dDS.pkl')

	# Pack

	# Initialize rasters
	N_Year=6
	vL=['ACTIVITY_TREATMENT_UNIT_ID','ID_SILV_FUND_SOURCE_CODE','SPH_Planted','ID_SILV_TECHNIQUE_CODE']
	zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])	
	zPac={'Year':{},'ACTIVITY_TREATMENT_UNIT_ID':{},'ID_SILV_FUND_SOURCE_CODE':{},'SPH_Planted':{},'ID_SILV_TECHNIQUE_CODE':{}}
	for iY in range(N_Year):
		for k in zPac.keys():
			if k=='ACTIVITY_TREATMENT_UNIT_ID':
				zPac[k][iY+1]=np.zeros(zRef['Data'].shape,dtype='int32')
			else:
				zPac[k][iY+1]=np.zeros(zRef['Data'].shape,dtype='int16')

	vNam='DS_All'
	for iT in range(tv.size):
		print(tv[iT])
		
		zYr={}
		for k in zPac.keys():
			zYr[k]=np.zeros(zRef['Data'].shape,dtype=float)

		# Add activities without spatial
		iS=0
		iA=indFC[0][dDS[tv[iT]][iS]['IndexToGrid']]
		iB=indFC[1][dDS[tv[iT]][iS]['IndexToGrid']]
		zYr['ACTIVITY_TREATMENT_UNIT_ID'][ iA,iB  ]=dDS[tv[iT]][iS]['ACTIVITY_TREATMENT_UNIT_ID']
		zYr['ID_SILV_FUND_SOURCE_CODE'][ iA,iB  ]=dDS[tv[iT]][iS]['ID_SILV_FUND_SOURCE_CODE']
		zYr['SPH_Planted'][ iA,iB  ]=dDS[tv[iT]][iS]['SPH_Planted']
		zYr['ID_SILV_TECHNIQUE_CODE'][ iA,iB  ]=dDS[tv[iT]][iS]['ID_SILV_TECHNIQUE_CODE']
		
		iS=1
		iA=indOP1[0][dDS[tv[iT]][iS]['IndexToGrid']]
		iB=indOP1[1][dDS[tv[iT]][iS]['IndexToGrid']]
		zYr['ACTIVITY_TREATMENT_UNIT_ID'][ iA,iB  ]=dDS[tv[iT]][iS]['ACTIVITY_TREATMENT_UNIT_ID']
		zYr['ID_SILV_FUND_SOURCE_CODE'][ iA,iB  ]=dDS[tv[iT]][iS]['ID_SILV_FUND_SOURCE_CODE']
		zYr['SPH_Planted'][ iA,iB  ]=dDS[tv[iT]][iS]['SPH_Planted']
		zYr['ID_SILV_TECHNIQUE_CODE'][ iA,iB  ]=dDS[tv[iT]][iS]['ID_SILV_TECHNIQUE_CODE']
		
		iS=2
		iA=indOP2[0][dDS[tv[iT]][iS]['IndexToGrid']]
		iB=indOP2[1][dDS[tv[iT]][iS]['IndexToGrid']]
		zYr['ACTIVITY_TREATMENT_UNIT_ID'][ iA,iB  ]=dDS[tv[iT]][iS]['ACTIVITY_TREATMENT_UNIT_ID']
		zYr['ID_SILV_FUND_SOURCE_CODE'][ iA,iB  ]=dDS[tv[iT]][iS]['ID_SILV_FUND_SOURCE_CODE']
		zYr['SPH_Planted'][ iA,iB  ]=dDS[tv[iT]][iS]['SPH_Planted']
		zYr['ID_SILV_TECHNIQUE_CODE'][ iA,iB  ]=dDS[tv[iT]][iS]['ID_SILV_TECHNIQUE_CODE']

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
		z1['Data']=zPac['ACTIVITY_TREATMENT_UNIT_ID'][iY+1].astype('int32')
		gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + vNam + '_' + str(iY+1) + '_ACTIVITY_TREATMENT_UNIT_ID.tif')
		z1=copy.deepcopy(zRef)
		z1['Data']=zPac['ID_SILV_FUND_SOURCE_CODE'][iY+1].astype('int16')
		gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + vNam + '_' + str(iY+1) + '_SILV_FUND_SOURCE_CODE.tif')
		z1=copy.deepcopy(zRef)
		z1['Data']=zPac['SPH_Planted'][iY+1].astype('int16')
		gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + vNam + '_' + str(iY+1) + '_SPH_Planted.tif')
		z1=copy.deepcopy(zRef)
		z1['Data']=zPac['ID_SILV_TECHNIQUE_CODE'][iY+1].astype('int16')
		gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + vNam + '_' + str(iY+1) + '_SILV_TECHNIQUE_CODE.tif')
	
	return

#%% Planting (Non-obligation by project type)
# Dabbled in doing it this way and then went back to doing it on the fly

def MaskPlantingNonOb(meta,zRef):
	zMask=copy.deepcopy(zRef)
	zMask['Data']=np.zeros(zRef['Data'].shape,dtype='int16')
	for iEY in range(6):
		zFSC=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(iEY+1) + '_SILV_FUND_SOURCE_CODE.tif')
		ind=np.where( (np.isin(zFSC['Data'],meta['Param']['BE']['FSC']['NO List ID'])==True) )
		zMask['Data'][ind]=1
	gis.SaveGeoTiff(zMask,meta['Paths']['bc1ha'] + '\\Management\\PL_NonOb_MaskAll.tif')
	return

#%%
def RasterizeSilviculture(meta,sbc,stc,smc,soc1,vNam):
	
	#sbc=meta,np.array(['SP']),np.array(['BU']),np.array(['BROAD']),np.array([]),np.array(['SP-BU-BROAD'])
	sbc=np.array(['SP'])
	stc=np.array(['ME'])
	smc=np.array(['DISC','EXCAV','MOUND','LRIP','RRIP','WING'])
	soc1=np.array([])
	vNam='SP-ME'
	
	tv=np.arange(1960,2024,1)
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
		AreaPlanted=ats['gdf']['ACTUAL_TREATMENT_AREA']		
		ats['gdf']=u1ha.CreateIdForCategoricalVariable(meta,'RSLT_ACTIVITY_TREATMENT_SVW','SILV_FUND_SOURCE_CODE',ats['gdf'])
		ats['gdf']=u1ha.CreateIdForCategoricalVariable(meta,'RSLT_ACTIVITY_TREATMENT_SVW','SILV_BASE_CODE',ats['gdf'])
		ats['gdf']=u1ha.CreateIdForCategoricalVariable(meta,'RSLT_ACTIVITY_TREATMENT_SVW','SILV_TECHNIQUE_CODE',ats['gdf'])
		ats['gdf']=u1ha.CreateIdForCategoricalVariable(meta,'RSLT_ACTIVITY_TREATMENT_SVW','SILV_METHOD_CODE',ats['gdf'])
		ats['gdf']=ats['gdf'][ (ats['gdf']['RESULTS_IND']=='Y') & (ats['gdf']['SILV_METHOD_CODE']!='LAYOT') ].copy()
		ats['gdf']=ats['gdf'][ats['gdf'].geometry!=None]
		ats['gdf']=ats['gdf'].reset_index()
		#ats['gdf'].to_file(r'C:\Data\BC1ha\RSLT_ACTIVITY_TREATMENT_SVW\ats' + vNam + '.geojson',driver="GeoJSON")
	#else:
	#	ats={}
	#	ats['gdf']=gpd.read_file(r'C:\Data\BC1ha\RSLT_ACTIVITY_TREATMENT_SVW\ats' + vNam + '.geojson')

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
		#gu.opickle(r'C:\Data\BC1ha\RSLT_ACTIVITY_TREATMENT_SVW\at_' + vNam + '.pkl',at)
	#else:
	#	at=gu.ipickle(r'C:\Data\BC1ha\RSLT_ACTIVITY_TREATMENT_SVW\at' + vNam + '.pkl')
	
	np.sum(ats['gdf']['ACTUAL_TREATMENT_AREA'])
	
	ind=np.where(at['gdf']['SILV_METHOD_CODE']=='BROAD')[0]
	np.sum(at['gdf']['ACTUAL_TREATMENT_AREA'])
	ind.size
	
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
	
	# Index to planting and year
	dP={}
	for iT in range(tv.size):
		dP[tv[iT]]={}
		for iS in range(4):
			dP[tv[iT]][iS]={}
			dP[tv[iT]][iS]={'IndexToGrid':np.array([],dtype=int),
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
	vL=['ACTIVITY_TREATMENT_UNIT_ID','ID_SILV_FUND_SOURCE_CODE','ID_SILV_BASE_CODE','ID_SILV_TECHNIQUE_CODE','ID_SILV_METHOD_CODE']
	zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])	
	zPac={'Year':{},'ACTIVITY_TREATMENT_UNIT_ID':{},'ID_SILV_FUND_SOURCE_CODE':{},'ID_SILV_BASE_CODE':{},'ID_SILV_TECHNIQUE_CODE':{},'ID_SILV_METHOD_CODE':{}}
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
		zYr['ACTIVITY_TREATMENT_UNIT_ID'][ iA,iB  ]=dP[tv[iT]][iS]['ACTIVITY_TREATMENT_UNIT_ID']
		zYr['ID_SILV_FUND_SOURCE_CODE'][ iA,iB  ]=dP[tv[iT]][iS]['ID_SILV_FUND_SOURCE_CODE']
		zYr['ID_SILV_BASE_CODE'][ iA,iB  ]=dP[tv[iT]][iS]['ID_SILV_BASE_CODE']
		zYr['ID_SILV_TECHNIQUE_CODE'][ iA,iB  ]=dP[tv[iT]][iS]['ID_SILV_TECHNIQUE_CODE']
		zYr['ID_SILV_METHOD_CODE'][ iA,iB  ]=dP[tv[iT]][iS]['ID_SILV_METHOD_CODE']
		
		iS=1
		iA=indOP1[0][dP[tv[iT]][iS]['IndexToGrid']]
		iB=indOP1[1][dP[tv[iT]][iS]['IndexToGrid']]
		zYr['ACTIVITY_TREATMENT_UNIT_ID'][ iA,iB  ]=dP[tv[iT]][iS]['ACTIVITY_TREATMENT_UNIT_ID']
		zYr['ID_SILV_FUND_SOURCE_CODE'][ iA,iB  ]=dP[tv[iT]][iS]['ID_SILV_FUND_SOURCE_CODE']
		zYr['ID_SILV_BASE_CODE'][ iA,iB  ]=dP[tv[iT]][iS]['ID_SILV_BASE_CODE']
		zYr['ID_SILV_TECHNIQUE_CODE'][ iA,iB  ]=dP[tv[iT]][iS]['ID_SILV_TECHNIQUE_CODE']
		zYr['ID_SILV_METHOD_CODE'][ iA,iB  ]=dP[tv[iT]][iS]['ID_SILV_METHOD_CODE']
		
		iS=2
		iA=indOP2[0][dP[tv[iT]][iS]['IndexToGrid']]
		iB=indOP2[1][dP[tv[iT]][iS]['IndexToGrid']]
		zYr['ACTIVITY_TREATMENT_UNIT_ID'][ iA,iB  ]=dP[tv[iT]][iS]['ACTIVITY_TREATMENT_UNIT_ID']
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
		z1['Data']=zPac['ACTIVITY_TREATMENT_UNIT_ID'][iY+1].astype('int32')
		gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + vNam + '_' + str(iY+1) + '_ACTIVITY_TREATMENT_UNIT_ID.tif')
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
	z1['Data']=np.zeros(zRef['Data'].shape,dtype='int16')
	for iY in range(N_Year):
		z=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + vNam + '_' + str(iY+1) + '_Year.tif')
		ind=np.where(z['Data']>0)
		z1['Data'][ind]=z['Data'][ind]
	gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + vNam + '_MaskAll.tif')

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

#%% 
def ClimateStatsByBGCZone(meta):
	# *** Needs updating ***
	zBGC=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\becz.tif')
	zBGC['Data']=zBGC['Data'].flatten()
	lutBGC=gu.ReadExcel(meta['Paths']['bc1ha'] + '\\VRI 2023\\becz_lut.xlsx')
	
	zMAT=gis.OpenGeoTiff(r'C:\Data\BC1ha\Climate\BC1ha_mat_norm_1971to2000_si_hist_v1.tif')
	zMAT['Data']=zMAT['Data'].flatten().astype(float)/10
	
	zWS=gis.OpenGeoTiff(r'C:\Data\BC1ha\Climate\BC1ha_ws_gs_norm_1971to2000_comp_hist_v1.tif')
	zWS['Data']=zWS['Data'].flatten().astype(float)
	
	zSI=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\si.tif')
	zSI['Data']=zSI['Data'].flatten().astype(float)
	
	zA=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\proj_age_1.tif')
	zA['Data']=zA['Data'].flatten().astype(float)
	
	lutBGC['MAT']=np.zeros(lutBGC['VALUE'].size)
	lutBGC['WS']=np.zeros(lutBGC['VALUE'].size)
	lutBGC['SI']=np.zeros(lutBGC['VALUE'].size)
	lutBGC['Age']=np.zeros(lutBGC['VALUE'].size)
	for i in range(lutBGC['VALUE'].size):
		ind=np.where( (zBGC['Data']==lutBGC['VALUE'][i]) & (zMAT['Data']>=-50) & (zWS['Data']>=0) & (zWS['Data']<=200) & (zSI['Data']>0) & (zSI['Data']<100) & (zA['Data']>=0) & (zA['Data']<1000) )[0]
		lutBGC['MAT'][i]=np.mean(zMAT['Data'][ind])
		lutBGC['WS'][i]=np.mean(zWS['Data'][ind])
		lutBGC['SI'][i]=np.mean(zSI['Data'][ind])
		lutBGC['Age'][i]=np.mean(zA['Data'][ind])
	
	df=pd.DataFrame(lutBGC)
	df.to_excel(r'C:\Data\BC1ha\Climate\tmp.xlsx')
	return

#%%
def RasterizeForestCoverInventory(meta):
	t0=time.time()

	zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
	zFID=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_FOREST_COVER_INV_SVW\\fcid.tif')['Data'].flatten()
	uFID=gu.IndicesFromUniqueArrayValues(zFID)
	
	fc={}
	fc['Path']=meta['Paths']['GDB']['Results']
	fc['Layer']='RSLT_FOREST_COVER_INV_SVW'; # fiona.listlayers(fc['Path'])
	fc['crs']=meta['Geos']['crs']
	fc['Keep Geom']='Off'
	fc['Select Openings']=np.array([])
	fc['SBC']=np.array([])
	fc['STC']=np.array([])
	fc['SMC']=np.array([])
	fc['FSC']=np.array([])
	fc['SOC1']=np.array([])
	fc['ROI']=[]
	fc['gdf']=qgdb.Query_Openings(fc,[])
	
	vL=['STOCKING_STATUS_CODE','STOCKING_TYPE_CODE','SILV_RESERVE_CODE','SILV_RESERVE_OBJECTIVE_CODE','TREE_COVER_PATTERN_CODE']
	for v in vL:
		fc['gdf']=u1ha.CreateIdForCategoricalVariable(meta,'RSLT_FOREST_COVER_INV_SVW',v,fc['gdf'])

	vL=['I_TOTAL_STEMS_PER_HA','I_TOTAL_WELL_SPACED_STEMS_HA','I_WELL_SPACED_STEMS_PER_HA','I_FREE_GROWING_STEMS_PER_HA','I_CROWN_CLOSURE_PERCENT']
	#vL=['REFERENCE_YEAR','OPENING_ID','ID_STOCKING_STATUS_CODE','ID_STOCKING_TYPE_CODE','ID_SILV_RESERVE_CODE','ID_SILV_RESERVE_OBJECTIVE_CODE','ID_TREE_COVER_PATTERN_CODE']
	z={}
	for v in vL:
		z[v]=np.zeros(zFID.size,dtype='int32')
		
	for i in range(fc['gdf']['OPENING_ID'].size):
		print(i)
		try:
			fid=fc['gdf']['FOREST_COVER_ID'][i]
			for v in vL:
				z[v][ uFID[fid] ]=fc['gdf'][v][i]
			#print('Working!')
		except:
			#print('Missing!')
			pass
	
	for v in vL:
		z1=copy.deepcopy(zRef)
		z1['Data']=np.reshape(z[v],zRef['Data'].shape)
		gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\RSLT_FOREST_COVER_INV_SVW\\' + v + '.tif')
		
	return   

#%% Clip rasters to standard grid
# *** This also compresses files that come out of Arc crazy big ***
def ClipToBC1ha(meta):

	# Forest Cover ID
	fin=r'C:\Data\BC1ha\RSLT_FOREST_COVER_INV_SVW\fcid.tif'
	fout=r'C:\Data\BC1ha\RSLT_FOREST_COVER_INV_SVW\fcid2.tif'
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

	# BTM
	fin=r'C:\Data\BC1ha\LandCoverUse\landuse.btm.tif'
	fout=r'C:\Data\BC1ha\LandCoverUse\landuse.btm.tif'
	gis.ClipToRaster_ByFile(fin,fout,meta['Paths']['bc1ha Ref Grid'])
	# Radiation
	fin=r'C:\Data\BC1ha\Climate\Seasonal\BC1ha_rswd_gs_norm_1971to2000_si_hist_v1_c.tif'
	fout=r'C:\Data\BC1ha\Climate\Seasonal\BC1ha_rswd_gs_norm_1971to2000_si_hist_v1_c.tif'
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
	d=gu.ReadExcel(path)
	vList=[]
	for i in range(d['Name'].size):
		if d['Included'][i]==1:
			vList.append(d['Name'][i])
	return vList

#%% Digitize the boundary of TSAs (the original is organized at the sub-TSA level)
def DigitizeTSABoundaries(meta):
	zTSA=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\FADM_TSA\\TSA_NUMBER_DESCRIPTION.tif')
	u=np.unique(zTSA['Data'])
	u=u[u>0]
	gdf=gpd.GeoDataFrame(data=[],columns=['Value','Name','geometry'])
	cnt=0
	for iU in range(u.size):
		ind=np.where(zTSA['Data']==u[iU])
		#z=np.zeros((zTSA['m'],zTSA['n'],3),dtype=np.uint8)
		z=np.zeros(zTSA['Data'].shape,dtype=np.uint8)
		z[ind]=255
		z=np.dstack([z]*3)
		z=cv2.cvtColor(z,cv2.COLOR_BGR2GRAY) # Convert to grey scale
		cont=cv2.findContours(image=z,mode=cv2.RETR_LIST,method=cv2.CHAIN_APPROX_SIMPLE)
		for j in range(len(cont[0])):
			cont_inner=cont[0][j].squeeze()
			if cont_inner.size==2:
				continue
			if cont_inner.shape[0]<3:
				continue
			pointList=[]
			for k in range(len(cont_inner)):
				c=cont_inner[k][0]
				r=cont_inner[k][1]
				x=int(zTSA['X'][0,c])
				y=int(zTSA['Y'][r,0])
				pointList.append(geometry.Point(x,y))
			gdf.loc[cnt,'Value']=int(u[iU])
			gdf.loc[cnt,'Name']=u1ha.lut_n2s(meta['LUT']['FADM_TSA']['TSA_NUMBER_DESCRIPTION'],u[iU])[0]
			gdf.loc[cnt,'geometry']=geometry.Polygon([[p.x,p.y] for p in pointList])
			cnt=cnt+1
	gdf.to_file(meta['Paths']['GDB']['GDB'] + '\\tsa.geojson',driver='GeoJSON')
	return

#%%
def DeriveHarvest_SILV_SYSTEM_CODE(meta):
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
	tv=np.arange(np.min(uYear),np.max(uYear),1)

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
def RasterizeBTKSpray(meta):
	gdf_spray=gpd.read_file(r'C:\Data\Aerial Btk Spray\Processed\btk_spray_comp.geojson')
	tv=np.arange(1950,2021,1)
	for iT in range(tv.size):
		print(tv[iT])
		zOut=copy.deepcopy(zRef)
		zOut['Data']=np.zeros(zRef['Data'].shape,dtype='int8')

		df0=gdf_spray[ (gdf_spray['Year']==tv[iT]) ].copy()
		df0=df0[df0.geometry!=None]
		df0['Dummy']=1

		if len(df0)>0:
			shapes=((geom,value) for geom, value in zip(df0.geometry,df0['Dummy']))
			z0=np.zeros(zRef['Data'].shape,dtype=float)
			burned=features.rasterize(shapes=shapes,fill=0,out=z0,transform=zRef['Transform'])
			zOut['Data']=burned.astype('int8')

		fout=meta['Paths']['bc1ha'] + '\\Management' + '\\btk_spray_' + str(tv[iT]) + '.tif'
		gis.SaveGeoTiff(zOut,fout)
	return

#%%
def RasterizeGFC_LossYear(meta):
	finL=['50N_120W','50N_130W','60N_120W','60N_130W','60N_140W']

	# Loss year
	z=copy.deepcopy(zRef)
	z['Data']=np.zeros(z['Data'].shape)
	for f in finL:
		fin=r'C:\Data\Global Forest Change\Hansen_GFC-2021-v1.9_lossyear_' + f + '.tif'
		fout=r'C:\Data\Global Forest Change\Hansen_GFC-2021-v1.9_lossyear_' + f + 'p.tif'
		gis.ReprojectRasterAndClipToRaster(fin,fout,meta['Paths']['bc1ha Ref Grid'],meta['Geos']['crs'])

		z0=gis.OpenGeoTiff(r'C:\Data\Global Forest Change\Hansen_GFC-2021-v1.9_lossyear_' + f + 'p.tif')
		ind=np.where(z0['Data']>0)
		z['Data'][ind]=z0['Data'][ind]

	ind=np.where(zRef['Data']==0)
	z['Data'][ind]=0

	z['Data']=z['Data'].astype('int16')
	ind=np.where(z['Data']>0)
	z['Data'][ind]=z['Data'][ind]+2000

	gis.SaveGeoTiff(z,meta['Paths']['bc1ha'] + '\\Disturbances\\GlobalForestChange_LossYear_2021.tif')
	return

#%%
def ImportBurnSeverityCanada(meta):
	# https://ftp.maps.canada.ca/pub/nrcan_rncan/Forest-fires_Incendie-de-foret/CanLaBS-Burned_Severity-Severite_des_feux/
	pthin=r'C:\Data\Wildfire\Burn Severity'
	fin=pthin + '\\bsn_yr1.tif'
	z=gis.OpenGeoTiff(fin); print(zRef['Data'].shape); print(z['Data'].shape)
	fout=pthin + '\\bsn_yr2.tif'
	gis.ReprojectRasterAndClipToRaster(fin,fout,meta['Paths']['bc1ha Ref Grid'],meta['Geos']['crs'])
	z=gis.OpenGeoTiff(fout); print(zRef['Data'].shape); print(z['Data'].shape) #z=gis.ClipToRaster(z,zRef)
	ind=np.where(zRef['Data']==0); z['Data'][ind]=0
	fout2=meta['Paths']['bc1ha'] + '\\Disturbances\\BurnSev_Can_Year.tif'
	gis.SaveGeoTiff(z,fout2)
	
	fin=pthin + '\\bsn_dnbr1.tif'
	z=gis.OpenGeoTiff(fin); print(zRef['Data'].shape); print(z['Data'].shape)
	fout=pthin + '\\bsn_dnbr2.tif'
	gis.ReprojectRasterAndClipToRaster(fin,fout,meta['Paths']['bc1ha Ref Grid'],meta['Geos']['crs'])
	z=gis.OpenGeoTiff(fout); print(zRef['Data'].shape); print(z['Data'].shape) #z=gis.ClipToRaster(z,zRef)
	ind=np.where(zRef['Data']==0); z['Data'][ind]=0
	fout2=meta['Paths']['bc1ha'] + '\\Disturbances\\BurnSev_Can_dNBR.tif'
	gis.SaveGeoTiff(z,fout2)
	return

#%%
def RasterizeBurnSeverityCurrentYear(meta,yr):
	zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
	# fiona.listlayers(pthin)
	# Current year
	yr=2023
	pthin=r'C:\Data\Wildfire\Wildfire.gdb'
	df=gpd.read_file(pthin,layer='VEG_BURN_SEVERITY_SAME_YR_SP')
	df['ID']=np.zeros(df['FIRE_YEAR'].size)
	for k in meta['LUT']['Derived']['burnsev_comp1'].keys():
		ind=np.where( (df['BURN_SEVERITY_RATING']==k) )
		df['ID'][ind]=meta['LUT']['Derived']['burnsev_comp1'][k]
	ind=np.where(df['FIRE_YEAR']==yr)[0]
	df=df.iloc[ind]
	df=df[df.geometry!=None]
	df=df.reset_index()
	shapes=((geom,value) for geom, value in zip(df['geometry'],df['ID']))
	z0=np.zeros(zRef['Data'].shape,dtype=float)
	burned=features.rasterize(shapes=shapes,fill=0,out=z0,transform=zRef['Transform'])
	z1=copy.deepcopy(zRef)
	z1['Data']=burned.astype('int16')
	plt.close(); plt.matshow(z1['Data'],clim=[0,5])
	gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\VEG_BURN_SEVERITY_SP\\BURN_SEVERITY_RATING_Current_' + str(yr) + '.tif')
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
	n,p=gu.CountByCategories(zBSH['Data'][zYearH['Data']>0],'Percent')
	plt.bar(n,p)

	# Import BC current year data
	zBSC=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VEG_BURN_SEVERITY_SP\\BURN_SEVERITY_RATING_Current_2023.tif')
	n,p=gu.CountByCategories(zBSC['Data'][zBSC['Data']>0],'Percent')
	plt.bar(n,p)

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
		zY=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\PROT_HISTORICAL_FIRE_POLYS_SP\\FIRE_YEAR_' + str(iY+1) + '_Year.tif')

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
		gis.SaveGeoTiff(zBS,meta['Paths']['bc1ha'] + '\\PROT_HISTORICAL_FIRE_POLYS_SP\\FIRE_YEAR_' + str(iY+1) + '_SevClass.tif')

	gis.SaveGeoTiff(zLastY,meta['Paths']['bc1ha'] + '\\Disturbances\\BurnSevComp_YearLast.tif')
	gis.SaveGeoTiff(zLastS,meta['Paths']['bc1ha'] + '\\Disturbances\\BurnSevComp_SevClassLast.tif')

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
		z=Import_Raster(meta,[],['refg','lc_comp1_2019','lc_vri_l5_02','harv_yr_con1'])
		z1=z['refg'].copy()
		z1['Data']=np.zeros(z['lc_vri_l5_02']['Data'].shape,dtype='int8')
		ind=np.where( (z['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (z['lc_vri_l5_02']['Data']==meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_5']['SP']) ); z1['Data'][ind]=meta['LUT']['Derived']['tdc']['Sparse']
		ind=np.where( (z['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (z['lc_vri_l5_02']['Data']==meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_5']['OP']) ); z1['Data'][ind]=meta['LUT']['Derived']['tdc']['Open']
		ind=np.where( (z['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (z['lc_vri_l5_02']['Data']==meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_5']['DE']) ); z1['Data'][ind]=meta['LUT']['Derived']['tdc']['Dense']
		#ind=np.where( (z['harv_yr_con1']['Data']==1) ); z1['Data'][ind]=meta['LUT']['Derived']['tdc']['Dense']
		ind=np.where( (z['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (z1['Data']==0) ); z1['Data'][ind]=meta['LUT']['Derived']['tdc']['Open']
		# plt.matshow(z1,clim=[0,3])
		gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\LandCoverUse\\TreeDensityClass_2002.tif')

	# 2015 VRI
	flg=1
	if flg==1:
		z=Import_Raster(meta,[],['refg','lc_comp1_2019','harv_yr_con1'])
		zL5=gis.OpenGeoTiff(r'C:\Data\BC1ha\VRI 2015\BCLCS_LEVEL_5.tif')
		z1=z['refg'].copy()
		z1['Data']=np.zeros(z['refg']['Data'].shape,dtype='int8')
		ind=np.where( (z['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (zL5['Data']==meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_5']['SP']) ); z1['Data'][ind]=meta['LUT']['Derived']['tdc']['Sparse']
		ind=np.where( (z['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (zL5['Data']==meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_5']['OP']) ); z1['Data'][ind]=meta['LUT']['Derived']['tdc']['Open']
		ind=np.where( (z['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (zL5['Data']==meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_5']['DE']) ); z1['Data'][ind]=meta['LUT']['Derived']['tdc']['Dense']
		#ind=np.where( (z['harv_yr_con1']['Data']==1) ); z1['Data'][ind]=meta['LUT']['Derived']['tdc']['Dense']
		ind=np.where( (z['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (z1['Data']==0) ); z1['Data'][ind]=meta['LUT']['Derived']['tdc']['Open']
		# plt.matshow(z1,clim=[0,3])
		gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\LandCoverUse\\TreeDensityClass_VRI2015.tif')

	# # 2023 VRI
	# z=Import_Raster(meta,[],['refg','lc_comp1_2019','lc_vri_l5','harv_yr_con1'])
		
	# z1=z['refg'].copy()
	# z1['Data']=np.zeros(z['lc_vri_l5']['Data'].shape,dtype='int8')
	# ind=np.where( (z['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (z['lc_vri_l5']['Data']==meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_5']['SP']) ); z1['Data'][ind]=meta['LUT']['Derived']['tdc']['Sparse']
	# ind=np.where( (z['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (z['lc_vri_l5']['Data']==meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_5']['OP']) ); z1['Data'][ind]=meta['LUT']['Derived']['tdc']['Open']
	# ind=np.where( (z['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (z['lc_vri_l5']['Data']==meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_5']['DE']) ); z1['Data'][ind]=meta['LUT']['Derived']['tdc']['Dense']
	# ind=np.where( (z['harv_yr_con1']['Data']==1) ); z1['Data'][ind]=meta['LUT']['Derived']['tdc']['Dense']
	# ind=np.where( (z['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (z1['Data']==0) ); z1['Data'][ind]=meta['LUT']['Derived']['tdc']['Open']
	# # plt.matshow(z1,clim=[0,3])
	# gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\LandCoverUse\\TreeDensityClass_Current.tif')
	
	# # With shrubs and grasses
	# z1=z['refg'].copy()
	# z1['Data']=np.zeros(z['lc_vri_l5']['Data'].shape,dtype='int8')
	# ind=np.where( (z['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (z['lc_vri_l5']['Data']==meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_5']['SP']) ); z1['Data'][ind]=1
	# ind=np.where( (z['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (z['lc_vri_l5']['Data']==meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_5']['OP']) ); z1['Data'][ind]=2
	# ind=np.where( (z['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (z['lc_vri_l5']['Data']==meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_5']['DE']) ); z1['Data'][ind]=3
	# ind=np.where( (z['harv_yr_con1']['Data']==1) ); z1['Data'][ind]=3
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

	vList=['bgcz','fire_yr','lc_comp1_2019','lu_comp1_2019','harv_yr_comp1','age_vri','d2road','d2fac']  # ,'d2fac','d2set'
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

#%% Planting layer (species and genetic worth)
def RasterizePlantingLayer(meta,zRef):
	pl={}
	pl['Path']=meta['Paths']['GDB']['Results']
	pl['Layer']='RSLT_PLANTING_SVW'; # fiona.listlayers(at['Path'])
	pl['crs']=meta['Geos']['crs']
	pl['Keep Geom']='Off'
	pl['Select Openings']=np.array([])
	pl['SBC']=np.array([])
	pl['STC']=np.array([])
	pl['SMC']=np.array([])
	pl['FSC']=np.array([])
	pl['SOC1']=np.array([])
	pl['ROI']=[]
	pl['gdf']=qgdb.Query_Openings(pl,[])
	pl['gdf']=CreateIdForCategoricalVariable(meta,'RSLT_PLANTING_SVW','SILV_TREE_SPECIES_CODE',pl['gdf'])
	
	dGW=gu.ReadExcel(meta['Paths']['Model']['Code'] + '\\Parameters\\Parameters_Seedlot_GW.xlsx')
	
	def round_retain_sum(x):
		N=np.round(np.sum(x)).astype(int)
		y=x.astype(int)
		M=np.sum(y)
		K=N-M
		z=y-x
		if K!=0:
			idx=np.argpartition(z,K)[:K]
			y[idx]+=1
		return y
	
	# Get information by ATU ID
	uAT=np.unique(pl['gdf']['ACTIVITY_TREATMENT_UNIT_ID'])
	
	plAT={}
	plAT['ACTIVITY_TREATMENT_UNIT_ID']=uAT
	vL=['PL_SPECIES_CD','PL_SPECIES_PCT','PL_SPECIES_GW']
	for v in vL:
		for i in range(6):
			plAT[v + str(i+1)]=np.zeros(uAT.size,dtype='int16')
	
	# Get weighted average genetic worth by AT ID
	N_no_spc=0
	N_no_num_planted=0
	N_sn_not_found=0
	for iU in range(uAT.size):
		iAT=np.where(pl['gdf']['ACTIVITY_TREATMENT_UNIT_ID']==uAT[iU])[0]
	
		cd0=pl['gdf']['ID_SILV_TREE_SPECIES_CODE'][iAT]
		pct0=np.nan_to_num(pl['gdf']['NUMBER_PLANTED'][iAT]/np.sum(pl['gdf']['NUMBER_PLANTED'][iAT])*100)
		sn0=np.nan_to_num(pl['gdf']['SEEDLOT_NUMBER'][iAT])
	
		if np.sum(cd0)==0:
			N_no_spc=N_no_spc+1
	
		if np.sum(pct0)==0:
			N_no_num_planted=N_no_num_planted+1
			pct0=(100/pct0.size)*np.ones(pct0.size)
	
		# Genetic worth
		gw0=np.zeros(pct0.size)
		for j in range(pct0.size):
			ind=np.where(dGW['SEEDLOT_NUMBER']==sn0[j])[0]
			if ind.size!=0:
				gw0[j]=dGW['GENETIC_WORTH_RTNG'][ind[0]]
			else:
				N_sn_not_found=N_sn_not_found+1
	
		# Dissolve to unique species combinations and calculate weighted average genetic worth
		cd1=np.unique(cd0)
		pct1=np.zeros(cd1.size)
		gw1=np.zeros(cd1.size)
		for iCD in range(cd1.size):
			ind=np.where(cd0==cd1[iCD])[0]
			pct1[iCD]=np.sum(pct0[ind])
			gw1[iCD]=np.sum(pct0[ind]*gw0[ind])/np.sum(pct0[ind])

		# Divi
		gw1=np.nan_to_num(gw1)

		# Sort
		iSort=np.argsort(-pct1)
		cd1=cd1[iSort]
		pct1=pct1[iSort]
		gw1=gw1[iSort]

		# Percents will be stored as integers, ensure they will sum to 100
		pct1=np.round(pct1,decimals=0).astype('int16')
		if (np.sum(pct1)!=100) & (pct1.size>1):
			pct1=round_retain_sum(pct1)

		# Populate structure
		if cd1.size>0:
			plAT['PL_SPECIES_CD1'][iU]=cd1[0]
			plAT['PL_SPECIES_PCT1'][iU]=pct1[0]
			plAT['PL_SPECIES_GW1'][iU]=gw1[0]
		if cd1.size>1:
			plAT['PL_SPECIES_CD2'][iU]=cd1[1]
			plAT['PL_SPECIES_PCT2'][iU]=pct1[1]
			plAT['PL_SPECIES_GW2'][iU]=gw1[1]
		if cd1.size>2:
			plAT['PL_SPECIES_CD3'][iU]=cd1[2]
			plAT['PL_SPECIES_PCT3'][iU]=pct1[2]
			plAT['PL_SPECIES_GW3'][iU]=gw1[2]
		if cd1.size>3:
			plAT['PL_SPECIES_CD4'][iU]=cd1[3]
			plAT['PL_SPECIES_PCT4'][iU]=pct1[3]
			plAT['PL_SPECIES_GW4'][iU]=gw1[3]
		if cd1.size>4:
			plAT['PL_SPECIES_CD5'][iU]=cd1[4]
			plAT['PL_SPECIES_PCT5'][iU]=pct1[4]
			plAT['PL_SPECIES_GW5'][iU]=gw1[4]
		if cd1.size>5:
			plAT['PL_SPECIES_CD6'][iU]=cd1[5]
			plAT['PL_SPECIES_PCT6'][iU]=pct1[5]
			plAT['PL_SPECIES_GW6'][iU]=gw1[5]
	
	# Populate packed layers
	N_Year=6
	lNam='RSLT_ACTIVITY_TREATMENT_SVW'
	vNam='PL_All'
	
	for iY in range(N_Year):
		t0=time.time()
		z=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\' + lNam + '\\' + vNam + '_' + str(iY+1) + '_ACTIVITY_TREATMENT_UNIT_ID.tif')
		uAT1=gu.IndicesFromUniqueArrayValues(z['Data'].flatten())
	
		z0={}
		for iS in range(6):
			z0['PL_SPECIES_CD' + str(iS+1)]=np.zeros(zRef['Data'].shape,dtype='int16').flatten()
			z0['PL_SPECIES_PCT' + str(iS+1)]=np.zeros(zRef['Data'].shape,dtype='int16').flatten()
			z0['PL_SPECIES_GW' + str(iS+1)]=np.zeros(zRef['Data'].shape,dtype='int16').flatten()
	
		for k in uAT1.keys():
			ind=np.where(plAT['ACTIVITY_TREATMENT_UNIT_ID']==k)[0]
			if ind.size==0:
				continue
			iIndex=uAT1[k]
			for iS in range(6):
				z0['PL_SPECIES_CD' + str(iS+1)][iIndex]=plAT['PL_SPECIES_CD' + str(iS+1)][ind]
				z0['PL_SPECIES_PCT' + str(iS+1)][iIndex]=plAT['PL_SPECIES_PCT' + str(iS+1)][ind]
				z0['PL_SPECIES_GW' + str(iS+1)][iIndex]=plAT['PL_SPECIES_GW' + str(iS+1)][ind]
	
		for iS in range(6):
			z=copy.deepcopy(zRef)
			z['Data']=np.reshape(z0['PL_SPECIES_CD' + str(iS+1)],(zRef['Data'].shape))
			gis.SaveGeoTiff(z,meta['Paths']['bc1ha'] + '\\' + lNam + '\\' + vNam + '_' + str(iY+1) + '_PL_SPECIES_CD' + str(iS+1) + '.tif')
			z=copy.deepcopy(zRef)
			z['Data']=np.reshape(z0['PL_SPECIES_PCT' + str(iS+1)],(zRef['Data'].shape))
			gis.SaveGeoTiff(z,meta['Paths']['bc1ha'] + '\\' + lNam + '\\' + vNam + '_' + str(iY+1) + '_PL_SPECIES_PCT' + str(iS+1) + '.tif')
			z=copy.deepcopy(zRef)
			z['Data']=np.reshape(z0['PL_SPECIES_GW' + str(iS+1)],(zRef['Data'].shape))
			gis.SaveGeoTiff(z,meta['Paths']['bc1ha'] + '\\' + lNam + '\\' + vNam + '_' + str(iY+1) + '_PL_SPECIES_GW' + str(iS+1) + '.tif')
	
		print((time.time()-t0)/60)
	return

#%% Global Forest Change Loss Year (adjusted to remove known disturbances)
def FilterGFC_LossYear(meta):
	zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
	vList=['gfcly','harv_yr_con1','fire_yr','ibm_yr']
	z=Import_Raster(meta,[],vList)
	ind=np.where( (np.abs(z['gfcly']['Data']-z['harv_yr_con1']['Data'])<3) | (np.abs(z['gfcly']['Data']-z['fire_yr']['Data'])<2) | (np.abs(z['gfcly']['Data']-z['ibm_yr']['Data'])<2) )
	z['gfcly']['Data'][ind]=0
	#plt.matshow(z['gfcly']['Data'])
	gis.SaveGeoTiff(z['gfcly'],meta['Paths']['bc1ha'] + '\\Disturbances\\GlobalForestChange_LossYear_2021_Filtered.tif')
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
	 
	z=Import_Raster(meta,[],['refg','lc_comp1_2019'])   
	
	if mask=='Province':
		zMask=z['refg']
	elif mask=='BCFCS_LUC':
		# Land use change
		zMask=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Disturbances\\LandCover_Comp1_DeforestationMaskAll.tif')
	elif mask=='NOSE':
		# Non-obligation stand establishment
		zMask=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Management\\PL_NonOb_MaskAll.tif')
	elif mask=='BCFCS_NMC':
		# Nutrient management
		zMask=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\FE-CA_MaskAll.tif')
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
	vNam='FIRE_YEAR'
	N_Year=6
	d=[None]*N_Year
	for iY in range(N_Year):
		z0=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\' + lNam + '\\' + vNam + '_' + str(iY+1) + '_Year.tif')
		z0=gis.UpdateGridCellsize(z0,rgsf)
		d[iY]=z0['Data'][iMask]
	gu.opickle(meta['Paths']['bc1ha'] + '\\Sparse\\RGSF' + str(rgsf) + '_Mask' + mask + '_' + lNam + '_Year.pkl',d)
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
	
	# Harvest
	N_Year=3
	d=[None]*N_Year
	for iY in range(N_Year):
		z0=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VEG_CONSOLIDATED_CUT_BLOCKS_SP\\HARVEST_YEAR_' + str(iY+1) + '_Year.tif')
		z0=gis.UpdateGridCellsize(z0,rgsf)
		d[iY]=z0['Data'][iMask]
	gu.opickle(meta['Paths']['bc1ha'] + '\\Sparse\\RGSF' + str(rgsf) + '_Mask' + mask + '_VEG_CONSOLIDATED_CUT_BLOCKS_SP_Year.pkl',d)
	
	# Planting	   
	N_Year=6
	d=[None]*N_Year
	for iY in range(N_Year):
		z0=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(iY+1) + '_Year.tif')
		z0=gis.UpdateGridCellsize(z0,rgsf)
		d[iY]=z0['Data'][iMask]
	gu.opickle(meta['Paths']['bc1ha'] + '\\Sparse\\RGSF' + str(rgsf) + '_Mask' + mask + '_PL_All_Year.pkl',d)
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
		z0=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(iY+1) + '_RegenType.tif')
		z0=gis.UpdateGridCellsize(z0,rgsf)
		d[iY]=z0['Data'][iMask]
	gu.opickle(meta['Paths']['bc1ha'] + '\\Sparse\\RGSF' + str(rgsf) + '_Mask' + mask + '_PL_All_RegenType.pkl',d)
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
	
	# Fertilization
	N_Year=3
	d=[None]*N_Year
	dFSC=[None]*N_Year
	for iY in range(N_Year):
		z0=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\FE-CA_' + str(iY+1) + '_Year.tif')
		z0=gis.UpdateGridCellsize(z0,rgsf)
		d[iY]=z0['Data'][iMask]
		z0=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\FE-CA_' + str(iY+1) + '_SILV_FUND_SOURCE_CODE.tif')
		z0=gis.UpdateGridCellsize(z0,rgsf)
		dFSC[iY]=z0['Data'][iMask]
	gu.opickle(meta['Paths']['bc1ha'] + '\\Sparse\\RGSF' + str(rgsf) + '_Mask' + mask + '_FE-CA_Year.pkl',d)
	gu.opickle(meta['Paths']['bc1ha'] + '\\Sparse\\RGSF' + str(rgsf) + '_Mask' + mask + '_FE-CA_SILV_FUND_SOURCE_CODE.pkl',dFSC)
	
	# Knockdown
	N_Year=3
	d=[None]*N_Year
	dFSC=[None]*N_Year
	for iY in range(N_Year):
		z0=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\SP_KD_' + str(iY+1) + '_Year.tif')
		z0=gis.UpdateGridCellsize(z0,rgsf)
		d[iY]=z0['Data'][iMask]
		z0=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\SP_KD_' + str(iY+1) + '_SILV_FUND_SOURCE_CODE.tif')
		z0=gis.UpdateGridCellsize(z0,rgsf)
		dFSC[iY]=z0['Data'][iMask]
	gu.opickle(meta['Paths']['bc1ha'] + '\\Sparse\\RGSF' + str(rgsf) + '_Mask' + mask + '_SP_KD_Year.pkl',d)
	gu.opickle(meta['Paths']['bc1ha'] + '\\Sparse\\RGSF' + str(rgsf) + '_Mask' + mask + '_SP_KD_SILV_FUND_SOURCE_CODE.pkl',dFSC)
	
	# Prescribed burning
	N_Year=3
	d=[None]*N_Year
	dFSC=[None]*N_Year
	for iY in range(N_Year):
		z0=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\SP-BU-BROAD_' + str(iY+1) + '_Year.tif')
		z0=gis.UpdateGridCellsize(z0,rgsf)
		d[iY]=z0['Data'][iMask]		
	gu.opickle(meta['Paths']['bc1ha'] + '\\Sparse\\RGSF' + str(rgsf) + '_Mask' + mask + '_SP-BU-BROAD_Year.pkl',d)
	#gu.opickle(meta['Paths']['bc1ha'] + '\\Sparse\\RGSF' + str(rgsf) + '_Mask' + mask + '_SP_KD_SILV_FUND_SOURCE_CODE.pkl',dFSC)
	
	return

#%% Rasterize OPENING ID from OPENING LAYER (1 hour)
def RasterizeOpeningID(meta):
	
	# Import opening ID with spatial
	zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
	RasterizeFromSource(meta,zRef,'RSLT_OPENING_SVW','OPENING_ID')
	
	# Open Opening ID 1
	zOP1=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_OPENING_SVW\\OPENING_ID.tif')['Data']
	iOP1=np.where(zOP1!=0)
	zOP2=zOP1[iOP1]
	
	# Import opening spatial (takes 26 min)
	t0=time.time()
	ops={}
	ops['Path']=meta['Paths']['GDB']['Results']
	ops['Layer']='RSLT_OPENING_SVW'; 
	ops['crs']=meta['Geos']['crs']
	ops['Keep Geom']='On'
	ops['Select Openings']=np.array([])
	ops['SBC']=np.array([])
	ops['STC']=np.array([])
	ops['SMC']=np.array([])
	ops['FSC']=np.array([])
	ops['SOC1']=np.array([])
	ops['ROI']=[]
	ops['gdf']=qgdb.Query_Openings(ops,[])
	ops['gdf']=ops['gdf'][ops['gdf'].geometry!=None]
	ops['gdf']=ops['gdf'].reset_index()   
	print((time.time()-t0)/60)
	
	flg=np.zeros(ops['gdf']['OPENING_ID'].size,dtype=int)
	for iOP in range(ops['gdf']['OPENING_ID'].size):
		ind=np.where(zOP2==ops['gdf']['OPENING_ID'][iOP])[0]
		if ind.size==0:
			flg[iOP]=1
	ikp=np.where(flg==1)[0]
	ops0=ops['gdf'].iloc[ikp].copy()
	ops0=ops0.reset_index()
	
	zOP0=np.zeros(zRef['Data'].shape,dtype=float)
	shapes=((geom,value) for geom, value in zip(ops0['geometry'],ops0['OPENING_ID']))
	burned=features.rasterize(shapes=shapes,fill=0,out=zOP0,transform=zRef['Transform'])
	
	ind=np.where(zOP0!=0)
	zOP2=copy.deepcopy(zRef)
	zOP2['Data']=np.zeros(zRef['Data'].shape,dtype='int32')
	zOP2['Data'][ind]=zOP0[ind]#.astype('int32')
	gis.SaveGeoTiff(zOP2,meta['Paths']['bc1ha'] + '\\RSLT_OPENING_SVW\\OPENING_ID_2.tif')	
	
	return

#%% Derive Artificial Stand Establishment Type (used for non-ob stand establishment)
def DeriveASETComp1(meta):
	meta=cbu.Load_LUTs_Modelling(meta)
	
	d={}
	d['tv']=np.arange(1960,2023,1)
	ptNam=np.array(list(meta['LUT']['Derived']['ASET'].keys()))
	ptID=np.array(list(meta['LUT']['Derived']['ASET'].values()))
	zD=Import_Raster(meta,[],['harv_yr_con1','kd_yr','pdead_cruise'])
	
	zP={}
	for iP in range(6):
		zP[iP]={}
		zP[iP]['Type']=np.zeros(zRef['Data'].shape,dtype='int8')
		zP[iP]['Year']=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(iP+1) + '_Year.tif')['Data']
		zP[iP]['STC']=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(iP+1) + '_SILV_TECHNIQUE_CODE.tif')['Data']
		zP[iP]['FSC']=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(iP+1) + '_SILV_FUND_SOURCE_CODE.tif')['Data']
	
	zDS={}
	for iDS in range(3):
		zDS[iDS]={}
		zDS[iDS]['Year']=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\DS_All_' + str(iP+1) + '_Year.tif')['Data']
	
	# Summarize frequency of each Stocking Type Code (All and NOSE)
	# d['STC Summary']={}
	# d['STC Summary']['stc']=['PL','RP','FP','CG','RO','RR','SE','SL']
	# d['STC Summary']['N']=np.zeros((len(d['STC Summary']['stc']),2))
	# for v in d['STC Summary']['stc']:
	#	 for iP in range(6):
	#		 ind=np.where( (zP[iP]['STC']==meta['LUT']['RSLT_ACTIVITY_TREATMENT_SVW']['SILV_TECHNIQUE_CODE'][v]) )
	#		 d['STC Summary']['N'][cnt,0]=d['STC Summary']['N'][cnt,0]+ind[0].size
	#		 ind=np.where( (zP[iP]['STC']==meta['LUT']['RSLT_ACTIVITY_TREATMENT_SVW']['SILV_TECHNIQUE_CODE'][v]) & (np.isin(zP[iP]['FSC'],meta['Param']['BE']['FSC']['NO List ID'])==True) )
	#		 d['STC Summary']['N'][cnt,1]=d['STC Summary']['N'][cnt,1]+ind[0].size
	
	cd=meta['LUT']['RSLT_ACTIVITY_TREATMENT_SVW']['SILV_TECHNIQUE_CODE']
	stcInclude=[cd['PL'],cd['RP'],cd['FP'],cd['RR']]
	
	d['A Tot']=np.zeros(d['tv'].size)
	d['A']=np.zeros((d['tv'].size,ptNam.size))
	d['A NOSE']=np.zeros((d['tv'].size,ptNam.size))
	DL=np.zeros(zRef['Data'].shape,dtype='int8')
	for iT in range(d['tv'].size):
		print(d['tv'][iT])
		zWF=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\PROT_HISTORICAL_FIRE_POLYS_SP\\FIRE_YEAR_' + str(d['tv'][iT]-1) + '.tif')['Data']
		zIBM=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\PEST_INFESTATION_POLY\\PEST_SEVERITY_CODE_IBM_' + str(d['tv'][iT]-1) + '.tif')['Data']
		ind=np.where(zIBM>=3); DL[ind]=meta['LUT']['Event']['IBM'] # Severity greater than light
		ind=np.where(zD['harv_yr_con1']['Data']==d['tv'][iT]-1); DL[ind]=meta['LUT']['Event']['Harvest']; #H_Yr[ind]=zD['harv_yr_con1']['Data'][ind]
		ind=np.where(zD['kd_yr']['Data']==d['tv'][iT]-1); DL[ind]=meta['LUT']['Event']['Knockdown']
		ind=np.where(zWF>0); DL[ind]=meta['LUT']['Event']['Wildfire']
		
		pl_oc=np.zeros(zRef['Data'].shape,dtype='int8')
		pl_stc=np.zeros(zRef['Data'].shape,dtype='int16')
		pl_fsc=np.zeros(zRef['Data'].shape,dtype='int16')
		for iP in range(6):
			ind=np.where( (zP[iP]['Year']==d['tv'][iT]) & (np.isin(zP[iP]['STC'],stcInclude)==True) )			
			pl_oc[ind]=1
			pl_stc[ind]=zP[iP]['STC'][ind]
			pl_fsc[ind]=zP[iP]['FSC'][ind]
		
		ind=np.where( (pl_oc==1) )
		d['A Tot'][iT]=ind[0].size
		
		zType0=np.zeros(zRef['Data'].shape,dtype='int8')
		
		# Replanting
		ind=np.where( (pl_oc==1) & (pl_stc==meta['LUT']['RSLT_ACTIVITY_TREATMENT_SVW']['SILV_TECHNIQUE_CODE']['RP']) & (np.isin(pl_fsc,meta['Param']['BE']['FSC']['NO List ID'])==True) )
		d['A NOSE'][iT,meta['LUT']['Derived']['ASET']['Replanting']-1]=ind[0].size
		ind=np.where( (pl_oc==1) & (pl_stc==meta['LUT']['RSLT_ACTIVITY_TREATMENT_SVW']['SILV_TECHNIQUE_CODE']['RP']) )
		d['A'][iT,meta['LUT']['Derived']['ASET']['Replanting']-1]=ind[0].size
		zType0[ind]=meta['LUT']['Derived']['ASET']['Replanting']
		
		# Fill Planting
		ind=np.where( (pl_oc==1) & (pl_stc==meta['LUT']['RSLT_ACTIVITY_TREATMENT_SVW']['SILV_TECHNIQUE_CODE']['FP']) & (np.isin(pl_fsc,meta['Param']['BE']['FSC']['NO List ID'])==True) )
		d['A NOSE'][iT,meta['LUT']['Derived']['ASET']['Fill Planting']-1]=ind[0].size
		ind=np.where( (pl_oc==1) & (pl_stc==meta['LUT']['RSLT_ACTIVITY_TREATMENT_SVW']['SILV_TECHNIQUE_CODE']['FP']) )		
		d['A'][iT,meta['LUT']['Derived']['ASET']['Fill Planting']-1]=ind[0].size
		zType0[ind]=meta['LUT']['Derived']['ASET']['Fill Planting']
		
		# Road rehab
		ind=np.where( (pl_oc==1) & (pl_stc==meta['LUT']['RSLT_ACTIVITY_TREATMENT_SVW']['SILV_TECHNIQUE_CODE']['RR']) & (np.isin(pl_fsc,meta['Param']['BE']['FSC']['NO List ID'])==True) )
		d['A NOSE'][iT,meta['LUT']['Derived']['ASET']['Road Rehabilitation']-1]=ind[0].size
		ind=np.where( (pl_oc==1) & (pl_stc==meta['LUT']['RSLT_ACTIVITY_TREATMENT_SVW']['SILV_TECHNIQUE_CODE']['RR']) )
		d['A'][iT,meta['LUT']['Derived']['ASET']['Road Rehabilitation']-1]=ind[0].size
		zType0[ind]=meta['LUT']['Derived']['ASET']['Road Rehabilitation']
		
		# Back to back planting
		ind=np.where( (pl_oc==1) & (zType0==0) & (DL==meta['LUT']['Event']['Planting']) & (np.isin(pl_fsc,meta['Param']['BE']['FSC']['NO List ID'])==True) )
		d['A NOSE'][iT,meta['LUT']['Derived']['ASET']['Back-to-back Planting']-1]=ind[0].size	
		ind=np.where( (pl_oc==1) & (zType0==0) & (DL==meta['LUT']['Event']['Planting']) )		
		d['A'][iT,meta['LUT']['Derived']['ASET']['Back-to-back Planting']-1]=ind[0].size	
		zType0[ind]=meta['LUT']['Derived']['ASET']['Back-to-back Planting']
		
		# Salvage
		ind=np.where( (pl_oc==1) & (d['tv'][iT]>=2004) & (zD['harv_yr_con1']['Data']>1987) & (zType0==0) & (DL==meta['LUT']['Event']['Harvest']) & (np.isin(pl_fsc,meta['Param']['BE']['FSC']['NO List ID'])==True) )
		d['A NOSE'][iT,meta['LUT']['Derived']['ASET']['Salvage and Planting']-1]=ind[0].size
		ind=np.where( (pl_oc==1) & (d['tv'][iT]>=2004) & (zD['harv_yr_con1']['Data']>1987) & (zType0==0) & (DL==meta['LUT']['Event']['Harvest']) & (np.isin(pl_fsc,meta['Param']['BE']['FSC']['NO List ID'])==True) | \
					 (pl_oc==1) & (zType0==0) & (DL==meta['LUT']['Event']['Harvest']) & (zD['pdead_cruise']['Data']>=50) )
		d['A'][iT,meta['LUT']['Derived']['ASET']['Salvage and Planting']-1]=ind[0].size		
		zType0[ind]=meta['LUT']['Derived']['ASET']['Salvage and Planting']
		
		# Knockdown
		ind=np.where( (pl_oc==1) & (zType0==0) & (DL==meta['LUT']['Event']['Knockdown']) & (np.isin(pl_fsc,meta['Param']['BE']['FSC']['NO List ID'])==True) )
		d['A NOSE'][iT,meta['LUT']['Derived']['ASET']['Knockdown and Planting']-1]=ind[0].size
		ind=np.where( (pl_oc==1) & (zType0==0) & (DL==meta['LUT']['Event']['Knockdown']) )
		d['A'][iT,meta['LUT']['Derived']['ASET']['Knockdown and Planting']-1]=ind[0].size		
		zType0[ind]=meta['LUT']['Derived']['ASET']['Knockdown and Planting']
		
		# Straight fire
		ind=np.where( (pl_oc==1) & (zType0==0) & (DL==meta['LUT']['Event']['Wildfire']) & (np.isin(pl_fsc,meta['Param']['BE']['FSC']['NO List ID'])==True) )
		d['A NOSE'][iT,meta['LUT']['Derived']['ASET']['Straight-to-planting Post Wildfire']-1]=ind[0].size
		ind=np.where( (pl_oc==1) & (zType0==0) & (DL==meta['LUT']['Event']['Wildfire']) )
		d['A'][iT,meta['LUT']['Derived']['ASET']['Straight-to-planting Post Wildfire']-1]=ind[0].size		
		zType0[ind]=meta['LUT']['Derived']['ASET']['Straight-to-planting Post Wildfire']		
		
		# Straight insect
		ind=np.where( (pl_oc==1) & (zType0==0) & (DL==meta['LUT']['Event']['IBM']) & (np.isin(pl_fsc,meta['Param']['BE']['FSC']['NO List ID'])==True) )
		d['A NOSE'][iT,meta['LUT']['Derived']['ASET']['Straight-to-planting Post Beetles']-1]=ind[0].size
		ind=np.where( (pl_oc==1) & (zType0==0) & (DL==meta['LUT']['Event']['IBM']) )
		d['A'][iT,meta['LUT']['Derived']['ASET']['Straight-to-planting Post Beetles']-1]=ind[0].size
		zType0[ind]=meta['LUT']['Derived']['ASET']['Straight-to-planting Post Beetles']
		
		# Straight other
		
		# Harvest and Planting NSR backlog
		ind=np.where( (pl_oc==1) & (zType0==0) & (DL==meta['LUT']['Event']['Harvest']) & (zD['harv_yr_con1']['Data']<=1987) & (np.isin(pl_fsc,meta['Param']['BE']['FSC']['NO List ID'])==True) )
		d['A NOSE'][iT,meta['LUT']['Derived']['ASET']['Harvest and Planting NSR Backlog']-1]=ind[0].size		
		d['A'][iT,meta['LUT']['Derived']['ASET']['Harvest and Planting NSR Backlog']-1]=ind[0].size		
		zType0[ind]=meta['LUT']['Derived']['ASET']['Harvest and Planting NSR Backlog']
		
		# Harvest and Planting
		ind=np.where( (pl_oc==1) & (zType0==0) & (DL==meta['LUT']['Event']['Harvest']) & (np.isin(pl_fsc,meta['Param']['BE']['FSC']['NO List ID'])==True) )
		d['A NOSE'][iT,meta['LUT']['Derived']['ASET']['Harvest and Planting']-1]=ind[0].size 
		ind=np.where( (pl_oc==1) & (zType0==0) & (DL==meta['LUT']['Event']['Harvest']) )
		d['A'][iT,meta['LUT']['Derived']['ASET']['Harvest and Planting']-1]=ind[0].size		 
		zType0[ind]=meta['LUT']['Derived']['ASET']['Harvest and Planting']
		
		# Unknown
		ind=np.where( (pl_oc==1) & (zType0==0) & (np.isin(pl_fsc,meta['Param']['BE']['FSC']['NO List ID'])==True) )
		d['A NOSE'][iT,meta['LUT']['Derived']['ASET']['Unknown']-1]=ind[0].size
		ind=np.where( (pl_oc==1) & (zType0==0) )
		d['A'][iT,meta['LUT']['Derived']['ASET']['Unknown']-1]=ind[0].size		
		zType0[ind]=meta['LUT']['Derived']['ASET']['Unknown']
		
		# Update last disturbance if planting occurs
		DL[(pl_oc==1)]=meta['LUT']['Event']['Planting']
		
		# Direct seeding
		for iDS in range(3):
			ind=np.where( (zDS[iDS]['Year']==d['tv'][iT]) & (np.isin(pl_fsc,meta['Param']['BE']['FSC']['NO List ID'])==True) )
			d['A NOSE'][iT,meta['LUT']['Derived']['ASET']['Direct Seeding']-1]=ind[0].size
			ind=np.where( (zDS[iDS]['Year']==d['tv'][iT]) )			
			d['A'][iT,meta['LUT']['Derived']['ASET']['Direct Seeding']-1]=ind[0].size
			zType0[ind]=meta['LUT']['Derived']['ASET']['Direct Seeding']
		
		# Pack
		for iP in range(6):
			ind=np.where( (zP[iP]['Year']==d['tv'][iT]) )
			zP[iP]['Type'][ind]=zType0[ind]
	
	# Save summary
	gu.opickle(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\ASET_Summary_All',d)
	
	# Save packed regen type
	for i in range(6):
		z1=copy.deepcopy(zRef)
		z1['Data']=zP[i]['Type']
		gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(i+1) + '_ASET.tif')
	
	return

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


#%%
def RasterizeInsectComp1(meta):
	zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
	#meta['LUT']['PEST_INFESTATION_POLY'].keys()
	
	# Import vector data
	lNam='PEST_INFESTATION_POLY'
	vNam='PEST_SEVERITY_CODE'
	if os.path.exists(meta['Paths']['bc1ha'] + '\\' + lNam)==False:
		os.mkdir(meta['Paths']['bc1ha'] + '\\' + lNam)	
	indI=np.where( (meta['Geos']['Variable Info']['Layer Name']==lNam) & (meta['Geos']['Variable Info']['Variable Name']==vNam) )[0]
	pthin=meta['Paths']['GDB'][ meta['Geos']['Variable Info']['Geodatabase Name'][indI[0]] ]
	df=gpd.read_file(pthin,layer=lNam)
	df=df[df.geometry!=None]
	df=df.reset_index()
	df=CreateIdForCategoricalVariable(meta,lNam,vNam,df)

	# Parameters
	pestL=['IBB','IBD','IBM','IBS','IDW','IDL']
	sevL=['M','S','V']
	tv=np.arange(1951,2023,1)

	# Initialize rasters
	N_Year=10
	z={'Year':{},'Type':{}}
	for iY in range(N_Year):
		z['Year'][iY+1]=np.zeros(zRef['Data'].shape,dtype='int16')
		z['Type'][iY+1]=np.zeros(zRef['Data'].shape,dtype='int8')

	for iT in range(tv.size):
		print(tv[iT])
		for pest in pestL:
			for sev in sevL:
				ind=np.where( (meta['Param']['BE']['InsectComp1']['Insect Code']==pest) & (meta['Param']['BE']['InsectComp1']['Severity Code']==sev) )[0][0]
				id=meta['Param']['BE']['InsectComp1']['ID'][ind]
				
				df0=df[ (df['PEST_SPECIES_CODE']==pest) & (df['PEST_SEVERITY_CODE']==sev) & (df['CAPTURE_YEAR']==tv[iT]) ].copy()
				df0=df0[df0.geometry!=None]
				df0=df0.reset_index()
				df0['ID']=np.ones(len(df0['PEST_SPECIES_CODE']))
				if len(df0)>0:
					shapes=((geom,value) for geom, value in zip(df0.geometry,df0['ID']))
					z0=np.zeros(zRef['Data'].shape,dtype=float)
					burned=features.rasterize(shapes=shapes,fill=0,out=z0,transform=zRef['Transform'])
				else:
					z0=np.zeros(zRef['Data'].shape,dtype=float)
				
				iY=1; ind=np.where((z['Year'][iY]==0) & (z0!=0))
				z['Year'][iY][ind]=tv[iT]
				z['Type'][iY][ind]=id
				for iY in range(2,N_Year): 
					ind=np.where((z['Year'][iY-1]!=0) & (z['Year'][iY-1]!=tv[iT]) & (z['Year'][iY]==0) & (z0!=0))
					z['Year'][iY][ind]=tv[iT]
					z['Type'][iY][ind]=id

	# Save
	for iY in range(N_Year):
		z1=copy.deepcopy(zRef)
		z1['Data']=z['Year'][iY+1].astype('int16')
		gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\' + lNam + '\\InsectComp1_' + str(iY+1) + '_Year.tif')
		z1=copy.deepcopy(zRef)
		z1['Data']=z['Type'][iY+1].astype('int8')
		gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\' + lNam + '\\InsectComp1_' + str(iY+1) + '_Type.tif')
	return

#%% Calculate time series summaries for insect comp 1
def CalcInsectComp1_TimeSeries(meta):
	zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
	tv=np.arange(1951,2024,1)
	id=meta['Param']['Raw']['InsectComp1']['ID']
	name=meta['Param']['Raw']['InsectComp1']['Insect Name']
	sev=meta['Param']['Raw']['InsectComp1']['Severity Code']
	uNam=np.unique(name)
	d={}
	for nam in uNam:
		d[nam]=np.zeros((tv.size,4))
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
			iSev=np.where( np.array(['L','M','S','V'])==sev[iID] )[0]
			ind1=np.where(id0==u)[0]
			yr1=yr0[ind1]
			idx=gu.IndicesFromUniqueArrayValues(yr1)
			for i in idx.keys():
				iT=np.where(tv==i)[0]
				d[name[iID[0]]][iT,iSev]=d[name[iID[0]]][iT,iSev]+idx[i].size
	gu.opickle(meta['Paths']['bc1ha'] + '\\PEST_INFESTATION_POLY\\InsectComp1_SummaryTS.pkl',d)
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
		gis.SaveGeoTiff(z1[i],meta['Paths']['bc1ha'] + '\\PEST_INFESTATION_POLY\\InsectComp1_' + i + '_YearLast_Severe.tif')

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



#%% Rasterize Watersheds (Named)
def RasterizeWatersheds(meta):
	# Import reference grid
	zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])

	pthin=meta['Paths']['GDB']['GDB'] + '\\WatershedsNamed.gdb'
	# fiona.listlayers(pthin)
	df=gpd.read_file(pthin,layer='FWA_NAMED_WATERSHEDS_POLY')

	df=df[df.geometry!=None]
	df=df.reset_index()

	#z1=copy.deepcopy(zRef)
	#z1['Data']=np.zeros(zRef['Data'].shape,dtype='int32')
	ord=np.flip(np.arange(1,11))
	for iO in ord:
		print(iO)
		df0=df[df['STREAM_ORDER']==iO]
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
	Name='Big Creek'
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

#%% Artificial Stand Establishment Type (ASET)
def DeriveLastASET(meta):
	zH=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VEG_CONSOLIDATED_CUT_BLOCKS_SP\\HARVEST_YEAR_YearLast.tif')
	zPD=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Disturbances\\HarvestPercentDead_FromCruise.tif')

	z1=copy.deepcopy(zRef)
	z1['Data']=np.zeros(zRef['Data'].shape,dtype='int16')
	for i in range(6):
		zY=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(i+1) + '_Year.tif')
		zT=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(i+1) + '_RegenType.tif')
		ind=np.where(zY['Data']>0)
		z1['Data'][ind]=zT['Data'][ind]
	# Add salvage
	ind=np.where( (zH['Data']>=2017) & (zPD['Data']>=15) & (z1['Data']!=meta['LUT']['Derived']['ASET']['Replanting']) )
	z1['Data'][ind]=meta['LUT']['Derived']['ASET']['Salvage and Planting']
	gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_ASET_All.tif')

	# Post 2017
	z1=copy.deepcopy(zRef)
	z1['Data']=np.zeros(zRef['Data'].shape,dtype='int16')
	for i in range(6):
		zY=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(i+1) + '_Year.tif')
		zT=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(i+1) + '_RegenType.tif')
		ind=np.where(zY['Data']>=2017)
		z1['Data'][ind]=zT['Data'][ind]
	# Add salvage
	ind=np.where( (zH['Data']>=2017) & (zPD['Data']>=15) & (z1['Data']!=meta['LUT']['Derived']['ASET']['Replanting']) )
	z1['Data'][ind]=meta['LUT']['Derived']['ASET']['Salvage and Planting']
	gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_ASET_Post2017.tif')
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

