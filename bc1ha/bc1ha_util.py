
#%% Import modules

import os
import numpy as np
import gc
import matplotlib.pyplot as plt
import matplotlib.colors
from matplotlib.colors import LightSource
import geopandas as gpd
import pandas as pd
from scipy.interpolate import griddata
import copy
from shapely.geometry import Polygon,Point,box,shape
from rasterio import features
import fiona
import time
import cv2
from perlin_noise import PerlinNoise
from numpy.random import default_rng
import fcgadgets.macgyver.util_general as gu
import fcgadgets.macgyver.util_gis as gis
import fcgadgets.macgyver.util_query_gdb as qgdb
import fcexplore.field_plots.Processing.psp_util as ugp
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
	meta['Paths']['GDB']={}
	meta['Paths']['GDB']['GDB']=r'C:\Data\Geodatabases'
	meta['Paths']['GDB']['LandCover']=r'C:\Data\Geodatabases\LandCover\20230607\LandCover.gdb'
	meta['Paths']['GDB']['LandUse']=r'C:\Data\Geodatabases\LandUse\20230501\LandUse.gdb'
	meta['Paths']['GDB']['Disturbance']=r'C:\Data\Geodatabases\Disturbances\20230501\Disturbances.gdb'
	meta['Paths']['GDB']['Results']=r'C:\Data\Geodatabases\Results\20230430\Results.gdb'
	meta['Paths']['GDB']['VRI 2002']=r'C:\Data\Geodatabases\VRI\2002\VRI2002_VEG_COMP_LYR_R1_POLY_FINAL_DELIVERYV4.gdb'
	meta['Paths']['GDB']['VRI 2023']=r'C:\Data\Geodatabases\VRI\20230401\VRI.gdb'
	meta['Paths']['Model']={}
	#meta['Paths']['Model']['Code']=r'C:\Code_Python\fcgadgets\cbrunner'
	meta['Paths']['Model']['Code']=r'G:\My Drive\Code_Python\fcgadgets\cbrunner'
	meta['Paths']['Model']['Parameters']=meta['Paths']['Model']['Code'] + '\\Parameters'
	meta['Paths']['Model']['Taz Datasets']=r'C:\Data\Taz Datasets'
	meta['Paths']['GP']={}
	meta['Paths']['GP']['DB']=r'C:\Data\GroundPlots\PSP-NADB2'
	meta['Paths']['GP']['Raw Data']={}
	meta['Paths']['GP']['Raw Data']['BC']=meta['Paths']['GP']['DB'] + '\\Given\BC\Received 2023-03-02'
	
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
	meta['Graphics']['Map']['Show TPF and Cities']='Off'
	meta['Graphics']['Map']['Show Symbol Labels']='Off'
	
	meta['Graphics']['Flowchart']={}
	meta['Graphics']['Flowchart']['Font Name']='Arial'
	meta['Graphics']['Flowchart']['Font Color']='#6b7d2a'
	meta['Graphics']['Flowchart']['Font Size']='10'
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
	gdf_bm=gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='TA_REGIONAL_DISTRICTS_SVW')
	meta['Geos']['crs']=gdf_bm.crs
	
	# Import LUTs
	if 'LUT' not in meta:
		meta=ImportLUTs(meta)

	if 'Include GPs' in argv:
		meta['GP'],data,soc=ugp.ImportGroundPlotData(meta,type='Just Parameters',include_soil='True')
		meta['GP']['Data']=data
		meta['GP']['SOC']=soc

	# Import model parameters
	meta=cbu.ImportParameters(meta)

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

    # Harvest retention compilation
    d=gu.ReadExcel(meta['Paths']['Model']['Parameters'] + '\\LUT_harv_retention_comp1.xlsx')
    vNam='harvret1'
    meta['LUT']['Raw'][vNam]=d
    meta['LUT'][lNam][vNam]={}
    for i in range(d['ID'].size):
        meta['LUT'][lNam][vNam][ d['Name'][i] ]=d['ID'][i]
    
    # Harvest regeneration type
    d=gu.ReadExcel(meta['Paths']['Model']['Parameters'] + '\\LUT_HarvestRegenType.xlsx')
    vNam='HarvestRegenType'
    meta['LUT']['Raw'][vNam]=d
    meta['LUT'][lNam][vNam]={}
    for i in range(d['ID'].size):
        meta['LUT'][lNam][vNam][ d['Name'][i] ]=d['ID'][i]

    # Planting regeneration type (used in NOSE modelling)
    d=gu.ReadExcel(meta['Paths']['Model']['Parameters'] + '\\LUT_RegenType.xlsx')
    vNam='RegenType'
    meta['LUT']['Raw'][vNam]=d
    meta['LUT'][lNam][vNam]={}
    for i in range(d['ID'].size):
        meta['LUT'][lNam][vNam][ d['Name'][i] ]=d['ID'][i]

    #d=gu.ReadExcel(meta['Paths']['Model']['Parameters'] + '\\LUT_EcozoneCanada.xlsx')
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

    # Forest Cover Stocking Type
    #d=gu.ReadExcel(r'C:\Data\BC1ha\Results\stocktype.tif.vat.xlsx')
    #d=gu.ReadExcel(r'C:\Code_Python\fcgadgets\cbrunner\Parameters\Parameters_ForestCover_StockingType.xlsx')
    # nam='fcst'
    # meta['LUT'][nam]={}
    # for i in range(d['VALUE'].size):
    #     meta['LUT'][nam][ d['STOCKING_T'][i] ]=d['VALUE'][i]

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
    t0=time.time()
    lNam='VEG_COMP_LYR_R1_POLY'

    # Import reference grid
    zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])

    # Import feature ID
    zFID=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2002\\featid.tif')
    #zFID=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\vri_feaid.tif')
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

    with fiona.open(meta['Paths']['VRI'],layer=lNam) as source:
        for feat in source:
            prp=dict(feat['properties'].items())
            # if prp['BCLCS_LEVEL_2']!='T':
            #     continue

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
                        d[vNam][ind]=meta['lut'][lNam][vNam][ prp[vNam] ]
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
            #     continue

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
    if roi['Type']=='ByTSA':
        roi['gdf']['bound within']=gdf['tsa']['gdf'].iloc[np.isin(gdf['tsa']['gdf'].Name,roi['List'])]
    elif roi['Type']=='ByRegDis':
        roi['gdf']['bound within']=gdf['regdis']['gdf'].iloc[ np.isin(gdf['regdis']['gdf']['REGIONAL_DISTRICT_NAME'],roi['List']) ]
    else:
        roi['gdf']['bound within']=roi['gdf']['bound']
    for k in gdf.keys():
        roi['gdf'][k]=gdf[k]['gdf'].cx[roi['grd']['xmin']:roi['grd']['xmax'],roi['grd']['ymin']:roi['grd']['ymax']]
        roi['gdf'][k]=roi['gdf'][k].reset_index(drop=True)
        if (roi['Type']=='ByTSA') | (roi['Type']=='ByRegDis'):
            roi['gdf'][k]=gpd.overlay(roi['gdf'][k],roi['gdf']['bound within'],how='intersection')

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
def Import_GDB_Over_ROI(meta,roi,vList):
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

        if roi['Type']=='ByTSA':
            roi['gdf'][nam]['gdf']=gpd.overlay(roi['gdf'][nam]['gdf'],roi['gdf']['bound within'],how='intersection')

    return roi

#%% Import variables for ROI
def Import_Raster(meta,roi,vList,*argv):
    d={}
    for v in vList:
        if roi!=[]:
            if 'grd' in roi.keys():
                if v in roi['grd'].keys():
                    continue
        if v=='age_ntems':
            d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Age\\Age_NTEM_2019.tif')
        elif v=='age_vri':
            d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\PROJ_AGE_1.tif')
        elif v=='age_vri02':
            d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2002\\PROJ_AGE_1.tif')
        elif v=='alr':
            d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\OATS_ALR_POLYS\\ALR_POLY_ID.tif')
        elif v=='BEC_ZONE_CODE':
            d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\BEC_ZONE_CODE.tif')
        elif v=='bgcz':
            d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\BEC_BIOGEOCLIMATIC_POLY\\ZONE_GapFilled.tif')
        elif v=='biomass_glob':
            d[v]=gis.OpenGeoTiff(r'C:\Data\Biomass\GlobBiomass\N80W140_agb pc8.tif')
            d[v]['Data']=np.squeeze(d[v]['Data'])
            d[v]['Data']=0.5*d[v]['Data']
        elif v=='bsr_yr':
            d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Disturbances\\BurnSevComp_YearLast.tif')
        elif v=='bsr_sc':
            d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Disturbances\\BurnSevComp_SevClassLast.tif')
        elif v=='btm':
            d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\landuse.btm.tif')
            d[v]['Data']=np.squeeze(d[v]['Data'])
            lut=pd.read_csv(r'C:\Data\BC1ha\LandCoverUse\landuse_btm_category_metadata.csv')
            cl=np.column_stack( (lut['C1'].values,lut['C2'].values,lut['C3'].values) )
            d[v]['Compressed']={}
            d[v]['Compressed']['Data'],d[v]['Compressed']['lab'],d[v]['Compressed']['cl1']=gis.CompressCats(d[v]['Data'],lut['Raster Value'].values,lut['PLU Label'].values,cl)
        elif v=='citym':
            d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\BC_MAJOR_CITIES_POINTS_500M\\NAME.tif')
        elif v=='crownc':
            d[v]=gis.OpenGeoTiff(r'C:\Data\BC1ha\VRI\crownc.tif')
        elif v=='d2road':
            d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Management\\DistanceFromRoads.tif')
        elif v=='d2fac':
            d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Management\\DistanceFromForestryFacility.tif')
        elif v=='d2set':
            d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Management\\DistanceFromSettlement.tif')
        elif v=='elev':
            d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Terrain\\elevation.tif')
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
        elif v=='height':
            d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\PROJ_HEIGHT_1.tif')
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
        elif v=='harvret1':
            d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Management\\HarvestRetentionComp1.tif')
        elif v=='harv_prob':
            d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Disturbances\\HarvestProbability.tif')
        elif v=='kd_yr':
            d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\SP_KD_YearLast.tif')
        elif v=='ibm_yr':
            d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\PEST_INFESTATION_POLY\\PEST_SEVERITY_CODE_IBM_YearLast.tif')
        elif v=='idw_mask':
            d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Disturbances\\IDW_Mask.tif')
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
        elif v=='plam':
            d[v]=gis.OpenGeoTiff(r'C:\Data\BC1ha\Results\Planting_FromRESULTS_MaskCount.tif')
        elif v=='popp':
            d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\NRC_POPULATED_PLACES_1M_SP\\NAME.tif')
        elif v=='prcp_ann_n':
            d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Climate\\Seasonal\\BC1ha_prcp_ann_norm_1971to2000_si_hist_v1.tif')
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
        elif v=='regentype':
            d[v]=gis.OpenGeoTiff(r'C:\Data\BC1ha\Results\Harvest_Regen_Type.tif')
        elif v=='rd':
            d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\TA_REGIONAL_DISTRICTS_SVW\\REGIONAL_DISTRICT_NAME.tif')
        elif v=='road_atl':
            d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\DRA_DGTL_ROAD_ATLAS_MPAR_SP\\DIGITAL_ROAD_ATLAS_LINE_ID.tif')
        elif v=='road_ften_s':
            d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\FTEN_ROAD_SEGMENT_LINES_SVW\\LIFE_CYCLE_STATUS_CODE.tif')
        elif v=='slope':
            d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Terrain\\Slope.tif')
        elif v=='si_vri23':
            d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\SITE_INDEX.tif')
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
        elif v=='si_spl_ntems':
            d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\SPL\\SI_From_NTEMS_Spc1_Filled_RD_CAPITAL_GF.tif')
        elif v=='si_spl_fd':
            d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\SPL\\Site_Prod_Fd.tif')
        elif v=='si_spl_hw':
            d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\SPL\\Site_Prod_Hw.tif')
        elif v=='tdc':
            d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\TreeDensityClass_Current.tif')
        elif v=='tdc_02':
            d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\TreeDensityClass_2002.tif')
        elif v=='tdc_wsg':
            d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\TreeDensityClass_WithShrubsGrasses.tif')
        elif v=='tmean_ann_n':
            d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Climate\\Seasonal\\BC1ha_tmean_ann_norm_1971to2000_si_hist_v1.tif')
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
            d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\BC_MAJOR_WATERSHEDS\\MAJOR_WATERSHED_SYSTEM.tif')
        elif v=='wbmm':
            d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\FWA_MANMADE_WATERBODIES_POLY\\WATERBODY_TYPE.tif')
        elif v=='ws_gs_n':
            d[v]=gis.OpenGeoTiff(r'C:\Data\BC1ha\Climate\Seasonal\wbm_ws_gs_norm_1971to2000.tif')
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

    # Get Current year
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
        tv,N=TimeSeriesOccurrenceFromPackedEvents(meta,lNam,vNam,6)
        plt.plot(tv,N,'-bo')
    return

#%% Current-year wildfire
# Current year + preveous year (sometimes missing)
def RasterizeWildfireCurrentYear(meta,yr):
    zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])    
    
    # Current year
    yr=2023
    pthin=r'C:\Data\Wildfire\Current Year\prot_current_fire_polys.shp'
    df=gpd.read_file(pthin)
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
    
    # Previous year when missing
    import geotable
    import pyproj
    srs=gis.ImportSRSs()
    crs=pyproj.CRS(srs['String']['Geographic'])
    
    a=geotable.load(r'C:\Data\Wildfire\BC Fire Perimeters 2020-2022.kmz')
    df=gpd.GeoDataFrame(data=a.Name,geometry=a.geometry_object)
    df['ID']=np.ones(len(df))
    df.crs=pyproj.CRS(srs['String']['Geographic'])
    df=df.to_crs({'init':'epsg:3005'})
    
    shapes=((geom,value) for geom, value in zip(df['geometry'],df['ID']))
    z0=np.zeros(zRef['Data'].shape,dtype=float)
    burned=features.rasterize(shapes=shapes,fill=0,out=z0,transform=zRef['Transform'])
    
    # Fix by removing other years
    tv=np.arange(2020,2021+1,1)
    for iT in range(tv.size):
        zF0=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\PROT_HISTORICAL_FIRE_POLYS_SP\\FIRE_YEAR_' + str(tv[iT]) + '.tif')
        ind=np.where(zF0['Data']>0)
        burned[ind]=0
    
    z1=copy.deepcopy(zRef)
    z1['Data']=np.zeros(zRef['Data'].shape,dtype='int16')
    ind=np.where(burned>0); z1['Data'][ind]=1
    print(np.sum(z1['Data']))
    gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\PROT_HISTORICAL_FIRE_POLYS_SP\\FIRE_YEAR_2022.tif')
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
        zFSC=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_ALL_' + str(iEY+1) + '_SILV_FUND_SOURCE_CODE.tif')
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
    #    ats={}
    #    ats['gdf']=gpd.read_file(r'C:\Data\BC1ha\RSLT_ACTIVITY_TREATMENT_SVW\ats' + vNam + '.geojson')

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
    #    at=gu.ipickle(r'C:\Data\BC1ha\RSLT_ACTIVITY_TREATMENT_SVW\at' + vNam + '.pkl')
    
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

#%% NALCMS CEC 2020 - Compress categories to remove tropics
def NALCMS_Compress(lut_in,zRef,zCEC10,zCEC20):

    cg=['Forest','Shrubland','Grassland','Lichen-Moss','Wetland','Cropland','Barren Ground','Urban','Water','Snow and Ice']
    d={'Value':(len(cg)+1)*np.ones(len(cg),dtype='uint8'),'Name':np.array(['empty' for _ in range(len(cg))],dtype=object)}
    lut_out={}
    cnt=1
    for k in cg:
        lut_out[k]=cnt
        d['Value'][cnt-1]=cnt
        d['Name'][cnt-1]=k
        cnt=cnt+1

    zCEC10c=copy.deepcopy(zRef)
    zCEC20c=copy.deepcopy(zRef)
    zCEC10c['Data']=(len(cg)+1)*np.ones(zRef['Data'].shape,dtype='uint8')
    zCEC20c['Data']=(len(cg)+1)*np.ones(zRef['Data'].shape,dtype='uint8')

    # Forest
    ind=np.where( (zCEC10['Data']==lut_in['Temperate or sub-polar needleleaf forest']) |
                 (zCEC10['Data']==lut_in['Sub-polar taiga needleleaf forest']) |
                 (zCEC10['Data']==lut_in['Tropical or sub-tropical broadleaf evergreen forest']) |
                 (zCEC10['Data']==lut_in['Tropical or sub-tropical broadleaf deciduous forest']) |
                 (zCEC10['Data']==lut_in['Temperate or sub-polar broadleaf deciduous forest']) |
                 (zCEC10['Data']==lut_in['Temperate or sub-polar needleleaf forest']) |
                 (zCEC10['Data']==lut_in['Mixed Forest']) )
    zCEC10c['Data'][ind]=lut_out['Forest']
    ind=np.where( (zCEC20['Data']==lut_in['Temperate or sub-polar needleleaf forest']) |
                 (zCEC20['Data']==lut_in['Sub-polar taiga needleleaf forest']) |
                 (zCEC20['Data']==lut_in['Tropical or sub-tropical broadleaf evergreen forest']) |
                 (zCEC20['Data']==lut_in['Tropical or sub-tropical broadleaf deciduous forest']) |
                 (zCEC20['Data']==lut_in['Temperate or sub-polar broadleaf deciduous forest']) |
                 (zCEC20['Data']==lut_in['Temperate or sub-polar needleleaf forest']) |
                 (zCEC20['Data']==lut_in['Mixed Forest']) )
    zCEC20c['Data'][ind]=lut_out['Forest']

    # Shrubland
    ind=np.where( (zCEC10['Data']==lut_in['Tropical or sub-tropical shrubland']) |
                 (zCEC10['Data']==lut_in['Temperate or sub-polar shrubland']) )
    zCEC10c['Data'][ind]=lut_out['Shrubland']
    ind=np.where( (zCEC20['Data']==lut_in['Tropical or sub-tropical shrubland']) |
                 (zCEC20['Data']==lut_in['Temperate or sub-polar shrubland']) )
    zCEC20c['Data'][ind]=lut_out['Shrubland']

    # Grassland
    ind=np.where( (zCEC10['Data']==lut_in['Tropical or sub-tropical grassland']) |
                 (zCEC10['Data']==lut_in['Temperate or sub-polar grassland']) )
    zCEC10c['Data'][ind]=lut_out['Grassland']
    ind=np.where( (zCEC20['Data']==lut_in['Tropical or sub-tropical grassland']) |
                 (zCEC20['Data']==lut_in['Temperate or sub-polar grassland']) )
    zCEC20c['Data'][ind]=lut_out['Grassland']

    # Liche-moss
    ind=np.where( (zCEC10['Data']==lut_in['Sub-polar or polar shrubland-lichen-moss']) |
                 (zCEC10['Data']==lut_in['Sub-polar or polar grassland-lichen-moss']) |
                 (zCEC10['Data']==lut_in['Sub-polar or polar barren-lichen-moss']) )
    zCEC10c['Data'][ind]=lut_out['Lichen-Moss']
    ind=np.where( (zCEC20['Data']==lut_in['Sub-polar or polar shrubland-lichen-moss']) |
                 (zCEC20['Data']==lut_in['Sub-polar or polar grassland-lichen-moss']) |
                 (zCEC20['Data']==lut_in['Sub-polar or polar barren-lichen-moss']) )
    zCEC20c['Data'][ind]=lut_out['Lichen-Moss']

    # Wetland
    ind=np.where( (zCEC10['Data']==lut_in['Wetland']) )
    zCEC10c['Data'][ind]=lut_out['Wetland']
    ind=np.where( (zCEC20['Data']==lut_in['Wetland']) )
    zCEC20c['Data'][ind]=lut_out['Wetland']
    cnt=cnt+1

    # Cropland
    ind=np.where( (zCEC10['Data']==lut_in['Cropland']) )
    zCEC10c['Data'][ind]=lut_out['Cropland']
    ind=np.where( (zCEC20['Data']==lut_in['Cropland']) )
    zCEC20c['Data'][ind]=lut_out['Cropland']

    # Barren Land
    ind=np.where( (zCEC10['Data']==lut_in['Barren lands']) )
    zCEC10c['Data'][ind]=lut_out['Barren Ground']
    ind=np.where( (zCEC20['Data']==lut_in['Barren lands']) )
    zCEC20c['Data'][ind]=lut_out['Barren Ground']

    # Urban
    ind=np.where( (zCEC10['Data']==lut_in['Urban']) )
    zCEC10c['Data'][ind]=lut_out['Urban']
    ind=np.where( (zCEC20['Data']==lut_in['Urban']) )
    zCEC20c['Data'][ind]=lut_out['Urban']

    # Water
    ind=np.where( (zCEC10['Data']==lut_in['Water']) )
    zCEC10c['Data'][ind]=lut_out['Water']
    ind=np.where( (zCEC20['Data']==lut_in['Water']) )
    zCEC20c['Data'][ind]=lut_out['Water']
    cnt=cnt+1

    # Snow and Ice
    ind=np.where( (zCEC10['Data']==lut_in['Snow and Ice']) )
    zCEC10c['Data'][ind]=lut_out['Snow and Ice']
    ind=np.where( (zCEC20['Data']==lut_in['Snow and Ice']) )
    zCEC20c['Data'][ind]=lut_out['Snow and Ice']
    cnt=cnt+1

    return lut_out,zCEC10c,zCEC20c

#%%
def PrepareLandUseCEC(meta):
    zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
    fin=r'C:\Data\Land Cover\NALCMS\NALCMS_lc_2010.tif'
    fout=r'C:\Data\Land Cover\NALCMS\NALCMS_lc_2010_bc1ha.tif'
    gis.ReprojectRasterAndClipToRaster(fin,fout,meta['Paths']['bc1ha Ref Grid'],meta['Geos']['crs'])
    zLCC_CEC10=gis.OpenGeoTiff(fout)
    zLCC_CEC10=gis.ClipToRaster(zLCC_CEC10,zRef)    
    ind=np.where(zRef['Data']==0)
    zLCC_CEC10['Data'][ind]=0 
    fin=r'C:\Data\Land Cover\NALCMS\NALCMS_lc_2020.tif'
    fout=r'C:\Data\Land Cover\NALCMS\NALCMS_lc_2020_bc1ha.tif'
    gis.ReprojectRasterAndClipToRaster(fin,fout,meta['Paths']['bc1ha Ref Grid'],meta['Geos']['crs'])    
    zLCC_CEC20=gis.OpenGeoTiff(fout)
    zLCC_CEC20=gis.ClipToRaster(zLCC_CEC20,zRef)    
    ind=np.where(zRef['Data']==0)
    zLCC_CEC20['Data'][ind]=0         
    # Compress categories to remove tropics    
    lut,zLCC_CEC10,zLCC_CEC20=NALCMS_Compress(meta['LUT']['Derived']['lcc_cec'],zRef,zLCC_CEC10,zLCC_CEC20)
    # Manually saved LUT to excell
    gis.SaveGeoTiff(zLCC_CEC10,meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCover_CEC_2010_Compressed.tif')
    gis.SaveGeoTiff(zLCC_CEC20,meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCover_CEC_2020_Compressed.tif')
    return

#%% Prepare NTEMS data
def PrepareNTEMS(meta):
    zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
    # Harvest year
    fin=r'C:\Data\Harvest\NTEMS\harv85to20.tif'
    a=gis.OpenGeoTiff(fin)
    fout=r'C:\Data\BC1ha\Disturbances\Harvest_NTEM_Year.tif'
    gis.ReprojectRasterAndClipToRaster(fin,fout,meta['Paths']['bc1ha Ref Grid'],meta['Geos']['crs'])
    z=gis.OpenGeoTiff(fout)
    z=gis.ClipToRaster(z,zRef)
    ind=np.where(zRef['Data']==0)
    z['Data'][ind]=0
    gis.SaveGeoTiff(z,fout)
    # Land cover 
    fin=r'C:\Data\Land Cover\NTEMS\vlce2_2019.tif'
    fout=meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCover_NTEMS_2019.tif'
    gis.ReprojectRasterAndClipToRaster(fin,fout,meta['Paths']['bc1ha Ref Grid'],meta['Geos']['crs'])
    z=gis.OpenGeoTiff(fout)
    z=gis.ClipToRaster(z,zRef)
    ind=np.where(zRef['Data']==0)
    z['Data'][ind]=0
    gis.SaveGeoTiff(z,fout)
    # Age
    fin=r'C:\Data\Age\NTEMS\age_ntem_c.tif'
    fout=meta['Paths']['bc1ha'] + '\\Age\\Age_NTEM_2019.tif'
    gis.ReprojectRasterAndClipToRaster(fin,fout,meta['Paths']['bc1ha Ref Grid'],meta['Geos']['crs'])
    return

#%% Reclassify NTEMS land cover into LCC Comp 1 classes

def ReclassifyNTEMS_LandCover(meta):

    zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
    vList=['lc_ntems_2019']
    z0=Import_Raster(meta,[],vList)['lc_ntems_2019']['Data'] 

    z1=copy.deepcopy(zRef)
    z1['Data']=(meta['LUT']['Derived']['lc_comp1']['Water']+1)*np.ones(zRef['Data'].shape,dtype='int8')

    ind=np.where( (np.isin(z0,[meta['LUT']['Derived']['lc_ntems_2019']['Coniferous'],meta['LUT']['Derived']['lc_ntems_2019']['Broadleaf'],meta['LUT']['Derived']['lc_ntems_2019']['Mixedwood'],meta['LUT']['Derived']['lc_ntems_2019']['Wetland-treed']])==True) )
    z1['Data'][ind]=meta['LUT']['Derived']['lc_comp1']['Forest']

    ind=np.where( (np.isin(z0,[meta['LUT']['Derived']['lc_ntems_2019']['Shrubs']])==True) )
    z1['Data'][ind]=meta['LUT']['Derived']['lc_comp1']['Shrub']

    ind=np.where( (np.isin(z0,[meta['LUT']['Derived']['lc_ntems_2019']['Herbs']])==True) )
    z1['Data'][ind]=meta['LUT']['Derived']['lc_comp1']['Herb']

    ind=np.where( (np.isin(z0,[meta['LUT']['Derived']['lc_ntems_2019']['Wetland']])==True) )
    z1['Data'][ind]=meta['LUT']['Derived']['lc_comp1']['Wetland']

    ind=np.where( (np.isin(z0,[meta['LUT']['Derived']['lc_ntems_2019']['Bryoids']])==True) )
    z1['Data'][ind]=meta['LUT']['Derived']['lc_comp1']['Bryoid']

    ind=np.where( (np.isin(z0,[meta['LUT']['Derived']['lc_ntems_2019']['Rock and rubble'],meta['LUT']['Derived']['lc_ntems_2019']['Exposed barren land']])==True) )
    z1['Data'][ind]=meta['LUT']['Derived']['lc_comp1']['Earth and Rock']

    ind=np.where( (np.isin(z0,[meta['LUT']['Derived']['lc_ntems_2019']['Snow and ice']])==True) )
    z1['Data'][ind]=meta['LUT']['Derived']['lc_comp1']['Snow and Ice']

    ind=np.where( (np.isin(z0,[meta['LUT']['Derived']['lc_ntems_2019']['Water']])==True) )
    z1['Data'][ind]=meta['LUT']['Derived']['lc_comp1']['Water']
    
    plt.close('all'); plt.matshow(z1['Data'])
    gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCover_NTEMS_2019_ReclassAsComp1.tif')
    return

#%% Reclassify VRI land cover into LCC Comp 1 classes

def ReclassifyVRI_LandCover(meta):
    vList=['refg','lc_vri_l2','lc_vri_l3','lc_vri_l4']
    z=Import_Raster(meta,[],vList)    
    z1=z['refg'].copy()
    z1['Data']=np.zeros(z['refg']['Data'].shape,dtype='int8')
    ind=np.where( (np.isin(z['lc_vri_l4']['Data'],[meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_4']['TC'],meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_4']['TM'],meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_4']['TB']])==True) )
    z1['Data'][ind]=meta['LUT']['Derived']['lc_comp1']['Forest']
    ind=np.where( (np.isin(z['lc_vri_l4']['Data'],[meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_4']['SI'],meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_4']['SL'],meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_4']['ST']])==True) )
    z1['Data'][ind]=meta['LUT']['Derived']['lc_comp1']['Shrub']
    ind=np.where( (np.isin(z['lc_vri_l4']['Data'],[meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_4']['HE'],meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_4']['HF'],meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_4']['HG']])==True) )
    z1['Data'][ind]=meta['LUT']['Derived']['lc_comp1']['Herb']
    ind=np.where( (np.isin(z['lc_vri_l4']['Data'],[meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_4']['BL'],meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_4']['BM'],meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_4']['BY']])==True) )
    z1['Data'][ind]=meta['LUT']['Derived']['lc_comp1']['Bryoid']    
    ind=np.where( (np.isin(z['lc_vri_l3']['Data'],[meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_3']['W']])==True) )
    z1['Data'][ind]=meta['LUT']['Derived']['lc_comp1']['Wetland']    
    ind=np.where( (np.isin(z['lc_vri_l4']['Data'],[meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_4']['EL'],meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_4']['RO']])==True) )
    z1['Data'][ind]=meta['LUT']['Derived']['lc_comp1']['Earth and Rock']
    ind=np.where( (np.isin(z['lc_vri_l4']['Data'],[meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_4']['SI']])==True) )
    z1['Data'][ind]=meta['LUT']['Derived']['lc_comp1']['Snow and Ice']
    ind=np.where( (np.isin(z['lc_vri_l2']['Data'],[meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_2']['W']])==True) )
    z1['Data'][ind]=meta['LUT']['Derived']['lc_comp1']['Water']
    ind=np.where( (z1['Data']==0) )
    z1['Data'][ind]=meta['LUT']['Derived']['lc_comp1']['Water']+1
    gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCover_VRI_ReclassAsComp1.tif')
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
    gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\Disturbances\\Harvest_PctDead_FromCruise.tif')
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
def DeriveBurnSeverityCompilation(meta):
    zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
        
    # Import provincial data and re-classify according to comp1
    zYb=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VEG_BURN_SEVERITY_SP\\FIRE_YEAR.tif')
    zBb=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VEG_BURN_SEVERITY_SP\\BURN_SEVERITY_RATING.tif')
    zBb2=copy.deepcopy(zBb)
    for k in meta['LUT']['Derived']['burnsev_comp1'].keys():
        ind=np.where( (zBb['Data']==meta['LUT']['VEG_BURN_SEVERITY_SP']['BURN_SEVERITY_RATING'][k]) )
        zBb2['Data'][ind]=meta['LUT']['Derived']['burnsev_comp1'][k]
    n,p=gu.CountByCategories(zBb2['Data'][zYb['Data']>0],'Percent')
    plt.bar(n,p)
    
    # Import national data, scale and classify according to (https://catalogue.data.gov.bc.ca/dataset/fire-burn-severity-same-year)
    zYn=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Disturbances\\BurnSev_Can_Year.tif')
    zBn=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Disturbances\\BurnSev_Can_dNBR.tif')
    #plt.hist(zBn['Data'][zBn['Data']>0].flatten()[0::50])
    zBn2=copy.deepcopy(zBn)
    zBn2['Data']=(zBn2['Data'].astype('float')+275)/5
    zBn2['Data']=zBn2['Data'].astype('int16')
    #plt.hist(zBn2['Data'][zBn2['Data']>0].flatten()[0::50])
    zBn3=copy.deepcopy(zRef)
    zBn3['Data']=np.zeros(zRef['Data'].shape,'int8')
    ind=np.where( (zBn2['Data']<76) ); zBn3['Data'][ind]=meta['LUT']['Derived']['burnsev_comp1']['Unburned']
    ind=np.where( (zBn2['Data']>=76) & (zBn2['Data']<110) ); zBn3['Data'][ind]=meta['LUT']['Derived']['burnsev_comp1']['Low']
    ind=np.where( (zBn2['Data']>=110) & (zBn2['Data']<187) ); zBn3['Data'][ind]=meta['LUT']['Derived']['burnsev_comp1']['Medium']
    ind=np.where( (zBn2['Data']>=187) ); zBn3['Data'][ind]=meta['LUT']['Derived']['burnsev_comp1']['High']
    
    #n,p=gu.CountByCategories(zBn3['Data'][zYn['Data']>0],'Percent')
    #plt.bar(n,p)
    
    # # Use distribution of BC classes to derive classes from dNBR for Canada (exclude unburned)
    # p=gu.CountByCategories(zBb['Data'][zBb['Data']>0],'Percent')
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
        zBS=copy.deepcopy(zRef)
        zBS['Data']=5*np.ones(zRef['Data'].shape,dtype='int8')
        
        # Add national data (converted from dNBR to severity class)
        uY=np.unique(zYn['Data'][zYn['Data']>0])
        for iU in range(uY.size):    
            ind=np.where( (zY['Data']==uY[iU]) & (zYn['Data']==uY[iU]) )
            zBS['Data'][ind]=zBn3['Data'][ind]
            zLastY['Data'][ind]=uY[iU]
            zLastS['Data'][ind]=zBn3['Data'][ind]
            
            # for iC in range(cpn.size):
            #     ind=np.where( (zY['Data']==uY[iU]) & (zYn['Data']==uY[iU]) & (zBn['Data']<cpn[iC]) )
            #     zBS['Data'][ind]=list(meta['LUT']['Derived']['burnsev_comp1'].values())[iC]
            #     zLastY['Data'][ind]=uY[iU]
            #     zLastS['Data'][ind]=list(meta['LUT']['Derived']['burnsev_comp1'].values())[iC]
        
        # Add BC data
        uY=np.unique(zYb['Data'][zYb['Data']>0])
        for iU in range(uY.size):
            ind=np.where( (zY['Data']==uY[iU]) & (zYb['Data']==uY[iU]) )
            zBS['Data'][ind]=zBb2['Data'][ind]
            zLastY['Data'][ind]=uY[iU]
            zLastS['Data'][ind]=zBb2['Data'][ind]
        #plt.matshow(zBS['Data'],clim=[0,4])
        # Save
        gis.SaveGeoTiff(zBS,meta['Paths']['bc1ha'] + '\\PROT_HISTORICAL_FIRE_POLYS_SP\\FIRE_YEAR_' + str(iY+1) + '_SevClass.tif')
    
    gis.SaveGeoTiff(zLastY,meta['Paths']['bc1ha'] + '\\Disturbances\\BurnSevComp_YearLast.tif')
    gis.SaveGeoTiff(zLastS,meta['Paths']['bc1ha'] + '\\Disturbances\\BurnSevComp_SevClassLast.tif')
    return
    
#%% 
def DeriveHarvestRetentionComp1(meta):
    zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
    lut1=meta['LUT']['Derived']['harvret1']
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
    gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\Management\\HarvestRetentionComp1.tif')
    return

#%% Harvest regeneration type
def DeriveHarvestRegenType(meta):
    zH=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Disturbances\\VEG_CONSOLIDATED_CUT_BLOCKS_SP_MaskAll.tif')
    zFCR=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Results\\RSLT_FOREST_COVER_RESERVE_SVW.tif')
    zPL=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Results\\Planting_FromRESULTS_MaskCount.tif')
    
    tv=np.arange(1960,2022,1)
    HY=0*zRef['Data'].astype('int16')
    for iT in range(tv.size):
        zH0=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Disturbances\\VEG_CONSOLIDATED_CUT_BLOCKS_SP_' + str(tv[iT]) + '.tif')
        ind=np.where( (zH0['Data']>0) & (HY==0) )
        HY[ind]=tv[iT]
    
    z=copy.deepcopy(zRef)
    z['Data']=np.zeros(zRef['Data'].shape,dtype='int8')
    # Harvested and planted
    ind=np.where( (zH['Data']>0) & (zFCR['Data']==0) & (zPL['Data']>0) )
    z['Data'][ind]=1
    # Harvested but not planted (before FRPA)
    ind=np.where( (zH['Data']>0) & (HY<=1987) & (zFCR['Data']==0) & (zPL['Data']==0) )
    z['Data'][ind]=2
    # Harvested but not planted (after FRPA)
    ind=np.where( (zH['Data']>0) & (HY>1987) & (HY<=2018) & (zFCR['Data']==0) & (zPL['Data']==0) )
    z['Data'][ind]=3
    # Harvested but not planted (after FRPA)
    ind=np.where( (zH['Data']>0) & (HY>2018) & (zFCR['Data']==0) & (zPL['Data']==0) )
    z['Data'][ind]=4
    # Not harvested, but planted
    ind=np.where( (zH['Data']==0) & (zPL['Data']>0) | (zFCR['Data']>0) & (zPL['Data']>0) )
    z['Data'][ind]=5
    # Not harvested, not planted
    ind=np.where( (zH['Data']==0) & (zFCR['Data']==0) & (zPL['Data']==0) )
    z['Data'][ind]=6
    np.unique(z['Data'])
    
    bin=np.arange(1,6,1)
    n=np.zeros(bin.size)
    for i in range(bin.size):
        ind=np.where(z['Data']==bin[i])
        n[i]=ind[0].size
    plt.bar(bin,n/np.sum(n)*100)
    
    # plt.matshow(z['Data'])
    gis.SaveGeoTiff(z,meta['Paths']['bc1ha'] + '\\Results\\Harvest_Regen_Type.tif')
    
    # # Just current era = 79%
    # n1=np.where( (HY>1987) & (HY<2018) & (zFCR['Data']==0) & (zPL['Data']>0) )[0].size
    # n2=np.where( (HY>1987) & (HY<2018) & (zFCR['Data']==0) )[0].size
    # print(n1/n2)
    
    # Time series
    yH=np.zeros(tv.size)
    yHP=np.zeros(tv.size)
    for iT in range(tv.size):
        zH=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Disturbances\\VEG_CONSOLIDATED_CUT_BLOCKS_SP_' + str(tv[iT]) + '.tif')
        ind=np.where( (zH['Data']>0) & (zFCR['Data']==0) )
        yH[iT]=ind[0].size
        ind=np.where( (zH['Data']>0) & (zFCR['Data']==0) & (zPL['Data']>0) )
        yHP[iT]=ind[0].size
    
    plt.plot(tv,(yH-yHP)/yH,'ob-')
    return

#%% Consolidate Tree Density Class
# *** Ensure harvested areas are listed as dense ***
def DeriveTreeDensityClass(meta):
    # 2002 VRI
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
def PlantingLayer(meta,zRef):
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

#%% Salvage logging mask
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

#%% Derive Regen Type (used for non-ob stand establishment)
def DeriveRegenTypeCompilation(meta):
    meta=cbu.Load_LUTs_Modelling(meta)
    
    d={}
    d['tv']=np.arange(1960,2023,1)
    ptNam=np.array(list(meta['LUT']['Derived']['RegenType'].keys()))
    ptID=np.array(list(meta['LUT']['Derived']['RegenType'].values()))
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
    #     for iP in range(6):
    #         ind=np.where( (zP[iP]['STC']==meta['LUT']['RSLT_ACTIVITY_TREATMENT_SVW']['SILV_TECHNIQUE_CODE'][v]) )
    #         d['STC Summary']['N'][cnt,0]=d['STC Summary']['N'][cnt,0]+ind[0].size
    #         ind=np.where( (zP[iP]['STC']==meta['LUT']['RSLT_ACTIVITY_TREATMENT_SVW']['SILV_TECHNIQUE_CODE'][v]) & (np.isin(zP[iP]['FSC'],meta['Param']['BE']['FSC']['NO List ID'])==True) )
    #         d['STC Summary']['N'][cnt,1]=d['STC Summary']['N'][cnt,1]+ind[0].size
    
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
        d['A NOSE'][iT,meta['LUT']['Derived']['RegenType']['Replanting']-1]=ind[0].size
        ind=np.where( (pl_oc==1) & (pl_stc==meta['LUT']['RSLT_ACTIVITY_TREATMENT_SVW']['SILV_TECHNIQUE_CODE']['RP']) )
        d['A'][iT,meta['LUT']['Derived']['RegenType']['Replanting']-1]=ind[0].size
        zType0[ind]=meta['LUT']['Derived']['RegenType']['Replanting']
        
        # Fill Planting        
        ind=np.where( (pl_oc==1) & (pl_stc==meta['LUT']['RSLT_ACTIVITY_TREATMENT_SVW']['SILV_TECHNIQUE_CODE']['FP']) & (np.isin(pl_fsc,meta['Param']['BE']['FSC']['NO List ID'])==True) )
        d['A NOSE'][iT,meta['LUT']['Derived']['RegenType']['Fill Planting']-1]=ind[0].size
        ind=np.where( (pl_oc==1) & (pl_stc==meta['LUT']['RSLT_ACTIVITY_TREATMENT_SVW']['SILV_TECHNIQUE_CODE']['FP']) )        
        d['A'][iT,meta['LUT']['Derived']['RegenType']['Fill Planting']-1]=ind[0].size
        zType0[ind]=meta['LUT']['Derived']['RegenType']['Fill Planting']
        
        # Road rehab
        ind=np.where( (pl_oc==1) & (pl_stc==meta['LUT']['RSLT_ACTIVITY_TREATMENT_SVW']['SILV_TECHNIQUE_CODE']['RR']) & (np.isin(pl_fsc,meta['Param']['BE']['FSC']['NO List ID'])==True) )
        d['A NOSE'][iT,meta['LUT']['Derived']['RegenType']['Road Rehabilitation']-1]=ind[0].size
        ind=np.where( (pl_oc==1) & (pl_stc==meta['LUT']['RSLT_ACTIVITY_TREATMENT_SVW']['SILV_TECHNIQUE_CODE']['RR']) )
        d['A'][iT,meta['LUT']['Derived']['RegenType']['Road Rehabilitation']-1]=ind[0].size
        zType0[ind]=meta['LUT']['Derived']['RegenType']['Road Rehabilitation']
        
        # Back to back planting
        ind=np.where( (pl_oc==1) & (zType0==0) & (DL==meta['LUT']['Event']['Planting']) & (np.isin(pl_fsc,meta['Param']['BE']['FSC']['NO List ID'])==True) )
        d['A NOSE'][iT,meta['LUT']['Derived']['RegenType']['Back-to-back Planting']-1]=ind[0].size    
        ind=np.where( (pl_oc==1) & (zType0==0) & (DL==meta['LUT']['Event']['Planting']) )        
        d['A'][iT,meta['LUT']['Derived']['RegenType']['Back-to-back Planting']-1]=ind[0].size    
        zType0[ind]=meta['LUT']['Derived']['RegenType']['Back-to-back Planting']
        
        # Salvage
        ind=np.where( (pl_oc==1) & (d['tv'][iT]>=2004) & (zD['harv_yr_con1']['Data']>1987) & (zType0==0) & (DL==meta['LUT']['Event']['Harvest']) & (np.isin(pl_fsc,meta['Param']['BE']['FSC']['NO List ID'])==True) )
        d['A NOSE'][iT,meta['LUT']['Derived']['RegenType']['Salvage and Planting']-1]=ind[0].size
        ind=np.where( (pl_oc==1) & (d['tv'][iT]>=2004) & (zD['harv_yr_con1']['Data']>1987) & (zType0==0) & (DL==meta['LUT']['Event']['Harvest']) & (np.isin(pl_fsc,meta['Param']['BE']['FSC']['NO List ID'])==True) | \
                     (pl_oc==1) & (zType0==0) & (DL==meta['LUT']['Event']['Harvest']) & (zD['pdead_cruise']['Data']>=50) )
        d['A'][iT,meta['LUT']['Derived']['RegenType']['Salvage and Planting']-1]=ind[0].size        
        zType0[ind]=meta['LUT']['Derived']['RegenType']['Salvage and Planting']
        
        # Knockdown
        ind=np.where( (pl_oc==1) & (zType0==0) & (DL==meta['LUT']['Event']['Knockdown']) & (np.isin(pl_fsc,meta['Param']['BE']['FSC']['NO List ID'])==True) )
        d['A NOSE'][iT,meta['LUT']['Derived']['RegenType']['Knockdown and Planting']-1]=ind[0].size
        ind=np.where( (pl_oc==1) & (zType0==0) & (DL==meta['LUT']['Event']['Knockdown']) )
        d['A'][iT,meta['LUT']['Derived']['RegenType']['Knockdown and Planting']-1]=ind[0].size        
        zType0[ind]=meta['LUT']['Derived']['RegenType']['Knockdown and Planting']
        
        # Straight fire
        ind=np.where( (pl_oc==1) & (zType0==0) & (DL==meta['LUT']['Event']['Wildfire']) & (np.isin(pl_fsc,meta['Param']['BE']['FSC']['NO List ID'])==True) )
        d['A NOSE'][iT,meta['LUT']['Derived']['RegenType']['Straight-to-planting Post Wildfire']-1]=ind[0].size
        ind=np.where( (pl_oc==1) & (zType0==0) & (DL==meta['LUT']['Event']['Wildfire']) )
        d['A'][iT,meta['LUT']['Derived']['RegenType']['Straight-to-planting Post Wildfire']-1]=ind[0].size        
        zType0[ind]=meta['LUT']['Derived']['RegenType']['Straight-to-planting Post Wildfire']        
        
        # Straight insect
        ind=np.where( (pl_oc==1) & (zType0==0) & (DL==meta['LUT']['Event']['IBM']) & (np.isin(pl_fsc,meta['Param']['BE']['FSC']['NO List ID'])==True) )
        d['A NOSE'][iT,meta['LUT']['Derived']['RegenType']['Straight-to-planting Post Beetles']-1]=ind[0].size
        ind=np.where( (pl_oc==1) & (zType0==0) & (DL==meta['LUT']['Event']['IBM']) )
        d['A'][iT,meta['LUT']['Derived']['RegenType']['Straight-to-planting Post Beetles']-1]=ind[0].size
        zType0[ind]=meta['LUT']['Derived']['RegenType']['Straight-to-planting Post Beetles']
        
        # Straight other
        
        # Harvest and Planting NSR backlog
        ind=np.where( (pl_oc==1) & (zType0==0) & (DL==meta['LUT']['Event']['Harvest']) & (zD['harv_yr_con1']['Data']<=1987) & (np.isin(pl_fsc,meta['Param']['BE']['FSC']['NO List ID'])==True) )
        d['A NOSE'][iT,meta['LUT']['Derived']['RegenType']['Harvest and Planting NSR Backlog']-1]=ind[0].size        
        d['A'][iT,meta['LUT']['Derived']['RegenType']['Harvest and Planting NSR Backlog']-1]=ind[0].size        
        zType0[ind]=meta['LUT']['Derived']['RegenType']['Harvest and Planting NSR Backlog']
        
        # Harvest and Planting
        ind=np.where( (pl_oc==1) & (zType0==0) & (DL==meta['LUT']['Event']['Harvest']) & (np.isin(pl_fsc,meta['Param']['BE']['FSC']['NO List ID'])==True) )
        d['A NOSE'][iT,meta['LUT']['Derived']['RegenType']['Harvest and Planting']-1]=ind[0].size 
        ind=np.where( (pl_oc==1) & (zType0==0) & (DL==meta['LUT']['Event']['Harvest']) )
        d['A'][iT,meta['LUT']['Derived']['RegenType']['Harvest and Planting']-1]=ind[0].size         
        zType0[ind]=meta['LUT']['Derived']['RegenType']['Harvest and Planting']
        
        # Unknown
        ind=np.where( (pl_oc==1) & (zType0==0) & (np.isin(pl_fsc,meta['Param']['BE']['FSC']['NO List ID'])==True) )
        d['A NOSE'][iT,meta['LUT']['Derived']['RegenType']['Unknown']-1]=ind[0].size
        ind=np.where( (pl_oc==1) & (zType0==0) )
        d['A'][iT,meta['LUT']['Derived']['RegenType']['Unknown']-1]=ind[0].size        
        zType0[ind]=meta['LUT']['Derived']['RegenType']['Unknown']
        
        # Update last disturbance if planting occurs
        DL[(pl_oc==1)]=meta['LUT']['Event']['Planting']
        
        # Direct seeding
        for iDS in range(3):
            ind=np.where( (zDS[iDS]['Year']==d['tv'][iT]) & (np.isin(pl_fsc,meta['Param']['BE']['FSC']['NO List ID'])==True) )
            d['A NOSE'][iT,meta['LUT']['Derived']['RegenType']['Direct Seeding']-1]=ind[0].size
            ind=np.where( (zDS[iDS]['Year']==d['tv'][iT]) )            
            d['A'][iT,meta['LUT']['Derived']['RegenType']['Direct Seeding']-1]=ind[0].size            
            zType0[ind]=meta['LUT']['Derived']['RegenType']['Direct Seeding']
        
        # Pack
        for iP in range(6):
            ind=np.where( (zP[iP]['Year']==d['tv'][iT]) )
            zP[iP]['Type'][ind]=zType0[ind]
    
    # Save summary
    gu.opickle(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\RegenTypeSummary_All',d)
    
    # Save packed regen type
    for i in range(6):
        z1=copy.deepcopy(zRef)
        z1['Data']=zP[i]['Type']
        gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(i+1) + '_RegenType.tif')
    
    return

#%%
# def DeriveConservationCompolation(meta):
#     vList=['refg','lc_comp1_2019','protected','park','ogma','ogdef','parknat']
#     z=Import_Raster(meta,[],vList)
    
#     #zLCC1=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCoverClass1.tif')
#     #zProt=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\TA_PROTECTED_LANDS_SV\\PROTECTED_LANDS_DESIGNATION.tif')
#     #zPark=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\TA_PARK_ECORES_PA_SVW\\ADMIN_AREA_SID.tif')
#     #zOGMA=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RMP_OGMA_LEGAL_ALL_SVW\\LEGAL_OGMA_PROVID.tif')
#     #zOGD=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\OGSR_TAP_PRIORITY_DEF_AREA_SP\\PRIORITY_DEFERRAL_ID.tif')
#     #zNP=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\CLAB_NATIONAL_PARKS\\ENGLISH_NAME.tif')
#     #zFCR=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\RSLT_FOREST_COVER_RESERVE_SVW.tif')
    
#     # ListCon=['No Timber Harvesting Areas','Landscape Corridors',
#     #   'Critical Deer Winter Range','Sensitive Watershed','Water Management Units','High Value Wetlands for Moose','Telkwa Caribou Recovery Area',
#     #   'Caribou Migration Corridor','High Biodiversity Emphasis Areas','Scenic Corridors']
    
#     # Generate random new protected areas
#     flg=0
#     if flg==1:
#         gdf=u1ha.Import_GDBs_ProvinceWide()
#         N=2500 # Number within bounding box
#         Dbuf=6200 # metres
    
#         x=np.random.uniform(zRef['xmin'],zRef['xmax'],size=N)
#         y=np.random.uniform(zRef['ymin'],zRef['ymax'],size=N)
#         points=[]
#         for k in range(x.size):
#             points.append(Point(x[k],y[k]))
#         gdf_xy=gpd.GeoDataFrame({'geometry':points,'ID':1})
#         gdf_xy.crs=gdf['bc_bound']['gdf'].crs
    
#         gdf_xy=gpd.sjoin(gdf_xy,gdf['bc_bound']['gdf'],op='within')
#         #gdf_xy.plot(ax=ax[0],markersize=8)
    
#         gdf_xyb=gdf_xy.geometry.buffer(Dbuf)
#         gdf_xyb=gpd.GeoDataFrame({'geometry':gdf_xyb,'ID':np.arange(0,gdf_xyb.size,1)})
#         gdf_xyb.crs=gdf['bc_bound']['gdf'].crs
#         #gdf_xyb.plot(ax=ax[0],color='r')
    
#         shapes=((geom,value) for geom, value in zip(gdf_xyb['geometry'],gdf_xyb['ID']))
#         z=np.zeros(zRef['Data'].shape,dtype=float)
#         burned=features.rasterize(shapes=shapes,fill=0,out=z,transform=zRef['Transform'])
    
#         zRn=copy.deepcopy(zRef)
#         zRn['Data']=z.astype('int32')
#         gis.SaveGeoTiff(zRn,meta['Paths']['bc1ha'] + '\\LandCoverUse\\REARs_RandomAreas.tif')
    
#     else:
#         zRn=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\REARs_RandomAreas.tif')
    
#     # # Protected areas with random areas
    
#     # Area treed
#     ind0=np.where( (zLCC1['Data']==meta['LUT']['Derived']['lc_comp1']['Forest']) );
#     A_treed=ind0[0].size
    
#     # Everything completed+proposed without random additions (Comp=1, Prop=2)
#     zCP=copy.deepcopy(zRef)
#     zCP['Data']=np.zeros(zRef['Data'].shape,dtype='int16')
#     ind=np.where( (zLCC1['Data']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (zPark['Data']>0) |
#                   (zLCC1['Data']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (zOGMA['Data']>0) |
#                   (zLCC1['Data']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (zNP['Data']>0) |
#                   (zLCC1['Data']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (zProt['Data']>0))
#     zCP['Data'][ind]=1
#     ind=np.where( (zCP['Data']==0) & (zLCC1['Data']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (zOGD['Data']>0) )
#     zCP['Data'][ind]=2
#     ind=np.where( (zLCC1['Data']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (zCP['Data']==0) ); zCP['Data'][ind]=3
#     ind=np.where( (zLCC1['Data']!=meta['LUT']['Derived']['lc_comp1']['Forest']) & (zCP['Data']==0) ); zCP['Data'][ind]=4
#     ind=np.where( (zRef['Data']==0) ); zCP['Data'][ind]=5
    
#     ind=np.where( (zLCC1['Data']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (zCP['Data']>=1) & (zCP['Data']<=2) ); A_CPsep=ind[0].size
#     print('Area Protected from completed+proposed (Sep) (%): ' + str(A_CPsep/A_treed*100))
    
#     zCP['Data']=zCP['Data'].astype('int16')
#     gis.SaveGeoTiff(zCP,meta['Paths']['bc1ha'] + '\\LandCoverUse\\REARs_CompPlusProp.tif')
    
#     # Everything completed+proposed with random additions (Comp=1, Prop=2)
#     zCP=copy.deepcopy(zRef)
#     zCP['Data']=np.zeros(zRef['Data'].shape,dtype='int16')
#     ind=np.where( (zLCC1['Data']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (zPark['Data']>0) |
#                   (zLCC1['Data']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (zOGMA['Data']>0) |
#                   (zLCC1['Data']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (zNP['Data']>0) |
#                   (zLCC1['Data']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (zProt['Data']>0) )
#     zCP['Data'][ind]=1
#     ind=np.where( (zCP['Data']==0) & (zLCC1['Data']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (zOGD['Data']>0) | (zCP['Data']==0) & (zLCC1['Data']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (zRn['Data']>0) )
#     zCP['Data'][ind]=2
#     ind=np.where( (zLCC1['Data']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (zCP['Data']==0) ); zCP['Data'][ind]=3
#     ind=np.where( (zLCC1['Data']!=meta['LUT']['Derived']['lc_comp1']['Forest']) & (zCP['Data']==0) ); zCP['Data'][ind]=4
#     ind=np.where( (zRef['Data']==0) ); zCP['Data'][ind]=5
    
#     ind=np.where( (zLCC1['Data']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (zCP['Data']>=1) & (zCP['Data']<=2) ); A_CPsep=ind[0].size
#     print('Area Protected from completed+proposed (Sep) (%): ' + str(A_CPsep/A_treed*100))
    
#     zCP['Data']=zCP['Data'].astype('int32')
#     gis.SaveGeoTiff(zCP,meta['Paths']['bc1ha'] + '\\LandCoverUse\\REARs_CompPlusPropWithRandomAdditions.tif')
#     return

#%%
def DeriveLandCoverComp1(meta):    
    
    zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
    vList=['lc_vri_l3','lc_vri_l4','lc_cec_2020','lc_ntems_2019','harv_yr_cc','harv_yr_ntems']
    z=Import_Raster(meta,[],vList,'Extract Grid')
    lut_comp1=meta['LUT']['Derived']['lc_comp1']
    lut_ntems=meta['LUT']['Derived']['lc_ntems_2019']    
    lut_lc_vri_l3=meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_3']
    lut_lc_vri_l4=meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_4']
    #--------------------------------------------------------------------------
    # 2019 Land Cover
    # In the past, lc_vri_l2 treed area has differed from that of lc_vri_l4 (treatment of TFLs and private land), but that appears to be fixed in 2023
    #--------------------------------------------------------------------------    
    z1=copy.deepcopy(zRef)
    z1['Data']=np.zeros(zRef['Data'].shape,dtype='int8')
    ind=np.where( (np.isin(z['lc_ntems_2019'],[lut_ntems['Coniferous'],lut_ntems['Broadleaf'],lut_ntems['Mixedwood'],lut_ntems['Wetland-treed']])==True) )
    z1['Data'][ind]=lut_comp1['Forest']
    ind=np.where( (np.isin(z['lc_ntems_2019'],[lut_ntems['Shrubs']])==True) )
    z1['Data'][ind]=lut_comp1['Shrub']
    ind=np.where( (np.isin(z['lc_ntems_2019'],[lut_ntems['Herbs']])==True) )
    z1['Data'][ind]=lut_comp1['Herb']    
    ind=np.where( (np.isin(z['lc_ntems_2019'],[lut_ntems['Wetland']])==True) )
    z1['Data'][ind]=lut_comp1['Wetland']
    ind=np.where( (np.isin(z['lc_ntems_2019'],[lut_ntems['Bryoids']])==True) )
    z1['Data'][ind]=lut_comp1['Bryoid']
    ind=np.where( (np.isin(z['lc_ntems_2019'],[lut_ntems['Rock and rubble'],lut_ntems['Exposed barren land']])==True) )
    z1['Data'][ind]=lut_comp1['Earth and Rock']
    ind=np.where( (np.isin(z['lc_ntems_2019'],[lut_ntems['Snow and ice']])==True) )
    z1['Data'][ind]=lut_comp1['Snow and Ice']
    ind=np.where( (np.isin(z['lc_ntems_2019'],[lut_ntems['Water']])==True) )
    z1['Data'][ind]=lut_comp1['Water']
    
    # If NTEM says harvest and CEC says not settlement, assume forest
    ind=np.where( (z['harv_yr_ntems']>0) & (z['lc_cec_2020']!=meta['LUT']['Derived']['lc_cec_c']['Urban']) & (z['lc_cec_2020']!=meta['LUT']['Derived']['lc_cec_c']['Cropland']) )
    z1['Data'][ind]=lut_comp1['Forest']
    
    # Override with forest from VRI lc_vri_l4
    # *** Stopped doing this because it leads to unrealistically high forest area ***
    #ind=np.where( (np.isin(z['lc_vri_l4'],[lut_lc_vri_l4['TC'],lut_lc_vri_l4['TM'],lut_lc_vri_l4['TB']])==True) )
    #z1['Data'][ind]=lut_comp1['Forest']
    
    # Convert areas  where CEC map says urban
    ind=np.where( (z['lc_cec_2020']==meta['LUT']['Derived']['lc_cec_c']['Urban']) )
    z1['Data'][ind]=lut_comp1['Built-up']
    
    # CEC is more reliable indicator of cropland -> convert to shrubland
    ind=np.where( (z1['Data']==lut_comp1['Herb']) & (z['lc_cec_2020']!=meta['LUT']['Derived']['lc_cec_c']['Cropland']) )
    z1['Data'][ind]=lut_comp1['Shrub']
    
    # Check for any cells left unclassified
    # Very small number, likely water, reclassify as water
    ind=np.where( (zRef['Data']==1) & (z1['Data']==0) )
    print(ind[0].size)
    z1['Data'][ind]=lut_comp1['Water']

    # Check forest area
    ind=np.where( (z1['Data']==lut_comp1['Forest']) )
    A['For 2019']=ind[0].size/1e6

    # Patches get left as earth and rock, possibly associated with wildfire and
    # fire guards. Fix by converting to forest outside alpine areas
    zE=Import_Raster(meta,[],['elev'],'Extract Grid')
    
    MaskBuf=np.zeros(zRef['Data'].shape)
    ind=np.where(z['lc_vri_l3']==lut_lc_vri_l3['A'])
    MaskBuf[ind]=1
    MaskBuf=gis.BufferRasterMask(MaskBuf,5)
    ind=np.where( (z1['Data']==lut_comp1['Earth and Rock']) & (MaskBuf==0) & (zE['elev']<1450) )
    A_ToAdd=ind[0].size/1e6
    print(A_ToAdd)
    # a=np.zeros(zRef['Data'].shape); a[ind]=1
    # plt.matshow(a)
    z1['Data'][ind]=lut_comp1['Forest']    

    # Reclassify outside land mask as water
    ind=np.where( (zRef['Data']==0) )
    z1['Data'][ind]=lut_comp1['Water']
    
    # Save
    gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCover_Comp1_2019.tif')

    # Look at result
    flg=0
    if flg==1:
        plt.close('all')
        plt.matshow(z1['Data'])
        
        a=zRef['Data'].copy()
        ind=np.where(z1['Data']==lut_comp1['Herb'])
        a[ind]=2
        plt.matshow(a)    
    
        a=zRef['Data'].copy()
        ind=np.where(z['lc_cec_2020']['Data']==meta['LUT']['Derived']['lc_cec_c']['Cropland'])
        a[ind]=2
        plt.matshow(a)
    
    return

#%%
def DeriveLandUseComp1(meta):
    zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
    
    lut_comp1=meta['LUT']['Derived']['lc_comp1']
    lut_lu=meta['LUT']['Derived']['lu_comp1']
    lut_ntems=meta['LUT']['Derived']['lc_ntems_2019']
    lut_lc_vri_l3=meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_3']
    lut_cec=meta['LUT']['Derived']['lc_cec_c']
    
    vList=['lc_comp1_2019','lc_vri_l3','munic'] 
    z=Import_Raster(meta,[],vList,'Extract Grid')
    
    z1=copy.deepcopy(zRef)
    z1['Data']=np.zeros(zRef['Data'].shape,dtype='int8')
    
    #--------------------------------------------------------------------------
    # Conservation
    #--------------------------------------------------------------------------
    vList=['protected','park','ogma','ogdef','parknat','harv_yr_comp2','own']
    zC=Import_Raster(meta,[],vList,'Extract Grid')
    ind=np.where( (zC['harv_yr_comp2']==0) & (zC['protected']>0) | 
                 (zC['harv_yr_comp2']==0) & (zC['park']>0) |
                 (zC['harv_yr_comp2']==0) & (zC['parknat']>0) |
                 (zC['harv_yr_comp2']==0) & (zC['ogma']>0) |
                 (zC['harv_yr_comp2']==0) & (zC['own']==meta['LUT']['F_OWN']['OWNERSHIP_DESCRIPTION']['Crown - Conservancy Area, Ecological Reserve, Protected Area, Provincial Park']) |
                 (zC['harv_yr_comp2']==0) & (zC['own']==meta['LUT']['F_OWN']['OWNERSHIP_DESCRIPTION']['Crown - Local/Regional Park']) )
    z1['Data'][ind]=lut_lu['Conservation Natural']
    
    ind=np.where( (zC['harv_yr_comp2']>0) & (zC['protected']>0) | 
                 (zC['harv_yr_comp2']>0) & (zC['park']>0) |
                 (zC['harv_yr_comp2']>0) & (zC['parknat']>0) |
                 (zC['harv_yr_comp2']>0) & (zC['ogma']>0) |
                 (zC['harv_yr_comp2']>0) & (zC['own']==meta['LUT']['F_OWN']['OWNERSHIP_DESCRIPTION']['Crown - Conservancy Area, Ecological Reserve, Protected Area, Provincial Park']) |
                 (zC['harv_yr_comp2']>0) & (zC['own']==meta['LUT']['F_OWN']['OWNERSHIP_DESCRIPTION']['Crown - Local/Regional Park']) )
    z1['Data'][ind]=lut_lu['Conservation Consistent']
    del zC
    
    #--------------------------------------------------------------------------
    # Cropland
    #--------------------------------------------------------------------------
    gdf=gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='OATS_ALR_POLYS')
    gdf['geometry']=gdf.geometry.buffer(10000)
    gdf=gdf[gdf.geometry!=None]
    gdf=gdf.reset_index()
    gdf['ID']=np.ones(len(gdf))
    shapes=((geom,value) for geom, value in zip(gdf['geometry'],gdf['ID']))
    zBuf=np.zeros(zRef['Data'].shape,dtype=float)
    burned=features.rasterize(shapes=shapes,fill=0,out=zBuf,transform=zRef['Transform'])    
    ind=np.where( (np.isin(z['lc_comp1_2019'],[lut_comp1['Herb']])==True) & (zBuf>0) )
    z1['Data'][ind]=lut_lu['Cropland']
    
    gdf=gpd.read_file(meta['Paths']['GDB']['GDB'] + '\\crop_per.shp')
    #gdf['geometry']=gdf.geometry.buffer(10000)
    gdf=gdf[gdf.geometry!=None]
    gdf=gdf.reset_index()
    gdf['ID']=np.ones(len(gdf))
    shapes=((geom,value) for geom, value in zip(gdf['geometry'],gdf['ID']))
    zBuf=np.zeros(zRef['Data'].shape,dtype=float)
    burned=features.rasterize(shapes=shapes,fill=0,out=zBuf,transform=zRef['Transform'])    
    ind=np.where( (np.isin(z['lc_comp1_2019'],[lut_comp1['Herb'],lut_comp1['Shrub']])==True) & (zBuf>0) )
    z1['Data'][ind]=lut_lu['Cropland']
    
    #--------------------------------------------------------------------------
    # Pasture
    #--------------------------------------------------------------------------
    gdf=gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='RANGE_PASTURE_POLY_SVW')
    gdf=gdf[(gdf['geometry'].area<1e8)]
    gdf['geometry']=gdf.geometry.buffer(5000)
    gdf=gdf[gdf.geometry!=None]
    gdf=gdf.reset_index()
    gdf['ID']=np.ones(len(gdf))
    shapes=((geom,value) for geom, value in zip(gdf['geometry'],gdf['ID']))
    zBuf=np.zeros(zRef['Data'].shape,dtype=float)
    burned=features.rasterize(shapes=shapes,fill=0,out=zBuf,transform=zRef['Transform'])    
    ind=np.where( (z1['Data']==0) & (np.isin(z['lc_comp1_2019'],[lut_comp1['Shrub']])==True) & (zBuf>0) )
    z1['Data'][ind]=lut_lu['Pasture']
    
    #--------------------------------------------------------------------------
    # Residential and commercial
    #--------------------------------------------------------------------------
    zRC=np.zeros(zRef['Data'].shape,dtype='int8') # Keep track of all residential and commercial search areas
    z['munic'][0:3500,:]=0 # Remove rediculously big municipality in Fort Nelson
    #ind=np.where(z['munic']>0); z['munic'][ind]=1
    #z['munic']=gis.BufferRasterMask(z['munic'],1)
    #plt.matshow(z['munic'])
    ind=np.where( (np.isin(z['lc_comp1_2019'],[lut_comp1['Built-up']])==True) & (z['munic']>0) )
    z1['Data'][ind]=lut_lu['Residential and Commercial']
    zRC[z['munic']>0]=1
    
    gdf=gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='BC_MAJOR_CITIES_POINTS_500M')
    gdf['geometry']=gdf.geometry.buffer(2000)
    #plt.close('all'); gdf.plot()
    gdf=gdf[gdf.geometry!=None]
    gdf=gdf.reset_index()
    gdf['ID']=np.ones(len(gdf))
    shapes=((geom,value) for geom, value in zip(gdf['geometry'],gdf['ID']))
    zBuf=np.zeros(zRef['Data'].shape,dtype=float)
    burned=features.rasterize(shapes=shapes,fill=0,out=zBuf,transform=zRef['Transform'])    
    ind=np.where( (np.isin(z['lc_comp1_2019'],[lut_comp1['Built-up'],lut_comp1['Earth and Rock']])==True) & (zBuf>0) )
    z1['Data'][ind]=lut_lu['Residential and Commercial']
    zRC[zBuf>0]=1
    
    gdf=gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='NRC_POPULATED_PLACES_1M_SP')
    gdf['geometry']=gdf.geometry.buffer(2000)
    #plt.close('all'); gdf.plot()
    gdf=gdf[gdf.geometry!=None]
    gdf=gdf.reset_index()
    gdf['ID']=np.ones(len(gdf))
    shapes=((geom,value) for geom, value in zip(gdf['geometry'],gdf['ID']))
    zBuf=np.zeros(zRef['Data'].shape,dtype=float)
    burned=features.rasterize(shapes=shapes,fill=0,out=zBuf,transform=zRef['Transform']) 
    ind=np.where( (np.isin(z['lc_comp1_2019'],[lut_comp1['Built-up'],lut_comp1['Earth and Rock']])==True) & (zBuf>0) )
    z1['Data'][ind]=lut_lu['Residential and Commercial']
    zRC[zBuf>0]=1
    
    gdf=gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='GSR_TMBR_PRCSSING_FAC_SV')
    gdf['geometry']=gdf.geometry.buffer(500)
    #plt.close('all'); gdf.plot()
    gdf=gdf[gdf.geometry!=None]
    gdf=gdf.reset_index()
    gdf['ID']=np.ones(len(gdf))
    shapes=((geom,value) for geom, value in zip(gdf['geometry'],gdf['ID']))
    zBuf=np.zeros(zRef['Data'].shape,dtype=float)
    burned=features.rasterize(shapes=shapes,fill=0,out=zBuf,transform=zRef['Transform'])    
    ind=np.where( (np.isin(z['lc_comp1_2019'],[lut_comp1['Built-up'],lut_comp1['Earth and Rock']])==True) & (zBuf>0) )
    z1['Data'][ind]=lut_lu['Residential and Commercial']
    zRC[zBuf>0]=1
    
    gdf=gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='REC_TENURE_ALPINE_SKI_AREAS_SP')
    gdf['geometry']=gdf.geometry.buffer(50)
    gdf=gdf[gdf.geometry!=None]
    gdf=gdf.reset_index()
    gdf['ID']=np.ones(len(gdf))
    shapes=((geom,value) for geom, value in zip(gdf['geometry'],gdf['ID']))
    zBuf=np.zeros(zRef['Data'].shape,dtype=float)
    burned=features.rasterize(shapes=shapes,fill=0,out=zBuf,transform=zRef['Transform']) 
    ind=np.where( (np.isin(z['lc_comp1_2019'],[lut_comp1['Built-up'],lut_comp1['Earth and Rock']])==True) & (zBuf>0) )
    z1['Data'][ind]=lut_lu['Residential and Commercial']
    zRC[zBuf>0]=1
    
    gdf=gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='CLAB_INDIAN_RESERVES')
    gdf['geometry']=gdf.geometry.buffer(1)
    gdf=gdf[gdf.geometry!=None]
    gdf=gdf.reset_index()
    gdf['ID']=np.ones(len(gdf))
    shapes=((geom,value) for geom, value in zip(gdf['geometry'],gdf['ID']))
    zBuf=np.zeros(zRef['Data'].shape,dtype=float)
    burned=features.rasterize(shapes=shapes,fill=0,out=zBuf,transform=zRef['Transform']) 
    ind=np.where( (np.isin(z['lc_comp1_2019'],[lut_comp1['Built-up'],lut_comp1['Earth and Rock']])==True) & (zBuf>0) )
    z1['Data'][ind]=lut_lu['Residential and Commercial']
    zRC[zBuf>0]=1
    
    # Golf courses
    gdf=gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='GSR_GOLF_COURSES_SV')
    gdf['geometry']=gdf.geometry.buffer(500)
    gdf=gdf[gdf.geometry!=None]
    gdf=gdf.reset_index()
    gdf['ID']=np.ones(len(gdf))
    shapes=((geom,value) for geom, value in zip(gdf['geometry'],gdf['ID']))
    zBuf=np.zeros(zRef['Data'].shape,dtype=float)
    burned=features.rasterize(shapes=shapes,fill=0,out=zBuf,transform=zRef['Transform']) 
    ind=np.where( (zRC==0) & (np.isin(z['lc_comp1_2019'],[lut_comp1['Built-up'],lut_comp1['Earth and Rock'],lut_comp1['Bryoid'],lut_comp1['Shrub']])==True) & (zBuf>0) )
    z1['Data'][ind]=lut_lu['Residential and Commercial']
    
    #--------------------------------------------------------------------------
    # Energy and mines
    #--------------------------------------------------------------------------
    
    # Mines
    gdf=gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='HSP_MJR_MINES_PERMTTD_AREAS_SP')
    gdf['geometry']=gdf.geometry.buffer(500)
    gdf=gdf[gdf.geometry!=None]
    gdf=gdf.reset_index()
    gdf['ID']=np.ones(len(gdf))
    shapes=((geom,value) for geom, value in zip(gdf['geometry'],gdf['ID']))
    zBuf=np.zeros(zRef['Data'].shape,dtype=float)
    burned=features.rasterize(shapes=shapes,fill=0,out=zBuf,transform=zRef['Transform']) 
    ind=np.where( (zRC==0) & (np.isin(z['lc_comp1_2019'],[lut_comp1['Built-up'],lut_comp1['Earth and Rock'],lut_comp1['Bryoid'],lut_comp1['Shrub']])==True) & (zBuf>0) )
    z1['Data'][ind]=lut_lu['Energy and Mines']
    
    # Transmission lines
    gdf=gpd.read_file(meta['Paths']['GDB']['GDB'] + '\\GBA_TRANSMISSION_LINES_SP.geojson')
    gdf['geometry']=gdf.geometry.buffer(200)
    gdf=gdf[gdf.geometry!=None]
    gdf=gdf.reset_index()
    gdf['ID']=np.ones(len(gdf))
    shapes=((geom,value) for geom, value in zip(gdf['geometry'],gdf['ID']))
    zBuf=np.zeros(zRef['Data'].shape,dtype=float)
    burned=features.rasterize(shapes=shapes,fill=0,out=zBuf,transform=zRef['Transform']) 
    ind=np.where( (zRC==0) & (np.isin(z['lc_comp1_2019'],[lut_comp1['Built-up'],lut_comp1['Earth and Rock'],lut_comp1['Bryoid'],lut_comp1['Shrub']])==True) & (zBuf>0) )
    z1['Data'][ind]=lut_lu['Energy and Mines']
    
    # BC Hydro facilities
    gdf=gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='DRP_BCHYDRO_FACILITY_PNTS_SP')
    gdf['geometry']=gdf.geometry.buffer(500)
    gdf=gdf[gdf.geometry!=None]
    gdf=gdf.reset_index()
    gdf['ID']=np.ones(len(gdf))
    shapes=((geom,value) for geom, value in zip(gdf['geometry'],gdf['ID']))
    zBuf=np.zeros(zRef['Data'].shape,dtype=float)
    burned=features.rasterize(shapes=shapes,fill=0,out=zBuf,transform=zRef['Transform']) 
    ind=np.where( (zRC==0) & (np.isin(z['lc_comp1_2019'],[lut_comp1['Built-up'],lut_comp1['Earth and Rock'],lut_comp1['Bryoid'],lut_comp1['Shrub']])==True) & (zBuf>0) )
    z1['Data'][ind]=lut_lu['Energy and Mines']
    
    # Oil and gas facilities
    gdf=gpd.read_file(meta['Paths']['GDB']['GDB'] + '\\ogf.geojson')
    gdf['geometry']=gdf.geometry.buffer(200)
    gdf=gdf[gdf.geometry!=None]
    gdf=gdf.reset_index()
    gdf['ID']=np.ones(len(gdf))
    shapes=((geom,value) for geom, value in zip(gdf['geometry'],gdf['ID']))
    zBuf=np.zeros(zRef['Data'].shape,dtype=float)
    burned=features.rasterize(shapes=shapes,fill=0,out=zBuf,transform=zRef['Transform']) 
    ind=np.where( (np.isin(z['lc_comp1_2019'],[lut_comp1['Built-up'],lut_comp1['Earth and Rock'],lut_comp1['Bryoid'],lut_comp1['Shrub']])==True) &  (zBuf>0) )
    z1['Data'][ind]=lut_lu['Energy and Mines']
    
    #--------------------------------------------------------------------------        
    # Transportation
    #--------------------------------------------------------------------------   
    buffer_width=55 # metres
    
    # Major roads
    gdf=gpd.read_file(meta['Paths']['GDB']['GDB'] + '\\road.geojson')
    gdf['geometry']=gdf.geometry.buffer(buffer_width)
    gdf=gdf[gdf.geometry!=None]
    gdf=gdf.reset_index()
    gdf['ID']=np.ones(len(gdf))
    shapes=((geom,value) for geom, value in zip(gdf['geometry'],gdf['ID']))
    zBuf=np.zeros(zRef['Data'].shape,dtype=float)
    burned=features.rasterize(shapes=shapes,fill=0,out=zBuf,transform=zRef['Transform']) 
    ind=np.where( (zRC==0) & (np.isin(z['lc_comp1_2019'],[lut_comp1['Built-up'],lut_comp1['Earth and Rock'],lut_comp1['Bryoid'],lut_comp1['Shrub']])==True) & (zBuf>0) )
    z1['Data'][ind]=lut_lu['Transportation']
    
    # Forestry roads
    gdf=gpd.read_file(meta['Paths']['GDB']['GDB'] + '\\FTEN_ROAD_SEGMENT_LINES_SVW.geojson')
    gdf['geometry']=gdf.geometry.buffer(buffer_width)
    gdf=gdf[gdf.geometry!=None]
    gdf=gdf.reset_index()
    gdf['ID']=np.ones(len(gdf))
    shapes=((geom,value) for geom, value in zip(gdf['geometry'],gdf['ID']))
    zBuf=np.zeros(zRef['Data'].shape,dtype=float)
    burned=features.rasterize(shapes=shapes,fill=0,out=zBuf,transform=zRef['Transform']) 
    ind=np.where( (zRC==0) & (np.isin(z['lc_comp1_2019'],[lut_comp1['Built-up'],lut_comp1['Earth and Rock'],lut_comp1['Bryoid'],lut_comp1['Shrub']])==True) & (zBuf>0) )
    z1['Data'][ind]=lut_lu['Transportation']
    
    # Oil and gas roads
    gdf=gpd.read_file(meta['Paths']['GDB']['GDB'] + '\\OG_ROAD_SEGMENT_PERMIT_SP.geojson')
    gdf['geometry']=gdf.geometry.buffer(buffer_width)
    gdf=gdf[gdf.geometry!=None]
    gdf=gdf.reset_index()
    gdf['ID']=np.ones(len(gdf))
    shapes=((geom,value) for geom, value in zip(gdf['geometry'],gdf['ID']))
    zBuf=np.zeros(zRef['Data'].shape,dtype=float)
    burned=features.rasterize(shapes=shapes,fill=0,out=zBuf,transform=zRef['Transform']) 
    ind=np.where( (zRC==0) & (np.isin(z['lc_comp1_2019'],[lut_comp1['Built-up'],lut_comp1['Earth and Rock'],lut_comp1['Bryoid'],lut_comp1['Shrub']])==True) & (zBuf>0) )
    z1['Data'][ind]=lut_lu['Transportation']
    
    # Oil and gas pipelines
    gdf=gpd.read_file(meta['Paths']['GDB']['GDB'] + '\\DRP_OIL_GAS_PIPELINES_BC_SP.geojson')
    gdf['geometry']=gdf.geometry.buffer(buffer_width)
    gdf=gdf[gdf.geometry!=None]
    gdf=gdf.reset_index()
    gdf['ID']=np.ones(len(gdf))
    shapes=((geom,value) for geom, value in zip(gdf['geometry'],gdf['ID']))
    zBuf=np.zeros(zRef['Data'].shape,dtype=float)
    burned=features.rasterize(shapes=shapes,fill=0,out=zBuf,transform=zRef['Transform']) 
    ind=np.where( (zRC==0) & (np.isin(z['lc_comp1_2019'],[lut_comp1['Built-up'],lut_comp1['Earth and Rock'],lut_comp1['Bryoid'],lut_comp1['Shrub']])==True) & (zBuf>0) )
    z1['Data'][ind]=lut_lu['Transportation']    
    
    # Airports
    gdf=gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='GSR_AIRPORTS_SVW')
    gdf['geometry']=gdf.geometry.buffer(500)
    gdf=gdf[gdf.geometry!=None]
    gdf=gdf.reset_index()
    gdf['ID']=np.ones(len(gdf))
    shapes=((geom,value) for geom, value in zip(gdf['geometry'],gdf['ID']))
    zBuf=np.zeros(zRef['Data'].shape,dtype=float)
    burned=features.rasterize(shapes=shapes,fill=0,out=zBuf,transform=zRef['Transform']) 
    ind=np.where( (zRC==0) & (np.isin(z['lc_comp1_2019'],[lut_comp1['Built-up'],lut_comp1['Earth and Rock'],lut_comp1['Bryoid'],lut_comp1['Shrub']])==True) & (zBuf>0) )
    z1['Data'][ind]=lut_lu['Transportation']
    
    # Ports
    gdf=gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='GSR_PORTS_TERMINALS_SVW')
    gdf['geometry']=gdf.geometry.buffer(500)
    gdf=gdf[gdf.geometry!=None]
    gdf=gdf.reset_index()
    gdf['ID']=np.ones(len(gdf))
    shapes=((geom,value) for geom, value in zip(gdf['geometry'],gdf['ID']))
    zBuf=np.zeros(zRef['Data'].shape,dtype=float)
    burned=features.rasterize(shapes=shapes,fill=0,out=zBuf,transform=zRef['Transform']) 
    ind=np.where( (zRC==0) & (np.isin(z['lc_comp1_2019'],[lut_comp1['Built-up'],lut_comp1['Earth and Rock'],lut_comp1['Bryoid'],lut_comp1['Shrub']])==True) & (zBuf>0) )
    z1['Data'][ind]=lut_lu['Transportation']
    
    #--------------------------------------------------------------------------
    # Water Management
    #--------------------------------------------------------------------------    
    gdf=gpd.read_file(meta['Paths']['GDB']['LandCover'],layer='FWA_MANMADE_WATERBODIES_POLY')
    gdf['geometry']=gdf.geometry.buffer(50)
    gdf=gdf[gdf.geometry!=None]
    gdf=gdf.reset_index()
    gdf['ID']=np.ones(len(gdf))
    shapes=((geom,value) for geom, value in zip(gdf['geometry'],gdf['ID']))
    zBuf=np.zeros(zRef['Data'].shape,dtype=float)
    burned=features.rasterize(shapes=shapes,fill=0,out=zBuf,transform=zRef['Transform']) 
    ind=np.where( (zRC==0) & (np.isin(z['lc_comp1_2019'],[lut_comp1['Built-up'],lut_comp1['Earth and Rock'],lut_comp1['Bryoid'],lut_comp1['Shrub'],lut_comp1['Water']])==True) & (zBuf>0) )
    z1['Data'][ind]=lut_lu['Water Management']
    ind=np.where( (np.isin(z['lc_comp1_2019'],[lut_comp1['Water']])==True) )
    z1['Data'][ind]=lut_lu['Water Management']
    
    #--------------------------------------------------------------------------
    # Timber production
    #--------------------------------------------------------------------------
    ind=np.where( (z['lc_comp1_2019']==lut_comp1['Forest']) & (z1['Data']==0) )
    z1['Data'][ind]=lut_lu['Timber Production Passive']
    
    #--------------------------------------------------------------------------
    # Transporation 2
    # There's a lot of logging roads that don't get picked up. We can take advanage
    # of the change detected by CEC between 2010 and 2020, applying to areas
    # where no other designation is found above.
    #--------------------------------------------------------------------------
    vList=['lc_cec_2010','lc_cec_2020','lc_vri_l3'] 
    z2=Import_Raster(meta,[],vList,'Extract Grid')
    
    MaskA=np.zeros(zRef['Data'].shape,dtype='int8')
    ind=np.where(z2['lc_vri_l3']==lut_lc_vri_l3['A'])
    MaskA[ind]=1
    #plt.matshow(MaskA)
    MaskA=gis.BufferRasterMask(MaskA,5)
    #plt.matshow(MaskA2)
    #ind=np.where( (MaskA==0) & (MaskA2==1) )
    #MaskA[ind]=2
    
    ind=np.where( (z2['lc_cec_2010']==lut_cec['Forest']) & (z2['lc_cec_2020']!=lut_cec['Forest']) & (np.isin(z['lc_comp1_2019'],[lut_comp1['Shrub'],lut_comp1['Wetland'],lut_comp1['Herb']])==False) & (MaskA==0) & (z1['Data']==0) )
    z1['Data'][ind]=lut_lu['Transportation']
    
    #--------------------------------------------------------------------------
    # No designation
    #--------------------------------------------------------------------------
    ind=np.where( (z1['Data']==0) & (zRef['Data']==1) )
    z1['Data'][ind]=lut_lu['No Designation']
    
    #--------------------------------------------------------------------------
    # Non-land
    #--------------------------------------------------------------------------
    ind=np.where( (zRef['Data']==0) )
    z1['Data'][ind]=lut_lu['No Designation']+1
    
    #--------------------------------------------------------------------------
    # Save
    #--------------------------------------------------------------------------
    # plt.close('all'); plt.matshow(z1['Data'],clim=[0,12])  
    gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandUse_Comp1_2019.tif')
    
    return

#%% This needs to be run after the LU comp!
def DeriveLandCoverComp1_1800(meta):
    zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
    vList=['lc_comp1_2019','lu_comp1_2019']
    z=u1ha.Import_Raster(meta,[],vList,'Extract Grid')
    lc_comp1=meta['LUT']['Derived']['lc_comp1']
    lu_comp1=meta['LUT']['Derived']['lu_comp1']
    
    z1=copy.deepcopy(zRef)
    z1['Data']=z['lc_comp1_2019'].copy()
    ind=np.where( (z['lc_comp1_2019']!=lc_comp1['Forest']) & (np.isin(z['lu_comp1_2019'],[lu_comp1['Residential and Commercial'],lu_comp1['Energy and Mines'],lu_comp1['Transportation']])==True) )
    z1['Data'][ind]=lc_comp1['Forest']
    
    ind=np.where( (z['lc_comp1_2019']!=lc_comp1['Forest']) & (np.isin(z['lu_comp1_2019'],[lu_comp1['Cropland'],lu_comp1['Pasture']])==True) )
    rn=np.random.random(ind[0].size)
    indF=np.where(rn<0.7)[0]
    z1['Data'][ ind[0][indF],ind[1][indF] ]=lc_comp1['Forest']
    indF=np.where( (rn>=0.7) & (rn<0.9) )[0]
    z1['Data'][ ind[0][indF],ind[1][indF] ]=lc_comp1['Shrub']
    indF=np.where( (rn>=0.9) )[0]
    z1['Data'][ ind[0][indF],ind[1][indF] ]=lc_comp1['Herb']
    gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCover_Comp1_1800.tif')
    return

#%% Derive land use compilation 2049 scenarios
def DeriveLandCoverLandUseComp1_2020to2049_Scenarios(meta):
    
    zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])        
    lut_lc=meta['LUT']['Derived']['lc_comp1']
    lut_lu=meta['LUT']['Derived']['lu_comp1']        
    lut_chng=meta['LUT']['Derived']['lclu_chng_comp1']
    vList=['lc_comp1_2019','lu_comp1_2019','harv_yr_comp2','ogdef','d2road','d2fac','d2set'] 
    z0=u1ha.Import_Raster(meta,[],vList,'Extract Grid')
    tv=np.arange(2020,2050,1)
    
    # Annual prob deforestation
    P_Def0=0.00065
    
    # Track area
    A={}
    A['For 2019']=np.where(z0['lc_comp1_2019']==lut_lc['Forest'])[0].size/1e6
    ind=np.where( (zRef['Data']==1) & (z0['lu_comp1_2019']==lut_lu['Conservation Natural']) | (zRef['Data']==1) & (z0['lu_comp1_2019']==lut_lu['Conservation Consistent']) )
    A['Tot Con 2019']=ind[0].size/1e6
    ind=np.where( (z0['lc_comp1_2019']==lut_lc['Forest']) & (zRef['Data']==1) & (z0['lu_comp1_2019']==lut_lu['Conservation Natural']) | (z0['lc_comp1_2019']==lut_lc['Forest']) & (zRef['Data']==1) & (z0['lu_comp1_2019']==lut_lu['Conservation Consistent']) )
    A['For Con 2019']=ind[0].size/1e6
    print('Current fraction of forest conserved:' + str(A['For Con 2019']/A['For 2019']))
    #A['Con 2049']=2.0*A['Con 2019']
    A['For Con 2049']=(0.3*A['For 2019'])
    A['OGD']=np.where( (z0['ogdef']>0) )[0].size/1e6
        
    #==========================================================================
    # LCLU-Comp1 2020to2049 Scenario 1: LUC-BAU / C-Hist / WM-None / FM-None
    #==========================================================================
    
    #ind=np.where( (z0['lu_comp1_2019']==lut_lu['Conservation Consistent']) & (z1['LU 2049']==lut_lu['Conservation Consistent']) )
    
    z1={}    
    z1['LC 2049']=z0['lc_comp1_2019'].copy()
    z1['LU 2049']=z0['lu_comp1_2019'].copy()
    z1['Year Change']=np.zeros(zRef['Data'].shape,dtype='int16')
    z1['Type Change']=np.zeros(zRef['Data'].shape,dtype='int8')
    for iT in range(tv.size):
        print(tv[iT])
        #--------------------------------------------------------------------------
        # Add LUC
        # Average NIR/PIR annual forest conversion = 7955ha/yr (0.01%/yr)
        # Modelling: x=np.arange(0,100,1); y=(1/(1+np.exp(0.14*(x-25)))); plt.close('all'); plt.plot(x,y,'r-')
        #--------------------------------------------------------------------------                    
        rn0=np.random.random(zRef['Data'].shape)
        P_Def=P_Def0*(1/(1+np.exp(0.14*(z0['d2road']-25))))*(1/(1+np.exp(0.14*(z0['d2set']-25))))
        #plt.close('all'); plt.matshow(P_Def); plt.colorbar()    
        #ind=np.where( (z1['LC 2049']==lut_lc['Forest']) & (rn0<P_Def) )
        #A_Def=ind[0].size
        #print(A_Def/A['For 2019']*100)
        #rn=-1*np.ones(zRef['Data'].shape)
        rn1=np.random.random(zRef['Data'].shape)
        
        ind=np.where( (z1['LC 2049']==lut_lc['Forest']) & (z0['lu_comp1_2019']!=lut_lu['Conservation Natural']) & (z0['lu_comp1_2019']!=lut_lu['Conservation Consistent']) & (rn0<P_Def) & (rn1>=0) & (rn1<0.15) )
        z1['LC 2049'][ind]=lut_lc['Herb']
        z1['LU 2049'][ind]=lut_lu['Cropland']
        z1['Year Change'][ind]=tv[iT]
        z1['Type Change'][ind]=lut_chng['FL-CL']
        
        ind=np.where( (z1['LC 2049']==lut_lc['Forest']) & (z0['lu_comp1_2019']!=lut_lu['Conservation Natural']) & (z0['lu_comp1_2019']!=lut_lu['Conservation Consistent']) & (rn0<P_Def) & (rn1>=0.15) & (rn1<0.45) )
        z1['LC 2049'][ind]=lut_lc['Shrub']
        z1['LU 2049'][ind]=lut_lu['Pasture']
        z1['Year Change'][ind]=tv[iT]
        z1['Type Change'][ind]=lut_chng['FL-PA']
        
        ind=np.where( (z1['LC 2049']==lut_lc['Forest']) & (z0['lu_comp1_2019']!=lut_lu['Conservation Natural']) & (z0['lu_comp1_2019']!=lut_lu['Conservation Consistent']) & (rn0<P_Def) & (rn1>=0.45) & (rn1<0.51) )
        z1['LC 2049'][ind]=lut_lc['Built-up']
        z1['LU 2049'][ind]=lut_lu['Residential and Commercial']        
        z1['Year Change'][ind]=tv[iT]
        z1['Type Change'][ind]=lut_chng['FL-RC']
        
        ind=np.where( (z1['LC 2049']==lut_lc['Forest']) & (z0['lu_comp1_2019']!=lut_lu['Conservation Natural']) & (z0['lu_comp1_2019']!=lut_lu['Conservation Consistent']) & (rn0<P_Def) & (rn1>=0.51) & (rn1<0.85) )
        z1['LC 2049'][ind]=lut_lc['Shrub']
        z1['LU 2049'][ind]=lut_lu['Transportation']
        z1['Year Change'][ind]=tv[iT]
        z1['Type Change'][ind]=lut_chng['FL-TR']
        
        ind=np.where( (z1['LC 2049']==lut_lc['Forest']) & (z0['lu_comp1_2019']!=lut_lu['Conservation Natural']) & (z0['lu_comp1_2019']!=lut_lu['Conservation Consistent']) & (rn0<P_Def) & (rn1>=0.85) )
        z1['LC 2049'][ind]=lut_lc['Built-up']
        z1['LU 2049'][ind]=lut_lu['Energy and Mines']
        z1['Year Change'][ind]=tv[iT]
        z1['Type Change'][ind]=lut_chng['FL-EM']
    
    #plt.close('all'); plt.matshow(z1['LU 2049']); plt.colorbar()    
    z2=copy.deepcopy(zRef); z2['Data']=z1['LC 2049']
    gis.SaveGeoTiff(z2,meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCover_Comp1_2049_Scn1.tif')
    z2=copy.deepcopy(zRef); z2['Data']=z1['LU 2049']
    gis.SaveGeoTiff(z2,meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandUse_Comp1_2049_Scn1.tif')
    z2=copy.deepcopy(zRef); z2['Data']=z1['Type Change']
    gis.SaveGeoTiff(z2,meta['Paths']['bc1ha'] + '\\Disturbances\\LandUseChange_Comp1_2020to2049_Scn1_Type.tif')
    z2=copy.deepcopy(zRef); z2['Data']=z1['Year Change']
    gis.SaveGeoTiff(z2,meta['Paths']['bc1ha'] + '\\Disturbances\\LandUseChange_Comp1_2020to2049_Scn1_Year.tif')
    
    #==========================================================================
    # LCLU-Comp1 2020to2049 Scenario 2: LUC-BAU / C-CAP / FM-Hist
    #==========================================================================
    
    z1={}    
    z1['LC 2049']=z0['lc_comp1_2019'].copy()
    z1['LU 2049']=z0['lu_comp1_2019'].copy()
    z1['Year Change']=np.zeros(zRef['Data'].shape,dtype='int16')
    z1['Type Change']=np.zeros(zRef['Data'].shape,dtype='int8')
    
    # Add old growth deferrals
    ind1=np.where( (z0['harv_yr_comp2']==0) & (z0['ogdef']>0) )
    z1['LU 2049'][ind1]=lut_lu['Conservation Natural']
    z1['Year Change'][ind1]=2020
    ind2=np.where( (z0['harv_yr_comp2']>0) & (z0['ogdef']>0) )
    z1['LU 2049'][ind2]=lut_lu['Conservation Consistent']
    z1['Year Change'][ind2]=2020
       
    for iT in range(tv.size):
        print(tv[iT])
        #--------------------------------------------------------------------------
        # Add LUC
        # Average NIR/PIR annual forest conversion = 7955ha/yr (0.01%/yr)
        # Modelling: x=np.arange(0,100,1); y=(1/(1+np.exp(0.14*(x-25)))); plt.close('all'); plt.plot(x,y,'r-')
        #--------------------------------------------------------------------------                    
        rn0=np.random.random(zRef['Data'].shape)
        P_Def=P_Def0*(1/(1+np.exp(0.14*(z0['d2road']-25))))*(1/(1+np.exp(0.14*(z0['d2set']-25))))
        #plt.close('all'); plt.matshow(P_Def); plt.colorbar()    
        #ind=np.where( (z1['LC 2049']==lut_lc['Forest']) & (rn0<P_Def) )
        #A_Def=ind[0].size
        #print(A_Def/A['For 2019']*100)
        #rn=-1*np.ones(zRef['Data'].shape)
        rn1=np.random.random(zRef['Data'].shape)
        
        ind=np.where( (z1['LC 2049']==lut_lc['Forest']) & \
                     (z0['lu_comp1_2019']!=lut_lu['Conservation Natural']) & (z0['lu_comp1_2019']!=lut_lu['Conservation Consistent']) & \
                     (z1['LU 2049']!=lut_lu['Conservation Natural']) & (z1['LU 2049']!=lut_lu['Conservation Consistent']) & \
                     (rn0<P_Def) & (rn1>=0) & (rn1<0.15) )
        z1['LC 2049'][ind]=lut_lc['Herb']
        z1['LU 2049'][ind]=lut_lu['Cropland']
        z1['Year Change'][ind]=tv[iT]
        z1['Type Change'][ind]=lut_chng['FL-CL']
        
        ind=np.where( (z1['LC 2049']==lut_lc['Forest']) & \
                     (z0['lu_comp1_2019']!=lut_lu['Conservation Natural']) & (z0['lu_comp1_2019']!=lut_lu['Conservation Consistent']) & \
                     (z1['LU 2049']!=lut_lu['Conservation Natural']) & (z1['LU 2049']!=lut_lu['Conservation Consistent']) & \
                     (rn0<P_Def) & (rn1>=0.15) & (rn1<0.45) )
        z1['LC 2049'][ind]=lut_lc['Shrub']
        z1['LU 2049'][ind]=lut_lu['Pasture']
        z1['Year Change'][ind]=tv[iT]
        z1['Type Change'][ind]=lut_chng['FL-PA']
        
        ind=np.where( (z1['LC 2049']==lut_lc['Forest']) & \
                     (z0['lu_comp1_2019']!=lut_lu['Conservation Natural']) & (z0['lu_comp1_2019']!=lut_lu['Conservation Consistent']) & \
                     (z1['LU 2049']!=lut_lu['Conservation Natural']) & (z1['LU 2049']!=lut_lu['Conservation Consistent']) & \
                     (rn0<P_Def) & (rn1>=0.45) & (rn1<0.51) )
        z1['LC 2049'][ind]=lut_lc['Built-up']
        z1['LU 2049'][ind]=lut_lu['Residential and Commercial']
        z1['Year Change'][ind]=tv[iT]
        z1['Type Change'][ind]=lut_chng['FL-RC']
        
        ind=np.where( (z1['LC 2049']==lut_lc['Forest']) & \
                     (z0['lu_comp1_2019']!=lut_lu['Conservation Natural']) & (z0['lu_comp1_2019']!=lut_lu['Conservation Consistent']) & \
                     (z1['LU 2049']!=lut_lu['Conservation Natural']) & (z1['LU 2049']!=lut_lu['Conservation Consistent']) & \
                     (rn0<P_Def) & (rn1>=0.51) & (rn1<0.85) )
        z1['LC 2049'][ind]=lut_lc['Shrub']
        z1['LU 2049'][ind]=lut_lu['Transportation']
        z1['Year Change'][ind]=tv[iT]
        z1['Type Change'][ind]=lut_chng['FL-TR']
        
        ind=np.where( (z1['LC 2049']==lut_lc['Forest']) & \
                     (z0['lu_comp1_2019']!=lut_lu['Conservation Natural']) & (z0['lu_comp1_2019']!=lut_lu['Conservation Consistent']) & \
                     (z1['LU 2049']!=lut_lu['Conservation Natural']) & (z1['LU 2049']!=lut_lu['Conservation Consistent']) & \
                     (rn0<P_Def) & (rn1>=0.85) )
        z1['LC 2049'][ind]=lut_lc['Built-up']
        z1['LU 2049'][ind]=lut_lu['Energy and Mines']
        z1['Year Change'][ind]=tv[iT]
        z1['Type Change'][ind]=lut_chng['FL-EM']
    
    #--------------------------------------------------------------------------
    # Add Conservation
    #--------------------------------------------------------------------------        
    
    # Add small areas
    seedval=10
    rng=default_rng(seed=seedval)
    sig=21
    noiseS=rng.integers(0,1000,zRef['Data'].shape,'int16',True)
    blurS=cv2.GaussianBlur(noiseS,(0,0),sigmaX=sig,sigmaY=sig,borderType=cv2.BORDER_DEFAULT)
    # plt.close('all'); plt.matshow(blurS)
    u=np.unique(blurS)
    #print(u)
    
    MaskS=np.zeros(blurS.shape,dtype='int8')
    ind=np.where(blurS<=u[15]); MaskS[ind]=1
    plt.close('all'); plt.matshow(MaskS)
    ind=np.where( (MaskS==1) & (zRef['Data']==1) & (np.isin(z1['LU 2049'],[lut_lu['Timber Production Passive'],lut_lu['Timber Production Intense']])==True) ); #Mask2[ind]=1
    A_PropS=ind[0].size/1e6
    print(A_PropS)
    print( (A['For Con 2019']+A['OGD']+A_PropS)/A['For Con 2049'] )
    
    # Add new areas
    ind=np.where( (MaskS==1) & (zRef['Data']==1) & (np.isin(z0['lu_comp1_2019'],[lut_lu['Timber Production Passive'],lut_lu['Timber Production Intense']])==True) & (z0['harv_yr_comp2']==0) )    
    z1['LU 2049'][ind]=lut_lu['Conservation Natural']
    ind=np.where( (MaskS==1) & (zRef['Data']==1) & (np.isin(z0['lu_comp1_2019'],[lut_lu['Timber Production Passive'],lut_lu['Timber Production Intense']])==True) & (z0['harv_yr_comp2']>0) )
    z1['LU 2049'][ind]=lut_lu['Conservation Consistent']
    
    # Save
    z2=copy.deepcopy(zRef); z2['Data']=z1['LC 2049']
    gis.SaveGeoTiff(z2,meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCover_Comp1_2049_Scn2.tif')
    z2=copy.deepcopy(zRef); z2['Data']=z1['LU 2049']
    gis.SaveGeoTiff(z2,meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandUse_Comp1_2049_Scn2.tif')
    z2=copy.deepcopy(zRef); z2['Data']=z1['Type Change']
    gis.SaveGeoTiff(z2,meta['Paths']['bc1ha'] + '\\Disturbances\\LandUseChange_Comp1_2020to2049_Scn2_Type.tif')
    z2=copy.deepcopy(zRef); z2['Data']=z1['Year Change']
    gis.SaveGeoTiff(z2,meta['Paths']['bc1ha'] + '\\Disturbances\\LandUseChange_Comp1_2020to2049_Scn2_Year.tif')
   
    #==========================================================================
    # LCLU-Comp1 2020to2049 Scenario 3: LUC-BAU / C-CAP / FM-CAP1
    #==========================================================================    
    # Import required data
    zLC=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCover_Comp1_2049_Scn2.tif')
    zLU=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandUse_Comp1_2049_Scn2.tif')
    zSlope=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Terrain\\slope.tif')['Data']
    #z1['Year Change']=np.zeros(zRef['Data'].shape,dtype='int16')
    #z1['Type Change']=np.zeros(zRef['Data'].shape,dtype='int8')
    
    # Create Perlin noise
    noise=PerlinNoise(octaves=10,seed=3)
    xpix,ypix=1000,1000
    zNoise0=np.zeros((xpix,ypix))
    for i in range(xpix):
        for j in range(ypix):
            zNoise0[i,j]=noise([i/xpix, j/ypix])            
    zNoise2=np.zeros(zNoise0.shape); 
    ind=np.where( (zNoise0<-0.12) | (zNoise0>0.12) ); 
    zNoise2[ind]=1
    zNoise2=np.tile(zNoise2,(int(np.ceil(zRef['m']/1000)),int(np.ceil(zRef['n']/1000))))
    zNoise2=zNoise2[0:zRef['m'],0:zRef['n']]  
    # plt.close('all'); plt.matshow(zNoise2)

    #--------------------------------------------------------------------------
    # Add Bands of intense forestry zone
    gdf=gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='BC_MAJOR_CITIES_POINTS_500M')
    gdf['geometry']=gdf.geometry.buffer(7*1000)
    gdf=gdf[gdf.geometry!=None]
    gdf=gdf.reset_index()
    gdf['ID']=np.ones(len(gdf))
    shapes=((geom,value) for geom, value in zip(gdf['geometry'],gdf['ID']))
    zInner=np.zeros(zRef['Data'].shape,dtype=float)
    burned=features.rasterize(shapes=shapes,fill=0,out=zInner,transform=zRef['Transform'])
    
    gdf=gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='BC_MAJOR_CITIES_POINTS_500M')
    gdf['geometry']=gdf.geometry.buffer(55*1000)
    gdf=gdf[gdf.geometry!=None]
    gdf=gdf.reset_index()
    gdf['ID']=np.ones(len(gdf))
    shapes=((geom,value) for geom, value in zip(gdf['geometry'],gdf['ID']))
    zOuter=np.zeros(zRef['Data'].shape,dtype=float)
    burned=features.rasterize(shapes=shapes,fill=0,out=zOuter,transform=zRef['Transform'])
    
    # gdf=gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='BC_MAJOR_CITIES_POINTS_500M')
    # gdf['geometry']=gdf.geometry.buffer(56*1000)
    # gdf=gdf[gdf.geometry!=None]
    # gdf=gdf.reset_index()
    # gdf['ID']=np.ones(len(gdf))
    # shapes=((geom,value) for geom, value in zip(gdf['geometry'],gdf['ID']))
    # zOuter2=np.zeros(zRef['Data'].shape,dtype=float)
    # burned=features.rasterize(shapes=shapes,fill=0,out=zOuter2,transform=zRef['Transform'])
    
    # gdf=gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='BC_MAJOR_CITIES_POINTS_500M')
    # gdf['geometry']=gdf.geometry.buffer(6*1000)
    # gdf=gdf[gdf.geometry!=None]
    # gdf=gdf.reset_index()
    # gdf['ID']=np.ones(len(gdf))
    # shapes=((geom,value) for geom, value in zip(gdf['geometry'],gdf['ID']))
    # zInner2=np.zeros(zRef['Data'].shape,dtype=float)
    # burned=features.rasterize(shapes=shapes,fill=0,out=zInner2,transform=zRef['Transform'])
    
    # zCb=np.zeros(zRef['Data'].shape,dtype='int8')
    # ind=np.where( (zF==1) & (zOuter==1) )
    # zCb[ind]=1
    # plt.close('all'); plt.matshow(zCb)
    
    gdf=gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='GSR_TMBR_PRCSSING_FAC_SV')
    gdf['geometry']=gdf.geometry.buffer(75*1000)
    gdf=gdf[gdf.geometry!=None]
    gdf=gdf.reset_index()
    gdf['ID']=np.ones(len(gdf))
    shapes=((geom,value) for geom, value in zip(gdf['geometry'],gdf['ID']))
    zF=np.zeros(zRef['Data'].shape,dtype=float)
    burned=features.rasterize(shapes=shapes,fill=0,out=zF,transform=zRef['Transform'])
    
    zMaskIFM=np.zeros(zRef['Data'].shape,dtype='int8')
    ind=np.where( (zInner==0) & (zOuter==1) & (zF==1) )
    zMaskIFM[ind]=1
    #plt.close('all'); plt.matshow(zMaskIFM)
        
    # Add intensee forest management
    ListLU=[lut_lu['Timber Production Passive'],lut_lu['No Designation']]
    ind=np.where( (zLC['Data']==lut_lc['Forest']) & (np.isin(zLU['Data'],ListLU)==True) & (zMaskIFM==1) & (zNoise2==1) )
    zLU['Data'][ind]=lut_lu['Timber Production Intense']
    
    # Add wildfire management zones
    ind=np.where( (zLC['Data']==lut_lc['Forest']) & (np.isin(zLU['Data'],ListLU)==True) & (zInner==1) )
    zLU['Data'][ind]=lut_lu['Wildfire Risk Management']
    
    # Save    
    gis.SaveGeoTiff(zLU,meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandUse_Comp1_2049_Scn3.tif')
   
    #==========================================================================
    # LCLU-Comp1 2020to2049 Scenario 4: LUC-BAU / C-CAP / FM-CAP2
    #==========================================================================    
    # Import required data
    zLC=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCover_Comp1_2049_Scn2.tif')
    zLU=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandUse_Comp1_2049_Scn2.tif')
    zSlope=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Terrain\\slope.tif')['Data']
    
    #--------------------------------------------------------------------------
    # Add Bands of intense forestry zone
    noise=PerlinNoise(octaves=10,seed=3)
    xpix,ypix=1000,1000
    zNoise0=np.zeros((xpix,ypix))
    for i in range(xpix):
        for j in range(ypix):
            zNoise0[i,j]=noise([i/xpix, j/ypix])            
    zNoise2=np.zeros(zNoise0.shape); 
    ind=np.where( (zNoise0<-0.12) | (zNoise0>0.12) ); 
    zNoise2[ind]=1
    zNoise2=np.tile(zNoise2,(int(np.ceil(zRef['m']/1000)),int(np.ceil(zRef['n']/1000))))
    zNoise2=zNoise2[0:zRef['m'],0:zRef['n']]  
    # plt.close('all'); plt.matshow(zNoise2)

    gdf=gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='BC_MAJOR_CITIES_POINTS_500M')
    gdf['geometry']=gdf.geometry.buffer(7*1000)
    gdf=gdf[gdf.geometry!=None]
    gdf=gdf.reset_index()
    gdf['ID']=np.ones(len(gdf))
    shapes=((geom,value) for geom, value in zip(gdf['geometry'],gdf['ID']))
    zInner=np.zeros(zRef['Data'].shape,dtype=float)
    burned=features.rasterize(shapes=shapes,fill=0,out=zInner,transform=zRef['Transform'])
    
    gdf=gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='BC_MAJOR_CITIES_POINTS_500M')
    gdf['geometry']=gdf.geometry.buffer(55*1000)
    gdf=gdf[gdf.geometry!=None]
    gdf=gdf.reset_index()
    gdf['ID']=np.ones(len(gdf))
    shapes=((geom,value) for geom, value in zip(gdf['geometry'],gdf['ID']))
    zOuter=np.zeros(zRef['Data'].shape,dtype=float)
    burned=features.rasterize(shapes=shapes,fill=0,out=zOuter,transform=zRef['Transform'])
    
    gdf=gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='GSR_TMBR_PRCSSING_FAC_SV')
    gdf['geometry']=gdf.geometry.buffer(75*1000)
    gdf=gdf[gdf.geometry!=None]
    gdf=gdf.reset_index()
    gdf['ID']=np.ones(len(gdf))
    shapes=((geom,value) for geom, value in zip(gdf['geometry'],gdf['ID']))
    zF=np.zeros(zRef['Data'].shape,dtype=float)
    burned=features.rasterize(shapes=shapes,fill=0,out=zF,transform=zRef['Transform'])
    
    zMaskIFM=np.zeros(zRef['Data'].shape,dtype='int8')
    ind=np.where( (zInner==0) & (zOuter==1) & (zF==1) )
    zMaskIFM[ind]=1
    #plt.close('all'); plt.matshow(zMaskIFM)
        
    # Add intensee forest management
    ListLU=[lut_lu['Timber Production Passive'],lut_lu['No Designation']]
    ind=np.where( (zLC['Data']==lut_lc['Forest']) & (np.isin(zLU['Data'],ListLU)==True) & (zMaskIFM==1) & (zNoise2==1) )
    zLU['Data'][ind]=lut_lu['Timber Production Intense']
    
    # Add wildfire management zones
    ind=np.where( (zLC['Data']==lut_lc['Forest']) & (np.isin(zLU['Data'],ListLU)==True) & (zInner==1) )
    zLU['Data'][ind]=lut_lu['Wildfire Risk Management']
    
    #--------------------------------------------------------------------------
    # Add bioenergy production
    
    # Fibrous pattern for bioenergy 
    noise=PerlinNoise(octaves=35,seed=3)
    xpix,ypix=1000,1000
    zNoise0=np.zeros((xpix,ypix))
    for i in range(xpix):
        for j in range(ypix):
            zNoise0[i,j]=noise([i/xpix, j/ypix])
    zNoise1=np.zeros(zNoise0.shape); 
    ind=np.where(np.abs(zNoise0)<0.05); 
    zNoise1[ind]=1    
    zNoise1=np.tile(zNoise1,(int(np.ceil(zRef['m']/1000)),int(np.ceil(zRef['n']/1000))))
    zNoise1=zNoise1[0:zRef['m'],0:zRef['n']]  
    #plt.close('all'); plt.matshow(zNoise1)    
    
    gdf=gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='BC_MAJOR_CITIES_POINTS_500M')
    gdf['geometry']=gdf.geometry.buffer(60*1000)
    gdf=gdf[gdf.geometry!=None]
    gdf=gdf.reset_index()
    gdf['ID']=np.ones(len(gdf))
    shapes=((geom,value) for geom, value in zip(gdf['geometry'],gdf['ID']))
    zOuterBE=np.zeros(zRef['Data'].shape,dtype=float)
    burned=features.rasterize(shapes=shapes,fill=0,out=zOuterBE,transform=zRef['Transform'])
    
    Mask1=np.zeros(zRef['Data'].shape,dtype='int8')
    ListLU=[lut_lu['Timber Production Intense'],lut_lu['Timber Production Passive'],lut_lu['No Designation']]
    ind=np.where( (zLC['Data']==lut_lc['Forest']) & (np.isin(zLU['Data'],ListLU)==True) & (zSlope<20) & (zNoise1==1) & (zOuterBE==1) )
    Mask1[ind]=1
    plt.close('all'); plt.matshow(Mask1)
    print(ind[0].size/1e6)
    zLU['Data'][ind]=lut_lu['Bioenergy Production']
    
    # Save    
    gis.SaveGeoTiff(zLU,meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandUse_Comp1_2049_Scn4.tif')
    
    
    
    
    # # Import required data
    # zLC=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCover_Comp1_2049_Scn2.tif')
    # zLU=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandUse_Comp1_2049_Scn2.tif')
    # zSlope=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Terrain\\slope.tif')['Data']
    # #z1['Year Change']=np.zeros(zRef['Data'].shape,dtype='int16')
    # #z1['Type Change']=np.zeros(zRef['Data'].shape,dtype='int8')
    
    # # Create Perlin noise
    # noise=PerlinNoise(octaves=35,seed=3)
    # xpix,ypix=1000,1000
    # zNoise0=np.zeros((xpix,ypix))
    # for i in range(xpix):
    #     for j in range(ypix):
    #         zNoise0[i,j]=noise([i/xpix, j/ypix])
    
    # # Fibrous pattern for bioenergy 
    # zNoise1=np.zeros(zNoise0.shape); 
    # ind=np.where(np.abs(zNoise0)<0.05); 
    # zNoise1[ind]=1    
    # zNoise1=np.tile(zNoise1,(int(np.ceil(zRef['m']/1000)),int(np.ceil(zRef['n']/1000))))
    # zNoise1=zNoise1[0:zRef['m'],0:zRef['n']]  
    # #plt.close('all'); plt.matshow(zNoise1)
    
    # # Exterior
    # noise=PerlinNoise(octaves=10,seed=3)
    # xpix,ypix=1000,1000
    # zNoise0=np.zeros((xpix,ypix))
    # for i in range(xpix):
    #     for j in range(ypix):
    #         zNoise0[i,j]=noise([i/xpix, j/ypix])            
    # zNoise2=np.zeros(zNoise0.shape); 
    # ind=np.where( (zNoise0<-0.12) | (zNoise0>0.12) ); 
    # zNoise2[ind]=1
    # zNoise2=np.tile(zNoise2,(int(np.ceil(zRef['m']/1000)),int(np.ceil(zRef['n']/1000))))
    # zNoise2=zNoise2[0:zRef['m'],0:zRef['n']]  
    # # plt.close('all'); plt.matshow(zNoise2)
    
    # #--------------------------------------------------------------------------
    # # Add bioenergy production    
    # Mask1=np.zeros(zRef['Data'].shape,dtype='int8')
    # ListLU=[lut_lu['Timber Production Passive'],lut_lu['No Designation']]
    # ind=np.where( (zLC['Data']==lut_lc['Forest']) & (np.isin(zLU['Data'],ListLU)==True) & (zSlope<20) & (zNoise1==1) )
    # Mask1[ind]=1
    # #plt.close('all'); plt.matshow(Mask1)
    # print(ind[0].size/1e6)
    # zLU['Data'][ind]=lut_lu['Bioenergy Production']

    # # #--------------------------------------------------------------------------
    # # # Wildfire management
    # # zMaskWM=np.zeros(zRef['Data'].shape,dtype=float)
    # # # Add FN communities    
    # # gdf=gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='FN_COMMUNITY_LOCATIONS_SP')
    # # gdf['geometry']=gdf.geometry.buffer(2*1000)    
    # # gdf=gdf[gdf.geometry!=None]
    # # gdf=gdf.reset_index()
    # # gdf['ID']=np.ones(len(gdf))
    # # shapes=((geom,value) for geom, value in zip(gdf['geometry'],gdf['ID']))    
    # # burned=features.rasterize(shapes=shapes,fill=0,out=zMaskWM,transform=zRef['Transform'])        
    # # ind=np.where( (zMaskWM>0) )        
    # # zMaskWM[ind]=1    
    # # # Add major cities
    # # gdf=gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='BC_MAJOR_CITIES_POINTS_500M')
    # # gdf['geometry']=gdf.geometry.buffer(2*1000)
    # # gdf=gdf[gdf.geometry!=None]
    # # gdf=gdf.reset_index()
    # # gdf['ID']=np.ones(len(gdf))
    # # shapes=((geom,value) for geom, value in zip(gdf['geometry'],gdf['ID']))    
    # # burned=features.rasterize(shapes=shapes,fill=0,out=zMaskWM,transform=zRef['Transform'])        
    # # # Add to LU map
    # # ListLU=[lut_lu['Timber Production Passive'],lut_lu['No Designation']]
    # # ind=np.where( (zLC['Data']==lut_lc['Forest']) & (np.isin(zLU['Data'],ListLU)==True) & (zSlope<30) & (zMaskWM==1) & (zNoise2==1) )
    # # zLU['Data'][ind]=lut_lu['Wildfire Risk Management']

    # #--------------------------------------------------------------------------
    # # Add Bands of intense forestry zone
    # gdf=gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='BC_MAJOR_CITIES_POINTS_500M')
    # gdf['geometry']=gdf.geometry.buffer(7*1000)
    # gdf=gdf[gdf.geometry!=None]
    # gdf=gdf.reset_index()
    # gdf['ID']=np.ones(len(gdf))
    # shapes=((geom,value) for geom, value in zip(gdf['geometry'],gdf['ID']))
    # zInner=np.zeros(zRef['Data'].shape,dtype=float)
    # burned=features.rasterize(shapes=shapes,fill=0,out=zInner,transform=zRef['Transform'])
    
    # gdf=gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='BC_MAJOR_CITIES_POINTS_500M')
    # gdf['geometry']=gdf.geometry.buffer(55*1000)
    # gdf=gdf[gdf.geometry!=None]
    # gdf=gdf.reset_index()
    # gdf['ID']=np.ones(len(gdf))
    # shapes=((geom,value) for geom, value in zip(gdf['geometry'],gdf['ID']))
    # zOuter=np.zeros(zRef['Data'].shape,dtype=float)
    # burned=features.rasterize(shapes=shapes,fill=0,out=zOuter,transform=zRef['Transform'])
    
    # gdf=gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='BC_MAJOR_CITIES_POINTS_500M')
    # gdf['geometry']=gdf.geometry.buffer(56*1000)
    # gdf=gdf[gdf.geometry!=None]
    # gdf=gdf.reset_index()
    # gdf['ID']=np.ones(len(gdf))
    # shapes=((geom,value) for geom, value in zip(gdf['geometry'],gdf['ID']))
    # zOuter2=np.zeros(zRef['Data'].shape,dtype=float)
    # burned=features.rasterize(shapes=shapes,fill=0,out=zOuter2,transform=zRef['Transform'])
    
    # gdf=gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='BC_MAJOR_CITIES_POINTS_500M')
    # gdf['geometry']=gdf.geometry.buffer(6*1000)
    # gdf=gdf[gdf.geometry!=None]
    # gdf=gdf.reset_index()
    # gdf['ID']=np.ones(len(gdf))
    # shapes=((geom,value) for geom, value in zip(gdf['geometry'],gdf['ID']))
    # zInner2=np.zeros(zRef['Data'].shape,dtype=float)
    # burned=features.rasterize(shapes=shapes,fill=0,out=zInner2,transform=zRef['Transform'])
    
    # zCb=np.zeros(zRef['Data'].shape,dtype='int8')
    # ind=np.where( (zF==1) & (zOuter==1) )
    # zCb[ind]=1
    # plt.close('all'); plt.matshow(zCb)
    
    # gdf=gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='GSR_TMBR_PRCSSING_FAC_SV')
    # gdf['geometry']=gdf.geometry.buffer(75*1000)
    # gdf=gdf[gdf.geometry!=None]
    # gdf=gdf.reset_index()
    # gdf['ID']=np.ones(len(gdf))
    # shapes=((geom,value) for geom, value in zip(gdf['geometry'],gdf['ID']))
    # zF=np.zeros(zRef['Data'].shape,dtype=float)
    # burned=features.rasterize(shapes=shapes,fill=0,out=zF,transform=zRef['Transform'])
    
    # zMaskIFM=np.zeros(zRef['Data'].shape,dtype='int8')
    # ind=np.where( (zInner==0) & (zOuter==1) & (zF==1) )
    # zMaskIFM[ind]=1
    # #plt.close('all'); plt.matshow(zMaskIFM)
    
    # zMaskWMRing=np.zeros(zRef['Data'].shape,dtype='int8')
    # ind=np.where( (zInner==0) & (zOuter==0) & (zOuter2==1) )
    # zMaskWMRing[ind]=1
    # ind=np.where( (zInner2==0) & (zInner==1) )
    # zMaskWMRing[ind]=1
    # #plt.close('all'); plt.matshow(zMaskWMRing)
        
    # # Add intensee forest management
    # ListLU=[lut_lu['Timber Production Passive'],lut_lu['No Designation']]
    # ind=np.where( (zLC['Data']==lut_lc['Forest']) & (np.isin(zLU['Data'],ListLU)==True) & (zMaskIFM==1) & (zNoise2==1) )
    # zLU['Data'][ind]=lut_lu['Timber Production Intense']
    
    # # # Add outter protective WM ring
    # # ListLU=[lut_lu['Timber Production Passive'],lut_lu['No Designation']]
    # # ind=np.where( (zLC['Data']==lut_lc['Forest']) & (np.isin(zLU['Data'],ListLU)==True) & (zMaskWMRing==1) )
    # # zLU['Data'][ind]=lut_lu['Wildfire Risk Management']
    
    # # Save    
    # gis.SaveGeoTiff(zLU,meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandUse_Comp1_2049_Scn3.tif') 
   
    return

#%% Derive deforestation mask from land use compilation 1
def DeriveLandCoverLandUseComp1_DeforstationMask(meta):    
    zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])        
    lut_lc=meta['LUT']['Derived']['lc_comp1']
    vList=['lc_comp1_1800','lc_comp1_2019','lc_comp1_2049s1','lc_comp1_2049s2']
    z0=u1ha.Import_Raster(meta,[],vList,'Extract Grid')
    z1=copy.deepcopy(zRef); 
    z1['Data']=np.zeros(zRef['Data'].shape,dtype='int8')
    ind=np.where( (z0['lc_comp1_2049s1']!=lut_lc['Forest']) & (z0['lc_comp1_1800']==lut_lc['Forest']) | (z0['lc_comp1_2049s2']!=lut_lc['Forest']) & (z0['lc_comp1_1800']==lut_lc['Forest']) )
    z1['Data'][ind]=1
    # plt.matshow(z1['Data'])
    gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\Disturbances\\LandCover_Comp1_DeforestationMaskAll.tif')
    return

#%%
def ExportSummaryLCLU(meta):
    
    zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
    
    lut_lc=meta['LUT']['Derived']['lc_comp1']
    lut_lu=meta['LUT']['Derived']['lu_comp1']
    
    vList=['lc_comp1_1800','lc_comp1_2019','lu_comp1_2019'] 
    z=u1ha.Import_Raster(meta,[],vList)
    for k in z.keys(): 
        z[k]=z[k]['Data']
    
    bc=np.array(list(lut_lc.values()))
    bu=np.array(list(lut_lu.values()))
    
    # Matrix of 1800 land cover / 2019 land cover
    M=np.zeros((bc.size,bc.size) )
    for i in range(bc.size):
        for j in range(bc.size):
            ind=np.where( (zRef['Data']==1) & (z['lc_comp1_1800']==bc[i]) & (z['lc_comp1_2019']==bu[j]) )        
            if ind[0].size>0:
                M[i,j]=ind[0].size    
    df=pd.DataFrame(M,columns=lut_lc.keys())
    df['LC']=lut_lc.keys()
    df.to_excel(meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCover1800_LandCover2019.xlsx')
    
    # Matrix of 2019 land cover / 2019 land use
    M=np.zeros((bu.size,bc.size) )
    for i in range(bu.size):
        for j in range(bc.size):
            ind=np.where( (z['lc_comp1_2019']==bc[j]) & (z['lu_comp1_2019']==bu[i]) )
            if ind[0].size>0:
                M[i,j]=ind[0].size    
    df=pd.DataFrame(M,columns=lut_lc.keys())
    df['LC']=lut_lu.keys()
    df.to_excel(meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCover2019_LandUse2019.xlsx')
    
    # Matrix of 1800 land cover / 2019 land use
    M=np.zeros((bu.size,bc.size) )
    for i in range(bu.size):
        for j in range(bc.size):
            ind=np.where( (z['lc_comp1_1800']==bc[j]) & (z['lu_comp1_2019']==bu[i]) )        
            if ind[0].size>0:
                M[i,j]=ind[0].size    
    df=pd.DataFrame(M,columns=lut_lc.keys())
    df['LC']=lut_lu.keys()
    df.to_excel(meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCover1800_LandUse2019.xlsx')    
    
    return

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
	id=meta['Param']['BE']['InsectComp1']['ID']
	name=meta['Param']['BE']['InsectComp1']['Insect Name']
	sev=meta['Param']['BE']['InsectComp1']['Severity Code']
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
def DeriveLandUseChangeYear_1800to2019(meta):
    
    zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
    vList=['lc_comp1_1800','lc_comp1_2019','lu_comp1_2019','lc_cec_2010','lc_cec_2020']
    z0=u1ha.Import_Raster(meta,[],vList,'Extract Grid')
    lut_lc=meta['LUT']['Derived']['lc_comp1']
    lut_lu=meta['LUT']['Derived']['lu_comp1']
    lut_cec=meta['LUT']['Derived']['lc_cec_c']
    lut_chng=meta['LUT']['Derived']['lclu_chng_comp1']
    
    cL1=['Cropland','Pasture','Residential and Commercial','Transportation','Energy and Mines','Water Management','No Designation']
    cL2=['FL-CL','FL-PA','FL-RC','FL-TR','FL-EM','FL-WM','FL-ND']
    
    A_For_1800=np.where(z0['lc_comp1_1800']==lut_lc['Forest'])[0].size/1e6
    A_For_2019=np.where(z0['lc_comp1_2019']==lut_lc['Forest'])[0].size/1e6
    A_For_Delta=A_For_1800-A_For_2019
    
    # Create rings from populated locations
    df=gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='BC_MAJOR_CITIES_POINTS_500M')
    df=df[df['CITY_TYPE']!='DM']
    df=df.reset_index()
    df['Seed Year']=1865*np.ones(len(df['BCMJ_TAG']))
    dSet=gu.ReadExcel(r'C:\Data\Cities\BC_MAJOR_CITIES_POINTS_500M.xlsx')    
    for i in range(len(dSet['BCMJ_TAG'])):
        ind=np.where(df['BCMJ_TAG']==dSet['BCMJ_TAG'][i])[0]; 
        df['Seed Year'][ind]=dSet['Year Settled'][i]
    df=df[df['Seed Year']>0]
    df=df[df.geometry!=None]
    df=df.reset_index()
    
    tv=np.arange(1851,2020,1)
    binA=9000*(1/(1+(np.exp(-0.16*(tv-1945)))))
    binGA=np.diff(binA)
    plt.close('all'); plt.plot(tv,binA,'bo')
    plt.close('all'); plt.plot(tv[1:],binGA,'bo')
    binR=np.sqrt(binA/np.pi) #A=pi*r^2
    binRG=np.diff(binR)
    plt.close('all'); plt.plot(tv,binR,'bo')
    #plt.close('all'); plt.plot(binRG,'ko')   
    
    zR=np.zeros(zRef['Data'].shape,dtype='float')
    zY=np.zeros(zRef['Data'].shape,dtype='int16')
    for iT in range(tv.size):
        print(binR[iT])
        ind=np.where(df['Seed Year']<=tv[iT])[0]
        df0=df.loc[ind].copy()
        df0['geometry']=df0.geometry.buffer(1000*binR[iT])
        z2=np.zeros(zRef['Data'].shape,dtype=float)
        shapes=((geom,value) for geom, value in zip(df0['geometry'],df0['BCMJ_TAG']))
        burned=features.rasterize(shapes=shapes,fill=0,out=z2,transform=zRef['Transform'])
        ind2=np.where( (burned>0) & (zR==0) )
        zR[ind2]=binR[iT]
        zY[ind2]=tv[iT]
    # indNZ=np.where(zY>0)
    # idx=gu.IndicesFromUniqueArrayValues(zY[indNZ])
    
    # #plt.close('all'); plt.matshow(zY,clim=[1850,2010])       
    # A=np.zeros(tv.size)
    # for i in idx.keys():
    #     iT=np.where(tv==i)[0]
    #     A[iT]=idx[i].size
    # plt.close('all'); plt.plot(tv,A,'bo')
    
    # LUC over 1800-2019    
    zYear=copy.deepcopy(zRef); zYear['Data']=np.zeros(zRef['Data'].shape,dtype='int16')
    zType=copy.deepcopy(zRef); zType['Data']=np.zeros(zRef['Data'].shape,dtype='int16')
    for iT in range(tv.size):        
        print(tv[iT])
        ind=np.where( (zY==tv[iT]) & (z0['lc_comp1_1800']==lut_lc['Forest']) & (z0['lu_comp1_2019']==lut_lu['Cropland']) )
        zYear['Data'][ind]=tv[iT]; zType['Data'][ind]=lut_chng['FL-CL']
        ind=np.where( (zY==tv[iT]) & (z0['lc_comp1_1800']==lut_lc['Forest']) & (z0['lu_comp1_2019']==lut_lu['Pasture']) )
        zYear['Data'][ind]=tv[iT]; zType['Data'][ind]=lut_chng['FL-PA']
        ind=np.where( (zY==tv[iT]) & (z0['lc_comp1_1800']==lut_lc['Forest']) & (z0['lu_comp1_2019']==lut_lu['Residential and Commercial']) )
        zYear['Data'][ind]=tv[iT]; zType['Data'][ind]=lut_chng['FL-RC']
        ind=np.where( (zY==tv[iT]) & (z0['lc_comp1_1800']==lut_lc['Forest']) & (z0['lu_comp1_2019']==lut_lu['Transportation']) )
        zYear['Data'][ind]=tv[iT]; zType['Data'][ind]=lut_chng['FL-TR']
        ind=np.where( (zY==tv[iT]) & (z0['lc_comp1_1800']==lut_lc['Forest']) & (z0['lu_comp1_2019']==lut_lu['Energy and Mines']) )
        zYear['Data'][ind]=tv[iT]; zType['Data'][ind]=lut_chng['FL-EM']
        ind=np.where( (zY==tv[iT]) & (z0['lc_comp1_1800']==lut_lc['Forest']) & (z0['lu_comp1_2019']==lut_lu['Water Management']) )
        zYear['Data'][ind]=tv[iT]; zType['Data'][ind]=lut_chng['FL-WM']
        ind=np.where( (zY==tv[iT]) & (z0['lc_comp1_1800']==lut_lc['Forest']) & (z0['lu_comp1_2019']==lut_lu['No Designation']) )
        zYear['Data'][ind]=tv[iT]; zType['Data'][ind]=lut_chng['FL-ND']
    A_For_Delta_Prime=np.where(zYear['Data']>0)[0].size/1e6
    zYearBK=copy.deepcopy(zYear)
    zTypeBK=copy.deepcopy(zType)
    #zYear=copy.deepcopy(zYearBK)
    #zType=copy.deepcopy(zTypeBK)    
    
    # Fill unclassified areas, but not areas with deforestaion over 2010-2019
    # Assume random uniform distribution over post-war period    
    for j in range(len(cL1)):
        ind=np.where( (zYear['Data']==0) &                  
                 (z0['lc_comp1_1800']==lut_lc['Forest']) & 
                 (z0['lc_comp1_2019']!=lut_lc['Forest']) &
                 (z0['lu_comp1_2019']==lut_lu[cL1[j]]) )
        zYear['Data'][ind]=np.random.randint(1955,2020,ind[0].size)
        zType['Data'][ind]=zType['Data'][ind]=lut_chng[cL2[j]]
    A_For_Delta_Prime=np.where(zYear['Data']>0)[0].size/1e6    
    #(zDef1019_Year['Data']==0) &
    
    idx=gu.IndicesFromUniqueArrayValues(zYear['Data'][zYear['Data']>0])
    tv=np.arange(1851,2020,1)
    A=np.zeros(tv.size)
    for i in idx.keys():
        iT=np.where(tv==i)[0]
        A[iT]=idx[i].size
    plt.close('all'); plt.plot(tv,A/1000,'-bo')
    
    # Deforestation indicated over 2010-2019
    zDef1019_Year=copy.deepcopy(zRef); zDef1019_Year['Data']=np.zeros(zRef['Data'].shape,dtype='int16')
    zDef1019_Type=copy.deepcopy(zRef); zDef1019_Type['Data']=np.zeros(zRef['Data'].shape,dtype='int16')
    ind=np.where( (z0['lc_cec_2010']==meta['LUT']['Derived']['lc_cec_c']['Forest']) & 
                 (np.isin(z0['lu_comp1_2019'],[lut_lu['Cropland'],lut_lu['Pasture'],lut_lu['Residential and Commercial'],lut_lu['Energy and Mines'],lut_lu['Water Management'],lut_lu['Transportation'],lut_lu['No Designation']])==True) & 
                 (z0['lc_comp1_1800']==lut_lc['Forest']) & 
                 (z0['lc_comp1_2019']!=lut_lc['Forest']) )
    zDef1019_Year['Data'][ind]=np.random.randint(2010,2019,ind[0].size)
    zDef1019_Type['Data'][ind]=zType['Data'][ind]
    A_Def1019=np.where(zDef1019_Year['Data']>0)[0].size/1e6
    
    print(A_For_Delta)
    print(A_For_Delta_Prime)
    print(A_Def1019)    
    
    # Add year of transition to conservation
    zYear=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Disturbances\\LandUseChange_Comp1_1800to2019_Year.tif')
    ind=np.where(z0['lu_comp1_2019']==lut_lu['Conservation Natural']); zYear['Data'][ind]=np.random.randint(1975,high=2019,size=ind[0].size)
    ind=np.where(z0['lu_comp1_2019']==lut_lu['Conservation Consistent']); zYear['Data'][ind]=np.random.randint(1975,high=2019,size=ind[0].size)
    
    # Save    
    gis.SaveGeoTiff(zYear,meta['Paths']['bc1ha'] + '\\Disturbances\\LandUseChange_Comp1_1800to2019_Year.tif')
    gis.SaveGeoTiff(zType,meta['Paths']['bc1ha'] + '\\Disturbances\\LandUseChange_Comp1_1800to2019_Type.tif')
    gis.SaveGeoTiff(zDef1019_Year,meta['Paths']['bc1ha'] + '\\Disturbances\\LandUseChange_Comp1_2010to2019_Year.tif')
    gis.SaveGeoTiff(zDef1019_Type,meta['Paths']['bc1ha'] + '\\Disturbances\\LandUseChange_Comp1_2010to2019_Type.tif')
        
    return
