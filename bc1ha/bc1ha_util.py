
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
import fiona
import time
import cv2
import fcgadgets.macgyver.util_general as gu
import fcgadgets.macgyver.util_gis as gis
import fcgadgets.macgyver.util_query_gdb as qgdb
import fcexplore.psp.Processing.psp_util as ugp
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
    meta['Paths']['bc1ha']=r'C:\Users\rhember\Documents\Data\BC1ha'
    meta['Paths']['bc1ha Ref Grid']=meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandMask.tif'
    meta['Paths']['GDB']={}
    meta['Paths']['GDB']['GDB']=r'C:\Users\rhember\Documents\Data\Geodatabases'
    meta['Paths']['GDB']['LandCover']=r'C:\Users\rhember\Documents\Data\Geodatabases\LandCover\20230607\LandCover.gdb'
    meta['Paths']['GDB']['LandUse']=r'C:\Users\rhember\Documents\Data\Geodatabases\LandUse\20230501\LandUse.gdb'
    meta['Paths']['GDB']['Disturbance']=r'C:\Users\rhember\Documents\Data\Geodatabases\Disturbances\20230501\Disturbances.gdb'
    meta['Paths']['GDB']['Results']=r'C:\Users\rhember\Documents\Data\Geodatabases\Results\20230430\Results.gdb'
    meta['Paths']['GDB']['VRI']=r'C:\Users\rhember\Documents\Data\Geodatabases\VRI\20230401\VRI.gdb'
    meta['Paths']['Model']={}
    meta['Paths']['Model']['Code']=r'C:\Users\rhember\Documents\Code_Python\fcgadgets\cbrunner'
    meta['Paths']['Model']['Parameters']=meta['Paths']['Model']['Code'] + '\\Parameters'
    meta['Paths']['Model']['Taz Datasets']=r'C:\Users\rhember\Documents\Data\Taz Datasets'
    meta['Paths']['GP']={}
    meta['Paths']['GP']['DB']=r'C:\Users\rhember\Documents\Data\GroundPlots\PSP-NADB2'

    meta['Graphics']={'Plot Style':{},'Map':{},'Flowchart':{}}
    meta['Graphics']['Plot Style']='Web' # Manuscript
    meta['Graphics']['gp']=gu.SetGraphics(meta['Graphics']['Plot Style'])
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

    meta['Graphics']['Flowchart']={}
    meta['Graphics']['Flowchart']['Font Name']='Arial'
    meta['Graphics']['Flowchart']['Font Color']='#6b7d2a'
    meta['Graphics']['Flowchart']['Font Size']='10'
    meta['Graphics']['Flowchart']['Node Background Color']='#f0fca2'
    meta['Graphics']['Flowchart']['Cluster Background Color']='#f4f5f2'
    meta['Graphics']['Flowchart']['Cluster Background Color 2']='#d1d1d1'
    
    meta['Graphics']['GP Comp']={'cl':np.array([[0.55,0.75,1],[0.85,1,0.65],[0.27,0.49,0.77],[0.45,0.85,0],[0.8,0.6,0.4],[0.4,0.3,0.2]])}    

    # Tracking parameters
    meta['Graphics']['Fig Count']=1
    meta['Graphics']['Tab Count']=1

    # Initiate geospatial info
    meta['Geos']={}

    # Import variable info
    meta['Geos']['Variable Info']=gu.ReadExcel(meta['Paths']['Model']['Code'] + '\\Parameters\\Parameters_BCFCS_BC1haRasterVariableList.xlsx')

    # Import coordinate reference system
    gdf_bm=gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='TA_REGIONAL_DISTRICTS_SVW')
    meta['Geos']['crs']=gdf_bm.crs

    # Import LUTs
    if 'LUT' not in meta:
        meta=ImportLUTs(meta)

    if 'Include GPs' in argv:
        meta['GP'],data=ugp.ImportPlotData(meta,type='Just Parameters')
        meta['GP']['Data']=data

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
        meta['LUT'][uL[iL]]=gu.ipickle(meta['Paths']['bc1ha'] + '\\LUTs\\LUTs_' + uL[iL] + '.pkl')

    # Override pest severity so that it is in order
    meta['LUT']['PEST_INFESTATION_POLY']['PEST_SEVERITY_CODE']={'T':1,'L':2,'M':3,'S':4,'V':5,'G':6}

    # Raw and Derived layer
    lNam='Derived'
    meta['LUT'][lNam]={}
    meta['LUT']['Raw']={}

    # Burn severity class compelation 1
    d=gu.ReadExcel(meta['Paths']['bc1ha'] + '\\LUTs\\LUT_burnsev_comp1.xlsx')
    vNam='burnsev_comp1'
    meta['LUT']['Raw'][vNam]=d
    meta['LUT'][lNam][vNam]={}
    for i in range(d['Value'].size):
        meta['LUT'][lNam][vNam][ d['Name'][i] ]=d['Value'][i]

    # Land cover - level 2
    vNam='lc2'
    meta['LUT'][lNam][vNam]={}
    meta['LUT'][lNam][vNam]['Water']=1
    meta['LUT'][lNam][vNam]['Land']=2
    meta['LUT'][lNam][vNam]['Non-treed']=3
    meta['LUT'][lNam][vNam]['Treed']=4

    # Land Cover Class 1
    d=gu.ReadExcel(meta['Paths']['bc1ha'] + '\\LUTs\\LUT_lcc1.xlsx')
    vNam='lcc1'
    meta['LUT']['Raw'][vNam]=d
    meta['LUT'][lNam][vNam]={}
    for i in range(d['Value'].size):
        meta['LUT'][lNam][vNam][ d['Name'][i] ]=d['Value'][i]

    # Land Cover Class - NTEMS
    d=gu.ReadExcel(meta['Paths']['bc1ha'] + '\\LUTs\\LUT_lcc_NTEMS.xlsx')
    vNam='lcc_ntem'
    meta['LUT']['Raw'][vNam]=d
    meta['LUT'][lNam][vNam]={}
    for i in range(d['Value'].size):
        meta['LUT'][lNam][vNam][ d['Name'][i] ]=d['Value'][i]

    # Land cover 2020 (CEC Given)
    d=gu.ReadExcel(meta['Paths']['bc1ha'] + '\\LUTs\\LUT_lcc_cec.xlsx')
    vNam='lcc_cec'
    meta['LUT']['Raw'][vNam]=d
    meta['LUT'][lNam][vNam]={}
    for i in range(d['Value'].size):
        meta['LUT'][lNam][vNam][ d['Name'][i] ]=d['Value'][i]

    # Land cover 2020 (CEC Compressed)
    d=gu.ReadExcel(meta['Paths']['bc1ha'] + '\\LUTs\\LUT_lcc_cec_Compressed.xlsx')
    vNam='lcc_cec_c'
    meta['LUT']['Raw'][vNam]=d
    meta['LUT'][lNam][vNam]={}
    for i in range(d['Value'].size):
        meta['LUT'][lNam][vNam][ d['Name'][i] ]=d['Value'][i]

    # Species leading - NTEMS
    d=gu.ReadExcel(meta['Paths']['bc1ha'] + '\\LUTs\\LUT_spc1_NTEMS.xlsx')
    vNam='spc1_ntems'
    meta['LUT']['Raw'][vNam]=d
    meta['LUT'][lNam][vNam]={}
    for i in range(d['Value'].size):
        meta['LUT'][lNam][vNam][ d['Name'][i] ]=d['Value'][i]

    # BGC Zone / NDT Zone Combo
    d=gu.ReadExcel(meta['Paths']['bc1ha'] + '\\LUTs\\LUT_bgc_zone.xlsx')
    vNam='bgc_zone'
    meta['LUT']['Raw'][vNam]=d

    # BGC Zone / NDT Zone Combo
    d=gu.ReadExcel(meta['Paths']['bc1ha'] + '\\LUTs\\LUT_bgcz_ndt_combo.xlsx')
    vNam='bgc-ndt'
    meta['LUT']['Raw'][vNam]=d
    meta['LUT'][lNam][vNam]={}
    for i in range(d['ID'].size):
        meta['LUT'][lNam][vNam][ d['BGC-NDT'][i] ]=d['ID'][i]

    # Tree Density Class
    vNam='tdc'
    meta['LUT'][lNam][vNam]={'Sparse':1,'Open':2,'Dense':3}

    # Reserves
    d=gu.ReadExcel(meta['Paths']['bc1ha'] + '\\LUTs\\LUT_HarvestRegenType.xlsx')
    vNam='HarvestRegenType'
    meta['LUT']['Raw'][vNam]=d
    meta['LUT'][lNam][vNam]={}
    for i in range(d['ID'].size):
        meta['LUT'][lNam][vNam][ d['Name'][i] ]=d['ID'][i]

    d=gu.ReadExcel(meta['Paths']['bc1ha'] + '\\LUTs\\LUT_regentype.xlsx')
    vNam='RegenType'
    meta['LUT']['Raw'][vNam]=d
    meta['LUT'][lNam][vNam]={}
    for i in range(d['ID'].size):
        meta['LUT'][lNam][vNam][ d['Name'][i] ]=d['ID'][i]

    #d=gu.ReadExcel(meta['Paths']['bc1ha'] + '\\LUTs\\LUT_EcozoneCanada.xlsx')
    vNam='ezcan'
    meta['LUT']['Raw'][vNam]=d
    meta['LUT'][lNam][vNam]={}
    for i in range(d['ID'].size):
        meta['LUT'][lNam][vNam][ d['Name'][i] ]=d['ID'][i]

    # Forest Cover Stocking Type
    #d=gu.ReadExcel(r'C:\Users\rhember\Documents\Data\BC1ha\Results\stocktype.tif.vat.xlsx')
    #d=gu.ReadExcel(r'C:\Users\rhember\Documents\Code_Python\fcgadgets\cbrunner\Parameters\Parameters_ForestCover_StockingType.xlsx')
    # nam='fcst'
    # meta['LUT'][nam]={}
    # for i in range(d['VALUE'].size):
    #     meta['LUT'][nam][ d['STOCKING_T'][i] ]=d['VALUE'][i]

    return meta

#%% Build LUTs for variables with categorical data from source geodatabases
# Takes 7.5 hours
#fiona.listlayers(r'C:\\Users\\rhember\\Documents\\Data\\Geodatabases\\LandUse\\20230501\\LandUse.gdb')
#fiona.listlayers(r'C:\Users\rhember\Documents\Data\Geodatabases\VRI\20230401\VRI.gdb')

def BuildLUTsFromSourceGDBs(meta):

    # Unique layers
    uL=np.unique(meta['Geos']['Variable Info']['Layer Name'])

    d={}
    for iL in range(uL.size):

        if uL[iL]=='FTEN_CUT_BLOCK_POLY_SVW': # 'BEC_NATURAL_DISTURBANCE_SV':
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
    else:
        shapes=((geom,value) for geom, value in zip(df['geometry'],df[vNam]))

    z0=np.zeros(zRef['Data'].shape,dtype=float)
    burned=features.rasterize(shapes=shapes,fill=0,out=z0,transform=zRef['Transform'])

    z1=zRef.copy()
    z1['Data']=z0.astype(meta['Geos']['Variable Info']['Precision'][ind[0]])
    gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\' + lNam + '\\' + vNam + '.tif')

    return

#%% Get list of variable labels from GDB

def GetVariablesFromGDB(meta,lNam):
    ind=np.where( (meta['Geos']['Variable Info']['Layer Name']==lNam) )[0]
    if ind.size>1:
        ind=ind[0]
    pthin=meta['Paths']['GDB'][ meta['Geos']['Variable Info']['Geodatabase Name'][ind] ] # fiona.listlayers(pthin)
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
    zFID=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\vri_feaid.tif')
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
        z=zRef.copy()
        z['Data']=np.reshape(d[vNam],zRef['Data'].shape)
        gis.SaveGeoTiff(z,meta['Paths']['bc1ha'] + '\\VRI ' + str(year) + '\\' + vNam + '.tif')

    print((time.time()-t0)/60/60)

    return

#%% Region of interest

def DefineROI(meta,roi,gdf):

    if roi['Type']=='Prov':
        meta['Graphics']['Map']['Legend X']=0.7
        meta['Graphics']['Map']['Show Lakes']='Off'
        meta['Graphics']['Map']['Show Rivers']='Off'

    # Import TSA grid
    tsa=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\FADM_TSA\\TSA_NUMBER_DESCRIPTION.tif')
    tsa=gis.UpdateGridCellsize(tsa,meta['Graphics']['Map']['RGSF'])

    roi['crs']=gdf['bc_bound']['gdf'].crs

    if roi['Type']=='ByTSA':

        # Index to TSAs in query
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

        # Clip mask
        roi['grd']=gis.ClipRasterByXYLimits(roi['grd'],xlim,ylim)
        gc.collect()

    elif roi['Type']=='ByRegDis':

        df=gdf['regdis']['gdf'][np.isin(gdf['regdis']['gdf']['REGIONAL_DISTRICT_NAME'],roi['List'])==True]
        df['ID']=1.0
        shapes=((geom,value) for geom, value in zip(df.geometry,df.ID))
        z0=np.zeros(tsa['Data'].shape,dtype=float)
        burned=features.rasterize(shapes=shapes,fill=0,out=z0,transform=tsa['Transform'])

        roi['grd']=tsa.copy()
        roi['grd']['Data']=burned.astype('int8')
        ind=np.where(burned>0)
        xlim=[np.min(tsa['X'][ind])-5000,np.max(tsa['X'][ind])+5000]
        ylim=[np.min(tsa['Y'][ind])-5000,np.max(tsa['Y'][ind])+5000]

        # Clip mask
        roi['grd']=gis.ClipRasterByXYLimits(roi['grd'],xlim,ylim)
        gc.collect()

    elif roi['Type']=='ByLatLon':

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

        roi['grd']=gis.ClipRasterByXYLimits(tsa,xlim,ylim)
        roi['grd']['Data']=0*roi['grd']['Data']+1
        roi['grd']['yxrat']=np.diff(ylim)[0]/np.diff(xlim)[0]

    elif roi['Type']=='Prov':
        roi['grd']=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
        roi['grd']=gis.UpdateGridCellsize(roi['grd'],meta['Graphics']['Map']['RGSF'])

    # Vector layers
    roi['gdf']={}

    #box(W, S, E, N)
    geom=box(roi['grd']['Extent'][0],roi['grd']['Extent'][2],roi['grd']['Extent'][1],roi['grd']['Extent'][3])
    roi['gdf']['bound']=gpd.GeoDataFrame({"id":1,"geometry":[geom]})

    #roi['gdf']['bound within']=gpd.overlay(gdf['tsa']['gdf'],roi['gdf']['bound'],how='intersection')
    roi['gdf']['tsa']=gdf['tsa']['gdf'].cx[roi['grd']['xmin']:roi['grd']['xmax'],roi['grd']['ymin']:roi['grd']['ymax']]
    roi['gdf']['tsa']=roi['gdf']['tsa'].reset_index(drop=True)
    roi['gdf']['tsa']=gpd.sjoin(roi['gdf']['tsa'],roi['gdf']['bound'],how='left')
    roi['gdf']['tsa']=gpd.overlay(roi['gdf']['tsa'],roi['gdf']['bound'],how='intersection')

    if roi['Type']=='ByTSA':
        roi['gdf']['bound within']=gdf['tsa']['gdf'].iloc[np.isin(gdf['tsa']['gdf'].Name,roi['List'])]
    elif roi['Type']=='ByRegDis':
        roi['gdf']['bound within']=gdf['regdis']['gdf'].iloc[ np.isin(gdf['regdis']['gdf']['REGIONAL_DISTRICT_NAME'],roi['List']) ]

    #roi['gdf']['lakes']=gpd.overlay(gdf['bc_bound']['gdf'][(gdf['bc_bound']['gdf']['TAG']=='lake')],roi['gdf']['bound'],how='intersection')

    roi['gdf']['lakes']=gdf['lakes']['gdf'].cx[roi['grd']['xmin']:roi['grd']['xmax'],roi['grd']['ymin']:roi['grd']['ymax']]
    roi['gdf']['lakes']=roi['gdf']['lakes'].reset_index(drop=True)

    roi['gdf']['rivers']=gdf['rivers']['gdf'].cx[roi['grd']['xmin']:roi['grd']['xmax'],roi['grd']['ymin']:roi['grd']['ymax']]
    roi['gdf']['rivers']=roi['gdf']['rivers'].reset_index(drop=True)

    roi['gdf']['tpf']=gdf['tpf']['gdf'].cx[roi['grd']['xmin']:roi['grd']['xmax'],roi['grd']['ymin']:roi['grd']['ymax']]
    roi['gdf']['tpf']=roi['gdf']['tpf'].reset_index(drop=True)

    roi['gdf']['fnc']=gdf['fnc']['gdf'].cx[roi['grd']['xmin']:roi['grd']['xmax'],roi['grd']['ymin']:roi['grd']['ymax']]
    roi['gdf']['fnc']=roi['gdf']['fnc'].reset_index(drop=True)

    roi['gdf']['bc_bound']=gdf['bc_bound']['gdf'].cx[roi['grd']['xmin']:roi['grd']['xmax'],roi['grd']['ymin']:roi['grd']['ymax']]
    roi['gdf']['bc_bound']=roi['gdf']['bc_bound'].reset_index(drop=True)

    roi['gdf']['popp']=gdf['popp']['gdf'].cx[roi['grd']['xmin']:roi['grd']['xmax'],roi['grd']['ymin']:roi['grd']['ymax']]
    roi['gdf']['popp']=roi['gdf']['popp'].reset_index(drop=True)

    if (roi['Type']=='ByTSA') | (roi['Type']=='ByRegDis'):
        roi['gdf']['lakes']=gpd.overlay(roi['gdf']['lakes'],roi['gdf']['bound within'],how='intersection')
        roi['gdf']['rivers']=gpd.overlay(roi['gdf']['rivers'],roi['gdf']['bound within'],how='intersection')
        roi['gdf']['tpf']=gpd.overlay(roi['gdf']['tpf'],roi['gdf']['bound within'],how='intersection')
        roi['gdf']['fnc']=gpd.overlay(roi['gdf']['fnc'],roi['gdf']['bound within'],how='intersection')
        roi['gdf']['popp']=gpd.overlay(roi['gdf']['popp'],roi['gdf']['bound within'],how='intersection')
        roi['gdf']['bc_bound']=gpd.overlay(roi['gdf']['bc_bound'],roi['gdf']['bound within'],how='intersection')

    try:
        roi['gdf']['road']=gdf['road']['gdf'].cx[roi['grd']['xmin']:roi['grd']['xmax'],roi['grd']['ymin']:roi['grd']['ymax']]
        roi['gdf']['road']=roi['gdf']['road'].reset_index(drop=True)
        if (roi['Type']=='ByTSA') | (roi['Type']=='ByRegDis'):
            roi['gdf']['road']=gpd.overlay(roi['gdf']['road'],roi['gdf']['bound within'],how='intersection')
    except:
        pass

    return roi

#%% Import province-wide basemaps

def Import_GDBs_ProvinceWide(meta):

    flg=0
    if flg==1:
        fiona.listlayers(meta['Paths']['GDB']['LandCover'])

    gdf={}

    # Regional districts
    gdf['regdis']={}
    gdf['regdis']['gdf']=gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='TA_REGIONAL_DISTRICTS_SVW')

    # Politial boundary
    gdf['bc_bound']={}
    gdf['bc_bound']['gdf']=gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='NRC_POLITICAL_BOUNDARIES_1M_SP')

    gdf['rivers']={}
    #gdf['rivers']['gdf']=gpd.read_file(meta['Paths']['GDB']['LandCover'],layer='NRC_RIVERS_1M_SP') # FWA_RIVERS_POLY too slow!
    gdf['rivers']['gdf']=gpd.read_file(meta['Paths']['GDB']['GDB'] + '\\rivers.geojson') 

    # This takes a long time
    gdf['lakes']={}
    #gdf['lakes']['gdf']=gpd.read_file(meta['Paths']['GDB']['LandCover'],layer='NRC_WATERBODIES_1M_SP') # FWA_LAKES_POLY too slow
    gdf['lakes']['gdf']=gpd.read_file(meta['Paths']['GDB']['GDB'] + '\\lakes.geojson') 

    # Import TSA info
    gdf['tsa']={}
    gdf['tsa']['gdf']=gpd.read_file(r'C:\Users\rhember\Documents\Data\Geodatabases\LandUse\tsa.geojson')

    # Roads
    gdf['road']={}
    #gdf['road']['gdf']=gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='MOT_ROAD_FEATURES_INVNTRY_SP')
    gdf['road']['gdf']=gpd.read_file(meta['Paths']['GDB']['GDB'] + '\\road.geojson') 

    # Timber processing facilities
    gdf['tpf']={}
    gdf['tpf']['gdf']=gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='GSR_TMBR_PRCSSING_FAC_SV')

    # First nations communities
    gdf['fnc']={}
    gdf['fnc']['gdf']=gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='FN_COMMUNITY_LOCATIONS_SP')

    # Oil and gas facilities
    gdf['ogf']={}
    #gdf['ogf']['gdf']=gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='DRP_OIL_GAS_FACILITIES_BC_SP')
    #gdf['ogf']['gdf']=gpd.read_file(meta['Paths']['GDB']['GDB'] + '\\ogf.geojson') 

    # Oil and gas pipeline
    gdf['ogp']={}
    #gdf['ogp']['gdf']=gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='DRP_OIL_GAS_PIPELINES_BC_SP')
    #gdf['ogp']['gdf']=gpd.read_file(meta['Paths']['GDB']['GDB'] + '\\ogp.geojson') 

    # Populated places
    gdf['popp']={}
    gdf['popp']['gdf']=gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='NRC_POPULATED_PLACES_1M_SP')

    # Cities
    gdf['cities']=gis.ImportCities(r'C:\Users\rhember\Documents\Data\Cities\Cities.xlsx','GDF')
    gdf['cities'].crs=meta['Geos']['crs']

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

def Import_GDB_Over_ROI(meta_bc1ha,roi,vList):

    for nam in vList:

        if nam=='op':
            roi['gdf'][nam]={}
            roi['gdf'][nam]['Path']=r'C:\Users\rhember\Documents\Data\ForestInventory\Results\20220422\Results.gdb'
            roi['gdf'][nam]['Layer']='RSLT_OPENING_SVW'; # fiona.listlayers(op['Path'])
            roi['gdf'][nam]['crs']=roi['crs']
            roi['gdf'][nam]['Keep Geom']='On'
            roi['gdf'][nam]['Select Openings']=np.array([])
            roi['gdf'][nam]['SBC']=np.array([])
            roi['gdf'][nam]['FSC']=np.array([])
            roi['gdf'][nam]['gdf']=qgdb.Query_Openings(roi['gdf'][nam],roi)

        elif nam=='fcinv':
            roi['gdf'][nam]={}
            roi['gdf'][nam]['Path']=r'C:\Users\rhember\Documents\Data\ForestInventory\Results\20220422\Results.gdb'
            roi['gdf'][nam]['Layer']='RSLT_FOREST_COVER_INV_SVW'; # fiona.listlayers(fcinv['Path'])
            roi['gdf'][nam]['crs']=roi['crs']
            roi['gdf'][nam]['Keep Geom']='On'
            roi['gdf'][nam]['Select Openings']=np.array([])
            roi['gdf'][nam]['SBC']=np.array([])
            roi['gdf'][nam]['FSC']=np.array([])
            roi['gdf'][nam]['gdf']=qgdb.Query_Openings(roi['gdf'][nam],roi)

        elif nam=='fcres':
            roi['gdf'][nam]={}
            roi['gdf'][nam]['Path']=r'C:\Users\rhember\Documents\Data\ForestInventory\Results\20220422\Results.gdb'
            roi['gdf'][nam]['Layer']='RSLT_FOREST_COVER_RESERVE_SVW'; # fiona.listlayers(fcres['Path'])
            roi['gdf'][nam]['crs']=roi['crs']
            roi['gdf'][nam]['Keep Geom']='On'
            roi['gdf'][nam]['Select Openings']=np.array([])
            roi['gdf'][nam]['SBC']=np.array([])
            roi['gdf'][nam]['FSC']=np.array([])
            roi['gdf'][nam]['gdf']=qgdb.Query_Openings(roi['gdf'][nam],roi)

        # Import atu within ROI
        elif nam=='pl':
            roi['gdf'][nam]={}
            roi['gdf'][nam]['Path']=r'C:\Users\rhember\Documents\Data\ForestInventory\Results\20220422\Results.gdb'
            roi['gdf'][nam]['Layer']='RSLT_ACTIVITY_TREATMENT_SVW'; # fiona.listlayers(atu['Path'])
            roi['gdf'][nam]['crs']=roi['crs']
            roi['gdf'][nam]['Keep Geom']='On'
            roi['gdf'][nam]['Select Openings']=np.array([])
            roi['gdf'][nam]['SBC']=np.array(['PL'])
            roi['gdf'][nam]['FSC']=np.array([])
            roi['gdf'][nam]['gdf']=qgdb.Query_Openings(roi['gdf'][nam],roi)

        elif nam=='wf':
            roi['gdf'][nam]={}
            roi['gdf'][nam]['Path']=r'C:\Users\rhember\Documents\Data\ForestInventory\Disturbances\20220422\Disturbances.gdb'
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
            roi['gdf'][nam]['Path']=r'C:\Users\rhember\Documents\Data\ForestInventory\Disturbances\20220422\Disturbances.gdb'
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
            roi['gdf'][nam]['Path']=r'C:\Users\rhember\Documents\Data\ForestInventory\VRI\20220404\VRI.gdb'
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
            roi['gdf'][nam]['Path']=r'C:\Users\rhember\Documents\Data\ForestInventory\LandUse\20220422\LandUse.gdb'
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
        if v=='age_ntem':
            d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Age\\Age_NTEM_2019.tif')
        elif v=='BEC_ZONE_CODE':
            d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\BEC_ZONE_CODE.tif')
        elif v=='bgcz':
            d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\BEC_BIOGEOCLIMATIC_POLY\\ZONE_GapFilled.tif')
        elif v=='biomass_glob':
            d[v]=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\Biomass\GlobBiomass\N80W140_agb pc8.tif')
            d[v]['Data']=np.squeeze(d[v]['Data'])
            d[v]['Data']=0.5*d[v]['Data']
        elif v=='bsr_yr':
            d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VEG_BURN_SEVERITY_SP\\FIRE_YEAR.tif')
        elif v=='bsr_sc':
            d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VEG_BURN_SEVERITY_SP\\BURN_SEVERITY_RATING.tif')
        elif v=='btm':
            d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\landuse.btm.tif')
            d[v]['Data']=np.squeeze(d[v]['Data'])
            lut=pd.read_csv(r'C:\Users\rhember\Documents\Data\BC1ha\LandCoverUse\landuse_btm_category_metadata.csv')
            cl=np.column_stack( (lut['C1'].values,lut['C2'].values,lut['C3'].values) )
            d[v]['Compressed']={}
            d[v]['Compressed']['Data'],d[v]['Compressed']['lab'],d[v]['Compressed']['cl1']=gis.CompressCats(d[v]['Data'],lut['Raster Value'].values,lut['PLU Label'].values,cl)
        elif v=='bsr':
            d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Disturbances\\VEG_BURN_SEVERITY_SP_2017.tif')
            d[v]['key']=gu.ReadExcel(meta['Paths']['bc1ha'] + '\\Disturbances\\VEG_BURN_SEVERITY_SP.xlsx')
        elif v=='crownc':
            d[v]=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\VRI\crownc.tif')
        elif v=='d2road':
            d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Terrain\\DistanceFromRoads.tif')
        elif v=='d2fac':
            d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Management\\DistanceFromForestryFacility.tif')
        elif v=='elev':
            d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Terrain\\elevation.tif')
        elif v=='ezcan':
            d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Ecozones_Canada\\Ecozones_Canada.tif')
        elif v=='feca_yr':
           d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\FECA_YearLast.tif')
        elif v=='fire_yr':
            d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\PROT_HISTORICAL_FIRE_POLYS_SP\\FIRE_YEAR_YearLast.tif')
        elif v=='fire_2023':
            d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\PROT_HISTORICAL_FIRE_POLYS_SP\\FIRE_YEAR_2023.tif')
        elif v=='gfcly':
            d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Disturbances\\GlobalForestChange_LossYear_2021.tif')
        elif v=='gfcly_filt':
            d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Disturbances\\GlobalForestChange_LossYear_2021_Filtered.tif')
        elif v=='gsoc':
            d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Soil\\gsoc2010_bc1ha.tif')
        elif v=='harv_yr_cc':
            d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VEG_CONSOLIDATED_CUT_BLOCKS_SP\\HARVEST_YEAR_YearLast.tif')
        elif v=='harv_yr_ntem':
            d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Disturbances\\Harvest_NTEM_Year.tif')
        elif v=='harv_yr_con1':
            d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Disturbances\\Harvest_Consol1_Year.tif')
        elif v=='harv_yr_con2':
            d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Disturbances\\Harvest_Consol2_Year.tif')
        elif v=='harv_salv':
            d[v]=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Disturbances\HarvestSalvageMask_FromCruise.tif')
        elif v=='harv_prob':
            d[v]=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Disturbances\HarvestProbability.tif')
        elif v=='kd_yr':
            d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\SP_KD_YearLast.tif')
        elif v=='ibm_yr':
            d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\PEST_INFESTATION_POLY\\PEST_SEVERITY_CODE_IBM_YearLast.tif')
        elif v=='idw_mask':
            d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Disturbances\\IDW_Mask.tif')
        elif v=='lc2':
            d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\BCLCS_LEVEL_2.tif')
        elif v=='lc4':
            d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\BCLCS_LEVEL_4.tif')
        elif v=='lc5':
            d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\BCLCS_LEVEL_5.tif')
        elif v=='lcc1_c':
            d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCoverClass1_Current.tif')
        elif v=='lcc1_pc':
            d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCoverClass1_PreContact.tif')
        elif v=='lcc_cec_2020':
            d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCoverClass_CEC_2020_Compressed.tif')
        elif v=='lcc_cec_2010':
            d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCoverClass_CEC_2010_Compressed.tif')
        elif v=='lcc_ntem_2019':
            d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCoverClass_NTEMS_2019.tif')
        elif v=='luc1_yr':
            d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Disturbances\\LandUseChange1_Year.tif')
        elif v=='luc_def_cec':
            d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Disturbances\\LandUseChange_Deforestation_10to20_CEC.tif')
        elif v=='luc_aff_cec':
            d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Disturbances\\LandUseChange_Afforestation_10to20_CEC.tif')
        elif v=='PROJ_AGE_1':
            d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\PROJ_AGE_1.tif')
        elif v=='pfi_c':
            d[v]=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\ForestInventory\LiDAR\Boundary TSA\From GQ 20230302\BoundaryTSA_PFI_LnTotalWSV_bc1ha.tif')
            d[v]['Data']=np.squeeze(d[v]['Data'])
            d[v]['Data']=0.5*0.5*d[v]['Data']
        elif v=='pdead_cruise':
            d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Disturbances\\HarvestPercentDead_FromCruise.tif')
        elif v=='plam':
            d[v]=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Results\Planting_FromRESULTS_MaskCount.tif')
        elif v=='prcp_ann_n':
            d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Climate\\Seasonal\\BC1ha_prcp_ann_norm_1971to2000_si_hist_v1.tif')
        elif v=='protected':
            d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\PROTECTED_LANDS_DESIGNATION.tif')
        elif v=='rangecon':
            d[v]=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\LandCoverUse\rangeten_consol.tif')
        elif v=='rears':
            d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\REARs_CompPlusProp.tif')
        elif v=='refg':
            d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
        elif v=='regentype':
            d[v]=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Results\Harvest_Regen_Type.tif')
        elif v=='rescon':
            d[v]=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Results\SILV_RESERVE_CODE_Consolidated.tif')
        elif v=='rd':
            d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\TA_REGIONAL_DISTRICTS_SVW\\REGIONAL_DISTRICT_NAME.tif')
        elif v=='si_vri23':
            d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\SITE_INDEX.tif')
        elif v=='sphl_vri23':
            d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI\\sphlive.tif')
        elif v=='sphd_vri23':
            d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI\\sphdead.tif')
        elif v=='soc':
            d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Soil\\soc_tot_forest_Shawetal2018.tif')
        elif v=='spc1_ntems':
            d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Species\\NTEMS_Spc1_Filled_RD_CAPITAL.tif')
        elif v=='tdc':
            d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\TreeDensityClass_Current.tif')
        elif v=='tdc_wgs':
            d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\TreeDensityClass_WithShrubsGrasses.tif')
        elif v=='tmean_ann_n':
            d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Climate\\Seasonal\\BC1ha_tmean_ann_norm_1971to2000_si_hist_v1.tif')
        elif v=='si_spl_ntems':
            d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\SPL\\SI_From_NTEMS_Spc1_Filled_RD_CAPITAL_GF.tif')
        elif v=='si_spl_fd':
            d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\SPL\\Site_Prod_Fd.tif')
        elif v=='si_spl_hw':
            d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\SPL\\Site_Prod_Hw.tif')
        elif v=='wbt':
            d[v]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\FWA_WETLANDS_POLY\\WATERBODY_TYPE.tif')
        elif v=='ws_gs_n':
            tmp=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
            #z_tmp=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Climate\Seasonal\wbm_ws_gs_norm_1971to2000.tif')
            #d[v]=tmp.copy()
            #d[v]['Data']=z_tmp['Data']
            #del z_tmp
            #d[v]['Data']=d[v]['Data'].astype('float')
            #d[v]['cm']=gu.ReadExcel(r'C:\Users\rhember\Documents\Data\Colormaps\colormap_ws.xlsx')
        else:
            pass

        if meta['Graphics']['Map']['RGSF']!=1:
            d[v]=gis.UpdateGridCellsize(d[v],meta['Graphics']['Map']['RGSF'])

        if roi!=[]:
            if 'points' in roi.keys():
                if v==vList[0]:
                    iPoints=gis.GetGridIndexToPoints(d[vList[0]],roi['points']['x'],roi['points']['y'])
                d[v]=d[v]['Data'][iPoints]

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
def GapFill_BGCZ(meta):
    zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
    z=u1ha.Import_Raster(meta,[],['lcc1_c','bgcz'])
    ivl=5
    iGap=np.where( (zRef['Data']==1) & (z['bgcz']['Data']==0) )
    iCal=np.where( (zRef['Data'][0::ivl,0::ivl]==1) & (z['bgcz']['Data'][0::ivl,0::ivl]>0) )
    xy=np.column_stack([zRef['X'][0::ivl,0::ivl][iCal],zRef['Y'][0::ivl,0::ivl][iCal]])
    vals=z['bgcz']['Data'][0::ivl,0::ivl][iCal]
    zFill=griddata(xy,vals,(zRef['X'][iGap],zRef['Y'][iGap]),method='nearest')
    z1=zRef.copy()
    z1['Data']=np.zeros(zRef['Data'].shape,dtype='int16')
    z1['Data']=z['bgcz']['Data']
    z1['Data'][iGap]=zFill
    plt.close('all'); plt.matshow(z1['Data'])
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
    zYearLast=zRef.copy()
    zYearLast['Data']=np.zeros(zRef['Data'].shape,dtype='int16')
    uYear=df[vNam].unique()
    tv=np.arange(np.min(uYear),np.max(uYear),1)
    for iT in range(tv.size):

        df0=df[df[vNam]==tv[iT]].copy()
        shapes=((geom,value) for geom, value in zip(df0.geometry,df0[vNam]))

        z0=np.zeros(zRef['Data'].shape,dtype=float)
        if len(df0)>0:
            burned=features.rasterize(shapes=shapes,fill=0,out=z0,transform=zRef['Transform'])

        z1=zRef.copy()
        z1['Data']=z0.astype(meta['Geos']['Variable Info']['Precision'][ind[0]])
        gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\' + lNam + '\\' + vNam + '_' + str(tv[iT]) + '.tif')

        # Update by year grid
        zYearLast['Data'][burned>0]=tv[iT]

    # Year of last occurrence
    z1=zRef.copy()
    z1['Data']=zYearLast['Data'].astype('int16')
    gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\' + lNam + '\\' + vNam + '_YearLast.tif')

    # Mask of occurrence
    z1=zRef.copy()
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
        z1=zRef.copy()
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

        zYearLast=zRef.copy()
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

            z1=zRef.copy()
            z1['Data']=z0.astype(meta['Geos']['Variable Info']['Precision'][indI[0]])
            gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\' + lNam + '\\' + vNam + '_' + pest + '_' + str(tv[iT]) + '.tif')

            # Update by year grid
            zYearLast['Data'][z0>0]=tv[iT]

        # Year of last occurrence
        z1=zRef.copy()
        z1['Data']=zYearLast['Data'].astype('int16')
        gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\' + lNam + '\\' + vNam + '_' + pest + '_YearLast.tif')

        # Mask of occurrence
        z1=zRef.copy()
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
            z1=zRef.copy()
            z1['Data']=z['Year'][iY+1].astype('int16')
            gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\' + lNam + '\\' + vNam + '_' + pest + '_' + str(iY+1) + '_Year.tif')
            z1=zRef.copy()
            z1['Data']=z['Severity'][iY+1].astype('int16')
            gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\' + lNam + '\\' + vNam + '_' + pest + '_' + str(iY+1) + '_Severity.tif')
    return

#%% Define insect outbreaks
def DefineOutbreaks(meta):
    # pestL=['IBM','IBS','IBB','IBD','IDW','IDL']

    # tv=np.arange(1951,2023,1)

    # for pest in pestL:

    #     zMask=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\PEST_INFESTATION_POLY\\PEST_SEVERITY_CODE_' + pest + '_MaskAll.tif')
    #     iMask=np.where(zMask['Data']==1)

    #     zO=[{}]*3
    #     for iO in range(3):
    #         zO[iO]['Year Start']=np.zeros(zMask['Data'].shape,dtype='int16')
    #         zO[iO]['Year End']=np.zeros(zMask['Data'].shape,dtype='int16')
    #         zO[iO]['Duration']=np.zeros(zMask['Data'].shape,dtype='int16')
    #         zO[iO]['Max Severity']=np.zeros(zMask['Data'].shape,dtype='int16')

    #     OutbreakNum=np.zeros(zMask['Data'].shape,dtype='int16')
    #     Active=np.zeros(zMask['Data'].shape,dtype='int16')
    #     Duration=np.zeros(zMask['Data'].shape,dtype='int16')
    #     for iT in range(tv.size):
    #         print(tv[iT])
    #         zS=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\PEST_INFESTATION_POLY\\PEST_SEVERITY_CODE_' + pest + '_' + str(tv[iT]) + '.tif')
    #         for iO in range(3):
    #             # Onset of outbreak
    #             ind=np.where( (OutbreakNum==iO) & (Active==0) & (zS['Data']>0) & (zS['Data']!=meta['LUT']['PEST_INFESTATION_POLY']['PEST_SEVERITY_CODE']['T']) )
    #             Active[ind]=1
    #             zO[iO]['Year Start'][ind]=tv[iT]
    #             zO[iO]['Max Severity'][ind]=np.maximum(zO[iO]['Max Severity'][ind],zS['Data'][ind])
    #             zO[iO]['Duration'][ind]=zO[iO]['Duration'][ind]+1

    #             # Continuation of outbreak
    #             ind=np.where( (OutbreakNum==iO) & (Active==1) & (zS['Data']>0) & (zS['Data']!=meta['LUT']['PEST_INFESTATION_POLY']['PEST_SEVERITY_CODE']['T']) )
    #             zO[iO]['Max Severity'][ind]=np.maximum(zO[iO]['Max Severity'][ind],zS['Data'][ind])
    #             zO[iO]['Duration'][ind]=zO[iO]['Duration'][ind]+1

    #             # End of outbreak
    #             ind=np.where( (OutbreakNum==iO) & (Active==1) & (zS['Data']==0) | (Active==1) & (zS['Data']!=meta['LUT']['PEST_INFESTATION_POLY']['PEST_SEVERITY_CODE']['T']) )
    #             zO[iO]['Year End'][ind]=tv[iT]-1
    #             OutbreakNum[ind]=OutbreakNum[ind]+1

    #     for i in range(3):
    #         tmp=zO[i]
    #         gu.opickle(meta['Paths']['bc1ha'] + '\\PEST_INFESTATION_POLY\\Outbreak_' + pest + '_' + str(i+1) + '.pkl',tmp)
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
        #ats['gdf'].to_file(r'C:\Users\rhember\Documents\Data\BC1ha\RSLT_ACTIVITY_TREATMENT_SVW\ats.geojson',driver="GeoJSON")
    else:
        ats={}
        ats['gdf']=gpd.read_file(r'C:\Users\rhember\Documents\Data\BC1ha\RSLT_ACTIVITY_TREATMENT_SVW\ats.geojson')

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
        gu.opickle(r'C:\Users\rhember\Documents\Data\BC1ha\RSLT_ACTIVITY_TREATMENT_SVW\at.pkl',at)
    else:
        at=gu.ipickle(r'C:\Users\rhember\Documents\Data\BC1ha\RSLT_ACTIVITY_TREATMENT_SVW\at.pkl')

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
    gu.opickle(r'C:\Users\rhember\Documents\Data\BC1ha\RSLT_ACTIVITY_TREATMENT_SVW\dP.pkl',dP)
    #dP=gu.ipickle(r'C:\Users\rhember\Documents\Data\BC1ha\RSLT_ACTIVITY_TREATMENT_SVW\dP.pkl')

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
        z1=zRef.copy()
        z1['Data']=zPac['Year'][iY+1].astype('int16')
        gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + vNam + '_' + str(iY+1) + '_Year.tif')
        z1=zRef.copy()
        z1['Data']=zPac['ACTIVITY_TREATMENT_UNIT_ID'][iY+1].astype('int32')
        gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + vNam + '_' + str(iY+1) + '_ACTIVITY_TREATMENT_UNIT_ID.tif')
        z1=zRef.copy()
        z1['Data']=zPac['ID_SILV_FUND_SOURCE_CODE'][iY+1].astype('int16')
        gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + vNam + '_' + str(iY+1) + '_SILV_FUND_SOURCE_CODE.tif')
        z1=zRef.copy()
        z1['Data']=zPac['SPH_Planted'][iY+1].astype('int16')
        gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + vNam + '_' + str(iY+1) + '_SPH_Planted.tif')
        z1=zRef.copy()
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
        #ats['gdf'].to_file(r'C:\Users\rhember\Documents\Data\BC1ha\RSLT_ACTIVITY_TREATMENT_SVW\ats_ds.geojson',driver="GeoJSON")
    else:
        ats={}
        ats['gdf']=gpd.read_file(r'C:\Users\rhember\Documents\Data\BC1ha\RSLT_ACTIVITY_TREATMENT_SVW\ats_ds.geojson')

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
        gu.opickle(r'C:\Users\rhember\Documents\Data\BC1ha\RSLT_ACTIVITY_TREATMENT_SVW\at_ds.pkl',at)
    else:
        at=gu.ipickle(r'C:\Users\rhember\Documents\Data\BC1ha\RSLT_ACTIVITY_TREATMENT_SVW\at_ds.pkl')

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
    gu.opickle(r'C:\Users\rhember\Documents\Data\BC1ha\RSLT_ACTIVITY_TREATMENT_SVW\dDS.pkl',dDS)
    #dDS=gu.ipickle(r'C:\Users\rhember\Documents\Data\BC1ha\RSLT_ACTIVITY_TREATMENT_SVW\dDS.pkl')

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
        z1=zRef.copy()
        z1['Data']=zPac['Year'][iY+1].astype('int16')
        gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + vNam + '_' + str(iY+1) + '_Year.tif')
        z1=zRef.copy()
        z1['Data']=zPac['ACTIVITY_TREATMENT_UNIT_ID'][iY+1].astype('int32')
        gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + vNam + '_' + str(iY+1) + '_ACTIVITY_TREATMENT_UNIT_ID.tif')
        z1=zRef.copy()
        z1['Data']=zPac['ID_SILV_FUND_SOURCE_CODE'][iY+1].astype('int16')
        gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + vNam + '_' + str(iY+1) + '_SILV_FUND_SOURCE_CODE.tif')
        z1=zRef.copy()
        z1['Data']=zPac['SPH_Planted'][iY+1].astype('int16')
        gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + vNam + '_' + str(iY+1) + '_SPH_Planted.tif')
        z1=zRef.copy()
        z1['Data']=zPac['ID_SILV_TECHNIQUE_CODE'][iY+1].astype('int16')
        gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + vNam + '_' + str(iY+1) + '_SILV_TECHNIQUE_CODE.tif')
    
    return

#%% Planting (Non-obligation by project type)
# Dabbled in doing it this way and then went back to doing it on the fly

def Planting_NonOb_Mask(meta,zRef):
    zMask=zRef.copy()
    zMask['Data']=np.zeros(zRef['Data'].shape,dtype='int16')
    for iEY in range(6):
        zFSC=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_ALL_' + str(iEY+1) + '_SILV_FUND_SOURCE_CODE.tif')
        ind=np.where( (np.isin(zFSC['Data'],meta['Param']['BE']['FSC']['NO List ID'])==True) )
        zMask['Data'][ind]=1
    gis.SaveGeoTiff(zMask,meta['Paths']['bc1ha'] + '\\Management\\PL_NonOb_MaskAll.tif')
    return

#%%

# For knockdown: ats['SMC']=np.array(['CABLE','GUARD','HARV','MDOWN','PUSH'])
#sbc='FE'; stc='CA'
sbc='SP'; stc='BU'; smc=None; vNam='SP-BU'

def RasterizeActivities(meta,sbc,stc,smc,vNam):
    
    tv=np.arange(1960,2024,1)
    zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
    
    # Start with planting with spatial from RESULTS (takes 15 min)
    flg=1
    if flg==1:
        t0=time.time()
        ats={}
        ats['Path']=meta['Paths']['GDB']['Results']
        ats['Layer']='RSLT_ACTIVITY_TREATMENT_SVW'; # fiona.listlayers(ats['Path'])
        ats['crs']=meta['Geos']['crs']
        ats['Keep Geom']='On'
        ats['Select Openings']=np.array([])
        ats['SBC']=np.array([sbc])
        ats['STC']=np.array([stc])
        ats['SMC']=np.array([smc])
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
        ats['gdf']=CreateIdForCategoricalVariable(meta,'RSLT_ACTIVITY_TREATMENT_SVW','SILV_METHOD_CODE',ats['gdf'])
        ats['gdf']=ats['gdf'][ (ats['gdf']['RESULTS_IND']=='Y') & (ats['gdf']['SILV_METHOD_CODE']!='LAYOT') ].copy()
        ats['gdf']=ats['gdf'][ats['gdf'].geometry!=None]
        ats['gdf']=ats['gdf'].reset_index()
        ats['gdf'].to_file(r'C:\Users\rhember\Documents\Data\BC1ha\RSLT_ACTIVITY_TREATMENT_SVW\ats' + vNam + '.geojson',driver="GeoJSON")
    else:
        ats={}
        ats['gdf']=gpd.read_file(r'C:\Users\rhember\Documents\Data\BC1ha\RSLT_ACTIVITY_TREATMENT_SVW\ats' + vNam + '.geojson')

    # Add areas where FC is artificial
    flg=1
    if flg==1:
        at={}
        at['Path']=meta['Paths']['GDB']['Results']
        at['Layer']='RSLT_ACTIVITY_TREATMENT_SVW'; # fiona.listlayers(at['Path'])
        at['crs']=meta['Geos']['crs']
        at['Keep Geom']='Off'
        at['Select Openings']=np.array([])
        at['SBC']=np.array([sbc])
        at['STC']=np.array([stc])
        at['SMC']=np.array([smc])
        at['FSC']=np.array([])
        at['SOC1']=np.array([])
        at['ROI']=[]
        at['gdf']=qgdb.Query_Openings(at,[])
        at['gdf']=CreateIdForCategoricalVariable(meta,'RSLT_ACTIVITY_TREATMENT_SVW','SILV_FUND_SOURCE_CODE',at['gdf'])
        at['gdf']=CreateIdForCategoricalVariable(meta,'RSLT_ACTIVITY_TREATMENT_SVW','SILV_TECHNIQUE_CODE',at['gdf'])
        at['gdf']=CreateIdForCategoricalVariable(meta,'RSLT_ACTIVITY_TREATMENT_SVW','SILV_METHOD_CODE',at['gdf'])
        # Make sure to remove entries that we know did not occur (planned or layout)
        ikp=np.where(  (at['gdf']['RESULTS_IND']=='Y') & (at['gdf']['SILV_METHOD_CODE']!='LAYOT') )[0]
        for k in at['gdf'].keys():
            at['gdf'][k]=at['gdf'][k][ikp]
        at['gdf']['Year']=np.zeros(at['gdf']['ACTIVITY_TREATMENT_UNIT_ID'].size)
        for i in range(at['gdf']['Year'].size):
            at['gdf']['Year'][i]=int(at['gdf']['ATU_COMPLETION_DATE'][i][0:4])
        gu.opickle(r'C:\Users\rhember\Documents\Data\BC1ha\RSLT_ACTIVITY_TREATMENT_SVW\at_' + vNam + '.pkl',at)
    else:
        at=gu.ipickle(r'C:\Users\rhember\Documents\Data\BC1ha\RSLT_ACTIVITY_TREATMENT_SVW\at' + vNam + '.pkl')
       
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
                       'ID_SILV_FUND_SOURCE_CODE':np.array([]),
                       'ACTIVITY_TREATMENT_UNIT_ID':np.array([]),
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
            dP[Year][iS]['ID_SILV_TECHNIQUE_CODE']=np.append(dP[Year][iS]['ID_SILV_TECHNIQUE_CODE'],STC*np.ones(ind2.size))
            dP[Year][iS]['ID_SILV_METHOD_CODE']=np.append(dP[Year][iS]['ID_SILV_METHOD_CODE'],STC*np.ones(ind2.size))
            print('From VRI')
     
        #print('Missing')
    #gu.opickle(r'C:\Users\rhember\Documents\Data\BC1ha\RSLT_ACTIVITY_TREATMENT_SVW\dP.pkl',dP)
    #dP=gu.ipickle(r'C:\Users\rhember\Documents\Data\BC1ha\RSLT_ACTIVITY_TREATMENT_SVW\dP.pkl')

    # Pack

    # Initialize rasters
    N_Year=6
    vL=['ACTIVITY_TREATMENT_UNIT_ID','ID_SILV_FUND_SOURCE_CODE','ID_SILV_TECHNIQUE_CODE','ID_SILV_METHOD_CODE']
    zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])    
    zPac={'Year':{},'ACTIVITY_TREATMENT_UNIT_ID':{},'ID_SILV_FUND_SOURCE_CODE':{},'ID_SILV_TECHNIQUE_CODE':{},'ID_SILV_METHOD_CODE':{}}
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
        zYr['ID_SILV_TECHNIQUE_CODE'][ iA,iB  ]=dP[tv[iT]][iS]['ID_SILV_TECHNIQUE_CODE']
        zYr['ID_SILV_METHOD_CODE'][ iA,iB  ]=dP[tv[iT]][iS]['ID_SILV_METHOD_CODE']
        
        iS=1
        iA=indOP1[0][dP[tv[iT]][iS]['IndexToGrid']]
        iB=indOP1[1][dP[tv[iT]][iS]['IndexToGrid']]
        zYr['ACTIVITY_TREATMENT_UNIT_ID'][ iA,iB  ]=dP[tv[iT]][iS]['ACTIVITY_TREATMENT_UNIT_ID']
        zYr['ID_SILV_FUND_SOURCE_CODE'][ iA,iB  ]=dP[tv[iT]][iS]['ID_SILV_FUND_SOURCE_CODE']
        zYr['ID_SILV_TECHNIQUE_CODE'][ iA,iB  ]=dP[tv[iT]][iS]['ID_SILV_TECHNIQUE_CODE']
        zYr['ID_SILV_METHOD_CODE'][ iA,iB  ]=dP[tv[iT]][iS]['ID_SILV_METHOD_CODE']
        
        iS=2
        iA=indOP2[0][dP[tv[iT]][iS]['IndexToGrid']]
        iB=indOP2[1][dP[tv[iT]][iS]['IndexToGrid']]
        zYr['ACTIVITY_TREATMENT_UNIT_ID'][ iA,iB  ]=dP[tv[iT]][iS]['ACTIVITY_TREATMENT_UNIT_ID']
        zYr['ID_SILV_FUND_SOURCE_CODE'][ iA,iB  ]=dP[tv[iT]][iS]['ID_SILV_FUND_SOURCE_CODE']
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
        z1=zRef.copy()
        z1['Data']=zPac['Year'][iY+1].astype('int16')
        gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + vNam + '_' + str(iY+1) + '_Year.tif')
        z1=zRef.copy()
        z1['Data']=zPac['ACTIVITY_TREATMENT_UNIT_ID'][iY+1].astype('int32')
        gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + vNam + '_' + str(iY+1) + '_ACTIVITY_TREATMENT_UNIT_ID.tif')
        z1=zRef.copy()
        z1['Data']=zPac['ID_SILV_FUND_SOURCE_CODE'][iY+1].astype('int16')
        gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + vNam + '_' + str(iY+1) + '_SILV_FUND_SOURCE_CODE.tif')        
        z1=zRef.copy()
        z1['Data']=zPac['ID_SILV_TECHNIQUE_CODE'][iY+1].astype('int16')
        gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + vNam + '_' + str(iY+1) + '_SILV_TECHNIQUE_CODE.tif')

    return

#lNam='RSLT_ACTIVITY_TREATMENT_SVW'
#vNam='FE-CA'
#tv,N=TimeSeriesOccurrenceFromPackedEvents(meta,lNam,vNam,6)
#plt.plot(tv,N,'-bo')

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

    z1=zRef.copy()
    z1['Data']=z
    gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\VRI 2023\\bgcz_ndt_combo.tif')
    return

#%%
def DeriveDistanceFromRoads(meta):
    # Open the shapefile
    df=gpd.read_file(r'C:\Users\rhember\Documents\Data\ForestInventory\Infrastructure.gdb',layer='MOT_ROAD_FEATURES_INVNTRY_SP')
    
    # Remove features with no geometry
    df=df[df.geometry!=None]    
    z=np.zeros(zRef['Data'].shape,dtype=np.int16)    
    bwD=5; binD=np.arange(bwD,200,bwD)
    for iD in range(binD.size):
        print(binD[iD])
        df0=df.copy()
        df0['geometry']=df0.geometry.buffer(1000*binD[iD])
        z0=np.zeros(zRef['Data'].shape,dtype=float)
        shapes=((geom,value) for geom, value in zip(df0['geometry'],df0['ROAD_FEATURE_INVNTRY_ID']))
        burned=features.rasterize(shapes=shapes,fill=0,out=z0,transform=zRef['Transform'])
        z[(burned>0) & (z==0)]=binD[iD]    
    #plt.matshow(z)    
    z1=zRef.copy()
    z1['Data']=z.astype(np.int16)
    gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\Terrain\DistanceFromRoads.tif')    
    return

#%%
def DeriveDistanceFromFacility(meta):
    #fiona.listlayers(r'C:\Users\rhember\Documents\Data\Geodatabases\LandCover\20230607\LandCover.gdb')
    df=gpd.read_file(r'C:\Users\rhember\Documents\Data\ForestInventory\Infrastructure.gdb',layer='GSR_TMBR_PRCSSING_FAC_SV')
    
    # Remove features with no geometry
    df=df[df.geometry!=None]
    z=np.zeros(zRef['Data'].shape,dtype=np.int16)
    bwD=5; binD=np.arange(bwD,850,bwD)
    for iD in range(binD.size):
        print(binD[iD])
        df0=df.copy()
        df0['geometry']=df0.geometry.buffer(1000*binD[iD])
        df0['ID']=1
        z0=np.zeros(zRef['Data'].shape,dtype=float)
        shapes=((geom,value) for geom, value in zip(df0['geometry'],df0['ID']))
        burned=features.rasterize(shapes=shapes,fill=0,out=z0,transform=zRef['Transform'])
        z[(burned>0) & (z==0)]=binD[iD]
    #plt.matshow(z)
    z1=zRef.copy()
    z1['Data']=z.astype(np.int16)
    gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\Management\\DistanceFromForestryFacility.tif')
    return

#%% Current-year wildfire
# Current year + preveous year (sometimes missing)
def RasterizeWildfireCurrentYear(meta,yr):
    zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])    
    
    # Current year
    yr=2023
    pthin=r'C:\Users\rhember\Documents\Data\Wildfire\Current Year\prot_current_fire_polys.shp'
    df=gpd.read_file(pthin)
    df=df[df.geometry!=None]
    df=df.reset_index()
    shapes=((geom,value) for geom, value in zip(df['geometry'],df['FIRE_YEAR']))
    z0=np.zeros(zRef['Data'].shape,dtype=float)
    burned=features.rasterize(shapes=shapes,fill=0,out=z0,transform=zRef['Transform'])
    z1=zRef.copy()
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
    
    a=geotable.load(r'C:\Users\rhember\Documents\Data\Wildfire\BC Fire Perimeters 2020-2022.kmz')
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
    
    z1=zRef.copy()
    z1['Data']=np.zeros(zRef['Data'].shape,dtype='int16')
    ind=np.where(burned>0); z1['Data'][ind]=1
    print(np.sum(z1['Data']))
    gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\PROT_HISTORICAL_FIRE_POLYS_SP\\FIRE_YEAR_2022.tif')
    return

#%% NALCMS CEC 2020 - Compress categories to remove tropics
def NALCMS_Compress(lut_in,zRef,zCEC):

    cg=['Forest','Shrubland','Grassland','Lichen-Moss','Wetland','Cropland','Barren Ground','Urban','Water','Snow and Ice']
    d={'Value':np.zeros(len(cg),dtype='uint8'),'Name':np.array(['empty' for _ in range(len(cg))],dtype=object)}
    lut_out={}
    cnt=1
    for k in cg:
        lut_out[k]=cnt
        d['Value'][cnt-1]=cnt
        d['Name'][cnt-1]=k
        cnt=cnt+1

    zCEC10c=zRef.copy()
    zCEC20c=zRef.copy()
    zCEC10c['Data']=np.zeros(zRef['Data'].shape,dtype='uint8')
    zCEC20c['Data']=np.zeros(zRef['Data'].shape,dtype='uint8')

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
def ReprojectAndCompressCECLandUse(meta):
    fin=r'C:\Users\rhember\Documents\Data\Land Cover\NALCMS\NALCMS_lc_2010.tif'
    fout=r'C:\Users\rhember\Documents\Data\Land Cover\NALCMS\NALCMS_lc_2010_bc1ha.tif'
    gis.ReprojectRasterAndClipToRaster(fin,fout,meta['Paths']['bc1ha Ref Grid'],meta['Geos']['crs'])    
    zLCC_CEC10=gis.OpenGeoTiff(fout)
    zLCC_CEC10=gis.ClipToRaster(zLCC_CEC10,zRef)    
    ind=np.where(zRef['Data']==0)
    zLCC_CEC10['Data'][ind]=0 
    fin=r'C:\Users\rhember\Documents\Data\Land Cover\NALCMS\NALCMS_lc_2020.tif'
    fout=r'C:\Users\rhember\Documents\Data\Land Cover\NALCMS\NALCMS_lc_2020_bc1ha.tif'
    gis.ReprojectRasterAndClipToRaster(fin,fout,meta['Paths']['bc1ha Ref Grid'],meta['Geos']['crs'])    
    zLCC_CEC20=gis.OpenGeoTiff(fout)
    zLCC_CEC20=gis.ClipToRaster(zLCC_CEC20,zRef)    
    ind=np.where(zRef['Data']==0)
    zLCC_CEC20['Data'][ind]=0         
    # Compress categories to remove tropics
    #zLCC_CEC10=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCoverClass_CEC_' + str(yr) + '.tif')
    #zLCC_CEC20=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCoverClass_CEC_' + str(yr) + '.tif')
    lut,zLCC_CEC10,zLCC_CEC20=u1ha.NALCMS_Compress(meta['LUT']['Derived']['lcc_cec'],zRef,zLCC_CEC10,zLCC_CEC20)
    # Manually saved LUT to excell
    gis.SaveGeoTiff(zLCC_CEC10,meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCoverClass_CEC_' + str(yr) + '_Compressed.tif')
    gis.SaveGeoTiff(zLCC_CEC20,meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCoverClass_CEC_' + str(yr) + '_Compressed.tif')
    return

#%% Reproject Harvest Year from NTEMS
def ReprojectDataFromNTEMS(meta):
    # Harvest year
    fin=r'C:\Users\rhember\Documents\Data\Harvest\NTEMS\harv85to20.tif'
    a=gis.OpenGeoTiff(fin)
    fout=r'C:\Users\rhember\Documents\Data\BC1ha\Disturbances\Harvest_NTEM_Year.tif'
    gis.ReprojectRasterAndClipToRaster(fin,fout,meta['Paths']['bc1ha Ref Grid'],meta['Geos']['crs'])
    z=gis.OpenGeoTiff(fout)
    z=gis.ClipToRaster(z,zRef)
    ind=np.where(zRef['Data']==0)
    z['Data'][ind]=0
    gis.SaveGeoTiff(z,fout)
    # Land cover 
    fin=r'C:\Users\rhember\Documents\Data\Land Cover\NTEMS\vlce2_2019.tif'
    fout=meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCoverClass_NTEMS_2019.tif'
    gis.ReprojectRasterAndClipToRaster(fin,fout,meta['Paths']['bc1ha Ref Grid'],meta['Geos']['crs'])
    z=gis.OpenGeoTiff(fout)
    z=gis.ClipToRaster(z,zRef)
    ind=np.where(zRef['Data']==0)
    z['Data'][ind]=0
    gis.SaveGeoTiff(z,fout)
    # Age
    fin=r'C:\Users\rhember\Documents\Data\Age\NTEMS\age_ntem_c.tif'
    fout=meta['Paths']['bc1ha'] + '\\Age\\Age_NTEM_2019.tif'
    gis.ReprojectRasterAndClipToRaster(fin,fout,meta['Paths']['bc1ha Ref Grid'],meta['Geos']['crs'])
    return

#%% 
def ClimateStatsByBGCZone(meta):
    # *** Needs updating ***
    zBGC=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\becz.tif')
    zBGC['Data']=zBGC['Data'].flatten()
    lutBGC=gu.ReadExcel(meta['Paths']['bc1ha'] + '\\VRI 2023\\becz_lut.xlsx')
    
    zMAT=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Climate\BC1ha_mat_norm_1971to2000_si_hist_v1.tif')
    zMAT['Data']=zMAT['Data'].flatten().astype(float)/10
    
    zWS=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Climate\BC1ha_ws_gs_norm_1971to2000_comp_hist_v1.tif')
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
    df.to_excel(r'C:\Users\rhember\Documents\Data\BC1ha\Climate\tmp.xlsx')
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
        z1=zRef.copy()
        z1['Data']=np.reshape(z[v],zRef['Data'].shape)
        gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\RSLT_FOREST_COVER_INV_SVW\\' + v + '.tif')
        
    return   

#%% Clip rasters to standard grid
# *** This also compresses files that come out of Arc crazy big ***
def ClipToBC1ha(meta):

    # Forest Cover ID
    fin=r'C:\Users\rhember\Documents\Data\BC1ha\RSLT_FOREST_COVER_INV_SVW\fcid.tif'
    fout=r'C:\Users\rhember\Documents\Data\BC1ha\RSLT_FOREST_COVER_INV_SVW\fcid2.tif'
    gis.ClipToRaster_ByFile(fin,fout,meta['Paths']['bc1ha Ref Grid'])

    # VRI Feature ID
    fin=r'C:\Users\rhember\Documents\Data\BC1ha\VRI 2023\vri_feaid.tif'
    fout=r'C:\Users\rhember\Documents\Data\BC1ha\VRI 2023\vri_feaid.tif'
    gis.ClipToRaster_ByFile(fin,fout,meta['Paths']['bc1ha Ref Grid'])

    # DEM
    fin=r'C:\Users\rhember\Documents\Data\BC1ha\Terrain\elevation.tif'
    fout=r'C:\Users\rhember\Documents\Data\BC1ha\Terrain\elevation.tif'
    gis.ClipToRaster_ByFile(fin,fout,meta['Paths']['bc1ha Ref Grid'])

    # Aspect
    fin=r'C:\Users\rhember\Documents\Data\BC1ha\Terrain\aspect.tif'
    fout=r'C:\Users\rhember\Documents\Data\BC1ha\Terrain\aspect.tif'
    gis.ClipToRaster_ByFile(fin,fout,meta['Paths']['bc1ha Ref Grid'])

    # Slope
    fin=r'C:\Users\rhember\Documents\Data\BC1ha\Terrain\slope.tif'
    fout=r'C:\Users\rhember\Documents\Data\BC1ha\Terrain\slope.tif'
    gis.ClipToRaster_ByFile(fin,fout,meta['Paths']['bc1ha Ref Grid'])

    # BTM
    fin=r'C:\Users\rhember\Documents\Data\BC1ha\LandCoverUse\landuse.btm.tif'
    fout=r'C:\Users\rhember\Documents\Data\BC1ha\LandCoverUse\landuse.btm.tif'
    gis.ClipToRaster_ByFile(fin,fout,meta['Paths']['bc1ha Ref Grid'])
    # Radiation
    fin=r'C:\Users\rhember\Documents\Data\BC1ha\Climate\Seasonal\BC1ha_rswd_gs_norm_1971to2000_si_hist_v1_c.tif'
    fout=r'C:\Users\rhember\Documents\Data\BC1ha\Climate\Seasonal\BC1ha_rswd_gs_norm_1971to2000_si_hist_v1_c.tif'
    gis.ClipToRaster_ByFile(fin,fout,meta['Paths']['bc1ha Ref Grid'])

    sL=['At','Ba','Bl','Cw','Ep','Fd','Hm','Hw','Lt','Lw','Pa','Pl','Pw','Py','Sb','Se','Ss','Sx','Sw','Yc']
    for s in sL:
        fin=r'C:\Users\rhember\Documents\Data\SiteProductivityLayer\Site_Prod_BC_Geotiffs\Site_Prod_' + s + '.tif'
        fout=r'C:\Users\rhember\Documents\Data\BC1ha\SPL\Site_Prod_' + s + '.tif'
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
    gdf.to_file(r'C:\Users\rhember\Documents\Data\Geodatabases\LandUse\tsa.geojson',driver='GeoJSON')

    return

#%%
def Plot_LandCoverClass1(meta,roi):

    if 'lcc1_c' in roi['grd'].keys():
        z0=roi['grd']['lcc1_c']['Data']
    else:
        z0=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCoverClass1_Current.tif')
        z0=gis.UpdateGridCellsize(z0,meta['Graphics']['Map']['RGSF'])['Data']

    lab=list(meta['LUT']['Derived']['lcc1'].keys())

    z1=len(lab)*np.ones(z0.shape,dtype='int8')
    for k in meta['LUT']['Derived']['lcc1'].keys():
        ind=np.where(z0==meta['LUT']['Derived']['lcc1'][k])
        if ind[0].size>0:
            z1[ind]=meta['LUT']['Derived']['lcc1'][k]
        else:
            z1[0,meta['LUT']['Derived']['lcc1'][k]]=meta['LUT']['Derived']['lcc1'][k]
    z1[roi['grd']['Data']==0]=meta['LUT']['Derived']['lcc1'][k]+1

    N_vis=len(lab)
    N_hidden=1
    N_tot=N_vis+N_hidden

    cm=np.vstack( ((0.3,0.6,0,1),(0.65,1,0,1),(1,1,0.5,1),(0.75,0.4,1,1),(0.75,0.5,0.5,1),(0.75,0.75,0.75,1),(0.8,0,0,1),(0.95,0.95,0.95,1),(0.75,0.85,0.95,1),(1,1,1,1)) )
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
    if meta['Graphics']['Print Figures']=='On':
        gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_lcc1_c','png',900)
    plt.show()
    return fig,ax

#%%
def Plot_REARs(meta,roi):
    if 'rears' in roi['grd'].keys():
        z1=roi['grd']['rears']['Data']
    else:
        z1=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\REARs_CompPlusProp.tif')['Data']
        z1=gis.UpdateGridCellsize(z1,meta['Graphics']['Map']['RGSF'])['Data']
    z1[(roi['grd']['Data']==1) & (z1==0) & (roi['grd']['lcc1_c']['Data']==meta['LUT']['Derived']['lcc1']['Forest'])]=3
    z1[(roi['grd']['Data']==1) & (roi['grd']['lcc1_c']['Data']!=meta['LUT']['Derived']['lcc1']['Forest'])]=4
    z1[(roi['grd']['Data']==0) | (roi['grd']['lcc1_c']['Data']==meta['LUT']['Derived']['lcc1']['Water'])]=5

    lab=['Protected','Protected (proposed)','Unprotected forest','Non-forest land']

    N_vis=4
    N_hidden=1
    N_tot=N_vis+N_hidden

    cm=np.vstack( ((0.7,0.6,1,1),(0.25,0.85,0.9,1),(0,0.4,0,1),(0.83,0.86,0.9,1),(1,1,1,1)) )
    cm=matplotlib.colors.ListedColormap(cm)

    plt.close('all')
    fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
    im=ax[0].matshow(z1[0::1,0::1],extent=roi['grd']['Extent'],cmap=cm)
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

    if meta['Graphics']['Print Figures']=='On':
        gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_REARs','png',900)

    return fig,ax

#%%
def Plot_Harvest_Po(meta,roi):
    z=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Disturbances\\HarvestProbability.tif')
    z=gis.UpdateGridCellsize(z,meta['Graphics']['Map']['RGSF'])
    sf=1000
    # Apply scale factor to convert to (%/yr)
    #z['Data']=z['Data'].astype('float')/1000
    bw=0.1; bin=np.arange(0,0.5+bw,bw)
    z1=np.ones( z['Data'].shape)
    for i in range(bin.size):
        ind=np.where(np.abs( z['Data']-bin[i]*sf)<=bw*sf/2)
        if ind[0].size==0:
            z1[ind[0][i],ind[1][i]]=i
        else:
            z1[ind]=i
    ind=np.where( z['Data']>=bin[i]*sf ); z1[ind]=i
    ind=np.where( roi['grd']['Data']!=1 ); z1[ind]=i+1

    lab=["%.2f" % x for x in bin]
    lab=np.append(lab,np.array(['Water']))

    N_vis=bin.size
    N_hidden=1
    N_tot=N_vis+N_hidden

    cm=plt.cm.get_cmap('viridis',bin.size)
    cm=np.vstack( (cm.colors,(1,1,1,1)) )
    cm=matplotlib.colors.ListedColormap(cm)

    plt.close('all'); fig,ax=plt.subplots(1,3,figsize=gu.cm2inch(18,18*0.5*roi['grd']['yxrat']))
    im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
    if meta['Graphics']['Map']['Show Bound Land Mask']=='On':
        roi['gdf']['bc_bound'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
    ax[0].set(position=[0,0,0.5,1],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'],aspect='auto')
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
    pos2=[0.39,0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],0.02,N_vis*meta['Graphics']['Map']['Legend Text Space']]
    ax[1].set(position=pos2)

    # Model
    V_Merch=np.arange(1,1200)
    beta=[0.005,-0.04,400]
    Po=beta[0]*(1/(1+np.exp(beta[1]*(V_Merch-beta[2]))))
    ax[2].plot(V_Merch,Po*100,'k-',linewidth=1,label='Harvest on-the-fly model 1')
    beta=[0.005,-0.04,500]
    Po=beta[0]*(1/(1+np.exp(beta[1]*(V_Merch-beta[2]))))
    ax[2].plot(V_Merch,Po*100,'g--',linewidth=1,label='Harvest on-the-fly model 2')
    ax[2].set(position=[0.55,0.12,0.4,0.86],xticks=np.arange(0,1300,100),xlabel='Merchantable volume (m$^3$ ha$^-$$^1$)',ylabel='Annual probability of harvest (%)',xlim=[0,800],ylim=[0,0.7])
    ax[2].legend(loc='upper left',bbox_to_anchor=(0.06,0.92),frameon=False,facecolor='w')
    ax[2].yaxis.set_ticks_position('both'); ax[2].xaxis.set_ticks_position('both'); ax[2].tick_params(length=meta['Graphics']['gp']['tickl'])

    if meta['Graphics']['Print Figures']=='On':
        gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\taz_ann_prob_harvest','png',900)
    plt.show()
    return

#%%
def Plot_BGC_Zone(meta,roi):

    z0=roi['grd']['bgcz']['Data']
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
    z1[(roi['grd']['Data']==0) | (roi['grd']['lcc1_c']['Data']==meta['LUT']['Derived']['lcc1']['Water'])]=N_vis+2

    lab1=np.append(lab1,[''])

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
    if meta['Graphics']['Print Figures']=='On':
        gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_bgcz','png',900)
    plt.show()
    return fig,ax

#%%
def Plot_SalvageLogging(meta,roi):

    lab=['Salvage (>50% dead)','Salvage (>10% dead)','Harvest (<10% dead)','No harvesting','Non-forest land','Outside']
    lab=lab[0:-1]
    z1=roi['grd']['harv_salv']['Data']
    z1[roi['grd']['Data']==0]=6

    # Number of colours and number of colours excluded from colorbar
    N_vis=4
    N_hidden=2
    N_tot=N_vis+N_hidden

    # Colormap
    cm=np.vstack( ((0.9,0,0,1),(1,0.85,0,1),(0.85,1,0.65,1),(0.6,0.6,0.6,1), (0.93,0.93,0.93,1),(1,1,1,1)) )
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

    cb.ax.set(yticklabels=lab)
    cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
    cb.outline.set_edgecolor('w')
    for i in range(cb_bnd.size):
        ax[1].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)

    pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
    ax[1].set(position=pos2)

    gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_harv_salv','png',900)

    return fig,ax

#%%
def Plot_MAP(meta,roi):

    z0=roi['grd']['prcp_ann_n']['Data']

    bw=200; bin=np.arange(0,2400+bw,bw)

    N_vis=bin.size
    N_hidden=2
    N_tot=N_vis+N_hidden

    z1=N_vis*np.ones(z0.shape)
    for i in range(N_vis):
        ind=np.where(np.abs(z0-bin[i])<=bw/2)
        if ind[0].size>0:
            z1[ind]=i+1
    ind=np.where(z0>=bin[i]); z1[ind]=i+1
    z1[1,1]=i+2
    z1[(roi['grd']['Data']==0) | (roi['grd']['lcc1_c']['Data']==0) | (roi['grd']['lcc1_c']['Data']==meta['LUT']['Derived']['lcc1']['Water'])]=i+3

    for i in range(N_vis):
        z1[0,i]=i+1

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

    gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_prcp_ann_n','png',900)

    return fig,ax

#%%
def Plot_MAT(meta,roi):

    z0=roi['grd']['tmean_ann_n']['Data']

    bw=10; bin=np.arange(-30,120+bw,bw)

    N_vis=bin.size
    N_hidden=2
    N_tot=N_vis+N_hidden

    z1=N_vis*np.ones(z0.shape)
    for i in range(N_vis):
        ind=np.where(np.abs(z0-bin[i])<=bw/2)
        if ind[0].size>0:
            z1[ind]=i+1
    ind=np.where(z0>=bin[i]); z1[ind]=i+1
    z1[1,1]=i+2
    z1[(roi['grd']['Data']==0) | (roi['grd']['lcc1_c']['Data']==0) | (roi['grd']['lcc1_c']['Data']==meta['LUT']['Derived']['lcc1']['Water'])]=i+3

    for i in range(N_vis):
        z1[0,i]=i+1

    lab=['']*(N_tot-1)
    lab[0:N_vis]=np.array(bin/10).astype(str)

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

    gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_tmean_ann_n','png',900)
    plt.show()
    return fig,ax

#%%
def Plot_Fire2023(meta,roi):

    z1=4*np.ones(roi['grd']['Data'].shape,dtype='int8')

    ind=np.where( (roi['grd']['lcc1_c']['Data']==1) ); z1[ind]=1
    ind=np.where( (roi['grd']['lcc1_c']['Data']==1) & (roi['grd']['fire_2023']['Data']>0) ); z1[ind]=2
    ind=np.where( (roi['grd']['lcc1_c']['Data']!=1) ); z1[ind]=3
    ind=np.where( (roi['grd']['lcc1_c']['Data']==meta['LUT']['Derived']['lcc1']['Water']) ); z1[ind]=4

    lab=['Forest land','Forest land (affected)','Non-forest land']

    N_vis=len(lab)
    N_hidden=1
    N_tot=N_vis+N_hidden

    cm=np.vstack( ((0.8,0.8,0.8,1),(0.85,0,0,1),(0.9,0.9,0.9,1),(1,1,1,1)) )
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

    gdf['cities'][gdf['cities']['Territory']=='BC'].plot(ax=ax[0],marker='s',edgecolor=[0,0,0],facecolor=[1,1,0],lw=0.25,markersize=9,alpha=1,zorder=2)
    for x,y,label in zip(gdf['cities'].geometry.x,gdf['cities'].geometry.y,gdf['cities'].Name):
        ax[0].annotate(label,xy=(x,y),xytext=(3,2),textcoords="offset points",color=[0,0,0],fontsize=4)

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
    if meta['Graphics']['Print Figures']=='On':
        gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_fire2023','png',900)

    return fig,ax

#%% Plot major timber producing facilities
def PlotMills(ax,gdf,mtypeL,labels):
    lw=0.75
    for mtype in mtypeL:
        #mtype='PLP'
        ind=np.where( (gdf['PRODUCT_CODE']==mtype) )[0]
        if mtype=='LBR':
            y=gdf.iloc[ind]['EST_AN_CAP_MLN_BOARD_FT']/453
            ms=500*y
            gdf.iloc[ind].plot(ax=ax,marker='o',edgecolor='g',facecolor='g',lw=lw,markersize=ms,alpha=0.35,zorder=2)
        elif (mtype=='PLY') | (mtype=='VNR') | (mtype=='OSB') | (mtype=='PNL'):
            # One PNL (panel) mill WestPine MDF in Quesnell - it is MDF
            y=gdf.iloc[ind]['EST_AN_CAP_MLN_SQ_FT']/885
            ms=300*y
            gdf.iloc[ind].plot(ax=ax,marker='o',edgecolor='c',facecolor='c',lw=lw,markersize=ms,alpha=0.35,zorder=2)
        elif mtype=='PLT':
            y=gdf.iloc[ind]['EST_AN_CAP_000_TONNES']*2
            ms=1.0*y
            gdf.iloc[ind].plot(ax=ax,marker='o',edgecolor='r',facecolor='r',lw=lw,markersize=ms,alpha=0.35,zorder=2)
        elif mtype=='CHP':
            y=gdf.iloc[ind]['EST_AN_CAP_000_BDUS']*2
            ms=0.75*y
            gdf.iloc[ind].plot(ax=ax,marker='o',edgecolor='y',facecolor='y',lw=lw,markersize=ms,alpha=0.35,zorder=2)
        elif mtype=='PLP':
            y=gdf.iloc[ind]['EST_AN_CAP_000_TONNES']*2
            ms=0.75*y
            gdf.iloc[ind].plot(ax=ax,marker='o',edgecolor='k',facecolor='k',lw=lw,markersize=ms,alpha=0.35,zorder=2)
        elif mtype=='LVL':
            # Laminated veneer lumber
            y=gdf.iloc[ind]['EST_AN_CAP_MLN_CUBIC_FT']#/885
            ms=200*y
            gdf.iloc[ind].plot(ax=ax,marker='o',edgecolor='k',facecolor='g',lw=lw,markersize=ms,alpha=0.35,zorder=2)

        if labels=='On':
            for x,y,label in zip(gdf.iloc[ind].geometry.x,gdf.iloc[ind].geometry.y,gdf.iloc[ind].COMPANY_NAME):
                ax.annotate(label,xy=(x,y),xytext=(5,4),textcoords="offset points")
    return ax

#%%
def RasterizeConsolidatedCutblocks(meta,zRef):
    lNam='VEG_CONSOLIDATED_CUT_BLOCKS_SP'
    vNam='HARVEST_YEAR'

    if os.path.exists(meta['Paths']['bc1ha'] + '\\' + lNam)==False:
        os.mkdir(meta['Paths']['bc1ha'] + '\\' + lNam)

    ind=np.where( (meta['Geos']['Variable Info']['Layer Name']==lNam) & (meta['Geos']['Variable Info']['Variable Name']==vNam) )[0]

    pthin=meta['Paths']['GDB'][ meta['Geos']['Variable Info']['Geodatabase Name'][ind[0]] ]

    df=gpd.read_file(pthin,layer=lNam)
    df=df[df.geometry!=None]
    df=df.reset_index()

    zYearLast=zRef.copy()
    zYearLast['Data']=np.zeros(zRef['Data'].shape,dtype='int16')

    uYear=df[vNam].unique()
    tv=np.arange(np.min(uYear),np.max(uYear),1)

    for iT in range(tv.size):

        df0=df[df[vNam]==tv[iT]].copy()
        shapes=((geom,value) for geom, value in zip(df0.geometry,df0[vNam]))

        z0=np.zeros(zRef['Data'].shape,dtype=float)
        if len(df0)>0:
            burned=features.rasterize(shapes=shapes,fill=0,out=z0,transform=zRef['Transform'])

        z1=zRef.copy()
        z1['Data']=z0.astype(meta['Geos']['Variable Info']['Precision'][ind[0]])
        gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\' + lNam + '\\' + vNam + '_' + str(tv[iT]) + '.tif')

        # Update by year grid
        zYearLast['Data'][burned>0]=tv[iT]

    # Year of last occurrence
    z1=zRef.copy()
    z1['Data']=zYearLast['Data'].astype('int16')
    gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\' + lNam + '\\' + vNam + '_YearLast.tif')

    # Mask of occurrence
    z1=zRef.copy()
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
        z1=zRef.copy()
        z1['Data']=z['Year'][iY+1].astype('int16')
        gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\' + lNam + '\\' + vNam + '_' + str(iY+1) + '_Year.tif')
    return

#%%
def RasterizeCruisePercentDead(meta):
    dC=gu.ipickle(r'C:\Users\rhember\Documents\Data\ECAS\Received 2023-04-04\rh_COMP_DATA_CleanCombined.pkl')

    vTM=list(meta['LUT']['FTEN_CUT_BLOCK_POLY_SVW']['TIMBER_MARK'].values())
    kTM=list(meta['LUT']['FTEN_CUT_BLOCK_POLY_SVW']['TIMBER_MARK'].keys())

    zTM=gis.OpenGeoTiff(meta['Paths'] + '\\FTEN_CUT_BLOCK_POLY_SVW\\TIMBER_MARK.tif')

    u=np.unique(dC['PRIMARY_MARK'])
    iZ=gu.IndicesFromUniqueArrayValues(zTM['Data'].flatten())

    z1=zRef.copy()
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

    z1=zRef.copy()
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
    z1=zRef.copy()
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

    pthin=r'C:\Users\rhember\Documents\Data\Ecozones\nef_ca_ter_ecozone_v2_2.geojson'
    df=gpd.read_file(pthin)
    df=df[df.geometry!=None]
    df=df.reset_index()
    df.crs=pyproj.CRS(srs['String']['Geographic'])
    df=df.to_crs({'init':'epsg:3005'})
    # Used to create LUT: df.drop(columns='geometry').to_excel(r'C:\Users\rhember\Documents\Data\Ecozones\table.xlsx')

    shapes=((geom,value) for geom, value in zip(df['geometry'],df['ECOZONE_ID']))
    z0=np.zeros(zRef['Data'].shape,dtype=float)
    burned=features.rasterize(shapes=shapes,fill=0,out=z0,transform=zRef['Transform'])
    z1=zRef.copy()
    z1['Data']=np.zeros(zRef['Data'].shape,dtype='int16')
    ind=np.where( (burned>0) &(zRef['Data']==1) ); z1['Data'][ind]=burned[ind]
    plt.close('all'); plt.matshow(z1['Data'],clim=[0,15])
    gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\Ecozones_Canada\\Ecozones_Canada.tif')
    return

#%%
def RasterizeBTKSpray(meta):
    gdf_spray=gpd.read_file(r'C:\Users\rhember\Documents\Data\Aerial Btk Spray\Processed\btk_spray_comp.geojson')
    tv=np.arange(1950,2021,1)
    for iT in range(tv.size):
        print(tv[iT])
        zOut=zRef.copy()
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
    z=zRef.copy()
    z['Data']=np.zeros(z['Data'].shape)
    for f in finL:
        fin=r'C:\Users\rhember\Documents\Data\Global Forest Change\Hansen_GFC-2021-v1.9_lossyear_' + f + '.tif'
        fout=r'C:\Users\rhember\Documents\Data\Global Forest Change\Hansen_GFC-2021-v1.9_lossyear_' + f + 'p.tif'
        gis.ReprojectRasterAndClipToRaster(fin,fout,meta['Paths']['bc1ha Ref Grid'],meta['Geos']['crs'])

        z0=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\Global Forest Change\Hansen_GFC-2021-v1.9_lossyear_' + f + 'p.tif')
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
def Import_BurnSeverityCanada(meta):
    pthin=r'C:\Users\rhember\Documents\Data\Wildfire\Burn Severity'
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
def DeriveBurnSeverityConsolidated(meta):
    # meta['LUT']['VEG_BURN_SEVERITY_SP']['BURN_SEVERITY_RATING']
    # meta['LUT']['Derived']['burnsev_comp1']
    zYn=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Disturbances\\BurnSev_Can_Year.tif')
    zBn=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Disturbances\\BurnSev_Can_dNBR.tif')
    zYb=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VEG_BURN_SEVERITY_SP\\FIRE_YEAR.tif')
    zBb=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VEG_BURN_SEVERITY_SP\\BURN_SEVERITY_RATING.tif')    
    # Use distribution of BC classes to derive classes from dNBR for Canada (exclude unburned)
    p=gu.CountByCategories(zBb['Data'][zBb['Data']>0],'Percent')
    p_ord=np.array([23.83,19.68,38.49,15.23])
    p_ord=p_ord/np.sum(p_ord)*100 # normalize without unknown class    
    # Get percentile values for national dNBR
    cpn=np.flip(np.percentile(zBn['Data'][zBn['Data']>0],np.cumsum(p_ord)))
    N_Year=6
    for iY in range(1,N_Year):
        zY=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\PROT_HISTORICAL_FIRE_POLYS_SP\\FIRE_YEAR_' + str(iY+1) + '_Year.tif')        
        zBS=zRef.copy()
        zBS['Data']=5*np.ones(zRef['Data'].shape,dtype='int8')
        
        # Add national data (converted from dNBR to severity class)
        uY=np.unique(zYn['Data'][zYn['Data']>0])
        for iU in range(uY.size):        
            for iC in range(cpn.size):
                ind=np.where( (zY['Data']==uY[iU]) & (zYn['Data']==uY[iU]) & (zBn['Data']<cpn[iC]) )
                zBS['Data'][ind]=list(meta['LUT']['Derived']['burnsev_comp1'].values())[iC]
        # Add BC data
        uY=np.unique(zYb['Data'][zYb['Data']>0])
        for iU in range(uY.size):
            for k in meta['LUT']['Derived']['burnsev_comp1'].keys():
                ind=np.where( (zY['Data']==uY[iU]) & (zYb['Data']==uY[iU]) & (zBb['Data']==meta['LUT']['VEG_BURN_SEVERITY_SP']['BURN_SEVERITY_RATING'][k]) )
                zBS['Data'][ind]=meta['LUT']['Derived']['burnsev_comp1'][k]
        #plt.matshow(zBS['Data'],clim=[0,4])
        # Save
        gis.SaveGeoTiff(zBS,meta['Paths']['bc1ha'] + '\\PROT_HISTORICAL_FIRE_POLYS_SP\\FIRE_YEAR_' + str(iY+1) + '_SevClass.tif')
        return

#%%
def ConsolidateProtectedLand(meta):
    zLCC1=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCoverClass1.tif')
    zProt=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\TA_PROTECTED_LANDS_SV\\PROTECTED_LANDS_DESIGNATION.tif')
    zPark=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\TA_PARK_ECORES_PA_SVW\\ADMIN_AREA_SID.tif')
    zOGMA=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RMP_OGMA_LEGAL_ALL_SVW\\LEGAL_OGMA_PROVID.tif')
    zOGD=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\OGSR_TAP_PRIORITY_DEF_AREA_SP\\PRIORITY_DEFERRAL_ID.tif')
    zNP=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\CLAB_NATIONAL_PARKS\\ENGLISH_NAME.tif')
    #zFCR=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\RSLT_FOREST_COVER_RESERVE_SVW.tif')
    
    # ListCon=['No Timber Harvesting Areas','Landscape Corridors',
    #   'Critical Deer Winter Range','Sensitive Watershed','Water Management Units','High Value Wetlands for Moose','Telkwa Caribou Recovery Area',
    #   'Caribou Migration Corridor','High Biodiversity Emphasis Areas','Scenic Corridors']
    
    # Generate random new protected areas
    flg=0
    if flg==1:
        gdf=u1ha.Import_GDBs_ProvinceWide()
        N=2500 # Number within bounding box
        Dbuf=6200 # metres
    
        x=np.random.uniform(zRef['xmin'],zRef['xmax'],size=N)
        y=np.random.uniform(zRef['ymin'],zRef['ymax'],size=N)
        points=[]
        for k in range(x.size):
            points.append(Point(x[k],y[k]))
        gdf_xy=gpd.GeoDataFrame({'geometry':points,'ID':1})
        gdf_xy.crs=gdf['bc_bound']['gdf'].crs
    
        gdf_xy=gpd.sjoin(gdf_xy,gdf['bc_bound']['gdf'],op='within')
        #gdf_xy.plot(ax=ax[0],markersize=8)
    
        gdf_xyb=gdf_xy.geometry.buffer(Dbuf)
        gdf_xyb=gpd.GeoDataFrame({'geometry':gdf_xyb,'ID':np.arange(0,gdf_xyb.size,1)})
        gdf_xyb.crs=gdf['bc_bound']['gdf'].crs
        #gdf_xyb.plot(ax=ax[0],color='r')
    
        shapes=((geom,value) for geom, value in zip(gdf_xyb['geometry'],gdf_xyb['ID']))
        z=np.zeros(zRef['Data'].shape,dtype=float)
        burned=features.rasterize(shapes=shapes,fill=0,out=z,transform=zRef['Transform'])
    
        zRn=zRef.copy()
        zRn['Data']=z.astype('int32')
        gis.SaveGeoTiff(zRn,meta['Paths']['bc1ha'] + '\\LandCoverUse\\REARs_RandomAreas.tif')
    
    else:
        zRn=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\REARs_RandomAreas.tif')
    
    # # Protected areas with random areas
    
    # Area treed
    ind0=np.where( (zLCC1['Data']==meta['LUT']['Derived']['lcc1']['Forest']) );
    A_treed=ind0[0].size
    
    # Everything completed+proposed without random additions (Comp=1, Prop=2)
    zCP=zRef.copy()
    zCP['Data']=np.zeros(zRef['Data'].shape,dtype='int16')
    ind=np.where( (zLCC1['Data']==meta['LUT']['Derived']['lcc1']['Forest']) & (zPark['Data']>0) |
                  (zLCC1['Data']==meta['LUT']['Derived']['lcc1']['Forest']) & (zOGMA['Data']>0) |
                  (zLCC1['Data']==meta['LUT']['Derived']['lcc1']['Forest']) & (zNP['Data']>0) |
                  (zLCC1['Data']==meta['LUT']['Derived']['lcc1']['Forest']) & (zProt['Data']>0))
    zCP['Data'][ind]=1
    ind=np.where( (zCP['Data']==0) & (zLCC1['Data']==meta['LUT']['Derived']['lcc1']['Forest']) & (zOGD['Data']>0) )
    zCP['Data'][ind]=2
    ind=np.where( (zLCC1['Data']==meta['LUT']['Derived']['lcc1']['Forest']) & (zCP['Data']==0) ); zCP['Data'][ind]=3
    ind=np.where( (zLCC1['Data']!=meta['LUT']['Derived']['lcc1']['Forest']) & (zCP['Data']==0) ); zCP['Data'][ind]=4
    ind=np.where( (zRef['Data']==0) ); zCP['Data'][ind]=5
    
    ind=np.where( (zLCC1['Data']==meta['LUT']['Derived']['lcc1']['Forest']) & (zCP['Data']>=1) & (zCP['Data']<=2) ); A_CPsep=ind[0].size
    print('Area Protected from completed+proposed (Sep) (%): ' + str(A_CPsep/A_treed*100))
    
    zCP['Data']=zCP['Data'].astype('int16')
    gis.SaveGeoTiff(zCP,meta['Paths']['bc1ha'] + '\\LandCoverUse\\REARs_CompPlusProp.tif')
    
    # Everything completed+proposed with random additions (Comp=1, Prop=2)
    zCP=zRef.copy()
    zCP['Data']=np.zeros(zRef['Data'].shape,dtype='int16')
    ind=np.where( (zLCC1['Data']==meta['LUT']['Derived']['lcc1']['Forest']) & (zPark['Data']>0) |
                  (zLCC1['Data']==meta['LUT']['Derived']['lcc1']['Forest']) & (zOGMA['Data']>0) |
                  (zLCC1['Data']==meta['LUT']['Derived']['lcc1']['Forest']) & (zNP['Data']>0) |
                  (zLCC1['Data']==meta['LUT']['Derived']['lcc1']['Forest']) & (zProt['Data']>0) )
    zCP['Data'][ind]=1
    ind=np.where( (zCP['Data']==0) & (zLCC1['Data']==meta['LUT']['Derived']['lcc1']['Forest']) & (zOGD['Data']>0) | (zCP['Data']==0) & (zLCC1['Data']==meta['LUT']['Derived']['lcc1']['Forest']) & (zRn['Data']>0) )
    zCP['Data'][ind]=2
    ind=np.where( (zLCC1['Data']==meta['LUT']['Derived']['lcc1']['Forest']) & (zCP['Data']==0) ); zCP['Data'][ind]=3
    ind=np.where( (zLCC1['Data']!=meta['LUT']['Derived']['lcc1']['Forest']) & (zCP['Data']==0) ); zCP['Data'][ind]=4
    ind=np.where( (zRef['Data']==0) ); zCP['Data'][ind]=5
    
    ind=np.where( (zLCC1['Data']==meta['LUT']['Derived']['lcc1']['Forest']) & (zCP['Data']>=1) & (zCP['Data']<=2) ); A_CPsep=ind[0].size
    print('Area Protected from completed+proposed (Sep) (%): ' + str(A_CPsep/A_treed*100))
    
    zCP['Data']=zCP['Data'].astype('int32')
    gis.SaveGeoTiff(zCP,meta['Paths']['bc1ha'] + '\\LandCoverUse\\REARs_CompPlusPropWithRandomAdditions.tif')
    return

#%% Consolidated retention layer
def ConsolidateRetention(meta):
    d=gu.ReadExcel(meta['Paths']['bc1ha'] + '\\Results\\SILV_RESERVE_CODE.xlsx')
    zLCC1=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCoverClass1.tif')
    zH=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Disturbances\\VEG_CONSOLIDATED_CUT_BLOCKS_SP_MaskAll.tif')
    zFCR=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Results\\SILV_RESERVE_CODE.tif')
    zIBR=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Results\\SILV_RESERVE_CODE_InBlock.tif')
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
    
    d1={1:['Dispersed'],2:'Group',3:'Riparian',4:'Wildlife trees',5:'Other',6:'No reserves',7:'Forest',8:'Non-forest'}
    
    z1=zRef.copy()
    z1['Data']=np.zeros(zRef['Data'].shape,dtype='int8')
    ind=np.where( (zC>1) ); z1['Data'][ind]=5
    ind=np.where( (zC==7) | (zC==4) ); z1['Data'][ind]=1
    ind=np.where( (zC==8) ); z1['Data'][ind]=2
    ind=np.where( (zC==9) ); z1['Data'][ind]=3
    ind=np.where( (zC==2) ); z1['Data'][ind]=4
    ind=np.where( (z1['Data']==0) & (zH['Data']>0) ); z1['Data'][ind]=6
    ind=np.where( (z1['Data']==0) & (zH['Data']==0) & (zLCC1['Data']==meta['LUT']['Derived']['lcc1']['Forest']) ); z1['Data'][ind]=7
    ind=np.where( (z1['Data']==0) & (zH['Data']==0) & (zLCC1['Data']!=meta['LUT']['Derived']['lcc1']['Forest']) ); z1['Data'][ind]=8
    gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\Results\\SILV_RESERVE_CODE_Consolidated.tif')
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
    
    z=zRef.copy()
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
def ConsolidateTreeDensityClass(meta):
    z=Import_Raster(meta,[],['refg','lcc1_c','lc5','harv_yr_con1'])
        
    z1=z['refg'].copy()
    z1['Data']=np.zeros(z['lc5']['Data'].shape,dtype='int8')
    ind=np.where( (z['lcc1_c']['Data']==meta['LUT']['Derived']['lcc1']['Forest']) & (z['lc5']['Data']==meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_5']['SP']) ); z1['Data'][ind]=meta['LUT']['Derived']['tdc']['Sparse']
    ind=np.where( (z['lcc1_c']['Data']==meta['LUT']['Derived']['lcc1']['Forest']) & (z['lc5']['Data']==meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_5']['OP']) ); z1['Data'][ind]=meta['LUT']['Derived']['tdc']['Open']
    ind=np.where( (z['lcc1_c']['Data']==meta['LUT']['Derived']['lcc1']['Forest']) & (z['lc5']['Data']==meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_5']['DE']) ); z1['Data'][ind]=meta['LUT']['Derived']['tdc']['Dense']
    ind=np.where( (z['harv_yr_con1']['Data']==1) ); z1['Data'][ind]=meta['LUT']['Derived']['tdc']['Dense']
    ind=np.where( (z['lcc1_c']['Data']==meta['LUT']['Derived']['lcc1']['Forest']) & (z1['Data']==0) ); z1['Data'][ind]=meta['LUT']['Derived']['tdc']['Open']
    # plt.matshow(z1,clim=[0,3])
    gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\LandCoverUse\\TreeDensityClass_Current.tif')
    
    # # With shrubs and grasses
    # z=np.zeros(zLC4['Data'].shape,dtype='int8')
    # ind=np.where( (z['lcc1_c']['Data']==lut['lc4']['Shrub Low']) | (zLC4['Data']==lut['lc4']['Shrub Tall']) ); z[ind]=4
    # ind=np.where( (z['lcc1_c']['Data']==lut['lc4']['Herb Gramanoid']) | (zLC4['Data']==lut['lc4']['Herb Forbs']) ); z[ind]=5
    # ind=np.where( (z['lcc1_c']['Data']==1) & (z['lc5']['Data']==lut['lc5']['Sparse']) ); z[ind]=1
    # ind=np.where( (z['lcc1_c']['Data']==1) & (z['lc5']['Data']==lut['lc5']['Open']) ); z[ind]=2
    # ind=np.where( (z['lcc1_c']['Data']==1) & (z['lc5']['Data']==lut['lc5']['Dense']) ); z[ind]=3
    # ind=np.where( (z['harv_yr_con1']['Data']==1) ); z[ind]=meta['LUT']['Derived']['tdc']['Dense']
    # ind=np.where( (zRef['Data']==1) & (z==0) ); z[ind]=6
    # ind=np.where( (zRef['Data']==0) ); z[ind]=7
    # plt.matshow(z[0::3,0::3])
    
    # z1=zRef.copy()
    # z1['Data']=z
    # gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\LandCoverUse\\TreeDensityClass_WithShrubsGrasses.tif')
    return

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
    zOut=zRef.copy()
    zOut['Data']=np.zeros(zRef['Data'].shape,dtype='int8')
    zOut['Data'][burned>0]=1
    zOut['Data']=zOut['Data'].astype('int8')
    gis.SaveGeoTiff(zOut,meta['Paths']['bc1ha'] + '\\LandCoverUse\LandMask.tif')
    return

#%% Harvest early reconstruction
def DeriveEarlyHarvest(meta):
    zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
    zLCC1=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCoverClass1.tif')
    zAge=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\PROJ_AGE_1.tif')
    zH_CC=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VEG_CONSOLIDATED_CUT_BLOCKS_SP\\HARVEST_YEAR_YearLast.tif')
    zH_NTEM=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Disturbances\\Harvest_NTEM_Year.tif')
    
    z0=np.zeros(zRef['Data'].shape,dtype='int16')
    ind=np.where( (zRef['Data']==1) & (zH_CC['Data']>0) ); z0[ind]=1
    ind=np.where( (zRef['Data']==1) & (z0==0) & (zH_NTEM['Data']>0) ); z0[ind]=2
    ind=np.where( (zRef['Data']==1) & (z0==0) & (zAge['Data']>0) & (zAge['Data']<150) & (zLCC1['Data']==meta['LUT']['Derived']['lcc1']['Forest']) ); z0[ind]=3
    ind=np.where( (zRef['Data']==1) & (z0==0) & (zAge['Data']>=150) & (zLCC1['Data']==meta['LUT']['Derived']['lcc1']['Forest']) ); z0[ind]=4
    ind=np.where( (zRef['Data']==1) & (z0==0) & (zAge['Data']==0) & (zLCC1['Data']==meta['LUT']['Derived']['lcc1']['Forest']) ); z0[ind]=4
    plt.close('all')
    plt.matshow(z0)
    
    # Create rings from populated locations
    df=gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='NRC_POPULATED_PLACES_1M_SP')
    df=df[df.geometry!=None]
    zD=np.zeros(zRef['Data'].shape,dtype='int16')
    #binD=np.arange(1,110,1)
    binD=np.array([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.25,2.5,2.75,3,3.25,3.5,3.75,4,4.5,5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10,10.5,11,11.5,12,12.5,13,13.5,14,14.5,15,16,17,18,19,20,22,24,26,28,30,40,50,60,120])
    for iD in range(binD.size):
        print(binD[iD])
        df0=df.copy()
        df0['geometry']=df0.geometry.buffer(1000*binD[iD])
        z1=np.zeros(zRef['Data'].shape,dtype=float)
        shapes=((geom,value) for geom, value in zip(df0['geometry'],df0['NRC_PP1M_SYSID']))
        burned=features.rasterize(shapes=shapes,fill=0,out=z1,transform=zRef['Transform'])
        zD[(burned>0) & (zD==0)]=binD[iD]
    
    # Year for each ring
    yr=np.linspace(1855,1960,binD.size).astype('int16')
    
    zH=zRef.copy()
    zH['Data']=np.zeros(zRef['Data'].shape,dtype='int16')
    for iD in range(binD.size):
        ind=np.where( (z0==3) & (zD==binD[iD]) )
        zH['Data'][ind]=yr[iD]
    plt.close('all'); plt.matshow(zH['Data'],clim=[1855,1955])
    
    gis.SaveGeoTiff(zH,meta['Paths']['bc1ha'] + '\\Disturbances\\Harvest_Early_Reconstruction_Year.tif')
    return

#%% Rasterize early land use change year
def RasterizeEaryLandUseChangeYear(meta):
    zLCC1=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCoverClass1.tif')
    
    # Create rings from populated locations
    df=gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='NRC_POPULATED_PLACES_1M_SP')
    df=df[df.geometry!=None]
    zD=np.zeros(zRef['Data'].shape,dtype='int16')
    binD=np.array([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.25,2.5,2.75,3,3.25,3.5,3.75,4,4.5,5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10,11,12,13,14,15,20,25,30,60,120,240,480])
    for iD in range(binD.size):
        print(binD[iD])
        df0=df.copy()
        df0['geometry']=df0.geometry.buffer(1000*binD[iD])
        z0=np.zeros(zRef['Data'].shape,dtype=float)
        shapes=((geom,value) for geom, value in zip(df0['geometry'],df0['NRC_PP1M_SYSID']))
        burned=features.rasterize(shapes=shapes,fill=0,out=z0,transform=zRef['Transform'])
        zD[(burned>0) & (zD==0)]=binD[iD]
    
    # Check that the rings are big enough
    ind=np.where( (zLCC1['Data']==meta['LUT']['Derived']['lcc1']['Settlement']) & (zD==0) | \
                 (zLCC1['Data']==meta['LUT']['Derived']['lcc1']['Herb']) & (zD==0))
    print(ind[0].size)
    
    # Year for each ring
    yr=np.linspace(1855,2010,binD.size).astype('int16')
    
    zLUC=zRef.copy()
    zLUC['Data']=np.zeros(zRef['Data'].shape,dtype='int16')
    for iD in range(binD.size):
        ind=np.where( (zLCC1['Data']==meta['LUT']['Derived']['lcc1']['Settlement']) & (zD==binD[iD]) | \
                     (zLCC1['Data']==meta['LUT']['Derived']['lcc1']['Herb']) & (zD==binD[iD]))
        zLUC['Data'][ind]=yr[iD]
    #plt.matshow(zLUC['Data'],clim=[1845,2020])
    
    gis.SaveGeoTiff(zLUC,meta['Paths']['bc1ha'] + '\\Disturbances\\LandUseChange1_Year.tif')
    
    #z=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\Terrain\DistanceFromRoads.tif')
    return

#%% Planting layer (species and genetic worth)
def lantingLayer(meta,zRef):
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
            z=zRef.copy()
            z['Data']=np.reshape(z0['PL_SPECIES_CD' + str(iS+1)],(zRef['Data'].shape))
            gis.SaveGeoTiff(z,meta['Paths']['bc1ha'] + '\\' + lNam + '\\' + vNam + '_' + str(iY+1) + '_PL_SPECIES_CD' + str(iS+1) + '.tif')
            z=zRef.copy()
            z['Data']=np.reshape(z0['PL_SPECIES_PCT' + str(iS+1)],(zRef['Data'].shape))
            gis.SaveGeoTiff(z,meta['Paths']['bc1ha'] + '\\' + lNam + '\\' + vNam + '_' + str(iY+1) + '_PL_SPECIES_PCT' + str(iS+1) + '.tif')
            z=zRef.copy()
            z['Data']=np.reshape(z0['PL_SPECIES_GW' + str(iS+1)],(zRef['Data'].shape))
            gis.SaveGeoTiff(z,meta['Paths']['bc1ha'] + '\\' + lNam + '\\' + vNam + '_' + str(iY+1) + '_PL_SPECIES_GW' + str(iS+1) + '.tif')
    
        print((time.time()-t0)/60)
    return

#%%
def ConsolidateHarvesting(meta):
    zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
    zLCC1=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCoverClass1.tif')
    zH_CC=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VEG_CONSOLIDATED_CUT_BLOCKS_SP\\HARVEST_YEAR_YearLast.tif')
    zH_NTEM=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Disturbances\\Harvest_NTEM_Year.tif')
    zH_Early=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Disturbances\\Harvest_Early_Reconstruction_Year.tif')
    
    # Harvest consolidated 1 (w/o early reconstruction)
    z=zRef.copy()
    z['Data']=np.zeros(zRef['Data'].shape,dtype='int16')
    ind=np.where( (zLCC1['Data']==meta['LUT']['Derived']['lcc1']['Forest']) & (zH_CC['Data']>0) ); z['Data'][ind]=zH_CC['Data'][ind]
    ind=np.where( (zLCC1['Data']==meta['LUT']['Derived']['lcc1']['Forest']) & (zH_NTEM['Data']>0) ); z['Data'][ind]=zH_NTEM['Data'][ind]
    gis.SaveGeoTiff(z,meta['Paths']['bc1ha'] + '\\Disturbances\\Harvest_Consol1_Year.tif')

    # Harvest consolidated 2 (with early reconstruction
    z=zRef.copy()
    z['Data']=np.zeros(zRef['Data'].shape,dtype='int16')
    ind=np.where( (zLCC1['Data']==meta['LUT']['Derived']['lcc1']['Forest']) & (zH_Early['Data']>0) ); z['Data'][ind]=zH_Early['Data'][ind]
    ind=np.where( (zLCC1['Data']==meta['LUT']['Derived']['lcc1']['Forest']) & (zH_CC['Data']>0) ); z['Data'][ind]=zH_CC['Data'][ind]
    ind=np.where( (zLCC1['Data']==meta['LUT']['Derived']['lcc1']['Forest']) & (zH_NTEM['Data']>0) ); z['Data'][ind]=zH_NTEM['Data'][ind]
    gis.SaveGeoTiff(z,meta['Paths']['bc1ha'] + '\\Disturbances\\Harvest_Consol2_Year.tif')
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
    zH=gis.OpenGeoTiff(meta['Paths'] + '\\Disturbances\\VEG_CONSOLIDATED_CUT_BLOCKS_SP_MaskAll.tif')
    zLCC1=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCoverClass1.tif')
    zPD=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Disturbances\\HarvestPercentDead_FromCruise.tif')
    
    z=zRef.copy()
    z['Data']=6*np.ones(zRef['Data'].shape,dtype='int8')
    ind=np.where( (zPD['Data']>=50) ); z['Data'][ind]=1 # Salvage high
    ind=np.where( (zPD['Data']>=10) & (zPD['Data']<50) ); z['Data'][ind]=2 # Salvage low
    ind=np.where( (zPD['Data']<10) & (zH['Data']>0) ); z['Data'][ind]=3 # Non-salvage harvested forest
    ind=np.where( (zH['Data']==0) & (zLCC1['Data']==1) ); z['Data'][ind]=4 # Unharvested forest
    ind=np.where( (zRef['Data']==1) & (zLCC1['Data']!=1) ); z['Data'][ind]=5 # Non-forest land
    ind=np.where( (zRef['Data']==0) ); z['Data'][ind]=6 # Non land
    gis.SaveGeoTiff(z,meta['Paths'] + '\\Disturbances\\HarvestSalvageMask_FromCruise.tif')
    return

#%% Generate sparse inputs
# This should speed import by an order of magnitude

def GenerateSparseInputs(meta,rgsf,mask):
    #rgsf=100
    
    z=Import_Raster(meta,[],['refg','lcc1_c'])
    #zRef=gis.UpdateGridCellsize(z['refg'],rgsf)
    zLCC1=gis.UpdateGridCellsize(z['lcc1_c'],rgsf)
    
    if mask=='All':
        zMask=z['refg']
    elif mask=='NOSE':
        # Non-obligation stand establishment
        zMask=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Management\\PL_NonOb_MaskAll.tif')
    else:
        pass
    zMask=gis.UpdateGridCellsize(zMask,rgsf)
    
    iMask=np.where( (zMask['Data']==1) & (zLCC1['Data']==meta['LUT']['Derived']['lcc1']['Forest']) )
        
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
    
    # Beetles
    N_Year=10
    d=[None]*N_Year
    for iY in range(N_Year):
        z0=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\PEST_INFESTATION_POLY\\PEST_SEVERITY_CODE_IBM_' + str(iY+1) + '_Year.tif')
        z0=gis.UpdateGridCellsize(z0,rgsf)
        d[iY]=z0['Data'][iMask]
    gu.opickle(meta['Paths']['bc1ha'] + '\\Sparse\\RGSF' + str(rgsf) + '_Mask' + mask + '_PEST_SEVERITY_CODE_IBM_Year.pkl',d)
    d=[None]*N_Year
    for iY in range(N_Year):
        z0=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\PEST_INFESTATION_POLY\\PEST_SEVERITY_CODE_IBM_' + str(iY+1) + '_Severity.tif')
        z0=gis.UpdateGridCellsize(z0,rgsf)
        d[iY]=z0['Data'][iMask]
    gu.opickle(meta['Paths']['bc1ha'] + '\\Sparse\\RGSF' + str(rgsf) + '_Mask' + mask + '_PEST_SEVERITY_CODE_IBM_SevClass.pkl',d)
    
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
    for iY in range(N_Year):
        z0=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\FE-CA_' + str(iY+1) + '_Year.tif')
        z0=gis.UpdateGridCellsize(z0,rgsf)
        d[iY]=z0['Data'][iMask]
    gu.opickle(meta['Paths']['bc1ha'] + '\\Sparse\\RGSF' + str(rgsf) + '_Mask' + mask + '_FE-CA_Year.pkl',d)
    
    # Knockdown
    N_Year=3
    d=[None]*N_Year
    for iY in range(N_Year):
        z0=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\SP_KD_' + str(iY+1) + '_Year.tif')
        z0=gis.UpdateGridCellsize(z0,rgsf)
        d[iY]=z0['Data'][iMask]
    gu.opickle(meta['Paths']['bc1ha'] + '\\Sparse\\RGSF' + str(rgsf) + '_Mask' + mask + '_SP_KD_Year.pkl',d)
    
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
    zOP2=zRef.copy()
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
    zD=u1ha.Import_Raster(meta,[],['harv_yr_con1','kd_yr'])
    
    zP={}
    for iP in range(6):
        zP[iP]={}
        zP[iP]['Type']=np.zeros(zRef['Data'].shape,dtype='int8')
        zP[iP]['Year']=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(iP+1) + '_Year.tif')['Data']
        zP[iP]['STC']=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(iP+1) + '_SILV_TECHNIQUE_CODE.tif')['Data']
        zP[iP]['FSC']=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(iP+1) + '_SILV_FUND_SOURCE_CODE.tif')['Data']
    
    # Summarize frequency of each Stocking Type Code (All and NOSE)
    d['STC Summary']={}
    d['STC Summary']['stc']=['PL','RP','FP','CG','RO','RR','SE','SL']
    d['STC Summary']['N']=np.zeros((len(d['STC Summary']['stc']),2))
    for v in d['STC Summary']['stc']:
        for iP in range(6):
            ind=np.where( (zP[iP]['STC']==meta['LUT']['RSLT_ACTIVITY_TREATMENT_SVW']['SILV_TECHNIQUE_CODE'][v]) )
            d['STC Summary']['N'][cnt,0]=d['STC Summary']['N'][cnt,0]+ind[0].size
            ind=np.where( (zP[iP]['STC']==meta['LUT']['RSLT_ACTIVITY_TREATMENT_SVW']['SILV_TECHNIQUE_CODE'][v]) & (np.isin(zP[iP]['FSC'],meta['Param']['BE']['FSC']['NO List ID'])==True) )
            d['STC Summary']['N'][cnt,1]=d['STC Summary']['N'][cnt,1]+ind[0].size
    
    cd=meta['LUT']['RSLT_ACTIVITY_TREATMENT_SVW']['SILV_TECHNIQUE_CODE']
    stcInclude=[cd['PL'],cd['RP'],cd['FP'],cd['RR']]
    
    d['A Tot']=np.zeros(d['tv'].size)
    d['A']=np.zeros((d['tv'].size,ptNam.size))
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
            #ind=np.where( (zP[iP]['Year']==d['tv'][iT]) & (np.isin(zP[iP]['FSC'],meta['Param']['BE']['FSC']['NO List ID'])==True) & (np.isin(zP[iP]['STC'],stcInclude)==True) )
            pl_oc[ind]=1
            pl_stc[ind]=zP[iP]['STC'][ind]
            pl_fsc[ind]=zP[iP]['FSC'][ind]
        
        ind=np.where( (pl_oc==1) )
        d['A Tot'][iT]=ind[0].size
        
        zType0=np.zeros(zRef['Data'].shape,dtype='int8')
        
        # Replanting
        ind=np.where( (pl_oc==1) & (pl_stc==meta['LUT']['RSLT_ACTIVITY_TREATMENT_SVW']['SILV_TECHNIQUE_CODE']['RP']) )
        zType0[ind]=meta['LUT']['Derived']['RegenType']['Replanting']
        d['A'][iT,meta['LUT']['Derived']['RegenType']['Replanting']-1]=ind[0].size
        # Fill Planting
        ind=np.where( (pl_oc==1) & (pl_stc==meta['LUT']['RSLT_ACTIVITY_TREATMENT_SVW']['SILV_TECHNIQUE_CODE']['FP']) )
        zType0[ind]=meta['LUT']['Derived']['RegenType']['Fill Planting']
        d['A'][iT,meta['LUT']['Derived']['RegenType']['Fill Planting']-1]=ind[0].size
        # Road rehab
        ind=np.where( (pl_oc==1) & (pl_stc==meta['LUT']['RSLT_ACTIVITY_TREATMENT_SVW']['SILV_TECHNIQUE_CODE']['RR']) )
        zType0[ind]=meta['LUT']['Derived']['RegenType']['Road Rehabilitation']
        d['A'][iT,meta['LUT']['Derived']['RegenType']['Road Rehabilitation']-1]=ind[0].size
        # Back to back planting
        ind=np.where( (pl_oc==1) & (zType0==0) & (DL==meta['LUT']['Event']['Planting']) )
        zType0[ind]=meta['LUT']['Derived']['RegenType']['Back-to-back Planting']
        d['A'][iT,meta['LUT']['Derived']['RegenType']['Back-to-back Planting']-1]=ind[0].size    
        # Salvage
        ind=np.where( (pl_oc==1) & (d['tv'][iT]>=2004) & (zD['harv_yr_con1']['Data']>1987) & (zType0==0) & (DL==meta['LUT']['Event']['Harvest']) )
        zType0[ind]=meta['LUT']['Derived']['RegenType']['Salvage and Planting']
        d['A'][iT,meta['LUT']['Derived']['RegenType']['Salvage and Planting']-1]=ind[0].size
        # Knockdown
        ind=np.where( (pl_oc==1) & (zType0==0) & (DL==meta['LUT']['Event']['Knockdown']) )
        zType0[ind]=meta['LUT']['Derived']['RegenType']['Knockdown and Planting']
        d['A'][iT,meta['LUT']['Derived']['RegenType']['Knockdown and Planting']-1]=ind[0].size
        # Straight fire
        ind=np.where( (pl_oc==1) & (zType0==0) & (DL==meta['LUT']['Event']['Wildfire']) )
        zType0[ind]=meta['LUT']['Derived']['RegenType']['Straight-to-planting Post Wildfire']
        d['A'][iT,meta['LUT']['Derived']['RegenType']['Straight-to-planting Post Wildfire']-1]=ind[0].size
        # Straight insect
        ind=np.where( (pl_oc==1) & (zType0==0) & (DL==meta['LUT']['Event']['IBM']) )
        zType0[ind]=meta['LUT']['Derived']['RegenType']['Straight-to-planting Post Insect Outbreak']
        d['A'][iT,meta['LUT']['Derived']['RegenType']['Straight-to-planting Post Insect Outbreak']-1]=ind[0].size
        # NSR backlog
        ind=np.where( (pl_oc==1) & (DL==meta['LUT']['Event']['Harvest']) & (zD['harv_yr_con1']['Data']<=1987) & (zType0==0) )
        zType0[ind]=meta['LUT']['Derived']['RegenType']['NSR Backlog']
        d['A'][iT,meta['LUT']['Derived']['RegenType']['NSR Backlog']-1]=ind[0].size        
        # Unknown
        ind=np.where( (pl_oc==1) & (zType0==0) )
        zType0[ind]=meta['LUT']['Derived']['RegenType']['Unknown']
        d['A'][iT,meta['LUT']['Derived']['RegenType']['Unknown']-1]=ind[0].size
        
        DL[(pl_oc==1)]=meta['LUT']['Event']['Planting']
        
        # Pack
        for iP in range(6):
            ind=np.where( (zP[iP]['Year']==d['tv'][iT]) )
            zP[iP]['Type'][ind]=zType0[ind]
    
    # Save summary
    gu.opickle(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\RegenTypeSummary_All',d)
    
    # Save packed regen type
    for i in range(6):
        z1=zRef.copy()
        z1['Data']=zP[i]['Type']
        gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(i+1) + '_RegenType.tif')
    
    return