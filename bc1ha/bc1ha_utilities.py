
#%% Import modules

import os
import numpy as np
import gc
import matplotlib.pyplot as plt
import geopandas as gpd
import pandas as pd
from shapely.geometry import Polygon,Point,box,shape
from rasterio import features
import fiona
import time
import fcgadgets.macgyver.utilities_general as gu
import fcgadgets.macgyver.utilities_gis as gis
import fcgadgets.macgyver.utilities_query_gdb as qgdb

#%% Initialize project

def Init():

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
    meta['Paths']['Model']['Taz Datasets']=r'C:\Users\rhember\Documents\Data\Taz Datasets'

    # Initiate geospatial info
    meta['Geos']={}

    # Import variable info
    meta['Geos']['Variable Info']=gu.ReadExcel(meta['Paths']['Model']['Code'] + '\\Parameters\\BC1ha Raster Variable List.xlsx')

    # Import coordinate reference system
    gdf_bm=gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='TA_REGIONAL_DISTRICTS_SVW')
    meta['Geos']['crs']=gdf_bm.crs

    return meta

#%% Build LUTs for variables with categorical data from source geodatabases
# Takes 7.5 hours
#fiona.listlayers(r'C:\\Users\\rhember\\Documents\\Data\\Geodatabases\\LandUse\\20230501\\LandUse.gdb')
#fiona.listlayers(r'C:\Users\rhember\Documents\Data\Geodatabases\VRI\20230401\VRI.gdb')

def BuildLUTsFromSourceDBs(meta):

    # Unique layers
    uL=np.unique(meta['Geos']['Variable Info']['Layer Name'])

    d={}
    for iL in range(uL.size):

        if uL[iL]=='VEG_COMP_LYR_R1_POLY': # 'BEC_NATURAL_DISTURBANCE_SV':
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
    vNam='lcc_ntems'
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

    # Forest Cover Stocking Type
    #d=gu.ReadExcel(r'C:\Users\rhember\Documents\Data\BC1ha\Results\stocktype.tif.vat.xlsx')
    #d=gu.ReadExcel(r'C:\Users\rhember\Documents\Code_Python\fcgadgets\cbrunner\Parameters\Parameters_ForestCover_StockingType.xlsx')
    # nam='fcst'
    # meta['LUT'][nam]={}
    # for i in range(d['VALUE'].size):
    #     meta['LUT'][nam][ d['STOCKING_T'][i] ]=d['VALUE'][i]

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

#%% Create ID for categorical variable

def CreateIdForCategoricalVariable(meta,lNam,vNam,df):
    v=list(df.keys())[0]
    df['ID_' + vNam]=np.zeros(df[v].size)
    for k in meta['LUT'][lNam][vNam].keys():
        ind=np.where(df[vNam]==k)[0]
        df['ID_' + vNam][ind]=meta['LUT'][lNam][vNam][k]
    return df

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

    # Import TSA grid
    tsa=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\FADM_TSA\\TSA_NUMBER_DESCRIPTION.tif')
    tsa['Data']=np.squeeze(tsa['Data'])

    roi['crs']=gdf['bc_bound']['gdf'].crs

    # Vector:
    roi['gdf']={}

    if roi['Type']=='ByTSA':

        # Index to TSAs in query
        iROI=gdf['tsa']['key'].VALUE[np.isin(gdf['tsa']['key'].Name,roi['List'])].values

        roi['grd']=tsa.copy()
        roi['grd']['Data']=np.zeros(tsa['Data'].shape)
        ind=(np.isin(tsa['Data'],iROI))
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
    gdf['rivers']['gdf']=gpd.read_file(meta['Paths']['GDB']['LandCover'],layer='NRC_RIVERS_1M_SP') # FWA_RIVERS_POLY too slow!

    # This takes a long time
    gdf['lakes']={}
    gdf['lakes']['gdf']=gpd.read_file(meta['Paths']['GDB']['LandCover'],layer='NRC_WATERBODIES_1M_SP') # FWA_LAKES_POLY too slow

    # Import TSA info
    gdf['tsa']={}
    gdf['tsa']['gdf']=gpd.read_file(r'C:\Users\rhember\Documents\Data\Geodatabases\LandUse\tsa.geojson')

    # Roads
    gdf['road']={}
    gdf['road']['gdf']=gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='MOT_ROAD_FEATURES_INVNTRY_SP')
    #gpd.read_file(r'C:\Users\rhember\Documents\Data\Basemaps\roads.shp')

    # Timber processing facilities
    gdf['tpf']={}
    gdf['tpf']['gdf']=gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='GSR_TMBR_PRCSSING_FAC_SV')

    # First nations communities
    gdf['fnc']={}
    gdf['fnc']['gdf']=gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='FN_COMMUNITY_LOCATIONS_SP')

    # Oil and gas facilities
    gdf['ogf']={}
    gdf['ogf']['gdf']=gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='DRP_OIL_GAS_FACILITIES_BC_SP')

    # Oil and gas pipeline
    gdf['ogp']={}
    gdf['ogp']['gdf']=gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='DRP_OIL_GAS_PIPELINES_BC_SP')

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

def Import_Raster_Over_ROI(meta,roi,vList):

    xlim=(roi['grd']['xmin'],roi['grd']['xmax']+roi['grd']['Cellsize'])
    ylim=(roi['grd']['ymin'],roi['grd']['ymax']+roi['grd']['Cellsize'])

    for nam in vList:

        if nam in roi['grd'].keys():
            continue
        if nam=='age_ntem':
            roi['grd'][nam]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Age\\Age_NTEM_2019.tif')
        elif nam=='aff_cec':
            roi['grd'][nam]=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\LandCoverUse\afforestation10to20_FromCEC.tif')
        elif nam=='BEC_ZONE_CODE':
            roi['grd'][nam]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\BEC_ZONE_CODE.tif')
        elif nam=='bgc_zone':
            roi['grd'][nam]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\BEC_BIOGEOCLIMATIC_POLY\\ZONE.tif')
        elif nam=='btm':
            roi['grd'][nam]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\landuse.btm.tif')
            roi['grd'][nam]['Data']=np.squeeze(roi['grd'][nam]['Data'])
            roi['grd'][nam]=gis.ClipToRaster(roi['grd'][nam],roi['grd'])
            lut=pd.read_csv(r'C:\Users\rhember\Documents\Data\BC1ha\LandCoverUse\landuse_btm_category_metadata.csv')
            cl=np.column_stack( (lut['C1'].values,lut['C2'].values,lut['C3'].values) )
            roi['grd'][nam]['Compressed']={}
            roi['grd'][nam]['Compressed']['Data'],roi['grd'][nam]['Compressed']['lab'],roi['grd'][nam]['Compressed']['cl1']=gis.CompressCats(roi['grd'][nam]['Data'],lut['Raster Value'].values,lut['PLU Label'].values,cl)
        elif nam=='crownc':
            roi['grd'][nam]=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\VRI\crownc.tif')
        elif nam=='cut_yr':
            roi['grd'][nam]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VEG_CONSOLIDATED_CUT_BLOCKS_SP\\HARVEST_YEAR_YearLast.tif')
        elif nam=='cut_yr_ntems':
            roi['grd'][nam]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Disturbances\\Harvest_NTEM_Year.tif')
        elif nam=='cut_yr_con1':
            roi['grd'][nam]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Disturbances\\Harvest_Consilidated1_Year.tif')
        elif nam=='denseclass':
            roi['grd'][nam]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\TreeDensityClass_WithShrubsGrasses.tif')
        elif nam=='d2road':
            roi['grd'][nam]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Terrain\\DistanceFromRoads.tif')
        elif nam=='def_cec':
            roi['grd'][nam]=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\LandCoverUse\deforestation10to20_FromCEC.tif')
        elif nam=='d2fac':
            roi['grd'][nam]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Management\\DistanceFromForestryFacility.tif')
        elif nam=='elev':
            roi['grd'][nam]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Terrain\\elevation.tif')
        elif nam=='fire_yr':
            roi['grd'][nam]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\PROT_HISTORICAL_FIRE_POLYS_SP\\FIRE_YEAR_YearLast.tif')
        elif nam=='gfcly':
            roi['grd'][nam]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Disturbances\\GlobalForestChange_LossYear_2021_KnownEventsRemoved.tif')
        elif nam=='globbio':
            roi['grd'][nam]=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\Biomass\GlobBiomass\N80W140_agb pc8.tif')
            roi['grd'][nam]['Data']=np.squeeze(roi['grd'][nam]['Data'])
            roi['grd'][nam]['Data']=0.5*roi['grd'][nam]['Data']
        elif nam=='gsoc':
            roi['grd'][nam]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Soil\\gsoc2010_bc1ha.tif')
            roi['grd'][nam]=gis.ClipToRaster(roi['grd'][nam],roi['grd'])
        elif nam=='idw_mask':
            roi['grd'][nam]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Disturbances\\IDW_Mask.tif')
        elif nam=='lcc1':
            roi['grd'][nam]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCoverClass1.tif')
        elif nam=='lc2':
            roi['grd'][nam]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI\\lc2.tif')
        elif nam=='lc20_cec':
            roi['grd'][nam]=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\LandCoverUse\LandCoverClass_CEC_2020_Compressed.tif')
        elif nam=='PROJ_AGE_1':
            roi['grd'][nam]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\PROJ_AGE_1.tif')
            roi['grd'][nam]=gis.ClipToRaster(roi['grd'][nam],roi['grd'])
        elif nam=='rears':
            roi['grd'][nam]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\REARs.tif')
        elif nam=='SITE_INDEX':
            roi['grd'][nam]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\SITE_INDEX.tif')
        elif nam=='sphlive':
            roi['grd'][nam]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI\\sphlive.tif')
        elif nam=='sphdead':
            roi['grd'][nam]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI\\sphdead.tif')
        elif nam=='bsr':
            roi['grd'][nam]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Disturbances\\VEG_BURN_SEVERITY_SP_2017.tif')
            roi['grd'][nam]['Data']=np.squeeze(roi['grd'][nam]['Data'])
            roi['grd'][nam]=gis.ClipToRaster(roi['grd'][nam],roi['grd'])
            roi['grd'][nam]['key']=gu.ReadExcel(meta['Paths']['bc1ha'] + '\\Disturbances\\VEG_BURN_SEVERITY_SP.xlsx')
        elif nam=='pfi_c':
            roi['grd'][nam]=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\ForestInventory\LiDAR\Boundary TSA\From GQ 20230302\BoundaryTSA_PFI_LnTotalWSV_bc1ha.tif')
            roi['grd'][nam]['Data']=np.squeeze(roi['grd'][nam]['Data'])
            roi['grd'][nam]['Data']=0.5*0.5*roi['grd'][nam]['Data']
            roi['grd'][nam]=gis.ClipToRaster(roi['grd'][nam],roi['grd'])
        elif nam=='plam':
            roi['grd'][nam]=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Results\Planting_FromRESULTS_MaskCount.tif')
        elif nam=='prcp_ann_n':
            roi['grd'][nam]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Climate\\Seasonal\\BC1ha_prcp_ann_norm_1971to2000_si_hist_v1.tif')
            roi['grd'][nam]['Data']=np.squeeze(roi['grd'][nam]['Data'])
            roi['grd'][nam]=gis.ClipToRaster(roi['grd'][nam],roi['grd'])
        elif nam=='protected':
            roi['grd'][nam]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\PROTECTED_LANDS_DESIGNATION.tif')
            roi['grd'][nam]['Data']=np.squeeze(roi['grd'][nam]['Data'])
            roi['grd'][nam]=gis.ClipToRaster(roi['grd'][nam],roi['grd'])
        elif nam=='soc':
            roi['grd'][nam]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Soil\\soc_tot_forest_Shawetal2018.tif')
            roi['grd'][nam]=gis.ClipToRaster(roi['grd'][nam],roi['grd'])
        elif nam=='spc1_ntems':
            roi['grd'][nam]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Species\\NTEMS_Spc1_Filled_RD_CAPITAL.tif')
        elif nam=='tmean_ann_n':
            roi['grd'][nam]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Climate\\Seasonal\\BC1ha_tmean_ann_norm_1971to2000_si_hist_v1.tif')
        elif nam=='regentype':
            roi['grd'][nam]=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Results\Harvest_Regen_Type.tif')
        elif nam=='rescon':
            roi['grd'][nam]=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Results\SILV_RESERVE_CODE_Consolidated.tif')
        elif nam=='rangecon':
            roi['grd'][nam]=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\LandCoverUse\rangeten_consol.tif')
        elif nam=='salvage':
            roi['grd'][nam]=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Disturbances\HarvestSalvageMask_FromCruise.tif')
        elif nam=='si_spl_ntems':
            roi['grd'][nam]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\SPL\\SI_From_NTEMS_Spc1_Filled_RD_CAPITAL_GF.tif')
        elif nam=='si_spl_fd':
            roi['grd'][nam]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\SPL\\Site_Prod_Fd.tif')
        elif nam=='si_spl_hw':
            roi['grd'][nam]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\SPL\\Site_Prod_Hw.tif')
        elif nam=='ws_gs_n':
            tmp=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
            z_tmp=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Climate\Seasonal\wbm_ws_gs_norm_1971to2000.tif')
            roi['grd'][nam]=tmp.copy()
            roi['grd'][nam]['Data']=z_tmp['Data']
            del z_tmp
            roi['grd'][nam]=gis.ClipToRaster(roi['grd'][nam],roi['grd'])
            roi['grd'][nam]['Data']=roi['grd'][nam]['Data'].astype('float')
            roi['grd'][nam]['cm']=gu.ReadExcel(r'C:\Users\rhember\Documents\Data\Colormaps\colormap_ws.xlsx')
        elif nam=='wf':
            roi['grd'][nam]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Disturbances\\PROT_HISTORICAL_FIRE_POLYS_SP_2017.tif')
            roi['grd'][nam]=gis.ClipToRaster(roi['grd'][nam],roi['grd'])

        # Clip
        roi['grd'][nam]=gis.ClipToRaster(roi['grd'][nam],roi['grd'])

    return roi

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

#%% NALCMS CEC 2020 - Compress categories to remove tropics

def NALCMS_Compress(lut_in,zRef,zCEC10,zCEC20):

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

#%% Clip rasters to standard grid
# *** This also compresses files that come out of Arc crazy big ***

def ClipToBC1ha(meta):

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

    # BGC zone
    fin=r'C:\Users\rhember\Documents\Data\BC1ha\VRI\becz.tif'
    fout=r'C:\Users\rhember\Documents\Data\BC1ha\VRI\becz.tif'
    gis.ClipToRaster_ByFile(fin,fout,meta['Paths']['bc1ha Ref Grid'])

    # BGC subzone
    fin=r'C:\Users\rhember\Documents\Data\BC1ha\VRI\becsz.tif'
    fout=r'C:\Users\rhember\Documents\Data\BC1ha\VRI\becsz.tif'
    gis.ClipToRaster_ByFile(fin,fout,meta['Paths']['bc1ha Ref Grid'])

    # VRI age
    fin=r'C:\Users\rhember\Documents\Data\BC1ha\VRI\proj_age_1.tif'
    fout=r'C:\Users\rhember\Documents\Data\BC1ha\VRI\proj_age_1.tif'
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
    tv=np.arange(1900,2025,1)
    N=np.zeros(tv.size)
    for i in range(nPack):
        z=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\' + lNam + '\\' + vNam + '_' + str(i+1) + '_Year.tif')['Data']
        u,c=np.unique(z,return_counts=True)
        for j in range(u.size):
            ind=np.where(tv==u[j])[0]
            N[ind]=N[ind]+c[j]
    return tv,N