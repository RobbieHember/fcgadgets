
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
import fcgadgets.macgyver.util_general as gu
import fcgadgets.macgyver.util_gis as gis
import fcgadgets.macgyver.util_query_gdb as qgdb
import fcexplore.psp.Processing.psp_util as ugp
import fcgadgets.cbrunner.cbrun_util as cbu

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

    # Tracking parameters
    meta['Graphics']['Fig Count']=1
    meta['Graphics']['Tab Count']=1

    # Initiate geospatial info
    meta['Geos']={}

    # Import variable info
    meta['Geos']['Variable Info']=gu.ReadExcel(meta['Paths']['Model']['Code'] + '\\Parameters\\Parameters_BC1haRasterVariableList.xlsx')

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

    d=gu.ReadExcel(meta['Paths']['bc1ha'] + '\\LUTs\\LUT_regentype_no.xlsx')
    vNam='RegenTypeNO'
    meta['LUT']['Raw'][vNam]=d
    meta['LUT'][lNam][vNam]={}
    for i in range(d['ID'].size):
        meta['LUT'][lNam][vNam][ d['Name'][i] ]=d['ID'][i]

    d=gu.ReadExcel(meta['Paths']['bc1ha'] + '\\LUTs\\LUT_EcozoneCanada.xlsx')
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

def BuildLUTsFromSourceDBs(meta):

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
