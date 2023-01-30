
#%% Import modules

import os
import numpy as np
import gc
from osgeo import gdal
import matplotlib.pyplot as plt
import matplotlib.colors
import geopandas as gpd
import pandas as pd
import fiona
from shapely.geometry import Polygon,Point,box,shape
import fcgadgets.macgyver.utilities_general as gu
import fcgadgets.macgyver.utilities_gis as gis
import fcgadgets.macgyver.utilities_query_gdb as qgdb
import fcgadgets.cbrunner.cbrun_utilities as cbu

#%% Region of interest

def DefineROI(roi,gdf):

    # Import TSA grid
    tsa=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Admin\tsa.tif')
    tsa['Data']=np.squeeze(tsa['Data'])

    roi['crs']=gdf['bc_land']['gdf'].crs

    # Vector:
    roi['gdf']={}

    if roi['Type']=='ByTSA':

        # Index to TSAs in query
        iROI=gdf['tsa']['key'].VALUE[np.isin(gdf['tsa']['key'].Name,roi['TSA List'])].values

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

    #roi['gdf']['tsa within']=gpd.overlay(gdf['tsa']['gdf'],roi['gdf']['bound'],how='intersection')
    roi['gdf']['tsa']=gdf['tsa']['gdf'].cx[roi['grd']['xmin']:roi['grd']['xmax'],roi['grd']['ymin']:roi['grd']['ymax']]
    roi['gdf']['tsa']=roi['gdf']['tsa'].reset_index(drop=True)
    roi['gdf']['tsa']=gpd.sjoin(roi['gdf']['tsa'],roi['gdf']['bound'],how='left')
    roi['gdf']['tsa']=gpd.overlay(roi['gdf']['tsa'],roi['gdf']['bound'],how='intersection')

    if roi['Type']=='ByTSA':
        roi['gdf']['tsa within']=gdf['tsa']['gdf'].iloc[np.isin(gdf['tsa']['gdf'].Name,roi['TSA List'])]

    roi['gdf']['lakes']=gpd.overlay(gdf['bc_land']['gdf'][(gdf['bc_land']['gdf']['TAG']=='lake')],roi['gdf']['bound'],how='intersection')

    roi['gdf']['rivers']=gpd.overlay(gdf['bc_land']['gdf'][(gdf['bc_land']['gdf']['TAG']=='river')],roi['gdf']['bound'],how='intersection')

    roi['gdf']['tpf']=gdf['tpf']['gdf'].cx[roi['grd']['xmin']:roi['grd']['xmax'],roi['grd']['ymin']:roi['grd']['ymax']]
    roi['gdf']['tpf']=roi['gdf']['tpf'].reset_index(drop=True)

    roi['gdf']['fnc']=gdf['fnc']['gdf'].cx[roi['grd']['xmin']:roi['grd']['xmax'],roi['grd']['ymin']:roi['grd']['ymax']]
    roi['gdf']['fnc']=roi['gdf']['fnc'].reset_index(drop=True)

    try:
        roi['gdf']['road']=gdf['road']['gdf'].cx[roi['grd']['xmin']:roi['grd']['xmax'],roi['grd']['ymin']:roi['grd']['ymax']]
        roi['gdf']['road']=roi['gdf']['road'].reset_index(drop=True)
    except:
        pass

    return roi

#%% Import basemaps

def Import_GDBs_ProvinceWide():

    path_inf=r'C:\Users\rhember\Documents\Data\ForestInventory\Infrastructure.gdb'
    path_basemaps=r'C:\Users\rhember\Documents\Data\Basemaps\Basemaps.gdb'

    flg=0
    if flg==1:
        fiona.listlayers(path_inf)
        fiona.listlayers(path_basemaps)

    gdf={}

    # Politial boundary
    gdf['bc_bound']={}
    gdf['bc_bound']['gdf']=gpd.read_file(path_basemaps,layer='NRC_POLITICAL_BOUNDARIES_1M_SP')

    # Not using
    gdf['bc_land']={}
    gdf['bc_land']['gdf']=gpd.read_file(r'C:\Users\rhember\Documents\Data\Basemaps\bc_land.shp')

    # Import TSA info
    gdf['tsa']={}
    gdf['tsa']['key']=pd.read_excel(r'C:\Users\rhember\Documents\Data\BC1ha\Admin\lut_tsa.xlsx')
    gdf['tsa']['gdf']=gpd.read_file(r'C:\Users\rhember\Documents\Data\TSA\tsa_boundaries.shp')

    # Roads
    gdf['road']={}
    gdf['road']['gdf']=gpd.read_file(path_inf,layer='MOT_ROAD_FEATURES_INVNTRY_SP')
    #gpd.read_file(r'C:\Users\rhember\Documents\Data\Basemaps\roads.shp')

    # Timber processing facilities
    gdf['tpf']={}
    gdf['tpf']['gdf']=gpd.read_file(path_inf,layer='GSR_TMBR_PRCSSING_FAC_SV')

    # First nations communities
    gdf['fnc']={}
    gdf['fnc']['gdf']=gpd.read_file(path_inf,layer='FN_COMMUNITY_LOCATIONS_SP')

    # Oil and gas facilities
    gdf['ogf']={}
    gdf['ogf']['gdf']=gpd.read_file(path_inf,layer='DRP_OIL_GAS_FACILITIES_BC_SP')

    # Oil and gas pipeline
    gdf['ogp']={}
    gdf['ogp']['gdf']=gpd.read_file(path_inf,layer='DRP_OIL_GAS_PIPELINES_BC_SP')

    # Cities
    gdf['cities']=gis.ImportCities(r'C:\Users\rhember\Documents\Data\Cities\Cities.xlsx','GDF')
    gdf['cities'].crs=gdf['bc_bound']['gdf'].crs

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

    return roi

#%% Import variables for ROI

def Import_Raster_Over_ROI(meta_bc1ha,roi,vList):

    xlim=(roi['grd']['xmin'],roi['grd']['xmax']+roi['grd']['Cellsize'])
    ylim=(roi['grd']['ymin'],roi['grd']['ymax']+roi['grd']['Cellsize'])

    for nam in vList:

        if nam in roi['grd'].keys():
            continue

        if nam=='lc2':
            #land cover scheme level 2 (Treed=4)
            roi['grd'][nam]=gis.OpenGeoTiff(meta_bc1ha['Paths']['BC1ha'] + '\\VRI\\lc2.tif')
            roi['grd'][nam]['Data']=np.squeeze(roi['grd'][nam]['Data'])
            roi['grd'][nam]=gis.ClipToRaster(roi['grd'][nam],roi['grd'])
            gc.collect()

        elif nam=='btm':
            # Import Base Thematic Map
            roi['grd'][nam]=gis.OpenGeoTiff(meta_bc1ha['Paths']['BC1ha'] + '\\LandUseLandCover\\landuse.btm.tif')
            roi['grd'][nam]['Data']=np.squeeze(roi['grd'][nam]['Data'])
            roi['grd'][nam]=gis.ClipToRaster(roi['grd'][nam],roi['grd'])
            lut=pd.read_csv(r'C:\Users\rhember\Documents\Data\BC1ha\LandUseLandCover\landuse_btm_category_metadata.csv')
            cl=np.column_stack( (lut['C1'].values,lut['C2'].values,lut['C3'].values) )
            roi['grd'][nam]['Compressed']={}
            roi['grd'][nam]['Compressed']['Data'],roi['grd'][nam]['Compressed']['lab'],roi['grd'][nam]['Compressed']['cl1']=gis.CompressCats(roi['grd'][nam]['Data'],lut['Raster Value'].values,lut['PLU Label'].values,cl)

        elif nam=='becz':
            # BGC classification
            roi['grd'][nam]=gis.OpenGeoTiff(meta_bc1ha['Paths']['BC1ha'] + '\\VRI\\becz.tif')
            roi['grd'][nam]['Data']=np.squeeze(roi['grd'][nam]['Data'])
            roi['grd'][nam]=gis.ClipToRaster(roi['grd'][nam],roi['grd'])
            roi['grd'][nam]['key']=pd.read_excel(meta_bc1ha['Paths']['BC1ha'] + '\\VRI\\becz_lut.xlsx')

        elif nam=='age1':
            # Age 1
            roi['grd'][nam]=gis.OpenGeoTiff(meta_bc1ha['Paths']['BC1ha'] + '\\VRI\\proj_age_1.tif')
            roi['grd'][nam]['Data']=np.squeeze(roi['grd'][nam]['Data'])
            roi['grd'][nam]=gis.ClipToRaster(roi['grd'][nam],roi['grd'])

        elif nam=='si':
            # Site index
            roi['grd'][nam]=gis.OpenGeoTiff(meta_bc1ha['Paths']['BC1ha'] + '\\VRI\\si.tif')
            roi['grd'][nam]['Data']=np.squeeze(roi['grd'][nam]['Data'])
            roi['grd'][nam]=gis.ClipToRaster(roi['grd'][nam],roi['grd'])

        elif nam=='sphlive':
            # SPH
            roi['grd'][nam]=gis.OpenGeoTiff(meta_bc1ha['Paths']['BC1ha'] + '\\VRI\\sphlive.tif')
            roi['grd'][nam]['Data']=np.squeeze(roi['grd'][nam]['Data'])
            roi['grd'][nam]=gis.ClipToRaster(roi['grd'][nam],roi['grd'])

        elif nam=='sphdead':
            # SPH
            roi['grd'][nam]=gis.OpenGeoTiff(meta_bc1ha['Paths']['BC1ha'] + '\\VRI\\sphdead.tif')
            roi['grd'][nam]['Data']=np.squeeze(roi['grd'][nam]['Data'])
            roi['grd'][nam]=gis.ClipToRaster(roi['grd'][nam],roi['grd'])

        elif nam=='cut_yr':

            roi['grd'][nam]=roi['grd'].copy()
            roi['grd'][nam]['Data']=0*roi['grd'][nam]['Data']
            for i in range(1,5):
                tmp=gis.OpenGeoTiff(meta_bc1ha['Paths']['BC1ha'] + '\\Disturbances\\VEG_CONSOLIDATED_CUT_BLOCKS_SP_year' + str(i) + 'a.tif')
                tmp['Data']=np.squeeze(tmp['Data'])
                tmp=gis.ClipToRaster(roi['grd'][nam],roi['grd'])
                roi['grd'][nam]['Data']=np.maximum(roi['grd'][nam]['Data'],tmp['Data'])
                del tmp
                gc.collect()

        elif nam=='bsr':

            roi['grd'][nam]=gis.OpenGeoTiff(meta_bc1ha['Paths']['BC1ha'] + '\\Disturbances\\VEG_BURN_SEVERITY_SP_2017.tif')
            roi['grd'][nam]['Data']=np.squeeze(roi['grd'][nam]['Data'])
            roi['grd'][nam]=gis.ClipToRaster(roi['grd'][nam],roi['grd'])
            roi['grd'][nam]['key']=gu.ReadExcel(meta_bc1ha['Paths']['BC1ha'] + '\\Disturbances\\VEG_BURN_SEVERITY_SP.xlsx')

        elif nam=='wf':

            roi['grd'][nam]=gis.OpenGeoTiff(meta_bc1ha['Paths']['BC1ha'] + '\\Disturbances\\PROT_HISTORICAL_FIRE_POLYS_SP_2017.tif')
            roi['grd'][nam]['Data']=np.squeeze(roi['grd'][nam]['Data'])
            roi['grd'][nam]=gis.ClipToRaster(roi['grd'][nam],roi['grd'])

        elif nam=='elev':
            roi['grd'][nam]=gis.OpenGeoTiff(meta_bc1ha['Paths']['BC1ha'] + '\\Terrain\\elevation.tif')
            roi['grd'][nam]['Data']=np.squeeze(roi['grd'][nam]['Data'])
            roi['grd'][nam]=gis.ClipToRaster(roi['grd'][nam],roi['grd'])

        elif nam=='soc':
            roi['grd'][nam]=gis.OpenGeoTiff(meta_bc1ha['Paths']['BC1ha'] + '\\Soil\\soc_tot_forest_Shawetal2018.tif')
            roi['grd'][nam]['Data']=np.squeeze(roi['grd'][nam]['Data'])
            roi['grd'][nam]=gis.ClipToRaster(roi['grd'][nam],roi['grd'])

        elif nam=='temp_norm':
            tsa=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Admin\tsa.tif')
            tsa['Data']=np.squeeze(tsa['Data'])
            z_tmp=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Climate\BC1ha_tmin_ann_norm_1971to2000_si_hist_v1.tif')
            roi['grd'][nam]=tsa.copy()
            roi['grd'][nam]['Data']=np.squeeze(z_tmp['Data'])
            roi['grd'][nam]=gis.ClipToRaster(roi['grd'][nam],roi['grd'])
            roi['grd'][nam]['Data']=roi['grd'][nam]['Data'].astype('float')/10
            gc.collect()

        elif nam=='ws_norm':
            tsa=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Admin\tsa.tif')
            tsa['Data']=np.squeeze(tsa['Data'])
            z_tmp=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Climate\BC1ha_ws_gs_norm_1971to2000_comp_hist_v1.tif')
            roi['grd'][nam]=tsa.copy()
            roi['grd'][nam]['Data']=z_tmp['Data']
            del z_tmp
            roi['grd'][nam]=gis.ClipToRaster(roi['grd'][nam],roi['grd'])
            roi['grd'][nam]['Data']=roi['grd'][nam]['Data'].astype('float')
            roi['grd'][nam]['cm']=gu.ReadExcel(r'C:\Users\rhember\Documents\Data\Colormaps\colormap_ws.xlsx')
            gc.collect()
            #plt.matshow(zW['Data']);plt.colorbar()
        elif nam=='d2road':
            roi['grd'][nam]=gis.OpenGeoTiff(meta_bc1ha['Paths']['BC1ha'] + '\\Terrain\\DistanceFromRoads.tif')
            roi['grd'][nam]['Data']=np.squeeze(roi['grd'][nam]['Data'])
            roi['grd'][nam]=gis.ClipToRaster(roi['grd'][nam],roi['grd'])
            gc.collect()
        elif nam=='d2fac':
            roi['grd'][nam]=gis.OpenGeoTiff(meta_bc1ha['Paths']['BC1ha'] + '\\Management\\DistanceFromForestryFacility.tif')
            roi['grd'][nam]['Data']=np.squeeze(roi['grd'][nam]['Data'])
            roi['grd'][nam]=gis.ClipToRaster(roi['grd'][nam],roi['grd'])
            gc.collect()
        elif nam=='idw_mask':
            roi['grd'][nam]=gis.OpenGeoTiff(meta_bc1ha['Paths']['BC1ha'] + '\\Disturbances\\IDW_Mask.tif')
            roi['grd'][nam]['Data']=np.squeeze(roi['grd'][nam]['Data'])
            roi['grd'][nam]=gis.ClipToRaster(roi['grd'][nam],roi['grd'])
        elif nam=='protected':
            roi['grd'][nam]=gis.OpenGeoTiff(meta_bc1ha['Paths']['BC1ha'] + '\\LandUseLandCover\\PROTECTED_LANDS_DESIGNATION.tif')
            roi['grd'][nam]['Data']=np.squeeze(roi['grd'][nam]['Data'])
            roi['grd'][nam]=gis.ClipToRaster(roi['grd'][nam],roi['grd'])

    return roi

#%% PLOT ROI mask

def Plot_ROI_Mask(meta_bc1ha,roi,gdf):

    # Grid and labels
    lab=[]
    z1=np.ones(roi['grd']['lc2']['Data'].shape)
    z1[(roi['grd']['Data']==1) & (roi['grd']['lc2']['Data']==4)]=0; lab.append('Treed')
    z1[(roi['grd']['Data']==1) & (roi['grd']['lc2']['Data']!=4)]=1; lab.append('Non-treed')
    z1[(roi['grd']['Data']!=1) & (roi['grd']['lc2']['Data']!=1)]=2; lab.append('hidden')
    z1[(roi['grd']['Data']!=1) & (roi['grd']['lc2']['Data']==1)]=3; lab.append('hidden')

    # Number of colours and number of colours excluded from colorbar
    N_color=4
    N_hidden=2

    # Colormap
    cm=np.vstack( ((0.7,0.7,0.7,1),(0.8,0.8,0.8,1),(0.93,0.93,0.93,1),(1,1,1,1)) )
    cm=matplotlib.colors.ListedColormap(cm)

    fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta_bc1ha['Graphics']['figwidth'],(1-meta_bc1ha['Graphics']['sidespace'])*meta_bc1ha['Graphics']['figwidth']*roi['grd']['yxrat']))
    im=ax[0].matshow(z1[0::1,0::1],clim=(0,N_color),extent=roi['grd']['lc2']['Extent'],cmap=cm)
    gdf['bc_bound']['gdf'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
    roi['gdf']['bound'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
    roi['gdf']['tsa'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)

    #roi['gdf_lakes'].plot(ax=ax[0],facecolor=[0.82,0.88,1],edgecolor=[0.7*0.82,0.7*0.88,0.7*1],linewidth=0.25,label='Water')
#    try:
#        roi['gdf_rivers'].plot(ax=ax[0],linecolor=[0,0,0.7],label='Water',linewidth=3.25)
#    except:
#        pass
#    roi['gdf']['bound'].plot(ax=ax[0],color=None,edgecolor=[0,0,0],facecolor='none')
#    try:
#        roi['gdf_roads'].plot(ax=ax[0],facecolor='none',edgecolor=[0,0,0],label='Roads',linewidth=0.75,alpha=1,zorder=1)
#    except:
#        # No roads exist in areas
#        pass
    ax[0].set(position=meta_bc1ha['Graphics']['ax1 pos'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'],aspect='auto')
    ax[0].yaxis.set_ticks_position('both');
    ax[0].xaxis.set_ticks_position('both')
    ax[0].grid(meta_bc1ha['Graphics']['ax1 gridvis'])
    ax[0].axis(meta_bc1ha['Graphics']['ax1 vis'])

    cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_color-1,1),ticks=np.arange(0.5,N_color+1.5,1))
    #cb.ax.set(yticklabels=lab)
    cb.ax.tick_params(labelsize=6,length=0)
    for i in range(0,N_color):
        ax[1].plot([0,100],[i,i],'k-',linewidth=0.5)
    pos2=meta_bc1ha['Graphics']['ax2 pos']
    pos2[1]=0.8
    pos2[3]=0.1
    ax[1].set(position=pos2)

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
            gdf.iloc[ind].plot(ax=ax[0],marker='o',edgecolor='g',facecolor='g',lw=lw,markersize=ms,alpha=0.35,zorder=2)
        elif (mtype=='PLY') | (mtype=='VNR') | (mtype=='OSB') | (mtype=='PNL'):
            # One PNL (panel) mill WestPine MDF in Quesnell - it is MDF
            y=gdf.iloc[ind]['EST_AN_CAP_MLN_SQ_FT']/885
            ms=300*y
            gdf.iloc[ind].plot(ax=ax[0],marker='o',edgecolor='c',facecolor='c',lw=lw,markersize=ms,alpha=0.35,zorder=2)
        elif mtype=='PLT':
            y=gdf.iloc[ind]['EST_AN_CAP_000_TONNES']*2
            ms=1.0*y
            gdf.iloc[ind].plot(ax=ax[0],marker='o',edgecolor='r',facecolor='r',lw=lw,markersize=ms,alpha=0.35,zorder=2)
        elif mtype=='CHP':
            y=gdf.iloc[ind]['EST_AN_CAP_000_BDUS']*2
            ms=0.75*y
            gdf.iloc[ind].plot(ax=ax[0],marker='o',edgecolor='y',facecolor='y',lw=lw,markersize=ms,alpha=0.35,zorder=2)
        elif mtype=='PLP':
            y=gdf.iloc[ind]['EST_AN_CAP_000_TONNES']*2
            ms=0.75*y
            gdf.iloc[ind].plot(ax=ax[0],marker='o',edgecolor='k',facecolor='k',lw=lw,markersize=ms,alpha=0.35,zorder=2)
        elif mtype=='LVL':
            # Laminated veneer lumber
            y=gdf.iloc[ind]['EST_AN_CAP_MLN_CUBIC_FT']#/885
            ms=200*y
            gdf.iloc[ind].plot(ax=ax[0],marker='o',edgecolor='k',facecolor='g',lw=lw,markersize=ms,alpha=0.35,zorder=2)

        if labels=='On':
            for x,y,label in zip(gdf.iloc[ind].geometry.x,gdf.iloc[ind].geometry.y,gdf.iloc[ind].COMPANY_NAME):
                ax[0].annotate(label,xy=(x,y),xytext=(5,4),textcoords="offset points")
    return ax

