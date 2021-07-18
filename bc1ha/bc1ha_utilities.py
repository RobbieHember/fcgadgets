
#%% Import modules

import os
import numpy as np
import gc
import gdal
import matplotlib.pyplot as plt
import matplotlib.colors
import geopandas as gpd
import pandas as pd
import fiona
from shapely.geometry import Polygon,Point,box,shape
import fcgadgets.macgyver.utilities_general as gu
import fcgadgets.macgyver.utilities_gis as gis
from fcgadgets.cbrunner import cbrun_utilities as cbu

#%% Import basemaps

def Import_BaseMaps():

    # BC basemap info
    bm={}
    bm['gdf_bc_bound']=gpd.read_file(r'C:\Users\rhember\Documents\Data\Basemaps\Basemaps.gdb',layer='NRC_POLITICAL_BOUNDARIES_1M_SP')
    bm['gdf_bm']=gpd.read_file(r'C:\Users\rhember\Documents\Data\Basemaps\bc_land.shp')

    # Import TSA info
    tsa={}
    tsa['grd']=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Admin\tsa.tif')
    tsa['key']=pd.read_excel(r'C:\Users\rhember\Documents\Data\BC1ha\Admin\lut_tsa.xlsx')
    tsa['gdf_bound']=gpd.read_file(r'C:\Users\rhember\Documents\Data\TSA\tsa_boundaries.shp')

    # Roads
    road={}
    road['gdf']=gpd.read_file(r'C:\Users\rhember\Documents\Data\Basemaps\roads.shp')

    # Districts
    district={}
    try:
        district['gdf']=gpd.read_file(r'Z:\!Workgrp\Forest Carbon\Data\Districts\district.shp')
    except:
        pass

    return bm,tsa,road,district

#%% Region of interest

def DefineROI(roi,tsa,bm,road):

    if roi['Type']=='ByTSA':
        
        # Index to TSAs in query
        iROI=tsa['key'].VALUE[np.isin(tsa['key'].Name,roi['TSA List'])].values

        roi['Mask']=tsa['grd'].copy()
        roi['Mask']['Data']=np.zeros(tsa['grd']['Data'].shape)
        ind=(np.isin(tsa['grd']['Data'],iROI))
        roi['Mask']['Data'][ind]=1

        # Define extent based on mask
        roi['xlim']=[np.min(tsa['grd']['X'][ind])-5000,np.max(tsa['grd']['X'][ind])+5000]
        roi['ylim']=[np.min(tsa['grd']['Y'][ind])-5000,np.max(tsa['grd']['Y'][ind])+5000]

        # Clip mask
        roi['Mask']=gis.ClipRaster(roi['Mask'],roi['xlim'],roi['ylim'])
        gc.collect()

        roi['gdf_bound']=tsa['gdf_bound'][np.isin(tsa['gdf_bound'].Name,roi['TSA List'])]
        roi['gdf_bound']=roi['gdf_bound'].reset_index(drop=True)
        
        roi['gdf_lakes']=gpd.overlay(bm['gdf_bm'][(bm['gdf_bm']['TAG']=='lake')],tsa['gdf_bound'][np.isin(tsa['gdf_bound'].Name,roi['TSA List'])],how='intersection')
        
        roi['gdf_rivers']=gpd.overlay(bm['gdf_bm'][(bm['gdf_bm']['TAG']=='river')],tsa['gdf_bound'][np.isin(tsa['gdf_bound'].Name,roi['TSA List'])],how='intersection')
    
        try:
            roi['gdf_roads']=road['gdf'].cx[roi['Mask']['xmin']:roi['Mask']['xmax'],roi['Mask']['ymin']:roi['Mask']['ymax']]
            roi['gdf_roads']=roi['gdf_roads'].reset_index(drop=True)
            roi['gdf_roads']=gpd.sjoin(roi['gdf_roads'],roi['gdf_bound'],how='left')
            roi['gdf_roads']=roi['gdf_roads'].groupby('index_right')
        except:
            roi['gdf_roads']=[]
    
    elif roi['Type']=='ByLatLon':
        
        srs=gis.ImportSRSs()
        xc,yc=srs['Proj']['BC1ha'](roi['Centre'][0],roi['Centre'][1])
        xc=np.round(xc/100)*100
        yc=np.round(yc/100)*100
        
        # Define extent based on mask
        roi['xlim']=[xc-roi['Radius'],xc+roi['Radius']]
        roi['ylim']=[yc-roi['Radius'],yc+roi['Radius']]
    
        roi['Mask']=gis.ClipRaster(tsa['grd'],roi['xlim'],roi['ylim'])
        roi['Mask']['Data']=0*roi['Mask']['Data']+1
        
        #box(W, S, E, N)
        geom=box(roi['Mask']['Extent'][0],roi['Mask']['Extent'][2],roi['Mask']['Extent'][1],roi['Mask']['Extent'][3])
        roi['gdf_bound']=gpd.GeoDataFrame({"id":1,"geometry":[geom]})
        
        roi['gdf_lakes']=gpd.overlay(bm['gdf_bm'][(bm['gdf_bm']['TAG']=='lake')],roi['gdf_bound'],how='intersection')
        
        roi['gdf_rivers']=gpd.overlay(bm['gdf_bm'][(bm['gdf_bm']['TAG']=='river')],roi['gdf_bound'],how='intersection')
    
        try:
            roi['gdf_roads']=road['gdf'].cx[roi['Mask']['xmin']:roi['Mask']['xmax'],roi['Mask']['ymin']:roi['Mask']['ymax']]
            roi['gdf_roads']=roi['gdf_roads'].reset_index(drop=True)
            roi['gdf_roads']=gpd.sjoin(roi['gdf_roads'],roi['gdf_bound'],how='left')
            roi['gdf_roads']=roi['gdf_roads'].groupby('index_right')
        except:
            # If there are no roads in the ROI, the sjoin fails
            roi['gdf_roads']=[]
        
    return roi

#%% Clip geodataframe to ROI

def ClipGDF_ByROI(gdf_in,roi):
    
    gdf_out=gdf_in.cx[roi['Mask']['xmin']:roi['Mask']['xmax'],roi['Mask']['ymin']:roi['Mask']['ymax']]
    gdf_out=gdf_out.reset_index(drop=True)
    gdf_out=gpd.sjoin(gdf_out,roi['gdf_bound'],how='left')
    gdf_out=gdf_out.groupby('index_right')
    
    return gdf_out

#%% Import variables for ROI

def Import_Raster_Over_ROI(v_cd,roi):

    d_out={}
    
    if v_cd=='lc2':
            
        #land cover scheme level 2 (Treed=4)        
        d_out['grd']=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\VRI\lc2.tif')
        d_out['grd']=gis.ClipRaster(d_out['grd'],roi['xlim'],roi['ylim'])        
        #d_out['grd']=tsa['grd'].copy()
        #tmp=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\VRI\lc2.tif')
        #d_out['grd']['Data']=tmp['Data']
        #del tmp
        #d_out['grd']=gis.ClipRaster(d_out['grd'],roi['xlim'],roi['ylim'])
        gc.collect()
            
    elif v_cd=='btm':

        # Import Base Thematic Map
        d_out['grd']=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\LandUseLandCover\landuse.btm.tif')
        d_out['grd']=gis.ClipRaster(d_out['grd'],roi['xlim'],roi['ylim'])        
        lut=pd.read_csv(r'C:\Users\rhember\Documents\Data\BC1ha\LandUseLandCover\landuse_btm_category_metadata.csv')
        cl=np.column_stack( (lut['C1'].values,lut['C2'].values,lut['C3'].values) )
        d_out['Data1'],d_out['lab1'],d_out['cl1']=gis.CompressCats(d_out['grd']['Data'],lut['Raster Value'].values,
             lut['PLU Label'].values,cl)
    
    elif v_cd=='bgcz':
        
        # BGC classification
        d_out['grd']=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\VRI\becz.tif')
        d_out['grd']=gis.ClipRaster(d_out['grd'],roi['xlim'],roi['ylim'])  
        d_out['key']=pd.read_excel(r'C:\Users\rhember\Documents\Data\BC1ha\VRI\becz.tif.vat.xlsx')
        #zBGC=tsa['grd'].copy()
        #zBGC['Data']=zBGC_tmp['Data']
        #del zBGC_tmp        
        #zBGC=gis.ClipRaster(zBGC,xlim,ylim)
        #gc.collect()
        
    elif v_cd=='cut_yr':
        
        d_out['grd']=roi['Mask'].copy()
        d_out['grd']['Data']=0*d_out['grd']['Data']
        for i in range(1,5):
            tmp=gis.OpenGeoTiff(r'Z:\!Workgrp\Forest Carbon\Data\BC1ha\Disturbances\VEG_CONSOLIDATED_CUT_BLOCKS_SP_year' + str(i) + '.tif')
            tmp=gis.ClipRaster(tmp,roi['xlim'],roi['ylim'])    
            d_out['grd']['Data']=np.maximum(d_out['grd']['Data'],tmp['Data'])
            del tmp
            gc.collect()        
    
    elif v_cd=='bsr':
    
        d_out['grd']=gis.OpenGeoTiff(r'Z:\!Workgrp\Forest Carbon\Data\BC1ha\Disturbances\VEG_BURN_SEVERITY_SP_2017.tif')    
        d_out['grd']=gis.ClipRaster(d_out['grd'],roi['xlim'],roi['ylim'])
        d_out['key']=gu.ReadExcel(r'Z:\!Workgrp\Forest Carbon\Data\BC1ha\Disturbances\VEG_BURN_SEVERITY_SP.xlsx')  
    
    return d_out

#%% Query openings
 
def Get_Vectors_For_ROI(roi,var_nam,year_start,year_end):
    
    # Load dataset with CRS
    gdf_bm=gpd.read_file(r'C:\Users\rhember\Documents\Data\Basemaps\bc_land.shp')
    crs=gdf_bm.crs
    del gdf_bm
    
    # Convert to shapely object
    roi_poly=roi['gdf_bound'].iloc[0]['geometry']
    
    if var_nam=='Openings':
    
        fin=r'C:\Users\rhember\Documents\Data\ForestInventory\Results\20210401\Results.gdb'    
        List=[None]*int(1e5)
        cnt=0
        with fiona.open(fin,layer='RSLT_OPENING_SVW') as source:
            for feat in source:
                if feat['geometry']==None:
                    continue                        
                s=shape(feat['geometry'])            
                if (s.overlaps(roi_poly)==False) & (s.within(roi_poly)==False):
                    continue            
                List[cnt]=feat
                cnt=cnt+1
        List=List[0:cnt-1]
    
    elif var_nam=='Planting':
        
        fin=r'C:\Users\rhember\Documents\Data\ForestInventory\Results\20210401\Results.gdb'
        List=[None]*int(1e5)
        cnt=0
        with fiona.open(fin,layer='RSLT_ACTIVITY_TREATMENT_SVW') as source:
            for feat in source:
                if feat['geometry']==None:
                    continue
                if feat['properties']['SILV_BASE_CODE']!='PL':
                    continue
                if feat['properties']['ATU_COMPLETION_DATE']==None:
                    continue
                Year=int(feat['properties']['ATU_COMPLETION_DATE'][0:4])
                if (Year<year_start) | (Year>year_end):
                    continue
                s=shape(feat['geometry'])
                if (s.overlaps(roi_poly)==False) & (s.within(roi_poly)==False):
                    continue            
                List[cnt]=feat
                cnt=cnt+1   
        List=List[0:cnt-1]
    
    elif var_nam=='vri':
        
        fin=r'C:\Users\rhember\Documents\Data\ForestInventory\VRI\20210401\VRI.gdb'
        List=[None]*int(1e5)
        cnt=0
        with fiona.open(fin,layer='VEG_COMP_LYR_R1_POLY') as source:
            for feat in source:
                if feat['geometry']==None:
                    continue
                
                shp=shape(feat['geometry'])
                bnd=shp.bounds # (W, S, E, N)
                
                # Check for overlap
                if (bnd[2]<roi['Mask']['Extent'][0]) | (bnd[0]>roi['Mask']['Extent'][1]) | (bnd[1]>roi['Mask']['Extent'][3]) | (bnd[3]<roi['Mask']['Extent'][2]):
                    continue
                
                #if (shp.overlaps(roi_poly)==False) & (shp.within(roi_poly)==False):
                #    continue 
                
                List[cnt]=feat
                cnt=cnt+1
                print('working')
        List=List[0:cnt-1]
        
    # Add to geodataframe
    gdf=gpd.GeoDataFrame.from_features(List,crs=crs)
    
    return gdf

#%% Get Survesy within ROI

def GetPlantingWithinROI(Year1,Year2,roi):
    
    fin=r'C:\Users\rhember\Documents\Data\ForestInventory\Results\20210401\Results.gdb'
    
    gdf_bm=gpd.read_file(r'C:\Users\rhember\Documents\Data\Basemaps\bc_land.shp')
    crs=gdf_bm.crs
    del gdf_bm
    
    # Convert to shapely object
    roi_poly=roi['gdf_bound'].iloc[0]['geometry']
    
    List=[None]*int(1e5)
    cnt=0
    with fiona.open(fin,layer='RSLT_ACTIVITY_TREATMENT_SVW') as source:
        for feat in source:
            if feat['geometry']==None:
                continue
            if feat['properties']['SILV_BASE_CODE']!='PL':
                continue
            if feat['properties']['ATU_COMPLETION_DATE']==None:
                continue
            Year=int(feat['properties']['ATU_COMPLETION_DATE'][0:4])
            if (Year<Year1) | (Year>Year2):
                continue
            s=shape(feat['geometry'])
            if (s.overlaps(roi_poly)==False) & (s.within(roi_poly)==False):
                continue
            
            List[cnt]=feat
            cnt=cnt+1   
    List=List[0:cnt-1]
    
    gdf=gpd.GeoDataFrame.from_features(List,crs=crs)
    
    #pls={}
    #pls['gdf']=gdf
    #gdf.plot()
    #gdf_fcinv2=gdf_fcinv.to_crs({'init':'epsg:4326'})
    # Save
    #gdf_fcinv.to_file(filename=meta['Paths']['Project'] + '\\Geospatial\\fcinv.geojson',driver='GeoJSON')
    return gdf

#%% Get Survesy within ROI

def GetSurveyWithinROI(Year1,Year2,roi):
    
    fin=r'C:\Users\rhember\Documents\Data\ForestInventory\Results\20210401\Results.gdb'
    
    gdf_bm=gpd.read_file(r'C:\Users\rhember\Documents\Data\Basemaps\bc_land.shp')
    crs=gdf_bm.crs
    del gdf_bm
    
    # Convert to shapely object
    roi_poly=roi['gdf_bound'].iloc[0]['geometry']
    
    List=[None]*int(1e5)
    cnt=0
    with fiona.open(fin,layer='RSLT_ACTIVITY_TREATMENT_SVW') as source:
        for feat in source:
            if feat['geometry']==None:
                continue
            if feat['properties']['SILV_BASE_CODE']!='SU':
                continue
            if feat['properties']['ATU_COMPLETION_DATE']==None:
                continue
            Year=int(feat['properties']['ATU_COMPLETION_DATE'][0:4])
            if (Year<Year1) | (Year>Year2):
                continue
            if feat['properties']['GEOMETRY_Area']/10000>2000:
                continue
            s=shape(feat['geometry'])
            if (s.overlaps(roi_poly)==False) & (s.within(roi_poly)==False):
                continue
            
            List[cnt]=feat
            cnt=cnt+1   
    List=List[0:cnt-1]
    
    gdf=gpd.GeoDataFrame.from_features(List,crs=crs)
    
    #su={}
    #su['gdf']=gdf
    #gdf.plot()
    #gdf_fcinv2=gdf_fcinv.to_crs({'init':'epsg:4326'})
    # Save
    #gdf_fcinv.to_file(filename=meta['Paths']['Project'] + '\\Geospatial\\fcinv.geojson',driver='GeoJSON')
    return gdf

#%% Get normal climate data

def GetClimateNormalFromGeoTiff(x,y,pthin,varl):

    d={}
    d['x']=x
    d['y']=y
    
    for iv in range(len(varl)):
    
        # Specify file name for each variable code and scale factor
        if varl[iv]=='tmean_gs':
            fin=pthin + '\\BC1ha_tmean_gs_norm_1971to2000_si_hist_v1.tif'
            sf=0.1
        elif varl[iv]=='tmin_ann':
            fin=pthin + '\\BC1ha_tmin_ann_norm_1971to2000_si_hist_v1.tif'
            sf=0.1
        elif varl[iv]=='etp_gs':
            fin=pthin + '\\BC1ha_etp_gs_norm_1971to2000_comp_hist_v1.tif'
            sf=0.1
        elif varl[iv]=='eta_gs':
            fin=pthin + '\\BC1ha_eta_gs_norm_1971to2000_comp_hist_v1.tif'
            sf=0.1            
        elif varl[iv]=='ws_gs':
            fin=pthin + '\\BC1ha_ws_gs_norm_1971to2000_comp_hist_v1.tif'
            sf=1
        
        # Load geotiff
        z=gis.OpenGeoTiff(fin)
    
        if iv==0:
            # Extract vector x and y
            xv=z['X'][0,:]
            yv=z['Y'][:,0]
        
            # Generate indices to x and y
            ix=np.zeros(x.size,dtype='int32')
            iy=np.zeros(x.size,dtype='int32')
            for i in range(x.size):
                dx=np.abs(xv-x[i])
                dy=np.abs(yv-y[i])
                ix[i]=np.where(dx==np.min(dx))[0][0]
                iy[i]=np.where(dy==np.min(dy))[0][0]
        
        # Extract nearest grid cell data
        d[varl[iv]]=np.nan*np.ones(x.size)
        for i in range(x.size):
            d[varl[iv]][i]=z['Data'][iy[i],ix[i]]

        # Covernt to float and apply scale factor
        d[varl[iv]]=d[varl[iv]].astype('float')*sf
    
    return d