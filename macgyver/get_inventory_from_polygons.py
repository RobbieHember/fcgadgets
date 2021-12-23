'''
PREPARE INVENTORY FROM POLYGONS
'''

#%% Import modules

import numpy as np
import pandas as pd
import geopandas as gpd
import fiona
import matplotlib.pyplot as plt
import numpy.matlib as ml
from shapely.geometry import Polygon,Point
import time
import gc as garc
from fcgadgets.macgyver import utilities_general as gu
from fcgadgets.macgyver import utilities_gis as gis
from fcgadgets.macgyver import utilities_inventory as invu
from fcgadgets.cbrunner import cbrun_utilities as cbu

#%% Project name

#project_name='FCI_RollupFCI_Inv'
#project_name='SummaryNutrientManagement'
#project_name='SummaryNutrientManagementFull'
#project_name='SummaryReforestationNonOb'
project_name='SummaryReforestation'
#project_name='SummarySurvey'

#%% Define paths

meta={}
meta['Paths']={}
meta['Paths']['Model Code']=r'C:\Users\rhember\Documents\Code_Python\fcgadgets\cbrunner'
meta['Paths']['Project']=r'D:\Data\FCI_Projects' + '\\' + project_name
#meta['Paths']['Project']=r'C:\Users\rhember\Documents\Data\FCI_Projects' + '\\' + project_name
meta['Paths']['Geospatial']=meta['Paths']['Project'] + '\\Geospatial'
meta['Paths']['Results']=r'C:\Users\rhember\Documents\Data\ForestInventory\Results\20210930'
meta['Paths']['VRI']=r'C:\Users\rhember\Documents\Data\ForestInventory\VRI\20210930'
meta['Paths']['Disturbances']=r'C:\Users\rhember\Documents\Data\ForestInventory\Disturbances\20210930'
meta['Paths']['LandUse']=r'C:\Users\rhember\Documents\Data\ForestInventory\LandUse\20210930'
meta['Paths']['Taz Datasets']=r'C:\Users\rhember\Documents\Data\Taz Datasets'

#%% Define subsampling frequency
# Some projects are way too big to collect 1-hectare coverage - subsample randomly

if project_name=='FCI_RollupFCI_Inv':    
    
    meta['subsampling_frequency']=0.25
    
    #FES recipients sometimes use the funding source code, "FES", for FCI-funded
    #projects. To include them in the query, import the FCI project list.    
    meta['Paths']['FCI DB File']=r'C:\Users\rhember\Documents\Data\FCI_Projects\FCI_RollupFCI_Inv\FCI_RollupProjects_01_Admin.xlsx'
    dAdmin=gu.ReadExcel(meta['Paths']['FCI DB File'])

    # Unique PP numbers in the database
    uPP=np.unique(dAdmin['PP Number'])
    
elif project_name=='SummaryNutrientManagement':
    
    # Import AIL- used to subsample certain multipolygons
    ail=gu.ipickle(r'D:\Data\FCI_Projects\SummaryNutrientManagement\Inputs\AnnualImplementationLevel.pkl')
    
    # Sparse grid subsampling rate
    meta['subsampling_frequency']=0.1
    
elif project_name=='SummaryNutrientManagementFull':
    
    # Import AIL- used to subsample certain multipolygons
    #ail=gu.ipickle(r'D:\Data\FCI_Projects\NutrientManagementSummary\Inputs\AnnualImplementationLevel.pkl')
    meta['subsampling_frequency']=0.05

elif project_name=='SummaryReforestationNonOb':
    
    # Import AIL- used to subsample certain multipolygons
    # *** Not subsampling - its too hard to partition into different project types afterwards ***
    #ail=gu.ipickle(r'D:\Data\FCI_Projects\SummaryReforestationNonOb\Inputs\AnnualImplementationLevel.pkl')
    
    # Sparse grid subsampling rate
    meta['subsampling_frequency']=0.05

elif project_name=='SummaryReforestation':  
    
    # Get random subsample
    ail=gu.ipickle(r'D:\Data\FCI_Projects\SummaryReforestation\Inputs\AnnualImplementationLevel.pkl')
    
    # Sparse grid subsampling rate
    meta['subsampling_frequency']=0.01
    
elif project_name=='SummarySurvey':  
    
    # Sparse grid subsampling rate
    meta['subsampling_frequency']=0.01

else:    
    
    meta['subsampling_frequency']=0.05

#%% Save metadata

gu.opickle(meta['Paths']['Project'] + '\\Inputs\\Metadata.pkl',meta)

#%% Define sparse grid based on geotiff from BC1ha database

# Import TSA maps
zTSA=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Admin\tsa.tif')
lut_tsa=gu.ReadExcel(r'C:\Users\rhember\Documents\Data\BC1ha\Admin\lut_tsa.xlsx')

# Grid interval (m)
ivl=zTSA['gt'][1]

# Domain (spans all of BC)
x_bc_all=zTSA['X'][0,:]
y_bc_all=zTSA['Y'][:,0]

# Load basemap (for CRS)
gdf_bm=gpd.read_file(r'C:\Users\rhember\Documents\Data\Basemaps\bc_land.shp')

#%% Compile inventory data at sparse grid locations

# Get layer names and variables
InvLyrInfo=invu.DefineInventoryLayersAndVariables()

#%% Import look-up tables for each layer

for iLyr in range(len(InvLyrInfo)):
    lut=gu.ipickle(InvLyrInfo[iLyr]['Path'] + '\\LUTs_' + InvLyrInfo[iLyr]['Layer Name'] +'.pkl')
    for key in lut.keys():
        InvLyrInfo[iLyr]['LUT'][key]=lut[key]
    del lut

#%% Open crosswalk between missing AT geometries and opening geometries
# If this doesn't work, you need to run the script that creates the crosswalk

missing_geo_atu_list=gu.ipickle(meta['Paths']['Results'] + '\\missing_geo_atu_list.pkl')
missing_geo_op_geos=gu.ipickle(meta['Paths']['Results'] + '\\missing_geo_op_geos.pkl')
missing_geo_fc_geos=gu.ipickle(meta['Paths']['Results'] + '\\missing_geo_fc_geos.pkl')

#%% Open missing FC layer geometries that were retreived from VRI
# Forest cover spatial reporting didn’t get turned on until 2004ish when RESULTS got turned on.  
# Previous forest cover reporting was delivered to the government as paper or as a PDF then 
# digitized directly into the VRI (or predecessors).  

dMisFC=gu.ipickle(meta['Paths']['Results'] + '\\missing_geo_fc_list.pkl')

#%% Get ATU lyaer path and layer name

# Get index
for iLyr in range(len(InvLyrInfo)):
    if InvLyrInfo[iLyr]['Layer Name']=='RSLT_ACTIVITY_TREATMENT_SVW':
        break

# Define path
path=InvLyrInfo[iLyr]['Path'] + '\\' + InvLyrInfo[iLyr]['File Name']

# Define layer name
lyr_nam=InvLyrInfo[iLyr]['Layer Name']
print(lyr_nam)
    
#%% Initialize list for multipolygons
# Keep multipolygons (the only way to get accurate area due to mismatches between
# opening area and treatment area)

atu_multipolygons=[None]*5000000
cnt_atu_multipolygons=0

# Track those that with no spatial as well
atu_mp_miss=[None]*5000000
cnt_atu_mp_miss=0

#%% Initialize a geodataframe that will store each polygon

cnt_atu_polygons=0
cn_atu_polygons=['ID_atu_polygons','ID_atu_multipolygons','MissStatus','Year','geometry']
list_atu_polygons=[None]*5000000
    
#%% Initialize sparse grid coordinate dictionary

L=20000000
sxy={}
sxy['x']=-999*np.ones(L)
sxy['y']=-999*np.ones(L)
sxy['ID_atu_multipolygons']=-999*np.ones(L,dtype=int)
sxy['ID_atu_polygons']=-999*np.ones(L,dtype=int)
cnt_sxy=0

#%% Get polygons from file

# Loop through features in layer
t0=time.time()

# Track amount of missing spatial
N_Missing_Spatial=0

trip=0

cnt=0

with fiona.open(path,layer=lyr_nam) as source:
        
    for feat in source:
            
        # Extract attributes and geometry
        prp=feat['properties']
        
        #if (prp['OPENING_ID']==1462976) & (prp['SILV_BASE_CODE']=='PL'):
        #    break
        
        geom=feat['geometry']
        
        # In rare instances, there are no dates - add dummy numbers so that it 
        # does not crash - it appears to be so rare that it should not be a big
        # problem
        if prp['ATU_COMPLETION_DATE']==None:
            prp['ATU_COMPLETION_DATE']='99990000'
        
        Year=int(prp['ATU_COMPLETION_DATE'][0:4]) 
        
        #----------------------------------------------------------------------
        # Project-specific query
        #----------------------------------------------------------------------
        
        if project_name=='SummaryNutrientManagement':
            
            flg_stop=1
            if (prp['SILV_BASE_CODE']=='FE') & (prp['SILV_TECHNIQUE_CODE']=='CA') & (prp['SILV_METHOD_CODE']=='HELI') & (prp['RESULTS_IND']=='Y') & (prp['ACTUAL_TREATMENT_AREA']!=None) & (prp['ATU_COMPLETION_DATE']!=None) & (prp['SILV_FUND_SOURCE_CODE']!=None):
                flg_stop=0
            
            if flg_stop==1:
                continue
            
            # Do subsampling to save time
            if np.isin(prp['ACTIVITY_TREATMENT_UNIT_ID'],ail['id_atu_subsample'])==False:
                continue
        
        elif project_name=='SummaryNutrientManagementFull':
            
            flg_stop=1
            if (prp['SILV_BASE_CODE']=='FE') & (prp['SILV_TECHNIQUE_CODE']=='CA') & (prp['SILV_METHOD_CODE']=='HELI') & (prp['RESULTS_IND']=='Y') & (prp['ACTUAL_TREATMENT_AREA']!=None) & (prp['ATU_COMPLETION_DATE']!=None) & (prp['SILV_FUND_SOURCE_CODE']!=None):
                flg_stop=0
            
            if flg_stop==1:
                continue
            
        elif project_name=='SummaryReforestationNonOb':
            
            # Define which funding source codes are licensee vs. non-ob
            ListOfNonObFSC=['FTL','FTM','RBM','RBL','FR','VG','FIL','FID','FIM','S', \
                    'FRP','XXX','O','GFS','IR','FES','FCE','FCM']
            ListOfLicenseeFSC=['BCT','LFP''IA','IR','VOI','SBF']
            
            # Only include certain years
            if Year<1990:
                continue
            
            # To get actual planting or direct seeding, exclude:
            #   seedling protection ('SILV_TECHNIQUE_CODE'=='SE')
            #   layout ('SILV_METHOD_CODE'=='LAYOT')
            #   fertilization ('SILV_TECHNIQUE_CODE'=='CG')            
            flg=0
            if (prp['SILV_BASE_CODE']=='PL') & (prp['SILV_TECHNIQUE_CODE']!='SE') & (prp['SILV_TECHNIQUE_CODE']!='CG') & (prp['SILV_METHOD_CODE']!='LAYOT') & (prp['RESULTS_IND']=='Y') & (np.isin(prp['SILV_FUND_SOURCE_CODE'],ListOfNonObFSC)==True):
                flg=1
            #if (prp['SILV_BASE_CODE']=='DS') & (prp['SILV_METHOD_CODE']!='LAYOT') & (prp['RESULTS_IND']=='Y') & (np.isin(prp['SILV_FUND_SOURCE_CODE'],ListOfNonObFSC)==True):
            #    flg=1
            
            # Do subsampling to save time
            #if np.isin(prp['ACTIVITY_TREATMENT_UNIT_ID'],ail['id_atu_subsample'])==False:
            #    continue
            
            if flg==0:
                continue
        
        elif project_name=='SummaryReforestation':
            
            flg=0
            if (prp['SILV_BASE_CODE']=='PL') & (prp['SILV_TECHNIQUE_CODE']!='SE') & (prp['SILV_TECHNIQUE_CODE']!='CG') & (prp['SILV_METHOD_CODE']!='LAYOT') & (prp['RESULTS_IND']=='Y') & (prp['ACTUAL_TREATMENT_AREA']!=None) & (prp['ATU_COMPLETION_DATE']!=None):
                flg=1            
            if flg==0:
                continue
            
            # Only include certain years
            if Year<1990:
                continue
            
            # Do subsampling to save time
            if np.isin(prp['ACTIVITY_TREATMENT_UNIT_ID'],ail['id_atu_subsample'])==False:
                continue
        
        elif project_name=='FESBC':
            
            if (prp['RESULTS_IND']!='Y') | (prp['SILV_FUND_SOURCE_CODE']!='FES') | (prp['SILV_BASE_CODE']=='SU'):
                continue   
        
        elif project_name=='SummarySurvey':
            
            Year=int(prp['ATU_COMPLETION_DATE'][0:4])   
            
            if (prp['SILV_BASE_CODE']!='SU'):
                continue
            if (Year<2018):
                continue
        
        elif project_name=='FCI_RollupFCI_Inv':
            
            # FCI query
            if prp['RESULTS_IND']=='N':
                continue
            if (prp['SILV_FUND_SOURCE_CODE']!='FCE') & (prp['SILV_FUND_SOURCE_CODE']!='FCM') & (np.isin(prp['FIA_PROJECT_ID'],uPP)==False):
                continue
            if (prp['FIA_PROJECT_ID']=='RB0000320') | (prp['FIA_PROJECT_ID']=='FC0000313') | (prp['FIA_PROJECT_ID']=='WR0000011'):
                continue
            if prp['SILV_BASE_CODE']=='SU':
                continue
        else:
            print('Project name unspecified!')
            continue
        
        #cnt=cnt+1
        print(cnt_atu_multipolygons)
        
        # Populate missing ATU layer geometry with geometry from 
        # OPENING or FC layer where possible.
        flg_geom_from_op=0
        flg_geom_from_fc=0
        flg_geom_from_cut=0
        if (geom==None):               
                    
            # Check to see if the opening is listed in the AT missing dictionary
            indMis=np.where( (missing_geo_atu_list['ACTIVITY_TREATMENT_UNIT_ID']==prp['ACTIVITY_TREATMENT_UNIT_ID']) )[0]
            
            if indMis.size>0:
                
                idx2fc=missing_geo_atu_list['IdxToFC'][indMis[0]]
                
                if len(idx2fc)>0:
                    
                    # Use forest cover geometries
                
                    geom={}
                    geom['coordinates']=[]
                    for i in range(len(idx2fc)):
                        geo0=missing_geo_fc_geos[prp['OPENING_ID']][idx2fc[i]]
                        if type(geo0)==dict:
                            geo1=geo0['coordinates']
                            geom['coordinates'].append(geo1[0])
                        else:
                            for j in range(len(geo0)):
                                geo1=geo0[j]['coordinates']
                                geom['coordinates'].append(geo1[0])
                    
                    flg_geom_from_fc=1
                    print('Missing spatial recovered from forest cover layer')
                    
                    # Plot (not working)
                    #flg=0
                    #if flg==1:                        
                    #    plt.close('all')
                    #    fig,ax=plt.subplots(1)
                    #    gdf_fc=gpd.GeoDataFrame.from_features(feat_fc)
                    #    gdf_fc.plot(ax=ax,facecolor='None',edgecolor='r',linewidth=1.25,linestyle='--')          
                
                elif len(missing_geo_op_geos[prp['OPENING_ID']])>0:
                    
                    # Use opening geometry
                
                    geom={}
                    geom['coordinates']=[]
                    geo0=missing_geo_op_geos[prp['OPENING_ID']]
                    if type(geo0)==dict:
                        geo1=geo0['coordinates']
                        geom['coordinates'].append(geo1[0])
                    else:
                        for j in range(len(geo0)):
                            geo1=geo0[j]['coordinates']
                            geom['coordinates'].append(geo1[0])
                            
                    flg_geom_from_op=1
                    print('Missing spatial recovered from opening layer')
                    
            else:
                    
                # Could not use either FC or opening layer
                print('Missing spatial could not be recovered')
            
        # Don't conitnue if no spatial data
        if (geom==None): 
            
            N_Missing_Spatial=N_Missing_Spatial+1
            
            # Update counter for missing multipolygon
            prp['Year']=int(prp['ATU_COMPLETION_DATE'][0:4])
            atu_mp_miss[cnt_atu_mp_miss]=prp
            cnt_atu_mp_miss=cnt_atu_mp_miss+1
            
            continue
        
        prp['ID_atu_multipolygons']=cnt_atu_multipolygons
        prp['Year']=int(prp['ATU_COMPLETION_DATE'][0:4])
        prp['geometry']=geom
        prp['GeomFromOpLyr']=flg_geom_from_op
        prp['GeomFromFcLyr']=flg_geom_from_fc
        prp['GeomFromCutLyr']=flg_geom_from_cut
        
        # Extract vector geometries and store as a list of dictionaries.
        # The list contains each polygon ("block"). Within each block, a
        # dictionary stores the x and y coordinates and the area.
        polys=[]
        polys_inner=[]
        coords0=geom['coordinates']
        for iPoly in range(len(coords0)):
            coords1=coords0[iPoly]
            
            x=[]; 
            y=[];
            for k in range(len(coords1[0])):
                x.append(coords1[0][k][0])
                y.append(coords1[0][k][1])
            
            l=[]
            for k in range(len(x)):
                l.append([x[k],y[k]])
            lp=Polygon(l)             
                    
            A=gis.PolyArea(x,y)/10000
                
            dct={}
            dct['Hectares']=A.copy()
            dct['x']=np.array(x).copy()
            dct['y']=np.array(y).copy()
            dct['x_Centroid']=lp.centroid.x
            dct['y_Centroid']=lp.centroid.y   
            dct['ID_atu_polygons']=cnt_atu_polygons
            polys.append(dct.copy())
            
            # Add to geodataframe
            dp0={}
            dp0['ID_atu_polygons']=cnt_atu_polygons
            dp0['ID_atu_multipolygons']=cnt_atu_multipolygons            
            dp0['OPENING_ID']=prp['OPENING_ID']
            dp0['ACTIVITY_TREATMENT_UNIT_ID']=prp['ACTIVITY_TREATMENT_UNIT_ID']
            dp0['Year']=prp['Year']
            dp0['SILV_FUND_SOURCE_CODE']=prp['SILV_FUND_SOURCE_CODE']
            dp0['FIA_PROJECT_ID']=prp['FIA_PROJECT_ID'] 
            dp0['SILV_BASE_CODE']=prp['SILV_BASE_CODE']
            dp0['SILV_TECHNIQUE_CODE']=prp['SILV_TECHNIQUE_CODE']
            dp0['SILV_METHOD_CODE']=prp['SILV_METHOD_CODE']
            dp0['GeomFromOpLyr']=prp['GeomFromOpLyr']
            dp0['GeomFromFcLyr']=prp['GeomFromFcLyr']
            dp0['GeomFromCutLyr']=prp['GeomFromCutLyr']
#            dp0['MissStatus']=0
#            if prp['GeomFromFcLyr']==1:
#                dp0['MissStatus']=1
#            if prp['GeomFromOpLyr']==1:
#                dp0['MissStatus']=2    
            #dp0['ACTUAL_TREATMENT_AREA']=prp['ACTUAL_TREATMENT_AREA']            
            #dp0['Year']=prp['Year']
            
            gdf_outer=gpd.GeoDataFrame(crs=gdf_bm.crs)
            gdf_outer.loc[0,'geometry']=Polygon(coords1[0])
            if len(coords1)>1:
                # Adjust polygons if they have an inner ring
                for j in range(1,len(coords1)):
                    gdf_inner=gpd.GeoDataFrame(crs=gdf_bm.crs)
                    gdf_inner.loc[0,'geometry']=Polygon(coords1[j])                
                    try:
                        # This very rarely receives "'NoneType' object has no attribute 'intersection'"
                        # Not sure why
                        gdf_outer=gpd.overlay(gdf_outer,gdf_inner,how='difference')
                    
                        x=[]; 
                        y=[];
                        for k in range(len(coords1[j])):
                            x.append(coords1[j][k][0])
                            y.append(coords1[j][k][1])
                        dct={}
                        dct['x']=np.array(x).copy()
                        dct['y']=np.array(y).copy()
                        polys_inner.append(dct.copy())
                    except:
                        pass
            
            dp0['geometry']=gdf_outer.loc[0,'geometry']
            list_atu_polygons[cnt_atu_polygons]=dp0
            cnt_atu_polygons=cnt_atu_polygons+1        
            
        # Add polygon info to atu query polygon list
        prp['polys']=polys
        
        #for i in range(len(polys)):
        #    plt.plot(polys[i]['x'],polys[i]['y'],'b--')
        
        # Extract grid cells that intersect the AT polygons
        for iPoly in range(len(polys)):
            
            if polys[iPoly]['Hectares']==None:
                print('Missing geometry')
                continue
            
            x_feat=polys[iPoly]['x']
            y_feat=polys[iPoly]['y']
        
            # Index to a bounding box around the feature geometry
            ix=np.where((x_bc_all>=np.min(x_feat)-2*ivl) & (x_bc_all<=np.max(x_feat)+2*ivl))[0]
            iy=np.where((y_bc_all>=np.min(y_feat)-2*ivl) & (y_bc_all<=np.max(y_feat)+2*ivl))[0]
        
            if (ix.size==0) | (iy.size==0):
                continue
    
            x_bb=np.reshape(x_bc_all[ix],(1,ix.shape[0]))
            y_bb=np.reshape(y_bc_all[iy],(iy.shape[0],1))
            x_bb=np.tile(x_bb,(y_bb.shape[0],1))
            y_bb=np.tile(y_bb,(1,x_bb.shape[1]))
            y_bb=np.flip(y_bb,0)
        
            # Indicator for whether a cell is within the geometry
            InPol=gis.InPolygon(x_bb,y_bb,x_feat,y_feat)
          
            # Index to cells within feature geometry
            ikp=np.where(InPol==1)
            
            # If grid cell(s) fall within the project area, use those grid cells
            if ikp[0].size!=0:
                
                # Subsampling
                if meta['subsampling_frequency']<1:
                    Nss=np.maximum(1,np.ceil(meta['subsampling_frequency']*ikp[0].size)).astype(int)
                    iSS=np.random.randint(ikp[0].size,size=Nss)   
                    ikp=(ikp[0][iSS],ikp[1][iSS])
                
                x_ikp=np.atleast_1d(x_bb[ikp])
                y_ikp=np.atleast_1d(y_bb[ikp])
                
                # Remove cells within inner rings
                for iPoly_inner in range(len(polys_inner)):                    
                    InPol=gis.InPolygon(x_ikp,y_ikp,polys_inner[iPoly_inner]['x'],polys_inner[iPoly_inner]['y'])
                    ind=np.where(InPol==0)[0]
                    if ind.size>0:
                        x_ikp=x_ikp[ind]
                        y_ikp=y_ikp[ind]
                
                for i in range(x_ikp.size):
                    sxy['x'][cnt_sxy]=x_ikp[i]
                    sxy['y'][cnt_sxy]=y_ikp[i]                                
                    sxy['ID_atu_multipolygons'][cnt_sxy]=cnt_atu_multipolygons
                    sxy['ID_atu_polygons'][cnt_sxy]=prp['polys'][iPoly]['ID_atu_polygons']
                    cnt_sxy=cnt_sxy+1
            else:
                # If no grid cells fall within the project area, don't use the 
                # polygon centriod (doesn't work) - make a higher res grid and
                # re-search for interior cells.
                # (This occurs in road plantings)
                # Go all the way down to 5 m, 10 m misses some small features
                bin=[10,5,1,0.5]
                for iBin in range(len(bin)):                
                    #print('Going down to higher res!')
                    x_bb2=np.arange(np.min(x_feat),np.max(x_feat),bin[iBin])
                    y_bb2=np.arange(np.min(y_feat),np.max(y_feat),bin[iBin])
                    x_bb2=np.matlib.repmat(np.reshape(x_bb2,(1,x_bb2.size)),y_bb2.size,1)
                    y_bb2=np.matlib.repmat(np.reshape(y_bb2,(y_bb2.size,1)),1,x_bb2.size)
                    InPol=gis.InPolygon(x_bb2,y_bb2,x_feat,y_feat)
                    ikp=np.where(InPol==1)
                    if ikp[0].size>0:
                        break
                
                r=np.random.permutation(ikp[0].size)
                x_ikp=x_bb2[ikp]
                y_ikp=y_bb2[ikp]
                for iR in range(np.minimum(ikp[0].size,5)):
                    x_ikp0=np.atleast_1d(x_ikp[r[iR]])
                    y_ikp0=np.atleast_1d(y_ikp[r[iR]])
                  
                    # Remove cells within inner rings
                    for iPoly_inner in range(len(polys_inner)):
                        InPol=gis.InPolygon(x_ikp0,y_ikp0,polys_inner[iPoly_inner]['x'],polys_inner[iPoly_inner]['y'])
                        ind=np.where(InPol==0)[0]
                        if ind.size>0:
                            x_ikp0=x_ikp0[ind]
                            y_ikp0=y_ikp0[ind]
                    
                    sxy['x'][cnt_sxy]=x_ikp0
                    sxy['y'][cnt_sxy]=y_ikp0
                    sxy['ID_atu_multipolygons'][cnt_sxy]=cnt_atu_multipolygons
                    sxy['ID_atu_polygons'][cnt_sxy]=prp['polys'][iPoly]['ID_atu_polygons']
                    cnt_sxy=cnt_sxy+1
                
        # Update counter for multipolygon
        atu_multipolygons[cnt_atu_multipolygons]=prp
        cnt_atu_multipolygons=cnt_atu_multipolygons+1

#%% Truncate data

# Vector geometries
atu_multipolygons=atu_multipolygons[0:cnt_atu_multipolygons]

list_atu_polygons=list_atu_polygons[0:cnt_atu_polygons]

# Sparse sample points
for k in sxy.keys():
    sxy[k]=sxy[k][0:cnt_sxy]

#%% Remove duplicate cells
# Notes: This means that some values of ID_atu_multipolygons may not exist in
# the sxy dictionaries - they were only retained for one of the IDs of the 
# overlapping multipolygons. 
              
# Unique cells (ind_xy is the index to the first instance of each unique row)
uxy,ind_xy,inv_xy=np.unique(np.column_stack((sxy['x'],sxy['y'])),return_index=True,return_inverse=True,axis=0)

for k in sxy.keys():
    sxy[k]=sxy[k][ind_xy]

#%% Check lengths and make note of (expected) discrepency in length between unique
# IDs from sxy dictionary and unique IDs from list of multipolygons

flg=0
if flg==1:
    atu_multipolygons=gu.ipickle(meta['Paths']['Geospatial'] + '\\atu_multipolygons.pkl')
    print(len(atu_multipolygons))
    a=[]
    for i in range(len(atu_multipolygons)):
        if atu_multipolygons[i]!=None:
            a.append(atu_multipolygons[i]['ID_atu_multipolygons'])
    u_mp=np.unique(np.array(a))
    print(u_mp.size)

    u_sxy=np.unique(sxy['ID_atu_multipolygons'])
    print(u_sxy.size)

    c=0
    for i in range(len(atu_multipolygons)):
        if atu_multipolygons[i]['OPENING_ID']==1702534:        
            for j in range(len(atu_multipolygons[i]['polys'])):
                plt.plot(atu_multipolygons[i]['polys'][j]['x'],atu_multipolygons[i]['polys'][j]['y'],'b-')
                #c=c+1
            ind=np.where(sxy['ID_atu_multipolygons']==atu_multipolygons[i]['ID_atu_multipolygons'])[0]
            c=c+ind.size
            plt.plot(sxy['x'][ind],sxy['y'][ind],'.')    

#%% Add TSA ID

sxy['ID_TSA']=0*sxy['x']
for i in range(sxy['x'].size):
    ix=np.where(zTSA['X'][0,:]==sxy['x'][i])[0]
    iy=np.where(zTSA['Y'][:,0]==sxy['y'][i])[0]
    id_tsa=zTSA['Data'][iy,ix]
    if id_tsa.size>0:
        sxy['ID_TSA'][i]=id_tsa
    #ind=np.where(lut_tsa['VALUE']==id_tsa)[0]
    #lut_tsa['Name'][ind]

#ind=np.where(lut_tsa['Name']=='Okanagan TSA')[0]

#%% Save multipolygon to file
    
# Save multipolygon file
gu.opickle(meta['Paths']['Project'] + '\\Geospatial\\atu_multipolygons.pkl',atu_multipolygons)

# Save multipolygon to geojson (this crashes - hard to fix)

flg=0
if flg==1:

    atu_multipolygons=gu.ipickle(meta['Paths']['Project'] + '\\Geospatial\\atu_multipolygons.pkl')
    
    List=[None]*len(atu_multipolygons)
    for i in range(len(atu_multipolygons)):
        if atu_multipolygons[i]==None:
            continue
        d={}
        for k in atu_multipolygons[i].keys():
            if (k!='polys'):
                d[k]=atu_multipolygons[i][k]
        List[i]=d
    
    df=pd.DataFrame(List[0])
    #gdf=gpd.GeoDataFrame(df,geometry=gpd.points_from_xy(df.Longitude, df.Latitude))
    gdf=gpd.GeoDataFrame(df,geometry=df.geometry)
    
    
    
    List=[None]*len(atu_multipolygons)
    for i in range(len(atu_multipolygons)):
        if atu_multipolygons[i]==None:
            continue
        d={}
        d['geometry']=atu_multipolygons[i]['geometry']
        prp={}
        for k in atu_multipolygons[i].keys():
            if (k!='geometry') & (k!='polys'):
                prp[k]=atu_multipolygons[i][k]
        d['properties']=prp
        List[i]=d
    
    gdf=gpd.GeoDataFrame.from_features(List,crs=gdf_bm.crs)

#%% Save polygons to file

# Save to pickle
gu.opickle(meta['Paths']['Project'] + '\\Geospatial\\atu_polygons.pkl',list_atu_polygons)

# Save to geojson
gdf_atu_polygons=gpd.GeoDataFrame(list_atu_polygons,crs=gdf_bm.crs)
gdf_atu_polygons=gdf_atu_polygons.set_geometry('geometry')
#gdf_atu_polygons=gdf_atu_polygons.to_crs({'init':'epsg:4326'})
gdf_atu_polygons.to_file(filename=meta['Paths']['Project'] + '\\Geospatial\\atu_polygons.geojson',driver='GeoJSON')
#gdf=gpd.read_file(meta['Paths']['Project'] + '\\Geospatial\\atu_polygons.geojson')

# Save shape file
fiona.supported_drivers
fnam=meta['Paths']['Project'] + '\\Geospatial\\atu_polygons.shp'
gdf_atu_polygons.to_file(filename=fnam,driver='ESRI Shapefile')

# Convert to geographic coordinates and then save to shapefile
#meta=gu.ipickle(r'C:\Users\rhember\Documents\Data\FCI_Projects\FCI_RollupFCI_Inv\Inputs\Metadata.pkl')
#gdf=gpd.read_file(meta['Paths']['Project'] + '\\Geospatial\\atu_polygons.geojson')
#gdf_atu_polygons.to_file(filename=meta['Paths']['Project'] + '\\Geospatial\\atu_polygons.geojson')
#gdf=gdf.to_crs({'init':'epsg:4326'})
#fiona.supported_drivers
#fnam=meta['Paths']['Project'] + '\\Geospatial\\atu_polygons_latlon.shp'
#gdf.to_file(filename=fnam,driver='ESRI Shapefile')

#%% Save sparse grid sample to file

# Add opening id to SXY dictionary before saving
sxy['OPENING_ID']=np.zeros(sxy['ID_atu_multipolygons'].size)
for iMP in range(len(atu_multipolygons)):
    ind=np.where(sxy['ID_atu_multipolygons']==iMP)[0]
    sxy['OPENING_ID'][ind]=atu_multipolygons[iMP]['OPENING_ID']

# Save sparse sample file
gu.opickle(meta['Paths']['Project'] + '\\Geospatial\\sxy.pkl',sxy)

# Save sparse grid as shapefile
flg=1
if flg==1:

    points=[]
    for k in range(sxy['x'].size):
        points.append(Point(sxy['x'][k],sxy['y'][k]))
        
    gdf_sxy=gpd.GeoDataFrame( {
        'geometry':points,
        'ID_atu_multipolygons':sxy['ID_atu_multipolygons'],
        'ID_atu_polygons':sxy['ID_atu_polygons'],
        'OPENING_ID':sxy['OPENING_ID'],
        'ID_TSA':sxy['ID_TSA']} )
    
    gdf_sxy=gdf_sxy.set_geometry('geometry')
    gdf_sxy.crs=gdf_bm.crs   
    gdf_sxy.to_file(meta['Paths']['Project'] + '\\Geospatial\\sxy.geojson',driver='GeoJSON')

#gdf_atu_polygons=gdf_atu_polygons.set_geometry('geometry')
# #gdf_atu_polygons=gdf_atu_polygons.to_crs({'init':'epsg:4326'})
#gdf_atu_polygons.to_file(filename=meta['Paths']['Project'] + '\\Geospatial\\atu_polygons.geojson',driver='GeoJSON')

print((time.time()-t0)/60)

#%% Extract inventory information for sparse grid sample over polygons
# *** See new removal of SU and incomplete entries from AT layer ***

garc.collect()

atu_multipolygons=gu.ipickle(meta['Paths']['Project'] + '\\Geospatial\\atu_multipolygons.pkl')
sxy=gu.ipickle(meta['Paths']['Project'] + '\\Geospatial\\sxy.pkl')

for iLyr in range(len(InvLyrInfo)):
    
    # Define path
    path=InvLyrInfo[iLyr]['Path'] + '\\' + InvLyrInfo[iLyr]['File Name']
    
    # Define layer
    lyr_nam=InvLyrInfo[iLyr]['Layer Name']
    print(lyr_nam)
    
    # Don't run this for planting - planting is done seperately below
    if lyr_nam=='RSLT_PLANTING_SVW':
        continue
     
    # Initialize index to inventory
    IdxToInv=[None]*sxy['x'].size

    # Initialize inventory dictionary
    L=20*sxy['x'].size
    data={}
    data['IdxToSXY']=np.zeros(L,dtype=int)
    for fnam,flag,dtype in InvLyrInfo[iLyr]['Field List']:
        if dtype=='<U20':
            data[fnam]=np.zeros(L,dtype=dtype)
        else:
            data[fnam]=-999*np.ones(L,dtype=dtype)
    cnt_inventory=0
    
    # Loop through features in layer
    t_start=time.time()

    # Scan through layer file to extract selected variables, and convert string
    # variables to numeric based on LUTs
    with fiona.open(path,layer=lyr_nam) as source:
        
        for feat in source:
            
            # Extract attributes and geometry
            prp=feat['properties']
            geom=feat['geometry']
            
            # If AT layer, may not be much need in retaining planned or SU
            if (lyr_nam=='RSLT_ACTIVITY_TREATMENT_SVW'):
                if (prp['RESULTS_IND']=='N') | (prp['SILV_BASE_CODE']=='SU'):
                    continue
            
            # Populate missing ATU layer geometry with geometry from 
            # OPENING or FC layer where possible.
            flg_geom_from_op=0
            flg_geom_from_fc=0
            if (lyr_nam=='RSLT_ACTIVITY_TREATMENT_SVW') & (geom==None):
                    
                # Check to see if the opening is listed in the AT missing dictionary
                indMis=np.where( (missing_geo_atu_list['ACTIVITY_TREATMENT_UNIT_ID']==prp['ACTIVITY_TREATMENT_UNIT_ID']) )[0]
            
                if indMis.size>0:
                
                    idx2fc=missing_geo_atu_list['IdxToFC'][indMis[0]]
                
                    if len(idx2fc)>0:
        
                        # Use forest cover geometries
                
                        geom={}
                        geom['coordinates']=[]
                        for i in range(len(idx2fc)):
                            geo0=missing_geo_fc_geos[prp['OPENING_ID']][idx2fc[i]]
                            if type(geo0)==dict:
                                geo1=geo0['coordinates']
                                geom['coordinates'].append(geo1[0])
                            else:
                                for j in range(len(geo0)):
                                    geo1=geo0[j]['coordinates']
                                    geom['coordinates'].append(geo1[0])
                    
                        flg_geom_from_fc=1
                    
                    #elif prp['OPENING_ID'] in missing_geo_op_geos==True:
                    
                    elif len(missing_geo_op_geos[prp['OPENING_ID']])>0:
                    
                        # Use opening geometry
                
                        geom={}
                        geom['coordinates']=[]
                        geo0=missing_geo_op_geos[prp['OPENING_ID']]
                        if type(geo0)==dict:
                            geo1=geo0['coordinates']
                            geom['coordinates'].append(geo1[0])
                        else:
                            for j in range(len(geo0)):
                                geo1=geo0[j]['coordinates']
                                geom['coordinates'].append(geo1[0])

                        flg_geom_from_op=1
                    
                else:
                
                    # Could not use either FC or openign layer
                    #print('Missing spatial could not be recovered')
                    pass
            
            # Populate missing FC layer geometry with geometry from 
            # OPENING or VRI where possible.
            if (lyr_nam=='RSLT_FOREST_COVER_INV_SVW') & (geom==None) | (lyr_nam=='RSLT_FOREST_COVER_SILV_SVW') & (geom==None):
                
                iMis_fc=np.where( (dMisFC['Unique Openings with Missing FC Geom']==prp['OPENING_ID']) )[0]
                iMis_fc=iMis_fc[0]
                
                if dMisFC['Geom from VRI'][iMis_fc]!=None:
                
                    D=np.zeros(len(dMisFC['Geom from VRI'][iMis_fc]))
                    for iV in range(len(dMisFC['Geom from VRI'][iMis_fc])):
                        D[iV]=np.abs(prp['SILV_POLYGON_AREA']-dMisFC['Geom from VRI'][iMis_fc][iV]['Hectares'])
                    iMinD=np.where(D==np.min(D))[0]
                    iMinD=iMinD[0]
            
                    geom=dMisFC['Geom from VRI'][iMis_fc][iMinD]
                
                    # QA: Look at polygon
                    #geom1=invu.GetPolygonsFromFionaFeature(geom)                
                    #gdf=gpd.GeoDataFrame(geom1,crs=gdf_bm.crs)            
            
            # Only continue if spatial info exists
            if (geom==None) | (geom==[]):
                continue
            
            # Extract multipolygon
            coords0=geom['coordinates']
            #plt.close('all')
            # loop through multipolygon
            for i in range(len(coords0)):
                
                # Extract multipolygon
                coords1=coords0[i]
                
                # loop through multipolygon
                for j in range(len(coords1)):
                    
                    # Extract polygon
                    coords2=np.asarray(coords1[j])
                    x_feat=coords2[:,0]
                    y_feat=coords2[:,1]
                    #plt.plot(x_feat,y_feat,'-',linewidth=5)
                    
                    # This should speed it up a bit
                    if np.max(x_feat)<np.min(sxy['x']):
                        continue
                    if np.min(x_feat)>np.max(sxy['x']):
                        continue
                    if np.max(y_feat)<np.min(sxy['y']): 
                        continue
                    if np.min(y_feat)>np.max(sxy['y']): 
                        continue
        
                    # Isolate cells within bounding box of the feature polygon
                    iBB=np.where( (sxy['x']>=np.min(x_feat)-1000) & (sxy['x']<=np.max(x_feat)+1000) & (sxy['y']>=np.min(y_feat)-1000) & (sxy['y']<=np.max(y_feat)+1000) )[0]
                    
                    # Only continue if overlaop
                    if iBB.size==0:
                        continue
        
                    # Isolate cells within the polygon
                    InPoly=gis.InPolygon(sxy['x'][iBB],sxy['y'][iBB],x_feat,y_feat)
                    
                    iInPoly=np.where(InPoly==1)[0]
        
                    # Index to cells within polygon
                    iKeep=iBB[iInPoly]
                    
                    # Only continue if overlap
                    if iKeep.size==0:     
                        continue
                    
                    #print('working...')
                    
                    # Add attributes of overlapping grid cells to list
                    for k in range(iKeep.size):
                        
                        iKeepK=iKeep[k]
                        
                        # Update index to inventory
                        if IdxToInv[iKeepK]==None:
                            IdxToInv[iKeepK]={}
                            IdxToInv[iKeepK]['Index']=np.array([cnt_inventory])
                        else:
                            IdxToInv[iKeepK]['Index']=np.append(IdxToInv[iKeepK]['Index'],cnt_inventory)
                        
                        data['IdxToSXY'][cnt_inventory]=iKeepK
                        for fnam,flag,dtype in InvLyrInfo[iLyr]['Field List']:
                            val=feat['properties'][fnam]
                            if val!=None:
                                if flag==0:
                                    # Numeric variable, leave as is
                                    data[fnam][cnt_inventory]=val
                                elif flag==1:
                                    # Convert string variable to numeric variable 
                                    # based on LUT
                                    data[fnam][cnt_inventory]=InvLyrInfo[iLyr]['LUT'][fnam][val]
                
                        cnt_inventory=cnt_inventory+1
    
    # Time
    t_ela=time.time()-t_start
    print(t_ela)
    
    #--------------------------------------------------------------------------
    # Truncate data
    #--------------------------------------------------------------------------
    
    # Truncate index to inventory
    IdxToInv=IdxToInv[0:sxy['x'].size+1]
            
    # Truncate variables at cnt_inventory
    data['IdxToSXY']=data['IdxToSXY'][0:cnt_inventory]
    for fnam,flag,dtype in InvLyrInfo[iLyr]['Field List']:
        data[fnam]=data[fnam][0:cnt_inventory]
    
    # Replace time string with numeric variables
    # Note: Don't do this for planting - planting will be re-compiled below
    # to include planting info for projects that reported no geometries in the
    # planting layer. The Strings will be fixed then.
    if lyr_nam!='RSLT_PLANTING_SVW':
        data=invu.ExtractDateStringsFromRESULTS(lyr_nam,data)      
       
    #--------------------------------------------------------------------------    
    # Save to file
    #--------------------------------------------------------------------------
    
    gu.opickle(meta['Paths']['Project'] + '\\Geospatial\\' + InvLyrInfo[iLyr]['Layer Name'] + '.pkl',data)
    gu.opickle(meta['Paths']['Project'] + '\\Geospatial\\' + InvLyrInfo[iLyr]['Layer Name'] + '_IdxToInv.pkl',IdxToInv)
        
#%% Retrieve planting information 
# Some projects did not report spatial planting info. Without the spatial info
# in the planting layer, the initial import of the PL layer (above) will miss
# planting info in some AT polygons. Use the ACTIVITY TREATMENT UNIT ID as
# a crosswalk to retrieve all planting info for each activity.
# (10 min)

t0=time.time()

# Import coordinates   
sxy=gu.ipickle(meta['Paths']['Geospatial'] + '\\sxy.pkl')
 
# Get planting lyaer id
for iLyr in range(len(InvLyrInfo)):
    if InvLyrInfo[iLyr]['Layer Name']=='RSLT_PLANTING_SVW':
        break
    
# Import geodataframe, drop geometry and convert to dict
path=InvLyrInfo[iLyr]['Path'] + '\\' + InvLyrInfo[iLyr]['File Name']
lyr_nam=InvLyrInfo[iLyr]['Layer Name']
gdf_pl=gpd.read_file(path,layer=lyr_nam)        
d_pl=gu.DataFrameToDict(gdf_pl.drop(columns='geometry'))

#ind=np.where( d_pl['OPENING_ID']==np.float(-603550000) )[0]
#ind=np.where( d_pl['OPENING_ID']==-603550000 )[0]

# Get keys for planting layer
key_pl=[]
for fnam,flag,dtype in InvLyrInfo[iLyr]['Field List']:
    key_pl.append(fnam)
            
# Open AT sparse grid
atu=gu.ipickle(meta['Paths']['Geospatial'] + '\\RSLT_ACTIVITY_TREATMENT_SVW.pkl')

pl_code=InvLyrInfo[0]['LUT']['SILV_BASE_CODE']['PL']
    
# Only proceed if planting occurs
ind_at=np.where(atu['SILV_BASE_CODE']==pl_code)[0]

# Initialize dictionary
L=20*sxy['x'].size
pl={}
pl['IdxToSXY']=np.zeros(L,dtype=int)
for fnam,flag,dtype in InvLyrInfo[iLyr]['Field List']:
    if dtype=='<U20':
        pl[fnam]=np.zeros(L,dtype=dtype)
    else:
        pl[fnam]=-999*np.ones(L,dtype=dtype)

# Initialize index to inventory
IdxToInv=[None]*sxy['x'].size
    
# Populate planting layer
cnt_inventory=0
for i in range(ind_at.size):
    
    # Matches are not found if the format is not the same so you need to convert the pl dictionary value to float32
    # to find matches
    ind_pl=np.where(np.float32(d_pl['ACTIVITY_TREATMENT_UNIT_ID'])==atu['ACTIVITY_TREATMENT_UNIT_ID'][ind_at[i]])[0]
    
    for j in range(ind_pl.size):
        
        idx=atu['IdxToSXY'][ind_at[i]]
        
        pl['IdxToSXY'][cnt_inventory]=idx
        
        # Update index to inventory
        if IdxToInv[idx]==None:
            IdxToInv[idx]={}
            IdxToInv[idx]['Index']=np.array([cnt_inventory])
        else:
            IdxToInv[idx]['Index']=np.append(IdxToInv[idx]['Index'],cnt_inventory)
        
        for fnam,flag,dtype in InvLyrInfo[iLyr]['Field List']:
            val=d_pl[fnam][ind_pl[j]]
            if val==None:
                continue
            if flag==0:
                # Numeric variable, leave as is
                pl[fnam][cnt_inventory]=val
            elif flag==1:
                # Convert string variable to numeric variable 
                # based on LUT
                pl[fnam][cnt_inventory]=InvLyrInfo[iLyr]['LUT'][fnam][val]
        # Update counter
        cnt_inventory=cnt_inventory+1
    
# Truncate variables at cnt_inventory
pl['IdxToSXY']=pl['IdxToSXY'][0:cnt_inventory]        
for fnam,flag,dtype in InvLyrInfo[iLyr]['Field List']:
    pl[fnam]=pl[fnam][0:cnt_inventory]
    
# Convert date string to numeric
pl=invu.ExtractDateStringsFromRESULTS(lyr_nam,pl)
       
# Save    
gu.opickle(meta['Paths']['Geospatial'] + '\\RSLT_PLANTING_SVW.pkl',pl)
gu.opickle(meta['Paths']['Geospatial'] + '\\RSLT_PLANTING_SVW_IdxToInv.pkl',IdxToInv)

t1=time.time()
print(t1-t0)


#%% Add FC Archive to FC Inventory dictionary

# Import inputs
#meta=gu.ipickle(r'D:\Data\FCI_Projects\SummaryReforestation\Inputs\Metadata.pkl')
#meta['Paths']['Model Code']=r'C:\Users\rhember\Documents\Code_Python\fcgadgets\cbrunner'
fcinv=gu.ipickle(meta['Paths']['Project'] + '\\Geospatial\\RSLT_FOREST_COVER_INV_SVW.pkl')
fcsilv=gu.ipickle(meta['Paths']['Project'] + '\\Geospatial\\RSLT_FOREST_COVER_SILV_SVW.pkl')

# Save old versions as alterntive names in case you need to redo this/troubleshoot a problem
flg=1
if flg==1:
    gu.opickle(meta['Paths']['Project'] + '\\Geospatial\\RSLT_FOREST_COVER_INV_SVW_before_adding_archive.pkl',fcinv)
    gu.opickle(meta['Paths']['Project'] + '\\Geospatial\\RSLT_FOREST_COVER_SILV_SVW_before_adding_archive.pkl',fcsilv)

# Load initial layers
flg=0
if flg==1:
    fcinv=gu.ipickle(meta['Paths']['Project'] + '\\Geospatial\\RSLT_FOREST_COVER_INV_SVW_before_adding_archive.pkl')
    fcsilv=gu.ipickle(meta['Paths']['Project'] + '\\Geospatial\\RSLT_FOREST_COVER_SILV_SVW_before_adding_archive.pkl')

# Run script
fcinv,fcsilv=invu.ForestCover_AddArchive(meta,fcinv,fcsilv)

# Save
gu.opickle(meta['Paths']['Project'] + '\\Geospatial\\RSLT_FOREST_COVER_INV_SVW.pkl',fcinv)
gu.opickle(meta['Paths']['Project'] + '\\Geospatial\\RSLT_FOREST_COVER_SILV_SVW.pkl',fcsilv)

# Add forest health to forest cover
fcinv=ForestCover_AddForestHealth(meta,fcinv)

# Save
gu.opickle(meta['Paths']['Project'] + '\\Geospatial\\RSLT_FOREST_COVER_INV_SVW.pkl',fcinv)


#%% Extract forest cover polygons for the AT sample
# Used by interactive map

def Get_FC_Polygons_ForMap():

    atu_multipolygons=gu.ipickle(meta['Paths']['Project'] + '\\Geospatial\\atu_multipolygons.pkl')

    # Unique list of opening ID
    opid_atu=np.zeros(len(atu_multipolygons))
    for i in range(len(atu_multipolygons)):
        opid_atu[i]=atu_multipolygons[i]['OPENING_ID']
    opid_atu=np.unique(opid_atu)

    fin=r'C:\Users\rhember\Documents\Data\ForestInventory\Results\20210401\Results.gdb'
    a=[None]*int(1e5)
    cnt=0
    with fiona.open(fin,layer='RSLT_FOREST_COVER_INV_SVW') as source:
        for feat in source:
            if feat['geometry']==None:
                continue
            if np.isin(feat['properties']['OPENING_ID'],opid_atu)==True:
                a[cnt]=feat
                cnt=cnt+1
                print(cnt)
    a=a[0:cnt-1]

    # Give it the BC spatial reference system
    gdf_bm=gpd.read_file(r'C:\Users\rhember\Documents\Data\Basemaps\bc_land.shp')
    gdf_fcinv=gpd.GeoDataFrame.from_features(a,crs=gdf_bm.crs)

    #gdf_fcinv2=gdf_fcinv.to_crs({'init':'epsg:4326'})

    # Save
    gdf_fcinv.to_file(filename=meta['Paths']['Project'] + '\\Geospatial\\fcinv.geojson',driver='GeoJSON')

    return

#%% Get harvest - FCI project intersections for mapping

def Get_HarvestProject_Intersection():

    # Import project polygons
    gdf=gpd.read_file(r'C:\Users\rhember\Documents\Data\FCI_Projects\FCI_RollupFCI_Inv\Geospatial\atu_polygons.geojson')

    # Import harvest after 2017
    gdf_cb=gpd.read_file(r'C:\Users\rhember\Documents\Data\ForestInventory\Disturbances\20210401\Disturbances.gdb',layer='VEG_CONSOLIDATED_CUT_BLOCKS_SP')
    gdf_cb1=gdf_cb[gdf_cb['HARVEST_YEAR']>=2018]

    # Overlay
    gdf_int=gpd.overlay(gdf,gdf_cb1,how='intersection')

    # Save
    gdf_int.to_file(filename=meta['Paths']['Project'] + '\\Geospatial\\harvest_intersection.geojson',driver='GeoJSON')

    return



