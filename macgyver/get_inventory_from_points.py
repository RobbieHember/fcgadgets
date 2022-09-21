'''

INVENTORY FOR SPARSE REGULAR GRID SAMPLE

'''

#%% Import modules

import os
import numpy as np
import pandas as pd
import geopandas as gpd
import fiona
import gc as garc
import matplotlib.pyplot as plt
import numpy.matlib as ml
import time
from shapely.geometry import Polygon,Point
from fcgadgets.macgyver import utilities_general as gu
from fcgadgets.macgyver import utilities_gis as gis
from fcgadgets.macgyver import utilities_inventory as invu
from fcgadgets.cbrunner import cbrun_utilities

#%% Project name

#name='SummaryBC5k'
#name='SummaryBC20k_H'
#name='SummaryBC20k_HAR'
#name='SummaryBC20k_NA'
#name='LICS Hanceville'
#name='SummaryQuesnel'
#name='SummaryBC_NOSE'
name='SummaryBC_OHS'

#%% Define paths

meta={}
meta['Paths']={}
meta['Paths']['Project']=r'D:\Data\FCI_Projects' + '\\' + name 
meta['Paths']['Geospatial']=meta['Paths']['Project'] + '\\Geospatial'
meta['Paths']['Results']=r'C:\Users\rhember\Documents\Data\ForestInventory\Results\20220422'
meta['Paths']['VRI']=r'C:\Users\rhember\Documents\Data\ForestInventory\VRI\20220404'
meta['Paths']['Disturbances']=r'C:\Users\rhember\Documents\Data\ForestInventory\Disturbances\20220422'
meta['Paths']['LandUse']=r'C:\Users\rhember\Documents\Data\ForestInventory\LandUse\20220422'
meta['Paths']['Model Code']=r'C:\Users\rhember\Documents\Code_Python\fcgadgets\cbrunner'
meta['Paths']['Taz Datasets']=r'C:\Users\rhember\Documents\Data\Taz Datasets'

# Save
gu.opickle(meta['Paths']['Project'] + '\\Inputs\\Metadata.pkl',meta)

#%% Import maps

# Import raster grids
zTSA=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Admin\tsa.tif')
lut_tsa=pd.read_excel(r'C:\Users\rhember\Documents\Data\BC1ha\Admin\lut_tsa.xlsx')
tsa_boundaries=gpd.read_file(r'C:\Users\rhember\Documents\Data\TSA\tsa_boundaries.shp')
zLC2=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\VRI\lc2.tif')
bc_grid_size=zTSA['Data'].size

#%% Define geospatial variables

# Initialize geospatial info structure
geos={}

# Define regular grid sampling frequency

#geos['rgsf']=1 # 100 m
#geos['rgsf']=5 # 500 m
#geos['rgsf']=10 # 1 km
#geos['rgsf']=20 # 2 km
#geos['rgsf']=40 # 4 km
#geos['rgsf']=50 # 5 km
#geos['rgsf']=100 # 10 km
geos['rgsf']=200 # 20 km

# Extract subgrid
zTSA['Data']=zTSA['Data'][0::geos['rgsf'],0::geos['rgsf']]
zLC2['Data']=zLC2['Data'][0::geos['rgsf'],0::geos['rgsf']]

# Populate geospatial structure with grid attributes
geos['X']=zTSA['X'][0::geos['rgsf'],0::geos['rgsf']]
geos['Y']=zTSA['Y'][0::geos['rgsf'],0::geos['rgsf']]
geos['m'],geos['n']=geos['X'].shape
geos['Mask']=np.zeros((geos['m'],geos['n']),dtype=np.int8)

# Area expansion factor
geos['AEF']=bc_grid_size/(geos['m']*geos['n'])

# Define additional sampling criteria

if (name=='SummaryBC5k') | (name=='SummaryBC_OHS'):
    
    # Treed, global
    geos['iMask']=np.where( (zLC2['Data']==4) )
    
elif (name=='SummbaryBC_Quesnel'):
    
    # Treed, Williams Lake TSA only
    iTSA=lut_tsa.loc[lut_tsa.Name=='Quesnel TSA','VALUE'].values
    geos['iMask']=np.where( (zLC2['Data']==4) & (zTSA['Data']==iTSA) )

elif (name=='BCTS'):
    
    zBCTS=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\LandUseLandCover\bcts_op_area.tif')
    zBCTS['Data']=zBCTS['Data'][0::geos['rgsf'],0::geos['rgsf']]
    geos['iMask']=np.where( (zLC2['Data']==4) & (zBCTS['Data']>0) )

elif (name=='LICS Hanceville'):
    
    lics=gu.ipickle(r'D:\Data\FCI_Projects\LICS Site List\LICS Sites.pkl')
    id=5
    width=30000
    zMask=np.zeros(geos['Mask'].shape,dtype=int)
    ind=np.where( (np.abs(geos['X']-lics[id]['X'])<=width/2) & (np.abs(geos['Y']-lics[id]['Y'])<=width/2) )
    zMask[ind]=1
    #geos['iMask']=np.where( (zLC2['Data']==4) & (zMask==1) )
    geos['iMask']=np.where( (zMask==1) )
    
elif (name=='Waste'):
    
    zW=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\Waste Wood\openings.tif')
    zW['Data']=zW['Data'][0::geos['rgsf'],0::geos['rgsf']]
    geos['iMask']=np.where( (zW['Data']>0) )
  
elif (name=='20k_SS'):
    
    zBGC=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\VRI\becz.tif')
    zBGC['Data']=zBGC['Data'][0::geos['rgsf'],0::geos['rgsf']]
    zMask=np.zeros(zBGC['Data'].shape)
    u=np.unique(zBGC['Data'])
    u=u[u!=255]
    th=16
    for iU in range(u.size):
        ind=np.where(zBGC['Data']==u[iU])
        rny=np.argsort(np.random.random(ind[0].size))
        rnx=np.argsort(np.random.random(ind[1].size))
        ind_y=ind[0][rny[0:th]].astype(int)
        ind_x=ind[1][rnx[0:th]].astype(int)
        zMask[ind_y,ind_x]=1
    
    geos['iMask']=np.where( (zLC2['Data']==4) & (zMask==1) )

elif (name=='SummaryBC_NOSE'):
    
    x=geos['X'].flatten()
    y=geos['Y'].flatten()
    points=[]
    for k in range(x.size):
        points.append(Point(x[k],y[k]))
    gdf_xy=gpd.GeoDataFrame({'geometry':points,'ID_TSA':1})
    gdf_xy.crs=tsa_boundaries.crs 
    
    # Import polygons of non-obligatoin stand establishment (from query script)
    gdf_nose=gpd.read_file(r'D:\Data\FCI_Projects\SummaryBC_NOSE\Geospatial\NOSE_query.geojson')

    # Get points within polygons
    gdf_sxy=gpd.sjoin(gdf_xy,gdf_nose,op='within')
    
    geos['Sparse']={}
    geos['Sparse']['X']=gdf_sxy['geometry'].x.values
    geos['Sparse']['Y']=gdf_sxy['geometry'].y.values

    for i in range(geos['Sparse']['X'].size):
        ind=np.where( (geos['X']==geos['Sparse']['X'][i]) & (geos['Y']==geos['Sparse']['Y'][i]) )
        geos['Mask'][ind]=1
        
    geos['iMask']=np.where( (geos['Mask']==1) )

# Revise mask
geos['Mask'][geos['iMask']]=1

# Generate sparse grid
geos['Sparse']={}
geos['Sparse']['X']=geos['X'][geos['iMask']]
geos['Sparse']['Y']=geos['Y'][geos['iMask']]
geos['Sparse']['ID_TSA']=zTSA['Data'][geos['iMask']]

# Save to pickle file
geos['X']=geos['X'][0,:]
geos['Y']=geos['Y'][:,0]
gu.opickle(meta['Paths']['Geospatial'] + '\\geos.pkl',geos)

# Save mask as geotiff
#gis.SaveGeoTiff(z,meta['Paths']['Project'] + '\\Geospatial\\geos_grid.tiff')
#plt.matshow(geos['Mask'])

# Save sparse points to geojson
flg=1
if flg==1:
    points=[]
    for k in range(geos['Sparse']['X'].size):
        points.append(Point(geos['Sparse']['X'][k],geos['Sparse']['Y'][k]))
    gdf_sxy=gpd.GeoDataFrame({'geometry':points,'ID_TSA':geos['Sparse']['ID_TSA']})
    gdf_sxy.crs=tsa_boundaries.crs  
    gdf_sxy.to_file(meta['Paths']['Geospatial'] + '\\geos.geojson',driver='GeoJSON')

#%% Plot

# Load basemap
gdf_bm=gpd.read_file(r'C:\Users\rhember\Documents\Data\Basemaps\Basemaps.gdb',layer='NRC_POLITICAL_BOUNDARIES_1M_SP')

#gdf_sxy=gpd.read_file(meta['Paths']['Geospatial'] + '\\geos.geojson')

plt.close('all')
fig,ax=plt.subplots(figsize=gu.cm2inch(7.8,6.6)); ms=3
#mngr=plt.get_current_fig_manager() 
#mngr.window.setGeometry(700,20,620,600)
gdf_bm.plot(ax=ax,facecolor=[0.8,0.8,0.8],edgecolor=[0,0,0],label='Political Boundary',linewidth=0.25,alpha=1)
tsa_boundaries.plot(ax=ax,facecolor='none',edgecolor=[0,0,0],linewidth=0.25)
gdf_sxy.plot(ax=ax,markersize=ms,facecolor=[0.75,0,0],edgecolor=None,linewidth=0.75,alpha=1)
ax.grid(color='k',linestyle='-',linewidth=0.25)
ax.set(position=[0.01,0.01,0.98,0.98],xticks=[],yticks=[])
#plt.savefig(PathProject + '\\SparseGrid_Map.png',format='png',dpi=900)
plt.close('all')

garc.collect()

#%% Open crosswalk between missing AT geometries and opening geometries
# If this doesn't work, you need to run the script that creates the crosswalk

missing_geo_atu_list=gu.ipickle(meta['Paths']['Results'] + '\\missing_geo_atu_list.pkl')
missing_geo_op_geos=gu.ipickle(meta['Paths']['Results'] + '\\missing_geo_op_geos.pkl')
missing_geo_fc_geos=gu.ipickle(meta['Paths']['Results'] + '\\missing_geo_fc_geos.pkl')

#%% Import inventory layer information (names, variables, LUTs)

InvLyrInfo=invu.DefineInventoryLayersAndVariables()

for iLyr in range(len(InvLyrInfo)):
    lut=gu.ipickle(InvLyrInfo[iLyr]['Path'] + '\\LUTs_' + InvLyrInfo[iLyr]['Layer Name'] +'.pkl')
    for key in lut.keys():
        InvLyrInfo[iLyr]['LUT'][key]=lut[key]
    del lut

#%% Compile inventory layers

for iLyr in range(len(InvLyrInfo)):

    # Loop through features in layer
    t_start=time.time()
    
    # Define path
    path=InvLyrInfo[iLyr]['Path'] + '\\' + InvLyrInfo[iLyr]['File Name']
    
    # Define layer
    lyr_nam=InvLyrInfo[iLyr]['Layer Name']
    print(lyr_nam)
        
    # Don't run this for planting - planting is done seperately below
    if lyr_nam=='RSLT_PLANTING_SVW':
        continue
        
    # Initialize index to inventory
    IdxToInv=[None]*geos['Sparse']['X'].size
        
    # Initialize inventory dictionary
    L=5000000
    data={}
    data['IdxToSXY']=np.zeros(L,dtype=int)
    for fnam,flag,dtype in InvLyrInfo[iLyr]['Field List']:
        if dtype=='<U20':
            data[fnam]=np.zeros(L,dtype=dtype)
        else:
            data[fnam]=-999*np.ones(L,dtype=dtype)
    cnt_inventory=0
        
    # Scan through layer file to extract selected variables, and convert string
    # variables to numeric based on LUTs
    cc=0
    with fiona.open(path,layer=lyr_nam) as source:    
        
        for feat in source:
            
            # Extract attributes and geometry
            prp=feat['properties']
            geom=feat['geometry']
            
            # Fill missing AT spatial with OP spatial            
            if (lyr_nam=='RSLT_ACTIVITY_TREATMENT_SVW'):
                
                # No need to record surveys
                #if prp['SILV_BASE_CODE']=='SU':
                #    continue
                
                # Populate missing ATU layer geometry with geometry from 
                # OPENING or FC layer where possible.
                flg_geom_from_op=0
                flg_geom_from_fc=0
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
                    
                            # Plot (not working)
                            #flg=0
                            #if flg==1:                        
                            #    plt.close('all')
                            #    fig,ax=plt.subplots(1)
                            #    gdf_fc=gpd.GeoDataFrame.from_features(feat_fc)
                            #    gdf_fc.plot(ax=ax,facecolor='None',edgecolor='r',linewidth=1.25,linestyle='--')          
                
                        if prp['OPENING_ID'] in missing_geo_op_geos:
                            
                            if len(missing_geo_op_geos[prp['OPENING_ID']])>0:
                    
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
                
                            # Could not use either FC or opening layer
                            print('Missing spatial could not be recovered')
                        
            # Only continue if spatial info exists
            if (geom==None) | (geom==[]):
                continue
            
            # Extract multipolygon
            coords0=geom['coordinates']
        
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
                    
                    # This should speed it up a bit
                    if np.max(x_feat)<np.min(geos['Sparse']['X']):
                        continue
                    if np.min(x_feat)>np.max(geos['Sparse']['X']):
                        continue
                    if np.max(y_feat)<np.min(geos['Sparse']['Y']): 
                        continue
                    if np.min(y_feat)>np.max(geos['Sparse']['Y']): 
                        continue
        
                    # Isolate cells within bounding box of the feature polygon
                    iBB=np.where( (geos['Sparse']['X']>=np.min(x_feat)-1000) & (geos['Sparse']['X']<=np.max(x_feat)+1000) & (geos['Sparse']['Y']>=np.min(y_feat)-1000) & (geos['Sparse']['Y']<=np.max(y_feat)+1000) )[0]
                    
                    # Only continue if overlaop
                    if iBB.size==0:
                        continue
        
                    # Isolate cells within the polygon
                    InPoly=gis.InPolygon(geos['Sparse']['X'][iBB],geos['Sparse']['Y'][iBB],x_feat,y_feat)
                    
                    iInPoly=np.where(InPoly==1)[0]
        
                    # Index to cells within polygon
                    iKeep=iBB[iInPoly]
                    
                    # Only continue if overlap
                    if iKeep.size==0:     
                        continue
                                        
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
    
    gu.opickle(meta['Paths']['Geospatial'] + '\\' + InvLyrInfo[iLyr]['Layer Name'] + '.pkl',data)
    gu.opickle(meta['Paths']['Geospatial'] + '\\' + InvLyrInfo[iLyr]['Layer Name'] + '_IdxToInv.pkl',IdxToInv)

#%% Retrieve planting information 
# Some projects did not report spatial planting info. Without the spatial info
# in the planting layer, the initial import of the PL layer (above) will miss
# planting info in some AT polygons. Use the ACTIVITY TREATMENT UNIT ID as
# a crosswalk to retrieve all planting info for each activity.
# (10 min)
    
t0=time.time()
    
# Get planting lyaer id
for iLyr in range(len(InvLyrInfo)):
    if InvLyrInfo[iLyr]['Layer Name']=='RSLT_PLANTING_SVW':
        break
    
# Import geodataframe, drop geometry and convert to dict
path=InvLyrInfo[iLyr]['Path'] + '\\' + InvLyrInfo[iLyr]['File Name']
lyr_nam=InvLyrInfo[iLyr]['Layer Name']
gdf_pl=gpd.read_file(path,layer=lyr_nam)        
d_pl=gu.DataFrameToDict(gdf_pl.drop(columns='geometry'))
        
# Open the planting sparse grid
#pl=gu.ipickle(meta['Paths']['TileGeospatial'] + '\\RSLT_PLANTING_SVW.pkl')
    
# Get keys for planting layer
key_pl=[]
for fnam,flag,dtype in InvLyrInfo[iLyr]['Field List']:
    key_pl.append(fnam)
            
# Open AT sparse grid
atu=gu.ipickle(meta['Paths']['Project'] + '\\Geospatial\\RSLT_ACTIVITY_TREATMENT_SVW.pkl')

pl_code=InvLyrInfo[0]['LUT']['SILV_BASE_CODE']['PL']
    
# Only proceed if planting occurs
ind_at=np.where(atu['SILV_BASE_CODE']==pl_code)[0]

# Initialize dictionary
L=1000000
pl={}
pl['IdxToSXY']=np.zeros(L,dtype=int)
for fnam,flag,dtype in InvLyrInfo[iLyr]['Field List']:
    if dtype=='<U20':
        pl[fnam]=np.zeros(L,dtype=dtype)
    else:
        pl[fnam]=-999*np.ones(L,dtype=dtype)

# Initialize index to inventory
IdxToInv=[None]*geos['Sparse']['X'].size
    
# Populate planting layer
cnt_inventory=0
for i in range(ind_at.size):
    ind_pl=np.where(d_pl['ACTIVITY_TREATMENT_UNIT_ID']==atu['ACTIVITY_TREATMENT_UNIT_ID'][ind_at[i]])[0]
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
        
t1=time.time()
print(t1-t0)
       
# Save    
gu.opickle(meta['Paths']['Geospatial'] + '\\RSLT_PLANTING_SVW.pkl',pl)
gu.opickle(meta['Paths']['Geospatial'] + '\\RSLT_PLANTING_SVW_IdxToInv.pkl',IdxToInv)

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

# Save revised versions
gu.opickle(meta['Paths']['Project'] + '\\Geospatial\\RSLT_FOREST_COVER_INV_SVW.pkl',fcinv)
gu.opickle(meta['Paths']['Project'] + '\\Geospatial\\RSLT_FOREST_COVER_SILV_SVW.pkl',fcsilv)