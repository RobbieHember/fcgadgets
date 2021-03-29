'''

FOREST CARBON BALANCE BY SPARSE GRID SAMPLING

'''

#%% Import modules

import os
import numpy as np
import pandas as pd
import geopandas as gpd
import fiona
import matplotlib.pyplot as plt
import numpy.matlib as ml
import matplotlib.pyplot as plt
import time
from shapely.geometry import Polygon,Point
from fcgadgets.utilities import utilities_general as gu
from fcgadgets.utilities import utilities_gis as gis
from fcgadgets.utilities import utilities_inventory as invu
from fcgadgets.cbrunner import cbrun_utilities

#%% Define paths

meta={}
meta['Paths']={}
#meta['Paths']['Project']=r'C:\Users\rhember\Documents\Data\FCI_Projects\FCI_SparseGrid'
meta['Paths']['Project']=r'D:\Data\FCI_Projects\SparseGrid_HighRes'
meta['Paths']['Results']=r'C:\Users\rhember\Documents\Data\ForestInventory\Results\20210208'
meta['Paths']['VRI']=r'C:\Users\rhember\Documents\Data\ForestInventory\VRI\20200430'
meta['Paths']['Disturbances']=r'C:\Users\rhember\Documents\Data\ForestInventory\Disturbances\20200430'
meta['Paths']['LandUse']=r'C:\Users\rhember\Documents\Data\ForestInventory\LandUse\20200706'
meta['Paths']['Geospatial']=meta['Paths']['Project'] + '\\Geospatial'
meta['Paths']['Taz Datasets']=r'C:\Users\rhember\Documents\Data\Taz Datasets'

# Save
gu.opickle(meta['Paths']['Project'] + '\\Inputs\\Metadata.pkl',meta)

#%% Define sparse grid sample

# Import TSA maps
zTSA=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Admin\tsa.tif')
lut_tsa=pd.read_excel(r'C:\Users\rhember\Documents\Data\BC1ha\Admin\lut_tsa.xlsx')
tsa_boundaries=gpd.read_file(r'C:\Users\rhember\Documents\Data\TSA\tsa_boundaries.shp')

# Import land cover
zLC2=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\VRI\lc2.tif')

# Define regular grid sampling frequency
sfreq=20

# Extract subgrid
zTSA['X']=zTSA['X'][0::sfreq,0::sfreq]
zTSA['Y']=zTSA['Y'][0::sfreq,0::sfreq]
zTSA['m'],zTSA['n']=zTSA['X'].shape
zTSA['Data']=zTSA['Data'][0::sfreq,0::sfreq]
zLC2['Data']=zLC2['Data'][0::sfreq,0::sfreq]

# Define additional inclusion criteria

# Treed, province-wide
iIreg=np.where( (zLC2['Data']==4) )

# Treed, Williams Lake TSA only
#iTSA=lut_tsa.loc[lut_tsa.Name=='Williams Lake TSA','VALUE'].values
#ind=np.where( (zLC2.Data==4) & (zTSA.Data==iTSA) )

# Save grid
flg=0
if flg==1:
    z=zTSA.copy()
    z['Data']=np.zeros(zTSA['Data'].shape,dtype='int16')
    z['Data'][iIreg]=1 # Treed = 1
    gis.SaveGeoTiff(z,meta['Paths']['Project'] + '\\Geospatial\\GridSXY.tiff')
    plt.matshow(z['Data'])

#%% Generate sparse grid

# Apply filters to BC1ha grid 
sxy={}
sxy['x']=zTSA['X'][iIreg]
sxy['y']=zTSA['Y'][iIreg]
sxy['ID_TSA']=zTSA['Data'][iIreg]

# Save to pickle file
gu.opickle(meta['Paths']['Geospatial'] + '\\sxy.pkl',sxy)

# Save as shapefile
flg=1
if flg==1:
    points=[]
    for k in range(sxy['x'].size):
        points.append(Point(sxy['x'][k],sxy['y'][k]))
    gdf_sxy=gpd.GeoDataFrame({'geometry':points,'ID_TSA':sxy['ID_TSA']})
    gdf_sxy.crs=tsa_boundaries.crs  
    gdf_sxy.to_file(meta['Paths']['Geospatial'] + '\\sxy.shp')

#%% Plot

# Load basemap
gdf_bm=gpd.read_file(r'C:\Users\rhember\Documents\Data\Basemaps\Basemaps.gdb',layer='NRC_POLITICAL_BOUNDARIES_1M_SP')

plt.close('all')
fig,ax=plt.subplots(figsize=gu.cm2inch(7.8,6.6))
#mngr=plt.get_current_fig_manager() 
#mngr.window.setGeometry(700,20,620,600)
gdf_bm.plot(ax=ax,facecolor=[0.8,0.8,0.8],edgecolor=[0,0,0],label='Political Boundary',linewidth=0.25,alpha=1)
tsa_boundaries.plot(ax=ax,facecolor='none',edgecolor=[0,0,0],linewidth=0.25)
gdf_sxy.plot(ax=ax,markersize=1,facecolor=[0,0,0.75],edgecolor=None,linewidth=0.75,alpha=1)
ax.grid(color='k',linestyle='-',linewidth=0.25)

#iP=2000
#for iD in range(len(nddat[iP])):
#    x,y=nddat[iP][iD]['Geometry'].exterior.xy
#    plt.plot(x,y,'r-')

ax.set(position=[0.01,0.01,0.98,0.98],xticks=[],yticks=[])
#plt.savefig(PathProject + '\\SparseGrid_Map.png',format='png',dpi=900)

##%% Get N deposition time series for each sparse grid cell
#
#tv=np.arange(1971,2021,1)
#ndep=np.zeros((tv.size,len(gdf_sxy)))
#for iMP in range(len(nddat)):
#    print(iMP)
#    it=np.where(tv==nddat[iMP][0]['Year'])[0]
#    if it.size==0:
#        continue
#    FlagDone=np.zeros(len(gdf_sxy))
#    for iD in range(len(nddat[iMP])):
#        InPoly=gdf_sxy.within(nddat[iMP][iD]['Geometry'])
#        ind=np.where( (InPoly==True) & (FlagDone==0) )[0]
#        #gdf_sxy.loc[InPoly].plot(ax=ax,markersize=1,facecolor=[1,0,0.25],edgecolor=None,linewidth=0.75,alpha=1)
#        if ind.size>0:
#            ndep[it,ind]=ndep[it,ind]+nddat[iMP][iD]['N deposition']
#            FlagDone[ind]=1
#
#plt.plot(tv,np.prctile(ndep,axis=1),'-k.')
#
#gu.opickle(r'C:\Users\rhember\Documents\Data\FCI_Projects\FertilizationSummaryNdep\Geospatial\ndep.pkl',ndep)


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
    IdxToInv=[None]*sxy['x'].size
        
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
                
                            # Could not use either FC or openign layer
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
IdxToInv=[None]*sxy['x'].size
    
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
