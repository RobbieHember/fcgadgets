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

Paths={}
Paths['Project']=r'C:\Users\rhember\Documents\Data\FCI_Projects\FCI_SparseGrid'
Paths['Results']=r'C:\Users\rhember\Documents\Data\ForestInventory\Results\20200430'
Paths['VRI']=r'C:\Users\rhember\Documents\Data\ForestInventory\VRI\20200430'
Paths['Disturbances']=r'C:\Users\rhember\Documents\Data\ForestInventory\Disturbances\20200430'
Paths['Geospatial']=Paths['Project'] + '\\Geospatial'
Paths['Mapping']=Paths['Project'] + '\\Outputs\Mapping'
Paths['Figures']=r'G:\My Drive\Figures\SparseGrid'
Paths['Results File']=Paths['Results'] + '\\Results.gdb'
Paths['VRI File']=Paths['VRI'] + '\\VRI.gdb'
Paths['Disturbances File']=Paths['Disturbances'] + '\\Disturbances.gdb'
Paths['Wildfire Stats and Scenarios File']=r'G:\My Drive\Data\Wildfire\Wildfire_Stats_Scenarios_By_BGCZ\Wildfire_Stats_Scenarios_By_BGCZ.pkl'
Paths['IBM Stats and Scenarios File']=r'G:\My Drive\Data\G:\My Drive\Data\Beetle_Stats_Scenarios_By_BGCZ\IBM_Stats_Scenarios_By_BGCZ.pkl'

# Save
gu.opickle(Paths['Project'] + '\\Inputs\Paths.pkl',Paths)

#%% Define sparse grid sample

# Import TSA maps
zTSA=gis.OpenGeoTiff(r'Z:\!Workgrp\Forest Carbon\Data\BC1ha\Admin\tsa.tif')
lut_tsa=pd.read_excel(r'Z:\!Workgrp\Forest Carbon\Data\BC1ha\Admin\lut_tsa.xlsx')
tsa_boundaries=gpd.read_file(r'C:\Users\rhember\Documents\Data\TSA\tsa_boundaries.shp')

# Import land cover
zLC2=gis.OpenGeoTiff(r'Z:\!Workgrp\Forest Carbon\Data\BC1ha\VRI\lc2.tif')

# Define regular grid sampling frequency
sfreq=100

# Extract subgrid
zTSA.X=zTSA.X[0::sfreq,0::sfreq]
zTSA.Y=zTSA.Y[0::sfreq,0::sfreq]
zTSA.Data=zTSA.Data[0::sfreq,0::sfreq]
zLC2.Data=zLC2.Data[0::sfreq,0::sfreq]

# Define additional inclusion criteria

# Treed, province-wide
iIreg=np.where( (zLC2.Data==4) )

# Treed, Williams Lake TSA only
#iTSA=lut_tsa.loc[lut_tsa.Name=='Williams Lake TSA','VALUE'].values
#ind=np.where( (zLC2.Data==4) & (zTSA.Data==iTSA) )

# Apply filters to BC1ha grid 
sxy={}
sxy['x']=zTSA.X[iIreg]
sxy['y']=zTSA.Y[iIreg]
sxy['ID_TSA']=zTSA.Data[iIreg]

# Save to pickle file
gu.opickle(Paths['Geospatial'] + '\\sxy.pkl',sxy)

# Save as shapefile
flg=1
if flg==1:
    points=[]
    for k in range(sxy['x'].size):
        points.append(Point(sxy['x'][k],sxy['y'][k]))
    gdf_sxy=gpd.GeoDataFrame({'geometry':points,'ID_TSA':sxy['ID_TSA']})
    gdf_sxy.crs=tsa_boundaries.crs  
    gdf_sxy.to_file(Paths['Geospatial'] + '\\sxy.shp')

#%% Plot

# Load basemap
gdf_bm=gpd.read_file(r'Z:\!Workgrp\Forest Carbon\Data\Basemaps\bcbound.shp')

plt.close('all')
fig,ax=plt.subplots(figsize=gu.cm2inch(7.8,6.6))
#mngr=plt.get_current_fig_manager() 
#mngr.window.setGeometry(700,20,620,600)
gdf_bm.plot(ax=ax,facecolor=[0.8,0.8,0.8],edgecolor=[0,0,0],label='Political Boundary',linewidth=0.25,alpha=1)
tsa_boundaries.plot(ax=ax,facecolor='none',edgecolor=[0,0,0],linewidth=0.25)
gdf_sxy.plot(ax=ax,markersize=1,facecolor=[0,0,0.75],edgecolor=None,linewidth=0.75,alpha=1)
ax.grid(color='k',linestyle='-',linewidth=0.25)

iP=2000
for iD in range(len(nddat[iP])):
    x,y=nddat[iP][iD]['Geometry'].exterior.xy
    plt.plot(x,y,'r-')

ax.set(position=[0.01,0.01,0.98,0.98],xticks=[],yticks=[])
#plt.savefig(PathProject + '\\SparseGrid_Map.png',format='png',dpi=900)

#%% Get N deposition time series for each sparse grid cell

tv=np.arange(1971,2021,1)
ndep=np.zeros((tv.size,len(gdf_sxy)))
for iMP in range(len(nddat)):
    print(iMP)
    it=np.where(tv==nddat[iMP][0]['Year'])[0]
    if it.size==0:
        continue
    FlagDone=np.zeros(len(gdf_sxy))
    for iD in range(len(nddat[iMP])):
        InPoly=gdf_sxy.within(nddat[iMP][iD]['Geometry'])
        ind=np.where( (InPoly==True) & (FlagDone==0) )[0]
        #gdf_sxy.loc[InPoly].plot(ax=ax,markersize=1,facecolor=[1,0,0.25],edgecolor=None,linewidth=0.75,alpha=1)
        if ind.size>0:
            ndep[it,ind]=ndep[it,ind]+nddat[iMP][iD]['N deposition']
            FlagDone[ind]=1

plt.plot(tv,np.prctile(ndep,axis=1),'-k.')

gu.opickle(r'C:\Users\rhember\Documents\Data\FCI_Projects\FertilizationSummaryNdep\Geospatial\ndep.pkl',ndep)


#%% Open crosswalk between missing AT geometries and opening geometries (if already exists)

atu_mis=gu.ipickle(Paths['Results'] + '\\atu_mis.pkl')
at_geo_from_op=gu.ipickle(Paths['Results'] + '\\at_geo_from_op.pkl')

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
    with fiona.open(path,layer=lyr_nam) as source:    
        
        for feat in source:
            
            # Extract attributes and geometry
            prp=feat['properties']
            geom=feat['geometry']
            
            # Fill missing AT spatial with OP spatial            
            if (lyr_nam=='RSLT_ACTIVITY_TREATMENT_SVW'):
                
                # Populate missing ATU layer geometry with geometry from 
                # OPENING layer where possible.
                flg_geom_from_op=0
                if (geom==None):
                    ind=np.where(atu_mis['OPENING_ID']==prp['OPENING_ID'])[0]
                    if ind.size>0:
                        geom=at_geo_from_op[ind[0]]
                        flg_geom_from_op=1
                    else:
                        print('Checked for opening spatial, but no match')
            
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
    
    gu.opickle(Paths['Geospatial'] + '\\' + InvLyrInfo[iLyr]['Layer Name'] + '.pkl',data)
    gu.opickle(Paths['Geospatial'] + '\\' + InvLyrInfo[iLyr]['Layer Name'] + '_IdxToInv.pkl',IdxToInv)

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
#pl=gu.ipickle(Paths['TileGeospatial'] + '\\RSLT_PLANTING_SVW.pkl')
    
# Get keys for planting layer
key_pl=[]
for fnam,flag,dtype in InvLyrInfo[iLyr]['Field List']:
    key_pl.append(fnam)
            
# Open AT sparse grid
atu=gu.ipickle(Paths['Project'] + '\\Inputs\\Geospatial\\' + '\\RSLT_ACTIVITY_TREATMENT_SVW.pkl')

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
gu.opickle(Paths['Geospatial'] + '\\RSLT_PLANTING_SVW.pkl',pl)
gu.opickle(Paths['Geospatial'] + '\\RSLT_PLANTING_SVW_IdxToInv.pkl',IdxToInv)
