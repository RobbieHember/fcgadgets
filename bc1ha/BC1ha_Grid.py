
#%% IMPORT MODULES

import sys
import numpy as np
from osgeo import gdal
from osgeo import osr
import matplotlib.pyplot as plt
import geopandas as gpd
from rasterio.transform import from_origin
import pandas as pd
import fiona
import rasterio
from rasterio import features
from shapely.geometry import Point, Polygon
from shapely import geometry
import cv2

import fcgadgets.macgyver.utilities_general as gu
import fcgadgets.macgyver.utilities_gis as gis
from fcgadgets.macgyver import utilities_inventory as invu
import fcgadgets.macgyver.utilities_query_gdb as qgdb

#%% Create a grid that is consistent with what FAIB uses

flg=0
if flg==1:
    # Load BC basemap to get projection information
    bm=gpd.read_file(r'C:\Users\rhember\Documents\Data\Basemaps\bc_land.shp')

    # Inputs (west,north,xsize,ysize) - received from Tyler Muhley
    transform=from_origin(159587.5,1748187.5,100,100)

    fin_Ref=r'C:\Users\rhember\Documents\Data\BC1ha\Admin\FAIB_Standard.tif'
    nd=rasterio.open(fin_Ref,'w',driver='GTiff',height=15744,width=17216,count=1,dtype=rasterio.float32,compress='lzw',crs=bm.crs,transform=transform)
    nd.write(nd,1) # This crashes but then it seems to work if you then close it.
    nd.close()
    #z=gis.OpenGeoTiff(pth)
else:
    fref=r'C:\Users\rhember\Documents\Data\BC1ha\Admin\FAIB_Standard.tif'
    zRef=gis.OpenGeoTiff(fref)

#%% Prepare timber supply area

# Create TSA Key
# This creates a key between the raster values for each TSA, and the TSA names
# from a shapefile of TSA downloaded from the BC Data Cat.

# Import the shapefile with the TSA names and get the unique list of TSA Numbers
# and each corresponding TSA name.
gdf_tsa=gpd.read_file(r'C:\Users\rhember\Documents\Data\TSA\tsa.shp')
u=gdf_tsa['TSA_NUMBER'].unique()
nam=[]
for i in range(len(u)):
    a=gdf_tsa[gdf_tsa.TSA_NUMBER==u[i]].TSA_NUMB_1.unique()
    nam.append(a[0])

# Key between raster values and TSA name (this fixes the original file created by ESRI the first time)

# Import raster metadata
db=dbfread.DBF(r'C:\Users\rhember\Documents\Data\BC1ha\Admin\tsa.tif.vat.dbf',load=True)

# Convert to dataframe for key
df=pd.DataFrame([])
for key in db.records[0].keys():
    tmp=np.array([])
    for i in range(len(db.records)):
        tmp=np.append(tmp,db.records[i][key])
    df[key]=tmp

# Add TSA names to key
df['Name']=df['TSA_NUMBER']
for i in range(len(u)):
    ind=np.where(np.array(df['TSA_NUMBER'])==np.array(u[i]))[0]
    df.loc[ind,'Name']=nam[i]
# Save to file
df.to_excel(r'C:\Users\rhember\Documents\Data\BC1ha\Admin\tsa_key.xlsx',index=False)

# Digitize the boundary of TSAs

# Open dataframe containing TSA names and values for raster grid of TSAs
zTSA=gis.OpenGeoTiff(r'Z:\!Workgrp\Forest Carbon\Data\BC1ha\Admin\tsa.tif')

# Import look up table
df_tsa=pd.read_excel(r'Z:\!Workgrp\Forest Carbon\Data\BC1ha\Admin\tsa.tif.vat.xlsx')

# Initialize geodatabase
gdf=gpd.GeoDataFrame(data=[],columns=['Value','Name','geometry'])
cnt=0
for i in range(len(df_tsa)):

    # Define the features (objects) that will be digitized
    id=df_tsa.loc[i,'VALUE']

    # Create binary image
    z=np.zeros((m,n,3),dtype=np.uint8)
    z[tsa==id,:]=255
    z=cv2.cvtColor(z,cv2.COLOR_BGR2GRAY) # Convert to grey scale

    # Calculate contour of object
    cont=cv2.findContours(image=z,mode=cv2.RETR_LIST,method=cv2.CHAIN_APPROX_SIMPLE)

    # Unpack silly tuple

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
            x=xG[c]
            y=yG[r]
            pointList.append(geometry.Point(x,y))
        gdf.loc[cnt,'Value']=df_tsa.loc[i,'VALUE'].astype(float)
        gdf.loc[cnt,'Name']=df_tsa.loc[i,'Name']
        gdf.loc[cnt,'geometry']=geometry.Polygon([[p.x,p.y] for p in pointList])
        cnt=cnt+1

gdf.to_file(filename=r'C:\Users\rhember\Documents\Data\TSA\tsa_boundaries.shp')

#%% Clip rasters to standard grid

# DEM
fin=r'C:\Users\rhember\Documents\Data\BC1ha\Terrain\elevation_old.tif'
fout=r'C:\Users\rhember\Documents\Data\BC1ha\Terrain\elevation.tif'
gis.ClipToRaster_ByFile(fin,fout,fref)

# TSA
fin=r'C:\Users\rhember\Documents\Data\BC1ha\Admin\tsa_old.tif'
fout=r'C:\Users\rhember\Documents\Data\BC1ha\Admin\tsa.tif'
gis.ClipToRaster_ByFile(fin,fout,fref)

# BTM
fin=r'C:\Users\rhember\Documents\Data\BC1ha\LandUseLandCover\landuse.btm_old.tif'
fout=r'C:\Users\rhember\Documents\Data\BC1ha\LandUseLandCover\landuse.btm.tif'
gis.ClipToRaster_ByFile(fin,fout,fref)

# LC2
fin=r'C:\Users\rhember\Documents\Data\BC1ha\VRI\lc2_old.tif'
fout=r'C:\Users\rhember\Documents\Data\BC1ha\VRI\lc2.tif'
gis.ClipToRaster_ByFile(fin,fout,fref)

# VRI age
fin=r'C:\Users\rhember\Documents\Data\BC1ha\VRI\proj_age_1.tif'
fout=r'C:\Users\rhember\Documents\Data\BC1ha\VRI\proj_age_1.tif'
gis.ClipToRaster_ByFile(fin,fout,fref)

# BGC zone
fin=r'C:\Users\rhember\Documents\Data\BC1ha\VRI\becsz_old.tif'
fout=r'C:\Users\rhember\Documents\Data\BC1ha\VRI\becsz.tif'
gis.ClipToRaster_ByFile(fin,fout,fref)

# Mean annual temp
fin=r'C:\Users\rhember\Documents\Data\BC1ha\Climate\BC1ha_mat_norm_1971to2000_si_hist_v1_old.tif'
fout=r'C:\Users\rhember\Documents\Data\BC1ha\Climate\BC1ha_mat_norm_1971to2000_si_hist_v1.tif'
gis.ClipToRaster_ByFile(fin,fout,fref)

fin=r'C:\Users\rhember\Documents\Data\BC1ha\Climate\BC1ha_ws_gs_norm_1971to2000_comp_hist_v1.tif'
fout=r'C:\Users\rhember\Documents\Data\BC1ha\Climate\BC1ha_ws_gs_norm_1971to2000_comp_hist_v1.tif'
gis.ClipToRaster_ByFile(fin,fout,fref)

#%% Rasterize protected lands

# Define paths
meta={}
meta['Paths']={}
meta['Paths']['Project']=r''
meta['Paths']['LandUse']=r'C:\Users\rhember\Documents\Data\ForestInventory\LandUse\20220422\LandUse.gdb'
meta['Paths']['Model Code']=r'C:\Users\rhember\Documents\Code_Python\fcgadgets\cbrunner'

meta=invu.Load_LUTs(meta)

pthin=meta['Paths']['LandUse']
pthout=r'C:\Users\rhember\Documents\Data\BC1ha\LandUseLandCover'

fiona.listlayers(pthin)
lnam='TA_PROTECTED_LANDS_SV'

list(meta['LUT']['LU L']['LEGAL_FEAT_OBJECTIVE'].keys())

# Open the shapefile
df=gpd.read_file(pthin,layer=lnam)

# Remove features with no geometry
df=df[df.geometry!=None]

desig=df['PROTECTED_LANDS_DESIGNATION'].unique()

z=np.zeros(zRef['Data'].shape,dtype=float)

for iD in range(desig.size):
    df0=df[df['PROTECTED_LANDS_DESIGNATION']==desig[iD]].copy()
    df0['ID']=(iD+1)*np.ones(len(df0))
    shapes=((geom,value) for geom, value in zip(df0['geometry'],df0['ID']))
    burned=features.rasterize(shapes=shapes,fill=0,out=z,transform=zRef['Transform'])
plt.matshow(z,clim=[0,6])

z1=zRef.copy()
z1['Data']=z.astype(np.int8)
gis.SaveGeoTiff(z1,pthout + '\\PROTECTED_LANDS_DESIGNATION.tif')

#%% Rasterize select openings

# List of openings to rasterize
dOp=gu.ipickle(r'C:\Users\rhember\Documents\Data\Waste Wood\WasteSummary_UniqueOpenings2006On.pkl')

# Import TSA raster
zTSA=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Admin\tsa.tif')

# Input path to RESULTS database (downloaded from BC data catalogue)
pthin=r'C:\Users\rhember\Documents\Data\ForestInventory\Results\20220422\Results.gdb'

lnam='RSLT_OPENING_SVW'

# Open the shapefile
df=gpd.read_file(pthin,layer=lnam)

# Remove features with no geometry
df=df[df.geometry!=None]

# Keep select openings

op=df['OPENING_ID'].to_numpy()

#u,idx,inv=np.unique(ar,return_index=True,return_inverse=True)
c,ia,ib=np.intersect1d(op,dOp['OPENING_ID'],return_indices=True)

df['Flag']=np.zeros(df['OPENING_ID'].size)
df['Flag'].iloc[ia]=1

df0=df.copy()
df0=df0[df0.Flag==1]

#df0=df.copy()
df0['dummy']=np.ones(len(df0))
shapes=((geom,value) for geom, value in zip(df0['geometry'],df0['dummy']))

z0=np.zeros(zTSA['Data'].shape,dtype=float)
burned=features.rasterize(shapes=shapes,fill=0,out=z0,transform=zTSA['Transform'])

zOut=zTSA.copy()
zOut['Data']=np.zeros(zTSA['Data'].shape,dtype='int16')
zOut['Data'][burned>0]=1
gis.SaveGeoTiff(zOut,r'C:\Users\rhember\Documents\Data\Waste Wood\openings.tif')

#plt.matshow(burned)

#%% Rasterize land use legal

# Define paths
meta={}
meta['Paths']={}
meta['Paths']['Project']=r''
#meta['Paths']['Geospatial']=meta['Paths']['Project'] + '\\Geospatial'
meta['Paths']['Results']=r'C:\Users\rhember\Documents\Data\ForestInventory\Results\20210930'
meta['Paths']['VRI']=r'C:\Users\rhember\Documents\Data\ForestInventory\VRI\20210930'
meta['Paths']['Disturbances']=r'C:\Users\rhember\Documents\Data\ForestInventory\Disturbances\20210930'
meta['Paths']['LandUse']=r'C:\Users\rhember\Documents\Data\ForestInventory\LandUse\20210930'
meta['Paths']['Model Code']=r'C:\Users\rhember\Documents\Code_Python\fcgadgets\cbrunner'
meta['Paths']['Taz Datasets']=r'C:\Users\rhember\Documents\Data\Taz Datasets'

meta=invu.Load_LUTs(meta)

# Import TSA raster
zTSA=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Admin\tsa.tif')

# Input path to RESULTS database (downloaded from BC data catalogue)
pthin=r'C:\Users\rhember\Documents\Data\ForestInventory\LandUse\20210930\LandUse.gdb'
pthout=r'C:\Users\rhember\Documents\Data\BC1ha\LandUseLandCover'

fiona.listlayers(pthin)
lnam='RMP_PLAN_LEGAL_POLY_SVW'

list(meta['LUT']['LU L']['LEGAL_FEAT_OBJECTIVE'].keys())

# Open the shapefile
df=gpd.read_file(pthin,layer=lnam)

# Remove features with no geometry
df=df[df.geometry!=None]

fn='LEGAL_FEAT_OBJECTIVE'

for val in meta['LUT']['LU L']['LEGAL_FEAT_OBJECTIVE'].keys():
    print(val)
    #val='Conservation Lands'

    df0=df[df[fn]==val].copy()
    df0['dummy']=np.ones(len(df0))
    shapes=((geom,value) for geom, value in zip(df0['geometry'],df0['dummy']))

    z0=np.zeros(zTSA['Data'].shape,dtype=float)
    burned=features.rasterize(shapes=shapes,fill=0,out=z0,transform=zTSA['Transform'])

    val1=val.replace(":","_")
    val1=val1.replace("/","_")
    val1=val1.replace("<","Less Than")

    zOut=zTSA.copy()
    zOut['Data']=np.zeros(zTSA['Data'].shape,dtype='int16')
    zOut['Data'][burned>0]=1
    gis.SaveGeoTiff(zOut,pthout + '\\' + fn + '_' + val1 + '.tif')


#plt.matshow(zOut['Data'][::50,::50])

#fout=r'Z:\!Workgrp\Forest Carbon\Data\BC1ha\Disturbances' + '\\' + lnam + '_' + str(tv[iT]) + '.tif'
#gis.SaveGeoTiff(zOut,fout)

# Test
#plt.matshow(z['Data'])

#%% Rasterize forest cover reserves

# First get year of harvest from opening layer
metaO={}
metaO['Path']=r'C:\Users\rhember\Documents\Data\ForestInventory\Results\20220422\Results.gdb'
metaO['Layer']='RSLT_OPENING_SVW'; # fiona.listlayers(op['Path'])
metaO['crs']=[]
metaO['Keep Geom']='Off'
metaO['Select Openings']=np.array([])
metaO['SBC']=np.array([])
metaO['FSC']=np.array([])
metaO['ROI']=[]
metaO['gdf']=qgdb.Query_Openings(metaO,[])

# Import reserve layer
pthin=r'C:\Users\rhember\Documents\Data\ForestInventory\Results\20220422\Results.gdb'
fiona.listlayers(pthin)
lnam='RSLT_FOREST_COVER_RESERVE_SVW'

df=gpd.read_file(pthin,layer=lnam)
df=df[df.geometry!=None] # Remove features with no geometry
df['Year']=np.ones(len(df))
for i in range(len(df)):
    ind=np.where(metaO['gdf']['OPENING_ID']==df.loc[i,'OPENING_ID'])[0]
    if metaO['gdf']['DENUDATION_1_COMPLETION_DATE'][ind[0]]!=None:
        df['Year'][i]=np.array(metaO['gdf']['DENUDATION_1_COMPLETION_DATE'][ind[0]][0:4],dtype=float)
    elif metaO['gdf']['DENUDATION_2_COMPLETION_DATE'][ind[0]]!=None:
        df['Year'][i]=np.array(metaO['gdf']['DENUDATION_2_COMPLETION_DATE'][ind[0]][0:4],dtype=float)

pthout=r'C:\Users\rhember\Documents\Data\BC1ha\Results'
z=np.zeros(zRef['Data'].shape,dtype=float)
shapes=((geom,value) for geom,value in zip(df['geometry'],df['Year']))
burned=features.rasterize(shapes=shapes,fill=0,out=z,transform=zRef['Transform'])
u,cbc=gu.CountByCategories(z)

# Ones are reserves with no apparent harvest from RESULTS, check to see that it is not a big problem
ind=np.where(u==1)[0]; iAll=np.where(u>0)[0]; print(cbc[ind]/np.sum(cbc[iAll])*100)

z1=zRef.copy()
z1['Data']=z.astype(np.int8)
gis.SaveGeoTiff(z1,pthout + '\\' + lnam + '.tif')

#%% Rasterize consolidated cutblocks by year

# Input path to RESULTS database (downloaded from BC data catalogue)
pthin=r'C:\Users\rhember\Documents\Data\ForestInventory\Disturbances\20220422\Disturbances.gdb'
fiona.listlayers(pthin)
lnam='VEG_CONSOLIDATED_CUT_BLOCKS_SP'

# Open the shapefile
df=gpd.read_file(pthin,layer=lnam)

# Remove features with no geometry
df=df[df.geometry!=None]

tv=np.arange(1990,2022,1)

for iT in range(tv.size):
    df0=df[df.HARVEST_YEAR==tv[iT]].copy()
    shapes=((geom,value) for geom, value in zip(df0.geometry,df0.HARVEST_YEAR))

    z0=np.zeros(zRef['Data'].shape,dtype=float)
    burned=features.rasterize(shapes=shapes,fill=0,out=z0,transform=zRef['Transform'])

    zOut=zRef.copy()
    zOut['Data']=np.zeros(zRef['Data'].shape,dtype='int16')
    zOut['Data'][burned>0]=1
    zOut['Data']=zOut['Data'].astype('int8')
    fout=r'C:\Users\rhember\Documents\Data\BC1ha\Disturbances' + '\\' + lnam + '_' + str(tv[iT]) + '.tif'
    gis.SaveGeoTiff(zOut,fout)

#%% Rasterize wildfire occurrence by year

# Input path to RESULTS database (downloaded from BC data catalogue)
pthin=r'C:\Users\rhember\Documents\Data\ForestInventory\Disturbances\20211030\Disturbances.gdb'
fiona.listlayers(pthin)
lnam='PROT_HISTORICAL_FIRE_POLYS_SP'

# Open the shapefile
df=gpd.read_file(pthin,layer=lnam)

# Remove features with no geometry
df=df[df.geometry!=None]

zTSA=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Admin\tsa.tif')

tv=np.arange(2020,2020,1)

for iT in range(tv.size):
    df0=df[df.FIRE_YEAR==tv[iT]].copy()
    shapes=((geom,value) for geom, value in zip(df0.geometry,df0.FIRE_YEAR))

    z0=np.zeros(zTSA['Data'].shape,dtype=float)
    burned=features.rasterize(shapes=shapes,fill=0,out=z0,transform=zTSA.Transform)

    zOut=zTSA.copy()
    zOut['Data']=np.zeros(zTSA['Data'].shape,dtype='int16')
    zOut['Data'][burned>0]=1
    fout=r'Z:\!Workgrp\Forest Carbon\Data\BC1ha\Disturbances' + '\\' + lnam + '_' + str(tv[iT]) + '.tif'
    gis.SaveGeoTiff(zOut,fout)

# Test
#z=gis.OpenGeoTiff(r'Z:\!Workgrp\Forest Carbon\Data\BC1ha\Disturbances\test1.tif')
#plt.matshow(z['Data'])

#%% Rasterize AOS occurrence by year

# Input path to RESULTS database (downloaded from BC data catalogue)
pthin=r'C:\Users\rhember\Documents\Data\ForestInventory\Disturbances\20210930\Disturbances.gdb'
fiona.listlayers(pthin)
lnam='PEST_INFESTATION_POLY'
fnam='PEST_SPECIES_CODE'

# Open the shapefile
df=gpd.read_file(pthin,layer=lnam)

# Remove features with no geometry
#df=df[df.geometry!=None]
#df=df.reset_index()

# Save code list
uPest=np.unique(df.PEST_SPECIES_CODE)
#df0=pd.DataFrame(data=np.arange(1,u.size+1,1),columns=['ID'])
#df0['PEST_SPECIES_CODE']=u
#df0.to_excel(pthout + '\\' + lnam + '_PEST_SPECIES_CODE.xlsx')

# Save code list
#ind=np.where(df['PEST_SPECIES_CODE']=='IDW')[0]
u_Severity=np.unique(df.PEST_SEVERITY_CODE)
#u_Severity=np.array(['T','L','M','S','V','G'])
#df0=pd.DataFrame(data=np.arange(1,u_Severity.size+1,1),columns=['ID'])
#df0['PEST_SEVERITY_CODE']=u_Severity
#df0.to_excel(pthout + '\\' + lnam + '_PEST_SEVERITY_CODE.xlsx')

df['PSC']=np.zeros(len(df),dtype='float')
for i in range(len(u_Severity)):
    ind=np.where(df['PEST_SEVERITY_CODE']==u_Severity[i])[0]
    df.loc[ind,'PSC']=i+1

zTSA=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Admin\tsa.tif')

tv=np.arange(1950,2021,1)

#fcd=['IDL','IBM','IBD','IBS','IDW','DFL']
fcd=['IDW']

for iP in range(len(fcd)):
    for iT in range(tv.size):
        print(tv[iT])
        zOut=zTSA.copy()
        zOut['Data']=np.zeros(zTSA['Data'].shape,dtype='int16')

        df0=df[ (df[fnam]==fcd[iP]) & (df.CAPTURE_YEAR==tv[iT]) ].copy()
        df0=df0[df0.geometry!=None]

        if len(df0)>0:
            shapes=((geom,value) for geom, value in zip(df0.geometry,df0.PSC))
            z0=np.zeros(zTSA['Data'].shape,dtype=float)
            burned=features.rasterize(shapes=shapes,fill=0,out=z0,transform=zTSA['Transform'])
            zOut['Data']=burned.astype('int16')

        fout=r'C:\Users\rhember\Documents\Data\BC1ha\Disturbances' + '\\' + lnam + '_' + fcd[iP] + '_SeverityClass_' + str(tv[iT]) + '.tif'
        gis.SaveGeoTiff(zOut,fout)

# Test
#z=gis.OpenGeoTiff(fout)
#plt.matshow(z['Data'])

#%% Rasterize aerial spray treatment

gdf_spray=gpd.read_file(r'C:\Users\rhember\Documents\Data\Aerial Btk Spray\Processed\btk_spray_comp.geojson')

tv=np.arange(1950,2021,1)

zTSA=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Admin\tsa.tif')

for iT in range(tv.size):
    print(tv[iT])
    zOut=zTSA.copy()
    zOut['Data']=np.zeros(zTSA['Data'].shape,dtype='int16')

    df0=gdf_spray[ (gdf_spray['Year']==tv[iT]) ].copy()
    df0=df0[df0.geometry!=None]
    df0['Dummy']=1

    if len(df0)>0:
        shapes=((geom,value) for geom, value in zip(df0.geometry,df0['Dummy']))
        z0=np.zeros(zTSA['Data'].shape,dtype=float)
        burned=features.rasterize(shapes=shapes,fill=0,out=z0,transform=zTSA['Transform'])
        zOut['Data']=burned.astype('int16')

    fout=r'C:\Users\rhember\Documents\Data\BC1ha\Management' + '\\btk_spray_' + str(tv[iT]) + '.tif'
    gis.SaveGeoTiff(zOut,fout)

#%% Rasterize burn severity

# Input path to RESULTS database (downloaded from BC data catalogue)
pthin=r'C:\Users\rhember\Documents\Data\ForestInventory\Disturbances\20200430\Disturbances.gdb'
fiona.listlayers(pthin)

lnam='VEG_BURN_SEVERITY_SP'

# Open the shapefile
df=gpd.read_file(pthin,layer=lnam)

# Remove features with no geometry
df=df[df.geometry!=None]

# Output path
pthout=r'Z:\!Workgrp\Forest Carbon\Data\BC1ha\Disturbances'

# Save severity code list
#uBSR=np.unique(df.BURN_SEVERITY_RATING)
uBSR=np.array(['Unburned','Low','Medium','High','Unknown'])
df0=pd.DataFrame(data=np.arange(1,uBSR.size+1,1),columns=['ID'])
df0['Code']=uBSR
df0.to_excel(pthout + '\\' + lnam + '.xlsx')

# Add burn severity rating ID to dataframe
df['ID_BSR']=np.zeros(len(df),dtype='float')
for i in range(len(uBSR)):
    ind=np.where(df['BURN_SEVERITY_RATING']==uBSR[i])[0]
    df.loc[ind,'ID_BSR']=i+1

zTSA=gis.OpenGeoTiff(r'Z:\!Workgrp\Forest Carbon\Data\BC1ha\Admin\tsa.tif')

tv=np.arange(2017,2020,1)

for iT in range(tv.size):

    df0=df[df.FIRE_YEAR==tv[iT]].copy()
    shapes=((geom,value) for geom, value in zip(df0.geometry,df0.ID_BSR))

    z0=np.zeros(zTSA['Data'].shape,dtype=float)
    burned=features.rasterize(shapes=shapes,fill=0,out=z0,transform=zTSA.Transform)

    zOut=zTSA.copy()
    zOut['Data']=5*np.ones(zTSA['Data'].shape,dtype='int16')
    zOut['Data'][burned>0]=burned[burned>0]
    fout=r'Z:\!Workgrp\Forest Carbon\Data\BC1ha\Disturbances' + '\\' + lnam + '_' + str(tv[iT]) + '.tif'
    gis.SaveGeoTiff(zOut,fout)


# Open template raster file and copy it's metadata
rst=rasterio.open(r'Z:\!Workgrp\Forest Carbon\Data\BC1ha\VRI\age1.tif')
meta=rst.meta.copy()
meta.update(compress='lzw')

# Output path
pthout=r'Z:\!Workgrp\Forest Carbon\Data\BC1ha\Disturbances'

# Save severity code list
u=np.unique(df.BURN_SEVERITY_RATING)
df0=pd.DataFrame(data=np.arange(1,u.size+1,1),columns=['ID'])
df0['BURN_SEVERITY_RATING']=u
df0.to_excel(pthout + '\\' + lnam + '_BURN_SEVERITY_RATING.xlsx')

tv=np.arange(1950,2020,1)

#------------------------------------------------------------------------------
# Year
#------------------------------------------------------------------------------

BR=['High','Medium','Low']

for iBR in range(len(BR)):

    out=rasterio.open(pthout + '\\test.tif', 'w', **meta)

    zy1=np.zeros(rst.shape,dtype=float)
    zy2=np.zeros(rst.shape,dtype=float)
    zy3=np.zeros(rst.shape,dtype=float)
    zy4=np.zeros(rst.shape,dtype=float)
    #zy5=np.zeros(rst.shape,dtype=float)
    #zy6=np.zeros(rst.shape,dtype=float)
    #zy7=np.zeros(rst.shape,dtype=float)
    #zy8=np.zeros(rst.shape,dtype=float)
    for i in range(tv.size):
        print(tv[i])

        # Isolate BR of interest in ith year
        df0=df[(df.BURN_SEVERITY_RATING==BR[iBR]) & (df.FIRE_YEAR==tv[i])].copy()
        df0=df0.reset_index(drop=True)
        if len(df0)==0:
            continue

        # Get index to entries with geometry
        ind=[]
        for j in range(len(df0)):
            if df0.loc[j,'geometry']!=None:
                ind.append(j)
        if len(ind)==0:
            continue
        df0=df0.loc[ind]

        # Year
        shapes=((geom,value) for geom, value in zip(df0.geometry,df0.FIRE_YEAR))
        z0=np.zeros(rst.shape,dtype=float)
        burned=features.rasterize(shapes=shapes,fill=0,out=z0,transform=out.transform)

        ind=np.where((zy1[:]==0) & (burned[:]!=0))
        zy1[ind]=tv[i]
        ind=np.where((zy1[:]!=0) & (zy1[:]!=tv[i]) & (zy2[:]==0) & (burned[:]!=0))
        zy2[ind]=tv[i]
        ind=np.where((zy2[:]!=0) & (zy2[:]!=tv[i]) & (zy3[:]==0) & (burned[:]!=0))
        zy3[ind]=tv[i]
        ind=np.where((zy3[:]!=0) & (zy3[:]!=tv[i]) & (zy4[:]==0) & (burned[:]!=0))
        zy4[ind]=tv[i]
        #ind=np.where((zy4[:]!=0) & (zy4[:]!=tv[i]) & (zy5[:]==0) & (burned[:]!=0))
        #zy5[ind]=tv[i]
        #ind=np.where((zy5[:]!=0) & (zy5[:]!=tv[i]) & (zy6[:]==0) & (burned[:]!=0))
        #zy6[ind]=tv[i]
        #ind=np.where((zy6[:]!=0) & (zy6[:]!=tv[i]) & (zy7[:]==0) & (burned[:]!=0))
        #zy7[ind]=tv[i]
        #ind=np.where((zy7[:]!=0) & (zy7[:]!=tv[i]) & (zy8[:]==0) & (burned[:]!=0))
        #zy8[ind]=tv[i]
    out.close()

    for i in range(4):
        out=rasterio.open(pthout + '\\' + lnam + '_' + BR[iBR] + '_year' + str(i+1) + '.tif', 'w', **meta)
        exec('a=zy' + str(i+1) + '.astype(np.int16).copy()')
        out.write_band(1,a); out.close()

#------------------------------------------------------------------------------
# Severity
#------------------------------------------------------------------------------

out=rasterio.open(pthout + '\\test.tif', 'w', **meta)

zs1=np.zeros(rst.shape,dtype=float)
zs2=np.zeros(rst.shape,dtype=float)
zs3=np.zeros(rst.shape,dtype=float)
zs4=np.zeros(rst.shape,dtype=float)
zs5=np.zeros(rst.shape,dtype=float)
zs6=np.zeros(rst.shape,dtype=float)
zs7=np.zeros(rst.shape,dtype=float)
zs8=np.zeros(rst.shape,dtype=float)
for i in range(tv.size):
    print(tv[i])

    # Isolate pest of interest in ith year
    df0=df[(df[fnam]==fcd) & (df.CAPTURE_YEAR==tv[i])].copy()
    df0=df0.reset_index(drop=True)
    if len(df0)==0:
        continue

    # Get index to entries with geometry
    ind=[]
    for j in range(len(df0)):
        if df0.loc[j,'geometry']!=None:
            ind.append(j)
    if len(ind)==0:
        continue
    df0=df0.loc[ind]

    # Severity
    shapes=((geom,value) for geom, value in zip(df0.geometry,df0.PSC))
    z0=np.zeros(rst.shape,dtype=float)
    burned=features.rasterize(shapes=shapes,fill=0,out=z0,transform=out.transform)

    ind=np.where((zs1[:]==0) & (burned[:]!=0))
    zs1[ind]=burned[ind]
    ind=np.where((zs1[:]!=0) & (zs1[:]!=tv[i]) & (zs2[:]==0) & (burned[:]!=0))
    zs2[ind]=burned[ind]
    ind=np.where((zs2[:]!=0) & (zs2[:]!=tv[i]) & (zs3[:]==0) & (burned[:]!=0))
    zs3[ind]=burned[ind]
    ind=np.where((zs3[:]!=0) & (zs3[:]!=tv[i]) & (zs4[:]==0) & (burned[:]!=0))
    zs4[ind]=burned[ind]
    ind=np.where((zs4[:]!=0) & (zs4[:]!=tv[i]) & (zs5[:]==0) & (burned[:]!=0))
    zs5[ind]=burned[ind]
    ind=np.where((zs5[:]!=0) & (zs5[:]!=tv[i]) & (zs6[:]==0) & (burned[:]!=0))
    zs6[ind]=burned[ind]
    ind=np.where((zs6[:]!=0) & (zs6[:]!=tv[i]) & (zs7[:]==0) & (burned[:]!=0))
    zs7[ind]=burned[ind]
    ind=np.where((zs7[:]!=0) & (zs7[:]!=tv[i]) & (zs8[:]==0) & (burned[:]!=0))
    zs8[ind]=burned[ind]
out.close()

for i in range(8):
    out=rasterio.open(pthout + '\\' + fnam + '_' + fcd + '_sev' + str(i+1) + '.tif', 'w', **meta)
    exec('a=zs' + str(i+1) + '.astype(np.int16).copy()')
    out.write_band(1,a); out.close()

#%% Rasterize BCTS

# Import TSA raster
zTSA=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Admin\tsa.tif')

# Input path to RESULTS database (downloaded from BC data catalogue)
pthin=r'C:\Users\rhember\Documents\Data\ForestInventory\LandUse\20220422\LandUse.gdb'

lnam='BCTS_OPERATING_AREAS_SP'

# Open the shapefile
df=gpd.read_file(pthin,layer=lnam)

# Remove features with no geometry
df=df[df.geometry!=None]

df0=df.copy()
df0['dummy']=np.ones(len(df0))
shapes=((geom,value) for geom, value in zip(df0['geometry'],df0['dummy']))

z0=np.zeros(zTSA['Data'].shape,dtype=float)
burned=features.rasterize(shapes=shapes,fill=0,out=z0,transform=zTSA['Transform'])

zOut=zTSA.copy()
zOut['Data']=np.zeros(zTSA['Data'].shape,dtype='int16')
zOut['Data'][burned>0]=1
gis.SaveGeoTiff(zOut,r'C:\Users\rhember\Documents\Data\BC1ha\LandUseLandCover\bcts_op_area.tif')

#%% RASTERIZE VRI VARIABLES

# # *** Gave up on this- takes too long ***

# # Open template raster file and copy it's metadata
# rst=rasterio.open(r'Z:\!Workgrp\Forest Carbon\Data\BC1ha\VRI\age1.tif')
# meta=rst.meta.copy()
# meta.update(compress='lzw')

# # Initialize raster
# z=np.zeros(rst.shape,dtype=float)

# # Input path to RESULTS database (downloaded from BC data catalogue)
# pthin=r'Z:\!Workgrp\Forest Carbon\Data\VRI\20190501\VRI.gdb'

# # Look at the layers in geodatabase
# fiona.listlayers(pthin)

# lnam='VEG_COMP_LYR_R1_POLY'

# # Output path
# pthout=r'Z:\!Workgrp\Forest Carbon\Data\BC1ha\VRI'

# # Initialize new raster file
# out=rasterio.open(pthout + '\\test.tif', 'w', **meta)

# # Get unique values
# s=['0']*10000
# cnt=0
# with fiona.open(pthin,layer=lnam) as source:
#     for feat in source:
#         tmp=feat['properties']['BCLCS_LEVEL_2']
#         if tmp!=None:
#             s[cnt]=tmp
#         cnt=cnt+1
# us=np.unique(s)
# Code=['L','W']
# ID=[1,2]

# df=gpd.GeoDataFrame()
# with fiona.open(pthin,layer=lnam) as source:
#     for feat in source:
#         df0=gpd.GeoDataFrame.from_features([feat])
#         df1=gpd.GeoDataFrame()
#         df1.geometry=df0.geometry
#         df1['BCLCS_LEVEL_2']=df0['BCLCS_LEVEL_2']
#         df=pd.concat([df,df1])

#         iCode=int(np.where(Code==df.BCLCS_LEVEL_2.values)[0])
#         df.BCLCS_LEVEL_2=ID[iCode]
#         shapes=((geom,value) for geom, value in zip(df.geometry,df.BCLCS_LEVEL_2))
#         burned=features.rasterize(shapes=shapes,fill=0,out=z,transform=out.transform)
#         ind=np.where((burned[:]!=0))
#         z[ind]=ID[iCode]
#         break

#         coords0=feat['geometry']['coordinates']
#         for i in range(len(coords0)):
#             coords1=coords0[i]
#             for j in range(len(coords1)):
#                 coords2=coords1[j]

#                 poly=np.asarray(coords2)

# # Output path
# pthout=r'Z:\!Workgrp\Forest Carbon\Data\BC1ha\VRI'

# # Initialize new raster file
# out=rasterio.open(pthout + '\\test.tif', 'w', **meta)

# # Set what variable will be done
# flg=1

# z1=np.zeros(rst.shape,dtype=float)

# shapes=((geom,value) for geom, value in zip(df.geometry,df.FIRE_YEAR))

# z0=np.zeros(rst.shape,dtype=float)
# burned=features.rasterize(shapes=shapes,fill=0,out=z0,transform=out.transform)

# #ind=np.where((z1[:]==0) & (burned[:]!=0))
# #z1[ind]=tv[i]

# out.close()

# out=rasterio.open(pthout + '\\' + lnam + '_year1.tif', 'w', **meta)
# out.write_band(1,z1.astype(np.int16))
# out.close()


#%% Import climate from old BC1ha project in matlab

# Annual summary
z=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Terrain\BC_1ha_twi.tif')
z=gis.ClipRaster(z,[zS['xmin'],zS['xmax']],[zS['ymin'],zS['ymax']])
gis.SaveGeoTiff(z,r'C:\Users\rhember\Documents\Data\BC1ha\Terrain\bc1ha_twi.tif')
del z
garc.collect()

# Annual summary
z=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Climate\BC1ha_tmin_ann_norm_1971to2000_si_hist_v1_c.tif')
z=gis.ClipRaster(z,[zS['xmin'],zS['xmax']],[zS['ymin'],zS['ymax']])
gis.SaveGeoTiff(z,r'C:\Users\rhember\Documents\Data\BC1ha\Climate\BC1ha_tmin_ann_norm_1971to2000_si_hist_v1.tif')
del z; garc.collect()

# Monthly
for mo in range(12,13):
    z=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Climate\BC1ha_dwf_ann_norm_1971to2000_si_hist_v1_c.tif')
    z=gis.ClipRaster(z,[zS['minx'],zS['maxx']],[zS['miny'],zS['maxy']])
    gis.SaveGeoTiff(z,r'C:\Users\rhember\Documents\Data\BC1ha_rswd_mon_norm_1971to2000_si_hist_v1_' + str(mo) + '.tif')
    del z
    gc.collect()

plt.matshow(z.Data)




#%% Extract mean climate data by BGC zone

zBGC=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\VRI\becsz.tif')
lutBGC=gu.ReadExcel(r'C:\Users\rhember\Documents\Data\BC1ha\VRI\becz_lut.xlsx')

zMAT=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Climate\BC1ha_mat_norm_1971to2000_si_hist_v1.tif')
zMAT['Data']=zMAT['Data'].astype(float)/10

zWS=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Climate\BC1ha_ws_gs_norm_1971to2000_comp_hist_v1.tif')
zWS['Data']=zWS['Data'].astype(float)


lutBGC['MAT']=np.zeros(lutBGC['VALUE'].size)
lutBGC['WS']=np.zeros(lutBGC['VALUE'].size)
for i in range(lutBGC['VALUE'].size):
    ta=zMAT['Data'].flatten()
    ws=zWS['Data'].flatten()
    ind=np.where( (zBGC['Data'].flatten()==lutBGC['VALUE'][i]) & (ta>=-50) & (ws>=0) & (ws<=200) )[0]
    lutBGC['MAT'][i]=np.mean(ta[ind])
    lutBGC['WS'][i]=np.mean(ws[ind])

df=pd.DataFrame(lutBGC)
df.to_excel(r'C:\Users\rhember\Documents\Data\BC1ha\Climate\tmp.xlsx')


