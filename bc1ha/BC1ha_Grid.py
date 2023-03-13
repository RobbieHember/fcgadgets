
#%% IMPORT MODULES

import sys
import numpy as np
from osgeo import gdal
from osgeo import osr
import matplotlib.pyplot as plt
import geopandas as gpd
from rasterio.transform import from_origin
import pandas as pd
import scipy.io as spio
import fiona
import rasterio
from rasterio import features
from shapely.geometry import Point,Polygon
from shapely import geometry
import cv2

import fcgadgets.macgyver.utilities_general as gu
import fcgadgets.macgyver.utilities_gis as gis
from fcgadgets.macgyver import utilities_inventory as invu
import fcgadgets.macgyver.utilities_query_gdb as qgdb
import fcgadgets.bc1ha.bc1ha_utilities as u1ha

#%% Import look-up-tables

lut=u1ha.Import_BC1ha_LUTs()
gp=gu.SetGraphics('Manuscript')

#%% Define reference grid

fref=r'C:\Users\rhember\Documents\Data\BC1ha\Admin\BC_Land_Mask.tif'
zRef=gis.OpenGeoTiff(fref)

# # Create a grid that is consistent with what FAIB uses
# # *** Stopped using this because there is a bunch of unnecessary space ***
# # Load BC basemap to get projection information
# bm=gpd.read_file(r'C:\Users\rhember\Documents\Data\Basemaps\bc_land.shp')

# # Inputs (west,north,xsize,ysize) - received from Tyler Muhley
# transform=from_origin(159587.5,1748187.5,100,100)

# fin_Ref=r'C:\Users\rhember\Documents\Data\BC1ha\Admin\FAIB_Standard.tif'
# nd=rasterio.open(fin_Ref,'w',driver='GTiff',height=15744,width=17216,count=1,dtype=rasterio.float32,compress='lzw',crs=bm.crs,transform=transform)
# nd.write(nd,1) # This crashes but then it seems to work if you then close it.
# nd.close()
# #z=gis.OpenGeoTiff(pth)

#%% Prepare timber supply area

# Create TSA Key
# This creates a key between the raster values for each TSA, and the TSA names
# from a shapefile of TSA downloaded from the BC Data Cat.

def TSA_Key():

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

    return

# Digitize the boundary of TSAs
def TSA_DigitizeBoundaries():
    # Open dataframe containing TSA names and values for raster grid of TSAs
    zTSA=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Admin\tsa.tif')

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

    return

#%% Clip rasters to standard grid
# *** This also compresses files that come out of Arc crazy big ***

# DEM
fin=r'C:\Users\rhember\Documents\Data\BC1ha\Terrain\elevation.tif'
fout=r'C:\Users\rhember\Documents\Data\BC1ha\Terrain\elevation.tif'
gis.ClipToRaster_ByFile(fin,fout,fref)

# Aspect
fin=r'C:\Users\rhember\Documents\Data\BC1ha\Terrain\aspect.tif'
fout=r'C:\Users\rhember\Documents\Data\BC1ha\Terrain\aspect.tif'
gis.ClipToRaster_ByFile(fin,fout,fref)

# Slope
fin=r'C:\Users\rhember\Documents\Data\BC1ha\Terrain\slope.tif'
fout=r'C:\Users\rhember\Documents\Data\BC1ha\Terrain\slope.tif'
gis.ClipToRaster_ByFile(fin,fout,fref)

# BTM
fin=r'C:\Users\rhember\Documents\Data\BC1ha\LandUseLandCover\landuse.btm.tif'
fout=r'C:\Users\rhember\Documents\Data\BC1ha\LandUseLandCover\landuse.btm.tif'
gis.ClipToRaster_ByFile(fin,fout,fref)

# LC2
fin=r'C:\Users\rhember\Documents\Data\BC1ha\VRI\lc2.tif'
fout=r'C:\Users\rhember\Documents\Data\BC1ha\VRI\lc2.tif'
gis.ClipToRaster_ByFile(fin,fout,fref)

# BGC zone
fin=r'C:\Users\rhember\Documents\Data\BC1ha\VRI\becz.tif'
fout=r'C:\Users\rhember\Documents\Data\BC1ha\VRI\becz.tif'
gis.ClipToRaster_ByFile(fin,fout,fref)

# BGC subzone
fin=r'C:\Users\rhember\Documents\Data\BC1ha\VRI\becsz.tif'
fout=r'C:\Users\rhember\Documents\Data\BC1ha\VRI\becsz.tif'
gis.ClipToRaster_ByFile(fin,fout,fref)

# VRI age
fin=r'C:\Users\rhember\Documents\Data\BC1ha\VRI\proj_age_1.tif'
fout=r'C:\Users\rhember\Documents\Data\BC1ha\VRI\proj_age_1.tif'
gis.ClipToRaster_ByFile(fin,fout,fref)

# VRI SI
fin=r'C:\Users\rhember\Documents\Data\BC1ha\VRI\si.tif'
fout=r'C:\Users\rhember\Documents\Data\BC1ha\VRI\si.tif'
gis.ClipToRaster_ByFile(fin,fout,fref)

fin=r'C:\Users\rhember\Documents\Data\BC1ha\VRI\sphdead.tif'
fout=r'C:\Users\rhember\Documents\Data\BC1ha\VRI\sphdead.tif'
gis.ClipToRaster_ByFile(fin,fout,fref)

fin=r'C:\Users\rhember\Documents\Data\BC1ha\VRI\crownc.tif'
fout=r'C:\Users\rhember\Documents\Data\BC1ha\VRI\crownc.tif'
gis.ClipToRaster_ByFile(fin,fout,fref)

fin=r'C:\Users\rhember\Documents\Data\BC1ha\SPL\Site_Prod_Pl.tif'
fout=r'C:\Users\rhember\Documents\Data\BC1ha\SPL\Site_Prod_Pl.tif'
gis.ClipToRaster_ByFile(fin,fout,fref)

# Mean annual temp
fin=r'C:\Users\rhember\Documents\Data\BC1ha\Climate\BC1ha_mat_norm_1971to2000_si_hist_v1.tif'
fout=r'C:\Users\rhember\Documents\Data\BC1ha\Climate\BC1ha_mat_norm_1971to2000_si_hist_v1.tif'
gis.ClipToRaster_ByFile(fin,fout,fref)

# Soil water content
fin=r'C:\Users\rhember\Documents\Data\BC1ha\Climate\BC1ha_ws_gs_norm_1971to2000_comp_hist_v1.tif'
fout=r'C:\Users\rhember\Documents\Data\BC1ha\Climate\BC1ha_ws_gs_norm_1971to2000_comp_hist_v1.tif'
gis.ClipToRaster_ByFile(fin,fout,fref)

fin=r'C:\Users\rhember\Documents\Data\BC1ha\Disturbances\PROT_HISTORICAL_FIRE_POLYS_SP_2017.tif'
fout=fin
gis.ClipToRaster_ByFile(fin,fout,fref)

fin=r'C:\Users\rhember\Documents\Data\BC1ha\Disturbances\VEG_CONSOLIDATED_CUT_BLOCKS_SP_year1a.tif'
fout=r'C:\Users\rhember\Documents\Data\BC1ha\Disturbances\VEG_CONSOLIDATED_CUT_BLOCKS_SP_year1.tif'
gis.ClipToRaster_ByFile(fin,fout,fref)

# Ownershp
fin=r'C:\Users\rhember\Documents\Data\BC1ha\LandUseLandCover\f_own.tif'
fout=r'C:\Users\rhember\Documents\Data\BC1ha\LandUseLandCover\f_own.tif'
gis.ClipToRaster_ByFile(fin,fout,fref)

# Stock type
fin=r'C:\Users\rhember\Documents\Data\BC1ha\Results\stocktype.tif'
fout=r'C:\Users\rhember\Documents\Data\BC1ha\Results\stocktype.tif'
gis.ClipToRaster_ByFile(fin,fout,fref)

#%% BC Land Mask

df=gpd.read_file(r'C:\Users\rhember\Documents\Data\Basemaps\Basemaps.gdb',layer='NRC_POLITICAL_BOUNDARIES_1M_SP')

# Remove features with no geometry
df=df[df.geometry!=None]

df['ID']=1
shapes=((geom,value) for geom, value in zip(df.geometry,df.ID))
z=np.zeros(zRef['Data'].shape,dtype=float)
burned=features.rasterize(shapes=shapes,fill=0,out=z,transform=zRef['Transform'])

plt.close('all')
plt.matshow(burned)

zOut=zRef.copy()
zOut['Data']=np.zeros(zRef['Data'].shape,dtype='int8')
zOut['Data'][burned>0]=1
zOut['Data']=zOut['Data'].astype('int8')

gis.SaveGeoTiff(zOut,r'C:\Users\rhember\Documents\Data\BC1ha\Admin\BC_Land_Mask.tif')

#%% Crown Land Mask

z0=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\LandUseLandCover\f_own.tif')
z1=zRef.copy()
z1['Data']=np.zeros(zRef['Data'].shape,dtype='int8')
ind=np.where( (zRef['Data']==1) & (z0['Data']>=9) )
z1['Data'][ind]=1
plt.matshow(z1['Data'])
gis.SaveGeoTiff(z1,r'C:\Users\rhember\Documents\Data\BC1ha\Admin\CrownForestLandMask.tif')

#%% BGCzone / NDT Combination

zBGC=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\VRI\becz.tif')
zNDT=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\VRI\ndt.tif')

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
df.to_excel(r'C:\Users\rhember\Documents\Data\BC1ha\VRI\lut_bgcz_ndt_combo.xlsx')

z1=zRef.copy()
z1['Data']=z
gis.SaveGeoTiff(z1,r'C:\Users\rhember\Documents\Data\BC1ha\VRI\bgcz_ndt_combo.tif')

#%% Land Cover Class

zRef=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Admin\BC_Land_Mask.tif')
zLC2=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\VRI\lc2.tif')
zLC4=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\VRI\lc4.tif')

z=zRef.copy()
z['Data']=np.zeros(z['Data'].shape,dtype='int8')

lab=['Forest Land','Shrubland','Herbs','Bryoids','Other']

A=np.zeros(len(lab))
ind=np.where( (zLC4['Data']==lut['lc4']['Treed Conifer']) | (zLC4['Data']==lut['lc4']['Treed Mixed']) | (zLC4['Data']==lut['lc4']['Treed Broadleaf']) )
A[0]=ind[0].size/1e6; z['Data'][ind]=1

ind=np.where( (zLC4['Data']==lut['lc4']['Shrub Low']) | (zLC4['Data']==lut['lc4']['Shrub Tall']) )
A[1]=ind[0].size/1e6; z['Data'][ind]=2

ind=np.where( (zLC4['Data']==lut['lc4']['Herb']) | (zLC4['Data']==lut['lc4']['Herb Gramanoid']) | (zLC4['Data']==lut['lc4']['Herb Forbs']) )
A[2]=ind[0].size/1e6; z['Data'][ind]=3

ind=np.where( (zLC4['Data']==lut['lc4']['Bryoid']) | (zLC4['Data']==lut['lc4']['Bryoid Moss']) | (zLC4['Data']==lut['lc4']['Bryoid Lichen']) )
A[3]=ind[0].size/1e6; z['Data'][ind]=4

ind=np.where( (zLC4['Data']==lut['lc4']['Rock ']) | (zLC4['Data']==lut['lc4']['Exposed Land']) | (zLC4['Data']==lut['lc4']['Snow and Ice']) )
A[4]=ind[0].size/1e6; z['Data'][ind]=5

# plt.matshow(z['Data'])
gis.SaveGeoTiff(z,r'C:\Users\rhember\Documents\Data\BC1ha\LandUseLandCover\LandCoverClass1.tif')

# Plot bar chart of area
x=np.arange(0,A.size,1)
plt.close('all');
fig,ax=plt.subplots(1,figsize=gu.cm2inch(10,6.5));
ax.bar(x,A)
ax.set(xticks=x,xticklabels=lab,ylabel='Area (Mha)')
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=gp['tickl'])

# Compare area treed from L2 and L4
ind=np.where( (zLC2['Data']==lut['lc2']['Treed']) )
A_treed_lc2=ind[0].size
print(A_treed_lc2/1e6)
print(A[0])

#%% Tree Density Class
# *** Ensure harvested areas are listed as dense ***

zLC=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\LandUseLandCover\LandCoverClass1.tif')
zLC5=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\VRI\lc5.tif')
zH=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Disturbances\VEG_CONSOLIDATED_CUT_BLOCKS_SP_All.tif')
zH2=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Results\Harvest_Mask_FromRESULTS.tif')

ind=np.where( (zH['Data']==1) | (zH2['Data']==1) )
A_H=ind[0].size/1e6
ind=np.where( (zLC['Data']==1) )
A_F=ind[0].size/1e6

z=np.zeros(zLC5['Data'].shape,dtype='int8')
ind=np.where( (zLC['Data']==1) & (zLC5['Data']==lut['lc5']['Sparse']) ); z[ind]=1
ind=np.where( (zLC['Data']==1) & (zLC5['Data']==lut['lc5']['Open']) ); z[ind]=2
ind=np.where( (zLC['Data']==1) & (zLC5['Data']==lut['lc5']['Dense']) ); z[ind]=3
ind=np.where( (zH['Data']==1) | (zH2['Data']==1) ); z[ind]=3
ind=np.where( (zRef['Data']==1) & (z==0) ); z[ind]=4
ind=np.where( (zRef['Data']==0) ); z[ind]=5
plt.matshow(z)

z1=zRef.copy()
z1['Data']=z
gis.SaveGeoTiff(z1,r'C:\Users\rhember\Documents\Data\BC1ha\LandUseLandCover\TreeDensityClass.tif')

# With shrubs and grasses
z=np.zeros(zLC4['Data'].shape,dtype='int8')
ind=np.where( (zLC['Data']==lut['lc4']['Shrub Low']) | (zLC4['Data']==lut['lc4']['Shrub Tall']) ); z[ind]=4
ind=np.where( (zLC['Data']==lut['lc4']['Herb Gramanoid']) | (zLC4['Data']==lut['lc4']['Herb Forbs']) ); z[ind]=5
ind=np.where( (zLC['Data']==1) & (zLC5['Data']==lut['lc5']['Sparse']) ); z[ind]=1
ind=np.where( (zLC['Data']==1) & (zLC5['Data']==lut['lc5']['Open']) ); z[ind]=2
ind=np.where( (zLC['Data']==1) & (zLC5['Data']==lut['lc5']['Dense']) ); z[ind]=3
ind=np.where( (zH['Data']==1) ); z[ind]=3
ind=np.where( (zRef['Data']==1) & (z==0) ); z[ind]=6
ind=np.where( (zRef['Data']==0) ); z[ind]=7
plt.matshow(z[0::3,0::3])

z1=zRef.copy()
z1['Data']=z
gis.SaveGeoTiff(z1,r'C:\Users\rhember\Documents\Data\BC1ha\LandUseLandCover\TreeDensityClass_WithShrubsGrasses.tif')

#%% Rasterize select openings

def RasterizeSelectOpenings():

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

    return

#%% Rasterize natural disturbance type

# Open the shapefile
df=gpd.read_file(r'C:\Users\rhember\Documents\Data\ForestInventory\VRI\20220404\VRI.gdb',layer='BEC_NATURAL_DISTURBANCE_SV')

# Remove features with no geometry
df=df[df.geometry!=None]

df['ID']=np.zeros(df['NATURAL_DISTURBANCE_TYPE_CODE'].size)
for i in range(1,5,1):
    df['ID'][df['NATURAL_DISTURBANCE_TYPE_CODE']=='NDT' + str(i)]=i

shapes=((geom,value) for geom, value in zip(df['geometry'],df['ID']))
z=np.zeros(zRef['Data'].shape,dtype=float)
burned=features.rasterize(shapes=shapes,fill=0,out=z,transform=zRef['Transform'])

z1=zRef.copy()
z1['Data']=z.astype('int16')
gis.SaveGeoTiff(z1,r'C:\Users\rhember\Documents\Data\BC1ha\VRI\ndt.tif')


#%% Rasterize distance from roads

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
gis.SaveGeoTiff(z1,r'C:\Users\rhember\Documents\Data\BC1ha\Terrain\DistanceFromRoads.tif')

#z=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Terrain\DistanceFromRoads.tif')
#plt.matshow(z['Data'])

#%% Rasterize distance from timber facilities

#fiona.listlayers(r'C:\Users\rhember\Documents\Data\ForestInventory\Infrastructure.gdb')
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
gis.SaveGeoTiff(z1,r'C:\Users\rhember\Documents\Data\BC1ha\Management\DistanceFromForestryFacility.tif')

#%% Rasterize old growth deferrals

# Define paths
meta={}
meta['Paths']={}
meta['Paths']['Project']=r''
meta['Paths']['LandUse']=r'C:\Users\rhember\Documents\Data\ForestInventory\LandUse\20220422'
meta['Paths']['Model Code']=r'C:\Users\rhember\Documents\Code_Python\fcgadgets\cbrunner'

meta=invu.Load_LUTs(meta)

pthin=meta['Paths']['LandUse'] + '\\LandUse.gdb'
pthout=r'C:\Users\rhember\Documents\Data\BC1ha\LandUseLandCover'

fiona.listlayers(pthin)
lnam='OGSR_TAP_PRIORITY_DEF_AREA_SP'

list(lut['LU L']['LEGAL_FEAT_OBJECTIVE'].keys())

# Open the shapefile
df=gpd.read_file(pthin,layer=lnam)

# Remove features with no geometry
df=df[df.geometry!=None]

shapes=((geom,value) for geom, value in zip(df['geometry'],df['OGSR_TPDA_SYSID']))
z=np.zeros(zRef['Data'].shape,dtype=float)
burned=features.rasterize(shapes=shapes,fill=0,out=z,transform=zRef['Transform'])

z1=zRef.copy()
z1['Data']=z.astype('int32')
gis.SaveGeoTiff(z1,pthout + '\\' + lnam + '.tif')

#%% Rasterize parks

# Define paths
meta={}
meta['Paths']={}
meta['Paths']['Project']=r''
meta['Paths']['LandUse']=r'C:\Users\rhember\Documents\Data\ForestInventory\LandUse\20220422'
meta['Paths']['Model Code']=r'C:\Users\rhember\Documents\Code_Python\fcgadgets\cbrunner'

meta=invu.Load_LUTs(meta)

pthin=meta['Paths']['LandUse'] + '\\LandUse.gdb'
pthout=r'C:\Users\rhember\Documents\Data\BC1ha\LandUseLandCover'

fiona.listlayers(pthin)
lnam='TA_PARK_ECORES_PA_SVW'

# Open the shapefile
df=gpd.read_file(pthin,layer=lnam)

# Remove features with no geometry
df=df[df.geometry!=None]

shapes=((geom,value) for geom, value in zip(df['geometry'],df['ADMIN_AREA_SID']))
z=np.zeros(zRef['Data'].shape,dtype=float)
burned=features.rasterize(shapes=shapes,fill=0,out=z,transform=zRef['Transform'])

z1=zRef.copy()
z1['Data']=z.astype('int32')
gis.SaveGeoTiff(z1,pthout + '\\' + lnam + '.tif')

#%% Rasterize OGMAs

meta={}
meta['Paths']={}
meta['Paths']['Project']=r''
meta['Paths']['LandUse']=r'C:\Users\rhember\Documents\Data\ForestInventory\LandUse\20220422'
meta['Paths']['Model Code']=r'C:\Users\rhember\Documents\Code_Python\fcgadgets\cbrunner'

meta=invu.Load_LUTs(meta)

pthin=meta['Paths']['LandUse'] + '\\LandUse.gdb'
pthout=r'C:\Users\rhember\Documents\Data\BC1ha\LandUseLandCover'

fiona.listlayers(pthin)
lnam='RMP_OGMA_LEGAL_ALL_SVW'

# Open the shapefile
df=gpd.read_file(pthin,layer=lnam)

# Remove features with no geometry
df=df[df.geometry!=None]

shapes=((geom,value) for geom, value in zip(df['geometry'],df['LEGAL_OGMA_INTERNAL_ID']))
z=np.zeros(zRef['Data'].shape,dtype=float)
burned=features.rasterize(shapes=shapes,fill=0,out=z,transform=zRef['Transform'])

z1=zRef.copy()
z1['Data']=z.astype('int32')
gis.SaveGeoTiff(z1,pthout + '\\' + lnam + '.tif')

#%% Rasterize forest cover reserves

# Define paths
meta={}
meta['Paths']={}
meta['Paths']['Project']=r''
meta['Paths']['Results']=r'C:\Users\rhember\Documents\Data\ForestInventory\Results\20220422'
meta['Paths']['Model Code']=r'C:\Users\rhember\Documents\Code_Python\fcgadgets\cbrunner'

meta=invu.Load_LUTs(meta)

pthin=meta['Paths']['Results'] + '\\Results.gdb'
pthout=r'C:\Users\rhember\Documents\Data\BC1ha\LandUseLandCover'

fiona.listlayers(pthin)
lnam='RSLT_FOREST_COVER_RESERVE_SVW'

# Open the shapefile
df=gpd.read_file(pthin,layer=lnam)

# Remove features with no geometry
df=df[df.geometry!=None]
df['ID']=np.arange(0,len(df),1)

shapes=((geom,value) for geom, value in zip(df['geometry'],df['ID']))
z=np.zeros(zRef['Data'].shape,dtype=float)
burned=features.rasterize(shapes=shapes,fill=0,out=z,transform=zRef['Transform'])

z1=zRef.copy()
z1['Data']=z.astype('int32')
gis.SaveGeoTiff(z1,pthout + '\\' + lnam + '.tif')

#%% Rasterize protected lands

# Define paths
meta={}
meta['Paths']={}
meta['Paths']['Project']=r''
meta['Paths']['LandUse']=r'C:\Users\rhember\Documents\Data\ForestInventory\LandUse\20220422'
meta['Paths']['Model Code']=r'C:\Users\rhember\Documents\Code_Python\fcgadgets\cbrunner'

meta=invu.Load_LUTs(meta)

pthin=meta['Paths']['LandUse'] + '\\LandUse.gdb'
pthout=r'C:\Users\rhember\Documents\Data\BC1ha\LandUseLandCover'

fiona.listlayers(pthin)
lnam='TA_PROTECTED_LANDS_SV'

list(lut['LU L']['LEGAL_FEAT_OBJECTIVE'].keys())

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

# Input path to RESULTS database (downloaded from BC data catalogue)
pthin=r'C:\Users\rhember\Documents\Data\ForestInventory\LandUse\20210930\LandUse.gdb'
pthout=r'C:\Users\rhember\Documents\Data\BC1ha\LandUseLandCover'

fiona.listlayers(pthin)
lnam='RMP_PLAN_LEGAL_POLY_SVW'

list(lut['LU L']['LEGAL_FEAT_OBJECTIVE'].keys())

# Open the shapefile
df=gpd.read_file(pthin,layer=lnam)

# Remove features with no geometry
df=df[df.geometry!=None]

fn='LEGAL_FEAT_OBJECTIVE'

for val in lut['LU L']['LEGAL_FEAT_OBJECTIVE'].keys():
    print(val)
    #val='Conservation Lands'

    df0=df[df[fn]==val].copy()
    df0['dummy']=np.ones(len(df0))
    shapes=((geom,value) for geom, value in zip(df0['geometry'],df0['dummy']))

    z0=np.zeros(zRef['Data'].shape,dtype=float)
    burned=features.rasterize(shapes=shapes,fill=0,out=z0,transform=zRef['Transform'])

    val1=val.replace(":","_")
    val1=val1.replace("/","_")
    val1=val1.replace("<","Less Than")

    zOut=zRef.copy()
    zOut['Data']=np.zeros(zRef['Data'].shape,dtype='int16')
    zOut['Data'][burned>0]=1
    gis.SaveGeoTiff(zOut,pthout + '\\' + fn + '_' + val1 + '.tif')

#plt.matshow(zOut['Data'][::50,::50])

#fout=r'Z:\!Workgrp\Forest Carbon\Data\BC1ha\Disturbances' + '\\' + lnam + '_' + str(tv[iT]) + '.tif'
#gis.SaveGeoTiff(zOut,fout)

# Test
#plt.matshow(z['Data'])

#%% Harvest Mask from CC

# Open the shapefile
df=gpd.read_file(r'C:\Users\rhember\Documents\Data\ForestInventory\Disturbances\20220422\Disturbances.gdb',layer='VEG_CONSOLIDATED_CUT_BLOCKS_SP')

# Remove features with no geometry
df=df[df.geometry!=None]

shapes=((geom,value) for geom, value in zip(df.geometry,df.HARVEST_YEAR))
z=np.zeros(zRef['Data'].shape,dtype=float)
burned=features.rasterize(shapes=shapes,fill=0,out=z,transform=zRef['Transform'])

zH=zRef.copy()
zH['Data']=np.zeros(zRef['Data'].shape,dtype='int8')
zH['Data'][burned>0]=1
zH['Data']=zH['Data'].astype('int8')

gis.SaveGeoTiff(zH,r'C:\Users\rhember\Documents\Data\BC1ha\Disturbances\VEG_CONSOLIDATED_CUT_BLOCKS_SP_All.tif')

#%% Rasterize consolidated cutblocks (by year)

# Input path to RESULTS database (downloaded from BC data catalogue)
pthin=r'C:\Users\rhember\Documents\Data\ForestInventory\Disturbances\20220422\Disturbances.gdb'
fiona.listlayers(pthin)
lnam='VEG_CONSOLIDATED_CUT_BLOCKS_SP'

# Open the shapefile
df=gpd.read_file(pthin,layer=lnam)

# Remove features with no geometry
df=df[df.geometry!=None]

tv=np.arange(1950,2022,1)

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

#%% Rasterize harvest mask from RESULTS (Mask all)

# Import opening layer
gdf_bm=gpd.read_file(r'C:\Users\rhember\Documents\Data\Basemaps\Basemaps.gdb',layer='NRC_POLITICAL_BOUNDARIES_1M_SP')
metaOP={}
metaOP['Path']=r'C:\Users\rhember\Documents\Data\ForestInventory\Results\20220422\Results.gdb'
metaOP['Layer']='RSLT_OPENING_SVW'; # fiona.listlayers(op['Path'])
metaOP['crs']=gdf_bm.crs
metaOP['Keep Geom']='On'
metaOP['Select Openings']=np.array([])
metaOP['Denudation']=np.array(['L','S'])
metaOP['SBC']=np.array([])
metaOP['FSC']=np.array([])
metaOP['ROI']=[]
metaOP['gdf']=qgdb.Query_Openings(metaOP,[])

# Remove features with no geometry
metaOP['gdf']=metaOP['gdf'][metaOP['gdf'].geometry!=None]
metaOP['gdf']['ID']=np.ones(len(metaOP['gdf']))

shapes=((geom,value) for geom, value in zip(metaOP['gdf']['geometry'],metaOP['gdf']['ID']))
z=np.zeros(zRef['Data'].shape,dtype=float)
burned=features.rasterize(shapes=shapes,fill=0,out=z,transform=zRef['Transform'])

z1=zRef.copy()
z1['Data']=z.astype('int32')
gis.SaveGeoTiff(z1,r'C:\Users\rhember\Documents\Data\BC1ha\Disturbances\Harvest_FromRESULTS_Mask.tif')

#%% Rasterize harvest from RESULTS (by year)

# Import opening layer
gdf_bm=gpd.read_file(r'C:\Users\rhember\Documents\Data\Basemaps\Basemaps.gdb',layer='NRC_POLITICAL_BOUNDARIES_1M_SP')
metaOP={}
metaOP['Path']=r'C:\Users\rhember\Documents\Data\ForestInventory\Results\20220422\Results.gdb'
metaOP['Layer']='RSLT_OPENING_SVW'; # fiona.listlayers(op['Path'])
metaOP['crs']=gdf_bm.crs
metaOP['Keep Geom']='On'
metaOP['Select Openings']=np.array([])
metaOP['Denudation']=np.array(['L','S'])
metaOP['SBC']=np.array([])
metaOP['FSC']=np.array([])
metaOP['ROI']=[]
metaOP['gdf']=qgdb.Query_Openings(metaOP,[])

# Remove features with no geometry
metaOP['gdf']=metaOP['gdf'][metaOP['gdf'].geometry!=None]
metaOP['gdf']['ID']=np.ones(len(metaOP['gdf']))

tv=np.arange(1950,2022,1)
for iT in range(tv.size):
    df0=metaOP['gdf'][metaOP['gdf']['Year']==tv[iT]].copy()
    shapes=((geom,value) for geom, value in zip(df0['geometry'],df0['ID']))
    z=np.zeros(zRef['Data'].shape,dtype=float)
    burned=features.rasterize(shapes=shapes,fill=0,out=z,transform=zRef['Transform'])
    z1=zRef.copy()
    z1['Data']=z.astype('int8')
    gis.SaveGeoTiff(z1,r'C:\Users\rhember\Documents\Data\BC1ha\Disturbances\Harvest_FromRESULTS_' + str(tv[iT]) + '.tif')

#%% Rasterize forest cover stocking type from RESULTS

# Import opening layer
gdf_bm=gpd.read_file(r'C:\Users\rhember\Documents\Data\Basemaps\Basemaps.gdb',layer='NRC_POLITICAL_BOUNDARIES_1M_SP')
metaOP={}
metaOP['Path']=r'C:\Users\rhember\Documents\Data\ForestInventory\Results\20220422\Results.gdb'
metaOP['Layer']='RSLT_FOREST_COVER_INV_SVW'; # fiona.listlayers(op['Path'])
metaOP['crs']=gdf_bm.crs
metaOP['Keep Geom']='On'
metaOP['Select Openings']=np.array([])
metaOP['SBC']=np.array([])
metaOP['FSC']=np.array([])
metaOP['ROI']=[]
metaOP['gdf']=qgdb.Query_Openings(metaOP,[])

# Remove features with no geometry
metaOP['gdf']=metaOP['gdf'][metaOP['gdf'].geometry!=None]

# Add ID
d=gu.ReadExcel(r'C:\Users\rhember\Documents\Code_Python\fcgadgets\cbrunner\Parameters\Parameters_ForestCover_StockingType.xlsx')
metaOP['gdf']['ID']=np.ones(len(metaOP['gdf']))
for i in range(d['ID'].size):
    ind=np.where(metaOP['gdf']['STOCKING_TYPE_CODE']==d['Code'])[0]
    metaOP['gdf']['ID'][ind]=d['ID'][i]

shapes=((geom,value) for geom, value in zip(metaOP['gdf']['geometry'],metaOP['gdf']['ID']))
z=np.zeros(zRef['Data'].shape,dtype=float)
burned=features.rasterize(shapes=shapes,fill=0,out=z,transform=zRef['Transform'])

z1=zRef.copy()
z1['Data']=z.astype('int8')
#plt.matshow(z1['Data'])
gis.SaveGeoTiff(z1,r'C:\Users\rhember\Documents\Data\BC1ha\Results\StockingTypeArt_FromRESULTS_Mask.tif')

#%% Rasterize planting from RESULTS

gdf_bm=gpd.read_file(r'C:\Users\rhember\Documents\Data\Basemaps\Basemaps.gdb',layer='NRC_POLITICAL_BOUNDARIES_1M_SP')

# Start with planting with spatial from RESULTS
metaATS={}
metaATS['Path']=r'C:\Users\rhember\Documents\Data\ForestInventory\Results\20220422\Results.gdb'
metaATS['Layer']='RSLT_ACTIVITY_TREATMENT_SVW'; # fiona.listlayers(metaATS['Path'])
metaATS['crs']=gdf_bm.crs
metaATS['Keep Geom']='On'
metaATS['Select Openings']=np.array([])
metaATS['SBC']=np.array(['PL'])
metaATS['FSC']=np.array([])
metaATS['ROI']=[]
metaATS['gdf']=qgdb.Query_Openings(metaATS,[])

# Remove features with no geometry
metaATS['gdf']=metaATS['gdf'][metaATS['gdf'].geometry!=None]

metaATS['gdf']['ID']=np.ones(len(metaATS['gdf']))
tv=np.arange(1950,2022,1)
for iT in range(tv.size):
    z=np.zeros(zRef['Data'].shape,dtype=float)
    df0=metaATS['gdf'][metaATS['gdf']['Year']==tv[iT]].copy()
    if len(df0)>0:
        shapes=((geom,value) for geom, value in zip(df0['geometry'],df0['ID']))
        burned=features.rasterize(shapes=shapes,fill=0,out=z,transform=zRef['Transform'])
    z1=zRef.copy()
    z1['Data']=z.astype('int8')
    gis.SaveGeoTiff(z1,r'C:\Users\rhember\Documents\Data\BC1ha\Results\Planting_FromRESULTS_' + str(tv[iT]) + '.tif')

# Add areas where FC is artificial

metaAT={}
metaAT['Path']=r'C:\Users\rhember\Documents\Data\ForestInventory\Results\20220422\Results.gdb'
metaAT['Layer']='RSLT_ACTIVITY_TREATMENT_SVW'; # fiona.listlayers(metaAT['Path'])
metaAT['crs']=gdf_bm.crs
metaAT['Keep Geom']='Off'
metaAT['Select Openings']=np.array([])
metaAT['SBC']=np.array(['PL'])
metaAT['FSC']=np.array([])
metaAT['ROI']=[]
metaAT['gdf']=qgdb.Query_Openings(metaAT,[])

# Import all openings from opening layer
metaOP={}
metaOP['Path']=r'C:\Users\rhember\Documents\Data\ForestInventory\Results\20220422\Results.gdb'
metaOP['Layer']='RSLT_OPENING_SVW'; # fiona.listlayers(op['Path'])
metaOP['crs']=gdf_bm.crs
metaOP['Keep Geom']='Off'
metaOP['Select Openings']=np.array([])
metaOP['SBC']=np.array([])
metaOP['FSC']=np.array([])
metaOP['ROI']=[]
metaOP['gdf']=qgdb.Query_Openings(metaOP,[])

metaOPS={}
metaOPS['Path']=r'C:\Users\rhember\Documents\Data\ForestInventory\Results\20220422\Results.gdb'
metaOPS['Layer']='RSLT_OPENING_SVW'; # fiona.listlayers(op['Path'])
metaOPS['crs']=gdf_bm.crs
metaOPS['Keep Geom']='On'
metaOPS['Select Openings']=np.array([])
metaOPS['SBC']=np.array([])
metaOPS['FSC']=np.array([])
metaOPS['ROI']=[]
metaOPS['Drop Props']='On'
metaOPS['gdf']=qgdb.Query_Openings(metaOPS,[])

# Rasterize opening ID
zOP=np.zeros(zRef['Data'].shape,dtype=float)
shapes=((geom,value) for geom, value in zip(metaOPS['gdf']['geometry'],metaOPS['gdf']['OPENING_ID']))
burned=features.rasterize(shapes=shapes,fill=0,out=zOP,transform=zRef['Transform'])
#z1=zRef.copy()
#z1['Data']=z.astype('int8')
#gis.SaveGeoTiff(z1,r'C:\Users\rhember\Documents\Data\BC1ha\Results\Planting_FromRESULTS_' + str(tv[iT]) + '.tif')
#plt.matshow(zOP)

# Import FC layer
metaFC={}
metaFC['Path']=r'C:\Users\rhember\Documents\Data\ForestInventory\Results\20220422\Results.gdb'
metaFC['Layer']='RSLT_FOREST_COVER_INV_SVW'; # fiona.listlayers(op['Path'])
metaFC['crs']=gdf_bm.crs
metaFC['Keep Geom']='Off'
metaFC['Select Openings']=np.array([])
metaFC['SBC']=np.array([])
metaFC['FSC']=np.array([])
metaFC['ROI']=[]
metaFC['gdf']=qgdb.Query_Openings(metaFC,[])

zArt=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Results\StockingTypeArt_FromRESULTS_Mask.tif')

# Reduce the size of rasters
indO1=np.where(zOP!=0)
zOP1=zOP[indO1]
zArt1=zArt['Data'][indO1]

# Index to planting and year
tv=np.arange(1920,2023,1)
d={}
for iT in range(tv.size):
    d[tv[iT]]={'Idx':[],'SPH':[]}

for i in range(metaOPS['gdf']['OPENING_ID'].size):
    print(i)
    indO2=np.where(zOP1==metaOPS['gdf']['OPENING_ID'][i])[0]

    #indAT=np.where( (metaAT['gdf']['OPENING_ID']==metaOPS['gdf']['OPENING_ID'][i]) & (metaAT['gdf']['SILV_TECHNIQUE_CODE']=='PL') )[0]
    indAT=np.where( (metaAT['gdf']['OPENING_ID']==metaOPS['gdf']['OPENING_ID'][i]) )[0]
    if indAT.size==0:
        continue
    #if indAT.size<2:
    #    continue

    #list(metaFC['gdf'].keys())
    #indFC=np.where( (metaFC['gdf']['OPENING_ID']==metaOPS['gdf']['OPENING_ID'][i]) )[0]
    #metaFC['gdf']['SILV_POLYGON_NUMBER'][indFC]
    #metaFC['gdf']['SILV_POLYGON_AREA'][indFC]
    #metaFC['gdf']['STOCKING_STATUS_CODE'][indFC]
    #metaFC['gdf']['STOCKING_TYPE_CODE'][indFC]

    #metaAT['gdf']['SILV_FUND_SOURCE_CODE'][indAT]
    #metaAT['gdf']['SILV_TECHNIQUE_CODE'][indAT]
    Yr=metaAT['gdf']['Year'][indAT]
    A_PL=metaAT['gdf']['ACTUAL_TREATMENT_AREA'][indAT]
    N_PL=metaAT['gdf']['ACTUAL_PLANTED_NUMBER'][indAT]
    SPH=N_PL/A_PL
    #A_PL_Sum=np.sum(A_PL)

    A_OP_FromGDB=metaOP['gdf']['GEOMETRY_Area'][i]/10000
    A_OP=indO2.size
    ind_Art=np.where(zArt1[indO2]==1)[0]
    A_ART=ind_Art.size

    #if (A_ART>0) & (A_ART<A_PL_Sum):
    #    break
    fA=np.max(A_PL)/A_ART

    for iY in range(Yr.size):
        if (Yr[iY]<tv[0]) | (Yr[iY]>tv[-1]):
            continue
        if (A_ART==0) | (fA>1.1):
            d[Yr[iY]]['Idx'].append(indO2)
            d[Yr[iY]]['SPH'].append(SPH[iY]*np.ones(indO2.size))
        else:
            d[Yr[iY]]['Idx'].append(indO2[ind_Art])
            d[Yr[iY]]['SPH'].append(SPH[iY]*np.ones(ind_Art.size))

zM=zRef.copy()
zM['Data']=np.zeros(zRef['Data'].shape,dtype='int8')
for iT in range(tv.size):
    z=zRef.copy()
    z['Data']=np.zeros(zRef['Data'].shape,dtype='int16')
    for i in range(len(d[tv[iT]]['Idx'])):
        ind=d[tv[iT]]['Idx'][i]
        z['Data'][ indO1[0][ind],indO1[1][ind]  ]=d[tv[iT]]['SPH'][i]
    ind=np.where(z['Data']>0)
    zM['Data'][ind]=zM['Data'][ind]+1

    #z['Data']=z['Data'].astype('int16')
    #gis.SaveGeoTiff(z,r'C:\Users\rhember\Documents\Data\BC1ha\Results\Planting_FromRESULTS_' + str(tv[iT]) + '.tif')
gis.SaveGeoTiff(zM,r'C:\Users\rhember\Documents\Data\BC1ha\Results\Planting_FromRESULTS_MaskCount.tif')

plt.matshow(zM['Data'])

#%% Rasterize wildfire occurrence by year

# Input path to RESULTS database (downloaded from BC data catalogue)
pthin=r'C:\Users\rhember\Documents\Data\ForestInventory\Disturbances\20220422\Disturbances.gdb'
fiona.listlayers(pthin)
lnam='PROT_HISTORICAL_FIRE_POLYS_SP'

# Open the shapefile
df=gpd.read_file(pthin,layer=lnam)

# Remove features with no geometry
df=df[df.geometry!=None]

tv=np.arange(1920,2022,1)

for iT in range(tv.size):
    df0=df[df.FIRE_YEAR==tv[iT]].copy()
    shapes=((geom,value) for geom, value in zip(df0.geometry,df0.FIRE_YEAR))

    z0=np.zeros(zRef['Data'].shape,dtype=float)
    burned=features.rasterize(shapes=shapes,fill=0,out=z0,transform=zRef['Transform'])

    zOut=zRef.copy()
    zOut['Data']=np.zeros(zRef['Data'].shape,dtype='int8')
    zOut['Data'][burned>0]=1
    fout=r'C:\Users\rhember\Documents\Data\BC1ha\Disturbances\\' + lnam + '_' + str(tv[iT]) + '.tif'
    gis.SaveGeoTiff(zOut,fout)

# Test
#z=gis.OpenGeoTiff(fout)
#plt.matshow(z['Data'])

# Mask of where wildfire has occurred
z=zRef.copy()
z['Data']=np.zeros(zRef['Data'].shape,dtype='int8')
for iT in range(tv.size):
    d=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Disturbances\PROT_HISTORICAL_FIRE_POLYS_SP_' + str(tv[iT]) + '.tif')
    z['Data'][ (d['Data']>0) ]=1
plt.matshow(z['Data'])
gis.SaveGeoTiff(z,r'C:\Users\rhember\Documents\Data\BC1ha\Disturbances\PROT_HISTORICAL_FIRE_POLYS_SP_All.tif')

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

u_Severity_New=np.array(['T','L','M','S','V','G'])

#df0=pd.DataFrame(data=np.arange(1,u_Severity.size+1,1),columns=['ID'])
#df0['PEST_SEVERITY_CODE']=u_Severity
#df0.to_excel(pthout + '\\' + lnam + '_PEST_SEVERITY_CODE.xlsx')

df['ID_Severity']=np.zeros(len(df),dtype='float')
for i in range(len(u_Severity_New)):
    ind=np.where(df['PEST_SEVERITY_CODE']==u_Severity_New[i])[0]
    df.loc[ind,'ID_Severity']=i+1

tv=np.arange(1990,2022,1)

#fcd=['IDL','IBM','IBD','IBS','IDW','DFL']
#fcd=['IBM']
fcd=['IBS','IBB','IBD']

for iP in range(len(fcd)):
    for iT in range(tv.size):
        print(tv[iT])
        zOut=zRef.copy()
        zOut['Data']=np.zeros(zRef['Data'].shape,dtype='int8')

        df0=df[ (df[fnam]==fcd[iP]) & (df.CAPTURE_YEAR==tv[iT]) ].copy()
        df0=df0[df0.geometry!=None]

        if len(df0)>0:
            shapes=((geom,value) for geom, value in zip(df0.geometry,df0.ID_Severity))
            z0=np.zeros(zRef['Data'].shape,dtype=float)
            burned=features.rasterize(shapes=shapes,fill=0,out=z0,transform=zRef['Transform'])
            zOut['Data']=burned.astype('int8')

        fout=r'C:\Users\rhember\Documents\Data\BC1ha\Disturbances' + '\\' + lnam + '_' + fcd[iP] + '_SeverityClass_' + str(tv[iT]) + '.tif'
        gis.SaveGeoTiff(zOut,fout)

# Test
#z=gis.OpenGeoTiff(fout)
#plt.matshow(z['Data'])

#%% Rasterize aerial spray treatment

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

tv=np.arange(2017,2020,1)

for iT in range(tv.size):

    df0=df[df.FIRE_YEAR==tv[iT]].copy()
    shapes=((geom,value) for geom, value in zip(df0.geometry,df0.ID_BSR))

    z0=np.zeros(zRef['Data'].shape,dtype=float)
    burned=features.rasterize(shapes=shapes,fill=0,out=z0,transform=zRef.Transform)

    zOut=zRef.copy()
    zOut['Data']=5*np.ones(zRef['Data'].shape,dtype='int16')
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

z0=np.zeros(zRef['Data'].shape,dtype=float)
burned=features.rasterize(shapes=shapes,fill=0,out=z0,transform=zRef['Transform'])

zOut=zRef.copy()
zOut['Data']=np.zeros(zRef['Data'].shape,dtype='int16')
zOut['Data'][burned>0]=1
gis.SaveGeoTiff(zOut,r'C:\Users\rhember\Documents\Data\BC1ha\LandUseLandCover\bcts_op_area.tif')

#%% Prepare global forest change

bm=gpd.read_file(r'C:\Users\rhember\Documents\Data\Basemaps\Basemaps.gdb',layer='NRC_POLITICAL_BOUNDARIES_1M_SP')

fref=r'C:\Users\rhember\Documents\Data\BC1ha\Admin\BC_Land_Mask.tif'

finL=['50N_120W','50N_130W','60N_120W','60N_130W','60N_140W']

# Loss year
z=zRef.copy()
z['Data']=np.zeros(z['Data'].shape)
for fin in finL:
    fin=r'C:\Users\rhember\Documents\Data\Global Forest Change\Hansen_GFC-2021-v1.9_lossyear_' + fin + '.tif'
    fout=r'C:\Users\rhember\Documents\Data\Global Forest Change\Hansen_GFC-2021-v1.9_lossyear_' + fin + 'p.tif'
    gis.ReprojectRasterAndClipToRaster(fin,fout,fref,bm.crs)

    z0=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\Global Forest Change\Hansen_GFC-2021-v1.9_lossyear_' + fin + 'p.tif')
    ind=np.where(z0['Data']>0)
    z['Data'][ind]=z0['Data'][ind]

ind=np.where(zRef['Data']==0)
z['Data'][ind]=0

z['Data']=z['Data'].astype('int16')
gis.SaveGeoTiff(z,r'C:\Users\rhember\Documents\Data\BC1ha\Disturbances\GlobalForestChange_LossYear_2021.tif')
plt.matshow(z['Data'])

#ind=np.where(z['Data']>0)
#ind[0].size


# Gain
z=zRef.copy()
z['Data']=np.zeros(z['Data'].shape)
for fin in finL:
    fin0=r'C:\Users\rhember\Documents\Data\Global Forest Change\Hansen_GFC-2021-v1.9_gain_' + fin + '.tif'
    fout=r'C:\Users\rhember\Documents\Data\Global Forest Change\Hansen_GFC-2021-v1.9_gain_' + fin + 'p.tif'
    gis.ReprojectRasterAndClipToRaster(fin0,fout,fref,bm.crs)

    z0=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\Global Forest Change\Hansen_GFC-2021-v1.9_gain_' + fin + 'p.tif')
    ind=np.where(z0['Data']>0)
    z['Data'][ind]=z0['Data'][ind]

ind=np.where(zRef['Data']==0)
z['Data'][ind]=0

z['Data']=z['Data'].astype('int16')
gis.SaveGeoTiff(z,r'C:\Users\rhember\Documents\Data\BC1ha\Disturbances\GlobalForestChange_Gain_2021.tif')
plt.matshow(z['Data'])

#%% RASTERIZE VRI VARIABLES

# # *** Gave up on this- takes too long - do it in ArcGIS ***

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

zBGC=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\VRI\becz.tif')
zBGC['Data']=zBGC['Data'].flatten()
lutBGC=gu.ReadExcel(r'C:\Users\rhember\Documents\Data\BC1ha\VRI\becz_lut.xlsx')

zMAT=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Climate\BC1ha_mat_norm_1971to2000_si_hist_v1.tif')
zMAT['Data']=zMAT['Data'].flatten().astype(float)/10

zWS=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Climate\BC1ha_ws_gs_norm_1971to2000_comp_hist_v1.tif')
zWS['Data']=zWS['Data'].flatten().astype(float)

zSI=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\VRI\si.tif')
zSI['Data']=zSI['Data'].flatten().astype(float)

zA=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\VRI\proj_age_1.tif')
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


#%% Convert old BC1ha .mat files to new BC1ha geotiffs with standardized extent

def ConvertClimateNormals():

    #zRef=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Climate\Monthly\BC1ha_tmean_gs_norm_1971to2000_si_hist_v1_c.tif')
    zRef=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Admin\BC_Land_Mask.tif')

    for mo in range(12):
        print(mo+1)
        z=spio.loadmat(r'C:\Users\rhember\Documents\Data\BC1ha\Climate\Monthly\BC1ha_tmean_mon_norm_1971to2000_si_hist_v1\BC1ha_tmean_mon_norm_1971to2000_si_hist_v1_' + str(mo+1) + '.mat',squeeze_me=True)
        idat=np.where(np.asarray(z['z'].dtype.names)=='Data')[0][0]
        iSF=np.where(np.asarray(z['z'].dtype.names)=='ScaleFactor')[0][0]
        z0=np.flip(z['z'][()][idat].astype(float)*z['z'][()][iSF],axis=0)
        z1=zRef.copy()
        z1['Data']=z0
        z1=gis.ClipToRaster(z1,zRef)
        gis.SaveGeoTiff(z1,r'C:\Users\rhember\Documents\Data\BC1ha\Climate\Monthly\BC1ha_tmean_mon_norm_1971to2000_si_hist_v1\BC1ha_tmean_mon_norm_1971to2000_si_hist_v1_' + str(mo+1) + '.tif')

    for mo in range(12):
        print(mo+1)
        z=spio.loadmat(r'C:\Users\rhember\Documents\Data\BC1ha\Climate\Monthly\BC1ha_prcp_mon_norm_1971to2000_si_hist_v1\BC1ha_prcp_mon_norm_1971to2000_si_hist_v1_' + str(mo+1) + '.mat',squeeze_me=True)
        idat=np.where(np.asarray(z['z'].dtype.names)=='Data')[0][0]
        iSF=np.where(np.asarray(z['z'].dtype.names)=='ScaleFactor')[0][0]
        z0=np.flip(z['z'][()][idat].astype(float)*z['z'][()][iSF],axis=0)
        z1=zRef.copy()
        z1['Data']=z0
        z1=gis.ClipToRaster(z1,zRef)
        gis.SaveGeoTiff(z1,r'C:\Users\rhember\Documents\Data\BC1ha\Climate\Monthly\BC1ha_prcp_mon_norm_1971to2000_si_hist_v1\BC1ha_prcp_mon_norm_1971to2000_si_hist_v1_' + str(mo+1) + '.tif')

    for mo in range(12):
        print(mo+1)
        z=spio.loadmat(r'C:\Users\rhember\Documents\Data\BC1ha\Climate\Monthly\BC1ha_rswd_mon_norm_1971to2000_si_hist_v1\BC1ha_rswd_mon_norm_1971to2000_si_hist_v1_' + str(mo+1) + '.mat',squeeze_me=True)
        idat=np.where(np.asarray(z['z'].dtype.names)=='Data')[0][0]
        iSF=np.where(np.asarray(z['z'].dtype.names)=='ScaleFactor')[0][0]
        z0=np.flip(z['z'][()][idat].astype(float)*z['z'][()][iSF],axis=0)
        z1=zRef.copy()
        z1['Data']=z0
        z1=gis.ClipToRaster(z1,zRef)
        gis.SaveGeoTiff(z1,r'C:\Users\rhember\Documents\Data\BC1ha\Climate\Monthly\BC1ha_rswd_mon_norm_1971to2000_si_hist_v1\BC1ha_rswd_mon_norm_1971to2000_si_hist_v1_' + str(mo+1) + '.tif')

    for mo in range(12):
        print(mo+1)
        z=spio.loadmat(r'C:\Users\rhember\Documents\Data\BC1ha\Climate\Monthly\BC1ha_vpd_mon_norm_1971to2000_si_hist_v1\BC1ha_vpd_mon_norm_1971to2000_si_hist_v1_' + str(mo+1) + '.mat',squeeze_me=True)
        idat=np.where(np.asarray(z['z'].dtype.names)=='Data')[0][0]
        iSF=np.where(np.asarray(z['z'].dtype.names)=='ScaleFactor')[0][0]
        z0=np.flip(z['z'][()][idat].astype(float)*z['z'][()][iSF],axis=0)
        z1=zRef.copy()
        z1['Data']=z0
        z1=gis.ClipToRaster(z1,zRef)
        gis.SaveGeoTiff(z1,r'C:\Users\rhember\Documents\Data\BC1ha\Climate\Monthly\BC1ha_vpd_mon_norm_1971to2000_si_hist_v1\BC1ha_vpd_mon_norm_1971to2000_si_hist_v1_' + str(mo+1) + '.tif')

    for mo in range(12):
        print(mo+1)
        z=spio.loadmat(r'C:\Users\rhember\Documents\Data\BC1ha\Climate\Monthly\BC1ha_etp_tmw_norm_1971to2000_comp_hist_v1\BC1ha_etp_mon_norm_1971to2000_comp_hist_v1_' + str(mo+1) + '.mat',squeeze_me=True)
        idat=np.where(np.asarray(z['z'].dtype.names)=='Data')[0][0]
        iSF=np.where(np.asarray(z['z'].dtype.names)=='ScaleFactor')[0][0]
        z0=np.flip(z['z'][()][idat].astype(float)*z['z'][()][iSF],axis=0)
        z0=z0*100
        z1=zRef.copy()
        z1['Data']=z0.astype('int16')
        z1=gis.ClipToRaster(z1,zRef)
        gis.SaveGeoTiff(z1,r'C:\Users\rhember\Documents\Data\BC1ha\Climate\Monthly\BC1ha_etp_tmw_norm_1971to2000_comp_hist_v1\BC1ha_etp_mon_norm_1971to2000_comp_hist_v1_' + str(mo+1) + '.tif')

    return