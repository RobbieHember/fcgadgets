
'''============================================================================

BC1ha GRID

Notes:
    - These rasters were produced from geodatabases in ArcGIS. The polygon area 
    was used as a priority field. (see BC1ha_RasterizeWithArcpyToolbox.py)
    - Fire, cutblocks and forest health were done in Python (this script) because
    ov overlapping polygons

============================================================================'''



#%% IMPORT MODULES

import sys
import numpy as np
import gdal
import osr 
import matplotlib.pyplot as plt
import geopandas as gpd
import pandas as pd
import fiona
import rasterio
from rasterio import features
#import cv2
from shapely.geometry import Point, Polygon
from shapely import geometry
import dbfread

import fcgadgets.utilities.utilities_general as gu
import fcgadgets.utilities.utilities_gis as gis


#%% DEFINE STANDARD GRID 
# Select an existing file that will define the SRS and extent of the BC1ha grid.
# Currently, it is a tif TSA file that was rasterized in ArcGIS.

finS=r'Z:\!Workgrp\Forest Carbon\Data\BC1ha\Admin\tsa.tif'
zS=gis.OpenGeoTiff(finS)


#%% TIMBER SUPPLY AREA

#------------------------------------------------------------------------------
# Create TSA Key
# This creates a key between the raster values for each TSA, and the TSA names
# from a shapefile of TSA downloaded from the BC Data Cat.
#------------------------------------------------------------------------------

# Import the shapefile with the TSA names and get the unique list of TSA Numbers 
# and each corresponding TSA name.
gdf_tsa=gpd.read_file(r'Z:\!Workgrp\Forest Carbon\Data\TSA\tsa.shp')
u=gdf_tsa['TSA_NUMBER'].unique()
nam=[]
for i in range(len(u)):
    a=gdf_tsa[gdf_tsa.TSA_NUMBER==u[i]].TSA_NUMB_1.unique()
    nam.append(a[0])

# Key between raster values and TSA name (this fixes the original file created by ESRI the first time)

# Import raster metadata
db=dbfread.DBF(r'Z:\!Workgrp\Forest Carbon\Data\BC1ha\Admin\tsa.tif.vat.dbf',load=True)  

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
df.to_excel(r'Z:\!Workgrp\Forest Carbon\Data\BC1ha\Admin\tsa_key.xlsx',index=False)

#------------------------------------------------------------------------------
# Digitize the boundary of TSAs
#------------------------------------------------------------------------------

# Open dataframe containing TSA names and values for raster grid of TSAs
zTSA=gis.OpenGeoTiff(r'Z:\!Workgrp\Forest Carbon\Data\BC1ha\Admin\tsa.tif')
#ds=gdal.Open(r'Z:\!Workgrp\Forest Carbon\Data\BC1ha\Admin\tsa.tif');
#tsa=ds.ReadAsArray(); 
#gt=ds.GetGeoTransform(); del ds; 
#m,n=tsa.shape; 
#extent=(gt[0],gt[0]+n*gt[1],gt[3]-m*gt[1],gt[3])
#xG=np.arange(gt[0],gt[0]+gt[1]*n,gt[1])
#yG=np.flip(np.arange(gt[3]-m*gt[1],gt[3],gt[1]))

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


#%% REVISE EXTENT TO BE CONSISTENT WITH STANDARD RASTER

# Admin boundaries

fin1=r'Z:\!Workgrp\Forest Carbon\Data\BC1ha\AdminBoundaries\Processed\tsa.tif'
fin1_adj=r'Z:\!Workgrp\Forest Carbon\Data\BC1ha\AdminBoundaries\Processed\tsa.tif'
gis.ReviseRasterExtent(fin1,finS,fin1_adj)

# Soils

fin1=r'Z:\!Workgrp\Forest Carbon\Data\BC1ha\Soil\Given\soil_symb1.tif'
fin1_adj=r'Z:\!Workgrp\Forest Carbon\Data\BC1ha\Soil\Processed\soil_symb1_rev.tif'
gis.ReviseRasterExtent(fin1,finS,fin1_adj)

fin1=r'Z:\!Workgrp\Forest Carbon\Data\BC1ha\Soil\Given\soil_pct1.tif'
fin1_adj=r'Z:\!Workgrp\Forest Carbon\Data\BC1ha\Soil\Processed\soil_pct1_rev.tif'
gis.ReviseRasterExtent(fin1,finS,fin1_adj)

fin1=r'Z:\!Workgrp\Forest Carbon\Data\BC1ha\Soil\Given\soil_symb2.tif'
fin1_adj=r'Z:\!Workgrp\Forest Carbon\Data\BC1ha\Soil\Processed\soil_symb2_rev.tif'
gis.ReviseRasterExtent(fin1,finS,fin1_adj)

fin1=r'Z:\!Workgrp\Forest Carbon\Data\BC1ha\Soil\Given\soil_pct2.tif'
fin1_adj=r'Z:\!Workgrp\Forest Carbon\Data\BC1ha\Soil\Processed\soil_pct2_rev.tif'
gis.ReviseRasterExtent(fin1,finS,fin1_adj)

# Land use and land cover

fin1=r'Z:\!Workgrp\Forest Carbon\Data\LandUseLandCover\landuse.ecophysio.tif'
fin1_adj=r'Z:\!Workgrp\Forest Carbon\Data\BC1ha\LandUseLandCover\Processed\landuse.ecophysio_rev.tif'
gis.ExtentRevise.ReviseRasterExtent(fin1,finS,fin1_adj)

fin1=r'Z:\!Workgrp\Forest Carbon\Data\LandUseLandCover\landuse.btm.tif'
fin1_adj=r'Z:\!Workgrp\Forest Carbon\Data\BC1ha\LandUseLandCover\Processed\landuse.btm_rev.tif'
gis.ExtentRevise.ReviseRasterExtent(fin1,finS,fin1_adj)

# Grasslands GCC

fin1=r'Z:\!Workgrp\Forest Carbon\Data\BC1ha\LandUseLandCover\Given\grasslands.grasslands.tif'
fin1_adj=r'Z:\!Workgrp\Forest Carbon\Data\BC1ha\LandUseLandCover\Processed\grasslands.grasslands_rev.tif'
gis.ExtentRevise.ReviseRasterExtent(fin1,finS,fin1_adj)

# Disturbance: Cutblocks

# Consistent with VRI to start with

ds=gdal.Open(fin1_adj)
data=ds.ReadAsArray()

plt.figure(4)
plt.matshow(data[0::5,0::5],clim=(0,150))
plt.colorbar()


#%% Rasterize wildfire occurrence by year

# Input path to RESULTS database (downloaded from BC data catalogue)
pthin=r'C:\Users\rhember\Documents\Data\ForestInventory\Disturbances\20200430\Disturbances.gdb'
fiona.listlayers(pthin)
lnam='PROT_HISTORICAL_FIRE_POLYS_SP'

# Open the shapefile
df=gpd.read_file(pthin,layer=lnam)

# Remove features with no geometry
df=df[df.geometry!=None]

zTSA=gis.OpenGeoTiff(r'Z:\!Workgrp\Forest Carbon\Data\BC1ha\Admin\tsa.tif')

tv=np.arange(1920,2020,1)

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
pthin=r'C:\Users\rhember\Documents\Data\ForestInventory\Disturbances\20200430\Disturbances.gdb'
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
u_Severity=np.unique(df.PEST_SEVERITY_CODE)
#u_Severity=np.array(['T','L','M','S','V','G'])
#df0=pd.DataFrame(data=np.arange(1,u_Severity.size+1,1),columns=['ID'])
#df0['PEST_SEVERITY_CODE']=u_Severity
#df0.to_excel(pthout + '\\' + lnam + '_PEST_SEVERITY_CODE.xlsx')

df['PSC']=np.zeros(len(df),dtype='float')
for i in range(len(u_Severity)):
    ind=np.where(df['PEST_SEVERITY_CODE']==u_Severity[i])[0]
    df.loc[ind,'PSC']=i+1

zTSA=gis.OpenGeoTiff(r'Z:\!Workgrp\Forest Carbon\Data\BC1ha\Admin\tsa.tif')

tv=np.arange(1950,2020,1)

fcd=['IDL','IBM','IBD','IBS','IDW','DFL']
for iP in range(0,1):#len(fcd)):
    for iT in range(tv.size):
    
        zOut=zTSA.copy()
        zOut['Data']=np.zeros(zTSA['Data'].shape,dtype='int16')
    
        df0=df[ (df[fnam]==fcd[iP]) & (df.CAPTURE_YEAR==tv[iT]) ].copy() 
        df0=df0[df0.geometry!=None]
        
        if len(df0)>0:
            shapes=((geom,value) for geom, value in zip(df0.geometry,df0.PSC))    
            z0=np.zeros(zTSA['Data'].shape,dtype=float)
            burned=features.rasterize(shapes=shapes,fill=0,out=z0,transform=zTSA['Transform'])
            zOut['Data']=burned.astype('int16')
            
        fout=r'Z:\!Workgrp\Forest Carbon\Data\BC1ha\Disturbances' + '\\' + lnam + '_' + fcd[iP] + '_SeverityClass_' + str(tv[iT]) + '.tif'
        gis.SaveGeoTiff(zOut,fout)

# Test
#z=gis.OpenGeoTiff(r'Z:\!Workgrp\Forest Carbon\Data\BC1ha\Disturbances\test1.tif')
#plt.matshow(z['Data'])



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





#%% RASTERIZE VRI VARIABLES
# *** Gave up on this- takes too long ***

# Open template raster file and copy it's metadata
rst=rasterio.open(r'Z:\!Workgrp\Forest Carbon\Data\BC1ha\VRI\age1.tif')
meta=rst.meta.copy()
meta.update(compress='lzw')

# Initialize raster
z=np.zeros(rst.shape,dtype=float)

# Input path to RESULTS database (downloaded from BC data catalogue)
pthin=r'Z:\!Workgrp\Forest Carbon\Data\VRI\20190501\VRI.gdb'

# Look at the layers in geodatabase
fiona.listlayers(pthin)

lnam='VEG_COMP_LYR_R1_POLY'

# Output path
pthout=r'Z:\!Workgrp\Forest Carbon\Data\BC1ha\VRI'

# Initialize new raster file
out=rasterio.open(pthout + '\\test.tif', 'w', **meta)

# Get unique values
s=['0']*10000
cnt=0
with fiona.open(pthin,layer=lnam) as source:
    for feat in source:
        tmp=feat['properties']['BCLCS_LEVEL_2']
        if tmp!=None:
            s[cnt]=tmp
        cnt=cnt+1
us=np.unique(s)
Code=['L','W']
ID=[1,2]

df=gpd.GeoDataFrame()
with fiona.open(pthin,layer=lnam) as source:
    for feat in source:
        df0=gpd.GeoDataFrame.from_features([feat])
        df1=gpd.GeoDataFrame()
        df1.geometry=df0.geometry
        df1['BCLCS_LEVEL_2']=df0['BCLCS_LEVEL_2']
        df=pd.concat([df,df1])
        
        iCode=int(np.where(Code==df.BCLCS_LEVEL_2.values)[0])
        df.BCLCS_LEVEL_2=ID[iCode]
        shapes=((geom,value) for geom, value in zip(df.geometry,df.BCLCS_LEVEL_2))        
        burned=features.rasterize(shapes=shapes,fill=0,out=z,transform=out.transform)
        ind=np.where((burned[:]!=0))
        z[ind]=ID[iCode]
        break
        
        coords0=feat['geometry']['coordinates']        
        for i in range(len(coords0)):            
            coords1=coords0[i]
            for j in range(len(coords1)):
                coords2=coords1[j]
                
                poly=np.asarray(coords2)

# Output path
pthout=r'Z:\!Workgrp\Forest Carbon\Data\BC1ha\VRI'

# Initialize new raster file
out=rasterio.open(pthout + '\\test.tif', 'w', **meta)

# Set what variable will be done
flg=1

z1=np.zeros(rst.shape,dtype=float)
        
shapes=((geom,value) for geom, value in zip(df.geometry,df.FIRE_YEAR))
    
z0=np.zeros(rst.shape,dtype=float)
burned=features.rasterize(shapes=shapes,fill=0,out=z0,transform=out.transform)
    
#ind=np.where((z1[:]==0) & (burned[:]!=0))
#z1[ind]=tv[i]

out.close()

out=rasterio.open(pthout + '\\' + lnam + '_year1.tif', 'w', **meta)
out.write_band(1,z1.astype(np.int16))
out.close()




#%% Import climate from old BC1ha project in matlab

for mo in range(12,13):
    z=gis.OpenGeoTiff(r'D:\BC1ha_rswd_mon_norm_1971to2000_si_hist_v1_' + str(mo) + '.tif')
    z=gis.ClipRaster(z,[zS['minx'],zS['maxx']],[zS['miny'],zS['maxy']])
    gis.SaveGeoTiff(z,r'C:\Users\rhember\Documents\Data\BC1ha_rswd_mon_norm_1971to2000_si_hist_v1_' + str(mo) + '.tif')
    del z
    gc.collect()

plt.matshow(z.Data)




