
#%% Import modules

import os
import numpy as np
import time
import matplotlib.pyplot as plt
import geopandas as gpd
import pandas as pd
import pyproj
import scipy.io as spio
import fiona
import rasterio
from rasterio import features
from shapely import geometry
from scipy.interpolate import griddata
import cv2
#from rasterio.transform import from_origin
#from shapely.geometry import Point,Polygon
import fcgadgets.macgyver.util_general as gu
import fcgadgets.macgyver.util_gis as gis
import fcgadgets.macgyver.util_query_gdb as qgdb
import fcgadgets.bc1ha.bc1ha_util as u1ha

#%% Import paths and look-up-tables

meta=u1ha.Init()
gp=gu.SetGraphics('Manuscript')

# Build look up tables (only do this once a year, takes 8 hours)
flg=0
if flg==1:
    u1ha.BuildLUTsFromSourceDBs(meta)

#%% Define reference grid

zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])

#%% Create land mask for BC

df=gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='NRC_POLITICAL_BOUNDARIES_1M_SP')
df=df[df.geometry!=None]
df['ID']=1
shapes=((geom,value) for geom, value in zip(df.geometry,df.ID))
z=np.zeros(zRef['Data'].shape,dtype=float)
burned=features.rasterize(shapes=shapes,fill=0,out=z,transform=zRef['Transform'])
# plt.close('all'); plt.matshow(burned) # Confirm that it worked
zOut=zRef.copy()
zOut['Data']=np.zeros(zRef['Data'].shape,dtype='int8')
zOut['Data'][burned>0]=1
zOut['Data']=zOut['Data'].astype('int8')
gis.SaveGeoTiff(zOut,meta['Paths']['bc1ha'] + '\\LandCoverUse\LandMask.tif')

#%% Digitize the boundary of TSAs (the original is organized at the sub-TSA level)

def DigitizeTSABoundaries():
    zTSA=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\FADM_TSA\\TSA_NUMBER_DESCRIPTION.tif')
    u=np.unique(zTSA['Data'])
    u=u[u>0]
    gdf=gpd.GeoDataFrame(data=[],columns=['Value','Name','geometry'])
    cnt=0
    for iU in range(u.size):
        ind=np.where(zTSA['Data']==u[iU])
        #z=np.zeros((zTSA['m'],zTSA['n'],3),dtype=np.uint8)
        z=np.zeros(zTSA['Data'].shape,dtype=np.uint8)
        z[ind]=255
        z=np.dstack([z]*3)
        z=cv2.cvtColor(z,cv2.COLOR_BGR2GRAY) # Convert to grey scale
        cont=cv2.findContours(image=z,mode=cv2.RETR_LIST,method=cv2.CHAIN_APPROX_SIMPLE)
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
                x=int(zTSA['X'][0,c])
                y=int(zTSA['Y'][r,0])
                pointList.append(geometry.Point(x,y))
            gdf.loc[cnt,'Value']=int(u[iU])
            gdf.loc[cnt,'Name']=u1ha.lut_n2s(meta['LUT']['FADM_TSA']['TSA_NUMBER_DESCRIPTION'],u[iU])[0]
            gdf.loc[cnt,'geometry']=geometry.Polygon([[p.x,p.y] for p in pointList])
            cnt=cnt+1
    gdf.to_file(r'C:\Users\rhember\Documents\Data\Geodatabases\LandUse\tsa.geojson',driver='GeoJSON')

    return

#%% Rasterize variables from source
# *** VRI and Forest Cover Inventory area too big - crash ***
# prp=u1ha.GetVariablesFromGDB(meta,'BC_MAJOR_WATERSHEDS')

u1ha.RasterizeFromSource(meta,zRef,'BC_MAJOR_WATERSHEDS','MAJOR_WATERSHED_CODE')
u1ha.RasterizeFromSource(meta,zRef,'BEC_BIOGEOCLIMATIC_POLY','ZONE')
u1ha.RasterizeFromSource(meta,zRef,'BEC_BIOGEOCLIMATIC_POLY','SUBZONE')
u1ha.RasterizeFromSource(meta,zRef,'BEC_NATURAL_DISTURBANCE_SV','NATURAL_DISTURBANCE_TYPE_CODE')
u1ha.RasterizeFromSource(meta,zRef,'FTEN_CUT_BLOCK_POLY_SVW','TIMBER_MARK')
u1ha.RasterizeFromSource(meta,zRef,'FTEN_CUT_BLOCK_POLY_SVW','CUT_BLOCK_ID')
u1ha.RasterizeFromSource(meta,zRef,'FWA_WETLANDS_POLY','WATERBODY_TYPE')
u1ha.RasterizeFromSource(meta,zRef,'OGSR_TAP_PRIORITY_DEF_AREA_SP','OGSR_TPDA_SYSID')
u1ha.RasterizeFromSource(meta,zRef,'OGSR_TAP_PRIORITY_DEF_AREA_SP','PRIORITY_DEFERRAL_ID')
u1ha.RasterizeFromSource(meta,zRef,'RMP_OGMA_LEGAL_ALL_SVW','LEGAL_OGMA_PROVID')
u1ha.RasterizeFromSource(meta,zRef,'RSLT_OPENING_SVW','OPENING_ID')
u1ha.RasterizeFromSource(meta,zRef,'RSLT_ACTIVITY_TREATMENT_SVW','OPENING_ID')
u1ha.RasterizeFromSource(meta,zRef,'RSLT_FOREST_COVER_RESERVE_SVW','SILV_RESERVE_CODE')
u1ha.RasterizeFromSource(meta,zRef,'RSLT_FOREST_COVER_INV_SVW','SILV_RESERVE_CODE')
u1ha.RasterizeFromSource(meta,zRef,'RSLT_FOREST_COVER_INV_SVW','SILV_RESERVE_OBJECTIVE_CODE')
u1ha.RasterizeFromSource(meta,zRef,'RSLT_FOREST_COVER_INV_SVW','STOCKING_TYPE_CODE') # takes 100 min
u1ha.RasterizeFromSource(meta,zRef,'RSLT_FOREST_COVER_INV_SVW','TREE_COVER_PATTERN_CODE')
u1ha.RasterizeFromSource(meta,zRef,'RSLT_FOREST_COVER_INV_SVW','I_TOTAL_STEMS_PER_HA')
u1ha.RasterizeFromSource(meta,zRef,'RSLT_FOREST_COVER_INV_SVW','I_TOTAL_WELL_SPACED_STEMS_HA')
u1ha.RasterizeFromSource(meta,zRef,'RSLT_FOREST_COVER_INV_SVW','I_CROWN_CLOSURE_PERCENT')
u1ha.RasterizeFromSource(meta,zRef,'RSLT_FOREST_COVER_INV_SVW','REFERENCE_YEAR')
u1ha.RasterizeFromSource(meta,zRef,'TA_REGIONAL_DISTRICTS_SVW','REGIONAL_DISTRICT_NAME')
u1ha.RasterizeFromSource(meta,zRef,'TA_PROTECTED_LANDS_SV','PROTECTED_LANDS_DESIGNATION')
u1ha.RasterizeFromSource(meta,zRef,'VEG_BURN_SEVERITY_SP','FIRE_YEAR')
u1ha.RasterizeFromSource(meta,zRef,'VEG_BURN_SEVERITY_SP','BURN_SEVERITY_RATING')

#%%

#z=u1ha.Import_Raster(meta,[],['bsr_yr','bsr_sc'])
#plt.hist(z['bsr_yr']['Data'][z['bsr_yr']['Data']>0].flatten()[0::20])

with fiona.open(meta['Paths']['GDB']['LandUse'],layer='FTEN_SPEC_USE_PERMIT_POLY_SVW') as source:
    for feat in source:
        prop=dict(feat['properties'].items())
        list(prop)
        #print(prop['LICENCE_TO_CUT_CODE'])
        break

#%% Gap-fill BGC Zone

zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
z=u1ha.Import_Raster(meta,[],['lcc1_c','bgcz'])

ivl=5
iGap=np.where( (zRef['Data']==1) & (z['bgcz']['Data']==0) )
iCal=np.where( (zRef['Data'][0::ivl,0::ivl]==1) & (z['bgcz']['Data'][0::ivl,0::ivl]>0) )
xy=np.column_stack([zRef['X'][0::ivl,0::ivl][iCal],zRef['Y'][0::ivl,0::ivl][iCal]])
vals=z['bgcz']['Data'][0::ivl,0::ivl][iCal]
zFill=griddata(xy,vals,(zRef['X'][iGap],zRef['Y'][iGap]),method='nearest')

z1=zRef.copy()
z1['Data']=np.zeros(zRef['Data'].shape,dtype='int16')
z1['Data']=z['bgcz']['Data']
z1['Data'][iGap]=zFill
plt.close('all'); plt.matshow(z1['Data'])
gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\BEC_BIOGEOCLIMATIC_POLY\\ZONE_GapFilled.tif')

#%% Rasterize wildfire

lNam='PROT_HISTORICAL_FIRE_POLYS_SP'
vNam='FIRE_YEAR'

if os.path.exists(meta['Paths']['bc1ha'] + '\\' + lNam)==False:
    os.mkdir(meta['Paths']['bc1ha'] + '\\' + lNam)

ind=np.where( (meta['Geos']['Variable Info']['Layer Name']==lNam) & (meta['Geos']['Variable Info']['Variable Name']==vNam) )[0]

pthin=meta['Paths']['GDB'][ meta['Geos']['Variable Info']['Geodatabase Name'][ind[0]] ]

df=gpd.read_file(pthin,layer=lNam)
df=df[df.geometry!=None]
df=df.reset_index()

zYearLast=zRef.copy()
zYearLast['Data']=np.zeros(zRef['Data'].shape,dtype='int16')

uYear=df[vNam].unique()
tv=np.arange(np.min(uYear),np.max(uYear),1)

for iT in range(tv.size):

    df0=df[df[vNam]==tv[iT]].copy()
    shapes=((geom,value) for geom, value in zip(df0.geometry,df0[vNam]))

    z0=np.zeros(zRef['Data'].shape,dtype=float)
    if len(df0)>0:
        burned=features.rasterize(shapes=shapes,fill=0,out=z0,transform=zRef['Transform'])

    z1=zRef.copy()
    z1['Data']=z0.astype(meta['Geos']['Variable Info']['Precision'][ind[0]])
    gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\' + lNam + '\\' + vNam + '_' + str(tv[iT]) + '.tif')

    # Update by year grid
    zYearLast['Data'][burned>0]=tv[iT]

# Year of last occurrence
z1=zRef.copy()
z1['Data']=zYearLast['Data'].astype('int16')
gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\' + lNam + '\\' + vNam + '_YearLast.tif')

# Mask of occurrence
z1=zRef.copy()
z1['Data']=np.zeros(zRef['Data'].shape,dtype='int8')
ind0=np.where(zYearLast['Data']>0)
z1['Data'][ind0]=1
gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\' + lNam + '\\' + vNam + '_MaskAll.tif')

# Pack into smaller number of layers

tv=np.arange(1920,2025,1)

# Initialize rasters
N_Year=6
z={'Year':{}}
for iY in range(N_Year):
    z['Year'][iY+1]=np.zeros(zRef['Data'].shape,dtype='int16')

for iT in range(tv.size):
    print(tv[iT])
    try:
        z0=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\' + lNam + '\\' + vNam + '_' + str(tv[iT]) + '.tif')['Data']
    except:
        continue
    ind=np.where((z['Year'][1]==0) & (z0!=0))
    z['Year'][1][ind]=tv[iT];
    ind=np.where((z['Year'][1]!=0) & (z['Year'][1]!=tv[iT]) & (z['Year'][2]==0) & (z0!=0))
    z['Year'][2][ind]=tv[iT];
    ind=np.where((z['Year'][2]!=0) & (z['Year'][2]!=tv[iT]) & (z['Year'][3]==0) & (z0!=0))
    z['Year'][3][ind]=tv[iT];
    ind=np.where((z['Year'][3]!=0) & (z['Year'][3]!=tv[iT]) & (z['Year'][4]==0) & (z0!=0))
    z['Year'][4][ind]=tv[iT];
    ind=np.where((z['Year'][4]!=0) & (z['Year'][4]!=tv[iT]) & (z['Year'][5]==0) & (z0!=0))
    z['Year'][5][ind]=tv[iT];
    ind=np.where((z['Year'][5]!=0) & (z['Year'][5]!=tv[iT]) & (z0!=0))
    z['Year'][6][ind]=tv[iT];

for iY in range(N_Year):
    z1=zRef.copy()
    z1['Data']=z['Year'][iY+1].astype('int16')
    gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\' + lNam + '\\' + vNam + '_' + str(iY+1) + '_Year.tif')

# Plot time series to confirm it worked
lNam='PROT_HISTORICAL_FIRE_POLYS_SP'
vNam='FIRE_YEAR'
tv,N=u1ha.TimeSeriesOccurrenceFromPackedEvents(meta,lNam,vNam,6)
plt.plot(tv,N,'-bo')

#%% Rasterize AOS occurrence by year

lNam='PEST_INFESTATION_POLY'
vNam='PEST_SEVERITY_CODE'

if os.path.exists(meta['Paths']['bc1ha'] + '\\' + lNam)==False:
    os.mkdir(meta['Paths']['bc1ha'] + '\\' + lNam)

indI=np.where( (meta['Geos']['Variable Info']['Layer Name']==lNam) & (meta['Geos']['Variable Info']['Variable Name']==vNam) )[0]

pthin=meta['Paths']['GDB'][ meta['Geos']['Variable Info']['Geodatabase Name'][indI[0]] ]

df=gpd.read_file(pthin,layer=lNam)
df=df[df.geometry!=None]
df=df.reset_index()

df=u1ha.CreateIdForCategoricalVariable(meta,lNam,vNam,df)

pestL=['IBM','IBS','IBB','IBD','IDW','IDL']

tv=np.arange(1951,2023,1)

for pest in pestL:

    zYearLast=zRef.copy()
    zYearLast['Data']=np.zeros(zRef['Data'].shape,dtype='int16')

    for iT in range(tv.size):

        df0=df[ (df['PEST_SPECIES_CODE']==pest) & (df['CAPTURE_YEAR']==tv[iT]) ].copy()
        df0=df0[df0.geometry!=None]
        df0=df0.reset_index()

        if len(df0)>0:
            shapes=((geom,value) for geom, value in zip(df0.geometry,df0['ID_PEST_SEVERITY_CODE']))
            z0=np.zeros(zRef['Data'].shape,dtype=float)
            burned=features.rasterize(shapes=shapes,fill=0,out=z0,transform=zRef['Transform'])
        else:
            z0=np.zeros(zRef['Data'].shape,dtype=float)

        z1=zRef.copy()
        z1['Data']=z0.astype(meta['Geos']['Variable Info']['Precision'][indI[0]])
        gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\' + lNam + '\\' + vNam + '_' + pest + '_' + str(tv[iT]) + '.tif')

        # Update by year grid
        zYearLast['Data'][z0>0]=tv[iT]

    # Year of last occurrence
    z1=zRef.copy()
    z1['Data']=zYearLast['Data'].astype('int16')
    gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\' + lNam + '\\' + vNam + '_' + pest + '_YearLast.tif')

    # Mask of occurrence
    z1=zRef.copy()
    z1['Data']=np.zeros(zRef['Data'].shape,dtype='int8')
    ind0=np.where(zYearLast['Data']>0)
    z1['Data'][ind0]=1
    gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\' + lNam + '\\' + vNam + '_' + pest + '_MaskAll.tif')

    # Pack into smaller number of layers

    # Initialize rasters
    N_Year=10
    z={'Year':{},'Severity':{}}
    for iY in range(N_Year):
        z['Year'][iY+1]=np.zeros(zRef['Data'].shape,dtype='int16')
        z['Severity'][iY+1]=np.zeros(zRef['Data'].shape,dtype='int16')

    for iT in range(tv.size):
        print(tv[iT])
        z0=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\' + lNam + '\\' + vNam + '_' + pest + '_' + str(tv[iT]) + '.tif')['Data']
        ind=np.where((z['Year'][1]==0) & (z0!=0))
        z['Year'][1][ind]=tv[iT]
        z['Severity'][1][ind]=z0[ind]
        ind=np.where((z['Year'][1]!=0) & (z['Year'][1]!=tv[iT]) & (z['Year'][2]==0) & (z0!=0))
        z['Year'][2][ind]=tv[iT]
        z['Severity'][2][ind]=z0[ind]
        ind=np.where((z['Year'][2]!=0) & (z['Year'][2]!=tv[iT]) & (z['Year'][3]==0) & (z0!=0))
        z['Year'][3][ind]=tv[iT]
        z['Severity'][3][ind]=z0[ind]
        ind=np.where((z['Year'][3]!=0) & (z['Year'][3]!=tv[iT]) & (z['Year'][4]==0) & (z0!=0))
        z['Year'][4][ind]=tv[iT]
        z['Severity'][4][ind]=z0[ind]
        ind=np.where((z['Year'][4]!=0) & (z['Year'][4]!=tv[iT]) & (z['Year'][5]==0) & (z0!=0))
        z['Year'][5][ind]=tv[iT]
        z['Severity'][5][ind]=z0[ind]
        ind=np.where((z['Year'][5]!=0) & (z['Year'][5]!=tv[iT]) & (z['Year'][6]==0) & (z0!=0))
        z['Year'][6][ind]=tv[iT]
        z['Severity'][6][ind]=z0[ind]
        ind=np.where((z['Year'][6]!=0) & (z['Year'][6]!=tv[iT]) & (z['Year'][7]==0) & (z0!=0))
        z['Year'][7][ind]=tv[iT]
        z['Severity'][7][ind]=z0[ind]
        ind=np.where((z['Year'][7]!=0) & (z['Year'][7]!=tv[iT]) & (z['Year'][8]==0) & (z0!=0))
        z['Year'][8][ind]=tv[iT]
        z['Severity'][8][ind]=z0[ind]
        ind=np.where((z['Year'][8]!=0) & (z['Year'][8]!=tv[iT]) & (z['Year'][9]==0) & (z0!=0))
        z['Year'][9][ind]=tv[iT]
        z['Severity'][9][ind]=z0[ind];
        ind=np.where((z['Year'][9]!=0) & (z['Year'][9]!=tv[iT]) & (z0!=0))
        z['Year'][10][ind]=tv[iT]
        z['Severity'][10][ind]=z0[ind]

    for iY in range(N_Year):
        z1=zRef.copy()
        z1['Data']=z['Year'][iY+1].astype('int16')
        gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\' + lNam + '\\' + vNam + '_' + pest + '_' + str(iY+1) + '_Year.tif')
        z1=zRef.copy()
        z1['Data']=z['Severity'][iY+1].astype('int16')
        gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\' + lNam + '\\' + vNam + '_' + pest + '_' + str(iY+1) + '_Severity.tif')

# Plot time series to confirm it worked
lNam='PEST_INFESTATION_POLY'
vNam='PEST_SEVERITY_CODE_IBM'
nPack=10
tv,N=u1ha.TimeSeriesOccurrenceFromPackedEvents(meta,lNam,vNam,nPack)
plt.plot(tv,N,'-bo')

#%% Define insect outbreaks

# pestL=['IBM','IBS','IBB','IBD','IDW','IDL']

# tv=np.arange(1951,2023,1)

# for pest in pestL:

#     zMask=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\PEST_INFESTATION_POLY\\PEST_SEVERITY_CODE_' + pest + '_MaskAll.tif')
#     iMask=np.where(zMask['Data']==1)

#     zO=[{}]*3
#     for iO in range(3):
#         zO[iO]['Year Start']=np.zeros(zMask['Data'].shape,dtype='int16')
#         zO[iO]['Year End']=np.zeros(zMask['Data'].shape,dtype='int16')
#         zO[iO]['Duration']=np.zeros(zMask['Data'].shape,dtype='int16')
#         zO[iO]['Max Severity']=np.zeros(zMask['Data'].shape,dtype='int16')

#     OutbreakNum=np.zeros(zMask['Data'].shape,dtype='int16')
#     Active=np.zeros(zMask['Data'].shape,dtype='int16')
#     Duration=np.zeros(zMask['Data'].shape,dtype='int16')
#     for iT in range(tv.size):
#         print(tv[iT])
#         zS=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\PEST_INFESTATION_POLY\\PEST_SEVERITY_CODE_' + pest + '_' + str(tv[iT]) + '.tif')
#         for iO in range(3):
#             # Onset of outbreak
#             ind=np.where( (OutbreakNum==iO) & (Active==0) & (zS['Data']>0) & (zS['Data']!=meta['LUT']['PEST_INFESTATION_POLY']['PEST_SEVERITY_CODE']['T']) )
#             Active[ind]=1
#             zO[iO]['Year Start'][ind]=tv[iT]
#             zO[iO]['Max Severity'][ind]=np.maximum(zO[iO]['Max Severity'][ind],zS['Data'][ind])
#             zO[iO]['Duration'][ind]=zO[iO]['Duration'][ind]+1

#             # Continuation of outbreak
#             ind=np.where( (OutbreakNum==iO) & (Active==1) & (zS['Data']>0) & (zS['Data']!=meta['LUT']['PEST_INFESTATION_POLY']['PEST_SEVERITY_CODE']['T']) )
#             zO[iO]['Max Severity'][ind]=np.maximum(zO[iO]['Max Severity'][ind],zS['Data'][ind])
#             zO[iO]['Duration'][ind]=zO[iO]['Duration'][ind]+1

#             # End of outbreak
#             ind=np.where( (OutbreakNum==iO) & (Active==1) & (zS['Data']==0) | (Active==1) & (zS['Data']!=meta['LUT']['PEST_INFESTATION_POLY']['PEST_SEVERITY_CODE']['T']) )
#             zO[iO]['Year End'][ind]=tv[iT]-1
#             OutbreakNum[ind]=OutbreakNum[ind]+1

#     for i in range(3):
#         tmp=zO[i]
#         gu.opickle(meta['Paths']['bc1ha'] + '\\PEST_INFESTATION_POLY\\Outbreak_' + pest + '_' + str(i+1) + '.pkl',tmp)

#%% Rasterize planting from RESULTS

lNam='RSLT_ACTIVITY_TREATMENT_SVW'

if os.path.exists(meta['Paths']['bc1ha'] + '\\' + lNam)==False:
    os.mkdir(meta['Paths']['bc1ha'] + '\\' + lNam)

tv=np.arange(1910,2024,1)

# Start with planting with spatial from RESULTS (takes 15 min)
t0=time.time()
ats={}
ats['Path']=meta['Paths']['GDB']['Results']
ats['Layer']=lNam; # fiona.listlayers(ats['Path'])
ats['crs']=meta['Geos']['crs']
ats['Keep Geom']='On'
ats['Select Openings']=np.array([])
ats['SBC']=np.array(['PL'])
ats['FSC']=np.array([])
ats['SOC1']=np.array([])
ats['ROI']=[]
ats['gdf']=qgdb.Query_Openings(ats,[])
ats['gdf']=ats['gdf'][ats['gdf'].geometry!=None]
ats['gdf']=ats['gdf'].reset_index()
print( (time.time()-t0)/60 )

ats['gdf']['Year']=np.zeros(len(ats['gdf']))
for i in range(ats['gdf']['Year'].size):
    ats['gdf']['Year'][i]=int(ats['gdf']['ATU_COMPLETION_DATE'][i][0:4])

AreaPlanted=ats['gdf']['ACTUAL_TREATMENT_AREA']
NumTreesPlanted=ats['gdf']['ACTUAL_PLANTED_NUMBER']
ats['gdf']['SPH_Planted']=NumTreesPlanted/AreaPlanted

ats['gdf']=u1ha.CreateIdForCategoricalVariable(meta,'RSLT_ACTIVITY_TREATMENT_SVW','SILV_FUND_SOURCE_CODE',ats['gdf'])

# Add areas where FC is artificial

at={}
at['Path']=meta['Paths']['GDB']['Results']
at['Layer']='RSLT_ACTIVITY_TREATMENT_SVW'; # fiona.listlayers(at['Path'])
at['crs']=meta['Geos']['crs']
at['Keep Geom']='Off'
at['Select Openings']=np.array([])
at['SBC']=np.array(['PL'])
at['FSC']=np.array([])
at['SOC1']=np.array([])
at['ROI']=[]
at['gdf']=qgdb.Query_Openings(at,[])
at['gdf']=u1ha.CreateIdForCategoricalVariable(meta,'RSLT_ACTIVITY_TREATMENT_SVW','SILV_FUND_SOURCE_CODE',at['gdf'])

# Import all openings from opening layer
op={}
op['Path']=meta['Paths']['GDB']['Results']
op['Layer']='RSLT_OPENING_SVW'; # fiona.listlayers(op['Path'])
op['crs']=meta['Geos']['crs']
op['Keep Geom']='Off'
op['Select Openings']=np.array([])
op['SBC']=np.array([])
op['FSC']=np.array([])
op['SOC1']=np.array([])
op['ROI']=[]
op['gdf']=qgdb.Query_Openings(op,[])

ops={}
ops['Path']=meta['Paths']['GDB']['Results']
ops['Layer']='RSLT_OPENING_SVW'; # fiona.listlayers(op['Path'])
ops['crs']=meta['Geos']['crs']
ops['Keep Geom']='On'
ops['Select Openings']=np.array([])
ops['SBC']=np.array([])
ops['FSC']=np.array([])
ops['SOC1']=np.array([])
ops['ROI']=[]
ops['Drop Props']='On'
ops['gdf']=qgdb.Query_Openings(ops,[])

# Rasterize opening ID
shapes=((geom,value) for geom, value in zip(ops['gdf']['geometry'],ops['gdf']['OPENING_ID']))
zOP=np.zeros(zRef['Data'].shape,dtype=float)
burned=features.rasterize(shapes=shapes,fill=0,out=zOP,transform=zRef['Transform'])

# Import FC layer
fc={}
fc['Path']=meta['Paths']['GDB']['Results']
fc['Layer']='RSLT_FOREST_COVER_INV_SVW'; # fiona.listlayers(op['Path'])
fc['crs']=meta['Geos']['crs']
fc['Keep Geom']='Off'
fc['Select Openings']=np.array([])
fc['SBC']=np.array([])
fc['FSC']=np.array([])
fc['SOC1']=np.array([])
fc['ROI']=[]
fc['gdf']=qgdb.Query_Openings(fc,[])

# Create artificial stocking type mask
z0=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_FOREST_COVER_INV_SVW\\STOCKING_TYPE_CODE.tif')
zArt=zRef.copy()
zArt['Data']=np.zeros(zRef['Data'].shape,dtype='int8')
ind=np.where(z0['Data']==meta['LUT']['RSLT_FOREST_COVER_INV_SVW']['STOCKING_TYPE_CODE']['ART'])
zArt['Data'][ind]=1

# Reduce the size of rasters
indO1=np.where(zOP!=0)
zOP1=zOP[indO1]
zArt1=zArt['Data'][indO1]

# Index to planting and year
d={}
for iT in range(tv.size):
    d[tv[iT]]={'IndexToGrid':[],
               'SPH_Planted':[],
               'SILV_FUND_SOURCE_CODE':[],
               'ACTIVITY_TREATMENT_UNIT_ID':[]}

# Loop through opening layer spatial
for iOP in range(ops['gdf']['OPENING_ID'].size):
    print(iOP)

    indO2=np.where(zOP1==ops['gdf']['OPENING_ID'][iOP])[0]

    indAT=np.where( (at['gdf']['OPENING_ID']==ops['gdf']['OPENING_ID'][iOP]) & (at['gdf']['RESULTS_IND']=='Y') & (at['gdf']['SILV_METHOD_CODE']!='LAYOT') )[0]
    if indAT.size==0:
        continue

    Year=at['gdf']['Year'][indAT]
    FSC=at['gdf']['ID_SILV_FUND_SOURCE_CODE'][indAT]
    ATUID=at['gdf']['ACTIVITY_TREATMENT_UNIT_ID'][indAT]

    AreaPlanted=at['gdf']['ACTUAL_TREATMENT_AREA'][indAT]
    NumTreesPlanted=at['gdf']['ACTUAL_PLANTED_NUMBER'][indAT]
    SPH_Planted=NumTreesPlanted/AreaPlanted
    #AreaPlanted_Sum=np.sum(AreaPlanted)

    A_OP_FromGDB=op['gdf']['GEOMETRY_Area'][iOP]/10000
    A_OP=indO2.size
    ind_Art=np.where(zArt1[indO2]==1)[0]
    A_ART=ind_Art.size

    #if (A_ART>0) & (A_ART<AreaPlanted_Sum):
    #    break
    fA=np.max(AreaPlanted)/A_ART

    for iY in range(Year.size):
        if (Year[iY]<tv[0]) | (Year[iY]>tv[-1]):
            continue
        if (A_ART==0) | (fA>1.1):
            d[Year[iY]]['IndexToGrid'].append(indO2)
            d[Year[iY]]['SPH_Planted'].append(SPH_Planted[iY]*np.ones(indO2.size))
            d[Year[iY]]['SILV_FUND_SOURCE_CODE'].append(FSC[iY]*np.ones(indO2.size))
            d[Year[iY]]['ACTIVITY_TREATMENT_UNIT_ID'].append(ATUID[iY]*np.ones(indO2.size))
            #print('working 1')
        else:
            d[Year[iY]]['IndexToGrid'].append(indO2[ind_Art])
            d[Year[iY]]['SPH_Planted'].append(SPH_Planted[iY]*np.ones(ind_Art.size))
            d[Year[iY]]['SILV_FUND_SOURCE_CODE'].append(FSC[iY]*np.ones(ind_Art.size))
            d[Year[iY]]['ACTIVITY_TREATMENT_UNIT_ID'].append(ATUID[iY]*np.ones(ind_Art.size))
            #print('working 2')

#gu.opickle(r'C:\Users\rhember\Documents\Data\BC1ha\RSLT_ACTIVITY_TREATMENT_SVW\Process Planting\PlantingGapFillInfo.pkl',d)

# Pack

# Initialize rasters

N_Year=6

z0={'Year':{},'ACTIVITY_TREATMENT_UNIT_ID':{},'SILV_FUND_SOURCE_CODE':{},'SPH_Planted':{}}
for iY in range(N_Year):
    for k in z0.keys():
        if k=='ACTIVITY_TREATMENT_UNIT_ID':
            z0[k][iY+1]=np.zeros(zRef['Data'].shape,dtype='int32')
        else:
            z0[k][iY+1]=np.zeros(zRef['Data'].shape,dtype='int16')

vNam='PL_All'
#vNam='PL_FirstOnly'

# Initialize counter (to deal with repeat planting)
zCounter=np.zeros(zRef['Data'].shape,dtype=float)

for iT in range(tv.size):

    print(tv[iT])

    # Activity layer with spatial
    df0=ats['gdf'][ (ats['gdf']['Year']==tv[iT]) & (ats['gdf']['RESULTS_IND']=='Y') & (ats['gdf']['SILV_METHOD_CODE']!='LAYOT') ].copy()
    df0=df0[df0.geometry!=None]
    df0=df0.reset_index()

    z1={}
    for k in z0.keys():
        z1[k]=np.zeros(zRef['Data'].shape,dtype=float)

    if len(df0)>0:
        shapes=((geom,value) for geom, value in zip(df0['geometry'],df0['ACTIVITY_TREATMENT_UNIT_ID']))
        burnedATUID=features.rasterize(shapes=shapes,fill=0,out=z1['ACTIVITY_TREATMENT_UNIT_ID'],transform=zRef['Transform'])

        shapes=((geom,value) for geom, value in zip(df0['geometry'],df0['ID_SILV_FUND_SOURCE_CODE']))
        burnedFSC=features.rasterize(shapes=shapes,fill=0,out=z1['SILV_FUND_SOURCE_CODE'],transform=zRef['Transform'])

        shapes=((geom,value) for geom, value in zip(df0['geometry'],df0['SPH_Planted']))
        burnedSPH=features.rasterize(shapes=shapes,fill=0,out=z1['SPH_Planted'],transform=zRef['Transform'])

    # Add activities without spatial
    for iD in range(len(d[tv[iT]]['IndexToGrid'])):
        ind2=d[tv[iT]]['IndexToGrid'][iD]
        z1['ACTIVITY_TREATMENT_UNIT_ID'][ indO1[0][ind2],indO1[1][ind2]  ]=d[tv[iT]]['ACTIVITY_TREATMENT_UNIT_ID'][iD]
        z1['SILV_FUND_SOURCE_CODE'][ indO1[0][ind2],indO1[1][ind2]  ]=d[tv[iT]]['SILV_FUND_SOURCE_CODE'][iD]
        z1['SPH_Planted'][ indO1[0][ind2],indO1[1][ind2]  ]=d[tv[iT]]['SPH_Planted'][iD]

    # Update counter
    ind=np.where( (z1['SILV_FUND_SOURCE_CODE']>0) | (zCounter>0) )
    zCounter[ind]=zCounter[ind]+1

    ind=np.where( (zCounter>15) )
    zCounter[ind]=0

    # Populate final grids

    ind=np.where( (zCounter<=1) & (z0['Year'][1]==0) & (z1['SILV_FUND_SOURCE_CODE']!=0) )
    z0['Year'][1][ind]=tv[iT]
    for k in z1.keys():
        if k=='Year':
            continue
        z0[k][1][ind]=z1[k][ind]

    ind=np.where( (zCounter<=1) & (z0['Year'][1]!=0) & (z0['Year'][1]!=tv[iT]) & (z0['Year'][2]==0) & (z1['SILV_FUND_SOURCE_CODE']!=0) )
    z0['Year'][2][ind]=tv[iT]
    for k in z1.keys():
        if k=='Year':
            continue
        z0[k][2][ind]=z1[k][ind]

    ind=np.where( (zCounter<=1) & (z0['Year'][2]!=0) & (z0['Year'][2]!=tv[iT]) & (z0['Year'][3]==0) & (z1['SILV_FUND_SOURCE_CODE']!=0) )
    z0['Year'][3][ind]=tv[iT]
    for k in z1.keys():
        if k=='Year':
            continue
        z0[k][3][ind]=z1[k][ind]

    ind=np.where( (zCounter<=1) & (z0['Year'][3]!=0) & (z0['Year'][3]!=tv[iT]) & (z0['Year'][4]==0) & (z1['SILV_FUND_SOURCE_CODE']!=0) )
    z0['Year'][4][ind]=tv[iT]
    for k in z1.keys():
        if k=='Year':
            continue
        z0[k][4][ind]=z1[k][ind]

    ind=np.where( (zCounter<=1) & (z0['Year'][4]!=0) & (z0['Year'][4]!=tv[iT]) & (z0['Year'][5]==0) & (z1['SILV_FUND_SOURCE_CODE']!=0) )
    z0['Year'][5][ind]=tv[iT]
    for k in z1.keys():
        if k=='Year':
            continue
        z0[k][5][ind]=z1[k][ind]

    ind=np.where( (zCounter<=1) & (z0['Year'][5]!=0) & (z0['Year'][5]!=tv[iT]) & (z1['SILV_FUND_SOURCE_CODE']!=0) )
    z0['Year'][6][ind]=tv[iT]
    for k in z1.keys():
        if k=='Year':
            continue
        z0[k][6][ind]=z1[k][ind]

# Save to file
for iY in range(N_Year):
    z1=zRef.copy()
    z1['Data']=z0['Year'][iY+1].astype('int16')
    gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\' + lNam + '\\' + vNam + '_' + str(iY+1) + '_Year.tif')
    z1=zRef.copy()
    z1['Data']=z0['ACTIVITY_TREATMENT_UNIT_ID'][iY+1].astype('int32')
    gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\' + lNam + '\\' + vNam + '_' + str(iY+1) + '_ACTIVITY_TREATMENT_UNIT_ID.tif')
    z1=zRef.copy()
    z1['Data']=z0['SILV_FUND_SOURCE_CODE'][iY+1].astype('int16')
    gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\' + lNam + '\\' + vNam + '_' + str(iY+1) + '_SILV_FUND_SOURCE_CODE.tif')
    z1=zRef.copy()
    z1['Data']=z0['SPH_Planted'][iY+1].astype('int16')
    gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\' + lNam + '\\' + vNam + '_' + str(iY+1) + '_SPH_Planted.tif')

# Plot time series to confirm it worked
lNam='RSLT_ACTIVITY_TREATMENT_SVW'
vNam='PL_FirstOnly'
nPack=6
tv1,N=u1ha.TimeSeriesOccurrenceFromPackedEvents(meta,lNam,vNam,nPack)
plt.plot(tv1,N,'-bo')

#%% Planting layer (species and genetic worth)

pl={}
pl['Path']=meta['Paths']['GDB']['Results']
pl['Layer']='RSLT_PLANTING_SVW'; # fiona.listlayers(at['Path'])
pl['crs']=meta['Geos']['crs']
pl['Keep Geom']='Off'
pl['Select Openings']=np.array([])
pl['SBC']=np.array([])
pl['STC']=np.array([])
pl['SMC']=np.array([])
pl['FSC']=np.array([])
pl['SOC1']=np.array([])
pl['ROI']=[]
pl['gdf']=qgdb.Query_Openings(pl,[])
pl['gdf']=u1ha.CreateIdForCategoricalVariable(meta,'RSLT_PLANTING_SVW','SILV_TREE_SPECIES_CODE',pl['gdf'])

dGW=gu.ReadExcel(meta['Paths']['Model']['Code'] + '\\Parameters\\Parameters_Seedlot_GW.xlsx')

def round_retain_sum(x):
    N=np.round(np.sum(x)).astype(int)
    y=x.astype(int)
    M=np.sum(y)
    K=N-M
    z=y-x
    if K!=0:
        idx=np.argpartition(z,K)[:K]
        y[idx]+=1
    return y

# Get information by ATU ID
uAT=np.unique(pl['gdf']['ACTIVITY_TREATMENT_UNIT_ID'])

plAT={}
plAT['ACTIVITY_TREATMENT_UNIT_ID']=uAT
vL=['PL_SPECIES_CD','PL_SPECIES_PCT','PL_SPECIES_GW']
for v in vL:
    for i in range(6):
        plAT[v + str(i+1)]=np.zeros(uAT.size,dtype='int16')

# Get weighted average genetic worth by AT ID
N_no_spc=0
N_no_num_planted=0
N_sn_not_found=0
for iU in range(uAT.size):
    iAT=np.where(pl['gdf']['ACTIVITY_TREATMENT_UNIT_ID']==uAT[iU])[0]

    cd0=pl['gdf']['ID_SILV_TREE_SPECIES_CODE'][iAT]
    pct0=np.nan_to_num(pl['gdf']['NUMBER_PLANTED'][iAT]/np.sum(pl['gdf']['NUMBER_PLANTED'][iAT])*100)
    sn0=np.nan_to_num(pl['gdf']['SEEDLOT_NUMBER'][iAT])

    if np.sum(cd0)==0:
        N_no_spc=N_no_spc+1

    if np.sum(pct0)==0:
        N_no_num_planted=N_no_num_planted+1
        pct0=(100/pct0.size)*np.ones(pct0.size)

    # Genetic worth
    gw0=np.zeros(pct0.size)
    for j in range(pct0.size):
        ind=np.where(dGW['SEEDLOT_NUMBER']==sn0[j])[0]
        if ind.size!=0:
            gw0[j]=dGW['GENETIC_WORTH_RTNG'][ind[0]]
        else:
            N_sn_not_found=N_sn_not_found+1

    # Dissolve to unique species combinations and calculate weighted average genetic worth
    cd1=np.unique(cd0)
    pct1=np.zeros(cd1.size)
    gw1=np.zeros(cd1.size)
    for iCD in range(cd1.size):
        ind=np.where(cd0==cd1[iCD])[0]
        pct1[iCD]=np.sum(pct0[ind])
        gw1[iCD]=np.sum(pct0[ind]*gw0[ind])/np.sum(pct0[ind])

    # Divi
    gw1=np.nan_to_num(gw1)

    # Sort
    iSort=np.argsort(-pct1)
    cd1=cd1[iSort]
    pct1=pct1[iSort]
    gw1=gw1[iSort]

    # Percents will be stored as integers, ensure they will sum to 100
    pct1=np.round(pct1,decimals=0).astype('int16')
    if (np.sum(pct1)!=100) & (pct1.size>1):
        pct1=round_retain_sum(pct1)

    # Populate structure
    if cd1.size>0:
        plAT['PL_SPECIES_CD1'][iU]=cd1[0]
        plAT['PL_SPECIES_PCT1'][iU]=pct1[0]
        plAT['PL_SPECIES_GW1'][iU]=gw1[0]
    if cd1.size>1:
        plAT['PL_SPECIES_CD2'][iU]=cd1[1]
        plAT['PL_SPECIES_PCT2'][iU]=pct1[1]
        plAT['PL_SPECIES_GW2'][iU]=gw1[1]
    if cd1.size>2:
        plAT['PL_SPECIES_CD3'][iU]=cd1[2]
        plAT['PL_SPECIES_PCT3'][iU]=pct1[2]
        plAT['PL_SPECIES_GW3'][iU]=gw1[2]
    if cd1.size>3:
        plAT['PL_SPECIES_CD4'][iU]=cd1[3]
        plAT['PL_SPECIES_PCT4'][iU]=pct1[3]
        plAT['PL_SPECIES_GW4'][iU]=gw1[3]
    if cd1.size>4:
        plAT['PL_SPECIES_CD5'][iU]=cd1[4]
        plAT['PL_SPECIES_PCT5'][iU]=pct1[4]
        plAT['PL_SPECIES_GW5'][iU]=gw1[4]
    if cd1.size>5:
        plAT['PL_SPECIES_CD6'][iU]=cd1[5]
        plAT['PL_SPECIES_PCT6'][iU]=pct1[5]
        plAT['PL_SPECIES_GW6'][iU]=gw1[5]

# Populate packed layers
N_Year=6
lNam='RSLT_ACTIVITY_TREATMENT_SVW'
vNam='PL_All'

for iY in range(N_Year):
    t0=time.time()
    z=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\' + lNam + '\\' + vNam + '_' + str(iY+1) + '_ACTIVITY_TREATMENT_UNIT_ID.tif')
    uAT1=gu.IndicesFromUniqueArrayValues(z['Data'].flatten())

    z0={}
    for iS in range(6):
        z0['PL_SPECIES_CD' + str(iS+1)]=np.zeros(zRef['Data'].shape,dtype='int16').flatten()
        z0['PL_SPECIES_PCT' + str(iS+1)]=np.zeros(zRef['Data'].shape,dtype='int16').flatten()
        z0['PL_SPECIES_GW' + str(iS+1)]=np.zeros(zRef['Data'].shape,dtype='int16').flatten()

    for k in uAT1.keys():
        ind=np.where(plAT['ACTIVITY_TREATMENT_UNIT_ID']==k)[0]
        if ind.size==0:
            continue
        iIndex=uAT1[k]
        for iS in range(6):
            z0['PL_SPECIES_CD' + str(iS+1)][iIndex]=plAT['PL_SPECIES_CD' + str(iS+1)][ind]
            z0['PL_SPECIES_PCT' + str(iS+1)][iIndex]=plAT['PL_SPECIES_PCT' + str(iS+1)][ind]
            z0['PL_SPECIES_GW' + str(iS+1)][iIndex]=plAT['PL_SPECIES_GW' + str(iS+1)][ind]

    for iS in range(6):
        z=zRef.copy()
        z['Data']=np.reshape(z0['PL_SPECIES_CD' + str(iS+1)],(zRef['Data'].shape))
        gis.SaveGeoTiff(z,meta['Paths']['bc1ha'] + '\\' + lNam + '\\' + vNam + '_' + str(iY+1) + '_PL_SPECIES_CD' + str(iS+1) + '.tif')
        z=zRef.copy()
        z['Data']=np.reshape(z0['PL_SPECIES_PCT' + str(iS+1)],(zRef['Data'].shape))
        gis.SaveGeoTiff(z,meta['Paths']['bc1ha'] + '\\' + lNam + '\\' + vNam + '_' + str(iY+1) + '_PL_SPECIES_PCT' + str(iS+1) + '.tif')
        z=zRef.copy()
        z['Data']=np.reshape(z0['PL_SPECIES_GW' + str(iS+1)],(zRef['Data'].shape))
        gis.SaveGeoTiff(z,meta['Paths']['bc1ha'] + '\\' + lNam + '\\' + vNam + '_' + str(iY+1) + '_PL_SPECIES_GW' + str(iS+1) + '.tif')

    print((time.time()-t0)/60)

#%% Rasterize fertilization

lNam='RSLT_ACTIVITY_TREATMENT_SVW'

if os.path.exists(meta['Paths']['bc1ha'] + '\\' + lNam)==False:
    os.mkdir(meta['Paths']['bc1ha'] + '\\' + lNam)

tv=np.arange(1970,2024,1)

# Start with planting with spatial from RESULTS (takes 15 min)
t0=time.time()
ats={}
ats['Path']=meta['Paths']['GDB']['Results']
ats['Layer']=lNam; # fiona.listlayers(ats['Path'])
ats['crs']=meta['Geos']['crs']
ats['Keep Geom']='On'
ats['Select Openings']=np.array([])
ats['SBC']=np.array(['FE'])
ats['STC']=np.array(['CA'])
ats['FSC']=np.array([])
ats['SOC1']=np.array([])
ats['ROI']=[]
ats['gdf']=qgdb.Query_Openings(ats,[])
ats['gdf']=ats['gdf'][ats['gdf'].geometry!=None]
ats['gdf']=ats['gdf'].reset_index()
print( (time.time()-t0)/60 )

ats['gdf']=ats['gdf'][ (ats['gdf']['RESULTS_IND']=='Y') & (ats['gdf']['ATU_COMPLETION_DATE']!=None) ]
ats['gdf']=ats['gdf'].reset_index()

ats['gdf']['Year']=np.zeros(len(ats['gdf']))
for i in range(ats['gdf']['Year'].size):
    ats['gdf']['Year'][i]=int(ats['gdf']['ATU_COMPLETION_DATE'][i][0:4])

ats['gdf']=u1ha.CreateIdForCategoricalVariable(meta,'RSLT_ACTIVITY_TREATMENT_SVW','SILV_FUND_SOURCE_CODE',ats['gdf'])

# Add areas where FC is artificial

at={}
at['Path']=meta['Paths']['GDB']['Results']
at['Layer']='RSLT_ACTIVITY_TREATMENT_SVW'; # fiona.listlayers(at['Path'])
at['crs']=meta['Geos']['crs']
at['Keep Geom']='Off'
at['Select Openings']=np.array([])
at['SBC']=np.array(['FE'])
at['STC']=np.array(['CA'])
at['FSC']=np.array([])
at['SOC1']=np.array([])
at['ROI']=[]
at['gdf']=qgdb.Query_Openings(at,[])
at['gdf']=u1ha.CreateIdForCategoricalVariable(meta,'RSLT_ACTIVITY_TREATMENT_SVW','SILV_FUND_SOURCE_CODE',at['gdf'])

# Import all openings from opening layer
ops={}
ops['Path']=meta['Paths']['GDB']['Results']
ops['Layer']='RSLT_OPENING_SVW'; # fiona.listlayers(op['Path'])
ops['crs']=meta['Geos']['crs']
ops['Keep Geom']='On'
ops['Select Openings']=np.array([])
ops['SBC']=np.array([])
ops['STC']=np.array([])
ops['FSC']=np.array([])
ops['SOC1']=np.array([])
ops['ROI']=[]
ops['Drop Props']='On'
ops['gdf']=qgdb.Query_Openings(ops,[])

# Rasterize opening ID
shapes=((geom,value) for geom, value in zip(ops['gdf']['geometry'],ops['gdf']['OPENING_ID']))
zOP=np.zeros(zRef['Data'].shape,dtype=float)
burned=features.rasterize(shapes=shapes,fill=0,out=zOP,transform=zRef['Transform'])

# Reduce the size of rasters
indO1=np.where(zOP!=0)
zOP1=zOP[indO1]

# Index to planting and year
d={}
for iT in range(tv.size):
    d[tv[iT]]={'IndexToGrid':[],
               'SILV_FUND_SOURCE_CODE':[],
               'ACTIVITY_TREATMENT_UNIT_ID':[]}

for iOP in range(ops['gdf']['OPENING_ID'].size):
    print(iOP)

    indO2=np.where(zOP1==ops['gdf']['OPENING_ID'][iOP])[0]

    indAT=np.where( (at['gdf']['OPENING_ID']==ops['gdf']['OPENING_ID'][iOP]) )[0]
    if indAT.size==0:
        continue

    Year=at['gdf']['Year'][indAT]
    FSC=at['gdf']['ID_SILV_FUND_SOURCE_CODE'][indAT]
    ATUID=at['gdf']['ACTIVITY_TREATMENT_UNIT_ID'][indAT]
    A_OP=indO2.size
    for iY in range(Year.size):
        if (Year[iY]<tv[0]) | (Year[iY]>tv[-1]):
            continue
        d[Year[iY]]['IndexToGrid'].append(indO2)
        d[Year[iY]]['SILV_FUND_SOURCE_CODE'].append(FSC[iY]*np.ones(indO2.size))
        d[Year[iY]]['ACTIVITY_TREATMENT_UNIT_ID'].append(ATUID[iY]*np.ones(indO2.size))

# Pack

# Initialize rasters

N_Year=6

z0={'Year':{},'ACTIVITY_TREATMENT_UNIT_ID':{},'SILV_FUND_SOURCE_CODE':{}}
for iY in range(N_Year):
    for k in z0.keys():
        if k=='ACTIVITY_TREATMENT_UNIT_ID':
            z0[k][iY+1]=np.zeros(zRef['Data'].shape,dtype='int32')
        else:
            z0[k][iY+1]=np.zeros(zRef['Data'].shape,dtype='int16')

vNam='FECA'

for iT in range(tv.size):

    print(tv[iT])

    # Activity layer with spatial
    df0=ats['gdf'][ (ats['gdf']['Year']==tv[iT]) ].copy()
    df0=df0[df0.geometry!=None]
    #df0=df0.reset_index()

    z1={}
    for k in z0.keys():
        z1[k]=np.zeros(zRef['Data'].shape,dtype=float)

    if len(df0)>0:
        shapes=((geom,value) for geom, value in zip(df0['geometry'],df0['ACTIVITY_TREATMENT_UNIT_ID']))
        burnedATUID=features.rasterize(shapes=shapes,fill=0,out=z1['ACTIVITY_TREATMENT_UNIT_ID'],transform=zRef['Transform'])

        shapes=((geom,value) for geom, value in zip(df0['geometry'],df0['ID_SILV_FUND_SOURCE_CODE']))
        burnedFSC=features.rasterize(shapes=shapes,fill=0,out=z1['SILV_FUND_SOURCE_CODE'],transform=zRef['Transform'])

    # Add activities without spatial
    for iD in range(len(d[tv[iT]]['IndexToGrid'])):
        ind2=d[tv[iT]]['IndexToGrid'][iD]
        z1['ACTIVITY_TREATMENT_UNIT_ID'][ indO1[0][ind2],indO1[1][ind2]  ]=d[tv[iT]]['ACTIVITY_TREATMENT_UNIT_ID'][iD]
        z1['SILV_FUND_SOURCE_CODE'][ indO1[0][ind2],indO1[1][ind2]  ]=d[tv[iT]]['SILV_FUND_SOURCE_CODE'][iD]

    # Populate final grids

    ind=np.where( (z0['Year'][1]==0) & (z1['SILV_FUND_SOURCE_CODE']!=0) )
    z0['Year'][1][ind]=tv[iT]
    for k in z1.keys():
        if k=='Year':
            continue
        z0[k][1][ind]=z1[k][ind]

    ind=np.where( (z0['Year'][1]!=0) & (z0['Year'][1]!=tv[iT]) & (z0['Year'][2]==0) & (z1['SILV_FUND_SOURCE_CODE']!=0) )
    z0['Year'][2][ind]=tv[iT]
    for k in z1.keys():
        if k=='Year':
            continue
        z0[k][2][ind]=z1[k][ind]

    ind=np.where( (z0['Year'][2]!=0) & (z0['Year'][2]!=tv[iT]) & (z0['Year'][3]==0) & (z1['SILV_FUND_SOURCE_CODE']!=0) )
    z0['Year'][3][ind]=tv[iT]
    for k in z1.keys():
        if k=='Year':
            continue
        z0[k][3][ind]=z1[k][ind]

    ind=np.where( (z0['Year'][3]!=0) & (z0['Year'][3]!=tv[iT]) & (z0['Year'][4]==0) & (z1['SILV_FUND_SOURCE_CODE']!=0) )
    z0['Year'][4][ind]=tv[iT]
    for k in z1.keys():
        if k=='Year':
            continue
        z0[k][4][ind]=z1[k][ind]

    ind=np.where( (z0['Year'][4]!=0) & (z0['Year'][4]!=tv[iT]) & (z0['Year'][5]==0) & (z1['SILV_FUND_SOURCE_CODE']!=0) )
    z0['Year'][5][ind]=tv[iT]
    for k in z1.keys():
        if k=='Year':
            continue
        z0[k][5][ind]=z1[k][ind]

    ind=np.where( (z0['Year'][5]!=0) & (z0['Year'][5]!=tv[iT]) & (z1['SILV_FUND_SOURCE_CODE']!=0) )
    z0['Year'][6][ind]=tv[iT]
    for k in z1.keys():
        if k=='Year':
            continue
        z0[k][6][ind]=z1[k][ind]

for iY in range(N_Year):
    z1=zRef.copy()
    z1['Data']=z0['Year'][iY+1].astype('int16')
    gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\' + lNam + '\\' + vNam + '_' + str(iY+1) + '_Year.tif')
    z1=zRef.copy()
    z1['Data']=z0['ACTIVITY_TREATMENT_UNIT_ID'][iY+1].astype('int32')
    gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\' + lNam + '\\' + vNam + '_' + str(iY+1) + '_ACTIVITY_TREATMENT_UNIT_ID.tif')
    z1=zRef.copy()
    z1['Data']=z0['SILV_FUND_SOURCE_CODE'][iY+1].astype('int16')
    gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\' + lNam + '\\' + vNam + '_' + str(iY+1) + '_SILV_FUND_SOURCE_CODE.tif')

# Mask all
z1=zRef.copy()
z1['Data']=np.zeros(zRef['Data'].shape,dtype='int8')
for iY in range(N_Year):
    z0=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\' + lNam + '\\' + vNam + '_' + str(iY+1) + '_Year.tif')
    ind=np.where(z0['Data']>0)
    z1['Data'][ind]=1
gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\' + lNam + '\\' + vNam + '_MaskAll.tif')

# Year last
z1=zRef.copy()
z1['Data']=np.zeros(zRef['Data'].shape,dtype='int16')
for iY in range(N_Year):
    z0=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\' + lNam + '\\' + vNam + '_' + str(iY+1) + '_Year.tif')
    ind=np.where(z0['Data']>0)
    z1['Data'][ind]=z0['Data'][ind]
gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\' + lNam + '\\' + vNam + '_YearLast.tif')

# Plot time series to confirm it worked
lNam='RSLT_ACTIVITY_TREATMENT_SVW'
vNam='FECA'
nPack=6
tv1,N=u1ha.TimeSeriesOccurrenceFromPackedEvents(meta,lNam,vNam,nPack)
plt.plot(tv1,N,'-bo')

#%% Rasterize knockdown

lNam='RSLT_ACTIVITY_TREATMENT_SVW'

if os.path.exists(meta['Paths']['bc1ha'] + '\\' + lNam)==False:
    os.mkdir(meta['Paths']['bc1ha'] + '\\' + lNam)

tv=np.arange(1970,2024,1)

# Start with planting with spatial from RESULTS (takes 15 min)
t0=time.time()
ats={}
ats['Path']=meta['Paths']['GDB']['Results']
ats['Layer']=lNam; # fiona.listlayers(ats['Path'])
ats['crs']=meta['Geos']['crs']
ats['Keep Geom']='On'
ats['Select Openings']=np.array([])
ats['SBC']=np.array(['SP'])
ats['STC']=np.array(['ME'])
ats['SMC']=np.array(['CABLE','GUARD','HARV','MDOWN','PUSH'])
ats['FSC']=np.array(meta['Param']['BE']['FSC']['NO List Name'])
ats['SOC1']=np.array([])
ats['ROI']=[]
ats['gdf']=qgdb.Query_Openings(ats,[])
ats['gdf']=ats['gdf'][ats['gdf'].geometry!=None]
ats['gdf']=ats['gdf'].reset_index()
print( (time.time()-t0)/60 )

ats['gdf']=ats['gdf'][ (ats['gdf']['RESULTS_IND']=='Y') & (ats['gdf']['ATU_COMPLETION_DATE']!=None) ]
ats['gdf']=ats['gdf'].reset_index()

ats['gdf']['Year']=np.zeros(len(ats['gdf']))
for i in range(ats['gdf']['Year'].size):
    ats['gdf']['Year'][i]=int(ats['gdf']['ATU_COMPLETION_DATE'][i][0:4])

ats['gdf']=u1ha.CreateIdForCategoricalVariable(meta,'RSLT_ACTIVITY_TREATMENT_SVW','SILV_FUND_SOURCE_CODE',ats['gdf'])

# Reconstruct location of activities with no reported spatial

at={}
at['Path']=meta['Paths']['GDB']['Results']
at['Layer']='RSLT_ACTIVITY_TREATMENT_SVW'; # fiona.listlayers(at['Path'])
at['crs']=meta['Geos']['crs']
at['Keep Geom']='Off'
at['Select Openings']=np.array([])
at['SBC']=np.array(['SP'])
at['STC']=np.array(['ME'])
at['SMC']=np.array(['CABLE','GUARD','HARV','MDOWN','PUSH'])
at['FSC']=np.array(meta['Param']['BE']['FSC']['NO List Name'])
at['SOC1']=np.array([])
at['ROI']=[]
at['gdf']=qgdb.Query_Openings(at,[])
at['gdf']=u1ha.CreateIdForCategoricalVariable(meta,'RSLT_ACTIVITY_TREATMENT_SVW','SILV_FUND_SOURCE_CODE',at['gdf'])

# Import all openings from opening layer
ops={}
ops['Path']=meta['Paths']['GDB']['Results']
ops['Layer']='RSLT_OPENING_SVW'; # fiona.listlayers(op['Path'])
ops['crs']=meta['Geos']['crs']
ops['Keep Geom']='On'
ops['Select Openings']=np.array([])
ops['SBC']=np.array([])
ops['STC']=np.array([])
ops['SMC']=np.array([])
ops['FSC']=np.array([])
ops['SOC1']=np.array([])
ops['ROI']=[]
ops['Drop Props']='On'
ops['gdf']=qgdb.Query_Openings(ops,[])

# Rasterize opening ID
shapes=((geom,value) for geom, value in zip(ops['gdf']['geometry'],ops['gdf']['OPENING_ID']))
zOP=np.zeros(zRef['Data'].shape,dtype=float)
burned=features.rasterize(shapes=shapes,fill=0,out=zOP,transform=zRef['Transform'])

# Reduce the size of rasters
indO1=np.where(zOP!=0)
zOP1=zOP[indO1]

# Index to year
d={}
for iT in range(tv.size):
    d[tv[iT]]={'IndexToGrid':[],
               'SILV_FUND_SOURCE_CODE':[],
               'ACTIVITY_TREATMENT_UNIT_ID':[]}

for iOP in range(ops['gdf']['OPENING_ID'].size):
    print(iOP)
    indO2=np.where(zOP1==ops['gdf']['OPENING_ID'][iOP])[0]
    indAT=np.where( (at['gdf']['OPENING_ID']==ops['gdf']['OPENING_ID'][iOP]) )[0]
    if indAT.size==0:
        continue
    Year=at['gdf']['Year'][indAT]
    FSC=at['gdf']['ID_SILV_FUND_SOURCE_CODE'][indAT]
    ATUID=at['gdf']['ACTIVITY_TREATMENT_UNIT_ID'][indAT]
    A_OP=indO2.size
    for iY in range(Year.size):
        if (Year[iY]<tv[0]) | (Year[iY]>tv[-1]):
            continue
        d[Year[iY]]['IndexToGrid'].append(indO2)
        d[Year[iY]]['SILV_FUND_SOURCE_CODE'].append(FSC[iY]*np.ones(indO2.size))
        d[Year[iY]]['ACTIVITY_TREATMENT_UNIT_ID'].append(ATUID[iY]*np.ones(indO2.size))

# Pack

vNam='SP_KD'

N_Year=6

z0={'Year':{},'ACTIVITY_TREATMENT_UNIT_ID':{},'SILV_FUND_SOURCE_CODE':{}}
for iY in range(N_Year):
    for k in z0.keys():
        if k=='ACTIVITY_TREATMENT_UNIT_ID':
            z0[k][iY+1]=np.zeros(zRef['Data'].shape,dtype='int32')
        else:
            z0[k][iY+1]=np.zeros(zRef['Data'].shape,dtype='int16')

for iT in range(tv.size):
    print(tv[iT])

    # Activity layer with spatial
    df0=ats['gdf'][ (ats['gdf']['Year']==tv[iT]) ].copy()
    df0=df0[df0.geometry!=None]
    #df0=df0.reset_index()

    z1={}
    for k in z0.keys():
        z1[k]=np.zeros(zRef['Data'].shape,dtype=float)

    if len(df0)>0:
        shapes=((geom,value) for geom, value in zip(df0['geometry'],df0['ACTIVITY_TREATMENT_UNIT_ID']))
        burnedATUID=features.rasterize(shapes=shapes,fill=0,out=z1['ACTIVITY_TREATMENT_UNIT_ID'],transform=zRef['Transform'])

        shapes=((geom,value) for geom, value in zip(df0['geometry'],df0['ID_SILV_FUND_SOURCE_CODE']))
        burnedFSC=features.rasterize(shapes=shapes,fill=0,out=z1['SILV_FUND_SOURCE_CODE'],transform=zRef['Transform'])

    # Add activities without spatial
    for iD in range(len(d[tv[iT]]['IndexToGrid'])):
        ind2=d[tv[iT]]['IndexToGrid'][iD]
        z1['ACTIVITY_TREATMENT_UNIT_ID'][ indO1[0][ind2],indO1[1][ind2]  ]=d[tv[iT]]['ACTIVITY_TREATMENT_UNIT_ID'][iD]
        z1['SILV_FUND_SOURCE_CODE'][ indO1[0][ind2],indO1[1][ind2]  ]=d[tv[iT]]['SILV_FUND_SOURCE_CODE'][iD]

    # Populate final grids

    ind=np.where( (z0['Year'][1]==0) & (z1['SILV_FUND_SOURCE_CODE']!=0) )
    z0['Year'][1][ind]=tv[iT]
    for k in z1.keys():
        if k=='Year':
            continue
        z0[k][1][ind]=z1[k][ind]

    ind=np.where( (z0['Year'][1]!=0) & (z0['Year'][1]!=tv[iT]) & (z0['Year'][2]==0) & (z1['SILV_FUND_SOURCE_CODE']!=0) )
    z0['Year'][2][ind]=tv[iT]
    for k in z1.keys():
        if k=='Year':
            continue
        z0[k][2][ind]=z1[k][ind]

    ind=np.where( (z0['Year'][2]!=0) & (z0['Year'][2]!=tv[iT]) & (z0['Year'][3]==0) & (z1['SILV_FUND_SOURCE_CODE']!=0) )
    z0['Year'][3][ind]=tv[iT]
    for k in z1.keys():
        if k=='Year':
            continue
        z0[k][3][ind]=z1[k][ind]

    ind=np.where( (z0['Year'][3]!=0) & (z0['Year'][3]!=tv[iT]) & (z0['Year'][4]==0) & (z1['SILV_FUND_SOURCE_CODE']!=0) )
    z0['Year'][4][ind]=tv[iT]
    for k in z1.keys():
        if k=='Year':
            continue
        z0[k][4][ind]=z1[k][ind]

    ind=np.where( (z0['Year'][4]!=0) & (z0['Year'][4]!=tv[iT]) & (z0['Year'][5]==0) & (z1['SILV_FUND_SOURCE_CODE']!=0) )
    z0['Year'][5][ind]=tv[iT]
    for k in z1.keys():
        if k=='Year':
            continue
        z0[k][5][ind]=z1[k][ind]

    ind=np.where( (z0['Year'][5]!=0) & (z0['Year'][5]!=tv[iT]) & (z1['SILV_FUND_SOURCE_CODE']!=0) )
    z0['Year'][6][ind]=tv[iT]
    for k in z1.keys():
        if k=='Year':
            continue
        z0[k][6][ind]=z1[k][ind]

for iY in range(N_Year):
    z1=zRef.copy()
    z1['Data']=z0['Year'][iY+1].astype('int16')
    gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\' + lNam + '\\' + vNam + '_' + str(iY+1) + '_Year.tif')
    z1=zRef.copy()
    z1['Data']=z0['ACTIVITY_TREATMENT_UNIT_ID'][iY+1].astype('int32')
    gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\' + lNam + '\\' + vNam + '_' + str(iY+1) + '_ACTIVITY_TREATMENT_UNIT_ID.tif')
    z1=zRef.copy()
    z1['Data']=z0['SILV_FUND_SOURCE_CODE'][iY+1].astype('int16')
    gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\' + lNam + '\\' + vNam + '_' + str(iY+1) + '_SILV_FUND_SOURCE_CODE.tif')

# Mask all
z1=zRef.copy()
z1['Data']=np.zeros(zRef['Data'].shape,dtype='int8')
for iY in range(N_Year):
    z0=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\' + lNam + '\\' + vNam + '_' + str(iY+1) + '_Year.tif')
    ind=np.where(z0['Data']>0)
    z1['Data'][ind]=1
gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\' + lNam + '\\' + vNam + '_MaskAll.tif')

# Year last
z1=zRef.copy()
z1['Data']=np.zeros(zRef['Data'].shape,dtype='int16')
for iY in range(N_Year):
    z0=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\' + lNam + '\\' + vNam + '_' + str(iY+1) + '_Year.tif')
    ind=np.where(z0['Data']>0)
    z1['Data'][ind]=z0['Data'][ind]
gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\' + lNam + '\\' + vNam + '_YearLast.tif')

# Plot time series to confirm it worked
lNam='RSLT_ACTIVITY_TREATMENT_SVW'
vNam='SP_KD'
nPack=6
tv1,N=u1ha.TimeSeriesOccurrenceFromPackedEvents(meta,lNam,vNam,nPack)
plt.plot(tv1,N,'-bo')

#%% Planting (Non-obligation by project type)

zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
meta['Graphics']['Map']['RGSF']=1

# Mask
zMask=zRef.copy()
zMask['Data']=np.zeros(zRef['Data'].shape,dtype='int16')
for iEY in range(6):
    zFSC=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_ALL_' + str(iEY+1) + '_SILV_FUND_SOURCE_CODE.tif')
    ind=np.where( (np.isin(zFSC['Data'],meta['Param']['BE']['FSC']['NO List ID'])==True) )
    zMask['Data'][ind]=1
gis.SaveGeoTiff(zMask,meta['Paths']['bc1ha'] + '\\Management\\PL_NonOb_MaskAll.tif')

th_Fill=15
yr_th=8
tv=np.arange(1960,2023,1)
zD=u1ha.Import_Raster(meta,[],['ibm_yr','fire_yr','harv_yr_con1','kd_yr'])

# First and last instances of planting, plus compacted indices to each year
zPL_First=zRef.copy()
zPL_First['Data']=3000*np.ones(zRef['Data'].shape,dtype='int16')
zPL_Last=zRef.copy()
zPL_Last['Data']=np.zeros(zRef['Data'].shape,dtype='int16')
idxPL={}
for iT in range(tv.size):
    idxPL[tv[iT]]=(np.array([],dtype=int),np.array([],dtype=int))
for iEY in range(6):
    print(iEY)
    zYr=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_ALL_' + str(iEY+1) + '_Year.tif')
    zFSC=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_ALL_' + str(iEY+1) + '_SILV_FUND_SOURCE_CODE.tif')
    ind=np.where( (np.isin(zFSC['Data'],meta['Param']['BE']['FSC']['NO List ID'])==True) ) # (zSPH['Data']>750)
    zPL_First['Data'][ind]=np.minimum(zPL_First['Data'][ind],zYr['Data'][ind])
    zPL_Last['Data'][ind]=np.maximum(zPL_Last['Data'][ind],zYr['Data'][ind])
    for iT in range(tv.size):
        ind=np.where( (zYr['Data']==tv[iT]) & (np.isin(zFSC['Data'],meta['Param']['BE']['FSC']['NO List ID'])==True) )
        idxPL[tv[iT]]=( np.append(idxPL[tv[iT]][0],ind[0]),np.append(idxPL[tv[iT]][1],ind[1]) )

for iT in range(tv.size):
    print(tv[iT])

    zPL=zRef.copy()
    zPL['Data']=np.zeros(zRef['Data'].shape,dtype='int16')
    zPL['Data'][ (idxPL[tv[iT]][0].astype(int),idxPL[tv[iT]][1].astype(int)) ]=1

    # Determine what is planting vs. fill-planting
    ind_First=np.where( (zPL['Data']==1) & (zPL_First['Data']==tv[iT]) | (zPL['Data']==1) & (tv[iT]-zPL_First['Data']>=th_Fill) )
    zPL_First1=zRef.copy()
    zPL_First1['Data']=np.zeros(zRef['Data'].shape,dtype='int16')
    zPL_First1['Data'][ind_First]=1

    ind_Fill=np.where( (zPL['Data']==1) & (tv[iT]-zPL_First['Data']>0) & (tv[iT]-zPL_First['Data']<th_Fill) )
    zPL_Fill1=zRef.copy()
    zPL_Fill1['Data']=np.zeros(zRef['Data'].shape,dtype='int16')
    zPL_Fill1['Data'][ind_Fill]=1
    # plt.close('all'); plt.hist(tv[iT]-zPL_First['Data'][ind],np.arange(0,51,1))
    #plt.close('all'); plt.hist(zPL_First['Data'][ind],np.arange(1950,2024,1))

    # Time since harvest
    tsh=tv[iT]-zD['harv_yr_con1']['Data']

    # Time since knockdown
    tsk=tv[iT]-zD['kd_yr']['Data']

    # Time since fire
    tsf=tv[iT]-zD['fire_yr']['Data']

    # Time since beetle
    tsb=tv[iT]-zD['ibm_yr']['Data']

    zPE=zRef.copy()
    zPE['Data']=np.zeros(zRef['Data'].shape,dtype='int16')
    # Salvage
    ind=np.where( (zPL_First1['Data']==1) & (tsh<yr_th) & (tsh<tsf) & (tsh<tsb) )
    zPE['Data'][ind]=meta['LUT']['Derived']['RegenTypeNO']['Salvage']
    # Knockdown
    ind=np.where( (zPL_First1['Data']==1) & (tsk<yr_th) & (tsh>yr_th) & (tsk<tsf) & (tsk<tsb) )
    zPE['Data'][ind]=meta['LUT']['Derived']['RegenTypeNO']['Knockdown']
    # Straight planting fire
    ind=np.where( (zPL_First1['Data']==1) & (tsh>2) & (tsf<yr_th) & (tsf<tsh) & (tsf<tsb) )
    zPE['Data'][ind]=meta['LUT']['Derived']['RegenTypeNO']['Straight Fire']
    # Straight planting beetle
    ind=np.where( (zPL_First1['Data']==1) & (tsh>2) & (tsb<yr_th) & (tsb<tsh) & (tsb<tsf) )
    zPE['Data'][ind]=meta['LUT']['Derived']['RegenTypeNO']['Straight Insect']
    # Unknown
    ind=np.where( (zPL_First1['Data']==1) & (zPE['Data']==0) )
    zPE['Data'][ind]=meta['LUT']['Derived']['RegenTypeNO']['Unknown']
    # NSR backlock
    ind=np.where( (zPE['Data']==meta['LUT']['Derived']['RegenTypeNO']['Unknown']) & (zD['harv_yr_con1']['Data']>0) & (zD['harv_yr_con1']['Data']<=1987) | \
                (zPL_Fill1['Data']==1) & (zD['harv_yr_con1']['Data']>0) & (zD['harv_yr_con1']['Data']<=1987) )
    zPE['Data'][ind]=meta['LUT']['Derived']['RegenTypeNO']['NSR Backlog']
    # Road
    ind=np.where( (zPE['Data']==meta['LUT']['Derived']['RegenTypeNO']['Unknown']) )
    zPE['Data'][ind]=meta['LUT']['Derived']['RegenTypeNO']['Road Planting']
    # Fill-planting
    ind=np.where( (zPL_Fill1['Data']==1) & (zPE['Data']!=meta['LUT']['Derived']['RegenTypeNO']['NSR Backlog']) )
    zPE['Data'][ind]=meta['LUT']['Derived']['RegenTypeNO']['Fill Planting']

    # ind_All=np.where(zPE['Data']>0)
    # d=meta['LUT']['Derived']['RegenTypeNO'].copy()
    # for k in meta['LUT']['Derived']['RegenTypeNO'].keys():
    #     ind=np.where(zPE['Data']==meta['LUT']['Derived']['RegenTypeNO'][k])
    #     d[k]=ind[0].size/ind_All[0].size
    # d
    gis.SaveGeoTiff(zPE,meta['Paths']['bc1ha'] + '\\Management\\PL_NonOb_Type_' + str(tv[iT]) + '.tif')


z=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Management\\PL_NonOb_Type_2022.tif')
d=meta['LUT']['Derived']['RegenTypeNO'].copy()
for k in meta['LUT']['Derived']['RegenTypeNO'].keys():
    ind=np.where(z['Data']==meta['LUT']['Derived']['RegenTypeNO'][k])
    d[k]=ind[0].size

#%% Rasterize consolidated cutblocks

lNam='VEG_CONSOLIDATED_CUT_BLOCKS_SP'
vNam='HARVEST_YEAR'

if os.path.exists(meta['Paths']['bc1ha'] + '\\' + lNam)==False:
    os.mkdir(meta['Paths']['bc1ha'] + '\\' + lNam)

ind=np.where( (meta['Geos']['Variable Info']['Layer Name']==lNam) & (meta['Geos']['Variable Info']['Variable Name']==vNam) )[0]

pthin=meta['Paths']['GDB'][ meta['Geos']['Variable Info']['Geodatabase Name'][ind[0]] ]

df=gpd.read_file(pthin,layer=lNam)
df=df[df.geometry!=None]
df=df.reset_index()

zYearLast=zRef.copy()
zYearLast['Data']=np.zeros(zRef['Data'].shape,dtype='int16')

uYear=df[vNam].unique()
tv=np.arange(np.min(uYear),np.max(uYear),1)

for iT in range(tv.size):

    df0=df[df[vNam]==tv[iT]].copy()
    shapes=((geom,value) for geom, value in zip(df0.geometry,df0[vNam]))

    z0=np.zeros(zRef['Data'].shape,dtype=float)
    if len(df0)>0:
        burned=features.rasterize(shapes=shapes,fill=0,out=z0,transform=zRef['Transform'])

    z1=zRef.copy()
    z1['Data']=z0.astype(meta['Geos']['Variable Info']['Precision'][ind[0]])
    gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\' + lNam + '\\' + vNam + '_' + str(tv[iT]) + '.tif')

    # Update by year grid
    zYearLast['Data'][burned>0]=tv[iT]

# Year of last occurrence
z1=zRef.copy()
z1['Data']=zYearLast['Data'].astype('int16')
gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\' + lNam + '\\' + vNam + '_YearLast.tif')

# Mask of occurrence
z1=zRef.copy()
z1['Data']=np.zeros(zRef['Data'].shape,dtype='int8')
ind0=np.where(zYearLast['Data']>0)
z1['Data'][ind0]=1
gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\' + lNam + '\\' + vNam + '_MaskAll.tif')

# Pack into smaller number of layers

# Initialize rasters
N_Year=6
z={'Year':{}}
for iY in range(N_Year):
    z['Year'][iY+1]=np.zeros(zRef['Data'].shape,dtype='int16')

for iT in range(tv.size):
    print(tv[iT])
    z0=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\' + lNam + '\\' + vNam + '_' + str(tv[iT]) + '.tif')['Data']
    ind=np.where((z['Year'][1]==0) & (z0!=0))
    z['Year'][1][ind]=tv[iT];
    ind=np.where((z['Year'][1]!=0) & (z['Year'][1]!=tv[iT]) & (z['Year'][2]==0) & (z0!=0))
    z['Year'][2][ind]=tv[iT];
    ind=np.where((z['Year'][2]!=0) & (z['Year'][2]!=tv[iT]) & (z['Year'][3]==0) & (z0!=0))
    z['Year'][3][ind]=tv[iT];
    ind=np.where((z['Year'][3]!=0) & (z['Year'][3]!=tv[iT]) & (z['Year'][4]==0) & (z0!=0))
    z['Year'][4][ind]=tv[iT];
    ind=np.where((z['Year'][4]!=0) & (z['Year'][4]!=tv[iT]) & (z['Year'][5]==0) & (z0!=0))
    z['Year'][5][ind]=tv[iT];
    ind=np.where((z['Year'][5]!=0) & (z['Year'][5]!=tv[iT]) & (z0!=0))
    z['Year'][6][ind]=tv[iT];

for iY in range(N_Year):
    z1=zRef.copy()
    z1['Data']=z['Year'][iY+1].astype('int16')
    gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\' + lNam + '\\' + vNam + '_' + str(iY+1) + '_Year.tif')

# Plot time series to confirm it worked
lNam='VEG_CONSOLIDATED_CUT_BLOCKS_SP'
vNam='HARVEST_YEAR'
nPack=6
tv,N=u1ha.TimeSeriesOccurrenceFromPackedEvents(meta,lNam,vNam,nPack)
plt.plot(tv,N,'-bo')

#%% Rasterize harvest from RESULTS (mask all)

# Import opening layer
gdf_bm=gpd.read_file(r'C:\Users\rhember\Documents\Data\Basemaps\Basemaps.gdb',layer='NRC_POLITICAL_BOUNDARIES_1M_SP')
metaOP={}
metaOP['Path']=r'C:\Users\rhember\Documents\Data\ForestInventory\Results\20220422\Results.gdb'
metaOP['Layer']='RSLT_OPENING_SVW'; # fiona.listlayers(op['Path'])
metaOP['crs']=gdf_meta['Geos']['crs']
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
metaOP={}
metaOP['Path']=r'C:\Users\rhember\Documents\Data\ForestInventory\Results\20220422\Results.gdb'
metaOP['Layer']='RSLT_OPENING_SVW'; # fiona.listlayers(op['Path'])
metaOP['crs']=gdf_meta['Geos']['crs']
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

#%% Harvest early reconstruction

zLCC1=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCoverClass1.tif')
zAge=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\PROJ_AGE_1.tif')
zH_CC=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VEG_CONSOLIDATED_CUT_BLOCKS_SP\\HARVEST_YEAR_YearLast.tif')
zH_NTEM=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Disturbances\\Harvest_NTEM_Year.tif')

z0=np.zeros(zRef['Data'].shape,dtype='int16')
ind=np.where( (zRef['Data']==1) & (zH_CC['Data']>0) ); z0[ind]=1
ind=np.where( (zRef['Data']==1) & (z0==0) & (zH_NTEM['Data']>0) ); z0[ind]=2
ind=np.where( (zRef['Data']==1) & (z0==0) & (zAge['Data']>0) & (zAge['Data']<150) & (zLCC1['Data']==meta['LUT']['Derived']['lcc1']['Forest']) ); z0[ind]=3
ind=np.where( (zRef['Data']==1) & (z0==0) & (zAge['Data']>=150) & (zLCC1['Data']==meta['LUT']['Derived']['lcc1']['Forest']) ); z0[ind]=4
ind=np.where( (zRef['Data']==1) & (z0==0) & (zAge['Data']==0) & (zLCC1['Data']==meta['LUT']['Derived']['lcc1']['Forest']) ); z0[ind]=4
plt.close('all')
plt.matshow(z0)

# Create rings from populated locations
df=gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='NRC_POPULATED_PLACES_1M_SP')
df=df[df.geometry!=None]
zD=np.zeros(zRef['Data'].shape,dtype='int16')
#binD=np.arange(1,110,1)
binD=np.array([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.25,2.5,2.75,3,3.25,3.5,3.75,4,4.5,5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10,10.5,11,11.5,12,12.5,13,13.5,14,14.5,15,16,17,18,19,20,22,24,26,28,30,40,50,60,120])
for iD in range(binD.size):
    print(binD[iD])
    df0=df.copy()
    df0['geometry']=df0.geometry.buffer(1000*binD[iD])
    z1=np.zeros(zRef['Data'].shape,dtype=float)
    shapes=((geom,value) for geom, value in zip(df0['geometry'],df0['NRC_PP1M_SYSID']))
    burned=features.rasterize(shapes=shapes,fill=0,out=z1,transform=zRef['Transform'])
    zD[(burned>0) & (zD==0)]=binD[iD]

# Year for each ring
yr=np.linspace(1855,1960,binD.size).astype('int16')

zH=zRef.copy()
zH['Data']=np.zeros(zRef['Data'].shape,dtype='int16')
for iD in range(binD.size):
    ind=np.where( (z0==3) & (zD==binD[iD]) )
    zH['Data'][ind]=yr[iD]
plt.close('all'); plt.matshow(zH['Data'],clim=[1855,1955])

gis.SaveGeoTiff(zH,meta['Paths']['bc1ha'] + '\\Disturbances\\Harvest_Early_Reconstruction_Year.tif')

#%% Harvest consolidated 1 (w/o early reconstruction)

zLCC1=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCoverClass1.tif')
zH_CC=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VEG_CONSOLIDATED_CUT_BLOCKS_SP\\HARVEST_YEAR_YearLast.tif')
zH_NTEM=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Disturbances\\Harvest_NTEM_Year.tif')

z=zRef.copy()
z['Data']=np.zeros(zRef['Data'].shape,dtype='int16')
ind=np.where( (zLCC1['Data']==meta['LUT']['Derived']['lcc1']['Forest']) & (zH_CC['Data']>0) ); z['Data'][ind]=zH_CC['Data'][ind]
ind=np.where( (zLCC1['Data']==meta['LUT']['Derived']['lcc1']['Forest']) & (zH_NTEM['Data']>0) ); z['Data'][ind]=zH_NTEM['Data'][ind]
gis.SaveGeoTiff(z,meta['Paths']['bc1ha'] + '\\Disturbances\\Harvest_Consol1_Year.tif')

#%% Harvest consolidated 2 (with early reconstruction

zH_CC=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VEG_CONSOLIDATED_CUT_BLOCKS_SP\\HARVEST_YEAR_YearLast.tif')
zH_NTEM=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Disturbances\\Harvest_NTEM_Year.tif')
zH_Early=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Disturbances\\Harvest_Early_Reconstruction_Year.tif')

z=zRef.copy()
z['Data']=np.zeros(zRef['Data'].shape,dtype='int16')
ind=np.where( (zLCC1['Data']==meta['LUT']['Derived']['lcc1']['Forest']) & (zH_Early['Data']>0) ); z['Data'][ind]=zH_Early['Data'][ind]
ind=np.where( (zLCC1['Data']==meta['LUT']['Derived']['lcc1']['Forest']) & (zH_CC['Data']>0) ); z['Data'][ind]=zH_CC['Data'][ind]
ind=np.where( (zLCC1['Data']==meta['LUT']['Derived']['lcc1']['Forest']) & (zH_NTEM['Data']>0) ); z['Data'][ind]=zH_NTEM['Data'][ind]
gis.SaveGeoTiff(z,meta['Paths']['bc1ha'] + '\\Disturbances\\Harvest_Consol2_Year.tif')

#%% Percent dead from cruises

dC=gu.ipickle(r'C:\Users\rhember\Documents\Data\ECAS\Received 2023-04-04\rh_COMP_DATA_CleanCombined.pkl')

vTM=list(meta['LUT']['FTEN_CUT_BLOCK_POLY_SVW']['TIMBER_MARK'].values())
kTM=list(meta['LUT']['FTEN_CUT_BLOCK_POLY_SVW']['TIMBER_MARK'].keys())

zTM=gis.OpenGeoTiff(meta['Paths'] + '\\FTEN_CUT_BLOCK_POLY_SVW\\TIMBER_MARK.tif')

u=np.unique(dC['PRIMARY_MARK'])
iZ=gu.IndicesFromUniqueArrayValues(zTM['Data'].flatten())

z1=zRef.copy()
z1['Data']=np.zeros(zRef['Data'].shape,dtype='int16').flatten()
for k in iZ.keys():
    ind1=np.where(vTM==k)[0]
    if ind1.size==0:
        continue
    ind2=np.where( (dC['PRIMARY_MARK']==kTM[ind1[0]]) )[0]
    if ind2.size>0:
        z1['Data'][iZ[k]]=np.nanmean(d['Pct Dead Net'][ind2])
        #print(np.nanmean(d['Pct Dead Net'][ind2]))

z1['Data']=np.reshape(z1['Data'],zRef['Data'].shape)
#plt.matshow(z1['Data'],clim=[0,100])
#plt.hist(z1['Data'][0::20,0::20].flatten(),np.arange(0,105,5))
gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\Disturbances\\Harvest_PctDead_FromCruise.tif')

#%% Salvage logging mask

zH=gis.OpenGeoTiff(meta['Paths'] + '\\Disturbances\\VEG_CONSOLIDATED_CUT_BLOCKS_SP_MaskAll.tif')
zLCC1=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCoverClass1.tif')
zPD=gis.OpenGeoTiff(meta['Paths'] + '\\Disturbances\\HarvestPercentDead_FromCruise.tif')

z=zRef.copy()
z['Data']=6*np.ones(zRef['Data'].shape,dtype='int8')
ind=np.where( (zPD['Data']>=50) ); z['Data'][ind]=1 # Salvage high
ind=np.where( (zPD['Data']>=10) & (zPD['Data']<50) ); z['Data'][ind]=2 # Salvage low
ind=np.where( (zPD['Data']<10) & (zH['Data']>0) ); z['Data'][ind]=3 # Non-salvage harvested forest
ind=np.where( (zH['Data']==0) & (zLCC1['Data']==1) ); z['Data'][ind]=4 # Unharvested forest
ind=np.where( (zRef['Data']==1) & (zLCC1['Data']!=1) ); z['Data'][ind]=5 # Non-forest land
ind=np.where( (zRef['Data']==0) ); z['Data'][ind]=6 # Non land
gis.SaveGeoTiff(z,meta['Paths'] + '\\Disturbances\\HarvestSalvageMask_FromCruise.tif')

#%% Create range tenure consolidated

zLCC1=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCoverClass1.tif')
zR=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\rangeten.tif')

d={1:'Forest with grazing tenure',2:'Forest with haycutting tenure',3:'Forest with no range tenure',4:'Non-forest land',5:'Non land'}

z1=zRef.copy()
z1['Data']=np.zeros(zRef['Data'].shape,dtype='int8')
ind=np.where( (zLCC1['Data']==1) & (zR['Data']==1) | (zLCC1['Data']==1) & (zR['Data']==2) | (zLCC1['Data']==1) & (zR['Data']==3) ); z1['Data'][ind]=1
ind=np.where( (zLCC1['Data']==1) & (zR['Data']==4) | (zLCC1['Data']==1) & (zR['Data']==5) | (zLCC1['Data']==1) & (zR['Data']==6) ); z1['Data'][ind]=2
ind=np.where( (zLCC1['Data']==1) & (zR['Data']==0) ); z1['Data'][ind]=3
ind=np.where( (zLCC1['Data']!=1) ); z1['Data'][ind]=4
ind=np.where( (zRef['Data']!=1) ); z1['Data'][ind]=5
gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\LandCoverUse\\rangeten_consol.tif')

#%% Create Crown land mask

z0=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\f_own.tif')
z1=zRef.copy()
z1['Data']=np.zeros(zRef['Data'].shape,dtype='int8')
ind=np.where( (zRef['Data']==1) & (z0['Data']>=9) )
z1['Data'][ind]=1
# plt.close('all'); plt.matshow(z1['Data'])# Confirm that it worked
gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\LandCoverUse\\CrownForestMask.tif')

#%% Create BGC Zone / NDT Combination

zBGC=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\BEC_ZONE_CODE.tif')
zNDT=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\ndt.tif')

n=np.zeros(5)
for i in range(1,5):
    ind=np.where(zNDT['Data']==i)
    n[i]=ind[0].size/1e6

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
df.to_excel(meta['Paths']['bc1ha'] + '\\VRI 2023\\lut_bgcz_ndt_combo.xlsx')

z1=zRef.copy()
z1['Data']=z
gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\VRI 2023\\bgcz_ndt_combo.tif')

#%% Create Tree Density Class
# *** Ensure harvested areas are listed as dense ***

z=u1ha.Import_Raster(meta,[],['lcc1_c','lc5','harv_yr_con1'])
#zLC=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCoverClass1_Current.tif')
#zLC5=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\lc5.tif')
#zH=gis.OpenGeoTiff(Disturbances\Harvest_Consol1_Year.tif')

z1=zRef.copy()
z1['Data']=np.zeros(z['lc5']['Data'].shape,dtype='int8')
ind=np.where( (z['lcc1_c']['Data']==meta['LUT']['Derived']['lcc1']['Forest']) & (z['lc5']['Data']==meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_5']['SP']) ); z1['Data'][ind]=meta['LUT']['Derived']['tdc']['Sparse']
ind=np.where( (z['lcc1_c']['Data']==meta['LUT']['Derived']['lcc1']['Forest']) & (z['lc5']['Data']==meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_5']['OP']) ); z1['Data'][ind]=meta['LUT']['Derived']['tdc']['Open']
ind=np.where( (z['lcc1_c']['Data']==meta['LUT']['Derived']['lcc1']['Forest']) & (z['lc5']['Data']==meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_5']['DE']) ); z1['Data'][ind]=meta['LUT']['Derived']['tdc']['Dense']
ind=np.where( (z['harv_yr_con1']['Data']==1) ); z1['Data'][ind]=meta['LUT']['Derived']['tdc']['Dense']
ind=np.where( (z['lcc1_c']['Data']==meta['LUT']['Derived']['lcc1']['Forest']) & (z1['Data']==0) ); z1['Data'][ind]=meta['LUT']['Derived']['tdc']['Open']
# plt.matshow(z1,clim=[0,3])
gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\LandCoverUse\\TreeDensityClass_Current.tif')

# # With shrubs and grasses
# z=np.zeros(zLC4['Data'].shape,dtype='int8')
# ind=np.where( (z['lcc1_c']['Data']==lut['lc4']['Shrub Low']) | (zLC4['Data']==lut['lc4']['Shrub Tall']) ); z[ind]=4
# ind=np.where( (z['lcc1_c']['Data']==lut['lc4']['Herb Gramanoid']) | (zLC4['Data']==lut['lc4']['Herb Forbs']) ); z[ind]=5
# ind=np.where( (z['lcc1_c']['Data']==1) & (z['lc5']['Data']==lut['lc5']['Sparse']) ); z[ind]=1
# ind=np.where( (z['lcc1_c']['Data']==1) & (z['lc5']['Data']==lut['lc5']['Open']) ); z[ind]=2
# ind=np.where( (z['lcc1_c']['Data']==1) & (z['lc5']['Data']==lut['lc5']['Dense']) ); z[ind]=3
# ind=np.where( (z['harv_yr_con1']['Data']==1) ); z[ind]=meta['LUT']['Derived']['tdc']['Dense']
# ind=np.where( (zRef['Data']==1) & (z==0) ); z[ind]=6
# ind=np.where( (zRef['Data']==0) ); z[ind]=7
# plt.matshow(z[0::3,0::3])

# z1=zRef.copy()
# z1['Data']=z
# gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\LandCoverUse\\TreeDensityClass_WithShrubsGrasses.tif')

#%% Protected lands consolidated

zLCC1=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCoverClass1.tif')
zProt=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\TA_PROTECTED_LANDS_SV\\PROTECTED_LANDS_DESIGNATION.tif')
zPark=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\TA_PARK_ECORES_PA_SVW\\ADMIN_AREA_SID.tif')
zOGMA=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RMP_OGMA_LEGAL_ALL_SVW\\LEGAL_OGMA_PROVID.tif')
zOGD=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\OGSR_TAP_PRIORITY_DEF_AREA_SP\\PRIORITY_DEFERRAL_ID.tif')
zNP=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\CLAB_NATIONAL_PARKS\\ENGLISH_NAME.tif')
#zFCR=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\RSLT_FOREST_COVER_RESERVE_SVW.tif')

# ListCon=['No Timber Harvesting Areas','Landscape Corridors',
#   'Critical Deer Winter Range','Sensitive Watershed','Water Management Units','High Value Wetlands for Moose','Telkwa Caribou Recovery Area',
#   'Caribou Migration Corridor','High Biodiversity Emphasis Areas','Scenic Corridors']

# Generate random new protected areas
flg=0
if flg==1:
    gdf=u1ha.Import_GDBs_ProvinceWide()
    N=2500 # Number within bounding box
    Dbuf=6200 # metres

    x=np.random.uniform(zRef['xmin'],zRef['xmax'],size=N)
    y=np.random.uniform(zRef['ymin'],zRef['ymax'],size=N)
    points=[]
    for k in range(x.size):
        points.append(Point(x[k],y[k]))
    gdf_xy=gpd.GeoDataFrame({'geometry':points,'ID':1})
    gdf_xy.crs=gdf['bc_bound']['gdf'].crs

    gdf_xy=gpd.sjoin(gdf_xy,gdf['bc_bound']['gdf'],op='within')
    #gdf_xy.plot(ax=ax[0],markersize=8)

    gdf_xyb=gdf_xy.geometry.buffer(Dbuf)
    gdf_xyb=gpd.GeoDataFrame({'geometry':gdf_xyb,'ID':np.arange(0,gdf_xyb.size,1)})
    gdf_xyb.crs=gdf['bc_bound']['gdf'].crs
    #gdf_xyb.plot(ax=ax[0],color='r')

    shapes=((geom,value) for geom, value in zip(gdf_xyb['geometry'],gdf_xyb['ID']))
    z=np.zeros(zRef['Data'].shape,dtype=float)
    burned=features.rasterize(shapes=shapes,fill=0,out=z,transform=zRef['Transform'])

    zRn=zRef.copy()
    zRn['Data']=z.astype('int32')
    gis.SaveGeoTiff(zRn,meta['Paths']['bc1ha'] + '\\LandCoverUse\\REARs_RandomAreas.tif')

else:
    zRn=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\REARs_RandomAreas.tif')

# # Protected areas with random areas

# Area treed
ind0=np.where( (zLCC1['Data']==meta['LUT']['Derived']['lcc1']['Forest']) );
A_treed=ind0[0].size

# Everything completed+proposed without random additions (Comp=1, Prop=2)
zCP=zRef.copy()
zCP['Data']=np.zeros(zRef['Data'].shape,dtype='int16')
ind=np.where( (zLCC1['Data']==meta['LUT']['Derived']['lcc1']['Forest']) & (zPark['Data']>0) |
              (zLCC1['Data']==meta['LUT']['Derived']['lcc1']['Forest']) & (zOGMA['Data']>0) |
              (zLCC1['Data']==meta['LUT']['Derived']['lcc1']['Forest']) & (zNP['Data']>0) |
              (zLCC1['Data']==meta['LUT']['Derived']['lcc1']['Forest']) & (zProt['Data']>0))
zCP['Data'][ind]=1
ind=np.where( (zCP['Data']==0) & (zLCC1['Data']==meta['LUT']['Derived']['lcc1']['Forest']) & (zOGD['Data']>0) )
zCP['Data'][ind]=2
ind=np.where( (zLCC1['Data']==meta['LUT']['Derived']['lcc1']['Forest']) & (zCP['Data']==0) ); zCP['Data'][ind]=3
ind=np.where( (zLCC1['Data']!=meta['LUT']['Derived']['lcc1']['Forest']) & (zCP['Data']==0) ); zCP['Data'][ind]=4
ind=np.where( (zRef['Data']==0) ); zCP['Data'][ind]=5

ind=np.where( (zLCC1['Data']==meta['LUT']['Derived']['lcc1']['Forest']) & (zCP['Data']>=1) & (zCP['Data']<=2) ); A_CPsep=ind[0].size
print('Area Protected from completed+proposed (Sep) (%): ' + str(A_CPsep/A_treed*100))

zCP['Data']=zCP['Data'].astype('int16')
gis.SaveGeoTiff(zCP,meta['Paths']['bc1ha'] + '\\LandCoverUse\\REARs_CompPlusProp.tif')

# Everything completed+proposed with random additions (Comp=1, Prop=2)
zCP=zRef.copy()
zCP['Data']=np.zeros(zRef['Data'].shape,dtype='int16')
ind=np.where( (zLCC1['Data']==meta['LUT']['Derived']['lcc1']['Forest']) & (zPark['Data']>0) |
              (zLCC1['Data']==meta['LUT']['Derived']['lcc1']['Forest']) & (zOGMA['Data']>0) |
              (zLCC1['Data']==meta['LUT']['Derived']['lcc1']['Forest']) & (zNP['Data']>0) |
              (zLCC1['Data']==meta['LUT']['Derived']['lcc1']['Forest']) & (zProt['Data']>0) )
zCP['Data'][ind]=1
ind=np.where( (zCP['Data']==0) & (zLCC1['Data']==meta['LUT']['Derived']['lcc1']['Forest']) & (zOGD['Data']>0) | (zCP['Data']==0) & (zLCC1['Data']==meta['LUT']['Derived']['lcc1']['Forest']) & (zRn['Data']>0) )
zCP['Data'][ind]=2
ind=np.where( (zLCC1['Data']==meta['LUT']['Derived']['lcc1']['Forest']) & (zCP['Data']==0) ); zCP['Data'][ind]=3
ind=np.where( (zLCC1['Data']!=meta['LUT']['Derived']['lcc1']['Forest']) & (zCP['Data']==0) ); zCP['Data'][ind]=4
ind=np.where( (zRef['Data']==0) ); zCP['Data'][ind]=5

ind=np.where( (zLCC1['Data']==meta['LUT']['Derived']['lcc1']['Forest']) & (zCP['Data']>=1) & (zCP['Data']<=2) ); A_CPsep=ind[0].size
print('Area Protected from completed+proposed (Sep) (%): ' + str(A_CPsep/A_treed*100))

zCP['Data']=zCP['Data'].astype('int32')
gis.SaveGeoTiff(zCP,meta['Paths']['bc1ha'] + '\\LandCoverUse\\REARs_CompPlusPropWithRandomAdditions.tif')

#%% Rasterize early land use change year

zLCC1=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCoverClass1.tif')

# Create rings from populated locations
df=gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='NRC_POPULATED_PLACES_1M_SP')
df=df[df.geometry!=None]
zD=np.zeros(zRef['Data'].shape,dtype='int16')
binD=np.array([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.25,2.5,2.75,3,3.25,3.5,3.75,4,4.5,5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10,11,12,13,14,15,20,25,30,60,120,240,480])
for iD in range(binD.size):
    print(binD[iD])
    df0=df.copy()
    df0['geometry']=df0.geometry.buffer(1000*binD[iD])
    z0=np.zeros(zRef['Data'].shape,dtype=float)
    shapes=((geom,value) for geom, value in zip(df0['geometry'],df0['NRC_PP1M_SYSID']))
    burned=features.rasterize(shapes=shapes,fill=0,out=z0,transform=zRef['Transform'])
    zD[(burned>0) & (zD==0)]=binD[iD]

# Check that the rings are big enough
ind=np.where( (zLCC1['Data']==meta['LUT']['Derived']['lcc1']['Settlement']) & (zD==0) | \
             (zLCC1['Data']==meta['LUT']['Derived']['lcc1']['Herb']) & (zD==0))
print(ind[0].size)

# Year for each ring
yr=np.linspace(1855,2010,binD.size).astype('int16')

zLUC=zRef.copy()
zLUC['Data']=np.zeros(zRef['Data'].shape,dtype='int16')
for iD in range(binD.size):
    ind=np.where( (zLCC1['Data']==meta['LUT']['Derived']['lcc1']['Settlement']) & (zD==binD[iD]) | \
                 (zLCC1['Data']==meta['LUT']['Derived']['lcc1']['Herb']) & (zD==binD[iD]))
    zLUC['Data'][ind]=yr[iD]
#plt.matshow(zLUC['Data'],clim=[1845,2020])

gis.SaveGeoTiff(zLUC,meta['Paths']['bc1ha'] + '\\Disturbances\\LandUseChange1_Year.tif')

#z=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\Terrain\DistanceFromRoads.tif')

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
gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\Terrain\DistanceFromRoads.tif')

#z=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\Terrain\DistanceFromRoads.tif')
#plt.matshow(z['Data'])

#%% Rasterize distance from timber facilities

#fiona.listlayers(r'C:\Users\rhember\Documents\Data\Geodatabases\LandCover\20230607\LandCover.gdb')
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

#%% Rasterize forest cover

gdf_FC={}
gdf_FC['Path']=r'C:\Users\rhember\Documents\Data\ForestInventory\Results\20230430\Results.gdb'
gdf_FC['Layer']='RSLT_FOREST_COVER_INV_SVW'; # fiona.listlayers(op['Path'])
gdf_FC['crs']=gdf_meta['Geos']['crs']
gdf_FC['Keep Geom']='Off'
gdf_FC['Select Openings']=np.array([])
gdf_FC['SBC']=np.array([])
gdf_FC['FSC']=np.array([])
gdf_FC['ROI']=[]
gdf_FC['gdf']=qgdb.Query_Openings(gdf_FC,[])

gdf_FCs={}
gdf_FCs['Path']=r'C:\Users\rhember\Documents\Data\ForestInventory\Results\20230430\Results.gdb'
gdf_FCs['Layer']='RSLT_FOREST_COVER_INV_SVW'; # fiona.listlayers(op['Path'])
gdf_FCs['crs']=gdf_meta['Geos']['crs']
gdf_FCs['Keep Geom']='On'
gdf_FCs['Select Openings']=np.array([])
gdf_FCs['SBC']=np.array([])
gdf_FCs['FSC']=np.array([])
gdf_FCs['SOC1']=np.array([])
gdf_FCs['ROI']=[]
gdf_FCs['gdf']=qgdb.Query_Openings(gdf_FCs,[])

dMisFC=gu.ipickle(r'C:\Users\rhember\Documents\Data\ForestInventory\Results\20230430\missing_geo_fc_list.pkl')

# Open the shapefile
df=gpd.read_file(pthin,layer=lnam)

# Remove features with no geometry
df=df[df.geometry!=None]
df=df.reset_index()

#%% Consolidated retention layer

d=gu.ReadExcel(meta['Paths']['bc1ha'] + '\\Results\\SILV_RESERVE_CODE.xlsx')
zLCC1=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCoverClass1.tif')
zH=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Disturbances\\VEG_CONSOLIDATED_CUT_BLOCKS_SP_MaskAll.tif')
zFCR=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Results\\SILV_RESERVE_CODE.tif')
zIBR=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Results\\SILV_RESERVE_CODE_InBlock.tif')
zC=np.maximum(zFCR['Data'],zIBR['Data'])
u,N=gu.CountByCategories(zC,'Percent')

flg=0
if flg==1:
    # Look at overlap
    ind1=np.where(zFCR['Data']>0)
    ind2=np.where(zIBR['Data']>0)
    ind3=np.where( (zFCR['Data']>0) & (zIBR['Data']>0) )
    ind4=np.where( (zFCR['Data']>0) & (zIBR['Data']==0) )
    ind5=np.where( (zFCR['Data']==0) & (zIBR['Data']>0) )

    print(ind1[0].size)
    print(ind2[0].size)
    print(ind3[0].size)
    print(ind4[0].size)
    print(ind5[0].size)

d1={1:['Dispersed'],2:'Group',3:'Riparian',4:'Wildlife trees',5:'Other',6:'No reserves',7:'Forest',8:'Non-forest'}

z1=zRef.copy()
z1['Data']=np.zeros(zRef['Data'].shape,dtype='int8')
ind=np.where( (zC>1) ); z1['Data'][ind]=5
ind=np.where( (zC==7) | (zC==4) ); z1['Data'][ind]=1
ind=np.where( (zC==8) ); z1['Data'][ind]=2
ind=np.where( (zC==9) ); z1['Data'][ind]=3
ind=np.where( (zC==2) ); z1['Data'][ind]=4
ind=np.where( (z1['Data']==0) & (zH['Data']>0) ); z1['Data'][ind]=6
ind=np.where( (z1['Data']==0) & (zH['Data']==0) & (zLCC1['Data']==meta['LUT']['Derived']['lcc1']['Forest']) ); z1['Data'][ind]=7
ind=np.where( (z1['Data']==0) & (zH['Data']==0) & (zLCC1['Data']!=meta['LUT']['Derived']['lcc1']['Forest']) ); z1['Data'][ind]=8
gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\Results\\SILV_RESERVE_CODE_Consolidated.tif')

#%% Harvest regeneration type

zH=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Disturbances\\VEG_CONSOLIDATED_CUT_BLOCKS_SP_MaskAll.tif')
zFCR=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Results\\RSLT_FOREST_COVER_RESERVE_SVW.tif')
zPL=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Results\\Planting_FromRESULTS_MaskCount.tif')

tv=np.arange(1960,2022,1)
HY=0*zRef['Data'].astype('int16')
for iT in range(tv.size):
    zH0=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Disturbances\\VEG_CONSOLIDATED_CUT_BLOCKS_SP_' + str(tv[iT]) + '.tif')
    ind=np.where( (zH0['Data']>0) & (HY==0) )
    HY[ind]=tv[iT]

z=zRef.copy()
z['Data']=np.zeros(zRef['Data'].shape,dtype='int8')
# Harvested and planted
ind=np.where( (zH['Data']>0) & (zFCR['Data']==0) & (zPL['Data']>0) )
z['Data'][ind]=1
# Harvested but not planted (before FRPA)
ind=np.where( (zH['Data']>0) & (HY<=1987) & (zFCR['Data']==0) & (zPL['Data']==0) )
z['Data'][ind]=2
# Harvested but not planted (after FRPA)
ind=np.where( (zH['Data']>0) & (HY>1987) & (HY<=2018) & (zFCR['Data']==0) & (zPL['Data']==0) )
z['Data'][ind]=3
# Harvested but not planted (after FRPA)
ind=np.where( (zH['Data']>0) & (HY>2018) & (zFCR['Data']==0) & (zPL['Data']==0) )
z['Data'][ind]=4
# Not harvested, but planted
ind=np.where( (zH['Data']==0) & (zPL['Data']>0) | (zFCR['Data']>0) & (zPL['Data']>0) )
z['Data'][ind]=5
# Not harvested, not planted
ind=np.where( (zH['Data']==0) & (zFCR['Data']==0) & (zPL['Data']==0) )
z['Data'][ind]=6
np.unique(z['Data'])

bin=np.arange(1,6,1)
n=np.zeros(bin.size)
for i in range(bin.size):
    ind=np.where(z['Data']==bin[i])
    n[i]=ind[0].size
plt.bar(bin,n/np.sum(n)*100)

# plt.matshow(z['Data'])
gis.SaveGeoTiff(z,meta['Paths']['bc1ha'] + '\\Results\\Harvest_Regen_Type.tif')

# # Just current era = 79%
# n1=np.where( (HY>1987) & (HY<2018) & (zFCR['Data']==0) & (zPL['Data']>0) )[0].size
# n2=np.where( (HY>1987) & (HY<2018) & (zFCR['Data']==0) )[0].size
# print(n1/n2)

# Time series
yH=np.zeros(tv.size)
yHP=np.zeros(tv.size)
for iT in range(tv.size):
    zH=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Disturbances\\VEG_CONSOLIDATED_CUT_BLOCKS_SP_' + str(tv[iT]) + '.tif')
    ind=np.where( (zH['Data']>0) & (zFCR['Data']==0) )
    yH[iT]=ind[0].size
    ind=np.where( (zH['Data']>0) & (zFCR['Data']==0) & (zPL['Data']>0) )
    yHP[iT]=ind[0].size

plt.plot(tv,(yH-yHP)/yH,'ob-')

#%% Recent wildfire
# Current year + preveous year (sometimes missing)

zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])

# Current year
pthin=r'C:\Users\rhember\Documents\Data\Wildfire\Current\prot_current_fire_polys.shp'
df=gpd.read_file(pthin)
df=df[df.geometry!=None]
df=df.reset_index()
shapes=((geom,value) for geom, value in zip(df['geometry'],df['FIRE_YEAR']))
z0=np.zeros(zRef['Data'].shape,dtype=float)
burned=features.rasterize(shapes=shapes,fill=0,out=z0,transform=zRef['Transform'])

z1=zRef.copy()
z1['Data']=np.zeros(zRef['Data'].shape,dtype='int16')
ind=np.where(burned>0); z1['Data'][ind]=1
plt.close(); plt.matshow(z1['Data'])
print(np.sum(z1['Data']))
gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\prot_current_fire_polys\prot_current_fire_polys.tif')
gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\PROT_HISTORICAL_FIRE_POLYS_SP\\FIRE_YEAR_2023.tif')

# Previous year when missing
import geotable
import pyproj
srs=gis.ImportSRSs()
crs=pyproj.CRS(srs['String']['Geographic'])

a=geotable.load(r'C:\Users\rhember\Documents\Data\Wildfire\BC Fire Perimeters 2020-2022.kmz')
df=gpd.GeoDataFrame(data=a.Name,geometry=a.geometry_object)
df['ID']=np.ones(len(df))
df.crs=pyproj.CRS(srs['String']['Geographic'])
df=df.to_crs({'init':'epsg:3005'})

shapes=((geom,value) for geom, value in zip(df['geometry'],df['ID']))
z0=np.zeros(zRef['Data'].shape,dtype=float)
burned=features.rasterize(shapes=shapes,fill=0,out=z0,transform=zRef['Transform'])

# Fix by removing other years
tv=np.arange(2020,2021+1,1)
for iT in range(tv.size):
    zF0=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\PROT_HISTORICAL_FIRE_POLYS_SP\\FIRE_YEAR_' + str(tv[iT]) + '.tif')
    ind=np.where(zF0['Data']>0)
    burned[ind]=0

z1=zRef.copy()
z1['Data']=np.zeros(zRef['Data'].shape,dtype='int16')
ind=np.where(burned>0); z1['Data'][ind]=1
print(np.sum(z1['Data']))
gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\PROT_HISTORICAL_FIRE_POLYS_SP\\FIRE_YEAR_2022.tif')

#%% Ecozones of Canada

zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])

srs=gis.ImportSRSs()
crs=pyproj.CRS(srs['String']['Geographic'])

pthin=r'C:\Users\rhember\Documents\Data\Ecozones\nef_ca_ter_ecozone_v2_2.geojson'
df=gpd.read_file(pthin)
df=df[df.geometry!=None]
df=df.reset_index()
df.crs=pyproj.CRS(srs['String']['Geographic'])
df=df.to_crs({'init':'epsg:3005'})
# Used to create LUT: df.drop(columns='geometry').to_excel(r'C:\Users\rhember\Documents\Data\Ecozones\table.xlsx')

shapes=((geom,value) for geom, value in zip(df['geometry'],df['ECOZONE_ID']))
z0=np.zeros(zRef['Data'].shape,dtype=float)
burned=features.rasterize(shapes=shapes,fill=0,out=z0,transform=zRef['Transform'])
z1=zRef.copy()
z1['Data']=np.zeros(zRef['Data'].shape,dtype='int16')
ind=np.where( (burned>0) &(zRef['Data']==1) ); z1['Data'][ind]=burned[ind]
plt.close('all'); plt.matshow(z1['Data'],clim=[0,15])
gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\Ecozones_Canada\\Ecozones_Canada.tif')

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

    fout=meta['Paths']['bc1ha'] + '\\Management' + '\\btk_spray_' + str(tv[iT]) + '.tif'
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

# #%% Rasterize BCTS

# # Input path to RESULTS database (downloaded from BC data catalogue)
# pthin=r'C:\Users\rhember\Documents\Data\ForestInventory\LandUse\20220422\LandUse.gdb'

# lnam='BCTS_OPERATING_AREAS_SP'

# # Open the shapefile
# df=gpd.read_file(pthin,layer=lnam)

# # Remove features with no geometry
# df=df[df.geometry!=None]

# df0=df.copy()
# df0['dummy']=np.ones(len(df0))
# shapes=((geom,value) for geom, value in zip(df0['geometry'],df0['dummy']))

# z0=np.zeros(zRef['Data'].shape,dtype=float)
# burned=features.rasterize(shapes=shapes,fill=0,out=z0,transform=zRef['Transform'])

# zOut=zRef.copy()
# zOut['Data']=np.zeros(zRef['Data'].shape,dtype='int16')
# zOut['Data'][burned>0]=1
# gis.SaveGeoTiff(zOut,meta['Paths']['bc1ha'] + '\\LandCoverUse\\bcts_op_area.tif')

#%% Global Forest Change Loss Year (mosaic and reproject)

finL=['50N_120W','50N_130W','60N_120W','60N_130W','60N_140W']

# Loss year
z=zRef.copy()
z['Data']=np.zeros(z['Data'].shape)
for f in finL:
    fin=r'C:\Users\rhember\Documents\Data\Global Forest Change\Hansen_GFC-2021-v1.9_lossyear_' + f + '.tif'
    fout=r'C:\Users\rhember\Documents\Data\Global Forest Change\Hansen_GFC-2021-v1.9_lossyear_' + f + 'p.tif'
    gis.ReprojectRasterAndClipToRaster(fin,fout,meta['Paths']['bc1ha Ref Grid'],meta['Geos']['crs'])

    z0=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\Global Forest Change\Hansen_GFC-2021-v1.9_lossyear_' + f + 'p.tif')
    ind=np.where(z0['Data']>0)
    z['Data'][ind]=z0['Data'][ind]

ind=np.where(zRef['Data']==0)
z['Data'][ind]=0

z['Data']=z['Data'].astype('int16')
ind=np.where(z['Data']>0)
z['Data'][ind]=z['Data'][ind]+2000

gis.SaveGeoTiff(z,meta['Paths']['bc1ha'] + '\\Disturbances\\GlobalForestChange_LossYear_2021.tif')

#%% Global Forest Change Loss Year (adjusted to remove known disturbances)

zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
vList=['gfcly','harv_yr_con1','fire_yr','ibm_yr']
z=u1ha.Import_Raster(meta,[],vList)

ind=np.where( (np.abs(z['gfcly']['Data']-z['harv_yr_con1']['Data'])<3) | (np.abs(z['gfcly']['Data']-z['fire_yr']['Data'])<2) | (np.abs(z['gfcly']['Data']-z['ibm_yr']['Data'])<2) )
z['gfcly']['Data'][ind]=0
plt.matshow(z['gfcly']['Data'])

gis.SaveGeoTiff(z['gfcly'],meta['Paths']['bc1ha'] + '\\Disturbances\\GlobalForestChange_LossYear_2021_Filtered.tif')

#%% Reproject CEC Land Use Map 2020

fin=r'C:\Users\rhember\Documents\Data\Land Cover\NALCMS\NALCMS_lc_2010.tif'
fout=r'C:\Users\rhember\Documents\Data\Land Cover\NALCMS\NALCMS_lc_2010_bc1ha.tif'
gis.ReprojectRasterAndClipToRaster(fin,fout,fref,meta['Geos']['crs'])

z=gis.OpenGeoTiff(fout)
z=gis.ClipToRaster(z,zRef)

ind=np.where(zRef['Data']==0)
z['Data'][ind]=0

gis.SaveGeoTiff(z,meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCoverClass_CEC_2010.tif')

#z=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCoverClass_CEC_.tif')
#plt.matshow(z['Data'][0::5,0::5])

#z=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\Land Cover\NALCMS\NALCMS_lc_2010.tif')

# Compress categories to remove tropics
zLCC_CEC10=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCoverClass_CEC_2010.tif')
zLCC_CEC20=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCoverClass_CEC_2020.tif')
lut,zLCC_CEC10,zLCC_CEC20=u1ha.NALCMS_Compress(meta['LUT']['Derived']['lcc_cec'],zRef,zLCC_CEC10,zLCC_CEC20)
# Manually saved LUT to excell
gis.SaveGeoTiff(zLCC_CEC10,meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCoverClass_CEC_2010_Compressed.tif')
gis.SaveGeoTiff(zLCC_CEC20,meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCoverClass_CEC_2020_Compressed.tif')

#%% Reproject Harvest Year from NTEMS

fin=r'C:\Users\rhember\Documents\Data\Harvest\NTEMS\harv85to20.tif'
fout=r'C:\Users\rhember\Documents\Data\BC1ha\Disturbances\Harvest_NTEM_Year.tif'
gis.ReprojectRasterAndClipToRaster(fin,fout,fref,meta['Geos']['crs'])

z=gis.OpenGeoTiff(fout)
z=gis.ClipToRaster(z,zRef)

ind=np.where(zRef['Data']==0)
z['Data'][ind]=0

gis.SaveGeoTiff(z,fout)

#%% Reproject Land Cover Class from NTEMS

fin=r'C:\Users\rhember\Documents\Data\Land Cover\NTEMS\vlce2_2019.tif'
fout=meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCoverClass_NTEMS_2019.tif'
gis.ReprojectRasterAndClipToRaster(fin,fout,fref,meta['Geos']['crs'])

z=gis.OpenGeoTiff(fout)
z=gis.ClipToRaster(z,zRef)

ind=np.where(zRef['Data']==0)
z['Data'][ind]=0

gis.SaveGeoTiff(z,fout)

#%% Reproject Age from NTEMS

fin=r'C:\Users\rhember\Documents\Data\Age\NTEMS\age_ntem_c.tif'
fout=meta['Paths']['bc1ha'] + '\\Age\\Age_NTEM_2019.tif'
gis.ReprojectRasterAndClipToRaster(fin,fout,meta['Paths']['bc1ha Ref Grid'],meta['Geos']['crs'])

#%% Rasterize land use legal

# # Define paths
# meta={}
# meta['Paths']={}
# meta['Paths']['Project']=r''
# #meta['Paths']['Geospatial']=meta['Paths']['Project'] + '\\Geospatial'
# meta['Paths']['Results']=r'C:\Users\rhember\Documents\Data\ForestInventory\Results\20210930'
# meta['Paths']['VRI']=r'C:\Users\rhember\Documents\Data\ForestInventory\VRI\20210930'
# meta['Paths']['Disturbances']=r'C:\Users\rhember\Documents\Data\ForestInventory\Disturbances\20210930'
# meta['Paths']['LandUse']=r'C:\Users\rhember\Documents\Data\ForestInventory\LandUse\20210930'
# meta['Paths']['Model Code']=r'C:\Users\rhember\Documents\Code_Python\fcgadgets\cbrunner'
# meta['Paths']['Taz Datasets']=r'C:\Users\rhember\Documents\Data\Taz Datasets'

# meta=invu.Load_LUTs(meta)

# # Input path to RESULTS database (downloaded from BC data catalogue)
# pthin=r'C:\Users\rhember\Documents\Data\ForestInventory\LandUse\20210930\LandUse.gdb'
# pthout=r'C:\Users\rhember\Documents\Data\BC1ha\LandUseLandCover'

# fiona.listlayers(pthin)
# lnam='RMP_PLAN_LEGAL_POLY_SVW'

# list(lut['LU L']['LEGAL_FEAT_OBJECTIVE'].keys())

# # Open the shapefile
# df=gpd.read_file(pthin,layer=lnam)

# # Remove features with no geometry
# df=df[df.geometry!=None]

# fn='LEGAL_FEAT_OBJECTIVE'

# for val in lut['LU L']['LEGAL_FEAT_OBJECTIVE'].keys():
#     print(val)
#     #val='Conservation Lands'

#     df0=df[df[fn]==val].copy()
#     df0['dummy']=np.ones(len(df0))
#     shapes=((geom,value) for geom, value in zip(df0['geometry'],df0['dummy']))

#     z0=np.zeros(zRef['Data'].shape,dtype=float)
#     burned=features.rasterize(shapes=shapes,fill=0,out=z0,transform=zRef['Transform'])

#     val1=val.replace(":","_")
#     val1=val1.replace("/","_")
#     val1=val1.replace("<","Less Than")

#     zOut=zRef.copy()
#     zOut['Data']=np.zeros(zRef['Data'].shape,dtype='int16')
#     zOut['Data'][burned>0]=1
#     gis.SaveGeoTiff(zOut,pthout + '\\' + fn + '_' + val1 + '.tif')

# #plt.matshow(zOut['Data'][::50,::50])

# #fout=r'Z:\!Workgrp\Forest Carbon\Data\BC1ha\Disturbances' + '\\' + lnam + '_' + str(tv[iT]) + '.tif'
# #gis.SaveGeoTiff(zOut,fout)

# # Test
# #plt.matshow(z['Data'])

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

#%% Import climate from old BC1ha project in matlab

# Annual summary
z=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\Terrain\BC_1ha_twi.tif')
z=gis.ClipRaster(z,[zS['xmin'],zS['xmax']],[zS['ymin'],zS['ymax']])
gis.SaveGeoTiff(z,meta['Paths']['bc1ha'] + '\Terrain\bc1ha_twi.tif')
del z
garc.collect()

# Annual summary
z=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Climate\BC1ha_tmin_ann_norm_1971to2000_si_hist_v1_c.tif')
z=gis.ClipRaster(z,[zS['xmin'],zS['xmax']],[zS['ymin'],zS['ymax']])
gis.SaveGeoTiff(z,r'C:\Users\rhember\Documents\Data\BC1ha\Climate\BC1ha_tmin_ann_norm_1971to2000_si_hist_v1.tif')
del z; garc.collect()


#%% Extract mean climate data by BGC zone

zBGC=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\becz.tif')
zBGC['Data']=zBGC['Data'].flatten()
lutBGC=gu.ReadExcel(meta['Paths']['bc1ha'] + '\\VRI 2023\\becz_lut.xlsx')

zMAT=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Climate\BC1ha_mat_norm_1971to2000_si_hist_v1.tif')
zMAT['Data']=zMAT['Data'].flatten().astype(float)/10

zWS=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Climate\BC1ha_ws_gs_norm_1971to2000_comp_hist_v1.tif')
zWS['Data']=zWS['Data'].flatten().astype(float)

zSI=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\si.tif')
zSI['Data']=zSI['Data'].flatten().astype(float)

zA=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\proj_age_1.tif')
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

    z2=zRef.copy()
    z2['Data']=np.zeros(zRef['Data'].shape,dtype='float')
    for mo in range(12):
        print(mo+1)
        fin=r'E:\Data\Climate\Canada\BC\Grids\BC1ha_tmean_mon_norm_1971to2000_si_hist_v1\BC1ha_tmean_mon_norm_1971to2000_si_hist_v1_' + str(mo+1) + '.mat'
        z=spio.loadmat(fin,squeeze_me=True)
        idat=np.where(np.asarray(z['z'].dtype.names)=='Data')[0][0]
        iSF=np.where(np.asarray(z['z'].dtype.names)=='ScaleFactor')[0][0]
        z0=np.flip(z['z'][()][idat],axis=0) # #.astype(float)*z['z'][()][iSF]
        z1=zRef.copy()
        z1['Data']=z0.astype('int16')
        z1=gis.ClipToRaster(z1,zRef)
        gis.SaveGeoTiff(z1,r'C:\Users\rhember\Documents\Data\BC1ha\Climate\Monthly\BC1ha_tmean_mon_norm_1971to2000_si_hist_v1\BC1ha_tmean_mon_norm_1971to2000_si_hist_v1_' + str(mo+1) + '.tif')
        z2['Data']=z2['Data']+z1['Data'].astype('float')
    z2['Data']=z2['Data']/12
    z2['Data']=z2['Data'].astype('int16')
    gis.SaveGeoTiff(z2,r'C:\Users\rhember\Documents\Data\BC1ha\Climate\Seasonal\BC1ha_tmean_ann_norm_1971to2000_si_hist_v1.tif')

    z2=zRef.copy()
    z2['Data']=np.zeros(zRef['Data'].shape,dtype='int16')
    for mo in range(12):
        print(mo+1)
        fin=r'E:\Data\Climate\Canada\BC\Grids\BC1ha_prcp_mon_norm_1971to2000_si_hist_v1\BC1ha_prcp_mon_norm_1971to2000_si_hist_v1_' + str(mo+1) + '.mat'
        z=spio.loadmat(fin,squeeze_me=True)
        idat=np.where(np.asarray(z['z'].dtype.names)=='Data')[0][0]
        iSF=np.where(np.asarray(z['z'].dtype.names)=='ScaleFactor')[0][0]
        z0=np.flip(z['z'][()][idat].astype(float)*z['z'][()][iSF],axis=0)
        z1=zRef.copy()
        z1['Data']=z0.astype('int16')
        z1=gis.ClipToRaster(z1,zRef)
        gis.SaveGeoTiff(z1,r'C:\Users\rhember\Documents\Data\BC1ha\Climate\Monthly\BC1ha_prcp_mon_norm_1971to2000_si_hist_v1\BC1ha_prcp_mon_norm_1971to2000_si_hist_v1_' + str(mo+1) + '.tif')
        z2['Data']=z2['Data']+np.maximum(0,z1['Data'])
    gis.SaveGeoTiff(z2,r'C:\Users\rhember\Documents\Data\BC1ha\Climate\Seasonal\BC1ha_prcp_ann_norm_1971to2000_si_hist_v1.tif')

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

#%% Convert old NACID .mat files to new BC1ha geotiffs with standardized extent

def ConvertClimateNormals():
    zRef=gis.OpenGeoTiff(r'E:\Data\Climate\NACID\Geotiff\NACID\grid.tif')
    for mo in range(12):
        print(mo+1)
        z=spio.loadmat(r'E:\Data\Climate\NACID\Grids\NACID_rswd_mon_norm_1971to2000_si_hist_v1\NACID_rswd_mon_norm_1971to2000_si_hist_v1_' + str(mo+1) + '.mat',squeeze_me=True)
        idat=np.where(np.asarray(z['z'].dtype.names)=='Data')[0][0]
        iSF=np.where(np.asarray(z['z'].dtype.names)=='ScaleFactor')[0][0]
        z0=np.flip(z['z'][()][idat].astype(float)*z['z'][()][iSF],axis=0)
        z1=zRef.copy()
        z1['Data']=10*z0
        z1['Data']=z1['Data'].astype('int16')
        #z1=gis.ClipToRaster(z1,zRef)
        gis.SaveGeoTiff(z1,r'E:\Data\Climate\NACID\Geotiff\NACID\NACID_rswd_mon_norm_1971to2000_si_hist_v1\NACID_rswd_mon_norm_1971to2000_si_hist_v1_' + str(mo+1) + '.tif')

    return

z=gis.OpenGeoTiff(r'E:\Data\Climate\NACID\Geotiff\NACID\NACID_rswd_mon_norm_1971to2000_si_hist_v1\NACID_rswd_mon_norm_1971to2000_si_hist_v1_1.tif')

