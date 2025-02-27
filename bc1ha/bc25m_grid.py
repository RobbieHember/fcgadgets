
#%% Import modules
import numpy as np
import matplotlib.pyplot as plt
import geopandas as gpd
import pandas as pd
from shapely.geometry import Point,Polygon
import fcgadgets.macgyver.util_general as gu
import fcgadgets.macgyver.util_gis as gis
import fcgadgets.bc1ha.bc1ha_utils as u1ha
import fcgadgets.bc1ha.bc25m_util as u25m
gp=gu.SetGraphics('Manuscript')

#%% Import project info
meta=u1ha.Init()
meta=u25m.Init_bc25m(meta)

#%% Define reference grid
u25m.DefineReferenceGrid(meta)

zRef=gis.OpenGeoTiff(meta['Paths']['bc25m Ref Grid'])
del zRef['Data']
gdf=u1ha.Import_GDBs_ProvinceWide(meta)

#%% Define tiles

CellSize=25 # metres
srs=gis.ImportSRSs()
cnt=0
x=np.linspace(zRef['xmin'],zRef['xmax'],20)
y=np.linspace(zRef['ymin'],zRef['ymax'],20)
gdf_tile=gpd.GeoDataFrame()
dLL={}
for i in range(y.size-1):
	for j in range(x.size-1):

		# Vector bounary of tile
		x_lo=x[j]
		x_hi=x[j+1]
		xv=np.array([x_lo,x_hi,x_hi,x_lo,x_lo])
		y_lo=y[i]
		y_hi=y[i+1]
		yv=np.array([y_hi,y_hi,y_lo,y_lo,y_hi])

		points=[]
		for k in range(len(yv)):
			points.append(Point(xv[k],yv[k]))
		coords=[(p.x,p.y) for p in points]
		poly=gpd.GeoSeries(Polygon(coords))

		# Geographic
		lon0,lat0=pyproj.transform(srs['Proj']['BC1ha'],srs['Proj']['Geographic'],xv,yv)

		gdf_tile0=gpd.GeoDataFrame({'ID':cnt, \
			'i':i, \
			'j':j, \
			'x':x[j], \
			'y':y[i], \
			'Cellsize':CellSize, \
			'WithOverlap':0, \
			'geometry':poly},
			crs=meta['Geos']['crs'])
		a=gpd.overlay(gdf['bc_bound']['gdf'],gdf_tile0,how='intersection')
		if len(a)>0:
			gdf_tile0['WithOverlap']=1
			dLL[cnt]=[np.min(lon0)-0.02,np.min(lat0)-0.02,np.max(lon0)+0.02,np.max(lat0)+0.02]
		cnt=cnt+1
		#gdf_tile=gdf_tile.append(gdf_tile0,ignore_index=True)
		gdf_tile=pd.concat([gdf_tile,gdf_tile0],ignore_index=True)

# Save tiles
#gdf_tile.to_file(meta['Paths']['Data'] + '\\TilesBC.geojson',driver='GeoJSON')

plt.close('all'); fig,ax=plt.subplots(figsize=gu.cm2inch(14,14))
#roi['gdf']['bound'].plot(ax=ax,facecolor='None',edgecolor=[0,0,0],label='Area',linewidth=1,alpha=1)
gdf['bc_bound']['gdf'].plot(ax=ax,facecolor=[0.8,0.8,0.8],edgecolor=[0,0,0],label='ROI',linewidth=0.25,alpha=1)
gdf_tile[(gdf_tile['WithOverlap']==1)].plot(ax=ax,facecolor='None',edgecolor=[1,0,0],linewidth=1,alpha=0.4)
#gdf_tile[(gdf_tile['i']==9) & (gdf_tile['j']==9)].plot(ax=ax,facecolor='r',edgecolor=[0.25,0,0],linewidth=1.5,alpha=0.5)
ax.set(position=[0.01,0.01,0.98,0.98],xticks=[],yticks=[])

for i in range(len(gdf_tile)):
	if gdf_tile.loc[i]['WithOverlap']==0:
		continue
	x=gdf_tile.loc[i].geometry.centroid.x
	y=gdf_tile.loc[i].geometry.centroid.y
	ax.text(x,y,str(gdf_tile.loc[i]['ID']),ha='center',va='center')

# Index to cells with land
indWO=[]
for i in range(len(gdf_tile)):
	if gdf_tile.loc[i]['WithOverlap']==0:
		continue
	indWO.append(i)
len(indWO)

#%% Import canopy height
import ee
ee.Authenticate(force=True)
ee.Initialize(project='cdf-zone1')

#%% Create tasks
scale=25
crs='EPSG:3005'
tva=np.arange(2002,2025,1)
iTile=70
for yr in tva:
	mo=7
	if mo<10:
		smo='0' + str(mo)
	else:
		smo=str(mo)
	t0=str(yr) + '-' + smo + '-01'
	t1=str(yr) + '-' + smo + '-31'

	im0=ee.ImageCollection("LANDSAT/COMPOSITES/C02/T1_L2_8DAY_NDWI").filterDate(t0,t1).select('NDWI').mosaic()
	#for j in range(5999,8000):
		#i=indWO[j]

	# Clip
	geomPoly=ee.Geometry.BBox(dLL[iTile][0],dLL[iTile][1],dLL[iTile][2],dLL[iTile][3]) # East,south,west,north
	im1=im0.clip(geomPoly)
	# Create forest land mask based on height
	#flm0=ch1.gte(h_TH)
	# Reproject
	im2=im1.reproject(crs=crs,scale=scale)
	# Reducer
	#flm0=flm0.reduceResolution(reducer=ee.Reducer.mean(),maxPixels=65536)#.reproject(crs=crs2,scale=scale3))
	#flm0=flm0.toDouble()
	# Start task
	name='ndwi_id' + str(iTile) + '_' + str(yr) + '_' + str(mo)
	task=ee.batch.Export.image.toDrive(image=im2,description=name,scale=scale,maxPixels=1e13,region=geomPoly,fileFormat="GeoTIFF") #
	task.start()

#%% Status
task.status()

#%%
z=gis.OpenGeoTiff(r'C:\Data\BC25m\NDWI\ndwi_id69.tif')
plt.close('all'); plt.matshow(z['Data'],clim=[-1,1]); plt.colorbar()





#%% Reproject CEC Land Use Map 2010
fin=r'C:\Users\rhember\Documents\Data\Land Cover\NALCMS\NALCMS_lc_2010.tif'
fout=r'C:\Users\rhember\Documents\Data\Land Cover\NALCMS\NALCMS_lc_2010c.tif'
gis.ClipToRaster_ByFile(fin,fout,meta['Paths']['bc1ha Ref Grid'])

#z=gis.OpenGeoTiff(fin)

fin=r'C:\Users\rhember\Documents\Data\Land Cover\NALCMS\NALCMS_lc_2010.tif'
fout=meta['Paths']['bc25m'] + '\\LandCoverUse\\LandCover_CEC_2010.tif'
gis.ReprojectRasterAndClipToRaster(fin,fout,meta['Paths']['bc25m Ref Grid'],meta['Geos']['crs'])

z=gis.OpenGeoTiff(fout)
z=gis.ClipToRaster(z,zRef)

ind=np.where(zRef['Data']==0)
z['Data'][ind]=0

gis.SaveGeoTiff(z,)

#z=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\CEC_LandCover2010.tif')
#plt.matshow(z['Data'][0::5,0::5])

#z=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\Land Cover\NALCMS\NALCMS_lc_2010.tif')

