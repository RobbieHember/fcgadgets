#%% Import modules
import os
from osgeo import gdal
import geopandas as gpd
import numpy as np
from osgeo import osr
from matplotlib import path
import copy
from shapely.geometry import Polygon,Point,LineString
import pyproj
from scipy import ndimage
import rasterio
from rasterio.features import shapes
from rasterio.enums import Resampling
import cv2
from shapely.geometry import Point, Polygon
from shapely import geometry
from rasterio.warp import calculate_default_transform, reproject, Resampling
from rasterio.transform import from_origin
from fcgadgets.macgyver import util_general as gu

#%% Import common spatial reference systems
# srs=gis.ImportSRSs()
def ImportSRSs():
	# Example:
	#	srs=gis.ImportSRSs()
	#	x,y=srs['Proj']['BC1ha'](lon[:,1],lat[:,0])
	#   lon,lat=pyproj.transform(srs['Proj']['BC1ha'],srs['Proj']['Geographic'],pl['X BC'],pl['Y BC'])
	#   x_NA,y_NA=pyproj.transform(srs['Proj']['BC1ha'],srs['Proj']['NACID'],pl['X BC'],pl['Y BC'])
	srs={}
	srs['String']={}
	srs['Proj']={}

	# Geographic
	srs['String']['Geographic']='+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
	srs['Proj']['Geographic']=pyproj.Proj(srs['String']['Geographic'])
	#srs['Proj']['Geographic']=pyproj.Proj(init='epsg:4326')

	# NACID
	srs['String']['NACID']='+proj=lcc +lat_1=49 +lat_2=77 +lat_0=40 +lon_0=-100 +x_0=0 +y_0=0 +no_defs +a=6378137 +rf=298.257222101 +to_meter=1'
	srs['Proj']['NACID']=pyproj.Proj(srs['String']['NACID'])

	# BC1ha
	srs['String']['BC1ha']='+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +no_defs +a=6378137 +rf=298.257222101 +to_meter=1'
	srs['Proj']['BC1ha']=pyproj.Proj(srs['String']['BC1ha'])

	# Robinson
	srs['String']['Robinson']='+proj=merc +a=6378137 +b=6378137 +lat_ts=0 +lon_0=0 +x_0=0 +y_0=0 +k=1 +units=m +nadgrids=@null +wktext +no_defs'
	srs['Proj']['Robinson']=pyproj.Proj(srs['String']['Robinson'])
	return srs

#%% Resample raster
def ResampleRaster(fin,sf):
	with rasterio.open(fin) as dataset:
		data=dataset.read(out_shape=(dataset.count,int(dataset.height*sf),int(dataset.width*sf)),resampling=Resampling.nearest)
		#transform=dataset.transform*dataset.transform.scale((dataset.width/data.shape[-1]),(dataset.height/data.shape[-2]))
	z=OpenGeoTiff(fin)
	gt=list(z['gt'])
	gt[1]=gt[1]/sf
	gt[5]=gt[5]/sf
	gt=tuple(gt)
	if data.shape[0]==1:
		data=np.squeeze(data)
		m,n=data.shape
	else:
		n_bands,m,n=data.shape
	Projection=z['Projection']
	proj4_str=z['Proj4_String']
	xmin=gt[0]
	ymax=gt[3]
	xmax=xmin+gt[1]*n
	ymin=ymax+gt[5]*m
	extent=(xmin,xmax,ymin,ymax)
	try:
		x=np.arange(gt[0],gt[0]+gt[1]*n,gt[1])
	except:
		x=np.arange(gt[0],gt[0]+gt[1]*n-gt[1],gt[1])
	try:
		y=np.flip(np.arange(gt[3]-m*gt[1],gt[3],gt[1]))
	except:
		y=np.flip(np.arange(gt[3]-m*gt[1]+gt[1],gt[3],gt[1]))
	x=np.tile(x,(int(m),1))
	y=np.tile(np.reshape(y,(m,1)),(1,int(n)))
	Cellsize=gt[1]
	Transform=from_origin(xmin,ymax,Cellsize,Cellsize)# Transform
	yxrat=data.shape[0]/data.shape[1] # y-x ratio
	z={'gt':gt,
	   'Data':data,
	   'X':x,
	   'Y':y,
	   'm':data.shape[0],
	   'n':data.shape[1],
	   'yxrat':yxrat,
	   'xmin':xmin,
	   'xmax':xmax,
	   'ymin':ymin,
	   'ymax':ymax,
	   'xlim':[xmin,xmax],
	   'ylim':[ymin,ymax],
	   'Extent':extent,
	   'Cellsize':Cellsize,
	   'Transform':Transform,
	   'Projection':Projection,
	   'Proj4_String':proj4_str}
	SaveGeoTiff(z,fin)
	return

#%%
# srs=gis.ImportSRSs()
# gis.ReprojectCoordinates(srs['Proj']['Geographic'],srs['Proj']['NACID'],x_in,y_in)
def ReprojectCoordinates(proj_in,proj_out,x_in,y_in):
	x_out,y_out=pyproj.transform(proj_in,proj_out,x_in,y_in)
	return x_out,y_out

#%% Open Geotiff
def OpenGeoTiff(pthin):

	ds=gdal.Open(pthin)
	gt=ds.GetGeoTransform()

	raster=rasterio.open(pthin)
	data=raster.read()

	if data.shape[0]==1:
		data=np.squeeze(data)
		m,n=data.shape
	else:
		n_bands,m,n=data.shape

	Projection=ds.GetProjection()
	prj=osr.SpatialReference(wkt=ds.GetProjection())
	proj4_str=prj.GetAttrValue('AUTHORITY',1)

	xmin=gt[0]
	ymax=gt[3]
	xmax=xmin+gt[1]*n
	ymin=ymax+gt[5]*m
	extent=(xmin,xmax,ymin,ymax)

	try:
		x=np.arange(gt[0],gt[0]+gt[1]*n,gt[1])
	except:
		x=np.arange(gt[0],gt[0]+gt[1]*n-gt[1],gt[1])

	try:
		y=np.flip(np.arange(gt[3]-m*gt[1],gt[3],gt[1]))
	except:
		y=np.flip(np.arange(gt[3]-m*gt[1]+gt[1],gt[3],gt[1]))

	x=np.tile(x,(int(m),1))

	try:
		y=np.tile(np.reshape(y,(m,1)),(1,int(n)))
	except:
		#print(data.shape)
		#print(y.shape)
		#print(x.shape)
		try:
			y=np.flip(np.arange(gt[3]-m*gt[1]+gt[1],gt[3],gt[1]))
			#print(y.shape)
			y=np.tile(np.reshape(y,(m,1)),(1,int(n)))
		except:
			y=np.flip(np.arange(gt[3]-m*gt[1],gt[3],gt[1]))
			y=y[0:-1]#.flatten()
			#print(y.shape)
			y=np.tile(np.reshape(y,(m,1)),(1,int(n)))

	Cellsize=gt[1]

	# Transform
	Transform=from_origin(xmin,ymax,Cellsize,Cellsize)

	# y-x ratio
	yxrat=data.shape[0]/data.shape[1]

	z={'gt':gt,
	   'Data':data,
	   'X':x,
	   'Y':y,
	   'm':data.shape[0],
	   'n':data.shape[1],
	   'yxrat':yxrat,
	   'xmin':xmin,
	   'xmax':xmax,
	   'ymin':ymin,
	   'ymax':ymax,
	   'xlim':[xmin,xmax],
	   'ylim':[ymin,ymax],
	   'Extent':extent,
	   'Cellsize':Cellsize,
	   'Transform':Transform,
	   'Projection':Projection,
	   'Proj4_String':proj4_str,
	   'crs':raster.crs}

	#z=Bunch(z)

	return z

#%% Save geotiff
def SaveGeoTiff(z,fout):
	driver=gdal.GetDriverByName("GTiff")
	# Data type
	if (z['Data'].dtype=='int8') | (z['Data'].dtype=='uint8'):
		dtype=gdal.GDT_Byte
	elif z['Data'].dtype=='int16':
		dtype=gdal.GDT_Int16
	elif z['Data'].dtype=='uint16':
		dtype=gdal.GDT_Int16
	elif z['Data'].dtype=='int32':
		dtype=gdal.GDT_Int32
	elif z['Data'].dtype=='float32':
		dtype=gdal.GDT_Float32
	else:
		print('Did not work, input dtype is:' + z['Data'].dtype)

	N_band=1
	ds=driver.Create(fout,z['n'],z['m'],N_band,dtype,[ 'COMPRESS=LZW' ])
	ds.SetProjection(z['Projection'])
	ds.SetGeoTransform(z['gt'])
	ds.GetRasterBand(1).WriteArray(z['Data'])
	return

#%% CLIP RASTER
#The raster input structure must be from the OpenGdal function from this module.
def ClipRasterByXYLimits(z_in,xlim,ylim):
	z=z_in.copy()
	ix=np.where((z['X'][0,:]>=xlim[0]) & (z['X'][0,:]<=xlim[1]))[0]
	iy=np.where((z['Y'][:,0]>=ylim[0]) & (z['Y'][:,0]<=ylim[1]))[0]
	ind=np.ix_(iy,ix)
	z['Data']=z['Data'][ind]
	z['X']=z['X'][ind]
	z['Y']=z['Y'][ind]
	z['m']=iy.size
	z['n']=ix.size
	z['xmin']=np.min(z['X'])
	z['xmax']=np.max(z['X'])
	z['ymin']=np.min(z['Y'])
	z['ymax']=np.max(z['Y'])
	z['gt']=(z['xmin'],z['gt'][1],z['gt'][2],z['ymax'],z['gt'][4],z['gt'][5])
	z['Extent']=(z['xmin'],z['xmax'],z['ymin'],z['ymax'])
	z['xlim']=[np.min(z['X']),np.max(z['X'])]
	z['ylim']=[np.min(z['Y']),np.max(z['Y'])]
	z['Transform']=from_origin(z['xmin'],z['ymax'],z['Cellsize'],z['Cellsize'])
	return z

#%% Clip raster based on extent of other raster
def ClipToRaster(z_in0,z_ref0):

	z_in=z_in0.copy()
	z_ref=z_ref0.copy()

	ix_in=np.where((z_in['X'][0,:]>=z_ref['xlim'][0]) & (z_in['X'][0,:]<=z_ref['xlim'][1]))[0]
	iy_in=np.where((z_in['Y'][:,0]>=z_ref['ylim'][0]) & (z_in['Y'][:,0]<=z_ref['ylim'][1]))[0]
	xmin=z_in['X'][0,ix_in[0]]
	xmax=z_in['X'][0,ix_in[-1]]
	ymin=z_in['Y'][iy_in[-1],0]
	ymax=z_in['Y'][iy_in[0],0]

	ix_ref=np.where((z_ref['X'][0,:]>=xmin) & (z_ref['X'][0,:]<=xmax))[0]
	iy_ref=np.where((z_ref['Y'][:,0]>=ymin) & (z_ref['Y'][:,0]<=ymax))[0]

	# Fix if there is a mismatch
	dx=ix_ref.size-ix_in.size
	dy=iy_ref.size-iy_in.size
	if (dx==1) & (dy==1):
		ix_ref=np.where((z_ref['X'][0,:]>=xmin) & (z_ref['X'][0,:]<xmax))[0]
		iy_ref=np.where((z_ref['Y'][:,0]>=ymin) & (z_ref['Y'][:,0]<ymax))[0]
	elif (dx==1) & (dy==0):
		ix_ref=np.where((z_ref['X'][0,:]>=xmin) & (z_ref['X'][0,:]<xmax))[0]
	elif (dx==0) & (dy==1):
		iy_ref=np.where((z_ref['Y'][:,0]>=ymin) & (z_ref['Y'][:,0]<ymax))[0]
	elif (dx==-1) & (dy==-1):
		ix_in=np.where((z_in['X'][0,:]>=xmin) & (z_in['X'][0,:]<xmax))[0]
		iy_in=np.where((z_in['Y'][:,0]>=ymin) & (z_in['Y'][:,0]<ymax))[0]
	elif (dx==0) & (dy==-1):
		iy_in=np.where((z_in['Y'][:,0]>=ymin) & (z_in['Y'][:,0]<ymax))[0]
	else:
		pass

	ind_in=np.ix_(iy_in,ix_in)
	ind_ref=np.ix_(iy_ref,ix_ref)

	z={}
	z['Data']=0*z_ref['Data']
	z['Data']=z['Data'].astype(z_in['Data'].dtype)
	z['Data'][ind_ref]=z_in['Data'][ind_in]

	for k in z_ref.keys():
		if (k=='Data'):
			continue
		z[k]=z_ref[k]

	return z

#%% Clip raster based on extent of other raster (from files)
def ClipToRaster_ByFile(fin,fout,fref):

	z_in=OpenGeoTiff(fin)
	z_ref=OpenGeoTiff(fref)
	z_out=ClipToRaster(z_in,z_ref)

	# Force output to be the same datatype as input
	z_out['Data']=z_out['Data'].astype(z_in['Data'].dtype)

	SaveGeoTiff(z_out,fout)

	return

#%% Reproject and clip raster
def ReprojectRasterAndClipToRaster(fin,fout,fref,crs):
	# Only outputs integers

	# Grid to reproject
	ds=gdal.Open(fin,0)

	# Open reference grid
	z_ref=OpenGeoTiff(fref)

	# Adopt the data type of the input
	#z_ref['Data']=z_ref['Data'].astype(z_in0['Data'].dtype)
	z_ref['Data']=z_ref['Data'].astype('float32')

	#fout_tmp=fout + 'a.tif'

	# Reproject
	ds_p=gdal.Warp(fout,fin,
				   xRes=z_ref['Cellsize'],
				   yRes=z_ref['Cellsize'],
				   dstSRS=crs,
				   creationOptions=["COMPRESS=LZW"] )

	# Open reprojected raster
	z_out=OpenGeoTiff(fout)

	ds=None
	ds_p=None
	#os.remove(fout_tmp)

	# Clip to extent of reference grid
	z=ClipToRaster(z_out,z_ref)

	# Save clipped and reprojected grid
	if np.max(z['Data'])>32767:
		z['Data']=z['Data'].astype('int32')
	else:
		z['Data']=z['Data'].astype('int16')
	SaveGeoTiff(z,fout)

	return z

#%% Reproject geotiff
def ReprojectGeoTiff(pthin,pthout,crs_dst):
	cmp='lzw'
	with rasterio.open(pthin) as src:
		transform,width,height=calculate_default_transform(src.crs,crs_dst,src.width,src.height,*src.bounds)
		kwargs=src.meta.copy()
		kwargs.update({'crs':crs_dst,'transform':transform,'width':width,
					   'height':height,'dtype':rasterio.int16,
					   'compress':cmp})
		with rasterio.open(pthout, 'w',**kwargs) as dst:
			for i in range(1,src.count+1):
				reproject(
					source=rasterio.band(src,i),
					destination=rasterio.band(dst,i),
					src_transform=src.transform,
					src_crs=src.crs,
					dst_transform=transform,
					dst_crs=crs_dst,
					resampling=Resampling.nearest)
	return

#%% Polygon area
def PolyArea(x,y):
	return 0.5*np.abs(np.dot(x,np.roll(y,1))-np.dot(y,np.roll(x,1)))

#%% Inside polygon
def InPolygon(xg,yg,xv,yv):

	l=[]
	for i in range(xv.size):
		l.append([xv[i],yv[i]])
	poly=Polygon(l)

	InPol=np.zeros(xg.shape)

	if xg.ndim<=1:
		# Points in one dimension
		for i in range(xg.size):
			pnt=Point(xg[i],yg[i])
			if poly.contains(pnt)==True:
				InPol[i]=1
	else:
		# Two-dimensional grid
		for i in range(xg.shape[0]):
			for j in range(xg.shape[1]):
				pnt=Point(xg[i,j],yg[i,j])
				if poly.contains(pnt)==True:
					InPol[i,j]=1

	return InPol

#%% Get extent from gdal
def GetExtentFromGDAL(gt,cols,rows):
	''' Return list of corner coordinates from a geotransform

		@type gt:   C{tuple/list}
		@param gt: geotransform
		@type cols:   C{int}
		@param cols: number of columns in the dataset
		@type rows:   C{int}
		@param rows: number of rows in the dataset
		@rtype:	C{[float,...,float]}
		@return:   coordinates of each corner
	'''
	ext=[]
	xarr=[0,cols]
	yarr=[0,rows]

	for px in xarr:
		for py in yarr:
			x=gt[0]+(px*gt[1])+(py*gt[2])
			y=gt[3]+(px*gt[4])+(py*gt[5])
			ext.append([x,y])
		yarr.reverse()
	return ext

#%% COMPRESS CATEGORIES IN RASTER DATASET
# z1,lab1,cl1=gis.CompressCats(z0,id0,lab0,cl0)
def CompressCats(z0,id0,lab0,cl0):
	uc=np.unique(z0)
	z1=uc.size*np.ones(z0.shape)
	cl1=0*np.ones((uc.size,3))
	lab1=[None]*uc.size
	for i in range(uc.size):
		ind=np.where(z0==uc[i])
		z1[ind]=i
		ind=np.where(id0==uc[i])[0]
		if ind.size==0:
			lab1[i]='Unknown'
			cl1[i,:]=[0,0,0]
		else:
			lab1[i]=lab0[ind[0]]
			cl1[i,:]=cl0[ind[0],:]
	lab1=np.array(lab1)
	return z1,lab1,cl1

#%% Digitize raster binary mask
# gdf=gis.DigitizeBinaryMask(zIn)
def DigitizeBinaryMask(zIn):

	# Create binary image
	z=np.zeros((zIn['m'],zIn['n'],3),dtype='uint8')
	z[zIn['Data']==1,:]=255
	z=cv2.cvtColor(z,cv2.COLOR_BGR2GRAY) # Convert to grey scale

	xg=zIn['X'][0,:]
	yg=zIn['Y'][:,0]

	# Calculate contour of object
	cont=cv2.findContours(image=z,mode=cv2.RETR_LIST,method=cv2.CHAIN_APPROX_SIMPLE)

	# Unpack silly tuple
	gdf=gpd.GeoDataFrame(data=[],columns=['Value','geometry'])
	cnt=0
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
			x=xg[c]
			y=yg[r]
			pointList.append(geometry.Point(x,y))
		gdf.loc[cnt,'Value']=1
		#gdf.loc[cnt,'Name']=df_tsa.loc[i,'Name']
		gdf.loc[cnt,'geometry']=geometry.Polygon([[p.x,p.y] for p in pointList])
		cnt=cnt+1

	return gdf

# #%% OLD digitize:
# def Digitize(BinaryMask,xv,yv):

#	 #xv=xv.flatten()
#	 #yv=yv.flatten()

#	 s=shapes(BinaryMask.astype('int16'),mask=None,connectivity=4)

#	 xy=[]
#	 for i in range(10000):
#		 try:
#			 a=next(s)
#			 b=a[0]['coordinates'][0]
#			 d={}
#			 d['x']=np.array([])
#			 d['y']=np.array([])
#			 for j in range(len(b)):
#				 col=int(b[j][0])
#				 row=int(b[j][1])
#				 d['x']=np.append(d['x'],xv[col-1])
#				 d['y']=np.append(d['y'],yv[row-1])
#			 xy.append(d)
#		 except:
#			 break
#	 return xy

#%% Clip geodataframe to user-specified x and y limits
def ClipGDF(gdf_in,xlim,ylim):

	gdf=gdf_in.cx[xlim[0]:xlim[1],ylim[0]:ylim[1]]
	#gdf=gdf.reset_index(drop=True)
	#gdf=gpd.sjoin(gdf,roi['gdb']['bound'],how='left')

	# This grouped by may not be necessary - it prevents the file from working in overlays
	#gdf=gdf.groupby('index_right')

	return gdf

#%% Adjust grid cellsize based on regular subsampling
# z=gis.UpdateGridCellsize(z_in,scale_factor)
def UpdateGridCellsize(z_in,scale_factor):
	z=z_in.copy()
	z['Data']=z['Data'][0::scale_factor,0::scale_factor]
	z['X']=z['X'][0::scale_factor,0::scale_factor]
	z['Y']=z['Y'][0::scale_factor,0::scale_factor]
	z['m'],z['n']=z['Data'].shape
	z['Cellsize']=z['Cellsize']*scale_factor
	gt=list(z['gt'])
	gt[1]=gt[1]*scale_factor
	gt[5]=gt[5]*scale_factor
	z['gt']=tuple(gt)
	t=list(z['Transform'])
	t[1]=t[0]*scale_factor
	t[5]=t[4]*scale_factor
	z['Transform']=tuple(t)
	return z

#%% Get grid index of points from x and y
# idx=gis.GetGridIndexToPoints(z,x,y)
def GetGridIndexToPoints(z,x,y):
	Xg=z['X'][0,:]
	Yg=z['Y'][:,0]
	ix=np.zeros(x.size,dtype='int16')
	iy=np.zeros(x.size,dtype='int16')
	for iPoint in range(x.size):
		if np.isnan(x[iPoint])==True:
			continue
		adx=np.abs(x[iPoint]-Xg)
		ady=np.abs(y[iPoint]-Yg)
		ind=np.where(adx==np.min(adx))[0]
		if ind.size>1:
			ind=ind[0]
		ix[iPoint]=ind
		ind=np.where(ady==np.min(ady))[0]
		if ind.size>1:
			ind=ind[0]
		iy[iPoint]=ind
	ind=tuple([iy,ix])
	return ind

#%% Import geographic coordinates
def ConvertGeographicToGeojson(pthin,pthout,type):

	ll=gu.ReadExcel(pthin,sheet_name='Sheet1',skiprows=0)

	# Import spatial reference systems
	srs=ImportSRSs()
	ll['X']=np.zeros(ll['Lon'].size)
	ll['Y']=np.zeros(ll['Lat'].size)
	for i in range(ll['Lat'].size):
		ll['X'][i],ll['Y'][i]=srs['Proj']['BC1ha'](ll['Lon'][i],ll['Lat'][i])

	points=[]
	for i in range(ll['X'].size):
		points.append(Point(ll['X'][i],ll['Y'][i]))

	if type=='Polygon':
		poly=Polygon(points)
		out=gpd.GeoDataFrame({'geometry':[poly]})
	elif type=='Line':
		poly=LineString(points)
		out=gpd.GeoDataFrame({'geometry':[poly]})
	elif type=='Point':
		out=gpd.GeoDataFrame({'geometry':points})
	out.to_file(pthout + '.geojson',driver='GeoJSON')

	return out

#%% Import geographic coordinates
def PolygonRotate(pthin,ll,dx,dy,r):
	xy0=gu.ReadExcel(pthin,sheet_name='Sheet1',skiprows=0)

	srs=ImportSRSs()
	xy={}
	xy['X']=np.zeros(xy0['Lon'].size)
	xy['Y']=np.zeros(xy0['Lon'].size)
	for i in range(xy0['Lon'].size):
		xy['X'][i],xy['Y'][i]=srs['Proj']['BC1ha'](ll[1],ll[0])
		xy['X'][i]=xy['X'][i]+xy0['Lon'][i]+dx
		xy['Y'][i]=xy['Y'][i]+xy0['Lat'][i]+dy

	points=[]
	for i in range(xy['X'].size):
		points.append(Point(xy['X'][i],xy['Y'][i]))
	poly=Polygon(points)
	out=gpd.GeoDataFrame({'geometry':[poly]})
	out=out.rotate(r)
	#out.to_file(pthout + '.geojson',driver='GeoJSON')
	return out,xy

#%% Import Cities
def ImportCities(pthin,output_type):

	Cities=gu.ReadExcel(pthin,sheet_name='Sheet1',skiprows=0)

	# Import spatial reference systems
	srs=ImportSRSs()
	Cities['X']=np.zeros(Cities['Lat'].size)
	Cities['Y']=np.zeros(Cities['Lat'].size)
	for i in range(Cities['Lat'].size):
		Cities['X'][i],Cities['Y'][i]=srs['Proj']['BC1ha'](Cities['Lon'][i],Cities['Lat'][i])

	if output_type=='Dict':
		out=Cities
	elif output_type=='GDF':
		points=[]
		for i in range(Cities['X'].size):
			points.append(Point(Cities['X'][i],Cities['Y'][i]))
		out=gpd.GeoDataFrame({'geometry':points,'City Name':Cities['City Name'],'Territory':Cities['Territory'],'Lat':Cities['Lat'],'Lon':Cities['Lon'],'Level':Cities['Level']})

	return out

#%% Import Bioenergy
def ImportBioenergy(pthin,output_type):

	d=gu.ReadExcel(pthin,sheet_name='Sheet1',skiprows=0)

	# Import spatial reference systems
	srs=ImportSRSs()
	d['X']=np.zeros(d['Lat'].size)
	d['Y']=np.zeros(d['Lat'].size)
	for i in range(d['Lat'].size):
		d['X'][i],d['Y'][i]=srs['Proj']['BC1ha'](d['Lon'][i],d['Lat'][i])

	if output_type=='Dict':
		out=d
	elif output_type=='GDF':
		points=[]
		for i in range(d['X'].size):
			points.append(Point(d['X'][i],d['Y'][i]))
		d['geometry']=points
		out=gpd.GeoDataFrame(d)

	return out

#%% Shift input matrix by dx,dy
def imshift(In,dx,dy):
	m,n=In.shape
	Out=-999*np.ones(In.shape)
	if (dx<0) & (dy<0):
		Out[0:m-np.abs(dy),0:n-np.abs(dx)]=In[np.abs(dy):m,np.abs(dx):n]
	elif (dx<0) & (dy>0):
	  Out[dy:m,0:n-abs(dx)]=In[0:m-dy,np.abs(dx):n]
	elif (dx>0) & (dy>0):
	  Out[dy:m,dx:n]=In[0:m-dy,0:n-dx]
	elif (dx>0) & (dy<0):
	  Out[0:m-np.abs(dy),dx:n]=In[np.abs(dy):m,0:n-dx]
	elif (dx==0) & (dy<0):
	  Out[0:m-np.abs(dy),0:n]=In[np.abs(dy):m,0:n]
	elif (dx==0) & (dy>0):
	  Out[dy:m,0:n]=In[0:m-dy,0:n]
	elif (dx==0) & (dy==0):
	  Out=In
	elif dx<0 & dy==0:
	  Out[0:m,0:n-np.abs(dx)]=In[0:m,abs(dx):n]
	elif (dx>0) & (dy==0):
	  Out[0:m,dx:n]=In[0:m,0:n-dx]
	return Out

#%% Buffer raster mask (buffer radius around binary mask)
# *** Don't use this - use cv2.dilate - super easy and fast or gdf.geometry.buffer(-20) ***
#import cv2
#kernel=np.ones((2,2),np.uint8)
#a=cv2.dilate(pbg,kernel,iterations=1)
def BufferRasterMask(MaskIn,r):
	# Create circular kernel
	def CreateKernel(radius):
		kernel=np.zeros((2*radius+1, 2*radius+1))
		y,x=np.ogrid[-radius:radius+1, -radius:radius+1]
		mask=x**2 + y**2 <= radius**2
		kernel[mask]=1
		return kernel
	MaskOut=ndimage.morphology.binary_dilation(MaskIn==1,structure=CreateKernel(r))
	
	# Test
	flg=0
	if flg==1:
		A=np.zeros((120,320))
		A[60,40]=1
		A[2,140]=1
		r=10
		B=ndimage.morphology.binary_dilation(A==1,structure=CreateKernel(r))
		plt.close('all'); plt.matshow(B)
	
	return MaskOut

# def BufferRasterMask_OLD(MaskIn,wShuf):
#	 id0=1
#	 id1=2
#	 MaskOut=MaskIn.copy()
#	 bin_x=np.arange(-wShuf,wShuf+1,1)
#	 bin_y=np.arange(-wShuf,wShuf+1,1)
#	 bin_x=bin_x[bin_x!=0]
#	 bin_y=bin_y[bin_y!=0]
#	 bin_x=bin_x[np.flip(np.argsort(np.abs(bin_x)))]
#	 bin_y=bin_y[np.flip(np.argsort(np.abs(bin_y)))]
#	 for i in range(bin_y.size):
#		 for j in range(bin_x.size):
#			 MaskShuf=imshift(MaskIn,bin_x[j],bin_y[i])
#			 ind=np.where( (MaskIn!=id0) & (MaskShuf==id0) )
#			 MaskOut[ind]=id1
#	 return MaskOut

#%% Fill raster missing values with shuffle
def ShuffleFill(zIn,wShuf,idMissing):
	zOut=idMissing*np.ones(zIn.shape,dtype=zIn.dtype)
	bin_x=np.arange(-wShuf,wShuf+1,1)
	bin_y=np.arange(-wShuf,wShuf+1,1)
	bin_x=bin_x[bin_x!=0]
	bin_y=bin_y[bin_y!=0]
	bin_x=bin_x[np.flip(np.argsort(np.abs(bin_x)))]
	bin_y=bin_y[np.flip(np.argsort(np.abs(bin_y)))]
	for i in range(bin_y.size):
		for j in range(bin_x.size):
			zShuf=imshift(zIn,bin_x[j],bin_y[i])
			ind=np.where( (zIn==idMissing) & (zShuf!=idMissing) )
			zOut[ind]=zShuf[ind]
	ind=np.where(zIn!=idMissing)
	zOut[ind]=zIn[ind]
	return zOut

#%%
def CreateEmptyRaster(path,xv,yv,crs):
	z=np.random.randint(5,size=(yv.size,xv.size),dtype='int8')
	with rasterio.open(
		path,
		mode="w",
		driver="GTiff",
		height=z.shape[0],
		width=z.shape[1],
		count=1,
		dtype=z.dtype,
		crs=crs,
		transform=from_origin(xv[0],yv[-1],10,10), # Inputs (west,north,xsize,ysize)
		) as new_dataset:
		new_dataset.write(z,1)
	return