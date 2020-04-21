
import gdal
import numpy as np
import osr
from matplotlib import path
from shapely.geometry import Polygon,Point
import pyproj
import rasterio
from rasterio.warp import calculate_default_transform, reproject, Resampling

'''============================================================================
IMPORT SPATIAL REFERENCE SYSTEMS
============================================================================'''

def ImportSRSs():
    
    srs={}
    
    # NACID
    str='+proj=lcc +lat_1=49 +lat_2=77 +lat_0=40 +lon_0=-100 +x_0=0 +y_0=0 +no_defs +a=6378137 +rf=298.257222101 +to_meter=1'
    srs['NACID']=pyproj.Proj(str)
    
    # BC1ha
    str='+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +no_defs +a=6378137 +rf=298.257222101 +to_meter=1'
    srs['BC1ha']=pyproj.Proj(str)
    return srs

'''============================================================================
BUNCH VARIABLES FROM DICTIONARY
============================================================================'''

class Bunch(dict):
    def __init__(self, *args, **kwds):
        super(Bunch, self).__init__(*args, **kwds)
        self.__dict__ = self

'''============================================================================
OPEN RASTER USING GDAL
Forerly "OpenGdal"
============================================================================'''

def OpenGeoTiff(pthin):
        
    ds=gdal.Open(pthin)
    
    gt=ds.GetGeoTransform(); 
    
    data=ds.ReadAsArray();
    
    Projection=ds.GetProjection()
    prj=osr.SpatialReference(wkt=ds.GetProjection())
    proj4_str=prj.GetAttrValue('AUTHORITY',1)
    
    m,n=data.shape    
    minx=gt[0]
    maxy=gt[3]
    maxx=minx+gt[1]*n  
    miny=maxy+gt[5]*m
    extent=(minx,maxx,miny,maxy)
    
    try:
        x=np.reshape(np.arange(gt[0],gt[0]+gt[1]*n,gt[1]),(1,n))
    except:
        x=np.reshape(np.arange(gt[0],gt[0]+gt[1]*n-gt[1],gt[1]),(1,n))
    try:
        y=np.reshape(np.flip(np.arange(gt[3]-m*gt[1],gt[3],gt[1])),(m,1))
    except:
        y=np.reshape(np.flip(np.arange(gt[3]-m*gt[1]+gt[1],gt[3],gt[1])),(m,1))
    x=np.tile(x,(m,1))
    y=np.tile(y,(1,n))
    
    z={'gt':gt,
       'Data':data,
       'X':x,
       'Y':y,
       'm':data.shape[0],
       'n':data.shape[1],
       'minx':minx,
       'maxx':maxx,
       'miny':miny,
       'maxy':maxy,
       'Extent':extent,
       'Projection':Projection,
       'Proj4_String':proj4_str}
    
    z=Bunch(z)
    
    return z

'''============================================================================
SAVE TO GEOTIFF
============================================================================'''

def SaveGeoTiff(z,fout):
    
    driver=gdal.GetDriverByName("GTiff")
    
    # Data type
    if z['Data'].dtype=='int8':
        dtype=gdal.GDT_Int8
    elif z['Data'].dtype=='int16':
        dtype=gdal.GDT_Int16    
    elif z['Data'].dtype=='int32':
        dtype=gdal.GDT_Int32
    elif z['Data'].dtype=='float32':
        dtype=gdal.GDT_Float32
    
    N_band=1
    ds_out=driver.Create(fout,z['n'],z['m'],N_band,dtype,[ 'COMPRESS=LZW' ])
    ds_out.SetProjection(z['Projection'])
    ds_out.SetGeoTransform(z['gt'])
    ds_out.GetRasterBand(1).WriteArray(z['Data'])    
    

'''============================================================================
CLIP RASTER
The raster input structure must be from the OpenGdal function from this module.
============================================================================'''

def ClipRaster(z_in,xlim,ylim):
    z=z_in.copy()
    z=Bunch(z)
    ix=np.where((z.X[0,:]>=xlim[0]) & (z.X[0,:]<=xlim[1]))[0]
    iy=np.where((z.Y[:,0]>=ylim[0]) & (z.Y[:,0]<=ylim[1]))[0]
    ind=np.ix_(iy,ix)
    z.Data=z.Data[ind]
    z.X=z.X[ind]
    z.Y=z.Y[ind]
    z.m=iy.size
    z.n=ix.size
    z.minx=np.min(z.X)
    z.maxx=np.max(z.X)
    z.miny=np.min(z.Y)
    z.maxy=np.max(z.Y)
    z.gt=(z.minx,z.gt[1],z.gt[2],z.maxy,z.gt[4],z.gt[5])
    z.Extent=(z.minx,z.maxx,z.miny,z.maxy)
    return z

'''============================================================================
REPROJECT RASTER
============================================================================'''

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


'''============================================================================
POLYGON AREA
============================================================================'''

def PolyArea(x,y):
    return 0.5*np.abs(np.dot(x,np.roll(y,1))-np.dot(y,np.roll(x,1)))


'''============================================================================
INPOLYGON
============================================================================'''

def InPolygon(xg,yg,xv,yv):
    
    l=[]
    for i in range(xv.size):
        l.append([xv[i],yv[i]])
    poly=Polygon(l)
    
    InPol=np.zeros(xg.shape)
    
    if xg.ndim==1:
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


'''============================================================================
EXTENT FROM GDAL
============================================================================'''

def GetExtentFromGDAL(gt,cols,rows):
    ''' Return list of corner coordinates from a geotransform

        @type gt:   C{tuple/list}
        @param gt: geotransform
        @type cols:   C{int}
        @param cols: number of columns in the dataset
        @type rows:   C{int}
        @param rows: number of rows in the dataset
        @rtype:    C{[float,...,float]}
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

'''============================================================================
COMPRESS CATEGORIES IN RASTER DATASET
============================================================================'''

def CompressCats(z0,lab0):
    uc=np.unique(z0)
    z1=uc[0]*np.ones(z0.shape)
    lab1=[]
    for i in range(uc.size):
        ind=np.where(z0==uc[i])
        ind2=np.ix_(ind[0],ind[1])
        z1[ind2[0].T,ind2[1]]=i+1
        lab1.append(lab0[uc[i]-1])
    lab1=np.array(lab1)
    return z1,lab1