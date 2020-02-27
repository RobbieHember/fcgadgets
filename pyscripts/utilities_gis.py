
import gdal
import numpy as np
import osr
from matplotlib import path

'''============================================================================
BUNCH VARIABLES FROM DICTIONARY
============================================================================'''

class Bunch(dict):
    def __init__(self, *args, **kwds):
        super(Bunch, self).__init__(*args, **kwds)
        self.__dict__ = self

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

'''============================================================================
OPEN RASTER USING GDAL
============================================================================'''

def OpenGdal(pthin):
        
    ds=gdal.Open(pthin)
    
    gt=ds.GetGeoTransform(); 
    
    data=ds.ReadAsArray();
    
    Projection=ds.GetProjection()
    prj=osr.SpatialReference(wkt=ds.GetProjection())
    proj4_str=prj.GetAttrValue('AUTHORITY',1)
    
    m,n=data.shape
    
    minx=gt[0]; 
    
    maxy=gt[3]; 
    
    maxx=minx+gt[1]*n; 
    
    miny=maxy+gt[5]*m
    
    extent=(minx,maxx,miny,maxy)
    
    x=np.reshape(np.arange(gt[0],gt[0]+gt[1]*n,gt[1]),(1,n))
    y=np.reshape(np.flip(np.arange(gt[3]-m*gt[1],gt[3],gt[1])),(m,1))
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
CLIP RASTER
============================================================================'''

def ClipRaster(z,xlim,ylim):
    ix=np.where((z.X[0,:]>=xlim[0]) & (z.X[0,:]<=xlim[1]))[0]
    iy=np.where((z.Y[:,0]>=ylim[0]) & (z.Y[:,0]<=ylim[1]))[0]
    ind=np.ix_(ix,iy)
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
POLYGON AREA
============================================================================'''

def PolyArea(x,y):
    return 0.5*np.abs(np.dot(x,np.roll(y,1))-np.dot(y,np.roll(x,1)))

'''============================================================================
REVISE RASTER EXTENT
============================================================================'''

def ReviseRasterExtent(fin1,fin2,fout1):

    # Notes:
    # 1) No data value specifically refers to values less than zero, set to -99
    # 2) Set up to compress output
    
    # Raster 1 (with extent that will be revised)
    ds1=gdal.Open(fin1)
    data1=ds1.ReadAsArray()
    m1=ds1.RasterYSize
    n1=ds1.RasterXSize

    # Raster 2 (with target extent)
    ds2=gdal.Open(fin2)
    gt2=ds2.GetGeoTransform()
    proj2=ds2.GetProjection()
    data2=ds2.ReadAsArray()
    m2=ds2.RasterYSize
    n2=ds2.RasterXSize
    pw=[gt2[0],gt2[3],gt2[0]+gt2[1]*n2,gt2[3]-gt2[1]*m2] # Target extent (ulx uly lrx lry)

    # Adjusted raster 1 (with same extent as raster 2)
    ds1_adj=gdal.Open(fin1)
    ds1_adj=gdal.Translate('',ds1_adj,projWin=pw,format="MEM")
    ds1_adj.SetGeoTransform(ds2.GetGeoTransform())
    ds1_adj.SetProjection(proj2)
    data1_adj=ds1_adj.ReadAsArray()
    m1_adj=ds1_adj.RasterYSize
    n1_adj=ds1_adj.RasterXSize

    # Adjust the no data value
    nodatval=-99
    data1_adj=np.where((data1_adj<0),nodatval,data1_adj)

    if fout1!='NoSave':
        driver=gdal.GetDriverByName('GTiff')
        outdata=driver.Create(fout1,n2,m2,1,gdal.GDT_UInt16,options=['COMPRESS=LZW'])
        outdata.SetGeoTransform(ds1_adj.GetGeoTransform()) # sets same geotransform as input
        outdata.SetProjection(ds1_adj.GetProjection()) # sets same projection as input
        outdata.GetRasterBand(1).WriteArray(data1_adj)
        outdata.GetRasterBand(1).SetNoDataValue(nodatval) # if you want these values transparent
        outdata.FlushCache() # saves to disk
        outdata=None
    
    return


'''============================================================================
INPOLYGON
============================================================================'''

def InPolygon(xq, yq, xv, yv):
    shape = xq.shape
    xq = xq.reshape(-1)
    yq = yq.reshape(-1)
    xv = xv.reshape(-1)
    yv = yv.reshape(-1)
    q = [(xq[i], yq[i]) for i in range(xq.shape[0])]
    p = path.Path([(xv[i], yv[i]) for i in range(xv.shape[0])])
    return p.contains_points(q).reshape(shape)

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