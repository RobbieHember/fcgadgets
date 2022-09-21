
import os
from osgeo import gdal
import numpy as np
from osgeo import osr
from matplotlib import path
import copy
from shapely.geometry import Polygon,Point
import pyproj
import rasterio
from rasterio.features import shapes
from rasterio.warp import calculate_default_transform, reproject, Resampling
from rasterio.transform import from_origin
from fcgadgets.macgyver import utilities_general as gu

'''============================================================================
IMPORT SPATIAL REFERENCE SYSTEMS
============================================================================'''

def ImportSRSs():

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
    gt=ds.GetGeoTransform()

    # Gdal does not recognize negative values so using rasterio to get the actual
    # grid data
    #data=ds.ReadAsArray()
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
       'Proj4_String':proj4_str}

    #z=Bunch(z)

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
    ds=driver.Create(fout,z['n'],z['m'],N_band,dtype,[ 'COMPRESS=LZW' ])
    ds.SetProjection(z['Projection'])
    ds.SetGeoTransform(z['gt'])
    ds.GetRasterBand(1).WriteArray(z['Data'])

    return

'''============================================================================
CLIP RASTER
The raster input structure must be from the OpenGdal function from this module.
============================================================================'''

def ClipRasterByXYLimits(z_in,xlim,ylim):

    #xlim=z['xlim']
    #ylim=z['ylim']

    z=z_in.copy()

    ix=np.where((z['X'][0,:]>=xlim[0]) & (z['X'][0,:]<=xlim[1]))[0]
    iy=np.where((z['Y'][:,0]>=ylim[0]) & (z['Y'][:,0]<=ylim[1]))[0]

#    dx=ix.size-z['Data'].shape[1]
#    if dx==1:
#        ix=np.where((z['X'][0,:]>=xlim[0]) & (z['X'][0,:]<xlim[1]))[0]
#    elif dx==-1:
#        ix=np.where((z['X'][0,:]>=xlim[0]) & (z['X'][0,:]<=xlim[1]))[0]
#
#    dy=iy.size-z['Data'].shape[0]
#    if dy==1:
#        iy=np.where((z['Y'][:,0]>=ylim[0]) & (z['Y'][:,0]<ylim[1]))[0]
#    elif dy==-1:
#        iy=np.where((z['Y'][:,0]>=ylim[0]) & (z['Y'][:,0]<=ylim[1]+z['cellsize']))[0]

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

#%%

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
    else:
        pass

    ind_in=np.ix_(iy_in,ix_in)
    ind_ref=np.ix_(iy_ref,ix_ref)

    z={}
    z['Data']=0*z_ref['Data']
    z['Data'][ind_ref]=z_in['Data'][ind_in]

    for k in z_ref.keys():
        if (k=='Data'):
            continue
        z[k]=z_ref[k]

    return z

'''============================================================================
REPROJECT RASTER
============================================================================'''

def ReprojectRasterAndClipToRaster(fin,fout,fref,crs):

    # Grid to reproject
    ds=gdal.Open(fin,0)

    # Open reference grid
    z_ref=OpenGeoTiff(fref)

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

    z['Data']=z['Data'].astype('int16')
    SaveGeoTiff(z,fout)

    return z

#%%

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

def CompressCats(z0,id0,lab0,cl0):
    uc=np.unique(z0)
    z1=0*np.ones(z0.shape)
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

'''============================================================================
DIGITIZE (GET POLYGONS FROM BINARY MASK)
============================================================================'''

def Digitize(BinaryMask,xv,yv):

    #xv=xv.flatten()
    #yv=yv.flatten()

    s=shapes(BinaryMask.astype('int16'),mask=None,connectivity=4)

    xy=[]
    for i in range(10000):
        try:
            a=next(s)
            b=a[0]['coordinates'][0]
            d={}
            d['x']=np.array([])
            d['y']=np.array([])
            for j in range(len(b)):
                col=int(b[j][0])
                row=int(b[j][1])
                d['x']=np.append(d['x'],xv[col-1])
                d['y']=np.append(d['y'],yv[row-1])
            xy.append(d)
        except:
            break
    return xy

'''============================================================================
REVISE RASTER EXTENT
============================================================================'''

def ReviseRasterExtent(fin_ToAdjust,fin_Ref,fout):

    # specify input and output filenames
    #fin_ToAdjust=r'C:\Users\rhember\Documents\Data\BC1ha\Disturbances\VEG_CONSOLIDATED_CUT_BLOCKS_SP_year4.tif'
    #fout=r'C:\Users\rhember\Documents\Data\BC1ha\Disturbances\VEG_CONSOLIDATED_CUT_BLOCKS_SP_year4a.tif'
    #fin_Ref = r'C:\Users\rhember\Documents\Data\BC1ha\Admin\tsa.tif'

    # open reference file and get resolution
    ds=gdal.Open(fin_Ref,0)
    gt=ds.GetGeoTransform()
    data=ds.ReadAsArray()
    rows,cols=data.shape
    x_res=gt[1]
    y_res=-gt[5]  # make sure this value is positive
    extent=GetExtentFromGDAL(gt,cols,rows)

    # call gdal Warp
    kwargs={"format":"GTiff",
            "xRes":x_res,
            "yRes":y_res}
    #"outputBounds":extent

    ds=gdal.Warp(fout,fin_ToAdjust,
        creationOptions=["COMPRESS=LZW"],
        options=gdal.WarpOptions(options=['outputBounds'],outputBounds=extent)
        **kwargs)

    return

#def ReviseRasterExtent_old(fin1,fin2,fout1):
#
#    # Notes:
#    # 1) No data value specifically refers to values less than zero, set to -99
#    # 2) Set up to compress output
#
#    # Raster 1 (with extent that will be revised)
#    ds1=gdal.Open(fin1)
#    data1=ds1.ReadAsArray()
#    m1=ds1.RasterYSize
#    n1=ds1.RasterXSize
#
#    # Raster 2 (with target extent)
#    ds2=gdal.Open(fin2)
#    gt2=ds2.GetGeoTransform()
#    proj2=ds2.GetProjection()
#    data2=ds2.ReadAsArray()
#    m2=ds2.RasterYSize
#    n2=ds2.RasterXSize
#    pw=[gt2[0],gt2[3],gt2[0]+gt2[1]*n2,gt2[3]-gt2[1]*m2] # Target extent (ulx uly lrx lry)
#
#    # Adjusted raster 1 (with same extent as raster 2)
#    ds1_adj=gdal.Open(fin1)
#    ds1_adj=gdal.Translate('',ds1_adj,projWin=pw,format="MEM")
#    ds1_adj.SetGeoTransform(ds2.GetGeoTransform())
#    ds1_adj.SetProjection(proj2)
#    data1_adj=ds1_adj.ReadAsArray()
#    m1_adj=ds1_adj.RasterYSize
#    n1_adj=ds1_adj.RasterXSize
#
#    # Adjust the no data value
#    nodatval=-99
#    data1_adj=np.where((data1_adj<0),nodatval,data1_adj)
#
#    if fout1!='NoSave':
#        driver=gdal.GetDriverByName('GTiff')
#        outdata=driver.Create(fout1,n2,m2,1,gdal.GDT_UInt16,options=['COMPRESS=LZW'])
#        outdata.SetGeoTransform(ds1_adj.GetGeoTransform()) # sets same geotransform as input
#        outdata.SetProjection(ds1_adj.GetProjection()) # sets same projection as input
#        outdata.GetRasterBand(1).WriteArray(data1_adj)
#        outdata.GetRasterBand(1).SetNoDataValue(nodatval) # if you want these values transparent
#        outdata.FlushCache() # saves to disk
#        outdata=None
#
#    return

#%% Clip geodataframe to user-specified x and y limits

def ClipGDF(gdf_in,xlim,ylim):

    gdf=gdf_in.cx[xlim[0]:xlim[1],ylim[0]:ylim[1]]
    #gdf=gdf.reset_index(drop=True)
    #gdf=gpd.sjoin(gdf,roi['gdb']['bound'],how='left')

    # This grouped by may not be necessary - it prevents the file from working in overlays
    #gdf=gdf.groupby('index_right')

    return gdf

#%% Import Cities

def ImportCities(pthin):

    Cities=gu.ReadExcel(pthin)

    # Import spatial reference systems
    srs=ImportSRSs()
    Cities['X']=np.zeros(Cities['Lat'].size)
    Cities['Y']=np.zeros(Cities['Lat'].size)
    for i in range(Cities['Lat'].size):
        Cities['X'][i],Cities['Y'][i]=srs['Proj']['BC1ha'](Cities['Lon'][i],Cities['Lat'][i])

    return Cities