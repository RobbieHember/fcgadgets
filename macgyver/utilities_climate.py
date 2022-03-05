'''

UTILITIES - CLIMATE DATA

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
from scipy.interpolate import griddata
import pyproj
import netCDF4 as nc
import gzip
from scipy.io import loadmat
from shapely.geometry import Polygon,Point
from fcgadgets.macgyver import utilities_general as gu
from fcgadgets.macgyver import utilities_gis as gis
from fcgadgets.macgyver import utilities_inventory as invu
from fcgadgets.cbrunner import cbrun_utilities

#%% Import historical monthly normals

def Import_BC1ha_Historical_Normals(meta,geos):
    
    # Initialize data structure
    data={}
    
    # Path to data
    pth=r'C:\Users\rhember\Documents\Data\BC1ha\Climate\Monthly'
    
    # List of variables
    namV=['tmean','rswd','prcp','vpd']
    
    # Import coordinates of BC1ha grid
    zTSA=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Admin\tsa.tif')
    
    X=zTSA['X'].astype(float)[0,:]
    Y=zTSA['Y'].astype(float)[:,0]
    
    # Calculate indices
    idx2bc1ha=np.zeros((geos['Sparse']['X'].size,2),dtype=int)
    for iS in range(idx2bc1ha.shape[0]):
        adx=np.abs(geos['Sparse']['X'][iS]-X)
        ady=np.abs(geos['Sparse']['Y'][iS]-Y)
        idx2bc1ha[iS,0]=int(np.where(adx==np.min(adx))[0][0])
        idx2bc1ha[iS,1]=int(np.where(ady==np.min(ady))[0][0])
    
    for iV in range(len(namV)):
        
        data[namV[iV]]=np.zeros((12,geos['Sparse']['X'].size))
        
        # Load data
        for mo in range(12):
            
            fin=pth + '\\BC1ha_' + namV[iV] + '_mon_norm_1971to2000_si_hist_v1\\BC1ha_' + namV[iV] + '_mon_norm_1971to2000_si_hist_v1_' + str(mo+1) + '.mat'
                   
            z0=loadmat(fin,squeeze_me=True)
            
            iDat=np.where(np.asarray(z0['z'].dtype.names)=='Data')[0][0]
            iSF=np.where(np.asarray(z0['z'].dtype.names)=='ScaleFactor')[0][0]
            z1=z0['z'][()][iDat].astype(float)*z0['z'][()][iSF]
            
            for iS in range(idx2bc1ha.shape[0]):
                data[namV[iV]][mo,iS]=z1[idx2bc1ha[iS,1],idx2bc1ha[iS,0]]
    
            # Fix missing data
            iBad=np.where(data[namV[iV]][mo,:]<=-90)[0]
            iGood=np.where(data[namV[iV]][mo,:]>-90)[0]
            data[namV[iV]][mo,:][iBad]=np.mean(data[namV[iV]][mo,:][iGood])
    
    return data

#data=Import_BC1ha_Historical_Normals(meta,geos)

#plt.plot(np.mean(data['tmean'],axis=1))

#%% Import CRU TS

def Import_CRU(meta,geos):
    
    # Projection info
    srs=gis.ImportSRSs()
    
    # Initialize data structure
    data={}
    
    # Time vector
    tv=gu.tvec('m',1901,2020)
    
    # List of variables
    namV=['tmn','tmx','pre','vap']
    
    for iV in range(len(namV)):
    
        data[namV[iV]]=np.zeros((tv.shape[0],geos['Sparse']['X'].size),dtype=np.float)
        
        # Path to data
        fin=r'C:\Users\rhember\Documents\Data\CRU\Version405\cru_ts4.05.1901.2020.' +  namV[iV] + '.dat.nc.gz'
        
        # Import netcdf data
        z0={}
        with gzip.open(fin) as gz:
            with nc.Dataset('dummy',mode='r',memory=gz.read()) as ds:
                for k in ds.variables.keys():
                    z0[k]=np.array(ds.variables[k][:])
        
        # Isolate BC area  (so that interpolation is faster)
        indX=np.where( (z0['lon']>=-143) & (z0['lon']<=-112) )[0]
        indY=np.where( (z0['lat']>=47) & (z0['lat']<=61) )[0]
        lon,lat=np.meshgrid(z0['lon'][indX],z0['lat'][indY],sparse=False) # ,indexing='ij'    
        x1,y1=srs['Proj']['BC1ha'](lon.flatten(),lat.flatten())
        z1=z0[namV[iV]][:,indY,:][:,:,indX]
        
        # Extract time interval
        for iT in range(tv.shape[0]):
            
            print(namV[iV] + ' ' + str(tv[iT,0]))
            
            # Get data for year and month
            iT=np.where( (tv[:,0]==tv[iT,0]) & (tv[:,1]==tv[iT,1]) )[0]            
            z1t=np.squeeze(z1[iT,:,:]).flatten()
            
            # Remove bad data
            iGood=np.where( np.abs(z1t)<9000 )
            x2=x1[iGood]
            y2=y1[iGood]
            z1t=z1t[iGood]
        
            # Interpolate
            z2=griddata( (x2,y2) ,z1t, (geos['Sparse']['X'],geos['Sparse']['Y']) ,method='linear')
        
            # QA
            flg=0
            if flg==1:
                Mask=geos['Data'].copy()
                Mask[geos['iMask']]=z2    
                plt.close('all')
                plt.matshow(Mask,vmin=-20,vmax=10)
            
            # Add to data structure
            data[namV[iV]][iT,:]=z2
    
    return data

#data=Import_CRU(meta,geos)

#plt.close('all')
#plt.plot(np.mean(data['tmn'],axis=1))


#%% Import NCEP reanlysis global

def Import_NCEP_Reanalysis_Global(meta,geos):
    
    # Projection info
    srs=gis.ImportSRSs()
    
    # Initialize data structure
    data={}
    
    # Time vector
    tv=gu.tvec('m',1951,2021)
    
    # List of variables
    namV=['dswrf']
    
    for iV in range(len(namV)):
    
        data[namV[iV]]=np.zeros((tv.shape[0],geos['Sparse']['X'].size),dtype=np.float)
        
        # Path to data
        fin=r'C:\Users\rhember\Documents\Data\Reanlysis\NCEP Global\dswrf.sfc.mon.mean.nc'
        
        # Import netcdf data        
        ds=nc.Dataset(fin)
        z0={}
        for k in ds.variables.keys():
            z0[k]=np.array(ds.variables[k][:])
        
        # Fix longitutude
        #latR=np.arange(-88.542,88.542,94)'; 
        #lonR=linspace(0,358.125,192)'; ind=find(lonR>180); lonR(ind)=lonR(ind)-360;
        ind=np.where(z0['lon']>=180)[0]
        z0['lon'][ind]=-1*(np.max(z0['lon'])-z0['lon'][ind])
        
        # Isolate BC area  (so that interpolation is faster)
        indX=np.where( (z0['lon']>=-143-5) & (z0['lon']<=-112+5) )[0]
        indY=np.where( (z0['lat']>=47-5) & (z0['lat']<=61+5) )[0]
        lon,lat=np.meshgrid(z0['lon'][indX],z0['lat'][indY],sparse=False) # ,indexing='ij'    
        x1,y1=srs['Proj']['BC1ha'](lon.flatten(),lat.flatten())
        z1=z0[namV[iV]][:,indY,:][:,:,indX]
        
        # Convert from given W/m2 to MJ m-2 d-1 
        z1=z1*86400*1e-6
        
        # Extract time interval
        for iT in range(tv.shape[0]):
            
            print(namV[iV] + ' ' + str(tv[iT,0]))
            
            # Get data for year and month
            iT=np.where( (tv[:,0]==tv[iT,0]) & (tv[:,1]==tv[iT,1]) )[0]            
            z1t=np.squeeze(z1[iT,:,:]).flatten()
            
            # Remove bad data
            iGood=np.where( np.abs(z1t)<9000 )
            x2=x1[iGood]
            y2=y1[iGood]
            z1t=z1t[iGood]
        
            # Interpolate
            z2=griddata( (x2,y2) ,z1t, (geos['Sparse']['X'],geos['Sparse']['Y']) ,method='linear')
        
            # QA
            flg=0
            if flg==1:
                Mask=geos['Data'].copy()
                Mask[geos['iMask']]=z2    
                plt.close('all')
                plt.matshow(Mask,vmin=-20,vmax=10)
            
            # Add to data structure
            data[namV[iV]][iT,:]=z2
    
    return data

#data=Import_NCEP_Reanalysis_Global(meta,geos)

#plt.close('all')
#plt.plot(np.mean(data['dswrf'],axis=1))


#%% Import monthly CMIP6 native grids and put them in project spatial reference system

def Import_CMIP6(meta,geos):

    # Path to CMIP6 data
    PathCMIP6=r'C:\Users\rhember\Documents\Data\CMIP6'
    
    # Time vectors
    tv=gu.tvec('m',1850,2100)
    tv_h=gu.tvec('m',1850,2014)
    tv_f=gu.tvec('m',2015,2100)
    #tva=np.arange(1850,2101,1)
    tva_h=np.arange(1850,2015,1)
    tva_f=np.arange(2015,2101,1)
    
    # Days in month
    dim=np.repeat(np.array([31,28,31,30,31,30,31,31,30,31,30,31]),tv.shape[0]/12)
    
    # Model names
    ms=gu.ReadExcel(PathCMIP6 + '\\ModelSummary.xlsx')
    
    # Exclusion of models 
    namM=ms['Model'][(ms['Model']!='CESM2-WACCM') & (ms['Model']!='EC-Earth3')]
    fsM=ms['File structure'][(ms['Model']!='CESM2-WACCM') & (ms['Model']!='EC-Earth3')]
    
    # Future scenarios
    namS=['ssp245','ssp585']
    
    # Variable list
    namV=['tasmin','tasmax','pr','rsds']
    
    # Scale factors
    sfV=[10,10,1,10]
    
    # Conversion factors
    sec2days=86400
    dk2dc=-273.15
    
    # Projection info
    srs=gis.ImportSRSs()
    
    data={}
    
    for iV in range(len(namV)):
        #iV=0
        
        nV=namV[iV]
        
        data[nV]={}
        
        for iM in range(namM.size):
        #for iM in range(1):        
            #iM=0
            
            nM=namM[iM]
            
            print(nM + '_' + nV)
            
            data[nV][nM]={}
               
            #----------------------------------------------------------------------
            # Determine what runs exist
            #----------------------------------------------------------------------
            
            namR=np.array([],dtype=np.int16)
            cnt=0
            for iR in range(0,100):
                if fsM[iM]=='Grouped':
                    pin=PathCMIP6 + '\\' + nM + '\\' + nV + '_Amon_' + nM + '_ssp585_r' + str(iR+1) + 'i1p1f1_gn_201501-210012.nc'
                else:
                    pin=PathCMIP6 + '\\' + nM + '\\' + nV + '_Amon_' + nM + '_ssp585_r' + str(iR+1) + 'i1p1f1_gn_201501-201512.nc'
                if os.path.isfile(pin)==True:
                    namR=np.append(namR,iR+1)
                    cnt=cnt+1
                if cnt==2:
                    break
            
            #----------------------------------------------------------------------
            # Loop through
            #----------------------------------------------------------------------
            
            for iR in range(namR.size):
                #iR=0
                
                nR=namR[iR]
                
                data[nV][nM][nR]={}
                
                # Initialize for each future scenario
                for iS in range(len(namS)):
                    data[nV][nM][nR][namS[iS]]=np.zeros( (tv.shape[0],geos['Sparse']['X'].size) ,dtype=np.int16 )
                
                if ms['File structure'][iM]=='Grouped':
                    
                    #--------------------------------------------------------------
                    # Historical
                    #--------------------------------------------------------------
                    
                    z0=gis.OpenGeoTiff(PathCMIP6 + '\\' + nM + '\\' + nV + '_Amon_' + nM + '_historical_r' + str(nR) + 'i1p1f1_gn_185001-201412.nc')   
                    x0,y0=srs['Proj']['BC1ha'](z0['X'].flatten(),z0['Y'].flatten())
                    
                    for iT in range(z0['Data'].shape[0]):
                        #print(tv_h[iT,0])
                        indT=np.where( (tv[:,0]==tv_h[iT,0]) & (tv[:,1]==tv_h[iT,1]) )[0]
                        
                        z0t=z0['Data'][iT,:,:].flatten()
                        
                        z1=griddata( (x0,y0) ,z0t, (geos['Sparse']['X'],geos['Sparse']['Y']) ,method='linear')
                    
                        for iS in range(len(namS)):
                            
                            nS=namS[iS]
                            if nV=='pr':
                                data[nV][nM][nR][nS][indT,:]=z1*sec2days*dim[indT]*sfV[iV]
                            elif (nV=='tasmin') | (nV=='tasmax'):
                                data[nV][nM][nR][nS][indT,:]=(z1+dk2dc)*sfV[iV]
                            elif (nV=='rsds'):
                                data[nV][nM][nR][nS][indT,:]=(z1*sec2days*1e-6)*sfV[iV]
    
    #                mu=np.zeros(tva.size)
    #                for iTa in range(tva.size):
    #                    iTm=np.where(tv[:,0]==tva[iTa])[0]
    #                    tmp=np.mean(data[nV][nM][nR][nS][iTm,:,:].astype(float)/sfV[iV],axis=0)
    #                    mu[iTa]=np.mean(tmp)
    #                
    #                plt.close('all')
    #                plt.plot(tva,mu-273.15)
                    
                    #--------------------------------------------------------------
                    # Future scenarios
                    #--------------------------------------------------------------
                    
                    for iS in range(len(namS)):
                        
                        nS=namS[iS]
                        
                        z0=gis.OpenGeoTiff(PathCMIP6 + '\\' + nM + '\\' + nV + '_Amon_' + nM + '_' + nS + '_r' + str(nR) + 'i1p1f1_gn_201501-210012.nc')   
                        x0,y0=srs['Proj']['BC1ha'](z0['X'].flatten(),z0['Y'].flatten())
                        
                        for iT in range(z0['Data'].shape[0]):
                            #print(tv_f[iT,0])
                            
                            indT=np.where( (tv[:,0]==tv_f[iT,0]) & (tv[:,1]==tv_f[iT,1]) )[0]
                            
                            z0t=z0['Data'][iT,:,:].flatten()
                            
                            z1=griddata( (x0,y0) ,z0t, (geos['Sparse']['X'],geos['Sparse']['Y']) ,method='linear')
                        
                            if nV=='pr':
                                data[nV][nM][nR][nS][indT,:]=z1*sec2days*dim[indT]*sfV[iV]
                            elif (nV=='tasmin') | (nV=='tasmax'):
                                data[nV][nM][nR][nS][indT,:]=(z1+dk2dc)*sfV[iV]
                            elif (nV=='rsds'):
                                data[nV][nM][nR][nS][indT,:]=(z1*sec2days*1e-6)*sfV[iV]
                    
                else:
                    
                    # Annual file structure
                    
                    #--------------------------------------------------------------
                    # Historical
                    #--------------------------------------------------------------
                    
                    for iT in tva_h:
                        
                        z0=gis.OpenGeoTiff(PathCMIP6 + '\\' + nM + '\\' + nV + '_Amon_' + nM + '_historical_r' + str(nR) + 'i1p1f1_gn_' + str(iT) + '01-' + str(iT) + '12.nc')
                        if iT==tva_h[0]:
                            x0,y0=srs['Proj']['BC1ha'](z0['X'].flatten(),z0['Y'].flatten())
                        
                        for mo in range(12):
                            
                            indT=np.where( (tv[:,0]==iT) & (tv[:,1]==mo+1) )[0]
                            
                            z0t=z0['Data'][mo,:,:].flatten()
                        
                            z1=griddata( (x0,y0) ,z0t, (geos['Sparse']['X'],geos['Sparse']['Y']) ,method='linear')
                    
                            for iS in range(len(namS)):
                        
                                nS=namS[iS]
                                
                                if nV=='pr':
                                    data[nV][nM][nR][nS][indT,:]=z1*sec2days*dim[indT]*sfV[iV]
                                elif (nV=='tasmin') | (nV=='tasmax'):
                                    data[nV][nM][nR][nS][indT,:]=(z1+dk2dc)*sfV[iV]
                                elif (nV=='rsds'):
                                    data[nV][nM][nR][nS][indT,:]=(z1*sec2days*1e-6)*sfV[iV]
                    
                    #--------------------------------------------------------------
                    # Future
                    #--------------------------------------------------------------
                    
                    for iS in range(len(namS)):
                        
                        nS=namS[iS]
                    
                        for iT in tva_f:
                            
                            z0=gis.OpenGeoTiff(PathCMIP6 + '\\' + nM + '\\' + '\\' + nV + '_Amon_' + nM + '_' + nS + '_r' + str(nR) + 'i1p1f1_gn_' + str(iT) + '01-' + str(iT) + '12.nc')
                            if iT==tva_f[0]:
                                x0,y0=srs['Proj']['BC1ha'](z0['X'].flatten(),z0['Y'].flatten())
                            
                            for mo in range(12):
                                
                                indT=np.where( (tv[:,0]==iT) & (tv[:,1]==mo+1) )[0]
                                
                                z0t=z0['Data'][mo,:,:].flatten()
                            
                                z1=griddata( (x0,y0) ,z0t, (geos['Sparse']['X'],geos['Sparse']['Y']) ,method='linear')
                        
                                if nV=='pr':
                                    data[nV][nM][nR][nS][indT,:]=z1*sec2days*dim[indT]*sfV[iV]
                                elif (nV=='tasmin') | (nV=='tasmax'):
                                    data[nV][nM][nR][nS][indT,:]=(z1+dk2dc)*sfV[iV]
                                elif (nV=='rsds'):
                                    data[nV][nM][nR][nS][indT,:]=(z1*sec2days*1e-6)*sfV[iV]
    
        # Make directory
        #if not os.path.exists(meta['Paths']['Inputs'] + '\\Climate'):
        #    os.makedirs(meta['Paths']['Inputs'] + '\\Climate')
    
        # Save data
        #gu.opickle(meta['Paths']['Inputs'] + '\\Climate\\cmip6_' + nV + '.pkl',data)
    
    return data


# Import CMIP6 data
#data=Import_CMIP6(meta,geos)

# Calculate mean air temperature
#data1={}
#data1['tmean']

#%%

#nV='tasmax'
#nM=ms['Model'][1]
#nS='ssp585'
#nR=1
#
#mua=np.zeros(tva.size)
#for i in range(tva.size):
#    iT=np.where( (tv[:,0]==tva[i]) & (tv[:,1]>=5) & (tv[:,1]<=9) )[0]
#    z=data[nV][nM][nR][nS][iT,:,:].astype(float)/10+273.15
#    mua[i]=np.mean(z)
#    
#plt.close('all')
#plt.plot(tva,mua/10,'-ko')
#

