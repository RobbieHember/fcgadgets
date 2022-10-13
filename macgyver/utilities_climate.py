'''

UTILITIES - CLIMATE DATA

'''

#%% Import modules

import os
import numpy as np
import numpy.matlib as mb
import pandas as pd
import geopandas as gpd
import fiona
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
import time
from scipy.interpolate import griddata
import pyproj
import netCDF4 as nc
import gzip
from scipy.io import loadmat
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

'''============================================================================

POTENTIAL EVAPOTRANSPIRATION

Input variables:
    Time - time vector [Year,Month]
    ta - monthly mean air temperature (degC)
    rswn - monthly mean net short wave radiation (MJ/m2/d)
    u - mean monthly wind speed (m/s) (defaults to 2.0 m/s if not supplied)
    Gs - mean monthly surface conductance (m/s) (required only for Method =
    Penman-Monteith)

Methods of calculation:
    1) Components
    2) Equilibrium rate
    3) Priestley-Taylor (1972)
    4) Penman (1948)
    5) Penman-Moneith

============================================================================'''

def GetETp(vi,Method,TimeInterval):

    # Dimensions (months, locations)
    try:
        N_m,N_s=vi['rswn'].shape
    except:
        N_m=vi['rswn'].size
        N_s=1

    # Number of years
    N_yr=int(N_m/12)

    # Wind speed - if not supplied, use default of 2.0 m s-1.
    if 'u' not in vi:
        vi['u']=2.0

    # Day of year
    #doyv=np.array([15,46,74,105,135,166,196,227,258,288,319,349])
    #if vi['Time'].shape[0]==1:
    #    doy=doyv[vi['Time'][:,1]]
    #else:
    #    doy=np.matlib.repmat(doyv,N_yr,N_s)

    # Constants
    con=HydroMetCon()

    # Convert net shortwave radiation from MJ m-2 d-1 to W m-2
    Rswn_conv=vi['rswn']*1e6/con['DayLength']

    # Convert net shortwave radiation (W m-2) to net radiation (W m-2)
    # Parameters from this were fitted to data at DF49 and agree roughly with
    # Landsberg's textbook.
    Rn=con['Rswn2Rn_Slope']*Rswn_conv+con['Rswn2Rn_AddOffset']

    # Psychrometric term (hPa K-1)
    Psychro=0.01*GetPsychrometric(vi['ta'],'Pressure')

    # Saturation vapour pressure/temperature slope (hPa K-1)
    Svps=0.01*GetSVPSlope(vi['ta'],'Pressure')

    if Method=='Components':

        # Radiative component (mm d-1) (McMahon et al. 2013)
        Eeq=((Svps/(Svps+Psychro))*Rn)*con['DayLength']/con['Lam']

        # Aerodynamic component (mm d-1) (McMahon et al. 2013)
        Ea=(1.31+1.38*vi['u'])*(vi['vpd']/10)

        # Add to tuple
        ETp=(Eeq,Ea)

    elif Method=='Equilibrium':

        # Radiative component (mm d-1) (McMahon et al. 2013)
        ETp=((Svps/(Svps+Psychro))*Rn)*con['DayLength']/con['Lam']

    elif Method=='Priestley-Taylor':

        # Radiative component (mm d-1) (McMahon et al. 2013)
        Eeq=((Svps/(Svps+Psychro))*Rn)*con['DayLength']/con['Lam']

        # Potential evaporatnspiration (Priestley and Taylor 1972) (mm d-1)
        ETp=con['Alpha_PT']*Eeq

    elif Method=='Penman':

        # Radiative component (mm d-1) (McMahon et al. 2013)
        Eeq=((Svps/(Svps+Psychro))*Rn)*con['DayLength']/con['Lam']

        # Aerodynamic component (mm d-1) (McMahon et al. 2013)
        Ea=(1.31+1.38*vi['u'])*(vi['vpd']/10)

        # Potential evapotranspiration (mm d-1) (Penman 1948)
        ETp=Eeq+(Psychro/(Svps+Psychro))*Ea

    elif Method=='Penman-Monteith':

        # Aerodynamic conductance (m s-1)
        if 'Ga' not in vi:
            vi['Ga']=0.01*vi['u']

        # Potential evapotranspiration from Penman-Monteith combination model (mm d-1)
        ETp=((Svps*Rn + con['RhoAir']*con['CpAir']*vi['vpd']*vi['Ga'])/(Svps+Psychro*(1+vi['Ga']/vi['Gs'])))/con['Lam']*con['DayLength']

    else:
        print('Method not recognized, quiting.')
        return ()

    # Unit conversion
    if (TimeInterval=='Month') | (TimeInterval=='M') | (TimeInterval=='m'):
        DIM=mb.repmat(np.reshape(con['DIM'],(-1,1)),N_yr,N_s)
        ETp=ETp*DIM

    return ETp

'''============================================================================

TADPOLE - MONTHLY SURFACE WATER BALANCE MODEL

Inputs variables (field names in vi structure):

    LAI, leaf area index (m2 m-2)
    lat, latitude (deg)
    ta, mean air temperature (degC)
    prcp, total precipitation (mm month-1)
    rswn, net shortwave radiation (MJ m-2 d-1)
    vpd, vapour pressure deficit (hPa)
    Ws, water content of soil (mm) ** Optional **
    Wsp, water content of snow (mm) ** Optional **

Outputs variables (field names in ov structure):

    ETp, potential evapotranspiration (mm month-1)
    ETa, actual evapotranspiration (mm month-1)
    Ws, water content of soil (mm)
    Wsp, water content of snow (mm)
    M, snowmelt(mm month-1)
    R, runoff (mm month-1)
    RSF, rain to snow ratio (dim)

Mode of operation:

If you are working with a small number of locations, use
"Combined" which runs all spatial units simultaneously. If running with large
grids, use "Grid" which inputs and outputs one month at a time. In the latter
approach, the state variables need to be initialized for the first time step,
and included as input arguments (from the output of t-1) for subsequent time
steps.

Parameters required in meta structure:

    Ws_max: water storage capacity (200 mm)
    Tmin: minimum rain/snow fraction temperature (-3 degC)
    Tmax: maximum rain/snow fraction temperature (-3 degC)
    Daily Interval: daily soil water computaiton interval (5)

Warnings:

    1) This code will not tolerate NaNs.
    2) The model requires spin-up for at least one year.

============================================================================'''

def Tadpole(par,vi):

    #--------------------------------------------------------------------------
    # Set up output arguments
    #--------------------------------------------------------------------------

    vo={}

    # Size of input variables
    N_mo,N_s=vi['ta'].shape

    if par['Method']=='Combined':

        # Running all time steps for a set of n locations
        vo['Ws']=par['Ws_max']*np.ones((N_mo,N_s))
        vo['Wsp']=np.zeros((N_mo,N_s))
        vo['ETp']=np.zeros((N_mo,N_s))
        vo['ETa']=np.zeros((N_mo,N_s))
        vo['R']=np.zeros((N_mo,N_s))
        vo['M']=np.zeros((N_mo,N_s))

    elif par['Method']=='Grid':

        # If running one month grid at a time, set initial values beased on input
        # from last run

        # Initialize state variables with maximum levels and add previous time
        # step for time steps beyond the first one

        if 'Ws' in vi:
            vo['Ws'][0,:]=vi['Ws']
            vo['Wsp'][0,:]=vi['Wsp']
        else:
            vo['Ws']=par['Ws_max']*np.ones((1,N_mo))
            vo['Wsp']=np.zeros((1,N_mo))

        vo['ETp']=np.zeros((1,N_mo))
        vo['ETa']=np.zeros((1,N_mo))
        vo['R']=np.zeros((1,N_mo))
        vo['M']=np.zeros((1,N_mo))

    #--------------------------------------------------------------------------
    # Constants
    #--------------------------------------------------------------------------

    con=HydroMetCon()

    #--------------------------------------------------------------------------
    # Potential evapotranspiratiom (mm month-1)
    #--------------------------------------------------------------------------

    vo['ETp']=GetETp(vi,par['ETp Method'],'Month')

    #--------------------------------------------------------------------------
    # Canopy interception as a function of leaf area (Landsberg et al. 2005)
    # Maximum proportion of rainfall evaporated from canopy MaxIntcptn – 0.15
    # Leaf area index for maximum rainfall interception LAImaxIntcptn – 5
    #--------------------------------------------------------------------------

    # Fraction of intercepted precipitation
    FracPrecipInt=par['Ei_FracMax']*np.minimum(1.0,vi['LAI']/par['Ei_ALMax'])

    # Potential evaporation of intercepted precipiration (mm month-1)
    Ei_Potential=FracPrecipInt*vi['prcp']

    # Actual evaporation of intercepted precipiration, defined as the minimum
    # between the energy-limited rate, and the water-limited rate (mm month-1)
    Ei_Actual=np.minimum(vo['ETp'],Ei_Potential)

    # Transpiration at energy-limited rate, i.e. prior to adding surface
    # constraints (mm month-1)
    Et_EnergyLimited=vo['ETp']-Ei_Actual

    # Throughfall (mm month-1)
    P_Throughfall=vi['prcp']-Ei_Actual

    # Partition total precipitation into solid and liquid components
    fT=(vi['ta']-par['Tmin'])/(par['Tmax']-par['Tmin'])
    fT=np.minimum(np.maximum(0,fT),1)
    Pr=fT*P_Throughfall
    Ps=P_Throughfall-Pr

    #--------------------------------------------------------------------------
    # Loop through months
    #--------------------------------------------------------------------------

    for iT in range(N_mo):

        #----------------------------------------------------------------------
        # Set inititial daily water pools at level sfromt he end of the last
        # month
        #----------------------------------------------------------------------

        if par['Method']=='Combined':

            if iT==0:
                Ws_d=vo['Ws'][iT,:]
                Wsp_d=vo['Wsp'][iT,:]
            else:
                Ws_d=vo['Ws'][iT-1,:]
                Wsp_d=vo['Wsp'][iT-1,:]

        elif par['Method']=='Grid':

            Ws_d=vo['Ws']
            Wsp_d=vo['Wsp']

        #----------------------------------------------------------------------
        # Daily fluxes (mm d-1)
        #----------------------------------------------------------------------

        Ps_d=Ps[iT,:]/(30/par['Daily_Interval'])
        Pr_d=Pr[iT,:]/(30/par['Daily_Interval'])
        Ei_Actual_d=Ei_Actual[iT,:]/(30/par['Daily_Interval'])
        Et_EnergyLimited_d=Et_EnergyLimited[iT,:]/(30/par['Daily_Interval'])

        for iDay in range(0,30,par['Daily_Interval']):

            # Potential snowmelt (mm d-1), equation from Thornthwaite and Mather (1955)
            M_d=2.63+2.55*vi['ta'][iT,:]+0.0912*vi['ta'][iT,:]*Pr_d

            # Actual snowmelt (mm d-1)
            M_d=np.maximum(np.zeros((1,N_s)),np.minimum(M_d,Wsp_d+Ps_d))

            # Cumulative snowmelt (mm)
            vo['M'][iT,:]=vo['M'][iT,:]+M_d

            # Updated snowpack water content (mm)
            Wsp_d=Wsp_d+Ps_d-M_d

            # Update soil water content (mm)
            Ws_d=Ws_d+M_d+Pr_d

            # Supply function (Willmott et al. 1985)
            # x=[0:0.1:1]' plot(x,1-exp(-6.68*(x)),'ko')
            fWs=np.minimum(1,np.maximum(0,1-np.exp(-6.68*(Ws_d/par['Ws_max']))))

            # Actual evaporation (mm d-1)
            Et_Actual_d=fWs*Et_EnergyLimited_d

            # Cumulative actual evapotranspiration as the sum of wet-canopy
            # evaporation and transpiration (mm)
            vo['ETa'][iT,:]=vo['ETa'][iT,:]+Ei_Actual_d+Et_Actual_d

            # Remove transpiration from soil water pool (mm)
            Ws_d=Ws_d-Et_Actual_d

            # Find any spatial units where soil water exceeds capacity and add
            # "surplus" to monthly runoff and restrict soil water content
            # to capacity
            R_d=np.maximum(0,Ws_d-par['Ws_max'])
            Ws_d=np.minimum(Ws_d,par['Ws_max'])

            # Update monthly runoff (mm)
            vo['R'][iT,:]=vo['R'][iT,:]+R_d

        # Update water pools
        vo['Ws'][iT,:]=Ws_d
        vo['Wsp'][iT,:]=Wsp_d

    #--------------------------------------------------------------------------
    # Constrain pools to be positve
    #--------------------------------------------------------------------------

    vo['Ws']=np.maximum(0,vo['Ws'])
    vo['Wsp']=np.maximum(0,vo['Wsp'])

    #--------------------------------------------------------------------------
    # Include rainfall fraction
    #--------------------------------------------------------------------------

    if par['Include Rainfall Fraction']=='Yes':
        vo['RF']=np.minimum(1,np.maximum(0,Pr/Ps))

    return vo

'''============================================================================

HYDROMETEOROLOGY CONSTANTS

============================================================================'''

def HydroMetCon():

    con={}

    # Air density (kg m-3)
    con['RhoAir']=1.2

    # Specific heat capacity (J kg-1 K-1)
    con['CpAir']=1010

    # Latent heat of vaporization (J kg-1)
    con['Lam']=2460000

    # Preistly taylor coefficient
    con['Alpha_PT']=1.26

    # Daylength (s)
    con['DayLength']=86400

    # Conversion of downwelling short wave solar radiation to net radiation
    con['Rswd2Rn_Slope']=11.96
    con['Rswd2Rn_AddOffset']=3.46

    # Convert net shortwave radiation (W m-2) to net radiation (W m-2)
    # Parameters from this were fitted to data at DF49 and agree roughly with
    # Landsberg's textbook.
    con['Rswn2Rn_Slope']=0.837
    con['Rswn2Rn_AddOffset']=-23.58

    # Albedo
    Albedo={}
    Albedo['Forest Coniferous']=0.04
    Albedo['Forest Deciduous']=0.09
    con['Albedo']=Albedo

    # Days in month
    con['DIM']=np.array([31,28,31,30,31,30,31,31,30,31,30,31])

    return con

'''============================================================================

CALCULATE PSYCHROMETRIC TERM

============================================================================'''

def GetPsychrometric(ta,Units):

    # This method compares closely with Fernandes et al. (2007)

    if Units=='Pressure':
        # Units: Pa K=1
        b=np.array([-6.02240896358241e-005,0.0616092436974788,64.9608123249299])
    elif Units=='Density':
        # Units: kg m-3 K-1
        b=np.array([4.2507002801123e-009,-1.40523109243698e-006,0.000515354446778711])

    y=b[0]*ta**2 + b[1]*ta+b[2]

    return y

'''============================================================================

CALCULATE SATURATION VAPOUR PRESSURE VS. TEMPERATURE SLOPE

============================================================================'''

def GetSVPSlope(ta,Units):

    # Calculate slope of the saturation vapour pressure - temperature curve
    # using the method METH based on input of air temperature (deg C).

    if Units=='Pressure':
        # Units: Pa K=1
        b=np.array([0.000011482039374,0.001256498041862,0.078296395471144,2.846599268013546,44.494094675319538])
    elif Units=='Density':
        # Units: kg m-3 K-1
        b=np.array([4.5645788679845e-011,7.49675848747055e-009,5.23844963086449e-007,2.00663120848879e-005,0.000335075613241248])

    y=b[0]*ta**4 + b[1]*ta**3 + b[2]*ta**2 + b[3]*ta + b[4]

    return y