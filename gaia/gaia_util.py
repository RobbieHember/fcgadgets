import os
import numpy as np
import copy
import matplotlib.pyplot as plt
import numpy.matlib as mb
import pyproj
import rasterio
from rasterio import features
from rasterio.transform import from_origin
from scipy.interpolate import griddata
from scipy.optimize import curve_fit
import fcgadgets.macgyver.util_general as gu
import fcgadgets.macgyver.util_gis as gis
import fcgadgets.bc1ha.bc1ha_util as u1ha

#%% Chapman-Richards difference equation
def fA(x,a,b,c):
	y=a*b*c*np.exp(-b*x)*(1-np.exp(-b*x))**(c-1)
	return y

#%% Logarithmic function of atmospheric carbon dioxide concentration (ppm)
def fCO2(x,b):
	y=b[0]*(1+b[1]*np.log(x/b[2]))
	return y

#%% Low N response (based on consideration of deficiency-sufficience phase transition assumption
def fNDEP(x,b):
	y=np.max(0,1-np.exp(-b[0]*(x-b[1])))
	return y

#%% Age response of growth
def AgeResponseFitTIPSY(meta):
	xhat=np.arange(0,200)
	flg=0
	if flg==1:
		b0=[140,0.032,4.4]
		y=fA(xhat,b0[0],b0[1],b0[2])
		plt.close('all'); plt.plot(xhat,y,'b-')
	dGY=gu.ReadExcel(r'C:\Data\TIPSY\tipsy_sw_carbon_plfd_si20_sph1600_oafdef.xlsx')
	p,pcov=curve_fit(fA,dGY['Age'],dGY['Gsw'],b0)
	y=fA(xhat,p[0],p[1],p[2])
	plt.close('all')
	plt.plot(dGY['Age'],dGY['Gsw'],'b-')
	plt.plot(xhat,y,'g--')
	return

#%%
def CreateRefGrid_BC5k(meta):
	# Manually copy the Land Mask from the bc1ha DB and then rescale it to 5k
	fin=meta['Paths']['bc5k Ref Grid']
	gis.ResampleRaster(fin,0.02)
	#Check that it worked: z=gis.OpenGeoTiff(meta['Paths']['bc5k Ref Grid'])
	return

#%% Import environmental data
def ImportEnvironment(meta,tv):
	d={'co2':np.zeros(tv.size),'ndep':np.zeros(tv.size)}

	# Import CO2 (ppm), source: https://www.pik-potsdam.de/~mmalte/rcps/ (ppm)
	dC=gu.ReadExcel(meta['Paths']['DB']['CO2'])
	for iT in range(tv.size):
		ind=np.where(dC['Year']==tv[iT])[0]
		d['co2'][iT]=dC['CO2 RCP45'][ind]

	# Import N deposition
	zRef5=gis.OpenGeoTiff(meta['Paths']['bc5k Ref Grid'])
	iLand=np.where(zRef5['Data']==1)
	dN=gu.ipickle(meta['Paths']['DB']['NDEP'] + '\\ISIMIP\\ndep_60.pkl')
	dN['ndep']['Data']=dN['ndep']['Data'].astype('float')*dN['ndep']['Scale Factor']
	for iT in range(tv.size):
		indT=np.where(dN['ndep']['tv']==tv[iT])[0]
		if indT.size==0:
			continue
		d['ndep'][iT]=np.mean(dN['ndep']['Data'][indT[0],:,:][iLand])
	indT=np.where(tv<dN['ndep']['tv'][0])[0]
	d['ndep'][indT]=np.mean(dN['ndep']['Data'][0,:,:][iLand])
	indT=np.where(tv>dN['ndep']['tv'][-1])[0]
	d['ndep'][indT]=np.mean(dN['ndep']['Data'][-1,:,:][iLand])
	return d

#%% Import nitrogen deposition from ISIMIP
def Process_Ndep_ISIMIP(meta):
	tv=np.arange(1850,2101,1)
	zRef=gis.OpenGeoTiff(meta['Paths']['bc5k Ref Grid'])
	srs=gis.ImportSRSs()
	lon,lat=pyproj.transform(srs['Proj']['BC1ha'],srs['Proj']['Geographic'],zRef['X'],zRef['Y'])
	rcpL=['26','60','85']
	for rcp in rcpL:
		d={'ndep':{}}
		d['ndep']['tv']=tv
		d['ndep']['Data']=np.zeros((tv.size,zRef['m'],zRef['n']),dtype='float')
		d['ndep']['TS Mean']=np.zeros(tv.size)
		d['ndep']['Scale Factor']=0.01
	
		tv0=np.arange(1861,2006,1)
		d0=gu.ReadNC(meta['Paths']['DB']['NDEP'] + '\\ISIMIP\\ndep_nhx_histsoc_annual_1861_2005.nc4')
		print(np.mean(d0['nhx'][140,:,:]))
		d1=gu.ReadNC(meta['Paths']['DB']['NDEP'] + '\\ISIMIP\\ndep_noy_histsoc_annual_1861_2005.nc4')
		#print(np.mean(d1['noy'][140,:,:]))
		d0['ntot']=10*(d0['nhx']+d1['noy'])
		lat0,lon0=np.meshgrid(d0['lat'],d0['lon'])
		#print(np.mean(d0['ntot'][140,:,:]))
		#plt.close('all'); plt.matshow(d0['ntot'][140,:,:],clim=[0,8]); plt.colorbar()
	
		# Isolate region
		indS=np.where( (lat0>42) & (lat0<62) & (lon0>-140) & (lon0<-100) )
		lat1=lat0[indS]
		lon1=lon0[indS]
		iLand=np.where(zRef['Data']==1)
	
		indT=np.where(tv<tv0[0])[0]
		z0=d0['ntot'][0,:,:].T # transposed to match lat and lon
		z1=z0[indS]
		z2=griddata( (lon1,lat1),z1,(lon,lat),method='cubic')
		for j in indT:
			d['ndep']['Data'][j,:,:]=z2
			d['ndep']['TS Mean'][j]=np.mean(z2[iLand])
	
		for iT in range(tv0.size):
			indT=np.where(tv==tv0[iT])[0]
			z0=d0['ntot'][iT,:,:].T # transposed to match lat and lon
			z1=z0[indS]
			z2=griddata( (lon1,lat1),z1,(lon,lat),method='cubic')
			d['ndep']['Data'][indT,:,:]=z2
			#plt.matshow(z2,clim=[0,8]); plt.colorbar()
			d['ndep']['TS Mean'][indT]=np.mean(z2[iLand])
			#print(np.mean(z2[ind]))

		# Future
		tv0=np.arange(2006,2100,1)
		d0=gu.ReadNC(meta['Paths']['DB']['NDEP'] + '\\ISIMIP\\ndep_nhx_rcp' + rcp + 'soc_monthly_2006_2099.nc4')
		d0['nhx']=gu.BlockSum(d0['nhx'],12)
		d1=gu.ReadNC(meta['Paths']['DB']['NDEP'] + '\\ISIMIP\\ndep_noy_rcp' + rcp + 'soc_monthly_2006_2099.nc4')
		d1['noy']=gu.BlockSum(d1['noy'],12)
		d0['ntot']=10*(d0['nhx']+d1['noy'])
		for iT in range(tv0.size):
			indT=np.where(tv==tv0[iT])[0]
			z0=d0['ntot'][iT,:,:].T # transposed to match lat and lon
			z1=z0[indS]
			z2=griddata( (lon1,lat1),z1,(lon,lat),method='cubic')
			d['ndep']['Data'][indT,:,:]=z2
			#plt.matshow(z2,clim=[0,8]); plt.colorbar()
			d['ndep']['TS Mean'][indT]=np.mean(z2[iLand])
	
		indT=np.where(tv>tv0[-1])[0]
		z0=d0['ntot'][-1,:,:].T # transposed to match lat and lon
		z1=z0[indS]
		z2=griddata( (lon1,lat1),z1,(lon,lat),method='cubic')
		for j in indT:
			d['ndep']['Data'][j,:,:]=z2
			d['ndep']['TS Mean'][j]=np.mean(z2[iLand])
		plt.plot(tv,np.mean(np.mean(d['ndep']['Data'],axis=2),axis=1),'k-')
		#plt.close('all'); plt.plot(tv,d['ndep']['TS Mean'],'g-.')
		#plt.matshow(np.mean(d['ntot'],axis=0))
		d['ndep']['Data']=d['ndep']['Data']/d['ndep']['Scale Factor']
		d['ndep']['Data']=d['ndep']['Data'].astype('int16')
		gu.opickle(meta['Paths']['DB']['NDEP'] + '\\ISIMIP\\ndep_' + rcp + '.pkl',d)
	return

#%%
def Process_Ndep_FromAckerman2018(meta):
	# Ackerman et al. (2018) (https://conservancy.umn.edu/handle/11299/197613)
	data=pd.read_csv(r'C:\Users\rhember\Documents\Data\Nitrogen Deposition\Ackerman\inorganic_N_deposition.csv')
	
	zTSA=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Admin\tsa.tif')

	#nam='tot_1985'
	nam='tot_2016'
	
	lat=np.unique(data['latitude'].to_numpy())
	lon=np.unique(data['longitude'].to_numpy())
	z0=data[nam].to_numpy()/100
	m=lat.size
	n=lon.size
	z0=np.flip(z0.reshape(n,m).T,0)
	
	z=zTSA.copy()
	z['Data']=z0
	z['Data']=z['Data'].astype('float32')
	
	srs=gis.ImportSRSs()
	#z['Projection']=srs['Proj']['Geographic']
	z['Proj4_String']=srs['String']['Geographic']
	
	sr=osr.SpatialReference()
	sr.ImportFromEPSG(4326)
	z['Projection']=sr.ExportToWkt()
	
	z['Y']=np.tile(np.reshape(lat,(-1,1)),(1,n))
	z['X']=np.tile(np.reshape(lon,(1,-1)),(m,1))
	z['m']=m
	z['n']=n
	
	# extent labels (left, right, bottom, top)
	arra=np.array([np.min(lon),np.max(lat),np.min(lat),np.max(lat)])
	z['Extent']=tuple(arra)
	
	# Transform
	arra=np.array([np.min(lon),2.5,0.0,np.max(lat),0.0,-2.0])
	z['gt']=tuple(arra)
	
	# Transform
	z['Transform']=from_origin(np.min(lon),np.max(lat),2.5,2.0)
	
	fout=r'C:\Users\rhember\Documents\Data\Nitrogen Deposition\Ackerman\ndep_' + nam + '.tif'
	gis.SaveGeoTiff(z,fout)
	
	# Clip to NA
	fout=r'C:\Users\rhember\Documents\Data\Nitrogen Deposition\Ackerman\ndep_' + nam + '.tif'
	z=gis.OpenGeoTiff(fout)
	
	#plt.close('all')
	#fig,ax=plt.subplots(1,figsize=gu.cm2inch(16,16))
	#im=ax.matshow(z['Data'],clim=(0,12),extent=z['Extent'])
	
	zC=gis.ClipRaster(z,[-140,-40],[0,90])
	
	fout=r'C:\Users\rhember\Documents\Data\Nitrogen Deposition\Ackerman\ndep_' + nam + '_clip.tif'
	gis.SaveGeoTiff(zC,fout)
	
	plt.close('all')
	fig,ax=plt.subplots(1,figsize=gu.cm2inch(16,16))
	im=ax.matshow(zC['Data'],clim=(0,12),extent=zC['Extent'])

	# Reproject
	pthin=r'C:\Users\rhember\Documents\Data\Nitrogen Deposition\Ackerman\ndep_' + nam + '_clip.tif'
	pthout=r'C:\Users\rhember\Documents\Data\Nitrogen Deposition\Ackerman\ndep_' + nam + '_clip_projected.tif'
	gis.ReprojectGeoTiff(pthin,pthout,srs['String']['BC1ha'])
	#gis.ReprojectGeoTiff(pthin,pthout,srs['String']['NACID'])
	
	z=gis.OpenGeoTiff(pthout)
	z['X']=z['X'][:,0:-1]
	
	plt.close('all')
	fig,ax=plt.subplots(1,figsize=gu.cm2inch(16,16))
	im=ax.matshow(z['Data'],clim=(0,12),extent=z['Extent'])
	
	#
	zND=zTSA.copy()
	ivl=10
	zND['X']=zTSA['X'][0::ivl,0::ivl]
	zND['Y']=zTSA['Y'][0::ivl,0::ivl]
	zND['Data']=griddata( (z['X'].flatten(),z['Y'].flatten()) ,z['Data'].flatten(),(zND['X'],zND['Y']),method='cubic')
	
	plt.close('all')
	fig,ax=plt.subplots(1,figsize=gu.cm2inch(16,16))
	im=ax.matshow(zND['Data'],clim=(0,5),extent=z['Extent'])
	
	# Plot
	# See bc1ha_map_full for plots
	return

#%% Potential evpapotranspiration
'''============================================================================
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
def GetETp(vi,Dims,Method,TimeInterval):

	if Dims=='Months by sites':
		# Dimensions (months by sites)
		N_m,N_s=vi['rswn'].shape
	elif Dims=='Sparse grid':
		N_m=1
		N_s=vi['rswn'].size
	elif Dims=='Single site':
		N_m=vi['rswn'].size
		N_s=1

	# Wind speed - if not supplied, use default of 2.0 m s-1.
	if 'u' not in vi:
		vi['u']=2.0

	# Constants
	con=HydroMetCon()

	# Convert net shortwave radiation from MJ m-2 d-1 to W m-2
	Rswn_conv=vi['rswn']*1e6/con['DayLength']

	# Convert net shortwave radiation (W m-2) to net radiation (W m-2)
	# Parameters from this were fitted to data at DF49 and agree roughly with
	# Landsberg's textbook.
	Rn=con['Rswn2Rn_Slope']*Rswn_conv+con['Rswn2Rn_AddOffset']

	# Psychrometric term (hPa K-1)
	Psychro=0.01*GetPsychrometric(vi['tmean'],'Pressure')

	# Saturation vapour pressure/temperature slope (hPa K-1)
	Svps=0.01*GetSVPSlope(vi['tmean'],'Pressure')

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

		# Extremely high combinations of temperature, radiation and vpd cuase extreme anomalies, truncating
		ETp=np.minimum(14,ETp)

	else:
		print('Method not recognized, quiting.')
		return ()

	# Unit conversion
	if Dims=='Months by sites':
		if (TimeInterval=='Month') | (TimeInterval=='M') | (TimeInterval=='m'):
			DIM=mb.repmat(np.reshape(con['DIM'],(-1,1)),N_yr,N_s)
			ETp=ETp*DIM
	elif Dims=='Sparse grid':
		ETp=ETp*con['DIM'][vi['Month']-1]

	ETp=np.maximum(0,ETp)

	return ETp

#%% WBM - Monthly Surface Water Balance Model
'''============================================================================

Inputs variables (field names in vi structure):

	LAI, leaf area index (m2 m-2)
	lat, latitude (deg)
	ta, mean air temperature (degC)
	prcp, total precipitation (mm month-1)
	rswn, net shortwave radiation (MJ m-2 d-1)
	vpd, vapour pressure deficit (hPa)
	ws, water content of soil (mm) ** Optional **
	wsp, water content of snow (mm) ** Optional **

Outputs variables (field names in ov structure):

	etp, potential evapotranspiration (mm month-1)
	eta, actual evapotranspiration (mm month-1)
	ws, water content of soil (mm)
	wsp, water content of snow (mm)
	melt, snowmelt(mm month-1)
	runoff, runoff (mm month-1)
	rainf, rain to snow ratio (dim)

Mode of operation:

If you are working with a small number of locations, use
"Combined" which runs all spatial units simultaneously. If running with large
grids, use "Grid" which inputs and outputs one month at a time. In the latter
approach, the state variables need to be initialized for the first time step,
and included as input arguments (from the output of t-1) for subsequent time
steps.

Parameters required in meta structure:

	Ws_max: water storage capacity (200 mm)
	Tmin: minimum rain/snow fraction temperature (degC)
	Tmax: maximum rain/snow fraction temperature (degC)
	Daily Interval: daily soil water computation interval (5)

Warnings:

	1) This code will not tolerate NaNs.
	2) The model requires spin-up for at least one year.

============================================================================'''
def WBM_SparseGrid(par,vi):
	N_s=vi['tmean'].size
	vo={}

	if 'ws' in vi:
		# Use final soil water content and snowpack water content from end of 
		# previous month
		vo['ws']=vi['ws']
		vo['wsp']=vi['wsp']
	else:
		# Initialize for first time step
		vo['ws']=par['Ws_max']*np.ones(N_s)
		vo['wsp']=np.zeros(N_s)
   
	# Constants
	con=HydroMetCon()
	
	# Potential evapotranspiratiom (mm month-1)
	vo['etp']=GetETp(vi,'Sparse grid',par['ETp Method'],'Month')

	vo['eta']=np.zeros(N_s)
	vo['runoff']=np.zeros(N_s)
	vo['melt']=np.zeros(N_s)

	# Canopy interception as a function of leaf area (Landsberg et al. 2005)
	# Maximum proportion of rainfall evaporated from canopy MaxIntcptn – 0.15
	# Leaf area index for maximum rainfall interception LAImaxIntcptn – 5
	
	# Fraction of intercepted precipitation
	FracPrecipInt=par['Ei_FracMax']*np.minimum(1.0,vi['LAI']/par['Ei_ALMax'])

	# Potential evaporation of intercepted precipiration (mm month-1)
	Ei_Potential=FracPrecipInt*vi['prcp']

	# Actual evaporation of intercepted precipiration, defined as the minimum
	# between the energy-limited rate, and the water-limited rate (mm month-1)
	Ei_Actual=np.minimum(vo['etp'],Ei_Potential)

	# Transpiration at energy-limited rate, i.e. prior to adding surface
	# constraints (mm month-1)
	Et_EnergyLimited=vo['etp']-Ei_Actual

	# Throughfall (mm month-1)
	P_Throughfall=vi['prcp']-Ei_Actual

	# Partition total precipitation into solid and liquid components
	fT=(vi['tmean']-par['Tmin'])/(par['Tmax']-par['Tmin'])
	fT=np.minimum(np.maximum(0,fT),1)
	Pr=fT*P_Throughfall
	Ps=P_Throughfall-Pr

	# Set inititial daily water pools at levels from the end of the last
	# month
	ws_d=vo['ws']
	wsp_d=vo['wsp']

	# Daily fluxes (mm d-1)
	Ps_d=Ps/(30/par['Daily_Interval'])
	Pr_d=Pr/(30/par['Daily_Interval'])
	Ei_Actual_d=Ei_Actual/(30/par['Daily_Interval'])
	Et_EnergyLimited_d=Et_EnergyLimited/(30/par['Daily_Interval'])

	for iDay in range(0,30,par['Daily_Interval']):

		# Potential snowmelt (mm d-1), equation from Thornthwaite and Mather (1955)
		M_d=2.63+2.55*vi['tmean']+0.0912*vi['tmean']*Pr_d

		# Actual snowmelt (mm d-1)
		M_d=np.maximum(0,np.minimum(M_d,wsp_d+Ps_d))

		# Cumulative snowmelt (mm)
		vo['melt']=vo['melt']+M_d

		# Updated snowpack water content (mm)
		wsp_d=wsp_d+Ps_d-M_d

		# Update soil water content (mm)
		ws_d=ws_d+M_d+Pr_d

		# Supply function (Willmott et al. 1985)
		# x=[0:0.1:1]' plot(x,1-exp(-6.68*(x)),'ko')
		fws=np.minimum(1,np.maximum(0,1-np.exp(-6.68*(ws_d/par['Ws_max']))))

		# Actual evaporation (mm d-1)
		Et_Actual_d=fws*Et_EnergyLimited_d

		# Cumulative actual evapotranspiration as the sum of wet-canopy
		# evaporation and transpiration (mm)
		vo['eta']=vo['eta']+Ei_Actual_d+Et_Actual_d

		# Remove transpiration from soil water pool (mm)
		ws_d=ws_d-Et_Actual_d

		# Find any spatial units where soil water exceeds capacity and add
		# "surplus" to monthly runoff and restrict soil water content
		# to capacity
		runoff_d=np.maximum(0,ws_d-par['Ws_max'])
		ws_d=np.minimum(ws_d,par['Ws_max'])

		# Update monthly runoff (mm)
		vo['runoff']=vo['runoff']+runoff_d

	# Update water pools (mm)
	vo['ws']=ws_d
	vo['wsp']=wsp_d

	# Constrain pools to be positive
	vo['ws']=np.maximum(0,vo['ws'])
	vo['wsp']=np.maximum(0,vo['wsp'])

	# Include rainfall fraction
	if par['Include Rainfall Fraction']=='Yes':
		vo['rainf']=np.minimum(1,np.maximum(0,Pr/Ps))

	return vo

# def WBM(par,vi):
#	 #--------------------------------------------------------------------------
#	 # Set up output arguments
#	 #--------------------------------------------------------------------------

#	 vo={}

#	 # If it is just one site, the user may import variables with 1d - change to 2d
#	 if vi['tmean'].ndim==1:
#		 for nam in vi.keys():
#			 vi[nam]=np.reshape(vi[nam],(-1,1))

#	 # Size of input variables
#	 N_mo,N_s=vi['tmean'].shape

#	 if par['Method']=='Combined':

#		 # Running all time steps for a set of n locations
#		 vo['ws']=par['Ws_max']*np.ones((N_mo,N_s))
#		 vo['Wsp']=np.zeros((N_mo,N_s))
#		 vo['ETp']=np.zeros((N_mo,N_s))
#		 vo['ETa']=np.zeros((N_mo,N_s))
#		 vo['R']=np.zeros((N_mo,N_s))
#		 vo['M']=np.zeros((N_mo,N_s))

#	 elif par['Method']=='Grid':

#		 # If running one month grid at a time, set initial values beased on input
#		 # from last run

#		 # Initialize state variables with maximum levels and add previous time
#		 # step for time steps beyond the first one

#		 if 'Ws' in vi:
#			 vo['Ws'][0,:]=vi['Ws']
#			 vo['Wsp'][0,:]=vi['Wsp']
#		 else:
#			 vo['Ws']=par['Ws_max']*np.ones((1,N_mo))
#			 vo['Wsp']=np.zeros((1,N_mo))

#		 vo['ETp']=np.zeros((1,N_mo))
#		 vo['ETa']=np.zeros((1,N_mo))
#		 vo['R']=np.zeros((1,N_mo))
#		 vo['M']=np.zeros((1,N_mo))

#	 #--------------------------------------------------------------------------
#	 # Constants
#	 #--------------------------------------------------------------------------

#	 con=HydroMetCon()

#	 #--------------------------------------------------------------------------
#	 # Potential evapotranspiratiom (mm month-1)
#	 #--------------------------------------------------------------------------

#	 vo['ETp']=GetETp(vi,par['ETp Method'],'Month')

#	 #--------------------------------------------------------------------------
#	 # Canopy interception as a function of leaf area (Landsberg et al. 2005)
#	 # Maximum proportion of rainfall evaporated from canopy MaxIntcptn – 0.15
#	 # Leaf area index for maximum rainfall interception LAImaxIntcptn – 5
#	 #--------------------------------------------------------------------------

#	 # Fraction of intercepted precipitation
#	 FracPrecipInt=par['Ei_FracMax']*np.minimum(1.0,vi['LAI']/par['Ei_ALMax'])

#	 # Potential evaporation of intercepted precipiration (mm month-1)
#	 Ei_Potential=FracPrecipInt*vi['prcp']

#	 # Actual evaporation of intercepted precipiration, defined as the minimum
#	 # between the energy-limited rate, and the water-limited rate (mm month-1)
#	 Ei_Actual=np.minimum(vo['ETp'],Ei_Potential)

#	 # Transpiration at energy-limited rate, i.e. prior to adding surface
#	 # constraints (mm month-1)
#	 Et_EnergyLimited=vo['ETp']-Ei_Actual

#	 # Throughfall (mm month-1)
#	 P_Throughfall=vi['prcp']-Ei_Actual

#	 # Partition total precipitation into solid and liquid components
#	 fT=(vi['tmean']-par['Tmin'])/(par['Tmax']-par['Tmin'])
#	 fT=np.minimum(np.maximum(0,fT),1)
#	 Pr=fT*P_Throughfall
#	 Ps=P_Throughfall-Pr

#	 #--------------------------------------------------------------------------
#	 # Loop through months
#	 #--------------------------------------------------------------------------

#	 for iT in range(N_mo):

#		 #----------------------------------------------------------------------
#		 # Set inititial daily water pools at level sfromt he end of the last
#		 # month
#		 #----------------------------------------------------------------------

#		 if par['Method']=='Combined':

#			 if iT==0:
#				 Ws_d=vo['Ws'][iT,:]
#				 Wsp_d=vo['Wsp'][iT,:]
#			 else:
#				 Ws_d=vo['Ws'][iT-1,:]
#				 Wsp_d=vo['Wsp'][iT-1,:]

#		 elif par['Method']=='Grid':

#			 Ws_d=vo['Ws']
#			 Wsp_d=vo['Wsp']

#		 #----------------------------------------------------------------------
#		 # Daily fluxes (mm d-1)
#		 #----------------------------------------------------------------------

#		 Ps_d=Ps[iT,:]/(30/par['Daily_Interval'])
#		 Pr_d=Pr[iT,:]/(30/par['Daily_Interval'])
#		 Ei_Actual_d=Ei_Actual[iT,:]/(30/par['Daily_Interval'])
#		 Et_EnergyLimited_d=Et_EnergyLimited[iT,:]/(30/par['Daily_Interval'])

#		 for iDay in range(0,30,par['Daily_Interval']):

#			 # Potential snowmelt (mm d-1), equation from Thornthwaite and Mather (1955)
#			 M_d=2.63+2.55*vi['tmean'][iT,:]+0.0912*vi['tmean'][iT,:]*Pr_d

#			 # Actual snowmelt (mm d-1)
#			 M_d=np.maximum(np.zeros((1,N_s)),np.minimum(M_d,Wsp_d+Ps_d))

#			 # Cumulative snowmelt (mm)
#			 vo['M'][iT,:]=vo['M'][iT,:]+M_d

#			 # Updated snowpack water content (mm)
#			 Wsp_d=Wsp_d+Ps_d-M_d

#			 # Update soil water content (mm)
#			 Ws_d=Ws_d+M_d+Pr_d

#			 # Supply function (Willmott et al. 1985)
#			 # x=[0:0.1:1]' plot(x,1-exp(-6.68*(x)),'ko')
#			 fWs=np.minimum(1,np.maximum(0,1-np.exp(-6.68*(Ws_d/par['Ws_max']))))

#			 # Actual evaporation (mm d-1)
#			 Et_Actual_d=fWs*Et_EnergyLimited_d

#			 # Cumulative actual evapotranspiration as the sum of wet-canopy
#			 # evaporation and transpiration (mm)
#			 vo['ETa'][iT,:]=vo['ETa'][iT,:]+Ei_Actual_d+Et_Actual_d

#			 # Remove transpiration from soil water pool (mm)
#			 Ws_d=Ws_d-Et_Actual_d

#			 # Find any spatial units where soil water exceeds capacity and add
#			 # "surplus" to monthly runoff and restrict soil water content
#			 # to capacity
#			 R_d=np.maximum(0,Ws_d-par['Ws_max'])
#			 Ws_d=np.minimum(Ws_d,par['Ws_max'])

#			 # Update monthly runoff (mm)
#			 vo['R'][iT,:]=vo['R'][iT,:]+R_d

#		 # Update water pools
#		 vo['Ws'][iT,:]=Ws_d
#		 vo['Wsp'][iT,:]=Wsp_d

#	 #--------------------------------------------------------------------------
#	 # Constrain pools to be positive
#	 #--------------------------------------------------------------------------

#	 vo['Ws']=np.maximum(0,vo['Ws'])
#	 vo['Wsp']=np.maximum(0,vo['Wsp'])

#	 #--------------------------------------------------------------------------
#	 # Include rainfall fraction
#	 #--------------------------------------------------------------------------

#	 if par['Include Rainfall Fraction']=='Yes':
#		 vo['RF']=np.minimum(1,np.maximum(0,Pr/Ps))
#	 return vo

#%% Hydrometerology constants
# con=gaia.HydroMetCon()
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

#%% Calculate psychrometric constant
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

#%% Calculate saturation vapour pressure (hPa) based on temperature (deg C)
def GetEstar(ta):
	# Buck's equation gives kPa, convert to hPa
	es=10*(0.61121*np.exp( (18.678-(ta/234.5)) * (ta/(257.14+ta)) ))
	#es=596.9277+44.1324*ta**2+1.6207*ta**3+0.0281*ta**4
	return es

#%% Calculate saturatoin vapour pressure vs. temperature slope
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

#%%
def GetDensityAir(x):
	# Calculate density of air (kg m-3) based on temperature (deg C)
	y=1.2756-0.0046478*x+2.2829e-005*x**2
	
	# Derivation
	#temp=[-100 -50 0 20 40 60]'; %C
	#rho=[1.98 1.534 1.293 1.205 1.127 1.067]'; %kg m-3
	#b=polyfit(temp,rho,2);
	#[r,p]=corrcoef(temp,rho);
	return y

#%%
def GYModelModifier(meta,pNam,NetGrowth,vo,iT):

	# Demonstrate
	flg=0
	if flg==1:
		t=np.arange(1850,2050,1)
		m0=0.5
		m1_young=1.5
		m1_old=1.0
		Age_young=10
		Age_old=100
		t0=1900
		t1=2000

		plt.close('all')
		A=0
		m1=(m1_young+m1_old)-gu.Clamp(m0+(Age_young-A)/(Age_young-Age_old),m1_old,m1_young)
		fG=gu.Clamp(m0+(m1-m0)/(t1-t0)*(t-t0),m0,m1)
		plt.plot(t,fG,'g-')
		A=150
		m1=(m1_young+m1_old)-gu.Clamp(m0+(Age_young-A)/(Age_young-Age_old),m1_old,m1_young)
		fG=gu.Clamp(m0+(m1-m0)/(t1-t0)*(t-t0),m0,m1)
		plt.plot(t,fG,'b--')

	# Import parameters
	if 'Opt' in meta[pNam].keys():
		#t0=meta[pNam]['Opt']['GrowthModifier']['t0']
		#t1=meta[pNam]['Opt']['GrowthModifier']['t1']
		t0=meta['Param']['BEV']['ByBGCZ']['GM t0']
		t1=meta['Param']['BEV']['ByBGCZ']['GM t1']
		m0=meta[pNam]['Opt']['GrowthModifier']['m0']
		m1_young=meta[pNam]['Opt']['GrowthModifier']['m1_y']
		m1_old=meta[pNam]['Opt']['GrowthModifier']['m1_o']
		Age_young=meta[pNam]['Opt']['GrowthModifier']['Age_y']
		Age_old=meta[pNam]['Opt']['GrowthModifier']['Age_o']
	else:
		t0=meta['Param']['BEV']['ByBGCZ']['GM t0']
		t1=meta['Param']['BEV']['ByBGCZ']['GM t1']
		m0=meta['Param']['BEV']['ByBGCZ']['GM m0']
		m1_young=meta['Param']['BEV']['ByBGCZ']['GM m1_y']
		m1_old=meta['Param']['BEV']['ByBGCZ']['GM m1_o']
		Age_young=meta['Param']['BEV']['ByBGCZ']['GM Age_y']
		Age_old=meta['Param']['BEV']['ByBGCZ']['GM Age_o']

	m1=(m1_young+m1_old)-gu.Clamp(m0+(Age_young-vo['A'][iT,:])/(Age_young-Age_old),m1_old,m1_young)
	fG=gu.Clamp(m0+(m1-m0)/(t1-t0)*(meta[pNam]['Year'][iT]-t0),m0,m1)
	fG=fG*np.array(vo['A'][iT,:]>1).astype('float')
	NetGrowth=np.tile(fG,(NetGrowth.shape[1],1)).T*NetGrowth

	# Plot the modifiers in the "by BGC zone" parameters spreadsheet
	flg=0
	if flg==1:
		d=gu.ReadExcel(r'G:\My Drive\Code_Python\fcgadgets\cbrunner\Parameters\Parameters_ByBGCZ.xlsx')
		t=np.arange(1850,2050,1)
		for iz in range(d['Name'].size):
			#np.where(d['Name']=='CWH')[0][0]
			plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(10,6));
			A=5
			m1=(d['GM m1_y'][iz]+d['GM m1_o'][iz])-gu.Clamp(d['GM m0'][iz]+(d['GM Age_y'][iz]-A)/(d['GM Age_y'][iz]-d['GM Age_o'][iz]),d['GM m1_o'][iz],d['GM m1_y'][iz])
			fG=gu.Clamp(d['GM m0'][iz]+(m1-d['GM m0'][iz])/(d['GM t1'][iz]-d['GM t0'][iz])*(t-d['GM t0'][iz]),d['GM m0'][iz],m1)
			plt.plot(t,fG,'g-',lw=1.25,color=[0.6,0.9,0],label='5 years old')
			A=175
			m1=(d['GM m1_y'][iz]+d['GM m1_o'][iz])-gu.Clamp(d['GM m0'][iz]+(d['GM Age_y'][iz]-A)/(d['GM Age_y'][iz]-d['GM Age_o'][iz]),d['GM m1_o'][iz],d['GM m1_y'][iz])
			fG=gu.Clamp(d['GM m0'][iz]+(m1-d['GM m0'][iz])/(d['GM t1'][iz]-d['GM t0'][iz])*(t-d['GM t0'][iz]),d['GM m0'][iz],m1)
			plt.plot(t,fG,'b--',lw=1.25,label='175 years old')
			ax.set(xlabel='Time, years',ylabel='Multiplier',ylim=[0,1.2])
			ax.legend(loc='upper left',frameon=False,facecolor=None,edgecolor='w')
			plt.tight_layout()
			gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\fcgadgets\GY Modifiers\Modifier_' + d['Name'][iz],'png',900)

	return NetGrowth

#%% Import CMIP6 data
def Process_CMIP6():
	#scn='ssp245'
	scn='ssp585'
	meta=u1ha.Init()
	zRef=gis.OpenGeoTiff(meta['Paths']['bc5k Ref Grid'])
	iLand=np.where(zRef['Data']==1)
	pth=r'C:\Data\Climate\CMIP6\Monthly'
	srs=gis.ImportSRSs()
	lon,lat=pyproj.transform(srs['Proj']['BC1ha'],srs['Proj']['Geographic'],zRef['X'],zRef['Y'])
	tvm=gu.tvec('m',1950,2101)
	tvH0=gu.tvec('m',1950,2014)
	tvF0=gu.tvec('m',2015,2100)
	press=1013.25
	con=HydroMetCon()
	vL=['tmean','ea','prcp','vpd','rswd']
	mdL=['CESM2','CNRM-CM6-1','GFDL-ESM4','HadGEM3-GC31-LL','MIROC-ES2L','MPI-ESM1-2-LR']
	#mdL=['MIROC-ES2L','MPI-ESM1-2-LR']
	dC={}
	for md in mdL:
		print(md)
		flL=os.listdir(pth + '\\' + md)
		def GetFileName(flL,v,period):
			fl1=[]
			for fl in flL:
				if fl.startswith(v):
					if period in fl:
						fl1=fl
			return fl1

		# Initialize
		dC[md]={}
		dC[md]={}
		for v in vL:
			dC[md][v]=np.zeros((tvm.shape[0],zRef['m'],zRef['n']))

		# Temperature
		fl=GetFileName(flL,'tas','historical')
		d0=gu.ReadNC(pth + '\\' + md + '\\' + fl)
		lat0,lon0=np.meshgrid(d0['lat'],d0['lon'])
		ind=np.where(lon0>180)[0]; lon0[ind]=lon0[ind]-360
		indS=np.where( (lat0>40) & (lat0<65) & (lon0>-140) & (lon0<-100) )
		lat1=lat0[indS]; lon1=lon0[indS]
		z0=d0['tas']
		for iT in range(z0.shape[0]):
			z1=z0[iT,:,:].T
			z1=z1[indS]
			z1=griddata( (lon1,lat1),z1,(lon,lat),method='cubic')
			indT=np.where( (tvm[:,0]==tvH0[iT,0]) & (tvm[:,1]==tvH0[iT,1]) )[0]
			dC[md]['tmean'][indT[0],:,:]=z1-273.15
		fl=GetFileName(flL,'tas',scn)
		d0=gu.ReadNC(pth + '\\' + md + '\\' + fl)
		lat0,lon0=np.meshgrid(d0['lat'],d0['lon'])
		ind=np.where(lon0>180)[0]; lon0[ind]=lon0[ind]-360
		indS=np.where( (lat0>40) & (lat0<65) & (lon0>-140) & (lon0<-100) )
		lat1=lat0[indS]; lon1=lon0[indS]
		z0=d0['tas']
		for iT in range(z0.shape[0]):
			z1=z0[iT,:,:].T
			z1=z1[indS]
			z1=griddata( (lon1,lat1),z1,(lon,lat),method='cubic')
			indT=np.where( (tvm[:,0]==tvF0[iT,0]) & (tvm[:,1]==tvF0[iT,1]) )[0]
			dC[md]['tmean'][indT[0],:,:]=z1-273.15
		#plt.plot(dC[md]['tmean'][:,10,10])

		# Humidity
		fl=GetFileName(flL,'huss','historical')
		d0=gu.ReadNC(pth + '\\' + md + '\\' + fl)
		lat0,lon0=np.meshgrid(d0['lat'],d0['lon'])
		ind=np.where(lon0>180)[0]; lon0[ind]=lon0[ind]-360
		indS=np.where( (lat0>40) & (lat0<65) & (lon0>-140) & (lon0<-100) )
		lat1=lat0[indS]; lon1=lon0[indS]
		z0=d0['huss']*press/(0.378*d0['huss']+0.622) # Convert to hPa
		for iT in range(z0.shape[0]):
			z1=z0[iT,:,:].T
			z1=z1[indS]
			z1=griddata( (lon1,lat1),z1,(lon,lat),method='cubic')
			indT=np.where( (tvm[:,0]==tvH0[iT,0]) & (tvm[:,1]==tvH0[iT,1]) )[0]
			dC[md]['ea'][indT[0],:,:]=z1
		fl=GetFileName(flL,'huss',scn)
		d0=gu.ReadNC(pth + '\\' + md + '\\' + fl)
		lat0,lon0=np.meshgrid(d0['lat'],d0['lon'])
		ind=np.where(lon0>180)[0]; lon0[ind]=lon0[ind]-360
		indS=np.where( (lat0>40) & (lat0<65) & (lon0>-140) & (lon0<-100) )
		lat1=lat0[indS]; lon1=lon0[indS]
		z0=d0['huss']*press/(0.378*d0['huss']+0.622) # Convert to hPa
		for iT in range(z0.shape[0]):
			z1=z0[iT,:,:].T
			z1=z1[indS]
			z1=griddata( (lon1,lat1),z1,(lon,lat),method='cubic')
			indT=np.where( (tvm[:,0]==tvF0[iT,0]) & (tvm[:,1]==tvF0[iT,1]) )[0]
			dC[md]['ea'][indT[0],:,:]=z1
		#plt.plot(dC[md]['ea'][:,10,10])

		# Radiation
		fl=GetFileName(flL,'rsds','historical')
		d0=gu.ReadNC(pth + '\\' + md + '\\' + fl)
		lat0,lon0=np.meshgrid(d0['lat'],d0['lon'])
		ind=np.where(lon0>180)[0]; lon0[ind]=lon0[ind]-360
		indS=np.where( (lat0>40) & (lat0<65) & (lon0>-140) & (lon0<-100) )
		lat1=lat0[indS]; lon1=lon0[indS]
		z0=d0['rsds']*86400/1e6 # Converted from W/m2 to MJ/m2/d-1
		for iT in range(z0.shape[0]):
			z1=z0[iT,:,:].T
			z1=z1[indS]
			z1=griddata( (lon1,lat1),z1,(lon,lat),method='cubic')
			indT=np.where( (tvm[:,0]==tvH0[iT,0]) & (tvm[:,1]==tvH0[iT,1]) )[0]
			dC[md]['rswd'][indT[0],:,:]=z1
		fl=GetFileName(flL,'rsds',scn)
		d0=gu.ReadNC(pth + '\\' + md + '\\' + fl)
		lat0,lon0=np.meshgrid(d0['lat'],d0['lon'])
		ind=np.where(lon0>180)[0]; lon0[ind]=lon0[ind]-360
		indS=np.where( (lat0>40) & (lat0<65) & (lon0>-140) & (lon0<-100) )
		lat1=lat0[indS]; lon1=lon0[indS]
		z0=d0['rsds']*86400/1e6 # Converted from W/m2 to MJ/m2/d-1
		for iT in range(z0.shape[0]):
			z1=z0[iT,:,:].T
			z1=z1[indS]
			z1=griddata( (lon1,lat1),z1,(lon,lat),method='cubic')
			indT=np.where( (tvm[:,0]==tvF0[iT,0]) & (tvm[:,1]==tvF0[iT,1]) )[0]
			dC[md]['rswd'][indT[0],:,:]=z1
		#plt.plot(dC[md]['rswd'][:,10,10])

		# Precipitation
		fl=GetFileName(flL,'pr','historical')
		d0=gu.ReadNC(pth + '\\' + md + '\\' + fl)
		lat0,lon0=np.meshgrid(d0['lat'],d0['lon'])
		ind=np.where(lon0>180)[0]; lon0[ind]=lon0[ind]-360
		indS=np.where( (lat0>40) & (lat0<65) & (lon0>-140) & (lon0<-100) )
		lat1=lat0[indS]; lon1=lon0[indS]
		# Convert from kg m-2 s-1 to mm/mo
		z0=d0['pr']
		for mo in range(12):
			z0[mo::12,:,:]=86400*z0[mo::12,:,:]*con['DIM'][mo]
		for iT in range(z0.shape[0]):
			z1=z0[iT,:,:].T
			z1=z1[indS]
			z1=griddata( (lon1,lat1),z1,(lon,lat),method='cubic')
			indT=np.where( (tvm[:,0]==tvH0[iT,0]) & (tvm[:,1]==tvH0[iT,1]) )[0]
			dC[md]['prcp'][indT[0],:,:]=z1
		fl=GetFileName(flL,'pr',scn)
		d0=gu.ReadNC(pth + '\\' + md + '\\' + fl)
		lat0,lon0=np.meshgrid(d0['lat'],d0['lon'])
		ind=np.where(lon0>180)[0]; lon0[ind]=lon0[ind]-360
		indS=np.where( (lat0>40) & (lat0<65) & (lon0>-140) & (lon0<-100) )
		lat1=lat0[indS]; lon1=lon0[indS]
		# Convert from kg m-2 s-1 to mm/mo
		z0=d0['pr']
		for mo in range(12):
			z0[mo::12,:,:]=86400*z0[mo::12,:,:]*con['DIM'][mo]
		for iT in range(z0.shape[0]):
			z1=z0[iT,:,:].T
			z1=z1[indS]
			z1=griddata( (lon1,lat1),z1,(lon,lat),method='cubic')
			indT=np.where( (tvm[:,0]==tvF0[iT,0]) & (tvm[:,1]==tvF0[iT,1]) )[0]
			dC[md]['prcp'][indT[0],:,:]=z1
		#plt.plot(dC[md]['rswd'][:,10,10])

	# Add vpd
	for md in mdL:
		dC[md]['vpd']=np.maximum(0,GetEstar(dC[md]['tmean'])-dC[md]['ea'])

	# Add water balance
	for md in mdL:
		vwbL=['cwd','eta','etp','melt','runoff','ws','wsp']
		for v in vwbL:
			dC[md][v]=np.zeros((tvm.shape[0],zRef['m'],zRef['n']))
		par={}
		par['Method']='Grid'
		par['Ws_max']=200.0 # mm
		par['Tmin']=-3.0
		par['Tmax']=3.0
		par['Ei_FracMax']=0.15
		par['Ei_ALMax']=5.0
		par['ETp Method']='Penman-Monteith'
		par['Daily_Interval']=5
		par['Include Rainfall Fraction']='No'
		vi={}
		for iT in range(tvm.shape[0]):
			print('Year:' + str(tvm[iT,0]) + ', Month:' + str(tvm[iT,1]) )
			mo=tvm[iT,1]
			vi['Month']=mo
			vi['LAI']=5.0
			vi['Gs']=0.010
			vi['Ga']=0.058
			vi['tmean']=dC[md]['tmean'][iT,:,:][iLand]
			vi['prcp']=dC[md]['prcp'][iT,:,:][iLand]
			vi['vpd']=dC[md]['vpd'][iT,:,:][iLand]
			vi['rswn']=dC[md]['rswd'][iT,:,:][iLand]
			vo=WBM_SparseGrid(par,vi)
			vi['ws']=vo['ws']
			vi['wsp']=vo['wsp']
			vo['cwd']=vo['etp']-vo['eta']
			for v in vwbL:
				dC[md][v][iT,:,:][iLand]=vo[v]
	#plt.plot(dC[md]['ws'][:,150,150])

	# Calculate seasonal indicators
	tva=np.arange(1950,2101,1)
	iTmu=np.where( (tva>=1971) & (tva<=2000) )[0]
	dCS={}
	for md in mdL:
		dCS[md]={}
		dCS[md]={}
		dCS[md]['tmean']={'ann':np.zeros((tva.shape[0],zRef['m'],zRef['n']))}
		dCS[md]['ws']={'mjjas':np.zeros((tva.shape[0],zRef['m'],zRef['n']))}
		dCS[md]['cwd']={'mjjas':np.zeros((tva.shape[0],zRef['m'],zRef['n']))}
		dCS[md]['cmi']={'mjjas':np.zeros((tva.shape[0],zRef['m'],zRef['n']))}
		dCS[md]['prcp']={'mjjas':np.zeros((tva.shape[0],zRef['m'],zRef['n']))}
		dCS[md]['etp']={'mjjas':np.zeros((tva.shape[0],zRef['m'],zRef['n']))}
		for iY in range(tva.size):
			iT=np.where( (tvm[:,0]==tva[iY]) )[0]
			dCS[md]['tmean']['ann'][iY]=np.mean(dC[md]['tmean'][iT,:,:],axis=0)
			iT=np.where( (tvm[:,0]==tva[iY]) & (tvm[:,1]>=5) & (tvm[:,1]<=9) )[0]
			dCS[md]['ws']['mjjas'][iY]=np.mean(dC[md]['ws'][iT,:,:],axis=0)
			dCS[md]['cwd']['mjjas'][iY]=np.mean(dC[md]['cwd'][iT,:,:],axis=0)
			dCS[md]['cmi']['mjjas'][iY]=np.mean(dC[md]['prcp'][iT,:,:]-dC[md]['etp'][iT,:,:],axis=0)
			dCS[md]['prcp']['mjjas'][iY]=np.mean(dC[md]['prcp'][iT,:,:],axis=0)
			dCS[md]['etp']['mjjas'][iY]=np.mean(dC[md]['etp'][iT,:,:],axis=0)

	# Convert to anomalies
	dCSa=copy.deepcopy(dCS)
	for md in mdL:
		for k in dCSa[md].keys():
			for j in dCSa[md][k].keys():
				mu=np.mean(dCSa[md][k][j][iTmu,:,:],axis=0)
				for iY in range(tva.size):
					dCSa[md][k][j][iY,:,:]=dCSa[md][k][j][iY,:,:]-mu

	# Compress
	for md in mdL:
		for k in dCS[md].keys():
			for j in dCS[md][k].keys():
				dCS[md][k][j]=dCS[md][k][j]/meta['SF']['Climate'][k]
				dCS[md][k][j]=dCS[md][k][j].astype('int16')
				dCSa[md][k][j]=dCSa[md][k][j]/meta['SF']['Climate'][k]
				dCSa[md][k][j]=dCSa[md][k][j].astype('int16')

	# Save
	gu.opickle(pth + '\\cmip6_actuals_' + scn + '.pkl',dCS)
	gu.opickle(pth + '\\cmip6_anomalies' + scn + '.pkl',dCSa)
	return

#%%



