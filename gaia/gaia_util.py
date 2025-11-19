import os
import numpy as np
import pandas as pd
import copy
import matplotlib.pyplot as plt
import numpy.matlib as mb
import pyproj
import rasterio
import time
import shutil
from fair import FAIR
from fair.interface import fill,initialise
from fair.io import read_properties
from scipy.optimize import curve_fit
import fcgadgets.macgyver.util_general as gu
import fcgadgets.macgyver.util_gis as gis
import fcgadgets.bc1ha.bc1ha_utils as u1ha
import fcgadgets.cbrunner.cbrun_postprocess as post

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

#%% Import environmental data
# *** Replaced by bc5k compiler ***
# def ImportEnvironment(meta,tv):
# 	d={'co2':{},'ndep':{}}

# 	# Import CO2 (ppm), source: https://www.pik-potsdam.de/~mmalte/rcps/ (ppm)
# 	dC=gu.ReadExcel(meta['Paths']['DB']['CO2'],sheet_name='Sheet1',skiprows=0)
# 	d['co2']['ssp245']=np.zeros(tv.size)
# 	d['co2']['ssp585']=np.zeros(tv.size)
# 	for iT in range(tv.size):
# 		ind=np.where(dC['Year']==tv[iT])[0]
# 		d['co2']['ssp245'][iT]=dC['CO2 RCP45'][ind]
# 		d['co2']['ssp585'][iT]=dC['CO2 RCP85'][ind]

# 	# Import N deposition
# 	zRef5=gis.OpenGeoTiff(meta['Paths']['bc5k Ref Grid'])
# 	iLand=np.where(zRef5['Data']==1)
# 	rcpL=['26','60','85']
# 	for rcp in rcpL:
# 		d['ndep'][rcp]=np.zeros(tv.size)
# 		dN=gu.ipickle(meta['Paths']['DB']['NDEP'] + '\\ISIMIP\\ndep_' + rcp + '.pkl')
# 		dN['ndep']['Data']=dN['ndep']['Data'].astype('float')*meta['Climate']['SF']['ndep']
# 		for iT in range(tv.size):
# 			indT=np.where(dN['ndep']['tv']==tv[iT])[0]
# 			if indT.size==0:
# 				continue
# 			d['ndep'][rcp][iT]=np.mean(dN['ndep']['Data'][indT[0],:,:][iLand])
# 		indT=np.where(tv<dN['ndep']['tv'][0])[0]
# 		d['ndep'][rcp][indT]=np.mean(dN['ndep']['Data'][0,:,:][iLand])
# 		indT=np.where(tv>dN['ndep']['tv'][-1])[0]
# 		d['ndep'][rcp][indT]=np.mean(dN['ndep']['Data'][-1,:,:][iLand])
# 	return d

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
		N_m,N_s=vi['tmean'].shape
	elif Dims=='Sparse grid':
		N_m=1
		N_s=vi['tmean'].size
	elif Dims=='Single site':
		N_m=vi['tmean'].size
		N_s=1

	# Wind speed - if not supplied, use default of 2.0 m s-1.
	if 'u' not in vi:
		vi['u']=2.0

	# Constants
	con=HydroMetCon()

	if 'rn' not in vi:
		# Convert net shortwave radiation from MJ m-2 d-1 to W m-2
		Rswn_conv=vi['rswn']*1e6/con['DayLength']
	
		# Convert net shortwave radiation (W m-2) to net radiation (W m-2)
		# Parameters from this were fitted to data at DF49 and agree roughly with
		# Landsberg's textbook.
		Rn=con['Rswn2Rn_Slope']*Rswn_conv+con['Rswn2Rn_AddOffset']
	else:
		Rn=vi['rn']

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

	# Stefan-Boltzmann constant (kg s-3 K-4)
	con['Stefan-Boltzmann Constant']=5.67037442*10**-8

	# Emissivity (dimensionless)
	con['Emissivity Conifer Forest']=0.98

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

	# Karman's constant
	con['Karmans Constant']=0.4

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

#%%
def Calc_TOA_Reflected_RSW(rswd,a_surf):
	# Bird et al. (Eq 5)
	f_cloud=0.7
	S=1-f_cloud # Sunniness (% of possible sunshine
	k_cloud=0.15 # reflectivity of clouds (0=clouds are transparent)
	a_cloud=k_cloud*(1-S)
	t_cloud=0.52 # Haigh (2007)
	rswr=rswd*(a_cloud+( (t_cloud**2*a_surf)/(1-a_cloud*a_surf) ))
	return rswr

#%%
def ConvertEmissionToRF(E):
	# Input: E, emission (tCO2e/ha)

	# Rotenberg and Yakir (2010)
	co2_re=5.35 # Radiative efficiency of CO2 (W m2)
	co2_ref=385 # Base CO2 concentration ppmv
	co2_specific_density=2.13E+12 #	Specific density of CO2 (kg C ppmv-1)
	af=0.5 # CO2 airborne fraction
	#M_co2=44.0065 #Molecular mass of carbon (g/mol)
	#M_air=28.95 # Molecular mass of dry air (g/mol)
	#m_air=5.148e15 # mass of atmosphere (Mg)
	nanoWperW=1e9
	sqmperha=1e4
	kgpertonne=1e3

	# Convert units of input emission
	E_c=E/3.667*kgpertonne/sqmperha # kgC m-2
	#E_c=E/3.667*kgpertonne # kgC ha-1

	# Convert emission to change in atmospheric concentration
	co2_delta=E_c/co2_specific_density # ppmv CO2

	# Radiative forcing (W m2)
	RF=co2_re*np.log(1+(co2_delta/co2_ref))

	# Radiative forcing (nW m-2 ha-1)
	RF=RF*nanoWperW*sqmperha
	#RF=RF*nanoWperW

	return RF

#%%
def Run_FAIR(meta,dfConfig,dfVolcanic,path_emissions):

	f=FAIR() # ch4_method='thornhill2021'
	f.define_time(1750,2150,1) # Define time period

	# Define forcings
	#frcL=['CO2','CH4','N2O','CO','Land use','Volcanic','Solar','Aerosol-cloud interactions','Aerosol-radiation interactions','Ozone','Stratospheric water vapour']
	#meta['Modules']['FAIR']['Forcings']

	# Define scenarios
	sL=['ssp370']
	f.define_scenarios(sL)

	# Define configs
	models=dfConfig['model'].unique()
	configs=[]
	for imodel, model in enumerate(models):
		for run in dfConfig.loc[dfConfig['model']==model,'run']:
			configs.append(f"{model}_{run}")
	f.define_configs(configs)
	
	# Define species
	#species,properties=read_properties(path + '//species_configs_properties.csv')
	species,properties=read_properties()
	f.define_species(species,properties)
	
	# Create arrays
	f.allocate()
	
	# Populate inputs
	f.fill_species_configs()
	
	# Populate emissions
	f.fill_from_rcmip()
	f.fill_from_csv(emissions_file=path_emissions)
	
	# Overwrite volcanic emissions
	dfVolcanic[1750:].head()
	#volcanic_forcing=np.zeros(351)
	volcanic_forcing=np.zeros(401)
	volcanic_forcing[:271]=dfVolcanic[1749:].groupby(np.ceil(dfVolcanic[1749:].index) // 1).mean().squeeze().values
	fill(f.forcing,volcanic_forcing[:,None,None],specie="Volcanic") # sometimes need to expand the array
	
	# The concentration in the first timestep will be set to the baseline concentration, which are the IPCC AR6 1750 values.
	initialise(f.concentration,f.species_configs['baseline_concentration'])
	initialise(f.forcing,0)
	initialise(f.temperature,0)
	initialise(f.cumulative_emissions,0)
	initialise(f.airborne_emissions,0)
	
	# Populate climate configurations
	models=dfConfig['model'].unique()
	seed=1355763
	for config in configs:
		model, run = config.split('_')
		condition = (dfConfig['model']==model) & (dfConfig['run']==run)
		fill(f.climate_configs['ocean_heat_capacity'], dfConfig.loc[condition, 'C1':'C3'].values.squeeze(), config=config)
		fill(f.climate_configs['ocean_heat_transfer'], dfConfig.loc[condition, 'kappa1':'kappa3'].values.squeeze(), config=config)
		fill(f.climate_configs['deep_ocean_efficacy'], dfConfig.loc[condition, 'epsilon'].values[0], config=config)
		fill(f.climate_configs['gamma_autocorrelation'], dfConfig.loc[condition, 'gamma'].values[0], config=config)
		fill(f.climate_configs['sigma_eta'], dfConfig.loc[condition, 'sigma_eta'].values[0], config=config)
		fill(f.climate_configs['sigma_xi'], dfConfig.loc[condition, 'sigma_xi'].values[0], config=config)
		fill(f.climate_configs['stochastic_run'], True, config=config)
		fill(f.climate_configs['use_seed'], True, config=config)
		fill(f.climate_configs['seed'], seed, config=config)
		seed = seed + 399
	
	# Run
	f.run()

	d={}
	d['Time']=f.timebounds
	for scn in sL:
		d[scn]={}
		d[scn]['Emissions']={}
		for sp in species:
			d[scn]['Emissions'][sp]=np.mean(f.emissions.loc[dict(scenario=scn,specie=sp)].to_numpy(),axis=1)
		d[scn]['Concentrations']={}
		for sp in species:
			d[scn]['Concentrations'][sp]=np.mean(f.concentration.loc[dict(scenario=scn,specie=sp)].to_numpy(),axis=1)
		d[scn]['RF']={}
		for frc in meta['Modules']['FAIR']['Forcings']:
			d[scn]['RF'][frc]=np.mean(f.forcing.loc[dict(scenario=scn,specie=frc)].to_numpy(),axis=1)
		d[scn]['RF Total']=np.mean(f.forcing_sum.loc[dict(scenario=scn)].to_numpy(),axis=1)
		d[scn]['Temperature']=np.mean(f.temperature.loc[dict(scenario=scn,layer=0)].to_numpy(),axis=1)

	return d

#%%
def CalcPulseEmissionBenchmark_RF_FAIR(meta,pNam,yrPulse):

	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)

	# Modelling in the time period for CMIP6
	d0={}

	# Import FAIR data
	dfConfig=pd.read_csv(meta['Paths']['DB']['FAIR'] + '//Variables_FAIR_4xCO2_cummins_ebm3.csv')
	dfVolcanic=pd.read_csv(meta['Paths']['DB']['FAIR'] + '//Variables_FAIR_volcanic_ERF_monthly_175001-201912.csv',index_col='year')

	# Run the reference (CMIP6-SSP370)
	path_emissions_ref=meta['Paths']['DB']['FAIR'] + '//Variables_FAIR_Emissions_SSP370.csv'
	d0['Reference']=Run_FAIR(meta,dfConfig,dfVolcanic,path_emissions_ref)

	 # 1 GtCO2/year, 1 MtCH4/year, 1 MtN2O/year, 1 MtCO/year
	E_PulseCO2=1
	E_PulseCH4=1
	E_PulseN2O=1
	E_PulseCO=1
	scnName=['CO2 Pulse','CH4 Pulse','N2O Pulse','CO Pulse']

	# Run each scenario
	for iScn,scn in enumerate(scnName):

		# Create a copy of SSP370 emissions
		path_emissions_adj=meta['Paths']['DB']['FAIR'] + '//Variables_FAIR_Emissions_SSP370_adjusted.csv'
		shutil.copyfile(path_emissions_ref,path_emissions_adj)
		E=gu.ReadCSV(meta['Paths']['DB']['FAIR'] + '//Variables_FAIR_Emissions_SSP370_adjusted.csv',skip_rows=0)

		# Add CO2
		if iScn==0:
			indV=np.where( (E['Variable']=='CO2 AFOLU') )[0]
			for i,yr in enumerate(tv):
				E[str(yrPulse+0.5)][indV]=E[str(yrPulse+0.5)][indV]+E_PulseCO2
		# Add CH4
		if iScn==1:
			indV=np.where( (E['Variable']=='CH4') )[0]
			for i,yr in enumerate(tv):
				E[str(yrPulse+0.5)][indV]=E[str(yrPulse+0.5)][indV]+E_PulseCH4
		# Add N2O
		if iScn==2:
			indV=np.where( (E['Variable']=='N2O') )[0]
			for i,yr in enumerate(tv):
				E[str(yrPulse+0.5)][indV]=E[str(yrPulse+0.5)][indV]+E_PulseN2O
		pd.DataFrame(E).to_csv(path_emissions_adj)
		# Add CO
		if iScn==3:
			indV=np.where( (E['Variable']=='CO') )[0]
			for i,yr in enumerate(tv):
				E[str(yrPulse+0.5)][indV]=E[str(yrPulse+0.5)][indV]+E_PulseCO
		pd.DataFrame(E).to_csv(path_emissions_adj)

		# Run adjusted emissions scenario
		d0[scn]=Run_FAIR(meta,dfConfig,dfVolcanic,path_emissions_adj)

	# Fit within timeframe of modelling project
	d1={}
	ind1=np.where( (d0['Reference']['Time']>=tv[0]) & (d0['Reference']['Time']<=tv[-1]) )[0]
	ind2=np.where( (tv>=d0['Reference']['Time'][0]) & (tv<=d0['Reference']['Time'][-1]) )[0]

	# Reference RF (FAIR prediction without any adjustement)
	d1['Reference']={}
	d1['Reference']['RF']={}
	for v in d0['Reference']['ssp370']['RF'].keys():
		d1['Reference']['RF'][v]={}
		d1['Reference']['RF'][v]=np.zeros(tv.size)
		d1['Reference']['RF'][v][ind2]=np.nan_to_num(d0['Reference']['ssp370']['RF'][v][ind1])*meta['Param']['BE']['Biophysical']['Nano Watts per Watt']#*meta['Param']['BE']['Biophysical']['Square meters per Hectare']
	d1['Reference']['RF Total']=np.zeros(tv.size)
	for v in d0['Reference']['ssp370']['RF'].keys():
		d1['Reference']['RF Total']=d1['Reference']['RF Total']+np.nan_to_num(d1['Reference']['RF'][v])
	d1['Reference']['Temperature']=np.zeros(tv.size)
	d1['Reference']['Temperature'][ind2]=d0['Reference']['ssp370']['Temperature'][ind1]

	# Add scenario RF
	for iScn,scn in enumerate(scnName):
		d1[scn]={'RF':{}}
		for v in d0[scn]['ssp370']['RF'].keys():
			d1[scn]['RF'][v]={}
			d1[scn]['RF'][v]=np.zeros(tv.size)
			d1[scn]['RF'][v][ind2]=np.nan_to_num(d0[scn]['ssp370']['RF'][v][ind1])*meta['Param']['BE']['Biophysical']['Nano Watts per Watt']
		d1[scn]['RF Total']=np.zeros(tv.size)
		for v in d0[scn]['ssp370']['RF'].keys():
			d1[scn]['RF Total']=d1[scn]['RF Total']+np.nan_to_num(d1[scn]['RF'][v])
		d1[scn]['Temperature']=np.zeros(tv.size)
		d1[scn]['Temperature'][ind2]=d0[scn]['ssp370']['Temperature'][ind1]

	return d1

#%% Benchmark pulse from FAIR
def FAIR_BenchmarkPulse(meta,pNam,YearPulse):
	#YearPulse=2020
	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	dFAIR_Pulse=CalcPulseEmissionBenchmark_RF_FAIR(meta,pNam,YearPulse)

	cl=np.array([[0.5,0.9,0],[0.5,0,1],[0.24,0.47,0.74],[1,0,0],[0.25,0,0],[1,1,0],[0,1,1],[0,0,1],[0,1,1],[1,0.6,0.7],[0,0.25,0]])
	xlim=[-5,110]
	iT=np.where( (tv>=YearPulse-5) )[0]
	plt.close('all'); fig,ax=plt.subplots(3,2,figsize=gu.cm2inch(12,12)); lw=0.5; lw2=1
	for i,f in enumerate(meta['Modules']['FAIR']['Forcings']):
		y=(dFAIR_Pulse['CO2 Pulse']['RF'][f][iT]-dFAIR_Pulse['Reference']['RF'][f][iT])/1e9
		ax[0,0].plot(tv[iT]-YearPulse,y,'k-',color=cl[i,:],lw=lw,label=f)
		if f=='Land use':
			yLU=y
		if f=='CO2':
			yCO2=y
	print('CO2 fertilization = ' + str(np.round(np.mean(yLU)/np.mean(yCO2),decimals=2)))
	y=(dFAIR_Pulse['CO2 Pulse']['RF Total']-dFAIR_Pulse['Reference']['RF Total'])/1e9
	ax[0,0].plot(tv[iT]-YearPulse,y[iT],'k-',color='k',lw=lw2,label='Total')
	ax[0,0].legend(loc='upper right',facecolor=[1,1,1],frameon=True,fontsize=4,ncol=2)
	ax[0,0].set(ylabel='$\Delta$RF (nW m$^{-2}$ ha$^{-1}$)',xlabel='Time since pulse emission (years)',xlim=xlim)
	y=dFAIR_Pulse['CO2 Pulse']['Temperature']-dFAIR_Pulse['Reference']['Temperature']
	ax[0,1].plot(tv[iT]-YearPulse,y[iT],'k-',lw=lw2)
	ax[0,1].set(ylabel=r'$\Delta$ temperature (K)',xlabel='Time since pulse emission (years)',xlim=xlim)

	for i,f in enumerate(meta['Modules']['FAIR']['Forcings']):
		y=(dFAIR_Pulse['CH4 Pulse']['RF'][f][iT]-dFAIR_Pulse['Reference']['RF'][f][iT])/1e9
		ax[1,0].plot(tv[iT]-YearPulse,y,'k-',color=cl[i,:],lw=lw,label=f)
	y=(dFAIR_Pulse['CH4 Pulse']['RF Total']-dFAIR_Pulse['Reference']['RF Total'])/1e9
	ax[1,0].plot(tv[iT]-YearPulse,y[iT],'k-',color='k',lw=lw2,label='Total')
	ax[1,0].set(ylabel='$\Delta$RF (nW m$^{-2}$ ha$^{-1}$)',xlabel='Time since pulse emission (years)',xlim=xlim)

	y=dFAIR_Pulse['CH4 Pulse']['Temperature']-dFAIR_Pulse['Reference']['Temperature']
	ax[1,1].plot(tv[iT]-YearPulse,y[iT],'k-',lw=lw2)
	ax[1,1].set(ylabel=r'$\Delta$ temperature (K)',xlabel='Time since pulse emission (years)',xlim=xlim)

	for i,f in enumerate(meta['Modules']['FAIR']['Forcings']):
		y=(dFAIR_Pulse['N2O Pulse']['RF'][f][iT]-dFAIR_Pulse['Reference']['RF'][f][iT])/1e9
		ax[2,0].plot(tv[iT]-YearPulse,y,'k-',color=cl[i,:],lw=lw,label=f)
	y=(dFAIR_Pulse['N2O Pulse']['RF Total']-dFAIR_Pulse['Reference']['RF Total'])/1e9
	ax[2,0].plot(tv[iT]-YearPulse,y[iT],'k-',color='k',lw=lw2,label='Total')
	ax[2,0].set(ylabel='$\Delta$RF (nW m$^{-2}$ ha$^{-1}$)',xlabel='Time since pulse emission (years)',xlim=xlim)

	y=dFAIR_Pulse['N2O Pulse']['Temperature']-dFAIR_Pulse['Reference']['Temperature']
	ax[2,1].plot(tv[iT]-YearPulse,y[iT],'k-',lw=lw2)
	ax[2,1].set(ylabel=r'$\Delta$ temperature (K)',xlabel='Time since pulse emission (years)',xlim=xlim)

	gu.axletters(ax,plt,0.48,0.87,FontColor=meta['Graphics']['gp']['cla'],
			LetterStyle='Default',FontWeight='Bold',
			Labels=['1 tonne CO$_2$ pulse','1 tonne CO$_2$ pulse','1 tonne CH$_4$ pulse','1 tonne CH$_4$ pulse','1 tonne N$_2$O pulse','1 tonne N$_2$O pulse'],LabelSpacer=0.055)
	plt.tight_layout()
	gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\FAIR\FAIR_Benchmark_PulseEmissions','png',900)
	return

#%%
def Calc_RF_FAIR(meta,pNam,mos):

	# Only populating All and only populating mean (not applied to sum yet)
	iPS=0
	iSS=0
	iYS=0
	iOS=0
	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)

	# Multiplier that will be applied to project emissions because single stand
	# emissions are below the precicion. This method assumes FAIR climate sensitivity
	# is independent of magintude
	sf_fair=1#e9

	# Modelling in the time period for CMIP6
	d0={}

	# Import FAIR data
	dfConfig=pd.read_csv(meta['Paths']['DB']['FAIR'] + '//Variables_FAIR_4xCO2_cummins_ebm3.csv')
	dfVolcanic=pd.read_csv(meta['Paths']['DB']['FAIR'] + '//Variables_FAIR_volcanic_ERF_monthly_175001-201912.csv',index_col='year')

	# Run the reference (CMIP6-SSP370)
	path_emissions_ref=meta['Paths']['DB']['FAIR'] + '//Variables_FAIR_Emissions_SSP370.csv'
	d0['Reference']=Run_FAIR(meta,dfConfig,dfVolcanic,path_emissions_ref)

	# Run each scenario
	for iScn in range(len(mos[pNam]['Scenarios'])):

		# GHG emissions from cbrunner adjsuted to units expected by FAIR
		E_co2=post.GetMosScnVar(meta,pNam,mos,iScn,'E_CO2',iPS,iSS,iYS,iOS)['Mean']['Ensemble Mean']/1e9 # GtCO2/year
		E_ch4=post.GetMosScnVar(meta,pNam,mos,iScn,'E_CH4',iPS,iSS,iYS,iOS)['Mean']['Ensemble Mean']/1e6 # MtCH4/year
		E_n2o=post.GetMosScnVar(meta,pNam,mos,iScn,'E_N2O',iPS,iSS,iYS,iOS)['Mean']['Ensemble Mean']/1e6 # MtN2O/year

		# Adjust by scale factor so that it is above the precision of FAIR
		E_co2=E_co2*sf_fair
		E_ch4=E_ch4*sf_fair
		E_n2o=E_n2o*sf_fair

		path_emissions_adj=meta['Paths']['DB']['FAIR'] + '//Variables_FAIR_Emissions_SSP370_adjusted.csv'
		shutil.copyfile(path_emissions_ref,path_emissions_adj)
		E=gu.ReadCSV(meta['Paths']['DB']['FAIR'] + '//Variables_FAIR_Emissions_SSP370_adjusted.csv',skip_rows=0)
		#indV=np.where( (E['Variable']=='CO2') | (E['Variable']=='CO2 AFOLU') )[0]
		indV=np.where( (E['Variable']=='CO2 AFOLU') )[0]
		for i,yr in enumerate(tv):
			if yr>2149:
				continue
			E[str(yr+0.5)][indV]=E[str(yr+0.5)][indV]+E_co2[i]
		indV=np.where( (E['Variable']=='CH4') )[0]
		for i,yr in enumerate(tv):
			if yr>2149:
				continue
			E[str(yr+0.5)][indV]=E[str(yr+0.5)][indV]+E_ch4[i]
		indV=np.where( (E['Variable']=='N2O') )[0]
		for i,yr in enumerate(tv):
			if yr>2149:
				continue
			E[str(yr+0.5)][indV]=E[str(yr+0.5)][indV]+E_n2o[i]
		pd.DataFrame(E).to_csv(path_emissions_adj)

		# Run adjusted emissions scenario
		d0[iScn]=Run_FAIR(meta,dfConfig,dfVolcanic,path_emissions_adj)

	# Fit within timeframe of modelling project
	d1={}
	ind1=np.where( (d0['Reference']['Time']>=tv[0]) & (d0['Reference']['Time']<=tv[-1]) )[0]
	ind2=np.where( (tv>=d0['Reference']['Time'][0]) & (tv<=d0['Reference']['Time'][-1]) )[0]

	# RF (nW m-2 ha-1)
	# Reference RF (FAIR prediction without any adjustement)
	d1['Reference']={}
	d1['Reference']['RF']={}
	for v in d0['Reference']['ssp370']['RF'].keys():
		d1['Reference']['RF'][v]={}
		d1['Reference']['RF'][v]=np.zeros(tv.size)
		d1['Reference']['RF'][v][ind2]=np.nan_to_num(d0['Reference']['ssp370']['RF'][v][ind1])*meta['Param']['BE']['Biophysical']['Nano Watts per Watt']
		#d1['Reference']['RF'][v][ind3]=0
	d1['Reference']['RF Total']=np.zeros(tv.size)
	for v in d0['Reference']['ssp370']['RF'].keys():
		d1['Reference']['RF Total'][ind2]=d1['Reference']['RF Total'][ind2]+np.nan_to_num(d0['Reference']['ssp370']['RF'][v][ind1])*meta['Param']['BE']['Biophysical']['Nano Watts per Watt']

	# Add scenario RF
	for iScn in range(len(mos[pNam]['Scenarios'])):
		d1[iScn]={'RF':{}}
		for v in d0[iScn]['ssp370']['RF'].keys():
			d1[iScn]['RF'][v]={}
			d1[iScn]['RF'][v]=np.zeros(tv.size)
			d1[iScn]['RF'][v][ind2]=np.nan_to_num(d0[iScn]['ssp370']['RF'][v][ind1])*meta['Param']['BE']['Biophysical']['Nano Watts per Watt']
			#d1[iScn]['RF'][v][ind3]=0

		d1[iScn]['RF Total']=np.zeros(tv.size)
		for v in d0[iScn]['ssp370']['RF'].keys():
			d1[iScn]['RF Total']=d1[iScn]['RF Total']+np.nan_to_num(d0[iScn]['ssp370']['RF'][v][ind1])*meta['Param']['BE']['Biophysical']['Nano Watts per Watt']

		for v in d0[iScn]['ssp370']['RF'].keys():
			mos[pNam]['Scenarios'][iScn]['Mean']['RF_Chem_' + v]={'Ensemble Mean':d1[iScn]['RF'][v]}
		mos[pNam]['Scenarios'][iScn]['Mean']['RF_Chem_Total']={'Ensemble Mean':d1[iScn]['RF Total']}

	# Calculate delta (nW m-2 ha-1)
	d2={}
	for scnC in mos[pNam]['Delta'].keys():
		d2[scnC]={}
		iB=mos[pNam]['Delta'][scnC]['iB']
		iP=mos[pNam]['Delta'][scnC]['iP']

		for v in d0[iScn]['ssp370']['RF'].keys():
			mos[pNam]['Delta'][scnC]['Data']['Mean']['RF_Chem_' + v]={'Ensemble Mean':d1[iP]['RF'][v]-d1[iB]['RF'][v]}
		mos[pNam]['Delta'][scnC]['Data']['Mean']['RF_Chem_Total']={'Ensemble Mean':d1[iP]['RF Total']-d1[iB]['RF Total']}

		#d2[scnC]['RF Total']=d1[iP]['RF Total']-d1[iB]['RF Total']
		#d2[scnC]['RF Total']=d1[iP]['RF Total']-d1[iB]['RF Total']
		#d2[scnC]['RF Total']=d1[iP]['RF Total']-d1[iB]['RF Total']
		#d2[scnC]['RF Total']=d1[iP]['RF Total']-d1[iB]['RF Total']

	#d={}
	#d['Scenario']=d1
	#d['Delta']=d2
	return mos

#%% Function for predicting albedo response to harvest
def Update_AlbedoRF_FromHarvest(tv,t_harv,beta,zone,var):
	if var=='Albedo':
		yhat0=beta[zone][0]+beta[zone][1]*(tv-t_harv)
		yhat=beta[zone][2]*np.ones(tv.size)
		yhat[tv>=t_harv]=np.maximum(beta[zone][2],yhat0[tv>=t_harv])
		yhat=yhat-yhat[0]
	else:
		yhat0=beta[zone][0]+beta[zone][1]*(tv-t_harv)
		yhat=beta[zone][2]*np.ones(tv.size)
		yhat[tv>=t_harv]=np.minimum(beta[zone][2],yhat0[tv>=t_harv])
		yhat=yhat-yhat[0]

	return yhat

#%%
def Plot_AlbedoResponseToHarvest(meta):
	beta_A=gu.ipickle(meta['Paths']['Model']['Parameters'] + '\\Parameters_AlbedoSurfaceShortwave_HarvestResponseByBGCZone.pkl')
	beta_ARF=gu.ipickle(meta['Paths']['Model']['Parameters'] + '\\Parameters_AlbedoSurfaceShortwaveRF_HarvestResponseByBGCZone.pkl')
	tv=np.arange(2000,2141,1)
	t_harv=2022
	lw=1.25
	
	plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(20,6))
	# Albedo
	var='Albedo'
	ax[0].plot([-10,200],[0,0],'-',lw=2,color=[0.8,0.8,0.8])
	zone='CWH'
	ax[0].plot(tv-t_harv,Update_AlbedoRF_FromHarvest(tv,t_harv,beta_A,zone,var),'-g',lw=lw,color=[0.5,1,0],label=zone)
	zone='ICH'
	ax[0].plot(tv-t_harv,Update_AlbedoRF_FromHarvest(tv,t_harv,beta_A,zone,var),'--g',lw=lw,color=[0,0.85,0],label=zone)
	zone='SBS'
	ax[0].plot(tv-t_harv,Update_AlbedoRF_FromHarvest(tv,t_harv,beta_A,zone,var),'-.c',lw=lw,label=zone)
	zone='IDF'
	ax[0].plot(tv-t_harv,Update_AlbedoRF_FromHarvest(tv,t_harv,beta_A,zone,var),':r',lw=lw,label=zone)
	zone='BWBS'
	ax[0].plot(tv-t_harv,Update_AlbedoRF_FromHarvest(tv,t_harv,beta_A,zone,var),'-b',color=[0,0,0.5],lw=lw,label=zone)
	ax[0].set(ylabel='$\Delta$ surface shortwave albedo',xlabel='Time since harest, years',xlim=[-10,85],xticks=np.arange(-10,90,10))
	ax[0].legend(loc='upper right',frameon=False,facecolor=None,edgecolor='w',fontsize=7)
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].tick_params(length=meta['Graphics']['gp']['tickl'])
	# Albedo RF
	var='Albedo RF'
	ax[1].plot([-10,200],[0,0],'-',lw=2,color=[0.8,0.8,0.8])
	zone='CWH'
	ax[1].plot(tv-t_harv,Update_AlbedoRF_FromHarvest(tv,t_harv,beta_ARF,zone,var),'-g',lw=lw,color=[0.5,1,0],label=zone)
	zone='ICH'
	ax[1].plot(tv-t_harv,Update_AlbedoRF_FromHarvest(tv,t_harv,beta_ARF,zone,var),'--g',lw=lw,color=[0,0.85,0],label=zone)
	zone='SBS'
	ax[1].plot(tv-t_harv,Update_AlbedoRF_FromHarvest(tv,t_harv,beta_ARF,zone,var),'-.c',lw=lw,label=zone)
	zone='IDF'
	ax[1].plot(tv-t_harv,Update_AlbedoRF_FromHarvest(tv,t_harv,beta_ARF,zone,var),':r',lw=lw,label=zone)
	zone='BWBS'
	ax[1].plot(tv-t_harv,Update_AlbedoRF_FromHarvest(tv,t_harv,beta_ARF,zone,var),'-b',color=[0,0,0.5],lw=lw,label=zone)
	
	ax[1].set(ylabel='$\Delta$ albedo RF (nW m$^{-2}$ ha$^{-1}$)',xlabel='Time since harest, years',xlim=[-10,85],xticks=np.arange(-10,90,10))
	#ax[1].legend(loc='lower right',frameon=False,facecolor=None,edgecolor='w',fontsize=7)
	ax[1].yaxis.set_ticks_position('both'); ax[1].xaxis.set_ticks_position('both'); ax[1].tick_params(length=meta['Graphics']['gp']['tickl'])
	#gu.axletters(ax,plt,0.015,0.88,FontColor=meta['Graphics']['gp']['cla'],LetterStyle=meta['Graphics']['Modelling']['AxesLetterStyle'],FontWeight=meta['Graphics']['Modelling']['AxesFontWeight'],Labels=['Cold season','Warm season'],LabelSpacer=0.02)
	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\AlbedoRF_ResponseToHarvesting','png',900)
	return

#%%
def DownloadAlbedo_Monthly(tvm):
	for i in range(tvm.shape[0]):
		dt1=str(tvm[i,0]) + '-' + str(tvm[i,1]) + '-01'
		dt2=str(tvm[i,0]) + '-' + str(tvm[i,1]) + '-28'
		y0=ee.ImageCollection("MODIS/061/MCD43A3").select('Albedo_BSA_shortwave').filterDate(dt1,dt2)

		# Mean
		y1=y0.reduce(ee.Reducer.mean())
	
		# Clip
		aoi=ee.Geometry.BBox(-139,47.25,-113,60.1) # Province wide
		#aoi=ee.Geometry.BBox(-122.15182314190407,50.51206023665638,-120.0185403010095,51.6781353097244) # Elephant hill fire
		y2=y1.clip(aoi)

		# Reproject
		crs='EPSG:3005'
		scale=500
		y3=y2.reproject(crs=crs,scale=scale)

		# Define task
		task=ee.batch.Export.image.toDrive(image=y3,description='albedo_' + str(tvm[i,0]) + '_' + str(tvm[i,1]),scale=scale,maxPixels=1e13,fileFormat="GeoTIFF")
		task.start()
	
	# Status
	# task.status()
	return

#%% Reproject albedo actual monthly variability
def Reproject_Albedo():
	sf_given=0.001
	for i in range(tvm.shape[0]):
		z=gis.OpenGeoTiff(r'D:\Albedo\Province\albedo_' + str(tvm[i,0]) + '_' + str(tvm[i,1]) + '.tif')
		z['Data']=z['Data']*sf_given
		z['Data']=z['Data']/meta['Climate']['SF']['Albedo']
		z['Data']=z['Data'].astype('int16')
		##plt.matshow(z['Data'],clim=[0,100])
		pth2=r'D:\Albedo\Province\albedo_' + str(tvm[i,0]) + '_' + str(tvm[i,1]) + 'b.tif'
		gis.SaveGeoTiff(z,pth2)
		#pth3=r'C:\Data\BC1ha\Climate\Monthly\Normals\bc1ha_albedo_surf_' + str(tvm[i,0]) + str(tvm[i,1]) + '.tif'
		pth3=r'D:\Albedo\Province\bc1ha_albedo_' + str(tvm[i,0]) + '_' + str(tvm[i,1]) + '.tif'
		gis.ReprojectRasterAndClipToRaster(pth2,pth3,meta['Paths']['bc1ha Ref Grid'],meta['Geos']['crs'])
	#z=gis.OpenGeoTiff(pth3);plt.matshow(z['Data'],clim=[0,1000])
	return

#%%
def Calc_Albedo_AnnualSummary(meta):
	meta=u1ha.Init()
	zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
	iMask=np.where(zRef['Data']==0)
	moL=[10,11,12,1,2,3]
	for yr in range(2002,2024):
		Ac=np.zeros(zRef['Data'].shape,dtype='int16')
		Aw=np.zeros(zRef['Data'].shape,dtype='int16')
		Aa=np.zeros(zRef['Data'].shape,dtype='int16')
		for mo in range(12):
			print(yr)
			print(mo+1)
			z=gis.OpenGeoTiff(r'D:\Albedo\Province\Albedo\Processed\bc1ha_albedo_' + str(yr) + '_' + str(mo+1) + '.tif')['Data'].astype('float')*meta['Climate']['SF']['Albedo'] #;plt.matshow(z['Data'],clim=[0,1000])

			Aa=Aa+z
			if np.isin(mo,moL)==True:
				Ac=Ac+z
			else:
				Aw=Aw+z

		# Corrections
		Aa=Aa/12
		Ac=Ac/6
		Aw=Aw/6

		z=copy.deepcopy(zRef)
		z['Data']=Aa
		z['Data']=z['Data']/meta['Climate']['SF']['Albedo']
		z['Data']=z['Data'].astype('int16')
		z['Data'][iMask]=0
		gis.SaveGeoTiff(z,r'D:\Albedo\Province\Albedo\Processed\bc1ha_albedo_ann_' + str(yr) + '.tif')

		z=copy.deepcopy(zRef)
		z['Data']=Ac
		z['Data']=z['Data']/meta['Climate']['SF']['Albedo']
		z['Data']=z['Data'].astype('int16')
		z['Data'][iMask]=0
		gis.SaveGeoTiff(z,r'D:\Albedo\Province\Albedo\Processed\bc1ha_albedo_coldsea_' + str(yr) + '.tif')

		z=copy.deepcopy(zRef)
		z['Data']=Aw
		z['Data']=z['Data']/meta['Climate']['SF']['Albedo']
		z['Data']=z['Data'].astype('int16')
		z['Data'][iMask]=0
		gis.SaveGeoTiff(z,r'D:\Albedo\Province\Albedo\Processed\bc1ha_albedo_warmsea_' + str(yr) + '.tif')

	return

#%%
def Calc_RSW_Absorption_Annual(meta):
	zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
	con=uga.HydroMetCon()
	tvm=gu.tvec('m',2002,2023)
	moL=[10,11,12,1,2,3]
	for yr in range(2002,2024):
		ArswC=np.zeros(zRef['Data'].shape,dtype='int16')
		ArswW=np.zeros(zRef['Data'].shape,dtype='int16')
		ArswA=np.zeros(zRef['Data'].shape,dtype='int16')
		for mo in range(12):
			zA=gis.OpenGeoTiff(r'D:\Albedo\Province\Albedo\Processed\bc1ha_albedo_' + str(yr) + '_' + str(mo+1) + '.tif') #;plt.matshow(z['Data'],clim=[0,1000])
			zA['Data']=zA['Data'].astype('float')*meta['Climate']['SF']['Albedo']

			zRg=gis.OpenGeoTiff(r'C:\Data\BC1ha\Climate\Monthly\Normals\bc1ha_rswd_norm_1971to2000_' + str(mo+1) + '.tif') #;plt.matshow(z['Data'],clim=[0,1000])
			zRg['Data']=zRg['Data'].astype('float')*meta['Climate']['SF']['rswd']

			# Convert net shortwave radiation from MJ m-2 d-1 to W m-2
			zRg['Data']=zRg['Data']*1e6/con['DayLength']

			Beta=0.23 # Bernier et al. 2011, Montenegro et al. 2009
			Arsw=(1-zA['Data'])*zRg['Data'] + Beta*zRg['Data']*zA['Data']

			ArswA=ArswA+Arsw
			if np.isin(mo,moL)==True:
				ArswC=ArswC+Arsw
			else:
				ArswW=ArswW+Arsw
		ArswA=ArswA/12
		ArswC=ArswC/6
		ArswW=ArswW/6
	
		z=copy.deepcopy(zA)
		z['Data']=ArswA
		z['Data']=z['Data']/meta['Climate']['SF']['AbsorptionRSW']
		z['Data']=z['Data'].astype('int16')
		gis.SaveGeoTiff(z,r'D:\Albedo\Province\AbsorptionRSW\bc1ha_AbsorptionRSW_ann_' + str(yr) + '.tif')
	
		z=copy.deepcopy(zA)
		z['Data']=ArswC
		z['Data']=z['Data']/meta['Climate']['SF']['AbsorptionRSW']
		z['Data']=z['Data'].astype('int16')
		gis.SaveGeoTiff(z,r'D:\Albedo\Province\AbsorptionRSW\bc1ha_AbsorptionRSW_coldsea_' + str(yr) + '.tif')
	
		z=copy.deepcopy(zA)
		z['Data']=ArswW
		z['Data']=z['Data']/meta['Climate']['SF']['AbsorptionRSW']
		z['Data']=z['Data'].astype('int16')
		gis.SaveGeoTiff(z,r'D:\Albedo\Province\AbsorptionRSW\bc1ha_AbsorptionRSW_warmsea_' + str(yr) + '.tif')
	return

#%%
def Calc_AlbedoRF_Climatology_ByBGCZone(meta,zG,con,ivl):
	nanoWperW=1e9
	sqmperha=1e4

	ikp={}
	d={}
	for zone in meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE'].keys():
		ikp[zone]=np.where( (zG['lc_comp1_2019']==1) & (zG['bgcz']==meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE'][zone]) )
		d[zone]={}
		d[zone]['Albedo']=np.zeros(12)
		d[zone]['RSWD']=np.zeros(12)
		d[zone]['AbsorptionRSW']=np.zeros(12)
		d[zone]['Albedo RF']=np.zeros(12)

	for mo in range(12):
		zA=gis.OpenGeoTiff(r'C:\Data\BC1ha\Climate\Monthly\Normals\bc1ha_albedo_surf_norm_2002to2023_' + str(mo+1) + '.tif')['Data'][0::ivl,0::ivl].astype('float')*meta['Climate']['SF']['Albedo']

		# Convert net shortwave radiation from MJ m-2 d-1 to W m-2
		zRg=gis.OpenGeoTiff(r'C:\Data\BC1ha\Climate\Monthly\Normals\bc1ha_rswd_norm_1971to2000_' + str(mo+1) + '.tif')['Data'][0::ivl,0::ivl].astype('float')*meta['Climate']['SF']['rswd']
		zRg=zRg*1e6/con['DayLength']
	
		Beta=0.23 # Bernier et al. 2011, Montenegro et al. 2009
		Arsw=(1-zA)*zRg + Beta*zRg*zA

		for zone in meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE'].keys():
			d[zone]['Albedo'][mo]=np.mean(zA[ikp[zone]])
			d[zone]['RSWD'][mo]=np.mean(zRg[ikp[zone]])
			d[zone]['AbsorptionRSW'][mo]=np.mean(Arsw[ikp[zone]])
			d[zone]['Albedo RF'][mo]=np.mean(Arsw[ikp[zone]]/meta['Param']['BE']['Biophysical']['Area_EarthSurface']*nanoWperW*sqmperha)

	return d

#%%
def Plot_Albedo_Climatology_ByBGCZone(meta):

	d=gu.ipickle(meta['Paths']['DB']['SurfaceClimate'] + '\\Variables_AlbedoClimatology_ByBGCZone.pkl')

	plt.close('all'); fig,ax=plt.subplots(1,3,figsize=gu.cm2inch(22,6)); lw=1
	v='Albedo'
	zone='CWH'
	ax[0].plot(d[zone][v],'-ko',lw=lw,color=[0.5,1,0],label=zone)
	zone='SBS'
	ax[0].plot(d[zone][v],'-ks',lw=lw,color=[0,0,0.6],label=zone)
	ax[0].set(ylabel='Albedo, surface shortwave',xlabel='Month',xticks=np.arange(12),xticklabels=['Jan','','Mar','','May','','Jul','','Sep','','Nov',''],xlim=[-1,12])
	ax[0].legend(loc='lower left',frameon=False,facecolor=None,edgecolor='w',fontsize=7)
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].tick_params(length=meta['Graphics']['gp']['tickl'])
	
	v='RSWD'
	zone='CWH'
	ax[1].plot(d[zone][v],'-ko',lw=lw,color=[0.5,1,0],label=zone)
	zone='SBS'
	ax[1].plot(d[zone][v],'-ks',lw=lw,color=[0,0,0.6],label=zone)
	ax[1].set(ylabel='Radiation, shortwave down (W m$^{-2}$)',xlabel='Month',xticks=np.arange(12),xticklabels=['Jan','','Mar','','May','','Jul','','Sep','','Nov',''],xlim=[-1,12])
	ax[1].yaxis.set_ticks_position('both'); ax[1].xaxis.set_ticks_position('both'); ax[1].tick_params(length=meta['Graphics']['gp']['tickl'])
	
	v='Albedo RF'
	zone='CWH'
	ax[2].plot(d[zone][v],'-ko',lw=lw,color=[0.5,1,0],label=zone)
	zone='SBS'
	ax[2].plot(d[zone][v],'-ks',lw=lw,color=[0,0,0.6],label=zone)
	ax[2].set(ylabel='Radiative forcing (nW m$^{-2}$ ha$^{-1}$)',xlabel='Month',xticks=np.arange(12),xticklabels=['Jan','','Mar','','May','','Jul','','Sep','','Nov',''],xlim=[-1,12])
	ax[2].yaxis.set_ticks_position('both'); ax[2].xaxis.set_ticks_position('both'); ax[2].tick_params(length=meta['Graphics']['gp']['tickl'])
	gu.axletters(ax,plt,0.0175,0.94,FontColor=meta['Graphics']['gp']['cla'],LetterStyle=meta['Graphics']['Modelling']['AxesLetterStyle'],FontWeight=meta['Graphics']['Modelling']['AxesFontWeight'])
	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\AlbedoRF_ResponseToHarvesting','png',900)
	return

#%%
def wy(x):
	y=np.append(x[10:],x[0:10])
	return y

#%%
def AerodynamicConductance(H_Canopy,U):

	# Constants
	con=HydroMetCon()

	# Reference height (m)
	zRef=H_Canopy+2

	# Zero-plane displacement height (m)
	H_ZPD=0.75*H_Canopy

	# Roughness length for heat (m)
	RL_H_Min=0.01
	RL_H_Max=1.00
	RL_H=RL_H_Min+(RL_H_Max-RL_H_Min)*(1-np.exp(-2.5*(H_Canopy/50)))

	# Roughness length for momentum (m)
	RL_M=RL_H

	flg=0
	if flg==1:
		H_Canopy=np.arange(0,50.5,0.5)
		RL_H=0.1+(1.2-0.1)*(1-np.exp(-2.5*(H_Canopy/50)))
		plt.close('all'); plt.plot(H_Canopy,RL_H,'b-')

	# Stabilty correction factors
	scf_M=0.6
	scf_H=0.6

	# Monin-Obukhov length (m)
	L_MO=-50

	# Terms
	y1=con['Karmans Constant']*U
	y2=np.log((zRef-H_ZPD)/RL_H)-scf_H*((zRef-H_ZPD)/L_MO)
	y3=np.log((zRef-H_ZPD)/RL_M)-scf_M*((zRef-H_ZPD)/L_MO)

	# Aerodynamic conductance (m s-1)
	y=y1/y2/y3

	return y