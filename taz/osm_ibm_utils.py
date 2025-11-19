'''
ONSET AND SPREAD MODEL UTILITIES
'''

#%% Import modules
import numpy as np
import scipy.stats as stats
from matplotlib import colors
from matplotlib import animation
import time
import copy
import cv2
import fcgadgets.macgyver.util_general as gu
import fcgadgets.macgyver.util_gis as gis
import fcgadgets.cbrunner.cbrun_util as cbu

#%%
def PercentPine(meta,z):
	pL=[ meta['LUT']['VEG_COMP_LYR_R1_POLY']['SPECIES_CD_1']['PL'],
		meta['LUT']['VEG_COMP_LYR_R1_POLY']['SPECIES_CD_1']['PLI'] ]

	z['Pct Pine 23']=np.zeros(z['age_vri23'].shape)
	ind=np.where( (np.isin(z['spc1_vri23'],pL)==True) )
	z['Pct Pine 23'][ind]=z['Pct Pine'][ind]+z['spc1_pct_vri23'][ind]
	ind=np.where( (np.isin(z['spc2_vri23'],pL)==True) )
	z['Pct Pine 23'][ind]=z['Pct Pine'][ind]+z['spc2_pct_vri23'][ind]
	ind=np.where( (np.isin(z['spc3_vri23'],pL)==True) )
	z['Pct Pine 23'][ind]=z['Pct Pine'][ind]+z['spc3_pct_vri23'][ind]
	ind=np.where( (np.isin(z['spc4_vri23'],pL)==True) )
	z['Pct Pine 23'][ind]=z['Pct Pine'][ind]+z['spc4_pct_vri23'][ind]
	z['Pct Pine Smoothed 23']=cv2.blur(z['Pct Pine 23'],(100,100)) # 5x5 kernel size

	z['Pct Pine 02']=np.zeros(z['age_vri02'].shape)
	ind=np.where( (np.isin(z['spc1_vri02'],pL)==True) )
	z['Pct Pine 02'][ind]=z['Pct Pine'][ind]+z['spc1_pct_vri02'][ind]
	ind=np.where( (np.isin(z['spc2_vri02'],pL)==True) )
	z['Pct Pine 02'][ind]=z['Pct Pine'][ind]+z['spc2_pct_vri02'][ind]
	ind=np.where( (np.isin(z['spc3_vri02'],pL)==True) )
	z['Pct Pine 02'][ind]=z['Pct Pine'][ind]+z['spc3_pct_vri02'][ind]
	ind=np.where( (np.isin(z['spc4_vri02'],pL)==True) )
	z['Pct Pine 02'][ind]=z['Pct Pine'][ind]+z['spc4_pct_vri02'][ind]
	z['Pct Pine Smoothed 02']=cv2.blur(z['Pct Pine 02'],(100,100)) # 5x5 kernel size

	#plt.matshow(z['Pct Pine'][0::5,0::5],clim=[0,100])
	#plt.matshow(z['Pct Pine Smoothed'][0::5,0::5],clim=[0,100])
	return z

#%%
def ConfigureProject(meta,roi):

	meta['Project']={}
	meta['Project']['Year']=np.arange(2025,2101)
	meta['Project']['N Time']=meta['Project']['Year'].size
	meta['Project']['N Ensemble']=1
	meta['Project']['N Scenario']=1
	meta['Project']['Grid Shape']=roi['grd']['Data'].shape
	meta['Project']['Save Daily Status']='On'

	# Parameters
	meta['Par']={}
	meta['Par']['neighbourhood']=((-1,1),(0,1),(1,1),(1,0),(1, -1),(0,-1),(-1,-1),(-1,0))
	#direction=['NW','N','NE','E','SE','S','SW','W']
	meta['Par']['UNSUSEPTABLE']=0
	meta['Par']['SUSEPTABLE']=1
	meta['Par']['ONSET']=2
	meta['Par']['Just Occurred']=3
	meta['Par']['non']=lambda s: s if s<0 else None
	meta['Par']['mom']=lambda s: max(0,s)
	
	meta['Par']['Season Length']=1
	meta['Par']['Duration Limit']=9
	meta['Par']['Prob Onset']=0.0005
	meta['Par']['Prob Spread']=0.18
	meta['Par']['Years To Regenerate']=25

	# Colours for visualization: brown for EMPTY, dark green for TREE and orange
	# for FIRE. Note that for the colormap to work, this list and the bounds list
	# must be one larger than the number of different values in the array.
	meta['Par']['colors_list']=[(0,0,0),(0.3,0.3,0.3),(1,0,0),(1,0.5,0)]
	meta['Par']['cmap']=colors.ListedColormap(meta['Par']['colors_list'])
	meta['Par']['bounds']=[0,1,2,3]
	meta['Par']['norm']=colors.BoundaryNorm(meta['Par']['bounds'],meta['Par']['cmap'].N)

	# Time series
	meta['TS']={}
	meta['TS']['Area Occurrence']=np.zeros(meta['Project']['N Time'])
	meta['TS']['Area Suseptable']=np.zeros(meta['Project']['N Time'])

	# Grids
	meta['Grids']={}
	#meta['Grids'][iScn]['Grid']['Prob Onset']=np.zeros((m,n),dtype=float)
	#meta['Grids'][iScn]['Grid']['Prob Spread']=np.zeros((m,n),dtype=float)

	# Initialize stats
	#meta['Grids'][iScn]['Area burned']=np.zeros((meta['Time'].size,meta['N Ensemble']),dtype=int)
	#meta['Grids'][iScn]['Area suseptable']=np.zeros((meta['Time'].size,meta['N Ensemble']),dtype=int)

	return meta

#%% Annual loop
def AnnualLoop(meta,zI):

	t0=time.time()

	# Unpack parameters
	mom=meta['Par']['mom']
	non=meta['Par']['non']

	zO={}

	# Initialize
	grid_size=(meta['Project']['N Time'],meta['Project']['Grid Shape'][0],meta['Project']['Grid Shape'][1])
	zO['Occurrence']=np.zeros(grid_size,dtype='int8')
	zO['Status']=meta['Par']['SUSEPTABLE']*np.ones(meta['Project']['Grid Shape'],dtype='int8')
	zO['Years Since Occurrence']=(meta['Par']['Years To Regenerate']+1)*np.ones(meta['Project']['Grid Shape'],dtype='int8')
	zO['Duration']=np.zeros(meta['Project']['Grid Shape'] ,dtype='int8')

	zO['P Onset']=meta['Par']['Prob Onset']*np.ones(meta['Project']['Grid Shape'],dtype='float64')
	zO['P Spread']=meta['Par']['Prob Spread']*np.ones(meta['Project']['Grid Shape'],dtype='float64')

	# Annual loop
	for iYear in range(meta['Project']['Year'].size):
		print(iYear)

		# Onset
		rn_onset=np.random.random( meta['Project']['Grid Shape'] )
		ind=(zO['Status']==meta['Par']['SUSEPTABLE']) & (rn_onset<zO['P Onset'])
		zO['Status'][ind]=meta['Par']['ONSET']
		zO['Years Since Occurrence'][ind]=0

		# Spread
		rn_spread=np.random.random( meta['Project']['Grid Shape'] )
		daily_weather_factor=np.random.normal(loc=0,scale=0.1,size=1)
		for dx,dy in meta['Par']['neighbourhood']:
			z_shift=zO['Status'].copy()
			z_shift[mom(dy):non(dy),mom(dx):non(dx)]=zO['Status'][mom(-dy):non(-dy),mom(-dx):non(-dx)]
			
			ind=(zO['Status']==meta['Par']['SUSEPTABLE']) & (z_shift==meta['Par']['ONSET']) & (rn_spread<(zO['P Spread']+daily_weather_factor))
			zO['Status'][ind]=meta['Par']['ONSET']
			zO['Years Since Occurrence'][ind]=0

		# Track duration of active onset
		ind=(zO['Status']==meta['Par']['ONSET'])
		zO['Duration'][ind]=zO['Duration'][ind]+1

		# Convert ONSET to UNSUSCEPTABLE
		ind=(zO['Duration']>=meta['Par']['Duration Limit'])
		zO['Status'][ind]=meta['Par']['UNSUSEPTABLE']

		# Record occurrence
		ind=(zO['Duration']>0)
		zO['Occurrence'][iYear,ind]=1

		# Update the years since occurrence counter
		ind=(zO['Years Since Occurrence']<meta['Par']['Years To Regenerate'])
		zO['Years Since Occurrence'][ind]=zO['Years Since Occurrence'][ind]+1

		# Transfer affected stands back to the suseptable class upon regeneration
		ind=(zO['Years Since Occurrence']==meta['Par']['Years To Regenerate'])
		zO['Status'][ind]=meta['Par']['SUSEPTABLE']

		# Adjust counter to be just above threshold regeneration period so that
		# they will not be found in queries
		zO['Years Since Occurrence'][ind]=meta['Par']['Years To Regenerate']+1

		# Update stats
		meta['TS']['Area Occurrence'][iYear]=np.sum(zO['Occurrence'][iYear,:,:])
		meta['TS']['Area Suseptable'][iYear]=np.sum( (zO['Status']==meta['Par']['SUSEPTABLE']) )

		#print(meta['Time'][iYear],z['Aburn'][iYear])

		#if meta['Save Daily Status']=='On':
		#	return meta

	flg=0
	if flg==1:
		a=np.sum(zO['Occurrence'],axis=(1,2))
		plt.plot(a,'-bo')

	# Save
	#idx=np.where(zO['Occurrence']==1)
	#gu.opickle(meta['Paths']['Project'] + '\\Inputs\\Ensembles\\wf_sim_osm_Scn' + cbu.FixFileNum(iScn) + '_Ens' + cbu.FixFileNum(iEns) + '.pkl',idx)

	print( np.round((time.time()-t0)/60,decimals=1) )

	return meta

#%%
def DailyLoop(meta,z):

	if meta['Project']['Save Daily Status']=='On':
		meta['Project']['Daily Status']=np.zeros((meta['Par']['Season Length'],meta['Project']['Grid Shape'][0],meta['Project']['Grid Shape'][1]),dtype='int8')

	# At the beginning of the year, convert ONSET from the end of the previous year to UNSUSEPTABLE
	if meta['Par']['Prob Continue']>0:
		# Event can persist across years
		rn=np.random.random( meta['Project']['Grid Shape'] )
		ind=(z['Status']==meta['Par']['ONSET']) & (rn>=meta['Par']['Prob Continue'])
		z['Status'][ind]=meta['Par']['UNSUSEPTABLE']
	else:
		# No carryover across years
		ind=(z['Status']==meta['Par']['ONSET'])
		z['Status'][ind]=meta['Par']['UNSUSEPTABLE']

	# Unpack parameters
	mom=meta['Par']['mom']
	non=meta['Par']['non']

	# Initialize duration of active events
	z['Duration']=np.zeros( meta['Project']['Grid Shape'] ,dtype='int8' )

	for iDay in range(meta['Par']['Season Length']):

		# Onset
		rn_onset=np.random.random( meta['Project']['Grid Shape'] )
		ind=(z['Status']==meta['Par']['SUSEPTABLE']) & (rn_onset<z['P Onset'])
		z['Status'][ind]=meta['Par']['ONSET']
		z['Years Since Occurrence'][ind]=0

		# Spread
		rn_spread=np.random.random( meta['Project']['Grid Shape'] )
		daily_weather_factor=np.random.normal(loc=0,scale=0.1,size=1)
		for dx,dy in meta['Par']['neighbourhood']:
			z_shift=z['Status'].copy()
			z_shift[mom(dy):non(dy),mom(dx):non(dx)]=z['Status'][mom(-dy):non(-dy),mom(-dx):non(-dx)]
			
			ind=(z['Status']==meta['Par']['SUSEPTABLE']) & (z_shift==meta['Par']['ONSET']) & (rn_spread<(z['P Spread']+daily_weather_factor))
			z['Status'][ind]=meta['Par']['ONSET']
			z['Years Since Occurrence'][ind]=0

		# Track duration of active onset
		ind=(z['Status']==meta['Par']['ONSET'])
		z['Duration'][ind]=z['Duration'][ind]+1

		# Convert ONSET to UNSUSCEPTABLE (or JUST BURNED FOR GRAPHIC PURPOSES)
		ind=(z['Duration']>=meta['Par']['Duration Limit'])

		# Save daily data
		if meta['Project']['Save Daily Status']=='On':
			print('Working on day ' + str(iDay))
			meta['Project']['Daily Status'][iDay,:,:]=z['Status'].copy()
			meta['Project']['Daily Status'][iDay,ind]=meta['Par']['Just Occurred']

		# Convert ONSET to UNSUSCEPTABLE
		z['Status'][ind]=meta['Par']['UNSUSEPTABLE']

	return meta,z
