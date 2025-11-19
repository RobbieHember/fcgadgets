'''
Beetle occurrence modelling
'''

#%% Import python modules
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import time
from matplotlib import animation
from matplotlib import colors
import fcgadgets.macgyver.util_general as gu
import fcgadgets.macgyver.util_gis as gis
import fcgadgets.bc1ha.bc1ha_utils as u1ha
import fcgadgets.cbrunner.cbrun_util as cbu
import fcgadgets.bc5k.bc5k_util as u5k
import fcgadgets.taz.osm_ibm_utils as uosm

#%% Import paths and look-up-tables
# *** Run BC1ha Map ***
#meta=u1ha.Init()
#zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])

#%% Import spatial layers
vList=['ibm_yr',
	'age_vri02','spc1_vri02','spc2_vri02','spc3_vri02','spc4_vri02',
	'spc1_pct_vri02','spc2_pct_vri02','spc3_pct_vri02','spc4_pct_vri02',
	'sphl_vri02',
	'age_vri23','sphl_vri23','fire_yl',
	'spc1_vri23','spc2_vri23','spc3_vri23','spc4_vri23',
	'spc1_pct_vri23','spc2_pct_vri23','spc3_pct_vri23','spc4_pct_vri23']
#z=u1ha.Import_Raster(meta,[],vList,'Extract Grid')
roi=u1ha.Import_Raster(meta,roi,vList,'Extract Grid')

# Calculate percent pine
roi['grd']=uosm.PercentPine(meta,roi['grd'])

#%% Changes between 2002 and 2023
ind=np.where( (roi['grd']['ibm_yr']>1998) & (roi['grd']['Pct Pine 02']>75) & (roi['grd']['fire_yl']>=2015) )
#ind=np.where( (roi['grd']['ibm_yr']>=1998) & (roi['grd']['Pct Pine 02']>75) )
#ind=np.where( (roi['grd']['ibm_yr']>=1998) & (roi['grd']['Pct Pine 02']>=0) )

a=roi['grd']['Pct Pine 23']-roi['grd']['Pct Pine 02']
plt.close('all'); plt.hist(a[ind].flatten()[0::1000]) #plt.close('all'); plt.matshow(a)
print(np.mean(a[ind]))

a=roi['grd']['sphl_vri23']-roi['grd']['sphl_vri02']
plt.close('all'); plt.hist(a[ind].flatten()[0::1000])
print(np.mean(a[ind]))

a=roi['grd']['age_vri23']-roi['grd']['age_vri02']
plt.close('all'); plt.hist(a[ind].flatten()[0::100])
print(np.mean(a[ind]))

#%% Change in abundance of susceptable age

plt.close('all');

t=np.arange(2002,2101)
y=np.zeros(t.size)
for i,yr in enumerate(t):
	ind=np.where( (roi['grd']['Pct Pine 02']>50) & (roi['grd']['age_vri02']+i>=80) & (roi['grd']['age_vri02']+i<=160) )
	y[i]=ind[0].size
plt.plot(t,y/1e6,'-ko')

t=np.arange(2023,2101)
y=np.zeros(t.size)
for i,yr in enumerate(t):
	ind=np.where( (roi['grd']['Pct Pine 23']>50) & (roi['grd']['age_vri23']+i>=80) & (roi['grd']['age_vri23']+i<=160) )
	y[i]=ind[0].size
plt.plot(t,y/1e6,'-bo')


#%% Load environmental variables
env=u5k.LoadEnvironmentalVariables(meta)

#%% Temperature
plt.close('all')
TaH=np.mean(env['CRU']['Data']['tmean']['mjjas']*meta['Climate']['SF']['tmean'],axis=(1,2))
TaF=np.mean(env['ESM']['Data']['ssp245']['CESM2']['tmean']['mjjas']*meta['Climate']['SF']['tmean'],axis=(1,2))

plt.plot(env['CRU']['tv'],TaH,'b-')
iT=np.where(env['ESM']['tv']>1971)[0]
plt.plot(env['ESM']['tv'][iT],TaF[iT],'r-')
#y=np.mean(env['ESM BC']['Data']['ssp585']['CESM2']['tmean']['mjjas']*meta['Climate']['SF']['tmean'],axis=(1,2))
#plt.plot(env['ESM']['tv'][iT],y[iT],'r--')

#%% Temperature
def fT(x):
	y=np.exp(1.5*x)
	return y

# fT(1)

flg=0
if flg==1:
	x=np.arange(-1,2,0.01)
	plt.close('all'); plt.plot(x,fT(x),'b-')

#%% Configure project
#meta=uosm.ConfigureProject(meta,roi)
def Run(meta,roi):

	meta['Project']={}
	meta['Project']['Year']=np.arange(1500,2101)
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
	mom=meta['Par']['mom']
	non=meta['Par']['non']
	
	meta['Par']['Duration Limit']=4
	meta['Par']['Prob Onset']=0.000001*(roi['grd']['Pct Pine'].astype('float')/100)
	meta['Par']['Prob Spread']=0.18*(roi['grd']['Pct Pine Smoothed'].astype('float')/100)
	meta['Par']['Years Until Suceptable']=30

	meta['Par']['colors_list']=[(0,0,0),(0.3,0.3,0.3),(1,0,0),(1,0.5,0)]
	meta['Par']['cmap']=colors.ListedColormap(meta['Par']['colors_list'])
	meta['Par']['bounds']=[0,1,2,3]
	meta['Par']['norm']=colors.BoundaryNorm(meta['Par']['bounds'],meta['Par']['cmap'].N)

	zO={}
	zO['Occurrence']=np.zeros((meta['Project']['N Time'],meta['Project']['Grid Shape'][0],meta['Project']['Grid Shape'][1]),dtype='int8')
	zO['Status']=meta['Par']['SUSEPTABLE']*np.ones(meta['Project']['Grid Shape'],dtype='int8')
	zO['Years Since Occurrence']=(meta['Par']['Years Until Suceptable']+1)*np.ones(meta['Project']['Grid Shape'],dtype='int8')
	zO['Duration']=np.zeros(meta['Project']['Grid Shape'] ,dtype='int8')

	# Annual loop
	for iT,year in enumerate(meta['Project']['Year']):
		print(iT)

		# Temperature
		if year<1851:
			rn=np.random.randint(1851,1920)
			iT1=np.where( (env['CRU']['tv']==rn) )[0]
			Ta=TaH[iT1]
		elif (year>=1851) & (year<2023):
			iT1=np.where( (env['CRU']['tv']==year) )[0]
			Ta=TaH[iT1]
		elif (year>=2023):
			iT1=np.where( (env['ESM']['tv']==year) )[0]
			if iT1.size==0:
				Ta=TaF[-1]
			else:
				Ta=TaF[iT1]
				print(Ta)
		else:
			pass

		Po=meta['Par']['Prob Onset']*fT(Ta)
		Ps=meta['Par']['Prob Spread']*fT(Ta)

		zO['P Onset']=Po*np.ones(meta['Project']['Grid Shape'],dtype='float')
		zO['P Spread']=Ps*np.ones(meta['Project']['Grid Shape'],dtype='float')

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
		ind=np.where( (zO['Status']==meta['Par']['ONSET']) )
		zO['Duration'][ind]=zO['Duration'][ind]+1
		zO['Occurrence'][iT,ind]=1

		# Convert ONSET to UNSUSCEPTABLE
		ind=np.where( (zO['Duration']>=meta['Par']['Duration Limit']) )
		zO['Status'][ind]=meta['Par']['UNSUSEPTABLE']
		zO['Occurrence'][iT,ind]=0

		# Update the years since occurrence counter
		ind=np.where( (zO['Years Since Occurrence']<meta['Par']['Years Until Suceptable']) )
		zO['Years Since Occurrence'][ind]=zO['Years Since Occurrence'][ind]+1

		# Transfer affected stands back to the suseptable class upon regeneration
		ind=np.where( (zO['Years Since Occurrence']==meta['Par']['Years Until Suceptable']) )
		zO['Status'][ind]=meta['Par']['SUSEPTABLE']

		# Adjust counter to be just above threshold regeneration period so that
		# they will not be found in queries
		zO['Years Since Occurrence'][ind]=meta['Par']['Years Until Suceptable']+1

	flg=0
	if flg==1:
		ind=np.where( roi['grd']['Pct Pine']>0 )
		N_tot=ind[0].size
		a=np.sum(zO['Occurrence'],axis=(1,2))
		plt.close('all'); plt.plot(meta['Project']['Year'],a/N_tot*100,'-bo')

		plt.close('all'); plt.matshow(zO['Status'])
		plt.close('all'); plt.matshow(zO['Duration'])
		plt.close('all'); plt.matshow(zO['Years Since Occurrence'])

	# Save
	#idx=np.where(zO['Occurrence']==1)
	#gu.opickle(meta['Paths']['Project'] + '\\Inputs\\Ensembles\\wf_sim_osm_Scn' + cbu.FixFileNum(iScn) + '_Ens' + cbu.FixFileNum(iEns) + '.pkl',idx)

	print( np.round((time.time()-t0)/60,decimals=1) )

	return meta

#%% Run annual loop
meta=uosm.AnnualLoop(meta,z)

#%% Plot daily spread for one year

metaF['Par']['colors_list']=[(0,0,0),(0.8,0.8,0.8),(1,0,0),(0.2,0.2,0.2)]
metaF['Par']['cmap']=colors.ListedColormap(metaF['Par']['colors_list'])

plt.close('all')
fig=plt.figure(figsize=gu.cm2inch(12,12))
ax1=fig.add_axes([0,0,1,1])
ax1.set_axis_off()
im=ax1.matshow(metaF['Daily Status'][0,:,:].copy(),clim=(0,4),cmap=metaF['Par']['cmap'])
#ax1.set(position=[0,0,1,1],aspect='auto')
title=ax1.text(0.03,0.94,'Day ' + str(0),bbox={'facecolor':'none','alpha':0.75,'pad':5,'edgecolor':'none'},transform=ax1.transAxes,ha='left',fontsize=14,fontweight='bold')

def animate(i):
    im.set_array(metaF['Daily Status'][i,:,:])
    title.set_text('Day ' + str(i))
    return im,title

# call the animator. blit=True means only re-draw the parts that have changed.
anim=animation.FuncAnimation(fig,animate,frames=metaF['Par']['Season Length'],blit=True,interval=100,repeat=False)

# anim.save(meta['Paths']['Figures'] + '\\Wildfire_DailySpreadExample.mp4',fps=3,extra_args=['-vcodec','libx264'],dpi=200)


#%% Map mean occurrence

occ=[None]*metaF['N Scenario']
for iScn in range(metaF['N Scenario']):
    occ[iScn]=np.zeros( (metaF['N Time'],metaF['N Y'],metaF['N X']) ,dtype='int8')
    for iEns in range(metaF['N Ensemble']):
        idx=gu.ipickle(meta['Paths']['Project'] + '\\Inputs\\Ensembles\\wf_sim_osm_Scn' + cbu.FixFileNum(iScn) + '_Ens' + cbu.FixFileNum(iEns) + '.pkl')
        occ[iScn][idx]=occ[iScn][idx]+1
    occ[iScn]=occ[iScn]/metaF['N Ensemble']

iT1=np.where(metaF['Time']>2020)[0]

Mean1=np.mean(occ[0][iT1,:,:])*100
Mean2=np.mean(occ[2][iT1,:,:])*100

Sum=np.sum(occ[2][iT1,:,:])-np.sum(occ[0][iT1,:,:])
Mean=(np.mean(occ[2][iT1,:,:])-np.mean(occ[0][iT1,:,:]))*100

def Plot():
    
    plt.close('all')
    fig,ax=plt.subplots(1,3,figsize=gu.cm2inch(17,6))
    cm=plt.cm.get_cmap('viridis',20)
    
    # BaseCase    
    im=ax[0].matshow(np.mean(occ[0][iT1,:,:],axis=0)*100,clim=(0,1.2),extent=roi['grd']['Extent'],cmap=cm)
    ax[0].set(position=[0.01,0.12,0.31,0.87],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'],aspect='auto',visible='off')
    ax[0].axis('off')   
    cax=fig.add_axes([0.01,0.05,0.31,0.05])
    fig.colorbar(im,orientation='horizontal',cax=cax)
    
    # S1
    im=ax[1].matshow(np.mean(occ[2][iT1,:,:],axis=0)*100,clim=(0,1.2),extent=roi['grd']['Extent'],cmap=cm)
    ax[1].set(position=[0.335,0.12,0.31,0.87],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'],aspect='auto',visible='off')
    ax[1].axis('off')   
    cax=fig.add_axes([0.335,0.05,0.31,0.05])
    fig.colorbar(im,orientation='horizontal',cax=cax)
    
    # Difference
    d=(np.mean(occ[2][iT1,:,:],axis=0)-np.mean(occ[0][iT1,:,:],axis=0))*100
    im=ax[2].matshow(d,clim=(-0.7,0.7),extent=roi['grd']['Extent'],cmap=cm)
    ax[2].set(position=[0.66,0.12,0.31,0.87],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'],aspect='auto',visible='off')
    ax[2].axis('off')   
    cax=fig.add_axes([0.66,0.05,0.31,0.05])
    fig.colorbar(im,orientation='horizontal',cax=cax)

    gu.PrintFig(meta['Paths']['Figures'] + '\\Wildfire_CompMeanProb','png',900)


#%% Plot time series of summary stats

iScn=1
iEns=0

plt.close('all')
plt.bar(metaF['Time'],metaF['Scenarios'][iScn]['Area burned'][:,iEns]/metaF['Size']*100,1)

plt.close('all')
plt.plot(metaF['Time'],metaF['Scenarios'][iScn]['Area suseptable'][:,iEns]/metaF['Size']*100,'-ko')



#%% Animation yearly

iScn=2
iEns=0

occ=[None]*metaF['N Scenario']
for iScn in range(metaF['N Scenario']):
    occ[iScn]=np.zeros( (metaF['N Time'],metaF['N Y'],metaF['N X']) ,dtype='int8')
    idx=gu.ipickle(meta['Paths']['Project'] + '\\Inputs\\Ensembles\\wf_sim_osm_Scn' + cbu.FixFileNum(iScn) + '_Ens' + cbu.FixFileNum(iEns) + '.pkl')
    occ[iScn][idx]=occ[iScn][idx]+1

plt.close('all')
fig=plt.figure(figsize=gu.cm2inch(12,9))
ax1=fig.add_axes([0.01,0.01,0.6,0.82])
ax1.set_axis_off()
im=ax1.matshow(occ[iScn][0,:,:].copy(),clim=(0,2))
ax1.set(position=[0.01,0.01,0.7,0.96],aspect='auto')
title=ax1.text(0.12,0.92,'Year ' + str(0),
              bbox={'facecolor':'w','alpha':0.75,'pad':5},
              transform=ax1.transAxes,ha="center",fontsize=14,fontweight='bold')

def animate(i):
    im.set_array(occ[iScn][i,:,:])
    title.set_text('Year ' + str(i))
    return im,title

# call the animator. blit=True means only re-draw the parts that have changed.
anim=animation.FuncAnimation(fig,animate,frames=120,blit=True,interval=250,repeat=True)

#%% Run over a range of parameters to define relaitonship wiht annual prob occ

binO=np.linspace(0.0000001,0.0000004,10)
binS=np.linspace(0.02,0.14,10)
Pa=np.zeros((binO.size*binS.size,3))
cnt=0
for iO in range(binO.size):
    for iS in range(binS.size):
        metaF['p_Occurrence']=binO[iO]
        metaF['p_Spread_Forest']=binS[iS]
        ind=np.where(metaF['Grid_LCC']==1); 
        metaF['Grid_p_Spread'][ind]=np.random.normal(metaF['p_Spread_Forest'],0.3,size=ind[0].size)
        z={}
        z['Data']=metaF['Grid_LCC']
        z['Aburn']=np.zeros(metaF['Time'].size)
        z['Asus']=np.zeros(metaF['Time'].size)
        DataAll=np.zeros((metaF['Time'].size,z['Data'].shape[0],z['Data'].shape[1]),dtype='int8')
        metaF,z,DataAll=osm.OSM_Annual(metaF,z,DataAll)

        #shape,loc,scale=stats.pareto.fit(z['Aburn']/metaF['Size'])
        #b1=np.array([shape,loc,scale])
        #print(b1)
        Pa[cnt,0]=binO[iO]
        Pa[cnt,1]=binS[iS]
        Pa[cnt,2]=DataAll[DataAll==2].size/DataAll.size*100
        cnt=cnt+1
        print(cnt/(binO.size*binS.size))
        #print(Pa)

# Save
gu.opickle(r'C:\Users\rhember\Documents\Data\DisturbanceStatsByBGCZ\Wildfire_FireballVsAnnProb.pkl',Pa)



#%% Get Pareto parameters

shape,loc,scale=stats.pareto.fit(z['Aburn']/metaF['Size'])
b1=np.array([shape,loc,scale])
print(b1)

Pa=DataAll[DataAll==2].size/DataAll.size*100
print(Pa)



