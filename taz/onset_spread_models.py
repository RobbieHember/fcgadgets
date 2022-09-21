'''

ONSET AND SPREAD MODELS

'''

#%% Import modules

import numpy as np
import scipy.stats as stats
from matplotlib import colors
from matplotlib import animation
import time
import copy
from fcgadgets.macgyver import utilities_general as gu
from fcgadgets.cbrunner import cbrun_utilities as cbu
from fcgadgets.taz import aspatial_stat_models as gensm

#%% Spread modelling tools

def ConfigureProject(meta,lsc):
    
    metaF={}
    metaF['N Time'],metaF['N Y'],metaF['N X']=lsc['Scenarios'][0]['Use'].shape
    metaF['Shape']=(metaF['N Y'],metaF['N X'])
    metaF['Size']=metaF['N Y']*metaF['N X']
    metaF['N Scenario']=lsc['N Scenario']
    metaF['Time']=lsc['tv']
    metaF['N Ensemble']=meta['Project']['N Ensemble']
    
    if 'Save Daily Status' not in metaF:
        metaF['Save Daily Status']='Off'    
    
    #--------------------------------------------------------------------------
    # Parameters
    #--------------------------------------------------------------------------
    
    metaF['Par']={}
    
    # Displacements from a cell to its eight nearest neighbours
    metaF['Par']['neighbourhood']=((-1,1),(0,1),(1,1),(1,0),(1, -1),(0,-1),(-1,-1),(-1,0))
    #direction=['NW','N','NE','E','SE','S','SW','W']
    metaF['Par']['UNSUSEPTABLE']=0
    metaF['Par']['SUSEPTABLE']=1    
    metaF['Par']['ONSET']=2
    metaF['Par']['JUST BURNED']=3
    metaF['Par']['non']=lambda s: s if s<0 else None
    metaF['Par']['mom']=lambda s: max(0,s)
    
    metaF['Par']['Season Length']=101
    metaF['Par']['Duration Limit']=3
    
    metaF['Par']['Prob Continue']=0.0
    
    metaF['Par']['Prob Onset']={}
    metaF['Par']['Prob Onset']['Timber Production']=0.000000075 # 0.0000005
    metaF['Par']['Prob Onset']['Energy Production']=0.5*metaF['Par']['Prob Onset']['Timber Production']
    metaF['Par']['Prob Onset']['Conservation Natural']=1.0*metaF['Par']['Prob Onset']['Timber Production']
    metaF['Par']['Prob Onset']['Conservation Consistent']=1.0*metaF['Par']['Prob Onset']['Timber Production']
    metaF['Par']['Prob Onset']['Fuel Break']=0.5*metaF['Par']['Prob Onset']['Timber Production']
    metaF['Par']['Prob Onset']['Peat']=0.25*metaF['Par']['Prob Onset']['Timber Production']
    metaF['Par']['Prob Onset']['Crop']=0.5*metaF['Par']['Prob Onset']['Timber Production']
    
    #metaF['Par']['p_Regen']=0.25
    metaF['Par']['Prob Spread']={}
    metaF['Par']['Prob Spread']['Timber Production']=0.18
    metaF['Par']['Prob Spread']['Energy Production']=0.5*metaF['Par']['Prob Spread']['Timber Production']
    metaF['Par']['Prob Spread']['Conservation Natural']=1.0*metaF['Par']['Prob Spread']['Timber Production']
    metaF['Par']['Prob Spread']['Conservation Consistent']=1.0*metaF['Par']['Prob Spread']['Timber Production']
    metaF['Par']['Prob Spread']['Fuel Break']=0.5*metaF['Par']['Prob Spread']['Timber Production']
    metaF['Par']['Prob Spread']['Peat']=0.25*metaF['Par']['Prob Spread']['Timber Production']
    metaF['Par']['Prob Spread']['Crop']=0.5*metaF['Par']['Prob Spread']['Timber Production']
    metaF['Par']['Years To Regenerate']=15
    
    # Colours for visualization: brown for EMPTY, dark green for TREE and orange
    # for FIRE. Note that for the colormap to work, this list and the bounds list
    # must be one larger than the number of different values in the array.
    metaF['Par']['colors_list']=[(0,0,0),(0.3,0.3,0.3),(1,0,0),(1,0.5,0)]
    metaF['Par']['cmap']=colors.ListedColormap(metaF['Par']['colors_list'])
    metaF['Par']['bounds']=[0,1,2,3]
    metaF['Par']['norm']=colors.BoundaryNorm(metaF['Par']['bounds'],metaF['Par']['cmap'].N)
    
    #--------------------------------------------------------------------------
    # Grids
    #--------------------------------------------------------------------------
    
    metaF['Scenarios']=[None]*lsc['N Scenario']
    for iScn in range(lsc['N Scenario']):
        
        metaF['Scenarios'][iScn]={}
        
        metaF['Scenarios'][iScn]['Name']=lsc['Scenarios'][iScn]['Name']
             
        metaF['Scenarios'][iScn]['Grid']={}
        
        # Initialize probability of onset
        metaF['Scenarios'][iScn]['Grid']['Prob Onset']=np.zeros(lsc['Scenarios'][iScn]['Cover'].shape,dtype=float)
        
        ind=np.where(lsc['Scenarios'][iScn]['Use']==meta['LUT']['LSC']['Use']['Timber Production'])
        metaF['Scenarios'][iScn]['Grid']['Prob Onset'][ind]=metaF['Par']['Prob Onset']['Timber Production']
        
        ind=np.where(lsc['Scenarios'][iScn]['Use']==meta['LUT']['LSC']['Use']['Energy Production'])
        metaF['Scenarios'][iScn]['Grid']['Prob Onset'][ind]=metaF['Par']['Prob Onset']['Energy Production']
        
        ind=np.where(lsc['Scenarios'][iScn]['Use']==meta['LUT']['LSC']['Use']['Conservation Natural'])
        metaF['Scenarios'][iScn]['Grid']['Prob Onset'][ind]=metaF['Par']['Prob Onset']['Conservation Natural']
        
        ind=np.where(lsc['Scenarios'][iScn]['Use']==meta['LUT']['LSC']['Use']['Conservation Consistent'])
        metaF['Scenarios'][iScn]['Grid']['Prob Onset'][ind]=metaF['Par']['Prob Onset']['Conservation Consistent']
        
        ind=np.where(lsc['Scenarios'][iScn]['Cover']==meta['LUT']['LSC']['Cover']['Peat'])
        metaF['Scenarios'][iScn]['Grid']['Prob Onset'][ind]=metaF['Par']['Prob Onset']['Peat']
        
        ind=np.where(lsc['Scenarios'][iScn]['Cover']==meta['LUT']['LSC']['Cover']['Crop'])
        metaF['Scenarios'][iScn]['Grid']['Prob Onset'][ind]=metaF['Par']['Prob Onset']['Crop']
        
        # Initialize probability of spread
        metaF['Scenarios'][iScn]['Grid']['Prob Spread']=np.zeros(lsc['Scenarios'][iScn]['Cover'].shape,dtype=float)
        
        ind=np.where(lsc['Scenarios'][iScn]['Use']==meta['LUT']['LSC']['Use']['Timber Production'])
        metaF['Scenarios'][iScn]['Grid']['Prob Spread'][ind]=metaF['Par']['Prob Spread']['Timber Production']
        
        ind=np.where(lsc['Scenarios'][iScn]['Use']==meta['LUT']['LSC']['Use']['Energy Production'])
        metaF['Scenarios'][iScn]['Grid']['Prob Spread'][ind]=metaF['Par']['Prob Spread']['Energy Production']
        
        ind=np.where(lsc['Scenarios'][iScn]['Use']==meta['LUT']['LSC']['Use']['Conservation Natural'])
        metaF['Scenarios'][iScn]['Grid']['Prob Spread'][ind]=metaF['Par']['Prob Spread']['Conservation Natural']
        
        ind=np.where(lsc['Scenarios'][iScn]['Use']==meta['LUT']['LSC']['Use']['Conservation Consistent'])
        metaF['Scenarios'][iScn]['Grid']['Prob Spread'][ind]=metaF['Par']['Prob Spread']['Conservation Consistent']
        
        ind=np.where(lsc['Scenarios'][iScn]['Cover']==meta['LUT']['LSC']['Cover']['Peat'])
        metaF['Scenarios'][iScn]['Grid']['Prob Spread'][ind]=metaF['Par']['Prob Spread']['Peat']
        
        ind=np.where(lsc['Scenarios'][iScn]['Cover']==meta['LUT']['LSC']['Cover']['Crop'])
        metaF['Scenarios'][iScn]['Grid']['Prob Spread'][ind]=metaF['Par']['Prob Spread']['Crop']
    
        # Initialize stats
        metaF['Scenarios'][iScn]['Area burned']=np.zeros((metaF['Time'].size,metaF['N Ensemble']),dtype=int)
        metaF['Scenarios'][iScn]['Area suseptable']=np.zeros((metaF['Time'].size,metaF['N Ensemble']),dtype=int)
    
    return metaF

#%% Onset-spread model annual loop
    
def Fireball(meta,metaF):
    
    t0=time.time()
    
    for iEns in range(metaF['N Ensemble']):
        
        for iScn in range(metaF['N Scenario']):    
    
            print('Running OS model for LSC scenario:' + metaF['Scenarios'][iScn]['Name'] + ', Ensemble:' + str(iEns+1) )
            
            z={}
            
            # Initialize occurrence
            z['Occurrence']=np.zeros(metaF['Scenarios'][iScn]['Grid']['Prob Spread'].shape,dtype=int)
            
            # Initialize status
            z['Status']=metaF['Par']['SUSEPTABLE']*np.ones((metaF['N Y'],metaF['N X']),dtype='int8')
            
            # Initialize years since fire
            z['Years Since Occurrence']=(metaF['Par']['Years To Regenerate']+1)*np.ones((metaF['N Y'],metaF['N X']),dtype='int8')
            
            # Annual loop
            for iYear in range(metaF['Time'].size):    
                
                # Run daily loop
                metaF,z=OSM_DailyLoop(metaF,z,iScn,iYear)
                
                # Record occurrence
                ind=(z['Duration']>0)
                z['Occurrence'][iYear,ind]=1
                
                # Update the years since occurrence counter        
                ind=(z['Years Since Occurrence']<metaF['Par']['Years To Regenerate'])
                z['Years Since Occurrence'][ind]=z['Years Since Occurrence'][ind]+1
            
                # Transfer affected stands back to the suseptable class upon regeneration
                ind=(z['Years Since Occurrence']==metaF['Par']['Years To Regenerate'])
                z['Status'][ind]=metaF['Par']['SUSEPTABLE']
                
                # Adjust counter to be just above threshold regeneration period so that
                # they will not be found in queries
                z['Years Since Occurrence'][ind]=metaF['Par']['Years To Regenerate']+1
                
                # Update stats
                metaF['Scenarios'][iScn]['Area burned'][iYear,iEns]=np.sum(z['Occurrence'][iYear,:,:])
                metaF['Scenarios'][iScn]['Area suseptable'][iYear,iEns]=np.sum( (z['Status']==metaF['Par']['SUSEPTABLE']) )
    
                #print(metaF['Time'][iYear],z['Aburn'][iYear])
                
                if metaF['Save Daily Status']=='On':
                    return metaF
                    
            # Save
            idx=np.where(z['Occurrence']==1)
            gu.opickle(meta['Paths']['Project'] + '\\Inputs\\Ensembles\\wf_sim_osm_Scn' + cbu.FixFileNum(iScn) + '_Ens' + cbu.FixFileNum(iEns) + '.pkl',idx)
    
    print( np.round((time.time()-t0)/60,decimals=1) )
    
    return metaF

#%% Onset-spread model daily loop

def OSM_DailyLoop(metaF,z,iScn,iYear):
    
    # Probability of onset and spread in current year
    p_Onset=metaF['Scenarios'][iScn]['Grid']['Prob Onset'][iYear,:,:].copy()
    p_Spread=metaF['Scenarios'][iScn]['Grid']['Prob Spread'][iYear,:,:].copy()
    
    if metaF['Save Daily Status']=='On':
        metaF['Daily Status']=np.zeros((metaF['Par']['Season Length'],metaF['N Y'],metaF['N X']),dtype='int8')
    
    # At the beginning of the year, convert ONSET from the end of the previous year to UNSUSEPTABLE        
    if metaF['Par']['Prob Continue']>0:
        
        # Event can persist across years
        rn=np.random.random( (metaF['N Y'],metaF['N X']) )
        ind=(z['Status']==metaF['Par']['ONSET']) & (rn>=metaF['Par']['Prob Continue'])
        z['Status'][ind]=metaF['Par']['UNSUSEPTABLE']
        
    else:
        
        # No carryover across years
        ind=(z['Status']==metaF['Par']['ONSET'])
        z['Status'][ind]=metaF['Par']['UNSUSEPTABLE']
        
    # Unpack parameters
    mom=metaF['Par']['mom']
    non=metaF['Par']['non']
    
    # Initialize duration of active events
    z['Duration']=np.zeros( (metaF['N Y'],metaF['N X']) ,dtype='int8' )
    
    for iDay in range(metaF['Par']['Season Length']):
        
        # Onset
        rn_onset=np.random.random( (metaF['N Y'],metaF['N X']) )
        ind=(z['Status']==metaF['Par']['SUSEPTABLE']) & (rn_onset<p_Onset)
        z['Status'][ind]=metaF['Par']['ONSET']
        z['Years Since Occurrence'][ind]=0
        
        # Spread       
        rn_spread=np.random.random( (metaF['N Y'],metaF['N X']) )
        daily_weather_factor=np.random.normal(loc=0,scale=0.1,size=1)
        for dx,dy in metaF['Par']['neighbourhood']:
            z_shift=z['Status'].copy()
            z_shift[mom(dy):non(dy),mom(dx):non(dx)]=z['Status'][mom(-dy):non(-dy),mom(-dx):non(-dx)]
            
            ind=(z['Status']==metaF['Par']['SUSEPTABLE']) & (z_shift==metaF['Par']['ONSET']) & (rn_spread<(p_Spread+daily_weather_factor))
            z['Status'][ind]=metaF['Par']['ONSET']
            z['Years Since Occurrence'][ind]=0
        
        # Track duration of active onset
        ind=(z['Status']==metaF['Par']['ONSET'])
        z['Duration'][ind]=z['Duration'][ind]+1
        
        # Convert ONSET to UNSUSCEPTABLE (or JUST BURNED FOR GRAPHIC PURPOSES)
        ind=(z['Duration']>=metaF['Par']['Duration Limit'])
        
        # Save daily data
        if metaF['Save Daily Status']=='On':
            metaF['Daily Status'][iDay,:,:]=z['Status'].copy()
            metaF['Daily Status'][iDay,ind]=metaF['Par']['JUST BURNED']
        
        # Convert ONSET to UNSUSCEPTABLE
        z['Status'][ind]=metaF['Par']['UNSUSEPTABLE']
        
    return metaF,z
