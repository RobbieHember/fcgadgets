'''

WILDFIRE STATS MODELS

See documentation.

'''

import numpy as np
import scipy.stats as stats
from fcgadgets.utilities import utilities_general as gu
from fcgadgets.cbrunner import cbrun_utilities as cbu
from fcgadgets.taz import general_stat_models as gensm

#%% Generate disturbance ensembles with AAO models

def GenerateWildfireEnsembleFromAAO(meta,par,id_bgcz,method_occ):
    
    # Import wildfire stats (by BGC zone)  
    wfss=gu.ipickle(meta['Paths']['Taz Datasets'] + '\\Wildfire Stats and Scenarios\\Wildfire_Stats_Scenarios_By_BGCZ.pkl')    
    tv_wfss=np.arange(-2000,2201,1)
    
    # Prepare mortality probability coefficients
    beta_pi=np.cumsum([par['WF']['p_Unburned_pi'],par['WF']['p_Low_pi'],par['WF']['p_Medium_pi'],par['WF']['p_High_pi']])
    beta_obs=np.cumsum([par['WF']['p_Unburned_obs'],par['WF']['p_Low_obs'],par['WF']['p_Medium_obs'],par['WF']['p_High_obs']])
    
    wf_sim={}

    #--------------------------------------------------------------------------
    # Occurrence
    #--------------------------------------------------------------------------

    # Initialize annual probability of occurrence (final with deterministic and
    # random components)
    wf_sim['Occurrence']=np.zeros((meta['Year'].size,meta['N Stand Full']),dtype='int8')
    
    # Get unique BGC zones
    uZone=np.unique(id_bgcz)
    
    for iZone in range(uZone.size):
        
        namZone=cbu.lut_n2s(meta['LUT BGC Zone'],uZone[iZone])[0]
        
        indZone=np.where(id_bgcz==uZone[iZone])[0]
        
        if method_occ=='DirectFromParetoDraw':
            
            Po_Det=wfss[namZone]['Po_Det_WF_Scn' + str(int(par['WF']['Scenario ID']))]
        
            for iT in range(meta['Year'].size):
            
                # Adjust shape parameter to match specified annual probability of 
                # occurrence from the deterministic component
                ind_scn=np.where(tv_wfss==meta['Year'][iT])[0]
                #b0=wfss[namZone]['Beta_Pareto'].copy()
                #b_shape=wfss[namZone]['Pareto_shape_for_Po'].copy()
                #b0[0]=np.exp(b_shape[1]*np.log(Po_Det[ind_scn])+b_shape[0])
                b0=wfss[namZone]['Beta_Pareto_Cal'].copy()
                Scale=wfss[namZone]['Pareto_scale_to_match_Po_mu'][1]*Po_Det[ind_scn]+wfss[namZone]['Pareto_scale_to_match_Po_mu'][0]
                b0[1]=-Scale
                b0[2]=Scale
            
                if meta['Scenario Source']=='Spreadsheet':
                    # When run from spreadsheet, stands are swapped for ensembles so
                    # generate different records for each stand
                    wf_sim['Occurrence'][iT,:]=gensm.GenerateDisturbancesFromPareto(1,indZone.size,b0)
                else:
                    # All stands get populated with the same prediction
                    wf_sim['Occurrence'][iT,indZone]=gensm.GenerateDisturbancesFromPareto(1,indZone.size,b0)        
    
        elif method_occ=='PreRun':
            # Using this will ensure consistency across tiles
            
            Po=wfss[namZone]['PctOcc_DetPlusRand_WF_Scn' + str(int(par['WF']['Scenario ID']))]/100
            for iT in range(meta['Year'].size):
                rn=np.random.random(indZone.size)
                ind_wfss=np.where(tv_wfss==meta['Year'][iT])[0]
                indMort=np.where(rn<Po[ind_wfss])[0]
                wf_sim['Occurrence'][iT,indZone[indMort]]=1    
    
    # Exclude inventory period
    if par['WF']['Exclude simulations during modern era']=='On':
        ind=np.where( (meta['Year']>=1920) & (meta['Year']<=meta['Year Project']) )[0]
        wf_sim['Occurrence'][ind,:]=0
    
    #--------------------------------------------------------------------------
    # Get mortality from probability of burn severity rating
    #--------------------------------------------------------------------------

    wf_sim['Mortality']=np.zeros((meta['Year'].size,meta['N Stand Full']),dtype='int8')
    
    for iT in range(meta['Year'].size):
        ind=np.where( (wf_sim['Occurrence'][iT,:]>0) )[0]
        if meta['Year'][iT]<1920:
            Mort=GetMortalityFromBurnSeverityRating(ind.size,beta_pi)
        else:
            Mort=GetMortalityFromBurnSeverityRating(ind.size,beta_obs)
        wf_sim['Mortality'][iT,ind]=Mort
        
    return wf_sim

#%% Mortality from burn severity rating

def GetMortalityFromBurnSeverityRating(n,beta):
    
    rn=np.random.random(n)
    y=np.zeros(rn.size)
    
    for i in range(n):
        
        if rn[i]<beta[0]:
            # Unburned
            y[i]=5
        elif (rn[i]>=beta[0]) & (rn[i]<beta[1]):
            # Low
            y[i]=50
        elif (rn[i]>=beta[1]) & (rn[i]<beta[2]): 
            # Medium
            y[i]=90
        elif (rn[i]>=beta[2]): 
            # High
            y[i]=100
        
    return y

#%% Spread modelling tools

def Import_OSM_Tools():
    
    ost={}
    
    # Displacements from a cell to its eight nearest neighbours
    #neighbourhood=((-1,-1),(-1,0),(-1,1),(0,-1),(0, 1),(1,-1),(1,0),(1,1))
    ost['neighbourhood']=((-1,1),(0,1),(1,1),(1,0),(1, -1),(0,-1),(-1,-1),(-1,0))
    #direction=['NW','N','NE','E','SE','S','SW','W']
    ost['SUSEPTABLE']=1
    ost['UNSUSEPTABLE']=0
    ost['ONSET']=2

    ost['non']=lambda s: s if s<0 else None
    ost['mom']=lambda s: max(0,s)
    
    return ost

#%% Onset-spread model annual loop
    
def OSM_Annual(metaOS,z,DataAll):
    
    for iYear in range(metaOS['Time'].size):    
    
        # Run daily loop
        z=OSM_DailyLoop(z,metaOS,iYear)
    
        # Regenerate previously burnt stands
        ind=(metaOS['Grid_YearsSinceFire']<metaOS['YearsToRegenerate'])
        metaOS['Grid_YearsSinceFire'][ind]=metaOS['Grid_YearsSinceFire'][ind]+1
    
        ind=(metaOS['Grid_YearsSinceFire']==metaOS['YearsToRegenerate'])
        metaOS['Grid_YearsSinceFire'][ind]=metaOS['YearsToRegenerate']+1
        z['Data'][ind]=metaOS['ost']['SUSEPTABLE']
        
        z['Asus'][iYear]=np.sum(z['Data']==1)
        
        DataAll[iYear,:,:]=z['Data']    
        #print(metaOS['Time'][iYear],z['Aburn'][iYear])
    
    return metaOS,z,DataAll

#%% Onset-spread model daily loop

def OSM_DailyLoop(z,metaOS,iYear):
    
    if iYear==0:
        # Initialize conditions at start of year
        z['Data']=metaOS['Grid_LCC'].copy()
    else:
        # Probability of regeneration
        #rn=np.random.random(metaOS['Shape'])
        #z['Data'][(metaOS['Grid_LCC']==SUSEPTABLE) & (z['Data']==UNSUSEPTABLE) & (rn<metaOS['p_Regen'])]=SUSEPTABLE        
    
        # Transition ONSET to UNSUSEPTABLE
        z['Data'][(z['Data']==metaOS['ost']['ONSET'])]=metaOS['ost']['UNSUSEPTABLE']
    
    neighbourhood=metaOS['ost']['neighbourhood']
    mom=metaOS['ost']['mom']
    non=metaOS['ost']['non']
    ONSET=metaOS['ost']['ONSET']
    SUSEPTABLE=metaOS['ost']['SUSEPTABLE']
    #UNSUSEPTABLE=metaOS['ost']['UNSUSEPTABLE']    
    
    for iDay in range(metaOS['Season_length']):
                
        # Spread       
        rn=np.random.random(metaOS['Shape'])
        #wind=[0.9,1.1,1.25,1.15,1.15,0.75,0.8,0.9]
        #wind=[0.25,0.25,0.25,0.25,1,2,2,2]
        wind=np.ones(8)
        cnt_wind=0
        #np.sum(wind)
        for dx,dy in neighbourhood:
            rn1=rn*wind[cnt_wind]
            z_shift=z['Data'].copy()
            z_shift[mom(dy):non(dy),mom(dx):non(dx)]=z['Data'][mom(-dy):non(-dy),mom(-dx):non(-dx)]
            ind=(z['Data']==SUSEPTABLE) & (z_shift==ONSET) & (rn1<metaOS['Grid_p_Spread'])
            z['Data'][ind]=ONSET
            metaOS['Grid_YearsSinceFire'][ind]=0
            cnt_wind=cnt_wind+1
            
        # Onset
        rn=np.random.random(metaOS['Shape'])
        ind=(metaOS['Grid_LCC']==1) & (rn<metaOS['p_Occurrence'])
        z['Data'][ind]=ONSET
        metaOS['Grid_YearsSinceFire'][ind]=0
            
        # Summarize results    
        z['Aburn'][iYear]=z['Aburn'][iYear]+np.sum(z['Data']==2)
    
    return z


