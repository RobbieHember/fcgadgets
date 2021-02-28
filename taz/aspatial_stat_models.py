'''
Aspatial Statistical Models of Disturbance Events

'''

#%% Import modules

import numpy as np
import scipy.stats as stats
from fcgadgets.utilities import utilities_general as gu
from fcgadgets.cbrunner import cbrun_utilities as cbu

#%% Generate disturbances from Pareto distribution

def GenerateDisturbancesFromPareto(N_t,N_s,beta):

    # Initialize occurrence array
    oc=np.zeros((N_t,N_s),dtype='int8')
    
    # Draw a probability of area disturbed per time step
    po=stats.pareto.rvs(beta[0],loc=beta[1],scale=beta[2],size=N_t)
    po=np.reshape(po,(-1,1))
    po=np.tile(po,N_s)
    
    # Loop through time steps
    rn=np.random.random((N_t,N_s))
    
    # Populate occurrence
    ind=np.where(rn<po)    
    oc[ind[0],ind[1]]=1
    
    return oc

#%% Simulate probability of stand breakup based on age

def PredictStandBreakup_OnTheFly(meta,vi,iT,Age):

    beta=[-0.05,300]
    
    # Take a look at function:
    #Age=np.arange(1,500)
    #Po=1/(1+np.exp(beta[0]*(Age-beta[1])))
    #plt.close('all')
    #plt.plot(Age,Po*100,'b-')        
    
    Po=1/(1+np.exp(beta[0]*(Age-beta[1])))

    rn=np.random.random(Age.size)
    
    indS=np.where(rn<Po)[0]
    if indS.size>0:
        for i in range(indS.size):
            iAvailable=np.where(vi['EC']['ID_Type'][iT,indS[i],:]==0)[0]
            if iAvailable.size>0:
                iE=iAvailable[0]+0
                #iE=-1
                vi['EC']['ID_Type'][iT,indS[i],iE]=meta['LUT']['Dist']['Mechanical']
                vi['EC']['MortalityFactor'][iT,indS[i],iE]=1
                vi['EC']['ID_GrowthCurve'][iT,indS[i],iE]=1
    
    return vi

#%% Simulate probability of harvesting on the fly

def PredictHarvesting_OnTheFly(meta,vi,iT,V_Merch,Period):
    
    # Indicator of THLB (THLB=1, Non-THLB=0)
    flag_thlb=vi['Inv']['THLB'][iT,:]
    
    # Deterministic component    
    if Period=='Historical':
        #f1=0.0014*25**((meta['Year'][iT]-1900)/100)
        #f2=(1/(1+np.exp(0.12*(Year-1950))))
        f1=0.0011*35**((meta['Year'][iT]-1900)/100)
        f2=(1/(1+np.exp(0.3*(meta['Year'][iT]-1960))))        
        beta0=f1*f2
    else:
        # Future
        beta0=0.03
    
    beta=[beta0,-0.025,400]
    #V_Merch=np.arange(1,1200)
    #Po=beta[0]*(1/(1+np.exp(beta[1]*(V_Merch-beta[2]))))
    #plt.close('all')
    #plt.plot(V_Merch,Po*100,'b-')
        
    Po=beta[0]*(1/(1+np.exp(beta[1]*(V_Merch-beta[2]))))
    rn=np.random.random(V_Merch.size)
    
    Oc=flag_thlb*np.floor(np.minimum(1,Po/rn))
    
    indS=np.where(Oc==1)[0]    
    
    if indS.size>0:
        for i in range(indS.size):
            iAvailable=np.where(vi['EC']['ID_Type'][iT,indS[i],:]==0)[0]        
            if iAvailable.size>0:
                iE=iAvailable[0]+0
                vi['EC']['ID_Type'][iT,indS[i],iE]=meta['LUT']['Dist']['Harvest']
                vi['EC']['MortalityFactor'][iT,indS[i],iE]=1
                vi['EC']['ID_GrowthCurve'][iT,indS[i],iE]=1
                iE=iAvailable[0]+1 # changing this to zero will cause the harvest to be overwritten
                vi['EC']['ID_Type'][iT,indS[i],iE]=meta['LUT']['Dist']['Slashpile Burn']
                vi['EC']['MortalityFactor'][iT,indS[i],iE]=1
                vi['EC']['ID_GrowthCurve'][iT,indS[i],iE]=1
                iE=iAvailable[0]+2
                vi['EC']['ID_Type'][iT,indS[i],iE]=meta['LUT']['Dist']['Planting']
                vi['EC']['MortalityFactor'][iT,indS[i],iE]=1
                vi['EC']['ID_GrowthCurve'][iT,indS[i],iE]=1
    
    return vi

#%% Generate disturbance ensembles with AAO models

def GenerateIBMEnsembleFromAAO(meta,par,id_bgcz):
    
    # Import IBM stats    
    ibmss=gu.ipickle(meta['Paths']['Taz Datasets'] + '\\Beetle Stats and Scenarios\\IBM_Stats_Scenarios_By_BGCZ.pkl')
    tv_scn=np.arange(-2000,2201,1)
    
    # Prepare mortality probability coefficients
    beta_obs=np.cumsum([par['IBM']['p_Trace_obs'],
                       par['IBM']['p_Low_obs'],
                       par['IBM']['p_Medium_obs'],
                       par['IBM']['p_Severe_obs'],
                       par['IBM']['p_VerySevere_obs']])
    
    ibm_sim={}
    
    #--------------------------------------------------------------------------
    # Occurrence
    #--------------------------------------------------------------------------

    # Initialize annual probability of occurrence 
    ibm_sim['Occurrence']=np.zeros((meta['Year'].size,meta['N Stand Full']),dtype='int16')

    uZone=np.unique(id_bgcz)
    
    for iZone in range(uZone.size):
        
        indZone=np.where(id_bgcz==uZone[iZone])[0]
            
        namZone=cbu.lut_n2s(meta['LUT']['VRI']['BEC_ZONE_CODE'],uZone[iZone])[0]
            
        # Alternative model
        b0=ibmss[namZone]['Beta_Pareto_Alt'].copy()
        for iT in range(meta['Year'].size):
            ibm_sim['Occurrence'][iT,indZone]=GenerateDisturbancesFromPareto(1,indZone.size,b0)
        
    # Exclude inventory period
    if par['IBM']['Exclude simulations during modern era']=='On':
        ind=np.where( (meta['Year']>=1951) & (meta['Year']<=meta['Year Project']) )[0]
        ibm_sim['Occurrence'][ind,:]=0    
        
    #--------------------------------------------------------------------------
    # Severity / mortality
    #--------------------------------------------------------------------------
        
    # Get mortality from probability of burn severity rating
    ibm_sim['Mortality']=np.zeros((meta['Year'].size,meta['N Stand Full']),dtype='int16')
    ind=np.where( (ibm_sim['Occurrence']>0) )
    
    # Add pre-inventory severity
    Mort=np.zeros((meta['Year'].size,meta['N Stand Full']),dtype='int16')    
    Mort[ind[0],ind[1]]=GetMortalityFromIBMSeverity(ind[0].size,beta_obs)
    it=np.where(meta['Year']<1920)[0]
    ibm_sim['Mortality'][it,:]=Mort[it,:]
    
    # Add post-inventory severity
    Mort=np.zeros((meta['Year'].size,meta['N Stand Full']),dtype='int16')    
    Mort[ind[0],ind[1]]=GetMortalityFromIBMSeverity(ind[0].size,beta_obs)
    it=np.where(meta['Year']>meta['Year Project'])[0]
    ibm_sim['Mortality'][it,:]=Mort[it,:]
    
    return ibm_sim

#%% Mortality from burn severity rating
# Mortality numbers come from DisturbanceBySeverityClass spreadsheet

def GetMortalityFromIBMSeverity(n,beta):
    
    rn=np.random.random(n)
    y=np.zeros(rn.size)
    
    for i in range(n):
        
        if rn[i]<beta[0]:
            # Trace
            y[i]=1
        elif (rn[i]>=beta[0]) & (rn[i]<beta[1]):
            # Low
            y[i]=5
        elif (rn[i]>=beta[1]) & (rn[i]<beta[2]): 
            # Medium
            y[i]=25
        elif (rn[i]>=beta[2]) & (rn[i]<beta[3]): 
            # Severe
            y[i]=50
        elif (rn[i]>=beta[3]): 
            # Severe
            y[i]=75    
        
    return y

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
        
        namZone=cbu.lut_n2s(meta['LUT']['VRI']['BEC_ZONE_CODE'],uZone[iZone])[0]
        
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
                    wf_sim['Occurrence'][iT,:]=GenerateDisturbancesFromPareto(1,indZone.size,b0)
                else:
                    # All stands get populated with the same prediction
                    wf_sim['Occurrence'][iT,indZone]=GenerateDisturbancesFromPareto(1,indZone.size,b0)        
    
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
