'''

PROBABILISTIC MODELS OF BREAKUP

See documentation.

'''

#%% Import modules

import numpy as np
import scipy.stats as stats
from fcgadgets.utilities import utilities_general as gu
from fcgadgets.cbrunner import cbrun_utilities as cbu

#%% Generate disturbances from Pareto distribution

def GenerateDisturbancesFromPareto(m,n,beta):
    po=stats.pareto.rvs(beta[0],loc=beta[1],scale=beta[2],size=m)
    po=np.reshape(po,(-1,1))
    po=np.tile(po,n)
    r=np.random.random((m,n))
    ind=np.where(r<po)
    y=np.zeros((m,n),dtype='int16')
    y[ind[0],ind[1]]=1    
    return y

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
                vi['EC']['ID_Type'][iT,indS[i],iE]=meta['LUT Dist']['Mechanical']
                vi['EC']['MortalityFactor'][iT,indS[i],iE]=1
                vi['EC']['ID_GrowthCurve'][iT,indS[i],iE]=1
    
    return vi

#%% Simulate probability of harvesting on the fly

def PredictHarvesting_OnTheFly(meta,vi,iT,flag_thlb,V_Merch):
    
    # Deterministic component    
    if meta['Year'][iT]<2020:
        #f1=0.0014*25**((meta['Year'][iT]-1900)/100)
        #f2=(1/(1+np.exp(0.12*(Year-1950))))
        f1=0.0011*35**((meta['Year'][iT]-1900)/100)
        f2=(1/(1+np.exp(0.3*(meta['Year'][iT]-1960))))        
        beta0=f1*f2
    else:
        beta0=0.02
    
    beta=[beta0,-0.025,550]
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
                vi['EC']['ID_Type'][iT,indS[i],iE]=meta['LUT Dist']['Harvest']
                vi['EC']['MortalityFactor'][iT,indS[i],iE]=1
                vi['EC']['ID_GrowthCurve'][iT,indS[i],iE]=1
                iE=iAvailable[0]+1
                vi['EC']['ID_Type'][iT,indS[i],iE]=meta['LUT Dist']['Slashpile Burn']
                vi['EC']['MortalityFactor'][iT,indS[i],iE]=1
                vi['EC']['ID_GrowthCurve'][iT,indS[i],iE]=1
                iE=iAvailable[0]+2
                vi['EC']['ID_Type'][iT,indS[i],iE]=meta['LUT Dist']['Planting']
                vi['EC']['MortalityFactor'][iT,indS[i],iE]=1
                vi['EC']['ID_GrowthCurve'][iT,indS[i],iE]=1
    
    return vi

#%% Generate disturbance ensembles with AAO models

def GenerateIBMEnsembleFromAAO(meta,par,id_bgcz):
    
    # Import IBM stats    
    ibmss=gu.ipickle(meta['Paths']['Beetle Stats and Scenarios File'])
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
            
        namZone=cbu.lut_n2s(meta['LUT VRI']['BEC_ZONE_CODE'],uZone[iZone])[0]
            
        # Alternative model
        b0=ibmss[namZone]['Beta_Pareto_Alt'].copy()
        ibm_sim['Occurrence'][:,indZone]=GenerateDisturbancesFromPareto(meta['Year'].size,indZone.size,b0)
        
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

def GenerateIBMEnsembleFromAAO_old(meta,par,id_bgcz):
    
    # Import wildfire stats    
    ibmss=gu.ipickle(meta['Paths']['Beetle Stats and Scenarios File'])
    tv_scn=np.arange(-2000,2201,1)
    
    # Prepare mortality probability coefficients
    beta_obs=np.cumsum([par['IBM']['p_Trace_obs'],
                       par['IBM']['p_Low_obs'],
                       par['IBM']['p_Medium_obs'],
                       par['IBM']['p_Severe_obs'],
                       par['IBM']['p_VerySevere_obs']])
    
    ibm_sim=[]
    for iEns in range(meta['N Ensemble']):
    
        # Initialize annual probability of occurrence (final with deterministic and
        # random components)
        Oc_sim=np.zeros((meta['Year'].size,meta['N Stand Full']))
    
        uZone=np.unique(id_bgcz)
    
        for iZone in range(uZone.size):
        
            indZone=np.where(id_bgcz==uZone[iZone])[0]
            
            namZone=cbu.lut_n2s(meta['LUT VRI']['BEC_ZONE_CODE'],uZone[iZone])[0]
            
            # Alternative model
            b0=ibmss[namZone]['Beta_Pareto_Alt'].copy()
            Oc_sim[:,indZone]=GenerateDisturbancesFromPareto(meta['Year'].size,indZone.size,b0)
            
            flg=0
            if flg==1:
                Po_Det=ibmss[namZone]['Po_Det_Scn' + str(par['IBM']['Scenario ID'].astype(int))]        
                for iT in range(meta['Year'].size):
            
                    ind_scn=np.where(tv_scn==meta['Year'][iT])[0]
            
                    # Adjust shape parameter to match specified annual probability of 
                    # occurrence from the deterministic component
                    b0=ibmss[namZone]['Beta_Pareto'].copy()
                    b_shape=ibmss[namZone]['Pareto_shape_for_Po'].copy()
                    b0[0]=np.exp(b_shape[1]*np.log(Po_Det[ind_scn])+b_shape[0])
            
                    Oc_sim[iT,indZone]=GenerateDisturbancesFromPareto(1,indZone.size,b0)
    
        # Exclude inventory period
        if par['IBM']['Exclude simulations during modern era']=='On':
            ind=np.where( (meta['Year']>=1951) & (meta['Year']<=2020) )[0]
            Oc_sim[ind,:]=0
    
        # Put events in a stand list, to be nexted within ensemble list
        # Some may be empty!
        ibm_sim0=[None]*meta['N Stand Full']
        for iStand in range(meta['N Stand Full']):
            
            ind=np.where(Oc_sim[:,iStand]==1)[0]
            if ind.size==0:
                continue
            
            d={}
            d['Year']=meta['Year'][ind]
            
            # Get mortality from probability of burn severity rating
            mort=np.zeros(ind.size,dtype='int16')
            for iyr in range(ind.size):
                #yr=d['Year'][iyr]
                b=beta_obs.copy()
                mort[iyr]=GetMortalityFromIBMSeverity(1,b)
            d['Mortality']=mort
            
            ibm_sim0[iStand]=d    
        
        del Oc_sim
    
    # Add to list of ensembles
    ibm_sim.append(ibm_sim0)
    
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
