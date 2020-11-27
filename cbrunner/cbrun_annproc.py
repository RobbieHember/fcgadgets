
#%% Import python modules

import numpy as np
from fcgadgets.activities.nutrient_application import update_nutrient_status 
from fcgadgets.taz import general_stat_models as gensm

#%% Define class objects

# Used to convert species parameter dictionaries into parameter vectors
class Bunch(object):
    def __init__(self, adict):
        self.__dict__.update(adict)

#%% Biomass dynamics

def BiomassFromBatchTIPSY(iScn,iT,vi,vo,psl,meta,iEP):
    
    # Update stand age
    vo['A'][iT,:]=vo['A'][iT-1,:]+1
     
    #--------------------------------------------------------------------------
    # Net growth of aboveground biomass
    #--------------------------------------------------------------------------
    
    # Index to growth curves at current age of stand
    iAge=np.minimum(vo['A'][iT,:],meta['GC']['BatchTIPSY Maximum Age'])-1
    
    # Convert to integer
    iAge=iAge.astype(int)
    
    # Extract net growth from growth curves
    NetGrowth=np.zeros((meta['N Stand'],6))
    for iS in range(meta['N Stand']):
        NetGrowth[iS,:]=vi['GC']['Active'][iAge[iS],iS,:]
    
    # Add net growth to output variable structure and apply scale factor
    # Oddly, using meta['iEP']['BiomassAboveground'] will invert the dimensions 
    # of C_G_Net - don't change it.
    vo['C_G_Net'][iT,:,0:5]=meta['GC']['Scale Factor']*NetGrowth[:,0:5].astype(float)
    
    # Total net growth of root biomass (Li et al. 2003, Eq. 4)
    G_Net_Root_Total=0.22*np.sum(vo['C_G_Net'][iT,:,0:5],axis=1)
    
    # Fine root fraction should decline with total biomass, but estimating it based on
    # size causes a lot of problems. Think about creating new equation.        
    vo['C_G_Net'][iT,:,iEP['RootCoarse']]=(1-0.072)*G_Net_Root_Total
    vo['C_G_Net'][iT,:,iEP['RootFine']]=0.072*G_Net_Root_Total
    
    # Update stemwood merchantable volume
    vo['V_StemMerch'][iT,:]=vo['V_StemMerch'][iT-1,:]+NetGrowth[:,5].astype(float)*meta['GC']['Scale Factor']      

    #--------------------------------------------------------------------------
    # Nutrient application effects to net growth
    #--------------------------------------------------------------------------
    
    # Index to stands that are stimulated by nutrient application
    meta['NM']['iApplication']=np.where( meta['NM']['ResponseCounter']>0 )[0]
    
    # Adjust N application response counter
    if meta['NM']['iApplication'].size>0:
        vi,vo,meta=update_nutrient_status(vi,vo,iT,meta,psl,'UpdateCounter')
    
    # Adjust root net growth
    if meta['NM']['iApplication'].size>0:
        vi,vo,meta=update_nutrient_status(vi,vo,iT,meta,psl,'BelowgroundNetGrowth')

    #--------------------------------------------------------------------------
    # Add net growth to biomass pools
    #--------------------------------------------------------------------------
    
    # Stemwood    
    vo['C_Eco_Pools'][iT,:,iEP['StemMerch']]=np.maximum(0,vo['C_Eco_Pools'][iT-1,:,iEP['StemMerch']]+vo['C_G_Net'][iT,:,iEP['StemMerch']])
    vo['C_Eco_Pools'][iT,:,iEP['StemNonMerch']]=np.maximum(0,vo['C_Eco_Pools'][iT-1,:,iEP['StemNonMerch']]+vo['C_G_Net'][iT,:,iEP['StemNonMerch']])
     
    # Foliage    
    # *** Notes 1: TIPSY foliage appears to be really low, use alternative equations 
    # instead. As net growth of foliage has been redefined, we must recalculate
    # net growth. Adjustment is applied in utilities script. ***
    vo['C_Eco_Pools'][iT,:,iEP['Foliage']]=np.maximum(0,vo['C_Eco_Pools'][iT-1,:,iEP['Foliage']]+vo['C_G_Net'][iT,:,iEP['Foliage']])
    
    # Branches
    vo['C_Eco_Pools'][iT,:,iEP['Branch']]=np.maximum(0,vo['C_Eco_Pools'][iT-1,:,iEP['Branch']]+vo['C_G_Net'][iT,:,iEP['Branch']])
    
    # Bark
    vo['C_Eco_Pools'][iT,:,iEP['Bark']]=np.maximum(0,vo['C_Eco_Pools'][iT-1,:,iEP['Bark']]+vo['C_G_Net'][iT,:,iEP['Bark']])
    
    # Coarse roots
    vo['C_Eco_Pools'][iT,:,iEP['RootCoarse']]=np.maximum(0,vo['C_Eco_Pools'][iT-1,:,iEP['RootCoarse']]+vo['C_G_Net'][iT,:,iEP['RootCoarse']])
    
    # Fine roots
    vo['C_Eco_Pools'][iT,:,iEP['RootFine']]=np.maximum(0,vo['C_Eco_Pools'][iT-1,:,iEP['RootFine']]+vo['C_G_Net'][iT,:,iEP['RootFine']])
  
    # Old root process based on size instead of growth - caused many problems
    flg=0
    if flg==1:
        
        # Calculate aboveground biomass
        vo['C_BiomassAG'][iT,:]=np.sum(vo['C_Eco_Pools'][iT,:,0:5],axis=1)
    
        # Total root biomass (Li et al. 2003, Eq. 4)
        BiomassRootTotal=0.222*vo['C_BiomassAG'][iT,:]
        
        # Fine root biomass (Li et al. 2003, Eq. 6)
        Pf=0.072+0.354*np.exp(-0.06*(2*BiomassRootTotal))
        vo['C_Eco_Pools'][iT,:,iEP['RootFine']]=Pf*BiomassRootTotal
  
        # Coarse root biomass
        Pc=1-Pf
        vo['C_Eco_Pools'][iT,:,iEP['RootCoarse']]=np.maximum(0,Pc*BiomassRootTotal)
    
        # Calculate net growth of roots from change in pools
        #Gnet_RC=vo['C_Eco_Pools'][iT,:,iEP['RootCoarse']]-vo['C_Eco_Pools'][iT-1,:,iEP['RootCoarse']]
        #Gnet_RF=vo['C_Eco_Pools'][iT,:,iEP['RootFine']]-vo['C_Eco_Pools'][iT-1,:,iEP['RootFine']]    
    
    #--------------------------------------------------------------------------
    # Biomass turnover
    #--------------------------------------------------------------------------
    
    # Calculate biomass loss due to mortality
    vo['C_M_Reg'][iT,:,0]=psl['bTR_StemMerch'][0,:]*vo['C_Eco_Pools'][iT,:,0]
    vo['C_M_Reg'][iT,:,1]=psl['bTR_StemMerch'][0,:]*vo['C_Eco_Pools'][iT,:,1]
    vo['C_M_Reg'][iT,:,2]=psl['bTR_StemMerch'][0,:]*vo['C_Eco_Pools'][iT,:,2]
    vo['C_M_Reg'][iT,:,3]=psl['bTR_StemMerch'][0,:]*vo['C_Eco_Pools'][iT,:,3]
    vo['C_M_Reg'][iT,:,4]=psl['bTR_StemMerch'][0,:]*vo['C_Eco_Pools'][iT,:,4]
    vo['C_M_Reg'][iT,:,5]=psl['bTR_StemMerch'][0,:]*vo['C_Eco_Pools'][iT,:,5]
    vo['C_M_Reg'][iT,:,6]=psl['bTR_StemMerch'][0,:]*vo['C_Eco_Pools'][iT,:,6]
    
    # Adjust mortality to account for N application response
    if meta['NM']['iApplication'].size>0:
        vi,vo,meta=update_nutrient_status(vi,vo,iT,meta,psl,'Mortality')       
    
    # Apply optional factor to mortality
    #if meta['Scenario'][iScn]['Status Mortality Factor']=='On':
    #    mortality_factor=meta['Scenario'][iScn]['Mortality Factor'][iT]
    #    if mortality_factor!=0:        
    #        # Normal operation
    #        vo['C_M_Reg'][iT,:,:]=net_growth_factor*vo['C_M_Reg'][iT,:,:]
    
    #--------------------------------------------------------------------------
    # If TIPSY indicates negative net growth (e.g., alder break-up), revise so
    # that net growth during the period of negative values equals the value in the
    # previous timestep.
    #--------------------------------------------------------------------------
    
    # Keep track of what regular mortality would be before it is affected by 
    # correction for catastrophic mortality
    C_M_Reg=vo['C_M_Reg'][iT,:,:].copy()
    
    # Define a threshold level of negative growth so it isn't triggered by noise
    NegGrowthThreshold=-0.25
    #NegGrowthThreshold=-1e6
    
    # Find stands with negative net growth    
    iNegNetG=np.where(vo['C_G_Net'][iT,:,0]<NegGrowthThreshold)[0]
    
    if iNegNetG.size>0:
        
        # If it is the first instance of negative net growth: 91) record net growth
        # of the preceeding timestep and (2) set flag = 1.
        iSwitchFlag=np.where(meta['FlagNegNetGrowth'][iNegNetG]==0)[0]
        meta['FlagNegNetGrowth'][iNegNetG[iSwitchFlag]]=1
        meta['G_Net_PriorToBreakup'][iNegNetG[iSwitchFlag],0:7]=vo['C_G_Net'][iT-1,iNegNetG[iSwitchFlag],0:7]
                        
        d=vo['C_G_Net'][iT,iNegNetG,0:7]-meta['G_Net_PriorToBreakup'][iNegNetG,:]
        
        CToTransfer=np.zeros((meta['N Stand'],7))
        CToTransfer[iNegNetG,:]=-1*d
        
        vo['C_G_Net'][iT,:,0:7]=vo['C_G_Net'][iT,:,0:7]+CToTransfer
        vo['C_M_Reg'][iT,:,0:7]=vo['C_M_Reg'][iT,:,0:7]+CToTransfer
        
        # Logbook entry
        for i in range(iNegNetG.size):
            txt='Scenario:' + str(iScn) + ', Stand:' + str(iNegNetG[i]) + ', Time:' + str(iT) + ', Negative net growth received from GY model is being adjusted.'
            meta['Logbook'].append(txt)        
    
    #--------------------------------------------------------------------------
    # Litterfall
    #--------------------------------------------------------------------------
    
    # Calculate foliage biomass turnover due to litterfall
    vo['C_LF'][iT,:,iEP['Foliage']]=psl['bTR_Foliage'][0,:]*vo['C_Eco_Pools'][iT,:,iEP['Foliage']]
    
    # Calculate branch biomass turnover due to litterfall
    vo['C_LF'][iT,:,iEP['Branch']]=psl['bTR_Branch'][0,:]*vo['C_Eco_Pools'][iT,:,iEP['Branch']]
    
    # Calculate bark biomass turnover due to litterfall
    vo['C_LF'][iT,:,iEP['Bark']]=psl['bTR_Bark'][0,:]*vo['C_Eco_Pools'][iT,:,iEP['Bark']]
    
    # Calculate coarse root biomass turnover due to litterfall
    vo['C_LF'][iT,:,iEP['RootCoarse']]=psl['bTR_RootCoarse'][0,:]*vo['C_Eco_Pools'][iT,:,iEP['RootCoarse']]
    
    # Calculate fine root biomass turnover due to litterfall
    vo['C_LF'][iT,:,iEP['RootFine']]=psl['bTR_RootFine'][0,:]*vo['C_Eco_Pools'][iT,:,iEP['RootFine']]
        
    # Adjust litterfall to account for N application response
    if meta['NM']['iApplication'].size>0:
        vi,vo,meta=update_nutrient_status(vi,vo,iT,meta,psl,'Litterfall')
    
    #--------------------------------------------------------------------------
    # Update summary variables
    #--------------------------------------------------------------------------
    
    # Update gross growth
    vo['C_G_Gross'][iT,:,:]=vo['C_G_Net'][iT,:,:]+C_M_Reg
    
    # Update NPP
    vo['C_NPP'][iT,:,:]=vo['C_G_Net'][iT,:,:]+C_M_Reg+vo['C_LF'][iT,:,:]
    
    return vo

#%% Annual tree biomass dynamics from Sawtooth

def BiomassFromSawtooth(iScn,iS,vi,vo,meta,iEP):
    
    #--------------------------------------------------------------------------
    # Define demensions
    #--------------------------------------------------------------------------
    
    N_t=meta['N Time']
    N_tree=meta['Max Num Trees Per Stand']
        
    #--------------------------------------------------------------------------
    # A vector of zeros
    #--------------------------------------------------------------------------
    
    ZerosYear=np.zeros((1,N_tree))
    
    #--------------------------------------------------------------------------
    # Function for calculating the difference
    #--------------------------------------------------------------------------
    
    def diff(x):    
        d=x[1:,:]-x[0:-1,:]
        y=np.concatenate((ZerosYear,d))
        return y
    
    #--------------------------------------------------------------------------
    # Multiplier to convert kg C to Mg C
    #--------------------------------------------------------------------------
    
    cm=1000
    
    #--------------------------------------------------------------------------
    # Initialize tree-level variables
    #--------------------------------------------------------------------------
    
    tl_ID_Srs=np.zeros((N_t,N_tree),dtype=np.int)
    tl_ID_Decid=np.zeros((N_t,N_tree),dtype=np.int)
    tl_A=np.zeros((N_t,N_tree))
    tl_H=np.zeros((N_t,N_tree))
    tl_D=np.zeros((N_t,N_tree))
    tl_N_R=np.zeros((N_t,N_tree))
    tl_N_M_Inv_Fir=np.zeros((N_t,N_tree))
    tl_N_M_Inv_Ins=np.zeros((N_t,N_tree))
    tl_N_M_Inv_Pat=np.zeros((N_t,N_tree))
    tl_N_M_Inv_Har=np.zeros((N_t,N_tree))
    tl_N_M_Inv_Win=np.zeros((N_t,N_tree))
    tl_N_M_Sim_Reg=np.zeros((N_t,N_tree))
    tl_N_M_Sim_Fir=np.zeros((N_t,N_tree))
    tl_N_M_Sim_Ins=np.zeros((N_t,N_tree))
    tl_N_M_Sim_Pat=np.zeros((N_t,N_tree))
    tl_N_M_Sim_Har=np.zeros((N_t,N_tree))
    tl_N_M_Sim_Win=np.zeros((N_t,N_tree))             
    tl_Csw=np.nan*np.zeros((N_t,N_tree))
    tl_Csw_Larger=np.nan*np.zeros((N_t,N_tree))
    tl_Csw_G=np.nan*np.zeros((N_t,N_tree))
    tl_Cf=np.zeros((N_t,N_tree))
    tl_Cbk=np.zeros((N_t,N_tree))
    tl_Cbr=np.zeros((N_t,N_tree))
    tl_Crc=np.zeros((N_t,N_tree))
    tl_Crf=np.zeros((N_t,N_tree))
    
    #--------------------------------------------------------------------------
    # Populate species ID
    #--------------------------------------------------------------------------
    
    rp=np.random.permutation(N_tree).astype(int)
    
    # Species 1
    n1=np.ceil(vi['Inv']['Srs1_Pct'][0,iS]/100*N_tree).astype(int)
    tl_ID_Srs[:,rp[0:n1]]=vi['Inv']['Srs1_ID'][0,iS]
    # Species 2
    if vi['Inv']['Srs2_Pct'][0,iS]>0:
        n2=np.ceil(vi['Inv']['Srs2_Pct'][0,iS]/100*N_tree).astype(int)
        tl_ID_Srs[:,rp[n1:n2]]=vi['Inv']['Srs2_ID'][0,iS]
    # Species 3
    if vi['Inv']['Srs3_Pct'][0,iS]>0:
        tl_ID_Srs[:,rp[n2:]]=vi['Inv']['Srs3_ID'][0,iS]
    
    #--------------------------------------------------------------------------
    # Initialize parameter vectors
    #--------------------------------------------------------------------------
    
    # Extract stand-level parameters from metadata structure
    ptl=Bunch(meta['ptl'])
    
    # Allometry
    bA=Bunch(ptl.Allom[0])
    for k in bA.__dict__.keys():
        setattr(bA,k,np.zeros((1,N_tree))) 
    
    # Recruitment
    iEqR=np.where(np.asarray(ptl.R_Type)==meta['Scenario'][iScn]['EquationSet_R_CD'])[0][0]
    bR=Bunch(ptl.R_Coef[iEqR][0])
    for k in bR.__dict__.keys():
        setattr(bR,k,np.zeros((1,N_tree)))
     
    # Mortality
    iEqM=np.where(np.asarray(ptl.M_Type)==meta['Scenario'][iScn]['EquationSet_M_CD'])[0][0]
    bM=Bunch(ptl.M_Coef[iEqM][0])
    for k in bM.__dict__.keys():
        setattr(bM,k,np.zeros((1,N_tree)))
    
    # Growth
    iEqG=np.where(np.asarray(ptl.G_Type)==meta['Scenario'][iScn]['EquationSet_G_CD'])[0][0]
    bG=Bunch(ptl.G_Coef[iEqG][0])
    for k in bG.__dict__.keys():
        setattr(bG,k,np.zeros((1,N_tree)))
    
    #--------------------------------------------------------------------------
    # Populate parameter vectors based on species to each tree based on fractions 
    # from inventory.
    #--------------------------------------------------------------------------
    
    uS=np.unique(tl_ID_Srs)
    for i in range(0,len(uS)):
        
        ind=np.where(tl_ID_Srs[0,:]==uS[i])[0]
        
        # Allometry
        for k in bA.__dict__.keys():
            if k!='SRS_CD':
                tmp=getattr(bA,k)
                tmp[0,ind]=ptl.Allom[uS[i]-1][k]
                setattr(bA,k,tmp)
        
        # Recruitment
        for k in bR.__dict__.keys():
            if k!='SRS_CD':
                tmp=getattr(bR,k)
                tmp[0,ind]=ptl.R_Coef[iEqR][uS[i]-1][k]
                setattr(bR,k,tmp)
        
        # Mortality
        for k in bM.__dict__.keys():
            if k!='SRS_CD':
                tmp=getattr(bM,k)
                tmp[0,ind]=ptl.M_Coef[iEqM][uS[i]-1][k]
                setattr(bM,k,tmp)
        
        # Growth
        for k in bG.__dict__.keys():
            if k!='SRS_CD':
                tmp=getattr(bG,k)
                tmp[0,ind]=ptl.G_Coef[iEqG][uS[i]-1][k]
                setattr(bG,k,tmp)
        
    #--------------------------------------------------------------------------
    # Germination
    #--------------------------------------------------------------------------
    
    # In the first time step, assume 1...seed_n seedlings germinate with initial
    # biomass drawn from a normal distribution, with a very small mean
    # biomass and a coefficient of variation of 10 percent. Force values to be positive.
    # Once a few live trees have been Invscribed, recruitment rate should
    # rise quickly because stand biomass (the primary driver of recruitment)
    # is low. This should lead to repopulation of the stand assuming mortality 
    # does not exceed recruitment, and growth is viable.
  
    # Planted:
    # seed_n=1100;
    # seed_mu=0.06; % (kg C tree-1)
    
    # Natural:
    seed_n=200
    seed_mu=0.001 # (kg C tree-1)
    seed_cf=10 # (%)
    seed_mx=0.000001
    Csw_seed=np.random.normal(seed_mu,seed_mu*(10/100),(1,seed_n))    
    Csw_seed=np.maximum(seed_mx,Csw_seed)
    tl_Csw[0,0:seed_n]=Csw_seed
    
    # *** This needs updating ***
    tl_H[0,0:seed_n]=bA.Cag2H1[0,0:seed_n]*((1-np.exp(-bA.Cag2H2[0,0:seed_n]*tl_Csw[0,0:seed_n])))**(1/(1-bA.Cag2H3[0,0:seed_n]))
  
    #--------------------------------------------------------------------------
    # Populate tree matrix with environmental variables for stand, s  
    #--------------------------------------------------------------------------
    
    # Working on it...
    
        
    #--------------------------------------------------------------------------
    # Loop through time intervals (start in second time step)
    #--------------------------------------------------------------------------
    
    for iT in range(1,N_t):
        
        #----------------------------------------------------------------------
        # Apply first disturbance cycle to subsequent disturbance cycles during 
        # the spinup period to save time.
        #----------------------------------------------------------------------
        
        if (vi['tv'][iT]>meta['SpinupSpanFastTrack'][0]+1) & (vi['tv'][iT]<=meta['SpinupSpanFastTrack'][1]+1):
            
            iT0=iT-meta['Spinup Disturbance Return Inverval']
            
            tl_A[iT,:]=tl_A[iT0,:]
            tl_H[iT,:]=tl_H[iT0,:]
            tl_D[iT,:]=tl_D[iT0,:]
            tl_N_R[iT,:]=tl_N_R[iT0,:]
            tl_N_M_Inv_Fir[iT,:]=tl_N_M_Inv_Fir[iT0,:]
            tl_N_M_Inv_Ins[iT,:]=tl_N_M_Inv_Ins[iT0,:]
            tl_N_M_Inv_Pat[iT,:]=tl_N_M_Inv_Pat[iT0,:]
            tl_N_M_Inv_Har[iT,:]=tl_N_M_Inv_Har[iT0,:]
            tl_N_M_Inv_Win[iT,:]=tl_N_M_Inv_Win[iT0,:]
            tl_N_M_Sim_Reg[iT,:]=tl_N_M_Sim_Reg[iT0,:]
            tl_N_M_Sim_Fir[iT,:]=tl_N_M_Sim_Fir[iT0,:]
            tl_N_M_Sim_Ins[iT,:]=tl_N_M_Sim_Ins[iT0,:]
            tl_N_M_Sim_Pat[iT,:]=tl_N_M_Sim_Pat[iT0,:]
            tl_N_M_Sim_Har[iT,:]=tl_N_M_Sim_Har[iT0,:]
            tl_N_M_Sim_Win[iT,:]=tl_N_M_Sim_Win[iT0,:]
            tl_Csw[iT,:]=tl_Csw[iT0,:]
            tl_Csw_Larger[iT,:]=tl_Csw_Larger[iT0,:]
            tl_Csw_G[iT,:]=tl_Csw_G[iT0,:]         
            vo['A'][iT,iS]=vo['A'][iT0,iS]
            
            continue
        
        #----------------------------------------------------------------------
        
        # Index to live trees    
        iLive=np.array(np.where(~np.isnan(tl_Csw[iT-1,:]))).flatten()
        nLive=iLive.size
            
        # Index to dead trees    
        iDead=np.array(np.where(np.isnan(tl_Csw[iT-1,:])))
        nDead=iDead.size
        
        # Create random number vectors that will be used for calculation of 
        # annual probabilities of recruitment and mortality
        rLive=np.random.uniform(0,1,nLive)
        rDead=np.random.uniform(0,1,nDead)
   
        # Update stand age (i.e., time since stand-replacing disturbance)
        vo['A'][iT,iS]=vo['A'][iT-1,iS]+1
    
        # Update tree age    
        tl_A[iT,iLive]=tl_A[iT-1,iLive]+1
        
        #----------------------------------------------------------------------
        # Calculate predictor variables used in the equations of recruitment, 
        # mortality and growth
        #----------------------------------------------------------------------
        
        # Tree age
        A=tl_A[iT-1,:]
    
        # Aboveground biomass of individual trees from t-1 (kg C tree-1)
        Csw=tl_Csw[iT-1,:]
    
        # Tree height (m)
        H=tl_H[iT-1,:]
    
        # Stand age
        SA=np.nanmean(tl_A[iT-1,:])
    
        # Stand-level biomass from t-1 (Mg C ha-1)
        SCsw=np.nansum(tl_Csw[iT-1,:])/cm
        
        # Stand density from t-1 (stems ha-1)
        SN=nLive
          
        # Biomass of larger trees (Mg C ha-1)        
        ListToSort=np.zeros((N_tree,2))
        ListToSort[:,0]=np.arange(0,N_tree,1)
        ListToSort[:,1]=Csw/cm
        SortedList=ListToSort[ListToSort[:,1].argsort()]
        SortedList=np.flip(SortedList,0)
        tmp=np.reshape(np.nancumsum(SortedList[:,1])-SortedList[:,1]/1000,(N_tree,1))
        SortedList2=np.append(SortedList,tmp,axis=1)
        
        Csw_Larger=np.zeros((1,N_tree))
        Csw_Larger[0,SortedList2[:,0].astype(int)]=SortedList2[:,2]
        
        tl_Csw_Larger[iT,:]=Csw_Larger
    
        #----------------------------------------------------------------------
        # Probability of recruitment
        #----------------------------------------------------------------------
        
        # Only do this if there are living trees
        if iDead.size!=0:
                 
            if meta['Scenario'][iScn]['EquationSet_R_CD']=='Def1':                
                SCsw_z=(SCsw-bR.SB_mu)/bR.SB_sig
                lgit=bR.Int+bR.SB*SCsw_z
                Pr=(np.exp(lgit)/(1+np.exp(lgit)))
                
            # Establish trees based on annual probability of recruitment.
            # Initial values of biomass and height are set low, arbitrary until
            # they can be set according to actual observations 
            # at end of first year.
      
            iRec=np.where(Pr[0,iDead]>=rDead)[1]
            tl_A[iT:,iDead[0,iRec]]=1
            tl_H[iT:,iDead[0,iRec]]=0.1
            tl_Csw[iT:,iDead[0,iRec]]=0.05
            tl_N_R[iT,iDead[0,iRec]]=1
         
        #----------------------------------------------------------------------    
        # Growth of stemwood biomass (kg C tree-1 yr-1)
        #----------------------------------------------------------------------
        
        if meta['Scenario'][iScn]['EquationSet_G_CD']=='Def1':
            
            # Standardization
            LnCsw_z=(np.log(Csw)-bG.LnB_mu)/bG.LnB_sig
            Csw_z=(Csw-bG.B_mu)/bG.B_sig
            SA_z=(A-bG.SA_mu)/bG.SA_sig
            SCswLT_z=(Csw_Larger-bG.SBLT_mu)/bG.SBLT_sig
            SCsw_z=(SCsw-bG.SB_mu)/bG.SB_sig                
               
            # Add all effects to intercept
            yhat=bG.Int + bG.LnB*LnCsw_z + bG.B*Csw_z + bG.SA*SA_z + bG.SBLT*SCswLT_z + bG.SB*SCsw_z
          
            # Back-transform and apply log correction
            yhat=np.exp(yhat)
              
            # Cap unrealistic growth
            #yhat[yhat>G_max]=G_max
        
            # Populate tree level structure with growth predictions
            tl_Csw_G[iT,iLive]=yhat[0,iLive]
            
            #******************************************************************
            # SPECIAL ORDER
            # Growth modifier for EP703 study
            
            #if iScn==1:
            #    try:
            #        gf=meta['Scenario'][iScn]['Growth factor'][iT,1]*1000/iLive.size
            #    except:
            #        gf=1
            #    tl_Csw_G[iT,iLive]=tl_Csw_G[iT,iLive]+gf             
            
            # END SPECIAL ORDER
            #******************************************************************
            
        
        #----------------------------------------------------------------------    
        # Update state variables
        #----------------------------------------------------------------------    
        
        # Update stemwood biomass (kg C tree-1)
        tl_Csw[iT,iLive]=tl_Csw[iT-1,iLive]+tl_Csw_G[iT,iLive]
    
        # Update tree height (m)   
        # *** This needs updating ***
        tl_H[iT,:]=bA.Cag2H1*((1-np.exp(-bA.Cag2H2*tl_Csw[iT,:]))**(1/(1-bA.Cag2H3)))
        
        #----------------------------------------------------------------------    
        # Probability of tree mortality
        #----------------------------------------------------------------------    
        
        if meta['Scenario'][iScn]['EquationSet_M_CD']=='Def1':
            
            Csw_z=(Csw-bM.B_mu)/bM.B_sig
            Csw2_z=(Csw**2-bM.B2_mu)/bM.B2_sig
            SA_z=(A-bM.SA_mu)/bM.SA_sig
            SCswLT_z=(Csw_Larger-bM.SBLT_mu)/bM.SBLT_sig
            SCsw_z=(SCsw-bM.SB_mu)/bM.SB_sig

            lgit=bM.Int + bM.B*Csw_z + bM.B2*Csw2_z + bM.SA*SA_z + bM.SBLT*SCswLT_z + bM.SB*SCsw_z
            
            Pm_Sim_Reg=np.exp(lgit)/(1+np.exp(lgit))
            Pm_Sim_Ins=np.zeros((1,N_tree))
            Pm_Sim_Pat=np.zeros((1,N_tree))
        
        #---------------------------------------------------------------------- 
        # Remove biomass of trees that died
        #---------------------------------------------------------------------- 
        
        # Remove biomass of trees that died directly from competition and 
        # environmental conditions
        iKill=np.where(Pm_Sim_Reg[0,iLive]>=rLive)
        tl_A[iT:,iLive[iKill]]=np.nan
        tl_D[iT:,iLive[iKill]]=np.nan
        tl_H[iT:,iLive[iKill]]=np.nan
        tl_Csw[iT:,iLive[iKill]]=np.nan
        tl_N_M_Sim_Reg[iT,iLive[iKill]]=1   
        rLive[iKill]=10; # Update rLive to avoid double-counting mortality
        
        # Remove biomass of trees that die due to insect attack    
        iKill=np.where(Pm_Sim_Ins[0,iLive]>=rLive)
        tl_A[iT:,iLive[iKill]]=np.nan
        tl_D[iT:,iLive[iKill]]=np.nan
        tl_H[iT:,iLive[iKill]]=np.nan
        tl_Csw[iT:,iLive[iKill]]=np.nan
        tl_N_M_Sim_Ins[iT,iLive[iKill]]=1
        rLive[iKill]=10
    
        # Remove biomass of trees that die due to pathogen infection    
        iKill=np.where(Pm_Sim_Pat[0,iLive]>=rLive)
        tl_A[iT:,iLive[iKill]]=np.nan
        tl_D[iT:,iLive[iKill]]=np.nan
        tl_H[iT:,iLive[iKill]]=np.nan
        tl_Csw[iT:,iLive[iKill]]=np.nan
        tl_N_M_Sim_Pat[iT,iLive[iKill]]=1
        rLive[iKill]=10
        
        #---------------------------------------------------------------------- 
        # Mortality from inventory data sources
        #---------------------------------------------------------------------- 
        
        iDist=np.where((vi['EH'][iS]['Year']==vi['tv'][iT]) & (np.isnan(vi['EH'][iS]['MortalityFactor'])==False))[0]
        if iDist.size!=0:            
            
            iLive=np.where(np.isnan(tl_Csw[iT,:])==False)[0]
            
            nKill=int(vi['EH'][iS]['MortalityFactor'][iDist]/100*iLive.size)
                        
            if vi['EH'][iS]['ID_Type'][iDist]==meta['LUT Dist']['Wildfire']:
                # Wildfire
                iKill=np.arange(0,nKill)
                tl_N_M_Inv_Fir[iT,iLive[iKill]]=1
                vo['A'][iT,iS]=0
                
            elif (vi['EH'][iS]['ID_Type'][iDist]==meta['LUT Dist']['Harvest']) | (vi['EH'][iS]['ID_Type'][iDist]==meta['LUT Dist']['Harvest and Slashpile Burn']):
                iKill=np.arange(0,nKill)
                tl_N_M_Inv_Har[iT,iLive[iKill]]=1
                vo['A'][iT,iS]=0
                               
            elif (vi['EH'][iS]['ID_Type'][iDist]==meta['LUT Dist']['Beetles']) | (vi['EH'][iS]['ID_Type'][iDist]==meta['LUT Dist']['Defoliators']):
                # Insects
                iKill=np.arange(0,nKill)
                tl_N_M_Inv_Ins[iT,iLive[iKill]]=1
                
            elif (vi['EH'][iS]['ID_Type'][iDist]==meta['LUT Dist']['Planting']):
                # Do nothing
                pass
                
            elif (vi['EH'][iS]['ID_Type'][iDist]==meta['LUT Dist']['Fertilization Aerial']):
                # Do nothing
                pass
            
            #elif vi['EH'][iS]['ID_Type'][iDist]==meta['LUT Dist']['Wind']:
            #    # Wind
            #    iKill=np.arange(0,nKill)
            #    tl_N_M_Inv_Win[iT,iLive[iKill]]=1
            #    vo['A'][iT,iS]=0
                
            else:
                pass
            
            # Update tree variables
            tl_A[iT:,iLive[iKill]]=0
            tl_D[iT:,iLive[iKill]]=np.nan
            tl_H[iT:,iLive[iKill]]=np.nan
            tl_Csw[iT:,iLive[iKill]]=np.nan
        
    #--------------------------------------------------------------------------
    # Biomass of other tissues
    #--------------------------------------------------------------------------
    
    # Merchantable stemwood fraction
    msf=0.82
    
    # Merchantable stemwood
    tl_Csw_Merch=msf*tl_Csw
    
    # Non-merchantable stemwood
    tl_Csw_NonMerch=(1-msf)*tl_Csw
    
    # Foliage biomass
    tl_Cf=tl_Csw*np.tile(bA.Csw2Cf1,(N_t,1))*tl_Csw**np.tile(bA.Csw2Cf2,(N_t,1))
    
    # Branch biomass
    tl_Cbk=tl_Csw*np.tile(bA.Csw2Cbk1,(N_t,1))*tl_Csw**np.tile(bA.Csw2Cbk2,(N_t,1))
  
    # Branch biomass
    tl_Cbr=tl_Csw*np.tile(bA.Csw2Cbr1,(N_t,1))*tl_Csw**np.tile(bA.Csw2Cbr2,(N_t,1))
    
    # Total aboveground biomass
    tl_Cag=tl_Csw+tl_Cf+tl_Cbk+tl_Cbr
    
    # Total root biomass (Li et al. 2003)  
  
    # Conifer
    tl_Cr=0.222*tl_Cag
  
    # Deciduous (equation from Li et al. XXXX)
    ind=np.where(tl_ID_Decid==1)[0]
    if ind.size!=0:
        tl_Cr[ind]=1.576*tl_Cag[ind]**0.615  
  
    # Fine root biomass (equation from Li et al. XXXX)
    tl_Crf=tl_Cr*(0.072+0.354*np.exp(-0.06*(2*tl_Cr)))
  
    # Coarse root biomass
    tl_Crc=tl_Cr-tl_Crf
  
    # Total biomass
    tl_Ctot=tl_Cag+tl_Cr
  
    #--------------------------------------------------------------------------
    # Biomass loss due to tree mortality (by tissue)
    #--------------------------------------------------------------------------
    
    # Derive variables needed to calculate tree mortality  
    indM=np.where(np.isnan(tl_Ctot)==True)
        
    # Calculate total demographic mortality
    N_M=tl_N_M_Sim_Reg+tl_N_M_Sim_Fir+tl_N_M_Sim_Ins+tl_N_M_Sim_Pat+ \
        tl_N_M_Sim_Win+tl_N_M_Sim_Har+tl_N_M_Inv_Fir+tl_N_M_Inv_Ins+ \
        tl_N_M_Inv_Pat+tl_N_M_Inv_Win+tl_N_M_Inv_Har
    N_M=np.minimum(1,N_M)
    
    # Create a vector that contains the biomass that will be lost due to
    # mortality
    Cx=tl_Csw_Merch
    Cx[indM[0],indM[1]]=0
    #Cx[np.unravel_index(indM,Cx.shape,'C')]=0    
    C_M=np.concatenate((ZerosYear,np.diff(Cx,axis=0)),axis=0)
    C_M[C_M>0]=0
    C_M=-C_M
    tl_Csw_Merch_M=N_M*C_M
    
    Cx=tl_Csw_NonMerch
    Cx[indM[0],indM[1]]=0
    #Cx[np.unravel_index(indM,Cx.shape,'C')]=0    
    C_M=np.concatenate((ZerosYear,np.diff(Cx,axis=0)),axis=0)
    C_M[C_M>0]=0
    C_M=-C_M
    tl_Csw_NonMerch_M=N_M*C_M
    
    Cx=tl_Cf
    Cx[indM[0],indM[1]]=0
    #Cx[np.unravel_index(indM,Cx.shape,'C')]=0    
    C_M=np.concatenate((ZerosYear,np.diff(Cx,axis=0)),axis=0)
    C_M[C_M>0]=0
    C_M=-C_M
    tl_Cf_M=N_M*C_M
    
    Cx=tl_Cbr
    Cx[indM[0],indM[1]]=0
    #Cx[np.unravel_index(indM,Cx.shape,'C')]=0      
    C_M=np.concatenate((ZerosYear,np.diff(Cx,axis=0)),axis=0)
    C_M[C_M>0]=0
    C_M=-C_M
    tl_Cbr_M=N_M*C_M
    
    Cx=tl_Cbk
    Cx[indM[0],indM[1]]=0
    #Cx[np.unravel_index(indM,Cx.shape,'C')]=0      
    C_M=np.concatenate((ZerosYear,np.diff(Cx,axis=0)),axis=0)
    C_M[C_M>0]=0
    C_M=-C_M
    tl_Cbk_M=N_M*C_M
    
    Cx=tl_Crc
    Cx[indM[0],indM[1]]=0
    #Cx[np.unravel_index(indM,Cx.shape,'C')]=0    
    C_M=np.concatenate((ZerosYear,np.diff(Cx,axis=0)),axis=0)
    C_M[C_M>0]=0
    C_M=-C_M
    tl_Crc_M=N_M*C_M
    
    Cx=tl_Crf
    Cx[indM[0],indM[1]]=0
    #Cx[np.unravel_index(indM,Cx.shape,'C')]=0       
    C_M=np.concatenate((ZerosYear,np.diff(Cx,axis=0)),axis=0)
    C_M[C_M>0]=0
    C_M=-C_M
    tl_Crf_M=N_M*C_M
    
    #--------------------------------------------------------------------------
    # Biomass loss due to tree mortality (by causal factor)
    #--------------------------------------------------------------------------
       
    # Create a vector that contains the biomass that will be lost due to
    # mortality  
    Cx=tl_Ctot
    Cx[indM[0],indM[1]]=0
    #Cx[np.unravel_index(indM,Cx.shape,'C')]=0    
    C_M=np.concatenate((ZerosYear,np.diff(Cx,axis=0)),axis=0)
    C_M[C_M>0]=0
    C_M=-C_M
  
    tl_Ctot_M_Sim_Reg=tl_N_M_Sim_Reg*C_M
    tl_Ctot_M_Sim_Fir=tl_N_M_Sim_Fir*C_M
    tl_Ctot_M_Sim_Ins=tl_N_M_Sim_Ins*C_M
    tl_Ctot_M_Sim_Pat=tl_N_M_Sim_Pat*C_M
    tl_Ctot_M_Sim_Win=tl_N_M_Sim_Win*C_M
    tl_Ctot_M_Sim_Har=tl_N_M_Sim_Har*C_M
  
    tl_Ctot_M_Inv_Fir=tl_N_M_Inv_Fir*C_M
    tl_Ctot_M_Inv_Ins=tl_N_M_Inv_Ins*C_M
    tl_Ctot_M_Inv_Pat=tl_N_M_Inv_Pat*C_M
    tl_Ctot_M_Inv_Win=tl_N_M_Inv_Win*C_M
    tl_Ctot_M_Inv_Har=tl_N_M_Inv_Har*C_M
    
    #--------------------------------------------------------------------------
    # Biomass turnover (litterfall)
    #--------------------------------------------------------------------------
    
    # Foliage turnover
    tl_Cf_LF=0.100*tl_Cf
    
    # Branch turnover
    tl_Cbr_LF=0.035*tl_Cbr
    
    # Bark turnover
    tl_Cbk_LF=0.035*tl_Cbk
    
    # Fine root turnover
    tl_Crf_LF=0.641*tl_Crf
        
    # Coarse root turnover
    tl_Crc_LF=0.02*tl_Crc
            
    #--------------------------------------------------------------------------
    # Populate stand variables (individual-tree statistics)
    #--------------------------------------------------------------------------
  
    # Mean tree age
    vo['TreeMean_A'][:,iS]=np.nanmean(tl_A,axis=1) 
    
    # Mean tree diamter at breast height
    vo['TreeMean_D'][:,iS]=np.nanmean(tl_D,axis=1)
    
    # Mean tree height (m)
    vo['TreeMean_H'][:,iS]=np.nanmean(tl_H,axis=1)
    
    # Mean tree stemwood biomass (kg C)
    vo['TreeMean_Csw'][:,iS]=np.nanmean(tl_Csw,axis=1) 
    
    # Mean stemwood biomass growth (kg C yr-1)
    vo['TreeMean_Csw_G'][:,iS]=np.nanmean(tl_Csw_G,axis=1)
    
    #--------------------------------------------------------------------------
    # Populate stand demographics
    #--------------------------------------------------------------------------
    
    # Stand density (stems ha-1)    
    vo['N'][:,iS]=np.sum(~np.isnan(tl_Csw),axis=1)
    
    # Stand density at beginning of time step
    N_t0=np.append(vo['N'][0,iS],vo['N'][0:-1,iS])
    
    # Demographic recruitment rate (% yr-1)
    vo['N_R'][:,iS]=np.minimum(100,np.sum(tl_N_R,axis=1)/N_t0*100)
  
    # Demographic mortality rate (% yr-1)
          
    # Need to update N0 or else stand-replacing disturbances will have lower
    # relative mortality rate
    # *** The other solution is to move the prescribed disturbances above the
    # background mortality, but perhaps GCBM will do this anyways ***
    #N0=N0-Sum_N_M_Sim_Reg-Sum_N_M_Sim_Ins-Sum_N_M_Sim_Pat;
    
    N_M_Tot=tl_N_M_Sim_Reg+tl_N_M_Sim_Fir+tl_N_M_Sim_Ins+tl_N_M_Sim_Pat+ \
        tl_N_M_Sim_Win+tl_N_M_Sim_Har+tl_N_M_Inv_Fir+tl_N_M_Inv_Ins+ \
        tl_N_M_Inv_Pat+tl_N_M_Inv_Win+tl_N_M_Inv_Har
    
    vo['N_M_Tot'][:,iS]=np.minimum(100,np.sum(N_M_Tot,axis=1)/N_t0*100)
    
    vo['N_M_Sim_Reg'][:,iS]=np.maximum(0,np.sum(tl_N_M_Sim_Reg,axis=1)/N_t0*100)
    vo['N_M_Sim_Fir'][:,iS]=np.maximum(0,np.sum(tl_N_M_Sim_Fir,axis=1)/N_t0*100)
    vo['N_M_Sim_Ins'][:,iS]=np.maximum(0,np.sum(tl_N_M_Sim_Ins,axis=1)/N_t0*100)
    vo['N_M_Sim_Pat'][:,iS]=np.maximum(0,np.sum(tl_N_M_Sim_Pat,axis=1)/N_t0*100)
    vo['N_M_Sim_Win'][:,iS]=np.maximum(0,np.sum(tl_N_M_Sim_Win,axis=1)/N_t0*100)
    vo['N_M_Sim_Har'][:,iS]=np.maximum(0,np.sum(tl_N_M_Sim_Har,axis=1)/N_t0*100)
  
    vo['N_M_Inv_Fir'][:,iS]=np.maximum(0,np.sum(tl_N_M_Inv_Fir,axis=1)/N_t0*100)
    vo['N_M_Inv_Ins'][:,iS]=np.maximum(0,np.sum(tl_N_M_Inv_Ins,axis=1)/N_t0*100)
    vo['N_M_Inv_Pat'][:,iS]=np.maximum(0,np.sum(tl_N_M_Inv_Pat,axis=1)/N_t0*100)
    vo['N_M_Inv_Win'][:,iS]=np.maximum(0,np.sum(tl_N_M_Inv_Win,axis=1)/N_t0*100)
    vo['N_M_Inv_Har'][:,iS]=np.maximum(0,np.sum(tl_N_M_Inv_Har,axis=1)/N_t0*100)
    
    #--------------------------------------------------------------------------
    # Populate stand gross growth of biomass (Mg C ha-1 yr-1)
    #--------------------------------------------------------------------------
    
    vo['C_G_Gross'][:,iS,0]=np.nansum(np.maximum(0,diff(tl_Csw_Merch)),axis=1)/cm
    vo['C_G_Gross'][:,iS,1]=np.nansum(np.maximum(0,diff(tl_Csw_NonMerch)),axis=1)/cm
    vo['C_G_Gross'][:,iS,2]=np.nansum(np.maximum(0,diff(tl_Cf)),axis=1)/cm
    vo['C_G_Gross'][:,iS,3]=np.nansum(np.maximum(0,diff(tl_Cbr)),axis=1)/cm
    vo['C_G_Gross'][:,iS,4]=np.nansum(np.maximum(0,diff(tl_Cbk)),axis=1)/cm
    vo['C_G_Gross'][:,iS,5]=np.nansum(np.maximum(0,diff(tl_Crc)),axis=1)/cm
    vo['C_G_Gross'][:,iS,6]=np.nansum(np.maximum(0,diff(tl_Crf)),axis=1)/cm
    
    #--------------------------------------------------------------------------
    # Populate stand biomass turnover (Mg C ha-1 yr-1)
    #--------------------------------------------------------------------------
    
    vo['C_LF'][:,iS,2]=np.nansum(tl_Cf_LF,axis=1)/cm
    vo['C_LF'][:,iS,3]=np.nansum(tl_Cbr_LF,axis=1)/cm
    vo['C_LF'][:,iS,4]=np.nansum(tl_Cbk_LF,axis=1)/cm   
    vo['C_LF'][:,iS,5]=np.nansum(tl_Crc_LF,axis=1)/cm
    vo['C_LF'][:,iS,6]=np.nansum(tl_Crf_LF,axis=1)/cm
    
    #--------------------------------------------------------------------------
    # Populate stand biomass loss due to mortality (Mg C ha-1 yr-1)
    #--------------------------------------------------------------------------
    
    vo['C_M_Reg'][:,iS,0]=np.nansum(tl_Csw_Merch_M,axis=1)/cm
    vo['C_M_Reg'][:,iS,1]=np.nansum(tl_Csw_NonMerch_M,axis=1)/cm
    vo['C_M_Reg'][:,iS,2]=np.nansum(tl_Cf_M,axis=1)/cm
    vo['C_M_Reg'][:,iS,3]=np.nansum(tl_Cbr_M,axis=1)/cm
    vo['C_M_Reg'][:,iS,4]=np.nansum(tl_Cbk_M,axis=1)/cm
    vo['C_M_Reg'][:,iS,5]=np.nansum(tl_Crc_M,axis=1)/cm
    vo['C_M_Reg'][:,iS,6]=np.nansum(tl_Crf_M,axis=1)/cm
    
    vo['C_M_Sim_Reg'][:,iS]=np.nansum(tl_Ctot_M_Sim_Reg,axis=1)/cm
    vo['C_M_Sim_Fir'][:,iS]=np.nansum(tl_Ctot_M_Sim_Fir,axis=1)/cm
    vo['C_M_Sim_Ins'][:,iS]=np.nansum(tl_Ctot_M_Sim_Ins,axis=1)/cm
    vo['C_M_Sim_Pat'][:,iS]=np.nansum(tl_Ctot_M_Sim_Pat,axis=1)/cm
    vo['C_M_Sim_Win'][:,iS]=np.nansum(tl_Ctot_M_Sim_Win,axis=1)/cm
    vo['C_M_Sim_Har'][:,iS]=np.nansum(tl_Ctot_M_Sim_Har,axis=1)/cm
  
    vo['C_M_Inv_Fir'][:,iS]=np.nansum(tl_Ctot_M_Inv_Fir,axis=1)/cm
    vo['C_M_Inv_Ins'][:,iS]=np.nansum(tl_Ctot_M_Inv_Ins,axis=1)/cm
    vo['C_M_Inv_Pat'][:,iS]=np.nansum(tl_Ctot_M_Inv_Pat,axis=1)/cm
    vo['C_M_Inv_Win'][:,iS]=np.nansum(tl_Ctot_M_Inv_Win,axis=1)/cm
    vo['C_M_Inv_Har'][:,iS]=np.nansum(tl_Ctot_M_Inv_Har,axis=1)/cm
            
    #--------------------------------------------------------------------------
    # Populate stand NPP by tissue (Mg C ha-1 yr-1)
    #--------------------------------------------------------------------------
    
    vo['C_NPP'][:,iS,0:7]=vo['C_G_Gross'][:,iS,0:7]+vo['C_LF'][:,iS,0:7]
    
    #--------------------------------------------------------------------------
    # Populate stand net growth by tissue (Mg C ha-1 yr-1)
    #--------------------------------------------------------------------------
    
    vo['C_G_Net'][:,iS,0:7]=vo['C_G_Gross'][:,iS,0:7]-vo['C_M_Reg'][:,iS,0:7]
    
    #--------------------------------------------------------------------------
    # Populate ecosystem carbon pools by tissue (Mg C ha-1 yr-1)
    #--------------------------------------------------------------------------
    
    vo['C_Eco_Pools'][:,iS,0]=np.maximum(0,np.nansum(tl_Csw_Merch,axis=1))/cm
    vo['C_Eco_Pools'][:,iS,1]=np.maximum(0,np.nansum(tl_Csw_NonMerch,axis=1))/cm
    vo['C_Eco_Pools'][:,iS,2]=np.maximum(0,np.nansum(tl_Cf,axis=1))/cm
    vo['C_Eco_Pools'][:,iS,3]=np.maximum(0,np.nansum(tl_Cbr,axis=1))/cm
    vo['C_Eco_Pools'][:,iS,4]=np.maximum(0,np.nansum(tl_Cbk,axis=1))/cm
    vo['C_Eco_Pools'][:,iS,5]=np.maximum(0,np.nansum(tl_Crc,axis=1))/cm
    vo['C_Eco_Pools'][:,iS,6]=np.maximum(0,np.nansum(tl_Crf,axis=1))/cm
        
    return vo

#%% STAND DEAD ORGANIC MATTER DYNAMICS

def DOM_like_CBM08(iT,vi,vo,psl,iEP,meta):
     
    #--------------------------------------------------------------------------
    # Flux of carbon between biomass components and dead organic matter components
    #--------------------------------------------------------------------------
    
    # Transfer biomass turnover to very fast litter pool
    vo['C_Eco_Pools'][iT,:,iEP['LitterVF']]=vo['C_Eco_Pools'][iT-1,:,iEP['LitterVF']]+ \
        psl['bIPF_FoliageLitToLitterVF']*vo['C_LF'][iT,:,iEP['Foliage']]+ \
        psl['bIPF_RootFineLitToLitterVF']*vo['C_LF'][iT,:,iEP['RootFine']]+ \
        psl['bIPF_FoliageMorToLitterVF']*vo['C_M_Reg'][iT,:,iEP['Foliage']]
    
    # Transfer biomass turnover to fast litter pool
    vo['C_Eco_Pools'][iT,:,iEP['LitterF']]=vo['C_Eco_Pools'][iT-1,:,iEP['LitterF']]+ \
        psl['bIPF_BranchLitToLitterF']*vo['C_LF'][iT,:,iEP['Branch']]+ \
        psl['bIPF_BarkLitToLitterF']*vo['C_LF'][iT,:,iEP['Bark']]+ \
        psl['bIPF_RootCoarseLitToLitterF']*vo['C_LF'][iT,:,iEP['RootCoarse']]+ \
        psl['bIPF_BranchMorToLitterF']*vo['C_M_Reg'][iT,:,iEP['Branch']]+ \
        psl['bIPF_StemNonMerchLitToLitterF']*vo['C_LF'][iT,:,iEP['StemNonMerch']]
            
    # Transfer biomass turnover to medium litter pool
    vo['C_Eco_Pools'][iT,:,iEP['LitterM']]=vo['C_Eco_Pools'][iT-1,:,iEP['LitterM']]+ \
        psl['bIPF_BarkLitToLitterM']*vo['C_LF'][iT,:,iEP['Bark']]+ \
        psl['bIPF_BarkMorToLitterM']*vo['C_M_Reg'][iT,:,iEP['Bark']]

    # Transfer biomass turnover to slow litter pool
    vo['C_Eco_Pools'][iT,:,iEP['LitterS']]=vo['C_Eco_Pools'][iT-1,:,iEP['LitterS']]+ \
        psl['bIPF_BarkLitToLitterS']*vo['C_LF'][iT,:,iEP['Bark']]+ \
        psl['bIPF_BarkMorToLitterS']*vo['C_M_Reg'][iT,:,iEP['Bark']]
        
    # Transfer biomass turnover to snag stemwood pool
    vo['C_Eco_Pools'][iT,:,iEP['SnagStem']]=vo['C_Eco_Pools'][iT-1,:,iEP['SnagStem']]+ \
        psl['bIPF_StemMerchMorToSnagStem']*vo['C_M_Reg'][iT,:,iEP['StemMerch']]+ \
        psl['bIPF_StemNonMerchMorToSnagStem']*vo['C_M_Reg'][iT,:,iEP['StemNonMerch']]
    
    # Transfer biomass turnover to snag branch pool
    vo['C_Eco_Pools'][iT,:,iEP['SnagBranch']]=vo['C_Eco_Pools'][iT-1,:,iEP['SnagBranch']]+ \
        psl['bIPF_BranchLitToSnagBranch']*vo['C_LF'][iT,:,iEP['Branch']]+ \
        psl['bIPF_BranchMorToSnagBranch']*vo['C_M_Reg'][iT,:,iEP['Branch']]+ \
        psl['bIPF_BarkMorToSnagBranch']*vo['C_M_Reg'][iT,:,iEP['Bark']]+ \
        psl['bIPF_StemNonMerchLitToSnagBranch']*vo['C_LF'][iT,:,iEP['Bark']]

    # Transfer biomass turnover to very fast soil pool
    vo['C_Eco_Pools'][iT,:,iEP['SoilVF']]=vo['C_Eco_Pools'][iT-1,:,iEP['SoilVF']]+ \
        psl['bIPF_RootFineLitToSoilVF']*vo['C_LF'][iT,:,iEP['RootFine']] +\
        psl['bIPF_RootCoarseMorToSoilVF']*vo['C_M_Reg'][iT,:,iEP['RootCoarse']]+ \
        psl['bIPF_RootFineMorToSoilVF']*vo['C_M_Reg'][iT,:,iEP['RootFine']]

    # Transfer biomass turnover to fast soil pool
    vo['C_Eco_Pools'][iT,:,iEP['SoilF']]=vo['C_Eco_Pools'][iT-1,:,iEP['SoilF']]+ \
        psl['bIPF_RootCoarseLitToSoilF']*vo['C_LF'][iT,:,iEP['RootCoarse']]+ \
        psl['bIPF_RootCoarseMorToSoilF']*vo['C_M_Reg'][iT,:,iEP['RootCoarse']]
        
    # Transfer biomass turnover to slow soil pool
    vo['C_Eco_Pools'][iT,:,iEP['SoilS']]=vo['C_Eco_Pools'][iT-1,:,iEP['SoilS']]    
    
    #--------------------------------------------------------------------------
    # Decomposition
    # Vectorizing the parts of this that can be vectorized was attempted, but
    # runtime actually increased a bit.
    #--------------------------------------------------------------------------
    
    # Prepare air temperature for respiration calculation
    Tref=10
    fT=(vi['Inv']['MAT']-Tref)/10
       
    # Respiration rate - Note that these terms do not equal the
    # atmosphere-bound efflux from heterotrophic respiration as a fraction is
    # emitted to the atmosphere and the remaining fraction is reorganized
    # within the ecosystem.
    meta['R_LitterVF']=psl['bDec_LitterVF_Rten']*vo['C_Eco_Pools'][iT,:,iEP['LitterVF']].flatten()*psl['bDec_LitterVF_Qten']**fT
    meta['R_LitterF']=psl['bDec_LitterF_Rten']*vo['C_Eco_Pools'][iT,:,iEP['LitterF']].flatten()*psl['bDec_LitterF_Qten']**fT
    meta['R_LitterM']=psl['bDec_LitterM_Rten']*vo['C_Eco_Pools'][iT,:,iEP['LitterM']].flatten()*psl['bDec_LitterM_Qten']**fT
    meta['R_LitterS']=psl['bDec_LitterS_Rten']*vo['C_Eco_Pools'][iT,:,iEP['LitterS']].flatten()*psl['bDec_LitterS_Qten']**fT
    meta['R_SnagStem']=psl['bDec_SnagStem_Rten']*vo['C_Eco_Pools'][iT,:,iEP['SnagStem']].flatten()*psl['bDec_SnagStem_Qten']**fT
    meta['R_SnagBranch']=psl['bDec_SnagBranch_Rten']*vo['C_Eco_Pools'][iT,:,iEP['SnagBranch']].flatten()*psl['bDec_SnagBranch_Qten']**fT    
    meta['R_SoilVF']=psl['bDec_SoilVF_Rten']*vo['C_Eco_Pools'][iT,:,iEP['SoilVF']].flatten()*psl['bDec_SoilVF_Qten']**fT
    meta['R_SoilF']=psl['bDec_SoilF_Rten']*vo['C_Eco_Pools'][iT,:,iEP['SoilF']].flatten()*psl['bDec_SoilF_Qten']**fT
    meta['R_SoilS']=psl['bDec_SoilS_Rten']*vo['C_Eco_Pools'][iT,:,iEP['SoilS']].flatten()*psl['bDec_SoilS_Qten']**fT
    
    # Adjust decomposition to account for N application response
    if meta['NM']['iApplication'].size>0:
        vi,vo,meta=update_nutrient_status(vi,vo,iT,meta,psl,'HeterotrophicRespiration')
    
    # Remove respired carbon from source DOM pools
    vo['C_Eco_Pools'][iT,:,iEP['LitterVF']]=vo['C_Eco_Pools'][iT,:,iEP['LitterVF']]-meta['R_LitterVF']
    vo['C_Eco_Pools'][iT,:,iEP['LitterF']]=vo['C_Eco_Pools'][iT,:,iEP['LitterF']]-meta['R_LitterF']
    vo['C_Eco_Pools'][iT,:,iEP['LitterM']]=vo['C_Eco_Pools'][iT,:,iEP['LitterM']]-meta['R_LitterM']
    vo['C_Eco_Pools'][iT,:,iEP['LitterS']]=vo['C_Eco_Pools'][iT,:,iEP['LitterS']]-meta['R_LitterS']
    vo['C_Eco_Pools'][iT,:,iEP['SnagStem']]=vo['C_Eco_Pools'][iT,:,iEP['SnagStem']]-meta['R_SnagStem']
    vo['C_Eco_Pools'][iT,:,iEP['SnagBranch']]=vo['C_Eco_Pools'][iT,:,iEP['SnagBranch']]-meta['R_SnagBranch']
    vo['C_Eco_Pools'][iT,:,iEP['SoilVF']]=vo['C_Eco_Pools'][iT,:,iEP['SoilVF']]-meta['R_SoilVF']
    vo['C_Eco_Pools'][iT,:,iEP['SoilF']]=vo['C_Eco_Pools'][iT,:,iEP['SoilF']]-meta['R_SoilF']
    vo['C_Eco_Pools'][iT,:,iEP['SoilS']]=vo['C_Eco_Pools'][iT,:,iEP['SoilS']]-meta['R_SoilS']
                
    # Re-define decayed fast litter
    vo['C_Eco_Pools'][iT,:,iEP['LitterF']]=vo['C_Eco_Pools'][iT,:,iEP['LitterF']]+psl['bIPF_SnagStemToLitterF']*meta['R_SnagStem']
    
    # Re-define decayed medium litter
    vo['C_Eco_Pools'][iT,:,iEP['LitterM']]=vo['C_Eco_Pools'][iT,:,iEP['LitterM']]+psl['bIPF_SnagBranchToLitterM']*meta['R_SnagBranch']
    
    # Re-define decayed slow litter
    vo['C_Eco_Pools'][iT,:,iEP['LitterS']]=vo['C_Eco_Pools'][iT,:,iEP['LitterS']]+psl['bIPF_LitterVFToLitterS']*meta['R_LitterVF']
    vo['C_Eco_Pools'][iT,:,iEP['LitterS']]=vo['C_Eco_Pools'][iT,:,iEP['LitterS']]+psl['bIPF_LitterFToLitterS']*meta['R_LitterF']
    vo['C_Eco_Pools'][iT,:,iEP['LitterS']]=vo['C_Eco_Pools'][iT,:,iEP['LitterS']]+psl['bIPF_LitterMToLitterS']*meta['R_LitterM']
    
    # Re-define decayed slow soil
    vo['C_Eco_Pools'][iT,:,iEP['SoilS']]=vo['C_Eco_Pools'][iT,:,iEP['SoilS']]+psl['bIPF_SoilVFToSoilS']*meta['R_SoilVF']
    vo['C_Eco_Pools'][iT,:,iEP['SoilS']]=vo['C_Eco_Pools'][iT,:,iEP['SoilS']]+psl['bIPF_SoilFToSoilS']*meta['R_SoilF']
        
    # Heterotrophic respiration 
    vo['C_RH'][iT,:,iEP['LitterVF']]=psl['bIPF_LitterVFToCO2']*meta['R_LitterVF']
    vo['C_RH'][iT,:,iEP['LitterF']]=psl['bIPF_LitterFToCO2']*meta['R_LitterF']
    vo['C_RH'][iT,:,iEP['LitterM']]=psl['bIPF_LitterMToCO2']*meta['R_LitterM']
    vo['C_RH'][iT,:,iEP['LitterS']]=psl['bIPF_LitterSToCO2']*meta['R_LitterS']
    vo['C_RH'][iT,:,iEP['SnagStem']]=psl['bIPF_SnagStemToCO2']*meta['R_SnagStem']
    vo['C_RH'][iT,:,iEP['SnagBranch']]=psl['bIPF_SnagBranchToCO2']*meta['R_SnagBranch']
    vo['C_RH'][iT,:,iEP['SoilVF']]=psl['bIPF_SoilVFToCO2']*meta['R_SoilVF']
    vo['C_RH'][iT,:,iEP['SoilF']]=psl['bIPF_SoilFToCO2']*meta['R_SoilF']
    vo['C_RH'][iT,:,iEP['SoilS']]=psl['bIPF_SoilSToCO2']*meta['R_SoilS']   
    
    #--------------------------------------------------------------------------
    # Physical transfer
    #--------------------------------------------------------------------------
    
    # Felled pools need to be initialized
    vo['C_Eco_Pools'][iT,:,iEP['FelledStemMerch']]=vo['C_Eco_Pools'][iT-1,:,iEP['FelledStemMerch']]
    vo['C_Eco_Pools'][iT,:,iEP['FelledStemNonMerch']]=vo['C_Eco_Pools'][iT-1,:,iEP['FelledStemNonMerch']]
    vo['C_Eco_Pools'][iT,:,iEP['FelledBranch']]=vo['C_Eco_Pools'][iT-1,:,iEP['FelledBranch']]
    vo['C_Eco_Pools'][iT,:,iEP['FelledBark']]=vo['C_Eco_Pools'][iT-1,:,iEP['FelledBark']]
    vo['C_Eco_Pools'][iT,:,iEP['FelledSnagStem']]=vo['C_Eco_Pools'][iT-1,:,iEP['FelledSnagStem']]
    vo['C_Eco_Pools'][iT,:,iEP['FelledSnagBranch']]=vo['C_Eco_Pools'][iT-1,:,iEP['FelledSnagBranch']]
    
    # Physical transfer rate
    PT_FelledStemMerch=psl['bDec_FelledStemMerch_PhysTransRate']*vo['C_Eco_Pools'][iT,:,iEP['FelledStemMerch']]
    PT_FelledStemNonMerch=psl['bDec_FelledStemNonMerch_PhysTransRate']*vo['C_Eco_Pools'][iT,:,iEP['FelledStemNonMerch']]
    PT_FelledBranch=psl['bDec_FelledBranch_PhysTransRate']*vo['C_Eco_Pools'][iT,:,iEP['FelledBranch']]
    PT_FelledBark=psl['bDec_FelledBark_PhysTransRate']*vo['C_Eco_Pools'][iT,:,iEP['FelledBark']]
    PT_FelledSnagStem=psl['bDec_FelledSnagStem_PhysTransRate']*vo['C_Eco_Pools'][iT,:,iEP['FelledSnagStem']]
    PT_FelledSnagBranch=psl['bDec_FelledSnagBranch_PhysTransRate']*vo['C_Eco_Pools'][iT,:,iEP['FelledSnagBranch']]
    
    PT_LitterS=psl['bDec_LitterS_PhysTransRate']*vo['C_Eco_Pools'][iT,:,iEP['LitterS']]
    PT_SnagStem=psl['bDec_SnagStem_PhysTransRate']*vo['C_Eco_Pools'][iT,:,iEP['SnagStem']]
    PT_SnagBranch=psl['bDec_SnagBranch_PhysTransRate']*vo['C_Eco_Pools'][iT,:,iEP['SnagBranch']] 
    
    # Remove carbon that is physically transferred
    vo['C_Eco_Pools'][iT,:,iEP['FelledStemMerch']]=vo['C_Eco_Pools'][iT,:,iEP['FelledStemMerch']]-PT_FelledStemMerch
    vo['C_Eco_Pools'][iT,:,iEP['FelledStemNonMerch']]=vo['C_Eco_Pools'][iT,:,iEP['FelledStemNonMerch']]-PT_FelledStemNonMerch
    vo['C_Eco_Pools'][iT,:,iEP['FelledBranch']]=vo['C_Eco_Pools'][iT,:,iEP['FelledBranch']]-PT_FelledBranch
    vo['C_Eco_Pools'][iT,:,iEP['FelledBark']]=vo['C_Eco_Pools'][iT,:,iEP['FelledBark']]-PT_FelledBark
    vo['C_Eco_Pools'][iT,:,iEP['FelledSnagStem']]=vo['C_Eco_Pools'][iT,:,iEP['FelledSnagStem']]-PT_FelledSnagStem
    vo['C_Eco_Pools'][iT,:,iEP['FelledSnagBranch']]=vo['C_Eco_Pools'][iT,:,iEP['FelledSnagBranch']]-PT_FelledSnagBranch    
    
    vo['C_Eco_Pools'][iT,:,iEP['LitterS']]=vo['C_Eco_Pools'][iT,:,iEP['LitterS']]-PT_LitterS
    vo['C_Eco_Pools'][iT,:,iEP['SnagStem']]=vo['C_Eco_Pools'][iT,:,iEP['SnagStem']]-PT_SnagStem
    vo['C_Eco_Pools'][iT,:,iEP['SnagBranch']]=vo['C_Eco_Pools'][iT,:,iEP['SnagBranch']]-PT_SnagBranch
    
    # Add carbon transferred from felled pools to DOM pools
    vo['C_Eco_Pools'][iT,:,iEP['LitterM']]=vo['C_Eco_Pools'][iT,:,iEP['LitterM']]+PT_FelledStemMerch
    vo['C_Eco_Pools'][iT,:,iEP['LitterM']]=vo['C_Eco_Pools'][iT,:,iEP['LitterM']]+PT_FelledStemNonMerch
    vo['C_Eco_Pools'][iT,:,iEP['LitterF']]=vo['C_Eco_Pools'][iT,:,iEP['LitterF']]+PT_FelledBranch
    vo['C_Eco_Pools'][iT,:,iEP['LitterF']]=vo['C_Eco_Pools'][iT,:,iEP['LitterF']]+PT_FelledBark
    vo['C_Eco_Pools'][iT,:,iEP['LitterM']]=vo['C_Eco_Pools'][iT,:,iEP['LitterM']]+PT_FelledSnagStem
    vo['C_Eco_Pools'][iT,:,iEP['LitterM']]=vo['C_Eco_Pools'][iT,:,iEP['LitterM']]+PT_FelledSnagBranch
    
    # Add decomposed carbon to more decomposed DOM pools
    vo['C_Eco_Pools'][iT,:,iEP['SoilS']]=vo['C_Eco_Pools'][iT,:,iEP['SoilS']]+PT_LitterS
    vo['C_Eco_Pools'][iT,:,iEP['LitterF']]=vo['C_Eco_Pools'][iT,:,iEP['LitterF']]+PT_SnagBranch    
    vo['C_Eco_Pools'][iT,:,iEP['LitterM']]=vo['C_Eco_Pools'][iT,:,iEP['LitterM']]+PT_SnagStem
    
    #--------------------------------------------------------------------------
    # litter decomposition
    # *** Added for fertilization study ***
    #--------------------------------------------------------------------------
    
    #vo['C_Eco_Pools'][iT,:,28]=meta['R_LitterVF']+meta['R_LitterF']+meta['R_LitterM']+meta['R_LitterS']
    
    return vo

#%% Disturbance and management events

def DisturbanceAndManagementEvents(iT,vi,vo,psl,meta,iEP):
    
    # Predict stand breakup (on the fly)
    if meta['Simulate breakup on the fly']=='On':
        vi=gensm.PredictStandBreakup_OnTheFly(meta,vi,iT,vo['A'][iT,:])    
    
    # Predict harvesting (on the fly)
    if meta['Simulate harvesting on the fly']=='On':            
        if meta['Dist on Fly']['Harvesting'][meta['iScn']]==1:
            vi=gensm.PredictHarvesting_OnTheFly(meta,vi,iT,vi['Inv']['THLB'][iT,:],vo['V_StemMerch'][iT,:])
    
    # Initialize indicator of aerial nutrient application
    flag_nutrient_application=np.zeros(meta['N Stand'])
    
    # Loop through events in year
    for iE in range(meta['Max Events Per Year']):
        
        # Event type IDs for the iE'th event of the year
        ID_Type=vi['EC']['ID_Type'][iT,:,iE].copy()
        
        # Record stands with aerial nutrient application
        ind=np.where(ID_Type==meta['LUT Dist']['Fertilization Aerial'])[0]
        flag_nutrient_application[ind]=1
        
        # Total affected biomass carbon
        MortalityFactor=vi['EC']['MortalityFactor'][iT,:,iE].copy()
        
        # Only continue if there are disturbance fluxes
        if np.sum(MortalityFactor)==0:
            continue
        
        # Get event-specific parameters        
        u,idx,inv=np.unique(ID_Type,return_index=True,return_inverse=True)
        
        b={}
        for k in psl['Dist'][1].keys():
            bU=np.zeros(u.size)
            for iU in range(u.size):
                if u[iU]==0:
                    continue
                bU[iU]=psl['Dist'][u[iU]][k]
            b[k]=bU[inv]
            
        #----------------------------------------------------------------------
        # Define the amount of each pool that is affected by the event
        #----------------------------------------------------------------------

        BiomassMerchAffectedFrac=b['BiomassMerch_Affected']*MortalityFactor
        BiomassNonMerchAffectedFrac=b['BiomassNonMerch_Affected']*MortalityFactor
        SnagsAffectedFrac=b['Snags_Affected']*MortalityFactor

        # Biomass carbon
        Affected_StemMerch=BiomassMerchAffectedFrac*vo['C_Eco_Pools'][iT,:,iEP['StemMerch']]
        Affected_StemNonMerch=BiomassNonMerchAffectedFrac*vo['C_Eco_Pools'][iT,:,iEP['StemNonMerch']]
        Affected_Foliage=BiomassNonMerchAffectedFrac*vo['C_Eco_Pools'][iT,:,iEP['Foliage']]
        Affected_Branch=BiomassNonMerchAffectedFrac*vo['C_Eco_Pools'][iT,:,iEP['Branch']]
        Affected_Bark=BiomassNonMerchAffectedFrac*vo['C_Eco_Pools'][iT,:,iEP['Bark']]            
        Affected_RootCoarse=BiomassNonMerchAffectedFrac*vo['C_Eco_Pools'][iT,:,iEP['RootCoarse']]
        Affected_RootFine=BiomassNonMerchAffectedFrac*vo['C_Eco_Pools'][iT,:,iEP['RootFine']]
                        
        # Partition bark into merch and non-merch components
        Affected_BarkMerch=psl['bASL_MerchBarkFrac']*Affected_Bark
        Affected_BarkNonMerch=(1-psl['bASL_MerchBarkFrac'])*Affected_Bark
            
        # Sum up total affected non-merchantable biomass
        Affected_TotNonMerch=Affected_StemNonMerch+Affected_Branch+Affected_BarkNonMerch
            
        # Snags           
        Affected_SnagStem=SnagsAffectedFrac*vo['C_Eco_Pools'][iT,:,iEP['SnagStem']]
        Affected_SnagBranch=SnagsAffectedFrac*vo['C_Eco_Pools'][iT,:,iEP['SnagBranch']]
            
        # All carbon
        Affected_All=Affected_StemMerch+Affected_StemNonMerch+Affected_Foliage+Affected_Branch+Affected_Bark+Affected_RootCoarse+Affected_RootFine
        
        # Stemwood merch. volume
        Affected_VolumeStemMerch=BiomassMerchAffectedFrac*vo['V_StemMerch'][iT,:]
        
        #----------------------------------------------------------------------
        # Calculate mortality
        #----------------------------------------------------------------------
            
        vo['C_M_Dist'][iT,:]=vo['C_M_Dist'][iT,:]+Affected_All
            
        #----------------------------------------------------------------------
        # Remove affected amount from each pool
        #----------------------------------------------------------------------
            
        # Remove carbon from affected biomass pools            
        vo['C_Eco_Pools'][iT,:,iEP['StemMerch']]=vo['C_Eco_Pools'][iT,:,iEP['StemMerch']]-Affected_StemMerch
        vo['C_Eco_Pools'][iT,:,iEP['StemNonMerch']]=vo['C_Eco_Pools'][iT,:,iEP['StemNonMerch']]-Affected_StemNonMerch
        vo['C_Eco_Pools'][iT,:,iEP['Foliage']]=vo['C_Eco_Pools'][iT,:,iEP['Foliage']]-Affected_Foliage
        vo['C_Eco_Pools'][iT,:,iEP['Branch']]=vo['C_Eco_Pools'][iT,:,iEP['Branch']]-Affected_Branch
        vo['C_Eco_Pools'][iT,:,iEP['Bark']]=vo['C_Eco_Pools'][iT,:,iEP['Bark']]-Affected_Bark
        vo['C_Eco_Pools'][iT,:,iEP['RootCoarse']]=vo['C_Eco_Pools'][iT,:,iEP['RootCoarse']]-Affected_RootCoarse
        vo['C_Eco_Pools'][iT,:,iEP['RootFine']]=vo['C_Eco_Pools'][iT,:,iEP['RootFine']]-Affected_RootFine
               
        # Remove carbon from snag pools
        vo['C_Eco_Pools'][iT,:,iEP['SnagStem']]=vo['C_Eco_Pools'][iT,:,iEP['SnagStem']]-Affected_SnagStem
        vo['C_Eco_Pools'][iT,:,iEP['SnagBranch']]=vo['C_Eco_Pools'][iT,:,iEP['SnagBranch']]-Affected_SnagBranch
            
        # Remove stemwood merch volume
        vo['V_StemMerch'][iT,:]=vo['V_StemMerch'][iT,:]-Affected_VolumeStemMerch
            
        #----------------------------------------------------------------------
        # Carbon removed (sent to mill)
        #----------------------------------------------------------------------
            
        # Merch biomass to mill - of the total amount of biomass affected, 
        # Add bark
        Biomass_Merch_Removed=b['BiomassMerch_Removed']*(Affected_StemMerch+Affected_BarkMerch)
        vo['C_RemovedMerch'][iT,:]=vo['C_RemovedMerch'][iT,:]+Biomass_Merch_Removed
        
        # NonMerch biomass to mill - of the total amount of biomass affected, 
        # what fraction of non-merch biomass was sent to the mill?
        # - NonMerch = NonMerchStem + Foliage + Branch + Bark
        Biomass_NonMerch_Removed=b['BiomassNonMerch_Removed']*Affected_TotNonMerch
        vo['C_RemovedNonMerch'][iT,:]=vo['C_RemovedNonMerch'][iT,:]+Biomass_NonMerch_Removed
        
        # Snag stemwood to mill           
        SnagStem_Removed=b['Snags_Removed']*Affected_SnagStem
        vo['C_RemovedSnagStem'][iT,:]=vo['C_RemovedSnagStem'][iT,:]+SnagStem_Removed
        
        # Snag branches to mill
        SnagBranch_Removed=b['Snags_Removed']*Affected_SnagBranch
        vo['C_RemovedNonMerch'][iT,:]=vo['C_RemovedNonMerch'][iT,:]+SnagBranch_Removed
        
        #----------------------------------------------------------------------
        # Carbon that is left on site (after felling or wind storms)
        #----------------------------------------------------------------------
        
        # Stem, branch and bark carbon transfered to felled pools (piles)
        vo['C_Eco_Pools'][iT,:,iEP['FelledStemMerch']]=vo['C_Eco_Pools'][iT,:,iEP['FelledStemMerch']]+b['BiomassMerch_LeftOnSite']*Affected_StemMerch
        vo['C_Eco_Pools'][iT,:,iEP['FelledStemNonMerch']]=vo['C_Eco_Pools'][iT,:,iEP['FelledStemNonMerch']]+b['BiomassNonMerch_LeftOnSite']*Affected_StemNonMerch            
        vo['C_Eco_Pools'][iT,:,iEP['FelledBranch']]=vo['C_Eco_Pools'][iT,:,iEP['FelledBranch']]+b['BiomassMerch_LeftOnSite']*Affected_Branch
        vo['C_Eco_Pools'][iT,:,iEP['FelledBark']]=vo['C_Eco_Pools'][iT,:,iEP['FelledBark']]+b['BiomassMerch_LeftOnSite']*Affected_Bark
        
        # Foliage transferred directly to very fast litter
        vo['C_Eco_Pools'][iT,:,iEP['LitterVF']]=vo['C_Eco_Pools'][iT,:,iEP['LitterVF']]+Affected_Foliage
        
        # Roots transferred directly to DOM       
        vo['C_Eco_Pools'][iT,:,iEP['LitterF']]=vo['C_Eco_Pools'][iT,:,iEP['LitterF']]+0.5*Affected_RootCoarse
        vo['C_Eco_Pools'][iT,:,iEP['SoilF']]=vo['C_Eco_Pools'][iT,:,iEP['SoilF']]+0.5*Affected_RootCoarse            
        vo['C_Eco_Pools'][iT,:,iEP['LitterVF']]=vo['C_Eco_Pools'][iT,:,iEP['LitterVF']]+0.5*Affected_RootFine
        vo['C_Eco_Pools'][iT,:,iEP['SoilVF']]=vo['C_Eco_Pools'][iT,:,iEP['SoilVF']]+0.5*Affected_RootFine
        
        # Snags transferred to felled pools (piles)
        vo['C_Eco_Pools'][iT,:,iEP['FelledSnagStem']]=vo['C_Eco_Pools'][iT,:,iEP['FelledSnagStem']]+b['Snags_LeftOnSite']*Affected_SnagStem
        vo['C_Eco_Pools'][iT,:,iEP['FelledSnagBranch']]=vo['C_Eco_Pools'][iT,:,iEP['FelledSnagBranch']]+b['Snags_LeftOnSite']*Affected_SnagBranch
        
        #----------------------------------------------------------------------
        # Biomass and DOM burned
        #----------------------------------------------------------------------
        
        Biomass_Merch_Burned=b['BiomassMerch_Burned']*(Affected_StemMerch+Affected_BarkMerch)
        Biomass_NonMerch_Burned=b['BiomassNonMerch_Burned']*Affected_TotNonMerch
        SnagStem_Burned=b['Snags_Burned']*Affected_SnagStem
        SnagBranch_Burned=b['Snags_Burned']*Affected_SnagBranch
        Total_Burned=Biomass_Merch_Burned+Biomass_NonMerch_Burned+SnagStem_Burned+SnagBranch_Burned
        
        vo['C_E_FireAsCO2'][iT,:]=vo['C_E_FireAsCO2'][iT,:]+psl['bCombFrac_CO2']*Total_Burned        
        vo['C_E_FireAsCH4'][iT,:]=vo['C_E_FireAsCH4'][iT,:]+psl['bCombFrac_CH4']*Total_Burned           
        vo['C_E_FireAsCO'][iT,:]=vo['C_E_FireAsCO'][iT,:]+psl['bCombFrac_CO']*Total_Burned        
        
        #----------------------------------------------------------------------
        # Carbon that is moved from biomass to snags
        #----------------------------------------------------------------------
            
        # Merch biomass that is killed by wildfire
        vo['C_Eco_Pools'][iT,:,iEP['SnagStem']]=vo['C_Eco_Pools'][iT,:,iEP['SnagStem']]+b['BiomassMerch_ToSnag']*(Affected_StemMerch+Affected_BarkMerch)
        
        # Non-merch biomass that is killed by wildfire
        vo['C_Eco_Pools'][iT,:,iEP['SnagBranch']]=vo['C_Eco_Pools'][iT,:,iEP['SnagBranch']]+b['BiomassNonMerch_ToSnag']*(Affected_StemNonMerch+Affected_BarkNonMerch)
        
        #----------------------------------------------------------------------
        # Felled carbon that is burned
        #----------------------------------------------------------------------           
            
        # Affected carbon
        StemMerch=b['FelledStemMerch_Burned']*vo['C_Eco_Pools'][iT,:,iEP['FelledStemMerch']]
        StemNonMerch=b['FelledStemNonMerch_Burned']*vo['C_Eco_Pools'][iT,:,iEP['FelledStemNonMerch']]
        Branch=b['FelledBranch_Burned']*vo['C_Eco_Pools'][iT,:,iEP['FelledBranch']]
        Bark=b['FelledBark_Burned']*vo['C_Eco_Pools'][iT,:,iEP['FelledBark']]
        SnagStem=b['FelledSnagStem_Burned']*vo['C_Eco_Pools'][iT,:,iEP['FelledSnagStem']]
        SnagBranch=b['FelledSnagStem_Burned']*vo['C_Eco_Pools'][iT,:,iEP['FelledSnagBranch']]
        
        # Remove affected carbon
        vo['C_Eco_Pools'][iT,:,iEP['FelledStemMerch']]=vo['C_Eco_Pools'][iT,:,iEP['FelledStemMerch']]-StemMerch
        vo['C_Eco_Pools'][iT,:,iEP['FelledStemNonMerch']]=vo['C_Eco_Pools'][iT,:,iEP['FelledStemNonMerch']]-StemNonMerch
        vo['C_Eco_Pools'][iT,:,iEP['FelledBranch']]=vo['C_Eco_Pools'][iT,:,iEP['FelledBranch']]-Branch
        vo['C_Eco_Pools'][iT,:,iEP['FelledBark']]=vo['C_Eco_Pools'][iT,:,iEP['FelledBark']]-Bark
        vo['C_Eco_Pools'][iT,:,iEP['FelledSnagStem']]=vo['C_Eco_Pools'][iT,:,iEP['FelledSnagStem']]-SnagStem
        vo['C_Eco_Pools'][iT,:,iEP['FelledSnagBranch']]=vo['C_Eco_Pools'][iT,:,iEP['FelledSnagBranch']]-SnagBranch
        
        Total_Burned=StemMerch+StemNonMerch+Branch+Bark+SnagStem+SnagBranch
        
        # Add to fire emissions
        vo['C_E_FireAsCO2'][iT,:]=vo['C_E_FireAsCO2'][iT,:]+psl['bCombFrac_CO2']*Total_Burned        
        vo['C_E_FireAsCH4'][iT,:]=vo['C_E_FireAsCH4'][iT,:]+psl['bCombFrac_CH4']*Total_Burned           
        vo['C_E_FireAsCO'][iT,:]=vo['C_E_FireAsCO'][iT,:]+psl['bCombFrac_CO']*Total_Burned
        
        #----------------------------------------------------------------------
        # Update stand age
        #----------------------------------------------------------------------
            
        if meta['Biomass Module']!='Sawtooth':
                    
            if (meta['Partial Mortality Affects Age']=='On'):
                        
                # Assume oldest trees were most affected, reduce age in prportion
                # with mortality rate
                # Not always realistic, but see how net growth is affected.
                vo['A'][iT,:]=vo['A'][iT,:]*(1-MortalityFactor)
                        
            else:                   
                        
                # Assume no change
                vo['A'][iT,:]=vo['A'][iT,:]
                
        # Ensure planting resets to age 0                
        vo['A'][iT,(ID_Type==meta['LUT Dist']['Planting']) | (ID_Type==meta['LUT Dist']['Direct Seeding'])]=0
                    
        #----------------------------------------------------------------------
        # Transition to new growth curve
        #----------------------------------------------------------------------
            
        # Initialize a flag that indicates whether the growth curve changes
        #flg_gc_change=0
        #if vi['EH'][iS]['ID_GrowthCurveM'][indDist[iDist]]!=vi['ID_GCA'][iS]:
        #    flg_gc_change=1
            
        # Only applies to BatchTIPSY
        if meta['Biomass Module']=='BatchTIPSY':
            
            for iGC in range(meta['GC']['N Growth Curves']):
                ind=np.where(vi['EC']['ID_GrowthCurve'][iT,:,iE]==meta['GC']['ID GC'][iGC])[0]
                if ind.size>0:
                    vi['GC']['Active'][:,ind,:]=vi['GC'][meta['GC']['ID GC'][iGC]][:,ind,:]
                    vi['GC']['ID_GCA'][ind]=int(meta['GC']['ID GC'][iGC])
    
    #--------------------------------------------------------------------------
    # Aerial nutrient application events
    #--------------------------------------------------------------------------
    
    # Generate index to stands that were fertilized
    meta['NM']['iApplication']=np.where(flag_nutrient_application==1)[0]    
        
    if meta['NM']['iApplication'].size>0:        

        # Adjust net growth of aboveground biomass
        vi,vo,meta=update_nutrient_status(vi,vo,iT,meta,psl,'AbovegroundNetGrowth')
                
        # Adjust emissions
        vi,vo,meta=update_nutrient_status(vi,vo,iT,meta,psl,'Emissions')
        
    return vo,vi

#%% Harvested wood products sector (from Dymond 2012)

def HWP_From_Dymond12(iT,vi,vo,psl,meta):
    
    #--------------------------------------------------------------------------
    # If custom harvests occur, revise parameters
    # *** This is going to be slow, but thus far, custom harvests have only been
    # used in small projects. It will require reworking if applied in big 
    # projects. ***
    #--------------------------------------------------------------------------
    
    if 'Harvest Custom' in meta:
        
        for iE in range(meta['Max Events Per Year']):
            
            ID_Type=vi['EC']['ID_Type'][iT,:,iE]
        
            for iHC in range(10):
                
                ID_HC=int(iHC+1)
                ID_Dist_HC=meta['LUT Dist']['Harvest Custom ' + str(ID_HC)]
                
                if ID_HC in meta['Harvest Custom']:
        
                    iStands=np.where(ID_Type==ID_Dist_HC)[0]
                    
                    if iStands.size==0:
                        continue
                    
                    for k in meta['Harvest Custom'][ID_HC].keys():
                        
                        if k[0:7]=='Removed':
                            psl['HWP'][k][iStands]=meta['Harvest Custom'][ID_HC][k]/100
    
    #--------------------------------------------------------------------------
    # Ecosystems --> Mills
    # Note: Mills act as transient reservoirs that are not tracked over time
    #--------------------------------------------------------------------------
    
    Mill_FuelFromChips=psl['HWP']['RemovedMerchToFuel']*vo['C_RemovedMerch'][iT,:] + \
                 psl['HWP']['RemovedNonMerchToFuel']*vo['C_RemovedNonMerch'][iT,:] + \
                 psl['HWP']['RemovedSnagStemToFuel']*vo['C_RemovedSnagStem'][iT,:]    
    Mill_PulpFromChips=psl['HWP']['RemovedMerchToPulp']*vo['C_RemovedMerch'][iT,:] + \
                 psl['HWP']['RemovedNonMerchToPulp']*vo['C_RemovedNonMerch'][iT,:] + \
                 psl['HWP']['RemovedSnagStemToPulp']*vo['C_RemovedSnagStem'][iT,:]    
    Mill_Lumber=psl['HWP']['RemovedMerchToLumber']*vo['C_RemovedMerch'][iT,:] + \
                 psl['HWP']['RemovedNonMerchToLumber']*vo['C_RemovedNonMerch'][iT,:] + \
                 psl['HWP']['RemovedSnagStemToLumber']*vo['C_RemovedSnagStem'][iT,:]
    Mill_Plywood=psl['HWP']['RemovedMerchToPlywood']*vo['C_RemovedMerch'][iT,:] + \
                 psl['HWP']['RemovedNonMerchToPlywood']*vo['C_RemovedNonMerch'][iT,:] + \
                 psl['HWP']['RemovedSnagStemToPlywood']*vo['C_RemovedSnagStem'][iT,:]
    Mill_OSB=psl['HWP']['RemovedMerchToOSB']*vo['C_RemovedMerch'][iT,:] + \
                 psl['HWP']['RemovedNonMerchToOSB']*vo['C_RemovedNonMerch'][iT,:] + \
                 psl['HWP']['RemovedSnagStemToOSB']*vo['C_RemovedSnagStem'][iT,:]
    Mill_MDF=psl['HWP']['RemovedMerchToMDF']*vo['C_RemovedMerch'][iT,:] + \
                 psl['HWP']['RemovedNonMerchToMDF']*vo['C_RemovedNonMerch'][iT,:] + \
                 psl['HWP']['RemovedSnagStemToMDF']*vo['C_RemovedSnagStem'][iT,:]
    Mill_Firewood=psl['HWP']['RemovedMerchToFirewood']*vo['C_RemovedMerch'][iT,:] + \
                 psl['HWP']['RemovedNonMerchToFirewood']*vo['C_RemovedNonMerch'][iT,:] + \
                 psl['HWP']['RemovedSnagStemToFirewood']*vo['C_RemovedSnagStem'][iT,:]
    Mill_Cants=psl['HWP']['RemovedMerchToCants']*vo['C_RemovedMerch'][iT,:] + \
                 psl['HWP']['RemovedNonMerchToCants']*vo['C_RemovedNonMerch'][iT,:] + \
                 psl['HWP']['RemovedSnagStemToCants']*vo['C_RemovedSnagStem'][iT,:]
     
    #--------------------------------------------------------------------------
    # Mills --> In-use pools or other mills
    #--------------------------------------------------------------------------
        
    # Transfer mill fibre to single-family homes
    ip=meta['iPP']['SFH']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT-1,:,ip] + \
        psl['HWP']['LumberToSFH']*Mill_Lumber + \
        psl['HWP']['PlywoodToSFH']*Mill_Plywood + \
        psl['HWP']['OSBToSFH']*Mill_OSB + \
        psl['HWP']['MDFToSFH']*Mill_MDF
    
    # Transfer mill fibre to multi-family homes
    ip=meta['iPP']['MFH']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT-1,:,ip] + \
        psl['HWP']['LumberToMFH']*Mill_Lumber + \
        psl['HWP']['PlywoodToMFH']*Mill_Plywood + \
        psl['HWP']['OSBToMFH']*Mill_OSB + \
        psl['HWP']['MDFToMFH']*Mill_MDF
    
    # Transfer mill fibre to commercial
    ip=meta['iPP']['Comm']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT-1,:,ip] + \
        psl['HWP']['LumberToCom']*Mill_Lumber + \
        psl['HWP']['PlywoodToCom']*Mill_Plywood + \
        psl['HWP']['OSBToCom']*Mill_OSB + \
        psl['HWP']['MDFToCom']*Mill_MDF
    
    # Transfer mill fibre to furniture
    ip=meta['iPP']['Furn']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT-1,:,ip] + \
        psl['HWP']['LumberToFurn']*Mill_Lumber + \
        psl['HWP']['PlywoodToFurn']*Mill_Plywood + \
        psl['HWP']['OSBToFurn']*Mill_OSB + \
        psl['HWP']['MDFToFurn']*Mill_MDF
    
    # Transfer mill fibre to shipping
    ip=meta['iPP']['Ship']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT-1,:,ip] + \
        psl['HWP']['LumberToShip']*Mill_Lumber + \
        psl['HWP']['PlywoodToShip']*Mill_Plywood + \
        psl['HWP']['OSBToShip']*Mill_OSB + \
        psl['HWP']['MDFToShip']*Mill_MDF
    
    # Transfer mill fibre to repairs
    ip=meta['iPP']['Repairs']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT-1,:,ip] + \
        psl['HWP']['LumberToRepairs']*Mill_Lumber + \
        psl['HWP']['PlywoodToRepairs']*Mill_Plywood + \
        psl['HWP']['OSBToRepairs']*Mill_OSB + \
        psl['HWP']['MDFToRepairs']*Mill_MDF
    
    # Transfer mill fibre to other
    ip=meta['iPP']['Other']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT-1,:,ip] + \
        psl['HWP']['LumberToOther']*Mill_Lumber + \
        psl['HWP']['PlywoodToOther']*Mill_Plywood + \
        psl['HWP']['OSBToOther']*Mill_OSB + \
        psl['HWP']['MDFToOther']*Mill_MDF
    
    # Transfer mill fibre to pulp mill
    Mill_Pulp=Mill_PulpFromChips + \
                 psl['HWP']['LumberToPulp']*Mill_Lumber + \
                 psl['HWP']['PlywoodToPulp']*Mill_Plywood
    
    # Transfer pulp mill fibre to paper
    ip=meta['iPP']['Paper']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT-1,:,ip] + \
                 psl['HWP']['PulpToPaper']*Mill_Pulp
    
    # Transfer mill fibre to fuel 
    ip=meta['iPP']['Fuel']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT-1,:,ip] + \
                 Mill_FuelFromChips + \
                 psl['HWP']['LumberToFuel']*Mill_Lumber + \
                 psl['HWP']['PlywoodToFuel']*Mill_Plywood + \
                 psl['HWP']['OSBToFuel']*Mill_OSB + \
                 psl['HWP']['MDFToFuel']*Mill_MDF
    
    # Transfer firewood to firewood pool
    ip=meta['iPP']['Firewood']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT-1,:,ip] + Mill_Firewood
    
    # Transfer mill fibre to cants
    ip=meta['iPP']['Cants']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT-1,:,ip] + Mill_Cants
    
    # Transfer pulp mill carbon to pulp-mill effluent
    ip=meta['iPP']['EffluentPulp']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT-1,:,ip] + \
                 psl['HWP']['PulpToEffluent']*Mill_Pulp   
    
    #--------------------------------------------------------------------------
    # Update dump and landfill reservoirs
    #--------------------------------------------------------------------------

    vo['C_Pro_Pools'][iT,:,meta['iPP']['DumpWood']]=vo['C_Pro_Pools'][iT-1,:,meta['iPP']['DumpWood']]
    vo['C_Pro_Pools'][iT,:,meta['iPP']['DumpPaper']]=vo['C_Pro_Pools'][iT-1,:,meta['iPP']['DumpPaper']]
    vo['C_Pro_Pools'][iT,:,meta['iPP']['LandfillWoodDegradable']]=vo['C_Pro_Pools'][iT-1,:,meta['iPP']['LandfillWoodDegradable']]
    vo['C_Pro_Pools'][iT,:,meta['iPP']['LandfillWoodNonDegradable']]=vo['C_Pro_Pools'][iT-1,:,meta['iPP']['LandfillWoodNonDegradable']]    
    vo['C_Pro_Pools'][iT,:,meta['iPP']['LandfillPaperDegradable']]=vo['C_Pro_Pools'][iT-1,:,meta['iPP']['LandfillPaperDegradable']]
    vo['C_Pro_Pools'][iT,:,meta['iPP']['LandfillPaperNonDegradable']]=vo['C_Pro_Pools'][iT-1,:,meta['iPP']['LandfillPaperNonDegradable']]    
        
    #--------------------------------------------------------------------------
    # Single-family homes --> dump and landfill
    #--------------------------------------------------------------------------
           
    # Turnover
    ip=meta['iPP']['SFH']
    C_retired=psl['HWP']['SFH_tr']*vo['C_Pro_Pools'][iT-1,:,ip]
    
    # Remove carbon
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] - C_retired
    
    # Transfer carbon to dump wood
    ip=meta['iPP']['DumpWood']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + psl['HWP']['SFHToDumpWood']*C_retired
    
    # Transfer carbon to landfill (degradable)
    ip=meta['iPP']['LandfillWoodDegradable']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + psl['HWP']['SFHToLandfillWood']*psl['HWP']['ToLandfillWoodDegradableFrac']*C_retired
    
    # Transfer carbon to landfill (non-degradable)
    ip=meta['iPP']['LandfillWoodNonDegradable']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + psl['HWP']['SFHToLandfillWood']*(1-psl['HWP']['ToLandfillWoodDegradableFrac'])*C_retired
    
    #--------------------------------------------------------------------------
    # Multi-family homes --> dump and landfill
    #--------------------------------------------------------------------------
    
    # Turnover
    ip=meta['iPP']['MFH']
    C_retired=psl['HWP']['MFH_tr']*vo['C_Pro_Pools'][iT,:,ip]
    
    # Remove carbon
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] - C_retired
    
    # Transfer carbon to dump wood
    ip=meta['iPP']['DumpWood']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + psl['HWP']['MFHToDumpWood']*C_retired
    
    # Transfer carbon to landfill (degradable)
    ip=meta['iPP']['LandfillWoodDegradable']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + psl['HWP']['MFHToLandfillWood']*psl['HWP']['ToLandfillWoodDegradableFrac']*C_retired
    
    # Transfer carbon to landfill (non-degradable)
    ip=meta['iPP']['LandfillWoodNonDegradable']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + psl['HWP']['MFHToLandfillWood']*(1-psl['HWP']['ToLandfillWoodDegradableFrac'])*C_retired
    
    #--------------------------------------------------------------------------
    # Commercial building --> dump and landfill
    #--------------------------------------------------------------------------
    
    # Turnover
    ip=meta['iPP']['Comm']
    C_retired=psl['HWP']['Comm_tr']*vo['C_Pro_Pools'][iT-1,:,ip]
    
    # Remove carbon
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] - C_retired
    
    # Transfer carbon to dump wood
    ip=meta['iPP']['DumpWood']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + psl['HWP']['CommToDumpWood']*C_retired
    
    # Transfer carbon to landfill (degradable)
    ip=meta['iPP']['LandfillWoodDegradable']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + psl['HWP']['CommToLandfillWood']*psl['HWP']['ToLandfillWoodDegradableFrac']*C_retired
    
    # Transfer carbon to landfill (non-degradable)
    ip=meta['iPP']['LandfillWoodNonDegradable']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + psl['HWP']['CommToLandfillWood']*(1-psl['HWP']['ToLandfillWoodDegradableFrac'])*C_retired
    
    #--------------------------------------------------------------------------
    # Furniture --> dump and landfill
    #--------------------------------------------------------------------------
    
    # Turnover
    ip=meta['iPP']['Furn']
    C_retired=psl['HWP']['Furn_tr']*vo['C_Pro_Pools'][iT-1,:,ip]
    
    # Remove carbon
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] - C_retired
    
    # Transfer carbon to dump wood
    ip=meta['iPP']['DumpWood']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + psl['HWP']['FurnToDumpWood']*C_retired
    
    # Transfer carbon to landfill (degradable)
    ip=meta['iPP']['LandfillWoodDegradable']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + psl['HWP']['FurnToLandfillWood']*psl['HWP']['ToLandfillWoodDegradableFrac']*C_retired
    
    # Transfer carbon to landfill (non-degradable)
    ip=meta['iPP']['LandfillWoodNonDegradable']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + psl['HWP']['FurnToLandfillWood']*(1-psl['HWP']['ToLandfillWoodDegradableFrac'])*C_retired
    
    
    #--------------------------------------------------------------------------
    # Shipping --> dump and landfill
    #--------------------------------------------------------------------------
    
    # Turnover
    ip=meta['iPP']['Ship']
    C_retired=psl['HWP']['Ship_tr']*vo['C_Pro_Pools'][iT-1,:,ip]
    
    # Remove carbon
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] - C_retired
    
    # Transfer carbon to dump wood
    ip=meta['iPP']['DumpWood']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + psl['HWP']['ShipToDumpWood']*C_retired
    
    # Transfer carbon to landfill (degradable)
    ip=meta['iPP']['LandfillWoodDegradable']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + psl['HWP']['ShipToLandfillWood']*psl['HWP']['ToLandfillWoodDegradableFrac']*C_retired
    
    # Transfer carbon to landfill (non-degradable)
    ip=meta['iPP']['LandfillWoodNonDegradable']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + psl['HWP']['ShipToLandfillWood']*(1-psl['HWP']['ToLandfillWoodDegradableFrac'])*C_retired
    
    #--------------------------------------------------------------------------
    # Repairs --> dump and landfill
    #--------------------------------------------------------------------------
    
    # Turnover
    ip=meta['iPP']['Repairs']
    C_retired=psl['HWP']['Repairs_tr']*vo['C_Pro_Pools'][iT-1,:,ip]
    
    # Remove carbon
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] - C_retired
    
    # Transfer carbon to dump wood
    ip=meta['iPP']['DumpWood']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + psl['HWP']['RepairsToDumpWood']*C_retired
    
    # Transfer carbon to landfill (degradble)
    ip=meta['iPP']['LandfillWoodDegradable']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + psl['HWP']['RepairsToLandfillWood']*psl['HWP']['ToLandfillWoodDegradableFrac']*C_retired
    
    # Transfer carbon to landfill (non-degradable)
    ip=meta['iPP']['LandfillWoodNonDegradable']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + psl['HWP']['RepairsToLandfillWood']*(1-psl['HWP']['ToLandfillWoodDegradableFrac'])*C_retired
    
    #--------------------------------------------------------------------------
    # Other --> dump and landfill
    #--------------------------------------------------------------------------
    
    # Turnover
    ip=meta['iPP']['Other']
    C_retired=psl['HWP']['Other_tr']*vo['C_Pro_Pools'][iT-1,:,ip]
    
    # Remove carbon
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] - C_retired
    
    # Transfer carbon to dump wood
    ip=meta['iPP']['DumpWood']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + psl['HWP']['OtherToDumpWood']*C_retired
    
    # Transfer carbon to landfill (degradble)
    ip=meta['iPP']['LandfillWoodDegradable']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + psl['HWP']['OtherToLandfillWood']*psl['HWP']['ToLandfillWoodDegradableFrac']*C_retired
    
    # Transfer carbon to landfill (non-degradable)
    ip=meta['iPP']['LandfillWoodNonDegradable']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + psl['HWP']['OtherToLandfillWood']*(1-psl['HWP']['ToLandfillWoodDegradableFrac'])*C_retired
    
    #--------------------------------------------------------------------------
    # Cants --> dump and landfill
    #--------------------------------------------------------------------------
    
    # Turnover
    ip=meta['iPP']['Cants']
    C_retired=psl['HWP']['Cants_tr']*vo['C_Pro_Pools'][iT-1,:,ip]
    
    # Remove carbon
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] - C_retired
    
    # Transfer carbon to dump wood
    ip=meta['iPP']['DumpWood']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + psl['HWP']['CantsToDumpWood']*C_retired
    
    # Transfer carbon to landfill (degradble)
    ip=meta['iPP']['LandfillWoodDegradable']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + psl['HWP']['CantsToLandfillWood']*psl['HWP']['ToLandfillWoodDegradableFrac']*C_retired
    
    # Transfer carbon to landfill (non-degradable)
    ip=meta['iPP']['LandfillWoodNonDegradable']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + psl['HWP']['CantsToLandfillWood']*(1-psl['HWP']['ToLandfillWoodDegradableFrac'])*C_retired
    
    #--------------------------------------------------------------------------
    # Paper --> dump and landfill
    #--------------------------------------------------------------------------
    
    # Turnover (with adjustment for recycling)
    ip=meta['iPP']['Paper']
    C_retired=(1-psl['HWP']['PaperRecycleRate'])*psl['HWP']['Paper_tr']*vo['C_Pro_Pools'][iT,:,ip]
    
    # Remove carbon
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] - C_retired
    
    # Transfer to dump
    ip=meta['iPP']['DumpPaper']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + psl['HWP']['PaperToDumpPaper']*C_retired
    
    # Transfer to landfill (degradable)
    ip=meta['iPP']['LandfillWoodDegradable']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + psl['HWP']['PaperToLandfillPaper']*psl['HWP']['ToLandfillPaperDegradableFrac']*C_retired
    
    # Transfer to landfill (non-degradable)
    ip=meta['iPP']['LandfillWoodNonDegradable']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + psl['HWP']['PaperToLandfillPaper']*(1-psl['HWP']['ToLandfillPaperDegradableFrac'])*C_retired
        
    #--------------------------------------------------------------------------
    # Emissions from fuel combustion
    #--------------------------------------------------------------------------
    
    # Turnover
    ip=meta['iPP']['Fuel']
    C_retired=psl['HWP']['Fuel_tr']*vo['C_Pro_Pools'][iT,:,ip]
    
    # Remove carbon
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] - C_retired
    
    # Emissions of CO2 from fuel use
    ip=meta['iPP']['E_CO2']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + psl['HWP']['FuelCombustionFracEmitCO2']*C_retired
    
    # Emissions of CH4 from fuel use
    ip=meta['iPP']['E_CH4']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + (1-psl['HWP']['FuelCombustionFracEmitCO2'])*C_retired
    
    #--------------------------------------------------------------------------
    # Emissions from firewood combustion
    #--------------------------------------------------------------------------
    
    # Turnover
    ip=meta['iPP']['Firewood']
    C_retired=psl['HWP']['Firewood_tr']*vo['C_Pro_Pools'][iT,:,ip]
    
    # Remove carbon
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] - C_retired
    
    # Emissions of CO2
    ip=meta['iPP']['E_CO2']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + C_retired
    
    # Emissions of CH4
    #ip=meta['iPP']['E_CH4']
    #vo['C_Pro_Pools[iT,:,ip]=vo['C_Pro_Pools[iT,:,ip] + (1-psl['FuelCombustionFracEmitCO2'])*C_retired
       
    #--------------------------------------------------------------------------
    # Emissions from pulp effluent
    #--------------------------------------------------------------------------
         
    # Emissions from pulp effluent (CO2 from aerobic decomposition)
    ip=meta['iPP']['EffluentPulp']
    c_emitted=psl['HWP']['EffluentPulp_tr']*vo['C_Pro_Pools'][iT,:,ip]
    
    # Remove emitted carbon
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip]-c_emitted
    
    # Add emitted carbon to CO2 emission "pool"
    ip=meta['iPP']['E_CO2']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip]+c_emitted    
    
    #--------------------------------------------------------------------------
    # Emissions from dump wood
    #--------------------------------------------------------------------------
    
    # Turnover
    ip=meta['iPP']['DumpWood']
    c_emitted=psl['HWP']['DumpWood_tr']*vo['C_Pro_Pools'][iT,:,ip]
    
    # Removal
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip]-c_emitted
    
    # Add to emissions (CO2 emission from aerobic decomposition)
    ip=meta['iPP']['E_CO2']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip]+c_emitted
    
    #--------------------------------------------------------------------------
    # Emissions from dump paper
    #--------------------------------------------------------------------------
    
    # Turnover
    ip=meta['iPP']['DumpPaper']
    c_emitted=psl['HWP']['DumpPaper_tr']*vo['C_Pro_Pools'][iT,:,ip]
    
    # Removal
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip]-c_emitted
    
    # Add to emissions (CO2 emission from aerobic decomposition)
    ip=meta['iPP']['E_CO2']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip]+c_emitted
    
    #--------------------------------------------------------------------------
    # Emissions from landfill degradable wood
    #--------------------------------------------------------------------------
               
    # Turnover
    ip=meta['iPP']['LandfillWoodDegradable']
    c_emitted=psl['HWP']['LandfillWoodDegradable_tr']*vo['C_Pro_Pools'][iT,:,ip]
    
    # Removal
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip]-c_emitted
    
    # Add to emissions (50% CO2 emissions during anaerobic decomposition)
    ip=meta['iPP']['E_CO2']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip]+psl['HWP']['LandfillDegradableFracEmitCO2']*c_emitted
    
    # Add to emissions (50% "potential" CH4 emissions during anaerobic decomposition)
    #E_ch4_pot=(1-psl['HWP']['LandfillDegradableFracEmitCO2)*c_emitted
    
    # Adjustment for proportion of degradable landfills with gas collection systems, 
    # efficiency of system, and methane oxided to CO2 from the landfill cover
    ch4_emitted=c_emitted*((1-psl['HWP']['LandfillMethaneEmit_GasColSysProp'])-psl['HWP']['LandfillMethaneOxidizedToCO2']*(1-psl['HWP']['LandfillMethaneEmit_GasColSysProp'])) + \
        c_emitted*psl['HWP']['LandfillMethaneEmit_GasColSysProp']*((1-psl['HWP']['LandfillMethaneEmit_GasColSysEffic'])-psl['HWP']['LandfillMethaneOxidizedToCO2']*(1-psl['HWP']['LandfillMethaneEmit_GasColSysProp']))
        
    ip=meta['iPP']['E_CH4']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip]+ch4_emitted
    
    #--------------------------------------------------------------------------
    # Emissions from landfill degradable paper
    #--------------------------------------------------------------------------
    
    # Turnover
    ip=meta['iPP']['LandfillPaperDegradable']
    c_emitted=psl['HWP']['LandfillWoodDegradable_tr']*vo['C_Pro_Pools'][iT,:,ip]
    
    # Removal
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip]-c_emitted
    
    # Add to emissions (50% CO2 emissions during anaerobic decomposition)
    ip=meta['iPP']['E_CO2']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip]+psl['HWP']['LandfillDegradableFracEmitCO2']*c_emitted
    
    # Add to emissions (50% "potential" CH4 emissions during anaerobic decomposition)
    #E_ch4_pot=(1-psl['HWP']['LandfillDegradableFracEmitCO2)*c_emitted
    
    # Adjustment for proportion of degradable landfills with gas collection systems, 
    # efficiency of system, and methane oxided to CO2 from the landfill cover
    ch4_emitted=c_emitted*((1-psl['HWP']['LandfillMethaneEmit_GasColSysProp'])-psl['HWP']['LandfillMethaneOxidizedToCO2']*(1-psl['HWP']['LandfillMethaneEmit_GasColSysProp'])) + \
        c_emitted*psl['HWP']['LandfillMethaneEmit_GasColSysProp']*((1-psl['HWP']['LandfillMethaneEmit_GasColSysEffic'])-psl['HWP']['LandfillMethaneOxidizedToCO2']*(1-psl['HWP']['LandfillMethaneEmit_GasColSysProp']))
        
    ip=meta['iPP']['E_CH4']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip]+ch4_emitted
    
    return vo


#%% Disturbance and management events

def DisturbanceAndManagementEvents_old(iT,vi,vo,psl,meta,iEP):
    
    # Predict harvesting (on the fly)
    if meta['Simulate harvesting on the fly']=='On':            
        if meta['Dist on Fly']['Harvesting'][meta['iScn']]==1:
            Oc_har_sim,id_har_sim,mort_har_sim=gensm.PredictHarvesting_OnTheFly(meta,vi['Inv']['THLB'][0,:],vo['V_StemMerch'][iT,:],vi['tv'][iT])
        else:
            Oc_har_sim=np.zeros(meta['N Stand'])
    
    # Predict stand breakup (on the fly)
    if meta['Simulate breakup on the fly']=='On':
        Oc_bu_sim,id_bu_sim,mort_bu_sim=gensm.PredictStandBreakup_OnTheFly(meta,vo['A'][iT,:])
    
    # Update disturbance flag
    #meta['Flag_Disturbed']=np.zeros(meta['N Stand'])
    
    for iS in range(meta['N Stand']):

        # Add harvesting (on the fly)
        if meta['Simulate harvesting on the fly']=='On':
            if Oc_har_sim[iS]==1:
                Year=vi['tv'][iT]
                # Add harvest  
                c=vi['EH'][iS]['cnt_fly']
                vi['EH'][iS]['Year'][c]=Year
                vi['EH'][iS]['ID_Type'][c]=id_har_sim[iS]
                vi['EH'][iS]['MortalityFactor'][c]=mort_har_sim[iS]
                vi['EH'][iS]['GrowthFactor'][c]=0
                vi['EH'][iS]['ID_GrowthCurve'][c]=1
                vi['EH'][iS]['ID_GrowthCurveM'][c]=1
                for k in range(len(meta['Name CustDistVar'])):
                    vnam=meta['Name CustDistVar'][k]
                    vi['EH'][iS][vnam][c]=0
                vi['EH'][iS]['cnt_fly']=vi['EH'][iS]['cnt_fly']+1
                
                # Add slashpile burn                
                c=vi['EH'][iS]['cnt_fly']
                vi['EH'][iS]['Year'][c]=Year
                vi['EH'][iS]['ID_Type'][c]=meta['LUT Dist']['Slashpile Burn']
                vi['EH'][iS]['MortalityFactor'][c]=100
                vi['EH'][iS]['GrowthFactor'][c]=0
                vi['EH'][iS]['ID_GrowthCurve'][c]=1
                vi['EH'][iS]['ID_GrowthCurveM'][c]=1
                for k in range(len(meta['Name CustDistVar'])):
                    vnam=meta['Name CustDistVar'][k]
                    vi['EH'][iS][vnam][c]=0
                vi['EH'][iS]['cnt_fly']=vi['EH'][iS]['cnt_fly']+1
                
                # Add planting
                c=vi['EH'][iS]['cnt_fly']
                vi['EH'][iS]['Year'][c]=Year
                vi['EH'][iS]['ID_Type'][c]=meta['LUT Dist']['Planting']
                vi['EH'][iS]['MortalityFactor'][c]=100
                vi['EH'][iS]['GrowthFactor'][c]=0
                vi['EH'][iS]['ID_GrowthCurve'][c]=1
                vi['EH'][iS]['ID_GrowthCurveM'][c]=1
                for k in range(len(meta['Name CustDistVar'])):
                    vnam=meta['Name CustDistVar'][k]
                    vi['EH'][iS][vnam][c]=0
                vi['EH'][iS]['cnt_fly']=vi['EH'][iS]['cnt_fly']+1
                
        # Add stand breakup (on the fly)
        if meta['Simulate breakup on the fly']=='On':
            if Oc_bu_sim[iS]==1:
                Year=vi['tv'][iT]
                # Add breakup
                c=vi['EH'][iS]['cnt_fly']
                vi['EH'][iS]['Year'][c]=Year
                vi['EH'][iS]['ID_Type'][c]=id_bu_sim[iS]
                vi['EH'][iS]['MortalityFactor'][c]=mort_bu_sim[iS]
                vi['EH'][iS]['GrowthFactor'][c]=0
                vi['EH'][iS]['ID_GrowthCurve'][c]=1
                vi['EH'][iS]['ID_GrowthCurveM'][c]=1
                for k in range(len(meta['Name CustDistVar'])):
                    vnam=meta['Name CustDistVar'][k]
                    vi['EH'][iS][vnam][c]=0
                vi['EH'][iS]['cnt_fly']=vi['EH'][iS]['cnt_fly']+1
                    
        # Index to the disturbance
        indDist=np.where(np.floor(vi['EH'][iS]['Year'])==vi['tv'][iT])[0]
        
        # Only continue if there are disturbances
        if indDist.size==0:
            continue
        
        for iDist in range(indDist.size):
            
            # Disturbance type
            ID_Type=vi['EH'][iS]['ID_Type'][indDist[iDist]]
            
            if ID_Type==0:
                continue

            # Index to default parameters for this disturbance type
            #iDefPar=np.where(psl['bDist_ID==ID_Type)[0]

            # Update disturbance flag
            #if (ID_Type==meta['LUT Dist']['Wildfire']) | (ID_Type==meta['LUT Dist']['IBM']) | (ID_Type==meta['LUT Dist']['Harvest']) | (ID_Type==meta['LUT Dist']['Harvest and Slashpile Burn']):
            #    meta['Flag_Disturbed'][iS]==1
            
            # Define affected biomass fractions
            if vi['EH'][iS]['Biomass_Affected_Pct'][indDist[iDist]]<=0:
                
                #--------------------------------------------------------------
                # No custom harvesting disturbance parameters supplied, use default
                # assumptions
                #--------------------------------------------------------------
                
                # Use mortality from disturbance history
                MortalityFraction=vi['EH'][iS]['MortalityFactor'][indDist[iDist]].astype(float)/100
                
                # Force Aerial Fertilization and Planting to be zero severity
                if (ID_Type==meta['LUT Dist']['Fertilization Aerial']) | (ID_Type==meta['LUT Dist']['Planting']):
                    MortalityFraction=0
                
                # Apply severity to biomass and snags
                Biomass_Affected_Frac=MortalityFraction
                Snags_Affected_Frac=MortalityFraction
                
#                # Use default disturbance fluxes
#                Biomass_Merch_Removed_Frac=psl['bDist_Biomass_Merch_Removed[iDefPar]
#                Biomass_Merch_Burned_Frac=psl['bDist_Biomass_Merch_Burned[iDefPar]
#                Biomass_Merch_LeftOnSite_Frac=psl['bDist_Biomass_Merch_LeftOnSite[iDefPar]
#                Biomass_Merch_ToSnag_Frac=psl['bDist_Biomass_Merch_ToSnag[iDefPar]
#                Biomass_NonMerch_Removed_Frac=psl['bDist_Biomass_NonMerch_Removed[iDefPar]
#                Biomass_NonMerch_Burned_Frac=psl['bDist_Biomass_NonMerch_Burned[iDefPar]
#                Biomass_NonMerch_LeftOnSite_Frac=psl['bDist_Biomass_NonMerch_LeftOnSite[iDefPar]
#                Biomass_NonMerch_ToSnag_Frac=psl['bDist_Biomass_NonMerch_ToSnag[iDefPar]
#                Snags_Removed_Frac=psl['bDist_Snags_Removed[iDefPar]
#                Snags_Burned_Frac=psl['bDist_Snags_Burned[iDefPar]
#                Snags_LeftOnSite_Frac=psl['bDist_Snags_LeftOnSite[iDefPar]                
#                FelledStemMerch_Burned_Frac=psl['bDist_FelledStemMerch_Burned[iDefPar]
#                FelledStemNonMerch_Burned_Frac=psl['bDist_FelledStemNonMerch_Burned[iDefPar]
#                FelledBranch_Burned_Frac=psl['bDist_FelledBranch_Burned[iDefPar]
#                FelledBark_Burned_Frac=psl['bDist_FelledBark_Burned[iDefPar]
#                FelledSnagStem_Burned_Frac=psl['bDist_FelledSnagStem_Burned[iDefPar]
#                FelledSnagBranch_Burned_Frac=psl['bDist_FelledSnagBranch_Burned[iDefPar]
                
            else:
                
                #--------------------------------------------------------------
                # Use custom harvesting disturbance parameter table
                #--------------------------------------------------------------
                
                # Use mortality from disturbance history
                MortalityFraction=vi['EH'][iS]['Biomass_Affected_Pct'][indDist[iDist]].astype(float)/100  
                
                # These come from a custom user input spreadsheet, but they have been added to DH structure in pre-processing
                Biomass_Affected_Frac=vi['EH'][iS]['Biomass_Affected_Pct'][indDist[iDist]].astype(float)/100                
                Snags_Affected_Frac=vi['EH'][iS]['Snags_Affected_Pct'][indDist[iDist]].astype(float)/100
                
                # Proportions removed, burned, and left on site                        
                Biomass_Merch_Removed_Frac=vi['EH'][iS]['Biomass_Merch_Removed_Pct'][indDist[iDist]].astype(float)/100
                Biomass_Merch_Burned_Frac=vi['EH'][iS]['Biomass_Merch_Burned_Pct'][indDist[iDist]].astype(float)/100          
                Biomass_Merch_LeftOnSite_Frac=vi['EH'][iS]['Biomass_Merch_LeftOnSite_Pct'][indDist[iDist]].astype(float)/100
                Biomass_Merch_ToSnag_Frac=0         
                
                Biomass_NonMerch_Removed_Frac=vi['EH'][iS]['Biomass_NonMerch_Removed_Pct'][indDist[iDist]].astype(float)/100
                Biomass_NonMerch_Burned_Frac=vi['EH'][iS]['Biomass_NonMerch_Burned_Pct'][indDist[iDist]].astype(float)/100          
                Biomass_NonMerch_LeftOnSite_Frac=vi['EH'][iS]['Biomass_NonMerch_LeftOnSite_Pct'][indDist[iDist]].astype(float)/100
                Biomass_NonMerch_ToSnag_Frac=0
                
                Snags_Removed_Frac=vi['EH'][iS]['Snags_Removed_Pct'][indDist[iDist]].astype(float)/100
                Snags_Burned_Frac=vi['EH'][iS]['Snags_Burned_Pct'][indDist[iDist]].astype(float)/100
                Snags_LeftOnSite_Frac=vi['EH'][iS]['Snags_LeftOnSite_Pct'][indDist[iDist]].astype(float)/100
                  
                FelledStemMerch_Burned_Frac=0
                FelledStemNonMerch_Burned_Frac=0
                FelledBranch_Burned_Frac=0
                FelledBark_Burned_Frac=0
                FelledSnagStem_Burned_Frac=0
                FelledSnagBranch_Burned_Frac=0
                
            #------------------------------------------------------------------
            # Total amount of each pool that is affected
            #------------------------------------------------------------------
            
            Affected_StemMerch=Biomass_Affected_Frac*vo['C_Eco_Pools'][iT,iS,iEP['StemMerch']]
            Affected_StemNonMerch=Biomass_Affected_Frac*vo['C_Eco_Pools'][iT,iS,iEP['StemNonMerch']]
            Affected_Foliage=Biomass_Affected_Frac*vo['C_Eco_Pools'][iT,iS,iEP['Foliage']]
            Affected_Branch=Biomass_Affected_Frac*vo['C_Eco_Pools'][iT,iS,iEP['Branch']]
            Affected_Bark=Biomass_Affected_Frac*vo['C_Eco_Pools'][iT,iS,iEP['Bark']]            
            Affected_RootCoarse=Biomass_Affected_Frac*vo['C_Eco_Pools'][iT,iS,iEP['RootCoarse']]
            Affected_RootFine=Biomass_Affected_Frac*vo['C_Eco_Pools'][iT,iS,iEP['RootFine']]
             
            # Volume
            Affected_VolumeStemMerch=Biomass_Affected_Frac*vo['V_StemMerch'][iT,iS]
                        
            # Partition bark into merch and non-merch components
            Affected_BarkMerch=psl['bASL_MerchBarkFrac']*Affected_Bark
            Affected_BarkNonMerch=(1-psl['bASL_MerchBarkFrac'])*Affected_Bark
            
            # Sum up total affected non-merchantable biomass
            Affected_TotNonMerch=Affected_StemNonMerch+Affected_Branch+Affected_BarkNonMerch
            
            # Define fraction of snags that are affected            
            Affected_SnagStem=Snags_Affected_Frac*vo['C_Eco_Pools'][iT,iS,iEP['SnagStem']]
            Affected_SnagBranch=Snags_Affected_Frac*vo['C_Eco_Pools'][iT,iS,iEP['SnagBranch']]
            
            # All
            Affected_All=Affected_StemMerch+Affected_StemNonMerch+Affected_Foliage+Affected_Branch+Affected_Bark+Affected_RootCoarse+Affected_RootFine
            
            #------------------------------------------------------------------
            # Calculate mortality
            #------------------------------------------------------------------
            
            vo['C_M_Dist'][iT,iS]=vo['C_M_Dist'][iT,iS]+Affected_All
            
            DistType_i=vi['EH'][iS]['ID_Type'][indDist[iDist]]
                 
            if DistType_i==meta['LUT Dist']['Wildfire']:
                vo['C_M_Inv_Fir'][iT,iS]=Affected_All
            elif (DistType_i==meta['LUT Dist']['Beetles']) | (DistType_i==meta['LUT Dist']['IBM']) | \
                (DistType_i==meta['LUT Dist']['IBB']) | (DistType_i==meta['LUT Dist']['IBS']) | \
                (DistType_i==meta['LUT Dist']['IBD']):
                vo['C_M_Inv_Ins'][iT,iS]=Affected_All
            elif DistType_i==meta['LUT Dist']['Defoliators']:
                vo['C_M_Inv_Ins'][iT,iS]=Affected_All 
            elif (DistType_i==meta['LUT Dist']['Harvest']) | (DistType_i==meta['LUT Dist']['Harvest Custom']) | \
                 (DistType_i==meta['LUT Dist']['Salvage Logging']) | (DistType_i==meta['LUT Dist']['Harvest and Slashpile Burn']):
                vo['C_M_Inv_Har'][iT,iS]=Affected_All            
            
            #------------------------------------------------------------------
            # Remove carbon from each affected pool
            #------------------------------------------------------------------
            
            # Remove carbon from affected biomass pools            
            vo['C_Eco_Pools'][iT,iS,iEP['StemMerch']]=vo['C_Eco_Pools'][iT,iS,iEP['StemMerch']]-Affected_StemMerch
            vo['C_Eco_Pools'][iT,iS,iEP['StemNonMerch']]=vo['C_Eco_Pools'][iT,iS,iEP['StemNonMerch']]-Affected_StemNonMerch
            vo['C_Eco_Pools'][iT,iS,iEP['Foliage']]=vo['C_Eco_Pools'][iT,iS,iEP['Foliage']]-Affected_Foliage
            vo['C_Eco_Pools'][iT,iS,iEP['Branch']]=vo['C_Eco_Pools'][iT,iS,iEP['Branch']]-Affected_Branch
            vo['C_Eco_Pools'][iT,iS,iEP['Bark']]=vo['C_Eco_Pools'][iT,iS,iEP['Bark']]-Affected_Bark
            vo['C_Eco_Pools'][iT,iS,iEP['RootCoarse']]=vo['C_Eco_Pools'][iT,iS,iEP['RootCoarse']]-Affected_RootCoarse
            vo['C_Eco_Pools'][iT,iS,iEP['RootFine']]=vo['C_Eco_Pools'][iT,iS,iEP['RootFine']]-Affected_RootFine
               
            # Remove carbon from snag pools
            vo['C_Eco_Pools'][iT,iS,iEP['SnagStem']]=vo['C_Eco_Pools'][iT,iS,iEP['SnagStem']]-Affected_SnagStem
            vo['C_Eco_Pools'][iT,iS,iEP['SnagBranch']]=vo['C_Eco_Pools'][iT,iS,iEP['SnagBranch']]-Affected_SnagBranch
            
            # Remove stemwood merch volume
            vo['V_StemMerch'][iT,iS]=vo['V_StemMerch'][iT,iS]-Affected_VolumeStemMerch
            
            #------------------------------------------------------------------
            # Carbon removed (sent to mill)
            #------------------------------------------------------------------
            
            # Merch biomass to mill - of the total amount of biomass affected, 
            # Add bark
            Biomass_Merch_Removed=Biomass_Merch_Removed_Frac*(Affected_StemMerch+Affected_BarkMerch)
            vo['C_Eco_Pools'][iT,iS,iEP['RemovedMerch']]=vo['C_Eco_Pools'][iT,iS,iEP['RemovedMerch']]+Biomass_Merch_Removed
            
            # NonMerch biomass to mill - of the total amount of biomass affected, 
            # what fraction of non-merch biomass was sent to the mill?
            # - NonMerch = NonMerchStem + Foliage + Branch + Bark
            Biomass_NonMerch_Removed=Biomass_NonMerch_Removed_Frac*Affected_TotNonMerch
            vo['C_Eco_Pools'][iT,iS,iEP['RemovedNonMerch']]=vo['C_Eco_Pools'][iT,iS,iEP['RemovedNonMerch']]+Biomass_NonMerch_Removed
                                        
            # Snag stemwood to mill           
            SnagStem_Removed=Snags_Removed_Frac*Affected_SnagStem
            vo['C_Eco_Pools'][iT,iS,iEP['RemovedSnagStem']]=vo['C_Eco_Pools'][iT,iS,iEP['RemovedSnagStem']]+SnagStem_Removed
            
            # Snag branches to mill
            SnagBranch_Removed=Snags_Removed_Frac*Affected_SnagBranch
            vo['C_Eco_Pools'][iT,iS,iEP['RemovedSnagStem']]=vo['C_Eco_Pools'][iT,iS,iEP['RemovedSnagStem']]+SnagBranch_Removed
            
            #------------------------------------------------------------------
            # Carbon that is left on site (after felling or wind storms)
            #------------------------------------------------------------------
            
            # Stem, branch and bark carbon transfered to felled pools (piles)
            vo['C_Eco_Pools'][iT,iS,iEP['FelledStemMerch']]=vo['C_Eco_Pools'][iT,iS,iEP['FelledStemMerch']]+Biomass_Merch_LeftOnSite_Frac*Affected_StemMerch
            vo['C_Eco_Pools'][iT,iS,iEP['FelledStemNonMerch']]=vo['C_Eco_Pools'][iT,iS,iEP['FelledStemNonMerch']]+Biomass_NonMerch_LeftOnSite_Frac*Affected_StemNonMerch            
            vo['C_Eco_Pools'][iT,iS,iEP['FelledBranch']]=vo['C_Eco_Pools'][iT,iS,iEP['FelledBranch']]+Biomass_Merch_LeftOnSite_Frac*Affected_Branch
            vo['C_Eco_Pools'][iT,iS,iEP['FelledBark']]=vo['C_Eco_Pools'][iT,iS,iEP['FelledBark']]+Biomass_Merch_LeftOnSite_Frac*Affected_Bark
                        
            # Foliage transferred directly to very fast litter
            vo['C_Eco_Pools'][iT,iS,iEP['LitterVF']]=vo['C_Eco_Pools'][iT,iS,iEP['LitterVF']]+Affected_Foliage
            
            # Roots transferred directly to DOM       
            vo['C_Eco_Pools'][iT,iS,iEP['LitterF']]=vo['C_Eco_Pools'][iT,iS,iEP['LitterF']]+0.5*Affected_RootCoarse
            vo['C_Eco_Pools'][iT,iS,iEP['SoilF']]=vo['C_Eco_Pools'][iT,iS,iEP['SoilF']]+0.5*Affected_RootCoarse            
            vo['C_Eco_Pools'][iT,iS,iEP['LitterVF']]=vo['C_Eco_Pools'][iT,iS,iEP['LitterVF']]+0.5*Affected_RootFine
            vo['C_Eco_Pools'][iT,iS,iEP['SoilVF']]=vo['C_Eco_Pools'][iT,iS,iEP['SoilVF']]+0.5*Affected_RootFine
            
            # Snags transferred to felled pools (piles)
            vo['C_Eco_Pools'][iT,iS,iEP['FelledSnagStem']]=vo['C_Eco_Pools'][iT,iS,iEP['FelledSnagStem']]+Snags_LeftOnSite_Frac*Affected_SnagStem
            vo['C_Eco_Pools'][iT,iS,iEP['FelledSnagBranch']]=vo['C_Eco_Pools'][iT,iS,iEP['FelledSnagBranch']]+Snags_LeftOnSite_Frac*Affected_SnagBranch
                        
            #------------------------------------------------------------------
            # Biomass and DOM burned
            #------------------------------------------------------------------
            
            # Profile of emitted compounds (N2O emissions calculated upon export)
            psl['bCombFrac_CO2']=psl['bCombFrac_CO2']
            psl['bCombFrac_CH4']=psl['bCombFrac_CH4']
            psl['bCombFrac_CO']=psl['bCombFrac_CO']
            
            # Merch biomass that is burned              
            Biomass_Merch_Burned_co2=psl['bCombFrac_CO2']*Biomass_Merch_Burned_Frac*(Affected_StemMerch+Affected_BarkMerch)
            vo['C_Eco_Pools'][iT,iS,iEP['ECO2asC']]=vo['C_Eco_Pools'][iT,iS,iEP['ECO2asC']]+Biomass_Merch_Burned_co2
            
            Biomass_Merch_Burned_ch4=psl['bCombFrac_CH4']*Biomass_Merch_Burned_Frac*(Affected_StemMerch+Affected_BarkMerch)
            vo['C_Eco_Pools'][iT,iS,iEP['ECH4asC']]=vo['C_Eco_Pools'][iT,iS,iEP['ECH4asC']]+Biomass_Merch_Burned_ch4
            
            Biomass_Merch_Burned_co=psl['bCombFrac_CO']*Biomass_Merch_Burned_Frac*(Affected_StemMerch+Affected_BarkMerch)
            vo['C_Eco_Pools'][iT,iS,iEP['ECOasC']]=vo['C_Eco_Pools'][iT,iS,iEP['ECOasC']]+Biomass_Merch_Burned_co

            # NonMerch biomass that is burned            
            Biomass_NonMerch_Burned_co2=psl['bCombFrac_CO2']*Biomass_NonMerch_Burned_Frac*Affected_TotNonMerch
            vo['C_Eco_Pools'][iT,iS,iEP['ECO2asC']]=vo['C_Eco_Pools'][iT,iS,iEP['ECO2asC']]+Biomass_NonMerch_Burned_co2
            
            Biomass_NonMerch_Burned_ch4=psl['bCombFrac_CH4']*Biomass_NonMerch_Burned_Frac*Affected_TotNonMerch
            vo['C_Eco_Pools'][iT,iS,iEP['ECH4asC']]=vo['C_Eco_Pools'][iT,iS,iEP['ECH4asC']]+Biomass_NonMerch_Burned_ch4
            
            Biomass_NonMerch_Burned_co=psl['bCombFrac_CO']*Biomass_NonMerch_Burned_Frac*Affected_TotNonMerch
            vo['C_Eco_Pools'][iT,iS,iEP['ECOasC']]=vo['C_Eco_Pools'][iT,iS,iEP['ECOasC']]+Biomass_NonMerch_Burned_co
              
            # Snag stems that burn
            SnagStem_Burned_co2=psl['bCombFrac_CO2']*Snags_Burned_Frac*Affected_SnagStem
            vo['C_Eco_Pools'][iT,iS,iEP['ECO2asC']]=vo['C_Eco_Pools'][iT,iS,iEP['ECO2asC']]+SnagStem_Burned_co2
            
            SnagStem_Burned_ch4=psl['bCombFrac_CH4']*Snags_Burned_Frac*Affected_SnagStem
            vo['C_Eco_Pools'][iT,iS,iEP['ECH4asC']]=vo['C_Eco_Pools'][iT,iS,iEP['ECH4asC']]+SnagStem_Burned_ch4
            
            SnagStem_Burned_co=psl['bCombFrac_CO']*Snags_Burned_Frac*Affected_SnagStem
            vo['C_Eco_Pools'][iT,iS,iEP['ECOasC']]=vo['C_Eco_Pools'][iT,iS,iEP['ECOasC']]+SnagStem_Burned_co
            
            # Snag branches that burn
            SnagBranch_Burned_co2=psl['bCombFrac_CO2']*Snags_Burned_Frac*Affected_SnagBranch
            vo['C_Eco_Pools'][iT,iS,iEP['ECO2asC']]=vo['C_Eco_Pools'][iT,iS,iEP['ECO2asC']]+SnagBranch_Burned_co2
            
            SnagBranch_Burned_ch4=psl['bCombFrac_CH4']*Snags_Burned_Frac*Affected_SnagBranch
            vo['C_Eco_Pools'][iT,iS,iEP['ECH4asC']]=vo['C_Eco_Pools'][iT,iS,iEP['ECH4asC']]+SnagBranch_Burned_ch4
            
            SnagBranch_Burned_co=psl['bCombFrac_CO']*Snags_Burned_Frac*Affected_SnagBranch
            vo['C_Eco_Pools'][iT,iS,iEP['ECOasC']]=vo['C_Eco_Pools'][iT,iS,iEP['ECOasC']]+SnagBranch_Burned_co
       
            #------------------------------------------------------------------
            # Carbon that is moved from biomass to snags
            #------------------------------------------------------------------
            
            # Merch biomass that is killed by wildfire
            vo['C_Eco_Pools'][iT,iS,iEP['SnagStem']]=vo['C_Eco_Pools'][iT,iS,iEP['SnagStem']]+Biomass_Merch_ToSnag_Frac*(Affected_StemMerch+Affected_BarkMerch)
            
            # Non-merch biomass that is killed by wildfire
            vo['C_Eco_Pools'][iT,iS,iEP['SnagBranch']]=vo['C_Eco_Pools'][iT,iS,iEP['SnagBranch']]+Biomass_NonMerch_ToSnag_Frac*(Affected_StemNonMerch+Affected_BarkNonMerch)
            
            #------------------------------------------------------------------
            # Felled carbon that is burned
            #------------------------------------------------------------------            
            
            # Felled stem merch that burn
            Affected=FelledStemMerch_Burned_Frac*vo['C_Eco_Pools'][iT,iS,iEP['FelledStemMerch']]
            vo['C_Eco_Pools'][iT,iS,iEP['FelledStemMerch']]=vo['C_Eco_Pools'][iT,iS,iEP['FelledStemMerch']]-Affected
            Burned_As_CO2=psl['bCombFrac_CO2']*Affected
            vo['C_Eco_Pools'][iT,iS,iEP['ECO2asC']]=vo['C_Eco_Pools'][iT,iS,iEP['ECO2asC']]+Burned_As_CO2            
            Burned_As_CH4=psl['bCombFrac_CH4']*Affected
            vo['C_Eco_Pools'][iT,iS,iEP['ECH4asC']]=vo['C_Eco_Pools'][iT,iS,iEP['ECH4asC']]+Burned_As_CH4            
            Burned_As_CO=psl['bCombFrac_CO']*Affected
            vo['C_Eco_Pools'][iT,iS,iEP['ECOasC']]=vo['C_Eco_Pools'][iT,iS,iEP['ECOasC']]+Burned_As_CO
          
            # Felled non merch stem that burn
            Affected=FelledStemNonMerch_Burned_Frac*vo['C_Eco_Pools'][iT,iS,iEP['FelledStemNonMerch']]
            vo['C_Eco_Pools'][iT,iS,iEP['FelledStemNonMerch']]=vo['C_Eco_Pools'][iT,iS,iEP['FelledStemNonMerch']]-Affected
            Burned_As_CO2=psl['bCombFrac_CO2']*Affected
            vo['C_Eco_Pools'][iT,iS,iEP['ECO2asC']]=vo['C_Eco_Pools'][iT,iS,iEP['ECO2asC']]+Burned_As_CO2            
            Burned_As_CH4=psl['bCombFrac_CH4']*Affected
            vo['C_Eco_Pools'][iT,iS,iEP['ECH4asC']]=vo['C_Eco_Pools'][iT,iS,iEP['ECH4asC']]+Burned_As_CH4            
            Burned_As_CO=psl['bCombFrac_CO']*Affected
            vo['C_Eco_Pools'][iT,iS,iEP['ECOasC']]=vo['C_Eco_Pools'][iT,iS,iEP['ECOasC']]+Burned_As_CO
          
            # Felled branch that burn
            Affected=FelledBranch_Burned_Frac*vo['C_Eco_Pools'][iT,iS,iEP['FelledBranch']]
            vo['C_Eco_Pools'][iT,iS,iEP['FelledBranch']]=vo['C_Eco_Pools'][iT,iS,iEP['FelledBranch']]-Affected
            Burned_As_CO2=psl['bCombFrac_CO2']*Affected
            vo['C_Eco_Pools'][iT,iS,iEP['ECO2asC']]=vo['C_Eco_Pools'][iT,iS,iEP['ECO2asC']]+Burned_As_CO2            
            Burned_As_CH4=psl['bCombFrac_CH4']*Affected
            vo['C_Eco_Pools'][iT,iS,iEP['ECH4asC']]=vo['C_Eco_Pools'][iT,iS,iEP['ECH4asC']]+Burned_As_CH4            
            Burned_As_CO=psl['bCombFrac_CO']*Affected
            vo['C_Eco_Pools'][iT,iS,iEP['ECOasC']]=vo['C_Eco_Pools'][iT,iS,iEP['ECOasC']]+Burned_As_CO
          
            # Felled bark that burn
            Affected=FelledBark_Burned_Frac*vo['C_Eco_Pools'][iT,iS,iEP['FelledBark']]
            vo['C_Eco_Pools'][iT,iS,iEP['FelledBark']]=vo['C_Eco_Pools'][iT,iS,iEP['FelledBark']]-Affected
            Burned_As_CO2=psl['bCombFrac_CO2']*Affected
            vo['C_Eco_Pools'][iT,iS,iEP['ECO2asC']]=vo['C_Eco_Pools'][iT,iS,iEP['ECO2asC']]+Burned_As_CO2            
            Burned_As_CH4=psl['bCombFrac_CH4']*Affected
            vo['C_Eco_Pools'][iT,iS,iEP['ECH4asC']]=vo['C_Eco_Pools'][iT,iS,iEP['ECH4asC']]+Burned_As_CH4            
            Burned_As_CO=psl['bCombFrac_CO']*Affected
            vo['C_Eco_Pools'][iT,iS,iEP['ECOasC']]=vo['C_Eco_Pools'][iT,iS,iEP['ECOasC']]+Burned_As_CO
          
            # Felled snag stem that burn
            Affected=FelledSnagStem_Burned_Frac*vo['C_Eco_Pools'][iT,iS,iEP['FelledSnagStem']]
            vo['C_Eco_Pools'][iT,iS,iEP['FelledSnagStem']]=vo['C_Eco_Pools'][iT,iS,iEP['FelledSnagStem']]-Affected
            Burned_As_CO2=psl['bCombFrac_CO2']*Affected
            vo['C_Eco_Pools'][iT,iS,iEP['ECO2asC']]=vo['C_Eco_Pools'][iT,iS,iEP['ECO2asC']]+Burned_As_CO2            
            Burned_As_CH4=psl['bCombFrac_CH4']*Affected
            vo['C_Eco_Pools'][iT,iS,iEP['ECH4asC']]=vo['C_Eco_Pools'][iT,iS,iEP['ECH4asC']]+Burned_As_CH4            
            Burned_As_CO=psl['bCombFrac_CO']*Affected
            vo['C_Eco_Pools'][iT,iS,iEP['ECOasC']]=vo['C_Eco_Pools'][iT,iS,iEP['ECOasC']]+Burned_As_CO
          
            # Felled snag branch that burn
            Affected=FelledSnagBranch_Burned_Frac*vo['C_Eco_Pools'][iT,iS,iEP['FelledSnagBranch']]
            vo['C_Eco_Pools'][iT,iS,iEP['FelledSnagBranch']]=vo['C_Eco_Pools'][iT,iS,iEP['FelledSnagBranch']]-Affected
            Burned_As_CO2=psl['bCombFrac_CO2']*Affected
            vo['C_Eco_Pools'][iT,iS,iEP['ECO2asC']]=vo['C_Eco_Pools'][iT,iS,iEP['ECO2asC']]+Burned_As_CO2            
            Burned_As_CH4=psl['bCombFrac_CH4']*Affected
            vo['C_Eco_Pools'][iT,iS,iEP['ECH4asC']]=vo['C_Eco_Pools'][iT,iS,iEP['ECH4asC']]+Burned_As_CH4            
            Burned_As_CO=psl['bCombFrac_CO']*Affected
            vo['C_Eco_Pools'][iT,iS,iEP['ECOasC']]=vo['C_Eco_Pools'][iT,iS,iEP['ECOasC']]+Burned_As_CO
                      
            #------------------------------------------------------------------
            # Update stand age
            #------------------------------------------------------------------
            
            if meta['Biomass Module']=='BatchTIPSY':
                
                if (ID_Type==meta['LUT Dist']['Planting']) | (ID_Type==meta['LUT Dist']['Direct Seeding']) | (MortalityFraction==1):
                    
                    vo['A'][iT,iS]=0
                
                else:                    
                    
                    if (meta['Partial Mortality Affects Age']=='On') | (meta['Partial Mortality Affects Age']=='on'):
                        
                        # Assume oldest trees were most affected, reduce age in prportion
                        # with mortality rate
                        # Not always realistic, but see how net growth is affected.
                        vo['A'][iT,iS]=vo['A'][iT,iS]*(1-Biomass_Affected_Frac)
                        
                    else:                   
                        
                        # Assume no change
                        vo['A'][iT,iS]=vo['A'][iT,iS]
                    
            #------------------------------------------------------------------
            # Transition to new growth curve
            #------------------------------------------------------------------
            
            # Initialize a flag that indicates whether the growth curve changes
            flg_gc_change=0
            if vi['EH'][iS]['ID_GrowthCurveM'][indDist[iDist]]!=vi['ID_GCA'][iS]:
                flg_gc_change=1
            
            # Only applies to TIPSY
            if meta['Biomass Module']!='Sawtooth':            
                
                # Does not apply to fertilization
                if (ID_Type!=meta['LUT Dist']['Fertilization Aerial']) & (ID_Type!=meta['LUT Dist']['Fertilization Hand']):
                    
                    if vi['EH'][iS]['ID_GrowthCurveM'][indDist[iDist]]==1:                        
                        vi['GCA'][:,iS,:]=vi['GC1'][:,iS,:]     
                        vi['ID_GCA'][iS]=1
                    elif vi['EH'][iS]['ID_GrowthCurveM'][indDist[iDist]]==2:
                        vi['GCA'][:,iS,:]=vi['GC2'][:,iS,:]
                        vi['ID_GCA'][iS]=2
                    elif vi['EH'][iS]['ID_GrowthCurveM'][indDist[iDist]]==3:
                        vi['GCA'][:,iS,:]=vi['GC3'][:,iS,:]
                        vi['ID_GCA'][iS]=3
                    elif vi['EH'][iS]['ID_GrowthCurveM'][indDist[iDist]]==4:
                        vi['GCA'][:,iS,:]=vi['GC4'][:,iS,:]
                        vi['ID_GCA'][iS]=4
                    elif vi['EH'][iS]['ID_GrowthCurveM'][indDist[iDist]]==5:
                        vi['GCA'][:,iS,:]=vi['GC5'][:,iS,:]
                        vi['ID_GCA'][iS]=5    
    
            #------------------------------------------------------------------
            # Update net growth (in response to lethal events)
            #
            # Partial disturbances likely have a lasting impact on gross growth 
            # by removing growing stock. Without ad hoc 
            # adjustment, this will not be reflected in the growth curve from TIPSY.
            #
            # This will apply to instances where Mortality=100% only if the
            # user does not specify a change in growth curve. 
            # With this functionality, the user has the options of:
            #  1) Changing the growth curve
            #  2) Assuming a gradual recovery to the previous curve, as defined
            #     by the default recovery half life for that disturbance type.
            #------------------------------------------------------------------
            flg=0
            if (flg==1) & (meta['Biomass Module']!='Sawtooth'):
                
                # Half life
                hl=psl['bDist_GrowthRecovery_HL'][iDefPar]
                
                # Only proceed if:
                #   1) the disturbance has a lasting growth impact/recovery
                #   2) the user is not specifying a change in growth curve
                if (hl>0) & (flg_gc_change==0):
                    
                    # Extract net growth for active growth curve
                    NetGrowth=vi['GCA'][:,iS,:].copy().astype(float)*meta['GC']['Scale Factor']
                
                    # Growth pre-event
                    G_pre=NetGrowth.copy()
                
                    # Growth post-event
                    G_post=(1-Biomass_Affected_Frac)*NetGrowth.copy()
                
                    # Difference in growth
                    dG=G_pre-G_post
                
                    # Age vector
                    A=np.arange(1,meta['GC']['BatchTIPSY Maximum Age']+1,1)
                
                    # Time since disturbance                    
                    TSD=np.tile(A-vo['A'][iT,iS],(G_pre.shape[1],1)).T
                
                    # Relative effect
                    fTSD=1/(1+np.exp(-0.5*(TSD-hl)))
                    fTSD[0:int(vo['A'][iT,iS]),:]=0
                
                    # Growth recovery
                    G_Recovery=fTSD*dG
                
                    # New estimate of net growth
                    NetGrowthNew=G_post+G_Recovery
                
                    # Add back to dictionary
                    NetGrowthNew=NetGrowthNew*meta['GC']['Scale Factor']
                    NetGrowthNew=NetGrowthNew.astype(np.int16)                
                    
                    #vi['GCA'][:,iS,:]=NetGrowthNew
    
            #------------------------------------------------------------------
            # Apply growth factors (in response to non-lethal events)
            #------------------------------------------------------------------

            flg=0
            if (flg==1) & (meta['Biomass Module']!='Sawtooth'):
                
                GF=vi['EH'][iS]['GrowthFactor'][indDist[iDist]]
                
                if (np.isnan(GF)==False) & (GF!=-999):
                    
                    GrowthFraction=1+GF.astype(float)/100               
                    NetGrowth=vi['GCA'][:,iS,:].copy().astype(float)/meta['GC']['Scale Factor']
                    
                    if ID_Type==meta['LUT Dist']['IDW']:
                        # Only a temporary change in growth (see severity class table)
                        ResponsePeriod=2  
                        A=np.arange(1,302,1)
                        iResponse=np.where( (A>=vo['A'][iT,iS]) & (A<=vo['A'][iT,iS]+ResponsePeriod) )[0]
                        NetGrowth[iResponse,:]=GrowthFraction*NetGrowth[iResponse,:]
                    else:
                        # A permanent change in growth
                        NetGrowth=GrowthFraction*NetGrowth                        
                    
                    NetGrowth=NetGrowth*meta['GC']['Scale Factor']
                    NetGrowth=NetGrowth.astype(np.int16)
                    vi['GCA'][:,iS,:]=NetGrowth
    
            #------------------------------------------------------------------
            # Aerial fertilization - growth response 
            #------------------------------------------------------------------
            
            if ID_Type==meta['LUT Dist']['Fertilization Aerial']:
                
                # Update N stimulus
                
                # Adjust net growth of aboveground biomass
                vi,vo,meta=update_nutrient_status(vi,vo,iT,iS,meta,psl,'AbovegroundNetGrowth')
                
                # Adjust emissions
                vi,vo,meta=update_nutrient_status(vi,vo,iT,iS,meta,psl,'Emissions')
    
    #--------------------------------------------------------------------------
    # Populate emission variables with values from ecosystem pool variable
    #--------------------------------------------------------------------------
    
    vo['C_E_FireAsCO2'][iT,:]=vo['C_Eco_Pools'][iT,:,iEP['ECO2asC']]        
    vo['C_E_FireAsCH4'][iT,:]=vo['C_Eco_Pools'][iT,:,iEP['ECH4asC']]
    vo['C_E_FireAsCO'][iT,:]=vo['C_Eco_Pools'][iT,:,iEP['ECOasC']]
    vo['C_E_FireAsN2O'][iT,:]=vo['C_Eco_Pools'][iT,:,iEP['EN2OasC']]
    vo['C_RemovedMerch'][iT,:]=vo['C_Eco_Pools'][iT,:,iEP['RemovedMerch']]
    vo['C_RemovedNonMerch'][iT,:]=vo['C_Eco_Pools'][iT,:,iEP['RemovedNonMerch']]
    vo['C_RemovedSnagStem'][iT,:]=vo['C_Eco_Pools'][iT,:,iEP['RemovedSnagStem']]
        
    return vo,vi