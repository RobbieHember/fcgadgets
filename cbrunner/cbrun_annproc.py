
#%% Import python modules

import numpy as np
from fcgadgets.cbrunner import cbrun_utilities as cbu
from fcgadgets.silviculture import nutrient_application as napp
from fcgadgets.taz import aspatial_stat_models as asm

#%% Biomass dynamics

def Biomass_FromTIPSYorTASS(iScn,iBat,iT,vi,vo,meta,iEP):
    
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
    NetGrowth=np.zeros((meta['Project']['Batch Size'][iBat],6))
    for iS in range(meta['Project']['Batch Size'][iBat]):
        NetGrowth[iS,:]=vi['GC']['Active'][iAge[iS],iS,:].copy().astype(float)
    
    # Apply growth factor
    NetGrowth=meta['GC']['Scale Factor']*NetGrowth
    
    # Net growth of total stemwood
    Gnet_Stem=NetGrowth[:,iEP['StemMerch']]+NetGrowth[:,iEP['StemNonMerch']]
    
    # Net growth of foliage    
    NetGrowth[:,iEP['Foliage']]=Gnet_Stem*(meta['Param']['BEV']['Biomass Allometry']['Gf1']+(meta['Param']['BEV']['Biomass Allometry']['Gf2']-meta['Param']['BEV']['Biomass Allometry']['Gf1'])*np.exp(-meta['Param']['BEV']['Biomass Allometry']['Gf3']*vo['A'][iT,:]))
    
    # Net growth of branches
    NetGrowth[:,iEP['Branch']]=Gnet_Stem*(meta['Param']['BEV']['Biomass Allometry']['Gbr1']+(meta['Param']['BEV']['Biomass Allometry']['Gbr2']-meta['Param']['BEV']['Biomass Allometry']['Gbr1'])*np.exp(-meta['Param']['BEV']['Biomass Allometry']['Gbr3']*vo['A'][iT,:]))
    
    # Net growth of bark
    NetGrowth[:,iEP['Bark']]=Gnet_Stem*(meta['Param']['BEV']['Biomass Allometry']['Gbk1']+(meta['Param']['BEV']['Biomass Allometry']['Gbk2']-meta['Param']['BEV']['Biomass Allometry']['Gbk1'])*np.exp(-meta['Param']['BEV']['Biomass Allometry']['Gbk3']*vo['A'][iT,:]))
    
    # Add net growth to output variable structure 
    # Oddly, using meta['iEP']['BiomassAboveground'] will invert the dimensions 
    # of C_G_Net - don't change it.
    vo['C_G_Net'][iT,:,0:5]=NetGrowth[:,0:5]
    
    # Total net growth of root biomass (Li et al. 2003, Eq. 4)
    G_Net_Root_Total=0.22*np.sum(vo['C_G_Net'][iT,:,0:5],axis=1)
    
    # Fine root fraction should decline with total biomass, but estimating it based on
    # size causes a lot of problems. Think about creating new equation.        
    vo['C_G_Net'][iT,:,iEP['RootCoarse']]=(1-0.072)*G_Net_Root_Total
    vo['C_G_Net'][iT,:,iEP['RootFine']]=0.072*G_Net_Root_Total
    
    # Update stemwood merchantable volume
    vo['V_StemMerch'][iT,:]=vo['V_StemMerch'][iT-1,:]+NetGrowth[:,5]#.astype(float)*meta['GC']['Scale Factor']      

    #--------------------------------------------------------------------------
    # Nutrient application effects to net growth
    #--------------------------------------------------------------------------
    
    # Index to stands that are stimulated by nutrient application
    meta['Nutrient Management']['iApplication']=np.where( meta['Nutrient Management']['ResponseCounter']>0 )[0]
    
    # Adjust N application response counter
    if meta['Nutrient Management']['iApplication'].size>0:
        vi,vo,meta=napp.UpdateStatus(vi,vo,iT,meta,'UpdateCounter')
    
    # Adjust root net growth
    if meta['Nutrient Management']['iApplication'].size>0:
        vi,vo,meta=napp.UpdateStatus(vi,vo,iT,meta,'BelowgroundNetGrowth')

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
    vo['C_M_Reg'][iT,:,0]=meta['Param']['BEV']['Biomass Turnover']['StemMerch'][0,:]*vo['C_Eco_Pools'][iT,:,0]
    vo['C_M_Reg'][iT,:,1]=meta['Param']['BEV']['Biomass Turnover']['StemMerch'][0,:]*vo['C_Eco_Pools'][iT,:,1]
    vo['C_M_Reg'][iT,:,2]=meta['Param']['BEV']['Biomass Turnover']['StemMerch'][0,:]*vo['C_Eco_Pools'][iT,:,2]
    vo['C_M_Reg'][iT,:,3]=meta['Param']['BEV']['Biomass Turnover']['StemMerch'][0,:]*vo['C_Eco_Pools'][iT,:,3]
    vo['C_M_Reg'][iT,:,4]=meta['Param']['BEV']['Biomass Turnover']['StemMerch'][0,:]*vo['C_Eco_Pools'][iT,:,4]
    vo['C_M_Reg'][iT,:,5]=meta['Param']['BEV']['Biomass Turnover']['StemMerch'][0,:]*vo['C_Eco_Pools'][iT,:,5]
    vo['C_M_Reg'][iT,:,6]=meta['Param']['BEV']['Biomass Turnover']['StemMerch'][0,:]*vo['C_Eco_Pools'][iT,:,6]
    
    # Adjust mortality to account for N application response
    if meta['Nutrient Management']['iApplication'].size>0:
        vi,vo,meta=napp.UpdateStatus(vi,vo,iT,meta,'Mortality')       
    
    # Apply optional factor to mortality
    #if meta['Scenario Switch']['Mortality Factor Status'][iScn]=='On':
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
        
        CToTransfer=np.zeros((meta['Project']['Batch Size'][iBat],7))
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
    
    # Setting turnover as a function of age will decouple NPP from net growth.
    fA=0
    #fA=-0.001    
    Aref=100
    
    # Calculate foliage biomass turnover due to litterfall
    tr=fA*(vo['A'][iT,:]-Aref)+meta['Param']['BEV']['Biomass Turnover']['Foliage'][0,:]
    vo['C_LF'][iT,:,iEP['Foliage']]=tr*vo['C_Eco_Pools'][iT,:,iEP['Foliage']]
    
    # Calculate branch biomass turnover due to litterfall
    tr=fA*(vo['A'][iT,:]-Aref)+meta['Param']['BEV']['Biomass Turnover']['Branch'][0,:]
    vo['C_LF'][iT,:,iEP['Branch']]=tr*vo['C_Eco_Pools'][iT,:,iEP['Branch']]
    
    # Calculate bark biomass turnover due to litterfall
    vo['C_LF'][iT,:,iEP['Bark']]=meta['Param']['BEV']['Biomass Turnover']['Bark'][0,:]*vo['C_Eco_Pools'][iT,:,iEP['Bark']]
    
    # Calculate coarse root biomass turnover due to litterfall
    tr=fA*(vo['A'][iT,:]-Aref)+meta['Param']['BEV']['Biomass Turnover']['RootCoarse'][0,:]
    vo['C_LF'][iT,:,iEP['RootCoarse']]=tr*vo['C_Eco_Pools'][iT,:,iEP['RootCoarse']]
    
    # Calculate fine root biomass turnover due to litterfall
    tr=fA*(vo['A'][iT,:]-Aref)+meta['Param']['BEV']['Biomass Turnover']['RootFine'][0,:]
    vo['C_LF'][iT,:,iEP['RootFine']]=tr*vo['C_Eco_Pools'][iT,:,iEP['RootFine']]
        
    # Adjust litterfall to account for N application response
    if meta['Nutrient Management']['iApplication'].size>0:
        vi,vo,meta=napp.UpdateStatus(vi,vo,iT,meta,'Litterfall')
    
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
    N_tree=meta['Project']['Max Num Trees Per Stand']
        
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
    #ptl=Bunch(meta['ptl'])
    
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
                        
            if vi['EH'][iS]['ID_Type'][iDist]==meta['LUT']['Dist']['Wildfire']:
                # Wildfire
                iKill=np.arange(0,nKill)
                tl_N_M_Inv_Fir[iT,iLive[iKill]]=1
                vo['A'][iT,iS]=0
                
            elif (vi['EH'][iS]['ID_Type'][iDist]==meta['LUT']['Dist']['Harvest']) | (vi['EH'][iS]['ID_Type'][iDist]==meta['LUT']['Dist']['Harvest and Slashpile Burn']):
                iKill=np.arange(0,nKill)
                tl_N_M_Inv_Har[iT,iLive[iKill]]=1
                vo['A'][iT,iS]=0
                               
            elif (vi['EH'][iS]['ID_Type'][iDist]==meta['LUT']['Dist']['Beetles']) | (vi['EH'][iS]['ID_Type'][iDist]==meta['LUT']['Dist']['Defoliators']):
                # Insects
                iKill=np.arange(0,nKill)
                tl_N_M_Inv_Ins[iT,iLive[iKill]]=1
                
            elif (vi['EH'][iS]['ID_Type'][iDist]==meta['LUT']['Dist']['Planting']):
                # Do nothing
                pass
                
            elif (vi['EH'][iS]['ID_Type'][iDist]==meta['LUT']['Dist']['Fertilization Aerial']):
                # Do nothing
                pass
            
            #elif vi['EH'][iS]['ID_Type'][iDist]==meta['LUT']['Dist']['Wind']:
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

def DOM_like_CBM08(iT,iBat,vi,vo,iEP,meta):
    
    # Extract parameters
    bIPF=meta['Param']['BEV']['Inter Pool Fluxes']
    bDec=meta['Param']['BEV']['Decomp']
    
    #--------------------------------------------------------------------------
    # Flux of carbon between biomass pools and dead organic matter pools
    #--------------------------------------------------------------------------
    
    # Transfer biomass turnover to very fast litter pool
    vo['C_Eco_Pools'][iT,:,iEP['LitterVF']]=vo['C_Eco_Pools'][iT-1,:,iEP['LitterVF']]+ \
        bIPF['FoliageLitToLitterVF']*vo['C_LF'][iT,:,iEP['Foliage']]+ \
        bIPF['RootFineLitToLitterVF']*vo['C_LF'][iT,:,iEP['RootFine']]+ \
        bIPF['FoliageMorToLitterVF']*vo['C_M_Reg'][iT,:,iEP['Foliage']]
    
    # Transfer biomass turnover to fast litter pool
    vo['C_Eco_Pools'][iT,:,iEP['LitterF']]=vo['C_Eco_Pools'][iT-1,:,iEP['LitterF']]+ \
        bIPF['BranchLitToLitterF']*vo['C_LF'][iT,:,iEP['Branch']]+ \
        bIPF['BarkLitToLitterF']*vo['C_LF'][iT,:,iEP['Bark']]+ \
        bIPF['RootCoarseLitToLitterF']*vo['C_LF'][iT,:,iEP['RootCoarse']]+ \
        bIPF['BranchMorToLitterF']*vo['C_M_Reg'][iT,:,iEP['Branch']]+ \
        bIPF['StemNonMerchLitToLitterF']*vo['C_LF'][iT,:,iEP['StemNonMerch']]
            
    # Transfer biomass turnover to medium litter pool
    vo['C_Eco_Pools'][iT,:,iEP['LitterM']]=vo['C_Eco_Pools'][iT-1,:,iEP['LitterM']]+ \
        bIPF['BarkLitToLitterM']*vo['C_LF'][iT,:,iEP['Bark']]+ \
        bIPF['BarkMorToLitterM']*vo['C_M_Reg'][iT,:,iEP['Bark']]

    # Transfer biomass turnover to slow litter pool
    vo['C_Eco_Pools'][iT,:,iEP['LitterS']]=vo['C_Eco_Pools'][iT-1,:,iEP['LitterS']]+ \
        bIPF['BarkLitToLitterS']*vo['C_LF'][iT,:,iEP['Bark']]+ \
        bIPF['BarkMorToLitterS']*vo['C_M_Reg'][iT,:,iEP['Bark']]
        
    # Transfer biomass turnover to snag stemwood pool
    vo['C_Eco_Pools'][iT,:,iEP['SnagStem']]=vo['C_Eco_Pools'][iT-1,:,iEP['SnagStem']]+ \
        bIPF['StemMerchMorToSnagStem']*vo['C_M_Reg'][iT,:,iEP['StemMerch']]+ \
        bIPF['StemNonMerchMorToSnagStem']*vo['C_M_Reg'][iT,:,iEP['StemNonMerch']]
    
    # Transfer biomass turnover to snag branch pool
    vo['C_Eco_Pools'][iT,:,iEP['SnagBranch']]=vo['C_Eco_Pools'][iT-1,:,iEP['SnagBranch']]+ \
        bIPF['BranchLitToSnagBranch']*vo['C_LF'][iT,:,iEP['Branch']]+ \
        bIPF['BranchMorToSnagBranch']*vo['C_M_Reg'][iT,:,iEP['Branch']]+ \
        bIPF['BarkMorToSnagBranch']*vo['C_M_Reg'][iT,:,iEP['Bark']]+ \
        bIPF['StemNonMerchLitToSnagBranch']*vo['C_LF'][iT,:,iEP['Bark']]

    # Transfer biomass turnover to very fast soil pool
    vo['C_Eco_Pools'][iT,:,iEP['SoilVF']]=vo['C_Eco_Pools'][iT-1,:,iEP['SoilVF']]+ \
        bIPF['RootFineLitToSoilVF']*vo['C_LF'][iT,:,iEP['RootFine']] +\
        bIPF['RootCoarseMorToSoilVF']*vo['C_M_Reg'][iT,:,iEP['RootCoarse']]+ \
        bIPF['RootFineMorToSoilVF']*vo['C_M_Reg'][iT,:,iEP['RootFine']]

    # Transfer biomass turnover to fast soil pool
    vo['C_Eco_Pools'][iT,:,iEP['SoilF']]=vo['C_Eco_Pools'][iT-1,:,iEP['SoilF']]+ \
        bIPF['RootCoarseLitToSoilF']*vo['C_LF'][iT,:,iEP['RootCoarse']]+ \
        bIPF['RootCoarseMorToSoilF']*vo['C_M_Reg'][iT,:,iEP['RootCoarse']]
        
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
    meta['R_LitterVF']=bDec['LitterVF_Rten']*vo['C_Eco_Pools'][iT,:,iEP['LitterVF']].flatten()*bDec['LitterVF_Qten']**fT
    meta['R_LitterF']=bDec['LitterF_Rten']*vo['C_Eco_Pools'][iT,:,iEP['LitterF']].flatten()*bDec['LitterF_Qten']**fT
    meta['R_LitterM']=bDec['LitterM_Rten']*vo['C_Eco_Pools'][iT,:,iEP['LitterM']].flatten()*bDec['LitterM_Qten']**fT
    meta['R_LitterS']=bDec['LitterS_Rten']*vo['C_Eco_Pools'][iT,:,iEP['LitterS']].flatten()*bDec['LitterS_Qten']**fT
    meta['R_SnagStem']=bDec['SnagStem_Rten']*vo['C_Eco_Pools'][iT,:,iEP['SnagStem']].flatten()*bDec['SnagStem_Qten']**fT
    meta['R_SnagBranch']=bDec['SnagBranch_Rten']*vo['C_Eco_Pools'][iT,:,iEP['SnagBranch']].flatten()*bDec['SnagBranch_Qten']**fT    
    meta['R_SoilVF']=bDec['SoilVF_Rten']*vo['C_Eco_Pools'][iT,:,iEP['SoilVF']].flatten()*bDec['SoilVF_Qten']**fT
    meta['R_SoilF']=bDec['SoilF_Rten']*vo['C_Eco_Pools'][iT,:,iEP['SoilF']].flatten()*bDec['SoilF_Qten']**fT
    meta['R_SoilS']=bDec['SoilS_Rten']*vo['C_Eco_Pools'][iT,:,iEP['SoilS']].flatten()*bDec['SoilS_Qten']**fT
    
    # Adjust decomposition to account for N application response
    if meta['Nutrient Management']['iApplication'].size>0:
        vi,vo,meta=napp.UpdateStatus(vi,vo,iT,meta,'HeterotrophicRespiration')
    
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
    vo['C_Eco_Pools'][iT,:,iEP['LitterF']]=vo['C_Eco_Pools'][iT,:,iEP['LitterF']]+bIPF['SnagBranchToLitterF']*meta['R_SnagBranch']
    
    # Re-define decayed medium litter
    vo['C_Eco_Pools'][iT,:,iEP['LitterM']]=vo['C_Eco_Pools'][iT,:,iEP['LitterM']]+bIPF['SnagStemToLitterM']*meta['R_SnagStem']
    
    # Re-define decayed slow litter
    vo['C_Eco_Pools'][iT,:,iEP['LitterS']]=vo['C_Eco_Pools'][iT,:,iEP['LitterS']]+bIPF['LitterVFToLitterS']*meta['R_LitterVF']
    vo['C_Eco_Pools'][iT,:,iEP['LitterS']]=vo['C_Eco_Pools'][iT,:,iEP['LitterS']]+bIPF['LitterFToLitterS']*meta['R_LitterF']
    vo['C_Eco_Pools'][iT,:,iEP['LitterS']]=vo['C_Eco_Pools'][iT,:,iEP['LitterS']]+bIPF['LitterMToLitterS']*meta['R_LitterM']
    
    # Re-define decayed slow soil
    vo['C_Eco_Pools'][iT,:,iEP['SoilS']]=vo['C_Eco_Pools'][iT,:,iEP['SoilS']]+bIPF['SoilVFToSoilS']*meta['R_SoilVF']
    vo['C_Eco_Pools'][iT,:,iEP['SoilS']]=vo['C_Eco_Pools'][iT,:,iEP['SoilS']]+bIPF['SoilFToSoilS']*meta['R_SoilF']
        
    # Heterotrophic respiration 
    vo['C_RH'][iT,:,iEP['LitterVF']]=bIPF['LitterVFToCO2']*meta['R_LitterVF']
    vo['C_RH'][iT,:,iEP['LitterF']]=bIPF['LitterFToCO2']*meta['R_LitterF']
    vo['C_RH'][iT,:,iEP['LitterM']]=bIPF['LitterMToCO2']*meta['R_LitterM']
    vo['C_RH'][iT,:,iEP['LitterS']]=bIPF['LitterSToCO2']*meta['R_LitterS']
    vo['C_RH'][iT,:,iEP['SnagStem']]=bIPF['SnagStemToCO2']*meta['R_SnagStem']
    vo['C_RH'][iT,:,iEP['SnagBranch']]=bIPF['SnagBranchToCO2']*meta['R_SnagBranch']
    vo['C_RH'][iT,:,iEP['SoilVF']]=bIPF['SoilVFToCO2']*meta['R_SoilVF']
    vo['C_RH'][iT,:,iEP['SoilF']]=bIPF['SoilFToCO2']*meta['R_SoilF']
    vo['C_RH'][iT,:,iEP['SoilS']]=bIPF['SoilSToCO2']*meta['R_SoilS']   
    
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
    PT_FelledStemMerch=bDec['FelledStemMerch_PhysTransRate']*vo['C_Eco_Pools'][iT,:,iEP['FelledStemMerch']]
    PT_FelledStemNonMerch=bDec['FelledStemNonMerch_PhysTransRate']*vo['C_Eco_Pools'][iT,:,iEP['FelledStemNonMerch']]
    PT_FelledBranch=bDec['FelledBranch_PhysTransRate']*vo['C_Eco_Pools'][iT,:,iEP['FelledBranch']]
    PT_FelledBark=bDec['FelledBark_PhysTransRate']*vo['C_Eco_Pools'][iT,:,iEP['FelledBark']]
    PT_FelledSnagStem=bDec['FelledSnagStem_PhysTransRate']*vo['C_Eco_Pools'][iT,:,iEP['FelledSnagStem']]
    PT_FelledSnagBranch=bDec['FelledSnagBranch_PhysTransRate']*vo['C_Eco_Pools'][iT,:,iEP['FelledSnagBranch']]
    
    PT_LitterS=bDec['LitterS_PhysTransRate']*vo['C_Eco_Pools'][iT,:,iEP['LitterS']]
    PT_SnagStem=bDec['SnagStem_PhysTransRate']*vo['C_Eco_Pools'][iT,:,iEP['SnagStem']]
    PT_SnagBranch=bDec['SnagBranch_PhysTransRate']*vo['C_Eco_Pools'][iT,:,iEP['SnagBranch']] 
    
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

#%% Disturbance and management events (from TAZ)

def Events_FromTaz(iT,iScn,iEns,iBat,vi,vo,meta,iEP):
    
    # Predict stand breakup (on the fly)
    if meta['Scenario'][iScn]['Breakup Status']=='On':
        vi=asm.PredictStandBreakup_OnTheFly(meta,vi,iT,iEns,vo['A'][iT,:])
    
    # Predict historical harvesting (on the fly)
    if meta['Scenario'][iScn]['Harvest Status Historical']=='On':
        if vi['tv'][iT]<meta['Scenario'][iScn]['Harvest Year Transition']:
            Period='Historical'
            Volume=vo['V_StemMerch'][iT,:]+2*(1/0.45)*vo['C_Eco_Pools'][iT,:,iEP['SnagStem']]
            vi=asm.PredictHarvesting_OnTheFly(meta,vi,iT,iScn,iEns,Volume,Period)
    
    # Predict future harvesting (on the fly)
    if meta['Scenario'][iScn]['Harvest Status Future']=='On':
        if vi['tv'][iT]>=meta['Scenario'][iScn]['Harvest Year Transition']:
            Period='Future'
            Volume=vo['V_StemMerch'][iT,:]+2*(1/0.45)*vo['C_Eco_Pools'][iT,:,iEP['SnagStem']]
            vi=asm.PredictHarvesting_OnTheFly(meta,vi,iT,iScn,iEns,Volume,Period) 
    
    # Predict future nutrient application (on the fly)
    if meta['Scenario'][iScn]['Nutrient Application Status']=='On':
        if vi['tv'][iT]>=meta['Project']['Year Project']:
            vi=napp.ScheduleApplication(meta,vi,vo,iT,iScn,iEns,iBat)
    
    # Initialize indicator of aerial nutrient application
    flag_nutrient_application=np.zeros(meta['Project']['Batch Size'][iBat])
    
    # Loop through events in year
    for iE in range(meta['Core']['Max Events Per Year']):
        
        # Event type IDs for the iE'th event of the year
        ID_Type=vi['EC']['ID_Type'][iT,:,iE].copy()
        
        # Record stands with aerial nutrient application
        iApp=np.where(ID_Type==meta['LUT']['Dist']['Fertilization Aerial'])[0]
        flag_nutrient_application[iApp]=1
        
        # Total affected biomass carbon
        MortalityFactor=vi['EC']['MortalityFactor'][iT,:,iE].copy()
        
        # Only continue if there are disturbance fluxes           
        #if (np.sum(MortalityFactor)==0) & (ID_Type!=meta['LUT']['Dist']['Planting']):
        #    # Don't do this 
        #    continue
        
        # Get event-specific parameters        
        u,idx,inv=np.unique(ID_Type,return_index=True,return_inverse=True)        
        b={}
        for k in meta['Param']['BEV']['Dist'][1].keys():
            bU=np.zeros(u.size)
            for iU in range(u.size):
                if u[iU]==0:
                    continue
                bU[iU]=meta['Param']['BEV']['Dist'][u[iU]][k]
            b[k]=bU[inv]
        
#        # *** Adjust utilization of merch stemwood on coast to reflect FPB study ***
#        # Need to work on regionalizing these parameters.
#        ind=np.where( (vi['Inv']['ID_BECZ']==meta['LUT']['VRI']['BEC_ZONE_CODE']['CWH']) | \
#                (vi['Inv']['ID_BECZ']==meta['LUT']['VRI']['BEC_ZONE_CODE']['CDF']) )[0]
#        if ind.size>0:
#            b['BiomassMerch_Removed'][ind]=0.88
#            b['BiomassMerch_LeftOnSite'][ind]=0.12
        
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
        Affected_BarkMerch=0.85*Affected_Bark
        Affected_BarkNonMerch=(1-0.85)*Affected_Bark
            
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
            
        # Mortality by category (lumping stands together)
        uType=np.unique(ID_Type)
        for iType in range(len(uType)):
            if uType[iType]==0:
                continue
            indType=np.where(ID_Type==uType[iType])[0]
            String_Type=cbu.lut_n2s(meta['LUT']['Dist'],uType[iType])[0]
            vo['C_M_ByAgent'][String_Type][iT,0]=vo['C_M_ByAgent'][String_Type][iT,0]+np.sum(Affected_All[indType])

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
        vo['V_StemMerch'][iT,:]=np.maximum(0,vo['V_StemMerch'][iT,:]-Affected_VolumeStemMerch)
        #vo['V_StemMerch'][iT,:]=vo['V_StemMerch'][iT,:]-Affected_VolumeStemMerch
            
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
        
        vo['C_E_FireAsCO2'][iT,:]=vo['C_E_FireAsCO2'][iT,:]+meta['Param']['BEV']['Biophysical']['CombFrac_CO2']*Total_Burned
        vo['C_E_FireAsCH4'][iT,:]=vo['C_E_FireAsCH4'][iT,:]+meta['Param']['BEV']['Biophysical']['CombFrac_CH4']*Total_Burned
        vo['C_E_FireAsCO'][iT,:]=vo['C_E_FireAsCO'][iT,:]+meta['Param']['BEV']['Biophysical']['CombFrac_CO']*Total_Burned
        
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
        Total_Burned=StemMerch+StemNonMerch+Branch+Bark+SnagStem+SnagBranch
        
        # Remove affected carbon
        vo['C_Eco_Pools'][iT,:,iEP['FelledStemMerch']]=vo['C_Eco_Pools'][iT,:,iEP['FelledStemMerch']]-StemMerch
        vo['C_Eco_Pools'][iT,:,iEP['FelledStemNonMerch']]=vo['C_Eco_Pools'][iT,:,iEP['FelledStemNonMerch']]-StemNonMerch
        vo['C_Eco_Pools'][iT,:,iEP['FelledBranch']]=vo['C_Eco_Pools'][iT,:,iEP['FelledBranch']]-Branch
        vo['C_Eco_Pools'][iT,:,iEP['FelledBark']]=vo['C_Eco_Pools'][iT,:,iEP['FelledBark']]-Bark
        vo['C_Eco_Pools'][iT,:,iEP['FelledSnagStem']]=vo['C_Eco_Pools'][iT,:,iEP['FelledSnagStem']]-SnagStem
        vo['C_Eco_Pools'][iT,:,iEP['FelledSnagBranch']]=vo['C_Eco_Pools'][iT,:,iEP['FelledSnagBranch']]-SnagBranch
        
        # Add to fire emissions
        vo['C_E_FireAsCO2'][iT,:]=vo['C_E_FireAsCO2'][iT,:]+meta['Param']['BEV']['Biophysical']['CombFrac_CO2']*Total_Burned
        vo['C_E_FireAsCH4'][iT,:]=vo['C_E_FireAsCH4'][iT,:]+meta['Param']['BEV']['Biophysical']['CombFrac_CH4']*Total_Burned
        vo['C_E_FireAsCO'][iT,:]=vo['C_E_FireAsCO'][iT,:]+meta['Param']['BEV']['Biophysical']['CombFrac_CO']*Total_Burned
        
        #----------------------------------------------------------------------
        # Update stand age
        #----------------------------------------------------------------------
            
        if (meta['Project']['Biomass Module']!='Sawtooth') & (meta['Project']['Partial Mortality Affects Age']=='On'):
            
            # List of exception event types (that will remain at the same age)
            Exceptions_to_Partial_Mortality=[  ]
            #v2[iB1]['A']

            # Index to types where age will change
            ind=np.where( (np.isin(ID_Type,Exceptions_to_Partial_Mortality)==False) )[0]
                
            # Assume oldest trees were most affected, reduce age in prportion
            # with mortality rate
            # Not always realistic, but see how net growth is affected.
            vo['A'][iT,ind]=vo['A'][iT,ind]*(1-MortalityFactor)
        
        # Ensure stand-replacing events reset stand age to 0
        vo['A'][iT,(MortalityFactor==1)]=0
        
        # Ensure planting resets to age 0                
        vo['A'][iT,(ID_Type==meta['LUT']['Dist']['Planting']) | (ID_Type==meta['LUT']['Dist']['Direct Seeding'])]=0
            
        #----------------------------------------------------------------------
        # Transition to new growth curve
        #----------------------------------------------------------------------
            
        # Initialize a flag that indicates whether the growth curve changes
        #flg_gc_change=0
        #if vi['EH'][iS]['ID_GrowthCurveM'][indDist[iDist]]!=vi['ID_GCA'][iS]:
        #    flg_gc_change=1
            
        # Only applies to BatchTIPSY
        if meta['Project']['Biomass Module']=='BatchTIPSY':
            
            for iGC in range(meta['GC']['N Growth Curves']):
                
                # Don't alter growth curve for fertilization
                ind=np.where( (vi['EC']['ID_GrowthCurve'][iT,:,iE]==meta['GC']['ID GC Unique'][iGC]) & (ID_Type!=meta['LUT']['Dist']['Fertilization Aerial']) )[0]
                
                if ind.size>0:
                    vi['GC']['Active'][:,ind,:]=vi['GC'][ meta['GC']['ID GC Unique'][iGC] ][:,ind,:]
                    vi['GC']['ID_GCA'][ind]=int(meta['GC']['ID GC Unique'][iGC])
    
        #----------------------------------------------------------------------
        # Impose regen failure
        #----------------------------------------------------------------------
            
        # Only applies to BatchTIPSY
        if meta['Project']['Biomass Module']=='BatchTIPSY':
            iFailure=np.where(ID_Type==meta['LUT']['Dist']['Regen Failure'])[0]
            if iFailure.size>0:
                vi['GC']['Active'][:,iFailure,:]=0
                vi['GC'][1][:,iFailure,:]=0
    
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
        
        # *** This is not up to date - will crash if attempted - needs update. ***
        flg=0
        
        if (flg==1) & (meta['Project']['Biomass Module']!='Sawtooth'):
                
            # Half life
            hl=meta['Param']['BEV']['bDist_GrowthRecovery_HL']
                
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
                    
                vi['GCA'][:,iS,:]=NetGrowthNew
    
        #----------------------------------------------------------------------
        # Apply growth factors (in response to non-lethal events)
        #----------------------------------------------------------------------

        flg=1
        if (flg==1) & (meta['Project']['Biomass Module']!='Sawtooth'):
                
            GF=vi['EC']['GrowthFactor'][iT,:,iE]
            
            indAdj=np.where(GF!=0)[0]
            
            if (indAdj.size>0):
                
                #print(vi['GC']['Active'].shape)
                
                GrowthFraction=1+GF[indAdj].astype(float)/100               
                NetGrowth=vi['GC']['Active'][:,indAdj,:].copy().astype(float)/meta['GC']['Scale Factor']
                    
                #if ID_Type==meta['LUT']['Dist']['IDW']:
                #    # Only a temporary change in growth (see severity class table)
                #    # *** I have not checked to see if this is working properly ***
                #    ResponsePeriod=2  
                #    A=np.arange(1,302,1)
                #    iResponse=np.where( (A>=vo['A'][iT,indAdj]) & (A<=vo['A'][iT,indAdj]+ResponsePeriod) )[0]
                #    NetGrowth[iResponse,:]=GrowthFraction*NetGrowth[iResponse,:]
                #else:
                
                # A permanent change in growth
                for iAdj in range(indAdj.size):
                    NetGrowth[:,iAdj,:]=GrowthFraction[iAdj]*NetGrowth[:,iAdj,:]                        
                    
                NetGrowth=NetGrowth*meta['GC']['Scale Factor']
                NetGrowth=NetGrowth.astype(np.int16)
                vi['GC']['Active'][:,indAdj,:]=NetGrowth

    #--------------------------------------------------------------------------
    # Aerial nutrient application events
    #--------------------------------------------------------------------------
    
    # Generate index to stands that were fertilized
    meta['Nutrient Management']['iApplication']=np.where(flag_nutrient_application==1)[0]    
        
    if meta['Nutrient Management']['iApplication'].size>0:
        
        # Adjust net growth of aboveground biomass
        vi,vo,meta=napp.UpdateStatus(vi,vo,iT,meta,'AbovegroundNetGrowth')
                
        # Adjust emissions
        vi,vo,meta=napp.UpdateStatus(vi,vo,iT,meta,'Emissions')
        
    return vo,vi

#%% Harvested wood products sector (from Dymond 2012)

def HWP_From_BCHWP12(iT,iBat,vi,vo,meta):
    
    #--------------------------------------------------------------------------
    # If custom harvests occur, revise parameters
    # *** This is slow, but thus far, custom harvests have only been used in
    # small projects. It will require improvements for speed if applied in big 
    # projects. It doesn't need to be in the annual loop, for starters. ***
    #--------------------------------------------------------------------------
    
    if 'Harvest Custom' in meta:
        
        for iE in range(meta['Core']['Max Events Per Year']):
            
            ID_Type=vi['EC']['ID_Type'][iT,:,iE]
        
            for iHC in range(10):
                
                ID_HC=int(iHC+1)
                ID_Dist_HC=meta['LUT']['Dist']['Harvest Custom ' + str(ID_HC)]
                
                if ID_HC in meta['Harvest Custom']:
        
                    iStands=np.where(ID_Type==ID_Dist_HC)[0]
                    
                    if iStands.size==0:
                        continue
                    
                    for k in meta['Harvest Custom'][ID_HC].keys():
                        
                        if k[0:7]=='Removed':
                            meta['Param']['BEV']['HWP'][k][iStands]=meta['Harvest Custom'][ID_HC][k]/100
    
    #--------------------------------------------------------------------------
    # Ecosystems --> Mills
    # Note: Mills act as transient reservoirs that are not tracked over time
    #--------------------------------------------------------------------------
    
    Mill_FuelFromChips=meta['Param']['BEV']['HWP']['RemovedMerchToFuel']*vo['C_RemovedMerch'][iT,:] + \
                 meta['Param']['BEV']['HWP']['RemovedNonMerchToFuel']*vo['C_RemovedNonMerch'][iT,:] + \
                 meta['Param']['BEV']['HWP']['RemovedSnagStemToFuel']*vo['C_RemovedSnagStem'][iT,:]    
    
    Mill_Lumber=meta['Param']['BEV']['HWP']['RemovedMerchToLumber']*vo['C_RemovedMerch'][iT,:] + \
                 meta['Param']['BEV']['HWP']['RemovedNonMerchToLumber']*vo['C_RemovedNonMerch'][iT,:] + \
                 meta['Param']['BEV']['HWP']['RemovedSnagStemToLumber']*vo['C_RemovedSnagStem'][iT,:]
    
    Mill_Plywood=meta['Param']['BEV']['HWP']['RemovedMerchToPlywood']*vo['C_RemovedMerch'][iT,:] + \
                 meta['Param']['BEV']['HWP']['RemovedNonMerchToPlywood']*vo['C_RemovedNonMerch'][iT,:] + \
                 meta['Param']['BEV']['HWP']['RemovedSnagStemToPlywood']*vo['C_RemovedSnagStem'][iT,:]
    
    Mill_PulpFromChips=meta['Param']['BEV']['HWP']['RemovedMerchToPulp']*vo['C_RemovedMerch'][iT,:] + \
                 meta['Param']['BEV']['HWP']['RemovedNonMerchToPulp']*vo['C_RemovedNonMerch'][iT,:] + \
                 meta['Param']['BEV']['HWP']['RemovedSnagStemToPulp']*vo['C_RemovedSnagStem'][iT,:]
    
    Mill_Pulp=Mill_PulpFromChips + \
                 meta['Param']['BEV']['HWP']['LumberToPulp']*Mill_Lumber + \
                 meta['Param']['BEV']['HWP']['PlywoodToPulp']*Mill_Plywood
    
    Mill_OSB=meta['Param']['BEV']['HWP']['RemovedMerchToOSB']*vo['C_RemovedMerch'][iT,:] + \
                 meta['Param']['BEV']['HWP']['RemovedNonMerchToOSB']*vo['C_RemovedNonMerch'][iT,:] + \
                 meta['Param']['BEV']['HWP']['RemovedSnagStemToOSB']*vo['C_RemovedSnagStem'][iT,:]
    
    Mill_MDF=meta['Param']['BEV']['HWP']['RemovedMerchToMDF']*vo['C_RemovedMerch'][iT,:] + \
                 meta['Param']['BEV']['HWP']['RemovedNonMerchToMDF']*vo['C_RemovedNonMerch'][iT,:] + \
                 meta['Param']['BEV']['HWP']['RemovedSnagStemToMDF']*vo['C_RemovedSnagStem'][iT,:]
    
    Mill_Firewood=meta['Param']['BEV']['HWP']['RemovedMerchToFirewood']*vo['C_RemovedMerch'][iT,:] + \
                 meta['Param']['BEV']['HWP']['RemovedNonMerchToFirewood']*vo['C_RemovedNonMerch'][iT,:] + \
                 meta['Param']['BEV']['HWP']['RemovedSnagStemToFirewood']*vo['C_RemovedSnagStem'][iT,:]
    
    Mill_Cants=meta['Param']['BEV']['HWP']['RemovedMerchToCants']*vo['C_RemovedMerch'][iT,:] + \
                 meta['Param']['BEV']['HWP']['RemovedNonMerchToCants']*vo['C_RemovedNonMerch'][iT,:] + \
                 meta['Param']['BEV']['HWP']['RemovedSnagStemToCants']*vo['C_RemovedSnagStem'][iT,:]    
    
    # Transfers of lumber from lumber mill
    Mill_LumberToSFH=Mill_Lumber*meta['Param']['BEV']['HWP']['LumberToSFH']
    Mill_LumberToMFH=Mill_Lumber*meta['Param']['BEV']['HWP']['LumberToMFH']
    Mill_LumberToCom=Mill_Lumber*meta['Param']['BEV']['HWP']['LumberToCom']
    Mill_LumberToFurn=Mill_Lumber*meta['Param']['BEV']['HWP']['LumberToFurn']
    Mill_LumberToShip=Mill_Lumber*meta['Param']['BEV']['HWP']['LumberToShip']
    Mill_LumberToRepairs=Mill_Lumber*meta['Param']['BEV']['HWP']['LumberToRepairs']
    Mill_LumberToOther=Mill_Lumber*meta['Param']['BEV']['HWP']['LumberToOther']
    #Mill_LumberToPulp=Mill_Lumber*meta['Param']['BEV']['HWP']['LumberToPulp']
    Mill_LumberToFuel=Mill_Lumber*meta['Param']['BEV']['HWP']['LumberToFuel']
    
    # Transfers of plywood from plywood mill
    Mill_PlywoodToSFH=Mill_Plywood*meta['Param']['BEV']['HWP']['PlywoodToSFH']
    Mill_PlywoodToMFH=Mill_Plywood*meta['Param']['BEV']['HWP']['PlywoodToMFH']
    Mill_PlywoodToCom=Mill_Plywood*meta['Param']['BEV']['HWP']['PlywoodToCom']
    Mill_PlywoodToFurn=Mill_Plywood*meta['Param']['BEV']['HWP']['PlywoodToFurn']
    Mill_PlywoodToShip=Mill_Plywood*meta['Param']['BEV']['HWP']['PlywoodToShip']
    Mill_PlywoodToRepairs=Mill_Plywood*meta['Param']['BEV']['HWP']['PlywoodToRepairs']
    Mill_PlywoodToOther=Mill_Plywood*meta['Param']['BEV']['HWP']['PlywoodToOther']
    #Mill_PlywoodToPulp=Mill_Plywood*meta['Param']['BEV']['HWP']['PlywoodToPulp']
    Mill_PlywoodToFuel=Mill_Plywood*meta['Param']['BEV']['HWP']['PlywoodToFuel']
    
    # Transfers of OSB from OSB mill
    Mill_OSBToSFH=Mill_OSB*meta['Param']['BEV']['HWP']['OSBToSFH']
    Mill_OSBToMFH=Mill_OSB*meta['Param']['BEV']['HWP']['OSBToMFH']
    Mill_OSBToCom=Mill_OSB*meta['Param']['BEV']['HWP']['OSBToCom']
    Mill_OSBToFurn=Mill_OSB*meta['Param']['BEV']['HWP']['OSBToFurn']
    Mill_OSBToShip=Mill_OSB*meta['Param']['BEV']['HWP']['OSBToShip']
    Mill_OSBToRepairs=Mill_OSB*meta['Param']['BEV']['HWP']['OSBToRepairs']
    Mill_OSBToOther=Mill_OSB*meta['Param']['BEV']['HWP']['OSBToOther']
    
    # Transfers of MDF from MDF mill
    Mill_MDFToSFH=Mill_MDF*meta['Param']['BEV']['HWP']['MDFToSFH']
    Mill_MDFToMFH=Mill_MDF*meta['Param']['BEV']['HWP']['MDFToMFH']
    Mill_MDFToCom=Mill_MDF*meta['Param']['BEV']['HWP']['MDFToCom']
    Mill_MDFToFurn=Mill_MDF*meta['Param']['BEV']['HWP']['MDFToFurn']
    Mill_MDFToShip=Mill_MDF*meta['Param']['BEV']['HWP']['MDFToShip']
    Mill_MDFToRepairs=Mill_MDF*meta['Param']['BEV']['HWP']['MDFToRepairs']
    Mill_MDFToOther=Mill_MDF*meta['Param']['BEV']['HWP']['MDFToOther']
    
    #--------------------------------------------------------------------------
    # Track sales (for economic modelling)
    #--------------------------------------------------------------------------
    
    Lumber_Total=Mill_LumberToSFH+Mill_LumberToMFH+Mill_LumberToCom+Mill_LumberToFurn+Mill_LumberToShip+Mill_LumberToRepairs+Mill_LumberToOther
    Plywood_Total=Mill_PlywoodToSFH+Mill_PlywoodToMFH+Mill_PlywoodToCom+Mill_PlywoodToFurn+Mill_PlywoodToShip+Mill_PlywoodToRepairs+Mill_PlywoodToOther
    OSB_Total=Mill_OSBToSFH+Mill_OSBToMFH+Mill_OSBToCom+Mill_OSBToFurn+Mill_OSBToShip+Mill_OSBToRepairs+Mill_OSBToOther
    MDF_Total=Mill_MDFToSFH+Mill_MDFToMFH+Mill_MDFToCom+Mill_MDFToFurn+Mill_MDFToShip+Mill_MDFToRepairs+Mill_MDFToOther
    
    Paper_Total=meta['Param']['BEV']['HWP']['PulpToPaper']*Mill_Pulp
    Fuel_Total=Mill_FuelFromChips+Mill_LumberToFuel+Mill_PlywoodToFuel  
    
    vo['C_Lumber'][iT,:]=Lumber_Total
    vo['C_Plywood'][iT,:]=Plywood_Total
    vo['C_OSB'][iT,:]=OSB_Total
    vo['C_MDF'][iT,:]=MDF_Total
    vo['C_Paper'][iT,:]=Paper_Total
    vo['C_Fuel'][iT,:]=Fuel_Total
    
    #--------------------------------------------------------------------------
    # Mills --> In-use pools or other mills
    #--------------------------------------------------------------------------
        
    # Transfer mill fibre to single-family homes
    ip=meta['Core']['iPP']['SFH']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT-1,:,ip] + \
        Mill_LumberToSFH + \
        Mill_PlywoodToSFH + \
        Mill_OSBToSFH + \
        Mill_MDFToSFH
    
    # Transfer mill fibre to multi-family homes
    ip=meta['Core']['iPP']['MFH']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT-1,:,ip] + \
        Mill_LumberToMFH + \
        Mill_PlywoodToMFH + \
        Mill_OSBToMFH + \
        Mill_MDFToMFH
    
    # Transfer mill fibre to commercial
    ip=meta['Core']['iPP']['Comm']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT-1,:,ip] + \
        Mill_LumberToCom + \
        Mill_PlywoodToCom + \
        Mill_OSBToCom + \
        Mill_MDFToCom
    
    # Transfer mill fibre to furniture
    ip=meta['Core']['iPP']['Furn']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT-1,:,ip] + \
        Mill_LumberToFurn + \
        Mill_PlywoodToFurn + \
        Mill_OSBToFurn + \
        Mill_MDFToFurn
    
    # Transfer mill fibre to shipping
    ip=meta['Core']['iPP']['Ship']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT-1,:,ip] + \
        Mill_LumberToShip + \
        Mill_PlywoodToShip + \
        Mill_OSBToShip + \
        Mill_MDFToShip
    
    # Transfer mill fibre to repairs
    ip=meta['Core']['iPP']['Repairs']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT-1,:,ip] + \
        Mill_LumberToRepairs + \
        Mill_PlywoodToRepairs + \
        Mill_OSBToRepairs + \
        Mill_MDFToRepairs
    
    # Transfer mill fibre to other
    ip=meta['Core']['iPP']['Other']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT-1,:,ip] + \
        Mill_LumberToOther + \
        Mill_PlywoodToOther + \
        Mill_OSBToOther + \
        Mill_MDFToOther
    
    # Transfer pulp mill fibre to paper
    ip=meta['Core']['iPP']['Paper']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT-1,:,ip] + Paper_Total
    
    # Transfer mill fibre to fuel 
    ip=meta['Core']['iPP']['Fuel']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT-1,:,ip] + Fuel_Total
    
    # Transfer firewood to firewood pool
    ip=meta['Core']['iPP']['Firewood']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT-1,:,ip] + Mill_Firewood
    
    # Transfer mill fibre to cants
    ip=meta['Core']['iPP']['Cants']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT-1,:,ip] + Mill_Cants
    
    # Transfer pulp mill carbon to pulp-mill effluent
    ip=meta['Core']['iPP']['EffluentPulp']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT-1,:,ip] + \
                 meta['Param']['BEV']['HWP']['PulpToEffluent']*Mill_Pulp   
    
    #--------------------------------------------------------------------------
    # Update dump and landfill reservoirs
    #--------------------------------------------------------------------------

    vo['C_Pro_Pools'][iT,:,meta['Core']['iPP']['DumpWood']]=vo['C_Pro_Pools'][iT-1,:,meta['Core']['iPP']['DumpWood']]
    vo['C_Pro_Pools'][iT,:,meta['Core']['iPP']['DumpPaper']]=vo['C_Pro_Pools'][iT-1,:,meta['Core']['iPP']['DumpPaper']]
    vo['C_Pro_Pools'][iT,:,meta['Core']['iPP']['LandfillWoodDegradable']]=vo['C_Pro_Pools'][iT-1,:,meta['Core']['iPP']['LandfillWoodDegradable']]
    vo['C_Pro_Pools'][iT,:,meta['Core']['iPP']['LandfillWoodNonDegradable']]=vo['C_Pro_Pools'][iT-1,:,meta['Core']['iPP']['LandfillWoodNonDegradable']]    
    vo['C_Pro_Pools'][iT,:,meta['Core']['iPP']['LandfillPaperDegradable']]=vo['C_Pro_Pools'][iT-1,:,meta['Core']['iPP']['LandfillPaperDegradable']]
    vo['C_Pro_Pools'][iT,:,meta['Core']['iPP']['LandfillPaperNonDegradable']]=vo['C_Pro_Pools'][iT-1,:,meta['Core']['iPP']['LandfillPaperNonDegradable']]    
        
    #--------------------------------------------------------------------------
    # Single-family homes --> dump and landfill
    #--------------------------------------------------------------------------
           
    # Turnover
    ip=meta['Core']['iPP']['SFH']
    C_retired=meta['Param']['BEV']['HWP']['SFH_tr']*vo['C_Pro_Pools'][iT-1,:,ip]
    
    # Remove carbon
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] - C_retired
    
    # Transfer carbon to dump wood
    ip=meta['Core']['iPP']['DumpWood']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + meta['Param']['BEV']['HWP']['SFHToDumpWood']*C_retired
    
    # Transfer carbon to landfill (degradable)
    ip=meta['Core']['iPP']['LandfillWoodDegradable']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + meta['Param']['BEV']['HWP']['SFHToLandfillWood']*meta['Param']['BEV']['HWP']['ToLandfillWoodDegradableFrac']*C_retired
    
    # Transfer carbon to landfill (non-degradable)
    ip=meta['Core']['iPP']['LandfillWoodNonDegradable']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + meta['Param']['BEV']['HWP']['SFHToLandfillWood']*(1-meta['Param']['BEV']['HWP']['ToLandfillWoodDegradableFrac'])*C_retired
    
    #--------------------------------------------------------------------------
    # Multi-family homes --> dump and landfill
    #--------------------------------------------------------------------------
    
    # Turnover
    ip=meta['Core']['iPP']['MFH']
    C_retired=meta['Param']['BEV']['HWP']['MFH_tr']*vo['C_Pro_Pools'][iT,:,ip]
    
    # Remove carbon
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] - C_retired
    
    # Transfer carbon to dump wood
    ip=meta['Core']['iPP']['DumpWood']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + meta['Param']['BEV']['HWP']['MFHToDumpWood']*C_retired
    
    # Transfer carbon to landfill (degradable)
    ip=meta['Core']['iPP']['LandfillWoodDegradable']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + meta['Param']['BEV']['HWP']['MFHToLandfillWood']*meta['Param']['BEV']['HWP']['ToLandfillWoodDegradableFrac']*C_retired
    
    # Transfer carbon to landfill (non-degradable)
    ip=meta['Core']['iPP']['LandfillWoodNonDegradable']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + meta['Param']['BEV']['HWP']['MFHToLandfillWood']*(1-meta['Param']['BEV']['HWP']['ToLandfillWoodDegradableFrac'])*C_retired
    
    #--------------------------------------------------------------------------
    # Commercial building --> dump and landfill
    #--------------------------------------------------------------------------
    
    # Turnover
    ip=meta['Core']['iPP']['Comm']
    C_retired=meta['Param']['BEV']['HWP']['Comm_tr']*vo['C_Pro_Pools'][iT-1,:,ip]
    
    # Remove carbon
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] - C_retired
    
    # Transfer carbon to dump wood
    ip=meta['Core']['iPP']['DumpWood']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + meta['Param']['BEV']['HWP']['CommToDumpWood']*C_retired
    
    # Transfer carbon to landfill (degradable)
    ip=meta['Core']['iPP']['LandfillWoodDegradable']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + meta['Param']['BEV']['HWP']['CommToLandfillWood']*meta['Param']['BEV']['HWP']['ToLandfillWoodDegradableFrac']*C_retired
    
    # Transfer carbon to landfill (non-degradable)
    ip=meta['Core']['iPP']['LandfillWoodNonDegradable']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + meta['Param']['BEV']['HWP']['CommToLandfillWood']*(1-meta['Param']['BEV']['HWP']['ToLandfillWoodDegradableFrac'])*C_retired
    
    #--------------------------------------------------------------------------
    # Furniture --> dump and landfill
    #--------------------------------------------------------------------------
    
    # Turnover
    ip=meta['Core']['iPP']['Furn']
    C_retired=meta['Param']['BEV']['HWP']['Furn_tr']*vo['C_Pro_Pools'][iT-1,:,ip]
    
    # Remove carbon
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] - C_retired
    
    # Transfer carbon to dump wood
    ip=meta['Core']['iPP']['DumpWood']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + meta['Param']['BEV']['HWP']['FurnToDumpWood']*C_retired
    
    # Transfer carbon to landfill (degradable)
    ip=meta['Core']['iPP']['LandfillWoodDegradable']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + meta['Param']['BEV']['HWP']['FurnToLandfillWood']*meta['Param']['BEV']['HWP']['ToLandfillWoodDegradableFrac']*C_retired
    
    # Transfer carbon to landfill (non-degradable)
    ip=meta['Core']['iPP']['LandfillWoodNonDegradable']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + meta['Param']['BEV']['HWP']['FurnToLandfillWood']*(1-meta['Param']['BEV']['HWP']['ToLandfillWoodDegradableFrac'])*C_retired
    
    
    #--------------------------------------------------------------------------
    # Shipping --> dump and landfill
    #--------------------------------------------------------------------------
    
    # Turnover
    ip=meta['Core']['iPP']['Ship']
    C_retired=meta['Param']['BEV']['HWP']['Ship_tr']*vo['C_Pro_Pools'][iT-1,:,ip]
    
    # Remove carbon
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] - C_retired
    
    # Transfer carbon to dump wood
    ip=meta['Core']['iPP']['DumpWood']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + meta['Param']['BEV']['HWP']['ShipToDumpWood']*C_retired
    
    # Transfer carbon to landfill (degradable)
    ip=meta['Core']['iPP']['LandfillWoodDegradable']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + meta['Param']['BEV']['HWP']['ShipToLandfillWood']*meta['Param']['BEV']['HWP']['ToLandfillWoodDegradableFrac']*C_retired
    
    # Transfer carbon to landfill (non-degradable)
    ip=meta['Core']['iPP']['LandfillWoodNonDegradable']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + meta['Param']['BEV']['HWP']['ShipToLandfillWood']*(1-meta['Param']['BEV']['HWP']['ToLandfillWoodDegradableFrac'])*C_retired
    
    #--------------------------------------------------------------------------
    # Repairs --> dump and landfill
    #--------------------------------------------------------------------------
    
    # Turnover
    ip=meta['Core']['iPP']['Repairs']
    C_retired=meta['Param']['BEV']['HWP']['Repairs_tr']*vo['C_Pro_Pools'][iT-1,:,ip]
    
    # Remove carbon
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] - C_retired
    
    # Transfer carbon to dump wood
    ip=meta['Core']['iPP']['DumpWood']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + meta['Param']['BEV']['HWP']['RepairsToDumpWood']*C_retired
    
    # Transfer carbon to landfill (degradble)
    ip=meta['Core']['iPP']['LandfillWoodDegradable']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + meta['Param']['BEV']['HWP']['RepairsToLandfillWood']*meta['Param']['BEV']['HWP']['ToLandfillWoodDegradableFrac']*C_retired
    
    # Transfer carbon to landfill (non-degradable)
    ip=meta['Core']['iPP']['LandfillWoodNonDegradable']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + meta['Param']['BEV']['HWP']['RepairsToLandfillWood']*(1-meta['Param']['BEV']['HWP']['ToLandfillWoodDegradableFrac'])*C_retired
    
    #--------------------------------------------------------------------------
    # Other --> dump and landfill
    #--------------------------------------------------------------------------
    
    # Turnover
    ip=meta['Core']['iPP']['Other']
    C_retired=meta['Param']['BEV']['HWP']['Other_tr']*vo['C_Pro_Pools'][iT-1,:,ip]
    
    # Remove carbon
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] - C_retired
    
    # Transfer carbon to dump wood
    ip=meta['Core']['iPP']['DumpWood']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + meta['Param']['BEV']['HWP']['OtherToDumpWood']*C_retired
    
    # Transfer carbon to landfill (degradble)
    ip=meta['Core']['iPP']['LandfillWoodDegradable']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + meta['Param']['BEV']['HWP']['OtherToLandfillWood']*meta['Param']['BEV']['HWP']['ToLandfillWoodDegradableFrac']*C_retired
    
    # Transfer carbon to landfill (non-degradable)
    ip=meta['Core']['iPP']['LandfillWoodNonDegradable']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + meta['Param']['BEV']['HWP']['OtherToLandfillWood']*(1-meta['Param']['BEV']['HWP']['ToLandfillWoodDegradableFrac'])*C_retired
    
    #--------------------------------------------------------------------------
    # Cants --> dump and landfill
    #--------------------------------------------------------------------------
    
    # Turnover
    ip=meta['Core']['iPP']['Cants']
    C_retired=meta['Param']['BEV']['HWP']['Cants_tr']*vo['C_Pro_Pools'][iT-1,:,ip]
    
    # Remove carbon
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] - C_retired
    
    # Transfer carbon to dump wood
    ip=meta['Core']['iPP']['DumpWood']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + meta['Param']['BEV']['HWP']['CantsToDumpWood']*C_retired
    
    # Transfer carbon to landfill (degradble)
    ip=meta['Core']['iPP']['LandfillWoodDegradable']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + meta['Param']['BEV']['HWP']['CantsToLandfillWood']*meta['Param']['BEV']['HWP']['ToLandfillWoodDegradableFrac']*C_retired
    
    # Transfer carbon to landfill (non-degradable)
    ip=meta['Core']['iPP']['LandfillWoodNonDegradable']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + meta['Param']['BEV']['HWP']['CantsToLandfillWood']*(1-meta['Param']['BEV']['HWP']['ToLandfillWoodDegradableFrac'])*C_retired
    
    #--------------------------------------------------------------------------
    # Paper --> dump and landfill
    #--------------------------------------------------------------------------
    
    # Turnover (with adjustment for recycling)
    ip=meta['Core']['iPP']['Paper']
    C_retired=(1-meta['Param']['BEV']['HWP']['PaperRecycleRate'])*meta['Param']['BEV']['HWP']['Paper_tr']*vo['C_Pro_Pools'][iT,:,ip]
    
    # Remove carbon
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] - C_retired
    
    # Transfer to dump
    ip=meta['Core']['iPP']['DumpPaper']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + meta['Param']['BEV']['HWP']['PaperToDumpPaper']*C_retired
    
    # Transfer to landfill (degradable)
    ip=meta['Core']['iPP']['LandfillWoodDegradable']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + meta['Param']['BEV']['HWP']['PaperToLandfillPaper']*meta['Param']['BEV']['HWP']['ToLandfillPaperDegradableFrac']*C_retired
    
    # Transfer to landfill (non-degradable)
    ip=meta['Core']['iPP']['LandfillWoodNonDegradable']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + meta['Param']['BEV']['HWP']['PaperToLandfillPaper']*(1-meta['Param']['BEV']['HWP']['ToLandfillPaperDegradableFrac'])*C_retired
        
    #--------------------------------------------------------------------------
    # Emissions from fuel combustion
    #--------------------------------------------------------------------------
    
    # Turnover
    ip=meta['Core']['iPP']['Fuel']
    C_retired=meta['Param']['BEV']['HWP']['Fuel_tr']*vo['C_Pro_Pools'][iT,:,ip]
    
    # Remove carbon
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] - C_retired
    
    # Emissions of CO2 from fuel use
    ip=meta['Core']['iPP']['E_CO2']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + meta['Param']['BEV']['HWP']['FuelCombustionFracEmitCO2']*C_retired
    
    # Emissions of CH4 from fuel use
    ip=meta['Core']['iPP']['E_CH4']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + (1-meta['Param']['BEV']['HWP']['FuelCombustionFracEmitCO2'])*C_retired
    
    #--------------------------------------------------------------------------
    # Emissions from firewood combustion
    #--------------------------------------------------------------------------
    
    # Turnover
    ip=meta['Core']['iPP']['Firewood']
    C_retired=meta['Param']['BEV']['HWP']['Firewood_tr']*vo['C_Pro_Pools'][iT,:,ip]
    
    # Remove carbon
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] - C_retired
    
    # Emissions of CO2
    ip=meta['Core']['iPP']['E_CO2']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + C_retired
    
    # Emissions of CH4
    #ip=meta['Core']['iPP']['E_CH4']
    #vo['C_Pro_Pools[iT,:,ip]=vo['C_Pro_Pools[iT,:,ip] + (1-psl['FuelCombustionFracEmitCO2'])*C_retired
       
    #--------------------------------------------------------------------------
    # Emissions from pulp effluent
    #--------------------------------------------------------------------------
         
    # Emissions from pulp effluent (CO2 from aerobic decomposition)
    ip=meta['Core']['iPP']['EffluentPulp']
    c_emitted=meta['Param']['BEV']['HWP']['EffluentPulp_tr']*vo['C_Pro_Pools'][iT,:,ip]
    
    # Remove emitted carbon
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip]-c_emitted
    
    # Add emitted carbon to CO2 emission "pool"
    ip=meta['Core']['iPP']['E_CO2']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip]+c_emitted    
    
    #--------------------------------------------------------------------------
    # Emissions from dump wood
    #--------------------------------------------------------------------------
    
    # Turnover
    ip=meta['Core']['iPP']['DumpWood']
    c_emitted=meta['Param']['BEV']['HWP']['DumpWood_tr']*vo['C_Pro_Pools'][iT,:,ip]
    
    # Removal
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip]-c_emitted
    
    # Add to emissions (CO2 emission from aerobic decomposition)
    ip=meta['Core']['iPP']['E_CO2']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip]+c_emitted
    
    #--------------------------------------------------------------------------
    # Emissions from dump paper
    #--------------------------------------------------------------------------
    
    # Turnover
    ip=meta['Core']['iPP']['DumpPaper']
    c_emitted=meta['Param']['BEV']['HWP']['DumpPaper_tr']*vo['C_Pro_Pools'][iT,:,ip]
    
    # Removal
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip]-c_emitted
    
    # Add to emissions (CO2 emission from aerobic decomposition)
    ip=meta['Core']['iPP']['E_CO2']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip]+c_emitted
    
    #--------------------------------------------------------------------------
    # Emissions from landfill degradable wood
    #--------------------------------------------------------------------------
               
    # Turnover
    ip=meta['Core']['iPP']['LandfillWoodDegradable']
    c_emitted=meta['Param']['BEV']['HWP']['LandfillWoodDegradable_tr']*vo['C_Pro_Pools'][iT,:,ip]
    
    # Removal
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip]-c_emitted
    
    # Add to emissions (50% CO2 emissions during anaerobic decomposition)
    ip=meta['Core']['iPP']['E_CO2']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip]+meta['Param']['BEV']['HWP']['LandfillDegradableFracEmitCO2']*c_emitted
    
    # Add to emissions (50% "potential" CH4 emissions during anaerobic decomposition)
    #E_ch4_pot=(1-meta['Param']['BEV']['HWP']['LandfillDegradableFracEmitCO2)*c_emitted
    
    # Adjustment for proportion of degradable landfills with gas collection systems, 
    # efficiency of system, and methane oxided to CO2 from the landfill cover
    ch4_emitted=c_emitted*((1-meta['Param']['BEV']['HWP']['LandfillMethaneEmit_GasColSysProp'])-meta['Param']['BEV']['HWP']['LandfillMethaneOxidizedToCO2']*(1-meta['Param']['BEV']['HWP']['LandfillMethaneEmit_GasColSysProp'])) + \
        c_emitted*meta['Param']['BEV']['HWP']['LandfillMethaneEmit_GasColSysProp']*((1-meta['Param']['BEV']['HWP']['LandfillMethaneEmit_GasColSysEffic'])-meta['Param']['BEV']['HWP']['LandfillMethaneOxidizedToCO2']*(1-meta['Param']['BEV']['HWP']['LandfillMethaneEmit_GasColSysProp']))
        
    ip=meta['Core']['iPP']['E_CH4']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip]+ch4_emitted
    
    #--------------------------------------------------------------------------
    # Emissions from landfill degradable paper
    #--------------------------------------------------------------------------
    
    # Turnover
    ip=meta['Core']['iPP']['LandfillPaperDegradable']
    c_emitted=meta['Param']['BEV']['HWP']['LandfillWoodDegradable_tr']*vo['C_Pro_Pools'][iT,:,ip]
    
    # Removal
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip]-c_emitted
    
    # Add to emissions (50% CO2 emissions during anaerobic decomposition)
    ip=meta['Core']['iPP']['E_CO2']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip]+meta['Param']['BEV']['HWP']['LandfillDegradableFracEmitCO2']*c_emitted
    
    # Add to emissions (50% "potential" CH4 emissions during anaerobic decomposition)
    #E_ch4_pot=(1-meta['Param']['BEV']['HWP']['LandfillDegradableFracEmitCO2)*c_emitted
    
    # Adjustment for proportion of degradable landfills with gas collection systems, 
    # efficiency of system, and methane oxided to CO2 from the landfill cover
    ch4_emitted=c_emitted*((1-meta['Param']['BEV']['HWP']['LandfillMethaneEmit_GasColSysProp'])-meta['Param']['BEV']['HWP']['LandfillMethaneOxidizedToCO2']*(1-meta['Param']['BEV']['HWP']['LandfillMethaneEmit_GasColSysProp'])) + \
        c_emitted*meta['Param']['BEV']['HWP']['LandfillMethaneEmit_GasColSysProp']*((1-meta['Param']['BEV']['HWP']['LandfillMethaneEmit_GasColSysEffic'])-meta['Param']['BEV']['HWP']['LandfillMethaneOxidizedToCO2']*(1-meta['Param']['BEV']['HWP']['LandfillMethaneEmit_GasColSysProp']))
        
    ip=meta['Core']['iPP']['E_CH4']
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip]+ch4_emitted
    
    return vo
