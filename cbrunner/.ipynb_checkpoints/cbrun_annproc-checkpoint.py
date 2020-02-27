# Import python modules
import os
import numpy as np
import pandas as pd

# Import custom modules
#from cbrun_init import *

# Used to convert species parameter dictionaries into parameter vectors
class Bunch(object):
    def __init__(self, adict):
        self.__dict__.update(adict)

'''============================================================================
BIOMASS DYNAMICS FROM TIPSY
============================================================================'''

def BiomassFromTIPSY(iScn,iT,vi,vo,psl,meta):
    
    # Update stand age
    vo.A[iT,:]=vo.A[iT-1,:]+1
     
    # BIOMASS GROWTH
    
    # Calculate net growth of aboveground biomass pools
    iAge=np.minimum(vo.A[iT,:],vi.GrowthCurveA['VStemMerch'].shape[0])-1
    iAge=iAge.astype(int)
    
    for i in range(meta['N Stand']):
        vo.V_StemMerch[iT,i]=vo.V_StemMerch[iT-1,i]+vi.GrowthCurveA['VStemMerch'][iAge[i],i]
        vo.C_G_Net[iT,i,0]=vi.GrowthCurveA['CStemMerch'][iAge[i],i].astype(float)/meta['Scale Factor GC']
        vo.C_G_Net[iT,i,1]=vi.GrowthCurveA['CStemNonMerch'][iAge[i],i].astype(float)/meta['Scale Factor GC']
        vo.C_G_Net[iT,i,2]=vi.GrowthCurveA['CFoliage'][iAge[i],i].astype(float)/meta['Scale Factor GC']
        vo.C_G_Net[iT,i,3]=vi.GrowthCurveA['CBranch'][iAge[i],i].astype(float)/meta['Scale Factor GC']
        vo.C_G_Net[iT,i,4]=vi.GrowthCurveA['CBark'][iAge[i],i].astype(float)/meta['Scale Factor GC']
        
    #======================================================================
    # *** SPECIAL ORDER ***
    #======================================================================
        
    # This optional component was written to explore the implications of 
    # historical trends in forest growth rate. It should be turned off 
    # unless accompanied with a clear explanation of why it is turned on.
        
    flg=0        
    if flg==1:        
        for i in range(meta['N_Stand']):
            mod=meta['Growth Modifier'][iT]
            vo.V_StemMerch[iT,i]=vo.V_StemMerch[iT-1,i]+mod*vi.GrowthCurveA['VStemMerch'][iAge[i],i]
            vo.C_G_Net[iT,i,0]=mod*vi.GrowthCurveA['CStemMerch'][iAge[i],i]
            vo.C_G_Net[iT,i,1]=mod*vi.GrowthCurveA['CStemNonMerch'][iAge[i],i]
            vo.C_G_Net[iT,i,2]=mod*vi.GrowthCurveA['CFoliage'][iAge[i],i]
            vo.C_G_Net[iT,i,3]=mod*vi.GrowthCurveA['CBranch'][iAge[i],i]
            vo.C_G_Net[iT,i,4]=mod*vi.GrowthCurveA['CBark'][iAge[i],i]
        
    #======================================================================
    # *** SPECIAL ORDER ***
    #======================================================================
                        
    # Add net growth to aboveground biomass pools
    vo.C_Eco_Pools[iT,:,0]=np.maximum(0,vo.C_Eco_Pools[iT-1,:,0]+vo.C_G_Net[iT,:,0])
    vo.C_Eco_Pools[iT,:,1]=np.maximum(0,vo.C_Eco_Pools[iT-1,:,1]+vo.C_G_Net[iT,:,1])
        
    # Foliage - TIPSY foliage appears to be really low, use alternative equations
    Stemwood=np.maximum(0.01,vo.C_Eco_Pools[iT,:,0]+0.2*vo.C_Eco_Pools[iT,:,1])
    vo.C_Eco_Pools[iT,:,2]=np.maximum(0,psl.bAllo_StemToF1*Stemwood**psl.bAllo_StemToF2)
        
    # We have redefined foliage growth, so make sure it is consistent
    vo.C_G_Net[iT,:,2]=vo.C_Eco_Pools[iT,:,2]-vo.C_Eco_Pools[iT-1,:,2]        
    
    vo.C_Eco_Pools[iT,:,3]=np.maximum(0,vo.C_Eco_Pools[iT-1,:,3]+vo.C_G_Net[iT,:,3])
    
    vo.C_Eco_Pools[iT,:,4]=np.maximum(0,vo.C_Eco_Pools[iT-1,:,4]+vo.C_G_Net[iT,:,4])
        
    # Calculate aboveground biomass
    vo.C_BiomassAG[iT,:]=np.sum(vo.C_Eco_Pools[iT,:,0:5],axis=1)
  
    # Total root biomass (Li et al. 2003, Eq. 4)
    BiomassRootTotal=0.222*vo.C_BiomassAG[iT,:]
        
    # Fine root biomass (Li et al. 2003, Eq. 6)
    Pf=0.072+0.354*np.exp(-0.06*(2*BiomassRootTotal))
    vo.C_Eco_Pools[iT,:,6]=Pf*BiomassRootTotal
  
    # Coarse root biomass
    Pc=1-Pf
    vo.C_Eco_Pools[iT,:,5]=np.maximum(0,Pc*BiomassRootTotal)
  
    # Calculate net growth of roots
    vo.C_G_Net[iT,:,5]=vo.C_Eco_Pools[iT,:,5]-vo.C_Eco_Pools[iT-1,:,5]
    vo.C_G_Net[iT,:,6]=vo.C_Eco_Pools[iT,:,6]-vo.C_Eco_Pools[iT-1,:,6]
    
    # BIOMASS TURNOVER
    
    # Calculate biomass loss due to mortality
    vo.C_M[iT,:,0]=psl.bTR_StemMerch[0,:]*vo.C_Eco_Pools[iT,:,0]
    vo.C_M[iT,:,1]=psl.bTR_StemMerch[0,:]*vo.C_Eco_Pools[iT,:,1]
    vo.C_M[iT,:,2]=psl.bTR_StemMerch[0,:]*vo.C_Eco_Pools[iT,:,2]
    vo.C_M[iT,:,3]=psl.bTR_StemMerch[0,:]*vo.C_Eco_Pools[iT,:,3]
    vo.C_M[iT,:,4]=psl.bTR_StemMerch[0,:]*vo.C_Eco_Pools[iT,:,4]
    vo.C_M[iT,:,5]=psl.bTR_StemMerch[0,:]*vo.C_Eco_Pools[iT,:,5]
    vo.C_M[iT,:,6]=psl.bTR_StemMerch[0,:]*vo.C_Eco_Pools[iT,:,6]
    
    # Calculate non-merchantable stemwood biomass turnover due to mortality (treated like litterfall)
    vo.C_LF[iT,:,1]=psl.bTR_StemNonMerch[0,:]*vo.C_Eco_Pools[iT,:,1]
    
    # Calculate foliage biomass turnover due to litterfall
    vo.C_LF[iT,:,2]=psl.bTR_Foliage[0,:]*vo.C_Eco_Pools[iT,:,2]
    
    # Calculate branch biomass turnover due to litterfall
    vo.C_LF[iT,:,3]=psl.bTR_Branch[0,:]*vo.C_Eco_Pools[iT,:,3]
    
    # Calculate bark biomass turnover due to litterfall
    vo.C_LF[iT,:,4]=psl.bTR_Bark[0,:]*vo.C_Eco_Pools[iT,:,4]
    
    # Calculate coarse root biomass turnover due to litterfall
    vo.C_LF[iT,:,5]=psl.bTR_RootCoarse[0,:]*vo.C_Eco_Pools[iT,:,5]
    
    # Calculate fine root biomass turnover due to litterfall
    vo.C_LF[iT,:,6]=psl.bTR_RootFine[0,:]*vo.C_Eco_Pools[iT,:,6]
            
    # Update gross growth
    vo.C_G_Gross[iT,:,:]=vo.C_G_Net[iT,:,:]+vo.C_M[iT,:,:]
    
    # Update NPP
    # NPP = G
    vo.C_NPP[iT,:,:]=vo.C_G_Net[iT,:,:]+vo.C_M[iT,:,:]+vo.C_LF[iT,:,:]
    
    return vo

'''============================================================================
ANNUAL TREE BIOMASS DYNAMICS
============================================================================'''

def BiomassFromSawtooth(iScn,iS,vi,vo,ptl,meta):
    
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
    tl_N_M_Sim_Dir=np.zeros((N_t,N_tree))
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
    n1=np.ceil(vi.inv.Srs1_Pct[0,iS]/100*N_tree).astype(int)
    tl_ID_Srs[:,rp[0:n1]]=vi.inv.Srs1_ID[0,iS]
    # Species 2
    if vi.inv.Srs2_Pct[0,iS]>0:
        n2=np.ceil(vi.inv.Srs2_Pct[0,iS]/100*N_tree).astype(int)
        tl_ID_Srs[:,rp[n1:n2]]=vi.inv.Srs2_ID[0,iS]
    # Species 3
    if vi.inv.Srs3_Pct[0,iS]>0:
        tl_ID_Srs[:,rp[n2:]]=vi.inv.Srs3_ID[0,iS]
    
    #--------------------------------------------------------------------------
    # Initialize parameter vectors
    #--------------------------------------------------------------------------
    
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
    
    u=np.unique(tl_ID_Srs)
    for i in range(0,len(u)):
        
        ind=np.where(tl_ID_Srs[0,:]==u[i])[0]
        
        # Allometry
        for k in bA.__dict__.keys():
            if k!='SRS_CD':
                tmp=getattr(bA,k)
                tmp[0,ind]=ptl.Allom[u[i]-1][k]
                setattr(bA,k,tmp)
        
        # Recruitment
        for k in bR.__dict__.keys():
            if k!='SRS_CD':
                tmp=getattr(bR,k)
                tmp[0,ind]=ptl.R_Coef[iEqR][u[i]-1][k]
                setattr(bR,k,tmp)
        
        # Mortality
        for k in bM.__dict__.keys():
            if k!='SRS_CD':
                tmp=getattr(bM,k)
                tmp[0,ind]=ptl.M_Coef[iEqM][u[i]-1][k]
                setattr(bM,k,tmp)
        
        # Growth
        for k in bG.__dict__.keys():
            if k!='SRS_CD':
                tmp=getattr(bG,k)
                tmp[0,ind]=ptl.G_Coef[iEqG][u[i]-1][k]
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
        
        if (vi.tv[iT]>meta['Spinup Disturbance Return Inverval']+1) & (vi.tv[iT]<=meta['Year Spinup Ends']+1):
            
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
            tl_N_M_Sim_Dir[iT,:]=tl_N_M_Sim_Dir[iT0,:]
            tl_N_M_Sim_Fir[iT,:]=tl_N_M_Sim_Fir[iT0,:]
            tl_N_M_Sim_Ins[iT,:]=tl_N_M_Sim_Ins[iT0,:]
            tl_N_M_Sim_Pat[iT,:]=tl_N_M_Sim_Pat[iT0,:]
            tl_N_M_Sim_Har[iT,:]=tl_N_M_Sim_Har[iT0,:]
            tl_N_M_Sim_Win[iT,:]=tl_N_M_Sim_Win[iT0,:]
            tl_Csw[iT,:]=tl_Csw[iT0,:]
            tl_Csw_Larger[iT,:]=tl_Csw_Larger[iT0,:]
            tl_Csw_G[iT,:]=tl_Csw_G[iT0,:]
         
            vo.A[iT,iS]=vo.A[iT0,iS]
            
            continue
        
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
        vo.A[iT,iS]=vo.A[iT-1,iS]+1
    
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
            
            Pm_Sim_Dir=np.exp(lgit)/(1+np.exp(lgit))
            Pm_Sim_Ins=np.zeros((1,N_tree))
            Pm_Sim_Pat=np.zeros((1,N_tree))
        
        #---------------------------------------------------------------------- 
        # Remove biomass of trees that died
        #---------------------------------------------------------------------- 
        
        # Remove biomass of trees that died directly from competition and 
        # environmental conditions
        iKill=np.where(Pm_Sim_Dir[0,iLive]>=rLive)
        tl_A[iT:,iLive[iKill]]=np.nan
        tl_D[iT:,iLive[iKill]]=np.nan
        tl_H[iT:,iLive[iKill]]=np.nan
        tl_Csw[iT:,iLive[iKill]]=np.nan
        tl_N_M_Sim_Dir[iT,iLive[iKill]]=1   
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
        
        iDist=np.where((vi.dh[iS].Year==vi.tv[iT]) & (np.isnan(vi.dh[iS].Severity)==False))[0]
        if iDist.size!=0:            
            
            iLive=np.where(np.isnan(tl_Csw[iT,:])==False)[0]
            
            nKill=int(vi.dh[iS].Severity[iDist]/100*iLive.size)
                        
            if vi.dh[iS].ID[iDist]==meta['LUT Dist']['Wildfire']:
                # Wildfire
                iKill=np.arange(0,nKill)
                tl_N_M_Inv_Fir[iT,iLive[iKill]]=1
                vo.A[iT,iS]=0
            elif vi.dh[iS].ID[iDist]==meta['LUT Dist']['Clearcut']:
                # Clearcut
                iKill=np.arange(0,nKill)
                tl_N_M_Inv_Har[iT,iLive[iKill]]=1
                vo.A[iT,iS]=0
            elif vi.dh[iS].ID[iDist]==meta['LUT Dist']['Clearcut and Slash Burn']:
                # Clearcut and slash burn
                iKill=np.arange(0,nKill)
                tl_N_M_Inv_Har[iT,iLive[iKill]]=1
                vo.A[iT,iS]=0
            elif vi.dh[iS].ID[iDist]==meta['LUT Dist']['Clearcut and Salvage Logging']:
                # Clearcut and salvage logging
                iKill=np.arange(0,nKill)
                tl_N_M_Inv_Har[iT,iLive[iKill]]=1
                vo.A[iT,iS]=0    
            elif vi.dh[iS].ID[iDist]==meta['LUT Dist']['Insects']:
                # Insects
                iKill=np.arange(0,nKill)
                tl_N_M_Inv_Ins[iT,iLive[iKill]]=1
            #elif vi.dh[iS].ID[iDist]==meta['LUT Dist']['Pre-commercial Thinning']:
                # Pre-commercial Thinning
            #    iKill=np.arange(0,nKill)
            #    tl_N_M_Inv_Ins[iT,iLive[iKill]]=1    
            else:
                pass
            
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
    N_M=tl_N_M_Sim_Dir+tl_N_M_Sim_Fir+tl_N_M_Sim_Ins+tl_N_M_Sim_Pat+ \
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
  
    tl_Ctot_M_Sim_Dir=tl_N_M_Sim_Dir*C_M
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
    vo.TreeMean_A[:,iS]=np.nanmean(tl_A,axis=1) 
    
    # Mean tree diamter at breast height
    vo.TreeMean_D[:,iS]=np.nanmean(tl_D,axis=1)
    
    # Mean tree height (m)
    vo.TreeMean_H[:,iS]=np.nanmean(tl_H,axis=1)
    
    # Mean tree stemwood biomass (kg C)
    vo.TreeMean_Csw[:,iS]=np.nanmean(tl_Csw,axis=1) 
    
    # Mean stemwood biomass growth (kg C yr-1)
    vo.TreeMean_Csw_G[:,iS]=np.nanmean(tl_Csw_G,axis=1)
    
    #--------------------------------------------------------------------------
    # Populate stand demographics
    #--------------------------------------------------------------------------
    
    # Stand density (stems ha-1)    
    vo.N[:,iS]=np.sum(~np.isnan(tl_Csw),axis=1)
    
    # Stand density at beginning of time step
    N_t0=np.append(vo.N[0,iS],vo.N[0:-1,iS])
    
    # Demographic recruitment rate (% yr-1)
    vo.N_R[:,iS]=np.minimum(100,np.sum(tl_N_R,axis=1)/N_t0*100)
  
    # Demographic mortality rate (% yr-1)
          
    # Need to update N0 or else stand-replacing disturbances will have lower
    # relative mortality rate
    # *** The other solution is to move the prescribed disturbances above the
    # background mortality, but perhaps GCBM will do this anyways ***
    #N0=N0-Sum_N_M_Sim_Dir-Sum_N_M_Sim_Ins-Sum_N_M_Sim_Pat;
    
    all=tl_N_M_Sim_Dir+tl_N_M_Sim_Fir+tl_N_M_Sim_Ins+tl_N_M_Sim_Pat+ \
        tl_N_M_Sim_Win+tl_N_M_Sim_Har+tl_N_M_Inv_Fir+tl_N_M_Inv_Ins+ \
        tl_N_M_Inv_Pat+tl_N_M_Inv_Win+tl_N_M_Inv_Har
    
    vo.N_M[:,iS]=np.minimum(100,np.sum(all,axis=1)/N_t0*100)
    
    vo.N_M_Sim_Dir[:,iS]=np.maximum(0,np.sum(tl_N_M_Sim_Dir,axis=1)/N_t0*100)
    vo.N_M_Sim_Fir[:,iS]=np.maximum(0,np.sum(tl_N_M_Sim_Fir,axis=1)/N_t0*100)
    vo.N_M_Sim_Ins[:,iS]=np.maximum(0,np.sum(tl_N_M_Sim_Ins,axis=1)/N_t0*100)
    vo.N_M_Sim_Pat[:,iS]=np.maximum(0,np.sum(tl_N_M_Sim_Pat,axis=1)/N_t0*100)
    vo.N_M_Sim_Win[:,iS]=np.maximum(0,np.sum(tl_N_M_Sim_Win,axis=1)/N_t0*100)
    vo.N_M_Sim_Har[:,iS]=np.maximum(0,np.sum(tl_N_M_Sim_Har,axis=1)/N_t0*100)
  
    vo.N_M_Inv_Fir[:,iS]=np.maximum(0,np.sum(tl_N_M_Inv_Fir,axis=1)/N_t0*100)
    vo.N_M_Inv_Ins[:,iS]=np.maximum(0,np.sum(tl_N_M_Inv_Ins,axis=1)/N_t0*100)
    vo.N_M_Inv_Pat[:,iS]=np.maximum(0,np.sum(tl_N_M_Inv_Pat,axis=1)/N_t0*100)
    vo.N_M_Inv_Win[:,iS]=np.maximum(0,np.sum(tl_N_M_Inv_Win,axis=1)/N_t0*100)
    vo.N_M_Inv_Har[:,iS]=np.maximum(0,np.sum(tl_N_M_Inv_Har,axis=1)/N_t0*100)
    
    #--------------------------------------------------------------------------
    # Populate stand gross growth of biomass (Mg C ha-1 yr-1)
    #--------------------------------------------------------------------------
    
    vo.C_G_Gross[:,iS,0]=np.nansum(np.maximum(0,diff(tl_Csw_Merch)),axis=1)/cm
    vo.C_G_Gross[:,iS,1]=np.nansum(np.maximum(0,diff(tl_Csw_NonMerch)),axis=1)/cm
    vo.C_G_Gross[:,iS,2]=np.nansum(np.maximum(0,diff(tl_Cf)),axis=1)/cm
    vo.C_G_Gross[:,iS,3]=np.nansum(np.maximum(0,diff(tl_Cbr)),axis=1)/cm
    vo.C_G_Gross[:,iS,4]=np.nansum(np.maximum(0,diff(tl_Cbk)),axis=1)/cm
    vo.C_G_Gross[:,iS,5]=np.nansum(np.maximum(0,diff(tl_Crc)),axis=1)/cm
    vo.C_G_Gross[:,iS,6]=np.nansum(np.maximum(0,diff(tl_Crf)),axis=1)/cm
    
    #--------------------------------------------------------------------------
    # Populate stand biomass turnover (Mg C ha-1 yr-1)
    #--------------------------------------------------------------------------
    
    vo.C_LF[:,iS,2]=np.nansum(tl_Cf_LF,axis=1)/cm
    vo.C_LF[:,iS,3]=np.nansum(tl_Cbr_LF,axis=1)/cm
    vo.C_LF[:,iS,4]=np.nansum(tl_Cbk_LF,axis=1)/cm   
    vo.C_LF[:,iS,5]=np.nansum(tl_Crc_LF,axis=1)/cm
    vo.C_LF[:,iS,6]=np.nansum(tl_Crf_LF,axis=1)/cm
    
    #--------------------------------------------------------------------------
    # Populate stand biomass loss due to mortality (Mg C ha-1 yr-1)
    #--------------------------------------------------------------------------
    
    vo.C_M[:,iS,0]=np.nansum(tl_Csw_Merch_M,axis=1)/cm
    vo.C_M[:,iS,1]=np.nansum(tl_Csw_NonMerch_M,axis=1)/cm
    vo.C_M[:,iS,2]=np.nansum(tl_Cf_M,axis=1)/cm
    vo.C_M[:,iS,3]=np.nansum(tl_Cbr_M,axis=1)/cm
    vo.C_M[:,iS,4]=np.nansum(tl_Cbk_M,axis=1)/cm
    vo.C_M[:,iS,5]=np.nansum(tl_Crc_M,axis=1)/cm
    vo.C_M[:,iS,6]=np.nansum(tl_Crf_M,axis=1)/cm
    
    vo.C_M_Sim_Dir[:,iS]=np.nansum(tl_Ctot_M_Sim_Dir,axis=1)/cm
    vo.C_M_Sim_Fir[:,iS]=np.nansum(tl_Ctot_M_Sim_Fir,axis=1)/cm
    vo.C_M_Sim_Ins[:,iS]=np.nansum(tl_Ctot_M_Sim_Ins,axis=1)/cm
    vo.C_M_Sim_Pat[:,iS]=np.nansum(tl_Ctot_M_Sim_Pat,axis=1)/cm
    vo.C_M_Sim_Win[:,iS]=np.nansum(tl_Ctot_M_Sim_Win,axis=1)/cm
    vo.C_M_Sim_Har[:,iS]=np.nansum(tl_Ctot_M_Sim_Har,axis=1)/cm
  
    vo.C_M_Inv_Fir[:,iS]=np.nansum(tl_Ctot_M_Inv_Fir,axis=1)/cm
    vo.C_M_Inv_Ins[:,iS]=np.nansum(tl_Ctot_M_Inv_Ins,axis=1)/cm
    vo.C_M_Inv_Pat[:,iS]=np.nansum(tl_Ctot_M_Inv_Pat,axis=1)/cm
    vo.C_M_Inv_Win[:,iS]=np.nansum(tl_Ctot_M_Inv_Win,axis=1)/cm
    vo.C_M_Inv_Har[:,iS]=np.nansum(tl_Ctot_M_Inv_Har,axis=1)/cm
            
    #--------------------------------------------------------------------------
    # Populate stand NPP by tissue (Mg C ha-1 yr-1)
    #--------------------------------------------------------------------------
    
    vo.C_NPP[:,iS,0:7]=vo.C_G_Gross[:,iS,0:7]+vo.C_LF[:,iS,0:7]
    
    #--------------------------------------------------------------------------
    # Populate stand net growth by tissue (Mg C ha-1 yr-1)
    #--------------------------------------------------------------------------
    
    vo.C_G_Net[:,iS,0:7]=vo.C_G_Gross[:,iS,0:7]-vo.C_M[:,iS,0:7]
    
    #--------------------------------------------------------------------------
    # Populate ecosystem carbon pools by tissue (Mg C ha-1 yr-1)
    #--------------------------------------------------------------------------
    
    vo.C_Eco_Pools[:,iS,0]=np.maximum(0,np.nansum(tl_Csw_Merch,axis=1))/cm
    vo.C_Eco_Pools[:,iS,1]=np.maximum(0,np.nansum(tl_Csw_NonMerch,axis=1))/cm
    vo.C_Eco_Pools[:,iS,2]=np.maximum(0,np.nansum(tl_Cf,axis=1))/cm
    vo.C_Eco_Pools[:,iS,3]=np.maximum(0,np.nansum(tl_Cbr,axis=1))/cm
    vo.C_Eco_Pools[:,iS,4]=np.maximum(0,np.nansum(tl_Cbk,axis=1))/cm
    vo.C_Eco_Pools[:,iS,5]=np.maximum(0,np.nansum(tl_Crc,axis=1))/cm
    vo.C_Eco_Pools[:,iS,6]=np.maximum(0,np.nansum(tl_Crf,axis=1))/cm
        
    return vo


flg=0
if flg==1:
    
    fs=8
    plt.close('all')
    fig,ax=plt.subplots(2,2,figsize=(14,10))
    plt.rcParams.update({'font.size':fs})
    plt.rc('xtick',labelsize=fs); plt.rc('ytick',labelsize=fs)
    ax[0,0].plot(vo.N[:,iS],'-',linewidth=0.75)
    #ax[0,0].set(ylabel=r'Tree age, years',xlabel=r'Time, calendar year');
    ax[0,1].plot(np.nanmean(tl_Csw_G,axis=1),'-',linewidth=0.75)
    ax[1,0].plot(np.nanmean(tl_Csw,axis=1),'-',linewidth=0.75)
    ax[1,1].plot(np.nansum(tl_Csw_Merch_M,axis=1),'-',linewidth=0.75)
    plt.tight_layout();


    fs=8
    plt.close('all')
    fig,ax=plt.subplots(3,2,figsize=(14,10))
    plt.rcParams.update({'font.size':fs})
    plt.rc('xtick',labelsize=fs); plt.rc('ytick',labelsize=fs)
    ax[0,0].plot(vo.N[:,iS],'-',linewidth=0.75)
    ax[0,1].plot(vo.C_Eco_Pools[:,iS,0],'-',linewidth=0.75)    
    ax[1,0].plot(vo.C_G_Net[:,iS,0],'-',linewidth=0.75)
    ax[1,1].plot(vo.C_M[:,iS,0],'-',linewidth=0.75)
    ax[2,0].plot(vo.C_LF[:,iS,2],'-',linewidth=0.75)
    plt.tight_layout();


'''============================================================================
STAND DEAD ORGANIC MATTER DYNAMICS
============================================================================'''    

def DeadOrganicMatterDynamics(iT,vi,vo,psl):
     
    # List of pools
    # Number 8          9         10        11        12         13           14      15       16      17           18     19     20    21   22    23          24                   25
    #Index:  7          8         9         10        11         12           13      14       15      16           17     18     19    20   21    22          23                   24
    #Name:  'LitterVF','LitterF','LitterM','LitterS','SnagStem','SnagBranch','SoilVF','SoilF','SoilS','BlackCarbon','Peat','CO2','CH4','CO','N2O','RemovedMerch','RemovedNonMerch','RemovedSnagStem'

    #--------------------------------------------------------------------------
    # FLUX FROM BIOMASS POOLS TO DEAD ORGANIC MATTER POOLS
    #--------------------------------------------------------------------------
    
    # Transfer biomass turnover to very fast litter pool
    vo.C_Eco_Pools[iT,:,7]=vo.C_Eco_Pools[iT-1,:,7]+ \
        psl.bIPF_FoliageLitToLitterVF*vo.C_LF[iT,:,2]+ \
        psl.bIPF_RootFineLitToLitterVF*vo.C_LF[iT,:,6]+ \
        psl.bIPF_FoliageMorToLitterVF*vo.C_M[iT,:,2]
    
    # Transfer biomass turnover to fast litter pool
    vo.C_Eco_Pools[iT,:,8]=vo.C_Eco_Pools[iT-1,:,8]+ \
        psl.bIPF_BranchLitToLitterF*vo.C_LF[iT,:,3]+ \
        psl.bIPF_BarkLitToLitterF*vo.C_LF[iT,:,4]+ \
        psl.bIPF_RootCoarseLitToLitterF*vo.C_LF[iT,:,5]+ \
        psl.bIPF_BranchMorToLitterF*vo.C_M[iT,:,3]+ \
        psl.bIPF_StemNonMerchLitToLitterF*vo.C_LF[iT,:,1]
            
    # Transfer biomass turnover to medium litter pool
    vo.C_Eco_Pools[iT,:,9]=vo.C_Eco_Pools[iT-1,:,9]+ \
        psl.bIPF_BarkLitToLitterM*vo.C_LF[iT,:,4]+ \
        psl.bIPF_BarkMorToLitterM*vo.C_M[iT,:,4]

    # Transfer biomass turnover to slow litter pool
    vo.C_Eco_Pools[iT,:,10]=vo.C_Eco_Pools[iT-1,:,10]+ \
        psl.bIPF_BarkLitToLitterS*vo.C_LF[iT,:,4]+ \
        psl.bIPF_BarkMorToLitterS*vo.C_M[iT,:,4]
        
    # Transfer biomass turnover to snag stemwood pool
    vo.C_Eco_Pools[iT,:,11]=vo.C_Eco_Pools[iT-1,:,11]+ \
        psl.bIPF_StemMerchMorToSnagStem*vo.C_M[iT,:,0]+ \
        psl.bIPF_StemNonMerchMorToSnagStem*vo.C_M[iT,:,1]
    
    # Transfer biomass turnvoer to snag branch pool
    vo.C_Eco_Pools[iT,:,12]=vo.C_Eco_Pools[iT-1,:,12]+ \
        psl.bIPF_BranchLitToSnagBranch*vo.C_LF[iT,:,3]+ \
        psl.bIPF_BranchMorToSnagBranch*vo.C_M[iT,:,3]+ \
        psl.bIPF_BarkMorToSnagBranch*vo.C_M[iT,:,4]+ \
        psl.bIPF_StemNonMerchLitToSnagBranch*vo.C_LF[iT,:,1]

    # Transfer biomass turnover to very fast soil pool
    vo.C_Eco_Pools[iT,:,13]=vo.C_Eco_Pools[iT-1,:,13]+ \
        psl.bIPF_RootFineLitToSoilVF*vo.C_LF[iT,:,6] +\
        psl.bIPF_RootCoarseMorToSoilVF*vo.C_M[iT,:,5]+ \
        psl.bIPF_RootFineMorToSoilVF*vo.C_M[iT,:,6]

    # Transfer biomass turnover to fast soil pool
    vo.C_Eco_Pools[iT,:,14]=vo.C_Eco_Pools[iT-1,:,14]+ \
        psl.bIPF_RootCoarseLitToSoilF*vo.C_LF[iT,:,5]+ \
        psl.bIPF_RootCoarseMorToSoilF*vo.C_M[iT,:,5]
        
    # Transfer biomass turnover to slow soil pool
    vo.C_Eco_Pools[iT,:,15]=vo.C_Eco_Pools[iT-1,:,15]    
    
    #--------------------------------------------------------------------------
    # DECOMPOSITION
    #--------------------------------------------------------------------------
    
    # Prepare air temperature for respiration calculation
    Tref=10
    fT=(vi.inv.MAT-Tref)/10    
       
    # Respiration rate - Note that these terms do not equal the
    # atmosphere-bound efflux from heterotrophic respiration as a fraction is
    # emitted to the atmosphere and the remaining fraction is reorganized
    # within the ecosystem.
    
    R_LitterVF=psl.bDec_LitterVF_Rten*vo.C_Eco_Pools[iT,:,7]*psl.bDec_LitterVF_Qten**fT
    R_LitterF=psl.bDec_LitterF_Rten*vo.C_Eco_Pools[iT,:,8]*psl.bDec_LitterF_Qten**fT
    R_LitterM=psl.bDec_LitterM_Rten*vo.C_Eco_Pools[iT,:,9]*psl.bDec_LitterM_Qten**fT
    R_LitterS=psl.bDec_LitterS_Rten*vo.C_Eco_Pools[iT,:,10]*psl.bDec_LitterS_Qten**fT
    R_SnagStem=psl.bDec_SnagStem_Rten*vo.C_Eco_Pools[iT,:,11]*psl.bDec_SnagStem_Qten**fT
    R_SnagBranch=psl.bDec_SnagBranch_Rten*vo.C_Eco_Pools[iT,:,12]*psl.bDec_SnagBranch_Qten**fT    
    R_SoilVF=psl.bDec_SoilVF_Rten*vo.C_Eco_Pools[iT,:,13]*psl.bDec_SoilVF_Qten**fT
    R_SoilF=psl.bDec_SoilF_Rten*vo.C_Eco_Pools[iT,:,14]*psl.bDec_SoilF_Qten**fT
    R_SoilS=psl.bDec_SoilS_Rten*vo.C_Eco_Pools[iT,:,15]*psl.bDec_SoilS_Qten**fT
    
    # These need to be flattened *** Is there a better way? ***
    R_LitterVF=R_LitterVF.flatten()
    R_LitterF=R_LitterF.flatten()
    R_LitterM=R_LitterM.flatten()
    R_LitterS=R_LitterS.flatten()
    R_SnagStem=R_SnagStem.flatten()
    R_SnagBranch=R_SnagBranch.flatten()
    R_SoilVF=R_SoilVF.flatten()
    R_SoilF=R_SoilF.flatten()
    R_SoilS=R_SoilS.flatten()
    
    # Remove respired carbon from source DOM pools
    vo.C_Eco_Pools[iT,:,7]=vo.C_Eco_Pools[iT,:,7]-R_LitterVF
    vo.C_Eco_Pools[iT,:,8]=vo.C_Eco_Pools[iT,:,8]-R_LitterF
    vo.C_Eco_Pools[iT,:,9]=vo.C_Eco_Pools[iT,:,9]-R_LitterM
    vo.C_Eco_Pools[iT,:,10]=vo.C_Eco_Pools[iT,:,10]-R_LitterS
    vo.C_Eco_Pools[iT,:,11]=vo.C_Eco_Pools[iT,:,11]-R_SnagStem
    vo.C_Eco_Pools[iT,:,12]=vo.C_Eco_Pools[iT,:,12]-R_SnagBranch
    vo.C_Eco_Pools[iT,:,13]=vo.C_Eco_Pools[iT,:,13]-R_SoilVF
    vo.C_Eco_Pools[iT,:,14]=vo.C_Eco_Pools[iT,:,14]-R_SoilF
    vo.C_Eco_Pools[iT,:,15]=vo.C_Eco_Pools[iT,:,15]-R_SoilS     
                
    # Re-define decayed fast litter
    vo.C_Eco_Pools[iT,:,8]=vo.C_Eco_Pools[iT,:,8]+psl.bIPF_SnagStemToLitterF*R_SnagStem
    
    # Re-define decayed medium litter
    vo.C_Eco_Pools[iT,:,9]=vo.C_Eco_Pools[iT,:,9]+psl.bIPF_SnagBranchToLitterM*R_SnagBranch
    
    # Re-define decayed slow litter
    vo.C_Eco_Pools[iT,:,10]=vo.C_Eco_Pools[iT,:,10]+psl.bIPF_LitterVFToLitterS*R_LitterVF
    vo.C_Eco_Pools[iT,:,10]=vo.C_Eco_Pools[iT,:,10]+psl.bIPF_LitterFToLitterS*R_LitterF
    vo.C_Eco_Pools[iT,:,10]=vo.C_Eco_Pools[iT,:,10]+psl.bIPF_LitterMToLitterS*R_LitterM    
    
    # Re-define decayed slow soil
    vo.C_Eco_Pools[iT,:,15]=vo.C_Eco_Pools[iT,:,15]+psl.bIPF_SoilVFToSoilS*R_SoilVF
    vo.C_Eco_Pools[iT,:,15]=vo.C_Eco_Pools[iT,:,15]+psl.bIPF_SoilFToSoilS*R_SoilF
        
    # Heterotrophic respiration 
    vo.C_RH[iT,:,7]=psl.bIPF_LitterVFToCO2*R_LitterVF
    vo.C_RH[iT,:,8]=psl.bIPF_LitterFToCO2*R_LitterF
    vo.C_RH[iT,:,9]=psl.bIPF_LitterMToCO2*R_LitterM
    vo.C_RH[iT,:,10]=psl.bIPF_LitterSToCO2*R_LitterS
    vo.C_RH[iT,:,11]=psl.bIPF_SnagStemToCO2*R_SnagStem
    vo.C_RH[iT,:,12]=psl.bIPF_SnagBranchToCO2*R_SnagBranch
    vo.C_RH[iT,:,13]=psl.bIPF_SoilVFToCO2*R_SoilVF
    vo.C_RH[iT,:,14]=psl.bIPF_SoilFToCO2*R_SoilF
    vo.C_RH[iT,:,15]=psl.bIPF_SoilSToCO2*R_SoilS    
    
    #--------------------------------------------------------------------------
    # PHYSICAL TRANSFER
    #--------------------------------------------------------------------------
    
    # Physical transfer rate
    PT_LitterS=psl.bDec_LitterS_PhysTransRate*vo.C_Eco_Pools[iT,:,10]
    PT_SnagStem=psl.bDec_SnagStem_PhysTransRate*vo.C_Eco_Pools[iT,:,11]
    PT_SnagBranch=psl.bDec_SnagBranch_PhysTransRate*vo.C_Eco_Pools[iT,:,12] 
    
    # Remove carbon that is physically transferred
    vo.C_Eco_Pools[iT,:,10]=vo.C_Eco_Pools[iT,:,10]-PT_LitterS
    vo.C_Eco_Pools[iT,:,11]=vo.C_Eco_Pools[iT,:,11]-PT_SnagStem
    vo.C_Eco_Pools[iT,:,12]=vo.C_Eco_Pools[iT,:,12]-PT_SnagBranch
    
    # Re-define pools after transfer
    vo.C_Eco_Pools[iT,:,15]=vo.C_Eco_Pools[iT,:,15]+PT_LitterS
    vo.C_Eco_Pools[iT,:,9]=vo.C_Eco_Pools[iT,:,9]+PT_SnagStem
    vo.C_Eco_Pools[iT,:,8]=vo.C_Eco_Pools[iT,:,8]+PT_SnagBranch    
    
    return vo

'''============================================================================
STAND DISTURBANCES
============================================================================'''

def Disturbances(iT,vi,vo,psl,meta):
        
    for iS in range(meta['N Stand']):
        
        # Only continue if there are disturbances
        iDist=np.where(vi.dh[iS].Year==vi.tv[iT])[0]
        
        if len(iDist)==0:
            continue
        
        for i in range(iDist.size):
            
            # Index to default parameters for this disturbance type
            iDefPar=np.where(psl.bDist_ID==vi.dh[iS].ID_Type[iDist[i]])[0]

            if vi.dh[iS].Biomass_Affected_Pct[iDist[i]]<=0:
                
                # No custom harvesting disturbance parameters have been supplied               
                
                # Use severity from disturbance history
                Severity=vi.dh[iS].Severity[iDist[i]].astype(float)/100
                Biomass_Affected_Pct=Severity
                Snags_Affected_Pct=Severity
                
                # Use default disturbance fluxes
                Biomass_Merch_Removed_Pct=psl.bDist_Biomass_Merch_Removed[iDefPar]
                Biomass_Merch_Burned_Pct=psl.bDist_Biomass_Merch_Burned[iDefPar]
                Biomass_Merch_LeftOnSite_Pct=psl.bDist_Biomass_Merch_LeftOnSite[iDefPar]
                Biomass_Merch_ToSnag_Pct=psl.bDist_Biomass_Merch_ToSnag[iDefPar]
                Biomass_NonMerch_Removed_Pct=psl.bDist_Biomass_NonMerch_Removed[iDefPar]
                Biomass_NonMerch_Burned_Pct=psl.bDist_Biomass_NonMerch_Burned[iDefPar]
                Biomass_NonMerch_LeftOnSite_Pct=psl.bDist_Biomass_NonMerch_LeftOnSite[iDefPar]
                Biomass_NonMerch_ToSnag_Pct=psl.bDist_Biomass_NonMerch_ToSnag[iDefPar]
                Snags_Removed_Pct=psl.bDist_Snags_Removed[iDefPar]
                Snags_Burned_Pct=psl.bDist_Snags_Burned[iDefPar]
                Snags_LeftOnSite_Pct=psl.bDist_Snags_LeftOnSite[iDefPar]
                
            else:
                
                # Use custom harvesting disturbance parameter table
                
                # These come from a custom user input spreadsheet, but they have been added to DH structure in pre-processing
                Biomass_Affected_Pct=vi.dh[iS].Biomass_Affected_Pct[iDist[i]].astype(float)/100
                Snags_Affected_Pct=vi.dh[iS].Snags_Affected_Pct[iDist[i]].astype(float)/100
                
                # Proportions removed, burned, and left on site                        
                Biomass_Merch_Removed_Pct=vi.dh[iS].Biomass_Merch_Removed_Pct[iDist[i]].astype(float)/100
                Biomass_Merch_Burned_Pct=vi.dh[iS].Biomass_Merch_Burned_Pct[iDist[i]].astype(float)/100          
                Biomass_Merch_LeftOnSite_Pct=vi.dh[iS].Biomass_Merch_LeftOnSite_Pct[iDist[i]].astype(float)/100
                Biomass_Merch_ToSnag_Pct=0
            
                Biomass_NonMerch_Removed_Pct=vi.dh[iS].Biomass_NonMerch_Removed_Pct[iDist[i]].astype(float)/100
                Biomass_NonMerch_Burned_Pct=vi.dh[iS].Biomass_NonMerch_Burned_Pct[iDist[i]].astype(float)/100          
                Biomass_NonMerch_LeftOnSite_Pct=vi.dh[iS].Biomass_NonMerch_LeftOnSite_Pct[iDist[i]].astype(float)/100
                Biomass_NonMerch_ToSnag_Pct=0
                
                Snags_Removed_Pct=vi.dh[iS].Snags_Removed_Pct[iDist[i]].astype(float)/100
                Snags_Burned_Pct=vi.dh[iS].Snags_Burned_Pct[iDist[i]].astype(float)/100
                Snags_LeftOnSite_Pct=vi.dh[iS].Snags_LeftOnSite_Pct[iDist[i]].astype(float)/100
                        
            #------------------------------------------------------------------
            # Total amount of each pool that is affected
            #------------------------------------------------------------------
            
            # Define biomass that is affected
            Affected_StemMerch=Biomass_Affected_Pct*vo.C_Eco_Pools[iT,iS,0]
            Affected_StemNonMerch=Biomass_Affected_Pct*vo.C_Eco_Pools[iT,iS,1]
            Affected_Foliage=Biomass_Affected_Pct*vo.C_Eco_Pools[iT,iS,2]
            Affected_Branch=Biomass_Affected_Pct*vo.C_Eco_Pools[iT,iS,3]
            Affected_Bark=Biomass_Affected_Pct*vo.C_Eco_Pools[iT,iS,4]
            Affected_RootCoarse=Biomass_Affected_Pct*vo.C_Eco_Pools[iT,iS,5]
            Affected_RootFine=Biomass_Affected_Pct*vo.C_Eco_Pools[iT,iS,6]
                                    
            # Assume 80% of bark is on merchantable-sized trees, 20% on non-merchantable sized trees
            Affected_Bark_Merch=0.85*Affected_Bark
            Affected_Bark_NonMerch=0.15*Affected_Bark
                
            # Sum up total affected non-merchantable biomass
            Affected_TotNonMerch=Affected_StemNonMerch+Affected_Branch+Affected_Bark_NonMerch
            
            # Define fraction of snags that are affected            
            Affected_SnagStem=Snags_Affected_Pct*vo.C_Eco_Pools[iT,iS,11]
            Affected_SnagBranch=Snags_Affected_Pct*vo.C_Eco_Pools[iT,iS,12]
            
            #------------------------------------------------------------------
            # Remove total amount of each pool affected
            #------------------------------------------------------------------
            
            # Remove affected biomass from biomass pools
            vo.C_Eco_Pools[iT,iS,0]=vo.C_Eco_Pools[iT,iS,0]-Affected_StemMerch
            vo.C_Eco_Pools[iT,iS,1]=vo.C_Eco_Pools[iT,iS,1]-Affected_StemNonMerch
            vo.C_Eco_Pools[iT,iS,2]=vo.C_Eco_Pools[iT,iS,2]-Affected_Foliage
            vo.C_Eco_Pools[iT,iS,3]=vo.C_Eco_Pools[iT,iS,3]-Affected_Branch
            vo.C_Eco_Pools[iT,iS,4]=vo.C_Eco_Pools[iT,iS,4]-Affected_Bark
            vo.C_Eco_Pools[iT,iS,5]=vo.C_Eco_Pools[iT,iS,5]-Affected_RootCoarse
            vo.C_Eco_Pools[iT,iS,6]=vo.C_Eco_Pools[iT,iS,6]-Affected_RootFine
               
            # Remove carbon from snag pools
            vo.C_Eco_Pools[iT,iS,11]=vo.C_Eco_Pools[iT,iS,11]-Affected_SnagStem
            vo.C_Eco_Pools[iT,iS,12]=vo.C_Eco_Pools[iT,iS,12]-Affected_SnagBranch
            
            # Remove stemwood merch volume
            vo.V_StemMerch[iT,iS]=(1-Biomass_Affected_Pct)*vo.V_StemMerch[iT,iS]
            
            #------------------------------------------------------------------
            # Carbon removed (sent to mill)
            #------------------------------------------------------------------
            
            # Merch biomass to mill - of the total amount of biomass affected, 
            # Add bark                        
            Biomass_Merch_Removed=Biomass_Merch_Removed_Pct*(Affected_StemMerch+Affected_Bark_Merch)
            vo.C_Eco_Pools[iT,iS,22]=vo.C_Eco_Pools[iT,iS,22]+Biomass_Merch_Removed
            
            # NonMerch biomass to mill - of the total amount of biomass affected, 
            # what fraction of non-merch biomass was sent to the mill?
            # - NonMerch = NonMerchStem + Foliage + Branch + Bark
            Biomass_NonMerch_Removed=Biomass_NonMerch_Removed_Pct*Affected_TotNonMerch
            vo.C_Eco_Pools[iT,iS,23]=vo.C_Eco_Pools[iT,iS,23]+Biomass_NonMerch_Removed
                                        
            # Snag stemwood to mill           
            SnagStem_Removed=Snags_Removed_Pct*Affected_SnagStem
            vo.C_Eco_Pools[iT,iS,24]=vo.C_Eco_Pools[iT,iS,24]+SnagStem_Removed
            
            # Snag branches to mill
            SnagBranch_Removed=Snags_Removed_Pct*Affected_SnagBranch
            vo.C_Eco_Pools[iT,iS,24]=vo.C_Eco_Pools[iT,iS,24]+SnagBranch_Removed
            
            #------------------------------------------------------------------
            # Carbon burned
            #------------------------------------------------------------------
            
            # Profile of emitted compounds (N2O emissions calculated upon export)
            tf_co2=psl.bCombFrac_CO2
            tf_ch4=psl.bCombFrac_CH4
            tf_co=psl.bCombFrac_CO
            
            # Merch biomass that is burned              
            Biomass_Merch_Burned_co2=tf_co2*Biomass_Merch_Burned_Pct*(Affected_StemMerch+Affected_Bark_Merch)
            vo.C_Eco_Pools[iT,iS,18]=vo.C_Eco_Pools[iT,iS,18]+Biomass_Merch_Burned_co2
            
            Biomass_Merch_Burned_ch4=tf_ch4*Biomass_Merch_Burned_Pct*(Affected_StemMerch+Affected_Bark_Merch)
            vo.C_Eco_Pools[iT,iS,19]=vo.C_Eco_Pools[iT,iS,19]+Biomass_Merch_Burned_ch4
            
            Biomass_Merch_Burned_co=tf_co*Biomass_Merch_Burned_Pct*(Affected_StemMerch+Affected_Bark_Merch)
            vo.C_Eco_Pools[iT,iS,20]=vo.C_Eco_Pools[iT,iS,20]+Biomass_Merch_Burned_co

            # NonMerch biomass that is burned            
            Biomass_NonMerch_Burned_co2=tf_co2*Biomass_NonMerch_Burned_Pct*Affected_TotNonMerch
            vo.C_Eco_Pools[iT,iS,18]=vo.C_Eco_Pools[iT,iS,18]+Biomass_NonMerch_Burned_co2
            
            Biomass_NonMerch_Burned_ch4=tf_ch4*Biomass_NonMerch_Burned_Pct*Affected_TotNonMerch
            vo.C_Eco_Pools[iT,iS,19]=vo.C_Eco_Pools[iT,iS,19]+Biomass_NonMerch_Burned_ch4
            
            Biomass_NonMerch_Burned_co=tf_co*Biomass_NonMerch_Burned_Pct*Affected_TotNonMerch
            vo.C_Eco_Pools[iT,iS,20]=vo.C_Eco_Pools[iT,iS,20]+Biomass_NonMerch_Burned_co
              
            # Snag stems that burn
            SnagStem_Burned_co2=tf_co2*Snags_Burned_Pct*Affected_SnagStem
            vo.C_Eco_Pools[iT,iS,18]=vo.C_Eco_Pools[iT,iS,18]+SnagStem_Burned_co2
            
            SnagStem_Burned_ch4=tf_ch4*Snags_Burned_Pct*Affected_SnagStem
            vo.C_Eco_Pools[iT,iS,19]=vo.C_Eco_Pools[iT,iS,19]+SnagStem_Burned_ch4
            
            SnagStem_Burned_co=tf_co*Snags_Burned_Pct*Affected_SnagStem
            vo.C_Eco_Pools[iT,iS,20]=vo.C_Eco_Pools[iT,iS,20]+SnagStem_Burned_co
            
            # Snag branches that burn
            SnagBranch_Burned_co2=tf_co2*Snags_Burned_Pct*Affected_SnagBranch
            vo.C_Eco_Pools[iT,iS,18]=vo.C_Eco_Pools[iT,iS,18]+SnagBranch_Burned_co2
            
            SnagBranch_Burned_ch4=tf_ch4*Snags_Burned_Pct*Affected_SnagBranch
            vo.C_Eco_Pools[iT,iS,19]=vo.C_Eco_Pools[iT,iS,19]+SnagBranch_Burned_ch4
            
            SnagBranch_Burned_co=tf_co*Snags_Burned_Pct*Affected_SnagBranch
            vo.C_Eco_Pools[iT,iS,20]=vo.C_Eco_Pools[iT,iS,20]+SnagBranch_Burned_co            
            
            #------------------------------------------------------------------
            # Check to confirm carbon in = carbon out
            #------------------------------------------------------------------
            
            # CarbonAffectedAG=affected_StemMerch+ \
            #    affected_Bark_Merch+ \
            #    affected_TotNonMerch+ \
            #    affected_SnagStem+ \
            #    affected_SnagBranch
            
            # CarbonRemoved=Biomass_Merch_Removed+Biomass_NonMerch_Removed+SnagStem_Removed+SnagBranch_Removed
            
            # CarbonBurned=Biomass_Merch_Burned_co2+Biomass_Merch_Burned_ch4+Biomass_Merch_Burned_co+ \
            #    Biomass_NonMerch_Burned_co2+Biomass_NonMerch_Burned_ch4+Biomass_NonMerch_Burned_co+ \
            #    SnagStem_Burned_co2+SnagStem_Burned_ch4+SnagStem_Burned_co+ \
            #    SnagBranch_Burned_co2+SnagBranch_Burned_ch4+SnagBranch_Burned_co
            
            # print([CarbonAffectedAG,CarbonRemoved,CarbonBurned,CarbonRemoved+CarbonBurned])
            
            #------------------------------------------------------------------
            # Carbon that is moved from biomass to snags
            #------------------------------------------------------------------
            
            # Merch biomass that is killed by wildfire
            vo.C_Eco_Pools[iT,iS,12]=vo.C_Eco_Pools[iT,iS,12]+Biomass_Merch_ToSnag_Pct*(Affected_StemMerch+Affected_Bark_Merch)
            
            # Non-merch biomass that is killed by wildfire
            vo.C_Eco_Pools[iT,iS,12]=vo.C_Eco_Pools[iT,iS,12]+Biomass_NonMerch_ToSnag_Pct*(Affected_StemNonMerch+Affected_Bark_NonMerch)
            
            #------------------------------------------------------------------
            # Carbon that is left on site (after felling)
            #------------------------------------------------------------------
                
            # Merch biomass that is cut, but left on site (transferred to DOM pools)            
            vo.C_Eco_Pools[iT,iS,9]=vo.C_Eco_Pools[iT,iS,9]+Biomass_Merch_LeftOnSite_Pct*(Affected_StemMerch+Affected_Bark_Merch)
            
            # Non-merch biomass that is left on site (transferred to DOM pools)
            vo.C_Eco_Pools[iT,iS,8]=vo.C_Eco_Pools[iT,iS,8]+Biomass_NonMerch_LeftOnSite_Pct*Affected_TotNonMerch
                                      
            # Foliage to very fast litter (not adjustable by user)
            vo.C_Eco_Pools[iT,iS,7]=vo.C_Eco_Pools[iT,iS,7]+1.0*Affected_Foliage
            
            # Transfer affected root biomass to DOM (not adjustable by user)        
            vo.C_Eco_Pools[iT,iS,8]=vo.C_Eco_Pools[iT,iS,8]+0.5*Affected_RootCoarse
            vo.C_Eco_Pools[iT,iS,14]=vo.C_Eco_Pools[iT,iS,14]+0.5*Affected_RootCoarse            
            vo.C_Eco_Pools[iT,iS,7]=vo.C_Eco_Pools[iT,iS,7]+0.5*Affected_RootFine
            vo.C_Eco_Pools[iT,iS,13]=vo.C_Eco_Pools[iT,iS,13]+0.5*Affected_RootFine            
            
            # Snag stemwood to medium litter pool
            vo.C_Eco_Pools[iT,iS,9]=vo.C_Eco_Pools[iT,iS,9]+Snags_LeftOnSite_Pct*Affected_SnagStem
            
            # Snag branches to fast litter
            vo.C_Eco_Pools[iT,iS,8]=vo.C_Eco_Pools[iT,iS,8]+Snags_LeftOnSite_Pct*Affected_SnagBranch            
              
            if meta['Biomass Module']=='TIPSY':
            
                # Update stand age
                if Biomass_Affected_Pct>0.75:
                    vo.A[iT,iS]=0            
               
                # Transition to new growth curve if there is a change in stand type
                if (vi.dh[iS].ID_Type[iDist[i]]!=meta['LUT Dist']['Aerial Fertilization']):
                    
                    if vi.dh[iS].ID_GrowthCurveM[iDist[i]]==1:
                        # Change to first growth curve
                        vi.GrowthCurveA['VStemMerch'][:,iS]=vi.GrowthCurve1['VStemMerch'][:,iS]
                        vi.GrowthCurveA['CStemMerch'][:,iS]=vi.GrowthCurve1['CStemMerch'][:,iS]
                        vi.GrowthCurveA['CStemNonMerch'][:,iS]=vi.GrowthCurve1['CStemNonMerch'][:,iS]
                        vi.GrowthCurveA['CFoliage'][:,iS]=vi.GrowthCurve1['CFoliage'][:,iS]
                        vi.GrowthCurveA['CBranch'][:,iS]=vi.GrowthCurve1['CBranch'][:,iS]
                        vi.GrowthCurveA['CBark'][:,iS]=vi.GrowthCurve1['CBark'][:,iS]
                    elif vi.dh[iS].ID_GrowthCurveM[iDist[i]]==2:
                        # Change to second growth curve
                        vi.GrowthCurveA['VStemMerch'][:,iS]=vi.GrowthCurve2['VStemMerch'][:,iS]
                        vi.GrowthCurveA['CStemMerch'][:,iS]=vi.GrowthCurve2['CStemMerch'][:,iS]
                        vi.GrowthCurveA['CStemNonMerch'][:,iS]=vi.GrowthCurve2['CStemNonMerch'][:,iS]
                        vi.GrowthCurveA['CFoliage'][:,iS]=vi.GrowthCurve2['CFoliage'][:,iS]
                        vi.GrowthCurveA['CBranch'][:,iS]=vi.GrowthCurve2['CBranch'][:,iS]
                        vi.GrowthCurveA['CBark'][:,iS]=vi.GrowthCurve2['CBark'][:,iS]
                    elif vi.dh[iS].ID_GrowthCurveM[iDist[i]]==3:
                        # Change to third growth curve
                        vi.GrowthCurveA['VStemMerch'][:,iS]=vi.GrowthCurve3['VStemMerch'][:,iS]
                        vi.GrowthCurveA['CStemMerch'][:,iS]=vi.GrowthCurve3['CStemMerch'][:,iS]
                        vi.GrowthCurveA['CStemNonMerch'][:,iS]=vi.GrowthCurve3['CStemNonMerch'][:,iS]
                        vi.GrowthCurveA['CFoliage'][:,iS]=vi.GrowthCurve3['CFoliage'][:,iS]
                        vi.GrowthCurveA['CBranch'][:,iS]=vi.GrowthCurve3['CBranch'][:,iS]
                        vi.GrowthCurveA['CBark'][:,iS]=vi.GrowthCurve3['CBark'][:,iS]
                    elif vi.dh[iS].ID_GrowthCurveM[iDist[i]]==4:
                        # Change to fourth growth curve
                        vi.GrowthCurveA['VStemMerch'][:,iS]=vi.GrowthCurve4['VStemMerch'][:,iS]
                        vi.GrowthCurveA['CStemMerch'][:,iS]=vi.GrowthCurve4['CStemMerch'][:,iS]
                        vi.GrowthCurveA['CStemNonMerch'][:,iS]=vi.GrowthCurve4['CStemNonMerch'][:,iS]
                        vi.GrowthCurveA['CFoliage'][:,iS]=vi.GrowthCurve4['CFoliage'][:,iS]
                        vi.GrowthCurveA['CBranch'][:,iS]=vi.GrowthCurve4['CBranch'][:,iS]
                        vi.GrowthCurveA['CBark'][:,iS]=vi.GrowthCurve4['CBark'][:,iS]
    
            #------------------------------------------------------------------
            # Calculate operational emissions from aerial fertilization
            #------------------------------------------------------------------
            
            if psl.bDist_Name[iDefPar]=='Aerial Fertilization':
                
                # Urea dose (kgUrea/ha)
                DoseUrea=psl.bDoseUrea
                
                # Nitrogen dose (kgN/ha)
                DoseN=psl.bUreaRatioThatIsNitrogen*DoseUrea

                # Emissions from manufacturing, transport prior to application 
                # and transport during application. This is derived from 
                # Sonne (2006) as 30% of the total emissions of 1.9 tCO2e/ha. 
                # It is higher than the work of Mel Scott (2010)
                E_Op=psl.bEmissionsFromManufacturingAndTransport

                # Volatilization: CO2 emissions following application (Tier 1 approach, IPCC 2006, 11.4.1)
                E_Vol=psl.bRatio_C_to_Urea*(DoseUrea/1000)*(1/psl.bRatio_C_to_CO2)

                # Denitrification: N2O emissions following application (Tier 1 approach, IPCC 2006, 11.4.1)
                E_Denit=psl.bFractionOfNAppliedThatIsEmittedAsN2OAsN*(DoseN/1000)* \
                    psl.bRatioN2OAsN_to_N2O*psl.bGWP_N2O_AR4
                
                # Total emissions (tCO2e/ha)
                E_Tot=E_Op+E_Vol+E_Denit
                
                # Convert to carbon (tC/ha)
                vo.C_E_Operations[iT,iS]=psl.bRatio_C_to_CO2*E_Tot
    
    #--------------------------------------------------------------------------
    # Populate emission variables with values from carbon structure
    #--------------------------------------------------------------------------
    
    vo.C_CombustCO2[iT,:]=vo.C_Eco_Pools[iT,:,18]        
    vo.C_CombustCH4[iT,:]=vo.C_Eco_Pools[iT,:,19]
    vo.C_CombustCO[iT,:]=vo.C_Eco_Pools[iT,:,20]
    vo.C_CombustN2O[iT,:]=vo.C_Eco_Pools[iT,:,21]
    vo.C_RemovedMerch[iT,:]=vo.C_Eco_Pools[iT,:,22]
    vo.C_RemovedNonMerch[iT,:]=vo.C_Eco_Pools[iT,:,23]
    vo.C_RemovedSnagStem[iT,:]=vo.C_Eco_Pools[iT,:,24]
        
    return vo

'''============================================================================
STAND PRODUCTS SECTOR
============================================================================'''

def ProductSector(iT,vi,vo,psl,meta):
    
    #             0     1     2      3      4      5         6       7       8      9          10             11         12          13                       14                          15                        16                           17      18     19
    # PoolNames=['SFH','MFH','Comm','Furn','Ship','Repairs','Other','Paper','Fuel','Firewood','EffluentPulp','DumpWood','DumpPaper','LandfillWoodDegradable','LandfillWoodNonDegradable','LandfillPaperDegradable','LandfillPaperNonDegradable','E_CO2','E_CH4','Cants']
            
    #--------------------------------------------------------------------------
    # Ecosystems --> Mills
    # Note: Mills act as transient reservoirs that are not tracked over time
    #--------------------------------------------------------------------------
    
    # Assume default parameters are applied, then update below if user has specified input variables
    Mill_Chips=psl.HWP_RemovedMerchToChipping*vo.C_RemovedMerch[iT,:] + \
                 psl.HWP_RemovedNonMerchToChipping*vo.C_RemovedNonMerch[iT,:] + \
                 psl.HWP_RemovedSnagStemToChipping*vo.C_RemovedSnagStem[iT,:]
    
    Mill_Chips_Pulp=psl.HWP_ChipsToPulp*Mill_Chips
    Mill_Chips_Fuel=psl.HWP_ChipsToFuel*Mill_Chips
    
    Mill_Lumber=psl.HWP_RemovedMerchToLumber*vo.C_RemovedMerch[iT,:] + \
                 psl.HWP_RemovedNonMerchToLumber*vo.C_RemovedNonMerch[iT,:] + \
                 psl.HWP_RemovedSnagStemToLumber*vo.C_RemovedSnagStem[iT,:]
    Mill_Plywood=psl.HWP_RemovedMerchToPlywood*vo.C_RemovedMerch[iT,:] + \
                 psl.HWP_RemovedNonMerchToPlywood*vo.C_RemovedNonMerch[iT,:] + \
                 psl.HWP_RemovedSnagStemToPlywood*vo.C_RemovedSnagStem[iT,:]
    Mill_OSB=psl.HWP_RemovedMerchToOSB*vo.C_RemovedMerch[iT,:] + \
                 psl.HWP_RemovedNonMerchToOSB*vo.C_RemovedNonMerch[iT,:] + \
                 psl.HWP_RemovedSnagStemToOSB*vo.C_RemovedSnagStem[iT,:]
    Mill_MDF=psl.HWP_RemovedMerchToMDF*vo.C_RemovedMerch[iT,:] + \
                 psl.HWP_RemovedNonMerchToMDF*vo.C_RemovedNonMerch[iT,:] + \
                 psl.HWP_RemovedSnagStemToMDF*vo.C_RemovedSnagStem[iT,:]
    Mill_Firewood=psl.HWP_RemovedMerchToFirewood*vo.C_RemovedMerch[iT,:] + \
                 psl.HWP_RemovedNonMerchToFirewood*vo.C_RemovedNonMerch[iT,:] + \
                 psl.HWP_RemovedSnagStemToFirewood*vo.C_RemovedSnagStem[iT,:]
    Mill_Cants=0*vo.C_RemovedMerch[iT,:] # Not considered among default in-use products
                 
    # Apply user-specified transfer percents
    for iS in range(meta['N Stand']):
        
        # Find custom disturbance inputs
        it=np.where((vi.dh[iS].Year==vi.tv[iT]) & (vi.dh[iS].Biomass_Affected_Pct>0) |
                (vi.dh[iS].Year==vi.tv[iT]) & (vi.dh[iS].Snags_Affected_Pct>0))[0]
        
        # Only continue if there are custom disturbance inputs
        if it.shape[0]==0: 
            continue
        
        # Run error if user inputs multiple custom disturbances in single year
        if it.shape[0]>1:
            print('Warning: Error - multiple custom disturbances in single year.')
            continue
                
        Mill_Chips_Pulp[iS]=vi.dh[iS].RemovedMerchToPulp_Pct[it]/100*vo.C_RemovedMerch[iT,iS] + \
                 vi.dh[iS].RemovedNonMerchToPulp_Pct[it]/100*vo.C_RemovedNonMerch[iT,iS] + \
                 vi.dh[iS].RemovedSnagStemToPulp_Pct[it]/100*vo.C_RemovedSnagStem[iT,iS]
        Mill_Chips_Fuel[iS]=vi.dh[iS].RemovedMerchToFuel_Pct[it]/100*vo.C_RemovedMerch[iT,iS] + \
                 vi.dh[iS].RemovedNonMerchToFuel_Pct[it]/100*vo.C_RemovedNonMerch[iT,iS] + \
                 vi.dh[iS].RemovedSnagStemToFuel_Pct[it]/100*vo.C_RemovedSnagStem[iT,iS]
        Mill_Lumber[iS]=vi.dh[iS].RemovedMerchToLumber_Pct[it]/100*vo.C_RemovedMerch[iT,iS] + \
                 vi.dh[iS].RemovedNonMerchToLumber_Pct[it]/100*vo.C_RemovedNonMerch[iT,iS] + \
                 vi.dh[iS].RemovedSnagStemToLumber_Pct[it]/100*vo.C_RemovedSnagStem[iT,iS]
        Mill_Plywood[iS]=vi.dh[iS].RemovedMerchToPlywood_Pct[it]/100*vo.C_RemovedMerch[iT,iS] + \
                 vi.dh[iS].RemovedNonMerchToPlywood_Pct[it]/100*vo.C_RemovedNonMerch[iT,iS] + \
                 vi.dh[iS].RemovedSnagStemToPlywood_Pct[it]/100*vo.C_RemovedSnagStem[iT,iS]
        Mill_OSB[iS]=vi.dh[iS].RemovedMerchToOSB_Pct[it]/100*vo.C_RemovedMerch[iT,iS] + \
                 vi.dh[iS].RemovedNonMerchToOSB_Pct[it]/100*vo.C_RemovedNonMerch[iT,iS] + \
                 vi.dh[iS].RemovedSnagStemToOSB_Pct[it]/100*vo.C_RemovedSnagStem[iT,iS]
        Mill_MDF[iS]=vi.dh[iS].RemovedMerchToMDF_Pct[it]/100*vo.C_RemovedMerch[iT,iS] + \
                 vi.dh[iS].RemovedNonMerchToMDF_Pct[it]/100*vo.C_RemovedNonMerch[iT,iS] + \
                 vi.dh[iS].RemovedSnagStemToMDF_Pct[it]/100*vo.C_RemovedSnagStem[iT,iS]
        Mill_Firewood[iS]=vi.dh[iS].RemovedMerchToFirewood_Pct[it]/100*vo.C_RemovedMerch[iT,iS] + \
                 vi.dh[iS].RemovedNonMerchToFirewood_Pct[it]/100*vo.C_RemovedNonMerch[iT,iS] + \
                 vi.dh[iS].RemovedSnagStemToFirewood_Pct[it]/100*vo.C_RemovedSnagStem[iT,iS]            
        Mill_Cants[iS]=vi.dh[iS].RemovedMerchToCants_Pct[it]/100*vo.C_RemovedMerch[iT,iS] + \
                 vi.dh[iS].RemovedNonMerchToCants_Pct[it]/100*vo.C_RemovedNonMerch[iT,iS] + \
                 vi.dh[iS].RemovedSnagStemToCants_Pct[it]/100*vo.C_RemovedSnagStem[iT,iS]                     
     
    #--------------------------------------------------------------------------
    # Mills --> In-use pools or other mills
    #--------------------------------------------------------------------------
            
    # Transfer mill fibre to single-family homes
    ip=0
    vo.C_Pro_Pools[iT,:,ip]=vo.C_Pro_Pools[iT-1,:,ip] + \
        psl.HWP_LumberToSFH*Mill_Lumber + \
        psl.HWP_PlywoodToSFH*Mill_Plywood + \
        psl.HWP_OSBToSFH*Mill_OSB + \
        psl.HWP_MDFToSFH*Mill_MDF
    
    # Transfer mill fibre to multi-family homes
    ip=1
    vo.C_Pro_Pools[iT,:,ip]=vo.C_Pro_Pools[iT-1,:,ip] + \
        psl.HWP_LumberToMFH*Mill_Lumber + \
        psl.HWP_PlywoodToMFH*Mill_Plywood + \
        psl.HWP_OSBToMFH*Mill_OSB + \
        psl.HWP_MDFToMFH*Mill_MDF
    
    # Transfer mill fibre to commercial
    ip=2
    vo.C_Pro_Pools[iT,:,ip]=vo.C_Pro_Pools[iT-1,:,ip] + \
        psl.HWP_LumberToCom*Mill_Lumber + \
        psl.HWP_PlywoodToCom*Mill_Plywood + \
        psl.HWP_OSBToCom*Mill_OSB + \
        psl.HWP_MDFToCom*Mill_MDF
    
    # Transfer mill fibre to furniture
    ip=3
    vo.C_Pro_Pools[iT,:,ip]=vo.C_Pro_Pools[iT-1,:,ip] + \
        psl.HWP_LumberToFurn*Mill_Lumber + \
        psl.HWP_PlywoodToFurn*Mill_Plywood + \
        psl.HWP_OSBToFurn*Mill_OSB + \
        psl.HWP_MDFToFurn*Mill_MDF
    
    # Transfer mill fibre to shipping
    ip=4
    vo.C_Pro_Pools[iT,:,ip]=vo.C_Pro_Pools[iT-1,:,ip] + \
        psl.HWP_LumberToShip*Mill_Lumber + \
        psl.HWP_PlywoodToShip*Mill_Plywood + \
        psl.HWP_OSBToShip*Mill_OSB + \
        psl.HWP_MDFToShip*Mill_MDF
    
    # Transfer mill fibre to repairs
    ip=5
    vo.C_Pro_Pools[iT,:,ip]=vo.C_Pro_Pools[iT-1,:,ip] + \
        psl.HWP_LumberToRepairs*Mill_Lumber + \
        psl.HWP_PlywoodToRepairs*Mill_Plywood + \
        psl.HWP_OSBToRepairs*Mill_OSB + \
        psl.HWP_MDFToRepairs*Mill_MDF
    
    # Transfer mill fibre to other
    ip=6
    vo.C_Pro_Pools[iT,:,ip]=vo.C_Pro_Pools[iT-1,:,ip] + \
        psl.HWP_LumberToOther*Mill_Lumber + \
        psl.HWP_PlywoodToOther*Mill_Plywood + \
        psl.HWP_OSBToOther*Mill_OSB + \
        psl.HWP_MDFToOther*Mill_MDF
    
    # Transfer mill fibre to pulp mill
    Mill_Pulp=Mill_Chips_Pulp + \
                 psl.HWP_LumberToPulp*Mill_Lumber + \
                 psl.HWP_PlywoodToPulp*Mill_Plywood
    
    # Transfer pulp mill fibre to paper
    ip=7
    vo.C_Pro_Pools[iT,:,ip]=vo.C_Pro_Pools[iT-1,:,ip] + \
                 psl.HWP_PulpToPaper*Mill_Pulp
    
    # Transfer mill fibre to fuel 
    ip=8
    vo.C_Pro_Pools[iT,:,ip]=vo.C_Pro_Pools[iT-1,:,ip] + \
                 Mill_Chips_Fuel + \
                 psl.HWP_LumberToFuel*Mill_Lumber + \
                 psl.HWP_PlywoodToFuel*Mill_Plywood + \
                 psl.HWP_OSBToFuel*Mill_OSB + \
                 psl.HWP_MDFToFuel*Mill_MDF
    
    # Transfer firewood to firewood pool
    ip=9
    vo.C_Pro_Pools[iT,:,ip]=vo.C_Pro_Pools[iT-1,:,ip] + Mill_Firewood
    
    # Transfer mill fibre to cants
    ip=19
    vo.C_Pro_Pools[iT,:,ip]=vo.C_Pro_Pools[iT-1,:,ip] + Mill_Cants
    
    # Transfer pulp mill carbon to pulp-mill effluent
    ip=10
    vo.C_Pro_Pools[iT,:,ip]=vo.C_Pro_Pools[iT-1,:,ip] + \
                 psl.HWP_PulpToEffluent*Mill_Pulp   
    
    #--------------------------------------------------------------------------
    # Update dump and landfill reservoirs
    #--------------------------------------------------------------------------

    ip=11
    vo.C_Pro_Pools[iT,:,ip]=vo.C_Pro_Pools[iT-1,:,ip]
    ip=12
    vo.C_Pro_Pools[iT,:,ip]=vo.C_Pro_Pools[iT-1,:,ip]
    ip=13
    vo.C_Pro_Pools[iT,:,ip]=vo.C_Pro_Pools[iT-1,:,ip]
    ip=14
    vo.C_Pro_Pools[iT,:,ip]=vo.C_Pro_Pools[iT-1,:,ip]    
    ip=15
    vo.C_Pro_Pools[iT,:,ip]=vo.C_Pro_Pools[iT-1,:,ip]
    ip=16
    vo.C_Pro_Pools[iT,:,ip]=vo.C_Pro_Pools[iT-1,:,ip]    
        
    #--------------------------------------------------------------------------
    # Single-family homes --> dump and landfill
    #--------------------------------------------------------------------------
        
    # Turnover
    ip=0
    Cretired=psl.HWP_SFH_tr*vo.C_Pro_Pools[iT-1,:,ip]
    
    # Remove carbon
    vo.C_Pro_Pools[iT,:,ip]=vo.C_Pro_Pools[iT,:,ip] - Cretired
    
    # Transfer carbon to dump wood
    ip=11
    vo.C_Pro_Pools[iT,:,ip]=vo.C_Pro_Pools[iT,:,ip] + psl.HWP_SFHToDumpWood*Cretired
    
    # Transfer carbon to landfill (degradable)
    ip=13
    vo.C_Pro_Pools[iT,:,ip]=vo.C_Pro_Pools[iT,:,ip] + psl.HWP_SFHToLandfillWood*psl.HWP_ToLandfillWoodDegradableFrac*Cretired
    
    # Transfer carbon to landfill (non-degradable)
    ip=14
    vo.C_Pro_Pools[iT,:,ip]=vo.C_Pro_Pools[iT,:,ip] + psl.HWP_SFHToLandfillWood*(1-psl.HWP_ToLandfillWoodDegradableFrac)*Cretired
    
    #--------------------------------------------------------------------------
    # Multi-family homes --> dump and landfill
    #--------------------------------------------------------------------------
    
    # Turnover
    ip=1
    Cretired=psl.HWP_MFH_tr*vo.C_Pro_Pools[iT,:,ip]
    
    # Remove carbon
    vo.C_Pro_Pools[iT,:,ip]=vo.C_Pro_Pools[iT,:,ip] - Cretired
    
    # Transfer carbon to dump wood
    ip=11
    vo.C_Pro_Pools[iT,:,ip]=vo.C_Pro_Pools[iT,:,ip] + psl.HWP_MFHToDumpWood*Cretired
    
    # Transfer carbon to landfill (degradable)
    ip=13
    vo.C_Pro_Pools[iT,:,ip]=vo.C_Pro_Pools[iT,:,ip] + psl.HWP_MFHToLandfillWood*psl.HWP_ToLandfillWoodDegradableFrac*Cretired
    
    # Transfer carbon to landfill (non-degradable)
    ip=14
    vo.C_Pro_Pools[iT,:,ip]=vo.C_Pro_Pools[iT,:,ip] + psl.HWP_MFHToLandfillWood*(1-psl.HWP_ToLandfillWoodDegradableFrac)*Cretired
    
    #--------------------------------------------------------------------------
    # Commercial building --> dump and landfill
    #--------------------------------------------------------------------------
    
    # Turnover
    ip=2
    Cretired=psl.HWP_Comm_tr*vo.C_Pro_Pools[iT-1,:,ip]
    
    # Remove carbon
    vo.C_Pro_Pools[iT,:,ip]=vo.C_Pro_Pools[iT,:,ip] - Cretired
    
    # Transfer carbon to dump wood
    ip=11
    vo.C_Pro_Pools[iT,:,ip]=vo.C_Pro_Pools[iT,:,ip] + psl.HWP_CommToDumpWood*Cretired
    
    # Transfer carbon to landfill (degradable)
    ip=13
    vo.C_Pro_Pools[iT,:,ip]=vo.C_Pro_Pools[iT,:,ip] + psl.HWP_CommToLandfillWood*psl.HWP_ToLandfillWoodDegradableFrac*Cretired
    
    # Transfer carbon to landfill (non-degradable)
    ip=14
    vo.C_Pro_Pools[iT,:,ip]=vo.C_Pro_Pools[iT,:,ip] + psl.HWP_CommToLandfillWood*(1-psl.HWP_ToLandfillWoodDegradableFrac)*Cretired
    
    #--------------------------------------------------------------------------
    # Furniture --> dump and landfill
    #--------------------------------------------------------------------------
    
    # Turnover
    ip=3
    Cretired=psl.HWP_Furn_tr*vo.C_Pro_Pools[iT-1,:,ip]
    
    # Remove carbon
    vo.C_Pro_Pools[iT,:,ip]=vo.C_Pro_Pools[iT,:,ip] - Cretired
    
    # Transfer carbon to dump wood
    ip=11
    vo.C_Pro_Pools[iT,:,ip]=vo.C_Pro_Pools[iT,:,ip] + psl.HWP_FurnToDumpWood*Cretired
    
    # Transfer carbon to landfill (degradable)
    ip=13
    vo.C_Pro_Pools[iT,:,ip]=vo.C_Pro_Pools[iT,:,ip] + psl.HWP_FurnToLandfillWood*psl.HWP_ToLandfillWoodDegradableFrac*Cretired
    
    # Transfer carbon to landfill (non-degradable)
    ip=14
    vo.C_Pro_Pools[iT,:,ip]=vo.C_Pro_Pools[iT,:,ip] + psl.HWP_FurnToLandfillWood*(1-psl.HWP_ToLandfillWoodDegradableFrac)*Cretired
    
    
    #--------------------------------------------------------------------------
    # Shipping --> dump and landfill
    #--------------------------------------------------------------------------
    
    # Turnover
    ip=4
    Cretired=psl.HWP_Ship_tr*vo.C_Pro_Pools[iT-1,:,ip]
    
    # Remove carbon
    vo.C_Pro_Pools[iT,:,ip]=vo.C_Pro_Pools[iT,:,ip] - Cretired
    
    # Transfer carbon to dump wood
    ip=11
    vo.C_Pro_Pools[iT,:,ip]=vo.C_Pro_Pools[iT,:,ip] + psl.HWP_ShipToDumpWood*Cretired
    
    # Transfer carbon to landfill (degradable)
    ip=13
    vo.C_Pro_Pools[iT,:,ip]=vo.C_Pro_Pools[iT,:,ip] + psl.HWP_ShipToLandfillWood*psl.HWP_ToLandfillWoodDegradableFrac*Cretired
    
    # Transfer carbon to landfill (non-degradable)
    ip=14
    vo.C_Pro_Pools[iT,:,ip]=vo.C_Pro_Pools[iT,:,ip] + psl.HWP_ShipToLandfillWood*(1-psl.HWP_ToLandfillWoodDegradableFrac)*Cretired
    
    #--------------------------------------------------------------------------
    # Repairs --> dump and landfill
    #--------------------------------------------------------------------------
    
    # Turnover
    ip=5
    Cretired=psl.HWP_Repairs_tr*vo.C_Pro_Pools[iT-1,:,ip]
    
    # Remove carbon
    vo.C_Pro_Pools[iT,:,ip]=vo.C_Pro_Pools[iT,:,ip] - Cretired
    
    # Transfer carbon to dump wood
    ip=11
    vo.C_Pro_Pools[iT,:,ip]=vo.C_Pro_Pools[iT,:,ip] + psl.HWP_RepairsToDumpWood*Cretired
    
    # Transfer carbon to landfill (degradble)
    ip=13
    vo.C_Pro_Pools[iT,:,ip]=vo.C_Pro_Pools[iT,:,ip] + psl.HWP_RepairsToLandfillWood*psl.HWP_ToLandfillWoodDegradableFrac*Cretired
    
    # Transfer carbon to landfill (non-degradable)
    ip=14
    vo.C_Pro_Pools[iT,:,ip]=vo.C_Pro_Pools[iT,:,ip] + psl.HWP_RepairsToLandfillWood*(1-psl.HWP_ToLandfillWoodDegradableFrac)*Cretired
    
    #--------------------------------------------------------------------------
    # Other --> dump and landfill
    #--------------------------------------------------------------------------
    
    # Turnover
    ip=6
    Cretired=psl.HWP_Other_tr*vo.C_Pro_Pools[iT-1,:,ip]
    
    # Remove carbon
    vo.C_Pro_Pools[iT,:,ip]=vo.C_Pro_Pools[iT,:,ip] - Cretired
    
    # Transfer carbon to dump wood
    ip=11
    vo.C_Pro_Pools[iT,:,ip]=vo.C_Pro_Pools[iT,:,ip] + psl.HWP_OtherToDumpWood*Cretired
    
    # Transfer carbon to landfill (degradble)
    ip=13
    vo.C_Pro_Pools[iT,:,ip]=vo.C_Pro_Pools[iT,:,ip] + psl.HWP_OtherToLandfillWood*psl.HWP_ToLandfillWoodDegradableFrac*Cretired
    
    # Transfer carbon to landfill (non-degradable)
    ip=14
    vo.C_Pro_Pools[iT,:,ip]=vo.C_Pro_Pools[iT,:,ip] + psl.HWP_OtherToLandfillWood*(1-psl.HWP_ToLandfillWoodDegradableFrac)*Cretired
    
    #--------------------------------------------------------------------------
    # Cants --> dump and landfill
    #--------------------------------------------------------------------------
    
    # Turnover
    ip=19
    Cretired=psl.HWP_Cants_tr*vo.C_Pro_Pools[iT-1,:,ip]
    
    # Remove carbon
    vo.C_Pro_Pools[iT,:,ip]=vo.C_Pro_Pools[iT,:,ip] - Cretired
    
    # Transfer carbon to dump wood
    ip=11
    vo.C_Pro_Pools[iT,:,ip]=vo.C_Pro_Pools[iT,:,ip] + psl.HWP_CantsToDumpWood*Cretired
    
    # Transfer carbon to landfill (degradble)
    ip=13
    vo.C_Pro_Pools[iT,:,ip]=vo.C_Pro_Pools[iT,:,ip] + psl.HWP_CantsToLandfillWood*psl.HWP_ToLandfillWoodDegradableFrac*Cretired
    
    # Transfer carbon to landfill (non-degradable)
    ip=14
    vo.C_Pro_Pools[iT,:,ip]=vo.C_Pro_Pools[iT,:,ip] + psl.HWP_CantsToLandfillWood*(1-psl.HWP_ToLandfillWoodDegradableFrac)*Cretired
    
    #--------------------------------------------------------------------------
    # Paper --> dump and landfill
    #--------------------------------------------------------------------------
    
    # Turnover (with adjustment for recycling)
    ip=7
    Cretired=(1-psl.HWP_PaperRecycleRate)*psl.HWP_Paper_tr*vo.C_Pro_Pools[iT,:,ip]
    
    # Remove carbon
    vo.C_Pro_Pools[iT,:,ip]=vo.C_Pro_Pools[iT,:,ip] - Cretired
    
    # Transfer to dump
    ip=12
    vo.C_Pro_Pools[iT,:,ip]=vo.C_Pro_Pools[iT,:,ip] + psl.HWP_PaperToDumpPaper*Cretired
    
    # Transfer to landfill (degradable)
    ip=13
    vo.C_Pro_Pools[iT,:,ip]=vo.C_Pro_Pools[iT,:,ip] + psl.HWP_PaperToLandfillPaper*psl.HWP_ToLandfillPaperDegradableFrac*Cretired
    
    # Transfer to landfill (non-degradable)
    ip=14
    vo.C_Pro_Pools[iT,:,ip]=vo.C_Pro_Pools[iT,:,ip] + psl.HWP_PaperToLandfillPaper*(1-psl.HWP_ToLandfillPaperDegradableFrac)*Cretired
        
    #--------------------------------------------------------------------------
    # Emissions from fuel combustion
    #--------------------------------------------------------------------------
    
    # Turnover
    ip=8
    Cretired=psl.HWP_Fuel_tr*vo.C_Pro_Pools[iT,:,ip]
    
    # Remove carbon
    vo.C_Pro_Pools[iT,:,ip]=vo.C_Pro_Pools[iT,:,ip] - Cretired
    
    # Emissions of CO2 from fuel use
    ip=17
    vo.C_Pro_Pools[iT,:,ip]=vo.C_Pro_Pools[iT,:,ip] + psl.HWP_FuelCombustionFracEmitCO2*Cretired
    
    # Emissions of CH4 from fuel use
    ip=18
    vo.C_Pro_Pools[iT,:,ip]=vo.C_Pro_Pools[iT,:,ip] + (1-psl.HWP_FuelCombustionFracEmitCO2)*Cretired
    
    #--------------------------------------------------------------------------
    # Emissions from firewood combustion
    #--------------------------------------------------------------------------
    
    # Turnover
    ip=9
    Cretired=psl.HWP_Firewood_tr*vo.C_Pro_Pools[iT,:,ip]
    
    # Remove carbon
    vo.C_Pro_Pools[iT,:,ip]=vo.C_Pro_Pools[iT,:,ip] - Cretired
    
    # Emissions of CO2
    ip=17
    vo.C_Pro_Pools[iT,:,ip]=vo.C_Pro_Pools[iT,:,ip] + Cretired
    
    # Emissions of CH4
    #ip=18
    #vo.C_Pro_Pools[iT,:,ip]=vo.C_Pro_Pools[iT,:,ip] + (1-psl.FuelCombustionFracEmitCO2)*Cretired
       
    #--------------------------------------------------------------------------
    # Emissions from pulp effluent
    #--------------------------------------------------------------------------
         
    # Emissions from pulp effluent (CO2 from aerobic decomposition)
    ip=10
    c_emitted=psl.HWP_EffluentPulp_tr*vo.C_Pro_Pools[iT,:,ip]
    
    # Remove emitted carbon
    vo.C_Pro_Pools[iT,:,ip]=vo.C_Pro_Pools[iT,:,ip]-c_emitted
    
    # Add emitted carbon to CO2 emission "pool"
    ip=17
    vo.C_Pro_Pools[iT,:,ip]=vo.C_Pro_Pools[iT,:,ip]+c_emitted    
    
    #--------------------------------------------------------------------------
    # Emissions from dump wood
    #--------------------------------------------------------------------------
    
    # Turnover
    ip=11
    c_emitted=psl.HWP_DumpWood_tr*vo.C_Pro_Pools[iT,:,ip]
    
    # Removal
    vo.C_Pro_Pools[iT,:,ip]=vo.C_Pro_Pools[iT,:,ip]-c_emitted
    
    # Add to emissions (CO2 emission from aerobic decomposition)
    ip=17
    vo.C_Pro_Pools[iT,:,ip]=vo.C_Pro_Pools[iT,:,ip]+c_emitted
    
    #--------------------------------------------------------------------------
    # Emissions from dump paper
    #--------------------------------------------------------------------------
    
    # Turnover
    ip=12
    c_emitted=psl.HWP_DumpPaper_tr*vo.C_Pro_Pools[iT,:,ip]
    
    # Removal
    vo.C_Pro_Pools[iT,:,ip]=vo.C_Pro_Pools[iT,:,ip]-c_emitted
    
    # Add to emissions (CO2 emission from aerobic decomposition)
    ip=17
    vo.C_Pro_Pools[iT,:,ip]=vo.C_Pro_Pools[iT,:,ip]+c_emitted
    
    #--------------------------------------------------------------------------
    # Emissions from landfill degradable wood
    #--------------------------------------------------------------------------
               
    # Turnover
    ip=13
    c_emitted=psl.HWP_LandfillWoodDegradable_tr*vo.C_Pro_Pools[iT,:,ip]
    
    # Removal
    vo.C_Pro_Pools[iT,:,ip]=vo.C_Pro_Pools[iT,:,ip]-c_emitted
    
    # Add to emissions (50% CO2 emissions during anaerobic decomposition)
    ip=17
    vo.C_Pro_Pools[iT,:,ip]=vo.C_Pro_Pools[iT,:,ip]+psl.HWP_LandfillDegradableFracEmitCO2*c_emitted
    
    # Add to emissions (50% "potential" CH4 emissions during anaerobic decomposition)
    E_ch4_pot=(1-psl.HWP_LandfillDegradableFracEmitCO2)*c_emitted
    
    # Adjustment for proportion of degradable landfills with gas collection systems, 
    # efficiency of system, and methane oxided to CO2 from the landfill cover
    ch4_emitted=c_emitted*((1-psl.HWP_LandfillMethaneEmit_GasColSysProp)-psl.HWP_LandfillMethaneOxidizedToCO2*(1-psl.HWP_LandfillMethaneEmit_GasColSysProp)) + \
        c_emitted*psl.HWP_LandfillMethaneEmit_GasColSysProp*((1-psl.HWP_LandfillMethaneEmit_GasColSysEffic)-psl.HWP_LandfillMethaneOxidizedToCO2*(1-psl.HWP_LandfillMethaneEmit_GasColSysProp))
        
    ip=18
    vo.C_Pro_Pools[iT,:,ip]=vo.C_Pro_Pools[iT,:,ip]+ch4_emitted
    
    #--------------------------------------------------------------------------
    # Emissions from landfill degradable paper
    #--------------------------------------------------------------------------
    
    # Turnover
    ip=15
    c_emitted=psl.HWP_LandfillWoodDegradable_tr*vo.C_Pro_Pools[iT,:,ip]
    
    # Removal
    vo.C_Pro_Pools[iT,:,ip]=vo.C_Pro_Pools[iT,:,ip]-c_emitted
    
    # Add to emissions (50% CO2 emissions during anaerobic decomposition)
    ip=17
    vo.C_Pro_Pools[iT,:,ip]=vo.C_Pro_Pools[iT,:,ip]+psl.HWP_LandfillDegradableFracEmitCO2*c_emitted
    
    # Add to emissions (50% "potential" CH4 emissions during anaerobic decomposition)
    E_ch4_pot=(1-psl.HWP_LandfillDegradableFracEmitCO2)*c_emitted
    
    # Adjustment for proportion of degradable landfills with gas collection systems, 
    # efficiency of system, and methane oxided to CO2 from the landfill cover
    ch4_emitted=c_emitted*((1-psl.HWP_LandfillMethaneEmit_GasColSysProp)-psl.HWP_LandfillMethaneOxidizedToCO2*(1-psl.HWP_LandfillMethaneEmit_GasColSysProp)) + \
        c_emitted*psl.HWP_LandfillMethaneEmit_GasColSysProp*((1-psl.HWP_LandfillMethaneEmit_GasColSysEffic)-psl.HWP_LandfillMethaneOxidizedToCO2*(1-psl.HWP_LandfillMethaneEmit_GasColSysProp))
        
    ip=18
    vo.C_Pro_Pools[iT,:,ip]=vo.C_Pro_Pools[iT,:,ip]+ch4_emitted
    
    return vo