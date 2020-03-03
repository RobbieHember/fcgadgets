
# Import python modules
import os
import time
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

def BiomassFromTIPSY(iScn,iT,vi,vo,psl,meta,iEP):
    
    # Update stand age
    vo['A'][iT,:]=vo['A'][iT-1,:]+1
     
    #--------------------------------------------------------------------------
    # Net growth of aboveground biomass
    #--------------------------------------------------------------------------
    
    iAge=np.minimum(vo['A'][iT,:],vi['GCA'].shape[0])-1
    iAge=iAge.astype(int)
       
    ng=np.zeros((meta['N Stand'],6))
    for i in range(meta['N Stand']):
        ng[i,:]=vi['GCA'][iAge[i],i,:]
    
    vo['C_G_Net'][iT,:,0:5]=ng[:,0:5].astype(float)/meta['Scale Factor GC']
    vo['V_StemMerch'][iT,:]=vo['V_StemMerch'][iT-1,:]+ng[:,5].astype(float)/meta['Scale Factor GC']
    
    # Apply optional factor to net growth        
    if meta['Scenario'][iScn]['Status Net Growth Factor']=='On':
        net_growth_factor=meta['Scenario'][iScn]['Net Growth Factor'][iT]
        if net_growth_factor!=0:
            # Normal operation
            vo['C_G_Net'][iT,:,:]=net_growth_factor*vo['C_G_Net'][iT,:,:]            

    # Add net growth to stemwood biomass pools
    vo['C_Eco_Pools'][iT,:,0]=np.maximum(0,vo['C_Eco_Pools'][iT-1,:,0]+vo['C_G_Net'][iT,:,0])
    vo['C_Eco_Pools'][iT,:,1]=np.maximum(0,vo['C_Eco_Pools'][iT-1,:,1]+vo['C_G_Net'][iT,:,1])
     
    #--------------------------------------------------------------------------
    # Problem with TIPSY foliage
    #--------------------------------------------------------------------------
    
    # Foliage - TIPSY foliage appears to be really low, use alternative equations
    Stemwood=np.maximum(0.01,vo['C_Eco_Pools'][iT,:,iEP['StemMerch']]+0.2*vo['C_Eco_Pools'][iT,:,iEP['StemNonMerch']])
    vo['C_Eco_Pools'][iT,:,iEP['Foliage']]=np.maximum(0,psl.bASL_StemToF1*Stemwood**psl.bASL_StemToF2)
        
    # We have redefined foliage growth, so make sure it is consistent
    vo['C_G_Net'][iT,:,iEP['Foliage']]=vo['C_Eco_Pools'][iT,:,iEP['Foliage']]-vo['C_Eco_Pools'][iT-1,:,iEP['Foliage']]        
    
    #--------------------------------------------------------------------------
    # End problem with TIPSY foliage
    #--------------------------------------------------------------------------
    
    # *** OLD: Add net growth to foliage ***
    #vo['C_Eco_Pools'][iT,:,iEP['Foliage']]=np.maximum(0,vo['C_Eco_Pools'][iT-1,:,iEP['Foliage']]+vo['C_G_Net'][iT,:,iEP['Foliage']])
    
    # Add net growth to branch
    vo['C_Eco_Pools'][iT,:,iEP['Branch']]=np.maximum(0,vo['C_Eco_Pools'][iT-1,:,iEP['Branch']]+vo['C_G_Net'][iT,:,iEP['Branch']])
    
    # Add net growth to bark
    vo['C_Eco_Pools'][iT,:,iEP['Bark']]=np.maximum(0,vo['C_Eco_Pools'][iT-1,:,iEP['Bark']]+vo['C_G_Net'][iT,:,iEP['Bark']])
        
    # Calculate aboveground biomass
    vo['C_BiomassAG'][iT,:]=np.sum(vo['C_Eco_Pools'][iT,:,0:5],axis=1)
  
    # Total root biomass (Li et al. 2003, Eq. 4)
    BiomassRootTotal=0.222*vo['C_BiomassAG'][iT,:]
        
    # Fine root biomass (Li et al. 2003, Eq. 6)
    Pf=0.072+0.354*np.exp(-0.06*(2*BiomassRootTotal))
    vo['C_Eco_Pools'][iT,:,6]=Pf*BiomassRootTotal
  
    # Coarse root biomass
    Pc=1-Pf
    vo['C_Eco_Pools'][iT,:,iEP['RootCoarse']]=np.maximum(0,Pc*BiomassRootTotal)
  
    # Calculate net growth of roots
    vo['C_G_Net'][iT,:,iEP['RootCoarse']]=vo['C_Eco_Pools'][iT,:,iEP['RootCoarse']]-vo['C_Eco_Pools'][iT-1,:,iEP['RootCoarse']]
    vo['C_G_Net'][iT,:,iEP['RootFine']]=vo['C_Eco_Pools'][iT,:,iEP['RootFine']]-vo['C_Eco_Pools'][iT-1,:,iEP['RootFine']]
    
    #--------------------------------------------------------------------------
    # Biomass turnover
    #--------------------------------------------------------------------------
    
    # Calculate biomass loss due to mortality
    vo['C_M_Reg'][iT,:,0]=psl.bTR_StemMerch[0,:]*vo['C_Eco_Pools'][iT,:,0]
    vo['C_M_Reg'][iT,:,1]=psl.bTR_StemMerch[0,:]*vo['C_Eco_Pools'][iT,:,1]
    vo['C_M_Reg'][iT,:,2]=psl.bTR_StemMerch[0,:]*vo['C_Eco_Pools'][iT,:,2]
    vo['C_M_Reg'][iT,:,3]=psl.bTR_StemMerch[0,:]*vo['C_Eco_Pools'][iT,:,3]
    vo['C_M_Reg'][iT,:,4]=psl.bTR_StemMerch[0,:]*vo['C_Eco_Pools'][iT,:,4]
    vo['C_M_Reg'][iT,:,5]=psl.bTR_StemMerch[0,:]*vo['C_Eco_Pools'][iT,:,5]
    vo['C_M_Reg'][iT,:,6]=psl.bTR_StemMerch[0,:]*vo['C_Eco_Pools'][iT,:,6]
    
    # Apply optional factor to mortality
    if meta['Scenario'][iScn]['Status Mortality Factor']=='On':
        mortality_factor=meta['Scenario'][iScn]['Mortality Factor'][iT]
        if mortality_factor!=0:        
            # Normal operation
            vo['C_M_Reg'][iT,:,:]=net_growth_factor*vo['C_M_Reg'][iT,:,:]
    
    #--------------------------------------------------------------------------
    # If TIPSY indicates negative net growth (e.g., alder break-up), revise so
    # that net growth during the period of negative values equals the value in the
    # previous timestep.
    #--------------------------------------------------------------------------
    
    # Keep track of what regular mortality would be before it is affected by 
    # correction for catastrophic mortality
    C_M_Reg=vo['C_M_Reg'][iT,:,:].copy()
    
    # Define a threshold level of negative growth so it isn't triggered by noise
    NegGrowthForRealz=-0.25
    
    # Find stands with negative net growth    
    iNegNetG=np.where(vo['C_G_Net'][iT,:,0]<NegGrowthForRealz)[0]
    
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
    
    # No littervall from stemwood
    #vo['C_LF'][iT,:,1]=psl.bTR_StemNonMerch[0,:]*vo['C_Eco_Pools'][iT,:,1]
    
    # Calculate foliage biomass turnover due to litterfall
    vo['C_LF'][iT,:,iEP['Foliage']]=psl.bTR_Foliage[0,:]*vo['C_Eco_Pools'][iT,:,iEP['Foliage']]
    
    # Calculate branch biomass turnover due to litterfall
    vo['C_LF'][iT,:,iEP['Branch']]=psl.bTR_Branch[0,:]*vo['C_Eco_Pools'][iT,:,iEP['Branch']]
    
    # Calculate bark biomass turnover due to litterfall
    vo['C_LF'][iT,:,iEP['Bark']]=psl.bTR_Bark[0,:]*vo['C_Eco_Pools'][iT,:,iEP['Bark']]
    
    # Calculate coarse root biomass turnover due to litterfall
    vo['C_LF'][iT,:,iEP['RootCoarse']]=psl.bTR_RootCoarse[0,:]*vo['C_Eco_Pools'][iT,:,iEP['RootCoarse']]
    
    # Calculate fine root biomass turnover due to litterfall
    vo['C_LF'][iT,:,iEP['RootFine']]=psl.bTR_RootFine[0,:]*vo['C_Eco_Pools'][iT,:,iEP['RootFine']]
        
    #--------------------------------------------------------------------------
    # Update summary variables
    #--------------------------------------------------------------------------
    
    # Update gross growth
    vo['C_G_Gross'][iT,:,:]=vo['C_G_Net'][iT,:,:]+C_M_Reg
    
    # Update NPP
    vo['C_NPP'][iT,:,:]=vo['C_G_Net'][iT,:,:]+C_M_Reg+vo['C_LF'][iT,:,:]
    
    return vo

'''============================================================================
BIOMASS DYNAMICS FROM TIPSY
*** THIS IS A MODIFIED MODULE VERSION USED FOR A SPECIFIC STUDY **
============================================================================'''

def BiomassFromTIPSY_SpecialAdjustments_EP703(iScn,iT,vi,vo,psl,meta,iEP):
    
    # Update stand age
    vo['A'][iT,:]=vo['A'][iT-1,:]+1
     
    #--------------------------------------------------------------------------
    # Net growth of aboveground biomass
    #--------------------------------------------------------------------------
    
    iAge=np.minimum(vo['A'][iT,:],vi['GCA'].shape[0])-1
    iAge=iAge.astype(int)
       
    ng=np.zeros((meta['N Stand'],6))
    for i in range(meta['N Stand']):
        ng[i,:]=vi['GCA'][iAge[i],i,:]
    
    vo['C_G_Net'][iT,:,0:5]=ng[:,0:5].astype(float)/meta['Scale Factor GC']
    vo['V_StemMerch'][iT,:]=vo['V_StemMerch'][iT-1,:]+ng[:,5].astype(float)/meta['Scale Factor GC']
    
    # Apply optional factor to net growth        
    if meta['Scenario'][iScn]['Status Net Growth Factor']=='On':
        net_growth_factor=meta['Scenario'][iScn]['Net Growth Factor'][iT]
        if net_growth_factor!=0:            
            
            # Normal operation
            # *** Don't use ***
            #vo['C_G_Net'][iT,:,:]=net_growth_factor*vo['C_G_Net'][iT,:,:]
        
            # *** For EP703 fertilization study ***
            
            # Apply a 8% reduction for wood density
            net_growth_factor=0.92*net_growth_factor
            
            response_ratio=net_growth_factor/vo['C_G_Net'][iT,0,0]
            
            # Assume 90% of total stemwood goes to merch
            fMerch=0.9; fNonMerch=0.1
            
            # Median ratios for Brix 1983, Mitchell et al. 1996, Gower et al. 1992
            fF=8.1; fBr=5.3; fBk=1.3
            
            # From broader sample of forests
            #fF=1.35; fBr=1.85; fBk=2.01
            
            vo['C_G_Net'][iT,:,0]=vo['C_G_Net'][iT,:,0]+fMerch*net_growth_factor
            vo['C_G_Net'][iT,:,1]=vo['C_G_Net'][iT,:,1]+fNonMerch*net_growth_factor
            vo['C_G_Net'][iT,:,2]=vo['C_G_Net'][iT,:,2]+vo['C_G_Net'][iT,:,2]*fF*response_ratio
            vo['C_G_Net'][iT,:,3]=vo['C_G_Net'][iT,:,3]+vo['C_G_Net'][iT,:,3]*fBr*response_ratio
            vo['C_G_Net'][iT,:,4]=vo['C_G_Net'][iT,:,4]+vo['C_G_Net'][iT,:,4]*fBk*response_ratio        

    # Add net growth to stemwood biomass pools
    vo['C_Eco_Pools'][iT,:,0]=np.maximum(0,vo['C_Eco_Pools'][iT-1,:,0]+vo['C_G_Net'][iT,:,0])
    vo['C_Eco_Pools'][iT,:,1]=np.maximum(0,vo['C_Eco_Pools'][iT-1,:,1]+vo['C_G_Net'][iT,:,1])
     
    #--------------------------------------------------------------------------
    # Problem with TIPSY foliage
    #--------------------------------------------------------------------------
    
    # Foliage - TIPSY foliage appears to be really low, use alternative equations
    Stemwood=np.maximum(0.01,vo['C_Eco_Pools'][iT,:,iEP['StemMerch']]+0.2*vo['C_Eco_Pools'][iT,:,iEP['StemNonMerch']])
    vo['C_Eco_Pools'][iT,:,iEP['Foliage']]=np.maximum(0,psl.bASL_StemToF1*Stemwood**psl.bASL_StemToF2)
        
    # We have redefined foliage growth, so make sure it is consistent
    vo['C_G_Net'][iT,:,iEP['Foliage']]=vo['C_Eco_Pools'][iT,:,iEP['Foliage']]-vo['C_Eco_Pools'][iT-1,:,iEP['Foliage']]        
    
    #--------------------------------------------------------------------------
    # End problem with TIPSY foliage
    #--------------------------------------------------------------------------
    
    # *** OLD: Add net growth to foliage ***
    #vo['C_Eco_Pools'][iT,:,iEP['Foliage']]=np.maximum(0,vo['C_Eco_Pools'][iT-1,:,iEP['Foliage']]+vo['C_G_Net'][iT,:,iEP['Foliage']])
    
    # Add net growth to branch
    vo['C_Eco_Pools'][iT,:,iEP['Branch']]=np.maximum(0,vo['C_Eco_Pools'][iT-1,:,iEP['Branch']]+vo['C_G_Net'][iT,:,iEP['Branch']])
    
    # Add net growth to bark
    vo['C_Eco_Pools'][iT,:,iEP['Bark']]=np.maximum(0,vo['C_Eco_Pools'][iT-1,:,iEP['Bark']]+vo['C_G_Net'][iT,:,iEP['Bark']])
        
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
  
    # Calculate net growth of roots
    vo['C_G_Net'][iT,:,iEP['RootCoarse']]=vo['C_Eco_Pools'][iT,:,iEP['RootCoarse']]-vo['C_Eco_Pools'][iT-1,:,iEP['RootCoarse']]
    vo['C_G_Net'][iT,:,iEP['RootFine']]=vo['C_Eco_Pools'][iT,:,iEP['RootFine']]-vo['C_Eco_Pools'][iT-1,:,iEP['RootFine']]
    
    # *** Apply a reduction in total root biomass growth
    Gnet_RC=0.61*(vo['C_Eco_Pools'][iT,:,iEP['RootCoarse']]-vo['C_Eco_Pools'][iT-1,:,iEP['RootCoarse']])
    #Gnet_RF=0.61*(vo['C_Eco_Pools'][iT,:,iEP['RootFine']]-vo['C_Eco_Pools'][iT-1,:,iEP['RootFine']])
    
    # *** Further reduce fine root growth by 12.8%. Shift the reduction to coarse
    # root growth so that the total root growth still equals the overall 68%
    # reduction to root growth. ***
    #ReductionOfFineRootGrowth=0.75*Gnet_RF
    #Gnet_RF=Gnet_RF-ReductionOfFineRootGrowth
    #Gnet_RC=Gnet_RC+ReductionOfFineRootGrowth
    
    vo['C_G_Net'][iT,:,iEP['RootCoarse']]=Gnet_RC
    #vo['C_G_Net'][iT,:,iEP['RootFine']]=Gnet_RF  
    
    # *** Recalculate pools and fluxes after reduction in growth ***
    vo['C_Eco_Pools'][iT,:,iEP['RootCoarse']]=vo['C_Eco_Pools'][iT-1,:,iEP['RootCoarse']]+Gnet_RC
    #vo['C_Eco_Pools'][iT,:,iEP['RootFine']]=vo['C_Eco_Pools'][iT-1,:,iEP['RootFine']]+Gnet_RF
        
    #--------------------------------------------------------------------------
    # Biomass turnover
    #--------------------------------------------------------------------------
    
    # Calculate biomass loss due to mortality
    vo['C_M_Reg'][iT,:,0]=psl.bTR_StemMerch[0,:]*vo['C_Eco_Pools'][iT,:,0]
    vo['C_M_Reg'][iT,:,1]=psl.bTR_StemMerch[0,:]*vo['C_Eco_Pools'][iT,:,1]
    vo['C_M_Reg'][iT,:,2]=psl.bTR_StemMerch[0,:]*vo['C_Eco_Pools'][iT,:,2]
    vo['C_M_Reg'][iT,:,3]=psl.bTR_StemMerch[0,:]*vo['C_Eco_Pools'][iT,:,3]
    vo['C_M_Reg'][iT,:,4]=psl.bTR_StemMerch[0,:]*vo['C_Eco_Pools'][iT,:,4]
    vo['C_M_Reg'][iT,:,5]=psl.bTR_StemMerch[0,:]*vo['C_Eco_Pools'][iT,:,5]
    vo['C_M_Reg'][iT,:,6]=psl.bTR_StemMerch[0,:]*vo['C_Eco_Pools'][iT,:,6]
    
    # Apply optional factor to mortality
    if meta['Scenario'][iScn]['Status Mortality Factor']=='On':
        mortality_factor=meta['Scenario'][iScn]['Mortality Factor'][iT]
        if mortality_factor!=0:            
        
            # Normal operation
            # *** Don't use ***
            #vo['C_M_Reg'][iT,:,:]=net_growth_factor*vo['C_M_Reg'][iT,:,:]
        
            # *** For EP703 fertilization study ***
            
            # Apply a 8% reduction for wood density
            net_growth_factor=0.92*net_growth_factor
            
            response_ratio=mortality_factor/vo['C_M_Reg'][iT,0,0]
            
            # Assume 90% of total stemwood goes to merch
            fMerch=0.9; fNonMerch=0.1
            
            # Median ratios for Brix 1983, Mitchell et al. 1996, Gower et al. 1992
            fF=8.1; fBr=5.3; fBk=1.3
            
            # From broader sample of forests
            #fF=1.35; fBr=1.85; fBk=2.01
            
            vo['C_M_Reg'][iT,:,0]=vo['C_M_Reg'][iT,:,0]+fMerch*net_growth_factor
            vo['C_M_Reg'][iT,:,1]=vo['C_M_Reg'][iT,:,1]+fNonMerch*net_growth_factor
            vo['C_M_Reg'][iT,:,2]=vo['C_M_Reg'][iT,:,2]+vo['C_M_Reg'][iT,:,2]*fF*response_ratio
            vo['C_M_Reg'][iT,:,3]=vo['C_M_Reg'][iT,:,3]+vo['C_M_Reg'][iT,:,3]*fBr*response_ratio
            vo['C_M_Reg'][iT,:,4]=vo['C_M_Reg'][iT,:,4]+vo['C_M_Reg'][iT,:,4]*fBk*response_ratio
            
            # Roots
            fR=0.61
            vo['C_M_Reg'][iT,:,5]=vo['C_M_Reg'][iT,:,5]+vo['C_M_Reg'][iT,:,5]*fR*response_ratio
            fR=1
            vo['C_M_Reg'][iT,:,6]=vo['C_M_Reg'][iT,:,6]+vo['C_M_Reg'][iT,:,6]*fR*response_ratio
    
    #--------------------------------------------------------------------------
    # If TIPSY indicates negative net growth (e.g., alder break-up), revise so
    # that net growth during the period of negative values equals the value in the
    # previous timestep.
    #--------------------------------------------------------------------------
    
    # Keep track of what regular mortality would be before it is affected by 
    # correction for catastrophic mortality
    C_M_Reg=vo['C_M_Reg'][iT,:,:].copy()
        
    # Define a threshold level of negative growth so it isn't triggered by noise
    NegGrowthForRealz=-0.25
    
    # Find stands with negative net growth    
    iNegNetG=np.where(vo['C_G_Net'][iT,:,0]<NegGrowthForRealz)[0]
    
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
    
    # No littervall from stemwood
    #vo['C_LF'][iT,:,1]=psl.bTR_StemNonMerch[0,:]*vo['C_Eco_Pools'][iT,:,1]
    
    # Modifier for EP703 study
    if (iScn==1) | (iScn==2) | (iScn==3) | (iScn==5) | (iScn==6) | (iScn==7):
        if meta['Time'][iT]>=2020:
            mod=1.0
        else:
            mod=1.0
    else:
        mod=1.0
    
    # Calculate foliage biomass turnover due to litterfall
    vo['C_LF'][iT,:,iEP['Foliage']]=mod*psl.bTR_Foliage[0,:]*vo['C_Eco_Pools'][iT,:,iEP['Foliage']]
    
    # Calculate branch biomass turnover due to litterfall
    vo['C_LF'][iT,:,iEP['Branch']]=mod*psl.bTR_Branch[0,:]*vo['C_Eco_Pools'][iT,:,iEP['Branch']]
    
    # Calculate bark biomass turnover due to litterfall
    vo['C_LF'][iT,:,iEP['Bark']]=mod*psl.bTR_Bark[0,:]*vo['C_Eco_Pools'][iT,:,iEP['Bark']]
    
    # Calculate coarse root biomass turnover due to litterfall
    vo['C_LF'][iT,:,iEP['RootCoarse']]=mod*psl.bTR_RootCoarse[0,:]*vo['C_Eco_Pools'][iT,:,iEP['RootCoarse']]
    
    # Calculate fine root biomass turnover due to litterfall
    vo['C_LF'][iT,:,iEP['RootFine']]=mod*psl.bTR_RootFine[0,:]*vo['C_Eco_Pools'][iT,:,iEP['RootFine']]
        
    #--------------------------------------------------------------------------
    # Update summary variables
    #--------------------------------------------------------------------------
    
    # Update gross growth
    vo['C_G_Gross'][iT,:,:]=vo['C_G_Net'][iT,:,:]+C_M_Reg
    
    # Update NPP
    vo['C_NPP'][iT,:,:]=vo['C_G_Net'][iT,:,:]+C_M_Reg+vo['C_LF'][iT,:,:]
    
    return vo

'''============================================================================
ANNUAL TREE BIOMASS DYNAMICS
============================================================================'''

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
        
        iDist=np.where((vi['DH'][iS]['Year']==vi['tv'][iT]) & (np.isnan(vi['DH'][iS]['Severity'])==False))[0]
        if iDist.size!=0:            
            
            iLive=np.where(np.isnan(tl_Csw[iT,:])==False)[0]
            
            nKill=int(vi['DH'][iS]['Severity'][iDist]/100*iLive.size)
                        
            if vi['DH'][iS]['ID_Type'][iDist]==meta['LUT Dist']['Wildfire']:
                # Wildfire
                iKill=np.arange(0,nKill)
                tl_N_M_Inv_Fir[iT,iLive[iKill]]=1
                vo['A'][iT,iS]=0
                
            elif vi['DH'][iS]['ID_Type'][iDist]==meta['LUT Dist']['Harvest']:
                iKill=np.arange(0,nKill)
                tl_N_M_Inv_Har[iT,iLive[iKill]]=1
                vo['A'][iT,iS]=0
                               
            elif (vi['DH'][iS]['ID_Type'][iDist]==meta['LUT Dist']['Beetles']) | (vi['DH'][iS]['ID_Type'][iDist]==meta['LUT Dist']['Defoliators']):
                # Insects
                iKill=np.arange(0,nKill)
                tl_N_M_Inv_Ins[iT,iLive[iKill]]=1
                
            elif (vi['DH'][iS]['ID_Type'][iDist]==meta['LUT Dist']['Planting']):
                # Do nothing
                pass
                
            elif (vi['DH'][iS]['ID_Type'][iDist]==meta['LUT Dist']['Fertilization Aerial']):
                # Do nothing
                pass
            
            elif vi['DH'][iS]['ID_Type'][iDist]==meta['LUT Dist']['Wind']:
                # Wind
                iKill=np.arange(0,nKill)
                tl_N_M_Inv_Win[iT,iLive[iKill]]=1
                vo['A'][iT,iS]=0
                
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

'''============================================================================
STAND DEAD ORGANIC MATTER DYNAMICS
============================================================================'''    

def DeadOrganicMatterDynamics(iT,vi,vo,psl,iEP):
     
    #--------------------------------------------------------------------------
    # Flux of carbon between biomass components and dead organic matter components
    #--------------------------------------------------------------------------
    
    # Transfer biomass turnover to very fast litter pool
    vo['C_Eco_Pools'][iT,:,iEP['LitterVF']]=vo['C_Eco_Pools'][iT-1,:,iEP['LitterVF']]+ \
        psl.bIPF_FoliageLitToLitterVF*vo['C_LF'][iT,:,iEP['Foliage']]+ \
        psl.bIPF_RootFineLitToLitterVF*vo['C_LF'][iT,:,iEP['RootFine']]+ \
        psl.bIPF_FoliageMorToLitterVF*vo['C_M_Reg'][iT,:,iEP['Foliage']]
    
    # Transfer biomass turnover to fast litter pool
    vo['C_Eco_Pools'][iT,:,iEP['LitterF']]=vo['C_Eco_Pools'][iT-1,:,iEP['LitterF']]+ \
        psl.bIPF_BranchLitToLitterF*vo['C_LF'][iT,:,iEP['Branch']]+ \
        psl.bIPF_BarkLitToLitterF*vo['C_LF'][iT,:,iEP['Bark']]+ \
        psl.bIPF_RootCoarseLitToLitterF*vo['C_LF'][iT,:,iEP['RootCoarse']]+ \
        psl.bIPF_BranchMorToLitterF*vo['C_M_Reg'][iT,:,iEP['Branch']]+ \
        psl.bIPF_StemNonMerchLitToLitterF*vo['C_LF'][iT,:,iEP['StemNonMerch']]
            
    # Transfer biomass turnover to medium litter pool
    vo['C_Eco_Pools'][iT,:,iEP['LitterM']]=vo['C_Eco_Pools'][iT-1,:,iEP['LitterM']]+ \
        psl.bIPF_BarkLitToLitterM*vo['C_LF'][iT,:,iEP['Bark']]+ \
        psl.bIPF_BarkMorToLitterM*vo['C_M_Reg'][iT,:,iEP['Bark']]

    # Transfer biomass turnover to slow litter pool
    vo['C_Eco_Pools'][iT,:,iEP['LitterS']]=vo['C_Eco_Pools'][iT-1,:,iEP['LitterS']]+ \
        psl.bIPF_BarkLitToLitterS*vo['C_LF'][iT,:,iEP['Bark']]+ \
        psl.bIPF_BarkMorToLitterS*vo['C_M_Reg'][iT,:,iEP['Bark']]
        
    # Transfer biomass turnover to snag stemwood pool
    vo['C_Eco_Pools'][iT,:,iEP['SnagStem']]=vo['C_Eco_Pools'][iT-1,:,iEP['SnagStem']]+ \
        psl.bIPF_StemMerchMorToSnagStem*vo['C_M_Reg'][iT,:,iEP['StemMerch']]+ \
        psl.bIPF_StemNonMerchMorToSnagStem*vo['C_M_Reg'][iT,:,iEP['StemNonMerch']]
    
    # Transfer biomass turnover to snag branch pool
    vo['C_Eco_Pools'][iT,:,iEP['SnagBranch']]=vo['C_Eco_Pools'][iT-1,:,iEP['SnagBranch']]+ \
        psl.bIPF_BranchLitToSnagBranch*vo['C_LF'][iT,:,iEP['Branch']]+ \
        psl.bIPF_BranchMorToSnagBranch*vo['C_M_Reg'][iT,:,iEP['Branch']]+ \
        psl.bIPF_BarkMorToSnagBranch*vo['C_M_Reg'][iT,:,iEP['Bark']]+ \
        psl.bIPF_StemNonMerchLitToSnagBranch*vo['C_LF'][iT,:,iEP['Bark']]

    # Transfer biomass turnover to very fast soil pool
    vo['C_Eco_Pools'][iT,:,iEP['SoilVF']]=vo['C_Eco_Pools'][iT-1,:,iEP['SoilVF']]+ \
        psl.bIPF_RootFineLitToSoilVF*vo['C_LF'][iT,:,iEP['RootFine']] +\
        psl.bIPF_RootCoarseMorToSoilVF*vo['C_M_Reg'][iT,:,iEP['RootCoarse']]+ \
        psl.bIPF_RootFineMorToSoilVF*vo['C_M_Reg'][iT,:,iEP['RootFine']]

    # Transfer biomass turnover to fast soil pool
    vo['C_Eco_Pools'][iT,:,iEP['SoilF']]=vo['C_Eco_Pools'][iT-1,:,iEP['SoilF']]+ \
        psl.bIPF_RootCoarseLitToSoilF*vo['C_LF'][iT,:,iEP['RootCoarse']]+ \
        psl.bIPF_RootCoarseMorToSoilF*vo['C_M_Reg'][iT,:,iEP['RootCoarse']]
        
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
    R_LitterVF=psl.bDec_LitterVF_Rten*vo['C_Eco_Pools'][iT,:,iEP['LitterVF']]*psl.bDec_LitterVF_Qten**fT
    R_LitterF=psl.bDec_LitterF_Rten*vo['C_Eco_Pools'][iT,:,iEP['LitterF']]*psl.bDec_LitterF_Qten**fT
    R_LitterM=psl.bDec_LitterM_Rten*vo['C_Eco_Pools'][iT,:,iEP['LitterM']]*psl.bDec_LitterM_Qten**fT
    R_LitterS=psl.bDec_LitterS_Rten*vo['C_Eco_Pools'][iT,:,iEP['LitterS']]*psl.bDec_LitterS_Qten**fT
    R_SnagStem=psl.bDec_SnagStem_Rten*vo['C_Eco_Pools'][iT,:,iEP['SnagStem']]*psl.bDec_SnagStem_Qten**fT
    R_SnagBranch=psl.bDec_SnagBranch_Rten*vo['C_Eco_Pools'][iT,:,iEP['SnagBranch']]*psl.bDec_SnagBranch_Qten**fT    
    R_SoilVF=psl.bDec_SoilVF_Rten*vo['C_Eco_Pools'][iT,:,iEP['SoilVF']]*psl.bDec_SoilVF_Qten**fT
    R_SoilF=psl.bDec_SoilF_Rten*vo['C_Eco_Pools'][iT,:,iEP['SoilF']]*psl.bDec_SoilF_Qten**fT
    R_SoilS=psl.bDec_SoilS_Rten*vo['C_Eco_Pools'][iT,:,iEP['SoilS']]*psl.bDec_SoilS_Qten**fT
    
    # These need to be flattened *** Is there a workaround for this? ***
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
    vo['C_Eco_Pools'][iT,:,iEP['LitterVF']]=vo['C_Eco_Pools'][iT,:,iEP['LitterVF']]-R_LitterVF
    vo['C_Eco_Pools'][iT,:,iEP['LitterF']]=vo['C_Eco_Pools'][iT,:,iEP['LitterF']]-R_LitterF
    vo['C_Eco_Pools'][iT,:,iEP['LitterM']]=vo['C_Eco_Pools'][iT,:,iEP['LitterM']]-R_LitterM
    vo['C_Eco_Pools'][iT,:,iEP['LitterS']]=vo['C_Eco_Pools'][iT,:,iEP['LitterS']]-R_LitterS
    vo['C_Eco_Pools'][iT,:,iEP['SnagStem']]=vo['C_Eco_Pools'][iT,:,iEP['SnagStem']]-R_SnagStem
    vo['C_Eco_Pools'][iT,:,iEP['SnagBranch']]=vo['C_Eco_Pools'][iT,:,iEP['SnagBranch']]-R_SnagBranch
    vo['C_Eco_Pools'][iT,:,iEP['SoilVF']]=vo['C_Eco_Pools'][iT,:,iEP['SoilVF']]-R_SoilVF
    vo['C_Eco_Pools'][iT,:,iEP['SoilF']]=vo['C_Eco_Pools'][iT,:,iEP['SoilF']]-R_SoilF
    vo['C_Eco_Pools'][iT,:,iEP['SoilS']]=vo['C_Eco_Pools'][iT,:,iEP['SoilS']]-R_SoilS     
                
    # Re-define decayed fast litter
    vo['C_Eco_Pools'][iT,:,iEP['LitterF']]=vo['C_Eco_Pools'][iT,:,iEP['LitterF']]+psl.bIPF_SnagStemToLitterF*R_SnagStem
    
    # Re-define decayed medium litter
    vo['C_Eco_Pools'][iT,:,iEP['LitterM']]=vo['C_Eco_Pools'][iT,:,iEP['LitterM']]+psl.bIPF_SnagBranchToLitterM*R_SnagBranch
    
    # Re-define decayed slow litter
    vo['C_Eco_Pools'][iT,:,iEP['LitterS']]=vo['C_Eco_Pools'][iT,:,iEP['LitterS']]+psl.bIPF_LitterVFToLitterS*R_LitterVF
    vo['C_Eco_Pools'][iT,:,iEP['LitterS']]=vo['C_Eco_Pools'][iT,:,iEP['LitterS']]+psl.bIPF_LitterFToLitterS*R_LitterF
    vo['C_Eco_Pools'][iT,:,iEP['LitterS']]=vo['C_Eco_Pools'][iT,:,iEP['LitterS']]+psl.bIPF_LitterMToLitterS*R_LitterM
    
    # Re-define decayed slow soil
    vo['C_Eco_Pools'][iT,:,iEP['SoilS']]=vo['C_Eco_Pools'][iT,:,iEP['SoilS']]+psl.bIPF_SoilVFToSoilS*R_SoilVF
    vo['C_Eco_Pools'][iT,:,iEP['SoilS']]=vo['C_Eco_Pools'][iT,:,iEP['SoilS']]+psl.bIPF_SoilFToSoilS*R_SoilF
        
    # Heterotrophic respiration 
    vo['C_RH'][iT,:,iEP['LitterVF']]=psl.bIPF_LitterVFToCO2*R_LitterVF
    vo['C_RH'][iT,:,iEP['LitterF']]=psl.bIPF_LitterFToCO2*R_LitterF
    vo['C_RH'][iT,:,iEP['LitterM']]=psl.bIPF_LitterMToCO2*R_LitterM
    vo['C_RH'][iT,:,iEP['LitterS']]=psl.bIPF_LitterSToCO2*R_LitterS
    vo['C_RH'][iT,:,iEP['SnagStem']]=psl.bIPF_SnagStemToCO2*R_SnagStem
    vo['C_RH'][iT,:,iEP['SnagBranch']]=psl.bIPF_SnagBranchToCO2*R_SnagBranch
    vo['C_RH'][iT,:,iEP['SoilVF']]=psl.bIPF_SoilVFToCO2*R_SoilVF
    vo['C_RH'][iT,:,iEP['SoilF']]=psl.bIPF_SoilFToCO2*R_SoilF
    vo['C_RH'][iT,:,iEP['SoilS']]=psl.bIPF_SoilSToCO2*R_SoilS    
    
    #--------------------------------------------------------------------------
    # Physical transfer
    #--------------------------------------------------------------------------
    
    # Physical transfer rate
    PT_FelledStemMerch=psl.bDec_FelledStemMerch_PhysTransRate*vo['C_Eco_Pools'][iT,:,iEP['FelledStemMerch']]
    PT_FelledStemNonMerch=psl.bDec_FelledStemNonMerch_PhysTransRate*vo['C_Eco_Pools'][iT,:,iEP['FelledStemNonMerch']]
    PT_FelledBranch=psl.bDec_FelledBranch_PhysTransRate*vo['C_Eco_Pools'][iT,:,iEP['FelledBranch']]
    PT_FelledBark=psl.bDec_FelledBark_PhysTransRate*vo['C_Eco_Pools'][iT,:,iEP['FelledBark']]
    PT_FelledSnagStem=psl.bDec_FelledSnagStem_PhysTransRate*vo['C_Eco_Pools'][iT,:,iEP['FelledSnagStem']]
    PT_FelledSnagBranch=psl.bDec_FelledSnagBranch_PhysTransRate*vo['C_Eco_Pools'][iT,:,iEP['FelledSnagBranch']]
    
    PT_LitterS=psl.bDec_LitterS_PhysTransRate*vo['C_Eco_Pools'][iT,:,iEP['LitterS']]
    PT_SnagStem=psl.bDec_SnagStem_PhysTransRate*vo['C_Eco_Pools'][iT,:,iEP['SnagStem']]
    PT_SnagBranch=psl.bDec_SnagBranch_PhysTransRate*vo['C_Eco_Pools'][iT,:,iEP['SnagBranch']] 
    
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
        
    return vo

'''============================================================================
STAND DISTURBANCES
============================================================================'''

def Disturbances(iT,vi,vo,psl,meta,iEP):
        
    for iS in range(meta['N Stand']):
        
        # Only continue if there are disturbances
        iDist=np.where(vi['DH'][iS]['Year']==vi['tv'][iT])[0]
        
        if len(iDist)==0:
            continue
        
        for i in range(iDist.size):
            
            # Disturbance type
            ID_Type=vi['DH'][iS]['ID_Type'][iDist[i]]
            
            # Index to default parameters for this disturbance type
            iDefPar=np.where(psl.bDist_ID==ID_Type)[0]

            if vi['DH'][iS]['Biomass_Affected_Pct'][iDist[i]]<=0:
                
                #--------------------------------------------------------------
                # No custom harvesting disturbance parameters supplied, use default
                # assumptions
                #--------------------------------------------------------------
                
                # Use severity from disturbance history
                Severity=vi['DH'][iS]['Severity'][iDist[i]].astype(float)/100
                
                # Force Aerial Fertilization and Planting to be zero severity
                if (ID_Type==meta['LUT Dist']['Fertilization Aerial']) | (ID_Type==meta['LUT Dist']['Planting']):
                    Severity=0
                
                # Apply severity to biomass and snags
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
                
                #--------------------------------------------------------------
                # Use custom harvesting disturbance parameter table
                #--------------------------------------------------------------
                
                # These come from a custom user input spreadsheet, but they have been added to DH structure in pre-processing
                Biomass_Affected_Pct=vi['DH'][iS]['Biomass_Affected_Pct'][iDist[i]].astype(float)/100                
                Snags_Affected_Pct=vi['DH'][iS]['Snags_Affected_Pct'][iDist[i]].astype(float)/100
                
                # Proportions removed, burned, and left on site                        
                Biomass_Merch_Removed_Pct=vi['DH'][iS]['Biomass_Merch_Removed_Pct'][iDist[i]].astype(float)/100
                Biomass_Merch_Burned_Pct=vi['DH'][iS]['Biomass_Merch_Burned_Pct'][iDist[i]].astype(float)/100          
                Biomass_Merch_LeftOnSite_Pct=vi['DH'][iS]['Biomass_Merch_LeftOnSite_Pct'][iDist[i]].astype(float)/100
                Biomass_Merch_ToSnag_Pct=0                
            
                Biomass_NonMerch_Removed_Pct=vi['DH'][iS]['Biomass_NonMerch_Removed_Pct'][iDist[i]].astype(float)/100
                Biomass_NonMerch_Burned_Pct=vi['DH'][iS]['Biomass_NonMerch_Burned_Pct'][iDist[i]].astype(float)/100          
                Biomass_NonMerch_LeftOnSite_Pct=vi['DH'][iS]['Biomass_NonMerch_LeftOnSite_Pct'][iDist[i]].astype(float)/100
                Biomass_NonMerch_ToSnag_Pct=0
                
                Snags_Removed_Pct=vi['DH'][iS]['Snags_Removed_Pct'][iDist[i]].astype(float)/100
                Snags_Burned_Pct=vi['DH'][iS]['Snags_Burned_Pct'][iDist[i]].astype(float)/100
                Snags_LeftOnSite_Pct=vi['DH'][iS]['Snags_LeftOnSite_Pct'][iDist[i]].astype(float)/100
                        
            #------------------------------------------------------------------
            # Total amount of each pool that is affected
            #------------------------------------------------------------------
            
            # Define biomass that is affected
            Affected_StemMerch=Biomass_Affected_Pct*vo['C_Eco_Pools'][iT,iS,iEP['StemMerch']]
            Affected_StemNonMerch=Biomass_Affected_Pct*vo['C_Eco_Pools'][iT,iS,iEP['StemNonMerch']]
            Affected_Foliage=Biomass_Affected_Pct*vo['C_Eco_Pools'][iT,iS,iEP['Foliage']]
            Affected_Branch=Biomass_Affected_Pct*vo['C_Eco_Pools'][iT,iS,iEP['Branch']]
            Affected_Bark=Biomass_Affected_Pct*vo['C_Eco_Pools'][iT,iS,iEP['Bark']]
            Affected_RootCoarse=Biomass_Affected_Pct*vo['C_Eco_Pools'][iT,iS,iEP['RootCoarse']]
            Affected_RootFine=Biomass_Affected_Pct*vo['C_Eco_Pools'][iT,iS,iEP['RootFine']]
            
            # Volume
            Affected_VolumeStemMerch=Biomass_Affected_Pct*vo['V_StemMerch'][iT,iS]
                        
            # Partition bark into merch and non-merch components
            Affected_Bark_Merch=psl.bASL_MerchBarkFrac*Affected_Bark
            Affected_Bark_NonMerch=(1-psl.bASL_MerchBarkFrac)*Affected_Bark
            
            # Sum up total affected non-merchantable biomass
            Affected_TotNonMerch=Affected_StemNonMerch+Affected_Branch+Affected_Bark_NonMerch
            
            # Define fraction of snags that are affected            
            Affected_SnagStem=Snags_Affected_Pct*vo['C_Eco_Pools'][iT,iS,iEP['SnagStem']]
            Affected_SnagBranch=Snags_Affected_Pct*vo['C_Eco_Pools'][iT,iS,iEP['SnagBranch']]
            
            # All
            Affected_All=Affected_StemMerch+Affected_StemNonMerch+Affected_Foliage+Affected_Branch+Affected_Bark+Affected_RootCoarse+Affected_RootFine
            
            #------------------------------------------------------------------
            # Calculate mortality
            #------------------------------------------------------------------
            
            DistType_i=vi['DH'][iS]['ID_Type'][iDist[i]]
            
            if DistType_i==meta['LUT Dist']['Wildfire']:
                vo['C_M_Inv_Fir'][iT,iS]=Affected_All
            elif DistType_i==meta['LUT Dist']['Beetles']:
                vo['C_M_Inv_Ins'][iT,iS]=Affected_All
            elif DistType_i==meta['LUT Dist']['Defoliators']:
                vo['C_M_Inv_Ins'][iT,iS]=Affected_All 
            elif (DistType_i==meta['LUT Dist']['Harvest']) | (DistType_i==meta['LUT Dist']['Harvest Custom']):
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
            Biomass_Merch_Removed=Biomass_Merch_Removed_Pct*(Affected_StemMerch+Affected_Bark_Merch)
            vo['C_Eco_Pools'][iT,iS,iEP['RemovedMerch']]=vo['C_Eco_Pools'][iT,iS,iEP['RemovedMerch']]+Biomass_Merch_Removed
            
            # NonMerch biomass to mill - of the total amount of biomass affected, 
            # what fraction of non-merch biomass was sent to the mill?
            # - NonMerch = NonMerchStem + Foliage + Branch + Bark
            Biomass_NonMerch_Removed=Biomass_NonMerch_Removed_Pct*Affected_TotNonMerch
            vo['C_Eco_Pools'][iT,iS,iEP['RemovedNonMerch']]=vo['C_Eco_Pools'][iT,iS,iEP['RemovedNonMerch']]+Biomass_NonMerch_Removed
                                        
            # Snag stemwood to mill           
            SnagStem_Removed=Snags_Removed_Pct*Affected_SnagStem
            vo['C_Eco_Pools'][iT,iS,iEP['RemovedSnagStem']]=vo['C_Eco_Pools'][iT,iS,iEP['RemovedSnagStem']]+SnagStem_Removed
            
            # Snag branches to mill
            SnagBranch_Removed=Snags_Removed_Pct*Affected_SnagBranch
            vo['C_Eco_Pools'][iT,iS,iEP['RemovedSnagStem']]=vo['C_Eco_Pools'][iT,iS,iEP['RemovedSnagStem']]+SnagBranch_Removed
            
            #------------------------------------------------------------------
            # Carbon burned
            #------------------------------------------------------------------
            
            # Profile of emitted compounds (N2O emissions calculated upon export)
            tf_co2=psl.bCombFrac_CO2
            tf_ch4=psl.bCombFrac_CH4
            tf_co=psl.bCombFrac_CO
            
            # Merch biomass that is burned              
            Biomass_Merch_Burned_co2=tf_co2*Biomass_Merch_Burned_Pct*(Affected_StemMerch+Affected_Bark_Merch)
            vo['C_Eco_Pools'][iT,iS,iEP['ECO2asC']]=vo['C_Eco_Pools'][iT,iS,iEP['ECO2asC']]+Biomass_Merch_Burned_co2
            
            Biomass_Merch_Burned_ch4=tf_ch4*Biomass_Merch_Burned_Pct*(Affected_StemMerch+Affected_Bark_Merch)
            vo['C_Eco_Pools'][iT,iS,iEP['ECH4asC']]=vo['C_Eco_Pools'][iT,iS,iEP['ECH4asC']]+Biomass_Merch_Burned_ch4
            
            Biomass_Merch_Burned_co=tf_co*Biomass_Merch_Burned_Pct*(Affected_StemMerch+Affected_Bark_Merch)
            vo['C_Eco_Pools'][iT,iS,iEP['ECOasC']]=vo['C_Eco_Pools'][iT,iS,iEP['ECOasC']]+Biomass_Merch_Burned_co

            # NonMerch biomass that is burned            
            Biomass_NonMerch_Burned_co2=tf_co2*Biomass_NonMerch_Burned_Pct*Affected_TotNonMerch
            vo['C_Eco_Pools'][iT,iS,iEP['ECO2asC']]=vo['C_Eco_Pools'][iT,iS,iEP['ECO2asC']]+Biomass_NonMerch_Burned_co2
            
            Biomass_NonMerch_Burned_ch4=tf_ch4*Biomass_NonMerch_Burned_Pct*Affected_TotNonMerch
            vo['C_Eco_Pools'][iT,iS,iEP['ECH4asC']]=vo['C_Eco_Pools'][iT,iS,iEP['ECH4asC']]+Biomass_NonMerch_Burned_ch4
            
            Biomass_NonMerch_Burned_co=tf_co*Biomass_NonMerch_Burned_Pct*Affected_TotNonMerch
            vo['C_Eco_Pools'][iT,iS,iEP['ECOasC']]=vo['C_Eco_Pools'][iT,iS,iEP['ECOasC']]+Biomass_NonMerch_Burned_co
              
            # Snag stems that burn
            SnagStem_Burned_co2=tf_co2*Snags_Burned_Pct*Affected_SnagStem
            vo['C_Eco_Pools'][iT,iS,iEP['ECO2asC']]=vo['C_Eco_Pools'][iT,iS,iEP['ECO2asC']]+SnagStem_Burned_co2
            
            SnagStem_Burned_ch4=tf_ch4*Snags_Burned_Pct*Affected_SnagStem
            vo['C_Eco_Pools'][iT,iS,iEP['ECH4asC']]=vo['C_Eco_Pools'][iT,iS,iEP['ECH4asC']]+SnagStem_Burned_ch4
            
            SnagStem_Burned_co=tf_co*Snags_Burned_Pct*Affected_SnagStem
            vo['C_Eco_Pools'][iT,iS,iEP['ECOasC']]=vo['C_Eco_Pools'][iT,iS,iEP['ECOasC']]+SnagStem_Burned_co
            
            # Snag branches that burn
            SnagBranch_Burned_co2=tf_co2*Snags_Burned_Pct*Affected_SnagBranch
            vo['C_Eco_Pools'][iT,iS,iEP['ECO2asC']]=vo['C_Eco_Pools'][iT,iS,iEP['ECO2asC']]+SnagBranch_Burned_co2
            
            SnagBranch_Burned_ch4=tf_ch4*Snags_Burned_Pct*Affected_SnagBranch
            vo['C_Eco_Pools'][iT,iS,iEP['ECH4asC']]=vo['C_Eco_Pools'][iT,iS,iEP['ECH4asC']]+SnagBranch_Burned_ch4
            
            SnagBranch_Burned_co=tf_co*Snags_Burned_Pct*Affected_SnagBranch
            vo['C_Eco_Pools'][iT,iS,iEP['ECOasC']]=vo['C_Eco_Pools'][iT,iS,iEP['ECOasC']]+SnagBranch_Burned_co
            
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
            vo['C_Eco_Pools'][iT,iS,iEP['SnagStem']]=vo['C_Eco_Pools'][iT,iS,iEP['SnagStem']]+Biomass_Merch_ToSnag_Pct*(Affected_StemMerch+Affected_Bark_Merch)
            
            # Non-merch biomass that is killed by wildfire
            vo['C_Eco_Pools'][iT,iS,iEP['SnagBranch']]=vo['C_Eco_Pools'][iT,iS,iEP['SnagBranch']]+Biomass_NonMerch_ToSnag_Pct*(Affected_StemNonMerch+Affected_Bark_NonMerch)
            
            #------------------------------------------------------------------
            # Carbon that is left on site (after felling or wind storms)
            #------------------------------------------------------------------
            
            # Stem, branch and bark carbon transferred to felled pools (piles)
            vo['C_Eco_Pools'][iT,iS,iEP['FelledStemMerch']]=vo['C_Eco_Pools'][iT,iS,iEP['FelledStemMerch']]+Biomass_Merch_LeftOnSite_Pct*Affected_StemMerch
            vo['C_Eco_Pools'][iT,iS,iEP['FelledStemNonMerch']]=vo['C_Eco_Pools'][iT,iS,iEP['FelledStemNonMerch']]+Biomass_NonMerch_LeftOnSite_Pct*Affected_StemNonMerch            
            vo['C_Eco_Pools'][iT,iS,iEP['FelledBranch']]=vo['C_Eco_Pools'][iT,iS,iEP['FelledBranch']]+Biomass_Merch_LeftOnSite_Pct*Affected_Branch
            vo['C_Eco_Pools'][iT,iS,iEP['FelledBark']]=vo['C_Eco_Pools'][iT,iS,iEP['FelledBark']]+Biomass_Merch_LeftOnSite_Pct*Affected_Bark
                        
            # Foliage transferred directly to very fast litter
            vo['C_Eco_Pools'][iT,iS,iEP['LitterVF']]=vo['C_Eco_Pools'][iT,iS,iEP['LitterVF']]+Affected_Foliage
            
            # Roots transferred directly to DOM       
            vo['C_Eco_Pools'][iT,iS,iEP['LitterF']]=vo['C_Eco_Pools'][iT,iS,iEP['LitterF']]+0.5*Affected_RootCoarse
            vo['C_Eco_Pools'][iT,iS,iEP['SoilF']]=vo['C_Eco_Pools'][iT,iS,iEP['SoilF']]+0.5*Affected_RootCoarse            
            vo['C_Eco_Pools'][iT,iS,iEP['LitterVF']]=vo['C_Eco_Pools'][iT,iS,iEP['LitterVF']]+0.5*Affected_RootFine
            vo['C_Eco_Pools'][iT,iS,iEP['SoilVF']]=vo['C_Eco_Pools'][iT,iS,iEP['SoilVF']]+0.5*Affected_RootFine
            
            # Snags transferred to felled pools (piles)
            vo['C_Eco_Pools'][iT,iS,iEP['FelledSnagStem']]=vo['C_Eco_Pools'][iT,iS,iEP['FelledSnagStem']]+Snags_LeftOnSite_Pct*Affected_SnagStem
            vo['C_Eco_Pools'][iT,iS,iEP['FelledSnagBranch']]=vo['C_Eco_Pools'][iT,iS,iEP['FelledSnagBranch']]+Snags_LeftOnSite_Pct*Affected_SnagBranch
            
            #------------------------------------------------------------------
            # Update stand age
            # *** Work on making this more flexible to users without hardwiring
            # how it is run here. ***
            #------------------------------------------------------------------
            
            if meta['Biomass Module']=='TIPSY':
                
                if (ID_Type==meta['LUT Dist']['Planting']) | (ID_Type==meta['LUT Dist']['Direct Seeding']):
                    
                    # If planting or direct seeding, force age to be 1
                    vo['A'][iT,iS]=1
                    
                else:
                
                    #age_mod_meth='CBM'
                    age_mod_meth='FunctionOfSeverity'
                
                    if age_mod_meth=='CBM':
                        
                        # I think this is how CBM works - something like this
                        if Biomass_Affected_Pct>0.9:
                            vo['A'][iT,iS]=1
                        
                    elif age_mod_meth=='FunctionOfSeverity':
                        vo['A'][iT,iS]=vo['A'][iT,iS]*(1-Biomass_Affected_Pct)
            
            #------------------------------------------------------------------
            # Update net growth (in response to lethal events)
            #------------------------------------------------------------------
            
            if meta['Biomass Module']=='TIPSY':
                
                G_pre=vi['GCA'][:,iS,:]
                G_post=(1-Biomass_Affected_Pct)*vi['GCA'][:,iS,:]
                dG=G_pre-G_post
                
                RecoveryHalfLife=15               
                TSD=np.tile(np.arange(0,301,1)-vo['A'][iT,iS],G_pre.shape[1])
                fTSD=1/(1+np.exp(-0.5*(TSD-RecoveryHalfLife))) 
                fTSD[0:vo['A'][iT,iS],:]=0
                G_Recovery=fTSD*dG
                
                vi['GCA'][:,iS,:]=G_post+G_Recovery
            
            #------------------------------------------------------------------
            # Transition to new growth curve
            #------------------------------------------------------------------
            
            # Only applies to TIPSY
            if meta['Biomass Module']=='TIPSY':            
                
                # Does not apply to fertilization
                if (vi['DH'][iS]['ID_Type'][iDist[i]]!=meta['LUT Dist']['Fertilization Aerial']):
                    
                    if vi['DH'][iS]['ID_GrowthCurveM'][iDist[i]]==1:                        
                        vi['GCA'][:,iS,:]=vi['GC1'][:,iS,:]                        
                    elif vi['DH'][iS]['ID_GrowthCurveM'][iDist[i]]==2:
                        vi['GCA'][:,iS,:]=vi['GC2'][:,iS,:] 
                    elif vi['DH'][iS]['ID_GrowthCurveM'][iDist[i]]==3:
                        vi['GCA'][:,iS,:]=vi['GC3'][:,iS,:]                        
                    elif vi['DH'][iS]['ID_GrowthCurveM'][iDist[i]]==4:
                        vi['GCA'][:,iS,:]=vi['GC4'][:,iS,:]                        
    
            #------------------------------------------------------------------
            # Apply growth factors (in response to non-lethal events)
            #------------------------------------------------------------------
            
            # Only applies to TIPSY
            if meta['Biomass Module']=='TIPSY':   
                vi['GCA'][:,iS,:]=psl.bDist_GrowthFactor[iDefPar]*vi['GCA'][:,iS,:]                
    
            #------------------------------------------------------------------
            # Aerial fertilization - growth response and
            #------------------------------------------------------------------
            
            if ID_Type==meta['LUT Dist']['Fertilization Aerial']:
                
                #--------------------------------------------------------------
                # Growth response
                # Affecting biomass pools equally
                #--------------------------------------------------------------
                
                if meta['Fertilization Source']=='CBRunner':
                
                    ResponseRatio=1.04
                
                    ResponsePeriod=10
                
                    iResponse=np.arange(vo['A'][iT,iS],vo['A'][iT,iS]+ResponsePeriod)
                
                    iP=meta['iEP']['BiomassTotal']
                
                    vi['GCA'][iResponse,iS,iP]=ResponseRatio*vi['GCA'][iResponse,iS,iP]
                
                #--------------------------------------------------------------
                # Operational emissions
                #--------------------------------------------------------------
                
                # Urea dose (kgUrea/ha)
                DoseUrea=psl.bDoseUrea
                
                # Nitrogen dose (kgN/ha)
                DoseN=psl.bUreaRatioThatIsNitrogen*DoseUrea

                # Emissions from manufacturing, transport prior to application 
                # and transport during application (tCO2e/ha). This is derived from 
                # Sonne (2006) as 30% of the total emissions of 1.9 tCO2e/ha. 
                # It is higher than the work of Mel Scott (2010)
                E_Op=psl.bEmissionsFromManufacturingAndTransport

                # Volatilization: CO2 emissions following application, Tier 1 
                # approach, IPCC 2006, 11.4.1 (tCO2e/ha)
                E_Vol=psl.bRatio_C_to_Urea*(DoseUrea/1000)*(1/psl.bRatio_C_to_CO2)

                # Denitrification: N2O emissions following application, Tier 1 
                # approach, IPCC 2006, 11.4.1 (tCO2e/ha)
                E_Denit=psl.bFractionOfNAppliedThatIsEmittedAsN2OAsN*(DoseN/1000)* \
                    psl.bRatioN2OAsN_to_N2O*psl.bGWP_N2O_AR4
                
                # Total emissions (tCO2e/ha)
                E_Tot=E_Op+E_Vol+E_Denit
                
                # Total emissions of carbon (MgC/ha) emitted as CO2e, converted
                # to carbon to be consistent with the rest of the variables in 
                # the vo object.
                vo['C_E_OperationsAsCO2'][iT,iS]=psl.bRatio_C_to_CO2*E_Tot
    
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
        
    return vo

'''============================================================================
PRODUCT SECTOR
============================================================================'''

def ProductSector(iT,vi,vo,psl,meta):
    
    #             0     1     2      3      4      5         6       7       8      9          10             11         12          13                       14                          15                        16                           17      18     19
    # PoolNames=['SFH','MFH','Comm','Furn','Ship','Repairs','Other','Paper','Fuel','Firewood','EffluentPulp','DumpWood','DumpPaper','LandfillWoodDegradable','LandfillWoodNonDegradable','LandfillPaperDegradable','LandfillPaperNonDegradable','E_CO2','E_CH4','Cants']
            
    #--------------------------------------------------------------------------
    # Ecosystems --> Mills
    # Note: Mills act as transient reservoirs that are not tracked over time
    #--------------------------------------------------------------------------
    
    # Assume default parameters are applied, then update below if user has specified input variables
    Mill_Chips=psl.HWP_RemovedMerchToChipping*vo['C_RemovedMerch'][iT,:] + \
                 psl.HWP_RemovedNonMerchToChipping*vo['C_RemovedNonMerch'][iT,:] + \
                 psl.HWP_RemovedSnagStemToChipping*vo['C_RemovedSnagStem'][iT,:]
    
    Mill_Chips_Pulp=psl.HWP_ChipsToPulp*Mill_Chips
    Mill_Chips_Fuel=psl.HWP_ChipsToFuel*Mill_Chips
    
    Mill_Lumber=psl.HWP_RemovedMerchToLumber*vo['C_RemovedMerch'][iT,:] + \
                 psl.HWP_RemovedNonMerchToLumber*vo['C_RemovedNonMerch'][iT,:] + \
                 psl.HWP_RemovedSnagStemToLumber*vo['C_RemovedSnagStem'][iT,:]
    Mill_Plywood=psl.HWP_RemovedMerchToPlywood*vo['C_RemovedMerch'][iT,:] + \
                 psl.HWP_RemovedNonMerchToPlywood*vo['C_RemovedNonMerch'][iT,:] + \
                 psl.HWP_RemovedSnagStemToPlywood*vo['C_RemovedSnagStem'][iT,:]
    Mill_OSB=psl.HWP_RemovedMerchToOSB*vo['C_RemovedMerch'][iT,:] + \
                 psl.HWP_RemovedNonMerchToOSB*vo['C_RemovedNonMerch'][iT,:] + \
                 psl.HWP_RemovedSnagStemToOSB*vo['C_RemovedSnagStem'][iT,:]
    Mill_MDF=psl.HWP_RemovedMerchToMDF*vo['C_RemovedMerch'][iT,:] + \
                 psl.HWP_RemovedNonMerchToMDF*vo['C_RemovedNonMerch'][iT,:] + \
                 psl.HWP_RemovedSnagStemToMDF*vo['C_RemovedSnagStem'][iT,:]
    Mill_Firewood=psl.HWP_RemovedMerchToFirewood*vo['C_RemovedMerch'][iT,:] + \
                 psl.HWP_RemovedNonMerchToFirewood*vo['C_RemovedNonMerch'][iT,:] + \
                 psl.HWP_RemovedSnagStemToFirewood*vo['C_RemovedSnagStem'][iT,:]
    Mill_Cants=0*vo['C_RemovedMerch'][iT,:] # Not considered among default in-use products
                 
    # Apply user-specified transfer percents
    for iS in range(meta['N Stand']):
        
        # Find custom disturbance inputs
        it=np.where((vi['DH'][iS]['Year']==vi['tv'][iT]) & (vi['DH'][iS]['Biomass_Affected_Pct']>0) |
                (vi['DH'][iS]['Year']==vi['tv'][iT]) & (vi['DH'][iS]['Snags_Affected_Pct']>0))[0]
        
        # Only continue if there are custom disturbance inputs
        if it.shape[0]==0: 
            continue
        
        # Run error if user inputs multiple custom disturbances in single year
        if it.shape[0]>1:
            print('Warning: Error - multiple custom disturbances in single year.')
            continue
                
        Mill_Chips_Pulp[iS]=vi['DH'][iS]['RemovedMerchToPulp_Pct'][it]/100*vo['C_RemovedMerch'][iT,iS] + \
                 vi['DH'][iS]['RemovedNonMerchToPulp_Pct'][it]/100*vo['C_RemovedNonMerch'][iT,iS] + \
                 vi['DH'][iS]['RemovedSnagStemToPulp_Pct'][it]/100*vo['C_RemovedSnagStem'][iT,iS]
        Mill_Chips_Fuel[iS]=vi['DH'][iS]['RemovedMerchToFuel_Pct'][it]/100*vo['C_RemovedMerch'][iT,iS] + \
                 vi['DH'][iS]['RemovedNonMerchToFuel_Pct'][it]/100*vo['C_RemovedNonMerch'][iT,iS] + \
                 vi['DH'][iS]['RemovedSnagStemToFuel_Pct'][it]/100*vo['C_RemovedSnagStem'][iT,iS]
        Mill_Lumber[iS]=vi['DH'][iS]['RemovedMerchToLumber_Pct'][it]/100*vo['C_RemovedMerch'][iT,iS] + \
                 vi['DH'][iS]['RemovedNonMerchToLumber_Pct'][it]/100*vo['C_RemovedNonMerch'][iT,iS] + \
                 vi['DH'][iS]['RemovedSnagStemToLumber_Pct'][it]/100*vo['C_RemovedSnagStem'][iT,iS]
        Mill_Plywood[iS]=vi['DH'][iS]['RemovedMerchToPlywood_Pct'][it]/100*vo['C_RemovedMerch'][iT,iS] + \
                 vi['DH'][iS]['RemovedNonMerchToPlywood_Pct'][it]/100*vo['C_RemovedNonMerch'][iT,iS] + \
                 vi['DH'][iS]['RemovedSnagStemToPlywood_Pct'][it]/100*vo['C_RemovedSnagStem'][iT,iS]
        Mill_OSB[iS]=vi['DH'][iS]['RemovedMerchToOSB_Pct'][it]/100*vo['C_RemovedMerch'][iT,iS] + \
                 vi['DH'][iS]['RemovedNonMerchToOSB_Pct'][it]/100*vo['C_RemovedNonMerch'][iT,iS] + \
                 vi['DH'][iS]['RemovedSnagStemToOSB_Pct'][it]/100*vo['C_RemovedSnagStem'][iT,iS]
        Mill_MDF[iS]=vi['DH'][iS]['RemovedMerchToMDF_Pct'][it]/100*vo['C_RemovedMerch'][iT,iS] + \
                 vi['DH'][iS]['RemovedNonMerchToMDF_Pct'][it]/100*vo['C_RemovedNonMerch'][iT,iS] + \
                 vi['DH'][iS]['RemovedSnagStemToMDF_Pct'][it]/100*vo['C_RemovedSnagStem'][iT,iS]
        Mill_Firewood[iS]=vi['DH'][iS]['RemovedMerchToFirewood_Pct'][it]/100*vo['C_RemovedMerch'][iT,iS] + \
                 vi['DH'][iS]['RemovedNonMerchToFirewood_Pct'][it]/100*vo['C_RemovedNonMerch'][iT,iS] + \
                 vi['DH'][iS]['RemovedSnagStemToFirewood_Pct'][it]/100*vo['C_RemovedSnagStem'][iT,iS]            
        Mill_Cants[iS]=vi['DH'][iS]['RemovedMerchToCants_Pct'][it]/100*vo['C_RemovedMerch'][iT,iS] + \
                 vi['DH'][iS]['RemovedNonMerchToCants_Pct'][it]/100*vo['C_RemovedNonMerch'][iT,iS] + \
                 vi['DH'][iS]['RemovedSnagStemToCants_Pct'][it]/100*vo['C_RemovedSnagStem'][iT,iS]
     
    #--------------------------------------------------------------------------
    # Mills --> In-use pools or other mills
    #--------------------------------------------------------------------------
            
    # Transfer mill fibre to single-family homes
    ip=0
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT-1,:,ip] + \
        psl.HWP_LumberToSFH*Mill_Lumber + \
        psl.HWP_PlywoodToSFH*Mill_Plywood + \
        psl.HWP_OSBToSFH*Mill_OSB + \
        psl.HWP_MDFToSFH*Mill_MDF
    
    # Transfer mill fibre to multi-family homes
    ip=1
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT-1,:,ip] + \
        psl.HWP_LumberToMFH*Mill_Lumber + \
        psl.HWP_PlywoodToMFH*Mill_Plywood + \
        psl.HWP_OSBToMFH*Mill_OSB + \
        psl.HWP_MDFToMFH*Mill_MDF
    
    # Transfer mill fibre to commercial
    ip=2
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT-1,:,ip] + \
        psl.HWP_LumberToCom*Mill_Lumber + \
        psl.HWP_PlywoodToCom*Mill_Plywood + \
        psl.HWP_OSBToCom*Mill_OSB + \
        psl.HWP_MDFToCom*Mill_MDF
    
    # Transfer mill fibre to furniture
    ip=3
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT-1,:,ip] + \
        psl.HWP_LumberToFurn*Mill_Lumber + \
        psl.HWP_PlywoodToFurn*Mill_Plywood + \
        psl.HWP_OSBToFurn*Mill_OSB + \
        psl.HWP_MDFToFurn*Mill_MDF
    
    # Transfer mill fibre to shipping
    ip=4
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT-1,:,ip] + \
        psl.HWP_LumberToShip*Mill_Lumber + \
        psl.HWP_PlywoodToShip*Mill_Plywood + \
        psl.HWP_OSBToShip*Mill_OSB + \
        psl.HWP_MDFToShip*Mill_MDF
    
    # Transfer mill fibre to repairs
    ip=5
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT-1,:,ip] + \
        psl.HWP_LumberToRepairs*Mill_Lumber + \
        psl.HWP_PlywoodToRepairs*Mill_Plywood + \
        psl.HWP_OSBToRepairs*Mill_OSB + \
        psl.HWP_MDFToRepairs*Mill_MDF
    
    # Transfer mill fibre to other
    ip=6
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT-1,:,ip] + \
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
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT-1,:,ip] + \
                 psl.HWP_PulpToPaper*Mill_Pulp
    
    # Transfer mill fibre to fuel 
    ip=8
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT-1,:,ip] + \
                 Mill_Chips_Fuel + \
                 psl.HWP_LumberToFuel*Mill_Lumber + \
                 psl.HWP_PlywoodToFuel*Mill_Plywood + \
                 psl.HWP_OSBToFuel*Mill_OSB + \
                 psl.HWP_MDFToFuel*Mill_MDF
    
    # Transfer firewood to firewood pool
    ip=9
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT-1,:,ip] + Mill_Firewood
    
    # Transfer mill fibre to cants
    ip=19
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT-1,:,ip] + Mill_Cants
    
    # Transfer pulp mill carbon to pulp-mill effluent
    ip=10
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT-1,:,ip] + \
                 psl.HWP_PulpToEffluent*Mill_Pulp   
    
    #--------------------------------------------------------------------------
    # Update dump and landfill reservoirs
    #--------------------------------------------------------------------------

    ip=11
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT-1,:,ip]
    ip=12
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT-1,:,ip]
    ip=13
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT-1,:,ip]
    ip=14
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT-1,:,ip]    
    ip=15
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT-1,:,ip]
    ip=16
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT-1,:,ip]    
        
    #--------------------------------------------------------------------------
    # Single-family homes --> dump and landfill
    #--------------------------------------------------------------------------
        
    # Turnover
    ip=0
    Cretired=psl.HWP_SFH_tr*vo['C_Pro_Pools'][iT-1,:,ip]
    
    # Remove carbon
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] - Cretired
    
    # Transfer carbon to dump wood
    ip=11
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + psl.HWP_SFHToDumpWood*Cretired
    
    # Transfer carbon to landfill (degradable)
    ip=13
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + psl.HWP_SFHToLandfillWood*psl.HWP_ToLandfillWoodDegradableFrac*Cretired
    
    # Transfer carbon to landfill (non-degradable)
    ip=14
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + psl.HWP_SFHToLandfillWood*(1-psl.HWP_ToLandfillWoodDegradableFrac)*Cretired
    
    #--------------------------------------------------------------------------
    # Multi-family homes --> dump and landfill
    #--------------------------------------------------------------------------
    
    # Turnover
    ip=1
    Cretired=psl.HWP_MFH_tr*vo['C_Pro_Pools'][iT,:,ip]
    
    # Remove carbon
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] - Cretired
    
    # Transfer carbon to dump wood
    ip=11
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + psl.HWP_MFHToDumpWood*Cretired
    
    # Transfer carbon to landfill (degradable)
    ip=13
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + psl.HWP_MFHToLandfillWood*psl.HWP_ToLandfillWoodDegradableFrac*Cretired
    
    # Transfer carbon to landfill (non-degradable)
    ip=14
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + psl.HWP_MFHToLandfillWood*(1-psl.HWP_ToLandfillWoodDegradableFrac)*Cretired
    
    #--------------------------------------------------------------------------
    # Commercial building --> dump and landfill
    #--------------------------------------------------------------------------
    
    # Turnover
    ip=2
    Cretired=psl.HWP_Comm_tr*vo['C_Pro_Pools'][iT-1,:,ip]
    
    # Remove carbon
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] - Cretired
    
    # Transfer carbon to dump wood
    ip=11
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + psl.HWP_CommToDumpWood*Cretired
    
    # Transfer carbon to landfill (degradable)
    ip=13
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + psl.HWP_CommToLandfillWood*psl.HWP_ToLandfillWoodDegradableFrac*Cretired
    
    # Transfer carbon to landfill (non-degradable)
    ip=14
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + psl.HWP_CommToLandfillWood*(1-psl.HWP_ToLandfillWoodDegradableFrac)*Cretired
    
    #--------------------------------------------------------------------------
    # Furniture --> dump and landfill
    #--------------------------------------------------------------------------
    
    # Turnover
    ip=3
    Cretired=psl.HWP_Furn_tr*vo['C_Pro_Pools'][iT-1,:,ip]
    
    # Remove carbon
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] - Cretired
    
    # Transfer carbon to dump wood
    ip=11
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + psl.HWP_FurnToDumpWood*Cretired
    
    # Transfer carbon to landfill (degradable)
    ip=13
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + psl.HWP_FurnToLandfillWood*psl.HWP_ToLandfillWoodDegradableFrac*Cretired
    
    # Transfer carbon to landfill (non-degradable)
    ip=14
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + psl.HWP_FurnToLandfillWood*(1-psl.HWP_ToLandfillWoodDegradableFrac)*Cretired
    
    
    #--------------------------------------------------------------------------
    # Shipping --> dump and landfill
    #--------------------------------------------------------------------------
    
    # Turnover
    ip=4
    Cretired=psl.HWP_Ship_tr*vo['C_Pro_Pools'][iT-1,:,ip]
    
    # Remove carbon
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] - Cretired
    
    # Transfer carbon to dump wood
    ip=11
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + psl.HWP_ShipToDumpWood*Cretired
    
    # Transfer carbon to landfill (degradable)
    ip=13
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + psl.HWP_ShipToLandfillWood*psl.HWP_ToLandfillWoodDegradableFrac*Cretired
    
    # Transfer carbon to landfill (non-degradable)
    ip=14
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + psl.HWP_ShipToLandfillWood*(1-psl.HWP_ToLandfillWoodDegradableFrac)*Cretired
    
    #--------------------------------------------------------------------------
    # Repairs --> dump and landfill
    #--------------------------------------------------------------------------
    
    # Turnover
    ip=5
    Cretired=psl.HWP_Repairs_tr*vo['C_Pro_Pools'][iT-1,:,ip]
    
    # Remove carbon
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] - Cretired
    
    # Transfer carbon to dump wood
    ip=11
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + psl.HWP_RepairsToDumpWood*Cretired
    
    # Transfer carbon to landfill (degradble)
    ip=13
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + psl.HWP_RepairsToLandfillWood*psl.HWP_ToLandfillWoodDegradableFrac*Cretired
    
    # Transfer carbon to landfill (non-degradable)
    ip=14
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + psl.HWP_RepairsToLandfillWood*(1-psl.HWP_ToLandfillWoodDegradableFrac)*Cretired
    
    #--------------------------------------------------------------------------
    # Other --> dump and landfill
    #--------------------------------------------------------------------------
    
    # Turnover
    ip=6
    Cretired=psl.HWP_Other_tr*vo['C_Pro_Pools'][iT-1,:,ip]
    
    # Remove carbon
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] - Cretired
    
    # Transfer carbon to dump wood
    ip=11
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + psl.HWP_OtherToDumpWood*Cretired
    
    # Transfer carbon to landfill (degradble)
    ip=13
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + psl.HWP_OtherToLandfillWood*psl.HWP_ToLandfillWoodDegradableFrac*Cretired
    
    # Transfer carbon to landfill (non-degradable)
    ip=14
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + psl.HWP_OtherToLandfillWood*(1-psl.HWP_ToLandfillWoodDegradableFrac)*Cretired
    
    #--------------------------------------------------------------------------
    # Cants --> dump and landfill
    #--------------------------------------------------------------------------
    
    # Turnover
    ip=19
    Cretired=psl.HWP_Cants_tr*vo['C_Pro_Pools'][iT-1,:,ip]
    
    # Remove carbon
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] - Cretired
    
    # Transfer carbon to dump wood
    ip=11
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + psl.HWP_CantsToDumpWood*Cretired
    
    # Transfer carbon to landfill (degradble)
    ip=13
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + psl.HWP_CantsToLandfillWood*psl.HWP_ToLandfillWoodDegradableFrac*Cretired
    
    # Transfer carbon to landfill (non-degradable)
    ip=14
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + psl.HWP_CantsToLandfillWood*(1-psl.HWP_ToLandfillWoodDegradableFrac)*Cretired
    
    #--------------------------------------------------------------------------
    # Paper --> dump and landfill
    #--------------------------------------------------------------------------
    
    # Turnover (with adjustment for recycling)
    ip=7
    Cretired=(1-psl.HWP_PaperRecycleRate)*psl.HWP_Paper_tr*vo['C_Pro_Pools'][iT,:,ip]
    
    # Remove carbon
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] - Cretired
    
    # Transfer to dump
    ip=12
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + psl.HWP_PaperToDumpPaper*Cretired
    
    # Transfer to landfill (degradable)
    ip=13
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + psl.HWP_PaperToLandfillPaper*psl.HWP_ToLandfillPaperDegradableFrac*Cretired
    
    # Transfer to landfill (non-degradable)
    ip=14
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + psl.HWP_PaperToLandfillPaper*(1-psl.HWP_ToLandfillPaperDegradableFrac)*Cretired
        
    #--------------------------------------------------------------------------
    # Emissions from fuel combustion
    #--------------------------------------------------------------------------
    
    # Turnover
    ip=8
    Cretired=psl.HWP_Fuel_tr*vo['C_Pro_Pools'][iT,:,ip]
    
    # Remove carbon
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] - Cretired
    
    # Emissions of CO2 from fuel use
    ip=17
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + psl.HWP_FuelCombustionFracEmitCO2*Cretired
    
    # Emissions of CH4 from fuel use
    ip=18
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + (1-psl.HWP_FuelCombustionFracEmitCO2)*Cretired
    
    #--------------------------------------------------------------------------
    # Emissions from firewood combustion
    #--------------------------------------------------------------------------
    
    # Turnover
    ip=9
    Cretired=psl.HWP_Firewood_tr*vo['C_Pro_Pools'][iT,:,ip]
    
    # Remove carbon
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] - Cretired
    
    # Emissions of CO2
    ip=17
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + Cretired
    
    # Emissions of CH4
    #ip=18
    #vo['C_Pro_Pools[iT,:,ip]=vo['C_Pro_Pools[iT,:,ip] + (1-psl.FuelCombustionFracEmitCO2)*Cretired
       
    #--------------------------------------------------------------------------
    # Emissions from pulp effluent
    #--------------------------------------------------------------------------
         
    # Emissions from pulp effluent (CO2 from aerobic decomposition)
    ip=10
    c_emitted=psl.HWP_EffluentPulp_tr*vo['C_Pro_Pools'][iT,:,ip]
    
    # Remove emitted carbon
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip]-c_emitted
    
    # Add emitted carbon to CO2 emission "pool"
    ip=17
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip]+c_emitted    
    
    #--------------------------------------------------------------------------
    # Emissions from dump wood
    #--------------------------------------------------------------------------
    
    # Turnover
    ip=11
    c_emitted=psl.HWP_DumpWood_tr*vo['C_Pro_Pools'][iT,:,ip]
    
    # Removal
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip]-c_emitted
    
    # Add to emissions (CO2 emission from aerobic decomposition)
    ip=17
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip]+c_emitted
    
    #--------------------------------------------------------------------------
    # Emissions from dump paper
    #--------------------------------------------------------------------------
    
    # Turnover
    ip=12
    c_emitted=psl.HWP_DumpPaper_tr*vo['C_Pro_Pools'][iT,:,ip]
    
    # Removal
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip]-c_emitted
    
    # Add to emissions (CO2 emission from aerobic decomposition)
    ip=17
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip]+c_emitted
    
    #--------------------------------------------------------------------------
    # Emissions from landfill degradable wood
    #--------------------------------------------------------------------------
               
    # Turnover
    ip=13
    c_emitted=psl.HWP_LandfillWoodDegradable_tr*vo['C_Pro_Pools'][iT,:,ip]
    
    # Removal
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip]-c_emitted
    
    # Add to emissions (50% CO2 emissions during anaerobic decomposition)
    ip=17
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip]+psl.HWP_LandfillDegradableFracEmitCO2*c_emitted
    
    # Add to emissions (50% "potential" CH4 emissions during anaerobic decomposition)
    E_ch4_pot=(1-psl.HWP_LandfillDegradableFracEmitCO2)*c_emitted
    
    # Adjustment for proportion of degradable landfills with gas collection systems, 
    # efficiency of system, and methane oxided to CO2 from the landfill cover
    ch4_emitted=c_emitted*((1-psl.HWP_LandfillMethaneEmit_GasColSysProp)-psl.HWP_LandfillMethaneOxidizedToCO2*(1-psl.HWP_LandfillMethaneEmit_GasColSysProp)) + \
        c_emitted*psl.HWP_LandfillMethaneEmit_GasColSysProp*((1-psl.HWP_LandfillMethaneEmit_GasColSysEffic)-psl.HWP_LandfillMethaneOxidizedToCO2*(1-psl.HWP_LandfillMethaneEmit_GasColSysProp))
        
    ip=18
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip]+ch4_emitted
    
    #--------------------------------------------------------------------------
    # Emissions from landfill degradable paper
    #--------------------------------------------------------------------------
    
    # Turnover
    ip=15
    c_emitted=psl.HWP_LandfillWoodDegradable_tr*vo['C_Pro_Pools'][iT,:,ip]
    
    # Removal
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip]-c_emitted
    
    # Add to emissions (50% CO2 emissions during anaerobic decomposition)
    ip=17
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip]+psl.HWP_LandfillDegradableFracEmitCO2*c_emitted
    
    # Add to emissions (50% "potential" CH4 emissions during anaerobic decomposition)
    E_ch4_pot=(1-psl.HWP_LandfillDegradableFracEmitCO2)*c_emitted
    
    # Adjustment for proportion of degradable landfills with gas collection systems, 
    # efficiency of system, and methane oxided to CO2 from the landfill cover
    ch4_emitted=c_emitted*((1-psl.HWP_LandfillMethaneEmit_GasColSysProp)-psl.HWP_LandfillMethaneOxidizedToCO2*(1-psl.HWP_LandfillMethaneEmit_GasColSysProp)) + \
        c_emitted*psl.HWP_LandfillMethaneEmit_GasColSysProp*((1-psl.HWP_LandfillMethaneEmit_GasColSysEffic)-psl.HWP_LandfillMethaneOxidizedToCO2*(1-psl.HWP_LandfillMethaneEmit_GasColSysProp))
        
    ip=18
    vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip]+ch4_emitted
    
    return vo