
'''
DEMO UTILITIES
'''

#%% Import Python modules

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from fcgadgets.macgyver import utilities_general as gu
from fcgadgets.cbrunner import cbrun_utilities as cbu

#%%
    
def CalculateAggregateVariables(meta,v1):
    
    for iScn in range(meta['Project']['N Scenario']):
        
        # Calculate carbon content of dead wood, organic and mineral soil horizons following Shaw et al. (2017)
        v1[iScn]['SoilOrgH']=v1[iScn]['C_Eco_Pools'][:,:,meta['Core']['iEP']['LitterVF']]+ \
            v1[iScn]['C_Eco_Pools'][:,:,meta['Core']['iEP']['LitterS']]
        v1[iScn]['SoilMinH']=v1[iScn]['C_Eco_Pools'][:,:,meta['Core']['iEP']['SoilVF']]+ \
            v1[iScn]['C_Eco_Pools'][:,:,meta['Core']['iEP']['SoilS']]+ \
            v1[iScn]['C_Eco_Pools'][:,:,meta['Core']['iEP']['SoilS']]
        v1[iScn]['DeadWood']=v1[iScn]['C_Eco_Pools'][:,:,meta['Core']['iEP']['SnagStem']]+ \
            v1[iScn]['C_Eco_Pools'][:,:,meta['Core']['iEP']['SnagBranch']]+ \
            v1[iScn]['C_Eco_Pools'][:,:,meta['Core']['iEP']['LitterM']]+ \
            v1[iScn]['C_Eco_Pools'][:,:,meta['Core']['iEP']['LitterF']]

        v1[iScn]['C_BiomassAG_Tot']=np.sum(v1[iScn]['C_Eco_Pools'][:,:,meta['Core']['iEP']['BiomassAboveground']],axis=2)
        v1[iScn]['C_BiomassBG_Tot']=np.sum(v1[iScn]['C_Eco_Pools'][:,:,meta['Core']['iEP']['BiomassBelowground']],axis=2)
        v1[iScn]['C_BiomassAG_Tot'][np.isnan(v1[iScn]['C_BiomassAG_Tot'])]=0
        
        v1[iScn]['C_Eco_Tot']=np.sum(v1[iScn]['C_Eco_Pools'],axis=2)
        v1[iScn]['C_Eco_Tot']=np.sum(v1[iScn]['C_Eco_Pools'],axis=2)
        v1[iScn]['C_Eco_Tot'][np.isnan(v1[iScn]['C_Eco_Tot'])]=0
    
    return v1

#%% Tabulate responses based on scenario comparison of specified time period

def CompareScenarios(meta,v1,iB,iP,iT):

    #--------------------------------------------------------------------------
    # Percent response of fluxes
    #--------------------------------------------------------------------------
    
    dFluxRel={}
    dFluxRel['NPP']=np.round(np.mean((v1[iP]['C_NPP_Tot'][iT]-v1[0]['C_NPP_Tot'][iT])/v1[iB]['C_NPP_Tot'][iT]*100))
    
    dFluxRel['Gross Growth']=np.round(np.mean((v1[iP]['C_G_Gross_Tot'][iT]-v1[0]['C_G_Gross_Tot'][iT])/v1[iB]['C_G_Gross_Tot'][iT]*100))
    
    cb0=np.mean(v1[iB]['C_G_Gross'][iT,0,meta['Core']['iEP']['Foliage']])
    cp0=np.mean(v1[iP]['C_G_Gross'][iT,0,meta['Core']['iEP']['Foliage']])
    dFluxRel['Foliage production']=np.round((cp0-cb0)/cb0*100)
    
    ind=np.where(np.array(meta['Core']['Name Pools Eco'])=='RootFine')[0]
    cb0=np.mean(v1[iB]['C_G_Gross'][iT,0,ind])
    cp0=np.mean(v1[iP]['C_G_Gross'][iT,0,ind])
    dFluxRel['Fine root production']=np.round((cp0-cb0)/cb0*100)
    
    dFluxRel['Net Growth']=np.round(np.mean((v1[iP]['C_G_Net_Tot'][iT]-v1[iB]['C_G_Net_Tot'][iT])/v1[iB]['C_G_Net_Tot'][iT]*100))
    
    dFluxRel['Tree Mortality']=np.round(np.mean((v1[iP]['C_M_Reg_Tot'][iT]-v1[iB]['C_M_Reg_Tot'][iT])/v1[iB]['C_M_Reg_Tot'][iT]*100))
    
    ind=np.where(np.array(meta['Core']['Name Pools Eco'])=='Foliage')[0]
    cb0=np.mean(v1[iB]['C_LF'][iT,0,ind]);cp0=np.mean(v1[iP]['C_LF'][iT,0,ind])
    ind=np.where(np.array(meta['Core']['Name Pools Eco'])=='Branch')[0]
    cb0=cb0+np.mean(v1[iB]['C_LF'][iT,0,ind]);cp0=cp0+np.mean(v1[iP]['C_LF'][iT,0,ind])
    dFluxRel['Foliage+Branch turnover']=np.round((cp0-cb0)/cb0*100)
    
    ind=np.where(np.array(meta['Core']['Name Pools Eco'])=='RootCoarse')[0]
    cb0=np.mean(v1[iB]['C_LF'][iT,0,ind]);cp0=np.mean(v1[iP]['C_LF'][iT,0,ind])
    dFluxRel['Coarse root turnover']=np.round((cp0-cb0)/cb0*100)
    
    ind=np.where(np.array(meta['Core']['Name Pools Eco'])=='RootFine')[0]
    cb0=np.mean(v1[iB]['C_LF'][iT,0,ind]);cp0=np.mean(v1[iP]['C_LF'][iT,0,ind])
    dFluxRel['Fine root turnover']=np.round((cp0-cb0)/cb0*100)
    
    dFluxRel['Litterfall']=np.round(np.mean((v1[1]['C_LF_Tot'][iT]-v1[0]['C_LF_Tot'][iT])/v1[0]['C_LF_Tot'][iT]*100))
    
    dFluxRel['RH']=np.round(np.mean((v1[1]['C_RH_Tot'][iT]-v1[0]['C_RH_Tot'][iT])/v1[0]['C_RH_Tot'][iT]*100))
    
    #cb0=np.mean(v1[1]['C_Eco_Pools'][iT,0,28]);cp0=np.mean(v1[0]['C_Eco_Pools'][iT,0,28])
    #dFluxRel['Litter Decomp']=np.round(np.mean((cp0-cb0)/cb0*100))
    
    #--------------------------------------------------------------------------
    # Percent respoonse of biomass pools
    #--------------------------------------------------------------------------
    
    s=['StemMerch','StemNonMerch','Foliage','Branch','Bark','RootCoarse','RootFine']
    dr0=np.zeros(len(s));cb=np.zeros(len(s));cp=np.zeros(len(s))
    for i in range(len(s)):
        ind=np.where(np.array(meta['Core']['Name Pools Eco'])==s[i])[0]
        cb0=np.mean(v1[iB]['C_G_Net'][iT,0,ind])
        cp0=np.mean(v1[iP]['C_G_Net'][iT,0,ind])
        dr0[i]=(cp0-cb0)/cb0*100
        cb[i]=cb0; cp[i]=cp0
    
    Stem_b=cb[0]+cb[1]+cb[4]; Stem_p=cp[0]+cp[1]+cp[4]
    
    dPoolRel={}
    dPoolRel['Stemwood']=np.round((Stem_p-Stem_b)/Stem_b*100)
    dPoolRel['Branch']=np.round(dr0[3])
    dPoolRel['Foliage']=np.round(dr0[2])
    dPoolRel['Coarse roots']=np.round(dr0[5])
    dPoolRel['Fine roots']=np.round(dr0[6])
        
    yP=v1[iP]['C_Biomass_Tot'][iT]
    yB=v1[iB]['C_Biomass_Tot'][iT]    
    dPoolRel['Biomass Total']=np.round(np.mean((yP-yB)/yB*100))
    
    dPoolRel['Dead Wood']=np.round(np.mean((v1[iP]['C_DeadWood_Tot'][iT]-v1[iB]['C_DeadWood_Tot'][iT])/v1[iB]['C_DeadWood_Tot'][iT]*100))
    dPoolRel['Litter']=np.round(np.mean(((v1[iP]['C_Litter_Tot'][iT])-(+v1[iB]['C_Litter_Tot'][iT]))/(v1[iB]['C_Litter_Tot'][iT])*100))
    dPoolRel['Soil organic horizon']=np.round(np.mean((v1[iP]['SoilOrgH'][iT]-v1[iB]['SoilOrgH'][iT])/v1[iB]['SoilOrgH'][iT]*100))
    dPoolRel['Soil mineral Horizon']=np.round(np.mean((v1[iP]['SoilMinH'][iT]-v1[iB]['SoilMinH'][iT])/v1[iB]['SoilMinH'][iT]*100))
    dPoolRel['Soil Total']=np.round(np.mean(((v1[iP]['C_Soil_Tot'][iT]+v1[iP]['C_Litter_Tot'][iT])-(v1[iB]['C_Soil_Tot'][iT]+v1[iB]['C_Litter_Tot'][iT]))/(v1[iB]['C_Soil_Tot'][iT]+v1[iB]['C_Litter_Tot'][iT])*100))
    
    #--------------------------------------------------------------------------
    # Actual respoonse of fluxes
    #--------------------------------------------------------------------------
    
    nam=['StemMerch','StemNonMerch','Foliage','Branch','Bark','RootCoarse','RootFine']
    dFluxAct={}
    for i in range(len(nam)):
        ind=np.where(np.array(meta['Core']['Name Pools Eco'])==nam[i])[0]
        dFluxAct[nam[i]]=np.mean(v1[iP]['C_G_Net'][iT,0,ind])-np.mean(v1[iB]['C_G_Net'][iT,0,ind])
    
    dFluxAct['Stem Total']=dFluxAct['StemMerch']+dFluxAct['StemNonMerch']
    
    y_b=v1[iB]['C_NPP_Tot'][iT,0]-v1[iB]['C_RH_Tot'][iT,0]
    y_p=v1[iP]['C_NPP_Tot'][iT,0]-v1[iP]['C_RH_Tot'][iT,0]
    dFluxAct['NEP']=np.mean(y_p-y_b)
    
    #--------------------------------------------------------------------------
    # Stemwood mortality
    #--------------------------------------------------------------------------
    
    # Actual and percent response of stemwood mortality for ten years following application
    ind=np.where(np.array(meta['Core']['Name Pools Eco'])=='StemMerch')[0]
    cb0=np.mean(v1[iB]['C_M_Reg'][iT,0,ind])
    cp0=np.mean(v1[iP]['C_M_Reg'][iT,0,ind])
    
    dStemMort={}
    dStemMort['Rel']=(cp0-cb0)/cb0*100
    dStemMort['Act']=cp0-cb0
    
    #--------------------------------------------------------------------------
    # Merch volume
    #--------------------------------------------------------------------------
    
    cb0=np.max(v1[iB]['V_StemMerch'][iT,0])
    cp0=np.max(v1[iP]['V_StemMerch'][iT,0])
    
    dMerchVolume={}
    dMerchVolume['Act']=cp0-cb0; 
    dMerchVolume['Rel']=(cp0-cb0)/cb0*100
    
    #--------------------------------------------------------------------------
    # Nitrogen use efficiency (applied)
    #--------------------------------------------------------------------------
    
    N=200 # NUE applied
    
    dcStemG=(np.mean(v1[iP]['C_G_Net'][iT,0,0]-v1[iB]['C_G_Net'][iT,0,0]))*1000
    dcStem=(np.mean(v1[iP]['C_Eco_Pools'][iT,0,0]-v1[iB]['C_Eco_Pools'][iT,0,0]))*1000
    dcStem=dcStem+(np.mean(v1[iP]['C_Eco_Pools'][iT,0,1]-v1[iB]['C_Eco_Pools'][iT,0,1]))*1000
    dcStem=dcStem+(np.mean(v1[iP]['C_Eco_Pools'][iT,0,4]-v1[iB]['C_Eco_Pools'][iT,0,4]))*1000
    
    dNUE_applied={}
    dNUE_applied['Stemwood']=np.round(dcStem/N)
    
    dcFoliage=(np.mean(v1[iP]['C_Eco_Pools'][iT,0,2]-v1[iB]['C_Eco_Pools'][iT,0,2]))*1000
    dNUE_applied['Foliage']=np.round(dcFoliage/N)
    
    dcBranch=(np.mean(v1[iP]['C_Eco_Pools'][iT,0,3]-v1[iB]['C_Eco_Pools'][iT,0,3]))*1000
    dNUE_applied['Branch']=np.round(dcBranch/N)
    
    dcRC=(np.mean(v1[iP]['C_Eco_Pools'][iT,0,5]-v1[iB]['C_Eco_Pools'][iT,0,5]))*1000
    dNUE_applied['Coarse root']=np.round(dcRC/N)
    
    dcRF=(np.mean(v1[iP]['C_Eco_Pools'][iT,0,6]-v1[iB]['C_Eco_Pools'][iT,0,6]))*1000
    dNUE_applied['Fine root']=np.round(dcRC/N)
    
    dcDW=(np.mean(v1[iP]['C_DeadWood_Tot'][iT]-v1[iB]['C_DeadWood_Tot'][iT]))*1000
    dNUE_applied['Dead Wood']=dcDW/N
    
    dcL=(np.mean(v1[iP]['C_Litter_Tot'][iT]-v1[iB]['C_Litter_Tot'][iT]))*1000
    dNUE_applied['Litter']=dcL/N
    
    dcS=(np.mean(v1[iP]['C_Soil_Tot'][iT]-v1[iB]['C_Soil_Tot'][iT]))*1000
    dNUE_applied['Soil']=dcS/N
    
    dcTot=(np.mean(v1[iP]['C_Eco_Tot'][iT]-v1[iB]['C_Eco_Tot'][iT]))*1000
    dNUE_applied['Total']=dcTot/N
    
    dNUE_applied['Biomass']=dNUE_applied['Stemwood']+dNUE_applied['Foliage']+dNUE_applied['Branch']+dNUE_applied['Coarse root']+dNUE_applied['Fine root']
    
    #--------------------------------------------------------------------------
    # Nitrogen use efficiency (utilized)
    #--------------------------------------------------------------------------
    
    N=40 # NUE utilized
    
    dcStemG=(np.mean(v1[iP]['C_G_Net'][iT,0,0]-v1[iB]['C_G_Net'][iT,0,0]))*1000
    dcStem=(np.mean(v1[iP]['C_Eco_Pools'][iT,0,0]-v1[iB]['C_Eco_Pools'][iT,0,0]))*1000
    dcStem=dcStem+(np.mean(v1[iP]['C_Eco_Pools'][iT,0,1]-v1[iB]['C_Eco_Pools'][iT,0,1]))*1000
    dcStem=dcStem+(np.mean(v1[iP]['C_Eco_Pools'][iT,0,4]-v1[iB]['C_Eco_Pools'][iT,0,4]))*1000
    
    dNUE_utilized={}
    dNUE_utilized['Stemwood']=np.round(dcStem/N)
    dcFoliage=(np.mean(v1[iP]['C_Eco_Pools'][iT,0,2]-v1[iB]['C_Eco_Pools'][iT,0,2]))*1000
    dNUE_utilized['Foliage']=np.round(dcFoliage/N)
    dcBranch=(np.mean(v1[iP]['C_Eco_Pools'][iT,0,3]-v1[iB]['C_Eco_Pools'][iT,0,3]))*1000
    dNUE_utilized['Branch']=np.round(dcBranch/N)
    dcRC=(np.mean(v1[iP]['C_Eco_Pools'][iT,0,5]-v1[iB]['C_Eco_Pools'][iT,0,5]))*1000
    dNUE_utilized['Coarse root']=np.round(dcRC/N)
    dcRF=(np.mean(v1[iP]['C_Eco_Pools'][iT,0,6]-v1[iB]['C_Eco_Pools'][iT,0,6]))*1000
    dNUE_utilized['Fine root']=np.round(dcRC/N)
    
    dcDW=(np.mean(v1[iP]['C_DeadWood_Tot'][iT]-v1[iB]['C_DeadWood_Tot'][iT]))*1000
    dNUE_utilized['Dead Wood']=dcDW/N
    
    dcL=(np.mean(v1[iP]['C_Litter_Tot'][iT]-v1[iB]['C_Litter_Tot'][iT]))*1000
    dNUE_utilized['Litter']=dcL/N
    
    dcS=(np.mean(v1[iP]['C_Soil_Tot'][iT]-v1[iB]['C_Soil_Tot'][iT]))*1000
    dNUE_utilized['Soil']=dcS/N
    
    dcTot=(np.mean(v1[iP]['C_Eco_Tot'][iT]-v1[iB]['C_Eco_Tot'][iT]))*1000
    dNUE_utilized['Total']=dcTot/N
    
    dNUE_utilized['Biomass']=dNUE_utilized['Stemwood']+dNUE_utilized['Foliage']+dNUE_utilized['Branch']+dNUE_utilized['Coarse root']+dNUE_utilized['Fine root']
    
    return dFluxRel,dFluxAct,dPoolRel,dStemMort,dMerchVolume,dNUE_applied,dNUE_utilized


#%% Plot cashflow

def Plot_Cashflow(meta,v1,econ,iB,iP,iT):
    
    lw=0.75
    
    fig,ax=plt.subplots(3,2,figsize=gu.cm2inch(15,11.5))
    
    ax[0,0].plot(v1[0]['Year'],econ[iB]['Cost Total']/1000,'-bo',lw=lw,ms=4,label='Baseline')
    ax[0,0].plot(v1[0]['Year'],econ[iP]['Cost Total']/1000,'--r^',lw=lw,ms=3,label='Project')
    ax[0,0].set(xlim=[v1[0]['Year'][iT[0]], v1[0]['Year'][iT[-1]]],ylabel='Cost (CDN$/000)')
    ax[0,0].legend(loc='upper right',frameon=False,facecolor=None)
    
    ax[0,1].plot(v1[0]['Year'],(econ[iP]['Cost Total']-econ[iB]['Cost Total'])/1000,'-g',lw=lw)
    ax[0,1].set(xlim=[v1[0]['Year'][iT[0]], v1[0]['Year'][iT[-1]]],ylabel='Cost (CDN$/000)')
    
    ax[1,0].plot(v1[0]['Year'],econ[iB]['Revenue Net']/1000,'-bo',ms=4,lw=lw)
    ax[1,0].plot(v1[0]['Year'],econ[iP]['Revenue Net']/1000,'--r^',ms=3,lw=lw)
    ax[1,0].set(xlim=[v1[0]['Year'][iT[0]], v1[0]['Year'][iT[-1]]],ylabel='Net revenue (CDN$/000)')
    
    ax[1,1].plot(v1[0]['Year'],(econ[iP]['Revenue Net']-econ[iB]['Revenue Net'])/1000,'-g',lw=lw)
    ax[1,1].set(xlim=[v1[0]['Year'][iT[0]], v1[0]['Year'][iT[-1]]],ylabel='Net revenue (CDN$/000)')
    
    ax[2,0].plot(v1[0]['Year'],np.cumsum(econ[iB]['Revenue Net'])/1000,'-b',ms=4,lw=lw)
    ax[2,0].plot(v1[0]['Year'],np.cumsum(econ[iP]['Revenue Net'])/1000,'--r',ms=3,lw=lw)
    ax[2,0].set(xlim=[v1[0]['Year'][iT[0]], v1[0]['Year'][iT[-1]]],ylabel='Cumulative net revenue (CDN$/000)')
    
    ax[2,1].plot(v1[0]['Year'],(np.cumsum(econ[iP]['Revenue Net'])-np.cumsum(econ[iB]['Revenue Net']))/1000,'-g',lw=lw)
    ax[2,1].set(xlim=[v1[0]['Year'][iT[0]], v1[0]['Year'][iT[-1]]],ylabel='Cumulative net revenue (CDN$/000)');
    
    return
