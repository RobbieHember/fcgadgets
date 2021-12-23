
'''
DEMO UTILITIES
'''

#%% Import Python modules

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from matplotlib.patches import Rectangle
from fcgadgets.macgyver import utilities_general as gu
from fcgadgets.cbrunner import cbrun_utilities as cbu
  
#%%

def GetSingleEnsembleResults(meta):

    v0=[]
    for iScn in range(meta['Project']['N Scenario']):    
        v0.append(cbu.LoadSingleOutputFile(meta,iScn,0,0))
    
    return v0

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
        
        #v1[iScn]['Sum']['C_Forest']=v1[iScn]['Sum']['C_Biomass_Tot']+v1[iScn]['Sum']['C_DeadWood_Tot']+v1[iScn]['Sum']['C_Litter_Tot']+v1[iScn]['Sum']['C_Soil_Tot']
    
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

#%%

def ExportSummariesByScenario(meta,tv,t_start,t_end,mu1):
    
    it=np.where( (tv>=t_start) & (tv<=t_end) )[0]
    
    for iScn in range(meta['Project']['N Scenario']):
        
        VL=['A','LogSizeEnhancement','V_StemMerch','V_StemMerchToMill',
            'C_Biomass_Tot',
            'C_DeadWood_Tot',
            'C_DumpLandfill_Tot',
            'C_Eco_Pools',
            'C_Piled_Tot',
            'C_G_Gross',
            'C_G_Gross_Tot',
            'C_G_Net',
            'C_G_Net_Tot',
            'C_InUse_Tot',
            'C_LF',
            'C_LF_Tot',
            'C_Litter_Tot',
            'C_M_Dist',
            'C_M_Reg',
            'C_M_Reg_Tot',
            'C_NPP_Tot',
            'C_Pro_Pools',
            'C_RH',
            'C_RH_Tot',
            'C_Soil_Tot',
            'C_ToFirewoodDom',
            'C_ToFirewoodFor',
            'C_ToLogExport',
            'C_ToLumber',
            'C_ToMDF',
            'C_ToMill',
            'C_ToMillMerch',
            'C_ToMillNonMerch',
            'C_ToMillSnagStem',
            'C_ToOSB',
            'C_ToPaper',
            'C_ToPellets',
            'C_ToPlywood',
            'C_ToPowerFacilityDom',
            'C_ToPowerFacilityFor',
            'C_ToPowerGrid',
            'C_ToSlashpileBurn',
            'E_CO2e_LULUCF_NEE',
            'E_CO2e_LULUCF_EcoOther',
            'E_CO2e_LULUCF_Fire',
            'E_CO2e_LULUCF_OpenBurning',
            'E_CO2e_LULUCF_Wildfire',
            'E_CO2e_LULUCF_HWP',
            'E_CO2e_ESC_Comb',
            'E_CO2e_ESC_SubBM',
            'E_CO2e_ESC_SubE',
            'E_CO2e_ET_Comb',
            'E_CO2e_IPPU_Comb',
            'E_CO2e_AGHGB_WOSub',
            'E_CO2e_AGHGB_WOSub_cumu',
            'E_CO2e_AGHGB_WSub',
            'E_CO2e_AGHGB_WSub_cumu',
            'Cost Roads',
            'Cost Knockdown',
            'Cost Ripping',
            'Cost Nutrient Management',
            'Cost PAS Deactivation',
            'Cost Harvest Felling and Piling',
            'Cost Harvest Hauling',
            'Cost Harvest Overhead',
            'Cost Harvest Residuals',
            'Cost Milling',
            'Cost Slashpile Burn',            
            'Cost Planting',
            'Cost Survey',
            'Cost Silviculture Total',
            'Cost Total',            
            'Cost Total Disc',
            'Cost Total Disc_cumu',
            'Revenue FirewoodDom',
            'Revenue LogExport',
            'Revenue Lumber',
            'Revenue MDF',
            'Revenue OSB',
            'Revenue Paper',
            'Revenue Pellets',
            'Revenue Plywood',
            'Revenue PowerFacilityDom',
            'Revenue PowerGrid',
            'Revenue Gross',
            'Revenue Gross Disc',
            'Revenue Gross Disc_cumu',
            'Revenue Net',
            'Revenue Net Disc',
            'Revenue Net Disc_cumu',
            'Yield FirewoodDom',
            'Yield LogExport',
            'Yield Lumber',
            'Yield MDF',
            'Yield OSB',
            'Yield Paper',
            'Yield Pellets',
            'Yield Plywood',
            'Yield PowerFacilityDom','Yield PowerGrid']
        
        d={}        
        for k in VL:
            d['Sum ' + k]=np.round(np.sum(mu1[iScn][k][it]),decimals=2)
        
        for k in VL:
            d['Mean ' + k]=np.round(np.mean(mu1[iScn][k][it]),decimals=2)
        
        if iScn==0:
            df=pd.DataFrame().from_dict(d,orient='index')
        else:
            df0=pd.DataFrame().from_dict(d,orient='index')
            df=pd.concat([df,df0],axis=1)
    
    #df.index.name='Variable'
    df.columns=[np.arange(1,df.columns.size+1)]
    #df=df.sort_index(axis=0)
    
    df.to_excel(meta['Paths']['Project'] + '\\Outputs\\TabularSummary_' + str(t_start) + '-' + str(t_end) + '.xlsx')
    
    return df

#%%
    
def PlotSchematicAtmoGHGBal(meta,mu1,iB,iP,t_start,t_end):
    
    tv=np.arange(meta['Project']['Year Start Saving'],meta['Project']['Year End']+1,1)
    
    iT=np.where( (tv>=t_start) & (tv<=t_end) )[0]
    
    # Calculate aggregate variables
    for iScn in range(len(mu1)):
        mu1[iScn]['C_Forest']=mu1[iScn]['C_Biomass_Tot']+mu1[iScn]['C_DeadWood_Tot']+mu1[iScn]['C_Litter_Tot']+mu1[iScn]['C_Soil_Tot']
        mu1[iScn]['C_HWP']=mu1[iScn]['C_InUse_Tot']+mu1[iScn]['C_DumpLandfill_Tot']
        mu1[iScn]['C_ToMill']=mu1[iScn]['C_ToMillMerch']+mu1[iScn]['C_ToMillNonMerch']+mu1[iScn]['C_ToMillSnagStem']
    
    y_b={}
    y_p={}
    y_d={}
    for k in mu1[0].keys():
        if (k=='C_Forest') | (k=='C_HWP') | (k=='C_ToMill'):
            y_b[k]=mu1[iB][k][iT[-1]]-mu1[iB][k][iT[0]]
            y_p[k]=mu1[iP][k][iT[-1]]-mu1[iP][k][iT[0]]
            y_d[k]=y_p[k]-y_b[k]
        else:
            y_b[k]=np.sum(mu1[iB][k][iT])
            y_p[k]=np.sum(mu1[iP][k][iT])
            y_d[k]=y_p[k]-y_b[k]
        
        # Round
        y_b[k]=y_b[k].astype(int)
        y_p[k]=y_p[k].astype(int)
        y_d[k]=y_d[k].astype(int)
    
    bx_ec='none'
    bx_fs=9
    bx_fc=[0.93,0.93,0.93]
    bx2_fc=[0.9,0.9,0.9]
    bx_lower_h=0.47
    bx_lulucf_w=0.48
    bx_esc_w=0.23
    bx_atmo_bottom=0.88
    arrow_head_w=0.007
    arrow_lw=0.05
    fs_flux=6.5
    decim=1
    
    def GetSign(y):
        if y>0:
            x='+'
        else:
            x=''
        return x
    
    fig,ax=plt.subplots(1,figsize=gu.cm2inch(18,10))
    
    # Background
    #ax.add_patch(Rectangle([0,1],0,1,fc='w',ec='k'))
    
    # Atmosphere
    ax.add_patch(Rectangle([0.01,bx_atmo_bottom],0.98,0.1,fc=bx_fc,ec=bx_ec))
    ax.text(0.5,0.935,'Atmosphere',size=bx_fs,ha='center')
    vr='E_CO2e_AGHGB_WSub'
    a1=np.round(y_b[vr],decimals=decim); a2=np.round(y_p[vr],decimals=decim); a3=np.round(y_d[vr],decimals=decim)
    txt=str(a1) + ',' + str(a2) + ' (' + GetSign(a3) + str(a3) + ')'
    ax.text(0.5,0.9,txt,size=fs_flux,ha='center')
    
    # LULUCF
    ax.add_patch(Rectangle([0.01,0.01],bx_lulucf_w,bx_lower_h,fc=bx_fc,ec=bx_ec))
    ax.text(0.25,0.04,'Land Use, Land Use Change and Forestry',size=bx_fs,ha='center')
    
    # Forest land
    ax.add_patch(Rectangle([0.02,0.1],bx_lulucf_w*0.53,bx_lower_h-0.11,fc=bx2_fc,ec=bx_ec))
    ax.text(0.15,0.28,'Forest Land',size=bx_fs,ha='center')    
    vr='C_Forest'
    a1=np.round(y_b[vr],decimals=decim); a2=np.round(y_p[vr],decimals=decim); a3=np.round(y_d[vr],decimals=decim)
    txt=str(a1) + ',' + str(a2) + ' (' + GetSign(a3) + str(a3) + ')'
    ax.text(0.15,0.24,txt,size=fs_flux,ha='center')
    
    # Harvested wood products
    ax.add_patch(Rectangle([0.36,0.1],0.12,bx_lower_h-0.11,fc=bx2_fc,ec=bx_ec))
    ax.text(0.42,0.24,'Harvested\nWood\nProducts',size=bx_fs,ha='center')
    vr='C_HWP'
    a1=np.round(y_b[vr],decimals=decim); a2=np.round(y_p[vr],decimals=decim); a3=np.round(y_d[vr],decimals=decim)
    txt=str(a1) + ',' + str(a2) + ' (' + GetSign(a3) + str(a3) + ')'
    ax.text(0.42,0.20,txt,size=fs_flux,ha='center')

    # Energy - Stationary Combustion
    ax.add_patch(Rectangle([bx_lulucf_w+0.02,0.01],bx_esc_w,bx_lower_h,fc=bx_fc,ec=bx_ec))
    ax.text(0.62,0.21,'Energy\nStationary\nCombustion',size=bx_fs,ha='center')
    
    # Energy - Transportation
    ax.add_patch(Rectangle([bx_lulucf_w+bx_esc_w+0.03,0.01],0.12,bx_lower_h,fc=bx_fc,ec=bx_ec))
    ax.text(0.8,0.23,'Energy\nTransportation',size=bx_fs,ha='center')
    
    # IPPU
    ax.add_patch(Rectangle([0.87,0.01],0.12,bx_lower_h,fc=bx_fc,ec=bx_ec))
    ax.text(0.93,0.18,'Industrial\nProducts\n&\nProduct\nUse',size=bx_fs,ha='center')
    
    # NEE
    vr='E_CO2e_LULUCF_NEE'
    a1=np.round(y_b[vr],decimals=decim); a2=np.round(y_p[vr],decimals=decim); a3=np.round(y_d[vr],decimals=decim)
    txt='Net ecosystem\nexchange\n' + str(a1) + ',' + str(a2) + ' (' + GetSign(a3) + str(a3) + ')'
    ax.text(0.02,0.79,txt,ha='left',size=fs_flux)
    ax.arrow(0.015,bx_atmo_bottom,0,-1*(bx_atmo_bottom-bx_lower_h-0.02),head_width=arrow_head_w,head_length=0.01,fc='k',ec='k',lw=arrow_lw)
    
    # Volatilization
    vr='E_CO2e_LULUCF_EcoOther'
    a1=np.round(y_b[vr],decimals=decim); a2=np.round(y_p[vr],decimals=decim); a3=np.round(y_d[vr],decimals=decim)
    txt='Volatilization\nand\ndenitrificaiton\n' + str(a1) + ',' + str(a2) + ' (' + GetSign(a3) + str(a3) + ')'
    ax.text(0.115,0.52,txt,ha='right',size=fs_flux)
    ax.arrow(0.12,bx_lower_h+0.01,0,bx_atmo_bottom-bx_lower_h-0.02,head_width=arrow_head_w,head_length=0.01,fc='k',ec='k',lw=arrow_lw)
    
    # Wildfire
    vr='E_CO2e_LULUCF_Wildfire'
    a1=np.round(y_b[vr],decimals=decim); a2=np.round(y_p[vr],decimals=decim); a3=np.round(y_d[vr],decimals=decim)
    txt='Wildfire\n' + str(a1) + ',' + str(a2) + ' (' + str(a3) + ')'
    ax.text(0.15,0.52,txt,ha='left',size=fs_flux)
    ax.arrow(0.145,bx_lower_h+0.01,0,bx_atmo_bottom-bx_lower_h-0.02,head_width=arrow_head_w,head_length=0.01,fc='k',ec='k',lw=arrow_lw)
    
    # Denitrification
    #a1=0.2; a2=0.5; a3=0.3
    #txt='Denitrification\n' + str(a1) + ',' + str(a2) + ' (' + str(a3) + ')'
    #ax.text(0.175,0.7,txt,ha='left',size=fs_flux)
    #ax.arrow(0.17,bx_lower_h+0.15,0,bx_atmo_bottom-bx_lower_h-0.15-0.01,head_width=arrow_head_w,head_length=0.01,fc='k',ec='k',lw=arrow_lw)
    
    # Open burning
    vr='E_CO2e_LULUCF_OpenBurning'
    a1=np.round(y_b[vr],decimals=decim); a2=np.round(y_p[vr],decimals=decim); a3=np.round(y_d[vr],decimals=decim)
    txt='Open burning\n' + str(a1) + ',' + str(a2) + ' (' + GetSign(a3) + str(a3) + ')'
    ax.text(0.27,0.52,txt,ha='left',size=fs_flux)
    ax.arrow(0.265,bx_lower_h+0.01,0,bx_atmo_bottom-bx_lower_h-0.02,head_width=arrow_head_w,head_length=0.01,fc='k',ec='k',lw=arrow_lw)
    
    # HWP fluxes
    vr='E_CO2e_LULUCF_HWP'
    a1=np.round(y_b[vr],decimals=decim); a2=np.round(y_p[vr],decimals=decim); a3=np.round(y_d[vr],decimals=decim)
    txt='Product decay and\ncombustion\n' + str(a1) + ',' + str(a2) + ' (' + GetSign(a3) + str(a3) + ')'
    ax.text(0.435,0.64,txt,ha='right',size=fs_flux)
    ax.arrow(0.44,bx_lower_h+0.01,0,bx_atmo_bottom-bx_lower_h-0.02,head_width=arrow_head_w,head_length=0.01,fc='k',ec='k',lw=arrow_lw)
    
    # Removals
    vr='C_ToMill'
    a1=np.round(y_b[vr],decimals=decim); a2=np.round(y_p[vr],decimals=decim); a3=np.round(y_d[vr],decimals=decim)
    #a1=0.2; a2=0.5; a3=0.3
    txt='Removals\n' + str(a1) + ',' + str(a2) + ' (' + GetSign(a3) + str(a3) + ')'
    ax.text(0.32,0.33,txt,ha='center',size=fs_flux)
    ax.arrow(0.28,0.3,0.07,0,head_width=0.01,head_length=arrow_head_w,fc='k',ec='k',lw=arrow_lw)
    
    # Bioenergy combustion
    vr='E_CO2e_ESC_Comb'
    a1=np.round(y_b[vr],decimals=decim); a2=np.round(y_p[vr],decimals=decim); a3=np.round(y_d[vr],decimals=decim)
    txt='Bioenergy\ncombustion\n' + str(a1) + ',' + str(a2) + ' (' + GetSign(a3) + str(a3) + ')'
    ax.text(0.515,0.5,txt,ha='left',size=fs_flux)
    ax.arrow(0.51,bx_lower_h+0.01,0,bx_atmo_bottom-bx_lower_h-0.02,head_width=arrow_head_w,head_length=0.01,fc='k',ec='k',lw=arrow_lw)
    
    # Displacement energy
    vr='E_CO2e_ESC_SubE'
    a1=np.round(y_b[vr],decimals=decim); a2=np.round(y_p[vr],decimals=decim); a3=np.round(y_d[vr],decimals=decim)
    txt='Displacement \nof\nnon-renewable\nenergy\nsources\n' + str(a1) + ',' + str(a2) + ' (' + GetSign(a3) + str(a3) + ')'
    ax.text(0.615,0.72,txt,ha='left',size=fs_flux)
    ax.arrow(0.61,bx_atmo_bottom,0,-1*(bx_atmo_bottom-bx_lower_h-0.02),head_width=arrow_head_w,head_length=0.01,fc='k',ec='k',lw=arrow_lw)
    
    # Displacement building materials
    vr='E_CO2e_ESC_SubBM'
    a1=np.round(y_b[vr],decimals=decim); a2=np.round(y_p[vr],decimals=decim); a3=np.round(y_d[vr],decimals=decim)
    txt='Displacement \nof\nnon-renewable\nenergy\nsources\nby solid\nwood\nproducts\n' + str(a1) + ',' + str(a2) + ' (' + GetSign(a3) + str(a3) + ')'
    ax.text(0.715,0.65,txt,ha='left',size=fs_flux)
    ax.arrow(0.71,bx_atmo_bottom,0,-1*(bx_atmo_bottom-bx_lower_h-0.02),head_width=arrow_head_w,head_length=0.01,fc='k',ec='k',lw=arrow_lw)
    
    # Combustion from transportation
    vr='E_CO2e_ET_Comb'
    a1=np.round(y_b[vr],decimals=decim); a2=np.round(y_p[vr],decimals=decim); a3=np.round(y_d[vr],decimals=decim)
    txt='Transportation\nfuel\ncombustion\n' + str(a1) + ',' + str(a2) + ' (' + GetSign(a3) + str(a3) + ')'
    ax.text(0.815,0.51,txt,ha='left',size=fs_flux)
    ax.arrow(0.81,bx_lower_h+0.01,0,bx_atmo_bottom-bx_lower_h-0.02,head_width=arrow_head_w,head_length=0.01,fc='k',ec='k',lw=arrow_lw)
    
    # Sequestration from IPPU
    vr='E_CO2e_IPPU_Comb'
    a1=np.round(y_b[vr],decimals=decim); a2=np.round(y_p[vr],decimals=decim); a3=np.round(y_d[vr],decimals=decim)
    txt='Urea\nsequestration\n' + str(a1) + ',' + str(a2) + ' (' + GetSign(a3) + str(a3) + ')'
    ax.text(0.915,0.79,txt,ha='left',size=fs_flux)
    ax.arrow(0.91,bx_atmo_bottom,0,-1*(bx_atmo_bottom-bx_lower_h-0.02),head_width=arrow_head_w,head_length=0.01,fc='k',ec='k',lw=arrow_lw)
    
    ax.set(position=[0,0,1,1],visible='Off',xticks=[],yticks=[])
    
    gu.PrintFig(meta['Paths']['Figures'] + '\\AGHGB Schematic_S' + str(iP) + 'minusS' + str(iB) + '_' + str(t_start) + 'to' + str(t_end),'png',900)
    
    return