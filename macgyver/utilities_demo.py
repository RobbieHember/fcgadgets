
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
            v1[iScn]['C_Eco_Pools'][:,:,meta['Core']['iEP']['SoilF']]+ \
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

    cb0=np.max(v1[iB]['V_MerchTotal'][iT,0])
    cp0=np.max(v1[iP]['V_MerchTotal'][iT,0])

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

def ExportSummariesByScenario(meta,mos,t_start,t_end,**kwargs):

    tv=np.arange(meta['Project']['Year Start Saving'],meta['Project']['Year End']+1,1)

    iT=np.where( (tv>=t_start) & (tv<=t_end) )[0]

    # Key word arguments
    if 'iSP' in kwargs.keys():
        iSP=kwargs['iSP']
    else:
        iSP=0

    if 'iSS' in kwargs.keys():
        iSS=kwargs['iSS']
    else:
        iSS=0

    if 'sum_mult' in kwargs.keys():
        sum_mult=kwargs['sum_mult']
    else:
        sum_mult=1.0

    if (meta['Project']['Scenario Source']!='Portfolio'):

        for iScn in range(meta['Project']['N Scenario']):

            d={}

            for k in meta['Core']['Output Variable List']:
                try:
                    d['Annual mean summed over area ' + k]=np.round(sum_mult*np.mean(mos['Scenarios'][iScn]['Sum'][k]['Ensemble Mean'][iT,iSP,iSS]),decimals=2)
                except:
                    pass

            for k in meta['Core']['Output Variable List']:
                try:
                    d['Per-hectare sum over time ' + k]=np.round(sum_mult*np.sum(mos['Scenarios'][iScn]['Mean'][k]['Ensemble Mean'][iT,iSP,iSS]),decimals=2)
                except:
                    pass

            for k in meta['Core']['Output Variable List']:
                try:
                    d['Per-hectare mean ' + k]=np.round(np.mean(mos['Scenarios'][iScn]['Mean'][k]['Ensemble Mean'][iT,iSP,iSS]),decimals=2)
                except:
                    pass

            if iScn==0:
                df=pd.DataFrame().from_dict(d,orient='index');
            else:
                df0=pd.DataFrame().from_dict(d,orient='index');
                df=pd.concat([df,df0],axis=1);

    else:

        for iPort in range(meta['Project']['N Portfolio']):

            for iScn in range(meta['Project']['N Scenario']):

                d={}

                for k in meta['Core']['Output Variable List']:
                    try:
                        d['Annual mean summed over area ' + k]=np.round(sum_mult*np.mean(mos[iPort][iScn]['Sum'][k]['Ensemble Mean'][iT]),decimals=2)
                    except:
                        pass

                for k in meta['Core']['Output Variable List']:
                    try:
                        d['Per-hectare sum over time ' + k]=np.round(sum_mult*np.sum(mos[iPort][iScn]['Mean'][k]['Ensemble Mean'][iT]),decimals=2)
                    except:
                        pass

                for k in meta['Core']['Output Variable List']:
                    try:
                        d['Per-hectare mean ' + k]=np.round(np.mean(mos[iPort][iScn]['Mean'][k]['Ensemble Mean'][iT]),decimals=2)
                    except:
                        pass

                if (iPort==0) & (iScn==0):
                    df=pd.DataFrame().from_dict(d,orient='index');
                else:
                    df0=pd.DataFrame().from_dict(d,orient='index');
                    df=pd.concat([df,df0],axis=1);

    df.columns=[np.arange(1,df.columns.size+1)];

    fout=meta['Paths']['Project'] + '\\Outputs\\TabularSummary_' + str(t_start) + '-' + str(t_end) + '_ProjectType' + str(iSP) + '_Region' + str(iSS) + '.xlsx'
    df.to_excel(fout);

    return df

#%% Custome scenario comparison tabular export

# Exmaple:
# tabNam='ScenarioComparisons1'
# thL=[50,100]
# vL=['Mean*IntSum*V_ToMillMerchTotal','Mean*Inst*E_CO2e_AGHGB_WSub_cumu_from_tref']
# udem.ExportDeltaTable(meta,mos,tabNam,thL,vL)

def ExportDeltaTable(meta,mos,tabNam,t_start,t_duration,vL,**kwargs):

    tv=np.arange(meta['Project']['Year Start Saving'],meta['Project']['Year End']+1,1)

    # Key word arguments
    if 'iSP' in kwargs.keys():
        iSP=kwargs['iSP']
    else:
        iSP=0

    if 'iSS' in kwargs.keys():
        iSS=kwargs['iSS']
    else:
        iSS=0

    if 'sum_mult' in kwargs.keys():
        sum_mult=kwargs['sum_mult']
    else:
        sum_mult=1.0

    df=pd.DataFrame()

    for sc in mos['Delta'].keys():

        for th in t_duration:

            iT=np.where( (tv>=t_start) & (tv<=t_start+th-1) )[0]

            d={}
            for v in vL:

                ind0=v.find('*')
                ind1=v.rfind('*')
                op0=v[0:ind0]
                op1=v[ind0+1:ind1]
                vnam=v[ind1+1:]

                if op1=='Inst':
                    # Instantaneous
                    y=sum_mult*mos['Delta'][sc]['ByStrata'][op0][vnam]['Ensemble Mean'][iT[-1],iSP,iSS]
                elif op1=='IntSum':
                    # Integrated sum
                    y=sum_mult*np.sum(mos['Delta'][sc]['ByStrata'][op0][vnam]['Ensemble Mean'][iT,iSP,iSS])
                elif op1=='MeanAnnualSumOverArea':
                    # Mean annual sum
                    y=sum_mult*np.mean(mos['Delta'][sc]['ByStrata'][op0][vnam]['Ensemble Mean'][iT,iSP,iSS])
                else:
                    # Integrated mean
                    y=np.mean(mos['Delta'][sc]['ByStrata'][op0][vnam]['Ensemble Mean'][iT,iSP,iSS])

                d[vnam]=np.round(y,decimals=2)

            df0=pd.DataFrame().from_dict(d,orient='index')
            df=pd.concat([df,df0],axis=1)

    #df.index.name='Variable'
    df.columns=[np.arange(1,df.columns.size+1)]
    #df=df.sort_index(axis=0)

    df.to_excel(meta['Paths']['Project'] + '\\Outputs\\TabularSummaryDelta_' + tabNam + '.xlsx')

    return df

#%% Plot mean fluxes and mean pools over a specified time horizon

def PlotSchematicAtmoGHGBal(meta,mos,iB,iP,t_start,t_end,**kwargs):

    # Key word arguments
    if 'iSP' in kwargs.keys():
        iSP=kwargs['iSP']
    else:
        iSP=0

    if 'iSS' in kwargs.keys():
        iSS=kwargs['iSS']
    else:
        iSS=0

    tv=np.arange(meta['Project']['Year Start Saving'],meta['Project']['Year End']+1,1)

    iT=np.where( (tv>=t_start) & (tv<=t_end) )[0]

    y_b={}
    y_p={}
    y_d={}
    for k in mos['Scenarios'][0]['Mean'].keys():
        if (k=='C_Forest_Tot') | (k=='C_HWP_Tot') | (k=='C_ToMill_Tot'):
            y_b[k]=mos['Scenarios'][iB]['Mean'][k]['Ensemble Mean'][iT[-1],iSP,iSS]-mos['Scenarios'][iB]['Mean'][k]['Ensemble Mean'][iT[0],iSP,iSS]
            y_p[k]=mos['Scenarios'][iP]['Mean'][k]['Ensemble Mean'][iT[-1],iSP,iSS]-mos['Scenarios'][iP]['Mean'][k]['Ensemble Mean'][iT[0],iSP,iSS]
            y_d[k]=y_p[k]-y_b[k]
        else:
            y_b[k]=np.sum(mos['Scenarios'][iB]['Mean'][k]['Ensemble Mean'][iT,iSP,iSS])
            y_p[k]=np.sum(mos['Scenarios'][iP]['Mean'][k]['Ensemble Mean'][iT,iSP,iSS])
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
    bx_esc_w=0.15
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
    ax.add_patch(Rectangle([0.01,bx_atmo_bottom],0.98,0.1,fc=[0.9,0.95,1],ec=bx_ec))
    ax.text(0.5,0.935,'Atmosphere',size=bx_fs,ha='center',fontweight='bold',color=[0.08,0.3,0.55])
    vr='E_CO2e_AGHGB_WSub'
    a1=np.round(y_b[vr],decimals=decim); a2=np.round(y_p[vr],decimals=decim); a3=np.round(y_d[vr],decimals=decim)
    txt='Change in storage (tCO$_2$e/ha): ' + str(a1) + ',' + str(a2) + ' (' + GetSign(a3) + str(a3) + ')'
    ax.text(0.5,0.9,txt,size=fs_flux+1,ha='center')

    # LULUCF
    ax.add_patch(Rectangle([0.01,0.01],bx_lulucf_w,bx_lower_h,fc=[0.85,0.9,0.85],ec=bx_ec))
    ax.text(0.25,0.04,'Land Use, Land Use Change and Forestry',size=bx_fs,ha='center',color=[0,0.5,0])

    # Forest land
    ax.add_patch(Rectangle([0.02,0.1],bx_lulucf_w*0.53,bx_lower_h-0.11,fc=[0.9,0.95,0.9],ec=bx_ec))
    ax.text(0.15,0.29,'Forest Land',size=bx_fs,ha='center',color=[0,0.5,0])
    ax.text(0.15,0.26,'Change in storage (tC/ha):',size=fs_flux,ha='center')
    vr='C_Forest_Tot'
    a1=np.round(y_b[vr],decimals=decim); a2=np.round(y_p[vr],decimals=decim); a3=np.round(y_d[vr],decimals=decim)
    txt=str(a1) + ',' + str(a2) + ' (' + GetSign(a3) + str(a3) + ')'
    ax.text(0.15,0.23,txt,size=fs_flux,ha='center')

    # Harvested wood products
    ax.add_patch(Rectangle([0.36,0.1],0.12,bx_lower_h-0.11,fc=[0.9,0.95,0.9],ec=bx_ec))
    ax.text(0.42,0.27,'Harvested\nWood\nProducts',size=bx_fs,ha='center',color=[0,0.5,0])
    vr='C_HWP_Tot'
    a1=np.round(y_b[vr],decimals=decim); a2=np.round(y_p[vr],decimals=decim); a3=np.round(y_d[vr],decimals=decim)
    txt='Change in\nstorage (tC/ha):\n' + str(a1) + ',' + str(a2) + ' (' + GetSign(a3) + str(a3) + ')'
    ax.text(0.42,0.185,txt,size=fs_flux,ha='center')

    # Lithosphere
    ax.add_patch( Rectangle([bx_lulucf_w+0.02,0.01],bx_esc_w*3+0.03+0.01,bx_lower_h,fc=[0.94,0.88,0.84],ec=bx_ec))
    ax.text(0.75,0.04,'Lithosphere (Fossil Fuels & Limestone)',size=bx_fs,ha='center',color=[0.5,0,0])

    # Energy - Stationary Combustion
    ax.add_patch(Rectangle([bx_lulucf_w+0.02+0.01,0.1],bx_esc_w,bx_lower_h-0.11,fc=[1,0.95,0.9],ec=bx_ec))
    ax.text(0.585,0.29,'Stationary\nCombustion',size=bx_fs,ha='center',color=[0.5,0,0])
    a1=np.round(y_d['E_CO2e_SUB_ESC']/3.667,decimals=decim);
    a2=np.round(-1*y_d['E_CO2e_ESC_OperFor']/3.667,decimals=decim);
    a3=-1*np.round((y_d['E_CO2e_SUB_ESC']+y_d['E_CO2e_ESC_OperFor'])/3.667,decimals=decim)
    txt='Change in\nstorage (tC/ha):\n' + GetSign(a1) + str(a1) + ',' + GetSign(a2) + str(a2) + ' (' + GetSign(a3) + str(a3) + ')'
    ax.text(0.585,0.21,txt,size=fs_flux,ha='center')

    # Energy - Transportation
    ax.add_patch(Rectangle([bx_lulucf_w+0.02+0.01+bx_esc_w+0.01,0.1],bx_esc_w,bx_lower_h-0.11,fc=[1,0.95,0.9],ec=bx_ec))
    ax.text(0.745,0.29,'Transportation',size=bx_fs,ha='center',color=[0.5,0,0])
    a1=np.round(y_d['E_CO2e_SUB_ET']/3.667,decimals=decim);
    a2=np.round(-1*y_d['E_CO2e_ET_OperFor']/3.667,decimals=decim);
    a3=-1*np.round((y_d['E_CO2e_SUB_ET']+y_d['E_CO2e_ET_OperFor'])/3.667,decimals=decim)
    txt='Change in\n storage (tC/ha):\n' + GetSign(a1) + str(a1) + ',' + GetSign(a2) + str(a2) + ' (' + GetSign(a3) + str(a3) + ')'
    ax.text(0.745,0.21,txt,size=fs_flux,ha='center')

    # IPPU
    ax.add_patch(Rectangle([bx_lulucf_w+0.02+0.01+bx_esc_w+0.01+bx_esc_w+0.01,0.1],bx_esc_w,bx_lower_h-0.11,fc=[1,0.95,0.9],ec=bx_ec))
    ax.text(0.905,0.27,'Industrial\nProcesses &\nProduct Use',size=bx_fs,ha='center',color=[0.5,0,0])
    a1=np.round(y_d['E_CO2e_SUB_IPPU']/3.667,decimals=decim);
    a2=np.round(-1*y_d['E_CO2e_IPPU_OperFor']/3.667,decimals=decim);
    a3=-1*np.round((y_d['E_CO2e_SUB_IPPU']+y_d['E_CO2e_IPPU_OperFor'])/3.667,decimals=decim)
    txt='Change in\n storage (tC/ha):\n' + GetSign(a1) + str(a1) + ',' + GetSign(a2) + str(a2) + ' (' + GetSign(a3) + str(a3) + ')'
    ax.text(0.905,0.185,txt,size=fs_flux,ha='center')

    #--------------------------------------------------------------------------
    # Fluxes
    #--------------------------------------------------------------------------

    # NEE
    vr='E_CO2e_LULUCF_NEE'
    a1=np.round(y_b[vr],decimals=decim); a2=np.round(y_p[vr],decimals=decim); a3=np.round(y_d[vr],decimals=decim)
    txt='Net ecosystem\nexchange\n' + str(a1) + ',' + str(a2) + ' (' + GetSign(a3) + str(a3) + ')'
    ax.text(0.025,0.79,txt,ha='left',size=fs_flux)
    ax.arrow(0.02,bx_atmo_bottom,0,-1*(bx_atmo_bottom-bx_lower_h-0.02),head_width=arrow_head_w,head_length=0.01,fc='k',ec='k',lw=arrow_lw)

    # Wildfire
    vr='E_CO2e_LULUCF_Wildfire'
    a1=np.round(y_b[vr],decimals=decim); a2=np.round(y_p[vr],decimals=decim); a3=np.round(y_d[vr],decimals=decim)
    txt='Wildfire\n' + str(a1) + ',' + str(a2) + ' (' + str(a3) + ')'
    ax.text(0.125,0.52,txt,ha='right',size=fs_flux)
    ax.arrow(0.13,bx_lower_h+0.01,0,bx_atmo_bottom-bx_lower_h-0.02,head_width=arrow_head_w,head_length=0.01,fc='k',ec='k',lw=arrow_lw)

    # Open burning
    vr='E_CO2e_LULUCF_OpenBurning'
    a1=np.round(y_b[vr],decimals=decim); a2=np.round(y_p[vr],decimals=decim); a3=np.round(y_d[vr],decimals=decim)
    txt='Open burning\n' + str(a1) + ',' + str(a2) + ' (' + GetSign(a3) + str(a3) + ')'
    ax.text(0.155,0.52,txt,ha='left',size=fs_flux)
    ax.arrow(0.15,bx_lower_h+0.01,0,bx_atmo_bottom-bx_lower_h-0.02,head_width=arrow_head_w,head_length=0.01,fc='k',ec='k',lw=arrow_lw)

    # Denitrification
    vr='E_CO2e_LULUCF_Denit'
    a1=np.round(y_b[vr],decimals=decim); a2=np.round(y_p[vr],decimals=decim); a3=np.round(y_d[vr],decimals=decim)
    txt='Denitrification\n' + str(a1) + ',' + str(a2) + ' (' + str(a3) + ')'
    ax.text(0.245,0.71,txt,ha='right',size=fs_flux)
    ax.arrow(0.25,bx_lower_h+0.01,0,bx_atmo_bottom-bx_lower_h-0.02,head_width=arrow_head_w,head_length=0.01,fc='k',ec='k',lw=arrow_lw)

    # Volatilization
    vr='E_CO2e_LULUCF_Other'
    a1=np.round(y_b[vr],decimals=decim); a2=np.round(y_p[vr],decimals=decim); a3=np.round(y_d[vr],decimals=decim)
    txt='Volatilization\n' + str(a1) + ',' + str(a2) + ' (' + GetSign(a3) + str(a3) + ')'
    ax.text(0.275,0.71,txt,ha='left',size=fs_flux)
    ax.arrow(0.27,bx_lower_h+0.01,0,bx_atmo_bottom-bx_lower_h-0.02,head_width=arrow_head_w,head_length=0.01,fc='k',ec='k',lw=arrow_lw)

    # HWP fluxes
    vr='E_CO2e_LULUCF_HWP'
    a1=np.round(y_b[vr],decimals=decim); a2=np.round(y_p[vr],decimals=decim); a3=np.round(y_d[vr],decimals=decim)
    txt='Product decay and\ncombustion\n' + str(a1) + ',' + str(a2) + ' (' + GetSign(a3) + str(a3) + ')'
    ax.text(0.415,0.51,txt,ha='right',size=fs_flux)
    ax.arrow(0.42,bx_lower_h+0.01,0,bx_atmo_bottom-bx_lower_h-0.02,head_width=arrow_head_w,head_length=0.01,fc='k',ec='k',lw=arrow_lw)

    # Removals
    vr='C_ToMill'
    a1=np.round(y_b[vr],decimals=decim); a2=np.round(y_p[vr],decimals=decim); a3=np.round(y_d[vr],decimals=decim)
    #a1=0.2; a2=0.5; a3=0.3
    txt='Removals\n' + str(a1) + ',' + str(a2) + ' (' + GetSign(a3) + str(a3) + ')'
    ax.text(0.315,0.31,txt,ha='center',size=fs_flux)
    ax.arrow(0.275,0.275,0.08,0,head_width=0.01,head_length=arrow_head_w,fc='k',ec='k',lw=arrow_lw)

    # Bioenergy combustion
    vr='E_CO2e_ESC_Bioenergy'
    a1=np.round(y_b[vr],decimals=decim); a2=np.round(y_p[vr],decimals=decim); a3=np.round(y_d[vr],decimals=decim)
    txt='Bioenergy\ncombustion\n' + str(a1) + ',' + str(a2) + ' (' + GetSign(a3) + str(a3) + ')'
    ax.text(0.515,0.64,txt,ha='right',va='top',size=fs_flux)
    ax.arrow(0.52,bx_lower_h+0.01,0,bx_atmo_bottom-bx_lower_h-0.02,head_width=arrow_head_w,head_length=0.01,fc='k',ec='k',lw=arrow_lw)

    # ESC operational emissions
    vr='E_CO2e_ESC_OperFor'
    a1=np.round(y_b[vr],decimals=decim); a2=np.round(y_p[vr],decimals=decim); a3=np.round(y_d[vr],decimals=decim)
    txt='Fossil fuel for\nstationary\ncombustion\n' + str(a1) + ',' + str(a2) + ' (' + GetSign(a3) + str(a3) + ')'
    ax.text(0.55,0.64,txt,ha='left',va='top',size=fs_flux)
    ax.arrow(0.545,bx_lower_h+0.01,0,bx_atmo_bottom-bx_lower_h-0.02,head_width=arrow_head_w,head_length=0.01,fc='k',ec='k',lw=arrow_lw)

    # Displacement energy
    vr='E_CO2e_SUB_E'
    a1=np.round(y_b[vr],decimals=decim); a2=np.round(y_p[vr],decimals=decim); a3=np.round(y_d[vr],decimals=decim)
    txt='Displacement \neffects of\nbioenergy\n' + str(a1) + ',' + str(a2) + ' (' + GetSign(a3) + str(a3) + ')'
    ax.text(0.655,0.82,txt,ha='right',va='top',size=fs_flux)
    ax.arrow(0.66,bx_atmo_bottom,0,-1*(bx_atmo_bottom-bx_lower_h-0.02),head_width=arrow_head_w,head_length=0.01,fc='k',ec='k',lw=arrow_lw)

    # Transportation
    vr='E_CO2e_ET_OperFor'
    a1=np.round(y_b[vr],decimals=decim); a2=np.round(y_p[vr],decimals=decim); a3=np.round(y_d[vr],decimals=decim)
    txt='Fossil fuel for\ntransportation\n' + str(a1) + ',' + str(a2) + ' (' + GetSign(a3) + str(a3) + ')'
    ax.text(0.695,0.64,txt,ha='left',va='top',size=fs_flux)
    ax.arrow(0.69,bx_lower_h+0.01,0,bx_atmo_bottom-bx_lower_h-0.02,head_width=arrow_head_w,head_length=0.01,fc='k',ec='k',lw=arrow_lw)

    # Displacement building materials
    vr='E_CO2e_SUB_M'
    a1=np.round(y_b[vr],decimals=decim); a2=np.round(y_p[vr],decimals=decim); a3=np.round(y_d[vr],decimals=decim)
    txt='Displacement \neffects of\nsolid wood\nmaterials\n' + str(a1) + ',' + str(a2) + ' (' + GetSign(a3) + str(a3) + ')'
    ax.text(0.805,0.82,txt,ha='left',va='top',size=fs_flux)
    ax.arrow(0.8,bx_atmo_bottom,0,-1*(bx_atmo_bottom-bx_lower_h-0.02),head_width=arrow_head_w,head_length=0.01,fc='k',ec='k',lw=arrow_lw)

    # IPPU
    vr='E_CO2e_IPPU_OperFor'
    a1=np.round(y_b[vr],decimals=decim); a2=np.round(y_p[vr],decimals=decim); a3=np.round(y_d[vr],decimals=decim)
    txt='Fossil fuel\ncombustion\nand\nurea\nsequestration\n' + str(a1) + ',' + str(a2) + ' (' + GetSign(a3) + str(a3) + ')'
    ax.text(0.915,0.64,txt,ha='left',va='top',size=fs_flux)
    ax.arrow(0.91,bx_lower_h+0.01,0,bx_atmo_bottom-bx_lower_h-0.02,head_width=arrow_head_w,head_length=0.01,fc='k',ec='k',lw=arrow_lw)

    ax.set(position=[0,0,1,1],visible='Off',xticks=[],yticks=[])

    try:
        if meta['Print Figures']=='On':
            gu.PrintFig(meta['Paths']['Figures'] + '\\AGHGB Schematic_S' + str(iP) + 'minusS' + str(iB) + '_' + str(t_start) + 'to' + str(t_end),'png',900)
    except:
        gu.PrintFig(meta['Paths']['Figures'] + '\\AGHGB Schematic_S' + str(iP) + 'minusS' + str(iB) + '_' + str(t_start) + 'to' + str(t_end),'png',900)

    return

#%% Plot Cashflow

def PlotCashflow(meta,mos,iB,iP,t_start,t_end):

    tv=np.arange(meta['Project']['Year Start Saving'],meta['Project']['Year End']+1,1)
    iT=np.where( (tv>=t_start) & (tv<=t_end) )[0]

    fig,ax=plt.subplots(3,2,figsize=gu.cm2inch(15,11.5)); lw=0.75
    ax[0,0].plot(tv,mos[iB]['Cashflow']['Mean']['Cost Total']['Ensemble Mean']/1000,'-bo',lw=lw,ms=4,label='Baseline')
    ax[0,0].plot(tv,mos[iP]['Cashflow']['Mean']['Cost Total']['Ensemble Mean']/1000,'--r^',lw=lw,ms=3,label='Project')
    ax[0,0].set(xlim=[tv[iT[0]], tv[iT[-1]]],ylabel='Cost (CDN$/000)')
    ax[0,0].legend(loc='upper right',frameon=False,facecolor=None)

    ax[0,1].plot(tv,(mos[iP]['Cashflow']['Mean']['Cost Total']['Ensemble Mean']-mos[iB]['Cashflow']['Mean']['Cost Total']['Ensemble Mean'])/1000,'-g',lw=lw)
    ax[0,1].set(xlim=[tv[iT[0]], tv[iT[-1]]],ylabel='Cost (CDN$/000)')

    ax[1,0].plot(tv,mos[iB]['Cashflow']['Mean']['Revenue Net']['Ensemble Mean']/1000,'-bo',ms=4,lw=lw)
    ax[1,0].plot(tv,mos[iP]['Cashflow']['Mean']['Revenue Net']['Ensemble Mean']/1000,'--r^',ms=3,lw=lw)
    ax[1,0].set(xlim=[tv[iT[0]], tv[iT[-1]]],ylabel='Net revenue (CDN$/000)')

    ax[1,1].plot(tv,(mos[iP]['Cashflow']['Mean']['Revenue Net']['Ensemble Mean']-mos[iB]['Cashflow']['Mean']['Revenue Net']['Ensemble Mean'])/1000,'-g',lw=lw)
    ax[1,1].set(xlim=[tv[iT[0]], tv[iT[-1]]],ylabel='Net revenue (CDN$/000)')

    ax[2,0].plot(tv,mos[iB]['Cashflow']['Mean']['Revenue Net_cumu']['Ensemble Mean']/1000,'-b',ms=4,lw=lw)
    ax[2,0].plot(tv,mos[iP]['Cashflow']['Mean']['Revenue Net_cumu']['Ensemble Mean']/1000,'--r',ms=3,lw=lw)
    ax[2,0].set(xlim=[tv[iT[0]], tv[iT[-1]]],ylabel='Cumulative net revenue (CDN$/000)')

    ax[2,1].plot(tv,(mos[iP]['Cashflow']['Mean']['Revenue Net_cumu']['Ensemble Mean']-mos[iB]['Cashflow']['Mean']['Revenue Net_cumu']['Ensemble Mean'])/1000,'-g',lw=lw)
    ax[2,1].set(xlim=[tv[iT[0]], tv[iT[-1]]],ylabel='Cumulative net revenue (CDN$/000)');

    if meta['Print Figures']=='On':
        gu.PrintFig(meta['Paths']['Figures'] + '\\Cashflow_S' + str(iP) + 'minusS' + str(iB) + '_' + str(t_start) + 'to' + str(t_end),'png',900)

    return

#%%

def PlotPools(meta,mos,tv,iT,**kwargs):

    iSP=0
    iSS=0

    vs=['C_Biomass_Tot','C_DeadWood_Tot','C_Litter_Tot','C_Soil_Tot','C_InUse_Tot','C_DumpLandfill_Tot']
    vs2=['Biomass','Dead Wood','Litter','Soil','In-use Products','Dump and Landfill']

    cl=np.array([[0,0.5,1],[0,0.6,0],[0.5,0,1],[0,1,1]])
    symb=['-','--','-.',':','-']

    if 'Custom Scenario List' not in kwargs.keys():

        # Generate one figure per scenario comparison

        for k in mos['Delta'].keys():
            cnt=0
            fig,ax=plt.subplots(3,2,figsize=gu.cm2inch(18,15)); Alpha=0.09
            sL=[mos['Delta'][k]['iB'],mos['Delta'][k]['iP']]
            for i in range(3):
                for j in range(2):

                    for iScn in range(len(sL)):
                        be=mos['Scenarios'][sL[iScn]]['Mean'][vs[cnt]]['Ensemble Mean'][iT,iSP,iSS]
                        lo=mos['Scenarios'][sL[iScn]]['Mean'][vs[cnt]]['Ensemble P025'][iT,iSP,iSS]
                        hi=mos['Scenarios'][sL[iScn]]['Mean'][vs[cnt]]['Ensemble P975'][iT,iSP,iSS]
                        lo2=mos['Scenarios'][sL[iScn]]['Mean'][vs[cnt]]['Ensemble P250'][iT,iSP,iSS]
                        hi2=mos['Scenarios'][sL[iScn]]['Mean'][vs[cnt]]['Ensemble P750'][iT,iSP,iSS]

                        ax[i,j].fill_between(tv[iT],lo,hi,color=cl[iScn,:],alpha=Alpha,linewidth=0)
                        ax[i,j].fill_between(tv[iT],lo2,hi2,color=cl[iScn,:],alpha=Alpha,linewidth=0)
                        ax[i,j].plot(tv[iT],be,symb[iScn],color=cl[iScn,:],lw=1,label='Scenario ' + str(iScn+1))

                    if (i==0) & (j==0):
                        ax[i,j].legend(loc="lower left")
                    ax[i,j].set(ylabel=vs2[cnt] + ' (MgC/ha)')
                    cnt=cnt+1

            gu.axletters(ax,plt,0.035,0.9)

            if meta['Print Figures']=='On':
                gu.PrintFig(meta['Paths']['Figures'] + '\\Pools_' + k,'png',900)

    else:

        # Generate one figure for a custom set of scenarios

        cnt=0
        fig,ax=plt.subplots(3,2,figsize=gu.cm2inch(18,15)); Alpha=0.09

        for iScn in range(len(kwargs['Custom Scenario List'])):
            s=kwargs['Custom Scenario List'][iScn]
            for i in range(3):
                for j in range(2):
                    be=mos['Scenarios'][s]['Mean'][vs[cnt]]['Ensemble Mean'][iT,iSP,iSS]
                    lo=mos['Scenarios'][s]['Mean'][vs[cnt]]['Ensemble P025'][iT,iSP,iSS]
                    hi=mos['Scenarios'][s]['Mean'][vs[cnt]]['Ensemble P975'][iT,iSP,iSS]
                    lo2=mos['Scenarios'][s]['Mean'][vs[cnt]]['Ensemble P250'][iT,iSP,iSS]
                    hi2=mos['Scenarios'][s]['Mean'][vs[cnt]]['Ensemble P750'][iT,iSP,iSS]

                    ax[i,j].fill_between(tv[iT],lo,hi,color=cl[iScn,:],alpha=Alpha,linewidth=0)
                    ax[i,j].fill_between(tv[iT],lo2,hi2,color=cl[iScn,:],alpha=Alpha,linewidth=0)
                    ax[i,j].plot(tv[iT],be,symb[iScn],color=cl[iScn,:],lw=1,label='Scenario ' + str(iScn+1))

                    if (i==0) & (j==0):
                        ax[i,j].legend(loc="lower left")
                    ax[i,j].set(ylabel=vs2[cnt] + ' (MgC/ha)')
            cnt=cnt+1

        gu.axletters(ax,plt,0.035,0.9)

        if meta['Print Figures']=='On':
            gu.PrintFig(meta['Paths']['Figures'] + '\\Pools_CustomScenarioList','png',900)

    return

#%%

def PlotFluxes(meta,mos,tv,iT,**kwargs):

    iSP=0
    iSS=0

    vs=['C_NPP_Tot','C_G_Net_Tot','C_RH_Tot','E_CO2e_LULUCF_OpenBurning','E_CO2e_LULUCF_Wildfire','E_CO2e_LULUCF_HWP','E_CO2e_SUB_Tot','E_CO2e_AGHGB_WSub']
    vs2=['NPP (tCO2e/ha/yr)','Net growth (tCO2e/ha/yr)','RH (tCO2e/ha/yr)','Open burning (tCO2e/ha/yr)','Wildfire (tCO2e/ha/yr)','HWP (tCO2e/ha/yr)','Substitutions (tCO2e/ha/yr)','GHG balance (tCO2e/ha/yr)']

    cl=np.array([[0,0.5,1],[0,0.6,0],[0.5,0,1],[0,1,1]])
    symb=['-','--','-.',':','-']

    if 'Custom Scenario List' not in kwargs.keys():

        # Generate one figure per scenario comparison

        for k in mos['Delta'].keys():
            cnt=0
            fig,ax=plt.subplots(3,2,figsize=gu.cm2inch(18,15)); Alpha=0.09
            sL=[mos['Delta'][k]['iB'],mos['Delta'][k]['iP']]
            for i in range(3):
                for j in range(2):

                    for iScn in range(len(sL)):
                        be=mos['Scenarios'][sL[iScn]]['Mean'][vs[cnt]]['Ensemble Mean'][iT,iSP,iSS]
                        lo=mos['Scenarios'][sL[iScn]]['Mean'][vs[cnt]]['Ensemble P025'][iT,iSP,iSS]
                        hi=mos['Scenarios'][sL[iScn]]['Mean'][vs[cnt]]['Ensemble P975'][iT,iSP,iSS]
                        lo2=mos['Scenarios'][sL[iScn]]['Mean'][vs[cnt]]['Ensemble P250'][iT,iSP,iSS]
                        hi2=mos['Scenarios'][sL[iScn]]['Mean'][vs[cnt]]['Ensemble P750'][iT,iSP,iSS]

                        ax[i,j].fill_between(tv[iT],lo,hi,color=cl[iScn,:],alpha=Alpha,linewidth=0)
                        ax[i,j].fill_between(tv[iT],lo2,hi2,color=cl[iScn,:],alpha=Alpha,linewidth=0)
                        ax[i,j].plot(tv[iT],be,symb[iScn],color=cl[iScn,:],lw=1,label='Scenario ' + str(iScn+1))

                    if (i==0) & (j==0):
                        ax[i,j].legend(loc="lower left")
                    ax[i,j].set(ylabel=vs2[cnt] + ' (MgC/ha)')
                    cnt=cnt+1

            gu.axletters(ax,plt,0.035,0.9)

            if meta['Print Figures']=='On':
                gu.PrintFig(meta['Paths']['Figures'] + '\\Fluxes_' + k,'png',900)

    else:

        # Generate one figure for a custom set of scenarios

        cnt=0
        fig,ax=plt.subplots(3,2,figsize=gu.cm2inch(18,15)); Alpha=0.09

        for iScn in range(len(kwargs['Custom Scenario List'])):
            s=kwargs['Custom Scenario List'][iScn]
            for i in range(3):
                for j in range(2):
                    be=mos['Scenarios'][s]['Mean'][vs[cnt]]['Ensemble Mean'][iT,iSP,iSS]
                    lo=mos['Scenarios'][s]['Mean'][vs[cnt]]['Ensemble P025'][iT,iSP,iSS]
                    hi=mos['Scenarios'][s]['Mean'][vs[cnt]]['Ensemble P975'][iT,iSP,iSS]
                    lo2=mos['Scenarios'][s]['Mean'][vs[cnt]]['Ensemble P250'][iT,iSP,iSS]
                    hi2=mos['Scenarios'][s]['Mean'][vs[cnt]]['Ensemble P750'][iT,iSP,iSS]

                    ax[i,j].fill_between(tv[iT],lo,hi,color=cl[iScn,:],alpha=Alpha,linewidth=0)
                    ax[i,j].fill_between(tv[iT],lo2,hi2,color=cl[iScn,:],alpha=Alpha,linewidth=0)
                    ax[i,j].plot(tv[iT],be,symb[iScn],color=cl[iScn,:],lw=1,label='Scenario ' + str(iScn+1))

                    if (i==0) & (j==0):
                        ax[i,j].legend(loc="lower left")
                    ax[i,j].set(ylabel=vs2[cnt])
            cnt=cnt+1

        gu.axletters(ax,plt,0.035,0.9)

        if meta['Print Figures']=='On':
            gu.PrintFig(meta['Paths']['Figures'] + '\\Fluxes_CustomScenarioList','png',900)

    return

#%% Plot NEP, RH, and NPP

def PlotNEP(meta,mos,tv,iT):

    iSP=0
    iSS=0

    for k in mos['Delta'].keys():

        iB=mos['Delta'][k]['iB']
        iP=mos['Delta'][k]['iP']

        fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(18,8)); Alpha=0.09

        for j in range(0,2):
            ax[j].yaxis.set_ticks_position('both');
            ax[j].xaxis.set_ticks_position('both')
            ax[j].plot(tv[iT],0*np.ones(tv[iT].shape),'-',lw=3,color=(0.8,0.8,0.8),label='')

        cl=np.array([[0.75,0,0],[0,0.85,0],[0.29,0.49,0.77]])
        ymin=0.0
        ymax=0.0

        vn='C_RH_Tot'
        lo=3.667*mos['Delta'][k]['ByStrata']['Mean'][vn]['Ensemble P025'][iT,iSP,iSS]
        hi=3.667*mos['Delta'][k]['ByStrata']['Mean'][vn]['Ensemble P975'][iT,iSP,iSS]
        lo2=3.667*mos['Delta'][k]['ByStrata']['Mean'][vn]['Ensemble P250'][iT,iSP,iSS]
        hi2=3.667*mos['Delta'][k]['ByStrata']['Mean'][vn]['Ensemble P750'][iT,iSP,iSS]
        #ax[0].fill_between(tv[iT],lo,hi,color=cl[0,:],alpha=Alpha,linewidth=0)
        ax[0].fill_between(tv[iT],lo2,hi2,color=cl[0,:],alpha=Alpha,linewidth=0)
        ax[0].plot(tv[iT],3.667*mos['Delta'][k]['ByStrata']['Mean'][vn]['Ensemble Mean'][iT,iSP,iSS],'--',color=cl[0,:],label='RH')
        ymin=np.minimum(ymin,np.min(lo))
        ymax=np.maximum(ymax,np.max(hi))

        vn='C_NPP_Tot'
        lo=3.667*mos['Delta'][k]['ByStrata']['Mean'][vn]['Ensemble P025'][iT,iSP,iSS]
        hi=3.667*mos['Delta'][k]['ByStrata']['Mean'][vn]['Ensemble P975'][iT,iSP,iSS]
        lo2=3.667*mos['Delta'][k]['ByStrata']['Mean'][vn]['Ensemble P250'][iT,iSP,iSS]
        hi2=3.667*mos['Delta'][k]['ByStrata']['Mean'][vn]['Ensemble P750'][iT,iSP,iSS]
        #ax[0].fill_between(tv[iT],lo,hi,color=cl[0,:],alpha=Alpha,linewidth=0)
        ax[0].fill_between(tv[iT],lo2,hi2,color=cl[1,:],alpha=Alpha,linewidth=0)
        ax[0].plot(tv[iT],3.667*mos['Delta'][k]['ByStrata']['Mean'][vn]['Ensemble Mean'][iT,iSP,iSS],'-.',color=cl[1,:],label='NPP')
        ymin=np.minimum(ymin,np.min(lo))
        ymax=np.maximum(ymax,np.max(hi))

        vn='E_CO2e_LULUCF_NEE'
        lo=-mos['Delta'][k]['ByStrata']['Mean'][vn]['Ensemble P025'][iT,iSP,iSS]
        hi=-mos['Delta'][k]['ByStrata']['Mean'][vn]['Ensemble P975'][iT,iSP,iSS]
        lo2=-mos['Delta'][k]['ByStrata']['Mean'][vn]['Ensemble P250'][iT,iSP,iSS]
        hi2=-mos['Delta'][k]['ByStrata']['Mean'][vn]['Ensemble P750'][iT,iSP,iSS]
        #ax[0].fill_between(tv[iT],lo,hi,color=cl[0,:],alpha=Alpha,linewidth=0)
        ax[0].fill_between(tv[iT],lo2,hi2,color=cl[2,:],alpha=Alpha,linewidth=0)
        ax[0].plot(tv[iT],-mos['Delta'][k]['ByStrata']['Mean'][vn]['Ensemble Mean'][iT,iSP,iSS],'-',color=cl[2,:],label='NEP')
        ymin=np.minimum(ymin,np.min(lo))
        ymax=np.maximum(ymax,np.max(hi))

        ax[0].legend(loc="lower right",frameon=0)
        ax[0].set(ylabel='Annual $\Delta$ (tCO$_2$e ha$^-$$^1$ yr$^-$$^1$)',xlabel='Time, years',ylim=[ymin,ymax],xlim=[tv[iT[0]],tv[iT[-1]]]);

        ymin=0.0
        ymax=0.0

        vn='C_RH_Tot'
        lo=np.cumsum(3.667*mos['Delta'][k]['ByStrata']['Mean'][vn]['Ensemble P025'][iT,iSP,iSS])
        hi=np.cumsum(3.667*mos['Delta'][k]['ByStrata']['Mean'][vn]['Ensemble P975'][iT,iSP,iSS])
        lo2=np.cumsum(3.667*mos['Delta'][k]['ByStrata']['Mean'][vn]['Ensemble P250'][iT,iSP,iSS])
        hi2=np.cumsum(3.667*mos['Delta'][k]['ByStrata']['Mean'][vn]['Ensemble P750'][iT,iSP,iSS])
        mu=np.cumsum(3.667*mos['Delta'][k]['ByStrata']['Mean'][vn]['Ensemble Mean'][iT,iSP,iSS])
        #ax[1].fill_between(tv[iT],lo,hi,color=[0.75,.5,1],alpha=Alpha,linewidth=0,label='95 C.I.')
        ax[1].fill_between(tv[iT],lo2,hi2,color=cl[0,:],alpha=Alpha,linewidth=0)
        ax[1].plot(tv[iT],mu,'--',color=cl[0,:])
        ax[1].set(ylabel='Cumulative $\Delta$ (tCO$_2$e ha$^-$$^1$)',xlabel='Time, years',xlim=[tv[iT[0]],tv[iT[-1]]]);
        ymin=np.minimum(ymin,np.min(mu))
        ymax=np.maximum(ymax,np.max(mu))

        vn='C_NPP_Tot'
        lo=np.cumsum(3.667*mos['Delta'][k]['ByStrata']['Mean'][vn]['Ensemble P025'][iT,iSP,iSS])
        hi=np.cumsum(3.667*mos['Delta'][k]['ByStrata']['Mean'][vn]['Ensemble P975'][iT,iSP,iSS])
        lo2=np.cumsum(3.667*mos['Delta'][k]['ByStrata']['Mean'][vn]['Ensemble P250'][iT,iSP,iSS])
        hi2=np.cumsum(3.667*mos['Delta'][k]['ByStrata']['Mean'][vn]['Ensemble P750'][iT,iSP,iSS])
        mu=np.cumsum(3.667*mos['Delta'][k]['ByStrata']['Mean'][vn]['Ensemble Mean'][iT,iSP,iSS])
        #ax[1].fill_between(tv[iT],lo,hi,color=[0.75,.5,1],alpha=Alpha,linewidth=0,label='95 C.I.')
        ax[1].fill_between(tv[iT],lo2,hi2,color=cl[1,:],alpha=Alpha,linewidth=0)
        ax[1].plot(tv[iT],mu,'-.',color=cl[1,:])
        ymin=np.minimum(ymin,np.min(mu))
        ymax=np.maximum(ymax,np.max(mu))

        vn='E_CO2e_LULUCF_NEE'
        lo=-np.cumsum(mos['Delta'][k]['ByStrata']['Mean'][vn]['Ensemble P025'][iT,iSP,iSS])
        hi=-np.cumsum(mos['Delta'][k]['ByStrata']['Mean'][vn]['Ensemble P975'][iT,iSP,iSS])
        lo2=-np.cumsum(mos['Delta'][k]['ByStrata']['Mean'][vn]['Ensemble P250'][iT,iSP,iSS])
        hi2=-np.cumsum(mos['Delta'][k]['ByStrata']['Mean'][vn]['Ensemble P750'][iT,iSP,iSS])
        mu=-np.cumsum(mos['Delta'][k]['ByStrata']['Mean'][vn]['Ensemble Mean'][iT,iSP,iSS])
        #ax[1].fill_between(tv[iT],lo,hi,color=[0.75,.5,1],alpha=Alpha,linewidth=0,label='95 C.I.')
        ax[1].fill_between(tv[iT],lo2,hi2,color=cl[2,:],alpha=Alpha,linewidth=0)
        ax[1].plot(tv[iT],mu,'-',color=cl[2,:])
        ymin=np.minimum(ymin,np.min(mu))
        ymax=np.maximum(ymax,np.max(mu))

        ax[1].set(ylabel='Cumulative $\Delta$ (tCO$_2$e ha$^-$$^1$)',xlabel='Time, years',ylim=[ymin,ymax],xlim=[tv[iT[0]],tv[iT[-1]]]);
        gu.axletters(ax,plt,0.03,0.89)
        if meta['Print Figures']=='On':
            gu.PrintFig(meta['Paths']['Figures'] + '\\NEE_Balance_' + k,'png',900)
        fig.suptitle(k)

    return

#%%

def PlotGHGB(meta,mos,tv,iT):

    iSP=0
    iSS=0

    for k in mos['Delta'].keys():

        iB=mos['Delta'][k]['iB']
        iP=mos['Delta'][k]['iP']

        fig,ax=plt.subplots(2,2,figsize=gu.cm2inch(18,10)); Alpha=0.09

        for i in range(0,2):
            for j in range(0,2):
                ax[i,j].yaxis.set_ticks_position('both'); ax[i,j].xaxis.set_ticks_position('both')
                ax[i,j].plot(tv[iT],0*np.ones(tv[iT].shape),'-',lw=3,color=(0.8,0.8,0.8),label='')

        ax[0,0].plot(tv[iT],mos['Scenarios'][iB]['Mean']['E_CO2e_AGHGB_WSub']['Ensemble Mean'][iT,iSP,iSS],'-',color=(0,0.5,1),label='Baseline')
        ax[0,0].plot(tv[iT],mos['Scenarios'][iP]['Mean']['E_CO2e_AGHGB_WSub']['Ensemble Mean'][iT,iSP,iSS],'--',color=(0,0.6,0),label='Project (With Subs.)')
        ax[0,0].legend(loc="upper right",frameon=0)
        ax[0,0].set(ylabel='AGHGB (tCO$_2$e ha$^-$$^1$ yr$^-$$^1$)',xlim=[tv[iT[0]],tv[iT[-1]]],xlabel='Time, years');

        ax[0,1].plot(tv[iT],mos['Scenarios'][iB]['Mean']['E_CO2e_AGHGB_WSub_cumu']['Ensemble Mean'][iT,iSP,iSS],'-',color=(0,0.5,1),label='Baseline SR')
        ax[0,1].plot(tv[iT],mos['Scenarios'][iP]['Mean']['E_CO2e_AGHGB_WSub_cumu']['Ensemble Mean'][iT,iSP,iSS],'--',color=(0,0.6,0),label='Baseline NSR')
        ax[0,1].set(ylabel='Cumulative AGHGB (tCO$_2$e ha$^-$$^1$)',xlim=[tv[iT[0]],tv[iT[-1]]],xlabel='Time, years');

        lo=mos['Delta'][k]['ByStrata']['Mean']['E_CO2e_AGHGB_WSub']['Ensemble P025'][iT,iSP,iSS]
        hi=mos['Delta'][k]['ByStrata']['Mean']['E_CO2e_AGHGB_WSub']['Ensemble P975'][iT,iSP,iSS]
        lo2=mos['Delta'][k]['ByStrata']['Mean']['E_CO2e_AGHGB_WSub']['Ensemble P250'][iT,iSP,iSS]
        hi2=mos['Delta'][k]['ByStrata']['Mean']['E_CO2e_AGHGB_WSub']['Ensemble P750'][iT,iSP,iSS]
        ax[1,0].fill_between(tv[iT],lo,hi,color=[0.15,0,0.75],alpha=Alpha,linewidth=0,label='95 C.I.')
        ax[1,0].fill_between(tv[iT],lo2,hi2,color=[0.05,0,0.6],alpha=Alpha,linewidth=0,label='50 C.I.')
        ax[1,0].plot(tv[iT],mos['Delta'][k]['ByStrata']['Mean']['E_CO2e_AGHGB_WSub']['Ensemble Mean'][iT,iSP,iSS],'-',color=(0.5,0,1),label='Best estimate (With Subs.)')
        ax[1,0].legend(loc="upper right",frameon=0)
        ax[1,0].set(ylabel='$\Delta$ AGHGB (tCO$_2$e ha$^-$$^1$ yr$^-$$^1$)',xlim=[tv[iT[0]],tv[iT[-1]]],xlabel='Time, years',ylim=[np.min(lo),np.max(hi)]);

        lo=mos['Delta'][k]['ByStrata']['Mean']['E_CO2e_AGHGB_WSub_cumu_from_tref']['Ensemble P025'][iT,iSP,iSS]
        hi=mos['Delta'][k]['ByStrata']['Mean']['E_CO2e_AGHGB_WSub_cumu_from_tref']['Ensemble P975'][iT,iSP,iSS]
        lo2=mos['Delta'][k]['ByStrata']['Mean']['E_CO2e_AGHGB_WSub_cumu_from_tref']['Ensemble P250'][iT,iSP,iSS]
        hi2=mos['Delta'][k]['ByStrata']['Mean']['E_CO2e_AGHGB_WSub_cumu_from_tref']['Ensemble P750'][iT,iSP,iSS]
        ax[1,1].fill_between(tv[iT],lo,hi,color=[0.75,.5,1],alpha=Alpha,linewidth=0,label='95 C.I.')
        ax[1,1].fill_between(tv[iT],lo2,hi2,color=[0.05,0,0.6],alpha=Alpha,linewidth=0,label='50 C.I.')
        ax[1,1].plot(tv[iT],mos['Delta'][k]['ByStrata']['Mean']['E_CO2e_AGHGB_WSub_cumu_from_tref']['Ensemble Mean'][iT,iSP,iSS],'-',color=(0.5,0,1),label='Best estimate (With Subs.)')
        ax[1,1].plot(tv[iT],mos['Delta'][k]['ByStrata']['Mean']['E_CO2e_AGHGB_WOSub_cumu_from_tref']['Ensemble Mean'][iT,iSP,iSS],'--',color=(0.7,0.2,1),label='Best estimate (W/O Subs.)')
        ax[1,1].set(ylabel='Cumulative $\Delta$ AGHGB (tCO$_2$e ha$^-$$^1$)',xlim=[tv[iT[0]],tv[iT[-1]]],xlabel='Time, years');

        gu.axletters(ax,plt,0.03,0.89)

        if meta['Print Figures']=='On':
            gu.PrintFig(meta['Paths']['Figures'] + '\\GHG_Balance_' + k,'png',900)

        fig.suptitle(k)

    return

#%%

def PlotGHGBenefit(meta,mos,tv,iT):

    iSP=0
    iSS=0

    cl=np.array([[0,0.5,1],[0,0.6,0],[0,1,1],[0.5,0,1]])
    symb=['-','--','-.',':','-']

    fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(16,6)); Alpha=0.09

    for i in range(0,2):
        ax[i].yaxis.set_ticks_position('both'); ax[i].xaxis.set_ticks_position('both')
        ax[i].plot(tv[iT],0*np.ones(tv[iT].shape),'-',lw=3,color=(0.8,0.8,0.8),label='')

    cnt=0
    for k in mos['Delta'].keys():
        lo=mos['Delta'][k]['ByStrata']['Mean']['E_CO2e_AGHGB_WSub']['Ensemble P025'][iT,iSP,iSS]
        hi=mos['Delta'][k]['ByStrata']['Mean']['E_CO2e_AGHGB_WSub']['Ensemble P975'][iT,iSP,iSS]
        lo2=mos['Delta'][k]['ByStrata']['Mean']['E_CO2e_AGHGB_WSub']['Ensemble P250'][iT,iSP,iSS]
        hi2=mos['Delta'][k]['ByStrata']['Mean']['E_CO2e_AGHGB_WSub']['Ensemble P750'][iT,iSP,iSS]
        ax[0].fill_between(tv[iT],lo,hi,color=cl[cnt,:],alpha=Alpha,linewidth=0)
        ax[0].fill_between(tv[iT],lo2,hi2,color=cl[cnt,:],alpha=Alpha,linewidth=0)
        ax[0].plot(tv[iT],mos['Delta'][k]['ByStrata']['Mean']['E_CO2e_AGHGB_WSub']['Ensemble Mean'][iT,iSP,iSS],symb[cnt],color=cl[cnt,:],label='SC ' + str(cnt+1) )
        ax[0].legend(loc="upper right",frameon=0)
        ax[0].set(ylabel='$\Delta$GHG (tCO$_2$e ha$^-$$^1$ yr$^-$$^1$)',xlim=[tv[iT[0]],tv[iT[-1]]], \
                  xlabel='Time, years',ylim=[np.min(lo),np.max(hi)]);

        lo=mos['Delta'][k]['ByStrata']['Mean']['E_CO2e_AGHGB_WSub_cumu_from_tref']['Ensemble P025'][iT,iSP,iSS]
        hi=mos['Delta'][k]['ByStrata']['Mean']['E_CO2e_AGHGB_WSub_cumu_from_tref']['Ensemble P975'][iT,iSP,iSS]
        lo2=mos['Delta'][k]['ByStrata']['Mean']['E_CO2e_AGHGB_WSub_cumu_from_tref']['Ensemble P250'][iT,iSP,iSS]
        hi2=mos['Delta'][k]['ByStrata']['Mean']['E_CO2e_AGHGB_WSub_cumu_from_tref']['Ensemble P750'][iT,iSP,iSS]
        ax[1].fill_between(tv[iT],lo,hi,color=cl[cnt,:],alpha=Alpha,linewidth=0)
        ax[1].fill_between(tv[iT],lo2,hi2,color=cl[cnt,:],alpha=Alpha,linewidth=0)
        ax[1].plot(tv[iT],mos['Delta'][k]['ByStrata']['Mean']['E_CO2e_AGHGB_WSub_cumu_from_tref']['Ensemble Mean'][iT,iSP,iSS],symb[cnt],color=cl[cnt,:],label='SC ' + str(cnt+1))
        ax[1].set(ylabel='Cumulative $\Delta$GHG (tCO$_2$e ha$^-$$^1$)',xlim=[tv[iT[0]],tv[iT[-1]]],xlabel='Time, years');
        cnt=cnt+1

    gu.axletters(ax,plt,0.035,0.92)

    if meta['Print Figures']=='On':
        gu.PrintFig(meta['Paths']['Figures'] + '\\GHG_Benefit','png',900)

    return

#%% Summary description

def SummaryDescription(meta,mos,sc,th):

    tv=np.arange(meta['Project']['Year Start Saving'],meta['Project']['Year End']+1,1)
    iT=np.where( (tv>=meta['Project']['Year Project']) & (tv<=meta['Project']['Year Project']+th) )[0]

    d={}
    d['EI PelletExport (tCO2e/ODT)']=np.sum(mos['Delta'][sc]['ByStrata']['Mean']['E_CO2e_ESC_BioenergyPelletExport']['Ensemble Mean'][iT,0,0])/np.sum(mos['Delta'][sc]['ByStrata']['Mean']['ODT PelletExport']['Ensemble Mean'][iT,0,0])
    d['EI PelletExport Boiler (tCO2e/GJ)']=np.sum(mos['Delta'][sc]['ByStrata']['Mean']['E_CO2e_ESC_BioenergyPelletExport']['Ensemble Mean'][iT,0,0])/np.sum(mos['Delta'][sc]['ByStrata']['Mean']['GJ PelletExport']['Ensemble Mean'][iT,0,0])
    d['EI PelletExport Boiler+Ops (tCO2e/GJ)']=np.sum( (mos['Delta'][sc]['ByStrata']['Mean']['E_CO2e_ESC_BioenergyPelletExport']['Ensemble Mean'][iT,0,0]+ \
                                                       mos['Delta'][sc]['ByStrata']['Mean']['E_CO2e_ESC_OperFor']['Ensemble Mean'][iT,0,0]+ \
                                                       mos['Delta'][sc]['ByStrata']['Mean']['E_CO2e_ET_OperFor']['Ensemble Mean'][iT,0,0]+ \
                                                       mos['Delta'][sc]['ByStrata']['Mean']['E_CO2e_IPPU_OperFor']['Ensemble Mean'][iT,0,0]+ \
                                                       mos['Delta'][sc]['ByStrata']['Mean']['E_CO2e_LULUCF_NEE']['Ensemble Mean'][iT,0,0]) )/np.sum(mos['Delta'][sc]['ByStrata']['Mean']['GJ PelletExport']['Ensemble Mean'][iT,0,0])
    d['EI Pellet Manufacture (tCO2e/ODT Pellets)']=np.sum( (mos['Delta'][sc]['ByStrata']['Mean']['E_CO2e_ESC_OperFor']['Ensemble Mean'][iT,0,0]+mos['Delta'][sc]['ByStrata']['Mean']['E_CO2e_IPPU_OperFor']['Ensemble Mean'][iT,0,0]) )/np.sum(mos['Delta'][sc]['ByStrata']['Mean']['ODT PelletExport']['Ensemble Mean'][iT,0,0])
    d['EI OperationForestry (tCO2e/ODT)']= \
            np.sum( (mos['Delta'][sc]['ByStrata']['Mean']['E_CO2e_ESC_OperFor']['Ensemble Mean'][iT,0,0]+mos['Delta'][sc]['ByStrata']['Mean']['E_CO2e_ET_OperFor']['Ensemble Mean'][iT,0,0]+mos['Delta'][sc]['ByStrata']['Mean']['E_CO2e_IPPU_OperFor']['Ensemble Mean'][iT,0,0]) ) / \
            np.sum( (mos['Delta'][sc]['ByStrata']['Mean']['ODT Lumber']['Ensemble Mean'][iT,0,0]+ \
            mos['Delta'][sc]['ByStrata']['Mean']['ODT LogExport']['Ensemble Mean'][iT,0,0]+ \
            mos['Delta'][sc]['ByStrata']['Mean']['ODT Plywood']['Ensemble Mean'][iT,0,0]+ \
            mos['Delta'][sc]['ByStrata']['Mean']['ODT OSB']['Ensemble Mean'][iT,0,0]+ \
            mos['Delta'][sc]['ByStrata']['Mean']['ODT MDF']['Ensemble Mean'][iT,0,0]+ \
            mos['Delta'][sc]['ByStrata']['Mean']['ODT Paper']['Ensemble Mean'][iT,0,0]+ \
            mos['Delta'][sc]['ByStrata']['Mean']['ODT PelletExport']['Ensemble Mean'][iT,0,0]+ \
            mos['Delta'][sc]['ByStrata']['Mean']['ODT PelletDomGrid']['Ensemble Mean'][iT,0,0]+ \
            mos['Delta'][sc]['ByStrata']['Mean']['ODT PelletDomRNG']['Ensemble Mean'][iT,0,0]+ \
            mos['Delta'][sc]['ByStrata']['Mean']['ODT PowerFacilityDom']['Ensemble Mean'][iT,0,0]+ \
            mos['Delta'][sc]['ByStrata']['Mean']['ODT PowerGrid']['Ensemble Mean'][iT,0,0]) )

    d['Energy efficiency (GJ/ODT)']=np.sum(mos['Delta'][sc]['ByStrata']['Mean']['GJ PelletExport']['Ensemble Mean'][iT,0,0])/np.sum(mos['Delta'][sc]['ByStrata']['Mean']['ODT PelletExport']['Ensemble Mean'][iT,0,0])

    # Displacement factor

    SubTot=-np.sum(mos['Delta'][sc]['ByStrata']['Mean']['E_CO2e_SUB_Tot']['Ensemble Mean'][iT,0,0])
    SubE=-np.sum(mos['Delta'][sc]['ByStrata']['Mean']['E_CO2e_SUB_E']['Ensemble Mean'][iT,0,0])
    SubM=-np.sum(mos['Delta'][sc]['ByStrata']['Mean']['E_CO2e_SUB_M']['Ensemble Mean'][iT,0,0])
    Wood=np.sum(mos['Delta'][sc]['ByStrata']['Mean']['C_ToLumber']['Ensemble Mean'][iT,0,0]+ \
        mos['Delta'][sc]['ByStrata']['Mean']['C_ToPlywood']['Ensemble Mean'][iT,0,0]+ \
        mos['Delta'][sc]['ByStrata']['Mean']['C_ToMDF']['Ensemble Mean'][iT,0,0]+ \
        mos['Delta'][sc]['ByStrata']['Mean']['C_ToOSB']['Ensemble Mean'][iT,0,0])
    Bioenergy=np.sum(mos['Delta'][sc]['ByStrata']['Mean']['E_CO2e_ESC_Bioenergy']['Ensemble Mean'][iT,0,0])
    BioenergyPlusOps=np.sum(mos['Delta'][sc]['ByStrata']['Mean']['E_CO2e_ESC_Bioenergy']['Ensemble Mean'][iT,0,0]+ \
                            mos['Delta'][sc]['ByStrata']['Mean']['E_CO2e_ESC_OperFor']['Ensemble Mean'][iT,0,0]+ \
                            mos['Delta'][sc]['ByStrata']['Mean']['E_CO2e_ET_OperFor']['Ensemble Mean'][iT,0,0]+ \
                            mos['Delta'][sc]['ByStrata']['Mean']['E_CO2e_IPPU_OperFor']['Ensemble Mean'][iT,0,0])
    BioenergyOpsAndHWPDecay=np.sum(mos['Delta'][sc]['ByStrata']['Mean']['E_CO2e_ESC_Bioenergy']['Ensemble Mean'][iT,0,0]+ \
                            mos['Delta'][sc]['ByStrata']['Mean']['E_CO2e_ESC_OperFor']['Ensemble Mean'][iT,0,0]+ \
                            mos['Delta'][sc]['ByStrata']['Mean']['E_CO2e_ET_OperFor']['Ensemble Mean'][iT,0,0]+ \
                            mos['Delta'][sc]['ByStrata']['Mean']['E_CO2e_IPPU_OperFor']['Ensemble Mean'][iT,0,0]+ \
                            mos['Delta'][sc]['ByStrata']['Mean']['E_CO2e_ESC_Bioenergy']['Ensemble Mean'][iT,0,0])

    d['Displacement Factor Total (tC/tC)']=SubTot/BioenergyOpsAndHWPDecay
    d['Displacement Factor Energy (tC/tC)']=SubE/BioenergyPlusOps
    d['Displacement Factor Materials (tC/tC)']=(SubM/3.667)/Wood

    return d