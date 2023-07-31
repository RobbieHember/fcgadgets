
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import geopandas as gpd
import os
import glob
import openpyxl
import copy
import gc as garc
import time
from shapely.geometry import Point
from fcgadgets.macgyver import util_general as gu
from fcgadgets.macgyver import util_gis as gis
from fcgadgets.macgyver import util_inventory as invu
from fcgadgets.hardhat import economics as econo
from fcgadgets.taz import aspatial_stat_models as asm

#%% CONVERT LUT NUMBER TO STRING NAME

def lut_n2s(dc,numb):
    if numb!=9999:
        vals=np.fromiter(dc.values(),dtype=float)
        keys=np.fromiter(dc.keys(),dtype='<U70')
        ind=np.where(vals==numb)[0]
        s=keys[ind]
    else:
        s=np.array(['Unidentified'],ndmin=1)
    return s

#%% Index to batch

def IndexToBatch(m,iBat):
    iStart=m['Project']['Batch Interval']*iBat
    iStop=np.minimum(m['Project']['N Stand'],iStart+m['Project']['Batch Interval'])
    indBat=np.arange(iStart,iStop,1)
    return indBat

#%% QUERY RESULTS CODES FOR MANAGEMENT ACTIVITY TYPES

def QueryResultsActivity(d):

    # Convert to arrays with at least 1d
    for key in d:
        d[key]=np.array(d[key],ndmin=1)

    Name=[]
    for i in range(d['SILV_BASE_CODE'].size):

        if (d['SILV_BASE_CODE'][i]=='FE') & (d['SILV_TECHNIQUE_CODE'][i]=='CA'):
            Name.append('Fertilization Aerial')

        elif (d['SILV_BASE_CODE'][i]=='FE') & (d['SILV_TECHNIQUE_CODE'][i]=='CG') & (d['SILV_METHOD_CODE'][i]!='BAGS'):
            Name.append('Fertilization Hand')

        elif (d['SILV_BASE_CODE'][i]=='FE') & (d['SILV_TECHNIQUE_CODE'][i]=='CG') & (d['SILV_METHOD_CODE'][i]=='BAGS'):
            Name.append('Fertilization Teabag')

        elif (d['SILV_BASE_CODE'][i]=='FE') & (d['SILV_TECHNIQUE_CODE'][i]=='OG'):
            Name.append('Fertilization Organic')

        elif (d['SILV_BASE_CODE'][i]=='PL') & (d['SILV_METHOD_CODE'][i]!='LAYOT'):
            # The planting in road rehab projects falls into this milestone type
            Name.append('Planting')

        elif (d['SILV_BASE_CODE'][i]=='DS') & (d['SILV_TECHNIQUE_CODE'][i]!='GS'):
            # Everything except grass seeding
            Name.append('Direct Seeding')

        elif (d['SILV_BASE_CODE'][i]=='PC') & (d['SILV_TECHNIQUE_CODE'][i]=='MA') & (d['SILV_OBJECTIVE_CODE_1'][i]=='DM'):
            # This will exclude SILV_TECHNIQUE_CODE=BI. Virtually all of it is mechanical.
            Name.append('Dwarf Mistletoe Control')

        elif (d['SILV_BASE_CODE'][i]=='PC') & (d['SILV_TECHNIQUE_CODE'][i]=='CA') & (d['SILV_OBJECTIVE_CODE_1'][i]=='ID'):
            Name.append('IDW Control')

        elif (d['SILV_BASE_CODE'][i]=='RD') & (d['SILV_BASE_CODE'][i]=='UP'):
            Name.append('Road Rehab')

        elif (d['SILV_BASE_CODE'][i]=='SP') & (d['SILV_TECHNIQUE_CODE'][i]=='BU') | \
            (d['SILV_BASE_CODE'][i]=='SP') & (d['SILV_TECHNIQUE_CODE'][i]=='ME') & (d['SILV_METHOD_CODE'][i]=='PBURN'):
            Name.append('Slashpile Burn')

        elif (d['SILV_BASE_CODE'][i]=='SP') & (d['SILV_TECHNIQUE_CODE'][i]=='MA') & (d['SILV_METHOD_CODE'][i]=='Unidentified') | \
            (d['SILV_BASE_CODE'][i]=='SP') & (d['SILV_TECHNIQUE_CODE'][i]=='MA') & (d['SILV_METHOD_CODE'][i]=='GUARD') | \
            (d['SILV_BASE_CODE'][i]=='SP') & (d['SILV_TECHNIQUE_CODE'][i]=='MA') & (d['SILV_METHOD_CODE'][i]=='HAND') | \
            (d['SILV_BASE_CODE'][i]=='SP') & (d['SILV_TECHNIQUE_CODE'][i]=='MA') & (d['SILV_METHOD_CODE'][i]=='KNOCK') | \
            (d['SILV_BASE_CODE'][i]=='SP') & (d['SILV_TECHNIQUE_CODE'][i]=='MA') & (d['SILV_METHOD_CODE'][i]=='POWER') | \
            (d['SILV_BASE_CODE'][i]=='SP') & (d['SILV_TECHNIQUE_CODE'][i]=='MA') & (d['SILV_METHOD_CODE'][i]=='MANCT') | \
            (d['SILV_BASE_CODE'][i]=='SP') & (d['SILV_TECHNIQUE_CODE'][i]=='MA') & (d['SILV_METHOD_CODE'][i]=='MDOWN') | \
            (d['SILV_BASE_CODE'][i]=='SP') & (d['SILV_TECHNIQUE_CODE'][i]=='MA') & (d['SILV_METHOD_CODE'][i]=='PILE') | \
            (d['SILV_BASE_CODE'][i]=='SP') & (d['SILV_TECHNIQUE_CODE'][i]=='MA') & (d['SILV_METHOD_CODE'][i]=='SNAG') | \
            (d['SILV_BASE_CODE'][i]=='SP') & (d['SILV_TECHNIQUE_CODE'][i]=='ME') & (d['SILV_METHOD_CODE'][i]=='CABLE') | \
            (d['SILV_BASE_CODE'][i]=='SP') & (d['SILV_TECHNIQUE_CODE'][i]=='ME') & (d['SILV_METHOD_CODE'][i]=='GUARD') | \
            (d['SILV_BASE_CODE'][i]=='SP') & (d['SILV_TECHNIQUE_CODE'][i]=='ME') & (d['SILV_METHOD_CODE'][i]=='MDOWN') | \
            (d['SILV_BASE_CODE'][i]=='SP') & (d['SILV_TECHNIQUE_CODE'][i]=='ME') & (d['SILV_METHOD_CODE'][i]=='PILE') | \
            (d['SILV_BASE_CODE'][i]=='SP') & (d['SILV_TECHNIQUE_CODE'][i]=='ME') & (d['SILV_METHOD_CODE'][i]=='PUSH'):
            Name.append('Knockdown')

        elif (d['SILV_BASE_CODE'][i]=='SP') & (d['SILV_TECHNIQUE_CODE'][i]=='ME') & (d['SILV_METHOD_CODE'][i]=='BRIP'):
            Name.append('Ripping')

        elif (d['SILV_BASE_CODE'][i]=='SP') & (d['SILV_TECHNIQUE_CODE'][i]=='ME') & (d['SILV_METHOD_CODE'][i]=='DISC'):
            Name.append('Disc Trenching')

        elif (d['SILV_BASE_CODE'][i]=='SP') & (d['SILV_TECHNIQUE_CODE'][i]=='ME') & (d['SILV_METHOD_CODE'][i]=='MULCH'):
            Name.append('Mulching')

        elif (d['SILV_BASE_CODE'][i]=='SP') & (d['SILV_TECHNIQUE_CODE'][i]=='ME') & (d['SILV_METHOD_CODE'][i]=='HARV'):
            Name.append('Harvest Salvage')

        elif (d['SILV_BASE_CODE'][i]=='LB') & (d['SILV_TECHNIQUE_CODE'][i]=='GR'):
            Name.append('LB-GR')

        elif (d['SILV_BASE_CODE'][i]=='SU'):
            Name.append('Surveys')

        else:
            Name.append('Undefined')

    return Name

#%% CONVERT DICTIONARY TO DATA STRUCTURE CLASS

class BunchDictionary(dict):
    def __init__(self, *args, **kwds):
        super(BunchDictionary, self).__init__(*args, **kwds)
        self.__dict__ = self


#%% Build event chronology from spreadsheet

def BuildEventChronologyFromSpreadsheet(meta,pNam):

    # Import inventory to get BGC zone
    iScn=0
    iBat=0
    inv=gu.ipickle(meta['Paths'][pNam]['Input Scenario'][iScn] + '\\Inventory_Bat' + FixFileNum(iBat) + '.pkl')

    # Simulate wildfires
    if (meta[pNam]['Scenario'][iScn]['Wildfire Status Pre-modern']=='On') | (meta[pNam]['Scenario'][iScn]['Wildfire Status Modern']=='On') | (meta[pNam]['Scenario'][iScn]['Wildfire Status Future']=='On'):
        asm.SimulateWildfireFromAAO(meta,inv)

    # Simulate MPB
    if (meta[pNam]['Scenario'][iScn]['MPB Status Pre-modern']=='On') | (meta[pNam]['Scenario'][iScn]['MPB Status Modern']=='On') | (meta[pNam]['Scenario'][iScn]['MPB Status Future']=='On'):
        asm.SimulateIBMFromAAO(meta,inv)

    for iScn in range(meta[pNam]['Project']['N Scenario']):

        for iEns in range(meta[pNam]['Project']['N Ensemble']):

            # Import wildfire simulations from Taz
            if (meta[pNam]['Scenario'][iScn]['Wildfire Status Pre-modern']=='On') | (meta[pNam]['Scenario'][iScn]['Wildfire Status Modern']=='On') | (meta[pNam]['Scenario'][iScn]['Wildfire Status Future']=='On'):

                wf_sim=gu.ipickle(meta['Paths'][pNam]['Data'] + '\\Inputs\\Ensembles\\wf_sim_Scn' + FixFileNum(iScn) + '_Ens' + FixFileNum(iEns) + '.pkl')
                if 'idx' in wf_sim:
                    idx=wf_sim['idx']
                    tmp=wf_sim.copy()
                    for v in ['Occurrence','Mortality']:
                        wf_sim[v]=np.zeros((meta[pNam]['Project']['N Time'],meta[pNam]['Project']['N Stand']),dtype='int16')
                        wf_sim[v][idx[0],idx[1]]=tmp[v]
                    del tmp

            # Import IBM
            if (meta[pNam]['Scenario'][iScn]['MPB Status Pre-modern']=='On') | (meta[pNam]['Scenario'][iScn]['MPB Status Modern']=='On') | (meta[pNam]['Scenario'][iScn]['MPB Status Future']=='On'):

                ibm_sim=gu.ipickle(meta['Paths'][pNam]['Data'] + '\\Inputs\\Ensembles\\ibm_sim_Scn' + FixFileNum(iScn) + '_Ens' + FixFileNum(iEns) + '.pkl')
                if 'idx' in ibm_sim:
                    idx=ibm_sim['idx']
                    tmp=ibm_sim.copy()
                    for v in ['Occurrence','Mortality']:
                        ibm_sim[v]=np.zeros((meta[pNam]['Project']['N Time'],meta[pNam]['Project']['N Stand']),dtype='int16')
                        ibm_sim[v][idx[0],idx[1]]=tmp[v]
                    del tmp

            for iBat in range(meta[pNam]['Project']['N Batch']):

                # Index to batch
                indBat=IndexToBatch(meta[pNam],iBat)

                # Always just one stand
                iS=0

                tv=np.arange(meta[pNam]['Project']['Year Start'],meta[pNam]['Project']['Year End']+1,1)

                # Initialize dictionary
                ec={}
                ec['ID Event Type']=np.zeros((meta[pNam]['Year'].size,indBat.size,meta['Core']['Max Events Per Year']),dtype='int16')
                ec['Mortality Factor']=np.zeros((meta[pNam]['Year'].size,indBat.size,meta['Core']['Max Events Per Year']),dtype='int16')
                ec['Growth Factor']=9999*np.ones((meta[pNam]['Year'].size,indBat.size,meta['Core']['Max Events Per Year']),dtype='int16')
                ec['ID Growth Curve']=np.zeros((meta[pNam]['Year'].size,indBat.size,meta['Core']['Max Events Per Year']),dtype='int16')

                #----------------------------------------------------------
                # Add spinup events
                #----------------------------------------------------------

                # Spinup interval
                if meta[pNam]['Project']['Return Interval Source']=='Custom':
                    # From custom input
                    ivl_spin=meta[pNam]['Project']['Custom Return Interval']
                elif meta[pNam]['Project']['Return Interval Source']=='BGC Zone':
                    # BGC Zone values
                    cd=meta[pNam]['Scenario'][iScn]['BGC Zone Code']
                    ind=np.where(meta['Param']['BE']['BGC Zone Averages']['Name']==cd)[0]
                    ivl_spin=meta['Param']['BE']['BGC Zone Averages']['Disturbance Return Interval'][ind]
                else:
                    print('Spin-up return interval source incorrect.')

                YearRef=meta[pNam]['Scenario'][iScn]['Year1_DisFromInv']
                AgeRef=meta[pNam]['Scenario'][iScn]['Age1_DisFromInv']
                if AgeRef>=0:
                    Year=np.arange(YearRef-AgeRef-100*ivl_spin,YearRef-AgeRef+ivl_spin,ivl_spin)
                else:
                    Year1=meta[pNam]['Project']['Year Start']+ivl_spin
                    Year2=meta[pNam]['Project']['Spinup Year End']
                    Year=np.arange(Year1,Year2+1,ivl_spin)

                for iYr in range(Year.size):
                    iT=np.where(tv==Year[iYr])[0]
                    ec['ID Event Type'][iT,:,0]=meta['LUT']['Event'][meta[pNam]['Project']['Spinup Disturbance Type']]
                    ec['Mortality Factor'][iT,:,0]=100
                    ec['Growth Factor'][iT,:,0]=9999
                    ec['ID Growth Curve'][iT,:,0]=meta[pNam]['Project']['Spinup Growth Curve ID']

                #----------------------------------------------------------
                # Add events from inventory
                #----------------------------------------------------------

                for iYr in range(1,10):

                    if ('Year' + str(iYr) + '_DisFromInv') not in meta[pNam]['Scenario'][iScn]:
                        continue

                    if np.isnan(meta[pNam]['Scenario'][iScn]['Year' + str(iYr) + '_DisFromInv'])==True:
                        continue

                    # If IDW, convert IDW class to growth and mortality factor
                    sc=np.array(['IDW-T','IDW-L','IDW-M','IDW-S','IDW-V','IDW-MM','IDW-MS','IDW-MV','IDW-SS','IDW-SV','IDW-VV'])
                    flg_i=0
                    indSc=np.where(sc==meta[pNam]['Scenario'][iScn]['Type' + str(iYr) + '_DisFromInv'])[0]
                    if indSc.size!=0:
                        if flg_i==0:
                            dfParDistBySC=pd.read_excel(meta['Paths']['Model']['Code'] + '\\Parameters\\Parameters_DisturbanceBySeverityClass.xlsx')
                            flg_i=1
                        indPar=np.where( (dfParDistBySC['Name']=='IDW') & (dfParDistBySC['SeverityCD']==sc[indSc[0]][4:]) )[0]
                        ID_TypeN=meta['LUT']['Event']['IDW']
                        MF=dfParDistBySC.loc[indPar,'Mortality Factor']
                        GF=dfParDistBySC.loc[indPar,'Growth Factor']
                    else:
                        ID_TypeS=meta[pNam]['Scenario'][iScn]['Type' + str(iYr) + '_DisFromInv']
                        try:
                            ID_TypeN=meta['LUT']['Event'][ID_TypeS]
                        except:
                            print(iScn)
                            print(iYr)
                            print(ID_TypeS)
                        MF=meta[pNam]['Scenario'][iScn]['Severity' + str(iYr) + '_DisFromInv']
                        GF=0

                    Year=meta[pNam]['Scenario'][iScn]['Year' + str(iYr) + '_DisFromInv']
                    iT=np.where(tv==Year)[0]

                    if iT.size==0:
                        print('Warning: An event was scheduled outside the timeframe of the simulation.')

                    iE=np.where(ec['ID Event Type'][iT,:,:]==0)[1]

                    ec['ID Event Type'][iT,:,iE[0]]=ID_TypeN
                    ec['Mortality Factor'][iT,:,iE[0]]=MF
                    ec['Growth Factor'][iT,:,iE[0]]=GF
                    ec['ID Growth Curve'][iT,:,iE[0]]=meta[pNam]['Scenario'][iScn]['GrowthCurve' + str(iYr) + '_DisFromInv']

                #----------------------------------------------------------
                # Add simulated wildfire from Taz
                #----------------------------------------------------------

                ind=np.array([],dtype=int)
                if meta[pNam]['Scenario'][iScn]['Wildfire Status Pre-modern']=='On':
                    ind0=np.where( (wf_sim['Occurrence'][:,iS]==1) & (meta[pNam]['Year']<1920) )[0]
                    ind=np.append(ind,ind0)
                if meta[pNam]['Scenario'][iScn]['Wildfire Status Modern']=='On':
                    ind0=np.where( (wf_sim['Occurrence'][:,iS]==1) & (meta[pNam]['Year']>=1920) & (meta[pNam]['Year']<meta[pNam]['Project']['Year Project']) )[0]
                    ind=np.append(ind,ind0)
                if meta[pNam]['Scenario'][iScn]['Wildfire Status Future']=='On':
                    ind0=np.where( (wf_sim['Occurrence'][:,iS]==1) & (meta[pNam]['Year']>=meta[pNam]['Project']['Year Project']) )[0]
                    ind=np.append(ind,ind0)

                if ind.size>0:

                    ID_Type=meta['LUT']['Event']['Wildfire']*np.ones(ind.size)
                    Year=tv[ind]
                    MortF=wf_sim['Mortality'][ind,iS]
                    GrowthF=9999*np.ones(ind.size)
                    ID_GrowthCurve=1*np.ones(ind.size)

                    for iYr in range(Year.size):
                        iT=np.where(tv==Year[iYr])[0]
                        ec['ID Event Type'][iT,:,0]=ID_Type[iYr]
                        ec['Mortality Factor'][iT,:,0]=MortF[iYr]
                        ec['Growth Factor'][iT,:,0]=GrowthF[iYr]
                        ec['ID Growth Curve'][iT,:,0]=ID_GrowthCurve[iYr]

                    #----------------------------------------------------------
                    # Add simulated MPB from Taz
                    #----------------------------------------------------------

                    ind=np.array([],dtype=int)
                    if meta[pNam]['Scenario'][iScn]['MPB Status Pre-modern']=='On':
                        ind0=np.where( (ibm_sim['Occurrence'][:,iS]==1) & (meta[pNam]['Year']<1920) )[0]
                        ind=np.append(ind,ind0)
                    if meta[pNam]['Scenario'][iScn]['MPB Status Modern']=='On':
                        ind0=np.where( (ibm_sim['Occurrence'][:,iS]==1) & (meta[pNam]['Year']>=1920) & (meta[pNam]['Year']<meta[pNam]['Project']['Year Project']) )[0]
                        ind=np.append(ind,ind0)
                    if meta[pNam]['Scenario'][iScn]['MPB Status Future']=='On':
                        ind0=np.where( (ibm_sim['Occurrence'][:,iS]==1) & (meta[pNam]['Year']>=meta[pNam]['Project']['Year Project']) )[0]
                        ind=np.append(ind,ind0)

                    if ind.size>0:

                        ID_Type=meta['LUT']['Event']['IBM']*np.ones(ind.size)
                        Year=tv[ind]
                        MortF=ibm_sim['Mortality'][ind,iS]
                        GrowthF=9999*np.ones(ind.size)
                        ID_GrowthCurve=1*np.ones(ind.size)

                        for iYr in range(Year.size):
                            iT=np.where(tv==Year[iYr])[0]
                            ec['ID Event Type'][iT,:,0]=ID_Type[iYr]
                            ec['Mortality Factor'][iT,:,0]=MortF[iYr]
                            ec['Growth Factor'][iT,:,0]=GrowthF[iYr]
                            ec['ID Growth Curve'][iT,:,0]=ID_GrowthCurve[iYr]

                #--------------------------------------------------------------
                # Save to file
                #--------------------------------------------------------------

                gu.opickle(meta['Paths'][pNam]['Input Scenario'][iScn] + '\\Events_Ens' + FixFileNum(iEns) + '_Bat' + FixFileNum(iBat) + '.pkl',ec)

    return

#%% Decompress event chronology

def EventChronologyDecompress(meta,pNam,ec,iScn,iEns,iBat):

    # Uncompress event chronology if it has been compressed
    if 'idx' in ec:
        idx=ec['idx']
        tmp=ec.copy()
        for v in ['ID Event Type','Mortality Factor','Growth Factor','ID Growth Curve']:
            ec[v]=np.zeros((meta[pNam]['Project']['N Time'],meta[pNam]['Project']['Batch Size'][iBat],meta['Core']['Max Events Per Year']),dtype='int16')
            ec[v][idx[0],idx[1],idx[2]]=tmp[v]
        del tmp

    return ec

#%% Fix ensemble name and numbering

def FixFileNum(ind):
    indStrFixed=str(ind+1)
    if len(indStrFixed)==1:
        indStrFixed='000' + indStrFixed
    elif len(indStrFixed)==2:
        indStrFixed='00' + indStrFixed
    elif len(indStrFixed)==3:
        indStrFixed='0' + indStrFixed
    return indStrFixed

#%% Configure project

def ImportProjectConfig(meta,pNam,**kwargs):

    #--------------------------------------------------------------------------
    # Initialize nested dictionaries
    #--------------------------------------------------------------------------

    if 'Core' not in meta:
        meta['Core']={}
    if 'Param' not in meta:
        meta['Param']={}
    if 'Modules' not in meta:
        meta['Modules']={}

    if pNam not in meta:
        meta[pNam]={}

    if 'Project' not in meta[pNam]:
        meta[pNam]['Project']={}

    if 'Scenario' not in meta[pNam]:
        meta[pNam]['Scenario']={}

    #--------------------------------------------------------------------------
    # Import project parameters from spreadsheet
    #--------------------------------------------------------------------------

    df=pd.read_excel(meta['Paths'][pNam]['Data'] + '\\Inputs\\ProjectConfig.xlsx',sheet_name='Project')
    for i in range(len(df)):
        Name=df['Name'].iloc[i]
        Value=df['Value'].iloc[i]
        if Name[-1]==':':
            # Exclude headers
            continue
        meta[pNam]['Project'][Name]=Value

    #--------------------------------------------------------------------------
    # Import look-up tables
    #--------------------------------------------------------------------------

    meta=Load_LUTs_Modelling(meta)

    #--------------------------------------------------------------------------
    # Define pool names
    #--------------------------------------------------------------------------

    # Pool names (ecosystem)
    # *** If you change this, you need to change the same list in "Update Parameters" ***
    meta['Core']['Name Pools Eco']=['StemMerch','StemNonMerch','Foliage','Branch','Bark','RootCoarse','RootFine', \
        'PiledStemMerch','PiledStemNonMerch','PiledBranch','PiledBark','PiledSnagStem','PiledSnagBranch', \
        'LitterVF','LitterF','LitterM','LitterS','SnagStem','SnagBranch','SoilVF','SoilF','SoilS']

    # Used in fert work: ,'LitterDecomp'

    # Number of ecosystem pools
    meta['Core']['N Pools Eco']=len(meta['Core']['Name Pools Eco'])

    # Pool names (products)
    meta['Core']['Name Pools Pro']=['SFH','MFH','Comm','Furn','Ship','Repairs', \
        'Other','Paper','EffluentPulp','PowerFacilityDom','PowerFacilityExport','PowerGrid', \
        'PelletExport','PelletDomGrid','PelletDomRNG','LogExport','FirewoodDom','FirewoodExport','DumpWood', \
        'DumpPaper','LandfillWoodDegradable','LandfillWoodNonDegradable', \
        'LandfillPaperDegradable','LandfillPaperNonDegradable','E_CO2','E_CH4']

    # Number of product pools
    meta['Core']['N Pools Pro']=len(meta['Core']['Name Pools Pro'])

    # List of output variables
    meta['Core']['Output Variable List']=[
        'A',
        'V_MerchLive',
        'V_MerchDead',
        'V_MerchTotal',
        'V_ToMillMerchLive',
        'V_ToMillMerchDead',
        'V_ToMillMerchTotal',
        'C_Forest_Tot',
        'C_Biomass_Tot',
        'C_Stemwood_Tot',
        'C_Foliage_Tot',
        'C_Branch_Tot',
        'C_Bark_Tot',
        'C_Root_Tot',
        'C_DeadWood_Tot',
        'C_DumpLandfill_Tot',
        'C_Piled_Tot',
        'C_G_Gross_Tot',
        'C_G_Net_Tot',
        'C_HWP_Tot',
        'C_InUse_Tot',
        'C_Buildings_Tot',
        'C_NonBuildings_Tot',
        'C_LF_Tot',
        'C_Litter_Tot',
        'C_M_Reg_Tot',
        'C_NPP_Tot',
        'C_RH_Tot',
        'C_Soil_Tot',
        'C_Soil_OHorizon',
        'C_Coal',
        'C_Oil',
        'C_Gas',
        'C_Limestone',
        'C_ToFirewoodDom',
        'C_ToFirewoodExport',
        'C_ToLogExport',
        'C_ToLumber',
        'C_ToMDF',
        'C_ToMill',
        'C_ToMillMerch',
        'C_ToMillNonMerch',
        'C_ToMillSnagStem',
        'C_ToOSB',
        'C_ToPaper',
        'C_ToPelletExport',
        'C_ToPelletDomGrid',
        'C_ToPelletDomRNG',
        'C_ToPlywood',
        'C_ToPowerFacilityDom',
        'C_ToPowerFacilityExport',
        'C_ToPowerGrid',
        'C_ToSlashpileBurnTot',
        'C_ToSlashpileBurnNonMerch',
        'E_CO2e_LULUCF_NEE',
        'E_CO2e_LULUCF_Denit',
        'E_CO2e_LULUCF_Other',
        'E_CO2e_LULUCF_OpenBurning',
        'E_CO2e_LULUCF_Wildfire',
        'E_CO2e_LULUCF_Fire',
        'E_CO2e_LULUCF_HWP',
        'E_CO2e_ESC_Bioenergy',
        'E_CO2e_ESC_BioenergyPowerFacilityDom',
        'E_CO2e_ESC_BioenergyPowerFacilityExport',
        'E_CO2e_ESC_BioenergyPowerGrid',
        'E_CO2e_ESC_BioenergyPelletExport',
        'E_CO2e_ESC_BioenergyPelletDomGrid',
        'E_CO2e_ESC_BioenergyPelletDomRNG',
        'E_CO2e_ESC_BioenergyFirewoodDom',
        'E_CO2e_ESC_BioenergyFirewoodExport',
        'E_CO2e_ESC_OperFor',
        'E_CO2e_ET_OperFor',
        'E_CO2e_IPPU_OperFor',
        'E_CO2e_SUB_E',
        'E_CO2e_SUB_M',
        'E_CO2e_SUB_Tot',
        'E_CO2e_SUB_Coal',
        'E_CO2e_SUB_Oil',
        'E_CO2e_SUB_Gas',
        'E_CO2e_SUB_ESC',
        'E_CO2e_SUB_ET',
        'E_CO2e_SUB_IPPU',
        'E_CO2e_SUB_Calcination',
        'E_CO2e_SUB_Sawnwood',
        'E_CO2e_SUB_Panel',
        'E_CO2e_SUB_PowerFacilityDom',
        'E_CO2e_SUB_PowerFacilityExport',
        'E_CO2e_SUB_PowerGrid',
        'E_CO2e_SUB_PelletExport',
        'E_CO2e_SUB_PelletDomGrid',
        'E_CO2e_SUB_PelletDomRNG',
        'E_CO2e_SUB_FirewoodDom',
        'E_CO2e_SUB_FirewoodExport',
        'E_CO2e_AGHGB_WOSub',
        'E_CO2e_AGHGB_WOSub_cumu',
        'E_CO2e_AGHGB_WSub',
        'E_CO2e_AGHGB_WSub_cumu',
        'ODT Sawnwood',
        'ODT Panel',
        'ODT Lumber',
        'ODT LogExport',
        'ODT Plywood',
        'ODT OSB',
        'ODT MDF',
        'ODT Paper',
        'ODT PelletExport',
        'ODT PelletDomGrid',
        'ODT PelletDomRNG',
        'ODT PowerFacilityDom',
        'ODT PowerGrid',
        'ODT FirewoodTot',
        'ODT FirewoodDom',
        'ODT Concrete',
        'ODT Steel',
        'ODT Aluminum',
        'ODT Plastic',
        'ODT Textile',
        'ODT Coal',
        'ODT Oil',
        'ODT Gas',
        'GJ PowerFacilityDom',
        'GJ PowerGrid',
        'GJ PelletExport',
        'GJ PelletDomGrid',
        'GJ PelletDomRNG',
        'GJ FirewoodDom',
        'Atm_CO2_In',
        'Atm_CH4_In',
        'Atm_N2O_In',
        'Atm_CO2_Out',
        'Atm_CH4_Out',
        'Atm_N2O_Out',
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
        'Revenue PelletExport',
        'Revenue PelletDom',
        'Revenue Plywood',
        'Revenue PowerFacilityDom',
        'Revenue PowerGrid',
        'Revenue Gross',
        'Revenue Gross Disc',
        'Revenue Gross Disc_cumu',
        'Revenue Net',
        'Revenue Net Disc',
        'Revenue Net Disc_cumu',
        'LogSizeEnhancement']

    #--------------------------------------------------------------------------
    # Define indices to each pool
    #--------------------------------------------------------------------------

    # Indices to ecosystem pools pools
    meta['Core']['iEP']={}; cnt=0
    for nam in meta['Core']['Name Pools Eco']:
        meta['Core']['iEP'][nam]=cnt
        cnt=cnt+1
    iEP=meta['Core']['iEP']
    meta['Core']['iEP']['BiomassTotal']=np.array([iEP['StemMerch'],iEP['StemNonMerch'],iEP['Foliage'],iEP['Branch'],iEP['Bark'],iEP['RootCoarse'],iEP['RootFine']])
    meta['Core']['iEP']['BiomassAboveground']=np.array([iEP['StemMerch'],iEP['StemNonMerch'],iEP['Foliage'],iEP['Branch'],iEP['Bark']])
    meta['Core']['iEP']['BiomassBelowground']=np.array([iEP['RootCoarse'],iEP['RootFine']])
    meta['Core']['iEP']['DeadWood']=np.array([iEP['PiledStemMerch'],iEP['PiledStemNonMerch'],iEP['PiledBranch'],iEP['PiledBark'],iEP['SnagStem'],iEP['SnagBranch']])
    meta['Core']['iEP']['Litter']=np.array([iEP['LitterVF'],iEP['LitterF'],iEP['LitterM'],iEP['LitterS']])
    meta['Core']['iEP']['Piled']=np.array([iEP['PiledStemMerch'],iEP['PiledStemNonMerch'],iEP['PiledBranch'],iEP['PiledBark'],iEP['PiledSnagStem'],iEP['PiledSnagBranch']])
    meta['Core']['iEP']['Soil']=np.array([iEP['SoilVF'],iEP['SoilF'],iEP['SoilS']])

    # Indices to produce pools pools
    meta['Core']['iPP']={}; cnt=0
    for nam in meta['Core']['Name Pools Pro']:
        meta['Core']['iPP'][nam]=cnt
        cnt=cnt+1
    iPP=meta['Core']['iPP']
    meta['Core']['iPP']['InUse']=np.array([ iPP['SFH'],iPP['MFH'],iPP['Comm'],iPP['Furn'],iPP['Ship'],iPP['Repairs'],iPP['Other'],iPP['Paper'] ])
    meta['Core']['iPP']['Buildings']=np.array([ iPP['SFH'],iPP['MFH'],iPP['Comm'] ])
    meta['Core']['iPP']['DumpLandfill']=np.array([ iPP['DumpWood'],iPP['DumpPaper'],iPP['LandfillWoodDegradable'],iPP['LandfillWoodNonDegradable'],iPP['LandfillPaperDegradable'],iPP['LandfillPaperNonDegradable'] ])

    #--------------------------------------------------------------------------
    # Maximum number of events per year
    # 8 appears to be sufficient but this may need to be changed for some
    # special projects
    #--------------------------------------------------------------------------

    meta['Core']['Max Events Per Year']=8

    #--------------------------------------------------------------------------
    # If Early Record Recycling is on, include a buffer before t_StartSaving
    #--------------------------------------------------------------------------

    meta['Core']['Recycle Early Record Buffer']=200

    #--------------------------------------------------------------------------
    # Define time
    #--------------------------------------------------------------------------

    # Calendar year
    meta[pNam]['Year']=np.arange(meta[pNam]['Project']['Year Start'],meta[pNam]['Project']['Year End']+1,1)

    meta[pNam]['Project']['N Time']=meta[pNam]['Year'].size

    #--------------------------------------------------------------------------
    # Define spatial domain for scripted projects
    #--------------------------------------------------------------------------

    if meta[pNam]['Project']['Scenario Source']=='Script':
        # Import land cover class
        zLCC1=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCoverClass1_Current.tif')

        #gdf_tsa=gpd.read_file(meta['Paths']['GDB']['GDB'] + '\\LandUse\\tsa.geojson')
        gdf_bcb=gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='NRC_POLITICAL_BOUNDARIES_1M_SP')

        # Import mask (0=excluded, 1=included)
        if meta[pNam]['Project']['ROI Source']=='Province':
            # Land mask for BC
            zMask=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])

        elif meta[pNam]['Project']['ROI Source']=='Regional District':
            # Mask from regional district
            zMask=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\TA_REGIONAL_DISTRICTS_SVW\\REGIONAL_DISTRICT_NAME.tif')
            ind=np.where( (zMask['Data']!=meta['LUT']['TA_REGIONAL_DISTRICTS_SVW']['REGIONAL_DISTRICT_NAME'][meta[pNam]['Project']['ROI Elements']]) )
            zMask['Data'][ind]=0
            ind=np.where( (zMask['Data']==meta['LUT']['TA_REGIONAL_DISTRICTS_SVW']['REGIONAL_DISTRICT_NAME'][meta[pNam]['Project']['ROI Elements']]) )
            zMask['Data'][ind]=1

        elif meta[pNam]['Project']['ROI Source']=='BCFCS_NMC':
            # Nutrient management
            zMask=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\FECA_MaskAll.tif')

        elif meta[pNam]['Project']['ROI Source']=='BCFCS_SPLC':
            # Non-obligation stand establishment
            zMask=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Management\\PL_NonOb_MaskAll.tif')

        elif meta[pNam]['Project']['ROI Source']=='EvalAtPlots':
            gplts=gu.ipickle(r'C:\Users\rhember\Documents\Data\GroundPlots\PSP-NADB2\Processed\L2\L2_BC.pkl')
            soils=gu.ipickle(r'C:\Users\rhember\Documents\Data\Soils\Shaw et al 2018 Database\SITES.pkl')
            x=np.append(gplts['sobs']['X'],soils['x'])
            y=np.append(gplts['sobs']['Y'],soils['y'])
            xy=np.unique(np.column_stack((x,y)),axis=0)
            zMask=zLCC1.copy()
            zMask['Data']=np.zeros(zLCC1['Data'].shape,dtype='int16')
            ind=gis.GetGridIndexToPoints(zLCC1,xy[:,0],xy[:,1])
            zMask['Data'][ind]=1
        else:
            print('ROI Source not recognized.')
            pass

        # Initialize geosptial info
        if 'Geos' not in meta.keys():
            meta['Geos']={}

        meta['Geos']['RGSF']=meta[pNam]['Project']['Regular grid sampling frequency (ha)']

        # Extract subgrid
        meta['Geos']['Grid']=zMask.copy()
        meta['Geos']['Grid']['Data']=np.zeros((zMask['Data'].shape),dtype='int8')
        meta['Geos']['Grid']=gis.UpdateGridCellsize(meta['Geos']['Grid'],meta['Geos']['RGSF'])

        # Resample required grids
        zLCC1_r=gis.UpdateGridCellsize(zLCC1,meta['Geos']['RGSF'])
        zMask_r=gis.UpdateGridCellsize(zMask,meta['Geos']['RGSF'])

        # Define additional sampling criteria
        if meta[pNam]['Project']['Land Cover Scope']=='Forest':

            iMask_Full=np.where( (zMask['Data']==1) & (zLCC1['Data']==meta['LUT']['Derived']['lcc1']['Forest']) )
            meta['Geos']['iMask']=np.where( (zMask_r['Data']==1) & (zLCC1_r['Data']==meta['LUT']['Derived']['lcc1']['Forest']) )

        elif meta[pNam]['Project']['Land Cover Scope']=='All Land':

            iMask_Full=np.where( (zMask['Data']==1) & (zLCC1['Data']>0) & (zLCC1['Data']!=meta['LUT']['Derived']['lcc1']['Water']) )
            meta['Geos']['iMask']=np.where( (zMask_r['Data']==1) & (zLCC1_r['Data']>0) & (zLCC1_r['Data']!=meta['LUT']['Derived']['lcc1']['Water']) )

        else:
            print('Land Cover Scope not recognized.')

        # Area expansion factor
        meta['Geos']['AEF']=iMask_Full[0].size/meta['Geos']['iMask'][0].size
        #meta['Geos']['AEF']=bc_grid_size/(meta['Geos']['Grid']['m']*meta['Geos']['Grid']['n'])
        meta[pNam]['Project']['AEF']=meta['Geos']['AEF']

        # Revise mask
        meta['Geos']['Grid']['Data'][meta['Geos']['iMask']]=1

        # Generate sparse grid
        meta['Geos']['Sparse']={}
        meta['Geos']['Sparse']['X']=meta['Geos']['Grid']['X'][meta['Geos']['iMask']]
        meta['Geos']['Sparse']['Y']=meta['Geos']['Grid']['Y'][meta['Geos']['iMask']]
        meta['Geos']['Sparse']['ID_Admin']=zMask_r['Data'][meta['Geos']['iMask']]

        # Save to pickle file
        # Flatten coordinate matrices first to save space
        meta['Geos']['Grid']['X']=meta['Geos']['Grid']['X'][0,:]
        meta['Geos']['Grid']['Y']=meta['Geos']['Grid']['Y'][:,0]
        gu.opickle(meta['Paths'][pNam]['Data'] + '\\geos.pkl',meta['Geos'])

        print('Number of stands = ' + str(meta['Geos']['Sparse']['X'].size))

        # Save mask as geotiff
        #gis.SaveGeoTiff(z,meta['Paths'][pNam]['Data'] + '\\geos_grid.tiff')
        #plt.matshow(meta['Geos']['Mask'])

        # Save sparse points to geojson
        flg=1
        if flg==1:
            points=[]
            for k in range(meta['Geos']['Sparse']['X'].size):
                points.append(Point(meta['Geos']['Sparse']['X'][k],meta['Geos']['Sparse']['Y'][k]))
            gdf_sxy=gpd.GeoDataFrame({'geometry':points,'ID_Admin':meta['Geos']['Sparse']['ID_Admin']})
            gdf_sxy.crs=meta['Geos']['crs']
            gdf_sxy.to_file(meta['Paths'][pNam]['Data'] + '\\geos.geojson',driver='GeoJSON')

        # Plot map
        if meta['Geos']['Sparse']['X'].size<25000:
            plt.close('all')
            fig,ax=plt.subplots(figsize=gu.cm2inch(14,14)); ms=2
            #mngr=plt.get_current_fig_manager()
            #mngr.window.setGeometry(700,20,620,600)
            gdf_bcb.plot(ax=ax,facecolor=[0.8,0.8,0.8],edgecolor=[0,0,0],label='Political Boundary',linewidth=0.25,alpha=1)
            #meta['Geos']['Boundary'].plot(ax=ax,facecolor='none',edgecolor=[0,0,1],label='Political Boundary',linewidth=0.75,alpha=1)

            #tsa_boundaries.plot(ax=ax,facecolor='none',edgecolor=[0,0,0],linewidth=0.25)
            gdf_sxy.plot(ax=ax,markersize=ms,facecolor=[0.75,0,0],edgecolor=None,linewidth=0.75,alpha=1)
            ax.grid(color='k',linestyle='-',linewidth=0.25)
            ax.set(position=[0.01,0.01,0.98,0.98],xticks=[],yticks=[])
            #gu.PrintFig(meta['Paths']['Figures'] + '\\Map_GroundPlots','png',900)
            #plt.savefig(PathProject + '\\SparseGrid_Map.png',format='png',dpi=900)
            #plt.close('all')
            #garc.collect()

    #--------------------------------------------------------------------------
    # Dimensions of simulation
    #--------------------------------------------------------------------------

    # Number of stands
    if meta[pNam]['Project']['Scenario Source']=='Spreadsheet':

        if 'Number of stands' in meta[pNam]['Project']:
            # Override!
            meta[pNam]['Project']['N Stand']=meta[pNam]['Project']['Number of stands']
        else:
            meta[pNam]['Project']['N Stand']=1

    elif meta[pNam]['Project']['Scenario Source']=='Portfolio':

        meta[pNam]['Project']['N Stand']=meta[pNam]['Project']['N Stand per Activity Type']*meta[pNam]['Project']['AIL']['N AT']*meta[pNam]['Project']['AIL']['N Years']

    elif meta[pNam]['Project']['Scenario Source']=='Script':

        try:
            meta[pNam]['Project']['N Stand']=meta['Geos']['Sparse']['X'].size
        except:
            # Tiled project
            meta[pNam]['Project']['N Stand']=meta[pNam]['Project']['iKeep'].size
            #meta[pNam]['Project']['N Stand']=kwargs['geos']['m']*kwargs['geos']['n']

    # Number of batches
    meta[pNam]['Project']['N Batch']=np.ceil(meta[pNam]['Project']['N Stand']/meta[pNam]['Project']['Batch Interval']).astype(int)

    # Initialize list that can keep track of batch sizes
    meta[pNam]['Project']['Batch Size']=[None]*meta[pNam]['Project']['N Batch']
    for iBat in range(meta[pNam]['Project']['N Batch']):
        meta[pNam]['Project']['Batch Size'][iBat]=IndexToBatch(meta[pNam],iBat).size

    #--------------------------------------------------------------------------
    # Import model parameters
    #--------------------------------------------------------------------------

    meta=ImportParameters(meta)

    #--------------------------------------------------------------------------
    # Define scenario parameters
    #--------------------------------------------------------------------------

    if meta[pNam]['Project']['Scenario Source']!='Portfolio':

        df=pd.read_excel(meta['Paths'][pNam]['Data'] + '\\Inputs\\ProjectConfig.xlsx',sheet_name='Scenarios') # ,usecols='A:OM'

        df=df.iloc[:,df.iloc[0,:].isnull().values==False]

        meta[pNam]['Scenario']=list()
        for i in range(1,df.shape[1]):
            pScn0={}
            for j in range(df.shape[0]):
                if df.iloc[j,0][-1]==':':
                    # Exclude headers
                    continue
                pScn0.update({df.iloc[j,0]:df.iat[j,i]})
            meta[pNam]['Scenario'].append(pScn0)

        # Number of scenarios
        meta[pNam]['Project']['N Scenario']=np.sum([i['Scenario Status']=='On' for i in meta[pNam]['Scenario']])

    #--------------------------------------------------------------------------
    # Number of Land Surface Scenarios
    #--------------------------------------------------------------------------

    meta[pNam]['Project']['LSC']={}
    if meta[pNam]['Project']['Scenario Source']=='Script':
        meta[pNam]['Project']['LSC']['Scenario Names Unique']=np.unique(np.array([meta[pNam]['Scenario'][i]['Land Surface Scenario'] for i in range(len(meta[pNam]['Scenario']))],dtype=object))
        meta[pNam]['Project']['LSC']['N Scenario']=meta[pNam]['Project']['LSC']['Scenario Names Unique'].size
    else:
        meta[pNam]['Project']['LSC']['Scenario Names Unique']=0
        meta[pNam]['Project']['LSC']['N Scenario']=0

    #--------------------------------------------------------------------------
    # Define strata for analyzing results (optional)
    #--------------------------------------------------------------------------
    meta[pNam]['Project']['Strata']={}
    meta[pNam]['Project']['Strata']['Project Type']={}
    meta[pNam]['Project']['Strata']['Project Type']['Unique CD']=np.array(['All'],dtype=object)
    meta[pNam]['Project']['Strata']['Project Type']['Unique ID']=np.array([77777],dtype='int32')
    meta[pNam]['Project']['Strata']['Project Type']['ID']=np.ones(meta[pNam]['Project']['N Stand'],dtype='int32')
    meta[pNam]['Project']['Strata']['Spatial']={}
    meta[pNam]['Project']['Strata']['Spatial']['Unique CD']=np.array(['All'],dtype=object)
    meta[pNam]['Project']['Strata']['Spatial']['Unique ID']=np.array([77777],dtype='int32')
    meta[pNam]['Project']['Strata']['Spatial']['ID']=np.ones(meta[pNam]['Project']['N Stand'],dtype='int32')
    meta[pNam]['Project']['Strata']['Year']={}
    meta[pNam]['Project']['Strata']['Year']['Unique CD']=np.array(['All'],dtype=object)
    meta[pNam]['Project']['Strata']['Year']['Unique ID']=np.array([77777],dtype='int32')
    meta[pNam]['Project']['Strata']['Year']['ID']=np.ones(meta[pNam]['Project']['N Stand'],dtype='int32')

    #--------------------------------------------------------------------------
    # Initialize project folders if they do not exist
    #--------------------------------------------------------------------------

    meta['Paths'][pNam]['Input Scenario']=[]
    meta['Paths'][pNam]['Output Scenario']=[]
    for iScn in range(0,meta[pNam]['Project']['N Scenario']):
        meta['Paths'][pNam]['Input Scenario'].append(meta['Paths'][pNam]['Data'] + '\\Inputs\\Scenario' + FixFileNum(iScn))
        if os.path.exists(meta['Paths'][pNam]['Input Scenario'][iScn])==False:
            os.mkdir(meta['Paths'][pNam]['Input Scenario'][iScn])
        meta['Paths'][pNam]['Output Scenario'].append(meta['Paths'][pNam]['Data'] + '\\Outputs\\Scenario' + FixFileNum(iScn))
        if os.path.exists(meta['Paths'][pNam]['Output Scenario'][iScn])==False:
            os.mkdir(meta['Paths'][pNam]['Output Scenario'][iScn])

    #--------------------------------------------------------------------------
    # Scale factors
    #--------------------------------------------------------------------------

    # *** Scale factor for saving results (this needs to be 100, 10 does not
    # capture carbon fluxes and it will affect GHG benefit estimates) ***
    # One variable ('CO2e_E_Products') requires the big one
    meta['Core']['Scale Factor Export Small']=0.001
    meta['Core']['Scale Factor Export Big']=0.001
    meta['Core']['Scale Factor C_M_ByAgent']=0.1

    #--------------------------------------------------------------------------
    # Define strings that frequently need to be populated with zeros
    #--------------------------------------------------------------------------

    meta['Core']['StringsToFill']=['Month','Day','SILV_FUND_SOURCE_CODE','FIA_PROJECT_ID','OPENING_ID','ACTUAL_TREATMENT_AREA','ACTUAL_PLANTED_NUMBER', \
               'PL_SPECIES_CD1','PL_SPECIES_PCT1','PL_SPECIES_GW1','PL_SPECIES_CD2','PL_SPECIES_PCT2','PL_SPECIES_GW2', \
               'PL_SPECIES_CD3','PL_SPECIES_PCT3','PL_SPECIES_GW3','PL_SPECIES_CD4','PL_SPECIES_PCT4','PL_SPECIES_GW4', \
               'PL_SPECIES_CD5','PL_SPECIES_PCT5','PL_SPECIES_GW5']

    #--------------------------------------------------------------------------
    # Growth curve information
    #--------------------------------------------------------------------------

    meta['Modules']['GYM']={}
    meta['Modules']['GYM']['N Growth Curves']=5
    meta['Modules']['GYM']['ID GC Unique']=np.array([1,2,3,4,5])
    meta['Modules']['GYM']['BatchTIPSY Maximum Age']=200
    meta['Modules']['GYM']['BatchTIPSY Column Names']=['Age','VolTot0','VolMerch125',
        'VolMerch175','ODT_Bark','ODT_Branch','ODT_Foliage','ODT_Roots',
        'ODT_Stem','MortalityVolumeTotal']

    # GC Input file variables
    meta['Modules']['GYM']['GC Input Indices']={'StemMerch':0,'StemNonMerch':1,'Bark':2,'Branch':3,'Foliage':4,'StemMerchV':5}

    # Scale factor for growth curves
    # Note: Do not change this to 0.1 - aerial Aerial BTK Spray response will not work properly at 0.1
    meta['Modules']['GYM']['Scale Factor']=0.001

    meta['Modules']['GYM']['GC_Variable_List']=['ID_Stand','ID_Scn','ID_GC','regeneration_method','s1','p1','i1','s2', \
        'p2','s3','p3','s4','p4','s5','p5','gain1','selage1','gain2','selage2', \
        'gain3','selage3','gain4','selage4','gain5','selage5', \
        'init_density','regen_delay','oaf1','oaf2','bec_zone','FIZ', \
        'fert_age1','fert_age2','fert_age3','fert_age4','fert_age5']

    # Operational adjustment factors (see TIPSY documentation)
    meta['Modules']['GYM']['OAF1 Default']=0.85
    meta['Modules']['GYM']['OAF2 Default']=0.95

    #--------------------------------------------------------------------------
    # Growth factor information
    # *** Not currently used ***
    #--------------------------------------------------------------------------

    #    # Default status of growth factors
    #    meta['Scenario Switch']['Net Growth Factor Status']=[None]*meta[pNam]['Project']['N Scenario']
    #    meta['Scenario Switch']['Mortality Factor Status']=[None]*meta[pNam]['Project']['N Scenario']
    #    for iScn in range(0,meta[pNam]['Project']['N Scenario']):
    #        meta['Scenario Switch']['Net Growth Factor Status'][iScn]='Off'
    #        meta['Scenario Switch']['Mortality Factor Status'][iScn]='Off'
    #        #meta['Scenario Switch'][iScn]['Status Net Growth Factor']='Off'
    #        #meta[pNam]['Scenario'][iScn]['Status Mortality Factor']='Off'

    #--------------------------------------------------------------------------
    # Harvested wood product information
    # Year to start calling annual HWP methods - running it before 1800 is a
    # waste of time.
    #--------------------------------------------------------------------------

    meta['Core']['HWP Year Start']=1850

    #--------------------------------------------------------------------------
    # Nutrient management information (for compatibility with "silviculture" module)
    #--------------------------------------------------------------------------

    # Initialize dictionary
    meta['Modules']['Nutrient Management']={}

    # Initialize index to stands affected by nutrient application
    # This needs to be populated with an empty array for when Sawtooth is used.
    meta['Modules']['Nutrient Management']['iApplication']=np.array([])

    # BGC zone exclusions (for on-the-fly application scheduler)
    meta['Modules']['Nutrient Management']['BGC Zone Exclusion CD']=['PP','MH','BAFA','BG','CMA','IMA']
    meta['Modules']['Nutrient Management']['BGC Zone Exclusion ID']=np.zeros(len(meta['Modules']['Nutrient Management']['BGC Zone Exclusion CD']))
    for iZ in range(len(meta['Modules']['Nutrient Management']['BGC Zone Exclusion CD'])):
        meta['Modules']['Nutrient Management']['BGC Zone Exclusion ID'][iZ]=meta['LUT']['VEG_COMP_LYR_R1_POLY']['BEC_ZONE_CODE'][ meta['Modules']['Nutrient Management']['BGC Zone Exclusion CD'][iZ] ]

    # Coastal zones used to make prob of occurrence region-specific
    meta['Modules']['Nutrient Management']['Coastal Zones CD']=['CWH','CDF']
    meta['Modules']['Nutrient Management']['Coastal Zones ID']=np.zeros(len(meta['Modules']['Nutrient Management']['Coastal Zones CD']))
    for iZ in range(len(meta['Modules']['Nutrient Management']['Coastal Zones CD'])):
        meta['Modules']['Nutrient Management']['Coastal Zones ID'][iZ]=meta['LUT']['VEG_COMP_LYR_R1_POLY']['BEC_ZONE_CODE'][ meta['Modules']['Nutrient Management']['Coastal Zones CD'][iZ] ]

    #--------------------------------------------------------------------------
    # Simulate random numbers that can be used for simulating harvest on the fly
    # The annual numbers will be the same among scenarios, but vary by ensemble
    #--------------------------------------------------------------------------

    # *** If you assign completely random numbers, random variation will occur among
    # scenarios, which can add considerable noise and demands many ensembles.
    # Conversely if you assign these pre-set sequeences, the random component will
    # vary among ensembles, but not among scenarios.
    meta[pNam]['Project']['On the Fly']={}
    meta[pNam]['Project']['On the Fly']['Random Numbers']={}
    meta[pNam]['Project']['On the Fly']['Random Numbers']['Scale Factor']=0.0001

    # Only create these files if they will be used

    # Not needed for portfolio projects

    if meta[pNam]['Project']['Scenario Source']!='Portfolio':

        flg_h=0
        for iScn in range(meta[pNam]['Project']['N Scenario']):
            if (meta[pNam]['Scenario'][iScn]['Harvest Status Historical']=='On') | (meta[pNam]['Scenario'][iScn]['Harvest Status Future']=='On'):
                flg_h=1
                break

        flg_b=0
        for iScn in range(meta[pNam]['Project']['N Scenario']):
            if (meta[pNam]['Scenario'][iScn]['Breakup Status Historical']=='On') | (meta[pNam]['Scenario'][iScn]['Breakup Status Future']=='On'):
                flg_b=1
                break

        # Create random numbers and save them
        if meta[pNam]['Project']['Frozen Ensembles Status']=='Off':
            for iEns in range(meta[pNam]['Project']['N Ensemble']):
                for iBat in range(meta[pNam]['Project']['N Batch']):

                    if flg_h==1:
                        rn=np.random.random( (meta[pNam]['Project']['N Time'],meta[pNam]['Project']['Batch Size'][iBat]) )
                        rn=rn/meta[pNam]['Project']['On the Fly']['Random Numbers']['Scale Factor']
                        rn=rn.astype('int16')
                        gu.opickle(meta['Paths'][pNam]['Data'] + '\\Inputs\\Ensembles\\RandomNumbers_Harvest_Ens' + FixFileNum(iEns) + '_Bat' + FixFileNum(iBat) + '.pkl',rn)

                    if flg_b==1:
                        rn=np.random.random( (meta[pNam]['Project']['N Time'],meta[pNam]['Project']['Batch Size'][iBat]) )
                        rn=rn/meta[pNam]['Project']['On the Fly']['Random Numbers']['Scale Factor']
                        rn=rn.astype('int16')
                        gu.opickle(meta['Paths'][pNam]['Data'] + '\\Inputs\\Ensembles\\RandomNumbers_Breakup_Ens' + FixFileNum(iEns) + '_Bat' + FixFileNum(iBat) + '.pkl',rn)

    #--------------------------------------------------------------------------
    # Parameter uncertainty by ensemble
    #--------------------------------------------------------------------------

    # Initialize list
    meta['Param']['By Ensemble']=[None]*meta[pNam]['Project']['N Ensemble']

    for iEns in range(meta[pNam]['Project']['N Ensemble']):

        # Initialize dictionary
        meta['Param']['By Ensemble'][iEns]={}

        #----------------------------------------------------------------------
        # Biomass turnover
        #----------------------------------------------------------------------

        if meta[pNam]['Project']['Uncertainty Status Biomass Turnover']=='On':

            meta['Param']['By Ensemble'][iEns]['Biomass Turnover']={}

            for k in meta['Param']['BE']['Biomass Turnover'].keys():

                mu=meta['Param']['BE']['Biomass Turnover'][k]
                sig=meta['Param']['Sigma']['Biomass Turnover'][k]
                bl=meta['Param']['BL']['Biomass Turnover'][k]
                bu=meta['Param']['BU']['Biomass Turnover'][k]

                r=np.random.normal(loc=mu,scale=mu*sig)

                if (bl!=-9999) & (bu!=-9999):
                    r=gu.Clamp(r,bl,bu)
                elif (bl!=-9999) & (bu==-9999):
                    r=np.maximum(bl,r)
                elif (bl==-9999) & (bu!=-9999):
                    r=np.minimum(bu,r)
                else:
                    r=r

                meta['Param']['By Ensemble'][iEns]['Biomass Turnover'][k]=r

        #----------------------------------------------------------------------
        # Decomposition
        #----------------------------------------------------------------------

        if meta[pNam]['Project']['Uncertainty Status Decomposition']=='On':

            meta['Param']['By Ensemble'][iEns]['Decomp']={}

            for k in meta['Param']['BE']['Decomp'].keys():

                mu=meta['Param']['BE']['Decomp'][k]
                sig=meta['Param']['Sigma']['Decomp'][k]
                bl=meta['Param']['BL']['Decomp'][k]
                bu=meta['Param']['BU']['Decomp'][k]

                r=np.random.normal(loc=mu,scale=mu*sig)

                if (bl!=-9999) & (bu!=-9999):
                    r=gu.Clamp(r,bl,bu)
                elif (bl!=-9999) & (bu==-9999):
                    r=np.maximum(bl,r)
                elif (bl==-9999) & (bu!=-9999):
                    r=np.minimum(bu,r)
                else:
                    r=r

                meta['Param']['By Ensemble'][iEns]['Decomp'][k]=r

        #----------------------------------------------------------------------
        # Harvesting
        #----------------------------------------------------------------------

        if meta[pNam]['Project']['Uncertainty Status Harvest Utilization']=='On':

            meta['Param']['By Ensemble'][iEns]['Event']={}

            EventList=['Harvest','Harvest Salvage']

            for Event in EventList:

                ID_Type=meta['LUT']['Event'][Event]
                meta['Param']['By Ensemble'][iEns]['Event'][ID_Type]={}

                #--------------------------------------------------------------
                # Biomass merch
                #--------------------------------------------------------------

                # Removed fraction
                mu=meta['Param']['BE']['Event'][ID_Type]['BiomassMerch_Removed']
                sig=np.array([0.1])
                bl=np.array([0.0])
                bu=np.array([1.0])
                r_Removed=np.random.normal(loc=mu,scale=mu*sig)
                r_Removed=gu.Clamp(r_Removed,bl,bu)
                meta['Param']['By Ensemble'][iEns]['Event'][ID_Type]['BiomassMerch_Removed']=r_Removed

                # Total fraction that is piled and dispersed
                r_PiledAndDispersed=1.0-r_Removed

                # Specific piled fraction
                mu=np.array([0.60]) # Specific fraction that is piled
                sig=np.array([0.1])
                bl=np.array([0.0])
                bu=np.array([1.0])
                rSpecific_Piled=np.random.normal(loc=mu,scale=mu*sig)
                rSpecific_Piled=gu.Clamp(rSpecific_Piled,bl,bu)

                # Piled fraction
                meta['Param']['By Ensemble'][iEns]['Event'][ID_Type]['BiomassMerch_Piled']=r_PiledAndDispersed*rSpecific_Piled

                # Specific dispersed fraction
                rSpecific_Dispersed=1.0-rSpecific_Piled

                # Dispersed fraction
                meta['Param']['By Ensemble'][iEns]['Event'][ID_Type]['BiomassMerch_LeftOnSite']=r_PiledAndDispersed*rSpecific_Dispersed

                #print(meta['Param']['By Ensemble'][iEns]['Event'][ID_Type]['BiomassMerch_Removed']+meta['Param']['By Ensemble'][iEns]['Event'][ID_Type]['BiomassMerch_Piled']+meta['Param']['By Ensemble'][iEns]['Event'][ID_Type]['BiomassMerch_LeftOnSite'])

                #--------------------------------------------------------------
                # Biomass non-merch
                #--------------------------------------------------------------

                # Removed fraction
                mu=meta['Param']['BE']['Event'][ID_Type]['BiomassNonMerch_Removed']
                sig=np.array([0.1])
                bl=np.array([0.0])
                bu=np.array([1.0])
                r_Removed=np.random.normal(loc=mu,scale=mu*sig)
                r_Removed=gu.Clamp(r_Removed,bl,bu)
                meta['Param']['By Ensemble'][iEns]['Event'][ID_Type]['BiomassNonMerch_Removed']=r_Removed

                # Total fraction that is piled and dispersed
                r_PiledAndDispersed=1.0-r_Removed

                # Specific piled fraction
                mu=np.array([0.60]) # Specific fraction that is piled
                sig=np.array([0.1])
                bl=np.array([0.0])
                bu=np.array([1.0])
                rSpecific_Piled=np.random.normal(loc=mu,scale=mu*sig)
                rSpecific_Piled=gu.Clamp(rSpecific_Piled,bl,bu)

                # Piled fraction
                meta['Param']['By Ensemble'][iEns]['Event'][ID_Type]['BiomassNonMerch_Piled']=r_PiledAndDispersed*rSpecific_Piled

                # Specific dispersed fraction
                rSpecific_Dispersed=1.0-rSpecific_Piled

                # Dispersed fraction
                meta['Param']['By Ensemble'][iEns]['Event'][ID_Type]['BiomassNonMerch_LeftOnSite']=r_PiledAndDispersed*rSpecific_Dispersed

                #print(meta['Param']['By Ensemble'][iEns]['Event'][ID_Type]['BiomassNonMerch_Removed']+meta['Param']['By Ensemble'][iEns]['Event'][ID_Type]['BiomassNonMerch_Piled']+meta['Param']['By Ensemble'][iEns]['Event'][ID_Type]['BiomassNonMerch_LeftOnSite'])

                #--------------------------------------------------------------
                # Snags
                #--------------------------------------------------------------

                # Removed fraction
                mu=meta['Param']['BE']['Event'][ID_Type]['Snags_Removed']
                sig=np.array([0.1])
                bl=np.array([0.0])
                bu=np.array([1.0])
                r_Removed=np.random.normal(loc=mu,scale=mu*sig)
                r_Removed=gu.Clamp(r_Removed,bl,bu)
                meta['Param']['By Ensemble'][iEns]['Event'][ID_Type]['Snags_Removed']=r_Removed

                # Total fraction that is piled and dispersed
                r_PiledAndDispersed=1.0-r_Removed

                # Specific piled fraction
                mu=np.array([0.60]) # Specific fraction that is piled
                sig=np.array([0.1])
                bl=np.array([0.0])
                bu=np.array([1.0])
                rSpecific_Piled=np.random.normal(loc=mu,scale=mu*sig)
                rSpecific_Piled=gu.Clamp(rSpecific_Piled,bl,bu)

                # Piled fraction
                meta['Param']['By Ensemble'][iEns]['Event'][ID_Type]['Snags_Piled']=r_PiledAndDispersed*rSpecific_Piled

                # Specific dispersed fraction
                rSpecific_Dispersed=1.0-rSpecific_Piled

                # Dispersed fraction
                meta['Param']['By Ensemble'][iEns]['Event'][ID_Type]['Snags_LeftOnSite']=r_PiledAndDispersed*rSpecific_Dispersed

                #print(meta['Param']['By Ensemble'][iEns]['Event'][ID_Type]['Snags_Removed']+meta['Param']['By Ensemble'][iEns]['Event'][ID_Type]['Snags_Piled']+meta['Param']['By Ensemble'][iEns]['Event'][ID_Type]['Snags_LeftOnSite'])

        #----------------------------------------------------------------------
        # Substitution effects
        #----------------------------------------------------------------------

        if meta[pNam]['Project']['Uncertainty Status Substitution']=='On':

            meta['Param']['By Ensemble'][iEns]['Substitution']={}

            #vL=['LumberDisplacementFactor','PanelDisplacementFactor']
            vL=['SawnwoodFracDisplacingConcrete','SawnwoodFracDisplacingSteel','SawnwoodFracDisplacingAluminum',
                'SawnwoodFracDisplacingPlastic','SawnwoodFracDisplacingTextile','PanelFracDisplacingConcrete','PanelFracDisplacingSteel',
                'PanelFracDisplacingAluminum','PanelFracDisplacingPlastic','PanelFracDisplacingTextile','ResidualsFracDisplacingConcrete',
                'ResidualsFracDisplacingSteel','ResidualsFracDisplacingAluminum','ResidualsFracDisplacingPlastic','ResidualsFracDisplacingTextile',
                'DisplacementRatio_ConcreteForSawnwood','DisplacementRatio_ConcreteForPanel','DisplacementRatio_ConcreteForResiduals',
                'DisplacementRatio_SteelForSawnwood','DisplacementRatio_SteelForPanel','DisplacementRatio_SteelForResiduals',
                'DisplacementRatio_AluminumForSawnwood','DisplacementRatio_AluminumForPanel','DisplacementRatio_AluminumForResiduals',
                'DisplacementRatio_PlasticForSawnwood','DisplacementRatio_PlasticForPanel','DisplacementRatio_PlasticForResiduals',
                'DisplacementRatio_TextileForSawnwood','DisplacementRatio_TextileForPanel','DisplacementRatio_TextileForResiduals']
            for k in vL:
                mu=meta['Param']['BE']['Substitution'][k]
                sig=meta['Param']['Sigma']['Substitution'][k]
                r=np.array(np.random.normal(loc=mu,scale=mu*sig),dtype=float)
                r=np.maximum(0,r)
                meta['Param']['By Ensemble'][iEns]['Substitution'][k]=r

            lv0=np.array(['PowerFacilityDomFracDisplacingRenewables','PowerFacilityDomFracDisplacingCoal','PowerFacilityDomFracDisplacingDiesel','PowerFacilityDomFracDisplacingNaturalGas','PowerFacilityDomFracDisplacingOil'])
            lv=lv0[np.argsort(np.random.random(len(lv0)))]
            r_remain=1.0
            for k in lv:
                mu=meta['Param']['BE']['Substitution'][k]
                sig=meta['Param']['Sigma']['Substitution'][k]
                r=np.array(np.random.normal(loc=mu,scale=mu*sig),dtype=float)
                r=gu.Clamp(r,0.0,r_remain)
                meta['Param']['By Ensemble'][iEns]['Substitution'][k]=r
                r_remain=r_remain-r

            lv0=np.array(['PowerFacilityExportFracDisplacingRenewables','PowerFacilityExportFracDisplacingCoal','PowerFacilityExportFracDisplacingDiesel','PowerFacilityExportFracDisplacingNaturalGas','PowerFacilityExportFracDisplacingOil'])
            lv=lv0[np.argsort(np.random.random(len(lv0)))]
            r_remain=1.0
            for k in lv:
                mu=meta['Param']['BE']['Substitution'][k]
                sig=meta['Param']['Sigma']['Substitution'][k]
                r=np.array(np.random.normal(loc=mu,scale=mu*sig),dtype=float)
                r=gu.Clamp(r,0.0,r_remain)
                meta['Param']['By Ensemble'][iEns]['Substitution'][k]=r
                r_remain=r_remain-r

            lv0=np.array(['PelletExportFracDisplacingRenewables','PelletExportFracDisplacingCoal','PelletExportFracDisplacingDiesel','PelletExportFracDisplacingNaturalGas','PelletExportFracDisplacingOil'])
            lv=lv0[np.argsort(np.random.random(len(lv0)))]
            r_remain=1.0
            for k in lv:
                mu=meta['Param']['BE']['Substitution'][k]
                sig=meta['Param']['Sigma']['Substitution'][k]
                r=np.array(np.random.normal(loc=mu,scale=mu*sig),dtype=float)
                r=gu.Clamp(r,0.0,r_remain)
                meta['Param']['By Ensemble'][iEns]['Substitution'][k]=r
                r_remain=r_remain-r

            lv0=np.array(['PelletDomGridFracDisplacingRenewables','PelletDomGridFracDisplacingCoal','PelletDomGridFracDisplacingDiesel','PelletDomGridFracDisplacingNaturalGas','PelletDomGridFracDisplacingOil'])
            lv=lv0[np.argsort(np.random.random(len(lv0)))]
            r_remain=1.0
            for k in lv:
                mu=meta['Param']['BE']['Substitution'][k]
                sig=meta['Param']['Sigma']['Substitution'][k]
                r=np.array(np.random.normal(loc=mu,scale=mu*sig),dtype=float)
                r=gu.Clamp(r,0.0,r_remain)
                meta['Param']['By Ensemble'][iEns]['Substitution'][k]=r
                r_remain=r_remain-r

            lv0=np.array(['PelletDomRNGFracDisplacingRenewables','PelletDomRNGFracDisplacingCoal','PelletDomRNGFracDisplacingDiesel','PelletDomRNGFracDisplacingNaturalGas','PelletDomRNGFracDisplacingOil'])
            lv=lv0[np.argsort(np.random.random(len(lv0)))]
            r_remain=1.0
            for k in lv:
                mu=meta['Param']['BE']['Substitution'][k]
                sig=meta['Param']['Sigma']['Substitution'][k]
                r=np.array(np.random.normal(loc=mu,scale=mu*sig),dtype=float)
                r=gu.Clamp(r,0.0,r_remain)
                meta['Param']['By Ensemble'][iEns]['Substitution'][k]=r
                r_remain=r_remain-r

            lv0=np.array(['FirewoodDomFracDisplacingRenewables','FirewoodDomFracDisplacingCoal','FirewoodDomFracDisplacingDiesel','FirewoodDomFracDisplacingNaturalGas','FirewoodDomFracDisplacingOil'])
            lv=lv0[np.argsort(np.random.random(len(lv0)))]
            r_remain=1.0
            for k in lv:
                mu=meta['Param']['BE']['Substitution'][k]
                sig=meta['Param']['Sigma']['Substitution'][k]
                r=np.array(np.random.normal(loc=mu,scale=mu*sig),dtype=float)
                r=gu.Clamp(r,0.0,r_remain)
                meta['Param']['By Ensemble'][iEns]['Substitution'][k]=r
                r_remain=r_remain-r

            lv0=np.array(['FirewoodExportFracDisplacingRenewables','FirewoodExportFracDisplacingCoal','FirewoodExportFracDisplacingDiesel','FirewoodExportFracDisplacingNaturalGas','FirewoodExportFracDisplacingOil'])
            lv=lv0[np.argsort(np.random.random(len(lv0)))]
            r_remain=1.0
            for k in lv:
                mu=meta['Param']['BE']['Substitution'][k]
                sig=meta['Param']['Sigma']['Substitution'][k]
                r=np.array(np.random.normal(loc=mu,scale=mu*sig),dtype=float)
                r=gu.Clamp(r,0.0,r_remain)
                meta['Param']['By Ensemble'][iEns]['Substitution'][k]=r
                r_remain=r_remain-r

            lv0=np.array(['PowerGridFracDisplacingRenewables','PowerGridFracDisplacingCoal','PowerGridFracDisplacingDiesel','PowerGridFracDisplacingNaturalGas','PowerGridFracDisplacingOil'])
            lv=lv0[np.argsort(np.random.random(len(lv0)))]
            r_remain=1.0
            for k in lv:
                mu=meta['Param']['BE']['Substitution'][k]
                sig=meta['Param']['Sigma']['Substitution'][k]
                r=np.array(np.random.normal(loc=mu,scale=mu*sig),dtype=float)
                r=gu.Clamp(r,0.0,r_remain)
                meta['Param']['By Ensemble'][iEns]['Substitution'][k]=r
                r_remain=r_remain-r

        #----------------------------------------------------------------------
        # Nutrient mmanagement
        #----------------------------------------------------------------------

        if meta[pNam]['Project']['Uncertainty Status Nutrient Application']=='On':

            meta['Param']['By Ensemble'][iEns]['Nutrient Management']={}

            for k in meta['Param']['BE']['Nutrient Management'].keys():
                mu=meta['Param']['BE']['Nutrient Management'][k]
                sig=meta['Param']['Sigma']['Nutrient Management'][k]
                r=np.random.normal(loc=mu,scale=mu*sig)
                r=np.maximum(0,r)
                meta['Param']['By Ensemble'][iEns]['Nutrient Management'][k]=r

        #----------------------------------------------------------------------
        # Mortality distribution
        #----------------------------------------------------------------------

        #meta[pNam]['Project']['Mortality Distribution']={}
        #meta[pNam]['Project']['Mortality Distribution']['Frequency (%)']=np.arange(0,101,1)
        #meta[pNam]['Project']['Mortality Distribution']['Severity (%)']=np.arange(0,101,1)
        #meta[pNam]['Project']['Mortality Distribution']['Variable List']=['Regular','Harvest','Wildfire','Beetles','Mechanical']
        #meta[pNam]['Project']['Mortality Distribution']['Data']={}
        #for v in meta[pNam]['Project']['Mortality Distribution']['Variable List']:
        #    meta[pNam]['Project']['Mortality Distribution']['Data'][v]=np.zeros((meta[pNam]['Project']['Mortality Distribution']['Frequency (%)'].size,meta[pNam]['Project']['Mortality Distribution']['Severity (%)'].size),dtype='int8')

    #--------------------------------------------------------------------------
    # Scenario info for portfolio projects
    #--------------------------------------------------------------------------

    if meta[pNam]['Project']['Scenario Source']=='Portfolio':

        # Index to rows with implementation
        indAT=np.where(np.sum(meta[pNam]['Project']['AIL']['Area'],axis=0)>0)[0]

        meta[pNam]['Project']['Portfolio']['ID Portfolio']=np.zeros(meta[pNam]['Project']['N Stand'],dtype=int)
        meta[pNam]['Project']['Portfolio']['ID AT']=np.zeros(meta[pNam]['Project']['N Stand'],dtype=int)
        meta[pNam]['Project']['Portfolio']['ID AT Unique']=np.zeros(meta[pNam]['Project']['N Stand'],dtype=int)
        meta[pNam]['Project']['Portfolio']['Area']=np.zeros(meta[pNam]['Project']['N Stand'])
        meta[pNam]['Project']['Portfolio']['Year']=np.zeros(meta[pNam]['Project']['N Stand'])
        meta[pNam]['Project']['Portfolio']['Region Code']=np.array(['' for _ in range(meta[pNam]['Project']['N Stand'])],dtype=object)
        meta[pNam]['Project']['Portfolio']['Felled Fate Scenario']=np.array(['' for _ in range(meta[pNam]['Project']['N Stand'])],dtype=object)
        meta[pNam]['Project']['Portfolio']['Removed Fate Scenario']=np.array(['' for _ in range(meta[pNam]['Project']['N Stand'])],dtype=object)
        meta[pNam]['Project']['Portfolio']['HWP End Use Scenario']=np.array(['' for _ in range(meta[pNam]['Project']['N Stand'])],dtype=object)
        cnt=0
        for iA in range(meta[pNam]['Project']['AIL']['N AT']):
            for iY in range(meta[pNam]['Project']['AIL']['N Years']):
                for iS in range(meta[pNam]['Project']['N Stand per Activity Type']):

                    ID_Portfolio=meta[pNam]['Project']['AIL']['ID Portfolio'][indAT[iA]]

                    meta[pNam]['Project']['Portfolio']['ID Portfolio'][cnt]=ID_Portfolio
                    meta[pNam]['Project']['Portfolio']['ID AT'][cnt]=meta[pNam]['Project']['AIL']['ID AT'][indAT[iA]]
                    meta[pNam]['Project']['Portfolio']['ID AT Unique'][cnt]=meta[pNam]['Project']['AIL']['ID AT Unique'][indAT[iA]]
                    meta[pNam]['Project']['Portfolio']['Area'][cnt]=meta[pNam]['Project']['AIL']['Area'][iY,indAT[iA]]
                    meta[pNam]['Project']['Portfolio']['Year'][cnt]=meta[pNam]['Project']['AIL']['Year'][iY]

                    # Region from Activities table
                    ind=np.where(meta[pNam]['Project']['Activities']['Activity ID']==meta[pNam]['Project']['Portfolio']['ID AT'][cnt])
                    meta[pNam]['Project']['Portfolio']['Region Code'][cnt]=meta[pNam]['Project']['Activities']['Region Code'][ind][0]

                    # Scenarios from Portfolio table
                    iPortfolio=np.where(meta[pNam]['Project']['Portfolio']['Raw']['ID_Portfolio']==ID_Portfolio)[0]
                    meta[pNam]['Project']['Portfolio']['Felled Fate Scenario'][cnt]=meta[pNam]['Project']['Portfolio']['Raw']['Felled Fate Scenario'][iPortfolio][0]
                    meta[pNam]['Project']['Portfolio']['Removed Fate Scenario'][cnt]=meta[pNam]['Project']['Portfolio']['Raw']['Removed Fate Scenario'][iPortfolio][0]
                    meta[pNam]['Project']['Portfolio']['HWP End Use Scenario'][cnt]=meta[pNam]['Project']['Portfolio']['Raw']['HWP End Use Scenario'][iPortfolio][0]
                    cnt=cnt+1

        # Scenario information
        # Will this tool ever be used with on-the-fly disturbances? Current set
        # to "off".
        meta[pNam]['Scenario']=[None]*meta[pNam]['Project']['N Scenario']
        for iScn in range(meta[pNam]['Project']['N Scenario']):

            meta[pNam]['Scenario'][iScn]={}

            # *** This is super awkward - simulations can't change between activities or scenarios!!!
            # Do we need that type of functionality for the PT? ***
            meta[pNam]['Scenario'][iScn]['Wildfire Scenario ID']=meta[pNam]['Project']['Activities']['Wildfire Scenario ID'][0]
            meta[pNam]['Scenario'][iScn]['Wildfire Status Pre-modern']=meta[pNam]['Project']['Activities']['Wildfire Status Pre-modern'][0]
            meta[pNam]['Scenario'][iScn]['Wildfire Status Modern']=meta[pNam]['Project']['Activities']['Wildfire Status Modern'][0]
            meta[pNam]['Scenario'][iScn]['Wildfire Status Future']=meta[pNam]['Project']['Activities']['Wildfire Status Future'][0]

            meta[pNam]['Scenario'][iScn]['Harvest Status Historical']='Off'
            meta[pNam]['Scenario'][iScn]['Harvest Status Future']='Off'
            meta[pNam]['Scenario'][iScn]['Breakup Status Historical']='Off'
            meta[pNam]['Scenario'][iScn]['Breakup Status Future']='Off'
            meta[pNam]['Scenario'][iScn]['Nutrient Application Status']='Off'

    #--------------------------------------------------------------------------
    # Grassland module
    #--------------------------------------------------------------------------

    for iScn in range(meta[pNam]['Project']['N Scenario']):
        if 'Grassland Module Status' not in meta[pNam]['Scenario'][iScn]:
            # Default is off
            meta[pNam]['Scenario'][iScn]['Grass Module Status']='Off'
            meta[pNam]['Scenario'][iScn]['Grass Module Year Start']=meta[pNam]['Project']['Year Project']+1

    #--------------------------------------------------------------------------
    # Initialize run time tracking
    #--------------------------------------------------------------------------

    meta[pNam]['Project']['Run Time Summary']={}

    return meta

#%% Load look-up-tables

def Load_LUTs_Modelling(meta):

    # Initialize LUTs dictionary
    if 'LUT' not in meta:
        meta['LUT']={}

    # Import distubance type
    p=gu.ReadExcel(meta['Paths']['Model']['Code'] + '\\Parameters\\' + 'Parameters_Disturbances.xlsx')
    meta['LUT']['Event']={}
    for i in range(p['Name'].size):
        meta['LUT']['Event'][p['Name'][i]]=p['ID'][i]

    # Region
    meta['LUT']['Region']={'Coast':1,'Interior':2,'GFS22':3}

    # # Added this to accommodate jupyter notebook demos - will need updating periodically
    # if 'Results' not in meta['Paths']:
    #     meta['Paths']['Results']=r'C:\Users\rhember\Documents\Data\ForestInventory\Results\20220422'
    # if 'VRI' not in meta['Paths']:
    #     meta['Paths']['VRI']=r'C:\Users\rhember\Documents\Data\ForestInventory\VRI\20220404'
    # if 'Disturbances' not in meta['Paths']:
    #     meta['Paths']['Disturbances']=r'C:\Users\rhember\Documents\Data\ForestInventory\Disturbances\20220422'
    # if 'LandUse' not in meta['Paths']:
    #     meta['Paths']['LandUse']=r'C:\Users\rhember\Documents\Data\ForestInventory\LandUse\20220422'

    # meta['LUT']['ATU']=gu.ipickle(meta['Paths']['Results'] + '\\LUTs_RSLT_ACTIVITY_TREATMENT_SVW.pkl')
    # meta['LUT']['OP']=gu.ipickle(meta['Paths']['Results'] + '\\LUTs_RSLT_OPENING_SVW.pkl')
    # meta['LUT']['PL']=gu.ipickle(meta['Paths']['Results'] + '\\LUTs_RSLT_PLANTING_SVW.pkl')
    # meta['LUT']['FC_I']=gu.ipickle(meta['Paths']['Results'] + '\\LUTs_RSLT_FOREST_COVER_INV_SVW.pkl')
    # meta['LUT']['FC_S']=gu.ipickle(meta['Paths']['Results'] + '\\LUTs_RSLT_FOREST_COVER_SILV_SVW.pkl')
    # meta['LUT']['VEG_COMP_LYR_R1_POLY']=gu.ipickle(meta['Paths']['VRI'] + '\\LUTs_VEG_COMP_LYR_R1_POLY.pkl')
    # meta['LUT']['BS']=gu.ipickle(meta['Paths']['Disturbances'] + '\\LUTs_VEG_BURN_SEVERITY_SP.pkl')
    # meta['LUT']['Pest']=gu.ipickle(meta['Paths']['Disturbances'] + '\\LUTs_PEST_INFESTATION_POLY.pkl')
    # meta['LUT']['FC_R']=gu.ipickle(meta['Paths']['Results'] + '\\LUTs_RSLT_FOREST_COVER_RESERVE_SVW.pkl')
    # #meta['LUT']['LU NL']=gu.ipickle(meta['Paths']['LandUse'] + '\\LUTs_RMP_PLAN_NON_LEGAL_POLY_SVW.pkl')
    # meta['LUT']['LU L']=gu.ipickle(meta['Paths']['LandUse'] + '\\LUTs_RMP_PLAN_LEGAL_POLY_SVW.pkl')
    # meta['LUT']['PARK']=gu.ipickle(meta['Paths']['LandUse'] + '\\LUTs_TA_PARK_ECORES_PA_SVW.pkl')
    # meta['LUT']['OGMA']=gu.ipickle(meta['Paths']['LandUse'] + '\\LUTs_RMP_OGMA_LEGAL_ALL_SVW.pkl')
    # meta['LUT']['OGSR']=gu.ipickle(meta['Paths']['LandUse'] + '\\LUTs_OGSR_TAP_PRIORITY_DEF_AREA_SP.pkl')
    # meta['LUT']['UWR']=gu.ipickle(meta['Paths']['LandUse'] + '\\LUTs_WCP_UNGULATE_WINTER_RANGE_SP.pkl')

    meta['LUT']['TIPSY']={}
    meta['LUT']['TIPSY']['FIZ']={'C':np.array(1,dtype=int),'I':np.array(2,dtype=int)}
    meta['LUT']['TIPSY']['regeneration_method']={'C':np.array(1,dtype=int),'N':np.array(2,dtype=int),'P':np.array(3,dtype=int)}

    # Land surface classification
    meta['LUT']['LSC']={}
    meta['LUT']['LSC']['Cover']={}
    data=gu.ReadExcel(meta['Paths']['Model']['Code'] + '\\Parameters\\Parameters_LSC_Cover.xlsx')
    for i in range(data['ID'].size):
        meta['LUT']['LSC']['Cover'][data['Name'][i]]=data['ID'][i]
    meta['LUT']['LSC']['Use']={}
    data=gu.ReadExcel(meta['Paths']['Model']['Code'] + '\\Parameters\\Parameters_LSC_Use.xlsx')
    for i in range(data['ID'].size):
        meta['LUT']['LSC']['Use'][data['Name'][i]]=data['ID'][i]

    # Species (for Sawtooth)
    #meta['LUT']['SRS']={}
    #for i in range(len(par['SRS']['SRS_CD'])):
    #    meta['LUT']['SRS'][par['SRS']['SRS_CD'][i]]=par['SRS']['SRS_ID'][i]

    return meta

#%% Load scenario results
# Return a list of dictionaries for each scenario. If multiple ensemble were run,
# the function will retun the average.

def LoadSingleOutputFile(meta,pNam,iScn,iEns,iBat):

    # Extract indices
    iEP=meta['Core']['iEP']

    # Extract biophysical parameters
    bB=meta['Param']['BE']['Biophysical']

    # Open batch results
    pth=meta['Paths'][pNam]['Data'] + '\\Outputs\\Scenario' + FixFileNum(iScn) + '\\Data_Scn' + FixFileNum(iScn) + '_Ens' + FixFileNum(iEns) + '_Bat' + FixFileNum(iBat) + '.pkl'
    v0=gu.ipickle(pth)

    # Convert to float and apply scale factor
    for k in v0.keys():

        # Skip mortality summary by agent
        if (k=='C_M_ByAgent') | (k=='C_M_Pct_ByAgent'):
            continue

        v0[k]=v0[k].astype(float)

        if (k=='E_CO2e_LULUCF_HWP') | (k=='E_CO2e_ESC_Comb') | (k=='E_CO2e_ET_Comb') | (k=='E_CO2e_IPPU_Comb'):
            v0[k]=v0[k]*meta['Core']['Scale Factor Export Big']
        else:
            v0[k]=v0[k]*meta['Core']['Scale Factor Export Small']

    #--------------------------------------------------------------------------
    # Add derived variables
    #--------------------------------------------------------------------------

    # Add year
    it=np.where(meta[pNam]['Year']>=meta[pNam]['Project']['Year Start Saving'])[0]
    v0['Year']=meta[pNam]['Year'][it]

    if meta[pNam]['Project']['Save Biomass Pools']=='On':

        # Aggregate variables not yet generated

        # Aggregate pools
        v0['C_Biomass_Tot']=np.sum(v0['C_Eco_Pools'][:,:,iEP['BiomassTotal']],axis=2)
        v0['C_Piled_Tot']=np.sum(v0['C_Eco_Pools'][:,:,iEP['Piled']],axis=2)
        v0['C_Litter_Tot']=np.sum(v0['C_Eco_Pools'][:,:,iEP['Litter']],axis=2)
        v0['C_DeadWood_Tot']=np.sum(v0['C_Eco_Pools'][:,:,iEP['DeadWood']],axis=2)
        v0['C_Soil_Tot']=np.sum(v0['C_Eco_Pools'][:,:,iEP['Soil']],axis=2)
        v0['C_InUse_Tot']=np.sum(v0['C_Pro_Pools'][:,:,meta['Core']['iPP']['InUse']],axis=2)
        v0['C_DumpLandfill_Tot']=np.sum(v0['C_Pro_Pools'][:,:,meta['Core']['iPP']['DumpLandfill']],axis=2)

        # Aggregate fluxes
        v0['C_G_Gross_Tot']=np.sum(v0['C_G_Gross'],axis=2)
        v0['C_G_Net_Tot']=np.sum(v0['C_G_Net'],axis=2)
        v0['C_M_Reg_Tot']=np.sum(v0['C_M_Reg'],axis=2)
        v0['C_LF_Tot']=np.sum(v0['C_LF'],axis=2)
        v0['C_RH_Tot']=np.sum(v0['C_RH'],axis=2)

    v0['C_NPP_Tot']=v0['C_G_Net_Tot']+v0['C_M_Reg_Tot']+v0['C_LF_Tot']
    v0['C_ToMill']=v0['C_ToMillMerch']+v0['C_ToMillNonMerch']+v0['C_ToMillSnagStem']
    v0['C_Forest_Tot']=v0['C_Biomass_Tot']+v0['C_DeadWood_Tot']+v0['C_Litter_Tot']+v0['C_Soil_Tot']
    v0['C_NonBuildings_Tot']=v0['C_InUse_Tot']-v0['C_Buildings_Tot']
    v0['C_HWP_Tot']=v0['C_InUse_Tot']+v0['C_DumpLandfill_Tot']
    v0['E_CO2e_LULUCF_NEE']=-1*bB['Ratio_CO2_to_C']*(v0['C_NPP_Tot']-v0['C_RH_Tot'])
    v0['E_CO2e_LULUCF_Fire']=v0['E_CO2e_LULUCF_Wildfire']+v0['E_CO2e_LULUCF_OpenBurning']

    v0['ODT Sawnwood']=v0['C_ToLumber']/bB['Carbon Content Wood']
    v0['ODT Panel']=(v0['C_ToPlywood']+v0['C_ToOSB']+v0['C_ToMDF'])/bB['Carbon Content Wood']
    v0['ODT Lumber']=v0['C_ToLumber']/bB['Carbon Content Wood']
    v0['ODT LogExport']=v0['C_ToLogExport']/bB['Carbon Content Wood']
    v0['ODT Plywood']=v0['C_ToPlywood']/bB['Carbon Content Wood']
    v0['ODT OSB']=v0['C_ToOSB']/bB['Carbon Content Wood']
    v0['ODT MDF']=v0['C_ToMDF']/bB['Carbon Content Wood']
    v0['ODT Paper']=v0['C_ToPaper']/bB['Carbon Content Wood']
    v0['ODT PelletExport']=v0['C_ToPelletExport']/bB['Carbon Content Wood']
    v0['ODT PelletDomGrid']=v0['C_ToPelletDomGrid']/bB['Carbon Content Wood']
    v0['ODT PelletDomRNG']=v0['C_ToPelletDomRNG']/bB['Carbon Content Wood']
    v0['ODT PowerGrid']=v0['C_ToPowerGrid']/bB['Carbon Content Wood']
    v0['ODT PowerFacilityDom']=v0['C_ToPowerFacilityDom']/bB['Carbon Content Wood']
    v0['ODT FirewoodDom']=v0['C_ToFirewoodDom']/bB['Carbon Content Wood']
    v0['ODT FirewoodTot']=(v0['C_ToFirewoodDom']+v0['C_ToFirewoodExport'])/bB['Carbon Content Wood']

    v0['ODT Coal']=v0['E_CO2e_SUB_CoalForBioenergy']/(bB['Emission Intensity Coal']/1000)/bB['Energy Content Coal']
    v0['ODT Oil']=v0['E_CO2e_SUB_OilForBioenergy']/(bB['Emission Intensity Oil']/1000)/bB['Energy Content Oil']
    v0['ODT Gas']=v0['E_CO2e_SUB_GasForBioenergy']/(bB['Emission Intensity Natural Gas']/1000)/bB['Energy Content Natural Gas']

    # Convert yield of bioenergy feedstock (ODT/ha) to energy (GJ/ha)
    v0['GJ PelletExport']=v0['ODT PelletExport']*bB['Energy Content Wood (0% moisture)']*bB['Electrical Conversion Efficiency of Pellet Electricity Plant (>25MW)']
    v0['GJ PelletDomGrid']=v0['ODT PelletDomGrid']*bB['Energy Content Wood (0% moisture)']*bB['Electrical Conversion Efficiency of Pellet Electricity Plant (>25MW)']
    v0['GJ PelletDomRNG']=v0['ODT PelletDomRNG']*bB['Energy Content Wood (0% moisture)']*bB['Electrical Conversion Efficiency of Pellet Electricity Plant (>25MW)']
    v0['GJ PowerGrid']=v0['ODT PowerGrid']*bB['Energy Content Wood (0% moisture)']*bB['Electrical Conversion Efficiency of Pellet Electricity Plant (>25MW)']
    v0['GJ PowerFacilityDom']=v0['ODT PowerFacilityDom']*bB['Energy Content Wood (0% moisture)']*bB['Electrical Conversion Efficiency of Pellet Electricity Plant (>25MW)']
    v0['GJ FirewoodDom']=v0['ODT FirewoodDom']*bB['Energy Content Wood (0% moisture)']*bB['Electrical Conversion Efficiency of Pellet Electricity Plant (>25MW)']
    v0['GJ FirewoodTot']=v0['ODT FirewoodTot']*bB['Energy Content Wood (0% moisture)']*bB['Electrical Conversion Efficiency of Pellet Electricity Plant (>25MW)']

    # Aggregate operational emissions
    v0['E_CO2e_ESC_OperFor']=v0['E_CO2e_ESC_OperForBurnCoal']+v0['E_CO2e_ESC_OperForBurnOil']+v0['E_CO2e_ESC_OperForBurnGas']
    v0['E_CO2e_ET_OperFor']=v0['E_CO2e_ET_OperForBurnCoal']+v0['E_CO2e_ET_OperForBurnOil']+v0['E_CO2e_ET_OperForBurnGas']
    v0['E_CO2e_IPPU_OperFor']=v0['E_CO2e_IPPU_OperForBurningCoal']+v0['E_CO2e_IPPU_OperForBurningOil']+v0['E_CO2e_IPPU_OperForBurningGas']

    v0['E_CO2e_OperForCoal']=v0['E_CO2e_ESC_OperForBurnCoal']+v0['E_CO2e_ET_OperForBurnCoal']+v0['E_CO2e_IPPU_OperForBurningCoal']
    v0['E_CO2e_OperForOil']=v0['E_CO2e_ESC_OperForBurnOil']+v0['E_CO2e_ET_OperForBurnOil']+v0['E_CO2e_IPPU_OperForBurningOil']
    v0['E_CO2e_OperForGas']=v0['E_CO2e_ESC_OperForBurnGas']+v0['E_CO2e_ET_OperForBurnGas']+v0['E_CO2e_IPPU_OperForBurningGas']
    v0['E_CO2e_OperForTot']=v0['E_CO2e_OperForCoal']+v0['E_CO2e_OperForOil']+v0['E_CO2e_OperForGas']

    # Revise the sign of substitution effects
    vL=['E_CO2e_SUB_CoalForBioenergy','E_CO2e_SUB_OilForBioenergy','E_CO2e_SUB_GasForBioenergy',
        'E_CO2e_SUB_PowerFacilityDom','E_CO2e_SUB_PowerFacilityExport','E_CO2e_SUB_PowerGrid',
        'E_CO2e_SUB_PelletExport','E_CO2e_SUB_PelletDomGrid','E_CO2e_SUB_PelletDomRNG',
        'E_CO2e_SUB_FirewoodDom','E_CO2e_SUB_FirewoodExport',
        'E_CO2e_SUB_CoalForWood','E_CO2e_SUB_OilForWood','E_CO2e_SUB_GasForWood',
        'E_CO2e_SUB_Sawnwood','E_CO2e_SUB_Panel',
        'E_CO2e_SUB_Concrete','E_CO2e_SUB_Steel','E_CO2e_SUB_Aluminum','E_CO2e_SUB_Plastic','E_CO2e_SUB_Textile',
        'E_CO2e_SUB_Calcination']
    for v in vL:
        v0[v]=-1*v0[v]

    # Reverse sign of building material production (saved as positive)
    vL=['ODT Concrete','ODT Steel','ODT Aluminum','ODT Plastic','ODT Textile','ODT Coal','ODT Oil','ODT Gas']
    for v in vL:
        v0[v]=-1*v0[v]

    v0['E_CO2e_SUB_Coal']=v0['E_CO2e_SUB_CoalForBioenergy']+v0['E_CO2e_SUB_CoalForWood']
    v0['E_CO2e_SUB_Oil']=v0['E_CO2e_SUB_OilForBioenergy']+v0['E_CO2e_SUB_OilForWood']
    v0['E_CO2e_SUB_Gas']=v0['E_CO2e_SUB_GasForBioenergy']+v0['E_CO2e_SUB_GasForWood']

    v0['E_CO2e_SUB_E']=v0['E_CO2e_SUB_CoalForBioenergy']+v0['E_CO2e_SUB_OilForBioenergy']+v0['E_CO2e_SUB_GasForBioenergy']
    v0['E_CO2e_SUB_M']=v0['E_CO2e_SUB_CoalForWood']+v0['E_CO2e_SUB_OilForWood']+v0['E_CO2e_SUB_GasForWood']+v0['E_CO2e_SUB_Calcination']
    v0['E_CO2e_SUB_Tot']=v0['E_CO2e_SUB_M']+v0['E_CO2e_SUB_E']

    v0['E_CO2e_SUB_ESC']=0.9*(v0['E_CO2e_SUB_Tot']-v0['E_CO2e_SUB_Calcination'])
    v0['E_CO2e_SUB_ET']=0.1*(v0['E_CO2e_SUB_Tot']-v0['E_CO2e_SUB_Calcination'])
    v0['E_CO2e_SUB_IPPU']=v0['E_CO2e_SUB_Calcination']

    v0['E_CO2e_Coal']=v0['E_CO2e_OperForCoal']+v0['E_CO2e_SUB_Coal']
    v0['E_CO2e_Oil']=v0['E_CO2e_OperForOil']+v0['E_CO2e_SUB_Oil']
    v0['E_CO2e_Gas']=v0['E_CO2e_OperForGas']+v0['E_CO2e_SUB_Gas']

    # Atmospheric GHG balance (tCO2e/ha/yr)
    v0['E_CO2e_AGHGB_WSub']=v0['E_CO2e_LULUCF_NEE']+v0['E_CO2e_LULUCF_Wildfire']+v0['E_CO2e_LULUCF_OpenBurning']+ \
        v0['E_CO2e_LULUCF_Denit']+v0['E_CO2e_LULUCF_Other']+v0['E_CO2e_LULUCF_HWP']+v0['E_CO2e_ESC_Bioenergy']+v0['E_CO2e_SUB_E']+v0['E_CO2e_SUB_M']+ \
        v0['E_CO2e_ESC_OperFor']+v0['E_CO2e_ET_OperFor']+v0['E_CO2e_IPPU_OperFor']

    v0['E_CO2e_AGHGB_WOSub']=v0['E_CO2e_LULUCF_NEE']+v0['E_CO2e_LULUCF_Wildfire']+v0['E_CO2e_LULUCF_OpenBurning']+ \
        v0['E_CO2e_LULUCF_Denit']+v0['E_CO2e_LULUCF_Other']+v0['E_CO2e_LULUCF_HWP']+v0['E_CO2e_ESC_Bioenergy']+ \
        v0['E_CO2e_ESC_OperFor']+v0['E_CO2e_ET_OperFor']+v0['E_CO2e_IPPU_OperFor']

    # Add cumulative
    v0['E_CO2e_AGHGB_WSub_cumu']=np.cumsum(v0['E_CO2e_AGHGB_WSub'],axis=0)
    v0['E_CO2e_AGHGB_WOSub_cumu']=np.cumsum(v0['E_CO2e_AGHGB_WOSub'],axis=0)

    # Add cumulative (starting from a specified start year)
    iT=np.where(v0['Year']>=meta[pNam]['Project']['Year Start Cumulative'])[0]
    v0['E_CO2e_AGHGB_WSub_cumu_from_tref']=np.zeros(v0['A'].shape)
    v0['E_CO2e_AGHGB_WSub_cumu_from_tref'][iT,:]=np.cumsum(v0['E_CO2e_AGHGB_WSub'][iT,:],axis=0)
    v0['E_CO2e_AGHGB_WOSub_cumu_from_tref']=np.zeros(v0['A'].shape)
    v0['E_CO2e_AGHGB_WOSub_cumu_from_tref'][iT,:]=np.cumsum(v0['E_CO2e_AGHGB_WOSub'][iT,:],axis=0)

    #--------------------------------------------------------------------------
    # Back-calculate production of fossil fuel consumption from operational use
    # and substitution effects (tonnesC)
    #--------------------------------------------------------------------------

    E_Op=v0['E_CO2e_ESC_OperForBurnCoal']+v0['E_CO2e_ET_OperForBurnCoal']+v0['E_CO2e_IPPU_OperForBurningCoal']
    E_Substitution=v0['E_CO2e_SUB_CoalForBioenergy']+v0['E_CO2e_SUB_CoalForWood']
    v0['C_Coal']=-1*(E_Op+E_Substitution)/(bB['Emission Intensity Coal']/1000)/bB['Energy Content Coal']*bB['Carbon Content Coal']

    E_Op=v0['E_CO2e_ESC_OperForBurnOil']+v0['E_CO2e_ET_OperForBurnOil']+v0['E_CO2e_IPPU_OperForBurningOil']
    E_Substitution=v0['E_CO2e_SUB_OilForBioenergy']+v0['E_CO2e_SUB_OilForWood']
    v0['C_Oil']=-1*(E_Op+E_Substitution)/(bB['Emission Intensity Oil']/1000)/bB['Energy Content Oil']*bB['Carbon Content Oil']

    E_Op=v0['E_CO2e_ESC_OperForBurnGas']+v0['E_CO2e_ET_OperForBurnGas']+v0['E_CO2e_IPPU_OperForBurningGas']
    E_Substitution=v0['E_CO2e_SUB_GasForBioenergy']+v0['E_CO2e_SUB_GasForWood']
    v0['C_Gas']=-1*(E_Op+E_Substitution)/(bB['Emission Intensity Natural Gas']/1000)/bB['Energy Content Natural Gas']*bB['Carbon Content Natural Gas']

    #--------------------------------------------------------------------------
    # Add and remove atmospheric gases
    # *** Fluxes from combustion and HWP decay already added ***
    #--------------------------------------------------------------------------

    # Add heterotrophic respiration
    v0['Atm_CO2_In']=v0['Atm_CO2_In']+bB['Ratio_CO2_to_C']*v0['C_RH_Tot']

    # Remove net primary productivity
    v0['Atm_CO2_Out']=v0['Atm_CO2_Out']+bB['Ratio_CO2_to_C']*v0['C_NPP_Tot']

    #v0['Atm_CH4_In']=v0['Atm_CH4_In'][iT,:]+E_CH4

    # Forestry operations
    v0['Atm_CO2_In']=v0['Atm_CO2_In']+v0['E_CO2e_OperForCoal']*bB['Emission Fraction Coal As CO2']
    v0['Atm_CH4_In']=v0['Atm_CH4_In']+v0['E_CO2e_OperForCoal']*bB['Emission Fraction Coal As CH4']
    v0['Atm_N2O_In']=v0['Atm_N2O_In']+v0['E_CO2e_OperForCoal']*bB['Emission Fraction Coal As N2O']

    v0['Atm_CO2_In']=v0['Atm_CO2_In']+v0['E_CO2e_OperForOil']*bB['Emission Fraction Oil As CO2']
    v0['Atm_CH4_In']=v0['Atm_CH4_In']+v0['E_CO2e_OperForOil']*bB['Emission Fraction Oil As CH4']
    v0['Atm_N2O_In']=v0['Atm_N2O_In']+v0['E_CO2e_OperForOil']*bB['Emission Fraction Oil As N2O']

    v0['Atm_CO2_In']=v0['Atm_CO2_In']+v0['E_CO2e_OperForGas']*bB['Emission Fraction Natural Gas As CO2']
    v0['Atm_CH4_In']=v0['Atm_CH4_In']+v0['E_CO2e_OperForGas']*bB['Emission Fraction Natural Gas As CH4']
    v0['Atm_N2O_In']=v0['Atm_N2O_In']+v0['E_CO2e_OperForGas']*bB['Emission Fraction Natural Gas As N2O']

    # Substitutions
    v0['Atm_CO2_Out']=v0['Atm_CO2_Out']-v0['E_CO2e_SUB_Coal']*bB['Emission Fraction Coal As CO2']
    v0['Atm_CH4_Out']=v0['Atm_CH4_Out']-v0['E_CO2e_SUB_Coal']*bB['Emission Fraction Coal As CH4']
    v0['Atm_N2O_Out']=v0['Atm_N2O_Out']-v0['E_CO2e_SUB_Coal']*bB['Emission Fraction Coal As N2O']

    v0['Atm_CO2_Out']=v0['Atm_CO2_Out']-v0['E_CO2e_SUB_Oil']*bB['Emission Fraction Oil As CO2']
    v0['Atm_CH4_Out']=v0['Atm_CH4_Out']-v0['E_CO2e_SUB_Oil']*bB['Emission Fraction Oil As CH4']
    v0['Atm_N2O_Out']=v0['Atm_N2O_Out']-v0['E_CO2e_SUB_Oil']*bB['Emission Fraction Oil As N2O']

    v0['Atm_CO2_Out']=v0['Atm_CO2_Out']-v0['E_CO2e_SUB_Gas']*bB['Emission Fraction Natural Gas As CO2']
    v0['Atm_CH4_Out']=v0['Atm_CH4_Out']-v0['E_CO2e_SUB_Gas']*bB['Emission Fraction Natural Gas As CH4']
    v0['Atm_N2O_Out']=v0['Atm_N2O_Out']-v0['E_CO2e_SUB_Gas']*bB['Emission Fraction Natural Gas As N2O']

    #--------------------------------------------------------------------------
    # Sawtooth variable adjustments
    #--------------------------------------------------------------------------

    if meta[pNam]['Project']['Biomass Module']=='Sawtooth':
        v0['N']=np.maximum(0,v0['N'])
        v0['N_R']=np.maximum(0,v0['N_R'])
        v0['N_M_Tot']=np.maximum(0,v0['N_M_Tot'])
        v0['TreeMean_D']=np.maximum(0,v0['TreeMean_D'])
        v0['TreeMean_Csw']=np.maximum(0,v0['TreeMean_Csw'])
        v0['TreeMean_Csw_G']=np.maximum(0,v0['TreeMean_Csw_G'])

    return v0

#%% LOAD SCENARIO RUSULTS
# Return a list of dictionaries for each scenario. If multiple ensemble were run,
# the function will retun the average.

def LoadScenarioResults(meta):

    # Initialize list that will contain scenarios
    v1=[]
    for iScn in range(meta[pNam]['Project']['N Scenario']):

        for iEns in range(meta[pNam]['Project']['N Ensemble']):

            for iBat in range(meta[pNam]['Project']['N Batch']):

                #--------------------------------------------------------------
                # Open batch results
                #--------------------------------------------------------------

                data_batch=LoadSingleOutputFile(meta,pNam,iScn,iEns,iBat)

                # Import event chronology
                if (meta[pNam]['Scenario'][iScn]['Harvest Status Future']=='On') | (meta[pNam]['Scenario'][iScn]['Breakup Status Historical']=='On') | (meta[pNam]['Scenario'][iScn]['Breakup Status Future']=='On'):
                    ec=gu.ipickle(meta['Paths'][pNam]['Input Scenario'][iScn] + '\\Modified_Events_Ens' + FixFileNum(iEns) + '_Bat' + FixFileNum(iBat) + '.pkl')
                else:
                    ec=gu.ipickle(meta['Paths'][pNam]['Input Scenario'][iScn] + '\\Events_Ens' + FixFileNum(iEns) + '_Bat' + FixFileNum(iBat) + '.pkl')

                # Uncompress event chronology if it has been compressed
                ec=EventChronologyDecompress(meta,ec,iScn,iEns,iBat)

                # Inventory
                inv=gu.ipickle(meta['Paths'][pNam]['Input Scenario'][iScn] + '\\Inventory_Bat' + FixFileNum(iBat) + '.pkl')

                # Cashflow
                econ=econo.CalculateNetRevenue(meta,iScn,iEns,iBat,inv,ec,data_batch)

                data_batch.update(econ)

                #--------------------------------------------------------------
                # Accumulate data in each batch
                #--------------------------------------------------------------

                if iBat==0:

                    data_all=data_batch

                else:

                    for key1 in data_batch.keys():

                        if key1=='Year':
                            # Only needed once
                            continue

                        elif (key1=='C_M_ByAgent') | (key1=='C_M_Pct_ByAgent'):
                            # Nested dictionary
                            for key2 in data_batch[key1].keys():
                                data_all[key1][key2]=np.append(data_all[key1][key2],data_batch[key1][key2],axis=1)

                        else:
                            # No nested dictionary
                            data_all[key1]=np.append(data_all[key1],data_batch[key1],axis=1)

            #------------------------------------------------------------------
            # Sum across ensembles
            #------------------------------------------------------------------

            if iEns==0:

                data_sum2ave=data_all

            else:

                for key1 in data_batch.keys():

                    if (key1=='C_M_ByAgent') | (key1=='C_M_Pct_ByAgent'):

                        # Nested dictionary
                        for key2 in data_batch[key1].keys():
                            data_sum2ave[key1][key2]=data_sum2ave[key1][key2]+data_all[key1][key2]

                    else:

                        # No nested dictionary
                        data_sum2ave[key1]=data_sum2ave[key1]+data_all[key1]

        #----------------------------------------------------------------------
        # If the simulation includes ensembles, calculate average
        #----------------------------------------------------------------------

        for key1 in data_batch.keys():

            # Skip mortality summary by agent
            if (key1=='C_M_ByAgent'):

                # Nested dictioanry
                for key2 in data_batch[key1].keys():
                    data_sum2ave[key1][key2]=data_sum2ave[key1][key2]/meta[pNam]['Project']['N Ensemble']

            else:

                # No nested dictionary
                data_sum2ave[key1]=data_sum2ave[key1]/meta[pNam]['Project']['N Ensemble']

        #----------------------------------------------------------------------
        # Add year
        #----------------------------------------------------------------------

        it=np.where(meta[pNam]['Year']>=meta[pNam]['Project']['Year Start Saving'])[0]
        data_sum2ave['Year']=meta[pNam]['Year'][it]

        #----------------------------------------------------------------------
        # Append to list
        #----------------------------------------------------------------------

        v1.append(data_sum2ave)

    return v1

#%% Calculate model output statistics for GHG variables (from points)

# Notes:

# You can't summarize the statistics in this script - the full set of ensemblers
# must be saved so that the statistics can be calculate for each scenario comparison.

# This also calculates age class distribution stats, but it does not include uncertainty
# (individual ensembles are not saved)

def Calc_MOS_FromPoints_GHG(meta,pNam,**kwargs):

    t0=time.time()
    print('Calculating model output statistics for GHG balance and age class distribution')

    # Key word arguments
    if 'StandsToInclude' in kwargs.keys():
        flag_stands_to_include=kwargs['StandsToInclude']
    else:
        flag_stands_to_include=[]
    if 'Scenarios' in kwargs.keys():
        scnL=kwargs['Scenarios']
    else:
        scnL=list(np.arange(0,meta[pNam]['Project']['N Scenario'],1).astype(int))

    # Time series of saved results
    tv_saving=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)

    # All non-economic values
    d1=LoadSingleOutputFile(meta,pNam,0,0,0)
    del d1['C_M_ByAgent']
    del d1['C_M_Pct_ByAgent']
    del d1['Year']

    if 'VariablesToKeep' not in kwargs.keys():
        v2include=list(d1.keys())
    else:
        v2include=kwargs['VariablesToKeep']

    # Add Sawtooth variables
    if meta[pNam]['Project']['Biomass Module']=='Sawtooth':
        v2include=v2include+['N','N_R','N_M_Tot','N_M_Reg','TreeMean_A','TreeMean_H','TreeMean_D','TreeMean_Csw','TreeMean_Csw_G']

    # Loop through scenarios
    for iScn0 in range(len(scnL)):

        iScn=scnL[iScn0]

        # Initialize data structures
        uPS_size=meta[pNam]['Project']['Strata']['Project Type']['Unique ID'].size
        uSS_size=meta[pNam]['Project']['Strata']['Spatial']['Unique ID'].size
        uYS_size=meta[pNam]['Project']['Strata']['Year']['Unique ID'].size

        # GHGs
        Data1={}
        for oper in ['Mean','Sum']:
            Data1[oper]={}
            for k in v2include:
                Data1[oper][k]=np.zeros( (tv_saving.size,meta[pNam]['Project']['N Ensemble'],uPS_size,uSS_size,uYS_size) ,dtype='float64')
                Data1[oper][k]=np.zeros( (tv_saving.size,meta[pNam]['Project']['N Ensemble'],uPS_size,uSS_size,uYS_size) ,dtype='float64')

        # Initialize age class distribution data
        acd={}
        acd['bwT']=10
        acd['binT']=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End'],acd['bwT'])
        acd['bwA']=5
        acd['binA']=np.arange(0,400,acd['bwA'])
        acd['Data']=np.zeros( (acd['binT'].size,acd['binA'].size,uPS_size,uSS_size,uYS_size) )

        # Loop through ensembles
        for iEns in range(meta[pNam]['Project']['N Ensemble']):

            # Initialize temporary data structure for full simulation
            Data0={}
            for k in v2include:
                Data0[k]=np.nan*np.empty( (tv_saving.size,meta[pNam]['Project']['N Stand']) ,dtype='float64')

            # Initialize age class distribution for this ensemble
            Age0=np.zeros( (acd['binT'].size,meta[pNam]['Project']['N Stand']) )

            for iBat in range(meta[pNam]['Project']['N Batch']):

                #print(str(iScn) + ' ' + str(iEns) + ' ' + str(iBat) )

                # Include specific subset of stands
                if len(flag_stands_to_include)==0:
                    iKeepStands=np.arange(0,meta[pNam]['Project']['Batch Size'][iBat],1,dtype=int)
                else:
                    iKeepStands=np.where(flag_stands_to_include[iScn][iEns][iBat]==1)[0]

                # Index to batch
                indBat=IndexToBatch(meta[pNam],iBat)

                d1=LoadSingleOutputFile(meta,pNam,iScn,iEns,iBat)

                for k in v2include:
                    if (d1[k].ndim>2):
                        # Skip C pools with more than 2 dims
                        continue
                    Data0[k][:,indBat[iKeepStands]]=d1[k][:,iKeepStands].copy()

                # Populate age class distribution
                for iT_acd in range(acd['binT'].size):
                    indT_acd=np.where(tv_saving==acd['binT'][iT_acd])[0]
                    Age0[iT_acd,indBat]=d1['A'][indT_acd,:]

            del d1
            garc.collect()

            # *** Fix bad RH - something wrong happens 1 in 500 simulations ***
            if 'C_RH_Tot' in Data0:
                iBad=np.where(Data0['C_RH_Tot']<0)
                if iBad[0].size>0:
                    Data0['C_RH_Tot'][iBad]=0

            # Summarize by project type, region, and time
            for iPS in range(uPS_size):
                for iSS in range(uSS_size):
                    for iYS in range(uYS_size):
                        if (meta[pNam]['Project']['Strata']['Project Type']['Unique CD'][iPS]=='All') & (meta[pNam]['Project']['Strata']['Spatial']['Unique CD'][iSS]=='All') & (meta[pNam]['Project']['Strata']['Year']['Unique CD'][iYS]=='All'):
                            # Index to all values for "All"
                            ind1=np.where(meta[pNam]['Project']['Strata']['Project Type']['ID']!=77777)[0]
                        elif (meta[pNam]['Project']['Strata']['Project Type']['Unique CD'][iPS]!='All') & (meta[pNam]['Project']['Strata']['Spatial']['Unique CD'][iSS]=='All') & (meta[pNam]['Project']['Strata']['Year']['Unique CD'][iYS]=='All'):
                            # Index to specific project type stratum
                            ind1=np.where( (meta[pNam]['Project']['Strata']['Project Type']['ID']==meta[pNam]['Project']['Strata']['Project Type']['Unique ID'][iPS]) )[0]
                        elif (meta[pNam]['Project']['Strata']['Project Type']['Unique CD'][iPS]=='All') & (meta[pNam]['Project']['Strata']['Spatial']['Unique CD'][iSS]!='All') & (meta[pNam]['Project']['Strata']['Year']['Unique CD'][iYS]=='All'):
                            # Index to specific spatial stratum
                            ind1=np.where( (meta[pNam]['Project']['Strata']['Spatial']['ID']==meta[pNam]['Project']['Strata']['Spatial']['Unique ID'][iSS]) )[0]
                        elif (meta[pNam]['Project']['Strata']['Project Type']['Unique CD'][iPS]=='All') & (meta[pNam]['Project']['Strata']['Spatial']['Unique CD'][iSS]=='All') & (meta[pNam]['Project']['Strata']['Year']['Unique CD'][iYS]!='All'):
                            # Index to specific year stratum
                            ind1=np.where( (meta[pNam]['Project']['Strata']['Year']['ID']==meta[pNam]['Project']['Strata']['Year']['Unique ID'][iSS]) )[0]
                        elif (meta[pNam]['Project']['Strata']['Project Type']['Unique CD'][iPS]!='All') & (meta[pNam]['Project']['Strata']['Year']['Unique CD'][iYS]=='All') | (meta[pNam]['Project']['Strata']['Spatial']['Unique CD'][iSS]!='All') & (meta[pNam]['Project']['Strata']['Year']['Unique CD'][iYS]=='All'):
                            # Index to specific project type and spatial stratum
                            ind1=np.where( (meta[pNam]['Project']['Strata']['Project Type']['ID']==meta[pNam]['Project']['Strata']['Project Type']['Unique ID'][iPS]) & \
                                      (meta[pNam]['Project']['Strata']['Spatial']['ID']==meta[pNam]['Project']['Strata']['Spatial']['Unique ID'][iSS]) )[0]
                        elif (meta[pNam]['Project']['Strata']['Project Type']['Unique CD'][iPS]!='All') & (meta[pNam]['Project']['Strata']['Year']['Unique CD'][iYS]!='All') & (meta[pNam]['Project']['Strata']['Spatial']['Unique CD'][iSS]=='All'):
                            # Index to specific project type and year stratum
                            ind1=np.where( (meta[pNam]['Project']['Strata']['Project Type']['ID']==meta[pNam]['Project']['Strata']['Project Type']['Unique ID'][iPS]) & \
                                      (meta[pNam]['Project']['Strata']['Year']['ID']==meta[pNam]['Project']['Strata']['Year']['Unique ID'][iYS]) )[0]
                        elif (meta[pNam]['Project']['Strata']['Project Type']['Unique CD'][iPS]=='All') & (meta[pNam]['Project']['Strata']['Year']['Unique CD'][iYS]!='All') & (meta[pNam]['Project']['Strata']['Spatial']['Unique CD'][iSS]!='All'):
                            # Index to specific spatial and year stratum
                            ind1=np.where( (meta[pNam]['Project']['Strata']['Spatial']['ID']==meta[pNam]['Project']['Strata']['Spatial']['Unique ID'][iSS]) & \
                                      (meta[pNam]['Project']['Strata']['Year']['ID']==meta[pNam]['Project']['Strata']['Year']['Unique ID'][iYS]) )[0]
                        else:
                            continue

                        for k in v2include:
                            Data1['Mean'][k][:,iEns,iPS,iSS,iYS]=np.mean(Data0[k][:,ind1],axis=1)
                            Data1['Sum'][k][:,iEns,iPS,iSS,iYS]=np.sum(Data0[k][:,ind1],axis=1)

                        # Add ensemble to age class distribution
                        for iT_acd in range(acd['binT'].size):
                            for iA_acd in range(acd['binA'].size):
                                ind_acd=np.where( np.abs(Age0[iT_acd,ind1]-acd['binA'][iA_acd])<=acd['bwA']/2 )[0]
                                acd['Data'][iT_acd,iA_acd,iPS,iSS]=acd['Data'][iT_acd,iA_acd,iPS,iSS,iYS]+ind_acd.size

        # Save GHGs
        gu.opickle(meta['Paths'][pNam]['Data'] + '\\Outputs\\MOS_ByStrata_GHGB_Scn' + str(iScn+1) + '.pkl',Data1)

        # Save age class distribution
        acd['Data']=acd['Data']/meta[pNam]['Project']['N Ensemble']
        gu.opickle(meta['Paths'][pNam]['Data'] + '\\Outputs\\MOS_ByStrata_AgeClassDist_Scn' + str(iScn+1) + '.pkl',acd)

    t1=time.time()
    print(str((t1-t0)/60) + ' min')
    return

#%% Calculate model output statistics for econ variables (from points)
# Notes: You can't summarize the statistics in this script - the full set of ensemblers
# must be saved so that the statistics can be calculate for each scenario comparison.

def Calc_MOS_FromPoints_Econ(meta,pNam,**kwargs):

    t0=time.time()
    print('Calculating model output statistics for economic variables')

    # Keyword argumenst
    if 'ScenariosToInclude' in kwargs.keys():
        ScenariosToInclude=kwargs['ScenariosToInclude']
    else:
        # Default is to not save individual ensembles
        ScenariosToInclude=np.arange(0,meta[pNam]['Project']['N Scenario'])
    if 'StandsToInclude' in kwargs.keys():
        flag_stands_to_include=kwargs['StandsToInclude']
    else:
        flag_stands_to_include=[]

    # Time series of saved results
    tv_saving=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)

    # All (excluding 'C_M_ByAgent' and 'Year','C_M_Dist')
    v2include=['Cost Roads','Cost Harvest Overhead','Cost Harvest Felling and Piling','Cost Harvest Hauling','Cost Harvest Residuals',
               'Cost Milling','Cost Nutrient Management','Cost Planting','Cost Survey','Cost Knockdown','Cost Ripping','Cost PAS Deactivation',
               'Cost Slashpile Burn','Cost Total','Cost Silviculture Total','Revenue Lumber','Revenue Plywood',
               'Revenue OSB','Revenue MDF','Revenue Paper','Revenue PowerFacilityDom','Revenue PowerGrid','Revenue PelletExport','Revenue PelletDom','Revenue FirewoodDom',
               'Revenue LogExport','Revenue Gross','Revenue Net','Revenue Net Disc','Revenue Gross Disc','Cost Total Disc','Cost Nutrient Management Disc',
               'Cost Silviculture Total Disc','Cost Total_cumu','Cost Silviculture Total_cumu','Cost Nutrient Management_cumu','Cost Total Disc_cumu',
               'Cost Silviculture Total Disc_cumu','Cost Nutrient Management Disc_cumu','Revenue Gross_cumu','Revenue Gross Disc_cumu','Revenue Net_cumu',
               'Revenue Net Disc_cumu']

    # Loop through scenarios
    for iS in range(ScenariosToInclude.size):

        iScn=ScenariosToInclude[iS]

        # Initialize data structures
        uPS_size=meta[pNam]['Project']['Strata']['Project Type']['Unique ID'].size
        uSS_size=meta[pNam]['Project']['Strata']['Spatial']['Unique ID'].size
        uYS_size=meta[pNam]['Project']['Strata']['Year']['Unique ID'].size

        Data1={}
        for oper in ['Mean','Sum']:
            Data1[oper]={}
            for k in v2include:
                Data1[oper][k]=np.zeros( (tv_saving.size,meta[pNam]['Project']['N Ensemble'],uPS_size,uSS_size,uYS_size) ,dtype='float64')
                Data1[oper][k]=np.zeros( (tv_saving.size,meta[pNam]['Project']['N Ensemble'],uPS_size,uSS_size,uYS_size) ,dtype='float64')

        # An option to skip economic calculations. The file will still be created, but all zeros
        if meta[pNam]['Project']['Skip Economics']=='Off':

            # Loop through ensembles
            for iEns in range(meta[pNam]['Project']['N Ensemble']):

                # Initialize temporary data structure for full simulation
                Data0={}
                for k in v2include:
                    Data0[k]=np.nan*np.empty( (tv_saving.size,meta[pNam]['Project']['N Stand']) ,dtype='float64')

                # Import batches
                for iBat in range(meta[pNam]['Project']['N Batch']):

                    # Include specific subset of stands
                    if len(flag_stands_to_include)==0:
                        iKeepStands=np.arange(0,meta[pNam]['Project']['Batch Size'][iBat],1,dtype=int)
                    else:
                        iKeepStands=np.where(flag_stands_to_include[iScn][iEns][iBat]==1)[0]

                    # Index to batch
                    indBat=IndexToBatch(meta[pNam],iBat)

                    d1=LoadSingleOutputFile(meta,pNam,iScn,iEns,iBat)

                    if (meta[pNam]['Scenario'][iScn]['Harvest Status Future']=='On') | (meta[pNam]['Scenario'][iScn]['Breakup Status Historical']=='On') | (meta[pNam]['Scenario'][iScn]['Breakup Status Future']=='On'):
                        ec=gu.ipickle(meta['Paths'][pNam]['Input Scenario'][iScn] + '\\Modified_Events_Ens' + FixFileNum(iEns) + '_Bat' + FixFileNum(iBat) + '.pkl')
                    else:
                        ec=gu.ipickle(meta['Paths'][pNam]['Input Scenario'][iScn] + '\\Events_Ens' + FixFileNum(iEns) + '_Bat' + FixFileNum(iBat) + '.pkl')

                    # Uncompress event chronology if it has been compressed
                    ec=EventChronologyDecompress(meta,pNam,ec,iScn,iEns,iBat)

                    # Inventory
                    inv=gu.ipickle(meta['Paths'][pNam]['Input Scenario'][iScn] + '\\Inventory_Bat' + FixFileNum(iBat) + '.pkl')

                    # Cashflow
                    econ1=econo.CalculateNetRevenue(meta,pNam,iScn,iEns,iBat,inv,ec,d1)

                    for k in v2include:
                        Data0[k][:,indBat[iKeepStands]]=econ1[k][:,iKeepStands].copy()

                del econ1,inv,ec
                #garc.collect()

                # Summarize by project type, region, and time
                for iPS in range(uPS_size):
                    for iSS in range(uSS_size):
                        for iYS in range(uYS_size):
                            if (meta[pNam]['Project']['Strata']['Project Type']['Unique CD'][iPS]=='All') & (meta[pNam]['Project']['Strata']['Spatial']['Unique CD'][iSS]=='All') & (meta[pNam]['Project']['Strata']['Year']['Unique CD'][iYS]=='All'):
                                # Index to all values for "All"
                                ind1=np.where(meta[pNam]['Project']['Strata']['Project Type']['ID']!=77777)[0]
                            elif (meta[pNam]['Project']['Strata']['Project Type']['Unique CD'][iPS]!='All') & (meta[pNam]['Project']['Strata']['Spatial']['Unique CD'][iSS]=='All') & (meta[pNam]['Project']['Strata']['Year']['Unique CD'][iYS]=='All'):
                                # Index to specific project type stratum
                                ind1=np.where( (meta[pNam]['Project']['Strata']['Project Type']['ID']==meta[pNam]['Project']['Strata']['Project Type']['Unique ID'][iPS]) )[0]
                            elif (meta[pNam]['Project']['Strata']['Project Type']['Unique CD'][iPS]=='All') & (meta[pNam]['Project']['Strata']['Spatial']['Unique CD'][iSS]!='All') & (meta[pNam]['Project']['Strata']['Year']['Unique CD'][iYS]=='All'):
                                # Index to specific spatial stratum
                                ind1=np.where( (meta[pNam]['Project']['Strata']['Spatial']['ID']==meta[pNam]['Project']['Strata']['Spatial']['Unique ID'][iSS]) )[0]
                            elif (meta[pNam]['Project']['Strata']['Project Type']['Unique CD'][iPS]=='All') & (meta[pNam]['Project']['Strata']['Spatial']['Unique CD'][iSS]=='All') & (meta[pNam]['Project']['Strata']['Year']['Unique CD'][iYS]!='All'):
                                # Index to specific year stratum
                                ind1=np.where( (meta[pNam]['Project']['Strata']['Year']['ID']==meta[pNam]['Project']['Strata']['Year']['Unique ID'][iSS]) )[0]
                            elif (meta[pNam]['Project']['Strata']['Project Type']['Unique CD'][iPS]!='All') & (meta[pNam]['Project']['Strata']['Year']['Unique CD'][iYS]=='All') | (meta[pNam]['Project']['Strata']['Spatial']['Unique CD'][iSS]!='All') & (meta[pNam]['Project']['Strata']['Year']['Unique CD'][iYS]=='All'):
                                # Index to specific project type and spatial stratum
                                ind1=np.where( (meta[pNam]['Project']['Strata']['Project Type']['ID']==meta[pNam]['Project']['Strata']['Project Type']['Unique ID'][iPS]) & \
                                          (meta[pNam]['Project']['Strata']['Spatial']['ID']==meta[pNam]['Project']['Strata']['Spatial']['Unique ID'][iSS]) )[0]
                            elif (meta[pNam]['Project']['Strata']['Project Type']['Unique CD'][iPS]!='All') & (meta[pNam]['Project']['Strata']['Year']['Unique CD'][iYS]!='All') & (meta[pNam]['Project']['Strata']['Spatial']['Unique CD'][iSS]=='All'):
                                # Index to specific project type and year stratum
                                ind1=np.where( (meta[pNam]['Project']['Strata']['Project Type']['ID']==meta[pNam]['Project']['Strata']['Project Type']['Unique ID'][iPS]) & \
                                          (meta[pNam]['Project']['Strata']['Year']['ID']==meta[pNam]['Project']['Strata']['Year']['Unique ID'][iYS]) )[0]
                            elif (meta[pNam]['Project']['Strata']['Project Type']['Unique CD'][iPS]=='All') & (meta[pNam]['Project']['Strata']['Year']['Unique CD'][iYS]!='All') & (meta[pNam]['Project']['Strata']['Spatial']['Unique CD'][iSS]!='All'):
                                # Index to specific spatial and year stratum
                                ind1=np.where( (meta[pNam]['Project']['Strata']['Spatial']['ID']==meta[pNam]['Project']['Strata']['Spatial']['Unique ID'][iSS]) & \
                                          (meta[pNam]['Project']['Strata']['Year']['ID']==meta[pNam]['Project']['Strata']['Year']['Unique ID'][iYS]) )[0]
                            else:
                                continue

                            for k in v2include:
                                Data1['Mean'][k][:,iEns,iPS,iSS,iYS]=np.mean(Data0[k][:,ind1],axis=1)
                                Data1['Sum'][k][:,iEns,iPS,iSS,iYS]=np.sum(Data0[k][:,ind1],axis=1)

        # Save
        gu.opickle(meta['Paths'][pNam]['Data'] + '\\Outputs\\MOS_ByStrata_Econ_Scn' + str(iScn+1) + '.pkl',Data1)

    t1=time.time()
    print(str((t1-t0)/60) + ' min')
    return

#%% Calculate model output statistics for area from points
# Notes: Unlike the GHG and Econ scripts, the area can calculate the statistics
# in this scripts (scenario comparisons don't need to include area)

def Calc_MOS_FromPoints_Area(meta,pNam,**kwargs):

    print('Calculating model output statistics for event areas')
    t0=time.time()

    # Key word arguments
    if 'StandsToInclude' in kwargs.keys():
        flag_stands_to_include=kwargs['StandsToInclude']
    else:
        flag_stands_to_include=[]

    tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)

    # Time series of saved results
    tv_saving=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)

    # Full time series (for event chronology)
    tv_full=np.arange(meta[pNam]['Project']['Year Start'],meta[pNam]['Project']['Year End']+1,1)

    # Loop through scenarios
    for iScn in range(meta[pNam]['Project']['N Scenario']):
        # Initialize data structures
        uPS_size=meta[pNam]['Project']['Strata']['Project Type']['Unique ID'].size
        uSS_size=meta[pNam]['Project']['Strata']['Spatial']['Unique ID'].size
        uYS_size=meta[pNam]['Project']['Strata']['Year']['Unique ID'].size

        Area1={}
        for oper in ['Mean','Sum']:
            Area1[oper]={}
            for k in meta['LUT']['Event'].keys():
                Area1[oper][k]=np.zeros( (tv_saving.size,meta[pNam]['Project']['N Ensemble'],uPS_size,uSS_size,uYS_size) ,dtype='float64')
                Area1[oper][k]=np.zeros( (tv_saving.size,meta[pNam]['Project']['N Ensemble'],uPS_size,uSS_size,uYS_size) ,dtype='float64')

        # An option to skip economic calculations. The file will still be created, but all zeros
        #if meta[pNam]['Project']['Skip Economics']=='Off':

        # Loop through ensembles
        for iEns in range(meta[pNam]['Project']['N Ensemble']):

            # Initialize temporary data structure for full simulation
            Area0={}
            for k in meta['LUT']['Event'].keys():
                Area0[k]=np.nan*np.empty( (tv_saving.size,meta[pNam]['Project']['N Stand']) ,dtype='float64')

            Area0_Full={}
            for k in meta['LUT']['Event'].keys():
                Area0_Full[k]=np.zeros( (tv_saving.size,meta[pNam]['Project']['N Stand']) ,dtype='float64')

            # Import batches
            for iBat in range(meta[pNam]['Project']['N Batch']):

                # Include specific subset of stands
                if len(flag_stands_to_include)==0:
                    iKeepStands=np.arange(0,meta[pNam]['Project']['Batch Size'][iBat],1,dtype=int)
                else:
                    iKeepStands=np.where(flag_stands_to_include[iScn][iEns][iBat]==1)[0]

                # Index to batches
                indBat=IndexToBatch(meta[pNam],iBat)

                # Import event chronology
                if (meta[pNam]['Scenario'][iScn]['Harvest Status Future']=='On') | (meta[pNam]['Scenario'][iScn]['Breakup Status Historical']=='On') | (meta[pNam]['Scenario'][iScn]['Breakup Status Future']=='On'):
                    ec=gu.ipickle(meta['Paths'][pNam]['Input Scenario'][iScn] + '\\Modified_Events_Ens' + FixFileNum(iEns) + '_Bat' + FixFileNum(iBat) + '.pkl')
                else:
                    ec=gu.ipickle(meta['Paths'][pNam]['Input Scenario'][iScn] + '\\Events_Ens' + FixFileNum(iEns) + '_Bat' + FixFileNum(iBat) + '.pkl')

                # Uncompress event chronology if it has been compressed
                ec=EventChronologyDecompress(meta,pNam,ec,iScn,iEns,iBat)
                for iT in range(tv_full.size):
                    indT=np.where(tv==tv_full[iT])[0]
                    if indT.size==0:
                        continue
                    idDist0=np.squeeze(ec['ID Event Type'][iT,:,:])
                    for iDistWithinYear in range(idDist0.shape[1]):
                        idDist1=idDist0[:,iDistWithinYear]
                        uidDist=np.unique(idDist1)
                        if np.sum(uidDist)==0:
                            continue
                        for iU in range(uidDist.size):
                            if uidDist[iU]==0:
                                continue
                            namDist=lut_n2s(meta['LUT']['Event'],uidDist[iU])[0]
                            indDist=np.where(idDist1==uidDist[iU])[0]

                            Area0_Full[namDist][indT,indBat[indDist]]=Area0_Full[namDist][indT,indBat[indDist]]+1

                # Pull results for subset of stands
                for k in meta['LUT']['Event'].keys():
                    Area0[k][:,indBat[iKeepStands]]=Area0_Full[k][:,indBat[iKeepStands]]

            del ec
            garc.collect()

            # Summarize by project type, region, and time
            for iPS in range(uPS_size):
                for iSS in range(uSS_size):
                    for iYS in range(uYS_size):
                        if (meta[pNam]['Project']['Strata']['Project Type']['Unique CD'][iPS]=='All') & (meta[pNam]['Project']['Strata']['Spatial']['Unique CD'][iSS]=='All') & (meta[pNam]['Project']['Strata']['Year']['Unique CD'][iYS]=='All'):
                            # Index to all values for "All"
                            ind1=np.where(meta[pNam]['Project']['Strata']['Project Type']['ID']!=77777)[0]
                        elif (meta[pNam]['Project']['Strata']['Project Type']['Unique CD'][iPS]!='All') & (meta[pNam]['Project']['Strata']['Spatial']['Unique CD'][iSS]=='All') & (meta[pNam]['Project']['Strata']['Year']['Unique CD'][iYS]=='All'):
                            # Index to specific project type stratum
                            ind1=np.where( (meta[pNam]['Project']['Strata']['Project Type']['ID']==meta[pNam]['Project']['Strata']['Project Type']['Unique ID'][iPS]) )[0]
                        elif (meta[pNam]['Project']['Strata']['Project Type']['Unique CD'][iPS]=='All') & (meta[pNam]['Project']['Strata']['Spatial']['Unique CD'][iSS]!='All') & (meta[pNam]['Project']['Strata']['Year']['Unique CD'][iYS]=='All'):
                            # Index to specific spatial stratum
                            ind1=np.where( (meta[pNam]['Project']['Strata']['Spatial']['ID']==meta[pNam]['Project']['Strata']['Spatial']['Unique ID'][iSS]) )[0]
                        elif (meta[pNam]['Project']['Strata']['Project Type']['Unique CD'][iPS]=='All') & (meta[pNam]['Project']['Strata']['Spatial']['Unique CD'][iSS]=='All') & (meta[pNam]['Project']['Strata']['Year']['Unique CD'][iYS]!='All'):
                            # Index to specific year stratum
                            ind1=np.where( (meta[pNam]['Project']['Strata']['Year']['ID']==meta[pNam]['Project']['Strata']['Year']['Unique ID'][iSS]) )[0]
                        elif (meta[pNam]['Project']['Strata']['Project Type']['Unique CD'][iPS]!='All') & (meta[pNam]['Project']['Strata']['Year']['Unique CD'][iYS]=='All') | (meta[pNam]['Project']['Strata']['Spatial']['Unique CD'][iSS]!='All') & (meta[pNam]['Project']['Strata']['Year']['Unique CD'][iYS]=='All'):
                            # Index to specific project type and spatial stratum
                            ind1=np.where( (meta[pNam]['Project']['Strata']['Project Type']['ID']==meta[pNam]['Project']['Strata']['Project Type']['Unique ID'][iPS]) & \
                                      (meta[pNam]['Project']['Strata']['Spatial']['ID']==meta[pNam]['Project']['Strata']['Spatial']['Unique ID'][iSS]) )[0]
                        elif (meta[pNam]['Project']['Strata']['Project Type']['Unique CD'][iPS]!='All') & (meta[pNam]['Project']['Strata']['Year']['Unique CD'][iYS]!='All') & (meta[pNam]['Project']['Strata']['Spatial']['Unique CD'][iSS]=='All'):
                            # Index to specific project type and year stratum
                            ind1=np.where( (meta[pNam]['Project']['Strata']['Project Type']['ID']==meta[pNam]['Project']['Strata']['Project Type']['Unique ID'][iPS]) & \
                                      (meta[pNam]['Project']['Strata']['Year']['ID']==meta[pNam]['Project']['Strata']['Year']['Unique ID'][iYS]) )[0]
                        elif (meta[pNam]['Project']['Strata']['Project Type']['Unique CD'][iPS]=='All') & (meta[pNam]['Project']['Strata']['Year']['Unique CD'][iYS]!='All') & (meta[pNam]['Project']['Strata']['Spatial']['Unique CD'][iSS]!='All'):
                            # Index to specific spatial and year stratum
                            ind1=np.where( (meta[pNam]['Project']['Strata']['Spatial']['ID']==meta[pNam]['Project']['Strata']['Spatial']['Unique ID'][iSS]) & \
                                      (meta[pNam]['Project']['Strata']['Year']['ID']==meta[pNam]['Project']['Strata']['Year']['Unique ID'][iYS]) )[0]
                        else:
                            continue

                        for k in meta['LUT']['Event'].keys():
                            Area1['Mean'][k][:,iEns,iPS,iSS,iYS]=np.mean(Area0[k][:,ind1],axis=1)
                            Area1['Sum'][k][:,iEns,iPS,iSS,iYS]=np.sum(Area0[k][:,ind1],axis=1)

        # Calculate statistics
        Area2={}
        for oper in ['Mean','Sum']:
            Area2[oper]={}
            for k in Area1[oper].keys():
                Area2[oper][k]={}
                Area2[oper][k]['Ensemble Mean']=np.mean(Area1[oper][k],axis=1)
                Area2[oper][k]['Ensemble SD']=np.std(Area1[oper][k],axis=1)
                Area2[oper][k]['Ensemble SE']=np.std(Area1[oper][k],axis=1)/np.sqrt(meta[pNam]['Project']['N Ensemble'])
                Area2[oper][k]['Ensemble P005']=np.percentile(Area1[oper][k],0.5,axis=1)
                Area2[oper][k]['Ensemble P025']=np.percentile(Area1[oper][k],2.5,axis=1)
                Area2[oper][k]['Ensemble P250']=np.percentile(Area1[oper][k],25,axis=1)
                Area2[oper][k]['Ensemble P750']=np.percentile(Area1[oper][k],75,axis=1)
                Area2[oper][k]['Ensemble P975']=np.percentile(Area1[oper][k],97.5,axis=1)
                Area2[oper][k]['Ensemble P995']=np.percentile(Area1[oper][k],99.5,axis=1)

        # Save
        gu.opickle(meta['Paths'][pNam]['Data'] + '\\Outputs\\MOS_ByStrata_Area_Scn' + str(iScn+1) + '.pkl',Area2)

        del Area1,Area2
        garc.collect()

    t1=time.time()
    print(str((t1-t0)/60) + ' min')
    return

#%% Calculate mortality by agent (from points)

def Calc_MOS_FromPoints_MortalityByAgent(meta,pNam,**kwargs):

    t0=time.time()
    print('Calculating model output statistics for mortality by agent')

    # Key word arguments

    if 'StandsToInclude' in kwargs.keys():
        flag_stands_to_include=kwargs['StandsToInclude']
    else:
        flag_stands_to_include=[]

    if 'Scenarios' in kwargs.keys():
        scnL=kwargs['Scenarios']
    else:
        scnL=list(np.arange(0,meta[pNam]['Project']['N Scenario'],1).astype(int))

    # Time series of saved results
    tv_saving=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)

    # All non-economic values
    d0=LoadSingleOutputFile(meta,pNam,0,0,0)

    # Loop through scenarios
    for iScn0 in range(len(scnL)):

        iScn=scnL[iScn0]

        #----------------------------------------------------------------------
        # Initialize data structures
        #----------------------------------------------------------------------

        uPS_size=meta[pNam]['Project']['Strata']['Project Type']['Unique ID'].size
        uSS_size=meta[pNam]['Project']['Strata']['Spatial']['Unique ID'].size
        uYS_size=meta[pNam]['Project']['Strata']['Year']['Unique ID'].size

        Data1={}
        for oper in ['Mean']:
            Data1[oper]={}
            for k in d0['C_M_ByAgent'].keys():
                Data1[oper][k]=np.zeros( (tv_saving.size,meta[pNam]['Project']['N Ensemble'],uPS_size,uSS_size,uYS_size) ,dtype='float32')
                Data1[oper][k]=np.zeros( (tv_saving.size,meta[pNam]['Project']['N Ensemble'],uPS_size,uSS_size,uYS_size) ,dtype='float32')

        #----------------------------------------------------------------------
        # Loop through ensembles
        #----------------------------------------------------------------------

        for iEns in range(meta[pNam]['Project']['N Ensemble']):

            # Initialize temporary data structure for full simulation
            Data0={}
            for k in d0['C_M_ByAgent'].keys():
                Data0[k]=np.nan*np.empty( (tv_saving.size,meta[pNam]['Project']['N Stand']) ,dtype='float32')

            for iBat in range(meta[pNam]['Project']['N Batch']):

                #print(str(iScn) + ' ' + str(iEns) + ' ' + str(iBat) )

                # Include specific subset of stands
                if len(flag_stands_to_include)==0:
                    iKeepStands=np.arange(0,meta[pNam]['Project']['Batch Size'][iBat],1,dtype=int)
                else:
                    iKeepStands=np.where(flag_stands_to_include[iScn][iEns][iBat]==1)[0]

                # Index to batch
                indBat=IndexToBatch(meta[pNam],iBat)

                d1=LoadSingleOutputFile(meta,pNam,iScn,iEns,iBat)

                for k in d1['C_M_ByAgent'].keys():
                    idx=d1['C_M_ByAgent'][k]['idx']
                    z=np.zeros( (tv_saving.size,indBat.size),dtype='float32')
                    z[idx[0],idx[1]]=meta['Core']['Scale Factor C_M_ByAgent']*d1['C_M_ByAgent'][k]['M'].astype('float32')
                    Data0[k][:,indBat[iKeepStands]]=z[:,iKeepStands].copy()
                    #Data0[k][:,indBat[iKeepStands]]=d1[k][:,iKeepStands].copy()

            del d1
            garc.collect()

            #------------------------------------------------------------------
            # Summarize by project type, region, and time
            #------------------------------------------------------------------
            for iPS in range(uPS_size):
                for iSS in range(uSS_size):
                    for iYS in range(uYS_size):
                        if (meta[pNam]['Project']['Strata']['Project Type']['Unique CD'][iPS]=='All') & (meta[pNam]['Project']['Strata']['Spatial']['Unique CD'][iSS]=='All') & (meta[pNam]['Project']['Strata']['Year']['Unique CD'][iYS]=='All'):
                            # Index to all values for "All"
                            ind1=np.where(meta[pNam]['Project']['Strata']['Project Type']['ID']!=77777)[0]
                        elif (meta[pNam]['Project']['Strata']['Project Type']['Unique CD'][iPS]!='All') & (meta[pNam]['Project']['Strata']['Spatial']['Unique CD'][iSS]=='All') & (meta[pNam]['Project']['Strata']['Year']['Unique CD'][iYS]=='All'):
                            # Index to specific project type stratum
                            ind1=np.where( (meta[pNam]['Project']['Strata']['Project Type']['ID']==meta[pNam]['Project']['Strata']['Project Type']['Unique ID'][iPS]) )[0]
                        elif (meta[pNam]['Project']['Strata']['Project Type']['Unique CD'][iPS]=='All') & (meta[pNam]['Project']['Strata']['Spatial']['Unique CD'][iSS]!='All') & (meta[pNam]['Project']['Strata']['Year']['Unique CD'][iYS]=='All'):
                            # Index to specific spatial stratum
                            ind1=np.where( (meta[pNam]['Project']['Strata']['Spatial']['ID']==meta[pNam]['Project']['Strata']['Spatial']['Unique ID'][iSS]) )[0]
                        elif (meta[pNam]['Project']['Strata']['Project Type']['Unique CD'][iPS]=='All') & (meta[pNam]['Project']['Strata']['Spatial']['Unique CD'][iSS]=='All') & (meta[pNam]['Project']['Strata']['Year']['Unique CD'][iYS]!='All'):
                            # Index to specific year stratum
                            ind1=np.where( (meta[pNam]['Project']['Strata']['Year']['ID']==meta[pNam]['Project']['Strata']['Year']['Unique ID'][iSS]) )[0]
                        elif (meta[pNam]['Project']['Strata']['Project Type']['Unique CD'][iPS]!='All') & (meta[pNam]['Project']['Strata']['Year']['Unique CD'][iYS]=='All') | (meta[pNam]['Project']['Strata']['Spatial']['Unique CD'][iSS]!='All') & (meta[pNam]['Project']['Strata']['Year']['Unique CD'][iYS]=='All'):
                            # Index to specific project type and spatial stratum
                            ind1=np.where( (meta[pNam]['Project']['Strata']['Project Type']['ID']==meta[pNam]['Project']['Strata']['Project Type']['Unique ID'][iPS]) & \
                                      (meta[pNam]['Project']['Strata']['Spatial']['ID']==meta[pNam]['Project']['Strata']['Spatial']['Unique ID'][iSS]) )[0]
                        elif (meta[pNam]['Project']['Strata']['Project Type']['Unique CD'][iPS]!='All') & (meta[pNam]['Project']['Strata']['Year']['Unique CD'][iYS]!='All') & (meta[pNam]['Project']['Strata']['Spatial']['Unique CD'][iSS]=='All'):
                            # Index to specific project type and year stratum
                            ind1=np.where( (meta[pNam]['Project']['Strata']['Project Type']['ID']==meta[pNam]['Project']['Strata']['Project Type']['Unique ID'][iPS]) & \
                                      (meta[pNam]['Project']['Strata']['Year']['ID']==meta[pNam]['Project']['Strata']['Year']['Unique ID'][iYS]) )[0]
                        elif (meta[pNam]['Project']['Strata']['Project Type']['Unique CD'][iPS]=='All') & (meta[pNam]['Project']['Strata']['Year']['Unique CD'][iYS]!='All') & (meta[pNam]['Project']['Strata']['Spatial']['Unique CD'][iSS]!='All'):
                            # Index to specific spatial and year stratum
                            ind1=np.where( (meta[pNam]['Project']['Strata']['Spatial']['ID']==meta[pNam]['Project']['Strata']['Spatial']['Unique ID'][iSS]) & \
                                      (meta[pNam]['Project']['Strata']['Year']['ID']==meta[pNam]['Project']['Strata']['Year']['Unique ID'][iYS]) )[0]
                        else:
                            continue

                        for k in Data1['Mean'].keys():
                            Data1['Mean'][k][:,iEns,iPS,iSS,iYS]=np.mean(Data0[k][:,ind1],axis=1)
                            #Data1['Sum'][k][:,iEns,iPS,iSS,iYS]=np.sum(Data0[k][:,ind1],axis=1)

        # Save
        gu.opickle(meta['Paths'][pNam]['Data'] + '\\Outputs\\MOS_ByStrata_MortByAgent_Scn' + str(iScn+1) + '.pkl',Data1)

    t1=time.time()
    print(str((t1-t0)/60) + ' min')

    return

#%% Calculate MOS for mortality from points

# *** This is not by project type and region because cbrunner is not currently
# set up to save mortality by stand (too much data) ***

def Calc_MOS_FromPoints_Mortality_WhenStandsCombined(meta,pNam):

    print('Calculating model output statistics for mortality')

    # Time series of saved results
    tv_saving=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)

    # Initialize list
    mos=[None]*meta[pNam]['Project']['N Scenario']

    # Loop through scenarios
    for iScn in range(meta[pNam]['Project']['N Scenario']):

        #--------------------------------------------------------------------------
        # Initialize data structures
        #--------------------------------------------------------------------------

        d0=LoadSingleOutputFile(meta,pNam,0,0,0)

        mos[iScn]={}
        for k in d0['C_M_ByAgent'].keys():
            mos[iScn][k]=np.zeros((tv_saving.size,meta[pNam]['Project']['N Ensemble']))

        # An option to skip economic calculations. The file will still be created, but all zeros
        #if meta[pNam]['Project']['Skip Economics']=='Off':

        # Loop through ensembles
        for iEns in range(meta[pNam]['Project']['N Ensemble']):

            # Initialize temporary data structure for full simulation
            Data={}
            for k in d0['C_M_ByAgent'].keys():
                Data[k]=np.zeros(tv_saving.size)

            # Loop through batches and add to Data structure
            for iBat in range(meta[pNam]['Project']['N Batch']):
                d1=LoadSingleOutputFile(meta,pNam,iScn,iEns,iBat)
                for k in d0['C_M_ByAgent'].keys():
                    Data[k]=Data[k]+d1['C_M_ByAgent'][k].flatten()

                del d1
                garc.collect()

            # Divide by N stand to get mean
            for k in d0['C_M_ByAgent'].keys():
                Data[k]=Data[k]/meta[pNam]['Project']['N Stand']

            # Populate ensembles
            for k in d0['C_M_ByAgent'].keys():
                mos[iScn][k][:,iEns]=Data[k].copy()

        # Average all ensembles
        for k in d0['C_M_ByAgent'].keys():
            mos[iScn][k]=np.nanmean(mos[iScn][k],axis=1)

    # Save
    gu.opickle(meta['Paths'][pNam]['Data'] + '\\Outputs\\MOS_ByAgent_Mortality.pkl',mos)

    return

#%%

def Calc_AgeClassDistribution(meta,pNam,acd):

    # Time series of saved results
    tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)

    acd['Data']=[None]*meta[pNam]['Project']['N Scenario']
    for iScn in range(meta[pNam]['Project']['N Scenario']):
        acd['Data'][iScn]=np.zeros((acd['binT'].size,acd['binA'].size))
        for iEns in range(meta[pNam]['Project']['N Ensemble']):
            A=np.zeros( (acd['binT'].size,meta[pNam]['Project']['N Stand']) )
            for iBat in range(meta[pNam]['Project']['N Batch']):
                indBat=IndexToBatch(meta[pNam],iBat)
                d1=LoadSingleOutputFile(meta,pNam,iScn,iEns,iBat)
                for iT in range(acd['binT'].size):
                    indT=np.where(tv==acd['binT'][iT])[0]
                    A[iT,indBat]=d1['A'][indT,:]
            del d1
            garc.collect()
            for iA in range(acd['binA'].size):
                for iT in range(acd['binT'].size):
                    ind=np.where( np.abs(A[iT,:]-acd['binA'][iA])<=acd['bwA']/2 )[0]
                    acd['Data'][iScn][iT,iA]=ind.size
        acd['Data'][iScn]=acd['Data'][iScn]/meta[pNam]['Project']['N Ensemble']

    gu.opickle(meta['Paths'][pNam]['Data'] + '\\Outputs\\MOS_Scn' + str(iScn+1) + '_AgeClassDistribution.pkl',acd)

    return acd

#%% Map mean of ensembles for specified time period

def Calc_MOS_MapMean(meta,pNam,iScn,tp,**kwargs):

    tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
    it_ref=np.where( (tv>=tp[0]) & (tv<=tp[1]) )[0]

    mu0={}
    for iEns in range(meta[pNam]['Project']['N Ensemble']):
        print(iEns)
        for iBat in range(meta[pNam]['Project']['N Batch']):
            indBat=IndexToBatch(meta[pNam],iBat)
            d0=LoadSingleOutputFile(meta,pNam,iScn,iEns,iBat)

            # Initialize
            if (iEns==0) & (iBat==0):

                if 'VariablesToKeep' not in kwargs.keys():
                    v2include=list(d0.keys())
                else:
                    v2include=kwargs['VariablesToKeep']

                for k in v2include:
                    mu0[k]=np.zeros(meta[pNam]['Project']['N Stand'])

            # Populate
            for k in v2include:
                if (type(d0[k])==dict) | (k=='Year'):
                    continue
                mu0[k][indBat]=mu0[k][indBat]+np.mean(d0[k][it_ref,:],axis=0)

    # Divide by number of ensembles to get ensemble mean
    for k in v2include:
        mu0[k]=mu0[k]/meta[pNam]['Project']['N Ensemble']

    return mu0

#%% Import scenario data from points
# *** You can't use the nanmean and nanpercentile - way too slow ***

def Import_MOS_ByScnAndStrata_GHGEcon_FromPoints(meta,pNam):

    mos={}
    mos[pNam]={}

    Data=[None]*meta[pNam]['Project']['N Scenario']

    for iScn in range(meta[pNam]['Project']['N Scenario']):

        Data[iScn]={}

        for oper in ['Mean','Sum']:

            Data[iScn][oper]={}

            #------------------------------------------------------------------
            # GHG
            #------------------------------------------------------------------

            d=gu.ipickle(meta['Paths'][pNam]['Data'] + '\\Outputs\\MOS_ByStrata_GHGB_Scn' + str(iScn+1) + '.pkl')
            for k in d[oper].keys():
                Data[iScn][oper][k]={}
                Data[iScn][oper][k]['Ensemble Mean']=np.mean(d[oper][k],axis=1)
                Data[iScn][oper][k]['Ensemble SD']=np.std(d[oper][k],axis=1)
                Data[iScn][oper][k]['Ensemble SE']=np.std(d[oper][k],axis=1)/np.sqrt(meta[pNam]['Project']['N Ensemble'])
                Data[iScn][oper][k]['Ensemble P005']=np.percentile(d[oper][k],0.5,axis=1)
                Data[iScn][oper][k]['Ensemble P025']=np.percentile(d[oper][k],2.5,axis=1)
                Data[iScn][oper][k]['Ensemble P250']=np.percentile(d[oper][k],25,axis=1)
                Data[iScn][oper][k]['Ensemble P750']=np.percentile(d[oper][k],75,axis=1)
                Data[iScn][oper][k]['Ensemble P975']=np.percentile(d[oper][k],97.5,axis=1)
                Data[iScn][oper][k]['Ensemble P995']=np.percentile(d[oper][k],99.5,axis=1)

            #------------------------------------------------------------------
            # Economics
            #------------------------------------------------------------------

            d=gu.ipickle(meta['Paths'][pNam]['Data'] + '\\Outputs\\MOS_ByStrata_Econ_Scn' + str(iScn+1) + '.pkl')
            for k in d[oper].keys():
                Data[iScn][oper][k]={}
                Data[iScn][oper][k]['Ensemble Mean']=np.mean(d[oper][k],axis=1)
                Data[iScn][oper][k]['Ensemble SD']=np.std(d[oper][k],axis=1)
                Data[iScn][oper][k]['Ensemble SE']=np.std(d[oper][k],axis=1)/np.sqrt(meta[pNam]['Project']['N Ensemble'])
                Data[iScn][oper][k]['Ensemble P005']=np.percentile(d[oper][k],0.5,axis=1)
                Data[iScn][oper][k]['Ensemble P025']=np.percentile(d[oper][k],2.5,axis=1)
                Data[iScn][oper][k]['Ensemble P250']=np.percentile(d[oper][k],25,axis=1)
                Data[iScn][oper][k]['Ensemble P750']=np.percentile(d[oper][k],75,axis=1)
                Data[iScn][oper][k]['Ensemble P975']=np.percentile(d[oper][k],97.5,axis=1)
                Data[iScn][oper][k]['Ensemble P995']=np.percentile(d[oper][k],99.5,axis=1)

            #------------------------------------------------------------------
            # Mortality by Agent
            #------------------------------------------------------------------

            if oper=='Mean':
                d=gu.ipickle(meta['Paths'][pNam]['Data'] + '\\Outputs\\MOS_ByStrata_MortbyAgent_Scn' + str(iScn+1) + '.pkl')
                for k in d[oper].keys():
                    Data[iScn][oper]['C_M_' + k]={}
                    Data[iScn][oper]['C_M_' + k]['Ensemble Mean']=np.mean(d[oper][k],axis=1)
                    Data[iScn][oper]['C_M_' + k]['Ensemble SD']=np.std(d[oper][k],axis=1)
                    Data[iScn][oper]['C_M_' + k]['Ensemble SE']=np.std(d[oper][k],axis=1)/np.sqrt(meta[pNam]['Project']['N Ensemble'])
                    Data[iScn][oper]['C_M_' + k]['Ensemble P005']=np.percentile(d[oper][k],0.5,axis=1)
                    Data[iScn][oper]['C_M_' + k]['Ensemble P025']=np.percentile(d[oper][k],2.5,axis=1)
                    Data[iScn][oper]['C_M_' + k]['Ensemble P250']=np.percentile(d[oper][k],25,axis=1)
                    Data[iScn][oper]['C_M_' + k]['Ensemble P750']=np.percentile(d[oper][k],75,axis=1)
                    Data[iScn][oper]['C_M_' + k]['Ensemble P975']=np.percentile(d[oper][k],97.5,axis=1)
                    Data[iScn][oper]['C_M_' + k]['Ensemble P995']=np.percentile(d[oper][k],99.5,axis=1)

    # Save to structure
    mos[pNam]['Scenarios']=Data

    return mos

#%% Import future scenario comparison

def Import_MOS_ByScnComparisonAndStrata_FromPoints(meta,pNam,mos):

    for cNam in mos[pNam]['Delta']:

        # Initialize
        mos[pNam]['Delta'][cNam]['ByStrata']={}

        dC={}
        for oper in ['Mean','Sum']:
            dC[oper]={}

        #----------------------------------------------------------------------
        # Add GHG
        #----------------------------------------------------------------------

        dB=gu.ipickle(meta['Paths'][pNam]['Data'] + '\\Outputs\\MOS_ByStrata_GHGB_Scn' + str(mos[pNam]['Delta'][cNam]['iB']+1) + '.pkl')
        dP=gu.ipickle(meta['Paths'][pNam]['Data'] + '\\Outputs\\MOS_ByStrata_GHGB_Scn' + str(mos[pNam]['Delta'][cNam]['iP']+1) + '.pkl')

        for oper in ['Mean','Sum']:
            for k in dB[oper].keys():
                dC[oper][k]={}
                dC[oper][k]['Ensemble Mean']=np.mean(dP[oper][k]-dB[oper][k],axis=1)
                dC[oper][k]['Ensemble SD']=np.std(dP[oper][k]-dB[oper][k],axis=1)
                dC[oper][k]['Ensemble SE']=np.std(dP[oper][k]-dB[oper][k],axis=1)/np.sqrt(meta[pNam]['Project']['N Ensemble'])
                dC[oper][k]['Ensemble P005']=np.percentile(dP[oper][k]-dB[oper][k],0.5,axis=1)
                dC[oper][k]['Ensemble P025']=np.percentile(dP[oper][k]-dB[oper][k],2.5,axis=1)
                dC[oper][k]['Ensemble P250']=np.percentile(dP[oper][k]-dB[oper][k],25,axis=1)
                dC[oper][k]['Ensemble P750']=np.percentile(dP[oper][k]-dB[oper][k],75,axis=1)
                dC[oper][k]['Ensemble P975']=np.percentile(dP[oper][k]-dB[oper][k],97.5,axis=1)
                dC[oper][k]['Ensemble P995']=np.percentile(dP[oper][k]-dB[oper][k],99.5,axis=1)

        #----------------------------------------------------------------------
        # Add Economics
        #----------------------------------------------------------------------

        dB=gu.ipickle(meta['Paths'][pNam]['Data'] + '\\Outputs\\MOS_ByStrata_Econ_Scn' + str(mos[pNam]['Delta'][cNam]['iB']+1) + '.pkl')
        dP=gu.ipickle(meta['Paths'][pNam]['Data'] + '\\Outputs\\MOS_ByStrata_Econ_Scn' + str(mos[pNam]['Delta'][cNam]['iP']+1) + '.pkl')

        for oper in ['Mean','Sum']:
            for k in dB[oper].keys():
                dC[oper][k]={}
                dC[oper][k]['Ensemble Mean']=np.mean(dP[oper][k]-dB[oper][k],axis=1)
                dC[oper][k]['Ensemble SD']=np.std(dP[oper][k]-dB[oper][k],axis=1)
                dC[oper][k]['Ensemble SE']=np.std(dP[oper][k]-dB[oper][k],axis=1)/np.sqrt(meta[pNam]['Project']['N Ensemble'])
                dC[oper][k]['Ensemble P005']=np.percentile(dP[oper][k]-dB[oper][k],0.5,axis=1)
                dC[oper][k]['Ensemble P025']=np.percentile(dP[oper][k]-dB[oper][k],2.5,axis=1)
                dC[oper][k]['Ensemble P250']=np.percentile(dP[oper][k]-dB[oper][k],25,axis=1)
                dC[oper][k]['Ensemble P750']=np.percentile(dP[oper][k]-dB[oper][k],75,axis=1)
                dC[oper][k]['Ensemble P975']=np.percentile(dP[oper][k]-dB[oper][k],97.5,axis=1)
                dC[oper][k]['Ensemble P995']=np.percentile(dP[oper][k]-dB[oper][k],99.5,axis=1)

        # Add to final structure
        mos[pNam]['Delta'][cNam]['ByStrata']=copy.deepcopy(dC)

        del dB,dP,dC
        garc.collect()

    return mos

#%% Import future scenario comparison

def Import_MOS_ByScnComparisonAndStrata_AcrossProjects(meta,mos,sc):

    # Initialize
    mos['Delta'][sc]['ByStrata']={}

    dC={}
    for oper in ['Mean','Sum']:
        dC[oper]={}

    # Add GHG
    dB=gu.ipickle(mos['Delta'][sc]['Source 1'] + '\\Outputs\\MOS_ByStrata_GHGB_Scn' + str(mos['Delta'][sc]['iB']+1) + '.pkl')
    dP=gu.ipickle(mos['Delta'][sc]['Source 2'] + '\\Outputs\\MOS_ByStrata_GHGB_Scn' + str(mos['Delta'][sc]['iP']+1) + '.pkl')

    for oper in ['Mean','Sum']:
        for k in dB[oper].keys():
            dC[oper][k]={}
            dC[oper][k]['Ensemble Mean']=np.mean(dP[oper][k]-dB[oper][k],axis=1)
            dC[oper][k]['Ensemble SD']=np.std(dP[oper][k]-dB[oper][k],axis=1)
            dC[oper][k]['Ensemble SE']=np.std(dP[oper][k]-dB[oper][k],axis=1)/np.sqrt(meta[pNam]['Project']['N Ensemble'])
            dC[oper][k]['Ensemble P005']=np.percentile(dP[oper][k]-dB[oper][k],0.5,axis=1)
            dC[oper][k]['Ensemble P025']=np.percentile(dP[oper][k]-dB[oper][k],2.5,axis=1)
            dC[oper][k]['Ensemble P250']=np.percentile(dP[oper][k]-dB[oper][k],25,axis=1)
            dC[oper][k]['Ensemble P750']=np.percentile(dP[oper][k]-dB[oper][k],75,axis=1)
            dC[oper][k]['Ensemble P975']=np.percentile(dP[oper][k]-dB[oper][k],97.5,axis=1)
            dC[oper][k]['Ensemble P995']=np.percentile(dP[oper][k]-dB[oper][k],99.5,axis=1)

    # Add Economics

    dB=gu.ipickle(mos['Delta'][sc]['Source 1'] + '\\Outputs\\MOS_ByStrata_Econ_Scn' + str(mos['Delta'][sc]['iB']+1) + '.pkl')
    dP=gu.ipickle(mos['Delta'][sc]['Source 2'] + '\\Outputs\\MOS_ByStrata_Econ_Scn' + str(mos['Delta'][sc]['iP']+1) + '.pkl')

    for oper in ['Mean','Sum']:
        for k in dB[oper].keys():
            dC[oper][k]={}
            dC[oper][k]['Ensemble Mean']=np.mean(dP[oper][k]-dB[oper][k],axis=1)
            dC[oper][k]['Ensemble SD']=np.std(dP[oper][k]-dB[oper][k],axis=1)
            dC[oper][k]['Ensemble SE']=np.std(dP[oper][k]-dB[oper][k],axis=1)/np.sqrt(meta[pNam]['Project']['N Ensemble'])
            dC[oper][k]['Ensemble P005']=np.percentile(dP[oper][k]-dB[oper][k],0.5,axis=1)
            dC[oper][k]['Ensemble P025']=np.percentile(dP[oper][k]-dB[oper][k],2.5,axis=1)
            dC[oper][k]['Ensemble P250']=np.percentile(dP[oper][k]-dB[oper][k],25,axis=1)
            dC[oper][k]['Ensemble P750']=np.percentile(dP[oper][k]-dB[oper][k],75,axis=1)
            dC[oper][k]['Ensemble P975']=np.percentile(dP[oper][k]-dB[oper][k],97.5,axis=1)
            dC[oper][k]['Ensemble P995']=np.percentile(dP[oper][k]-dB[oper][k],99.5,axis=1)

    # Add to final structure
    mos['Delta'][sc]['ByStrata']=copy.deepcopy(dC)

    del dB,dP,dC
    garc.collect()

    return mos

#%% Import model output statistics for area (from points)
# Area isn't considered in scenario comparisons so the model output statistics
# have already been calculated and just need to be imported.

def Import_MOS_ByScnAndStrata_Area_FromPoints(meta,pNam,mos):

    for iScn in range(meta[pNam]['Project']['N Scenario']):

        d=gu.ipickle(meta['Paths'][pNam]['Data'] + '\\Outputs\\MOS_ByStrata_Area_Scn' + str(iScn+1) + '.pkl')

        for oper in ['Mean','Sum']:
            for k in d[oper].keys():
                mos[pNam]['Scenarios'][iScn][oper]['Area_' + k]=d[oper][k]

    return mos

# #%% Save output variables by multipolygon

# # This is only designed for projects that apply inventory_from_polygons method.
# # Not only is this outputting by project, but the project flux sums are accurately
# # representing the given treatment area, rather than the area inferred from the
# # total number of sparse sample points within a project area, which can be wrong.

# # *** Use this if you want to keep a time series for each MP (e.g. LCELF). If
# # you are content with summaries by project type and region, then use other
# # functions below. ***

# def MosByMultipolygon(meta,switch_area,switch_cashflow):

#     # Import multipolygons
#     atu_multipolygons=gu.ipickle(meta['Paths']['Geospatial'] + '\\atu_multipolygons.pkl')

#     # Import sxy
#     geos=gu.ipickle(meta['Paths']['Geospatial'] + '\\geos.pkl')

#     # Unique MPs
#     uMP=np.unique(geos['Sparse']['ID_atu_multipolygons'])

#     # Create listed index (faster than indexing on the fly)
#     Crosswalk_sxy_to_mp=[None]*uMP.size
#     for iMP in range(uMP.size):
#         d={}
#         d['Index']=np.where(geos['Sparse']['ID_atu_multipolygons']==uMP[iMP])[0]
#         Crosswalk_sxy_to_mp[iMP]=d

#     # Time series of saved results
#     tv_full=np.arange(meta[pNam]['Project']['Year Start'],meta[pNam]['Project']['Year End']+1,1)
#     tv_saving=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
#     it=np.where( (tv_full>=tv_saving[0]) & (tv_full<=tv_saving[-1]) )[0]

#     # Variables to save
#     d1=LoadSingleOutputFile(meta,0,0,0)
#     nam1=list(d1.keys())
#     nam1.remove('Year')
#     nam1.remove('C_M_ByAgent')

#     #nam1=['A','V_StemMerch','C_Biomass_Tot','C_DeadWood_Tot','C_Litter_Tot','C_Soil_Tot', \
#     #      'C_InUse_Tot','C_DumpLandfill_Tot','C_ToMillMerch','C_ToMillNonMerch','C_ToMillSnagStem', \
#     #      'C_ToLumber','C_ToPlywood','C_ToOSB','C_ToMDF','C_ToPaper','C_ToPowerFacilityDom','C_ToPowerGrid','C_ToPellets','C_ToFirewoodDom','C_ToLogExport', \
#     #      'C_G_Net_Tot','C_LF_Tot','C_RH_Tot','E_CO2e_LULUCF_NEE','E_CO2e_LULUCF_Wildfire','E_CO2e_LULUCF_OpenBurning','CO2e_LULUCF_E_EcoOther', \
#     #      'CO2e_LULUCF_E_HWP','E_CO2e_ESC_Comb','E_CO2e_ESC_SubE','E_CO2e_ESC_SubBM','E_CO2e_ET_Comb','E_CO2e_IPPU_Comb','E_CO2e_LULUCF_Fire','E_CO2e_AGHGB']

#     nam_cashflow=['Cost Total','Revenue Gross']

#     # Scale factor used to temporarily store data
#     ScaleFactor=0.001

#     #--------------------------------------------------------------------------
#     # Initialize data by multipolygon structure
#     #--------------------------------------------------------------------------

#     MosByMP=[None]*meta[pNam]['Project']['N Scenario']
#     for iScn in range(meta[pNam]['Project']['N Scenario']):

#         d={}

#         d['v1']={}
#         d['v1']['Mean']={}
#         d['v1']['Sum']={}
#         for iV in range(len(nam1)):
#             d['v1']['Mean'][nam1[iV]]={}
#             d['v1']['Mean'][nam1[iV]]['Ensemble Mean']=np.zeros((tv_saving.size,uMP.size))
#             d['v1']['Mean'][nam1[iV]]['Ensemble SD']=np.zeros((tv_saving.size,uMP.size))
#             d['v1']['Mean'][nam1[iV]]['Ensemble P10']=np.zeros((tv_saving.size,uMP.size))
#             d['v1']['Mean'][nam1[iV]]['Ensemble P25']=np.zeros((tv_saving.size,uMP.size))
#             d['v1']['Mean'][nam1[iV]]['Ensemble P75']=np.zeros((tv_saving.size,uMP.size))
#             d['v1']['Mean'][nam1[iV]]['Ensemble P90']=np.zeros((tv_saving.size,uMP.size))
#             d['v1']['Sum'][nam1[iV]]={}
#             d['v1']['Sum'][nam1[iV]]['Ensemble Mean']=np.zeros((tv_saving.size,uMP.size))
#             d['v1']['Sum'][nam1[iV]]['Ensemble SD']=np.zeros((tv_saving.size,uMP.size))
#             d['v1']['Sum'][nam1[iV]]['Ensemble P10']=np.zeros((tv_saving.size,uMP.size))
#             d['v1']['Sum'][nam1[iV]]['Ensemble P25']=np.zeros((tv_saving.size,uMP.size))
#             d['v1']['Sum'][nam1[iV]]['Ensemble P75']=np.zeros((tv_saving.size,uMP.size))
#             d['v1']['Sum'][nam1[iV]]['Ensemble P90']=np.zeros((tv_saving.size,uMP.size))

#         if switch_area=='On':
#             d['Area']={}
#             for k in meta['LUT']['Event'].keys():
#                 d['Area'][k]={}
#                 d['Area'][k]['Ensemble Mean']=np.zeros((tv_saving.size,uMP.size))
#                 #d['Area'][k]['Ensemble SD']=np.zeros((tv_saving.size,uMP.size))

#         if switch_cashflow=='On':
#             d['Cashflow']={}
#             d['Cashflow']['Mean']={}
#             d['Cashflow']['Sum']={}
#             for iV in range(len(nam_cashflow)):
#                 d['Cashflow']['Mean'][nam_cashflow[iV]]={}
#                 d['Cashflow']['Mean'][nam_cashflow[iV]]['Ensemble Mean']=np.zeros((tv_saving.size,uMP.size))
#                 d['Cashflow']['Mean'][nam_cashflow[iV]]['Ensemble SD']=np.zeros((tv_saving.size,uMP.size))
#                 #d['Cashflow']['Mean'][nam_cashflow[iV]]['Ensemble P10']=np.zeros((tv_saving.size,uMP.size))
#                 #d['Cashflow']['Mean'][nam_cashflow[iV]]['Ensemble P25']=np.zeros((tv_saving.size,uMP.size))
#                 #d['Cashflow']['Mean'][nam_cashflow[iV]]['Ensemble P75']=np.zeros((tv_saving.size,uMP.size))
#                 #d['Cashflow']['Mean'][nam_cashflow[iV]]['Ensemble P90']=np.zeros((tv_saving.size,uMP.size))
#                 d['Cashflow']['Sum'][nam_cashflow[iV]]={}
#                 d['Cashflow']['Sum'][nam_cashflow[iV]]['Ensemble Mean']=np.zeros((tv_saving.size,uMP.size))
#                 d['Cashflow']['Sum'][nam_cashflow[iV]]['Ensemble SD']=np.zeros((tv_saving.size,uMP.size))
#                 d['Cashflow']['Sum'][nam_cashflow[iV]]['Ensemble P10']=np.zeros((tv_saving.size,uMP.size))
#                 d['Cashflow']['Sum'][nam_cashflow[iV]]['Ensemble P25']=np.zeros((tv_saving.size,uMP.size))
#                 d['Cashflow']['Sum'][nam_cashflow[iV]]['Ensemble P75']=np.zeros((tv_saving.size,uMP.size))
#                 d['Cashflow']['Sum'][nam_cashflow[iV]]['Ensemble P90']=np.zeros((tv_saving.size,uMP.size))

#         MosByMP[iScn]=d

#     #--------------------------------------------------------------------------
#     # Loop through scenarios
#     #--------------------------------------------------------------------------

#     for iScn in range(meta[pNam]['Project']['N Scenario']):

#         for iEns in range(meta[pNam]['Project']['N Ensemble']):

#             #------------------------------------------------------------------
#             # Initialize temporary data structure for full simulation
#             #------------------------------------------------------------------

#             Data={}

#             Data['v1']={}
#             for iV in range(len(nam1)):
#                 Data['v1'][nam1[iV]]=np.zeros((tv_saving.size,meta[pNam]['Project']['N Stand']),dtype=int)

#             if switch_area=='On':
#                 Data['Area']={}
#                 for k in MosByMP[iScn]['Area']:
#                     Data['Area'][k]=np.zeros((tv_saving.size,meta[pNam]['Project']['N Stand']),dtype=int)

#             if switch_cashflow=='On':
#                 Data['Cashflow']={}
#                 for iV in range(len(nam_cashflow)):
#                     nam=nam_cashflow[iV]
#                     Data['Cashflow'][nam]=np.zeros((tv_saving.size,meta[pNam]['Project']['N Stand']),dtype=int)

#             #------------------------------------------------------------------
#             # Populate full simulation results
#             #------------------------------------------------------------------

#             for iBat in range(meta[pNam]['Project']['N Batch']):

#                 indBat=IndexToBatch(meta[pNam],iBat)

#                 d1=LoadSingleOutputFile(meta,pNam,iScn,iEns,iBat)

#                 for iV in range(len(nam1)):
#                     tmp=d1[nam1[iV]]/ScaleFactor
#                     Data['v1'][nam1[iV]][:,indBat]=tmp.copy().astype(int)

#                 if (switch_area=='On' ) | (switch_cashflow=='On'):

#                     # Import event chronology
#                     if (meta[pNam]['Scenario'][iScn]['Harvest Status Future']=='On') | (meta[pNam]['Scenario'][iScn]['Breakup Status Historical']=='On') | (meta[pNam]['Scenario'][iScn]['Breakup Status Future']=='On'):
#                         ec=gu.ipickle(meta['Paths'][pNam]['Input Scenario'][iScn] + '\\Modified_Events_Ens' + FixFileNum(iEns) + '_Bat' + FixFileNum(iBat) + '.pkl')
#                     else:
#                         ec=gu.ipickle(meta['Paths'][pNam]['Input Scenario'][iScn] + '\\Events_Ens' + FixFileNum(iEns) + '_Bat' + FixFileNum(iBat) + '.pkl')

#                     # Uncompress event chronology if it has been compressed
#                     ec=EventChronologyDecompress(meta,pNam,ec,iScn,iEns,iBat)

#                 if switch_area=='On':

#                     for k in Data['Area'].keys():
#                         Data0=np.zeros((tv_saving.size,indBat.size))
#                         for iEY in range(meta['Core']['Max Events Per Year']):
#                             ind=np.where(ec['ID Event Type'][it,:,iEY]==meta['LUT']['Event'][k])[0]
#                             Data0[ind]=Data0[ind]+1
#                         Data['Area'][k][:,indBat]=Data0

#                 if switch_cashflow=='On':

#                     inv=gu.ipickle(meta['Paths'][pNam]['Input Scenario'][iScn] + '\\Inventory_Bat' + FixFileNum(iBat) + '.pkl')

#                     econ=econo.CalculateNetRevenue(meta,iScn,iEns,iBat,inv,ec,d1)

#                     for iV in range(len(nam_cashflow)):
#                         tmp=econ[nam_cashflow[iV]]/ScaleFactor
#                         Data['Cashflow'][nam_cashflow[iV]][:,indBat]=tmp.copy().astype(int)

#                     del econ

#                 del d1
#                 garc.collect()

#             if switch_area=='On':
#                 # Populating the final structure with area data is slow - get a flag
#                 # indicator of whether each event ID can be skipped because it has no
#                 # info
#                 flag_Area={}
#                 for k in Data['Area'].keys():
#                     if np.sum(Data['Area'][k])>0:
#                         flag_Area[k]=1
#                     else:
#                         flag_Area[k]=0

#             #------------------------------------------------------------------
#             # Calculate stats and populate results for each treatment area
#             #------------------------------------------------------------------

#             for iMP in range(uMP.size):

#                 ATA=atu_multipolygons[uMP[iMP]]['ACTUAL_TREATMENT_AREA']
#                 if ATA==None:
#                     print('Encounterd no area, using zero')
#                     ATA=0

#                 ind=Crosswalk_sxy_to_mp[iMP]['Index']

#                 for iV in range(len(nam1)):
#                     tmp=ScaleFactor*Data['v1'][nam1[iV]][:,ind].astype(float)
#                     MosByMP[iScn]['v1']['Mean'][nam1[iV]]['Ensemble Mean'][:,iMP]=MosByMP[iScn]['v1']['Mean'][nam1[iV]]['Ensemble Mean'][:,iMP]+np.mean(tmp,axis=1)
#                     MosByMP[iScn]['v1']['Mean'][nam1[iV]]['Ensemble SD'][:,iMP]=MosByMP[iScn]['v1']['Mean'][nam1[iV]]['Ensemble SD'][:,iMP]+np.std(tmp,axis=1)
#                     MosByMP[iScn]['v1']['Sum'][nam1[iV]]['Ensemble Mean'][:,iMP]=MosByMP[iScn]['v1']['Sum'][nam1[iV]]['Ensemble Mean'][:,iMP]+ATA*np.mean(tmp,axis=1)
#                     MosByMP[iScn]['v1']['Sum'][nam1[iV]]['Ensemble SD'][:,iMP]=MosByMP[iScn]['v1']['Sum'][nam1[iV]]['Ensemble SD'][:,iMP]+ATA*np.std(tmp,axis=1)
#                     MosByMP[iScn]['v1']['Sum'][nam1[iV]]['Ensemble P10'][:,iMP]=MosByMP[iScn]['v1']['Sum'][nam1[iV]]['Ensemble P10'][:,iMP]+ATA*np.percentile(tmp,10,axis=1)
#                     MosByMP[iScn]['v1']['Sum'][nam1[iV]]['Ensemble P25'][:,iMP]=MosByMP[iScn]['v1']['Sum'][nam1[iV]]['Ensemble P25'][:,iMP]+ATA*np.percentile(tmp,25,axis=1)
#                     MosByMP[iScn]['v1']['Sum'][nam1[iV]]['Ensemble P75'][:,iMP]=MosByMP[iScn]['v1']['Sum'][nam1[iV]]['Ensemble P75'][:,iMP]+ATA*np.percentile(tmp,75,axis=1)
#                     MosByMP[iScn]['v1']['Sum'][nam1[iV]]['Ensemble P90'][:,iMP]=MosByMP[iScn]['v1']['Sum'][nam1[iV]]['Ensemble P90'][:,iMP]+ATA*np.percentile(tmp,90,axis=1)

#                 if switch_cashflow=='On':
#                     for iV in range(len(nam_cashflow)):
#                         tmp=ScaleFactor*Data['Cashflow'][nam_cashflow[iV]][:,ind].astype(float)
#                         MosByMP[iScn]['Cashflow']['Mean'][nam_cashflow[iV]]['Ensemble Mean'][:,iMP]=MosByMP[iScn]['Cashflow']['Mean'][nam_cashflow[iV]]['Ensemble Mean'][:,iMP]+np.mean(tmp,axis=1)
#                         MosByMP[iScn]['Cashflow']['Mean'][nam_cashflow[iV]]['Ensemble SD'][:,iMP]=MosByMP[iScn]['Cashflow']['Mean'][nam_cashflow[iV]]['Ensemble SD'][:,iMP]+np.std(tmp,axis=1)
#                         MosByMP[iScn]['Cashflow']['Sum'][nam_cashflow[iV]]['Ensemble Mean'][:,iMP]=MosByMP[iScn]['Cashflow']['Sum'][nam_cashflow[iV]]['Ensemble Mean'][:,iMP]+ATA*np.mean(tmp,axis=1)
#                         MosByMP[iScn]['Cashflow']['Sum'][nam_cashflow[iV]]['Ensemble SD'][:,iMP]=MosByMP[iScn]['Cashflow']['Sum'][nam_cashflow[iV]]['Ensemble SD'][:,iMP]+ATA*np.std(tmp,axis=1)
#                         MosByMP[iScn]['Cashflow']['Sum'][nam_cashflow[iV]]['Ensemble P10'][:,iMP]=MosByMP[iScn]['Cashflow']['Sum'][nam_cashflow[iV]]['Ensemble P10'][:,iMP]+ATA*np.percentile(tmp,10,axis=1)
#                         MosByMP[iScn]['Cashflow']['Sum'][nam_cashflow[iV]]['Ensemble P25'][:,iMP]=MosByMP[iScn]['Cashflow']['Sum'][nam_cashflow[iV]]['Ensemble P25'][:,iMP]+ATA*np.percentile(tmp,25,axis=1)
#                         MosByMP[iScn]['Cashflow']['Sum'][nam_cashflow[iV]]['Ensemble P75'][:,iMP]=MosByMP[iScn]['Cashflow']['Sum'][nam_cashflow[iV]]['Ensemble P75'][:,iMP]+ATA*np.percentile(tmp,75,axis=1)
#                         MosByMP[iScn]['Cashflow']['Sum'][nam_cashflow[iV]]['Ensemble P90'][:,iMP]=MosByMP[iScn]['Cashflow']['Sum'][nam_cashflow[iV]]['Ensemble P90'][:,iMP]+ATA*np.percentile(tmp,90,axis=1)

#                 if switch_area=='On':
#                     for k in MosByMP[iScn]['Area']:
#                         if flag_Area[k]==1:
#                             # Only continue if there are some events
#                             MosByMP[iScn]['Area'][k]['Ensemble Mean'][:,iMP]=MosByMP[iScn]['Area'][k]['Ensemble Mean'][:,iMP]+np.sum(Data['Area'][k][:,ind],axis=1)

#     #--------------------------------------------------------------------------
#     # Divide by number of ensembles
#     #--------------------------------------------------------------------------

#     for iScn in range(meta[pNam]['Project']['N Scenario']):

#         for iV in range(len(nam1)):
#             MosByMP[iScn]['v1']['Mean'][nam1[iV]]['Ensemble Mean']=MosByMP[iScn]['v1']['Mean'][nam1[iV]]['Ensemble Mean']/meta[pNam]['Project']['N Ensemble']
#             MosByMP[iScn]['v1']['Mean'][nam1[iV]]['Ensemble SD']=MosByMP[iScn]['v1']['Mean'][nam1[iV]]['Ensemble SD']/meta[pNam]['Project']['N Ensemble']
#             MosByMP[iScn]['v1']['Sum'][nam1[iV]]['Ensemble Mean']=MosByMP[iScn]['v1']['Sum'][nam1[iV]]['Ensemble Mean']/meta[pNam]['Project']['N Ensemble']
#             MosByMP[iScn]['v1']['Sum'][nam1[iV]]['Ensemble SD']=MosByMP[iScn]['v1']['Sum'][nam1[iV]]['Ensemble SD']/meta[pNam]['Project']['N Ensemble']
#             MosByMP[iScn]['v1']['Sum'][nam1[iV]]['Ensemble P10']=MosByMP[iScn]['v1']['Sum'][nam1[iV]]['Ensemble P10']/meta[pNam]['Project']['N Ensemble']
#             MosByMP[iScn]['v1']['Sum'][nam1[iV]]['Ensemble P25']=MosByMP[iScn]['v1']['Sum'][nam1[iV]]['Ensemble P25']/meta[pNam]['Project']['N Ensemble']
#             MosByMP[iScn]['v1']['Sum'][nam1[iV]]['Ensemble P75']=MosByMP[iScn]['v1']['Sum'][nam1[iV]]['Ensemble P75']/meta[pNam]['Project']['N Ensemble']
#             MosByMP[iScn]['v1']['Sum'][nam1[iV]]['Ensemble P90']=MosByMP[iScn]['v1']['Sum'][nam1[iV]]['Ensemble P90']/meta[pNam]['Project']['N Ensemble']

#         if switch_cashflow=='On':
#             for iV in range(len(nam_cashflow)):
#                 MosByMP[iScn]['Cashflow']['Mean'][nam_cashflow[iV]]['Ensemble Mean']=MosByMP[iScn]['Cashflow']['Mean'][nam_cashflow[iV]]['Ensemble Mean']/meta[pNam]['Project']['N Ensemble']
#                 MosByMP[iScn]['Cashflow']['Mean'][nam_cashflow[iV]]['Ensemble SD']=MosByMP[iScn]['Cashflow']['Mean'][nam_cashflow[iV]]['Ensemble SD']/meta[pNam]['Project']['N Ensemble']
#                 MosByMP[iScn]['Cashflow']['Sum'][nam_cashflow[iV]]['Ensemble Mean']=MosByMP[iScn]['Cashflow']['Sum'][nam_cashflow[iV]]['Ensemble Mean']/meta[pNam]['Project']['N Ensemble']
#                 MosByMP[iScn]['Cashflow']['Sum'][nam_cashflow[iV]]['Ensemble SD']=MosByMP[iScn]['Cashflow']['Sum'][nam_cashflow[iV]]['Ensemble SD']/meta[pNam]['Project']['N Ensemble']
#                 MosByMP[iScn]['Cashflow']['Sum'][nam_cashflow[iV]]['Ensemble P10']=MosByMP[iScn]['Cashflow']['Sum'][nam_cashflow[iV]]['Ensemble P10']/meta[pNam]['Project']['N Ensemble']
#                 MosByMP[iScn]['Cashflow']['Sum'][nam_cashflow[iV]]['Ensemble P25']=MosByMP[iScn]['Cashflow']['Sum'][nam_cashflow[iV]]['Ensemble P25']/meta[pNam]['Project']['N Ensemble']
#                 MosByMP[iScn]['Cashflow']['Sum'][nam_cashflow[iV]]['Ensemble P75']=MosByMP[iScn]['Cashflow']['Sum'][nam_cashflow[iV]]['Ensemble P75']/meta[pNam]['Project']['N Ensemble']
#                 MosByMP[iScn]['Cashflow']['Sum'][nam_cashflow[iV]]['Ensemble P90']=MosByMP[iScn]['Cashflow']['Sum'][nam_cashflow[iV]]['Ensemble P90']/meta[pNam]['Project']['N Ensemble']

#         if switch_area=='On':
#             for iV in MosByMP[iScn]['Area']:
#                 MosByMP[iScn]['Area'][iV]['Ensemble Mean'][:,iMP]=MosByMP[iScn]['Area'][iV]['Ensemble Mean'][:,iMP]/meta[pNam]['Project']['N Ensemble']

#     #--------------------------------------------------------------------------
#     # Save
#     #--------------------------------------------------------------------------

#     gu.opickle(meta['Paths'][pNam]['Data'] + '\\Outputs\\MosByMultipolygon.pkl',MosByMP)

#     return

# #%% Summarize GHG for multipolygons

# def MOS_FromMPs_ByProjTypeRegAndYear_GHG(meta):

#     t0=time.time()

#     # Import geos
#     geos=gu.ipickle(meta['Paths']['Geospatial'] + '\\geos.pkl')

#     # Unique MPs
#     uMP=np.unique(geos['Sparse']['ID_atu_multipolygons'])

#     # Unique project types
#     uPT=np.unique(meta[pNam]['Project']['ByMP']['Project Type'])

#     # Unique regions
#     uReg=np.unique(meta[pNam]['Project']['ByMP']['Region ID'])

#     # Create listed index (faster than indexing on the fly)
#     Crosswalk_sxy_to_mp=[None]*uMP.size
#     for iMP in range(uMP.size):
#         d={}
#         d['Index']=np.where(geos['Sparse']['ID_atu_multipolygons']==uMP[iMP])[0]
#         Crosswalk_sxy_to_mp[iMP]=d

#     # Pull out project type and year for unique MPs
#     ProjectType=np.zeros(uMP.size)
#     ProjectArea=np.zeros(uMP.size)
#     ProjectYear=np.zeros(uMP.size)
#     Region=np.zeros(uMP.size)
#     for iMP in range(uMP.size):
#         ProjectType[iMP]=meta[pNam]['Project']['ByMP']['Project Type'][uMP[iMP]]
#         ProjectArea[iMP]=meta[pNam]['Project']['ByMP']['Project Area'][uMP[iMP]]
#         ProjectYear[iMP]=meta[pNam]['Project']['ByMP']['Project Year'][uMP[iMP]]
#         Region[iMP]=meta[pNam]['Project']['ByMP']['Region ID'][uMP[iMP]]

#     # Time series of saved results
#     tv_saving=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)

#     # Years of implementation
#     uT=np.unique(ProjectYear[ProjectYear!=0])
#     tv_Implementation=np.arange(np.min(uT),np.max(uT)+1,1)

#     # All (excluding 'C_M_ByAgent' and 'Year','C_M_Dist')
#     v2include=['A','V_MerchLive','V_MerchDead','V_MerchTotal','V_ToMillMerchLive','V_ToMillMerchDead','V_ToMillMerchTotal','V_ToMillNonMerch',
#                'LogSizeEnhancement','C_Forest_Tot','C_HWP_Tot','C_NPP_Tot','C_ToMill','C_Biomass_Tot','C_Piled_Tot','C_Litter_Tot','C_DeadWood_Tot',
#                'C_Soil_Tot','C_InUse_Tot','C_DumpLandfill_Tot','C_G_Gross_Tot','C_G_Net_Tot','C_M_Reg_Tot','C_LF_Tot','C_RH_Tot',
#                'C_ToMillMerch','C_ToMillNonMerch','C_ToMillSnagStem','C_ToSlashpileBurnTot','C_ToSlashpileBurnNonMerch','C_ToLumber','C_ToPlywood','C_ToOSB','C_ToMDF',
#                'C_ToPaper','C_ToPowerFacilityDom','C_ToPowerFacilityExport','C_ToPowerGrid','C_ToPellets','C_ToFirewoodDom',
#                'C_ToFirewoodExport','C_ToLogExport',
#                'E_CO2e_LULUCF_NEE','E_CO2e_LULUCF_Denit','E_CO2e_LULUCF_Other','E_CO2e_LULUCF_Wildfire','E_CO2e_LULUCF_OpenBurning',
#                'E_CO2e_LULUCF_Fire','E_CO2e_LULUCF_HWP','E_CO2e_ESC_Bioenergy',
#                'E_CO2e_ESC_OperFor','E_CO2e_ET_OperFor','E_CO2e_IPPU_OperFor',
#                'E_CO2e_Coal','E_CO2e_Oil','E_CO2e_Gas',
#                'E_CO2e_SUB_Coal','E_CO2e_SUB_Oil','E_CO2e_SUB_Gas','E_CO2e_SUB_Calcination','E_CO2e_SUB_E','E_CO2e_SUB_M','E_CO2e_SUB_Tot',
#                'E_CO2e_AGHGB_WSub','E_CO2e_AGHGB_WOSub','E_CO2e_AGHGB_WSub_cumu','E_CO2e_AGHGB_WOSub_cumu',
#                'E_CO2e_AGHGB_WSub_cumu_from_tref','E_CO2e_AGHGB_WOSub_cumu_from_tref',
#                'C_Coal','C_Oil','C_Gas','ODT Sawnwood','ODT Panel','ODT Concrete',
#                'ODT Steel','ODT Aluminum','ODT Plastic','ODT Textile']

#     # Scale factor used to temporarily store data
#     #ScaleFactor=1.0

#     # Loop through scenarios
#     for iScn in range(meta[pNam]['Project']['N Scenario']):

#         #--------------------------------------------------------------------------
#         # Initialize data structures
#         #--------------------------------------------------------------------------

#         MoMeanByPT={}
#         MoSumByPT={}
#         MoMeanByPTAndYr={}
#         MoSumByPTAndYr={}
#         for k in v2include:
#             MoMeanByPT[k]=np.zeros( (tv_saving.size,meta[pNam]['Project']['N Ensemble'],uPT.size,uReg.size) ,dtype='float64')
#             MoSumByPT[k]=np.zeros( (tv_saving.size,meta[pNam]['Project']['N Ensemble'],uPT.size,uReg.size) ,dtype='float64')

#             MoMeanByPTAndYr[k]={}
#             MoSumByPTAndYr[k]={}
#             for t in tv_Implementation:
#                 MoMeanByPTAndYr[k][t]=np.zeros( (tv_saving.size,meta[pNam]['Project']['N Ensemble'],uPT.size,uReg.size) ,dtype='float64')
#                 MoSumByPTAndYr[k][t]=np.zeros( (tv_saving.size,meta[pNam]['Project']['N Ensemble'],uPT.size,uReg.size) ,dtype='float64')

#         # Loop through ensembles
#         for iEns in range(meta[pNam]['Project']['N Ensemble']):

#             #--------------------------------------------------------------------------
#             # Initialize temporary data structure for full simulation
#             #--------------------------------------------------------------------------

#             Data0={}
#             for k in v2include:
#                 Data0[k]=np.zeros( (tv_saving.size,meta[pNam]['Project']['N Stand']) ,dtype='float64')

#             #--------------------------------------------------------------------------
#             # Import batches
#             #--------------------------------------------------------------------------

#             for iBat in range(meta[pNam]['Project']['N Batch']):

#                 indBat=IndexToBatch(meta[pNam],iBat)

#                 d1=LoadSingleOutputFile(meta,pNam,iScn,iEns,iBat)

#                 for k in v2include:

#                     Data0[k][:,indBat]=d1[k].copy()

#             del d1
#             garc.collect()

#             #--------------------------------------------------------------------------
#             # Convert from SXY results to MP averages
#             #--------------------------------------------------------------------------

#             MeanByMP={}
#             SumByMP={}
#             for k in v2include:
#                 MeanByMP[k]=np.zeros( (tv_saving.size,uMP.size) )
#                 SumByMP[k]=np.zeros( (tv_saving.size,uMP.size) )

#             for iMP in range(uMP.size):

#                 # Index to SXY for the ith multipolygon
#                 indSXY=Crosswalk_sxy_to_mp[iMP]['Index']

#                 for k in v2include:
#                     MeanByMP[k][:,iMP]=MeanByMP[k][:,iMP]+np.mean(Data0[k][:,indSXY],axis=1)
#                     SumByMP[k][:,iMP]=SumByMP[k][:,iMP]+np.sum(ProjectArea[iMP])*np.mean(Data0[k][:,indSXY],axis=1)

#             #--------------------------------------------------------------------------
#             # Summarize by project type, region, and time
#             #--------------------------------------------------------------------------

#             for iPT in range(uPT.size):

#                 for iReg in range(uReg.size):

#                     ind1=np.where( (ProjectType==uPT[iPT]) & (Region==uReg[iReg]) )[0]

#                     for k in v2include:

#                         MoMeanByPT[k][:,iEns,iPT,iReg]=np.mean(MeanByMP[k][:,ind1],axis=1)
#                         MoSumByPT[k][:,iEns,iPT,iReg]=np.sum(SumByMP[k][:,ind1],axis=1)

#                         # Summmarize by PT and year
#                         for t in tv_Implementation:

#                             ind2=np.where( (ProjectType==uPT[iPT]) & (Region==uReg[iReg]) & (ProjectYear==t) )[0]

#                             MoMeanByPTAndYr[k][t][:,iEns,iPT,iReg]=np.mean(MeanByMP[k][:,ind2],axis=1)
#                             MoSumByPTAndYr[k][t][:,iEns,iPT,iReg]=np.sum(SumByMP[k][:,ind2],axis=1)

#         #--------------------------------------------------------------------------
#         # Save
#         #--------------------------------------------------------------------------

#         gu.opickle(meta['Paths'][pNam]['Data'] + '\\Outputs\\GHGB_ByProjTypeReg_Mean_Scn' + str(iScn+1) + '_FromMPs.pkl',MoMeanByPT)
#         gu.opickle(meta['Paths'][pNam]['Data'] + '\\Outputs\\GHGB_ByProjTypeReg_Sum_Scn' + str(iScn+1) + '_FromMPs.pkl',MoSumByPT)

#         gu.opickle(meta['Paths'][pNam]['Data'] + '\\Outputs\\GHGB_ByProjTypeRegAndYear_Mean_Scn' + str(iScn+1) + '_FromMPs.pkl',MoMeanByPTAndYr)
#         gu.opickle(meta['Paths'][pNam]['Data'] + '\\Outputs\\GHGB_ByProjTypeRegAndYear_Sum_Scn' + str(iScn+1) + '_FromMPs.pkl',MoSumByPTAndYr)

#     t1=time.time()
#     print((t1-t0)/60)

#     return

# #%% Summarize economcis for multipolygons

# def MOS_FromMPs_ByProjTypeRegAndYear_Econ(meta):

#     t0=time.time()

#     # Import geos
#     geos=gu.ipickle(meta['Paths']['Geospatial'] + '\\geos.pkl')

#     # Unique MPs
#     uMP=np.unique(geos['Sparse']['ID_atu_multipolygons'])

#     # Unique project types
#     uPT=np.unique(meta[pNam]['Project']['ByMP']['Project Type'])

#     # Unique regions
#     uReg=np.unique(meta[pNam]['Project']['ByMP']['Region ID'])

#     # Create listed index (faster than indexing on the fly)
#     Crosswalk_sxy_to_mp=[None]*uMP.size
#     for iMP in range(uMP.size):
#         d={}
#         d['Index']=np.where(geos['Sparse']['ID_atu_multipolygons']==uMP[iMP])[0]
#         Crosswalk_sxy_to_mp[iMP]=d

#     # Pull out project type and year for unique MPs
#     ProjectType=np.zeros(uMP.size)
#     ProjectArea=np.zeros(uMP.size)
#     ProjectYear=np.zeros(uMP.size)
#     Region=np.zeros(uMP.size)
#     for iMP in range(uMP.size):
#         ProjectType[iMP]=meta[pNam]['Project']['ByMP']['Project Type'][uMP[iMP]]
#         ProjectArea[iMP]=meta[pNam]['Project']['ByMP']['Project Area'][uMP[iMP]]
#         ProjectYear[iMP]=meta[pNam]['Project']['ByMP']['Project Year'][uMP[iMP]]
#         Region[iMP]=meta[pNam]['Project']['ByMP']['Region ID'][uMP[iMP]]

#     # Time series of saved results
#     tv_saving=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)

#     # Years of implementation
#     uT=np.unique(ProjectYear[ProjectYear!=0])
#     tv_Implementation=np.arange(np.min(uT),np.max(uT)+1,1)

#     # Variables to include
#     v2include=['Yield Lumber','Yield Plywood','Yield OSB','Yield MDF','Yield Paper','Yield Pellets','Yield PowerGrid',
#                'Yield PowerFacilityDom','Yield FirewoodDom','Yield LogExport','Cost Roads','Cost Harvest Overhead',
#                'Cost Harvest Felling and Piling','Cost Harvest Hauling','Cost Harvest Residuals',
#                'Cost Milling','Cost Nutrient Management','Cost Planting','Cost Survey','Cost Knockdown','Cost Ripping','Cost PAS Deactivation',
#                'Cost Slashpile Burn','Cost Total','Cost Silviculture Total','Revenue Lumber','Revenue Plywood',
#                'Revenue OSB','Revenue MDF','Revenue Paper','Revenue PowerFacilityDom','Revenue PowerGrid','Revenue Pellets','Revenue FirewoodDom',
#                'Revenue LogExport','Revenue Gross','Revenue Net','Revenue Net Disc','Revenue Gross Disc','Cost Total Disc','Cost Nutrient Management Disc',
#                'Cost Silviculture Total Disc','Cost Total_cumu','Cost Silviculture Total_cumu','Cost Nutrient Management_cumu','Cost Total Disc_cumu',
#                'Cost Silviculture Total Disc_cumu','Cost Nutrient Management Disc_cumu','Revenue Gross_cumu','Revenue Gross Disc_cumu','Revenue Net_cumu',
#                'Revenue Net Disc_cumu']

#     # Scale factor used to temporarily store data
#     ScaleFactor=1.0

#     # Loop through scenarios
#     for iScn in range(meta[pNam]['Project']['N Scenario']):

#         #--------------------------------------------------------------------------
#         # Initialize data structures
#         #--------------------------------------------------------------------------

#         MoMeanByPT={}
#         MoSumByPT={}
#         MoMeanByPTAndYr={}
#         MoSumByPTAndYr={}
#         for k in v2include:
#             MoMeanByPT[k]=np.zeros( (tv_saving.size,meta[pNam]['Project']['N Ensemble'],uPT.size,uReg.size) ,dtype='float64')
#             MoSumByPT[k]=np.zeros( (tv_saving.size,meta[pNam]['Project']['N Ensemble'],uPT.size,uReg.size) ,dtype='float64')

#             MoMeanByPTAndYr[k]={}
#             MoSumByPTAndYr[k]={}
#             for t in tv_Implementation:
#                 MoMeanByPTAndYr[k][t]=np.zeros( (tv_saving.size,meta[pNam]['Project']['N Ensemble'],uPT.size,uReg.size) ,dtype='float64')
#                 MoSumByPTAndYr[k][t]=np.zeros( (tv_saving.size,meta[pNam]['Project']['N Ensemble'],uPT.size,uReg.size) ,dtype='float64')

#         # Loop through ensembles
#         for iEns in range(meta[pNam]['Project']['N Ensemble']):

#             #--------------------------------------------------------------------------
#             # Initialize temporary data structure for full simulation
#             #--------------------------------------------------------------------------

#             Data0={}
#             for k in v2include:
#                 Data0[k]=np.zeros( (tv_saving.size,meta[pNam]['Project']['N Stand']) ,dtype='float64')

#             #--------------------------------------------------------------------------
#             # Import batches
#             #--------------------------------------------------------------------------

#             for iBat in range(meta[pNam]['Project']['N Batch']):

#                 indBat=IndexToBatch(meta[pNam],iBat)

#                 d1=LoadSingleOutputFile(meta,pNam,iScn,iEns,iBat)

#                 if (meta[pNam]['Scenario'][iScn]['Harvest Status Future']=='On') | (meta[pNam]['Scenario'][iScn]['Breakup Status Historical']=='On') | (meta[pNam]['Scenario'][iScn]['Breakup Status Future']=='On'):
#                     ec=gu.ipickle(meta['Paths'][pNam]['Input Scenario'][iScn] + '\\Modified_Events_Ens' + FixFileNum(iEns) + '_Bat' + FixFileNum(iBat) + '.pkl')
#                 else:
#                     ec=gu.ipickle(meta['Paths'][pNam]['Input Scenario'][iScn] + '\\Events_Ens' + FixFileNum(iEns) + '_Bat' + FixFileNum(iBat) + '.pkl')

#                 # Uncompress event chronology if it has been compressed
#                 ec=EventChronologyDecompress(meta,pNam,ec,iScn,iEns,iBat)

#                 # Inventory
#                 inv=gu.ipickle(meta['Paths'][pNam]['Input Scenario'][iScn] + '\\Inventory_Bat' + FixFileNum(iBat) + '.pkl')

#                 # Cashflow
#                 econ1=econo.CalculateNetRevenue(meta,iScn,iEns,iBat,inv,ec,d1)

#                 for k in v2include:

#                     Data0[k][:,indBat]=econ1[k].copy()

#             del d1,ec,inv,econ1
#             #garc.collect()

#             #--------------------------------------------------------------------------
#             # Convert from SXY results to MP averages
#             #--------------------------------------------------------------------------

#             MeanByMP={}
#             SumByMP={}
#             for k in v2include:
#                 MeanByMP[k]=np.zeros( (tv_saving.size,uMP.size) )
#                 SumByMP[k]=np.zeros( (tv_saving.size,uMP.size) )

#             for iMP in range(uMP.size):

#                 # Index to SXY for the ith multipolygon
#                 indSXY=Crosswalk_sxy_to_mp[iMP]['Index']

#                 for k in v2include:
#                     MeanByMP[k][:,iMP]=MeanByMP[k][:,iMP]+np.mean(Data0[k][:,indSXY],axis=1)
#                     SumByMP[k][:,iMP]=SumByMP[k][:,iMP]+np.sum(ProjectArea[iMP])*np.mean(Data0[k][:,indSXY],axis=1)

#             #--------------------------------------------------------------------------
#             # Summarize by project type, region, and time
#             #--------------------------------------------------------------------------

#             for iPT in range(uPT.size):

#                 for iReg in range(uReg.size):

#                     ind1=np.where( (ProjectType==uPT[iPT]) & (Region==uReg[iReg]) )[0]

#                     for k in v2include:

#                         MoMeanByPT[k][:,iEns,iPT,iReg]=ScaleFactor*np.mean(MeanByMP[k][:,ind1],axis=1)
#                         MoSumByPT[k][:,iEns,iPT,iReg]=ScaleFactor*np.sum(SumByMP[k][:,ind1],axis=1)

#                         # Summmarize by PT and year
#                         for t in tv_Implementation:

#                             ind2=np.where( (ProjectType==uPT[iPT]) & (Region==uReg[iReg]) & (ProjectYear==t) )[0]

#                             MoMeanByPTAndYr[k][t][:,iEns,iPT,iReg]=ScaleFactor*np.mean(MeanByMP[k][:,ind2],axis=1)
#                             MoSumByPTAndYr[k][t][:,iEns,iPT,iReg]=ScaleFactor*np.sum(SumByMP[k][:,ind2],axis=1)

#         #--------------------------------------------------------------------------
#         # Save
#         #--------------------------------------------------------------------------

#         gu.opickle(meta['Paths'][pNam]['Data'] + '\\Outputs\\Econ_ByProjTypeReg_Mean_Scn' + str(iScn+1) + '_FromMPs.pkl',MoMeanByPT)
#         gu.opickle(meta['Paths'][pNam]['Data'] + '\\Outputs\\Econ_ByProjTypeReg_Sum_Scn' + str(iScn+1) + '_FromMPs.pkl',MoSumByPT)

#         gu.opickle(meta['Paths'][pNam]['Data'] + '\\Outputs\\Econ_ByProjTypeRegAndYear_Mean_Scn' + str(iScn+1) + '_FromMPs.pkl',MoMeanByPTAndYr)
#         gu.opickle(meta['Paths'][pNam]['Data'] + '\\Outputs\\Econ_ByProjTypeRegAndYear_Sum_Scn' + str(iScn+1) + '_FromMPs.pkl',MoSumByPTAndYr)

#     t1=time.time()
#     print((t1-t0)/60)

#     return

# #%% Import scenario data from MP data

# def Import_Scenario_Data_FromMPs(meta,mos):

#     Data=[None]*meta[pNam]['Project']['N Scenario']

#     for iS in range(meta[pNam]['Project']['N Scenario']):

#         Data[iS]={}

#         for oper in ['Mean','Sum']:

#             Data[iS][oper]={}

#             #------------------------------------------------------------------
#             # GHG
#             #------------------------------------------------------------------

#             d=gu.ipickle(meta['Paths'][pNam]['Data'] + '\\Outputs\\GHGB_ByProjTypeReg_' + oper + '_Scn' + str(iS+1) + '_FromMPs.pkl')

#             for k in d.keys():

#                 Data[iS][oper][k]={}
#                 Data[iS][oper][k]['Ensemble Mean']=np.mean(d[k],axis=1)
#                 Data[iS][oper][k]['Ensemble SD']=np.std(d[k],axis=1)
#                 Data[iS][oper][k]['Ensemble SE']=np.std(d[k],axis=1)/np.sqrt(meta[pNam]['Project']['N Ensemble'])
#                 Data[iS][oper][k]['Ensemble P005']=np.percentile(d[k],0.5,axis=1)
#                 Data[iS][oper][k]['Ensemble P025']=np.percentile(d[k],2.5,axis=1)
#                 Data[iS][oper][k]['Ensemble P250']=np.percentile(d[k],25,axis=1)
#                 Data[iS][oper][k]['Ensemble P750']=np.percentile(d[k],75,axis=1)
#                 Data[iS][oper][k]['Ensemble P975']=np.percentile(d[k],97.5,axis=1)
#                 Data[iS][oper][k]['Ensemble P995']=np.percentile(d[k],99.5,axis=1)

#             #------------------------------------------------------------------
#             # Economics
#             #------------------------------------------------------------------

#             d=gu.ipickle(meta['Paths'][pNam]['Data'] + '\\Outputs\\Econ_ByProjTypeReg_' + oper + '_Scn' + str(iS+1) + '_FromMPs.pkl')

#             for k in d.keys():

#                 Data[iS][oper][k]={}
#                 Data[iS][oper][k]['Ensemble Mean']=np.mean(d[k],axis=1)
#                 Data[iS][oper][k]['Ensemble SD']=np.std(d[k],axis=1)
#                 Data[iS][oper][k]['Ensemble SE']=np.std(d[k],axis=1)/np.sqrt(meta[pNam]['Project']['N Ensemble'])
#                 Data[iS][oper][k]['Ensemble P005']=np.percentile(d[k],0.5,axis=1)
#                 Data[iS][oper][k]['Ensemble P025']=np.percentile(d[k],2.5,axis=1)
#                 Data[iS][oper][k]['Ensemble P250']=np.percentile(d[k],25,axis=1)
#                 Data[iS][oper][k]['Ensemble P750']=np.percentile(d[k],75,axis=1)
#                 Data[iS][oper][k]['Ensemble P975']=np.percentile(d[k],97.5,axis=1)
#                 Data[iS][oper][k]['Ensemble P995']=np.percentile(d[k],99.5,axis=1)

#     # Add to MOS structure
#     mos['Scenarios']=Data

#     return mos

# #%% Import scenario comparisons from MP data

# def Import_ScenarioComparisons_FromMPs(meta,mos):

#     for sc in mos['Delta']:

#         mos['Delta'][sc]['ByStrata']={}
#         mos['Delta'][sc]['ByPTAndYear']={}

#         for oper in ['Mean','Sum']:

#             #------------------------------------------------------------------
#             # By projet type and region
#             #------------------------------------------------------------------

#             dC={}

#             # GHG

#             dB=gu.ipickle(meta['Paths'][pNam]['Data'] + '\\Outputs\\GHGB_ByProjTypeReg_' + oper + '_Scn' + str(mos['Delta'][sc]['iB']+1) + '_FromMPs.pkl')
#             dP=gu.ipickle(meta['Paths'][pNam]['Data'] + '\\Outputs\\GHGB_ByProjTypeReg_' + oper + '_Scn' + str(mos['Delta'][sc]['iP']+1) + '_FromMPs.pkl')

#             for k in dB.keys():

#                 dC[k]={}
#                 dC[k]['Ensemble Mean']=np.mean(dP[k]-dB[k],axis=1)
#                 dC[k]['Ensemble SD']=np.std(dP[k]-dB[k],axis=1)
#                 dC[k]['Ensemble SE']=np.std(dP[k]-dB[k],axis=1)/np.sqrt(meta[pNam]['Project']['N Ensemble'])
#                 dC[k]['Ensemble P005']=np.percentile(dP[k]-dB[k],0.5,axis=1)
#                 dC[k]['Ensemble P025']=np.percentile(dP[k]-dB[k],2.5,axis=1)
#                 dC[k]['Ensemble P250']=np.percentile(dP[k]-dB[k],25,axis=1)
#                 dC[k]['Ensemble P750']=np.percentile(dP[k]-dB[k],75,axis=1)
#                 dC[k]['Ensemble P975']=np.percentile(dP[k]-dB[k],97.5,axis=1)
#                 dC[k]['Ensemble P995']=np.percentile(dP[k]-dB[k],99.5,axis=1)

#             # Economics

#             dB=gu.ipickle(meta['Paths'][pNam]['Data'] + '\\Outputs\\Econ_ByProjTypeReg_' + oper + '_Scn' + str(mos['Delta'][sc]['iB']+1) + '_FromMPs.pkl')
#             dP=gu.ipickle(meta['Paths'][pNam]['Data'] + '\\Outputs\\Econ_ByProjTypeReg_' + oper + '_Scn' + str(mos['Delta'][sc]['iP']+1) + '_FromMPs.pkl')

#             for k in dB.keys():

#                 dC[k]={}
#                 dC[k]['Ensemble Mean']=np.mean(dP[k]-dB[k],axis=1)
#                 dC[k]['Ensemble SD']=np.std(dP[k]-dB[k],axis=1)
#                 dC[k]['Ensemble SE']=np.std(dP[k]-dB[k],axis=1)/np.sqrt(meta[pNam]['Project']['N Ensemble'])
#                 dC[k]['Ensemble P005']=np.percentile(dP[k]-dB[k],0.5,axis=1)
#                 dC[k]['Ensemble P025']=np.percentile(dP[k]-dB[k],2.5,axis=1)
#                 dC[k]['Ensemble P250']=np.percentile(dP[k]-dB[k],25,axis=1)
#                 dC[k]['Ensemble P750']=np.percentile(dP[k]-dB[k],75,axis=1)
#                 dC[k]['Ensemble P975']=np.percentile(dP[k]-dB[k],97.5,axis=1)
#                 dC[k]['Ensemble P995']=np.percentile(dP[k]-dB[k],99.5,axis=1)

#             mos['Delta'][sc]['ByStrata'][oper]=dC

#             #------------------------------------------------------------------
#             # By project type, region and year
#             #------------------------------------------------------------------

#             dC={}

#             # GHG

#             dB=gu.ipickle(meta['Paths'][pNam]['Data'] + '\\Outputs\\GHGB_ByProjTypeRegAndYear_' + oper + '_Scn' + str(mos['Delta'][sc]['iB']+1) + '_FromMPs.pkl')
#             dP=gu.ipickle(meta['Paths'][pNam]['Data'] + '\\Outputs\\GHGB_ByProjTypeRegAndYear_' + oper + '_Scn' + str(mos['Delta'][sc]['iP']+1) + '_FromMPs.pkl')

#             for k in dB.keys():

#                 dC[k]={}
#                 for t in dB[k].keys():

#                     dC[k][t]={}
#                     dC[k][t]['Ensemble Mean']=np.mean(dP[k][t]-dB[k][t],axis=1)
#                     dC[k][t]['Ensemble SD']=np.std(dP[k][t]-dB[k][t],axis=1)
#                     dC[k][t]['Ensemble SE']=np.std(dP[k][t]-dB[k][t],axis=1)/np.sqrt(meta[pNam]['Project']['N Ensemble'])
#                     dC[k][t]['Ensemble P005']=np.percentile(dP[k][t]-dB[k][t],0.5,axis=1)
#                     dC[k][t]['Ensemble P025']=np.percentile(dP[k][t]-dB[k][t],2.5,axis=1)
#                     dC[k][t]['Ensemble P250']=np.percentile(dP[k][t]-dB[k][t],25,axis=1)
#                     dC[k][t]['Ensemble P750']=np.percentile(dP[k][t]-dB[k][t],75,axis=1)
#                     dC[k][t]['Ensemble P975']=np.percentile(dP[k][t]-dB[k][t],97.5,axis=1)
#                     dC[k][t]['Ensemble P995']=np.percentile(dP[k][t]-dB[k][t],99.5,axis=1)

#             # Economics

#             dB=gu.ipickle(meta['Paths'][pNam]['Data'] + '\\Outputs\\Econ_ByProjTypeRegAndYear_' + oper + '_Scn' + str(mos['Delta'][sc]['iB']+1) + '_FromMPs.pkl')
#             dP=gu.ipickle(meta['Paths'][pNam]['Data'] + '\\Outputs\\Econ_ByProjTypeRegAndYear_' + oper + '_Scn' + str(mos['Delta'][sc]['iP']+1) + '_FromMPs.pkl')

#             for k in dB.keys():
#                 dC[k]={}
#                 for t in dB[k].keys():

#                     #dB[k][t]=dB[k][t].astype(float)/ScaleFactor
#                     #dP[k][t]=dP[k][t].astype(float)/ScaleFactor

#                     dC[k][t]={}
#                     dC[k][t]['Ensemble Mean']=np.mean(dP[k][t]-dB[k][t],axis=1)
#                     dC[k][t]['Ensemble SD']=np.std(dP[k][t]-dB[k][t],axis=1)
#                     dC[k][t]['Ensemble SE']=np.std(dP[k][t]-dB[k][t],axis=1)/np.sqrt(meta[pNam]['Project']['N Ensemble'])
#                     dC[k][t]['Ensemble P005']=np.percentile(dP[k][t]-dB[k][t],0.5,axis=1)
#                     dC[k][t]['Ensemble P025']=np.percentile(dP[k][t]-dB[k][t],2.5,axis=1)
#                     dC[k][t]['Ensemble P250']=np.percentile(dP[k][t]-dB[k][t],25,axis=1)
#                     dC[k][t]['Ensemble P750']=np.percentile(dP[k][t]-dB[k][t],75,axis=1)
#                     dC[k][t]['Ensemble P975']=np.percentile(dP[k][t]-dB[k][t],97.5,axis=1)
#                     dC[k][t]['Ensemble P995']=np.percentile(dP[k][t]-dB[k][t],99.5,axis=1)

#             mos['Delta'][sc]['ByPTAndYear'][oper]=dC

#     return mos

#%% Save MOS GHG output variables by multipolygon subset

#def MosByMPSubset_GHGB(meta,ListSubsetMP):
#
#    t0=time.time()
#
#    # Error multiplier
#    sigma_multiplier=1.0
#
#    # Import multipolygons
#    atu_multipolygons=gu.ipickle(meta['Paths']['Geospatial'] + '\\atu_multipolygons.pkl')
#
#    # Import sxy
#    sxy=gu.ipickle(meta['Paths']['Geospatial'] + '\\sxy.pkl')
#
#    # Create listed index (faster than indexing on the fly)
#    Crosswalk_sxy_to_mp=[None]*len(ListSubsetMP)
#    for iMP in range(len(ListSubsetMP)):
#        d={}
#        d['Index']=np.where(sxy['ID_atu_multipolygons']==ListSubsetMP[iMP])[0]
#        Crosswalk_sxy_to_mp[iMP]=d
#
#    # Area
#    Area=np.zeros( len(ListSubsetMP) )
#    for iMP in range(len(ListSubsetMP)):
#        Area0=atu_multipolygons[ListSubsetMP[iMP]]['ACTUAL_TREATMENT_AREA']
#        if Area0==None:
#            print('Encounterd no area, using zero')
#            Area0=0
#        Area[iMP]=Area0
#
#    # Time series of saved results
#    tv_saving=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
#
#    # Op
#    oper=['Mean','Sum']
#
#    # GHG balance variables
#    d0=LoadSingleOutputFile(meta,0,0,0)
#
#    # All (excluding 'C_M_ByAgent' and 'Year')
#    v2include=['A', 'V_Merch', 'V_MerchToMill', 'LogSizeEnhancement', 'C_Biomass_Tot', 'C_Piled_Tot', 'C_Litter_Tot', 'C_DeadWood_Tot', 'C_Soil_Tot',
#               'C_InUse_Tot', 'C_DumpLandfill_Tot', 'C_M_Dist', 'C_G_Gross_Tot', 'C_G_Net_Tot', 'C_M_Reg_Tot', 'C_LF_Tot', 'C_RH_Tot', 'C_ToMillMerch',
#               'C_ToMillNonMerch', 'C_ToMillSnagStem', 'C_ToSlashpileBurnTot', 'C_ToLumber', 'C_ToPlywood', 'C_ToOSB', 'C_ToMDF', 'C_ToPaper', 'C_ToPowerFacilityDom',
#               'C_ToPowerFacilityExport', 'C_ToPowerGrid', 'C_ToPellets', 'C_ToFirewoodDom', 'C_ToFirewoodExport', 'C_ToLogExport', 'E_CO2e_LULUCF_NEE', 'E_CO2e_LULUCF_Wildfire',
#               'E_CO2e_LULUCF_OpenBurning', 'E_CO2e_LULUCF_EcoOther', 'E_CO2e_LULUCF_HWP', 'E_CO2e_ESC_Comb', 'E_CO2e_ESC_SubE', 'E_CO2e_ESC_SubBM', 'E_CO2e_ET_Comb',
#               'E_CO2e_IPPU_Comb', 'C_NPP_Tot', 'C_ToMill', 'E_CO2e_LULUCF_Fire', 'E_CO2e_AGHGB_WSub', 'E_CO2e_AGHGB_WOSub', 'E_CO2e_AGHGB_WSub_cumu', 'E_CO2e_AGHGB_WOSub_cumu',
#               'E_CO2e_AGHGB_WSub_cumu_from_tref', 'E_CO2e_AGHGB_WOSub_cumu_from_tref']
#
#    #--------------------------------------------------------------------------
#    # Initialize data by multipolygon structure
#    #--------------------------------------------------------------------------
#
#    mos=[None]*meta[pNam]['Project']['N Scenario']
#
#    for iScn in range(meta[pNam]['Project']['N Scenario']):
#
#        d={}
#        for op in oper:
#            d[op]={}
#
#        for k in d0.keys():
#
#            if np.isin(k,v2include)==False:
#                continue
#
#            for op in oper:
#                d[op][k]={}
#                d[op][k]['Ensemble Mean']=np.zeros((tv_saving.size,len(ListSubsetMP)))
#                d[op][k]['Ensemble CIL']=np.zeros((tv_saving.size,len(ListSubsetMP)))
#                d[op][k]['Ensemble CIH']=np.zeros((tv_saving.size,len(ListSubsetMP)))
#
#        mos[iScn]=d
#
#    #--------------------------------------------------------------------------
#    # Loop through scenarios
#    #--------------------------------------------------------------------------
#
#    for iScn in range(meta[pNam]['Project']['N Scenario']):
#
#        # Initialize data matrix
#        Data={}
#        for op in oper:
#            Data[op]=np.zeros( (tv_saving.size,meta[pNam]['Project']['N Ensemble'],len(ListSubsetMP),len(v2include)) ,dtype='float64')
#
#        #----------------------------------------------------------------------
#        # Loop through ensembles
#        #----------------------------------------------------------------------
#
#        for iEns in range(meta[pNam]['Project']['N Ensemble']):
#            print(iEns)
#
#            #------------------------------------------------------------------
#            # Import batches
#            #------------------------------------------------------------------
#
#            # Initialize temporary data structure for full simulation
#            DataSXY=np.zeros( (tv_saving.size,meta[pNam]['Project']['N Stand'],len(v2include)) ,dtype='float64')
#
#            # Import batches
#            for iBat in range(meta[pNam]['Project']['N Batch']):
#
#                indBat=IndexToBatch(meta[pNam],iBat)
#
#                d1=LoadSingleOutputFile(meta,pNam,iScn,iEns,iBat)
#
#                cnt_k=0
#                for k in d0.keys():
#
#                    if np.isin(k,v2include)==False:
#                        continue
#
#                    DataSXY[:,indBat,cnt_k]=d1[k].copy()
#
#                    cnt_k=cnt_k+1
#
#            del d1
#            #garc.collect()
#
#            #------------------------------------------------------------------
#            # Calculate mean for multipolygon subset
#            #------------------------------------------------------------------
#
#            for iMP in range(len(ListSubsetMP)):
#
#                indMP=Crosswalk_sxy_to_mp[iMP]['Index']
#
#                mu=np.mean(DataSXY[:,indMP,:],axis=1)
#
#                for op in oper:
#
#                    if op=='Mean':
#                        Data[op][:,iEns,iMP,:]=mu
#                    else:
#                        Data[op][:,iEns,iMP,:]=np.sum(Area[iMP])*mu
#
#        #----------------------------------------------------------------------
#        # Add stats to MOS structure
#        #----------------------------------------------------------------------
#
#        cnt_k=0
#        for k in d0.keys():
#
#            if np.isin(k,v2include)==False:
#                continue
#
#            for op in oper:
#
#                mu=np.mean(Data[op][:,:,:,cnt_k],axis=1)
#                sd=np.std(Data[op][:,:,:,cnt_k],axis=1)
#                cil=mu-sigma_multiplier*sd/np.sqrt(meta[pNam]['Project']['N Ensemble'])
#                cih=mu+sigma_multiplier*sd/np.sqrt(meta[pNam]['Project']['N Ensemble'])
#
#                mos[iScn][op][k]['Ensemble Mean']=mu
#                mos[iScn][op][k]['Ensemble CIL']=cil
#                mos[iScn][op][k]['Ensemble CIH']=cih
#
#            cnt_k=cnt_k+1
#
#    #--------------------------------------------------------------------------
#    # Save
#    #--------------------------------------------------------------------------
#
#    gu.opickle(meta['Paths'][pNam]['Data'] + '\\Outputs\\MosByMPSubset_GHGB.pkl',mos)
#
#    t1=time.time()
#    print((t1-t0)/60)
#
#    return


#%% GET TASS GROWTH CURVES

def GetTASSCurves(meta,iScn,iGC,fny,fnc,fnm):

    # Assume just one stand per scenario (i.e., demo)
    iS=0

    dfY=pd.read_csv(fny,header=30)

    # Initialize age response of net growth
    G=np.zeros((N_Age,1,6),dtype='int16')

    G[:,iS,0]=G_StemMerch[:,indTIPSY[0]]/meta['Modules']['GYM']['Scale Factor']
    G[:,iS,1]=G_StemNonMerch[:,indTIPSY[0]]/meta['Modules']['GYM']['Scale Factor']
    G[:,iS,2]=G_Bark[:,indTIPSY[0]]/meta['Modules']['GYM']['Scale Factor']
    G[:,iS,3]=G_Branch[:,indTIPSY[0]]/meta['Modules']['GYM']['Scale Factor']
    G[:,iS,4]=G_Foliage[:,indTIPSY[0]]/meta['Modules']['GYM']['Scale Factor']
    G[:,iS,5]=G_VStemMerch[:,indTIPSY[0]]/meta['Modules']['GYM']['Scale Factor']

    # Save data to file in input variables folder of project
    gu.opickle(meta['Paths'][pNam]['Data'] + '\\Inputs\\Scenario' + FixFileNum(iScn) + '\\GrowthCurve' + str(iGC+1) + '_Bat' + FixFileNum(1) + '.pkl',G)

    return

#%%

def Import_BatchTIPSY_Output(meta,pNam,iScn,iGC):

    # Import unique growth curves
    ugc=gu.ipickle(meta['Paths'][pNam]['Data'] + '\\Inputs\\ugc.pkl')

    # TIPSY exports curves as MgDM/ha/yr, CBRunner expects inputs of MgC/ha/yr. Create
    # conversion factor.
    dm2c=0.5

    # Growth curve parameters and TIPSY outputs
    #dfPar=pd.read_excel(meta['Paths'][pNam]['Data'] + '\\Inputs\\GrowthCurvesTIPSY_Parameters.xlsx',sheet_name='Sheet1',skiprows=7)
    txtDat=np.loadtxt(meta['Paths'][pNam]['Data'] + '\\Inputs\\GrowthCurvesTIPSY_Output.out',skiprows=4)

    # TIPSY saves to text file -> convert to dataframe (column names must match TIPSY output file design)
    dfDat=pd.DataFrame(txtDat,columns=meta['Modules']['GYM']['BatchTIPSY Column Names'])

    del txtDat

    # Define age vector (must be consistent with how TIPSY was set up)
    Age=np.arange(0,meta['Modules']['GYM']['BatchTIPSY Maximum Age']+1,1)

    # Get dimensions of the TIPSY output file to reshape the data into Age x Stand
    N_Age=Age.size
    N_GC=int(dfDat.shape[0]/N_Age)

    # Define the fraction of merchantable stemwood
    fMerch=np.nan_to_num(np.reshape(dfDat['VolMerch125'].values,(N_Age,N_GC),order='F')/np.reshape(dfDat['VolTot0'].values,(N_Age,N_GC),order='F'))
    fNonMerch=1-fMerch

    # Merchantable stemwood volume
    V_Merch=np.reshape(dfDat['VolMerch125'].values,(N_Age,N_GC),order='F')
    G_VStemMerch=np.append(np.zeros((1,N_GC)),np.diff(V_Merch,axis=0),axis=0)

    # Extract age responses for each biomass pool
    C_Stem=dm2c*np.reshape(dfDat['ODT_Stem'].values,(N_Age,N_GC),order='F')
    C_StemMerch=fMerch*C_Stem
    C_StemNonMerch=fNonMerch*C_Stem
    C_Foliage=dm2c*np.reshape(dfDat['ODT_Foliage'].values,(N_Age,N_GC),order='F')
    C_Branch=dm2c*np.reshape(dfDat['ODT_Branch'].values,(N_Age,N_GC),order='F')
    C_Bark=dm2c*np.reshape(dfDat['ODT_Bark'].values,(N_Age,N_GC),order='F')

    # Calculate growth
    z=np.zeros((1,N_GC))
    G_StemMerch=np.append(z,np.diff(C_StemMerch,axis=0),axis=0)
    G_StemNonMerch=np.append(z,np.diff(C_StemNonMerch,axis=0),axis=0)
    G_Stem=G_StemMerch+G_StemNonMerch
    G_Foliage=np.append(z,np.diff(C_Foliage,axis=0),axis=0)
    G_Branch=np.append(z,np.diff(C_Branch,axis=0),axis=0)
    G_Bark=np.append(z,np.diff(C_Bark,axis=0),axis=0)

    # Fix growth of year zero
    G_Stem[0,:]=G_Stem[1,:]
    G_StemMerch[0,:]=G_StemMerch[1,:]
    G_StemNonMerch[0,:]=G_StemNonMerch[1,:]
    G_Foliage[0,:]=G_Foliage[1,:]
    G_Branch[0,:]=G_Branch[1,:]
    G_Bark[0,:]=G_Bark[1,:]

    #del C_Stem,C_StemMerch,C_StemNonMerch,C_Foliage,C_Branch,C_Bark,fMerch,fNonMerch

    # Index to the full set of growth curves for scenario iScn and growth curve iGC
    ind_ugc_ScnAndGc=np.where( (ugc['Full'][:,1]==iScn) & (ugc['Full'][:,2]==meta['Modules']['GYM']['ID GC Unique'][iGC]) )[0]

    # Extract the unique growth curve ID for scenario iScn and growth curve iGC
    ID_ugc_ScnAndGc=ugc['Full'][ind_ugc_ScnAndGc,0]

    # Extract the inverse index for scenario iScn and growth curve iGC
    Inverse_ugc_ScnAndGc=ugc['Inverse'][ind_ugc_ScnAndGc]

    # Intersect
    ind=np.arange(0,meta[pNam]['Project']['N Stand'],1,dtype=int)
    c,inda,indb=np.intersect1d(ID_ugc_ScnAndGc,ind,return_indices=True)

    Inverse_ugc_ScnAndGcAndBat=Inverse_ugc_ScnAndGc[inda]

    # Initialize array of growth data
    d={}
    d['Csw']=np.zeros((N_Age,ind.size))
    d['Gsw_Net']=np.zeros((N_Age,ind.size))
    for i in range(inda.size):
        iStand=indb[i]
        iGC_Unique=Inverse_ugc_ScnAndGcAndBat[i]
        d['Csw'][:,iStand]=C_Stem[:,iGC_Unique].T
        d['Gsw_Net'][:,iStand]=G_Stem[:,iGC_Unique].T

    return d

#%% POST-PROCESS TIPSY GROWTH CURVES
# Nested list, gc[Scenario][Stand][Growth Curve]
# *** This is problematic - I think it only works when total GCs = unique GCs ***

def Import_BatchTIPSY_Output_OLD(meta):

    # Growth curve parameters and TIPSY outputs
    dfPar=pd.read_excel(meta['Paths'][pNam]['Data'] + '\\Inputs\\GrowthCurvesTIPSY_Parameters.xlsx',sheet_name='Sheet1',skiprows=6)
    txtDat=np.loadtxt(meta['Paths'][pNam]['Data'] + '\\Inputs\\GrowthCurvesTIPSY_Output.out',skiprows=4)

    # TIPSY saves to text file -> convert to dataframe (column names must match TIPSY output file design)
    dfDat=pd.DataFrame(txtDat,columns=meta['Modules']['GYM']['BatchTIPSY Column Names'])

    # Define age vector (must be consistent with how TIPSY was set up)
    Age=np.arange(0,meta['Modules']['GYM']['BatchTIPSY Maximum Age']+1,1)

    # Get dimensions of the TIPSY output file to reshape the data into Age x Stand
    N_Age=Age.size
    N_GC=int(dfDat.shape[0]/N_Age)

    gc=[None]*N_GC
    for i in range(N_GC):
        gc[i]={}
    for i in range(len(meta['Modules']['GYM']['BatchTIPSY Column Names'])):
        data=np.reshape(dfDat[meta['Modules']['GYM']['BatchTIPSY Column Names'][i]].values,(N_Age,N_GC),order='F')
        for j in range(N_GC):
            gc[j][meta['Modules']['GYM']['BatchTIPSY Column Names'][i]]=data[:,j]

    gc2=[]
    uScn=np.unique(dfPar['ID_Scenario'])
    for iScn in range(uScn.size):
        ind=np.where(dfPar['ID_Scenario']==uScn[iScn])[0]
        uStand=np.unique(dfPar.loc[ind,'ID_Stand'])
        gc1=[]
        for iS in range(uStand.size):
            ind=np.where( (dfPar['ID_Scenario']==uScn[iScn]) & (dfPar['ID_Stand']==uStand[iS]) )[0]
            uGC=np.unique(dfPar.loc[ind,'ID_GC'])
            gc0=[]
            for iGC in range(uGC.size):
                ind=np.where( (dfPar['ID_Scenario']==uScn[iScn]) & (dfPar['ID_Stand']==uStand[iS]) & (dfPar['ID_GC']==uGC[iGC]) )[0]
                d={}
                for i in range(len(meta['Modules']['GYM']['BatchTIPSY Column Names'])):
                    data=np.reshape(dfDat[meta['Modules']['GYM']['BatchTIPSY Column Names'][i]].values,(N_Age,N_GC),order='F')
                    d[meta['Modules']['GYM']['BatchTIPSY Column Names'][i]]=data[:,ind]
                gc0.append(d)
            gc1.append(gc0)
        gc2.append(gc1)
    return gc2

#%% Import growth curves

def Import_CompiledGrowthCurves(meta,scn):

    gc=[]
    for iScn in range(len(scn)): #range(meta[pNam]['Project']['N Scenario']):
        gc0=[]

        gc1=[]
        for iBat in range(0,meta[pNam]['Project']['N Batch']):
            tmp=gu.ipickle(meta['Paths'][pNam]['Data'] + '\\Inputs\\Scenario' + FixFileNum(iScn) + '\\GrowthCurve1_Bat' + FixFileNum(iBat) + '.pkl')
            tmp=tmp[:,:,0].astype(float)
            for iS in range(tmp.shape[1]):
                gc1.append(tmp[:,iS].copy()*meta['Modules']['GYM']['Scale Factor'])
        gc0.append(gc1.copy())

        gc1=[]
        for iBat in range(0,meta[pNam]['Project']['N Batch']):
            tmp=gu.ipickle(meta['Paths'][pNam]['Data'] + '\\Inputs\\Scenario' + FixFileNum(iScn) + '\\GrowthCurve2_Bat' + FixFileNum(iBat) + '.pkl')
            tmp=tmp[:,:,0].astype(float)
            for iS in range(tmp.shape[1]):
                gc1.append(tmp[:,iS].copy()*meta['Modules']['GYM']['Scale Factor'])
        gc0.append(gc1.copy())

        gc1=[]
        for iBat in range(0,meta[pNam]['Project']['N Batch']):
            tmp=gu.ipickle(meta['Paths'][pNam]['Data'] + '\\Inputs\\Scenario' + FixFileNum(iScn) + '\\GrowthCurve3_Bat' + FixFileNum(iBat) + '.pkl')
            tmp=tmp[:,:,0].astype(float)
            for iS in range(tmp.shape[1]):
                gc1.append(tmp[:,iS].copy()*meta['Modules']['GYM']['Scale Factor'])
        gc0.append(gc1.copy())
        gc.append(gc0.copy())
    return gc

#%% WRITE SPREADSHEET OF BatchTIPSY PARAMTERS

def Write_BatchTIPSY_Input_Spreadsheet(meta,pNam,ugc):

    # Create a function that will return the column corresponding to a variable name
    fin=meta['Paths']['Model']['Code'] + '\\Parameters\\Parameters_TemplateBatchTIPSY.xlsx'
    df_frmt=pd.read_excel(fin,sheet_name='Sheet1')
    gy_labels=df_frmt.loc[5,:].values
    def GetColumn(lab):
        ind=np.where(gy_labels==lab)[0]
        return int(ind+1)

    # Open spreadsheet
    #PathGrowthCurveParameters=meta['Paths'][pNam]['Data'] + '\\Inputs\\GrowthCurvesTIPSY_Parameters.xlsx'
    PathGrowthCurveParameters=fin
    xfile=openpyxl.load_workbook(PathGrowthCurveParameters)
    sheet=xfile.get_sheet_by_name('Sheet1')
    N_headers=7

    # # Overwrite existing data entries with empty cells
    # # *** This is really important - failing to wipe it clean first will lead to
    # # weird parameters ***
    # for i in range(int(1.5*ugc['Unique'].shape[0])):
    #     for j in range(len(gy_labels)):
    #         sheet.cell(row=i+1+N_headers,column=j+1).value=''

    # Initialize counter
    cnt=1

    # Loop through unique stand types
    for iUGC in range(ugc['Unique'].shape[0]):

        # It isn't necessary to populate these, as we will use the crosswalk in
        # Python session instead
        sheet.cell(row=cnt+N_headers,column=1).value=cnt
        sheet.cell(row=cnt+N_headers,column=2).value=0
        sheet.cell(row=cnt+N_headers,column=3).value=0
        sheet.cell(row=cnt+N_headers,column=4).value=0

        # Regeneration type (N, C, P)
        vnam='regeneration_method'
        id=ugc['Unique'][iUGC,np.where(ugc['GC_Variable_List']==vnam)[0]]
        cd=lut_n2s(meta['LUT']['TIPSY'][vnam],id)[0]
        sheet.cell(row=cnt+N_headers,column=GetColumn(vnam)).value=cd

        # Species 1
        vnam='s1'
        id=ugc['Unique'][iUGC,np.where(ugc['GC_Variable_List']==vnam)[0]]
        cd=lut_n2s(meta['LUT']['VEG_COMP_LYR_R1_POLY']['SPECIES_CD_1'],id)[0]
        sheet.cell(row=cnt+N_headers,column=GetColumn(vnam)).value=cd

        vnam='p1'
        dat=ugc['Unique'][iUGC,np.where(ugc['GC_Variable_List']==vnam)[0]]
        sheet.cell(row=cnt+N_headers,column=GetColumn(vnam)).value=int(dat)

        vnam='i1'
        dat=ugc['Unique'][iUGC,np.where(ugc['GC_Variable_List']==vnam)[0]]
        sheet.cell(row=cnt+N_headers,column=GetColumn(vnam)).value=int(dat)

        vnam='gain1'
        dat=ugc['Unique'][iUGC,np.where(ugc['GC_Variable_List']==vnam)[0]]
        if dat!=9999:
            sheet.cell(row=cnt+N_headers,column=GetColumn(vnam)).value=int(dat)
            sheet.cell(row=cnt+N_headers,column=GetColumn(vnam)).value=int(10)

        # Species 2
        vnam='s2'
        id=ugc['Unique'][iUGC,np.where(ugc['GC_Variable_List']==vnam)[0]]
        if id!=9999:
            cd=lut_n2s(meta['LUT']['VEG_COMP_LYR_R1_POLY']['SPECIES_CD_1'],id)[0]
            sheet.cell(row=cnt+N_headers,column=GetColumn(vnam)).value=cd

            vnam='p2'
            dat=ugc['Unique'][iUGC,np.where(ugc['GC_Variable_List']==vnam)[0]]
            sheet.cell(row=cnt+N_headers,column=GetColumn(vnam)).value=int(dat)

            vnam='gain2'
            dat=ugc['Unique'][iUGC,np.where(ugc['GC_Variable_List']==vnam)[0]]
            if dat!=9999:
                sheet.cell(row=cnt+N_headers,column=GetColumn(vnam)).value=int(dat)
                sheet.cell(row=cnt+N_headers,column=GetColumn(vnam)).value=int(10)

        # Species 3
        vnam='s3'
        id=ugc['Unique'][iUGC,np.where(ugc['GC_Variable_List']==vnam)[0]]
        if id!=9999:
            cd=lut_n2s(meta['LUT']['VEG_COMP_LYR_R1_POLY']['SPECIES_CD_1'],id)[0]
            sheet.cell(row=cnt+N_headers,column=GetColumn(vnam)).value=cd

            vnam='p3'
            dat=ugc['Unique'][iUGC,np.where(ugc['GC_Variable_List']==vnam)[0]]
            sheet.cell(row=cnt+N_headers,column=GetColumn(vnam)).value=int(dat)

            vnam='gain3'
            dat=ugc['Unique'][iUGC,np.where(ugc['GC_Variable_List']==vnam)[0]]
            if dat!=9999:
                sheet.cell(row=cnt+N_headers,column=GetColumn(vnam)).value=int(dat)
                sheet.cell(row=cnt+N_headers,column=GetColumn(vnam)).value=int(10)

        # Species 4
        vnam='s4'
        id=ugc['Unique'][iUGC,np.where(ugc['GC_Variable_List']==vnam)[0]]
        if id!=9999:
            cd=lut_n2s(meta['LUT']['VEG_COMP_LYR_R1_POLY']['SPECIES_CD_1'],id)[0]
            sheet.cell(row=cnt+N_headers,column=GetColumn(vnam)).value=cd

            vnam='p4'
            dat=ugc['Unique'][iUGC,np.where(ugc['GC_Variable_List']==vnam)[0]]
            sheet.cell(row=cnt+N_headers,column=GetColumn(vnam)).value=int(dat)

            vnam='gain4'
            dat=ugc['Unique'][iUGC,np.where(ugc['GC_Variable_List']==vnam)[0]]
            if dat!=9999:
                sheet.cell(row=cnt+N_headers,column=GetColumn(vnam)).value=int(dat)
                sheet.cell(row=cnt+N_headers,column=GetColumn(vnam)).value=int(10)

        # Species 5
        vnam='s5'
        id=ugc['Unique'][iUGC,np.where(ugc['GC_Variable_List']==vnam)[0]]
        if id!=9999:
            cd=lut_n2s(meta['LUT']['VEG_COMP_LYR_R1_POLY']['SPECIES_CD_1'],id)[0]
            sheet.cell(row=cnt+N_headers,column=GetColumn(vnam)).value=cd

            vnam='p5'
            dat=ugc['Unique'][iUGC,np.where(ugc['GC_Variable_List']==vnam)[0]]
            sheet.cell(row=cnt+N_headers,column=GetColumn(vnam)).value=int(dat)

            vnam='gain5'
            dat=ugc['Unique'][iUGC,np.where(ugc['GC_Variable_List']==vnam)[0]]
            if dat!=9999:
                sheet.cell(row=cnt+N_headers,column=GetColumn(vnam)).value=int(dat)
                sheet.cell(row=cnt+N_headers,column=GetColumn(vnam)).value=int(10)

        # Planting density
        vnam='init_density'
        dat=ugc['Unique'][iUGC,np.where(ugc['GC_Variable_List']==vnam)[0]]
        sheet.cell(row=cnt+N_headers,column=GetColumn(vnam)).value=int(dat)

        # Regeneration delay
        vnam='regen_delay'
        dat=ugc['Unique'][iUGC,np.where(ugc['GC_Variable_List']==vnam)[0]]
        sheet.cell(row=cnt+N_headers,column=GetColumn(vnam)).value=int(dat)

        # OAF1
        vnam='oaf1'
        dat=ugc['Unique'][iUGC,np.where(ugc['GC_Variable_List']==vnam)[0]]
        sheet.cell(row=cnt+N_headers,column=GetColumn(vnam)).value=dat[0]

        # OAF2
        vnam='oaf2'
        dat=ugc['Unique'][iUGC,np.where(ugc['GC_Variable_List']==vnam)[0]]
        sheet.cell(row=cnt+N_headers,column=GetColumn(vnam)).value=dat[0]

        # BEC zone
        vnam='bec_zone'
        id=ugc['Unique'][iUGC,np.where(ugc['GC_Variable_List']==vnam)[0]]
        cd=lut_n2s(meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE'],id)[0]
        sheet.cell(row=cnt+N_headers,column=GetColumn(vnam)).value=cd

        # FIZ
        vnam='FIZ'
        id=ugc['Unique'][iUGC,np.where(ugc['GC_Variable_List']==vnam)[0]]
        cd=lut_n2s(meta['LUT']['TIPSY']['FIZ'],id)[0]
        sheet.cell(row=cnt+N_headers,column=GetColumn(vnam)).value=cd

        # Age at Aerial BTK Spray
        # vnam='fert_age1'
        #dat=ugc['Unique'][iUGC,np.where(ugc['GC_Variable_List']==vnam)[0]]
        #if dat!=9999:
        #    sheet.cell(row=cnt+N_headers,column=GetColumn(vnam)).value=int(dat[0])

        # Update counter
        cnt=cnt+1

    #------------------------------------------------------------------------------
    # Save to spreadsheet
    #------------------------------------------------------------------------------

    xfile.save(meta['Paths'][pNam]['Data'] + '\\Inputs\\GrowthCurvesTIPSY_Parameters.xlsx')

#%% Write BatchTIPSY input file
# Notes:
#    If the input spreadsheet has nan's, all data will be converted to float

def Write_BatchTIPSY_Input_File(meta,pNam):

    # Import format info (number of designated positions per variable)
    fin=meta['Paths']['Model']['Code'] + '\\Parameters\\Parameters_TemplateBatchTIPSY.xlsx'
    df_frmt=pd.read_excel(fin,sheet_name='Sheet1')

    # Format array
    nfrmt=df_frmt.iloc[1,4:53]

    # Import input data
    df=pd.read_excel(meta['Paths'][pNam]['Data'] + '\\Inputs\\GrowthCurvesTIPSY_Parameters.xlsx',sheet_name='Sheet1',skiprows=6)

    # Convert to dictionary
    Data=df.to_dict('list')
    for k in Data.keys():
        Data[k]=np.array(Data[k])

    varnams=list(df.columns[4:53])

    # Change to list
    DataList=[None]*Data['ID'].size
    for i in range(len(DataList)):
        d={}
        for k in Data.keys():
            d[k]=Data[k][i]
        DataList[i]=d

    # Create a spacer
    Spacer=" "

    # Open file
    fid=open(meta['Paths'][pNam]['Data'] + '\\Inputs\\GrowthCurvesTIPSY_InputVariables.dat','w')

    # Loop through rows of table
    for iStand in range(len(DataList)):

        LineString=""

        for iV in range(len(varnams)):

            s0=str(DataList[iStand][varnams[iV]])

            # Make sure NANs are blank for TIPSY
            if s0=='nan': s0=''
            if s0=='NA': s0=''
            if s0=='NaN': s0=''
            if s0=='': s0=''

            n=nfrmt[iV]
            if len(s0)==n:
                s1=s0
            else:
                s1=s0
                d=n-len(s0)
                for k in range(d):
                    s1=s1 + Spacer

            # Add individual column variables
            LineString=LineString + s1 + Spacer

        # Add end of line
        LineString=LineString + '\n'

        # Save
        fid.write(LineString)

    # Close
    fid.close()

#%% IMPORT DISTURBANCE HISTORY

def GetDisturbanceHistory(meta):
    dh=[]
    for iScn in range(meta[pNam]['Project']['N Scenario']):
        dhB=[]
        for iBat in range(meta[pNam]['Project']['N Batch']):
            iEns=0
            dh0=gu.ipickle(meta['Paths'][pNam]['Input Scenario'][iScn] + '\\Events_Ens' + FixFileNum(iEns) + '_Bat' + FixFileNum(iBat) + '.pkl')
            if iBat==0:
                dhB=dh0
            else:
                for i in range(len(dh0)):
                    dhB.append(dh0[i])
            dh.append(dhB)
    return dh

#%% GRAPHICS

# Look at variables in rcParams:
#plt.rcParams.keys()

def Import_GraphicsParameters(x):

    params={}

    if x=='FCI_Demo':

        fs1=6
        fs2=7

        params={'font.sans-serif':'Arial',
                'font.size':fs1,
                'figure.titlesize':fs2,
                'figure.dpi':150,
                'figure.constrained_layout.use':True,
                'axes.edgecolor':'black',
                'axes.labelsize':fs1,
                'axes.labelcolor':'black',
                'axes.titlesize':fs2,
                'axes.titlepad':2,
                'axes.linewidth':0.5,
                'lines.linewidth':1,
                'text.color':'black',
                'xtick.color':'black',
                'xtick.labelsize':fs1,
                'xtick.major.width':0.5,
                'xtick.major.size':3,
                'xtick.direction':'in',
                'ytick.color':'black',
                'ytick.labelsize':fs1,
                'ytick.major.width':0.5,
                'ytick.major.size':3,
                'ytick.direction':'in',
                'legend.fontsize':fs1,
                'savefig.dpi':900,
                'savefig.transparent':True,
                'savefig.format':'png',
                'savefig.pad_inches':0.1,
                'savefig.bbox':'tight'}

    elif x=='bc1ha_1':

        params={'font.sans-serif':'Arial',
                'font.size':7,
                'axes.labelsize':7,
                'axes.titlesize':14,
                'axes.linewidth':0.5,
                'xtick.labelsize':7,
                'xtick.major.width':0.5,
                'xtick.major.size':5,
                'xtick.direction':'in',
                'ytick.labelsize':7,
                'ytick.major.width':0.5,
                'ytick.major.size':5,
                'ytick.direction':'in',
                'legend.fontsize':10,
                'savefig.dpi':150}

    elif (x=='Article') | (x=='article'):

        fs=6

        params={'font.sans-serif':'Arial','font.size':fs,'axes.edgecolor':'black','axes.labelsize':fs,'axes.labelcolor':'black',
                'axes.titlesize':fs,'axes.linewidth':0.5,'lines.linewidth':0.5,'text.color':'black','xtick.color':'black',
                'xtick.labelsize':fs,'xtick.major.width':0.5,'xtick.major.size':3,'xtick.direction':'in','ytick.color':'black',
                'ytick.labelsize':fs,'ytick.major.width':0.5,'ytick.major.size':3,'ytick.direction':'in','legend.fontsize':fs}

    else:

        params={}

    return params

#%% Prepare inventory from spreadsheet

def PrepareInventoryFromSpreadsheet(meta,pNam):

    for iScn in range(0,meta[pNam]['Project']['N Scenario']):

        # Loop through batches, saving inventory to file
        for iBat in range(0,meta[pNam]['Project']['N Batch']):

            inv={}

            # Index to batch
            indBat=IndexToBatch(meta[pNam],iBat)
            N_StandsInBatch=len(indBat)

            # Initialize inventory variables
            inv['Lat']=np.zeros((1,N_StandsInBatch))
            inv['Lon']=np.zeros((1,N_StandsInBatch))
            inv['X']=inv['Lat']
            inv['Y']=inv['Lon']

            # BEC zone
            inv['ID_BGCZ']=np.zeros((1,N_StandsInBatch),dtype='int16')
            inv['ID_BGCZ'][0,:]=meta['LUT']['VEG_COMP_LYR_R1_POLY']['BEC_ZONE_CODE'][meta[pNam]['Scenario'][iScn]['BGC Zone Code']]

            # Timber harvesting landbase (1=yes, 0=no)
            inv['THLB']=meta[pNam]['Scenario'][iScn]['THLB Status']*np.ones((meta[pNam]['Year'].size,N_StandsInBatch))

            # Temperature will be updated automatically
            inv['MAT']=4.0*np.ones((1,N_StandsInBatch))
            #cd=meta[pNam]['Scenario'][iScn]['BGC Zone Code']
            #ind=np.where(meta['Param']['BE']['BGC Zone Averages']['Name']==cd)[0]
            #mat=meta['Param']['BE']['BGC Zone Averages']['MAT'][ind]
            #inv['MAT']=mat*np.ones((1,N_StandsInBatch))

            inv['Wood Density']=meta['Param']['BE']['Biophysical']['Density Wood']*np.ones((1,N_StandsInBatch))

            if meta[pNam]['Project']['Biomass Module']=='Sawtooth':

                ind=np.where( meta['Param']['BE']['Sawtooth']['SRS Key']['SRS_CD']==meta[pNam]['Scenario'][iScn]['SRS1_CD'])[0]

                id=meta['Param']['BE']['Sawtooth']['SRS Key']['SRS_ID'][ind]

                inv['SRS1_ID']=id*np.ones((1,N_StandsInBatch),dtype='int16')
                inv['SRS1_PCT']=100*np.ones((1,N_StandsInBatch),dtype='int16')

                inv['SRS2_ID']=np.ones((1,N_StandsInBatch),dtype='int16')
                inv['SRS2_PCT']=0*np.ones((1,N_StandsInBatch),dtype='int16')

                inv['SRS3_ID']=np.ones((1,N_StandsInBatch),dtype='int16')
                inv['SRS3_PCT']=0*np.ones((1,N_StandsInBatch),dtype='int16')

            # Save
            gu.opickle(meta['Paths'][pNam]['Input Scenario'][iScn] + '\\Inventory_Bat' + FixFileNum(iBat) + '.pkl',inv)

    return

#%% Compile events

def CompileEvents(ec,tv,iS,ID_Type,Year,MortalityFactor,GrowthFactor,ID_GrowthCurve):

    if Year.size>0:
        YearFloor=np.floor(Year)
        uYearFloor=np.unique(YearFloor)
        for iU in range(uYearFloor.size):
            iT=np.where(tv==uYearFloor[iU])[0]
            if iT.size==0:
                continue
            indYear=np.where(YearFloor==uYearFloor[iU])[0]
            for iY in range(indYear.size):
                indAvailable=np.where(ec['ID Event Type'][iT,iS,:].flatten()==0)[0]
                if indAvailable.size==0:
                    print('Warning, more events per year than can be handled!')
                else:
                    iE=indAvailable[0]
                    ec['ID Event Type'][iT,iS,iE]=ID_Type[indYear[iY]]
                    ec['Mortality Factor'][iT,iS,iE]=MortalityFactor[indYear[iY]]
                    ec['Growth Factor'][iT,iS,iE]=GrowthFactor[indYear[iY]]
                    ec['ID Growth Curve'][iT,iS,iE]=ID_GrowthCurve[indYear[iY]]
    return ec

#%% Mortality frequency distribution

def GetMortalityFrequencyDistribution(meta,pNam):

    iEns=0

    M=[None]*meta[pNam]['Project']['N Scenario']
    for iScn in range(meta[pNam]['Project']['N Scenario']):

        tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
        M[iScn]={}
        M[iScn]['Ma']={}
        M[iScn]['Mr']={}
        for k in meta['LUT']['Event'].keys():
            M[iScn]['Ma'][k]=np.zeros((tv.size,meta[pNam]['Project']['N Stand']))
            M[iScn]['Mr'][k]=np.zeros((tv.size,meta[pNam]['Project']['N Stand']))
        M[iScn]['Ma']['Reg']=np.zeros((tv.size,meta[pNam]['Project']['N Stand']))
        M[iScn]['Mr']['Reg']=np.zeros((tv.size,meta[pNam]['Project']['N Stand']))

        for iBat in range(meta[pNam]['Project']['N Batch']):

            indBat=IndexToBatch(meta[pNam],iBat)

            d1=LoadSingleOutputFile(meta,pNam,iScn,iEns,iBat)
            dh=gu.ipickle(meta['Paths'][pNam]['Input Scenario'][iScn] + '\\Events_Ens' + FixFileNum(iEns) + '_Bat' + FixFileNum(iBat) + '.pkl')

            for iStandInBat in range(len(indBat)):
                iStand=indBat[iStandInBat]
                for iYr in range(dh[iStandInBat]['Year'].size):
                    it=np.where(tv==int(dh[iStandInBat]['Year'][iYr]))[0]
                    if it.size==0:
                        continue
                    nam=lut_n2s(meta['LUT']['Event'],dh[iStandInBat]['ID Event Type'][iYr])[0]
                    M[iScn]['Ma'][nam][it,iStand]=d1[0]['C_M_Dist'][it,iStandInBat]
                    M[iScn]['Mr'][nam][it,iStand]=d1[0]['C_M_Dist'][it,iStandInBat]/np.maximum(0.0001,d1[0]['Eco_Biomass'][it-1,iStandInBat])
                    M[iScn]['Ma']['Reg'][it,iStand]=d1[0]['C_M_Reg'][it,iStandInBat]
                    M[iScn]['Mr']['Reg'][it,iStand]=d1[0]['C_M_Reg'][it,iStandInBat]/np.maximum(0.0001,d1[0]['Eco_Biomass'][it-1,iStandInBat])
            del d1,dh
            garc.collect()

    return M

#%% Summarize affected area due to natural disturbance and management

def SummarizeAreaAffected(meta,pNam,mos,tv,iScn,iPS,iSS,iYS,ivlT):

    A={}
    A['Nat Dist']=[None]*5; c=-1
    c=c+1; A['Nat Dist'][c]={}; A['Nat Dist'][c]['Name']='Wildfire'; A['Nat Dist'][c]['Color']=[0.75,0,0]; A['Nat Dist'][c]['Data']=mos[pNam]['Scenarios'][iScn]['Sum']['Area_Wildfire']['Ensemble Mean'][:,iPS,iSS,iYS]*meta[pNam]['Project']['AEF']
    c=c+1; A['Nat Dist'][c]={}; A['Nat Dist'][c]['Name']='Bark beetles'; A['Nat Dist'][c]['Color']=[0.3,0.8,0.2]; A['Nat Dist'][c]['Data']=(mos[pNam]['Scenarios'][iScn]['Sum']['Area_IBM']['Ensemble Mean'][:,iPS,iSS,iYS]+mos[pNam]['Scenarios'][iScn]['Sum']['Area_IBB']['Ensemble Mean'][:,iPS,iSS,iYS]+mos[pNam]['Scenarios'][iScn]['Sum']['Area_IBD']['Ensemble Mean'][:,iPS,iSS,iYS]+mos[pNam]['Scenarios'][iScn]['Sum']['Area_IBS']['Ensemble Mean'][:,iPS,iSS,iYS])*meta[pNam]['Project']['AEF']
    c=c+1; A['Nat Dist'][c]={}; A['Nat Dist'][c]['Name']='Disease'; A['Nat Dist'][c]['Color']=[1,0.5,0.25]; A['Nat Dist'][c]['Data']=(mos[pNam]['Scenarios'][iScn]['Sum']['Area_Disease Root']['Ensemble Mean'][:,iPS,iSS,iYS]+mos[pNam]['Scenarios'][iScn]['Sum']['Area_Disease Foliage']['Ensemble Mean'][:,iPS,iSS,iYS]+mos[pNam]['Scenarios'][iScn]['Sum']['Area_Disease Stem']['Ensemble Mean'][:,iPS,iSS,iYS])*meta[pNam]['Project']['AEF']
    c=c+1; A['Nat Dist'][c]={}; A['Nat Dist'][c]['Name']='Defoliators'; A['Nat Dist'][c]['Color']=[0,0.75,1]; A['Nat Dist'][c]['Data']=mos[pNam]['Scenarios'][iScn]['Sum']['Area_IDW']['Ensemble Mean'][:,iPS,iSS,iYS]*meta[pNam]['Project']['AEF']
    c=c+1; A['Nat Dist'][c]={}; A['Nat Dist'][c]['Name']='Wind, snow & ice'; A['Nat Dist'][c]['Color']=[0,0,0.6]; A['Nat Dist'][c]['Data']=mos[pNam]['Scenarios'][iScn]['Sum']['Area_Wind']['Ensemble Mean'][:,iPS,iSS,iYS]*meta[pNam]['Project']['AEF']
    A['Nat Dist']=A['Nat Dist'][0:c+1]

    A['Management']=[None]*8; c=-1
    c=c+1; A['Management'][c]={}; A['Management'][c]['Name']='Harvest'; A['Management'][c]['Color']=[0,0.75,1]; A['Management'][c]['Data']=(mos[pNam]['Scenarios'][iScn]['Sum']['Area_Harvest']['Ensemble Mean'][:,iPS,iSS,iYS]+mos[pNam]['Scenarios'][iScn]['Sum']['Area_Harvest Salvage']['Ensemble Mean'][:,iPS,iSS,iYS])*meta[pNam]['Project']['AEF']
    c=c+1; A['Management'][c]={}; A['Management'][c]['Name']='Slashpile burn'; A['Management'][c]['Color']=[0.75,0,0]; A['Management'][c]['Data']=(mos[pNam]['Scenarios'][iScn]['Sum']['Area_Slashpile Burn']['Ensemble Mean'][:,iPS,iSS,iYS]+mos[pNam]['Scenarios'][iScn]['Sum']['Area_Harvest Salvage']['Ensemble Mean'][:,iPS,iSS,iYS])*meta[pNam]['Project']['AEF']
    c=c+1; A['Management'][c]={}; A['Management'][c]['Name']='Thinning'; A['Management'][c]['Color']=[0.2,0.4,0.7]; A['Management'][c]['Data']=(mos[pNam]['Scenarios'][iScn]['Sum']['Area_Knockdown']['Ensemble Mean'][:,iPS,iSS,iYS]+mos[pNam]['Scenarios'][iScn]['Sum']['Area_Thinning']['Ensemble Mean'][:,iPS,iSS,iYS])*meta[pNam]['Project']['AEF']
    c=c+1; A['Management'][c]={}; A['Management'][c]['Name']='Site preparation'; A['Management'][c]['Color']=[1,0.7,0.7]; A['Management'][c]['Data']=(mos[pNam]['Scenarios'][iScn]['Sum']['Area_Mechanical Site Prep']['Ensemble Mean'][:,iPS,iSS,iYS])*meta[pNam]['Project']['AEF']
    c=c+1; A['Management'][c]={}; A['Management'][c]['Name']='Prescribed burn'; A['Management'][c]['Color']=[0.5,0,0]; A['Management'][c]['Data']=mos[pNam]['Scenarios'][iScn]['Sum']['Area_Prescribed Burn']['Ensemble Mean'][:,iPS,iSS,iYS]*meta[pNam]['Project']['AEF']
    c=c+1; A['Management'][c]={}; A['Management'][c]['Name']='Planting'; A['Management'][c]['Color']=[0.3,0.8,0.2]; A['Management'][c]['Data']=(mos[pNam]['Scenarios'][iScn]['Sum']['Area_Planting']['Ensemble Mean'][:,iPS,iSS,iYS]+mos[pNam]['Scenarios'][iScn]['Sum']['Area_Direct Seeding']['Ensemble Mean'][:,iPS,iSS,iYS])*meta[pNam]['Project']['AEF']
    c=c+1; A['Management'][c]={}; A['Management'][c]['Name']='Foliage protection'; A['Management'][c]['Color']=[1,0.7,0]; A['Management'][c]['Data']=mos[pNam]['Scenarios'][iScn]['Sum']['Area_Aerial BTK Spray']['Ensemble Mean'][:,iPS,iSS,iYS]*meta[pNam]['Project']['AEF']
    c=c+1; A['Management'][c]={}; A['Management'][c]['Name']='Fertilization'; A['Management'][c]['Color']=[0.75,0.55,1]; A['Management'][c]['Data']=mos[pNam]['Scenarios'][iScn]['Sum']['Area_Fertilization Aerial']['Ensemble Mean'][:,iPS,iSS,iYS]*meta[pNam]['Project']['AEF']
    #c=c+1; A['Management'][c]={}; A['Management'][c]['Name']='Dwarf Mistletoe control'; A['Management'][c]['Color']=[1,0.5,0]; A['Management'][c]['Data']=mos[pNam]['Scenarios'][iScn]['Sum']['Area_Dwarf Mistletoe Control']['Ensemble Mean'][:,iPS,iSS,iYS]*meta[pNam]['Project']['AEF']
    A['Management']=A['Management'][0:c+1]

    # Convert to x-year intervals
    A['tv']=gu.BlockMean(tv,ivlT)
    for i in range(len(A['Nat Dist'])):
        A['Nat Dist'][i]['Data']=gu.BlockMean(A['Nat Dist'][i]['Data'],ivlT)
    for i in range(len(A['Management'])):
        A['Management'][i]['Data']=gu.BlockMean(A['Management'][i]['Data'],ivlT)

    return A

#%% Area affected by individual multipolygon

def AreaAffectedInSingleMultipolygon(meta,iScn,ivlT,tv,MosByMP,iMP):

    A={}
    # Ensemble mean (currently not equipped to give individual ensembles)
    A['Nat Dist']=[None]*10;
    c=-1
    c=c+1; A['Nat Dist'][c]={}; A['Nat Dist'][c]['Name']='Wildfire'; A['Nat Dist'][c]['Color']=[0.75,0,0]; A['Nat Dist'][c]['Data']=MosByMP[iScn]['Area']['Wildfire']['Ensemble Mean'][:,iMP]
    c=c+1; A['Nat Dist'][c]={}; A['Nat Dist'][c]['Name']='Mountain Pine beetle'; A['Nat Dist'][c]['Color']=[0,0.8,0]; A['Nat Dist'][c]['Data']=MosByMP[iScn]['Area']['IBM']['Ensemble Mean'][:,iMP]
    c=c+1; A['Nat Dist'][c]={}; A['Nat Dist'][c]['Name']='Douglas-fir beetle'; A['Nat Dist'][c]['Color']=[0.6,1,0]; A['Nat Dist'][c]['Data']=MosByMP[iScn]['Area']['IBD']['Ensemble Mean'][:,iMP]
    c=c+1; A['Nat Dist'][c]={}; A['Nat Dist'][c]['Name']='Spruce beetle'; A['Nat Dist'][c]['Color']=[0.25,1,1]; A['Nat Dist'][c]['Data']=MosByMP[iScn]['Area']['IBS']['Ensemble Mean'][:,iMP]
    c=c+1; A['Nat Dist'][c]={}; A['Nat Dist'][c]['Name']='W. balsam beetle'; A['Nat Dist'][c]['Color']=[0,0.45,0]; A['Nat Dist'][c]['Data']=MosByMP[iScn]['Area']['IBB']['Ensemble Mean'][:,iMP]
    c=c+1; A['Nat Dist'][c]={}; A['Nat Dist'][c]['Name']='Other pests'; A['Nat Dist'][c]['Color']=[0.8,1,0]; A['Nat Dist'][c]['Data']=MosByMP[iScn]['Area']['Beetles']['Ensemble Mean'][:,iMP]
    c=c+1; A['Nat Dist'][c]={}; A['Nat Dist'][c]['Name']='W. spruce budworm'; A['Nat Dist'][c]['Color']=[0,0.75,1]; A['Nat Dist'][c]['Data']=MosByMP[iScn]['Area']['IDW']['Ensemble Mean'][:,iMP]
    c=c+1; A['Nat Dist'][c]={}; A['Nat Dist'][c]['Name']='Rust'; A['Nat Dist'][c]['Color']=[0.75,0.5,1]; A['Nat Dist'][c]['Data']=MosByMP[iScn]['Area']['Rust Onset']['Ensemble Mean'][:,iMP]
    c=c+1; A['Nat Dist'][c]={}; A['Nat Dist'][c]['Name']='Dwarf Mistletoe'; A['Nat Dist'][c]['Color']=[1,0.5,0.25]; A['Nat Dist'][c]['Data']=MosByMP[iScn]['Area']['Dwarf Mistletoe Onset']['Ensemble Mean'][:,iMP]
    c=c+1; A['Nat Dist'][c]={}; A['Nat Dist'][c]['Name']='Mechanical'; A['Nat Dist'][c]['Color']=[0,0,0.6]; A['Nat Dist'][c]['Data']=MosByMP[iScn]['Area']['Mechanical']['Ensemble Mean'][:,iMP]
    A['Nat Dist']=A['Nat Dist'][0:c+1]

    A['Management']=[None]*10;
    c=-1
    c=c+1; A['Management'][c]={}; A['Management'][c]['Name']='Harvest'; A['Management'][c]['Color']=[0,0.75,1]; A['Management'][c]['Data']=MosByMP[iScn]['Area']['Harvest']['Ensemble Mean'][:,iMP]+MosByMP[iScn]['Area']['Harvest Salvage']['Ensemble Mean'][:,iMP]
    c=c+1; A['Management'][c]={}; A['Management'][c]['Name']='Slashpile burn'; A['Management'][c]['Color']=[0.75,0,0]; A['Management'][c]['Data']=MosByMP[iScn]['Area']['Slashpile Burn']['Ensemble Mean'][:,iMP]+MosByMP[iScn]['Area']['Harvest Salvage']['Ensemble Mean'][:,iMP]
    c=c+1; A['Management'][c]={}; A['Management'][c]['Name']='Thinning and knockdown'; A['Management'][c]['Color']=[0.2,0.4,0.7]; A['Management'][c]['Data']=MosByMP[iScn]['Area']['Knockdown']['Ensemble Mean'][:,iMP]+MosByMP[iScn]['Area']['Thinning']['Ensemble Mean'][:,iMP]
    c=c+1; A['Management'][c]={}; A['Management'][c]['Name']='Site preparation'; A['Management'][c]['Color']=[1,0.7,0.7]; A['Management'][c]['Data']=MosByMP[iScn]['Area']['Disc Trenching']['Ensemble Mean'][:,iMP]+MosByMP[iScn]['Area']['Ripping']['Ensemble Mean'][:,iMP]
    c=c+1; A['Management'][c]={}; A['Management'][c]['Name']='Prescribed burn'; A['Management'][c]['Color']=[0.5,0,0]; A['Management'][c]['Data']=MosByMP[iScn]['Area']['Prescribed Burn']['Ensemble Mean'][:,iMP]
    c=c+1; A['Management'][c]={}; A['Management'][c]['Name']='Dwarf Mistletoe control'; A['Management'][c]['Color']=[1,0.5,0]; A['Management'][c]['Data']=MosByMP[iScn]['Area']['Dwarf Mistletoe Control']['Ensemble Mean'][:,iMP]
    c=c+1; A['Management'][c]={}; A['Management'][c]['Name']='Planting'; A['Management'][c]['Color']=[0.3,0.8,0.2]; A['Management'][c]['Data']=MosByMP[iScn]['Area']['Planting']['Ensemble Mean'][:,iMP]+MosByMP[iScn]['Area']['Direct Seeding']['Ensemble Mean'][:,iMP]
    c=c+1; A['Management'][c]={}; A['Management'][c]['Name']='Foliage protection'; A['Management'][c]['Color']=[1,0.7,0]; A['Management'][c]['Data']=MosByMP[iScn]['Area']['IDW Btk Spray']['Ensemble Mean'][:,iMP]
    c=c+1; A['Management'][c]={}; A['Management'][c]['Name']='Aerial nutrient application'; A['Management'][c]['Color']=[0.65,0,1]; A['Management'][c]['Data']=MosByMP[iScn]['Area']['Aerial BTK Spray Aerial']['Ensemble Mean'][:,iMP]
    A['Management']=A['Management'][0:c+1]

    # Convert to x-year intervals
    A['tv']=gu.BlockMean(tv,ivlT)
    for i in range(len(A['Nat Dist'])):
        A['Nat Dist'][i]['Data']=gu.BlockMean(A['Nat Dist'][i]['Data'],ivlT)
    for i in range(len(A['Management'])):
        A['Management'][i]['Data']=gu.BlockMean(A['Management'][i]['Data'],ivlT)

    return A

#%% Post process BatchTIPSY output

def PrepGrowthCurvesForCBR(meta,pNam):

    # Function used to smooth curves
    def smooth(y, box_pts):
        box = np.ones(box_pts)/box_pts
        y_smooth = np.convolve(y, box, mode='same')
        return y_smooth

    # TIPSY Export curves as MgDM/ha/yr, CBRunner expects inputs of MgC/ha/yr - create conversion factor
    dm2c=0.5

    # Growth curve parameters and TIPSY outputs
    dfPar=pd.read_excel(meta['Paths'][pNam]['Data'] + '\\Inputs\\GrowthCurvesTIPSY_Parameters.xlsx',sheet_name='Sheet1',skiprows=6)
    txtDat=np.loadtxt(meta['Paths'][pNam]['Data'] + '\\Inputs\\GrowthCurvesTIPSY_Output.out',skiprows=4)

    # TIPSY saves to text file -> convert to dataframe (column names must match TIPSY output file design)
    dfDat=pd.DataFrame(txtDat,columns=meta['Modules']['GYM']['BatchTIPSY Column Names'])

    # Define age vector (must be consistent with how TIPSY was set up)
    Age=np.arange(0,meta['Modules']['GYM']['BatchTIPSY Maximum Age']+1,1)

    # Get dimensions of the TIPSY output file to reshape the data into Age x Stand
    N_Age=Age.size
    N_GC=int(dfDat.shape[0]/N_Age)

    # Stemwood

    # Merchantable stemwood volume
    V_Merch=np.reshape(dfDat['VolMerch125'].values,(N_Age,N_GC),order='F')
    V_StemTot=np.reshape(dfDat['VolTot0'].values,(N_Age,N_GC),order='F')
    G_VStemMerch=np.append(np.zeros((1,N_GC)),np.diff(V_Merch,axis=0),axis=0)

    # Extract age responses for each biomass pool
    C_Stem=dm2c*np.reshape(dfDat['ODT_Stem'].values,(N_Age,N_GC),order='F')

    #import matplotlib.pyplot as plt
    #plt.close('all')
    #plt.plot(C_Stem[:,0])

    # Apply smoothing - it messes up the last ten years so don't smooth that part
    for j in range(C_Stem.shape[1]):
        a=smooth(C_Stem[:,j],10)
        C_Stem[:-10,j]=a[:-10]
        a=smooth(V_Merch[:,j],10)
        V_Merch[:-10,j]=a[:-10]
        a=smooth(V_StemTot[:,j],10)
        V_StemTot[:-10,j]=a[:-10]
    #plt.plot(C_Stem[:,0],'--')

    # Define the fraction of merchantable stemwood
    fMerch=np.nan_to_num(V_Merch/V_StemTot)
    fNonMerch=1-fMerch

    # Calculate growth
    z=np.zeros((1,N_GC))
    G_Stem=np.append(z,np.diff(C_Stem,axis=0),axis=0)

    # Adjust early net growth, but don't change the total stock change
    A_th=30
    ind=np.where(Age<=A_th)[0]
    bin=np.arange(0.005,0.15,0.005)
    x=np.arange(0,A_th+1,1)
    for j in range(N_GC):
        Gtot=C_Stem[ind[-1],j]
        y_th=G_Stem[ind[-1],j]
        Gtot_hat=1000*np.ones(bin.size)
        for k in range(bin.size):
            Gtot_hat[k]=np.sum(y_th*np.exp(bin[k]*(x-A_th)))
        ind1=np.where(np.abs(Gtot_hat-Gtot)==np.min(np.abs(Gtot_hat-Gtot)))[0]
        G_Stem[ind,j]=y_th*np.exp(bin[ind1[0]]*(x-A_th))

    # Update merch and nonmerch growth
    C_Stem=np.cumsum(G_Stem,axis=0)
    C_StemMerch=fMerch*C_Stem
    C_StemNonMerch=fNonMerch*C_Stem

    G_StemMerch=np.append(z,np.diff(C_StemMerch,axis=0),axis=0)
    G_StemNonMerch=np.append(z,np.diff(C_StemNonMerch,axis=0),axis=0)

    # Fix growth of year zero
    G_StemMerch[0,:]=G_StemMerch[1,:]
    G_StemNonMerch[0,:]=G_StemNonMerch[1,:]

    # Add negative nonmerch to merch
    ind=np.where(G_StemNonMerch<0)
    G_StemMerch[ind]=G_StemMerch[ind]+G_StemNonMerch[ind]
    G_StemNonMerch[ind]=0

    #plt.close('all')
    #plt.plot(np.maximum(-1,G_Stem[:,1]))
    #plt.plot(np.maximum(-1,G_StemMerch[:,1]),'--')
    #plt.plot(np.maximum(-1,G_StemNonMerch[:,1]),'-.')

    # Other pools - these are no longer being used. Net growth of non-stemwood
    # biomass is simulated in the annual loop based on allometric reatlionships
    # with stemwood net growth and stand age.

    # Foliage biomass is very low, revise
    #bF1=0.579
    #bF2=0.602
    #C_Foliage=np.maximum(0,bF1*C_Stem**bF2)

    C_Foliage=dm2c*np.reshape(dfDat['ODT_Foliage'].values,(N_Age,N_GC),order='F')
    C_Branch=dm2c*np.reshape(dfDat['ODT_Branch'].values,(N_Age,N_GC),order='F')
    C_Bark=dm2c*np.reshape(dfDat['ODT_Bark'].values,(N_Age,N_GC),order='F')

    G_Foliage=np.append(z,np.diff(C_Foliage,axis=0),axis=0)
    G_Branch=np.append(z,np.diff(C_Branch,axis=0),axis=0)
    G_Bark=np.append(z,np.diff(C_Bark,axis=0),axis=0)

    for iScn in range(meta[pNam]['Project']['N Scenario']):

        # Define growth curves (from TIPSY output)
        for iBat in range(meta[pNam]['Project']['N Batch']):

            # Index to batch
            indBat=IndexToBatch(meta[pNam],iBat)

            for iGC in range(3):

                # Initialize age response of net growth
                G=np.zeros((N_Age,indBat.size,6),dtype='int32')

                # Populate the growth curve
                for iS in range(indBat.size):

                    #u=np.unique(ec['ID Growth Curve'][:,iS,:])

                    if (meta[pNam]['Project']['Scenario Source']=='Spreadsheet'):

                        indTIPSY=np.where(
                                (dfPar['ID_Scenario']==iScn+1) &
                                (dfPar['ID_GC']==int(meta['Modules']['GYM']['ID GC Unique'][iGC])) )[0]

                    elif (meta[pNam]['Project']['Scenario Source']=='Script') | (meta[pNam]['Project']['Scenario Source']=='Portfolio'):

                        indTIPSY=np.where(
                            (dfPar['ID_Stand']==indBat[iS]+1) &
                            (dfPar['ID_Scenario']==iScn+1) &
                            (dfPar['ID_GC']==int(iGC+1)))[0]

                    if (indTIPSY.size==0):
                        # This can happen if only some stands have a third GC, for example
                        continue

                    G[:,iS,0]=G_StemMerch[:,indTIPSY[0]]/meta['Modules']['GYM']['Scale Factor']
                    G[:,iS,1]=G_StemNonMerch[:,indTIPSY[0]]/meta['Modules']['GYM']['Scale Factor']
                    G[:,iS,2]=G_Bark[:,indTIPSY[0]]/meta['Modules']['GYM']['Scale Factor']
                    G[:,iS,3]=G_Branch[:,indTIPSY[0]]/meta['Modules']['GYM']['Scale Factor']
                    G[:,iS,4]=G_Foliage[:,indTIPSY[0]]/meta['Modules']['GYM']['Scale Factor']
                    G[:,iS,5]=G_VStemMerch[:,indTIPSY[0]]/meta['Modules']['GYM']['Scale Factor']

                # Save data to file in input variables folder of project
                gu.opickle(meta['Paths'][pNam]['Data'] + '\\Inputs\\Scenario' + FixFileNum(iScn) + '\\GrowthCurve' + str(iGC+1) + '_Bat' + FixFileNum(iBat) + '.pkl',G)

    return

#%% Prepare growth curves (with early correction)

def PrepGrowthCurvesUniqueForCBR(meta,pNam,ugc):

    # *** This adjusts early net growth so that it is not zero. There is a
    # second version of this function that excludes the correction. ***

    # TIPSY Export curves as MgDM/ha/yr, CBRunner expects inputs of MgC/ha/yr. Create
    # conversion factor.
    dm2c=0.5

    # Growth curve parameters and TIPSY outputs
    #dfPar=pd.read_excel(meta['Paths'][pNam]['Data'] + '\\Inputs\\GrowthCurvesTIPSY_Parameters.xlsx',sheet_name='Sheet1',skiprows=7)
    txtDat=np.loadtxt(meta['Paths'][pNam]['Data'] + '\\Inputs\\GrowthCurvesTIPSY_Output.out',skiprows=4)

    # TIPSY saves to text file -> convert to dataframe (column names must match TIPSY output file design)
    dfDat=pd.DataFrame(txtDat,columns=meta['Modules']['GYM']['BatchTIPSY Column Names'])

    del txtDat

    # Define age vector (must be consistent with how TIPSY was set up)
    Age=np.arange(0,meta['Modules']['GYM']['BatchTIPSY Maximum Age']+1,1)

    # Get dimensions of the TIPSY output file to reshape the data into Age x Stand
    N_Age=Age.size
    N_GC=int(dfDat.shape[0]/N_Age)

    # Define the fraction of merchantable stemwood
    fMerch=np.nan_to_num(np.reshape(dfDat['VolMerch125'].values,(N_Age,N_GC),order='F')/np.reshape(dfDat['VolTot0'].values,(N_Age,N_GC),order='F'))
    fNonMerch=1-fMerch

    ## Merchantable stemwood volume
    #V_Merch=np.reshape(dfDat['VolMerch125'].values,(N_Age,N_GC),order='F')
    #G_VStemMerch=np.append(np.zeros((1,N_GC)),np.diff(V_Merch,axis=0),axis=0)
    #
    ## Extract age responses for each biomass pool
    #C_Stem=dm2c*np.reshape(dfDat['ODT_Stem'].values,(N_Age,N_GC),order='F')
    #C_StemMerch=fMerch*C_Stem
    #C_StemNonMerch=fNonMerch*C_Stem
    #C_Foliage=dm2c*np.reshape(dfDat['ODT_Foliage'].values,(N_Age,N_GC),order='F')
    #C_Branch=dm2c*np.reshape(dfDat['ODT_Branch'].values,(N_Age,N_GC),order='F')
    #C_Bark=dm2c*np.reshape(dfDat['ODT_Bark'].values,(N_Age,N_GC),order='F')
    #
    ## Calculate growth
    #z=np.zeros((1,N_GC))
    #G_StemMerch=np.append(z,np.diff(C_StemMerch,axis=0),axis=0)
    #G_StemNonMerch=np.append(z,np.diff(C_StemNonMerch,axis=0),axis=0)
    #G_Foliage=np.append(z,np.diff(C_Foliage,axis=0),axis=0)
    #G_Branch=np.append(z,np.diff(C_Branch,axis=0),axis=0)
    #G_Bark=np.append(z,np.diff(C_Bark,axis=0),axis=0)
    #
    #plt.plot(np.mean(G_StemMerch+G_StemNonMerch,axis=1))
    ##plt.plot(np.percentile(G_StemMerch+G_StemNonMerch,20,axis=1))
    ##plt.plot(np.percentile(G_StemMerch+G_StemNonMerch,80,axis=1))

    # Function used to smooth curves
    def smooth(y, box_pts):
        box = np.ones(box_pts)/box_pts
        y_smooth = np.convolve(y, box, mode='same')
        return y_smooth

    # Stemwood

    # Merchantable stemwood volume
    V_Merch=np.reshape(dfDat['VolMerch125'].values,(N_Age,N_GC),order='F')
    V_StemTot=np.reshape(dfDat['VolTot0'].values,(N_Age,N_GC),order='F')
    G_VStemMerch=np.append(np.zeros((1,N_GC)),np.diff(V_Merch,axis=0),axis=0)

    # Extract age responses for each biomass pool
    C_Stem=dm2c*np.reshape(dfDat['ODT_Stem'].values,(N_Age,N_GC),order='F')

    #import matplotlib.pyplot as plt
    #plt.close('all')
    #plt.plot(C_Stem[:,0])

    # Apply smoothing - it messes up the last ten years so don't smooth that part
    for j in range(C_Stem.shape[1]):
        a=smooth(C_Stem[:,j],10)
        C_Stem[:-10,j]=a[:-10]
        a=smooth(V_Merch[:,j],10)
        V_Merch[:-10,j]=a[:-10]
        a=smooth(V_StemTot[:,j],10)
        V_StemTot[:-10,j]=a[:-10]
        #plt.plot(C_Stem[:,0],'--')

    # Define the fraction of merchantable stemwood
    fMerch=np.nan_to_num(V_Merch/V_StemTot)
    fNonMerch=1-fMerch

    # Calculate growth
    z=np.zeros((1,N_GC))
    G_Stem=np.append(z,np.diff(C_Stem,axis=0),axis=0)

    # Adjust early net growth, but don't change the total stock change
    A_th=30
    ind=np.where(Age<=A_th)[0]
    bin=np.arange(0.005,0.15,0.005)
    x=np.arange(0,A_th+1,1)
    for j in range(N_GC):
        Gtot=C_Stem[ind[-1],j]
        y_th=G_Stem[ind[-1],j]
        Gtot_hat=1000*np.ones(bin.size)
        for k in range(bin.size):
            Gtot_hat[k]=np.sum(y_th*np.exp(bin[k]*(x-A_th)))
        ind1=np.where(np.abs(Gtot_hat-Gtot)==np.min(np.abs(Gtot_hat-Gtot)))[0]
        G_Stem[ind,j]=y_th*np.exp(bin[ind1[0]]*(x-A_th))

    # Update merch and nonmerch growth
    C_Stem=np.cumsum(G_Stem,axis=0)
    C_StemMerch=fMerch*C_Stem
    C_StemNonMerch=fNonMerch*C_Stem
    G_StemMerch=np.append(z,np.diff(C_StemMerch,axis=0),axis=0)
    G_StemNonMerch=np.append(z,np.diff(C_StemNonMerch,axis=0),axis=0)

    # Fix growth of year zero
    G_StemMerch[0,:]=G_StemMerch[1,:]
    G_StemNonMerch[0,:]=G_StemNonMerch[1,:]

    # Add negative nonmerch to merch
    ind=np.where(G_StemNonMerch<0)
    G_StemMerch[ind]=G_StemMerch[ind]+G_StemNonMerch[ind]
    G_StemNonMerch[ind]=0

    # Other pools - these are no longer being used. Net growth of non-stemwood
    # biomass is simulated in the annual loop based on allometric reatlionships
    # with stemwood net growth and stand age.

    # Foliage biomass is very low, revise
    #bF1=0.579
    #bF2=0.602
    #C_Foliage=np.maximum(0,bF1*C_Stem**bF2)

    C_Foliage=dm2c*np.reshape(dfDat['ODT_Foliage'].values,(N_Age,N_GC),order='F')
    C_Branch=dm2c*np.reshape(dfDat['ODT_Branch'].values,(N_Age,N_GC),order='F')
    C_Bark=dm2c*np.reshape(dfDat['ODT_Bark'].values,(N_Age,N_GC),order='F')

    G_Foliage=np.append(z,np.diff(C_Foliage,axis=0),axis=0)
    G_Branch=np.append(z,np.diff(C_Branch,axis=0),axis=0)
    G_Bark=np.append(z,np.diff(C_Bark,axis=0),axis=0)

    del C_Stem,C_StemMerch,C_StemNonMerch,C_Foliage,C_Branch,C_Bark,fMerch,fNonMerch

    for iScn in range(meta[pNam]['Project']['N Scenario']):

        for iGC in range(meta['Modules']['GYM']['N Growth Curves']):

            # Index to the full set of growth curves for scenario iScn and growth curve iGC
            ind_ugc_ScnAndGc=np.where( (ugc['Full'][:,1]==iScn) & (ugc['Full'][:,2]==meta['Modules']['GYM']['ID GC Unique'][iGC]) )[0]

            # Extract the unique growth curve ID for scenario iScn and growth curve iGC
            ID_ugc_ScnAndGc=ugc['Full'][ind_ugc_ScnAndGc,0]

            # Extract the inverse index for scenario iScn and growth curve iGC
            Inverse_ugc_ScnAndGc=ugc['Inverse'][ind_ugc_ScnAndGc]

            for iBat in range(0,meta[pNam]['Project']['N Batch']):

                # Index to batch
                indBat=IndexToBatch(meta[pNam],iBat)

                # Intersect
                c,inda,indb=np.intersect1d(ID_ugc_ScnAndGc,indBat,return_indices=True)

                Inverse_ugc_ScnAndGcAndBat=Inverse_ugc_ScnAndGc[inda]

                # Initialize array of growth data
                G=np.zeros((N_Age,indBat.size,6),dtype='int16')

                for i in range(inda.size):

                    iStand=indb[i]
                    iGC_Unique=Inverse_ugc_ScnAndGcAndBat[i]

                    try:
                        # This crashes often so leave it in try for troubleshooting
                        # *** It's a sign that the BatchTIPSY inputs are pointing to the wrong directory ***
                        G[:,iStand,0]=G_StemMerch[:,iGC_Unique].T/meta['Modules']['GYM']['Scale Factor']
                    except:
                        print(iStand)
                        print(G.shape)
                        print(iGC_Unique)
                        print(G_StemMerch.shape)
                    G[:,iStand,1]=G_StemNonMerch[:,iGC_Unique].T/meta['Modules']['GYM']['Scale Factor']
                    G[:,iStand,2]=G_Bark[:,iGC_Unique].T/meta['Modules']['GYM']['Scale Factor']
                    G[:,iStand,3]=G_Branch[:,iGC_Unique].T/meta['Modules']['GYM']['Scale Factor']
                    G[:,iStand,4]=G_Foliage[:,iGC_Unique].T/meta['Modules']['GYM']['Scale Factor']
                    G[:,iStand,5]=G_VStemMerch[:,iGC_Unique].T/meta['Modules']['GYM']['Scale Factor']

                # Save data to file in input variables folder of project
                gu.opickle(meta['Paths'][pNam]['Data'] + '\\Inputs\\Scenario' + FixFileNum(iScn) + '\\GrowthCurve' + str(meta['Modules']['GYM']['ID GC Unique'][iGC]) + '_Bat' + FixFileNum(iBat) + '.pkl',G)

    return

#%% Prepare growth curves (without early correction)

def PrepGrowthCurvesUniqueForCBR_WithoutEarlyCorrection(meta,pNam,ugc):

    # TIPSY exports curves as MgDM/ha/yr, CBRunner expects inputs of MgC/ha/yr. Create
    # conversion factor.
    dm2c=0.5

    # Growth curve parameters and TIPSY outputs
    #dfPar=pd.read_excel(meta['Paths'][pNam]['Data'] + '\\Inputs\\GrowthCurvesTIPSY_Parameters.xlsx',sheet_name='Sheet1',skiprows=7)
    txtDat=np.loadtxt(meta['Paths'][pNam]['Data'] + '\\Inputs\\GrowthCurvesTIPSY_Output.out',skiprows=4)

    # TIPSY saves to text file -> convert to dataframe (column names must match TIPSY output file design)
    dfDat=pd.DataFrame(txtDat,columns=meta['Modules']['GYM']['BatchTIPSY Column Names'])

    del txtDat

    # Define age vector (must be consistent with how TIPSY was set up)
    Age=np.arange(0,meta['Modules']['GYM']['BatchTIPSY Maximum Age']+1,1)

    # Get dimensions of the TIPSY output file to reshape the data into Age x Stand
    N_Age=Age.size
    N_GC=int(dfDat.shape[0]/N_Age)

    # Define the fraction of merchantable stemwood
    fMerch=np.nan_to_num(np.reshape(dfDat['VolMerch125'].values,(N_Age,N_GC),order='F')/np.reshape(dfDat['VolTot0'].values,(N_Age,N_GC),order='F'))
    fNonMerch=1-fMerch

    # Merchantable stemwood volume
    V_Merch=np.reshape(dfDat['VolMerch125'].values,(N_Age,N_GC),order='F')
    G_VStemMerch=np.append(np.zeros((1,N_GC)),np.diff(V_Merch,axis=0),axis=0)

    # Extract age responses for each biomass pool
    C_Stem=dm2c*np.reshape(dfDat['ODT_Stem'].values,(N_Age,N_GC),order='F')
    C_StemMerch=fMerch*C_Stem
    C_StemNonMerch=fNonMerch*C_Stem
    C_Foliage=dm2c*np.reshape(dfDat['ODT_Foliage'].values,(N_Age,N_GC),order='F')
    C_Branch=dm2c*np.reshape(dfDat['ODT_Branch'].values,(N_Age,N_GC),order='F')
    C_Bark=dm2c*np.reshape(dfDat['ODT_Bark'].values,(N_Age,N_GC),order='F')

    # Calculate growth
    z=np.zeros((1,N_GC))
    G_StemMerch=np.append(z,np.diff(C_StemMerch,axis=0),axis=0)
    G_StemNonMerch=np.append(z,np.diff(C_StemNonMerch,axis=0),axis=0)
    G_Foliage=np.append(z,np.diff(C_Foliage,axis=0),axis=0)
    G_Branch=np.append(z,np.diff(C_Branch,axis=0),axis=0)
    G_Bark=np.append(z,np.diff(C_Bark,axis=0),axis=0)

    # Fix growth of year zero
    G_StemMerch[0,:]=G_StemMerch[1,:]
    G_StemNonMerch[0,:]=G_StemNonMerch[1,:]
    G_Foliage[0,:]=G_Foliage[1,:]
    G_Branch[0,:]=G_Branch[1,:]
    G_Bark[0,:]=G_Bark[1,:]

    del C_Stem,C_StemMerch,C_StemNonMerch,C_Foliage,C_Branch,C_Bark,fMerch,fNonMerch

    for iScn in range(meta[pNam]['Project']['N Scenario']):
        for iGC in range(meta['Modules']['GYM']['N Growth Curves']):

            # Index to the full set of growth curves for scenario iScn and growth curve iGC
            ind_ugc_ScnAndGc=np.where( (ugc['Full'][:,1]==iScn) & (ugc['Full'][:,2]==meta['Modules']['GYM']['ID GC Unique'][iGC]) )[0]

            # Extract the unique growth curve ID for scenario iScn and growth curve iGC
            ID_ugc_ScnAndGc=ugc['Full'][ind_ugc_ScnAndGc,0]

            # Extract the inverse index for scenario iScn and growth curve iGC
            Inverse_ugc_ScnAndGc=ugc['Inverse'][ind_ugc_ScnAndGc]

            for iBat in range(0,meta[pNam]['Project']['N Batch']):

                # Index to batch
                indBat=IndexToBatch(meta[pNam],iBat)

                # Intersect
                c,inda,indb=np.intersect1d(ID_ugc_ScnAndGc,indBat,return_indices=True)

                Inverse_ugc_ScnAndGcAndBat=Inverse_ugc_ScnAndGc[inda]

                # Initialize array of growth data
                G=np.zeros((N_Age,indBat.size,6),dtype='int16')

                for i in range(inda.size):

                    iStand=indb[i]
                    iGC_Unique=Inverse_ugc_ScnAndGcAndBat[i]

                    G[:,iStand,0]=G_StemMerch[:,iGC_Unique].T/meta['Modules']['GYM']['Scale Factor']
                    G[:,iStand,1]=G_StemNonMerch[:,iGC_Unique].T/meta['Modules']['GYM']['Scale Factor']
                    G[:,iStand,2]=G_Bark[:,iGC_Unique].T/meta['Modules']['GYM']['Scale Factor']
                    G[:,iStand,3]=G_Branch[:,iGC_Unique].T/meta['Modules']['GYM']['Scale Factor']
                    G[:,iStand,4]=G_Foliage[:,iGC_Unique].T/meta['Modules']['GYM']['Scale Factor']
                    G[:,iStand,5]=G_VStemMerch[:,iGC_Unique].T/meta['Modules']['GYM']['Scale Factor']

                # Save data to file in input variables folder of project
                gu.opickle(meta['Paths'][pNam]['Data'] + '\\Inputs\\Scenario' + FixFileNum(iScn) + '\\GrowthCurve' + str(meta['Modules']['GYM']['ID GC Unique'][iGC]) + '_Bat' + FixFileNum(iBat) + '.pkl',G)

    return

#%% Import parameters

def ImportParameters(meta):

    # Path to parameters
    pthin=meta['Paths']['Model']['Code'] + '\\Parameters\\'

    # Initialize parameter structure
    if 'Param' not in meta:
        meta['Param']={}

    # Best estimates, variance, lower and upper confidence levels
    meta['Param']['BE']={}
    meta['Param']['Sigma']={}
    meta['Param']['BL']={}
    meta['Param']['BU']={}

    #--------------------------------------------------------------------------
    # Biophysical
    #--------------------------------------------------------------------------

    df=pd.read_excel(pthin + '\Parameters_Biophysical.xlsx',sheet_name='Sheet1')
    meta['Param']['BE']['Biophysical']={}
    for i in range(len(df)):
        Name=df['Name'][i]
        meta['Param']['BE']['Biophysical'][Name]=df['Value'][i]

    #--------------------------------------------------------------------------
    # Biomass - allometry
    # Values are location specific so no specific processing is done here (see cbrun.py)
    #--------------------------------------------------------------------------

    meta['Param']['BE']['Biomass Allometry']={}
    meta['Param']['BE']['Biomass Allometry']['Raw']=gu.ReadExcel(pthin + '\Parameters_BiomassAllometrySL.xlsx')

    #--------------------------------------------------------------------------
    # Biomass - turnover
    #--------------------------------------------------------------------------

    df=pd.read_excel(pthin + '\Parameters_BiomassTurnover.xlsx',sheet_name='Sheet1')
    meta['Param']['BE']['Biomass Turnover']={}
    meta['Param']['Sigma']['Biomass Turnover']={}
    meta['Param']['BL']['Biomass Turnover']={}
    meta['Param']['BU']['Biomass Turnover']={}
    for i in range(len(df)):
        Name=df['Name'][i]
        meta['Param']['BE']['Biomass Turnover'][Name]=df['Best Estimate'][i]
        meta['Param']['Sigma']['Biomass Turnover'][Name]=df['Sigma'][i]
        meta['Param']['BL']['Biomass Turnover'][Name]=df['BL'][i]
        meta['Param']['BU']['Biomass Turnover'][Name]=df['BU'][i]

    #--------------------------------------------------------------------------
    # By BGC zone
    #--------------------------------------------------------------------------

    meta['Param']['BE']['BGC Zone Averages']=gu.ReadExcel(pthin + '\Parameters_ByBGC.xlsx')

    #--------------------------------------------------------------------------
    # Decomposition parameters (Kurz et al. 2009)
    #--------------------------------------------------------------------------

    df=pd.read_excel(pthin + '\Parameters_Decomposition.xlsx',sheet_name='Sheet1')
    meta['Param']['BE']['Decomp']={}
    meta['Param']['Sigma']['Decomp']={}
    meta['Param']['BL']['Decomp']={}
    meta['Param']['BU']['Decomp']={}
    for i in range(len(df)):
        Name=df['Name'].iloc[i]
        meta['Param']['BE']['Decomp'][Name]=df['Best Estimate'].iloc[i]
        meta['Param']['Sigma']['Decomp'][Name]=df['Sigma'].iloc[i]
        meta['Param']['BL']['Decomp'][Name]=df['BL'].iloc[i]
        meta['Param']['BU']['Decomp'][Name]=df['BU'].iloc[i]

    #--------------------------------------------------------------------------
    # Disturbance - static parameters
    #--------------------------------------------------------------------------

    p=gu.ReadExcel(pthin + '\Parameters_Disturbances.xlsx')

    str_to_exclude=['ID','Name','Species_CD','MortalityOccurs','Growth Factor',
                    'GrowthFactor_Source','GrowthRecovery_HL','GrowthRecovery_HL_Source',
                    'QA1','QA2','QA3']

    meta['Param']['BE']['Event']={}
    for i in range(p['ID'].size):
        meta['Param']['BE']['Event'][p['ID'][i]]={}
        #meta['Param']['BE']['Event'][p['ID'][i]]['BiomassMerch_Affected']=1
        #meta['Param']['BE']['Event'][p['ID'][i]]['BiomassNonMerch_Affected']=1
        #meta['Param']['BE']['Event'][p['ID'][i]]['Snags_Affected']=1
        for k in p.keys():
            # Exclude some legacy variables that may be re-implemented
            if np.isin(k,str_to_exclude)==True:
                continue
            meta['Param']['BE']['Event'][p['ID'][i]][k]=np.nan_to_num(p[k][i])

    #--------------------------------------------------------------------------
    # Disutrbance - fate of felled material
    # Values are location specific so no specific processing is done here (see cbrun.py)
    #--------------------------------------------------------------------------

    meta['Param']['BE']['Felled Fate']={}
    meta['Param']['BE']['Felled Fate']=gu.ipickle(meta['Paths']['Model']['Code'] + '\\Parameters\\Variables_FelledFate.pkl')

    #--------------------------------------------------------------------------
    # Disturbance - Wildfire aspatial (Taz-AAO) parameters
    #--------------------------------------------------------------------------

    meta['Param']['BE']['Taz']={}

    wf=gu.ReadExcel(meta['Paths']['Model']['Code'] + '\\Parameters\\Parameters_WildfireStatsMod.xlsx')

    meta['Param']['BE']['Taz']['WF']={}
    for i in range(wf['Name'].size):
        try:
            meta['Param']['BE']['Taz']['WF'][wf['Name'][i]]=wf['Value'][i].astype(float)
        except:
            meta['Param']['BE']['Taz']['WF'][wf['Name'][i]]=wf['Value'][i]

    #--------------------------------------------------------------------------
    # Disturbance - Mountain Pine Beetle (Taz-AAO) parameters
    #--------------------------------------------------------------------------

    ibm=gu.ReadExcel(meta['Paths']['Model']['Code'] + '\\Parameters\\Parameters_IBMStatsMod.xlsx')

    meta['Param']['BE']['Taz']['IBM']={}
    for i in range(ibm['Name'].size):
        try:
            meta['Param']['BE']['Taz']['IBM'][ibm['Name'][i]]=ibm['Value'][i].astype(float)
        except:
            meta['Param']['BE']['Taz']['IBM'][ibm['Name'][i]]=ibm['Value'][i]

    #--------------------------------------------------------------------------
    # Disturbance - On the fly parameters
    #--------------------------------------------------------------------------

    df=pd.read_excel(pthin + '\Parameters_OnTheFly.xlsx',sheet_name='Sheet1')
    meta['Param']['BE']['On The Fly']={}
    for i in range(len(df)):
        Name=df['Name'].iloc[i]
        meta['Param']['BE']['On The Fly'][Name]=df['Value'].iloc[i]

    #--------------------------------------------------------------------------
    # Disturbance - harvesting - basic early historical reconstruction
    #--------------------------------------------------------------------------

    fin=meta['Paths']['Model']['Taz Datasets'] + '\\Harvest Stats and Scenarios\\HarvestHistoricalProbabilitySimple.xlsx'
    meta['Param']['BE']['Taz']['Ph_Simp']=gu.ReadExcel(fin)

    #--------------------------------------------------------------------------
    # Disturbance - By severity class
    #--------------------------------------------------------------------------

    meta['Param']['BE']['DistBySC']=gu.ReadExcel(meta['Paths']['Model']['Code'] + '\\Parameters\\Parameters_DisturbanceBySeverityClass.xlsx')

    #--------------------------------------------------------------------------
    # Economics
    #--------------------------------------------------------------------------

    df=pd.read_excel(pthin + '\Parameters_Economics.xlsx',sheet_name='Sheet1')
    meta['Param']['BE']['Econ']={}
    for i in range(len(df)):
        Name=df['Name'].iloc[i]
        meta['Param']['BE']['Econ'][Name]=df['Value'].iloc[i]

    #--------------------------------------------------------------------------
    # Funding source code
    #--------------------------------------------------------------------------

    meta['Param']['BE']['FSC']={}
    meta['Param']['BE']['FSC']['Raw']=gu.ReadExcel(meta['Paths']['Model']['Code'] + '\\Parameters\\Parameters_ByFundingSourceCode.xlsx')

    ind=np.where(meta['Param']['BE']['FSC']['Raw']['Legal Obligation of Licensee To Fund Milestone']=='No')[0]
    meta['Param']['BE']['FSC']['NO List Name']=meta['Param']['BE']['FSC']['Raw']['Name'][ind]

    meta['Param']['BE']['FSC']['NO List ID']=meta['LUT']['RSLT_ACTIVITY_TREATMENT_SVW']['SILV_FUND_SOURCE_CODE']['FTM']*np.ones(meta['Param']['BE']['FSC']['NO List Name'].size,dtype='int16')
    for i in range(meta['Param']['BE']['FSC']['NO List Name'].size):
        try:
            meta['Param']['BE']['FSC']['NO List ID'][i]=meta['LUT']['RSLT_ACTIVITY_TREATMENT_SVW']['SILV_FUND_SOURCE_CODE'][ meta['Param']['BE']['FSC']['NO List Name'][i] ]
        except:
            pass

    #--------------------------------------------------------------------------
    # HWP - static parameters
    # Values are location specific so no specific processing is done here (see cbrun.py)
    #--------------------------------------------------------------------------

    meta['Param']['BE']['HWP']={}
    meta['Param']['BE']['HWP']['Raw']=pd.read_excel(pthin + '\Parameters_HWP.xlsx',sheet_name='Main')

    #--------------------------------------------------------------------------
    # HWP - End Use parameters
    # Values are location specific so no specific processing is done here (see cbrun.py)
    #--------------------------------------------------------------------------

    meta['Param']['BE']['HWP End Use']={}
    meta['Param']['BE']['HWP End Use']=gu.ipickle(meta['Paths']['Model']['Code'] + '\\Parameters\\Variables_HWP_EndUse.pkl')

    #--------------------------------------------------------------------------
    # Genetic worth
    #--------------------------------------------------------------------------

    meta['Param']['BE']['Genetic Worth']=gu.ReadExcel(pthin + '\\Parameters_Seedlot_GW.xlsx')

    #--------------------------------------------------------------------------
    # Inter-pool transfer parameters (Kurz et al. 2009)
    #--------------------------------------------------------------------------

    df=pd.read_excel(pthin + '\Parameters_InterPoolFluxes.xlsx',sheet_name='Sheet1')
    meta['Param']['BE']['Inter Pool Fluxes']={}
    meta['Param']['Sigma']['Inter Pool Fluxes']={}
    for i in range(len(df)):
        Name=df['Name'][i]
        meta['Param']['BE']['Inter Pool Fluxes'][Name]=df['Best Estimate'][i]
        meta['Param']['Sigma']['Inter Pool Fluxes'][Name]=df['Sigma'][i]

    #--------------------------------------------------------------------------
    # Nutrient application paramaters
    #--------------------------------------------------------------------------

    df=pd.read_excel(pthin + '\Parameters_NutrientApplication.xlsx',sheet_name='Sheet1')
    meta['Param']['BE']['Nutrient Management']={}
    meta['Param']['Sigma']['Nutrient Management']={}
    for i in range(len(df)):
        Name=df['Name'].iloc[i]
        meta['Param']['BE']['Nutrient Management'][Name]=df['Value'].iloc[i]
        meta['Param']['Sigma']['Nutrient Management'][Name]=df['Sigma'].iloc[i]

    #--------------------------------------------------------------------------
    # Removed fate scenarios
    # Values are location specific so no specific processing is done here (see cbrun.py)
    #--------------------------------------------------------------------------

    meta['Param']['BE']['Removed Fate']=gu.ipickle(meta['Paths']['Model']['Code'] + '\\Parameters\\Variables_RemovedFate.pkl')

    #--------------------------------------------------------------------------
    # Substitution effects (use of basic displacement factors)
    #--------------------------------------------------------------------------

    df=pd.read_excel(pthin + '\\Parameters_Substitution.xlsx',sheet_name='Sheet1')
    meta['Param']['BE']['Substitution']={}
    meta['Param']['Sigma']['Substitution']={}
    for i in range(len(df)):
        Name=df['Name'].iloc[i]
        meta['Param']['BE']['Substitution'][Name]=df['Value'].iloc[i]
        meta['Param']['Sigma']['Substitution'][Name]=df['Sigma'].iloc[i]

    #--------------------------------------------------------------------------
    # By Tree Density Class
    #--------------------------------------------------------------------------

    meta['Param']['BE']['Tree Density Class Averages']=gu.ReadExcel(pthin + '\Parameters_TreeDensityClassStatistics.xlsx')

    #--------------------------------------------------------------------------
    # Wood density by species
    #--------------------------------------------------------------------------

    meta['Param']['BE']['Wood Density']=gu.ReadExcel(pthin + '\Parameters_WoodDensity.xlsx')

    #--------------------------------------------------------------------------
    # Sawtooth parameters
    #--------------------------------------------------------------------------

    #if meta[pNam]['Project']['Biomass Module']=='Sawtooth':

    ps={}
    #pthin=r'C:\Users\rhember\Documents\Code_Python\fcgadgets\cbrunner\Parameters'

    # Core sawtooth parameters
    d=gu.ReadExcel(pthin + '\Parameters_Sawtooth.xlsx','Core')
    ps['Core']={}
    for i in range(d['Name'].size):
        ps['Core'][d['Name'][i]]=d['Value'][i]

    # Sawtooth key between provincial species codes adn species-region sample
    # codes
    ps['SRS Key']=gu.ReadExcel(pthin + '\Parameters_Sawtooth.xlsx','SRS Key')

    # Tree-level allometry
    ps['Allom']=gu.ReadExcel(pthin + '\Parameters_Sawtooth.xlsx','Allometry')
    ps['Eq R']={}
    ps['Eq R']['Def1']=gu.ReadExcel(pthin + '\Parameters_Sawtooth.xlsx','R_Def1')

    ps['Eq M']={}
    ps['Eq M']['Def1']=gu.ReadExcel(pthin + '\Parameters_Sawtooth.xlsx','M_Def1')

    ps['Eq G']={}
    ps['Eq G']['Def1']=gu.ReadExcel(pthin + '\Parameters_Sawtooth.xlsx','G_Def1')

    meta['Param']['BE']['Sawtooth']=ps

    #--------------------------------------------------------------------------
    # Other stuff
    #--------------------------------------------------------------------------

    meta['Param']['BE']['Comparison_With_NIR']=gu.ReadExcel(meta['Paths']['Model']['Parameters'] + '\\Parameters_Comparison_With_NIR.xlsx')
    meta['Param']['BE']['Reporting_Version_Comparison']=gu.ReadExcel(meta['Paths']['Model']['Parameters'] + '\\Parameters_Reporting_Version_Comparison.xlsx')
    meta['Param']['BE']['Forcing_Categories']=gu.ReadExcel(meta['Paths']['Model']['Parameters'] + '\\Parameters_ForcingCategories.xlsx')
    meta['Param']['BE']['Level_4_Categories_Status']=gu.ReadExcel(meta['Paths']['Model']['Parameters'] + '\\Parameters_Level_4_Categories_Status.xlsx')

    return meta

#%% Delete all output files

def DeleteAllOutputFiles(meta,pNam):

    for iScn in range(meta[pNam]['Project']['N Scenario']):
        files=glob.glob(meta['Paths'][pNam]['Output Scenario'][iScn] + '\\*')
        for f in files:
            os.remove(f)

    for iBat in range(meta[pNam]['Project']['N Batch']):
        try:
            pth=meta['Paths'][pNam]['Data'] + '\\Outputs\\WorkingOnBatch_' + FixFileNum(iBat) + '.pkl'
            os.remove(pth)
        except:
            pass

    return

#%% Import custom harvest assumptions (optional)
# *** Decomissioned ***

# def ImportCustomHarvestAssumptions(pthin):

#     # Import parameters
#     df=pd.read_excel(pthin,sheet_name='Sheet1',skiprows=1)

#     d={}
#     d['BiomassMerch_Affected']=df.iloc[1,1]
#     d['BiomassMerch_Burned']=df.iloc[3,1]
#     d['BiomassMerch_LeftOnSite']=df.iloc[4,1]
#     d['BiomassMerch_Removed']=df.iloc[5,1]

#     d['BiomassNonMerch_Affected']=df.iloc[1,2]
#     d['BiomassNonMerch_Burned']=df.iloc[3,2]
#     d['BiomassNonMerch_LeftOnSite']=df.iloc[4,2]
#     d['BiomassNonMerch_Removed']=df.iloc[5,2]

#     d['Snags_Affected']=df.iloc[1,3]
#     d['Snags_Burned']=df.iloc[3,3]
#     d['Snags_LeftOnSite']=df.iloc[4,3]
#     d['Snags_Removed']=df.iloc[5,3]

#     d['RemovedMerchToLumberMill']=df.iloc[7,1]
#     d['RemovedMerchToPulpMill']=df.iloc[8,1]
#     d['RemovedMerchToPelletMill']=df.iloc[9,1]
#     d['RemovedMerchToPlywoodMill']=df.iloc[10,1]
#     d['RemovedMerchToOSBMill']=df.iloc[11,1]
#     d['RemovedMerchToMDFMill']=df.iloc[12,1]
#     d['RemovedMerchToFirewood']=df.iloc[13,1]
#     d['RemovedMerchToIPP']=df.iloc[14,1]
#     d['RemovedMerchToLogExport']=df.iloc[15,1]

#     d['RemovedNonMerchToLumberMill']=df.iloc[7,2]
#     d['RemovedNonMerchToPulpMill']=df.iloc[8,2]
#     d['RemovedNonMerchToPelletMill']=df.iloc[9,2]
#     d['RemovedNonMerchToPlywoodMill']=df.iloc[10,2]
#     d['RemovedNonMerchToOSBMill']=df.iloc[11,2]
#     d['RemovedNonMerchToMDFMill']=df.iloc[12,2]
#     d['RemovedNonMerchToFirewood']=df.iloc[13,2]
#     d['RemovedNonMerchToIPP']=df.iloc[14,2]
#     d['RemovedNonMerchToLogExport']=df.iloc[15,2]

#     d['RemovedSnagStemToLumberMill']=df.iloc[7,3]
#     d['RemovedSnagStemToPulpMill']=df.iloc[8,3]
#     d['RemovedSnagStemToPelletMill']=df.iloc[9,3]
#     d['RemovedSnagStemToPlywoodMill']=df.iloc[10,3]
#     d['RemovedSnagStemToOSBMill']=df.iloc[11,3]
#     d['RemovedSnagStemToMDFMill']=df.iloc[12,3]
#     d['RemovedSnagStemToFirewood']=df.iloc[13,3]
#     d['RemovedSnagStemToIPP']=df.iloc[14,3]
#     d['RemovedSnagStemToLogExport']=df.iloc[15,3]

#     return d

#%% Unpack ensemble stats from MOS

def UnpackEnsembleStatsFromMos(meta,pNam,mos):

    tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)

    #--------------------------------------------------------------------------
    # Unpack contents for easy use
    #--------------------------------------------------------------------------

    mu_mos=[]
    for iScn in range(meta[pNam]['Project']['N Scenario']):
        d={}
        for k in mos[iScn]['v1']['Mean'].keys():
            d[k]=mos[iScn]['v1']['Mean'][k]['Ensemble Mean']
        for k in mos[iScn]['Cashflow']['Mean'].keys():
            d[k]=mos[iScn]['Cashflow']['Mean'][k]['Ensemble Mean']
        mu_mos.append(d)

    cil_mos=[]
    for iScn in range(meta[pNam]['Project']['N Scenario']):
        d={}
        for k in mos[iScn]['v1']['Mean'].keys():
            d[k]=mos[iScn]['v1']['Mean'][k]['Ensemble CIL']
        for k in mos[iScn]['Cashflow']['Mean'].keys():
            d[k]=mos[iScn]['Cashflow']['Mean'][k]['Ensemble CIL']
        cil_mos.append(d)

    cih_mos=[]
    for iScn in range(meta[pNam]['Project']['N Scenario']):
        d={}
        for k in mos[iScn]['v1']['Mean'].keys():
            d[k]=mos[iScn]['v1']['Mean'][k]['Ensemble CIH']
        for k in mos[iScn]['Cashflow']['Mean'].keys():
            d[k]=mos[iScn]['Cashflow']['Mean'][k]['Ensemble CIH']
        cih_mos.append(d)

    sdl_mos=[]
    for iScn in range(meta[pNam]['Project']['N Scenario']):
        d={}
        for k in mos[iScn]['v1']['Mean'].keys():
            d[k]=mos[iScn]['v1']['Mean'][k]['Ensemble Mean']-mos[iScn]['v1']['Mean'][k]['Ensemble SD']
        for k in mos[iScn]['Cashflow']['Mean'].keys():
            d[k]=mos[iScn]['Cashflow']['Mean'][k]['Ensemble Mean']-mos[iScn]['Cashflow']['Mean'][k]['Ensemble SD']
        sdl_mos.append(d)

    sdh_mos=[]
    for iScn in range(meta[pNam]['Project']['N Scenario']):
        d={}
        for k in mos[iScn]['v1']['Mean'].keys():
            d[k]=mos[iScn]['v1']['Mean'][k]['Ensemble Mean']+mos[iScn]['v1']['Mean'][k]['Ensemble SD']
        for k in mos[iScn]['Cashflow']['Mean'].keys():
            d[k]=mos[iScn]['Cashflow']['Mean'][k]['Ensemble Mean']+mos[iScn]['Cashflow']['Mean'][k]['Ensemble SD']
        sdh_mos.append(d)

#    p1_mos=[]
#    for iScn in range(meta[pNam]['Project']['N Scenario']):
#        d={}
#        for k in mos[iScn]['v1']['Sum'].keys():
#            d[k]=mos[iScn]['v1']['Sum'][k]['Ensemble P1']
#        for k in mos[iScn]['Cashflow']['Sum'].keys():
#            d[k]=mos[iScn]['Cashflow']['Sum'][k]['Ensemble P1']
#        p1_mos.append(d)
#
#    p10_mos=[]
#    for iScn in range(meta[pNam]['Project']['N Scenario']):
#        d={}
#        for k in mos[iScn]['v1']['Sum'].keys():
#            d[k]=mos[iScn]['v1']['Sum'][k]['Ensemble P10']
#        for k in mos[iScn]['Cashflow']['Sum'].keys():
#            d[k]=mos[iScn]['Cashflow']['Sum'][k]['Ensemble P10']
#        p10_mos.append(d)
#
#
#    p90_mos=[]
#    for iScn in range(meta[pNam]['Project']['N Scenario']):
#        d={}
#        for k in mos[iScn]['v1']['Sum'].keys():
#            d[k]=mos[iScn]['v1']['Sum'][k]['Ensemble P90']
#        for k in mos[iScn]['Cashflow']['Sum'].keys():
#            d[k]=mos[iScn]['Cashflow']['Sum'][k]['Ensemble P90']
#        p90_mos.append(d)
#
#    p99_mos=[]
#    for iScn in range(meta[pNam]['Project']['N Scenario']):
#        d={}
#        for k in mos[iScn]['v1']['Sum'].keys():
#            d[k]=mos[iScn]['v1']['Sum'][k]['Ensemble P99']
#        for k in mos[iScn]['Cashflow']['Sum'].keys():
#            d[k]=mos[iScn]['Cashflow']['Sum'][k]['Ensemble P99']
#        p99_mos.append(d)

    return tv,mu_mos,cil_mos,cih_mos,sdl_mos,sdh_mos

#%% Combine completed and future
# projects must have the same time period saved
def CombineProjectMOSs(meta,mos,pNamL):
    mosCF=copy.deepcopy(mos[pNamL[0]])
    for iScn in range(meta[pNamL[0]]['Project']['N Scenario']):
        for k1 in mosCF['Scenarios'][iScn]['Sum'].keys():
            for k2 in mosCF['Scenarios'][iScn]['Sum'][k1].keys():
                mosCF['Scenarios'][iScn]['Sum'][k1][k2]=mosCF['Scenarios'][iScn]['Sum'][k1][k2]*meta[pNamL[0]]['Project']['AEF']+mos[pNamL[1]]['Scenarios'][iScn]['Sum'][k1][k2]*meta[pNamL[1]]['Project']['AEF']
                mosCF['Scenarios'][iScn]['Mean'][k1][k2]=mosCF['Scenarios'][iScn]['Mean'][k1][k2]+mos[pNamL[1]]['Scenarios'][iScn]['Mean'][k1][k2]
    for cNam in mosCF['Delta'].keys():
        for k1 in mosCF['Delta'][cNam]['ByStrata']['Sum'].keys():
            for k2 in mosCF['Delta'][cNam]['ByStrata']['Sum'][k1].keys():
                mosCF['Delta'][cNam]['ByStrata']['Sum'][k1][k2]=mosCF['Delta'][cNam]['ByStrata']['Sum'][k1][k2]*meta[pNamL[0]]['Project']['AEF']+mos[pNamL[1]]['Delta'][cNam]['ByStrata']['Sum'][k1][k2]*meta[pNamL[1]]['Project']['AEF']
                mosCF['Delta'][cNam]['ByStrata']['Mean'][k1][k2]=mosCF['Delta'][cNam]['ByStrata']['Mean'][k1][k2]+mos[pNamL[1]]['Delta'][cNam]['ByStrata']['Mean'][k1][k2]
    return mosCF

