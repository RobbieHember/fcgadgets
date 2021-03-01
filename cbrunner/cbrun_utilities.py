
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import openpyxl
import gc as garc
import time
from fcgadgets.utilities import utilities_general as gu
from fcgadgets.utilities import utilities_inventory as invu
from fcgadgets.taz import aspatial_stat_models as asm

#%% CONVERT LUT NUMBER TO STRING NAME

def lut_n2s(dc,numb):
    if numb!=-999:
        vals=np.fromiter(dc.values(),dtype=float)
        keys=np.fromiter(dc.keys(),dtype='<U70')
        ind=np.where(vals==numb)[0]
        s=keys[ind]
    else:
        s=np.array(['Unidentified'],ndmin=1)
    return s

#%% Index to batch

def IndexToBatch(meta,iBat):
    iStart=meta['Batch Interval']*iBat
    iStop=np.minimum(meta['N Stand Full'],iStart+meta['Batch Interval'])
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
            Name.append('Salvage Logging')
            
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

def BuildEventChronologyFromSpreadsheet(meta):       
    
    for iScn in range(meta['N Scenario']):    
        for iEns in range(meta['N Ensemble']):        
            for iBat in range(meta['N Batch']):
    
                # Index to batch
                indBat=IndexToBatch(meta,iBat)
                N_StandsInBatch=len(indBat)
    
                tv=np.arange(meta['Year Start'],meta['Year End']+1,1)
    
                # Initialize dictionary
                ec={}
                ec['ID_Type']=np.zeros((meta['Year'].size,indBat.size,meta['Max Events Per Year']),dtype='int16')
                ec['MortalityFactor']=np.zeros((meta['Year'].size,indBat.size,meta['Max Events Per Year']),dtype='int16')
                ec['GrowthFactor']=np.zeros((meta['Year'].size,indBat.size,meta['Max Events Per Year']),dtype='int16')
                ec['ID_GrowthCurve']=np.zeros((meta['Year'].size,indBat.size,meta['Max Events Per Year']),dtype='int16')
            
                #--------------------------------------------------------------
                # Simulate wildfire occurrence and severity from Taz
                #--------------------------------------------------------------
                    
                if meta['Scenario'][iScn]['AAO Wildfire Status']=='On':
                        
                    # Import inventory to get BGC zone
                    inv=gu.ipickle(meta['Paths']['Input Scenario'][iScn] + '\\Inventory_Bat' + FixFileNum(iBat) + '.pkl')
                        
                    # Prepare required parameters dictionary
                    # Import fcgadgets parameters
                    par=invu.Load_Params(meta)
                    par['WF']['Scenario ID']=meta['Scenario'][iScn]['AAO Wildfire Scenario ID']
                    par['WF']['Exclude simulations during modern era']='On'
                    
                    # Think about moving this into the par dictionary
                    method_occ='DirectFromParetoDraw'
                    
                    # Normally, this would be run here, but this approach assumes
                    # that stands have been swapped with ensembles (when running
                    # from spreadsheet). Instead, it is being re-run for each stand
                    wf_sim=asm.GenerateWildfireEnsembleFromAAO(meta,par,inv['ID_BECZ'],method_occ)  
            
                for iS in range(N_StandsInBatch): 
                    
                    #----------------------------------------------------------
                    # Add spinup events
                    #----------------------------------------------------------
                    
                    ivl_spin=meta['Spinup Disturbance Return Inverval']                
                    YearRef=meta['Scenario'][iScn]['Year1_DisFromInv']
                    AgeRef=meta['Scenario'][iScn]['Age1_DisFromInv']                    
                    if AgeRef>=0:
                        Year=np.arange(YearRef-AgeRef-100*ivl_spin,YearRef-AgeRef+ivl_spin,ivl_spin)        
                    else:
                        Year1=meta['Year Start']+ivl_spin
                        Year2=meta['Spinup Year End']
                        Year=np.arange(Year1,Year2+1,meta['Spinup Disturbance Return Inverval'])
                    
                    for iYr in range(Year.size):
                        iT=np.where(tv==Year[iYr])[0]
                        ec['ID_Type'][iT,iS,0]=meta['LUT']['Dist'][meta['Spinup Disturbance Type']]
                        ec['MortalityFactor'][iT,iS,0]=100
                        ec['GrowthFactor'][iT,iS,0]=0
                        ec['ID_GrowthCurve'][iT,iS,0]=meta['Spinup Growth Curve ID']
                    
                    #----------------------------------------------------------
                    # Add simulated constant disturbances
                    #----------------------------------------------------------
               
                    # Historical disturbance from simulation 1
                    ri=meta['Scenario'][iScn]['ReturnInterval1_Hist_DisFromSim']
                    if (ri!=0) & (np.isnan(ri)==False):
                    
                        p_Dist=1/ri
                        p_Rand=np.random.uniform(0,1,size=(meta['Year'].size))        
                        it=np.where((p_Rand<p_Dist) & (meta['Year']>meta['Spinup Year End']) & (meta['Year']<meta['Year Project']))[0]
                        Year=meta['Year'][it]                    
                        for iYr in range(Year.size):
                            iT=np.where(tv==Year[iYr])[0]
                            ec['ID_Type'][iT,iS,0]=meta['LUT']['Dist'][meta['Scenario'][iScn]['Type1_Hist_DisFromSim']]
                            ec['MortalityFactor'][iT,iS,0]=100
                            ec['GrowthFactor'][iT,iS,0]=0
                            ec['ID_GrowthCurve'][iT,iS,0]=meta['Spinup Growth Curve ID']
                    
                    # Historical disturbance from simulation 2
                    ri=meta['Scenario'][iScn]['ReturnInterval2_Hist_DisFromSim']
                    if (ri!=0) & (np.isnan(ri)==False):
                    
                        p_Dist=1/ri
                        p_Rand=np.random.uniform(0,1,size=(meta['Year'].size))        
                        it=np.where((p_Rand<p_Dist) & (meta['Year']>meta['Spinup Year End']) & (meta['Year']<meta['Year Project']))[0]
                        Year=meta['Year'][it]                    
                        for iYr in range(Year.size):
                            iT=np.where(tv==Year[iYr])[0]
                            ec['ID_Type'][iT,iS,0]=meta['LUT']['Dist'][meta['Scenario'][iScn]['Type2_Hist_DisFromSim']]
                            ec['MortalityFactor'][iT,iS,0]=100
                            ec['GrowthFactor'][iT,iS,0]=0
                            ec['ID_GrowthCurve'][iT,iS,0]=meta['Spinup Growth Curve ID']
      
                    #----------------------------------------------------------
                    # Add events from inventory
                    #----------------------------------------------------------
                    
                    for iYr in range(1,7):
                        
                        if np.isnan(meta['Scenario'][iScn]['Year' + str(iYr) + '_DisFromInv'])==True:
                            continue
            
                        # If IDW, convert IDW class to growth and mortality factor
                        sc=np.array(['IDW-T','IDW-L','IDW-M','IDW-S','IDW-V','IDW-MM','IDW-MS','IDW-MV','IDW-SS','IDW-SV','IDW-VV'])
                        flg_i=0
                        indSc=np.where(sc==meta['Scenario'][iScn]['Type' + str(iYr) + '_DisFromInv'])[0]
                        if indSc.size!=0:
                            if flg_i==0:
                                dfParDistBySC=pd.read_excel(meta['Paths']['Model Code'] + '\\Parameters\\Parameters_DisturbanceBySeverityClass.xlsx')
                                flg_i=1
                            indPar=np.where( (dfParDistBySC['Name']=='IDW') & (dfParDistBySC['SeverityCD']==sc[indSc[0]][4:]) )[0]
                            ID_TypeN=meta['LUT']['Dist']['IDW']
                            MF=dfParDistBySC.loc[indPar,'MortalityFactor']
                            GF=dfParDistBySC.loc[indPar,'GrowthFactor']
                        else:
                            ID_TypeS=meta['Scenario'][iScn]['Type' + str(iYr) + '_DisFromInv']
                            ID_TypeN=meta['LUT']['Dist'][ID_TypeS]
                            MF=meta['Scenario'][iScn]['Severity' + str(iYr) + '_DisFromInv']
                            GF=0
            
                        Year=meta['Scenario'][iScn]['Year' + str(iYr) + '_DisFromInv']
                        iT=np.where(tv==Year)[0]
                        
                        iE=np.where(ec['ID_Type'][iT,iS,:]==0)[1]

                        ec['ID_Type'][iT,iS,iE[0]]=ID_TypeN
                        ec['MortalityFactor'][iT,iS,iE[0]]=MF
                        ec['GrowthFactor'][iT,iS,iE[0]]=GF
                        ec['ID_GrowthCurve'][iT,iS,iE[0]]=meta['Scenario'][iScn]['GrowthCurve' + str(iYr) + '_DisFromInv']

                    #----------------------------------------------------------
                    # Add simulated constant future disturbances
                    #----------------------------------------------------------
               
                    # Future disturbance from simulation 1
                    ri=meta['Scenario'][iScn]['ReturnInterval1_Fut_DisFromSim']
                    if (ri!=0) & (np.isnan(ri)==False):
                    
                        p_Dist=1/ri
                        p_Rand=np.random.uniform(0,1,size=(meta['Year'].size))        
                        it=np.where((p_Rand<p_Dist) & (meta['Year']>meta['Year Project']))[0]
                        Year=meta['Year'][it]                    
                        for iYr in range(Year.size):
                            iT=np.where(tv==Year[iYr])[0]
                            ec['ID_Type'][iT,iS,0]=meta['LUT']['Dist'][meta['Scenario'][iScn]['Type1_Fut_DisFromSim']]
                            ec['MortalityFactor'][iT,iS,0]=100
                            ec['GrowthFactor'][iT,iS,0]=0
                            ec['ID_GrowthCurve'][iT,iS,0]=meta['Spinup Growth Curve ID']
                    
                    # Future disturbance from simulation 2
                    ri=meta['Scenario'][iScn]['ReturnInterval2_Fut_DisFromSim']
                    if (ri!=0) & (np.isnan(ri)==False):
                    
                        p_Dist=1/ri
                        p_Rand=np.random.uniform(0,1,size=(meta['Year'].size))        
                        it=np.where((p_Rand<p_Dist) & (meta['Year']>meta['Year Project']))[0]
                        Year=meta['Year'][it]                    
                        for iYr in range(Year.size):
                            iT=np.where(tv==Year[iYr])[0]
                            ec['ID_Type'][iT,iS,0]=meta['LUT']['Dist'][meta['Scenario'][iScn]['Type2_Fut_DisFromSim']]
                            ec['MortalityFactor'][iT,iS,0]=100
                            ec['GrowthFactor'][iT,iS,0]=0
                            ec['ID_GrowthCurve'][iT,iS,0]=meta['Spinup Growth Curve ID']    
            
                    #----------------------------------------------------------
                    # Add simulated wildfire from Taz
                    #----------------------------------------------------------
                    
                    if meta['Scenario'][iScn]['AAO Wildfire Status']=='On':
            
                        ind=np.where(wf_sim['Occurrence'][:,iS]==1)[0]
                        if ind.size==0:
                            continue
                        
                        ID_Type=meta['LUT']['Dist']['Wildfire']*np.ones(ind.size)
                        Year=tv[ind]
                        MortF=wf_sim['Mortality'][ind,iS]
                        GrowthF=0*np.ones(ind.size)
                        ID_GrowthCurve=1*np.ones(ind.size)
                            
                        for iYr in range(Year.size):
                            iT=np.where(tv==Year[iYr])[0]
                            ec['ID_Type'][iT,iS,0]=ID_Type[iYr]
                            ec['MortalityFactor'][iT,iS,0]=MortF[iYr]
                            ec['GrowthFactor'][iT,iS,0]=GrowthF[iYr]
                            ec['ID_GrowthCurve'][iT,iS,0]=ID_GrowthCurve[iYr]
    
                    #----------------------------------------------------------
                    # Add simulated MPB from Taz
                    #----------------------------------------------------------

            
                #--------------------------------------------------------------
                # Save to file            
                #--------------------------------------------------------------
                
                gu.opickle(meta['Paths']['Input Scenario'][iScn] + '\\Events_Ens' + FixFileNum(iEns) + '_Bat' + FixFileNum(iBat) + '.pkl',ec)

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

##%% Import look-up tables
#
#def ImportLUTs(pthin):
#    
#    # Open connection to parameter database    
#    par=gu.ipickle(pthin + '\\Parameters\\Parameters.pkl')
#    
#    # Import distubance type        
#    LUT_Dist={}
#    for i in range(len(par['Disturbances']['Name'])):
#        LUT_Dist[par['Disturbances']['Name'][i]]=par['Disturbances']['ID'][i]
#    
#    # BGC zone     
#    LUT_BGC_Zone={}
#    for i in range(len(par['BGC_ZONE']['CODE_BGC_ZONE'])):
#        LUT_BGC_Zone[par['BGC_ZONE']['CODE_BGC_ZONE'][i]]=par['BGC_ZONE']['ID_BGC_ZONE'][i]
#    
#    # Species
#    LUT_Spc={}
#    for i in range(len(par['SRS']['SRS_CD'])):
#        LUT_Spc[par['SRS']['SRS_CD'][i]]=par['SRS']['SRS_ID'][i]
#    
#    return LUT_Dist,LUT_Spc,LUT_BGC_Zone

#%% Configure project

def ImportProjectConfig(meta):
    
    #--------------------------------------------------------------------------
    # Import project parameters from spreadsheet
    #--------------------------------------------------------------------------
    
    df_p=pd.read_excel(meta['Paths']['Project'] + '\\Inputs\\ProjectConfig.xlsx',sheet_name='Inputs',skiprows=1,usecols='A:B')
    for i in range(df_p.shape[0]): 
        meta.update({df_p.iloc[i,0]:df_p.iat[i,1]})    
    
    #--------------------------------------------------------------------------
    # Import look-up tables
    #--------------------------------------------------------------------------
    
    meta=invu.Load_LUTs(meta)
    
#    meta['LUT']={}
#    
#    LUT_Dist,LUT_Spc,LUT_BGC_Zone=ImportLUTs(meta['Paths']['Model Code'])
#    
#    meta['LUT']['Dist']=LUT_Dist
#    meta['LUT']['Spc']=LUT_Spc
#    meta['LUT']['BGC Zone']=LUT_BGC_Zone

    #--------------------------------------------------------------------------
    # Define pool names
    #--------------------------------------------------------------------------
    
    # Pool names (ecosystem)
    # *** If you change this, you need to change the same list in "Update Parameters" ***
    meta['Name Pools Eco']=['StemMerch','StemNonMerch','Foliage','Branch','Bark','RootCoarse','RootFine', \
        'FelledStemMerch','FelledStemNonMerch','FelledBranch','FelledBark','FelledSnagStem','FelledSnagBranch', \
        'LitterVF','LitterF','LitterM','LitterS','SnagStem','SnagBranch','SoilVF','SoilF','SoilS']
    
    # Used in fert work: ,'LitterDecomp'
    
    # Number of ecosystem pools
    meta['N Pools Eco']=len(meta['Name Pools Eco'])
    
    # Pool names (products)
    meta['Name Pools Pro']=['SFH','MFH','Comm','Furn','Ship','Repairs','Other','Paper','Fuel','Firewood','EffluentPulp', \
        'DumpWood','DumpPaper','LandfillWoodDegradable','LandfillWoodNonDegradable', \
        'LandfillPaperDegradable','LandfillPaperNonDegradable','E_CO2','E_CH4','Cants']    
    
    # Number of product pools
    meta['N Pools Pro']=len(meta['Name Pools Pro'])
    
    #--------------------------------------------------------------------------
    # Define indices to each pool
    #--------------------------------------------------------------------------
    
    # Indices to ecosystem pools pools
    meta['iEP']={}; cnt=0
    for nam in meta['Name Pools Eco']:
        meta['iEP'][nam]=cnt
        cnt=cnt+1
    iEP=meta['iEP']
    meta['iEP']['BiomassTotal']=np.array([iEP['StemMerch'],iEP['StemNonMerch'],iEP['Foliage'],iEP['Branch'],iEP['Bark'],iEP['RootCoarse'],iEP['RootFine']])
    meta['iEP']['BiomassAboveground']=np.array([iEP['StemMerch'],iEP['StemNonMerch'],iEP['Foliage'],iEP['Branch'],iEP['Bark']])
    meta['iEP']['BiomassBelowground']=np.array([iEP['RootCoarse'],iEP['RootFine']])
    meta['iEP']['DeadWood']=np.array([iEP['FelledStemMerch'],iEP['FelledStemNonMerch'],iEP['FelledBranch'],iEP['FelledBark'],iEP['SnagStem'],iEP['SnagBranch']])
    meta['iEP']['Litter']=np.array([iEP['LitterVF'],iEP['LitterF'],iEP['LitterM'],iEP['LitterS']])
    meta['iEP']['Felled']=np.array([iEP['FelledStemMerch'],iEP['FelledStemNonMerch'],iEP['FelledBranch'],iEP['FelledBark'],iEP['FelledSnagStem'],iEP['FelledSnagBranch']])
    meta['iEP']['Soil']=np.array([iEP['SoilVF'],iEP['SoilF'],iEP['SoilS']])
    
    # Indices to produce pools pools
    meta['iPP']={}; cnt=0
    for nam in meta['Name Pools Pro']:
        meta['iPP'][nam]=cnt
        cnt=cnt+1
    
    #--------------------------------------------------------------------------
    # Define time
    #--------------------------------------------------------------------------
    
    # Calendar year
    meta['Year']=np.arange(meta['Year Start'],meta['Year End']+1,1)
    meta['N Time']=meta['Year'].size
    
    #--------------------------------------------------------------------------
    # Dimensions of simulation
    #--------------------------------------------------------------------------
    
    # Number of stands 
    if meta['Scenario Source']=='Spreadsheet' and meta['N Ensemble']==1:
        meta['N Stand']=1    
    elif meta['Scenario Source']=='Spreadsheet' and meta['N Ensemble']>1:
        # If running ensembles from spreadsheet, it is faster to run them as stands
        meta['N Stand']=meta['N Ensemble']
        meta['N Ensemble']=1
    
    # Number of batches
    meta['N Batch']=np.ceil(meta['N Stand']/meta['Batch Interval']).astype(int)

    # Projects with more than one batch change N Stand - create a full version
    meta['N Stand Full']=meta['N Stand']
    
    #--------------------------------------------------------------------------
    # Define scenario parameters
    #--------------------------------------------------------------------------
    
    if meta['Scenario Source']=='Spreadsheet':
        
        # Import from spreadsheet
        df_s=pd.read_excel(meta['Paths']['Project'] + '\\Inputs\\ProjectConfig.xlsx',sheet_name='Inputs',skiprows=1,usecols='D:AN')
        df_s=df_s.iloc[:,df_s.iloc[0,:].isnull().values==False] 
        meta['Scenario']=list()
        for i in range(1,df_s.shape[1]):
            pScn0={}
            for j in range(df_s.shape[0]):    
                pScn0.update({df_s.iloc[j,0]:df_s.iat[j,i]})
            meta['Scenario'].append(pScn0)

        # Number of scenarios
        meta['N Scenario']=np.sum([i['Scenario_OnSwitch']=='On' for i in meta['Scenario']])
        
    elif (meta['Scenario Source']=='Script') & (meta['Biomass Module']=='Sawtooth'):
        
        # Initialize scenario data structure
        meta['Scenario']=list()
        for i in range(meta['N Scenario']):
            pScn0={'SRS1_CD':' '}
            meta['Scenario'].append(pScn0)
            
    elif (meta['Scenario Source']=='Script') & (meta['Biomass Module']=='BatchTIPSY'):
        
        # Initialize scenario data structure
        meta['Scenario']=list()
        for i in range(meta['N Scenario']):
            pScn0={}
            meta['Scenario'].append(pScn0)

    #--------------------------------------------------------------------------
    # Initialize project folders if they do not exist
    #--------------------------------------------------------------------------
    
    meta['Paths']['Input Scenario']=[]
    meta['Paths']['Output Scenario']=[]
    for iScn in range(0,meta['N Scenario']):    
        meta['Paths']['Input Scenario'].append(meta['Paths']['Project'] + '\\Inputs\\Scenario' + FixFileNum(iScn))
        if os.path.exists(meta['Paths']['Input Scenario'][iScn])==False:
            os.mkdir(meta['Paths']['Input Scenario'][iScn])
        meta['Paths']['Output Scenario'].append(meta['Paths']['Project'] + '\\Outputs\\Scenario' + FixFileNum(iScn))
        if os.path.exists(meta['Paths']['Output Scenario'][iScn])==False:
            os.mkdir(meta['Paths']['Output Scenario'][iScn])
    
    #--------------------------------------------------------------------------
    # Scale factors
    #--------------------------------------------------------------------------
    
    # *** Scale factor for saving results (this needs to be 100, 10 does not 
    # capture carbon fluxes and it will affect GHG benefit estimates) ***
    meta['Scale Factor Export']=0.01
    
    #--------------------------------------------------------------------------
    # Disturbance information
    #--------------------------------------------------------------------------
    
    # Initialize flag indicating whether a disturbance occurred
    meta['Flag_Disturbed']=np.zeros(meta['N Stand'])
    
    #--------------------------------------------------------------------------
    # Growth curve information
    #--------------------------------------------------------------------------
    
    meta['GC']={}
    meta['GC']['N Growth Curves']=5
    meta['GC']['ID GC']=np.array([1,2,3,4,5])
    meta['GC']['BatchTIPSY Maximum Age']=200
    meta['GC']['BatchTIPSY Column Names']=['Age','VolTot0','VolMerch125',
        'VolMerch175','ODT_Bark','ODT_Branch','ODT_Foliage','ODT_Roots',
        'ODT_Stem','MortalityVolumeTotal']
    
    # Scale factor for growth curves
    # Note: Do not change this to 0.1 - aerial fertilization response will not work properly at 0.1
    meta['GC']['Scale Factor']=0.01
    
    #--------------------------------------------------------------------------
    # Initialize scenario switches
    #--------------------------------------------------------------------------
    
    meta['Scenario Switch']={}
    
    #--------------------------------------------------------------------------
    # Growth factor information
    # *** Not currently used ***
    #--------------------------------------------------------------------------
    
    # Default status of growth factors
    meta['Scenario Switch']['Net Growth Factor Status']=[None]*meta['N Scenario']
    meta['Scenario Switch']['Mortality Factor Status']=[None]*meta['N Scenario']
    for iScn in range(0,meta['N Scenario']): 
        meta['Scenario Switch']['Net Growth Factor Status'][iScn]='Off'
        meta['Scenario Switch']['Mortality Factor Status'][iScn]='Off'
        #meta['Scenario Switch'][iScn]['Status Net Growth Factor']='Off'
        #meta['Scenario'][iScn]['Status Mortality Factor']='Off'
    
    #--------------------------------------------------------------------------
    # Harvested wood product information
    #--------------------------------------------------------------------------
    
    # Year to start calling annual HWP methods
    meta['HWP Year Start']=1800
    
    #--------------------------------------------------------------------------
    # Nutrient application information (for compatibility with "actions" module)
    #--------------------------------------------------------------------------
    
    meta['NM']={}
    
    # Initialize index to stands affected by nutrient application
    meta['NM']['iApplication']=[]
    
    # Nutrient addition response yearly counter
    meta['NM']['ResponseCounter']=np.zeros(meta['N Stand'])
    
    return meta

#%% LOAD SINGLE OUTPUT FILE FOR SCENARIO A, ENSEMBLE, B AND BATCH C

def LoadSingleOutputFile(meta,iScn,iEns,iBat):
       
    # Open output data
    data=gu.ipickle(meta['Paths']['Project'] + '\\Outputs\\Scenario' + FixFileNum(iScn) + '\\Data_Scn' + FixFileNum(iScn) + '_Ens' + FixFileNum(iEns) + '_Bat' + FixFileNum(iBat) + '.pkl')
                
    # Convert to float and apply scale factor
    for k in data.keys():
        
        # Skip mortality summary by agent
        if k=='C_M_ByAgent':
            continue
        
        data[k]=data[k].astype(float)
        data[k]=data[k]*meta['Scale Factor Export']
    
    # Mortality summary by agent
    for k in data['C_M_ByAgent'].keys():
        data['C_M_ByAgent'][k]=data['C_M_ByAgent'][k].astype(float)
        #data['C_M_ByAgent'][k]=data['C_M_ByAgent'][k]*meta['Scale Factor Export']
    
    v=BunchDictionary(data)
        
    return v

#%% LOAD SCENARIO RUSULTS
# Return a list of dictionaries for each scenario. If multiple ensemble were run, 
# the function will retun the average.

def LoadScenarioResults(meta,scn):
    
    # Initialize list that will contain scenarios
    v=[]
    for iScn in scn:
        for iEns in range(0,meta['N Ensemble']):            
            for iBat in range(0,meta['N Batch']):
                
                # Open results 
                pth=meta['Paths']['Project'] + '\\Outputs\\Scenario' + FixFileNum(iScn) + '\\Data_Scn' + FixFileNum(iScn) + '_Ens' + FixFileNum(iEns) + '_Bat' + FixFileNum(iBat) + '.pkl'
                data_bat=gu.ipickle(pth)
                
                # Convert to float and apply scale factor
                for k in data_bat.keys():
                    
                    # Skip mortality summary by agent
                    if k=='C_M_ByAgent':
                        continue                      
                    
                    data_bat[k]=data_bat[k].astype(float)
                    data_bat[k]=data_bat[k]*meta['Scale Factor Export']
                
                # Accumulate data in each batch
                key=list(data_bat.keys())
                if iBat==0:
                    data_all=data_bat
                else:                    
                    for key in data_bat.keys():
                        
                        if key=='Year':
                            continue
                        
                        # Skip mortality summary by agent
                        if key=='C_M_ByAgent':
                            continue  
                    
                        data_all[key]=np.append(data_all[key],data_bat[key],axis=1)
            
            # Sum across ensembles
            if iEns==0:
                data_sum2ave=data_all
            else:                    
                for key in data_bat.keys():
                    
                    # Skip mortality summary by agent
                    if key=='C_M_ByAgent':
                        continue 
                    data_sum2ave[key]=data_sum2ave[key]+data_all[key]
        
        # If the simulation includes ensembles, calculate average
        for key in data_bat.keys():
            
            # Skip mortality summary by agent
            if key=='C_M_ByAgent':
                continue 
            
            data_sum2ave[key]=data_sum2ave[key]/meta['N Ensemble']        
        
        v.append(BunchDictionary(data_sum2ave))
        
    return v

#%% Calculate GHG balance

def CalculateGHGBalance(v1,meta):

    if type(v1)!=list:
        v1=[v1]
    
    # Global warming potential of CO2
    gwp_co2=1
        
    # Global warming potential for CH4 -> using IPCC 2007 (AR4) values to be 
    # consistent with Greenhouse Gas Reduction Targets Act Carbon Neutral Government Regulation (B.C. Reg. 193/2014)
    #gwp_ch4=28
    gwp_ch4=meta['psl']['bGWP_CH4_AR5']
    
    # CO not recognized as GHG in BC Reg 193/2014, so using most recent
    #gwp_co=3.3
    gwp_co=meta['psl']['bGWP_CO_AR5']
    
    # Global warming potential for CH4 -> using IPCC 2007 (AR4) values to be 
    # consistent with Greenhouse Gas Reduction Targets Act Carbon Neutral Government Regulation (B.C. Reg. 193/2014)
    #gwp_n2o=298
    gwp_n2o=meta['psl']['bGWP_N2O_AR5']
       
    v2=[]
    for i in range(len(v1)):
        Year=np.arange(meta['Year Start Saving'],meta['Year End']+1,1)
        A=v1[i].A.copy()
        
        if meta['Save Biomass Pools']=='On':
            Eco_G_Gross=meta['psl']['bRatio_CO2_to_C']*np.sum(v1[i].C_G_Gross,axis=2)
            Eco_G_Net=meta['psl']['bRatio_CO2_to_C']*np.sum(v1[i].C_G_Net,axis=2)
            Eco_M_Reg=meta['psl']['bRatio_CO2_to_C']*np.sum(v1[i].C_M_Reg,axis=2)
            Eco_NPP=meta['psl']['bRatio_CO2_to_C']*np.sum(v1[i].C_NPP,axis=2)
            Eco_RH=meta['psl']['bRatio_CO2_to_C']*np.sum(v1[i].C_RH,axis=2)
            Eco_LF=meta['psl']['bRatio_CO2_to_C']*np.sum(v1[i].C_LF,axis=2)
        elif meta['Save Biomass Pools']!='On': 
            Eco_G_Gross=meta['psl']['bRatio_CO2_to_C']*v1[i].C_G_Gross.copy()
            Eco_G_Net=meta['psl']['bRatio_CO2_to_C']*v1[i].C_G_Net.copy()
            Eco_M_Reg=meta['psl']['bRatio_CO2_to_C']*v1[i].C_M_Reg.copy()
            Eco_NPP=meta['psl']['bRatio_CO2_to_C']*v1[i].C_NPP.copy()
            Eco_RH=meta['psl']['bRatio_CO2_to_C']*v1[i].C_RH.copy()
            Eco_LF=meta['psl']['bRatio_CO2_to_C']*v1[i].C_LF.copy()
    
        # Carbon dioxide flux (tCO2/ha/yr)
        Eco_E_Fire=v1[i].CO2e_E_Fire.copy()
        
        Eco_E_Operations=meta['psl']['bRatio_CO2_to_C']*v1[i].C_E_Operations.copy()
        
        Eco_Removals=meta['psl']['bRatio_CO2_to_C']* \
                        (v1[i].C_RemovedMerch.copy()+ \
                         v1[i].C_RemovedNonMerch.copy()+ \
                         v1[i].C_RemovedSnagStem.copy())
        
        Eco_E_OpenBurning=np.zeros(Eco_E_Fire.shape)
        ind=np.where(Eco_Removals>0)
        if ind[0].size>0:
            Eco_E_OpenBurning[ind[0],ind[1]]=Eco_E_Fire[ind[0],ind[1]]
        
        Eco_E_Wildfire=np.zeros(Eco_E_Fire.shape)
        ind=np.where(Eco_Removals==0)
        if ind[0].size>0:
            Eco_E_Wildfire[ind[0],ind[1]]=Eco_E_Fire[ind[0],ind[1]]
        
        Eco_NGHGB=Eco_NPP-Eco_RH-Eco_E_Fire-Eco_E_Operations-Eco_Removals
    
        if meta['Save Biomass Pools']=='On':
            Eco_Biomass=meta['psl']['bRatio_CO2_to_C']*np.sum(v1[i].C_Eco_Pools[:,:,meta['iEP']['BiomassTotal']].copy(),axis=2)       
            Eco_BiomassAG=meta['psl']['bRatio_CO2_to_C']*np.nansum(v1[i].C_Eco_Pools[:,:,meta['iEP']['BiomassAboveground']].copy(),axis=2)
            Eco_Felled=meta['psl']['bRatio_CO2_to_C']*np.sum(v1[i].C_Eco_Pools[:,:,meta['iEP']['Felled']].copy(),axis=2)
            Eco_Litter=meta['psl']['bRatio_CO2_to_C']*np.sum(v1[i].C_Eco_Pools[:,:,meta['iEP']['Litter']].copy(),axis=2)
            Eco_DeadWood=meta['psl']['bRatio_CO2_to_C']*np.sum(v1[i].C_Eco_Pools[:,:,meta['iEP']['DeadWood']].copy(),axis=2)
            Eco_Soil=meta['psl']['bRatio_CO2_to_C']*np.sum(v1[i].C_Eco_Pools[:,:,meta['iEP']['Soil']].copy(),axis=2)
            Pro_InUse=meta['psl']['bRatio_CO2_to_C']*np.sum(v1[i].C_Pro_Pools[:,:,0:10].copy(),axis=2)
            Pro_DumpLandfill=meta['psl']['bRatio_CO2_to_C']*np.sum(v1[i].C_Pro_Pools[:,:,10:17].copy(),axis=2)
            Pro_Emissions=gwp_co2*meta['psl']['bRatio_CO2_to_C']*v1[i].C_Pro_Pools[:,:,17].copy()+gwp_ch4*meta['psl']['bRatio_CO2_to_C']*v1[i].C_Pro_Pools[:,:,18].copy()                 
        else: 
            Eco_Biomass=meta['psl']['bRatio_CO2_to_C']*v1[i].C_Biomass.copy()
            Eco_BiomassAG=meta['psl']['bRatio_CO2_to_C']*v1[i].C_BiomassAG.copy()
            Eco_Felled=meta['psl']['bRatio_CO2_to_C']*v1[i].C_Felled.copy()
            Eco_Litter=meta['psl']['bRatio_CO2_to_C']*v1[i].C_Litter.copy()
            Eco_DeadWood=meta['psl']['bRatio_CO2_to_C']*v1[i].C_DeadWood.copy()
            Eco_Soil=meta['psl']['bRatio_CO2_to_C']*v1[i].C_Soil.copy()
            Pro_InUse=meta['psl']['bRatio_CO2_to_C']*v1[i].C_InUse.copy()
            Pro_DumpLandfill=meta['psl']['bRatio_CO2_to_C']*v1[i].C_DumpLandfill.copy()
            Pro_Emissions=v1[i].CO2e_E_Products.copy() # Already converted to CO2e
        
        Eco_Total=Eco_Biomass+Eco_Felled+Eco_Litter+Eco_DeadWood+Eco_Soil
        Pro_Total=Pro_InUse+Pro_DumpLandfill        
        Sec_NGHGB=Eco_NPP-Eco_RH-Eco_E_Fire-Eco_E_Operations-Pro_Emissions
    
        # Add data to dictionary        
        d={'Year':Year,
              'A':A,
              'Eco_G_Gross':Eco_G_Gross,
              'Eco_G_Net':Eco_G_Net,
              'Eco_M_Reg':Eco_M_Reg,
              'Eco_NPP':Eco_NPP,
              'Eco_RH':Eco_RH,
              'Eco_LF':Eco_LF,
              'Eco_E_Wildfire':Eco_E_Wildfire,
              'Eco_E_OpenBurning':Eco_E_OpenBurning,              
              'Eco_E_Operations':Eco_E_Operations,
              'Eco_Removals':Eco_Removals,
              'Eco_NGHGB':Eco_NGHGB,
              'Eco_Biomass':Eco_Biomass,
              'Eco_Felled':Eco_Felled,
              'Eco_Litter':Eco_Litter,
              'Eco_DeadWood':Eco_DeadWood,
              'Eco_Soil':Eco_Soil,
              'Eco_Total':Eco_Total,
              'Pro_InUse':Pro_InUse,
              'Pro_DumpLandfill':Pro_DumpLandfill,
              'Pro_Total':Pro_Total,
              'Pro_Emissions':Pro_Emissions,
              'Sec_NGHGB':Sec_NGHGB}

        v2.append(BunchDictionary(d))        
        
        #----------------------------------------------------------------------
        # Keep labels for each variable for plotting
        #----------------------------------------------------------------------
        
        HandleLabels=['Time, years','Stand age (years)','Gross growth (tCO2e/ha/yr)',
             'Net growth (tCO2e/ha/yr)','Regular mortality (tCO2e/ha/yr)',
             'NPP (tCO2e/ha/yr)','RH (tCO2e/ha/yr)','LF (tCO2e/ha/yr)',
             'Wildfire emissions (tCO2e/ha/yr)',
             'Open burning emissions (tCO2e/ha/yr)','Operational emissions (tCO2e/ha/yr)',
             'Removals (tCO2e/ha/yr)','Net ecosystem GHG balance (tCO2e/ha/yr)',
             'Biomass (tCO2e/ha)','Felled (tCO2e/ha)','Litter (tCO2e/ha)','Dead wood (tCO2e/ha)',
             'Soil (tCO2e/ha)','Total ecosystem (tCO2e/ha)','In-use products (tCO2e/ha)',
             'Dump and landfill (tCO2e/ha)','Total product sector (tCO2e/ha)',
             'Product emissions (tCO2e/ha/yr)','Net sector GHG balance (tCO2e/ha/yr)']        
        Keys=np.array(list(d.keys()))
        d1={}
        for j in range(len(HandleLabels)):
            d1[Keys[j]]=HandleLabels[j]  
        meta['Labels GHG Balance']=d1
    
    # Not sure why, but the list is being nested in a second list
    #v2=v2[0]
    
    return v2,meta

#%% Post process BatchTIPSY output

def PostProcessBatchTIPSY(meta):

    # Function used to smooth curves
    def smooth(y, box_pts):
        box = np.ones(box_pts)/box_pts
        y_smooth = np.convolve(y, box, mode='same')
        return y_smooth
    
    # TIPSY exports curves as MgDM/ha/yr, CBRunner expects inputs of MgC/ha/yr - create conversion factor
    dm2c=0.5

    # Growth curve parameters and TIPSY outputs
    dfPar=pd.read_excel(meta['Paths']['Project'] + '\\Inputs\\GrowthCurvesTIPSY_Parameters.xlsx',sheet_name='Sheet1',skiprows=6)
    txtDat=np.loadtxt(meta['Paths']['Project'] + '\\Inputs\\GrowthCurvesTIPSY_Output.out',skiprows=4)

    # TIPSY saves to text file -> convert to dataframe (column names must match TIPSY output file design)
    dfDat=pd.DataFrame(txtDat,columns=meta['GC']['BatchTIPSY Column Names'])

    # Define age vector (must be consistent with how TIPSY was set up)
    Age=np.arange(0,meta['GC']['BatchTIPSY Maximum Age']+1,1)

    # Get dimensions of the TIPSY output file to reshape the data into Age x Stand
    N_Age=Age.size
    N_GC=int(dfDat.shape[0]/N_Age)

    # Stemwood

    # Merchantable stemwood volume
    V_StemMerch=np.reshape(dfDat['VolMerch125'].values,(N_Age,N_GC),order='F')
    V_StemTot=np.reshape(dfDat['VolTot0'].values,(N_Age,N_GC),order='F')
    G_VStemMerch=np.append(np.zeros((1,N_GC)),np.diff(V_StemMerch,axis=0),axis=0)

    # Extract age responses for each biomass pool
    C_Stem=dm2c*np.reshape(dfDat['ODT_Stem'].values,(N_Age,N_GC),order='F')

    #import matplotlib.pyplot as plt
    #plt.close('all')
    #plt.plot(C_Stem[:,0])

    # Apply smoothing - it messes up the last ten years so don't smooth that part
    for j in range(C_Stem.shape[1]):
        a=smooth(C_Stem[:,j],10)
        C_Stem[:-10,j]=a[:-10]
        a=smooth(V_StemMerch[:,j],10)
        V_StemMerch[:-10,j]=a[:-10]
        a=smooth(V_StemTot[:,j],10)
        V_StemTot[:-10,j]=a[:-10]
    #plt.plot(C_Stem[:,0],'--')
    
    # Define the fraction of merchantable stemwood
    fMerch=np.nan_to_num(V_StemMerch/V_StemTot)
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
    
    # Add negative nonmerch to merch 
    ind=np.where(G_StemNonMerch<0)
    G_StemMerch[ind]=G_StemMerch[ind]+G_StemNonMerch[ind]
    G_StemNonMerch[ind]=0    
    
    #import matplotlib.pyplot as plt
    #plt.close('all')
    #plt.plot(np.maximum(-1,G_Stem[:,0]))
    #plt.plot(np.maximum(-1,G_StemMerch[:,0]),'--')
    #plt.plot(np.maximum(-1,G_StemNonMerch[:,0]),'-.')

    # Other pools   

    # Foliage biomass is very low, revise
    bF1=0.579
    bF2=0.602
    C_Foliage=np.maximum(0,bF1*C_Stem**bF2)   
    #C_Foliage=dm2c*np.reshape(dfDat['ODT_Foliage'].values,(N_Age,N_GC),order='F')
    
    C_Branch=dm2c*np.reshape(dfDat['ODT_Branch'].values,(N_Age,N_GC),order='F')
    C_Bark=dm2c*np.reshape(dfDat['ODT_Bark'].values,(N_Age,N_GC),order='F')

    G_Foliage=np.append(z,np.diff(C_Foliage,axis=0),axis=0)
    G_Branch=np.append(z,np.diff(C_Branch,axis=0),axis=0)
    G_Bark=np.append(z,np.diff(C_Bark,axis=0),axis=0)

    for iScn in range(meta['N Scenario']):
    
        # Define growth curves (from TIPSY output)
        for iBat in range(meta['N Batch']):
                
            # Index to batch
            indBat=IndexToBatch(meta,iBat)
                
            # Import event chronology
            #ec=gu.ipickle(meta['Paths']['Input Scenario'][iScn] + '\\Events_Ens' + FixFileNum(0) + '_Bat' + FixFileNum(iBat) + '.pkl')
            
            for iGC in range(3):
                
                # Initialize age response of net growth
                G=np.zeros((N_Age,indBat.size,6),dtype=np.int32)
    
                # Populate the growth curve
                for iS in range(indBat.size):
               
                    #u=np.unique(ec['ID_GrowthCurve'][:,iS,:])
                    
                    if meta['Scenario Source']=='Spreadsheet':
                        
                        indTIPSY=np.where(
                                (dfPar['ID_Scenario']==iScn+1) &
                                (dfPar['ID_GC']==int(meta['GC']['ID GC'][iGC])) )[0]                    
                    
                    elif meta['Scenario Source']=='Script':                        
                        
                        indTIPSY=np.where(
                            (dfPar['ID_Stand']==indBat[iS]+1) & 
                            (dfPar['ID_Scenario']==iScn+1) &
                            (dfPar['ID_GC']==int(iGC+1)))[0]
                        # dh[iS].ID_GrowthCurve[indDH][0]
                        
                    
                    if indTIPSY.size==0:
                        # This can happen if only some stands have a third GC, for example
                        continue
                    
                    G[:,iS,0]=G_StemMerch[:,indTIPSY[0]]/meta['GC']['Scale Factor']
                    G[:,iS,1]=G_StemNonMerch[:,indTIPSY[0]]/meta['GC']['Scale Factor']
                    G[:,iS,2]=G_Bark[:,indTIPSY[0]]/meta['GC']['Scale Factor']
                    G[:,iS,3]=G_Branch[:,indTIPSY[0]]/meta['GC']['Scale Factor']
                    G[:,iS,4]=G_Foliage[:,indTIPSY[0]]/meta['GC']['Scale Factor']
                    G[:,iS,5]=G_VStemMerch[:,indTIPSY[0]]/meta['GC']['Scale Factor']                  
                    
                # Save data to file in input variables folder of project
                gu.opickle(meta['Paths']['Project'] + '\\Inputs\\Scenario' + FixFileNum(iScn) + '\\GrowthCurve' + str(iGC+1) + '_Bat' + FixFileNum(iBat) + '.pkl',G)
                    
    return

#%% GET TASS GROWTH CURVES

def GetTASSCurves(meta,iScn,iGC,fny,fnc,fnm):
    
    # Assume just one stand per scenario (i.e., demo)
    iS=0
    
    dfY=pd.read_csv(fny,header=30)
    
    # Initialize age response of net growth
    G=np.zeros((N_Age,1,6),dtype=np.int16)
    
    G[:,iS,0]=G_StemMerch[:,indTIPSY[0]]/meta['GC']['Scale Factor']
    G[:,iS,1]=G_StemNonMerch[:,indTIPSY[0]]/meta['GC']['Scale Factor']
    G[:,iS,2]=G_Bark[:,indTIPSY[0]]/meta['GC']['Scale Factor']
    G[:,iS,3]=G_Branch[:,indTIPSY[0]]/meta['GC']['Scale Factor']
    G[:,iS,4]=G_Foliage[:,indTIPSY[0]]/meta['GC']['Scale Factor']
    G[:,iS,5]=G_VStemMerch[:,indTIPSY[0]]/meta['GC']['Scale Factor']
                    
    # Save data to file in input variables folder of project
    gu.opickle(meta['Paths']['Project'] + '\\Inputs\\Scenario' + FixFileNum(iScn) + '\\GrowthCurve' + str(iGC+1) + '_Bat' + FixFileNum(1) + '.pkl',G)
    
    return    

#%% POST-PROCESS TIPSY GROWTH CURVES
# Nested list, gc[Scenario][Stand][Growth Curve]

def Import_BatchTIPSY_Output(meta):
        
    # Growth curve parameters and TIPSY outputs
    dfPar=pd.read_excel(meta['Paths']['Project'] + '\\Inputs\\GrowthCurvesTIPSY_Parameters.xlsx',sheet_name='Sheet1',skiprows=6)
    txtDat=np.loadtxt(meta['Paths']['Project'] + '\\Inputs\\GrowthCurvesTIPSY_Output.out',skiprows=4)
    
    # TIPSY saves to text file -> convert to dataframe (column names must match TIPSY output file design)
    dfDat=pd.DataFrame(txtDat,columns=meta['GC']['BatchTIPSY Column Names'])

    # Define age vector (must be consistent with how TIPSY was set up)
    Age=np.arange(0,meta['GC']['BatchTIPSY Maximum Age']+1,1)

    # Get dimensions of the TIPSY output file to reshape the data into Age x Stand
    N_Age=Age.size
    N_GC=int(dfDat.shape[0]/N_Age)

    gc=[None]*N_GC
    for i in range(N_GC): 
        gc[i]={}    
    for i in range(len(meta['GC']['BatchTIPSY Column Names'])):
        data=np.reshape(dfDat[meta['GC']['BatchTIPSY Column Names'][i]].values,(N_Age,N_GC),order='F')
        for j in range(N_GC):
            gc[j][meta['GC']['BatchTIPSY Column Names'][i]]=data[:,j]

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
                for i in range(len(meta['GC']['BatchTIPSY Column Names'])):
                    data=np.reshape(dfDat[meta['GC']['BatchTIPSY Column Names'][i]].values,(N_Age,N_GC),order='F')
                    d[meta['GC']['BatchTIPSY Column Names'][i]]=data[:,ind]
                gc0.append(d)
            gc1.append(gc0)
        gc2.append(gc1)
    return gc2

#%% Import growth curves

def Import_CompiledGrowthCurves(meta,scn):
    # *** This needs to be fixed - Vstem is on the end column ***
    # gc[iScn][iGC][iStand]
    gc=[]
    for iScn in range(len(scn)): #range(meta['N Scenario']):
        gc0=[]
        
        gc1=[]
        for iBat in range(0,meta['N Batch']):            
            tmp=gu.ipickle(meta['Paths']['Project'] + '\\Inputs\\Scenario' + FixFileNum(iScn) + '\\GrowthCurve1_Bat' + FixFileNum(iBat) + '.pkl')
            tmp=np.sum(tmp,axis=2).astype(float)
            for iS in range(tmp.shape[1]):
                gc1.append(tmp[:,iS].copy()/meta['GC']['Scale Factor'])
        gc0.append(gc1.copy())
        
        gc1=[]
        for iBat in range(0,meta['N Batch']):            
            tmp=gu.ipickle(meta['Paths']['Project'] + '\\Inputs\\Scenario' + FixFileNum(iScn) + '\\GrowthCurve2_Bat' + FixFileNum(iBat) + '.pkl')
            tmp=np.sum(tmp,axis=2).astype(float)
            for iS in range(tmp.shape[1]):
                gc1.append(tmp[:,iS].copy()/meta['GC']['Scale Factor'])
        gc0.append(gc1.copy())
        
        gc1=[]
        for iBat in range(0,meta['N Batch']):            
            tmp=gu.ipickle(meta['Paths']['Project'] + '\\Inputs\\Scenario' + FixFileNum(iScn) + '\\GrowthCurve3_Bat' + FixFileNum(iBat) + '.pkl')
            tmp=np.sum(tmp,axis=2).astype(float)
            for iS in range(tmp.shape[1]):
                gc1.append(tmp[:,iS].copy()/meta['GC']['Scale Factor'])
        gc0.append(gc1.copy())
        gc.append(gc0.copy())
    return gc

#%% Import custom harvest assumptions (optional)

def ImportCustomHarvestAssumptions(pthin):

    # Import parameters
    df=pd.read_excel(pthin,sheet_name='Sheet1',skiprows=1)
    
    d={}
    d['BiomassMerch_Affected']=df.iloc[1,1]    
    d['BiomassMerch_Removed']=df.iloc[5,1]
    d['BiomassMerch_Burned']=df.iloc[3,1]
    d['BiomassMerch_LeftOnSite']=df.iloc[4,1]    
    d['BiomassNonMerch_Affected']=df.iloc[1,2]
    d['BiomassNonMerch_Removed']=df.iloc[5,2]
    d['BiomassNonMerch_Burned']=df.iloc[3,2]
    d['BiomassNonMerch_LeftOnSite']=df.iloc[4,2]    
    d['Snags_Affected']=df.iloc[1,3]
    d['Snags_Removed']=df.iloc[5,3]
    d['Snags_Burned']=df.iloc[3,3]
    d['Snags_LeftOnSite']=df.iloc[4,3]    
    d['RemovedMerchToPulp']=df.iloc[11,1]
    d['RemovedMerchToFuel']=df.iloc[12,1]
    d['RemovedMerchToLumber']=df.iloc[7,1]
    d['RemovedMerchToPlywood']=df.iloc[8,1]
    d['RemovedMerchToOSB']=df.iloc[9,1]
    d['RemovedMerchToMDF']=df.iloc[10,1]    
    d['RemovedMerchToCants']=df.iloc[13,1]
    d['RemovedMerchToFirewood']=df.iloc[14,1]    
    d['RemovedNonMerchToFuel']=df.iloc[12,2]
    d['RemovedNonMerchToLumber']=df.iloc[7,2]
    d['RemovedNonMerchToPlywood']=df.iloc[8,2]
    d['RemovedNonMerchToOSB']=df.iloc[9,2]
    d['RemovedNonMerchToMDF']=df.iloc[10,2]
    d['RemovedNonMerchToPulp']=df.iloc[11,2]
    d['RemovedNonMerchToCants']=df.iloc[13,2]
    d['RemovedNonMerchToFirewood']=df.iloc[14,2]    
    d['RemovedSnagStemToFuel']=df.iloc[12,3]
    d['RemovedSnagStemToLumber']=df.iloc[7,3]
    d['RemovedSnagStemToPlywood']=df.iloc[8,3]
    d['RemovedSnagStemToOSB']=df.iloc[9,3]
    d['RemovedSnagStemToMDF']=df.iloc[10,3]
    d['RemovedSnagStemToPulp']=df.iloc[11,3]
    d['RemovedSnagStemToCants']=df.iloc[13,3]
    d['RemovedSnagStemToFirewood']=df.iloc[14,3]
    
    return d


#%% Update parameters

def UpdateParamaters(pthin):

    #------------------------------------------------------------------------------
    # Carbon Pools
    #------------------------------------------------------------------------------

    # Define pool names
    PoolNames=['StemMerch','StemNonMerch','Foliage','Branch','Bark','RootCoarse','RootFine', \
        'FelledStemMerch','FelledStemNonMerch','FelledBranch','FelledBark','FelledSnagStem','FelledSnagBranch', \
        'LitterVF','LitterF','LitterM','LitterS','SnagStem','SnagBranch','SoilVF','SoilF','SoilS', \
        'ECO2asC','ECH4asC','ECOasC','EN2OasC','RemovedMerch','RemovedNonMerch','RemovedSnagStem']   

    pPools={'ID_Pool':np.arange(1,len(PoolNames),1),'Name_Pool':PoolNames}

    #--------------------------------------------------------------------------
    # Biomass allometry (stand level)
    #--------------------------------------------------------------------------
    
    df=pd.read_excel(pthin + "\Parameters_BiomassAllometrySL.xlsx",sheet_name='Sheet1')
    m,n=df.shape
    
    pBiomassAllomSL=df.to_dict()
    
    #------------------------------------------------------------------------------
    # Biomass Turnover Rates (from Kurz et al. 2009)
    #------------------------------------------------------------------------------
    
    df=pd.read_excel(pthin + "\Parameters_BiomassTurnover.xlsx",sheet_name='Sheet1')
    m,n=df.shape

    # Pool identifiers
    idp=[0,1,2,3,4,5,6]

    ID_Species=np.zeros(m*len(idp),dtype=int)
    ID_Pool=np.zeros(m*len(idp),dtype=int)
    Value=np.zeros(m*len(idp),dtype=float)
    cnt=0
    for i in range(m):
        for j in range(len(idp)):
            ID_Species[cnt]=int(df['VALUE'][i])
            ID_Pool[cnt]=int(idp[j])
            Value[cnt]=df.iloc[i,j+2]
            cnt=cnt+1

    pBiomassTurnover={'ID_Species':ID_Species,'ID_Pool':ID_Pool,'Value':Value}

    #------------------------------------------------------------------------------
    # Inter-pool carbon fluxes (from Kurz et al. 2009)
    #------------------------------------------------------------------------------

    df=pd.read_excel(pthin + "\Parameters_InterPoolFluxes.xlsx",sheet_name='Sheet1',skiprows=2)
    m,n=df.shape

    Name=[]
    Value=[]
    for i in range(m):    
        for j in range(n-3):
            b=df.iloc[i,j+1]
            if (str(b)!='NaN') & (str(b)!='nan'):
                s=df.iloc[i,0] + 'To' + df.columns[j+1]            
                Name.append(s)
                Value.append(b)            
                
    Name=np.asarray(Name)
    Value=np.asarray(Value)

    pInterPoolFluxes={'Name':Name,'Value':Value}

    #------------------------------------------------------------------------------
    # Decomposition and physical transfer
    #------------------------------------------------------------------------------

    df=pd.read_excel(pthin + "\Parameters_Decomposition.xlsx",sheet_name='Sheet1',skiprows=1)
    m,n=df.shape
    
    Name_Pool=np.asarray(df['Name_Pool'])    
    Rten=np.asarray(df['Rten'])
    Qten=np.asarray(df['Qten'])
    PhysTransRate=np.asarray(df['PhysTransRate'])    

    pDecomp={'Name_Pool':Name_Pool,'Rten':Rten,'Qten':Qten,'PhysTransRate':PhysTransRate}

    #------------------------------------------------------------------------------
    # Distrubances
    #------------------------------------------------------------------------------
    
    # Full disturbance info:
    df=pd.read_excel(pthin + "\Parameters_Disturbances.xlsx",sheet_name='Sheet1')    
    pDistFull=df.to_dict()
    
    # Condensed dictionary
    
    p=gu.ReadExcel(pthin + "\Parameters_Disturbances.xlsx")
    
    str_to_exclude=['ID','Name','Species_CD','MortalityOccurs','GrowthFactor',
                    'GrowthFactor_Source','GrowthRecovery_HL','GrowthRecovery_HL_Source',
                    'QA1','QA2','QA3']    
    
    pDist={}
    for i in range(p['ID'].size):
        pDist[p['ID'][i]]={}
        pDist[p['ID'][i]]['BiomassMerch_Affected']=1
        pDist[p['ID'][i]]['BiomassNonMerch_Affected']=1
        pDist[p['ID'][i]]['Snags_Affected']=1
        for k in p.keys():
            if np.isin(k,str_to_exclude)==True:
                continue
            pDist[p['ID'][i]][k]=np.nan_to_num(p[k][i])

    #------------------------------------------------------------------------------
    # Harvested wood products (from Dymond 2012, Skog 2008)
    #------------------------------------------------------------------------------

    df=pd.read_excel(pthin + "\Parameters_HWP.xlsx",sheet_name='Default')
    m,n=df.shape

    Name=[]
    Value=[]
    for i in range(m):
        Name.append(df.iloc[i,0])
        Value.append(df.iloc[i,1])

    Name=np.asarray(Name)
    Value=np.asarray(Value)

    pHWP={'Name':Name,'Value':Value}

    #------------------------------------------------------------------------------
    # Biogeoclimatic Zone
    #------------------------------------------------------------------------------

    df=pd.read_excel(pthin + "\Parameters_BGC.xlsx",sheet_name='Zone')
    m,n=df.shape

    ID_BGC_ZONE=[]
    CODE_BGC_ZONE=[]
    for j in range(m):
        ID_BGC_ZONE.append(j+1)
        CODE_BGC_ZONE.append(df['Name'][j])

    ID_BGC_ZONE=np.asarray(ID_BGC_ZONE)
    CODE_BGC_ZONE=np.asarray(CODE_BGC_ZONE)

    pBGC_ZONE={'ID_BGC_ZONE':ID_BGC_ZONE,'CODE_BGC_ZONE':CODE_BGC_ZONE}

    #------------------------------------------------------------------------------
    # VRI - Species
    #------------------------------------------------------------------------------

    df=pd.read_excel(pthin + "\Parameters_VRI.xlsx",sheet_name='SPECIES')
    m,n=df.shape

    ID_SPECIES=[]
    SPECIES_CODE=[]
    Descript=[]
    for j in range(m-1):
        ID_SPECIES.append(j+1)
        SPECIES_CODE.append(df['Code'][j+1])
        Descript.append(df['Description'][j+1])

    ID_SPECIES=np.asarray(ID_SPECIES)
    SPECIES_CODE=np.asarray(SPECIES_CODE)
    Descript=np.asarray(Descript)

    pSpeciesVRI={'ID_SPECIES':ID_SPECIES,'SPECIES_CODE':SPECIES_CODE,'Description':Descript}

    #------------------------------------------------------------------------------
    # Biophysical
    #------------------------------------------------------------------------------

    df=pd.read_excel(pthin + "\Parameters_Biophysical.xlsx",sheet_name='Sheet1')
    m,n=df.shape

    HANDLE=[]
    VALUE=[]
    for j in range(m-1):
        HANDLE.append(df['Handle'][j+1])
        VALUE.append(df['Value'][j+1])

    HANDLE=np.asarray(HANDLE)
    VALUE=np.asarray(VALUE)
    
    pBiophysical={'Handle':HANDLE,'Value':VALUE}
    
    #------------------------------------------------------------------------------
    # Nutrient application
    #------------------------------------------------------------------------------

    df=pd.read_excel(pthin + "\Parameters_NutrientApplication.xlsx",sheet_name='Sheet1')
    m,n=df.shape

    HANDLE=[]
    VALUE=[]
    for j in range(m):
        HANDLE.append(df['Handle'][j])
        VALUE.append(df['Value'][j])

    HANDLE=np.asarray(HANDLE)
    VALUE=np.asarray(VALUE)
    
    pNutAdd={'Handle':HANDLE,'Value':VALUE}

    #------------------------------------------------------------------------------
    # SAWTOOTH
    #------------------------------------------------------------------------------
    
    # Species Region Sample
    df=pd.read_excel(pthin + "\Parameters_SRS.xlsx",'Sheet1')
    cn=df.columns
    pSRS={}
    for i in range(df.shape[1]):
        pSRS[cn[i]]=df[cn[i]].values

    # Allometry
    df=pd.read_excel(pthin + "\Parameters_TreeAllometry.xlsx",'Sheet1')
    cn=df.columns
    pTreeAllometry={}
    for i in range(df.shape[1]):
        pTreeAllometry[cn[i]]=df[cn[i]].values

    # Recruitment Default 1
    df=pd.read_excel(pthin + "\Parameters_TreeRecruitmentDef1.xlsx",'Sheet1')
    cn=df.columns
    pR_Def1={}
    for i in range(df.shape[1]):
        pR_Def1[cn[i]]=df[cn[i]].values
    
    # Mortality Default 1
    df=pd.read_excel(pthin + "\Parameters_TreeMortalityDef1.xlsx",'Sheet1')
    cn=df.columns
    pM_Def1={}
    for i in range(df.shape[1]):
        pM_Def1[cn[i]]=df[cn[i]].values

    # Growth Default 1
    df=pd.read_excel(pthin + "\Parameters_TreeGrowthDef1.xlsx",'Sheet1')
    cn=df.columns
    pG_Def1={}
    for i in range(df.shape[1]):
        pG_Def1[cn[i]]=df[cn[i]].values

    #------------------------------------------------------------------------------
    # Nest parameter dictionaries into a dictionary
    #------------------------------------------------------------------------------

    pts={'Pools':pPools, \
         'Biophysical':pBiophysical, \
         'NutrientAddition':pNutAdd, \
         'BiomassAllomSL':pBiomassAllomSL, \
         'BiomassTurnover':pBiomassTurnover, \
         'InterPoolFluxes':pInterPoolFluxes, \
         'Decomposition':pDecomp, \
         'Disturbances':pDistFull, \
         'Dist':pDist, \
         'BGC_ZONE':pBGC_ZONE, \
         'SpeciesVRI':pSpeciesVRI, \
         'HWP':pHWP, \
         'SRS':pSRS, \
         'TreeAllometry':pTreeAllometry, \
         'R_Def1':pR_Def1, \
         'M_Def1':pM_Def1, \
         'G_Def1':pG_Def1, \
         }

    # Save to file
    gu.opickle(pthin + '\Parameters.pkl',pts)
    
#%% Write BatchTIPSY input file
# Notes:
#    if the input spreadsheet has nan's, all data will be converted to float

def Write_BatchTIPSY_Input_File(meta):
    
    # Import format info (number of designated positions per variable)
    fin=meta['Paths']['Model Code'] + '\\Parameters\\GrowthCurvesTIPSY_Parameters_Template.xlsx'
    df_frmt=pd.read_excel(fin,sheet_name='Sheet1')

    # Format array
    nfrmt=df_frmt.iloc[1,4:53]

    # Import input data   
    df=pd.read_excel(meta['Paths']['Project'] + '\\Inputs\\GrowthCurvesTIPSY_Parameters.xlsx',sheet_name='Sheet1',skiprows=6)

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
    fid=open(meta['Paths']['Project'] + '\\Inputs\\GrowthCurvesTIPSY_InputVariables.dat','w')    
    
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

#%% WRITE SPREADSHEET OF BatchTIPSY PARAMTERS

def Write_BatchTIPSY_Input_Spreadsheet(meta,ugc):
    
    # Create a function that will return the column corresponding to a variable name
    fin=meta['Paths']['Model Code'] + '\\Parameters\\GrowthCurvesTIPSY_Parameters_Template.xlsx'
    df_frmt=pd.read_excel(fin,sheet_name='Sheet1')
    gy_labels=df_frmt.loc[5,:].values
    def GetColumn(lab):
        ind=np.where(gy_labels==lab)[0]
        return int(ind+1)

    # Open spreadsheet
    PathGrowthCurveParameters=meta['Paths']['Project'] + '\\Inputs\\GrowthCurvesTIPSY_Parameters.xlsx'
    xfile=openpyxl.load_workbook(PathGrowthCurveParameters)
    sheet=xfile.get_sheet_by_name('Sheet1')    
    N_headers=7

    # Overwrite existing data entries with empty cells
    for i in range(22000):
        for j in range(len(gy_labels)):
            sheet.cell(row=i+1+N_headers,column=j+1).value=''

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
        cd=lut_n2s(meta['LUT']['VRI']['SPECIES_CD_1'],id)[0]
        sheet.cell(row=cnt+N_headers,column=GetColumn(vnam)).value=cd
    
        vnam='p1'
        dat=ugc['Unique'][iUGC,np.where(ugc['GC_Variable_List']==vnam)[0]]
        sheet.cell(row=cnt+N_headers,column=GetColumn(vnam)).value=int(dat)
    
        vnam='i1'
        dat=ugc['Unique'][iUGC,np.where(ugc['GC_Variable_List']==vnam)[0]]
        sheet.cell(row=cnt+N_headers,column=GetColumn(vnam)).value=int(dat)
    
        vnam='gain1'
        dat=ugc['Unique'][iUGC,np.where(ugc['GC_Variable_List']==vnam)[0]]
        if dat!=-999:
            sheet.cell(row=cnt+N_headers,column=GetColumn(vnam)).value=int(dat)
            sheet.cell(row=cnt+N_headers,column=GetColumn(vnam)).value=int(10)
    
        # Species 2
        vnam='s2'
        id=ugc['Unique'][iUGC,np.where(ugc['GC_Variable_List']==vnam)[0]]
        if id!=-999:
            cd=lut_n2s(meta['LUT']['VRI']['SPECIES_CD_1'],id)[0]
            sheet.cell(row=cnt+N_headers,column=GetColumn(vnam)).value=cd
        
            vnam='p2'
            dat=ugc['Unique'][iUGC,np.where(ugc['GC_Variable_List']==vnam)[0]]        
            sheet.cell(row=cnt+N_headers,column=GetColumn(vnam)).value=int(dat)
        
            vnam='gain2'
            dat=ugc['Unique'][iUGC,np.where(ugc['GC_Variable_List']==vnam)[0]]
            if dat!=-999:
                sheet.cell(row=cnt+N_headers,column=GetColumn(vnam)).value=int(dat)
                sheet.cell(row=cnt+N_headers,column=GetColumn(vnam)).value=int(10)
    
        # Species 3
        vnam='s3'
        id=ugc['Unique'][iUGC,np.where(ugc['GC_Variable_List']==vnam)[0]]
        if id!=-999:
            cd=lut_n2s(meta['LUT']['VRI']['SPECIES_CD_1'],id)[0]
            sheet.cell(row=cnt+N_headers,column=GetColumn(vnam)).value=cd
        
            vnam='p3'
            dat=ugc['Unique'][iUGC,np.where(ugc['GC_Variable_List']==vnam)[0]]        
            sheet.cell(row=cnt+N_headers,column=GetColumn(vnam)).value=int(dat)
        
            vnam='gain3'
            dat=ugc['Unique'][iUGC,np.where(ugc['GC_Variable_List']==vnam)[0]]
            if dat!=-999:
                sheet.cell(row=cnt+N_headers,column=GetColumn(vnam)).value=int(dat)
                sheet.cell(row=cnt+N_headers,column=GetColumn(vnam)).value=int(10)
                
        # Species 4
        vnam='s4'
        id=ugc['Unique'][iUGC,np.where(ugc['GC_Variable_List']==vnam)[0]]
        if id!=-999:
            cd=lut_n2s(meta['LUT']['VRI']['SPECIES_CD_1'],id)[0]
            sheet.cell(row=cnt+N_headers,column=GetColumn(vnam)).value=cd
        
            vnam='p4'
            dat=ugc['Unique'][iUGC,np.where(ugc['GC_Variable_List']==vnam)[0]]        
            sheet.cell(row=cnt+N_headers,column=GetColumn(vnam)).value=int(dat)
        
            vnam='gain4'
            dat=ugc['Unique'][iUGC,np.where(ugc['GC_Variable_List']==vnam)[0]]
            if dat!=-999:
                sheet.cell(row=cnt+N_headers,column=GetColumn(vnam)).value=int(dat)
                sheet.cell(row=cnt+N_headers,column=GetColumn(vnam)).value=int(10)
            
        # Species 5
        vnam='s5'
        id=ugc['Unique'][iUGC,np.where(ugc['GC_Variable_List']==vnam)[0]]
        if id!=-999:
            cd=lut_n2s(meta['LUT']['VRI']['SPECIES_CD_1'],id)[0]
            sheet.cell(row=cnt+N_headers,column=GetColumn(vnam)).value=cd
        
            vnam='p5'
            dat=ugc['Unique'][iUGC,np.where(ugc['GC_Variable_List']==vnam)[0]]        
            sheet.cell(row=cnt+N_headers,column=GetColumn(vnam)).value=int(dat)
        
            vnam='gain5'
            dat=ugc['Unique'][iUGC,np.where(ugc['GC_Variable_List']==vnam)[0]]
            if dat!=-999:
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
        cd=lut_n2s(meta['LUT']['VRI']['BEC_ZONE_CODE'],id)[0]    
        sheet.cell(row=cnt+N_headers,column=GetColumn(vnam)).value=cd
    
        # FIZ
        vnam='FIZ'
        id=ugc['Unique'][iUGC,np.where(ugc['GC_Variable_List']==vnam)[0]]
        cd=lut_n2s(meta['LUT']['TIPSY']['FIZ'],id)[0]
        sheet.cell(row=cnt+N_headers,column=GetColumn(vnam)).value=cd
    
        # Age at fertilization
        # vnam='fert_age1'
        #dat=ugc['Unique'][iUGC,np.where(ugc['GC_Variable_List']==vnam)[0]]
        #if dat!=-999:
        #    sheet.cell(row=cnt+N_headers,column=GetColumn(vnam)).value=int(dat[0])
    
        # Update counter
        cnt=cnt+1
    
    #------------------------------------------------------------------------------
    # Save to spreadsheet
    #------------------------------------------------------------------------------
    
    xfile.save(meta['Paths']['Project'] + '\\Inputs\\GrowthCurvesTIPSY_Parameters.xlsx')

#%% IMPORT DISTURBANCE HISTORY

def GetDisturbanceHistory(meta):
    dh=[]
    for iScn in range(meta['N Scenario']):
        dhB=[]    
        for iBat in range(meta['N Batch']):
            iEns=0
            dh0=gu.ipickle(meta['Paths']['Input Scenario'][iScn] + '\\Events_Ens' + FixFileNum(iEns) + '_Bat' + FixFileNum(iBat) + '.pkl')
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

def Import_GraphicsParameters(type):
    if type=='FCI Demo':
        fs1=7
        fs2=9
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
    return params

#%% Prepare inventory from spreadsheet

def PrepareInventoryFromSpreadsheet(meta):
    
    for iScn in range(0,meta['N Scenario']):
    
        # Loop through batches, saving inventory to file
        for iBat in range(0,meta['N Batch']):
      
            inv={}
    
            # Index to batch
            indBat=IndexToBatch(meta,iBat)    
            N_StandsInBatch=len(indBat)
    
            # Initialize inventory variables
            inv['Lat']=np.zeros((1,N_StandsInBatch))
            inv['Lon']=np.zeros((1,N_StandsInBatch))
            inv['X']=inv['Lat']
            inv['Y']=inv['Lon']
        
            # BEC zone
            inv['ID_BECZ']=np.zeros((1,N_StandsInBatch),dtype=np.int)
            inv['ID_BECZ'][0,:]=meta['LUT']['VRI']['BEC_ZONE_CODE'][meta['Scenario'][iScn]['BGC Zone Code']]
    
            # Timber harvesting landbase (1=yes, 0=no)
            inv['THLB']=1*np.ones((1,N_StandsInBatch))
            inv['THLB'][0,:]=meta['Scenario'][iScn]['THLB Status']
        
            # Temperature will be updated automatically
            inv['MAT']=4*np.ones((1,N_StandsInBatch))
            
            if meta['Biomass Module']=='Sawtooth':
                inv['Srs1_ID']=meta['LUT']['Spc'][meta['Scenario'][iScn]['SRS1_CD']]*np.ones((1,N_StandsInBatch),dtype=np.int)
            else:
                inv['Srs1_ID']=9999*np.ones((1,N_StandsInBatch),dtype=np.int)
                inv['Srs1_Pct']=100*np.ones((1,N_StandsInBatch),dtype=np.int)
                inv['Srs2_ID']=0*np.ones((1,N_StandsInBatch),dtype=np.int)
                inv['Srs2_Pct']=0*np.ones((1,N_StandsInBatch),dtype=np.int)
                inv['Srs3_ID']=0*np.ones((1,N_StandsInBatch),dtype=np.int)
                inv['Srs3_Pct']=0*np.ones((1,N_StandsInBatch),dtype=np.int)
                inv['Spc1_ID']=0*np.ones((1,N_StandsInBatch),dtype=np.int)
                inv['Spc1_Pct']=0*np.ones((1,N_StandsInBatch),dtype=np.int)
                inv['Spc2_ID']=0*np.ones((1,N_StandsInBatch),dtype=np.int)
                inv['Spc2_Pct']=0*np.ones((1,N_StandsInBatch),dtype=np.int)
                inv['Spc3_ID']=0*np.ones((1,N_StandsInBatch),dtype=np.int)
                inv['Spc3_Pct']=0*np.ones((1,N_StandsInBatch),dtype=np.int)

            # Save
            gu.opickle(meta['Paths']['Input Scenario'][iScn] + '\\Inventory_Bat' + FixFileNum(iBat) + '.pkl',inv)
        
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
                indAvailable=np.where(ec['ID_Type'][iT,iS,:].flatten()==0)[0]
                if indAvailable.size==0:
                    print('Warning, more events per year than can be handled!')
                else:
                    iE=indAvailable[0]
                    ec['ID_Type'][iT,iS,iE]=ID_Type[indYear[iY]]
                    ec['MortalityFactor'][iT,iS,iE]=MortalityFactor[indYear[iY]]
                    ec['GrowthFactor'][iT,iS,iE]=GrowthFactor[indYear[iY]]
                    ec['ID_GrowthCurve'][iT,iS,iE]=ID_GrowthCurve[indYear[iY]]
    return ec

#%% Mortality frequency distribution

def GetMortalityFrequencyDistribution(meta):
    
    iEns=0
    
    M=[None]*meta['N Scenario']    
    for iScn in range(meta['N Scenario']):
        
        tv=np.arange(meta['Year Start Saving'],meta['Year End']+1,1)
        M[iScn]={}
        M[iScn]['Ma']={}
        M[iScn]['Mr']={}
        for k in meta['LUT']['Dist'].keys():
            M[iScn]['Ma'][k]=np.zeros((tv.size,meta['N Stand Full']))
            M[iScn]['Mr'][k]=np.zeros((tv.size,meta['N Stand Full']))
        M[iScn]['Ma']['Reg']=np.zeros((tv.size,meta['N Stand Full']))
        M[iScn]['Mr']['Reg']=np.zeros((tv.size,meta['N Stand Full']))
        
        for iBat in range(meta['N Batch']): 
            
            indBat=IndexToBatch(meta,iBat)
            
            d1=LoadSingleOutputFile(meta,iScn,iEns,iBat)        
            dh=gu.ipickle(meta['Paths']['Input Scenario'][iScn] + '\\Events_Ens' + FixFileNum(iEns) + '_Bat' + FixFileNum(iBat) + '.pkl')
            
            for iStandInBat in range(len(indBat)):
                iStand=indBat[iStandInBat]
                for iYr in range(dh[iStandInBat]['Year'].size):
                    it=np.where(tv==int(dh[iStandInBat]['Year'][iYr]))[0]
                    if it.size==0:
                        continue
                    nam=lut_n2s(meta['LUT']['Dist'],dh[iStandInBat]['ID_Type'][iYr])[0]                     
                    M[iScn]['Ma'][nam][it,iStand]=d1[0]['C_M_Dist'][it,iStandInBat]
                    M[iScn]['Mr'][nam][it,iStand]=d1[0]['C_M_Dist'][it,iStandInBat]/np.maximum(0.0001,d1[0]['Eco_Biomass'][it-1,iStandInBat])
                    M[iScn]['Ma']['Reg'][it,iStand]=d1[0]['C_M_Reg'][it,iStandInBat]
                    M[iScn]['Mr']['Reg'][it,iStand]=d1[0]['C_M_Reg'][it,iStandInBat]/np.maximum(0.0001,d1[0]['Eco_Biomass'][it-1,iStandInBat])
            del d1,dh
            garc.collect()
    
    return M

#%% Save output variables by multipolygon
# This is only designed for projects that apply inventory_from_polygons method.
# Not only is this outputting by project, but the project flux sums are accurately
# representing the given treatment area, rather than the area inferred from the 
# total number of sparse sample points within a project area, which can be wrong.

def MosByMultipolygon(meta,include_area):

    # Import multipolygons
    atu_multipolygons=gu.ipickle(meta['Paths']['Geospatial'] + '\\atu_multipolygons.pkl')

    # Import sxy
    sxy=gu.ipickle(meta['Paths']['Geospatial'] + '\\sxy.pkl')

    uMP=np.unique(sxy['ID_atu_multipolygons'])

    # Create listed index (faster than indexing on the fly)
    Crosswalk_sxy_to_mp=[None]*uMP.size
    for iMP in range(uMP.size):
        d={}
        d['Index']=np.where(sxy['ID_atu_multipolygons']==uMP[iMP])[0]
        Crosswalk_sxy_to_mp[iMP]=d

    # Time series of saved results
    tv_full=np.arange(meta['Year Start'],meta['Year End']+1,1)
    tv_saving=np.arange(meta['Year Start Saving'],meta['Year End']+1,1)
    it=np.where( (tv_full>=tv_saving[0]) & (tv_full<=tv_saving[-1]) )[0]
    
    # Variables to save
    nam1=['V_StemMerch','C_Ecosystem', 'C_InUse', 'C_DumpLandfill','C_RemovedMerch','C_RemovedNonMerch','C_RemovedSnagStem']
    nam2=['A','Eco_Biomass','Eco_DeadWood','Eco_Litter','Eco_Total','Eco_RH','Eco_E_Wildfire','Eco_E_OpenBurning','Eco_E_Operations','Eco_Removals','Eco_NGHGB','Sec_NGHGB']

    # Initialize data by multipolygon structure    
    MosByMP=[None]*meta['N Scenario']
    for iScn in range(meta['N Scenario']):
        d={}
        d['v1']={}
        d['v1']['Mean']={}
        d['v1']['Sum']={}
        for iV in range(len(nam1)):
            d['v1']['Mean'][nam1[iV]]={}
            d['v1']['Mean'][nam1[iV]]['Ensemble Mean']=np.zeros((tv_saving.size,uMP.size))
            d['v1']['Mean'][nam1[iV]]['Ensemble SD']=np.zeros((tv_saving.size,uMP.size))
            d['v1']['Sum'][nam1[iV]]={}
            d['v1']['Sum'][nam1[iV]]['Ensemble Mean']=np.zeros((tv_saving.size,uMP.size))
            d['v1']['Sum'][nam1[iV]]['Ensemble SD']=np.zeros((tv_saving.size,uMP.size))
        d['v2']={}
        d['v2']['Mean']={}
        d['v2']['Sum']={}
        for iV in range(len(nam2)):
            d['v2']['Mean'][nam2[iV]]={}
            d['v2']['Mean'][nam2[iV]]['Ensemble Mean']=np.zeros((tv_saving.size,uMP.size))
            d['v2']['Mean'][nam2[iV]]['Ensemble SD']=np.zeros((tv_saving.size,uMP.size))
            d['v2']['Sum'][nam2[iV]]={}
            d['v2']['Sum'][nam2[iV]]['Ensemble Mean']=np.zeros((tv_saving.size,uMP.size))
            d['v2']['Sum'][nam2[iV]]['Ensemble SD']=np.zeros((tv_saving.size,uMP.size))
        if include_area=='On':
            d['Area']={}
            for k in meta['LUT']['Dist'].keys():
                d['Area'][k]={}
                d['Area'][k]['Ensemble Mean']=np.zeros((tv_saving.size,uMP.size))
                #d['Area'][k]['Ensemble SD']=np.zeros((tv_saving.size,uMP.size))
        MosByMP[iScn]=d
    
    # Scale factor used to temporarily store data
    ScaleFactor=0.01
    
    # Loop through scenarios
    for iScn in range(meta['N Scenario']):
    
        for iEns in range(meta['N Ensemble']):
            
            # Initialize temporary data structure for full simulation
            # Give each variable two columns for instances where areas overlap multiple batches.
            # Values for areas with data from multiple batches will be averaged.        
            Data={}
            Data['v1']={}
            for iV in range(len(nam1)):
                Data['v1'][nam1[iV]]=np.zeros((tv_saving.size,meta['N Stand Full']),dtype=int)        
            Data['v2']={}
            for iV in range(len(nam2)):
                Data['v2'][nam2[iV]]=np.zeros((tv_saving.size,meta['N Stand Full']),dtype=int)
            if include_area=='On':
                Data['Area']={}
                for k in MosByMP[iScn]['Area']:
                    Data['Area'][k]=np.zeros((tv_saving.size,meta['N Stand Full']),dtype=int)

            # Populate full simulation results       
            for iBat in range(meta['N Batch']):
        
                indBat=IndexToBatch(meta,iBat)            
            
                d1=LoadSingleOutputFile(meta,iScn,iEns,iBat)
                d2=CalculateGHGBalance(d1,meta)        
                for iV in range(len(nam1)):
                    tmp=d1[nam1[iV]]/ScaleFactor
                    Data['v1'][nam1[iV]][:,indBat]=tmp.copy().astype(int)            
                for iV in range(len(nam2)):
                    tmp=d2[0][0][nam2[iV]]/ScaleFactor
                    Data['v2'][nam2[iV]][:,indBat]=tmp.copy().astype(int)
            
                if include_area=='On':
                    # Import event chronology
                    ec=gu.ipickle(meta['Paths']['Input Scenario'][iScn] + '\\Events_Ens' + FixFileNum(iEns) + '_Bat' + FixFileNum(iBat) + '.pkl')
            
                    # Uncompress event chronology if it has been compressed
                    if 'idx' in ec:
                        idx=ec['idx']
                        tmp=ec.copy()
                        for v in ['ID_Type','MortalityFactor','GrowthFactor','ID_GrowthCurve']:
                            ec[v]=np.zeros((tv_full.size,d1['A'].shape[1],meta['Max Events Per Year']),dtype='int16')
                            ec[v][idx[0],idx[1],idx[2]]=tmp[v]
                    del tmp            
                
                    for k in Data['Area'].keys():
                        Data0=np.zeros((tv_saving.size,indBat.size))
                        for iEY in range(meta['Max Events Per Year']):
                            ind=np.where(ec['ID_Type'][it,:,iEY]==meta['LUT']['Dist'][k])[0]
                            Data0[ind]=Data0[ind]+1
                        Data['Area'][k][:,indBat]=Data0
                    
                    del ec
            
                del d1,d2
                garc.collect()
        
            if include_area=='On':
                # Populating the final structure with area data is slow - get a flag
                # indicator of whether each event ID can be skipped because it has no
                # info
                flg_area={}
                for k in Data['Area'].keys():            
                    if np.sum(Data['Area'][k])>0:
                        flg_area[k]=1
                    else:
                        flg_area[k]=0
        
            # Calculate stats and populate results for each treatment area    
            for iMP in range(uMP.size):
                #print(iMP)
                ATA=atu_multipolygons[uMP[iMP]]['ACTUAL_TREATMENT_AREA']
                if ATA==None:
                    print('Encounterd no area, using zero')
                    ATA=0
                    
                ind=Crosswalk_sxy_to_mp[iMP]['Index']
                for iV in range(len(nam1)):
                    tmp=ScaleFactor*Data['v1'][nam1[iV]][:,ind].astype(float)
                    MosByMP[iScn]['v1']['Mean'][nam1[iV]]['Ensemble Mean'][:,iMP]=MosByMP[iScn]['v1']['Mean'][nam1[iV]]['Ensemble Mean'][:,iMP]+np.mean(tmp,axis=1)
                    MosByMP[iScn]['v1']['Mean'][nam1[iV]]['Ensemble SD'][:,iMP]=MosByMP[iScn]['v1']['Mean'][nam1[iV]]['Ensemble SD'][:,iMP]+np.std(tmp,axis=1)
                    MosByMP[iScn]['v1']['Sum'][nam1[iV]]['Ensemble Mean'][:,iMP]=MosByMP[iScn]['v1']['Sum'][nam1[iV]]['Ensemble Mean'][:,iMP]+ATA*np.mean(tmp,axis=1)
                    MosByMP[iScn]['v1']['Sum'][nam1[iV]]['Ensemble SD'][:,iMP]=MosByMP[iScn]['v1']['Sum'][nam1[iV]]['Ensemble SD'][:,iMP]+ATA*np.std(tmp,axis=1)
                for iV in range(len(nam2)):
                    tmp=ScaleFactor*Data['v2'][nam2[iV]][:,ind].astype(float)
                    MosByMP[iScn]['v2']['Mean'][nam2[iV]]['Ensemble Mean'][:,iMP]=MosByMP[iScn]['v2']['Mean'][nam2[iV]]['Ensemble Mean'][:,iMP]+np.mean(tmp,axis=1)
                    MosByMP[iScn]['v2']['Mean'][nam2[iV]]['Ensemble SD'][:,iMP]=MosByMP[iScn]['v2']['Mean'][nam2[iV]]['Ensemble SD'][:,iMP]+np.std(tmp,axis=1)
                    MosByMP[iScn]['v2']['Sum'][nam2[iV]]['Ensemble Mean'][:,iMP]=MosByMP[iScn]['v2']['Sum'][nam2[iV]]['Ensemble Mean'][:,iMP]+ATA*np.mean(tmp,axis=1)
                    MosByMP[iScn]['v2']['Sum'][nam2[iV]]['Ensemble SD'][:,iMP]=MosByMP[iScn]['v2']['Sum'][nam2[iV]]['Ensemble SD'][:,iMP]+ATA*np.std(tmp,axis=1)
                if include_area=='On':
                    for k in MosByMP[iScn]['Area']:
                        if flg_area[k]==1:
                            # Only continue if there are some events
                            MosByMP[iScn]['Area'][k]['Ensemble Mean'][:,iMP]=MosByMP[iScn]['Area'][k]['Ensemble Mean'][:,iMP]+np.sum(Data['Area'][k][:,ind],axis=1)

    # Divide by number of ensembles
    for iScn in range(meta['N Scenario']):
        for iV in range(len(nam1)):
            MosByMP[iScn]['v1']['Mean'][nam1[iV]]['Ensemble Mean']=MosByMP[iScn]['v1']['Mean'][nam1[iV]]['Ensemble Mean']/meta['N Ensemble']
            MosByMP[iScn]['v1']['Mean'][nam1[iV]]['Ensemble SD']=MosByMP[iScn]['v1']['Mean'][nam1[iV]]['Ensemble SD']/meta['N Ensemble']
            MosByMP[iScn]['v1']['Sum'][nam1[iV]]['Ensemble Mean']=MosByMP[iScn]['v1']['Sum'][nam1[iV]]['Ensemble Mean']/meta['N Ensemble']
            MosByMP[iScn]['v1']['Sum'][nam1[iV]]['Ensemble SD']=MosByMP[iScn]['v1']['Sum'][nam1[iV]]['Ensemble SD']/meta['N Ensemble']
        for iV in range(len(nam2)):
            MosByMP[iScn]['v2']['Mean'][nam2[iV]]['Ensemble Mean']=MosByMP[iScn]['v2']['Mean'][nam2[iV]]['Ensemble Mean']/meta['N Ensemble']
            MosByMP[iScn]['v2']['Mean'][nam2[iV]]['Ensemble SD']=MosByMP[iScn]['v2']['Mean'][nam2[iV]]['Ensemble SD']/meta['N Ensemble']
            MosByMP[iScn]['v2']['Sum'][nam2[iV]]['Ensemble Mean']=MosByMP[iScn]['v2']['Sum'][nam2[iV]]['Ensemble Mean']/meta['N Ensemble']
            MosByMP[iScn]['v2']['Sum'][nam2[iV]]['Ensemble SD']=MosByMP[iScn]['v2']['Sum'][nam2[iV]]['Ensemble SD']/meta['N Ensemble']
        if include_area=='On':
            for iV in MosByMP[iScn]['Area']:
                MosByMP[iScn]['Area'][iV]['Ensemble Mean'][:,iMP]=MosByMP[iScn]['Area'][iV]['Ensemble Mean'][:,iMP]/meta['N Ensemble']

    # Save
    gu.opickle(meta['Paths']['Project'] + '\\Outputs\\MosByMultipolygon.pkl',MosByMP)


#%% CALCULATE MODEL OUTPUT STATISTICS

def ModelOutputStats(meta,flag_save):

    mos=[None]*meta['N Scenario']
    
    for iScn in range(meta['N Scenario']):
        
        tv=np.arange(meta['Year Start Saving'],meta['Year End']+1,1)
        tv_full=np.arange(meta['Year Start'],meta['Year End']+1,1)
        
        # Initialize structure
        mos[iScn]={}
        
        mos[iScn]['v1']={}
        mos[iScn]['v1']['Mean']={}
        mos[iScn]['v1']['Sum']={}
        d1=LoadSingleOutputFile(meta,0,0,0)
        for k in d1.keys(): 
            if k=='Year':
                continue
            if k=='C_M_ByAgent':
                continue
            mos[iScn]['v1']['Mean'][k]={}
            mos[iScn]['v1']['Mean'][k]['Ensembles']=np.zeros((tv.size,meta['N Ensemble']))
            mos[iScn]['v1']['Mean'][k]['Ensemble Mean']=np.zeros(tv.size)
            mos[iScn]['v1']['Mean'][k]['Ensemble SD']=np.zeros(tv.size)
            mos[iScn]['v1']['Sum'][k]={}
            mos[iScn]['v1']['Sum'][k]['Ensembles']=np.zeros((tv.size,meta['N Ensemble']))
            mos[iScn]['v1']['Sum'][k]['Ensemble Mean']=np.zeros(tv.size)
            mos[iScn]['v1']['Sum'][k]['Ensemble SD']=np.zeros(tv.size)
            
        mos[iScn]['v2']={}
        mos[iScn]['v2']['Mean']={}
        mos[iScn]['v2']['Sum']={}
        d2=CalculateGHGBalance(d1,meta)        
        for k in d2[0][0].keys(): 
            if (k=='Year'):
                continue
            mos[iScn]['v2']['Mean'][k]={}
            mos[iScn]['v2']['Mean'][k]['Ensembles']=np.zeros((tv.size,meta['N Ensemble']))
            mos[iScn]['v2']['Mean'][k]['Ensemble Mean']=np.zeros(tv.size)
            mos[iScn]['v2']['Mean'][k]['Ensemble SD']=np.zeros(tv.size)
            mos[iScn]['v2']['Sum'][k]={}
            mos[iScn]['v2']['Sum'][k]['Ensembles']=np.zeros((tv.size,meta['N Ensemble']))
            mos[iScn]['v2']['Sum'][k]['Ensemble Mean']=np.zeros(tv.size)
            mos[iScn]['v2']['Sum'][k]['Ensemble SD']=np.zeros(tv.size)
        
        mos[iScn]['C_M_ByAgent']={}
        mos[iScn]['C_M_ByAgent']['Mean']={}
        mos[iScn]['C_M_ByAgent']['Sum']={}   
        for k in d1['C_M_ByAgent'].keys(): 
            mos[iScn]['C_M_ByAgent']['Mean'][k]={}
            mos[iScn]['C_M_ByAgent']['Mean'][k]['Ensembles']=np.zeros((tv.size,meta['N Ensemble']))
            mos[iScn]['C_M_ByAgent']['Mean'][k]['Ensemble Mean']=np.zeros(tv.size)
            mos[iScn]['C_M_ByAgent']['Mean'][k]['Ensemble SD']=np.zeros(tv.size)
            mos[iScn]['C_M_ByAgent']['Sum'][k]={}
            mos[iScn]['C_M_ByAgent']['Sum'][k]['Ensembles']=np.zeros((tv.size,meta['N Ensemble']))
            mos[iScn]['C_M_ByAgent']['Sum'][k]['Ensemble Mean']=np.zeros(tv.size)
            mos[iScn]['C_M_ByAgent']['Sum'][k]['Ensemble SD']=np.zeros(tv.size)
        
        mos[iScn]['Area']={}
        for k in meta['LUT']['Dist'].keys():
            mos[iScn]['Area'][k]={}
            mos[iScn]['Area'][k]['Ensembles']=np.zeros((tv.size,meta['N Ensemble']))
            mos[iScn]['Area'][k]['Ensemble Mean']=np.zeros(tv.size)
            mos[iScn]['Area'][k]['Ensemble SD']=np.zeros(tv.size)
            
        # Loop through ensembles
        for iEns in range(meta['N Ensemble']):
            
            v1={}
            v2={}          
            C_M_ByAgent={}
            for iBat in range(meta['N Batch']):
            
                # Basic output
                d1=LoadSingleOutputFile(meta,iScn,iEns,iBat)
            
                for k in d1.keys(): 
                    
                    if k=='Year':
                        continue
                    
                    # Skip mortality summary by agent
                    if k=='C_M_ByAgent':
                        continue
                    
                    if iBat==0:
                        v1[k]=np.sum(d1[k],axis=1)
                    else:
                        v1[k]=v1[k]+np.sum(d1[k],axis=1)
                
                # Mortality summary by agent                
                for k in d1['C_M_ByAgent'].keys():
                    if iBat==0:
                        C_M_ByAgent[k]=d1['C_M_ByAgent'][k].flatten()
                    else:
                        C_M_ByAgent[k]=C_M_ByAgent[k]+d1['C_M_ByAgent'][k].flatten()
                
                # CO2e fluxes and pools
                d2=CalculateGHGBalance(d1,meta)
                
                for k in d2[0][0].keys(): 
                    
                    if k=='Year':
                        continue
                    
                    if iBat==0:
                        v2[k]=np.sum(d2[0][0][k],axis=1)
                    else:
                        v2[k]=v2[k]+np.sum(d2[0][0][k],axis=1)    
            
                # Import event chronology
                ec=gu.ipickle(meta['Paths']['Input Scenario'][iScn] + '\\Events_Ens' + FixFileNum(iEns) + '_Bat' + FixFileNum(iBat) + '.pkl')
            
                # Uncompress event chronology if it has been compressed
                if 'idx' in ec:
                    idx=ec['idx']
                    tmp=ec.copy()
                    for v in ['ID_Type','MortalityFactor','GrowthFactor','ID_GrowthCurve']:
                        ec[v]=np.zeros((tv_full.size,d1['A'].shape[1],meta['Max Events Per Year']),dtype='int16')
                        ec[v][idx[0],idx[1],idx[2]]=tmp[v]
                del tmp
            
                for iYr in range(tv_full.size):
                    it=np.where(tv==tv_full[iYr])[0]
                    if it.size==0:
                        continue
                    ID_Type0=ec['ID_Type'][iYr,:,:].flatten()
                    u=np.unique(ID_Type0)
                    for iU in range(u.size):
                        if u[iU]==0:
                            continue
                        id=lut_n2s(meta['LUT']['Dist'],u[iU])[0]
                        ind=np.where(ID_Type0==u[iU])[0]
                        mos[iScn]['Area'][id]['Ensembles'][it,iEns]=mos[iScn]['Area'][id]['Ensembles'][it,iEns]+ind.size
                del d1,d2,ec
                garc.collect()
            
            # Populate mos for each scenario
            for k in v1.keys():
                
                if k=='Year':
                    continue
                
                # Skip mortality summary by agent
                if k=='C_M_ByAgent':
                    continue
                
                mos[iScn]['v1']['Sum'][k]['Ensembles'][:,iEns]=v1[k].copy()
                mos[iScn]['v1']['Mean'][k]['Ensembles'][:,iEns]=v1[k].copy()/meta['N Stand Full']
            
            for k in C_M_ByAgent.keys():                
                mos[iScn]['C_M_ByAgent']['Sum'][k]['Ensembles'][:,iEns]=C_M_ByAgent[k].copy()
                mos[iScn]['C_M_ByAgent']['Mean'][k]['Ensembles'][:,iEns]=C_M_ByAgent[k].copy()/meta['N Stand Full']
            
            for k in v2.keys():
                if k=='Year':
                    continue
                mos[iScn]['v2']['Sum'][k]['Ensembles'][:,iEns]=v2[k].copy()
                mos[iScn]['v2']['Mean'][k]['Ensembles'][:,iEns]=v2[k].copy()/meta['N Stand Full']                                   
        
        for k in v1.keys():
            if k=='Year':
                continue
            mos[iScn]['v1']['Sum'][k]['Ensemble Mean']=np.mean(mos[iScn]['v1']['Sum'][k]['Ensembles'],axis=1)
            mos[iScn]['v1']['Sum'][k]['Ensemble SD']=np.std(mos[iScn]['v1']['Sum'][k]['Ensembles'],axis=1)
            mos[iScn]['v1']['Mean'][k]['Ensemble Mean']=np.mean(mos[iScn]['v1']['Mean'][k]['Ensembles'],axis=1)
            mos[iScn]['v1']['Mean'][k]['Ensemble SD']=np.std(mos[iScn]['v1']['Mean'][k]['Ensembles'],axis=1)
        
        for k in C_M_ByAgent.keys():
            mos[iScn]['C_M_ByAgent']['Sum'][k]['Ensemble Mean']=np.mean(mos[iScn]['C_M_ByAgent']['Sum'][k]['Ensembles'],axis=1)
            mos[iScn]['C_M_ByAgent']['Sum'][k]['Ensemble SD']=np.std(mos[iScn]['C_M_ByAgent']['Sum'][k]['Ensembles'],axis=1)
            mos[iScn]['C_M_ByAgent']['Mean'][k]['Ensemble Mean']=np.mean(mos[iScn]['C_M_ByAgent']['Mean'][k]['Ensembles'],axis=1)
            mos[iScn]['C_M_ByAgent']['Mean'][k]['Ensemble SD']=np.std(mos[iScn]['C_M_ByAgent']['Mean'][k]['Ensembles'],axis=1)
        
        for k in v2.keys():
            if k=='Year':
                continue
            mos[iScn]['v2']['Sum'][k]['Ensemble Mean']=np.mean(mos[iScn]['v2']['Sum'][k]['Ensembles'],axis=1)
            mos[iScn]['v2']['Sum'][k]['Ensemble SD']=np.std(mos[iScn]['v2']['Sum'][k]['Ensembles'],axis=1)
            mos[iScn]['v2']['Mean'][k]['Ensemble Mean']=np.mean(mos[iScn]['v2']['Mean'][k]['Ensembles'],axis=1)
            mos[iScn]['v2']['Mean'][k]['Ensemble SD']=np.std(mos[iScn]['v2']['Mean'][k]['Ensembles'],axis=1)
        
        for k in mos[iScn]['Area'].keys():
            mos[iScn]['Area'][k]['Ensemble Mean']=np.mean(mos[iScn]['Area'][k]['Ensembles'],axis=1)
            mos[iScn]['Area'][k]['Ensemble SD']=np.std(mos[iScn]['Area'][k]['Ensembles'],axis=1)
            mos[iScn]['Area'][k]['Ensemble Mean']=np.mean(mos[iScn]['Area'][k]['Ensembles'],axis=1)
            mos[iScn]['Area'][k]['Ensemble SD']=np.std(mos[iScn]['Area'][k]['Ensembles'],axis=1)
        
    # Save
    if flag_save=='On':
        gu.opickle(meta['Paths']['Project'] + '\\Outputs\\MOS.pkl',mos)
    
    return mos

#%% Summarize affected area due to natural disturbance and management

def SummarizeAreaAffected(meta,iScn,iEns,AEF,ivlT,tv,mos):
    
    A={}
    if iEns>=0:
        # Individual ensemble        
        A['Nat Dist']=[None]*10; c=-1
        c=c+1; A['Nat Dist'][c]={}; A['Nat Dist'][c]['Name']='Wildfire'; A['Nat Dist'][c]['Color']=[0.75,0,0]; A['Nat Dist'][c]['Data']=mos[iScn]['Area']['Wildfire']['Ensembles'][:,iEns]
        c=c+1; A['Nat Dist'][c]={}; A['Nat Dist'][c]['Name']='Mountain Pine beetle'; A['Nat Dist'][c]['Color']=[0,0.8,0]; A['Nat Dist'][c]['Data']=mos[iScn]['Area']['IBM']['Ensembles'][:,iEns]
        c=c+1; A['Nat Dist'][c]={}; A['Nat Dist'][c]['Name']='Douglas-fir beetle'; A['Nat Dist'][c]['Color']=[0.6,1,0]; A['Nat Dist'][c]['Data']=mos[iScn]['Area']['IBD']['Ensembles'][:,iEns]
        c=c+1; A['Nat Dist'][c]={}; A['Nat Dist'][c]['Name']='Spruce beetle'; A['Nat Dist'][c]['Color']=[0.25,1,1]; A['Nat Dist'][c]['Data']=mos[iScn]['Area']['IBS']['Ensembles'][:,iEns]
        c=c+1; A['Nat Dist'][c]={}; A['Nat Dist'][c]['Name']='W. balsam beetle'; A['Nat Dist'][c]['Color']=[0,0.45,0]; A['Nat Dist'][c]['Data']=mos[iScn]['Area']['IBB']['Ensembles'][:,iEns]
        c=c+1; A['Nat Dist'][c]={}; A['Nat Dist'][c]['Name']='Other pests'; A['Nat Dist'][c]['Color']=[0.8,1,0]; A['Nat Dist'][c]['Data']=mos[iScn]['Area']['Beetles']['Ensembles'][:,iEns]
        c=c+1; A['Nat Dist'][c]={}; A['Nat Dist'][c]['Name']='W. spruce budworm'; A['Nat Dist'][c]['Color']=[0,0.75,1]; A['Nat Dist'][c]['Data']=mos[iScn]['Area']['IDW']['Ensembles'][:,iEns]
        c=c+1; A['Nat Dist'][c]={}; A['Nat Dist'][c]['Name']='Rust'; A['Nat Dist'][c]['Color']=[0.75,0.5,1]; A['Nat Dist'][c]['Data']=mos[iScn]['Area']['Rust Onset']['Ensembles'][:,iEns]
        c=c+1; A['Nat Dist'][c]={}; A['Nat Dist'][c]['Name']='Dwarf Mistletoe'; A['Nat Dist'][c]['Color']=[1,0.5,0.25]; A['Nat Dist'][c]['Data']=mos[iScn]['Area']['Dwarf Mistletoe Onset']['Ensembles'][:,iEns]
        c=c+1; A['Nat Dist'][c]={}; A['Nat Dist'][c]['Name']='Mechanical'; A['Nat Dist'][c]['Color']=[0,0,0.6]; A['Nat Dist'][c]['Data']=mos[iScn]['Area']['Mechanical']['Ensembles'][:,iEns]
        A['Nat Dist']=A['Nat Dist'][0:c+1]

        A['Management']=[None]*10; c=-1
        c=c+1; A['Management'][c]={}; A['Management'][c]['Name']='Harvest'; A['Management'][c]['Color']=[0,0.75,1]; A['Management'][c]['Data']=mos[iScn]['Area']['Harvest']['Ensembles'][:,iEns]+mos[iScn]['Area']['Salvage Logging']['Ensembles'][:,iEns]
        c=c+1; A['Management'][c]={}; A['Management'][c]['Name']='Slashpile burn'; A['Management'][c]['Color']=[0.75,0,0]; A['Management'][c]['Data']=mos[iScn]['Area']['Slashpile Burn']['Ensembles'][:,iEns]+mos[iScn]['Area']['Salvage Logging']['Ensembles'][:,iEns]
        c=c+1; A['Management'][c]={}; A['Management'][c]['Name']='Thinning and knockdown'; A['Management'][c]['Color']=[0.2,0.4,0.7]; A['Management'][c]['Data']=mos[iScn]['Area']['Knockdown']['Ensembles'][:,iEns]+mos[iScn]['Area']['Thinning']['Ensembles'][:,iEns]
        c=c+1; A['Management'][c]={}; A['Management'][c]['Name']='Site preparation'; A['Management'][c]['Color']=[1,0.7,0.7]; A['Management'][c]['Data']=mos[iScn]['Area']['Disc Trenching']['Ensembles'][:,iEns]+mos[iScn]['Area']['Ripping']['Ensembles'][:,iEns]
        c=c+1; A['Management'][c]={}; A['Management'][c]['Name']='Prescribed burn'; A['Management'][c]['Color']=[0.5,0,0]; A['Management'][c]['Data']=mos[iScn]['Area']['Prescribed Burn']['Ensembles'][:,iEns]
        c=c+1; A['Management'][c]={}; A['Management'][c]['Name']='Dwarf Mistletoe control'; A['Management'][c]['Color']=[1,0.5,0]; A['Management'][c]['Data']=mos[iScn]['Area']['Dwarf Mistletoe Control']['Ensembles'][:,iEns]
        c=c+1; A['Management'][c]={}; A['Management'][c]['Name']='Planting'; A['Management'][c]['Color']=[0.3,0.8,0.2]; A['Management'][c]['Data']=mos[iScn]['Area']['Planting']['Ensembles'][:,iEns]+mos[iScn]['Area']['Direct Seeding']['Ensembles'][:,iEns]
        c=c+1; A['Management'][c]={}; A['Management'][c]['Name']='Foliage protection'; A['Management'][c]['Color']=[1,0.7,0]; A['Management'][c]['Data']=mos[iScn]['Area']['IDW Btk Spray']['Ensembles'][:,iEns]
        c=c+1; A['Management'][c]={}; A['Management'][c]['Name']='Aerial nutrient application'; A['Management'][c]['Color']=[0.65,0,1]; A['Management'][c]['Data']=mos[iScn]['Area']['Fertilization Aerial']['Ensembles'][:,iEns]
        A['Management']=A['Management'][0:c+1]
        
    else:
        # Ensemble mean
        A['Nat Dist']=[None]*10; c=-1
        c=c+1; A['Nat Dist'][c]={}; A['Nat Dist'][c]['Name']='Wildfire'; A['Nat Dist'][c]['Color']=[0.75,0,0]; A['Nat Dist'][c]['Data']=mos[iScn]['Area']['Wildfire']['Ensemble Mean']
        c=c+1; A['Nat Dist'][c]={}; A['Nat Dist'][c]['Name']='Mountain Pine beetle'; A['Nat Dist'][c]['Color']=[0,0.8,0]; A['Nat Dist'][c]['Data']=mos[iScn]['Area']['IBM']['Ensemble Mean']
        c=c+1; A['Nat Dist'][c]={}; A['Nat Dist'][c]['Name']='Douglas-fir beetle'; A['Nat Dist'][c]['Color']=[0.6,1,0]; A['Nat Dist'][c]['Data']=mos[iScn]['Area']['IBD']['Ensemble Mean']
        c=c+1; A['Nat Dist'][c]={}; A['Nat Dist'][c]['Name']='Spruce beetle'; A['Nat Dist'][c]['Color']=[0.25,1,1]; A['Nat Dist'][c]['Data']=mos[iScn]['Area']['IBS']['Ensemble Mean']
        c=c+1; A['Nat Dist'][c]={}; A['Nat Dist'][c]['Name']='W. balsam beetle'; A['Nat Dist'][c]['Color']=[0,0.45,0]; A['Nat Dist'][c]['Data']=mos[iScn]['Area']['IBB']['Ensemble Mean']
        c=c+1; A['Nat Dist'][c]={}; A['Nat Dist'][c]['Name']='Other pests'; A['Nat Dist'][c]['Color']=[0.8,1,0]; A['Nat Dist'][c]['Data']=mos[iScn]['Area']['Beetles']['Ensemble Mean']
        c=c+1; A['Nat Dist'][c]={}; A['Nat Dist'][c]['Name']='W. spruce budworm'; A['Nat Dist'][c]['Color']=[0,0.75,1]; A['Nat Dist'][c]['Data']=mos[iScn]['Area']['IDW']['Ensemble Mean']
        c=c+1; A['Nat Dist'][c]={}; A['Nat Dist'][c]['Name']='Rust'; A['Nat Dist'][c]['Color']=[0.75,0.5,1]; A['Nat Dist'][c]['Data']=mos[iScn]['Area']['Rust Onset']['Ensemble Mean']
        c=c+1; A['Nat Dist'][c]={}; A['Nat Dist'][c]['Name']='Dwarf Mistletoe'; A['Nat Dist'][c]['Color']=[1,0.5,0.25]; A['Nat Dist'][c]['Data']=mos[iScn]['Area']['Dwarf Mistletoe Onset']['Ensemble Mean']
        c=c+1; A['Nat Dist'][c]={}; A['Nat Dist'][c]['Name']='Mechanical'; A['Nat Dist'][c]['Color']=[0,0,0.6]; A['Nat Dist'][c]['Data']=mos[iScn]['Area']['Mechanical']['Ensemble Mean']
        A['Nat Dist']=A['Nat Dist'][0:c+1]

        A['Management']=[None]*10; c=-1
        c=c+1; A['Management'][c]={}; A['Management'][c]['Name']='Harvest'; A['Management'][c]['Color']=[0,0.75,1]; A['Management'][c]['Data']=mos[iScn]['Area']['Harvest']['Ensemble Mean']+mos[iScn]['Area']['Salvage Logging']['Ensemble Mean']
        c=c+1; A['Management'][c]={}; A['Management'][c]['Name']='Slashpile burn'; A['Management'][c]['Color']=[0.75,0,0]; A['Management'][c]['Data']=mos[iScn]['Area']['Slashpile Burn']['Ensemble Mean']+mos[iScn]['Area']['Salvage Logging']['Ensemble Mean']
        c=c+1; A['Management'][c]={}; A['Management'][c]['Name']='Thinning and knockdown'; A['Management'][c]['Color']=[0.2,0.4,0.7]; A['Management'][c]['Data']=mos[iScn]['Area']['Knockdown']['Ensemble Mean']+mos[iScn]['Area']['Thinning']['Ensemble Mean']
        c=c+1; A['Management'][c]={}; A['Management'][c]['Name']='Site preparation'; A['Management'][c]['Color']=[1,0.7,0.7]; A['Management'][c]['Data']=mos[iScn]['Area']['Disc Trenching']['Ensemble Mean']+mos[iScn]['Area']['Ripping']['Ensemble Mean']
        c=c+1; A['Management'][c]={}; A['Management'][c]['Name']='Prescribed burn'; A['Management'][c]['Color']=[0.5,0,0]; A['Management'][c]['Data']=mos[iScn]['Area']['Prescribed Burn']['Ensemble Mean']
        c=c+1; A['Management'][c]={}; A['Management'][c]['Name']='Dwarf Mistletoe control'; A['Management'][c]['Color']=[1,0.5,0]; A['Management'][c]['Data']=mos[iScn]['Area']['Dwarf Mistletoe Control']['Ensemble Mean']
        c=c+1; A['Management'][c]={}; A['Management'][c]['Name']='Planting'; A['Management'][c]['Color']=[0.3,0.8,0.2]; A['Management'][c]['Data']=mos[iScn]['Area']['Planting']['Ensemble Mean']+mos[iScn]['Area']['Direct Seeding']['Ensemble Mean']
        c=c+1; A['Management'][c]={}; A['Management'][c]['Name']='Foliage protection'; A['Management'][c]['Color']=[1,0.7,0]; A['Management'][c]['Data']=mos[iScn]['Area']['IDW Btk Spray']['Ensemble Mean']
        c=c+1; A['Management'][c]={}; A['Management'][c]['Name']='Aerial nutrient application'; A['Management'][c]['Color']=[0.65,0,1]; A['Management'][c]['Data']=mos[iScn]['Area']['Fertilization Aerial']['Ensemble Mean']
        A['Management']=A['Management'][0:c+1]

    # Apply area expansion factor
    for i in range(len(A['Nat Dist'])):
        A['Nat Dist'][i]['Data']=AEF*A['Nat Dist'][i]['Data']
    for i in range(len(A['Management'])):
        A['Management'][i]['Data']=AEF*A['Management'][i]['Data']    
        
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
    c=c+1; A['Management'][c]={}; A['Management'][c]['Name']='Harvest'; A['Management'][c]['Color']=[0,0.75,1]; A['Management'][c]['Data']=MosByMP[iScn]['Area']['Harvest']['Ensemble Mean'][:,iMP]+MosByMP[iScn]['Area']['Salvage Logging']['Ensemble Mean'][:,iMP]
    c=c+1; A['Management'][c]={}; A['Management'][c]['Name']='Slashpile burn'; A['Management'][c]['Color']=[0.75,0,0]; A['Management'][c]['Data']=MosByMP[iScn]['Area']['Slashpile Burn']['Ensemble Mean'][:,iMP]+MosByMP[iScn]['Area']['Salvage Logging']['Ensemble Mean'][:,iMP]
    c=c+1; A['Management'][c]={}; A['Management'][c]['Name']='Thinning and knockdown'; A['Management'][c]['Color']=[0.2,0.4,0.7]; A['Management'][c]['Data']=MosByMP[iScn]['Area']['Knockdown']['Ensemble Mean'][:,iMP]+MosByMP[iScn]['Area']['Thinning']['Ensemble Mean'][:,iMP]
    c=c+1; A['Management'][c]={}; A['Management'][c]['Name']='Site preparation'; A['Management'][c]['Color']=[1,0.7,0.7]; A['Management'][c]['Data']=MosByMP[iScn]['Area']['Disc Trenching']['Ensemble Mean'][:,iMP]+MosByMP[iScn]['Area']['Ripping']['Ensemble Mean'][:,iMP]
    c=c+1; A['Management'][c]={}; A['Management'][c]['Name']='Prescribed burn'; A['Management'][c]['Color']=[0.5,0,0]; A['Management'][c]['Data']=MosByMP[iScn]['Area']['Prescribed Burn']['Ensemble Mean'][:,iMP]
    c=c+1; A['Management'][c]={}; A['Management'][c]['Name']='Dwarf Mistletoe control'; A['Management'][c]['Color']=[1,0.5,0]; A['Management'][c]['Data']=MosByMP[iScn]['Area']['Dwarf Mistletoe Control']['Ensemble Mean'][:,iMP]
    c=c+1; A['Management'][c]={}; A['Management'][c]['Name']='Planting'; A['Management'][c]['Color']=[0.3,0.8,0.2]; A['Management'][c]['Data']=MosByMP[iScn]['Area']['Planting']['Ensemble Mean'][:,iMP]+MosByMP[iScn]['Area']['Direct Seeding']['Ensemble Mean'][:,iMP]
    c=c+1; A['Management'][c]={}; A['Management'][c]['Name']='Foliage protection'; A['Management'][c]['Color']=[1,0.7,0]; A['Management'][c]['Data']=MosByMP[iScn]['Area']['IDW Btk Spray']['Ensemble Mean'][:,iMP]
    c=c+1; A['Management'][c]={}; A['Management'][c]['Name']='Aerial nutrient application'; A['Management'][c]['Color']=[0.65,0,1]; A['Management'][c]['Data']=MosByMP[iScn]['Area']['Fertilization Aerial']['Ensemble Mean'][:,iMP]
    A['Management']=A['Management'][0:c+1]  
        
    # Convert to x-year intervals
    A['tv']=gu.BlockMean(tv,ivlT)
    for i in range(len(A['Nat Dist'])):
        A['Nat Dist'][i]['Data']=gu.BlockMean(A['Nat Dist'][i]['Data'],ivlT)
    for i in range(len(A['Management'])):
        A['Management'][i]['Data']=gu.BlockMean(A['Management'][i]['Data'],ivlT)        
    
    return A

#%% QA plot
    
def QA_Plot_ByMultiPolygon(meta,uMP,ivlMP,iScnForArea,ivlT,tv,it,MosByMP):
    
    for iMP in range(0,uMP.size,ivlMP):
    
        # Get area affected for multipolygon    
        A=cbu.AreaAffectedInSingleMultipolygon(meta,iScnForArea,ivlT,tv,MosByMP,iMP)
    
        #atu_multipolygons[uMP[iMP]]
    
        lw1=1; cle1=[0,0,1]; cle2=[1,0,0]; ms=3;
        xlim=[tv[it[0]],tv[it[-1]]]; xticks=np.arange(1500,2200,10);
    
        plt.close('all'); fig,ax=plt.subplots(4,2,figsize=gu.cm2inch(28,20))
   
        pl_d=[None]*len(A['Nat Dist']); nams_d=[None]*len(A['Nat Dist']);
        for i in range(len(A['Nat Dist'])):
            bottom=0; 
            if i!=0:
                for j in range(i):
                    bottom=bottom+A['Nat Dist'][j]['Data']
            pl_d[i]=ax[0,0].bar(A['tv'],A['Nat Dist'][i]['Data'],ivlT,color=A['Nat Dist'][i]['Color'],bottom=bottom)
            nams_d[i]=A['Nat Dist'][i]['Name']
        ax[0,0].legend(pl_d,nams_d,loc='upper left',bbox_to_anchor=(0.05,0.99),labelspacing=0.12,facecolor=[1,1,1],frameon=False)
        ax[0,0].set(position=[0.04,0.77,0.44,0.22],xlim=xlim,xticks=xticks,ylabel='Area affected (ha)');

        pl_m=[None]*len(A['Management']); nams_m=[None]*len(A['Management']);
        for i in range(len(A['Management'])):
            bottom=0; 
            if i!=0:
                for j in range(i):
                    bottom=bottom+A['Management'][j]['Data']
            pl_m[i]=ax[0,1].bar(A['tv'],A['Management'][i]['Data'],ivlT,color=A['Management'][i]['Color'],bottom=bottom)
            nams_m[i]=A['Management'][i]['Name']
        ax[0,1].legend(pl_m,nams_m,loc='upper left',bbox_to_anchor=(0.05,0.99),labelspacing=0.12,facecolor=[1,1,1],frameon=False)
        ax[0,1].set(position=[0.54,0.77,0.44,0.22],xlim=xlim,xlabel='',ylabel='Area affected (ha)');

        ax[1,0].plot(tv,MosByMP[iB]['v2']['Mean']['A']['Ensemble Mean'][:,iMP],'-',color=cle1,linewidth=lw1)
        ax[1,0].plot(tv,MosByMP[iP]['v2']['Mean']['A']['Ensemble Mean'][:,iMP],'--',color=cle2,linewidth=lw1)
        ax[1,0].set(position=[0.04,0.53,0.44,0.22],xlim=xlim,xlabel='',ylabel='Age, years')
    
        ax[1,1].plot(tv,MosByMP[iB]['v2']['Mean']['Eco_Biomass']['Ensemble Mean'][:,iMP],'-',color=cle1,linewidth=lw1)
        ax[1,1].plot(tv,MosByMP[iP]['v2']['Mean']['Eco_Biomass']['Ensemble Mean'][:,iMP],'--',color=cle2,linewidth=lw1)
        ax[1,1].set(position=[0.54,0.53,0.44,0.22],xlim=xlim,xlabel='',ylabel='Biomass (MgC/ha)')
    
        ax[2,0].plot(tv,MosByMP[iB]['v2']['Mean']['Eco_DeadWood']['Ensemble Mean'][:,iMP],'-',color=cle1,linewidth=lw1)
        ax[2,0].plot(tv,MosByMP[iP]['v2']['Mean']['Eco_DeadWood']['Ensemble Mean'][:,iMP],'--',color=cle2,linewidth=lw1)
        ax[2,0].set(position=[0.04,0.285,0.44,0.22],xlim=xlim,xlabel='',ylabel='Dead wood (MgC/ha)')
    
        ax[2,1].plot(tv,MosByMP[iB]['v2']['Mean']['Eco_E_Wildfire']['Ensemble Mean'][:,iMP],'o-',color=cle1,linewidth=lw1,markersize=ms)
        ax[2,1].plot(tv,MosByMP[iP]['v2']['Mean']['Eco_E_Wildfire']['Ensemble Mean'][:,iMP],'s--',color=cle2,linewidth=lw1,markersize=ms)
        ax[2,1].set(position=[0.54,0.285,0.44,0.22],xlim=xlim,xlabel='',ylabel='Wildfire emissions (MgC/ha)')
    
        ax[3,0].plot(tv,MosByMP[iB]['v2']['Mean']['Eco_Removals']['Ensemble Mean'][:,iMP],'o-',color=cle1,linewidth=lw1,markersize=ms)
        ax[3,0].plot(tv,MosByMP[iP]['v2']['Mean']['Eco_Removals']['Ensemble Mean'][:,iMP],'s--',color=cle2,linewidth=lw1,markersize=ms)
        ax[3,0].set(position=[0.04,0.04,0.44,0.22],xlim=xlim,xlabel='',ylabel='Removals (MgC/ha)')
    
        ax[3,1].plot(tv,MosByMP[iB]['v2']['Mean']['Sec_NGHGB']['Ensemble Mean'][:,iMP],'-',color=cle1,linewidth=lw1)
        ax[3,1].plot(tv,MosByMP[iP]['v2']['Mean']['Sec_NGHGB']['Ensemble Mean'][:,iMP],'--',color=cle2,linewidth=lw1)
        ax[3,1].set(position=[0.54,0.04,0.44,0.22],xlim=xlim,xlabel='',ylabel='Sector GHG balance (MgC/ha)')
    
        gu.axletters(ax,plt,0.01,0.91)
        pt=meta['LUT']['ProjectType'][meta['ProjectType'][iMP]]
        gu.PrintFig(meta['Paths']['Figures'] + '\\BySparseGridSample\\MP' + str(iMP) + '_' + pt,'png',200)
    
    return
