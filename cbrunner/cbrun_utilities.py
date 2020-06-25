
import numpy as np
import pandas as pd
import pickle
import os
import matplotlib.pyplot as plt
from fcgadgets.pyscripts import utilities_general as gu

'''============================================================================
CONVERT LUT NUMBER TO STRING NAME
============================================================================'''

def lut_n2s(dc,numb):
    if numb!=-999:
        vals=np.fromiter(dc.values(),dtype=float)
        keys=np.fromiter(dc.keys(),dtype='<U30')
        ind=np.where(vals==numb)[0]
        s=keys[ind]
    else:
        s=np.array(['Unidentified'],ndmin=1)
    return s

'''============================================================================
QUERY RESULTS CODES FOR MANAGEMENT ACTIVITY TYPES
============================================================================'''

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
        
        elif (d['SILV_BASE_CODE'][i]=='SP') & (d['SILV_TECHNIQUE_CODE'][i]=='BU') & (d['SILV_METHOD_CODE'][i]=='PBURN') | \
            (d['SILV_BASE_CODE'][i]=='SP') & (d['SILV_TECHNIQUE_CODE'][i]=='ME') & (d['SILV_METHOD_CODE'][i]=='PBURN'):
            Name.append('Slashpile Burn')
            
        elif (d['SILV_BASE_CODE'][i]=='SP') & (d['SILV_TECHNIQUE_CODE'][i]=='MA') & (d['SILV_METHOD_CODE'][i]=='SNAG') | \
            (d['SILV_BASE_CODE'][i]=='SP') & (d['SILV_TECHNIQUE_CODE'][i]=='ME') & (d['SILV_METHOD_CODE'][i]=='MDOWN') | \
            (d['SILV_BASE_CODE'][i]=='SP') & (d['SILV_TECHNIQUE_CODE'][i]=='MA') & (d['SILV_METHOD_CODE'][i]=='HAND') | \
            (d['SILV_BASE_CODE'][i]=='SP') & (d['SILV_TECHNIQUE_CODE'][i]=='MA') & (d['SILV_METHOD_CODE'][i]=='MANCT'):
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


'''============================================================================
CONVERT DICTIONARY TO DATA STRUCTURE CLASS
============================================================================'''

class BunchDictionary(dict):
    def __init__(self, *args, **kwds):
        super(BunchDictionary, self).__init__(*args, **kwds)
        self.__dict__ = self

'''============================================================================
CONSTRUCT DISTURBANCES FOR SPINUP
============================================================================'''

def CompileDisturbancesForSpinup(meta,iScn,eh):

    ivl_spin=meta['Spinup Disturbance Return Inverval']
    yr1=meta['Year Start']+ivl_spin
    yr2=meta['Spinup Year End']
                
    YearRef=meta['Scenario'][iScn]['Year1_DisFromInv']
    AgeRef=meta['Scenario'][iScn]['Age1_DisFromInv']                
    if AgeRef>=0:
        Year=np.arange(YearRef-AgeRef-100*ivl_spin,YearRef-AgeRef+ivl_spin,ivl_spin)        
    else:
        Year=np.arange(yr1,yr2+1,meta['Spinup Disturbance Return Inverval'])
        
    eh['Year']=Year[np.where(Year>=meta['Time'][0])[0]]               
                
    m_Dist=eh['Year'].size
               
    eh['ID_Type']=meta['LUT Dist'][meta['Spinup Disturbance Type']]*np.ones(m_Dist)
                
    eh['MortalityFactor']=100*np.ones(m_Dist)
    
    eh['GrowthFactor']=0*np.ones(m_Dist)
                
    if (meta['Biomass Module']=='TIPSY') | (meta['Biomass Module']=='TIPSY_EP703'):
        eh['ID_GrowthCurve']=meta['Spinup Growth Curve ID']*np.ones(m_Dist)
    else:
        eh['ID_GrowthCurve']=-999*np.ones(m_Dist)
   
    return eh

'''============================================================================
CONSTRUCT DISTURBANCES FROM INVENTORY
============================================================================'''

def CompileDisturbancesFromInventory(meta,iScn,eh):
    
    for iYr in range(1,7):
        if np.isnan(meta['Scenario'][iScn]['Year' + str(iYr) + '_DisFromInv'])==False:  
            
            # If IDW, convert IDW class to growth and mortality factor
            sc=np.array(['IDW-T','IDW-L','IDW-M','IDW-S','IDW-V','IDW-MM','IDW-MS','IDW-MV','IDW-SS','IDW-SV','IDW-VV'])
            flg_i=0
            indSc=np.where(sc==meta['Scenario'][iScn]['Type' + str(iYr) + '_DisFromInv'])[0]
            if indSc.size!=0:                
                if flg_i==0:
                    dfParDistBySC=pd.read_excel(meta['Path Model Code'] + '\\Parameters\\Parameters_DisturbanceBySeverityClass.xlsx')
                    flg_i=1
                indPar=np.where( (dfParDistBySC['Name']=='IDW') & (dfParDistBySC['SeverityCD']==sc[indSc[0]][4:]) )[0]
                ID_TypeN=meta['LUT Dist']['IDW']
                MF=dfParDistBySC.loc[indPar,'MortalityFactor']
                GF=dfParDistBySC.loc[indPar,'GrowthFactor']
            else:
                ID_TypeS=meta['Scenario'][iScn]['Type' + str(iYr) + '_DisFromInv']
                ID_TypeN=meta['LUT Dist'][ID_TypeS]
                MF=meta['Scenario'][iScn]['Severity' + str(iYr) + '_DisFromInv']
                GF=0
                
            # Populate event history dictionary
            eh['Year']=np.append(eh['Year'],meta['Scenario'][iScn]['Year' + str(iYr) + '_DisFromInv'])
            eh['ID_Type']=np.append(eh['ID_Type'],ID_TypeN)
            eh['MortalityFactor']=np.append(eh['MortalityFactor'],MF)
            eh['GrowthFactor']=np.append(eh['GrowthFactor'],GF)
            if (meta['Biomass Module']=='TIPSY') | (meta['Biomass Module']=='TIPSY_EP703'):
                eh['ID_GrowthCurve']=np.append(eh['ID_GrowthCurve'],meta['Scenario'][iScn]['GrowthCurve' + str(iYr) + '_DisFromInv'])
            else:
                eh['ID_GrowthCurve']=np.append(eh['ID_GrowthCurve'],-999)
                
    return eh

'''============================================================================
CONSTRUCT DISTURBANCES FROM SIMULATION
============================================================================'''

def CompileHistoricalDisturbancesFromSimulation(meta,iScn,eh):

    # Historical disturbance from simulation 1
    ri=meta['Scenario'][iScn]['ReturnInterval1_Hist_DisFromSim']
    if (ri!=0) & (np.isnan(ri)!=True):
        p_Dist=1/ri
        p_Rand=np.random.uniform(0,1,size=(meta['Time'].size))        
        it=np.where((p_Rand<p_Dist) & (meta['Time']>meta['Spinup Year End']) & (meta['Time']<meta['Year Project']))[0]
        eh['Year']=np.append(eh['Year'],meta['Time'][it])
        ID_Type0=meta['LUT Dist'][meta['Scenario'][iScn]['Type1_Hist_DisFromSim']]
        eh['ID_Type']=np.append(eh['ID_Type'],ID_Type0*np.ones(it.size))
        MortalityFactor0=100
        eh['MortalityFactor']=np.append(eh['MortalityFactor'],MortalityFactor0*np.ones(it.size))
        eh['GrowthFactor']=np.append(eh['GrowthFactor'],0*np.ones(it.size))
        eh['ID_GrowthCurve']=np.append(eh['ID_GrowthCurve'],eh['ID_GrowthCurve'][-1]*np.ones(it.size))
    
    # Historical disturbance from simulation 2
    ri=meta['Scenario'][iScn]['ReturnInterval2_Hist_DisFromSim']
    if (ri!=0) & (np.isnan(ri)!=True):
        p_Dist=1/ri
        p_Rand=np.random.uniform(0,1,size=(meta['Time'].size))        
        it=np.where((p_Rand<p_Dist) & (meta['Time']>meta['Spinup Year End']) & (meta['Time']<meta['Year Project']))[0]
        eh['Year']=np.append(eh['Year'],meta['Time'][it])
        ID_Type0=meta['LUT Dist'][meta['Scenario'][iScn]['Type2_Hist_DisFromSim']]
        eh['ID_Type']=np.append(eh['ID_Type'],ID_Type0*np.ones(it.size))
        MortalityFactor0=100
        eh['MortalityFactor']=np.append(eh['MortalityFactor'],MortalityFactor0*np.ones(it.size))
        eh['GrowthFactor']=np.append(eh['GrowthFactor'],0*np.ones(it.size))
        eh['ID_GrowthCurve']=np.append(eh['ID_GrowthCurve'],eh['ID_GrowthCurve'][-1]*np.ones(it.size))

    return eh

def CompileFutureDisturbancesFromSimulation(meta,iScn,eh):

    # Future disturbance from simulation 1
    ri=meta['Scenario'][iScn]['ReturnInterval1_Fut_DisFromSim']
    if (ri!=0) & (np.isnan(ri)!=True):
        p_Dist=1/ri
        p_Rand=np.random.uniform(0,1,size=(meta['Time'].size))        
        it=np.where((p_Rand<p_Dist) & (meta['Time']>meta['Year Project']))[0]
        eh['Year']=np.append(eh['Year'],meta['Time'][it])
        ID_Type0=meta['LUT Dist'][meta['Scenario'][iScn]['Type1_Fut_DisFromSim']]
        eh['ID_Type']=np.append(eh['ID_Type'],ID_Type0*np.ones(it.size))
        MortalityFactor0=100
        eh['MortalityFactor']=np.append(eh['MortalityFactor'],MortalityFactor0*np.ones(it.size))
        eh['GrowthFactor']=np.append(eh['GrowthFactor'],0*np.ones(it.size)) 
        eh['ID_GrowthCurve']=np.append(eh['ID_GrowthCurve'],eh['ID_GrowthCurve'][-1]*np.ones(it.size))             
    
    # Future disturbance from simulation 2
    ri=meta['Scenario'][iScn]['ReturnInterval2_Fut_DisFromSim']
    if (ri!=0) & (np.isnan(ri)!=True):
        p_Dist=1/ri
        p_Rand=np.random.uniform(0,1,size=(meta['Time'].size))        
        it=np.where((p_Rand<p_Dist) & (meta['Time']>meta['Year Project']))[0]
        eh['Year']=np.append(eh['Year'],meta['Time'][it])
        ID_Type0=meta['LUT Dist'][meta['Scenario'][iScn]['Type2_Fut_DisFromSim']]
        eh['ID_Type']=np.append(eh['ID_Type'],ID_Type0*np.ones(it.size))
        MortalityFactor0=100
        eh['MortalityFactor']=np.append(eh['MortalityFactor'],MortalityFactor0*np.ones(it.size))
        eh['GrowthFactor']=np.append(eh['GrowthFactor'],0*np.ones(it.size)) 
        eh['ID_GrowthCurve']=np.append(eh['ID_GrowthCurve'],eh['ID_GrowthCurve'][-1]*np.ones(it.size))
    
    return eh

'''============================================================================
FIX ENSEMBLE NAME NUMBERING
============================================================================'''

def FixFileNum(ind):
    indStrFixed=str(ind+1)
    if len(indStrFixed)==1:
        indStrFixed='000' + indStrFixed
    elif len(indStrFixed)==2:
        indStrFixed='00' + indStrFixed    
    elif len(indStrFixed)==3:
        indStrFixed='0' + indStrFixed
    return indStrFixed

'''============================================================================
IMPORT LOOKUP TABLES
============================================================================'''

def ImportLUTs(pthin):
    
    # Open connection to parameter database    
    par=gu.ipickle(pthin + '\\Parameters\\Parameters.pkl')
    
    # Import distubance type        
    LUT_Dist={}
    for i in range(len(par['Disturbances']['Name'])):
        LUT_Dist[par['Disturbances']['Name'][i]]=par['Disturbances']['ID'][i]
    
    # BGC zone     
    LUT_BGC_Zone={}
    for i in range(len(par['BGC_ZONE']['CODE_BGC_ZONE'])):
        LUT_BGC_Zone[par['BGC_ZONE']['CODE_BGC_ZONE'][i]]=par['BGC_ZONE']['ID_BGC_ZONE'][i]
    
    # Species
    LUT_Spc={}
    for i in range(len(par['SRS']['SRS_CD'])):
        LUT_Spc[par['SRS']['SRS_CD'][i]]=par['SRS']['SRS_ID'][i]
    
    return LUT_Dist,LUT_Spc,LUT_BGC_Zone

'''============================================================================
IMPORT PROJECT CONFIGURATION
============================================================================'''

def ImportProjectConfig(meta):
    
    # Import project parameters from spreadsheet
    df_p=pd.read_excel(meta['Path Project'] + '\\Inputs\\ProjectConfig.xlsx',sheet_name='Inputs',skiprows=1,usecols='A:B')
    for i in range(df_p.shape[0]): 
        meta.update({df_p.iloc[i,0]:df_p.iat[i,1]})    
    
    # Import look-up tables
    LUT_Dist,LUT_Spc,LUT_BGC_Zone=ImportLUTs(meta['Path Model Code'])
    meta['LUT Dist']=LUT_Dist
    meta['LUT Spc']=LUT_Spc
    meta['LUT BGC Zone']=LUT_BGC_Zone

    # Pool names (ecosystem)
    # *** If you change this, you need to change the same list in "UpdateParameters" ***
    meta['Name Pools Eco']=['StemMerch','StemNonMerch','Foliage','Branch','Bark','RootCoarse','RootFine', \
        'FelledStemMerch','FelledStemNonMerch','FelledBranch','FelledBark','FelledSnagStem','FelledSnagBranch', \
        'LitterVF','LitterF','LitterM','LitterS','SnagStem','SnagBranch','SoilVF','SoilF','SoilS', \
        'ECO2asC','ECH4asC','ECOasC','EN2OasC','RemovedMerch','RemovedNonMerch','RemovedSnagStem']
    meta['N Pools Eco']=len(meta['Name Pools Eco'])
    
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
    
    # Pool names (products)
    meta['Name Pools Pro']=['SFH','MFH','Comm','Furn','Ship','Repairs','Other','Paper','Fuel','Firewood','EffluentPulp', \
        'DumpWood','DumpPaper','LandfillWoodDegradable','LandfillWoodNonDegradable', \
        'LandfillPaperDegradable','LandfillPaperNonDegradable','E_CO2','E_CH4','Cants']    
    meta['N Pools Pro']=len(meta['Name Pools Pro'])
    
    # Indices to produce pools pools
    meta['iPP']={}; cnt=0
    for nam in meta['Name Pools Pro']:
        meta['iPP'][nam]=cnt
        cnt=cnt+1
    iPP=meta['iPP']
    
    # Custom disturbance variable names
    meta['Name CustDistVar']=[
            'Biomass_Affected_Pct','Biomass_Merch_Removed_Pct','Biomass_Merch_Burned_Pct','Biomass_Merch_LeftOnSite_Pct', \
            'Biomass_NonMerch_Removed_Pct','Biomass_NonMerch_Burned_Pct','Biomass_NonMerch_LeftOnSite_Pct', \
            'Snags_Affected_Pct','Snags_Removed_Pct','Snags_Burned_Pct','Snags_LeftOnSite_Pct', \
            'RemovedMerchToFuel_Pct','RemovedMerchToLumber_Pct','RemovedMerchToPlywood_Pct', \
            'RemovedMerchToOSB_Pct','RemovedMerchToMDF_Pct','RemovedMerchToPulp_Pct', \
            'RemovedMerchToFirewood_Pct','RemovedMerchToCants_Pct','RemovedNonMerchToFuel_Pct','RemovedNonMerchToLumber_Pct', \
            'RemovedNonMerchToPlywood_Pct','RemovedNonMerchToOSB_Pct','RemovedNonMerchToMDF_Pct', \
            'RemovedNonMerchToPulp_Pct','RemovedNonMerchToFirewood_Pct','RemovedNonMerchToCants_Pct','RemovedSnagStemToFuel_Pct', \
            'RemovedSnagStemToLumber_Pct','RemovedSnagStemToPlywood_Pct','RemovedSnagStemToOSB_Pct', \
            'RemovedSnagStemToMDF_Pct','RemovedSnagStemToPulp_Pct','RemovedSnagStemToFirewood_Pct','RemovedSnagStemToCants_Pct']   
    
    # Calendar year
    meta['Time']=np.arange(meta['Year Start'],meta['Year End']+1,1)
    
    # Number of stands 
    if meta['Scenario Source']=='Spreadsheet' and meta['N Ensemble']==1:
        meta['N Stand']=1
    elif meta['Scenario Source']=='Spreadsheet' and meta['N Ensemble']>1:
        # If running ensembles from spreadsheet, it is faster to run them as stands
        meta['N Stand']=meta['N Ensemble']
        meta['N Ensemble']=1
    
    # Define scenario parameters
    if meta['Scenario Source']=='Spreadsheet':
        
        # Import from spreadsheet
        df_s=pd.read_excel(meta['Path Project'] + '\\Inputs\\ProjectConfig.xlsx',sheet_name='Inputs',skiprows=1,usecols='D:AN')
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
            
    elif (meta['Scenario Source']=='Script') & (meta['Biomass Module']=='TIPSY'):
        
        # Initialize scenario data structure
        meta['Scenario']=list()
        for i in range(meta['N Scenario']):
            pScn0={}
            meta['Scenario'].append(pScn0)

    # Number of batches
    meta['N Batch']=np.ceil(meta['N Stand']/meta['Batch Interval']).astype(int)

    # Projects with more than one batch change N Stand - create a full value
    meta['N Stand Full']=meta['N Stand']

    # Make all the folders
    meta['Path Input Scenario']=[]
    meta['Path Output Scenario']=[]
    for iScn in range(0,meta['N Scenario']):    
        meta['Path Input Scenario'].append(meta['Path Project'] + '\\Inputs\\Scenario' + FixFileNum(iScn))
        if os.path.exists(meta['Path Input Scenario'][iScn])==False:
            os.mkdir(meta['Path Input Scenario'][iScn])
        meta['Path Output Scenario'].append(meta['Path Project'] + '\\Outputs\\Scenario' + FixFileNum(iScn))
        if os.path.exists(meta['Path Output Scenario'][iScn])==False:
            os.mkdir(meta['Path Output Scenario'][iScn])
   
    # Scale factor for growth curves
    # Note: Do not change this to 10 - aerial fertilization response will not work properly at 10
    meta['Scale Factor GC']=100
    
    # Scale factor for saving results (this needs to be 100, 10 does not capture HWP fluxes)
    meta['Scale Factor Export']=1000
    
    # Default status of growth factors
    for iScn in range(0,meta['N Scenario']): 
        meta['Scenario'][iScn]['Status Net Growth Factor']='Off'
        meta['Scenario'][iScn]['Status Mortality Factor']='Off'
    
    # Initialize flag indicating whether a disturbance occurred
    meta['Flag_Disturbed']=np.zeros(meta['N Stand'])
    
    return meta

'''============================================================================
LOAD SCENARIO RUSULTS
Return a list of dictionaries for each scenario. If multiple ensemble were run, 
the function will retun the average.
============================================================================'''

def LoadScenarioResults(meta,scn):

    # Open connection to parameter database   
    meta=gu.ipickle(meta['Path Project'] + '\\Outputs\\Scenario' + FixFileNum(0) + '\\Metadata.pkl')
    
    # Initialize list that will contain scenarios
    v=[]
    for iScn in scn:
        for iEns in range(0,meta['N Ensemble']):            
            for iBat in range(0,meta['N Batch']):
                
                # Open results 
                fin=open(meta['Path Project'] + '\\Outputs\\Scenario' + FixFileNum(iScn) + '\\Data_Scn' + FixFileNum(iScn) + '_Ens' + FixFileNum(iEns) + '_Bat' + FixFileNum(iBat) + '.pkl','rb')
                data_bat=pickle.load(fin)
                fin.close()
                
                # Convert to float and apply scale factor
                for k in data_bat.keys():
                    if k!='ScaleFactor':                        
                        data_bat[k]=data_bat[k].astype(float)
                        data_bat[k]=data_bat[k]/data_bat['ScaleFactor']
                del data_bat['ScaleFactor']
                
                # Accumulate data in each batch
                key=list(data_bat.keys())
                if iBat==0:
                    data_all=data_bat
                else:                    
                    for key in data_bat.keys():
                        if key=='Year':
                            continue
                        data_all[key]=np.append(data_all[key],data_bat[key],axis=1)
            
            # Sum across ensembles
            if iEns==0:
                data_sum2ave=data_all
            else:                    
                for key in data_bat.keys():
                    data_sum2ave[key]=data_sum2ave[key]+data_all[key]
        
        # If the simulation includes ensembles, calculate average
        for key in data_bat.keys():
            data_sum2ave[key]=data_sum2ave[key]/meta['N Ensemble']        
        
        v.append(BunchDictionary(data_sum2ave))
        
    return v,meta

'''============================================================================
LOAD SINGLE OUTPUT FILE FOR SCENARIO A, ENSEMBLE, B AND BATCH C
============================================================================'''

def LoadSingleOutputFile(meta,iScn,iEns,iBat):
       
    # Open output data
    data=gu.ipickle(meta['Path Project'] + '\\Outputs\\Scenario' + FixFileNum(iScn) + '\\Data_Scn' + FixFileNum(iScn) + '_Ens' + FixFileNum(iEns) + '_Bat' + FixFileNum(iBat) + '.pkl')
                
    # Convert to float and apply scale factor
    for k in data.keys():
        if k!='ScaleFactor':                        
            data[k]=data[k].astype(float)
            data[k]=data[k]/data['ScaleFactor']
    del data['ScaleFactor']
    
    key=list(data.keys())
    
    if iEns==0:
        data_sum2ave=data
    else:                    
        for k in range(len(key)):
            data_sum2ave[key[k]]=data_sum2ave[key[k]]+data[key[k]]
    for k in range(len(key)):
        data_sum2ave[key[k]]=data_sum2ave[key[k]]/meta['N Ensemble']        
    
    v=[BunchDictionary(data_sum2ave)]
        
    return v


'''============================================================================
POST-PROCESS TIPSY GROWTH CURVES
============================================================================'''

def PostProcessTIPSY(meta):

    # TIPSY exports curves as MgDM/ha/yr, CBRunner expects inputs of MgC/ha/yr - create conversion factor
    dm2c=0.5

    for iScn in range(meta['N Scenario']):

        # Growth curve parameters and TIPSY outputs
        if meta['Growth Curves Lumped']=='Yes':
            dfPar=pd.read_excel(meta['Path Project'] + '\\Inputs\\GrowthCurvesTIPSY_Parameters.xlsx',sheet_name='Sheet1',skiprows=6)
            txtDat=np.loadtxt(meta['Path Project'] + '\\Inputs\\GrowthCurvesTIPSY_Output.out',skiprows=4)
        else:
            dfPar=pd.read_excel(meta['Path Project'] + '\\Inputs\\Scenario' + FixFileNum(iScn) + '\\GrowthCurvesTIPSY_Parameters.xlsx',sheet_name='Sheet1',skiprows=6)
            txtDat=np.loadtxt(meta['Path Project'] + '\\Inputs\\Scenario' + FixFileNum(iScn) + '\\GrowthCurvesTIPSY_Output.out',skiprows=4)
        
        # TIPSY saves to text file -> convert to dataframe (column names must match TIPSY output file design)
        #colnams=['Age','VolTot0','VolMerch75','VolMerch125','VolMerch175','ODT_Bark','ODT_Branch','ODT_Foliage','ODT_Roots','ODT_Stem','Lum1_2x4','Lum1_2x6','Lum1_2x8','Lum1_2x10','Lum2orBetter_2x4','Lum2orBetter_2x6','Lum2orBetter_2x8','Lum2orBetter_2x10']
        colnams=['Age','VolTot0','VolMerch125','VolMerch175','ODT_Bark','ODT_Branch','ODT_Foliage','ODT_Roots','ODT_Stem','MortalityVolumeTotal']
        dfDat=pd.DataFrame(txtDat,columns=colnams)

        # Define age vector (must be consistent with how TIPSY was set up)
        Age=np.arange(0,301,1)

        # Get dimensions of the TIPSY output file to reshape the data into Age x Stand
        N_Age=Age.size
        N_GC=int(dfDat.shape[0]/N_Age)
    
        # Define the fraction of merchantable stemwood
        fMerch=np.nan_to_num(np.reshape(dfDat['VolMerch125'].values,(N_Age,N_GC),order='F')/np.reshape(dfDat['VolTot0'].values,(N_Age,N_GC),order='F'))
        fNonMerch=1-fMerch

        # Merchantable stemwood volume
        V_StemMerch=np.reshape(dfDat['VolMerch125'].values,(N_Age,N_GC),order='F')
        G_VStemMerch=np.append(np.zeros((1,N_GC)),np.diff(V_StemMerch,axis=0),axis=0)

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
    
        # Define growth curves (from TIPSY output)
        for iGC in range(3):
    
            for iBat in range(0,meta['N Batch']):
            
                # Index to batch
                iStart=meta['Batch Interval']*iBat
                iStop=np.minimum(meta['N Stand'],iStart+meta['Batch Interval'])
                indBat=np.arange(iStart,iStop,1)
                
                # Import disturbance history
                dh=gu.ipickle(meta['Path Input Scenario'][iScn] + '\\Disturbance_Ens' + FixFileNum(0) + '_Bat' + FixFileNum(iBat) + '.pkl')
                
                # Initialize age response of net growth
                G=np.zeros((N_Age,indBat.size,6),dtype=np.int16)
    
                # Populate the growth curve
                for iS in range(len(dh)):
        
                    ID_GrowthCurveM=1*np.ones((dh[iS]['ID_GrowthCurve'].shape[0],))
                    d=np.diff(dh[iS]['ID_GrowthCurve'])
                    cnt=1
                    for j in range(0,d.shape[0]):
                        if d[j]!=0:
                            cnt=cnt+1
                            ID_GrowthCurveM[j+1]=cnt 
                        else:
                            ID_GrowthCurveM[j+1]=ID_GrowthCurveM[j]        
                    u=np.unique(ID_GrowthCurveM)
        
                    if iGC+1>u.size:
                        continue
        
                    indDH=np.where(ID_GrowthCurveM==u[iGC])[0]
                    if indDH.size==0:
                        continue
                    
                    if meta['Scenario Source']=='Spreadsheet':
                        
                        indTIPSY=np.where(
                                (dfPar['ID_Scenario']==iScn+1) &
                                (dfPar['ID_GC']==int(dh[iS]['ID_GrowthCurve'][indDH][0])))[0]                    
                    
                    elif meta['Scenario Source']=='Script':                        
                        
                        if meta['Growth Curves Lumped']=='Yes':
                            indTIPSY=np.where(
                                (dfPar['ID_Stand']==indBat[iS]+1) & 
                                (dfPar['ID_Scenario']==iScn+1) &
                                (dfPar['ID_GC']==int(iGC+1)))[0]
                            # dh[iS].ID_GrowthCurve[indDH][0]
                        else:
                            indTIPSY=np.where(
                                (dfPar['ID_Stand']==indBat[iS]+1) &
                                (dfPar['ID_GC']==int(dh[iS]['ID_GrowthCurve'][indDH][0])))[0]                    
                    
                    if indTIPSY.size==0:
                        # This can happen if only some stands have a third GC, for example
                        continue
                    
                    G[:,iS,0]=meta['Scale Factor GC']*G_StemMerch[:,indTIPSY[0]]
                    G[:,iS,1]=meta['Scale Factor GC']*G_StemNonMerch[:,indTIPSY[0]]
                    G[:,iS,2]=meta['Scale Factor GC']*G_Bark[:,indTIPSY[0]]
                    G[:,iS,3]=meta['Scale Factor GC']*G_Branch[:,indTIPSY[0]]
                    G[:,iS,4]=meta['Scale Factor GC']*G_Foliage[:,indTIPSY[0]]
                    G[:,iS,5]=meta['Scale Factor GC']*G_VStemMerch[:,indTIPSY[0]]                    
                    
                # Save data to file in input variables folder of project
                gu.opickle(meta['Path Project'] + '\\Inputs\\Scenario' + FixFileNum(iScn) + '\\GrowthCurve' + str(iGC+1) + '_Bat' + FixFileNum(iBat) + '.pkl',G)
                    
    return


'''============================================================================
POST-PROCESS TIPSY GROWTH CURVES
Nested list, gc[Scenario][Stand][Growth Curve]
============================================================================'''

def Import_BatchTIPSY_Output(meta):
        
    # Growth curve parameters and TIPSY outputs
    if meta['Growth Curves Lumped']!='Yes':
        print('This wont work')
    dfPar=pd.read_excel(meta['Path Project'] + '\\Inputs\\GrowthCurvesTIPSY_Parameters.xlsx',sheet_name='Sheet1',skiprows=6)
    txtDat=np.loadtxt(meta['Path Project'] + '\\Inputs\\GrowthCurvesTIPSY_Output.out',skiprows=4)
    
    # TIPSY saves to text file -> convert to dataframe (column names must match TIPSY output file design)
    #colnams=['Age','VolTot0','VolMerch75','VolMerch125','VolMerch175','ODT_Bark','ODT_Branch','ODT_Foliage','ODT_Roots','ODT_Stem','Lum1_2x4','Lum1_2x6','Lum1_2x8','Lum1_2x10','Lum2orBetter_2x4','Lum2orBetter_2x6','Lum2orBetter_2x8','Lum2orBetter_2x10']
    colnams=['Age','VolTot0','VolMerch125','VolMerch175','ODT_Bark','ODT_Branch','ODT_Foliage','ODT_Roots','ODT_Stem','MortalityVolumeTotal']
    dfDat=pd.DataFrame(txtDat,columns=colnams)

    # Define age vector (must be consistent with how TIPSY was set up)
    Age=np.arange(0,301,1)

    # Get dimensions of the TIPSY output file to reshape the data into Age x Stand
    N_Age=Age.size
    N_GC=int(dfDat.shape[0]/N_Age)

    gc=[None]*N_GC
    for i in range(N_GC): gc[i]={}    
    for i in range(len(colnams)):
        data=np.reshape(dfDat[colnams[i]].values,(N_Age,N_GC),order='F')
        for j in range(N_GC):
            gc[j][colnams[i]]=data[:,j]

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
                for i in range(len(colnams)):
                    data=np.reshape(dfDat[colnams[i]].values,(N_Age,N_GC),order='F')
                    d[colnams[i]]=data[:,ind]
                gc0.append(d)
            gc1.append(gc0)
        gc2.append(gc1)
    return gc2

def Import_CompiledGrowthCurves(meta,scn):
    # *** This needs to be fixed - Vstem is on the end column ***
    # gc[iScn][iGC][iStand]
    gc=[]
    for iScn in range(len(scn)): #range(meta['N Scenario']):
        gc0=[]
        gc1=[]
        for iBat in range(0,meta['N Batch']):            
            tmp=gu.ipickle(meta['Path Project'] + '\\Inputs\\Scenario' + FixFileNum(iScn) + '\\GrowthCurve1_Bat' + FixFileNum(iBat) + '.pkl')
            tmp=np.sum(tmp,axis=2).astype(float)
            for iS in range(tmp.shape[1]):
                gc1.append(tmp[:,iS].copy()/meta['Scale Factor GC'])
        gc0.append(gc1.copy())
        gc1=[]
        for iBat in range(0,meta['N Batch']):            
            tmp=gu.ipickle(meta['Path Project'] + '\\Inputs\\Scenario' + FixFileNum(iScn) + '\\GrowthCurve2_Bat' + FixFileNum(iBat) + '.pkl')
            tmp=np.sum(tmp,axis=2).astype(float)
            for iS in range(tmp.shape[1]):
                gc1.append(tmp[:,iS].copy()/meta['Scale Factor GC'])
        gc0.append(gc1.copy())
        gc1=[]
        for iBat in range(0,meta['N Batch']):            
            tmp=gu.ipickle(meta['Path Project'] + '\\Inputs\\Scenario' + FixFileNum(iScn) + '\\GrowthCurve3_Bat' + FixFileNum(iBat) + '.pkl')
            tmp=np.sum(tmp,axis=2).astype(float)
            for iS in range(tmp.shape[1]):
                gc1.append(tmp[:,iS].copy()/meta['Scale Factor GC'])
        gc0.append(gc1.copy())
        gc.append(gc0.copy())
    return gc

'''============================================================================
IMPORT CUSTOM HARVEST ASSUMPTIONS
============================================================================'''

def ImportCustomHarvestAssumptions(pthin):
    
    # Initialize dictionary
    d={}

    # Import table
    df=pd.read_excel(pthin,sheet_name='Sheet1',skiprows=1)
    
    d['Biomass_Affected_Pct']=df.iloc[1,1]    
    d['Biomass_Merch_Removed_Pct']=df.iloc[5,1]
    d['Biomass_Merch_Burned_Pct']=df.iloc[3,1]
    d['Biomass_Merch_LeftOnSite_Pct']=df.iloc[4,1]    
    d['Biomass_NonMerch_Removed_Pct']=df.iloc[5,2]
    d['Biomass_NonMerch_Burned_Pct']=df.iloc[3,2]
    d['Biomass_NonMerch_LeftOnSite_Pct']=df.iloc[4,2]    
    d['Snags_Affected_Pct']=df.iloc[1,3]
    d['Snags_Removed_Pct']=df.iloc[5,3]
    d['Snags_Burned_Pct']=df.iloc[3,3]
    d['Snags_LeftOnSite_Pct']=df.iloc[4,3]    
    d['RemovedMerchToFuel_Pct']=df.iloc[12,1]
    d['RemovedMerchToLumber_Pct']=df.iloc[7,1]
    d['RemovedMerchToPlywood_Pct']=df.iloc[8,1]
    d['RemovedMerchToOSB_Pct']=df.iloc[9,1]
    d['RemovedMerchToMDF_Pct']=df.iloc[10,1]
    d['RemovedMerchToPulp_Pct']=df.iloc[11,1]
    d['RemovedMerchToCants_Pct']=df.iloc[13,1]
    d['RemovedMerchToFirewood_Pct']=df.iloc[14,1]    
    d['RemovedNonMerchToFuel_Pct']=df.iloc[12,2]
    d['RemovedNonMerchToLumber_Pct']=df.iloc[7,2]
    d['RemovedNonMerchToPlywood_Pct']=df.iloc[8,2]
    d['RemovedNonMerchToOSB_Pct']=df.iloc[9,2]
    d['RemovedNonMerchToMDF_Pct']=df.iloc[10,2]
    d['RemovedNonMerchToPulp_Pct']=df.iloc[11,2]
    d['RemovedNonMerchToCants_Pct']=df.iloc[13,2]
    d['RemovedNonMerchToFirewood_Pct']=df.iloc[14,2]    
    d['RemovedSnagStemToFuel_Pct']=df.iloc[12,3]
    d['RemovedSnagStemToLumber_Pct']=df.iloc[7,3]
    d['RemovedSnagStemToPlywood_Pct']=df.iloc[8,3]
    d['RemovedSnagStemToOSB_Pct']=df.iloc[9,3]
    d['RemovedSnagStemToMDF_Pct']=df.iloc[10,3]
    d['RemovedSnagStemToPulp_Pct']=df.iloc[11,3]
    d['RemovedSnagStemToCants_Pct']=df.iloc[13,3]
    d['RemovedSnagStemToFirewood_Pct']=df.iloc[14,3]
    
    return d


'''============================================================================
CALCULATE NET SECTOR GREENHOUSE GAS BALANCE
============================================================================'''

def CalculateGHGBalance(v1,meta):
    
    # Conversion factors
    co2_to_c=meta['psl']['bRatio_CO2_to_C']
    co_to_c=meta['psl']['bRatio_CO_to_C']
    ch4_to_c=meta['psl']['bRatio_CH4_to_C']
    EF_N2O_fromCO2=meta['psl']['bEF_N2O_fromCO2']

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
    
        Year=v1[i].Year.copy()
        A=v1[i].A.copy()
        
        if meta['Save Biomass Pools']=='Yes':
            Eco_G_Gross=co2_to_c*np.sum(v1[i].C_G_Gross,axis=2)
            Eco_G_Net=co2_to_c*np.sum(v1[i].C_G_Net,axis=2)
            Eco_M_Reg=co2_to_c*np.sum(v1[i].C_M_Reg,axis=2)
            Eco_NPP=co2_to_c*np.sum(v1[i].C_NPP,axis=2)
            Eco_RH=co2_to_c*np.sum(v1[i].C_RH,axis=2)
            Eco_LF=co2_to_c*np.sum(v1[i].C_LF,axis=2)
        elif meta['Save Biomass Pools']=='No': 
            Eco_G_Gross=co2_to_c*v1[i].C_G_Gross.copy()
            Eco_G_Net=co2_to_c*v1[i].C_G_Net.copy()
            Eco_M_Reg=co2_to_c*v1[i].C_M_Reg.copy()
            Eco_NPP=co2_to_c*v1[i].C_NPP.copy()
            Eco_RH=co2_to_c*v1[i].C_RH.copy()
            Eco_LF=co2_to_c*v1[i].C_LF.copy()
    
        # Carbon dioxide flux (tCO2/ha/yr)
        Eco_E_CO2=co2_to_c*v1[i].C_E_FireAsCO2.copy()
    
        # Carbon monoxide flux (tCO/ha/yr)
        Eco_E_CO=co_to_c*v1[i].C_E_FireAsCO.copy()
    
        # Methan flux *(tCH4/ha/yr)
        Eco_E_CH4=ch4_to_c*v1[i].C_E_FireAsCH4.copy()
    
        # Nitrous oxide flux (tN2O/ha/yr)
        Eco_E_N2O=EF_N2O_fromCO2*Eco_E_CO2.copy()
    
        Eco_E_Fire=gwp_co2*Eco_E_CO2+gwp_co*Eco_E_CO+gwp_ch4*Eco_E_CH4+gwp_n2o*Eco_E_N2O
        
        Eco_E_Operations=gwp_co2*co2_to_c*v1[i].C_E_OperationsAsCO2.copy()
        
        Eco_Removals=co2_to_c*(v1[i].C_RemovedMerch.copy()+ \
                         v1[i].C_RemovedNonMerch.copy()+ \
                         v1[i].C_RemovedSnagStem.copy())        
        
        Eco_NGHGB=Eco_NPP-Eco_RH-Eco_E_Fire-Eco_E_Operations-Eco_Removals
        
        Eco_BiomassStemMerch=co2_to_c*v1[i].C_Eco_Pools.copy()
    
        if meta['Save Biomass Pools']=='Yes':
            Eco_Biomass=co2_to_c*np.sum(v1[i].C_Eco_Pools[:,:,meta['iEP']['BiomassTotal']].copy(),axis=2)       
            Eco_BiomassAG=co2_to_c*np.nansum(v1[i].C_Eco_Pools[:,:,meta['iEP']['BiomassAboveground']].copy(),axis=2)
            Eco_Felled=co2_to_c*np.sum(v1[i].C_Eco_Pools[:,:,meta['iEP']['Felled']].copy(),axis=2)
            Eco_Litter=co2_to_c*np.sum(v1[i].C_Eco_Pools[:,:,meta['iEP']['Litter']].copy(),axis=2)
            Eco_DeadWood=co2_to_c*np.sum(v1[i].C_Eco_Pools[:,:,meta['iEP']['DeadWood']].copy(),axis=2)
            Eco_Soil=co2_to_c*np.sum(v1[i].C_Eco_Pools[:,:,meta['iEP']['Soil']].copy(),axis=2)
            Pro_InUse=co2_to_c*np.sum(v1[i].C_Pro_Pools[:,:,0:10].copy(),axis=2)
            Pro_DumpLandfill=co2_to_c*np.sum(v1[i].C_Pro_Pools[:,:,10:17].copy(),axis=2)
            Pro_Emissions=gwp_co2*co2_to_c*v1[i].C_Pro_Pools[:,:,17].copy()+gwp_ch4*co2_to_c*v1[i].C_Pro_Pools[:,:,18].copy()
                      
        elif meta['Save Biomass Pools']=='No': 
            Eco_Biomass=co2_to_c*v1[i].Eco_Biomass.copy()
            Eco_BiomassAG=co2_to_c*v1[i].Eco_BiomassAG.copy()
            Eco_Felled=co2_to_c*v1[i].Eco_Felled.copy()
            Eco_Litter=co2_to_c*v1[i].Eco_Litter.copy()
            Eco_DeadWood=co2_to_c*v1[i].Eco_DeadWood.copy()
            Eco_Soil=co2_to_c*v1[i].Eco_Soil.copy()
            Pro_InUse=co2_to_c*v1[i].Pro_InUse.copy()
            Pro_DumpLandfill=co2_to_c*v1[i].Pro_DumpLandfill.copy()
            Pro_Emissions=v1[i].Pro_Emissions_co2e.copy() # Already converted to CO2e
        
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
              'Eco_E_Fire':Eco_E_Fire,
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
             'Fire emissions (tCO2e/ha/yr)','Operational emissions (tCO2e/ha/yr)',
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

'''============================================================================
IMPORT PARAMETERS
============================================================================'''

def UpdateParamaters(pthin):

    #------------------------------------------------------------------------------
    # Carbon Pools (from Kurz et al. 2009 with modifications)
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

    df=pd.read_excel(pthin + "\Parameters_Disturbances.xlsx",sheet_name='Sheet1')
    pDist=df.to_dict()    

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
         'BiomassAllomSL':pBiomassAllomSL, \
         'BiomassTurnover':pBiomassTurnover, \
         'InterPoolFluxes':pInterPoolFluxes, \
         'Decomposition':pDecomp, \
         'Disturbances':pDist, \
         'BGC_ZONE':pBGC_ZONE, \
         'SpeciesVRI':pSpeciesVRI, \
         'HWP':pHWP, \
         'SRS':pSRS, \
         'TreeAllometry':pTreeAllometry, \
         'R_Def1':pR_Def1, \
         'M_Def1':pM_Def1, \
         'G_Def1':pG_Def1, \
         }

    # Save to file (this had permission issues on Google Drive)
    gu.opickle(pthin + '\Parameters.pkl',pts)
    
'''============================================================================
TIPSY BUILD INPUT FILE

Notes:
    if the input spreadsheet has nan's, all data will be converted to float

============================================================================'''

def BuildTIPSYInputs(meta):
    
    # Import format info (number of designated positions per variable)
    fin=meta['Path Model Code'] + '\\Parameters\\GrowthCurvesTIPSY_Parameters_Template.xlsx'
    df_frmt=pd.read_excel(fin,sheet_name='Sheet1')

    # Format array
    nfrmt=df_frmt.iloc[1,4:53]

    # Are the growth curves stored all together, or by scenario
    if meta['Growth Curves Lumped']=='Yes':
        n=1
    else:
        n=meta['N Scenario']

    for iScn in range(n):    
        
        # Import input data   
        if meta['Growth Curves Lumped']=='Yes':
            df=pd.read_excel(meta['Path Project'] + '\\Inputs\\GrowthCurvesTIPSY_Parameters.xlsx',sheet_name='Sheet1',skiprows=6)
        else:
            df=pd.read_excel(meta['Path Project'] + '\\Inputs\\Scenario' + FixFileNum(iScn) + '\\GrowthCurvesTIPSY_Parameters.xlsx',sheet_name='Sheet1',skiprows=6)
    
        # Variable names to include
        varnams=df.columns[4:53]

        # Create a spacer
        spc=" " 

        # Number of rows to run
        N_Run=df.shape[0]

        # Do it in batches of no more than 5000 rows because - slow when not defining size upfront
        Ivl_Batch=5000
    
        N_Batch=int(np.ceil(N_Run/Ivl_Batch))
        
        BigString=""
    
        for h in range(0,N_Batch):
                        
            # Index to batch
            iStart=Ivl_Batch*h
            iStop=np.minimum(N_Run,iStart+Ivl_Batch)
            ikp=np.arange(iStart,iStop,1)
            n_Batch=len(ikp)
        
            LittleString="" 
        
            # Loop through rows of table
            for i in range(n_Batch):
            
                ikp0=ikp[i]
            
                for j in range(len(varnams)):
            
                    #print(h + i + j)
                
                    s0=str(df.loc[ikp0,varnams[j]])
    
                    # Make sure NANs are blank for TIPSY
                    if s0=='nan': s0=''
                    if s0=='NA': s0=''
                    if s0=='NaN': s0=''
                    if s0=='': s0=''
    
                    n=nfrmt[j]
                    if len(s0)==n: 
                        s1=s0
                    else: 
                        s1=s0
                        d=n-len(s0)
                        for k in range(d): 
                            s1=s1 + " "
                
                    # Add individual column variables
                    LittleString=LittleString + s1 + spc
                
                # Add end of line
                LittleString=LittleString + '\n'
        
            # Write to big string
            BigString=BigString+LittleString    

        # Save
        if meta['Growth Curves Lumped']=='Yes':
            fid=open(meta['Path Project'] + '\\Inputs\\GrowthCurvesTIPSY_InputVariables.dat','w')
        else:
            fid=open(meta['Path Project'] + '\\Inputs\\Scenario' + FixFileNum(iScn) + '\\GrowthCurvesTIPSY_InputVariables.dat','w')
        fid.write(BigString)
        fid.close()

'''============================================================================
SEVERITY STOCHASTIC PREDICTIONS BASED ON CONSTANT DISTRIBUTION ASSUMPTIONS
============================================================================'''

def SeverityGet(DistType,m):
    
    rn=np.random.random(m)
    y=np.zeros(rn.size)
    
    for i in range(rn.size):
        
        #--------------------------------------------------------------------------
        # Wildfire severity based on a function that predicts burn severity based 
        # on the histagram of burn severity class
        #--------------------------------------------------------------------------
        if DistType=='Wildfire':
        
            if rn[i]<0.24: y[i]=5
            elif (rn[i]>=0.24) & (rn[i]<0.48): y[i]=50
            elif (rn[i]>=0.48) & (rn[i]<0.92): y[i]=90
            elif (rn[i]>=0.92): y[i]=100
        
        elif DistType=='Beetles':
        
            if rn[i]<0.25: y[i]=5
            elif (rn[i]>=0.25) & (rn[i]<0.50): y[i]=25
            elif (rn[i]>=0.50) & (rn[i]<0.75): y[i]=75
            elif (rn[i]>=0.75): y[i]=100
    
        elif DistType=='Defoliators':
        
            y[i]=10
    
        elif DistType=='RootRot':
            
            y[i]=5
        
        elif DistType=='Wind':
            
            y[i]=75
        
        elif DistType=='SnowAndIce':
            
            y[i]=25

    return y

'''============================================================================
IMPORT DISTURBANCE HISTORY
============================================================================'''

def GetDisturbanceHistory(meta):
    dh=[]
    for iScn in range(meta['N Scenario']):
        dhB=[]    
        for iBat in range(meta['N Batch']):
            iEns=0
            dh0=gu.ipickle(meta['Path Input Scenario'][iScn] + '\\Disturbance_Ens' + FixFileNum(iEns) + '_Bat' + FixFileNum(iBat) + '.pkl')
            if iBat==0:
                dhB=dh0
            else:
                for i in range(len(dh0)):
                    dhB.append(dh0[i])
            dh.append(dhB)
    return dh

'''============================================================================
TIME SERIES OF AREA DISTURBED
============================================================================'''

def GetAreasAffected(meta):
    
    it=np.where(meta['Time']>=meta['Year Start Saving'])
    tv=meta['Time'][it]
    
    A=[]
    for iScn in range(meta['N Scenario']):
        A.append({})
        for key in meta['LUT Dist'].keys():
            A[iScn][key]=np.zeros((tv.size,meta['N Stand']))
        for iBat in range(meta['N Batch']):            
            iEns=0            
            dh=gu.ipickle(meta['Path Input Scenario'][iScn] + '\\Disturbance_Ens' + FixFileNum(iEns) + '_Bat' + FixFileNum(iBat) + '.pkl')    
            for iS in range(len(dh)):                
                for iY in range(dh[iS]['Year'].size):
                    it=np.where(tv==np.round(dh[iS]['Year'][iY]))[0]
                    if it.size==0: 
                        continue
                    for key in meta['LUT Dist'].keys():
                        if dh[iS]['ID_Type'][iY]==meta['LUT Dist'][key]:
                            A[iScn][key][it,iS]=A[iScn][key][it,iS]+1
        
        #----------------------------------------------------------------------
        # Summaries
        #----------------------------------------------------------------------
        
        A[iScn]['Harvest (all)']=A[iScn]['Harvest and Slashpile Burn']+ \
            A[iScn]['Harvest']+ \
            A[iScn]['Harvest Custom']
    
    return A

'''============================================================================
GRAPHICS

# Look at variables in rcParams:
plt.rcParams.keys()
============================================================================'''

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
