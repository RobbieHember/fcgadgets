
import os
import numpy as np
import pandas as pd
import copy
import gc as garc
import time
import datetime
from fcgadgets.macgyver import utilities_general as gu
from fcgadgets.cbrunner import cbrun_utilities as cbu
from fcgadgets.cbrunner import cbrun_annproc as annproc

#%% Run simulation

def MeepMeep(meta):
    
    # Identify scenarios to run
    # If 'Scenario Override' is absent, it defaults to running all scenarios. 
    # Specifying a list of scenarios in 'Scenario Override' facilitates running
    # the simulation in multiple instances of Python.
    if 'Scenario Override' in meta['Project']:
        ScenariosToRun=meta['Project']['Scenario Override']['Scenario List']
    else:
        ScenariosToRun=range(meta['Project']['N Scenario'])
    
    # Initialize run time tracking info
    rt_info={}
    
    # Loop through batches
    for iBat in range(meta['Project']['N Batch']):
        
        flag_WorkingOnBatch=0
        
        # Loop through scenarios    
        for iScn in ScenariosToRun: 
        
            # Track the current scenario index
            meta['iScn']=iScn
        
            # Loop through ensembles
            for iEns in range(meta['Project']['N Ensemble']):
                
                # Track time
                t0=time.time()
                rt_info['t1']=time.time()
                
                # Path to temporary "working on batch..." file -> tells other instances 
                pth_WorkingOnBatch=meta['Paths']['Project'] + '\\Outputs\\WorkingOnBatch_' + cbu.FixFileNum(iBat) + '.pkl'
                
                # Only proceed if the file does not exist running multiple instances
                if meta['Project']['Skip Completed Runs']=='On':                    
                    if os.path.exists(pth_WorkingOnBatch)==False:
                        gu.opickle(pth_WorkingOnBatch,[])
                        flag_WorkingOnBatch=1
                    else:
                        if flag_WorkingOnBatch==0:
                            continue
            
                # Report progress
                if (meta['Project']['Scenario Source']=='Spreadsheet'):
                    #print('Running Scenario ' + cbu.FixFileNum(iScn) )
                    pass
                else:
                    print('Running Scenario ' + cbu.FixFileNum(iScn) + ', Ensemble ' + cbu.FixFileNum(iEns) + ', Batch ' + cbu.FixFileNum(iBat))
                
                # Initialize stands
                meta,vi,vo=InitializeStands(meta,iScn,iEns,iBat)
                
                rt_info['t2']=time.time()
                
                # Set location-specific parameters
                meta,vi=PrepareParametersForBatch(meta,vi,iEns,iBat)
                
                rt_info['t3']=time.time()
              
                # Indices to ecosystem pools
                iEP=meta['Core']['iEP']
                
                # Biomass dynamics from Sawtooth
                if meta['Project']['Biomass Module']=='Sawtooth':                    
                    for iS in range(meta['Project']['Batch Size'][iBat]):                        
                        vo=annproc.BiomassFromSawtooth(iScn,iS,vi,vo,meta,iEP)
                
                # Try to start at a later date and adopt spinup from a previous run
                #t_start=1
                if (meta['Project']['N Ensemble']==1) | (iScn==0) & (iEns==0):
                    
                    t_start=1
                
                else:
                    
                    t_start=np.where(meta['Year']==meta['Project']['Year Start Saving'])[0][0]
                    
                    # Get output variables from first run of this batch
                    #vo=copy.deepcopy(vo_full)
                    it=np.where(meta['Year']==meta['Project']['Year Start Saving']-1)[0]                    
                    for k in vo.keys():
                        if (k=='C_M_ByAgent'):
                            continue
                        vo[k][it,:]=vo_full[k]
                    
                    # Set future periods to zero
                    it=np.where(meta['Year']>=meta['Project']['Year Start Saving'])[0]
                    for k in vo.keys():
                        if (k=='C_M_ByAgent'):
                            # Nested dictionaries
                            # *** These are already shortened time periods so set whole array to zero
                            for k2 in vo[k].keys():
                                vo[k][k2]=0*vo[k][k2]
                        else:
                            # Not a nested dictionary
                            vo[k][it,:]=0*vo[k][it,:]
                    
                rt_info['t4']=time.time()    
                
                # Loop through time intervals (start in second time step)
                for iT in range(t_start,meta['Project']['N Time']):
                    
                    # Biomass dynamics
                    if meta['Project']['Biomass Module']=='BatchTIPSY':
                        
                        # Biomass from BatchTIPSY.exe (or TASS)
                        vo=annproc.Biomass_FromTIPSYorTASS(iScn,iBat,iT,vi,vo,meta,iEP)
                    
                    # Calculate annual dead organic matter dynamics
                    vo=annproc.DOM_like_CBM08(iT,iBat,vi,vo,iEP,meta)
                    
                    # Calculate effects of disturbance and management                    
                    vo,vi=annproc.Events_FromTaz(iT,iScn,iEns,iBat,vi,vo,meta,iEP)
                    
                    # Calculate products sector                    
                    if meta['Year'][iT]>=meta['Core']['HWP Year Start']:
                        
                        # No need to run this before a certain date
                        vo=annproc.HWP_From_BCHWP12_Update21(iT,iBat,vi,vo,meta)
                
                rt_info['t5']=time.time()
                
                # Export simulation results to file
                if (meta['Project']['N Ensemble']==1) | (iScn==0) & (iEns==0):
                    # Only save output if it is the first instance of a batch
                    vo_full=ExportSimulation(meta,vi,vo,iScn,iEns,iBat,iEP)
                    
                else:                   
                    ExportSimulation(meta,vi,vo,iScn,iEns,iBat,iEP)
                
                rt_info['t6']=time.time()
                
                # Delete 'working on' file
                #if meta['Project']['Skip Completed Runs']=='On':   
                    #os.remove(pthWO)
                    
                # Delete variables
                del vi,vo
                garc.collect()
                
                # Track simulation time
                #print(rt_info['t2']-rt_info['t1'])
                #print(rt_info['t3']-rt_info['t2'])
                #print(rt_info['t4']-rt_info['t3'])
                #print(rt_info['t5']-rt_info['t4'])
                #print(rt_info['t6']-rt_info['t5'])
                t1=time.time()
                #print(t1-t0)
                
            # Calculate and save model output stats for this ensemble, delete
            # full model output data to save space            
            if meta['Project']['Save MOS on fly']=='On':                
                SaveOutputToMOS(meta,iScn,iEns)
            
#%% Initialize stands

def InitializeStands(meta,iScn,iEns,iBat):
    
    #--------------------------------------------------------------------------
    # Input variables
    #--------------------------------------------------------------------------

    vi={}
        
    vi['tv']=meta['Year']
    tv=meta['Year']
    meta['Project']['N Time']=len(tv)
    
    # Import inventory
    vi['Inv']=gu.ipickle(meta['Paths']['Input Scenario'][iScn] + '\\Inventory_Bat' + cbu.FixFileNum(iBat) + '.pkl')
        
    # Update number of stands for batch
    meta['Project']['N Stand Batch']=vi['Inv']['X'].shape[1]          
    
    # Import event chronology
    vi['EC']=gu.ipickle(meta['Paths']['Input Scenario'][iScn] + '\\Events_Ens' + cbu.FixFileNum(iEns) + '_Bat' + cbu.FixFileNum(iBat) + '.pkl')
    
    vi['EC']=cbu.EventChronologyDecompress(meta,vi['EC'],iScn,iEns,iBat)    
    
    # Convert mortality percent to fraction
    vi['EC']['MortalityFactor']=vi['EC']['MortalityFactor'].astype(float)/100
    
    # Import growth curves
    if (meta['Project']['Biomass Module']=='BatchTIPSY'):
        
        vi['GC']={}
        
        # Import growth curve 1
        vi['GC'][1]=gu.ipickle(meta['Paths']['Input Scenario'][iScn] + '\\GrowthCurve1_Bat' + cbu.FixFileNum(iBat) + '.pkl')
        
        # Set active growth curve to growth curve 1
        vi['GC']['Active']=vi['GC'][1].copy()
        
        # Initialize an indicator of the active growth curve ID
        vi['GC']['ID_GCA']=np.ones(meta['Project']['Batch Size'][iBat])
        
        # Import growth curve 2
        vi['GC'][2]=gu.ipickle(meta['Paths']['Input Scenario'][iScn] + '\\GrowthCurve2_Bat' + cbu.FixFileNum(iBat) + '.pkl')
        
        # Import growth curve 3
        vi['GC'][3]=gu.ipickle(meta['Paths']['Input Scenario'][iScn] + '\\GrowthCurve3_Bat' + cbu.FixFileNum(iBat) + '.pkl')
        
        # Import growth curve 4 (optional)
        try:
            vi['GC'][4]=gu.ipickle(meta['Paths']['Input Scenario'][iScn] + '\\GrowthCurve4_Bat' + cbu.FixFileNum(iBat) + '.pkl')
        except:
            vi['GC'][4]=0
        
        # Import growth curve 5 (optional)
        try:
            vi['GC'][5]=gu.ipickle(meta['Paths']['Input Scenario'][iScn] + '\\GrowthCurve5_Bat' + cbu.FixFileNum(iBat) + '.pkl')
        except:
            vi['GC'][5]=0
        
    else:        
        vi['GC']={}               
        vi['GC'][1]=0
        vi['GC'][2]=0
        vi['GC'][3]=0
        vi['GC'][4]=0
        vi['GC'][5]=0
        vi['GC']['Active']=0 
        vi['GC']['ID_GCA']=np.ones(meta['Project']['Batch Size'][iBat])
    
    #--------------------------------------------------------------------------
    # Initialize output variables
    #--------------------------------------------------------------------------
        
    # Output variables dictionary   
    vo={}
    
    # Get dimensions    
    m=meta['Project']['N Time']
    n=meta['Project']['Batch Size'][iBat]
    o=meta['Core']['N Pools Eco']
    
    # Stand age (i.e. time since stand-replacing disturbance)    
    vo['A']=np.zeros((m,n))
      
    # Species composition
    #vo['Spc_ID']=np.zeros((m,n,6),dtype='int8')
    #vo['Spc_Percent']=np.zeros((m,n,6),dtype='int8')
    
    # Stemwood merch volume
    vo['V_StemMerch']=np.zeros((m,n))
    
    # Carbon density of ecosystem (Mg C ha-1)
    # -> 3-D matrix: Time x Stand x Carbon pool
    vo['C_Eco_Pools']=np.zeros((m,n,o))
        
    # Carbon density of products sector (Mg C ha-1)
    # -> 3-D matrix: Time x Stand x Carbon pool
    vo['C_Pro_Pools']=np.zeros((m,n,meta['Core']['N Pools Pro']))  
    
    # Aggregate pools (Mg C ha-1) (polulated upon export)
    vo['C_Biomass_Tot']=np.array([])
    vo['C_Felled_Tot']=np.array([])
    vo['C_Litter_Tot']=np.array([])
    vo['C_DeadWood_Tot']=np.array([])
    vo['C_Soil_Tot']=np.array([])
    vo['C_InUse_Tot']=np.array([])
    vo['C_DumpLandfill_Tot']=np.array([])
    
    # Carbon flux densities (Mg C ha-1 yr-1)
    vo['C_NPP']=np.zeros((m,n,o))
    vo['C_G_Gross']=np.zeros((m,n,o))
    vo['C_G_Net']=np.zeros((m,n,o))    
    vo['C_M_Reg']=np.zeros((m,n,o))
    vo['C_M_Dist']=np.zeros((m,n))
    vo['C_M_ByAgent']={}
    for k in meta['LUT']['Dist']:
        vo['C_M_ByAgent'][k]=np.zeros((m,1))        
    vo['C_LF']=np.zeros((m,n,o))    
    vo['C_RH']=np.zeros((m,n,o))  
    
    # Aggregate pools (Mg C ha-1) (polulated upon export)
    vo['C_G_Gross_Tot']=np.array([])
    vo['C_G_Net_Tot']=np.array([])
    vo['C_M_Reg_Tot']=np.array([])
    vo['C_LF_Tot']=np.array([])
    vo['C_RH_Tot']=np.array([])
    
    # Keep track of carbon transfers
    vo['C_ToMillMerch']=np.zeros((m,n))
    vo['C_ToMillNonMerch']=np.zeros((m,n))
    vo['C_ToMillSnagStem']=np.zeros((m,n))
    vo['C_ToSlashpileBurn']=np.zeros((m,n))
    vo['C_ToLumber']=np.zeros((m,n))
    vo['C_ToPlywood']=np.zeros((m,n))
    vo['C_ToOSB']=np.zeros((m,n))
    vo['C_ToMDF']=np.zeros((m,n))
    vo['C_ToPolePost']=np.zeros((m,n))
    vo['C_ToPaper']=np.zeros((m,n))
    vo['C_ToPowerGeneration']=np.zeros((m,n))
    vo['C_ToPellets']=np.zeros((m,n))
    vo['C_ToFirewood']=np.zeros((m,n))
    
    # Emissions from wildfire and open burning (will be deleted upon export)
    vo['C_E_FireAsCO2']=np.zeros((m,n))
    vo['C_E_FireAsCH4']=np.zeros((m,n))
    vo['C_E_FireAsCO']=np.zeros((m,n))
    vo['C_E_FireAsN2O']=np.zeros((m,n))
    
    # LULUCF Sector: Net ecosystem production (polulated in load results)
    vo['CO2e_LULUCF_NEP']=np.array([])
    
    # LULUCF Sector: Net ecosystem production (polulated upon export)
    vo['CO2e_LULUCF_E_Wildfire']=np.array([])
    
    # LULUCF Sector: Net ecosystem production (polulated upon export)
    vo['CO2e_LULUCF_E_OpenBurning']=np.array([])
    
    # LULUCF Sector: Denitrification, valatilization
    vo['CO2e_LULUCF_E_EcoOther']=np.zeros((m,n))
    
    # LULUCF Sector: RH from dumps and landfills
    vo['CO2e_LULUCF_E_HWP']=np.zeros((m,n)) 
    
    # Energy - Stationary Combustion Sector
    vo['CO2e_StatComb_E']=np.zeros((m,n))
    
    # Energy - Transporation Sector
    vo['CO2e_Transp_E']=np.zeros((m,n))
    
    # Inudstrial Produciton and Product Use Sector
    vo['CO2e_IPPU_E']=np.zeros((m,n))
    
    #vo['C_M_Sim_Reg']=np.zeros((m,n))
    #vo['C_M_Sim_Fir']=np.zeros((m,n))
    #vo['C_M_Sim_Ins']=np.zeros((m,n))
    #vo['C_M_Sim_Pat']=np.zeros((m,n))
    #vo['C_M_Sim_Har']=np.zeros((m,n))
    #vo['C_M_Sim_Win']=np.zeros((m,n))   
    
    # Stand density
    #vo['N']=np.zeros((m,n))
    
    # Change in stand density (stems ha-1 yr-1)
    #vo['N_R']=np.zeros((m,n))
    #vo['N_M_Tot']=np.zeros((m,n))
    #vo['N_M_Reg']=np.zeros((m,n))
    #vo['N_M_Inv_Fir']=np.zeros((m,n))
    #vo['N_M_Inv_Ins']=np.zeros((m,n))
    #vo['N_M_Inv_Pat']=np.zeros((m,n))
    #vo['N_M_Inv_Har']=np.zeros((m,n))
    #vo['N_M_Inv_Win']=np.zeros((m,n))
    #vo['N_M_Sim_Reg']=np.zeros((m,n))
    #vo['N_M_Sim_Fir']=np.zeros((m,n))
    #vo['N_M_Sim_Ins']=np.zeros((m,n))
    #vo['N_M_Sim_Pat']=np.zeros((m,n))
    #vo['N_M_Sim_Har']=np.zeros((m,n))
    #vo['N_M_Sim_Win']=np.zeros((m,n))
    
    # Mean of tree attributes
    #vo['TreeMean_A']=np.zeros((m,n))
    #vo['TreeMean_H']=np.zeros((m,n))
    #vo['TreeMean_D']=np.zeros((m,n))
    #vo['TreeMean_Csw']=np.zeros((m,n))
    #vo['TreeMean_Csw_G']=np.zeros((m,n))
    
    # Lables for plotting
#    HandleLabels=['Stand age (years)','Gross growth (MgC/ha/yr)',
#             'Net growth (MgC/ha/yr)','Background mortality (MgC/ha/yr)',
#             'Litterfall (MgC/ha/yr)','NPP (MgC/ha/yr)','RH (MgC/ha/yr)',
#             'Disturbance mortality (MgC/ha/yr)',
#             'Removed merch. (MgC/ha/yr)',
#             'Removed non-merch. (MgC/ha/yr)','Removed dead (MgC/ha/yr)',
#             'Wildfire emissions as CH4 (MgC/ha/yr)','Wildfire emissions as CO (MgC/ha/yr)',
#             'Wildfire emissions as CO2 (MgC/ha/yr)','Wildfire emissions as N2O (MgC/ha/yr)',
#             'Operation emissions as CO2 (MgC/ha/yr)','Merch. volume (m3/ha)',
#             'Biomass (MgC/ha)','Aboveground biomass (MgC/ha)','Felled (MgC/ha)',
#             'Litter (MgC/ha)','Dead wood (MgC/ha)',
#             'Soil (MgC/ha)','Total ecosystem carbon (MgC/ha)','In-use products (MgC/ha)',
#             'Dump and landfill (MgC/ha)','Product sector E (MgC/ha/yr)']        
    #'Wildfire mortality (MgC/ha/yr)','Harvest mortality (MgC/ha/yr)',
    #         'Insect mortality (MgC/ha/yr)','Pathogen mortality (MgC/ha/yr)',
    #         'Wind mortality (MgC/ha/yr)',
    #Keys=np.array(list(vo.keys()))
    #d1={}
    #for j in range(len(HandleLabels)):
    #    d1[Keys[j]]=HandleLabels[j]  
    #meta['Labels Output Variables']=d1
    
    #--------------------------------------------------------------------------
    # Configure batch-specific setttings
    #--------------------------------------------------------------------------
                
    # Nutrient application response yearly counter
    meta['Nutrient Management']['ResponseCounter']=np.zeros(meta['Project']['Batch Size'][iBat])
    
    # Identify the start and end of a fixed spin up disturbance interval
    # This is used to fast-track through the spinup period.   
    if meta['Project']['Biomass Module']=='Sawtooth':        
        d=np.append(0,np.diff(vi['EH'][0]['Year']))
        ind=np.where(d==meta['Project']['Spinup Disturbance Return Inverval'])[0]
        meta['Project']['SpinupSpanFastTrack']=[np.min(vi['EH'][0]['Year'][ind])+meta['Project']['Spinup Disturbance Return Inverval'],np.max(vi['EH'][0]['Year'][ind])]
    
    # Initialize flag for fixing negative net growth. When TIPSY yields negative
    # net growth, the fluxes of gross growth and mortality need adjustment. 
    # This flag helps achieve that.    
    meta['FlagNegNetGrowth']=np.zeros(meta['Project']['Batch Size'][iBat])
    meta['G_Net_PriorToBreakup']=np.zeros((meta['Project']['Batch Size'][iBat],7))
    
    #--------------------------------------------------------------------------
    # Get random numbers for on-the-fly disturbance types
    #--------------------------------------------------------------------------
    
    if (meta['Scenario'][iScn]['Harvest Status Historical']=='On') | (meta['Scenario'][iScn]['Harvest Status Future']=='On'):
        rn=gu.ipickle(meta['Paths']['Project'] + '\\Inputs\\Ensembles\\RandomNumbers_Harvest_Ens' + cbu.FixFileNum(iEns) + '_Bat' + cbu.FixFileNum(iBat) + '.pkl')
        rn=rn.astype(float)
        rn=rn*meta['Project']['On the Fly']['Random Numbers']['Scale Factor']        
        meta['Project']['On the Fly']['Random Numbers']['Harvest']=rn.copy()
    
    if (meta['Scenario'][iScn]['Breakup Status']=='On'):
        rn=gu.ipickle(meta['Paths']['Project'] + '\\Inputs\\Ensembles\\RandomNumbers_Breakup_Ens' + cbu.FixFileNum(iEns) + '_Bat' + cbu.FixFileNum(iBat) + '.pkl')
        rn=rn.astype(float)
        rn=rn*meta['Project']['On the Fly']['Random Numbers']['Scale Factor']
        meta['Project']['On the Fly']['Random Numbers']['Breakup']=rn.copy()
    
    #--------------------------------------------------------------------------
    # Initialize a log book that will record various diagnostics, warnings, 
    # and error flags    
    #--------------------------------------------------------------------------
    
    meta['Logbook']=list()
    
    return meta,vi,vo

#%% Import parameters

def PrepareParametersForBatch(meta,vi,iEns,iBat):
      
    #--------------------------------------------------------------------------
    # Populate final parameters with initial best estimates
    #--------------------------------------------------------------------------
    
    meta['Param']['BEV']=copy.deepcopy(meta['Param']['BE'])
    
    #--------------------------------------------------------------------------
    # Add error variance
    #--------------------------------------------------------------------------
    
    if meta['Project']['Uncertainty NM Status']=='On':
        
        for k in meta['Param']['BE']['Nutrient Management'].keys():        
            meta['Param']['BEV']['Nutrient Management'][k]=meta['Param']['By Ensemble'][iEns]['Nutrient Management'][k]           
    
    #--------------------------------------------------------------------------
    # Biomass allometry (stand level)
    #--------------------------------------------------------------------------
    
    # Unique BGC zones
    u=np.unique(vi['Inv']['ID_BECZ'].flatten())
    
    MaritimeZones=['CDF','CWH','ICH']
    
    for k in meta['Param']['BEV']['Biomass Allometry']['Raw'].keys():
        
        if k=='Region':
            continue    
        
        meta['Param']['BEV']['Biomass Allometry'][k]=np.zeros(meta['Project']['Batch Size'][iBat])
        
        for iU in range(u.size):
            
            bgc_cd=cbu.lut_n2s(meta['LUT']['VRI']['BEC_ZONE_CODE'],u[iU])
            
            ind=np.where(vi['Inv']['ID_BECZ'].flatten()==u[iU])[0]
            
            if np.isin(bgc_cd,MaritimeZones)==True:
                meta['Param']['BEV']['Biomass Allometry'][k][ind]=meta['Param']['BEV']['Biomass Allometry']['Raw'][k][0]
            else:
                meta['Param']['BEV']['Biomass Allometry'][k][ind]=meta['Param']['BEV']['Biomass Allometry']['Raw'][k][1]
    
    #--------------------------------------------------------------------------
    # Biomass turnover
    #--------------------------------------------------------------------------
    
    PoolNames=np.array(meta['Core']['Name Pools Eco'])[meta['Core']['iEP']['BiomassTotal']]
    
    indPar=np.where( (meta['Param']['BEV']['Biomass Turnover']['Raw']['SPECIES_CD']=='PL') )[0]
    
    for pn in PoolNames:
        meta['Param']['BEV']['Biomass Turnover'][pn]=meta['Param']['BEV']['Biomass Turnover']['Raw'][pn][indPar]*np.ones((1,meta['Project']['Batch Size'][iBat]))
        
    #--------------------------------------------------------------------------
    # HWP
    #--------------------------------------------------------------------------
    
    for i in range(len(meta['Param']['BEV']['HWP']['Raw'])):
        
        Name=meta['Param']['BEV']['HWP']['Raw']['Name'].iloc[i]
        Value=meta['Param']['BEV']['HWP']['Raw']['Value'].iloc[i]
        
        # Convert half life to turnover rate
        if Name[-2:]=='hl':            
            Name=Name[0:-2] + 'tr'
            Value=1-np.exp(-np.log(2)/Value)        
        
        meta['Param']['BEV']['HWP'][Name]=Value*np.ones(meta['Project']['Batch Size'][iBat])
    
    #--------------------------------------------------------------------------
    # Populate custom harvest parameters with those supplied for project
    #--------------------------------------------------------------------------
    
    if 'Harvest Custom' in meta:
        
        for iHC in range(10):
            
            if int(iHC+1) in meta['Harvest Custom']:
                
                for k in meta['Harvest Custom'][int(iHC+1)].keys():
                    
                    if k[0:7]!='Removed':
                        id=meta['LUT']['Dist']['Harvest Custom ' + str(int(iHC+1))]
                        meta['Param']['BEV']['Dist'][id][k]=meta['Harvest Custom'][int(iHC+1)][k]/100 

    return meta,vi

#%% Export simulation results

def ExportSimulation(meta,vi,vo,iScn,iEns,iBat,iEP):
    
    #--------------------------------------------------------------------------
    # Save project metadata    
    #--------------------------------------------------------------------------
    
    if (iBat==0) & (iScn==0):
        gu.opickle(meta['Paths']['Project'] + '\\Inputs\\Metadata.pkl',meta)
    
    #--------------------------------------------------------------------------
    # Save a full copy of the first ensemble for re-use in big projects
    #--------------------------------------------------------------------------
    
    if (meta['Project']['N Ensemble']>1) & (iScn==0) & (iEns==0):        
        it=np.where(meta['Year']==meta['Project']['Year Start Saving']-1)[0]
        vo_full={}
        for k in vo.keys():
            if (k=='C_M_ByAgent'):
                continue
            vo_full[k]=vo[k][it,:]        
        #vo_full=copy.deepcopy(vo)
    else:
        vo_full=[]
    
    #--------------------------------------------------------------------------
    # Extract parameters
    #--------------------------------------------------------------------------
    
    bB=meta['Param']['BEV']['Biophysical']
    
    #--------------------------------------------------------------------------
    # Isolate time period that will be saved to file
    #--------------------------------------------------------------------------
    
    it=np.where(meta['Year']>=meta['Project']['Year Start Saving'])[0]
    
    # Main output variables
    for k in vo:
        
        # Skip mortality summary
        if type(vo[k])==dict:
            continue
        
        # Skip variables that were initiated but not yet populated
        if vo[k].size==0:
            continue

        vo[k]=vo[k][it,:]
    
    # Mortality summaries
    for k in vo['C_M_ByAgent']:
        vo['C_M_ByAgent'][k]=vo['C_M_ByAgent'][k][it,:]
    
    #--------------------------------------------------------------------------
    # Convert fire to CO2e
    #--------------------------------------------------------------------------
    
    # Carbon dioxide flux (tCO2/ha/yr)
    E_Fire_CO2=bB['Ratio_CO2_to_C']*vo['C_E_FireAsCO2']
    
    # Carbon monoxide flux (tCO/ha/yr)
    E_Fire_CO=bB['Ratio_CO_to_C']*vo['C_E_FireAsCO']
    
    # Methan flux *(tCH4/ha/yr)
    E_Fire_CH4=bB['Ratio_CH4_to_C']*vo['C_E_FireAsCH4']
    
    # Nitrous oxide flux (tN2O/ha/yr)
    E_Fire_N2O=bB['EF_N2O_fromCO2']*E_Fire_CO2
    
    # Convert fluxes to CO2e using global warming potential estimates
    CO2e_E_Fire_AsCO2=1*E_Fire_CO2    
    CO2e_E_Fire_AsCH4=bB['GWP_CH4_AR5']*E_Fire_CH4    
    CO2e_E_Fire_AsCO=bB['GWP_CO_AR5']*E_Fire_CO    
    CO2e_E_Fire_AsN2O=bB['GWP_N2O_AR5']*E_Fire_N2O
    
    CO2e_E_Fire_Total=CO2e_E_Fire_AsCO2+CO2e_E_Fire_AsCH4+CO2e_E_Fire_AsCO+CO2e_E_Fire_AsN2O
    
    #--------------------------------------------------------------------------
    # Emissions from wildfire
    #--------------------------------------------------------------------------
    
    vo['CO2e_LULUCF_E_Wildfire']=np.zeros(vo['C_ToMillMerch'].shape)
    ind=np.where( (vo['C_ToMillMerch']==0) & (vo['C_ToMillSnagStem']==0) )    
    if ind[0].size>0:
        vo['CO2e_LULUCF_E_Wildfire'][ind[0],ind[1]]=CO2e_E_Fire_Total[ind[0],ind[1]]

    #--------------------------------------------------------------------------
    # Emissions from open burning    
    #--------------------------------------------------------------------------
    
    vo['CO2e_LULUCF_E_OpenBurning']=np.zeros(vo['C_ToMillMerch'].shape)
    ind=np.where( (vo['C_ToMillMerch']>0) | (vo['C_ToMillSnagStem']>0) )    
    if ind[0].size>0:
        vo['CO2e_LULUCF_E_OpenBurning'][ind[0],ind[1]]=CO2e_E_Fire_Total[ind[0],ind[1]]
    
    #--------------------------------------------------------------------------
    # Remove unnecessary fire emission variables
    #--------------------------------------------------------------------------
    
    del vo['C_E_FireAsCO2']
    del vo['C_E_FireAsCO']
    del vo['C_E_FireAsCH4']
    del vo['C_E_FireAsN2O']
    
    #--------------------------------------------------------------------------
    # Substitution effects
    # This needs to be here so that parameter uncertainty can be considered.
    #--------------------------------------------------------------------------
    
    if meta['Scenario'][iScn]['Substitution Effects Status']=='On':
    
        # Extract parameters
        bS=meta['Param']['BEV']['Substitution']
        WoodDensity=0.42
        
        # Power generation (MgC/ha) to (green tonne/ha)    
        PowerGeneration=vo['C_ToPowerGeneration']/WoodDensity/bB['WoodMoistureContent']
        
        # Power generation (GJ/ha)
        GJ=bB['Wood Fuel Industrial (50% moisture) (GJ/kg)']*(PowerGeneration*1000)
        
        E_Coal=bS['PowerGenFracDisplacingCoal']*bB['Coal Emission (kgCO2e/GJ)']*GJ/1000
        E_Diesel=bS['PowerGenFracDisplacingDiesel']*bB['Diesel Fuel Emission (kgCO2e/GJ)']*GJ/1000
        E_NatGas=bS['PowerGenFracDisplacingNaturalGas']*bB['Natural gas Emission (kgCO2e/GJ)']*GJ/1000
        E_Oil=bS['PowerGenFracDisplacingOil']*bB['Oil Emission (kgCO2e/GJ)']*GJ/1000
        E_PowerGenSub=E_Coal+E_Diesel+E_NatGas+E_Oil
        
        # Pellets (MgC/ha) to (kiln dried tonne/ha)
        Pellets=vo['C_ToPellets']/WoodDensity
        
        # Pellets (GJ/ha)
        GJ=bB['Wood Fuel Kiln-dried (GJ/kg)']*(Pellets*1000)
        
        E_Coal=bS['PelletFracDisplacingCoal']*bB['Coal Emission (kgCO2e/GJ)']*GJ/1000
        E_Diesel=bS['PelletFracDisplacingDiesel']*bB['Diesel Fuel Emission (kgCO2e/GJ)']*GJ/1000
        E_NatGas=bS['PelletFracDisplacingNaturalGas']*bB['Natural gas Emission (kgCO2e/GJ)']*GJ/1000
        E_Oil=bS['PelletFracDisplacingOil']*bB['Oil Emission (kgCO2e/GJ)']*GJ/1000
        E_PelletSub=E_Coal+E_Diesel+E_NatGas+E_Oil
        
        # Firewood (MgC/ha) to (air dried tonne/ha)
        Firewood=vo['C_ToFirewood']/WoodDensity
        
        # Firewood (GJ/ha)
        GJ=bB['Wood Fuel Residential (0% moisture) (GJ/kg)']*(Firewood*1000)
        
        E_Coal=bS['FirewoodFracDisplacingCoal']*bB['Coal Emission (kgCO2e/GJ)']*GJ/1000
        E_Diesel=bS['FirewoodFracDisplacingDiesel']*bB['Diesel Fuel Emission (kgCO2e/GJ)']*GJ/1000
        E_NatGas=bS['FirewoodFracDisplacingNaturalGas']*bB['Natural gas Emission (kgCO2e/GJ)']*GJ/1000
        E_Oil=bS['FirewoodFracDisplacingOil']*bB['Oil Emission (kgCO2e/GJ)']*GJ/1000
        E_FirewoodSub=E_Coal+E_Diesel+E_NatGas+E_Oil
        
        # Sum substitution of power generation, pellets, and firewood
        vo['CO2e_StatComb_E']=vo['CO2e_StatComb_E']-1*(E_PowerGenSub+E_PelletSub+E_FirewoodSub)
        
        # Building materials (tCO2e/ha)
        E_LumberSub=bB['Ratio_CO2_to_C']*bS['LumberDisplacementFactor']*vo['C_ToLumber']
        E_PanelSub=bB['Ratio_CO2_to_C']*bS['PanelDisplacementFactor']*(vo['C_ToPlywood']+vo['C_ToOSB']+vo['C_ToMDF'])
        
        # Sum substitution of building materials
        vo['CO2e_StatComb_E']=vo['CO2e_StatComb_E']-1*(E_LumberSub+E_PanelSub)
    
    #--------------------------------------------------------------------------
    # Sum pools and fluxes if "Save Biomass Pools" is Off
    # These variables only need to be saved to file if the indiviudal pools are
    # not being stored.
    #--------------------------------------------------------------------------    
    
    if meta['Project']['Save Biomass Pools']!='On':
        
        # Aggregate pools
        vo['C_Biomass_Tot']=np.sum(vo['C_Eco_Pools'][:,:,iEP['BiomassTotal']],axis=2)
        vo['C_Felled_Tot']=np.sum(vo['C_Eco_Pools'][:,:,iEP['Felled']],axis=2)
        vo['C_Litter_Tot']=np.sum(vo['C_Eco_Pools'][:,:,iEP['Litter']],axis=2)
        vo['C_DeadWood_Tot']=np.sum(vo['C_Eco_Pools'][:,:,iEP['DeadWood']],axis=2)
        vo['C_Soil_Tot']=np.sum(vo['C_Eco_Pools'][:,:,iEP['Soil']],axis=2)
        vo['C_InUse_Tot']=np.sum(vo['C_Pro_Pools'][:,:,meta['Core']['iPP']['InUse']],axis=2)
        vo['C_DumpLandfill_Tot']=np.sum(vo['C_Pro_Pools'][:,:,meta['Core']['iPP']['DumpLandfill']],axis=2)
        
        # Remove pools
        del vo['C_Eco_Pools']
        del vo['C_Pro_Pools']
        
        # Aggregate fluxes
        vo['C_G_Gross_Tot']=np.sum(vo['C_G_Gross'],axis=2)
        vo['C_G_Net_Tot']=np.sum(vo['C_G_Net'],axis=2)
        vo['C_M_Reg_Tot']=np.sum(vo['C_M_Reg'],axis=2)
        vo['C_LF_Tot']=np.sum(vo['C_LF'],axis=2)
        vo['C_RH_Tot']=np.sum(vo['C_RH'],axis=2)
        
        # Remove fluxes
        del vo['C_G_Gross']
        del vo['C_G_Net']
        del vo['C_G_M_Reg']
        del vo['C_G_LF']
        del vo['C_G_RH']
        
    #--------------------------------------------------------------------------
    # Apply scale factor to data
    #--------------------------------------------------------------------------    
    
    for k in vo.keys():
        
        # Skip mortality summary by agent
        if (k=='C_M_ByAgent'):
            continue
        
        if vo[k].size==0:
            continue
        
        # This one variable needs a larger scale factor
        if (k=='CO2e_LULUCF_HWP'):
            vo[k]=vo[k]/meta['Core']['Scale Factor Export Big']
        else:
            vo[k]=vo[k]/meta['Core']['Scale Factor Export Small']
        
        if np.max(vo[k])<32767:
            vo[k]=vo[k].astype('int16')
        else:
            vo[k]=vo[k].astype(int)
    
    # Mortality summary by agent
    for k in vo['C_M_ByAgent'].keys():
        if np.max(vo['C_M_ByAgent'][k]<32767):
            vo['C_M_ByAgent'][k]=vo['C_M_ByAgent'][k].astype('int16')
        else:
            vo['C_M_ByAgent'][k]=vo['C_M_ByAgent'][k].astype(int)  
    
    #--------------------------------------------------------------------------  
    # Sawtooth variables
    #--------------------------------------------------------------------------  
    
#    if meta['Project']['Biomass Module']=='Sawtooth':
#        dat['N']=vo['N'][it,:]
#        dat['N_M']=vo['N_M_Tot'][it,:]
#        dat['N_R']=vo['N_R'][it,:]
#        dat['C_M_Sim_Reg']=vo['C_M_Sim_Reg'][it,:]
#        dat['C_M_Sim_Fir']=vo['C_M_Sim_Fir'][it,:]
#        dat['C_M_Sim_Har']=vo['C_M_Sim_Har'][it,:]
#        dat['C_M_Sim_Ins']=vo['C_M_Sim_Ins'][it,:]
#        dat['C_M_Sim_Pat']=vo['C_M_Sim_Pat'][it,:]
#        dat['C_M_Sim_Win']=vo['C_M_Sim_Win'][it,:]
#        dat['TreeMean_A']=vo['TreeMean_A'][it,:]
#        dat['TreeMean_Csw']=vo['TreeMean_Csw'][it,:]
#        dat['TreeMean_Csw_G']=vo['TreeMean_Csw_G'][it,:]
#        dat['TreeMean_D']=vo['TreeMean_D'][it,:]
#        dat['TreeMean_H']=vo['TreeMean_H'][it,:]
    
    #--------------------------------------------------------------------------  
    # Save data
    #--------------------------------------------------------------------------  
    
    fout=meta['Paths']['Output Scenario'][iScn] + '\\Data_Scn' + cbu.FixFileNum(iScn) + \
        '_Ens' + cbu.FixFileNum(iEns) + '_Bat' + cbu.FixFileNum(iBat) + '.pkl'
    
    gu.opickle(fout,vo)
    
    #--------------------------------------------------------------------------
    # Save disturbance/management event chronology
    # If events are added on the fly, they will only be accessable if resaved.
    #--------------------------------------------------------------------------
    
    if (meta['Scenario'][iScn]['Harvest Status Historical']=='On') | (meta['Scenario'][iScn]['Harvest Status Future']=='On') | (meta['Scenario'][iScn]['Breakup Status']=='On'):
        
        # If it was input as compressed, output as re-compressed
        if 'idx' in vi['EC']:
            # Revise index
            vi['EC']['idx']=np.where(vi['EC']['ID_Type']>0)
            vi['EC']['ID_Type']=vi['EC']['ID_Type'][vi['EC']['idx']]
            vi['EC']['MortalityFactor']=vi['EC']['MortalityFactor'][vi['EC']['idx']]
            vi['EC']['GrowthFactor']=vi['EC']['GrowthFactor'][vi['EC']['idx']]
            vi['EC']['ID_GrowthCurve']=vi['EC']['ID_GrowthCurve'][vi['EC']['idx']]        
        
        fout=meta['Paths']['Input Scenario'][iScn] + '\\Modified_Events_Ens' + cbu.FixFileNum(iEns) + '_Bat' + cbu.FixFileNum(iBat) + '.pkl'
        
        gu.opickle(fout,vi['EC'])
    
    #--------------------------------------------------------------------------
    # Save diagnostics (turned off for now)
    #--------------------------------------------------------------------------
    
    flg=0
    if flg==1:
        #if (iBat==0) & (iScn==0):
    
        cnams=['Value','Variable']
        df=pd.DataFrame(columns=cnams)
    
        # Model version
        df.loc[0]=[meta['Paths']['Model Code'],'Version']    
    
        # Date
        now=datetime.datetime.now()
        df.loc[1]=[now,'Run date']
    
        # Simulation time
        #meta['t_Sim']=time.time()-meta['t_SimStart']
        #df.loc[2]=[meta['t_Sim']/60,'Simulation time (min)']
    
        # Check conservation of mass (stock change = NEBP)
        # *** only done for the first stand ***
        NPP=np.sum(vo['C_NPP'][:,0,0:7],axis=1)
        RH=np.sum(vo['C_RH'][:,0,7:16],axis=1)
        E=vo['C_E_FireAsCO2'][:,0]+vo['C_E_FireAsCO'][:,0]+vo['C_E_FireAsCH4'][:,0]+vo['C_E_FireAsN2O'][:,0]
        R=vo['C_ToMillMerch'][:,0]+vo['C_ToMillNonMerch'][:,0]
        x=np.sum(vo['C_Eco_Pools'][:,0,0:16],axis=1)
        y=np.sum(vo['C_Eco_Pools'][0,0,0:16],axis=0)+np.cumsum(NPP-RH-R-E)
        D_abs=np.mean(np.abs(y-x))
        D_rel=np.mean(np.abs(y-x)/np.maximum(0.000001,x)*100)
        df.loc[3]=[D_abs,'Mean absolute difference between stock change and NECB (MgC/ha)']
        df.loc[4]=[D_rel,'Mean relative difference between stock change and NECB (%)']    
        #plt.close(6)
        #plt.figure(6)
        #plt.plot(x)
        #plt.plot(y)
    
        # Year that slow soil carbon pool reaches dynamic equilibrium        
        def runningMean(x, N):
            y=np.zeros((len(x),))
            for ctr in range(len(x)):
                y[ctr]=np.sum(x[ctr:(ctr+N)])
            return y/N
    
        #plt.plot(runningMean(vo['C_Eco_Pools[:,0,15],150))
        #x=np.diff(vo['C_Eco_Pools'][:,0,15])/vo['C_Eco_Pools'][0:-1,0,15]*100
        #rmx=runningMean(x,150)
        #plt.figure(33)
        #plt.plot(rmx)
        #iDE=np.min(np.where(rmx<0.001))
        #yrDE=vi['tv'][iDE]
        #tDE=meta.N_t-iDE
        #df.loc[5]=[yrDE,'Equilibrium in SoilS reached at (year)']
        #df.loc[6]=[tDE,'Equlibrium in Soil S reached at (time before end of simulation)']
        
        # Save
        pthoutD=meta['Paths']['Output Scenario'][iScn] + '\\Diagnostics.xlsx'
        writer=pd.ExcelWriter(pthoutD)
        df.to_excel(writer,'Sheet1')
        writer.save()
    
    return vo_full

#%% Save outputs to model output statistics on the fly
                
def SaveOutputToMOS(meta,iScn,iEns):
    
    tv=np.arange(meta['Project']['Year Start Saving'],meta['Project']['Year End']+1,1)
    tv_full=np.arange(meta['Year Start'],meta['Project']['Year End']+1,1)
    
    if (iEns==0):
        # Create file
        mos=[None]*meta['Project']['N Scenario']
    else:
        # Import existing file
        mos=gu.ipickle(meta['Paths']['Project'] + '\\Outputs\\MOS_' + cbu.FixFileNum(iScn) + '.pkl')
        
    if iEns==0:
        
        # Initialize dictionaries for current scenario
        mos[iScn]={}
        
        mos[iScn]['v1']={}
        mos[iScn]['v1']['Mean']={}
        mos[iScn]['v1']['Sum']={}
        d1=cbu.LoadSingleOutputFile(meta,iScn,iEns,0)
        for k in d1.keys(): 
            if (k=='Year'):
                continue
            mos[iScn]['v1']['Mean'][k]={}
            mos[iScn]['v1']['Mean'][k]['Ensembles']=np.zeros((tv.size,meta['Project']['N Ensemble']))
            mos[iScn]['v1']['Sum'][k]={}
            mos[iScn]['v1']['Sum'][k]['Ensembles']=np.zeros((tv.size,meta['Project']['N Ensemble']))
            
        mos[iScn]['v2']={}
        mos[iScn]['v2']['Mean']={}
        mos[iScn]['v2']['Sum']={}
        d2=cbu.CalculateGHGBalance(d1,meta)        
        for k in d2[0][0].keys(): 
            if (k=='Year'):
                continue
            mos[iScn]['v2']['Mean'][k]={}
            mos[iScn]['v2']['Mean'][k]['Ensembles']=np.zeros((tv.size,meta['Project']['N Ensemble']))
            mos[iScn]['v2']['Sum'][k]={}
            mos[iScn]['v2']['Sum'][k]['Ensembles']=np.zeros((tv.size,meta['Project']['N Ensemble']))
            
        mos[iScn]['Area']={}
        for k in meta['LUT']['Dist'].keys():
            mos[iScn]['Area'][k]={}
            mos[iScn]['Area'][k]['Ensembles']=np.zeros((tv.size,meta['Project']['N Ensemble']))
        
    v1={}
    v2={}
    for iBat in range(meta['Project']['N Batch']):
            
        d1=cbu.LoadSingleOutputFile(meta,iScn,iEns,iBat)
        d2=cbu.CalculateGHGBalance(d1,meta)
            
        for k in d1.keys(): 
            if (k=='Year'):
                continue
            if iBat==0:
                v1[k]=np.sum(d1[k],axis=1)
            else:
                v1[k]=v1[k]+np.sum(d1[k],axis=1)
            
        for k in d2[0][0].keys():
            if k=='Year':
                continue
            if iBat==0:
                v2[k]=np.sum(d2[0][0][k],axis=1)
            else:
                v2[k]=v2[k]+np.sum(d2[0][0][k],axis=1)    
            
        # Import event chronology
        ec=gu.ipickle(meta['Paths']['Input Scenario'][iScn] + '\\Events_Ens' + cbu.FixFileNum(iEns) + '_Bat' + cbu.FixFileNum(iBat) + '.pkl')
            
        # Uncompress event chronology if it has been compressed
        if 'idx' in ec:
            ec=cbu.EventHistoryDecompress(meta,ec,iScn,iEns,iBat) 
            
        for iYr in range(tv_full.size):
            it=np.where(tv==tv_full[iYr])[0]
            if it.size==0:
                continue
            ID_Type0=ec['ID_Type'][iYr,:,:].flatten()
            u=np.unique(ID_Type0)
            for iU in range(u.size):
                if u[iU]==0:
                    continue
                id=cbu.lut_n2s(meta['LUT']['Dist'],u[iU])[0]
                ind=np.where(ID_Type0==u[iU])[0]
                mos[iScn]['Area'][id]['Ensembles'][it,iEns]=mos[iScn]['Area'][id]['Ensembles'][it,iEns]+ind.size
        del d1,d2,ec
        garc.collect()
        
    # Populate mos for each scenario
    for k in v1.keys():
        if k=='Year':
            continue
        mos[iScn]['v1']['Sum'][k]['Ensembles'][:,iEns]=v1[k].copy()
        mos[iScn]['v1']['Mean'][k]['Ensembles'][:,iEns]=v1[k].copy()/meta['Project']['N Stand']
    for k in v2.keys():
        if k=='Year':
            continue
        mos[iScn]['v2']['Sum'][k]['Ensembles'][:,iEns]=v2[k].copy()
        mos[iScn]['v2']['Mean'][k]['Ensembles'][:,iEns]=v2[k].copy()/meta['Project']['N Stand']
        
    # Save MOS
    gu.opickle(meta['Paths']['Project'] + '\\Outputs\\MOS_' + cbu.FixFileNum(iScn) + '.pkl',mos)
    
    # Delete full output data
    fnams=os.listdir(meta['Paths']['Output Scenario'][iScn])
    for i in range(len(fnams)):
        os.remove(meta['Paths']['Output Scenario'][iScn] + '\\' + fnams[i])
    
    return mos