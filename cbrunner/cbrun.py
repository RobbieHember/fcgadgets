
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
from fcgadgets.hardhat import geological as geologic

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

                # Path to temporary "working on batch..." file -> tells other instances
                pth_WorkingOnBatch=meta['Paths']['Project'] + '\\Outputs\\WorkingOnBatch_' + cbu.FixFileNum(iBat) + '.pkl'

                # Only proceed if the file does not exist running multiple instances
                if meta['Project']['Skip Completed Runs']=='On':
                    if os.path.exists(pth_WorkingOnBatch)==False:
                        gu.opickle(pth_WorkingOnBatch,[])
                        flag_WorkingOnBatch=1
                    else:
                        if flag_WorkingOnBatch==0:
                            print(pth_WorkingOnBatch)
                            continue

                # Report progress
                if (meta['Project']['Scenario Source']=='Spreadsheet'):
                    pass
                else:
                    print('Running Scenario ' + cbu.FixFileNum(iScn) + ', Ensemble ' + cbu.FixFileNum(iEns) + ', Batch ' + cbu.FixFileNum(iBat))

                # Initialize stands
                meta,vi,vo=InitializeStands(meta,iScn,iEns,iBat)

                # Track time
                t1=time.time()
                meta['Project']['Run Time Summary']['Stand initialization']=t1-t0

                # Set location-specific parameters
                meta,vi=PrepareParametersForBatch(meta,vi,iEns,iBat,iScn)

                # Track time
                t2=time.time()
                meta['Project']['Run Time Summary']['Set location-specific parameters']=t2-t1

                # Indices to ecosystem pools
                iEP=meta['Core']['iEP']

                # Try to start at a later date and adopt spinup from a previous run
                if (meta['Project']['N Ensemble']==1) | (iScn==0) & (iEns==0):
                    #Start from beginning
                    t_start=1
                else:
                    # Start at an advanced date using previous data for spinup period
                    t_start=np.where(meta['Year']==meta['Project']['Year Start Saving'])[0][0]

                    # Get output variables from first run of this batch
                    it=np.where(meta['Year']==meta['Project']['Year Start Saving']-1)[0]
                    for k in vo.keys():
                        if (k=='C_M_ByAgent'):
                            continue
                        if vo[k].size==0:
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

                            if vo[k].size==0:
                                continue

                            # Not a nested dictionary
                            vo[k][it,:]=0*vo[k][it,:]

                # Track time
                t3=time.time()

                # Biomass dynamics from Sawtooth
                if meta['Project']['Biomass Module']=='Sawtooth':
                    for iS in range(meta['Project']['Batch Size'][iBat]):
                        vo=annproc.BiomassFromSawtooth(iScn,iS,vi,vo,meta,iEP)

                # Track time
                t4=time.time()
                meta['Project']['Run Time Summary']['Running biomass dynamics from sawtooth']=t4-t3

                meta['Project']['Run Time Summary']['Biomass from GY model']=0.0
                meta['Project']['Run Time Summary']['DOM']=0.0
                meta['Project']['Run Time Summary']['Events']=0.0
                meta['Project']['Run Time Summary']['HWP']=0.0

                meta['Project']['Run Time Summary']['Test1']=0.0
                meta['Project']['Run Time Summary']['Test2']=0.0
                meta['Project']['Run Time Summary']['Test3']=0.0
                meta['Project']['Run Time Summary']['Test4']=0.0
                meta['Project']['Run Time Summary']['Test5']=0.0
                meta['Project']['Run Time Summary']['Test6']=0.0
                meta['Project']['Run Time Summary']['Test7']=0.0
                meta['Project']['Run Time Summary']['Test8']=0.0
                meta['Project']['Run Time Summary']['Test9']=0.0

                # Loop through time intervals (start in second time step)
                for iT in range(t_start,meta['Project']['N Time']):

                    # Biomass dynamics
                    if meta['Project']['Biomass Module']=='BatchTIPSY':

                        t0=time.time()
                        # Biomass from growth and yield model
                        vo=annproc.Biomass_FromGYM(iScn,iBat,iT,vi,vo,meta,iEP)
                        t1=time.time()
                        meta['Project']['Run Time Summary']['Biomass from GY model']=meta['Project']['Run Time Summary']['Biomass from GY model']+t1-t0

                    # Grassland dynamics
                    if (meta['Scenario'][iScn]['Grass Module Status']=='On') & (meta['Year'][iT]>=meta['Scenario'][iScn]['Grass Module Year Start']):
                        vo=annproc.Biomass_FromGrasses(iScn,iBat,iT,vi,vo,meta,iEP)

                    # Calculate annual dead organic matter dynamics
                    t0=time.time()
                    vo=annproc.DOM_LikeKurzetal2009(iT,iBat,vi,vo,iEP,meta)
                    t1=time.time()
                    meta['Project']['Run Time Summary']['DOM']=meta['Project']['Run Time Summary']['DOM']+t1-t0

                    # Calculate effects of disturbance and management
                    t0=time.time()
                    vo,vi=annproc.Events_FromTaz(iT,iScn,iEns,iBat,vi,vo,meta,iEP)
                    t1=time.time()
                    meta['Project']['Run Time Summary']['Events']=meta['Project']['Run Time Summary']['Events']+t1-t0

                    # Calculate products sector
                    if meta['Year'][iT]>=meta['Core']['HWP Year Start']:

                        # No need to run this before a certain date
                        t0=time.time()
                        vo=annproc.HWP_Update21(iT,iBat,vi,vo,meta)
                        t1=time.time()
                        meta['Project']['Run Time Summary']['HWP']=meta['Project']['Run Time Summary']['HWP']+t1-t0

                # Track time
                t5=time.time()
                meta['Project']['Run Time Summary']['Full annual loop']=t5-t4

                # Calculate fossil fuel emissions from operations
                vo=geologic.FossilFuelEmissions(meta,vi,vo)

                # Calculate substitution effects
                vo=geologic.SubstitutionEffects(meta,vi,vo)

                # Track time
                t6=time.time()
                meta['Project']['Run Time Summary']['Run geological']=t6-t5

                # Export simulation results to file
                #ExportSimulation(meta,vi,vo,iScn,iEns,iBat,iEP)
                if (meta['Project']['N Ensemble']==1) | (iScn==0) & (iEns==0):
                    # Only save output if it is the first instance of a batch
                    vo_full=ExportSimulation(meta,vi,vo,iScn,iEns,iBat,iEP)
                else:
                    ExportSimulation(meta,vi,vo,iScn,iEns,iBat,iEP)

                # Track time
                t7=time.time()
                meta['Project']['Run Time Summary']['Export results to file']=t7-t6

                # Delete 'working on' file
                #if meta['Project']['Skip Completed Runs']=='On':
                    #os.remove(pthWO)

                # Delete variables
                del vi,vo
                garc.collect()

                # Track time
                t8=time.time()
                meta['Project']['Run Time Summary']['Collect garbage']=t8-t7

    return meta

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
    meta['Project']['N Stand Batch']=vi['Inv']['ID_BECZ'].shape[1]

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
        #vi['GC']['Active']=vi['GC'][1].copy()
        vi['GC']['Active']=vi['GC'][1].copy().astype(float)*meta['GC']['Scale Factor']

        # Initialize an indicator of the active growth curve ID
        vi['GC']['ID_GCA']=np.ones(meta['Project']['Batch Size'][iBat])

        # Import growth curve 2
        vi['GC'][2]=gu.ipickle(meta['Paths']['Input Scenario'][iScn] + '\\GrowthCurve2_Bat' + cbu.FixFileNum(iBat) + '.pkl')

        # Import growth curve 3
        try:
            vi['GC'][3]=gu.ipickle(meta['Paths']['Input Scenario'][iScn] + '\\GrowthCurve3_Bat' + cbu.FixFileNum(iBat) + '.pkl')
        except:
            vi['GC'][3]=0

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
    vo['V_MerchLive']=np.zeros((m,n))
    vo['V_MerchDead']=np.zeros((m,n))
    vo['V_MerchTotal']=np.zeros((m,n))

    # Stemwood merch volume sent to mill
    vo['V_ToMillMerchLive']=np.zeros((m,n))
    vo['V_ToMillMerchDead']=np.zeros((m,n))
    vo['V_ToMillMerchTotal']=np.zeros((m,n))
    vo['V_ToMillNonMerch']=np.zeros((m,n))

    # Piece size
    vo['LogSizeEnhancement']=np.zeros((m,n))

    # Carbon density of ecosystem (Mg C ha-1)
    # -> 3-D matrix: Time x Stand x Carbon pool
    vo['C_Eco_Pools']=np.zeros((m,n,o))

    # Carbon density of products sector (Mg C ha-1)
    # -> 3-D matrix: Time x Stand x Carbon pool
    vo['C_Pro_Pools']=np.zeros((m,n,meta['Core']['N Pools Pro']))

    # Aggregate pools (Mg C ha-1) (polulated upon export)
    vo['C_Biomass_Tot']=np.array([])
    vo['C_Piled_Tot']=np.array([])
    vo['C_Litter_Tot']=np.array([])
    vo['C_DeadWood_Tot']=np.array([])
    vo['C_Soil_Tot']=np.array([])
    vo['C_Soil_OHorizon']=np.array([])
    vo['C_InUse_Tot']=np.array([])
    vo['C_DumpLandfill_Tot']=np.array([])
    vo['C_Buildings_Tot']=np.zeros((m,n))

    # Carbon flux densities (Mg C ha-1 yr-1)
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

    # Keep track of carbon transfers (for economics)
    vo['C_ToMillMerch']=np.zeros((m,n))
    vo['C_ToMillNonMerch']=np.zeros((m,n))
    vo['C_ToMillSnagStem']=np.zeros((m,n))
    vo['C_ToSlashpileBurnTot']=np.zeros((m,n))
    vo['C_ToSlashpileBurnNonMerch']=np.zeros((m,n))
    vo['C_ToLumber']=np.zeros((m,n))
    vo['C_ToPlywood']=np.zeros((m,n))
    vo['C_ToOSB']=np.zeros((m,n))
    vo['C_ToMDF']=np.zeros((m,n))
    vo['C_ToPaper']=np.zeros((m,n))
    vo['C_ToPowerFacilityDom']=np.zeros((m,n))
    vo['C_ToPowerFacilityExport']=np.zeros((m,n))
    vo['C_ToPowerGrid']=np.zeros((m,n))
    vo['C_ToPelletExport']=np.zeros((m,n))
    vo['C_ToPelletDomRNG']=np.zeros((m,n))
    vo['C_ToPelletDomGrid']=np.zeros((m,n))
    vo['C_ToFirewoodDom']=np.zeros((m,n))
    vo['C_ToFirewoodExport']=np.zeros((m,n))
    vo['C_ToLogExport']=np.zeros((m,n))

    # Emissions from wildfire and open burning (will be deleted upon export)
    vo['C_E_OpenBurningAsCO2']=np.zeros((m,n))
    vo['C_E_OpenBurningAsCH4']=np.zeros((m,n))
    vo['C_E_OpenBurningAsCO']=np.zeros((m,n))
    vo['C_E_OpenBurningAsN2O']=np.zeros((m,n))

    vo['C_E_WildfireAsCO2']=np.zeros((m,n))
    vo['C_E_WildfireAsCH4']=np.zeros((m,n))
    vo['C_E_WildfireAsCO']=np.zeros((m,n))
    vo['C_E_WildfireAsN2O']=np.zeros((m,n))

    vo['C_Coal']=np.zeros((m,n))
    vo['C_Oil']=np.zeros((m,n))
    vo['C_Gas']=np.zeros((m,n))
    vo['C_Limestone']=np.zeros((m,n))

    # LULUCF Sector: Net ecosystem exchange (polulated in load results)
    vo['E_CO2e_LULUCF_NEE']=np.array([])

    # LULUCF Sector: Net ecosystem production (polulated upon export)
    vo['E_CO2e_LULUCF_Wildfire']=np.array([])

    # LULUCF Sector: Net ecosystem production (polulated upon export)
    vo['E_CO2e_LULUCF_OpenBurning']=np.array([])

    # LULUCF Sector: Denitrification
    vo['E_CO2e_LULUCF_Denit']=np.zeros((m,n))

    # LULUCF Sector: valatilization and deposition of N
    vo['E_CO2e_LULUCF_Other']=np.zeros((m,n))

    # LULUCF Sector: RH from dumps and landfills
    vo['E_CO2e_LULUCF_HWP']=np.zeros((m,n))

    # Energy - Stationary Combustion Sector
    vo['E_CO2e_ESC_Bioenergy']=np.zeros((m,n))
    vo['E_CO2e_ESC_BioenergyPowerFacilityDom']=np.zeros((m,n))
    vo['E_CO2e_ESC_BioenergyPowerFacilityExport']=np.zeros((m,n))
    vo['E_CO2e_ESC_BioenergyPowerGrid']=np.zeros((m,n))
    vo['E_CO2e_ESC_BioenergyPelletExport']=np.zeros((m,n))
    vo['E_CO2e_ESC_BioenergyPelletDomGrid']=np.zeros((m,n))
    vo['E_CO2e_ESC_BioenergyPelletDomRNG']=np.zeros((m,n))
    vo['E_CO2e_ESC_BioenergyFirewoodDom']=np.zeros((m,n))
    vo['E_CO2e_ESC_BioenergyFirewoodExport']=np.zeros((m,n))

    vo['E_CO2e_ESC_OperForBurnCoal']=np.zeros((m,n))
    vo['E_CO2e_ESC_OperForBurnOil']=np.zeros((m,n))
    vo['E_CO2e_ESC_OperForBurnGas']=np.zeros((m,n))

    # Energy - Transporation Sector
    vo['E_CO2e_ET_OperForBurnCoal']=np.zeros((m,n))
    vo['E_CO2e_ET_OperForBurnOil']=np.zeros((m,n))
    vo['E_CO2e_ET_OperForBurnGas']=np.zeros((m,n))

    # Industrial Produciton and Product Use Sector
    vo['E_CO2e_IPPU_OperForBurningCoal']=np.zeros((m,n))
    vo['E_CO2e_IPPU_OperForBurningOil']=np.zeros((m,n))
    vo['E_CO2e_IPPU_OperForBurningGas']=np.zeros((m,n))

    if meta['Project']['Biomass Module']=='Sawtooth':

        meta['Core']['Sawtooth']={}
        meta['Core']['Sawtooth']['DBH Classes']=np.linspace(0,100,50)
        vo['DBH_Class']=np.zeros((m,n,meta['Core']['Sawtooth']['DBH Classes'].size))

        # Stand density
        vo['N']=np.zeros((m,n))

        # Change in stand density (stems ha-1 yr-1)
        vo['N_R']=np.zeros((m,n))
        vo['N_M_Tot']=np.zeros((m,n))
        vo['N_M_Reg']=np.zeros((m,n))

        # Mean of tree attributes
        vo['TreeMean_A']=np.zeros((m,n))
        vo['TreeMean_H']=np.zeros((m,n))
        vo['TreeMean_D']=np.zeros((m,n))
        vo['TreeMean_Csw']=np.zeros((m,n))
        vo['TreeMean_Csw_G']=np.zeros((m,n))

        # Needed to supply affected carbon to disturbance module
        vo['C_M_Tot']=np.zeros((m,n,o))

        # Needed to track dead merch volume
        vo['V_Merch_M']=np.zeros((m,n))

    # Atmosphere
    vo['Atm_CO2_In']=np.zeros((m,n))
    vo['Atm_CO2_Out']=np.zeros((m,n))
    vo['Atm_CH4_In']=np.zeros((m,n))
    vo['Atm_CH4_Out']=np.zeros((m,n))
    vo['Atm_N2O_In']=np.zeros((m,n))
    vo['Atm_N2O_Out']=np.zeros((m,n))

    #--------------------------------------------------------------------------
    # Specify a fast-track period that spans the spin-up period for use of Sawtooth
    # module
    # *** Not working ***
    #--------------------------------------------------------------------------

#    if meta['Project']['Biomass Module']=='Sawtooth':
#
#        meta['Project']['SpinupSpanFastTrack']=[None]*meta['Project']['N Scenario']
#
#        for iScn in range(meta['Project']['N Scenario']):
#
#            # Get spinup event years
#            ivl_spin=meta['Project']['Spinup Disturbance Return Inverval']
#            YearRef=meta['Scenario'][iScn]['Year1_DisFromInv']
#            AgeRef=meta['Scenario'][iScn]['Age1_DisFromInv']
#            if AgeRef>=0:
#                Year=np.arange(YearRef-AgeRef-100*ivl_spin,YearRef-AgeRef+ivl_spin,ivl_spin)
#            else:
#                Year1=meta['Project']['Year Start']+ivl_spin
#                Year2=meta['Project']['Spinup Year End']
#                Year=np.arange(Year1,Year2+1,meta['Project']['Spinup Disturbance Return Inverval'])
#
#            meta['Project']['SpinupSpanFastTrack'][iScn]={}
#            meta['Project']['SpinupSpanFastTrack'][iScn]['Start']=Year[3]
#            meta['Project']['SpinupSpanFastTrack'][iScn]['End']=Year[-2]

    #--------------------------------------------------------------------------
    # Configure batch-specific setttings
    #--------------------------------------------------------------------------

    # Nutrient application response yearly counter
    meta['Nutrient Management']['ResponseCounter']=np.zeros(meta['Project']['Batch Size'][iBat])

    # Initialize flag for fixing negative net growth. When TIPSY yields negative
    # net growth, the fluxes of gross growth and mortality need adjustment.
    # This flag helps achieve that.
    meta['FlagNegNetGrowth']=np.zeros(meta['Project']['Batch Size'][iBat])
    meta['G_Net_PriorToBreakup']=np.zeros((meta['Project']['Batch Size'][iBat],7))

    #--------------------------------------------------------------------------
    # Generate random numbers for on-the-fly disturbance types
    #--------------------------------------------------------------------------

    if (meta['Scenario'][iScn]['Harvest Status Historical']=='On') | (meta['Scenario'][iScn]['Harvest Status Future']=='On'):
        if meta['Project']['Use Frozen Ensembles']=='Off':
            rn=gu.ipickle(meta['Paths']['Project'] + '\\Inputs\\Ensembles\\RandomNumbers_Harvest_Ens' + cbu.FixFileNum(iEns) + '_Bat' + cbu.FixFileNum(iBat) + '.pkl')
        else:
            rn=gu.ipickle(meta['Project']['Use Frozen Ensembles'] + '\\RandomNumbers_Harvest_Ens' + cbu.FixFileNum(iEns) + '_Bat' + cbu.FixFileNum(iBat) + '.pkl')
        rn=rn.astype(float)
        rn=rn*meta['Project']['On the Fly']['Random Numbers']['Scale Factor']
        meta['Project']['On the Fly']['Random Numbers']['Harvest']=rn.copy()

    if (meta['Scenario'][iScn]['Breakup Status']=='On'):
        if meta['Project']['Use Frozen Ensembles']=='Off':
            rn=gu.ipickle(meta['Paths']['Project'] + '\\Inputs\\Ensembles\\RandomNumbers_Breakup_Ens' + cbu.FixFileNum(iEns) + '_Bat' + cbu.FixFileNum(iBat) + '.pkl')
        else:
            rn=gu.ipickle(meta['Project']['Use Frozen Ensembles'] + '\\RandomNumbers_Breakup_Ens' + cbu.FixFileNum(iEns) + '_Bat' + cbu.FixFileNum(iBat) + '.pkl')
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

def PrepareParametersForBatch(meta,vi,iEns,iBat,iScn):

    #--------------------------------------------------------------------------
    # Populate final parameters with initial best estimates
    #--------------------------------------------------------------------------

    meta['Param']['BEV']=copy.deepcopy(meta['Param']['BE'])

    #--------------------------------------------------------------------------
    # Add error variance to parameters
    #--------------------------------------------------------------------------

    if meta['Project']['Uncertainty Status Biomass Turnover']=='On':
        for k in meta['Param']['BE']['Biomass Turnover'].keys():
            meta['Param']['BEV']['Biomass Turnover'][k]=meta['Param']['By Ensemble'][iEns]['Biomass Turnover'][k]

    if meta['Project']['Uncertainty Status Decomposition']=='On':
        for k in meta['Param']['BE']['Decomp'].keys():
            meta['Param']['BEV']['Decomp'][k]=meta['Param']['By Ensemble'][iEns]['Decomp'][k]

    if meta['Project']['Uncertainty Status Inter Pool Fluxes']=='On':
        for k in meta['Param']['BE']['Inter Pool Fluxes'].keys():
            meta['Param']['BEV']['Inter Pool Fluxes'][k]=meta['Param']['By Ensemble'][iEns]['Inter Pool Fluxes'][k]

    if meta['Project']['Uncertainty Status Harvest Utilization']=='On':
        EventList=['Harvest','Harvest Salvage']
        VariableList=['BiomassMerch_Removed','BiomassNonMerch_Removed','Snags_Removed', \
                      'BiomassMerch_Piled','BiomassNonMerch_Piled','Snags_Piled', \
                      'BiomassMerch_LeftOnSite','BiomassNonMerch_LeftOnSite','Snags_LeftOnSite']
        for Event in EventList:
            ID_Type=meta['LUT']['Dist'][Event]
            for Variable in VariableList:
                meta['Param']['BE']['Dist'][ID_Type][Variable]=meta['Param']['By Ensemble'][iEns]['Dist'][ID_Type][Variable]

    if meta['Project']['Uncertainty Status Substitution']=='On':
        for k in meta['Param']['By Ensemble'][iEns]['Substitution'].keys():
            meta['Param']['BEV']['Substitution'][k]=meta['Param']['By Ensemble'][iEns]['Substitution'][k]

    if meta['Project']['Uncertainty Status Nutrient Application']=='On':
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
    # Fate of Felled Material
    #--------------------------------------------------------------------------

    # Initialize
    meta['Param']['BEV']['Felled Fate']={}
    for k in meta['Param']['BE']['Felled Fate']['BaseCase']['Coast'].keys():
        meta['Param']['BEV']['Felled Fate'][k]=np.zeros( (meta['Param']['BE']['Felled Fate']['Year'].size,meta['Project']['Batch Size'][iBat]) )

    # Papulate batch-specific parameters
    if meta['Project']['Scenario Source']=='Portfolio':

        indBat=cbu.IndexToBatch(meta,iBat)

        # Isolate felled fate scenario names within this batch
        Scenario=meta['Project']['Portfolio']['Felled Fate Scenario'][indBat]

        # Isolate region
        Region=meta['Project']['Portfolio']['Region Code'][indBat]

        # Unique scenario and region (must be converted to string)
        SR=np.column_stack((Scenario,Region))
        SR=SR.astype(str)
        u=np.unique(SR,axis=0)

        for iU in range(u.shape[0]):

            # Index to each scenario
            scn=u[iU,0]
            reg=u[iU,1]
            ind=np.where( (Scenario==scn) & (Region==reg) )[0]

            for k in meta['Param']['BE']['Felled Fate'][scn][reg].keys():
                x=meta['Param']['BE']['Felled Fate'][scn][reg][k]
                for i in range(ind.size):
                    meta['Param']['BEV']['Felled Fate'][k][:,ind[i]]=x

    elif meta['Project']['Scenario Source']=='Spreadsheet':

        Scenario=meta['Scenario'][meta['iScn']]['Felled Fate Scenario']
        Region=meta['Scenario'][meta['iScn']]['Region Code']
        for k in meta['Param']['BE']['Felled Fate'][Scenario][Region].keys():
            x=meta['Param']['BE']['Felled Fate'][Scenario][Region][k]
            x=np.reshape(x,(-1,1))
            x=np.tile(x,(1,meta['Project']['Batch Size'][iBat]))
            meta['Param']['BEV']['Felled Fate'][k]=x

    elif meta['Project']['Scenario Source']=='Script':

        # Index to batch
        indBat=cbu.IndexToBatch(meta,iBat)

        # Isolate felled fate scenario names within this batch
        Scenario=meta['Scenario'][meta['iScn']]['Felled Fate Scenario']

        for reg in meta['LUT']['Region'].keys():
            ind=np.where( vi['Inv']['Region Code'][0,:]==meta['LUT']['Region'][reg] )[0]
            for k in meta['Param']['BE']['Felled Fate'][Scenario][reg].keys():
                x=meta['Param']['BE']['Felled Fate'][Scenario][reg][k]
                for i in range(ind.size):
                    meta['Param']['BEV']['Felled Fate'][k][:,ind[i]]=x

        # Override regional parameters for stands that have land use = energy production
        if meta['Scenario'][iScn]['Land Surface Scenario']!='None':
            iEnergy=np.where(vi['Inv']['LSC']['Use']==meta['LUT']['LSC']['Use']['Energy Production'])
            if iEnergy[0].size>0:
                for k in meta['Param']['BEV']['Felled Fate'].keys():
                    meta['Param']['BEV']['Felled Fate'][k][iEnergy]=meta['Param']['BE']['Felled Fate'][Scenario]['Energy Production'][k][0]

    #--------------------------------------------------------------------------
    # Removed Fate
    #--------------------------------------------------------------------------

    # Initialize
    meta['Param']['BEV']['Removed Fate']={}
    for k in meta['Param']['BE']['Removed Fate']['BaseCase']['Coast'].keys():
        meta['Param']['BEV']['Removed Fate'][k]=np.zeros( (meta['Param']['BE']['Removed Fate']['Year'].size,meta['Project']['Batch Size'][iBat]) )

    # Papulate batch-specific parameters
    if meta['Project']['Scenario Source']=='Portfolio':

        indBat=cbu.IndexToBatch(meta,iBat)

        # Isolate Removal Fate scenario names within this batch
        Scenario=meta['Project']['Portfolio']['Removed Fate Scenario'][indBat]

        # Isolate region
        Region=meta['Project']['Portfolio']['Region Code'][indBat]

        # Unique scenario and region (must be converted to string)
        SR=np.column_stack( (Scenario,Region) )
        SR=SR.astype(str)
        u=np.unique(SR,axis=0)

        for iU in range(u.shape[0]):

            # Index to each scenario
            scn=u[iU,0]
            reg=u[iU,1]
            ind=np.where( (Scenario==scn) & (Region==reg) )[0]

            for k in meta['Param']['BE']['Removed Fate'][scn][reg].keys():
                x=meta['Param']['BE']['Removed Fate'][scn][reg][k]
                for i in range(ind.size):
                    meta['Param']['BEV']['Removed Fate'][k][:,ind[i]]=x

    elif meta['Project']['Scenario Source']=='Spreadsheet':

        Scenario=meta['Scenario'][meta['iScn']]['Removed Fate Scenario']
        Region=meta['Scenario'][meta['iScn']]['Region Code']
        for k in meta['Param']['BE']['Removed Fate'][Scenario][Region].keys():
            x=meta['Param']['BE']['Removed Fate'][Scenario][Region][k]
            x=np.reshape(x,(-1,1))
            x=np.tile(x,(1,meta['Project']['Batch Size'][iBat]))
            meta['Param']['BEV']['Removed Fate'][k]=x

    elif meta['Project']['Scenario Source']=='Script':

        # Index to batch
        indBat=cbu.IndexToBatch(meta,iBat)

        # Isolate Removal Fate scenario names within this batch
        Scenario=meta['Scenario'][meta['iScn']]['Removed Fate Scenario']

        for reg in meta['LUT']['Region'].keys():
            ind=np.where( vi['Inv']['Region Code'][0,:]==meta['LUT']['Region'][reg] )[0]
            for k in meta['Param']['BE']['Removed Fate'][Scenario][reg].keys():
                x=meta['Param']['BE']['Removed Fate'][Scenario][reg][k]
                for i in range(ind.size):
                    meta['Param']['BEV']['Removed Fate'][k][:,ind[i]]=x

        # Override regional parameters for stands that have land use = energy production
        if meta['Scenario'][iScn]['Land Surface Scenario']!='None':
            iEnergy=np.where(vi['Inv']['LSC']['Use']==meta['LUT']['LSC']['Use']['Energy Production'])
            if iEnergy[0].size>0:
                for k in meta['Param']['BEV']['Removed Fate'].keys():
                    meta['Param']['BEV']['Removed Fate'][k][iEnergy]=meta['Param']['BE']['Removed Fate'][Scenario]['Energy Production'][k][0]

    #--------------------------------------------------------------------------
    # Harvested Wood Products - static
    #--------------------------------------------------------------------------

    for i in range(len(meta['Param']['BEV']['HWP']['Raw'])):

        Name=meta['Param']['BEV']['HWP']['Raw']['Name'].iloc[i]
        Value=meta['Param']['BEV']['HWP']['Raw']['Coast'].iloc[i]

        # Convert half life to turnover rate
        if Name[-2:]=='hl':
            Name=Name[0:-2] + 'tr'
            Value=1-np.exp(-np.log(2)/Value)

        meta['Param']['BEV']['HWP'][Name]=Value*np.ones(meta['Project']['Batch Size'][iBat])

    #--------------------------------------------------------------------------
    # Override mill transfers (for salvage logging)
    # Shift transfers to pulp in salvage logging projects
    #--------------------------------------------------------------------------

    if 'Salvage Mill Transfers' in meta['Scenario'][iScn]:
        if meta['Scenario'][iScn]['Salvage Mill Transfers']=='On':
            y=meta['Param']['BEV']['HWP']['RemovedMerchToSawMill'].copy()
            meta['Param']['BEV']['HWP']['RemovedMerchToSawMill']=0.0*y
            meta['Param']['BEV']['HWP']['RemovedMerchToPulpMill']=meta['Param']['BEV']['HWP']['RemovedMerchToPulpMill']+1.0*y

            y=meta['Param']['BEV']['HWP']['RemovedSnagStemToSawMill'].copy()
            meta['Param']['BEV']['HWP']['RemovedSnagStemToSawMill']=0.0*y
            meta['Param']['BEV']['HWP']['RemovedSnagStemToPulpMill']=meta['Param']['BEV']['HWP']['RemovedSnagStemToPulpMill']+1.0*y

    #--------------------------------------------------------------------------
    # Harvested Wood Products - End Uses
    #--------------------------------------------------------------------------

    # Initialize
    meta['Param']['BEV']['HWP End Use']={}
    for k in meta['Param']['BE']['HWP End Use']['BaseCase']['Coast'].keys():
        meta['Param']['BEV']['HWP End Use'][k]=np.zeros( (meta['Param']['BE']['HWP End Use']['Year'].size,meta['Project']['Batch Size'][iBat]) )

    # Papulate batch-specific parameters
    if meta['Project']['Scenario Source']=='Portfolio':

        indBat=cbu.IndexToBatch(meta,iBat)

        # Isolate scenario names within this batch
        Scenario=meta['Project']['Portfolio']['HWP End Use Scenario'][indBat]

        # Isolate region
        Region=meta['Project']['Portfolio']['Region Code'][indBat]

        # Unique scenario and region (must be converted to string)
        SR=np.column_stack( (Scenario,Region) )
        SR=SR.astype(str)
        u=np.unique(SR,axis=0)

        for iU in range(u.shape[0]):

            # Index to each scenario
            scn=u[iU,0]
            reg=u[iU,1]
            ind=np.where( (Scenario==scn) & (Region==reg) )[0]

            for k in meta['Param']['BE']['HWP End Use'][scn][reg].keys():
                x=meta['Param']['BE']['HWP End Use'][scn][reg][k]
                x=np.reshape(x,(-1,1))
                x=np.tile(x,(1,ind.size))
                meta['Param']['BEV']['HWP End Use'][k][:,ind]=x

    elif meta['Project']['Scenario Source']=='Spreadsheet':

        Scenario=meta['Scenario'][meta['iScn']]['HWP End Use Scenario']
        Region=meta['Scenario'][meta['iScn']]['Region Code']
        for k in meta['Param']['BE']['HWP End Use'][Scenario][Region].keys():
            x=meta['Param']['BE']['HWP End Use'][Scenario][Region][k]
            x=np.reshape(x,(-1,1))
            x=np.tile(x,(1,meta['Project']['Batch Size'][iBat]))
            meta['Param']['BEV']['HWP End Use'][k]=x

    else:

        # Index to batch
        indBat=cbu.IndexToBatch(meta,iBat)

        # Isolate Removal Fate scenario names within this batch
        Scenario=meta['Scenario'][meta['iScn']]['HWP End Use Scenario']

        for reg in meta['LUT']['Region'].keys():
            ind=np.where( vi['Inv']['Region Code'][0,:]==meta['LUT']['Region'][reg] )[0]
            for k in meta['Param']['BE']['HWP End Use'][Scenario][reg].keys():
                x=meta['Param']['BE']['HWP End Use'][Scenario][reg][k]
                for i in range(ind.size):
                    meta['Param']['BEV']['HWP End Use'][k][:,ind[i]]=x

    #--------------------------------------------------------------------------
    # Populate custom harvest parameters with those supplied for project
    # *** Retired ***
    #--------------------------------------------------------------------------

    #    if 'Harvest Custom' in meta:
    #
    #        for iHC in range(10):
    #
    #            if int(iHC+1) in meta['Harvest Custom']:
    #
    #                for k in meta['Harvest Custom'][int(iHC+1)].keys():
    #
    #                    if k[0:7]!='Removed':
    #                        id=meta['LUT']['Dist']['Harvest Custom ' + str(int(iHC+1))]
    #                        meta['Param']['BEV']['Dist'][id][k]=meta['Harvest Custom'][int(iHC+1)][k]/100

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

            if vo[k].size==0:
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
    # Emissions from wildfire
    #--------------------------------------------------------------------------

    # Carbon dioxide flux (tCO2/ha/yr)
    E_CO2=bB['Ratio_CO2_to_C']*vo['C_E_WildfireAsCO2']

    # Carbon monoxide flux (tCO/ha/yr)
    E_CO=bB['Ratio_CO_to_C']*vo['C_E_WildfireAsCO']

    # Methan flux *(tCH4/ha/yr)
    E_CH4=bB['Ratio_CH4_to_C']*vo['C_E_WildfireAsCH4']

    # Nitrous oxide flux (tN2O/ha/yr)
    E_N2O=bB['EF_N2O_fromCO2']*E_CO2

    # Convert fluxes to CO2e using global warming potential estimates
    CO2e_E_AsCO2=1*E_CO2
    CO2e_E_AsCH4=bB['GWP_CH4_AR5']*E_CH4
    CO2e_E_AsCO=bB['GWP_CO_AR5']*E_CO
    CO2e_E_AsN2O=bB['GWP_N2O_AR5']*E_N2O

    vo['E_CO2e_LULUCF_Wildfire']=CO2e_E_AsCO2+CO2e_E_AsCH4+CO2e_E_AsCO+CO2e_E_AsN2O

    # Add to atmosphere
    vo['Atm_CO2_In']=vo['Atm_CO2_In']+E_CO2
    vo['Atm_CH4_In']=vo['Atm_CH4_In']+E_CH4
    vo['Atm_N2O_In']=vo['Atm_N2O_In']+E_N2O

    #--------------------------------------------------------------------------
    # Emissions from open burning
    #--------------------------------------------------------------------------

    # Carbon dioxide flux (tCO2/ha/yr)
    E_CO2=bB['Ratio_CO2_to_C']*vo['C_E_OpenBurningAsCO2']

    # Carbon monoxide flux (tCO/ha/yr)
    E_CO=bB['Ratio_CO_to_C']*vo['C_E_OpenBurningAsCO']

    # Methan flux *(tCH4/ha/yr)
    E_CH4=bB['Ratio_CH4_to_C']*vo['C_E_OpenBurningAsCH4']

    # Nitrous oxide flux (tN2O/ha/yr)
    E_N2O=bB['EF_N2O_fromCO2']*E_CO2

    # Convert fluxes to CO2e using global warming potential estimates
    CO2e_E_AsCO2=1*E_CO2
    CO2e_E_AsCH4=bB['GWP_CH4_AR5']*E_CH4
    CO2e_E_AsCO=bB['GWP_CO_AR5']*E_CO
    CO2e_E_AsN2O=bB['GWP_N2O_AR5']*E_N2O

    vo['E_CO2e_LULUCF_OpenBurning']=CO2e_E_AsCO2+CO2e_E_AsCH4+CO2e_E_AsCO+CO2e_E_AsN2O

    # Add to atmosphere
    vo['Atm_CO2_In']=vo['Atm_CO2_In']+E_CO2
    vo['Atm_CH4_In']=vo['Atm_CH4_In']+E_CH4
    vo['Atm_N2O_In']=vo['Atm_N2O_In']+E_N2O

    #--------------------------------------------------------------------------
    # Delete unnecessary fire emission variables
    #--------------------------------------------------------------------------

    del vo['C_E_WildfireAsCO2']
    del vo['C_E_WildfireAsCO']
    del vo['C_E_WildfireAsCH4']
    del vo['C_E_WildfireAsN2O']

    del vo['C_E_OpenBurningAsCO2']
    del vo['C_E_OpenBurningAsCO']
    del vo['C_E_OpenBurningAsCH4']
    del vo['C_E_OpenBurningAsN2O']

    #--------------------------------------------------------------------------
    # Sum pools and fluxes if "Save Biomass Pools" is Off
    # These variables only need to be saved to file if the indiviudal pools are
    # not being stored.
    #--------------------------------------------------------------------------

    # Calculate carbon content of dead wood, organic and mineral soil horizons following Shaw et al. (2017)
    vo['C_Soil_OHorizon']=vo['C_Eco_Pools'][:,:,iEP['LitterVF']]+vo['C_Eco_Pools'][:,:,iEP['LitterS']]

    if meta['Project']['Save Biomass Pools']!='On':

        # Aggregate pools
        vo['C_Biomass_Tot']=np.sum(vo['C_Eco_Pools'][:,:,iEP['BiomassTotal']],axis=2)
        vo['C_Piled_Tot']=np.sum(vo['C_Eco_Pools'][:,:,iEP['Piled']],axis=2)
        vo['C_Litter_Tot']=np.sum(vo['C_Eco_Pools'][:,:,iEP['Litter']],axis=2)
        vo['C_DeadWood_Tot']=np.sum(vo['C_Eco_Pools'][:,:,iEP['DeadWood']],axis=2)
        vo['C_Soil_Tot']=np.sum(vo['C_Eco_Pools'][:,:,iEP['Soil']],axis=2)
        vo['C_InUse_Tot']=np.sum(vo['C_Pro_Pools'][:,:,meta['Core']['iPP']['InUse']],axis=2)
        vo['C_Buildings_Tot']=np.sum(vo['C_Pro_Pools'][:,:,meta['Core']['iPP']['Buildings'] ],axis=2)
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
        del vo['C_M_Reg']
        del vo['C_LF']
        del vo['C_RH']

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
        if (k=='E_CO2e_LULUCF_HWP') | (k=='E_CO2e_ESC_Comb') | (k=='E_CO2e_ET_Comb') | (k=='E_CO2e_IPPU_Comb'):
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
