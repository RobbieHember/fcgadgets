
import numpy as np
from fcgadgets.cbrunner import cbrun_utilities as cbu

#%% Nutrient application effects

def UpdateStatus(meta,pNam,vi,vo,iT,comp):

    # Exctract parameters
    bNA=meta['Param']['BEV']['Nutrient Management']

    # Pull out index to applications
    iApp=meta['Modules']['Nutrient Management']['iApplication']

    if comp=='UpdateCounter':

        meta['Modules']['Nutrient Management']['ResponseCounter'][iApp]=meta['Modules']['Nutrient Management']['ResponseCounter'][iApp]+1

    elif (comp=='AbovegroundNetGrowth') & (meta[pNam]['Project']['Nutrient Application Module']=='cbrunner'):

        # Think about revising to be based on N
        # r_NonUniform_Spatial_Distribution

        #----------------------------------------------------------------------
        # Start response counter
        #----------------------------------------------------------------------

        meta['Modules']['Nutrient Management']['ResponseCounter'][iApp]=meta['Modules']['Nutrient Management']['ResponseCounter'][iApp]+1

        #----------------------------------------------------------------------
        # Update log-size enhancement factor
        #----------------------------------------------------------------------

        vo['LogSizeEnhancement'][iT:,iApp]=vo['LogSizeEnhancement'][iT:,iApp]+1

        #----------------------------------------------------------------------
        # Adjust net growth on treatment area
        #----------------------------------------------------------------------

        # Response ratio of stemwood (before wood density effect)
        rrS=bNA['r_Stemwood']

        # Wood density effect
        rr_wd=bNA['r_WoodDensity']
        f_wd=rr_wd-1
        f_wd=np.array([f_wd,f_wd,0,0,0])

        # Response relative to stemwood
        # Can't do roots here because they are not tracked in growth curve dict
        rXS=np.array([1,
                1,
                bNA['Ratio_Foliage_to_Stemwood'],
                bNA['Ratio_Branch_to_Stemwood'],
                bNA['Ratio_Bark_to_Stemwood']])

        # Biomass response ratios (after wood density effect accounted for)
        rr=1+(rrS-1)*rXS+f_wd

        # Append biomass response with volume response
        rr=np.append(rr,rrS)

        # Get age vector for GC's
        A_gc=np.arange(0,meta['Modules']['GYM']['BatchTIPSY Maximum Age']+1,1)

        # Extract active growth curves, apply scale factor
        #GCA_SP=vi['GC']['Active'][:,iApp,:].copy().astype(float)*meta['Modules']['GYM']['Scale Factor']
        GCA_SP=vi['GC']['Active'][:,iApp,:].copy()#.astype(float)*meta['Modules']['GYM']['Scale Factor']

        for iStand in range(iApp.size):

            AgeAtApplication=int(vo['A'][iT,iApp[iStand]])

            # Index to the response period that will be alterred
            iResponse=np.where( (A_gc>=AgeAtApplication) & (A_gc<AgeAtApplication+bNA['ResponseDuration']) )[0]

            # This will crash if an application occurs within 10 years of the maximum
            # TIPSY age curve (eg 200 years) -> adjust response period so that it
            # does not crash
            if iResponse.size==int(bNA['ResponseDuration']):

                for iRR in range(rr.size):
                    GCA_SP[iResponse,iStand,iRR]=rr[iRR]*GCA_SP[iResponse,iStand,iRR]

            else:
                em='Error: Nutrient application not implemented - stand age exceeds max age of growth curves.'
                print(em)
                #print(AgeAtApplication)

        # Re-applly scalefactor
        #GCA_SP=GCA_SP/meta['Modules']['GYM']['Scale Factor']

        # Convert to 16-bit integers
        #GCA_SP=GCA_SP.astype('int16')

        # Repopulate in input variable dictionary
        vi['GC']['Active'][:,iApp,:]=GCA_SP

    elif (comp=='BelowgroundNetGrowth') & (meta[pNam]['Project']['Nutrient Application Module']=='cbrunner'):

        #----------------------------------------------------------------------
        # Adjust root net growth
        #----------------------------------------------------------------------

        iEP=meta['Core']['iEP']

        # Responses
        rrRC=1+(bNA['r_Stemwood']-1)*bNA['Ratio_RootC_to_Stemwood']
        rrRF=1+(bNA['r_Stemwood']-1)*bNA['Ratio_RootF_to_Stemwood']

        # Calculate net growth of roots from change in pools
        #Gnet_RC=vo['C_Eco_Pools'][iT,:,iEP['RootCoarse']]-vo['C_Eco_Pools'][iT-1,:,iEP['RootCoarse']]
        #Gnet_RF=vo['C_Eco_Pools'][iT,:,iEP['RootFine']]-vo['C_Eco_Pools'][iT-1,:,iEP['RootFine']]

        # Remove the false stimulus that arises from estimating roots from AG
        # biomass
        vo['C_G_Net'][iT,iApp,iEP['RootCoarse']]=0.84*vo['C_G_Net'][iT,iApp,iEP['RootCoarse']]
        vo['C_G_Net'][iT,iApp,iEP['RootFine']]=0.84*vo['C_G_Net'][iT,iApp,iEP['RootFine']]
        #Gnet_RC[iApp]=0.65*(2-bNA['r_Stemwood'])*Gnet_RC[iApp]
        #Gnet_RF[iApp]=0.65*(2-bNA['r_Stemwood'])*Gnet_RF[iApp]

        # Stimulate root growth from N application
        #Gnet_RC[iApp]=rrRC*Gnet_RC[iApp]
        #Gnet_RF[iApp]=rrRF*Gnet_RF[iApp]

        # Recalculate current-year root biomass where stimulated by N application
        #vo['C_Eco_Pools'][iT,iApp,iEP['RootCoarse']]=vo['C_Eco_Pools'][iT-1,iApp,iEP['RootCoarse']]+Gnet_RC[iApp]
        #vo['C_Eco_Pools'][iT,iApp,iEP['RootFine']]=vo['C_Eco_Pools'][iT-1,iApp,iEP['RootFine']]+Gnet_RF[iApp]

        # Populate net growth of root biomass
        vo['C_G_Net'][iT,iApp,iEP['RootCoarse']]=rrRC*vo['C_G_Net'][iT,iApp,iEP['RootCoarse']]
        vo['C_G_Net'][iT,iApp,iEP['RootFine']]=rrRF*vo['C_G_Net'][iT,iApp,iEP['RootFine']]
        #vo['C_G_Net'][iT,iApp,iEP['RootCoarse']]=Gnet_RC[iApp]
        #vo['C_G_Net'][iT,iApp,iEP['RootFine']]=Gnet_RF[iApp]

        #----------------------------------------------------------------------
        # Stop stimulation counter when it passes the response duration
        #----------------------------------------------------------------------

        iStop=np.where( meta['Modules']['Nutrient Management']['ResponseCounter']>bNA['ResponseDuration'] )[0]

        if iStop.size>0:

            meta['Modules']['Nutrient Management']['ResponseCounter'][iStop]=0

            #uGC=np.unique(meta['Modules']['GYM']['ID GC'][iStop])
            #for iGC in range(uGC.size):
            #    indGC=np.where(meta['Modules']['GYM']['ID GC'][iStop]==uGC[iGC])[0]
            #    vi['GC']['Active'][:,iStop[indGC],:]=vi['GC'][ uGC[iGC] ][:,iStop[indGC],:]

    elif (comp=='Mortality') & (meta[pNam]['Project']['Nutrient Application Module']=='cbrunner'):

        #----------------------------------------------------------------------
        # Adjust mortality
        #----------------------------------------------------------------------

        vo['C_M_Reg'][iT,iApp,0:7]=bNA['rPrime_TreeMortality']*vo['C_M_Reg'][iT,iApp,0:7]

    elif (comp=='Litterfall') & (meta[pNam]['Project']['Nutrient Application Module']=='cbrunner'):

        #----------------------------------------------------------------------
        # Adjust biomass turnover
        #----------------------------------------------------------------------

        iEP=meta['Core']['iEP']

        vo['C_LF'][iT,iApp,iEP['Foliage']]=bNA['rPrime_Litterfall']*vo['C_LF'][iT,iApp,iEP['Foliage']]
        vo['C_LF'][iT,iApp,iEP['Branch']]=bNA['rPrime_Litterfall']*vo['C_LF'][iT,iApp,iEP['Branch']]
        vo['C_LF'][iT,iApp,iEP['Bark']]=bNA['rPrime_Litterfall']*vo['C_LF'][iT,iApp,iEP['Bark']]

        vo['C_LF'][iT,iApp,iEP['RootCoarse']]=bNA['rPrime_TurnoverRootCoarse']*vo['C_LF'][iT,iApp,iEP['RootCoarse']]
        vo['C_LF'][iT,iApp,iEP['RootFine']]=bNA['rPrime_TurnoverRootFine']*vo['C_LF'][iT,iApp,iEP['RootFine']]

    elif comp=='Emissions':

        #----------------------------------------------------------------------
        # Emissions from manufacture of ammonia and urea
        #----------------------------------------------------------------------

        # Urea dose (kgUrea/ha)
        DoseUrea=bNA['DoseUrea_Standard']

        # Nitrogen dose (kgN/ha)
        DoseN=bNA['Ratio_N_to_Urea']*DoseUrea

        # Emissions from production of ammonia (tCO2e/ha)
        E_ProdNH3_per_t=bNA['EmissionFromAmmoniaProduction_NRCAN']*bNA['TonneNH3PerTonneUrea']
        E_ProdNH3=E_ProdNH3_per_t*(DoseUrea/1000)

        # Emissions from production of Urea from ammonia (tCO2e/ha)
        MMBtu_per_app=(DoseUrea/1000)*bNA['UreaEnergyConsumption']
        therm_per_app=MMBtu_per_app/bNA['MMBtu_per_therm']
        E_ProdUrea=bNA['EmissionFromUreaProduction_per_therm']*therm_per_app

        vo['E_CO2e_ESC_OperForBurnGas'][iT,iApp]=vo['E_CO2e_ESC_OperForBurnGas'][iT,iApp] + \
            (E_ProdNH3+E_ProdUrea)

        #----------------------------------------------------------------------
        # Emissions from transportation (tCO2e/ha)
        #----------------------------------------------------------------------

        vo['E_CO2e_ET_OperForBurnOil'][iT,iApp]=vo['E_CO2e_ET_OperForBurnOil'][iT,iApp] + \
            (bNA['EmissionFromRailBargeTruck_Workbook']+bNA['EmissionFromHelicopter_SP10'])

        #----------------------------------------------------------------------
        # Denitrification: N2O emissions following application, Tier 1
        # approach, IPCC 2006, 11.4.1 (tCO2e/ha)
        #----------------------------------------------------------------------

        vo['E_CO2e_LULUCF_Denit'][iT,iApp]=vo['E_CO2e_LULUCF_Denit'][iT,iApp] + \
            (bNA['EmissionFactor_N2O_Jassaletal2008']*(DoseN/1000)*bNA['Ratio_N2OAsN_to_N2O']*meta['Param']['BEV']['Biophysical']['GWP_N2O_AR4'])

        #----------------------------------------------------------------------
        # Volatilization
        #----------------------------------------------------------------------

        # CO2 emissions following application, Tier 1 approach, IPCC 2006, 11.4.1 (tCO2e/ha)
        # 0.2*(430/1000)*(1/0.27) = 0.32 tCO2e/ha
        #E_Vol=bNA['Ratio_C_to_Urea']*(DoseUrea/1000)*(1/meta['Param']['BEV']['Biophysical']['Ratio_C_to_CO2'])
        # Assume sequestration and volatilization of CO2 cancel out

        E_vol=0.3

        vo['E_CO2e_IPPU_OperForBurningGas'][iT,iApp]=vo['E_CO2e_IPPU_OperForBurningGas'][iT,iApp] - E_vol

        vo['E_CO2e_LULUCF_Other'][iT,iApp]=vo['E_CO2e_LULUCF_Other'][iT,iApp] + E_vol

        #----------------------------------------------------------------------
        # Exterior area (volatilization/deposition effects)
        #----------------------------------------------------------------------

        flg=1

        if 'Nutrient Application Footprint Status' in meta[pNam]['Scenario'][meta[pNam]['iScn']]:
            if meta[pNam]['Scenario'][meta[pNam]['iScn']]['Nutrient Application Footprint Status']!='On':
                flg=0

        if (meta[pNam]['Project']['External Footprint Effect Status']=='On') & (flg==1):

            # Dose (kgN/ha)
            DoseN=bNA['DoseUrea_Standard']*bNA['Ratio_N_to_Urea']

            # Emissions (NH3-N ha-1)
            EmissionNH3_as_N=bNA['EA Fraction volatilized']*DoseN

            # Canopy uptake (kgN)
            CanopyUptakeN=EmissionNH3_as_N*bNA['EA Forest deposition fraction']*bNA['EA Leaf uptake fraction']

            # Root uptake (kgN)
            RootUptakeN=EmissionNH3_as_N*bNA['EA Forest deposition fraction']*bNA['EA Throughfall fraction']*bNA['EA RootUptakeFraction']

            # GHG benefit from canopy uptake (MgC/ha)
            GHG_Benefit_CanupyUptake_Coast=bNA['EA_NUEu_NGTT_Coast']*CanopyUptakeN/1000
            GHG_Benefit_CanupyUptake_Interior=CanopyUptakeN/1000*bNA['EA_NUEu_NGTT_Interior']

            # GHG benefit from root uptake (MgC/ha)
            GHG_Benefit_RootUptake_Coast=RootUptakeN/1000*bNA['EA_NUEa_NGTTD_Coast']
            GHG_Benefit_RootUptake_Interior=RootUptakeN/1000*bNA['EA_NUEa_NGTTD_Interior']

            # Total GHG benefit (MgC/ha)
            GHG_Benefit_Tot_Coast=GHG_Benefit_CanupyUptake_Coast+GHG_Benefit_RootUptake_Coast
            GHG_Benefit_Tot_Interior=GHG_Benefit_CanupyUptake_Interior+GHG_Benefit_RootUptake_Interior

            EA_GHG_Benefit=GHG_Benefit_Tot_Coast*bNA['EA Fraction of footprint coast']+GHG_Benefit_Tot_Interior*bNA['EA Fraction of footprint interior']

            # Annaul GHG benefit over response duration (tCO2e/ha/yr)
            #EA_GHG_Benefit=EA_GHG_Benefit/bNA['ResponseDuration']

            # Convert to CO2e (tCO2e/ha/yr)
            EA_GHG_Benefit=meta['Param']['BEV']['Biophysical']['Ratio_CO2_to_C']*EA_GHG_Benefit

            # Subtract from ecosystem LULUCF emissions
            # *** it will crash in the last time step ***
            try:
                vo['E_CO2e_LULUCF_Other'][iT+1,iApp]=vo['E_CO2e_LULUCF_Other'][iT+1,iApp]-EA_GHG_Benefit
            except:
                pass

    elif (comp=='HeterotrophicRespiration') & (meta[pNam]['Project']['Nutrient Application Module']=='cbrunner'):

        #----------------------------------------------------------------------
        # Adjust rate of heterotrophic respiration
        #----------------------------------------------------------------------

        rr=bNA['r_Decomp']
        meta[pNam]['Project']['R_LitterVF'][0,iApp]=rr*meta[pNam]['Project']['R_LitterVF'][0,iApp]
        meta[pNam]['Project']['R_LitterF'][0,iApp]=rr*meta[pNam]['Project']['R_LitterF'][0,iApp]
        meta[pNam]['Project']['R_LitterM'][0,iApp]=rr*meta[pNam]['Project']['R_LitterM'][0,iApp]
        meta[pNam]['Project']['R_LitterS'][0,iApp]=rr*meta[pNam]['Project']['R_LitterS'][0,iApp]
        meta[pNam]['Project']['R_SoilVF'][0,iApp]=rr*meta[pNam]['Project']['R_SoilVF'][0,iApp]
        meta[pNam]['Project']['R_SoilF'][0,iApp]=rr*meta[pNam]['Project']['R_SoilF'][0,iApp]
        meta[pNam]['Project']['R_SoilS'][0,iApp]=rr*meta[pNam]['Project']['R_SoilS'][0,iApp]

    return vi,vo,meta

#%% Simulate probability of harvesting on the fly

def ScheduleNutrientApplication(meta,pNam,vi,vo,iT,iScn,iEns,iBat):

    # Create a random number
    rn=np.random.random(meta[pNam]['Project']['Batch Size'][iBat])

    # Regional probabilities
    Po_Sat_Coast=0.00775
    Po_Sat_Interior=0.00375

    flg=0
    if flg==1:
        vi={}
        vi['tv']=tv

        Po_Coast=np.maximum(Po_Sat_Coast,Po_Sat_Coast+0.002*np.maximum(1,vi['tv']-2021) )
        plt.close('all')
        plt.plot(vi['tv'],Po_Coast,'r-',lw=1.5 )

        Po_Interior=np.maximum(Po_Sat_Interior,Po_Sat_Interior+0.0012*np.maximum(1,vi['tv']-2045) )
        plt.close('all')
        plt.plot(vi['tv'],Po_Interior,'r-',lw=1.5 )

    # Contstant
    Po_Coast=Po_Sat_Coast*np.ones(vi['tv'].size)
    Po_Interior=Po_Sat_Interior*np.ones(vi['tv'].size)

    # Time-dependent models to compensate for aging forests (paper)
    #Po_Coast=np.maximum(Po_Sat_Coast,Po_Sat_Coast+0.0005*np.maximum(1,vi['tv']-2021) )
    #Po_Interior=np.maximum(Po_Sat_Interior,Po_Sat_Interior+0.0005*np.maximum(1,vi['tv']-2021) )

    # Find eligible coastal stands to fertilize
    indS_Coast=np.where( (meta['Modules']['Nutrient Management']['ResponseCounter']==0) & \
            (vo['A'][iT,:]>=9) & \
            (vo['A'][iT,:]<=71) & \
            (rn<Po_Coast[iT]) & \
            (vo['V_MerchLive'][iT,:]>10) & \
            (np.isin(vi['Inv']['ID_BECZ'][0,:],meta['Modules']['Nutrient Management']['BGC Zone Exclusion ID'])==False) & \
            (np.isin(vi['Inv']['ID_BECZ'][0,:],meta['Modules']['Nutrient Management']['Coastal Zones ID'])==True) )[0]

    indS_Interior=np.where( (meta['Modules']['Nutrient Management']['ResponseCounter']==0) & \
            (vo['A'][iT,:]>=9) & \
            (vo['A'][iT,:]<=71) & \
            (rn<Po_Interior[iT]) & \
            (vo['V_MerchLive'][iT,:]>10) & \
            (np.isin(vi['Inv']['ID_BECZ'][0,:],meta['Modules']['Nutrient Management']['BGC Zone Exclusion ID'])==False) & \
            (np.isin(vi['Inv']['ID_BECZ'][0,:],meta['Modules']['Nutrient Management']['Coastal Zones ID'])==False) )[0]

    indS=np.append(indS_Coast,indS_Interior)

    if indS.size>0:
        for i in range(indS.size):
            iAvailable=np.where(vi['EC']['ID Event Type'][iT,indS[i],:]==0)[0]
            if iAvailable.size>0:
                iE=iAvailable[0]
                vi['EC']['ID Event Type'][iT,indS[i],iE]=meta['LUT']['Event']['Fertilization Aerial']
                vi['EC']['MortalityFactor'][iT,indS[i],iE]=np.array(0,dtype='int16')
                #vi['EC']['ID Growth Curve'][iT,indS[i],iE]=np.max(vi['EC']['ID Growth Curve'][0:iT,indS[i],:])

    return vi