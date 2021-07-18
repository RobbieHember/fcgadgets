
import numpy as np

#%% Nutrient application effects

def update_nutrient_status(vi,vo,iT,meta,comp):
    
    # Exctract parameters
    bNA=meta['Param']['Nutrients']
    
    if comp=='UpdateCounter':
        
        meta['Nutrient Management']['ResponseCounter'][meta['Nutrient Management']['iApplication']]=meta['Nutrient Management']['ResponseCounter'][meta['Nutrient Management']['iApplication']]+1
    
    elif (comp=='AbovegroundNetGrowth') & (meta['Project']['Nutrient Application Module']=='cbrunner'):

        # Need to revise to be based on N
        # r_NonUniform_Spatial_Distribution
        
        #----------------------------------------------------------------------
        # Adjust net growth
        #----------------------------------------------------------------------
        
        # Trigger flag indicating stand is being affectd by N application
        # This is used to specify duration of the stimulus 
        meta['Nutrient Management']['ResponseCounter'][meta['Nutrient Management']['iApplication']]=meta['Nutrient Management']['ResponseCounter'][meta['Nutrient Management']['iApplication']]+1
        
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
        A_gc=np.arange(0,meta['GC']['BatchTIPSY Maximum Age']+1,1)
        
        # Extract active growth curves, apply scale factor
        GCA_SP=vi['GC']['Active'][:,meta['Nutrient Management']['iApplication'],:].copy().astype(float)*meta['GC']['Scale Factor']

        for iStand in range(meta['Nutrient Management']['iApplication'].size):
            
            A_app=int(vo['A'][iT,meta['Nutrient Management']['iApplication'][iStand]])
            
            # Index to the response period that will be alterred
            iResponse=np.where( (A_gc>=A_app) & (A_gc<A_app+bNA['ResponseDuration']) )[0]
            
            # This will crash if an application occurs within 10 years of the maximum
            # TIPSY age curve (eg 200 years) -> adjust response period so that it
            # does not crash
            if iResponse.size==int(bNA['ResponseDuration']):
                
                for iRR in range(rr.size):
                    GCA_SP[iResponse,iStand,iRR]=rr[iRR]*GCA_SP[iResponse,iStand,iRR]
        
            else:
                em='Error: Nutrient application not implemented - stand age exceeds max age of growth curves.'
                print(em)
                print(A_app)
                
        # Re-applly scalefactor
        GCA_SP=GCA_SP/meta['GC']['Scale Factor']
            
        # Convert to 16-bit integers
        GCA_SP=GCA_SP.astype(np.int16)
            
        # Repopulate in input variable dictionary
        vi['GC']['Active'][:,meta['Nutrient Management']['iApplication'],:]=GCA_SP
    
    elif (comp=='BelowgroundNetGrowth') & (meta['Project']['Nutrient Application Module']=='cbrunner'):
        
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
        vo['C_G_Net'][iT,meta['Nutrient Management']['iApplication'],iEP['RootCoarse']]=0.84*vo['C_G_Net'][iT,meta['Nutrient Management']['iApplication'],iEP['RootCoarse']]
        vo['C_G_Net'][iT,meta['Nutrient Management']['iApplication'],iEP['RootFine']]=0.84*vo['C_G_Net'][iT,meta['Nutrient Management']['iApplication'],iEP['RootFine']]
        #Gnet_RC[meta['Nutrient Management']['iApplication']]=0.65*(2-bNA['r_Stemwood'])*Gnet_RC[meta['Nutrient Management']['iApplication']]
        #Gnet_RF[meta['Nutrient Management']['iApplication']]=0.65*(2-bNA['r_Stemwood'])*Gnet_RF[meta['Nutrient Management']['iApplication']] 
        
        # Stimulate root growth from N application
        #Gnet_RC[meta['Nutrient Management']['iApplication']]=rrRC*Gnet_RC[meta['Nutrient Management']['iApplication']]
        #Gnet_RF[meta['Nutrient Management']['iApplication']]=rrRF*Gnet_RF[meta['Nutrient Management']['iApplication']]        
        
        # Recalculate current-year root biomass where stimulated by N application
        #vo['C_Eco_Pools'][iT,meta['Nutrient Management']['iApplication'],iEP['RootCoarse']]=vo['C_Eco_Pools'][iT-1,meta['Nutrient Management']['iApplication'],iEP['RootCoarse']]+Gnet_RC[meta['Nutrient Management']['iApplication']]
        #vo['C_Eco_Pools'][iT,meta['Nutrient Management']['iApplication'],iEP['RootFine']]=vo['C_Eco_Pools'][iT-1,meta['Nutrient Management']['iApplication'],iEP['RootFine']]+Gnet_RF[meta['Nutrient Management']['iApplication']]
        
        # Populate net growth of root biomass
        vo['C_G_Net'][iT,meta['Nutrient Management']['iApplication'],iEP['RootCoarse']]=rrRC*vo['C_G_Net'][iT,meta['Nutrient Management']['iApplication'],iEP['RootCoarse']]
        vo['C_G_Net'][iT,meta['Nutrient Management']['iApplication'],iEP['RootFine']]=rrRF*vo['C_G_Net'][iT,meta['Nutrient Management']['iApplication'],iEP['RootFine']]
        #vo['C_G_Net'][iT,meta['Nutrient Management']['iApplication'],iEP['RootCoarse']]=Gnet_RC[meta['Nutrient Management']['iApplication']]
        #vo['C_G_Net'][iT,meta['Nutrient Management']['iApplication'],iEP['RootFine']]=Gnet_RF[meta['Nutrient Management']['iApplication']]    
        
        #----------------------------------------------------------------------
        # Stop stimulation counter when it passes the response duration
        #----------------------------------------------------------------------
    
        iStop=np.where( meta['Nutrient Management']['ResponseCounter']>bNA['ResponseDuration'] )[0]
        if iStop.size>0:
            meta['Nutrient Management']['ResponseCounter'][iStop]=0
        
    elif (comp=='Mortality') & (meta['Project']['Nutrient Application Module']=='cbrunner'):
        
        #----------------------------------------------------------------------
        # Adjust mortality
        #----------------------------------------------------------------------
        
        vo['C_M_Reg'][iT,meta['Nutrient Management']['iApplication'],0:7]=bNA['rPrime_TreeMortality']*vo['C_M_Reg'][iT,meta['Nutrient Management']['iApplication'],0:7] 
     
    elif (comp=='Litterfall') & (meta['Project']['Nutrient Application Module']=='cbrunner'):
        
        #----------------------------------------------------------------------
        # Adjust biomass turnover
        #----------------------------------------------------------------------
        
        iEP=meta['Core']['iEP']
        
        vo['C_LF'][iT,meta['Nutrient Management']['iApplication'],iEP['Foliage']]=bNA['rPrime_Litterfall']*vo['C_LF'][iT,meta['Nutrient Management']['iApplication'],iEP['Foliage']]
        vo['C_LF'][iT,meta['Nutrient Management']['iApplication'],iEP['Branch']]=bNA['rPrime_Litterfall']*vo['C_LF'][iT,meta['Nutrient Management']['iApplication'],iEP['Branch']]
        vo['C_LF'][iT,meta['Nutrient Management']['iApplication'],iEP['Bark']]=bNA['rPrime_Litterfall']*vo['C_LF'][iT,meta['Nutrient Management']['iApplication'],iEP['Bark']]
        
        vo['C_LF'][iT,meta['Nutrient Management']['iApplication'],iEP['RootCoarse']]=bNA['rPrime_TurnoverRootCoarse']*vo['C_LF'][iT,meta['Nutrient Management']['iApplication'],iEP['RootCoarse']]
        vo['C_LF'][iT,meta['Nutrient Management']['iApplication'],iEP['RootFine']]=bNA['rPrime_TurnoverRootFine']*vo['C_LF'][iT,meta['Nutrient Management']['iApplication'],iEP['RootFine']]
    
    elif comp=='Emissions':
        
        #----------------------------------------------------------------------
        # Adjust emissions
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
            
        # Emissions from operations (tCO2e/ha)
        E_Ops=bNA['EmissionFromRailBargeTruck_Workbook']+ \
            bNA['EmissionFromHelicopter_SP10']

        # Volatilization: CO2 emissions following application, Tier 1 
        # approach, IPCC 2006, 11.4.1 (tCO2e/ha)
        # 0.2*(430/1000)*(1/0.27) = 0.32 tCO2e/ha
        #E_Vol=bNA['Ratio_C_to_Urea']*(DoseUrea/1000)*(1/meta['Param']['Biophysical']['bRatio_C_to_CO2'])
        # Assume sequestration and volatilization of CO2 cancel out 
        E_Vol=0

        # Denitrification: N2O emissions following application, Tier 1 
        # approach, IPCC 2006, 11.4.1 (tCO2e/ha)
        E_Denit=bNA['EmissionFactor_N2O_Jassaletal2008']*(DoseN/1000)* \
            bNA['Ratio_N2OAsN_to_N2O']*meta['Param']['Biophysical']['GWP_N2O_AR4']
                
        # Total emissions (tCO2e/ha)
        E_Tot=E_ProdNH3+E_ProdUrea+E_Ops+E_Vol+E_Denit
                
        # Total emissions of carbon (MgC/ha) emitted as CO2e, converted
        # to carbon to be consistent with the rest of the variables in 
        # the vo object.
        vo['C_E_Operations'][iT,meta['Nutrient Management']['iApplication']]=meta['Param']['Biophysical']['Ratio_C_to_CO2']*E_Tot
    
    elif (comp=='HeterotrophicRespiration') & (meta['Project']['Nutrient Application Module']=='cbrunner'):
        
        #----------------------------------------------------------------------
        # Adjust rate of heterotrophic respiration
        #----------------------------------------------------------------------
        
        rr=bNA['r_Decomp']
        meta['R_LitterVF'][0,meta['Nutrient Management']['iApplication']]=rr*meta['R_LitterVF'][0,meta['Nutrient Management']['iApplication']]
        meta['R_LitterF'][0,meta['Nutrient Management']['iApplication']]=rr*meta['R_LitterF'][0,meta['Nutrient Management']['iApplication']]
        meta['R_LitterM'][0,meta['Nutrient Management']['iApplication']]=rr*meta['R_LitterM'][0,meta['Nutrient Management']['iApplication']]
        meta['R_LitterS'][0,meta['Nutrient Management']['iApplication']]=rr*meta['R_LitterS'][0,meta['Nutrient Management']['iApplication']]
        meta['R_SoilVF'][0,meta['Nutrient Management']['iApplication']]=rr*meta['R_SoilVF'][0,meta['Nutrient Management']['iApplication']]
        meta['R_SoilF'][0,meta['Nutrient Management']['iApplication']]=rr*meta['R_SoilF'][0,meta['Nutrient Management']['iApplication']]
        meta['R_SoilS'][0,meta['Nutrient Management']['iApplication']]=rr*meta['R_SoilS'][0,meta['Nutrient Management']['iApplication']]
    
    return vi,vo,meta
