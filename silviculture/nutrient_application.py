
import numpy as np
from fcgadgets.cbrunner import cbrun_utilities as cbu

#%% Nutrient application effects

def update_nutrient_status(vi,vo,iT,meta,comp):
    
    # Exctract parameters
    bNA=meta['Param']['Nutrients']
    
    if comp=='UpdateCounter':
        
        iApplication=meta['Nutrient Management']['iApplication']
        
        meta['Nutrient Management']['ResponseCounter'][iApplication]=meta['Nutrient Management']['ResponseCounter'][iApplication]+1
    
    elif (comp=='AbovegroundNetGrowth') & (meta['Project']['Nutrient Application Module']=='cbrunner'):

        # Need to revise to be based on N
        # r_NonUniform_Spatial_Distribution
        
        # Index to stands where application occurs
        iApplication=meta['Nutrient Management']['iApplication']

        #----------------------------------------------------------------------
        # Start response counter
        #----------------------------------------------------------------------
        
        meta['Nutrient Management']['ResponseCounter'][iApplication]=meta['Nutrient Management']['ResponseCounter'][iApplication]+1
        
        #----------------------------------------------------------------------
        # Adjust net growth
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
        A_gc=np.arange(0,meta['GC']['BatchTIPSY Maximum Age']+1,1)
        
        # Extract active growth curves, apply scale factor
        GCA_SP=vi['GC']['Active'][:,iApplication,:].copy().astype(float)*meta['GC']['Scale Factor']
        
        for iStand in range(iApplication.size):
            
            AgeAtApplication=int(vo['A'][iT,iApplication[iStand]])
            
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
                print(AgeAtApplication)
            
        # Re-applly scalefactor
        GCA_SP=GCA_SP/meta['GC']['Scale Factor']
            
        # Convert to 16-bit integers
        GCA_SP=GCA_SP.astype(np.int16)
            
        # Repopulate in input variable dictionary
        vi['GC']['Active'][:,iApplication,:]=GCA_SP        
    
    elif (comp=='BelowgroundNetGrowth') & (meta['Project']['Nutrient Application Module']=='cbrunner'):
        
        #----------------------------------------------------------------------
        # Adjust root net growth
        #----------------------------------------------------------------------
        
        # Index to stands where application occurs
        iApplication=meta['Nutrient Management']['iApplication']
        
        iEP=meta['Core']['iEP']
        
        # Responses
        rrRC=1+(bNA['r_Stemwood']-1)*bNA['Ratio_RootC_to_Stemwood']
        rrRF=1+(bNA['r_Stemwood']-1)*bNA['Ratio_RootF_to_Stemwood']
        
        # Calculate net growth of roots from change in pools
        #Gnet_RC=vo['C_Eco_Pools'][iT,:,iEP['RootCoarse']]-vo['C_Eco_Pools'][iT-1,:,iEP['RootCoarse']]
        #Gnet_RF=vo['C_Eco_Pools'][iT,:,iEP['RootFine']]-vo['C_Eco_Pools'][iT-1,:,iEP['RootFine']] 
        
        # Remove the false stimulus that arises from estimating roots from AG
        # biomass  
        vo['C_G_Net'][iT,iApplication,iEP['RootCoarse']]=0.84*vo['C_G_Net'][iT,iApplication,iEP['RootCoarse']]
        vo['C_G_Net'][iT,iApplication,iEP['RootFine']]=0.84*vo['C_G_Net'][iT,iApplication,iEP['RootFine']]
        #Gnet_RC[iApplication]=0.65*(2-bNA['r_Stemwood'])*Gnet_RC[iApplication]
        #Gnet_RF[iApplication]=0.65*(2-bNA['r_Stemwood'])*Gnet_RF[iApplication] 
        
        # Stimulate root growth from N application
        #Gnet_RC[iApplication]=rrRC*Gnet_RC[iApplication]
        #Gnet_RF[iApplication]=rrRF*Gnet_RF[iApplication]        
        
        # Recalculate current-year root biomass where stimulated by N application
        #vo['C_Eco_Pools'][iT,iApplication,iEP['RootCoarse']]=vo['C_Eco_Pools'][iT-1,iApplication,iEP['RootCoarse']]+Gnet_RC[iApplication]
        #vo['C_Eco_Pools'][iT,iApplication,iEP['RootFine']]=vo['C_Eco_Pools'][iT-1,iApplication,iEP['RootFine']]+Gnet_RF[iApplication]
        
        # Populate net growth of root biomass
        vo['C_G_Net'][iT,iApplication,iEP['RootCoarse']]=rrRC*vo['C_G_Net'][iT,iApplication,iEP['RootCoarse']]
        vo['C_G_Net'][iT,iApplication,iEP['RootFine']]=rrRF*vo['C_G_Net'][iT,iApplication,iEP['RootFine']]
        #vo['C_G_Net'][iT,iApplication,iEP['RootCoarse']]=Gnet_RC[iApplication]
        #vo['C_G_Net'][iT,iApplication,iEP['RootFine']]=Gnet_RF[iApplication]    
        
        #----------------------------------------------------------------------
        # Stop stimulation counter when it passes the response duration
        #----------------------------------------------------------------------
    
        iStop=np.where( meta['Nutrient Management']['ResponseCounter']>bNA['ResponseDuration'] )[0]
        
        if iStop.size>0:
            
            meta['Nutrient Management']['ResponseCounter'][iStop]=0
            
            #uGC=np.unique(meta['GC']['ID GC'][iStop])
            #for iGC in range(uGC.size):
            #    indGC=np.where(meta['GC']['ID GC'][iStop]==uGC[iGC])[0]
            #    vi['GC']['Active'][:,iStop[indGC],:]=vi['GC'][ uGC[iGC] ][:,iStop[indGC],:]
        
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

#%% Simulate probability of harvesting on the fly

def ScheduleNutrientApplication_OnTheFly(meta,vi,vo,iT,iScn,iEns):
    
   rn=np.random.random(vo['A'][iT,:].size) 
    
   indS=np.where( (meta['Nutrient Management']['ResponseCounter']==0) & \
                 (vo['A'][iT,:]>=10) & (vo['A'][iT,:]<=71) & \
                 (rn<meta['Scenario'][iScn]['Nutrient Application Prob']) & \
                 (vo['V_StemMerch'][iT,:]>1) )[0]

   if indS.size>0:
       for i in range(indS.size):
           iAvailable=np.where(vi['EC']['ID_Type'][iT,indS[i],:]==0)[0]        
           if iAvailable.size>0:
               iE=iAvailable[0]
               vi['EC']['ID_Type'][iT,indS[i],iE]=meta['LUT']['Dist']['Fertilization Aerial']
               vi['EC']['MortalityFactor'][iT,indS[i],iE]=np.array(0,dtype='int16')
               vi['EC']['ID_GrowthCurve'][iT,indS[i],iE]=2
    
   return vi