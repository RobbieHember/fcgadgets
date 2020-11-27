
import numpy as np

#%% Nutrient application effects

def update_nutrient_status(vi,vo,iT,meta,psl,comp):
    
    if comp=='UpdateCounter':
        
        meta['NM']['ResponseCounter'][meta['NM']['iApplication']]=meta['NM']['ResponseCounter'][meta['NM']['iApplication']]+1
    
    elif (comp=='AbovegroundNetGrowth') & (meta['Fertilization Source']=='CBRunner'):
        
        # Need to revise to be based on N
        # r_NonUniform_Spatial_Distribution
        
        #----------------------------------------------------------------------
        # Adjust net growth
        #----------------------------------------------------------------------
        
        # Trigger flag indicating stand is being affectd by N application
        # This is used to specify duration of the stimulus 
        meta['NM']['ResponseCounter'][meta['NM']['iApplication']]=meta['NM']['ResponseCounter'][meta['NM']['iApplication']]+1
        
        # Duration
        Duration=int(psl['bNA_ResponseDuration'])
        
        # Response ratio of stemwood (before wood density effect)
        rrS=psl['bNA_r_Stemwood']
        
        # Wood density effect
        rr_wd=psl['bNA_r_WoodDensity']
        f_wd=rr_wd-1
        f_wd=np.array([f_wd,f_wd,0,0,0])
        
        # Response relative to stemwood
        # Can't do roots here because they are not tracked in growth curve dict
        rXS=np.array([1,
                1,
                psl['bNA_Ratio_Foliage_to_Stemwood'],
                psl['bNA_Ratio_Branch_to_Stemwood'],
                psl['bNA_Ratio_Bark_to_Stemwood']])
        
        # Biomass response ratios (after wood density effect accounted for) 
        rr=1+(rrS-1)*rXS+f_wd
            
        # Append biomass response with volume response
        rr=np.append(rr,rrS)
        
        A=np.arange(1,meta['GC']['BatchTIPSY Maximum Age'],1)
        
        iResponse=np.where( (A>=vo['A'][iT,meta['NM']['iApplication']]) & (A<vo['A'][iT,meta['NM']['iApplication']]+Duration) )[0]
        
        if iResponse.size==Duration:        
            
            GCA_SP=vi['GC']['Active'][:,meta['NM']['iApplication'],:].copy().astype(float)*meta['GC']['Scale Factor']
            
            print(GCA_SP.shape)
            print(np.tile(rr,(Duration,1)).shape)
            print(rr)
            
            GCA_SP[iResponse,0,:]=np.tile(rr,(Duration,1))*GCA_SP[iResponse,0,:]
            #GCA_SP=GCA_SP*meta['GC']['Scale Factor']
            #GCA_SP=GCA_SP.astype(np.int16)
            vi['GC']['Active'][:,meta['NM']['iApplication'],:]=GCA_SP
    
    elif (comp=='BelowgroundNetGrowth') & (meta['Fertilization Source']=='CBRunner'):
        
        #----------------------------------------------------------------------
        # Adjust root net growth
        #----------------------------------------------------------------------
        
        iEP=meta['iEP']
        
        # Responses
        rrRC=1+(psl['bNA_r_Stemwood']-1)*psl['bNA_Ratio_RootC_to_Stemwood']
        rrRF=1+(psl['bNA_r_Stemwood']-1)*psl['bNA_Ratio_RootF_to_Stemwood']
        
        # Calculate net growth of roots from change in pools
        #Gnet_RC=vo['C_Eco_Pools'][iT,:,iEP['RootCoarse']]-vo['C_Eco_Pools'][iT-1,:,iEP['RootCoarse']]
        #Gnet_RF=vo['C_Eco_Pools'][iT,:,iEP['RootFine']]-vo['C_Eco_Pools'][iT-1,:,iEP['RootFine']] 
        
        # Remove the false stimulus that arises from estimating roots from AG
        # biomass  
        vo['C_G_Net'][iT,meta['NM']['iApplication'],iEP['RootCoarse']]=0.84*vo['C_G_Net'][iT,meta['NM']['iApplication'],iEP['RootCoarse']]
        vo['C_G_Net'][iT,meta['NM']['iApplication'],iEP['RootFine']]=0.84*vo['C_G_Net'][iT,meta['NM']['iApplication'],iEP['RootFine']]
        #Gnet_RC[meta['NM']['iApplication']]=0.65*(2-psl['bNA_r_Stemwood'])*Gnet_RC[meta['NM']['iApplication']]
        #Gnet_RF[meta['NM']['iApplication']]=0.65*(2-psl['bNA_r_Stemwood'])*Gnet_RF[meta['NM']['iApplication']] 
        
        # Stimulate root growth from N application
        #Gnet_RC[meta['NM']['iApplication']]=rrRC*Gnet_RC[meta['NM']['iApplication']]
        #Gnet_RF[meta['NM']['iApplication']]=rrRF*Gnet_RF[meta['NM']['iApplication']]        
        
        # Recalculate current-year root biomass where stimulated by N application
        #vo['C_Eco_Pools'][iT,meta['NM']['iApplication'],iEP['RootCoarse']]=vo['C_Eco_Pools'][iT-1,meta['NM']['iApplication'],iEP['RootCoarse']]+Gnet_RC[meta['NM']['iApplication']]
        #vo['C_Eco_Pools'][iT,meta['NM']['iApplication'],iEP['RootFine']]=vo['C_Eco_Pools'][iT-1,meta['NM']['iApplication'],iEP['RootFine']]+Gnet_RF[meta['NM']['iApplication']]
        
        # Populate net growth of root biomass
        vo['C_G_Net'][iT,meta['NM']['iApplication'],iEP['RootCoarse']]=rrRC*vo['C_G_Net'][iT,meta['NM']['iApplication'],iEP['RootCoarse']]
        vo['C_G_Net'][iT,meta['NM']['iApplication'],iEP['RootFine']]=rrRF*vo['C_G_Net'][iT,meta['NM']['iApplication'],iEP['RootFine']]
        #vo['C_G_Net'][iT,meta['NM']['iApplication'],iEP['RootCoarse']]=Gnet_RC[meta['NM']['iApplication']]
        #vo['C_G_Net'][iT,meta['NM']['iApplication'],iEP['RootFine']]=Gnet_RF[meta['NM']['iApplication']]    
        
        #----------------------------------------------------------------------
        # Stop stimulation counter when it passes the response duration
        #----------------------------------------------------------------------
    
        iStop=np.where( meta['NM']['ResponseCounter']>10 )[0]
        if iStop.size>0:
            meta['NM']['ResponseCounter'][iStop]=0
        
    elif (comp=='Mortality') & (meta['Fertilization Source']=='CBRunner'):
        
        #----------------------------------------------------------------------
        # Adjust mortality
        #----------------------------------------------------------------------
        
        vo['C_M_Reg'][iT,meta['NM']['iApplication'],0:7]=psl['bNA_rPrime_TreeMortality']*vo['C_M_Reg'][iT,meta['NM']['iApplication'],0:7] 
     
    elif (comp=='Litterfall') & (meta['Fertilization Source']=='CBRunner'):
        
        #----------------------------------------------------------------------
        # Adjust biomass turnover
        #----------------------------------------------------------------------
        
        iEP=meta['iEP']
        
        vo['C_LF'][iT,meta['NM']['iApplication'],iEP['Foliage']]=psl['bNA_rPrime_Litterfall']*vo['C_LF'][iT,meta['NM']['iApplication'],iEP['Foliage']]
        vo['C_LF'][iT,meta['NM']['iApplication'],iEP['Branch']]=psl['bNA_rPrime_Litterfall']*vo['C_LF'][iT,meta['NM']['iApplication'],iEP['Branch']]
        vo['C_LF'][iT,meta['NM']['iApplication'],iEP['Bark']]=psl['bNA_rPrime_Litterfall']*vo['C_LF'][iT,meta['NM']['iApplication'],iEP['Bark']]
        
        vo['C_LF'][iT,meta['NM']['iApplication'],iEP['RootCoarse']]=psl['bNA_rPrime_TurnoverRootCoarse']*vo['C_LF'][iT,meta['NM']['iApplication'],iEP['RootCoarse']]
        vo['C_LF'][iT,meta['NM']['iApplication'],iEP['RootFine']]=psl['bNA_rPrime_TurnoverRootFine']*vo['C_LF'][iT,meta['NM']['iApplication'],iEP['RootFine']]
    
    elif comp=='Emissions':
        
        #----------------------------------------------------------------------
        # Adjust emissions
        #----------------------------------------------------------------------         
        
        # Urea dose (kgUrea/ha)
        DoseUrea=psl['bNA_DoseUrea_Standard']
                   
        # Nitrogen dose (kgN/ha)
        DoseN=psl['bNA_Ratio_N_to_Urea']*DoseUrea
            
        # Emissions from production of ammonia (tCO2e/ha)
        Ratio_UreaToAmmonia=psl['bNA_MolecularWeight_Urea']/psl['bNA_MolecularWeight_NH3']
        E_ProdNH3_per_t=psl['bNA_EmissionFromAmmoniaProduction_NRCAN']/Ratio_UreaToAmmonia
        E_ProdNH3=E_ProdNH3_per_t*(DoseUrea/1000)
            
        # Emissions from production of Urea from ammonia (tCO2e/ha)    
        MMBtu_per_app=(DoseUrea/1000)*psl['bNA_UreaEnergyConsumption']
        therm_per_app=MMBtu_per_app/psl['bNA_MMBtu_per_therm']
        E_ProdUrea=psl['bNA_EmissionFromUreaProduction_per_therm']*therm_per_app
            
        # Emissions from operations (tCO2e/ha)
        E_Ops=psl['bNA_EmissionFromLogisticsPreApplication_SP10']+ \
        psl['bNA_EmissionFromLogisticsApplication_SP10']

        # Volatilization: CO2 emissions following application, Tier 1 
        # approach, IPCC 2006, 11.4.1 (tCO2e/ha)
        # 0.2*(430/1000)*(1/0.27) = 0.32 tCO2e/ha
        E_Vol=psl['bNA_Ratio_C_to_Urea']*(DoseUrea/1000)*(1/psl['bRatio_C_to_CO2'])

        # Denitrification: N2O emissions following application, Tier 1 
        # approach, IPCC 2006, 11.4.1 (tCO2e/ha)
        E_Denit=psl['bNA_EmissionFactor_N2O_Jassaletal2008']*(DoseN/1000)* \
            psl['bNA_Ratio_N2OAsN_to_N2O']*psl['bGWP_N2O_AR4']
                
        # Total emissions (tCO2e/ha)
        E_Tot=E_ProdNH3+E_ProdUrea+E_Ops+E_Vol+E_Denit
                
        # Total emissions of carbon (MgC/ha) emitted as CO2e, converted
        # to carbon to be consistent with the rest of the variables in 
        # the vo object.
        vo['C_E_Operations'][iT,meta['NM']['iApplication']]=psl['bRatio_C_to_CO2']*E_Tot
    
    elif (comp=='HeterotrophicRespiration') & (meta['Fertilization Source']=='CBRunner'):
        
        #----------------------------------------------------------------------
        # Adjust rate of heterotrophic respiration
        #----------------------------------------------------------------------
        
        rr=psl['bNA_r_Decomp']
        meta['R_LitterVF'][0,meta['NM']['iApplication']]=rr*meta['R_LitterVF'][0,meta['NM']['iApplication']]
        meta['R_LitterF'][0,meta['NM']['iApplication']]=rr*meta['R_LitterF'][0,meta['NM']['iApplication']]
        meta['R_LitterM'][0,meta['NM']['iApplication']]=rr*meta['R_LitterM'][0,meta['NM']['iApplication']]
        meta['R_LitterS'][0,meta['NM']['iApplication']]=rr*meta['R_LitterS'][0,meta['NM']['iApplication']]
        meta['R_SoilVF'][0,meta['NM']['iApplication']]=rr*meta['R_SoilVF'][0,meta['NM']['iApplication']]
        meta['R_SoilF'][0,meta['NM']['iApplication']]=rr*meta['R_SoilF'][0,meta['NM']['iApplication']]
        meta['R_SoilS'][0,meta['NM']['iApplication']]=rr*meta['R_SoilS'][0,meta['NM']['iApplication']]
    
    return vi,vo,meta
