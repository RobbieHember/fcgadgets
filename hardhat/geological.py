
#%% Import modules

import numpy as np
import time
from fcgadgets.macgyver import utilities_general as gu
from fcgadgets.cbrunner import cbrun_utilities as cbu

#%% Calculate fossil fuel GHG emissions during operations

def FossilFuelEmissions(meta,pNam,vi,vo):

    # Extract parameters
    bB=meta['Param']['BEV']['Biophysical']

    # Total removals (m3/ha)
    V_Tot=vo['V_ToMillMerchTotal']+vo['V_ToMillNonMerch']

    # Total removals (ODT/ha)
    ODT_Tot=V_Tot*bB['Density Wood']

    #--------------------------------------------------------------------------
    # Resource extraction
    #--------------------------------------------------------------------------

    # Construction and maintenance of roads
    vo['E_CO2e_ET_OperForBurnOil']=vo['E_CO2e_ET_OperForBurnOil']+bB['EI Road Construction']*V_Tot

    # Cruising and reconnaissance
    vo['E_CO2e_ET_OperForBurnOil']=vo['E_CO2e_ET_OperForBurnOil']+bB['EI Cruise And Recon']*V_Tot

    # Felling and processing logs
    vo['E_CO2e_ET_OperForBurnOil']=vo['E_CO2e_ET_OperForBurnOil']+bB['EI Felling Process Logs']*V_Tot

    # Skidding trees to landing
    vo['E_CO2e_ET_OperForBurnOil']=vo['E_CO2e_ET_OperForBurnOil']+bB['EI Skidding To Landing']*V_Tot

    # Piling and sorting
    vo['E_CO2e_ET_OperForBurnOil']=vo['E_CO2e_ET_OperForBurnOil']+bB['EI Piling And Sorting Logs']*V_Tot

    # Loading logs
    vo['E_CO2e_ET_OperForBurnOil']=vo['E_CO2e_ET_OperForBurnOil']+bB['EI Loading Logs At Landing']*V_Tot

    # Chipping (non-merch)
    vo['E_CO2e_ET_OperForBurnOil']=vo['E_CO2e_ET_OperForBurnOil']+bB['EI Chipping']*vo['V_ToMillNonMerch']

    # Site preparation
    vo['E_CO2e_ET_OperForBurnOil']=vo['E_CO2e_ET_OperForBurnOil']+bB['EI Site Prep']*V_Tot

    # Sowing seeds
    vo['E_CO2e_ET_OperForBurnOil']=vo['E_CO2e_ET_OperForBurnOil']+bB['EI Sowing']*V_Tot

    # Planting
    vo['E_CO2e_ET_OperForBurnOil']=vo['E_CO2e_ET_OperForBurnOil']+bB['EI Planting']*V_Tot

    # Surveying
    vo['E_CO2e_ET_OperForBurnOil']=vo['E_CO2e_ET_OperForBurnOil']+bB['EI Surveying']*V_Tot

    # Hauling (forest ecosystem to mill)
    vo['E_CO2e_ET_OperForBurnOil']=vo['E_CO2e_ET_OperForBurnOil']+ODT_Tot/bB['Moisture Content Wood']*bB['Emission Intensity Transport Truck']*bB['Distance Forest To Mill (One Way)']

    #--------------------------------------------------------------------------
    # Mill operations
    #--------------------------------------------------------------------------

    # Unloading
    vo['E_CO2e_ET_OperForBurnOil']=vo['E_CO2e_ET_OperForBurnOil']+bB['EI Unloading At Mill']*V_Tot

    # Sawing and processing lumber
    vo['E_CO2e_ESC_OperForBurnOil']=vo['E_CO2e_ESC_OperForBurnOil']+bB['EI Sawing Processing Lumber']*vo['C_ToLumber']/bB['Carbon Content Wood']

    # Sawing and processing plywood
    vo['E_CO2e_ESC_OperForBurnOil']=vo['E_CO2e_ESC_OperForBurnOil']+bB['EI Sawing Processing Plywood']*vo['C_ToPlywood']/bB['Carbon Content Wood']

    # Sawing and processing OSB
    vo['E_CO2e_ESC_OperForBurnOil']=vo['E_CO2e_ESC_OperForBurnOil']+bB['EI Sawing Processing OSB']*vo['C_ToOSB']/bB['Carbon Content Wood']

    # Sawing and processing MDF
    vo['E_CO2e_ESC_OperForBurnOil']=vo['E_CO2e_ESC_OperForBurnOil']+bB['EI Sawing Processing MDF']*vo['C_ToMDF']/bB['Carbon Content Wood']

    # Sawing and processing pulp
    vo['E_CO2e_ESC_OperForBurnOil']=vo['E_CO2e_ESC_OperForBurnOil']+bB['EI Processing Pulp']*vo['C_ToPaper']/bB['Carbon Content Wood']

    # Pellets
    C_Pellet=vo['C_ToPelletExport']+vo['C_ToPelletDomGrid']+vo['C_ToPelletDomRNG']

    # Size reduction of pellets
    vo['E_CO2e_ESC_OperForBurnOil']=vo['E_CO2e_ESC_OperForBurnOil']+bB['EI Pellet Size Reduction']*C_Pellet/bB['Carbon Content Wood']

    # Drying pellets
    vo['E_CO2e_ESC_OperForBurnOil']=vo['E_CO2e_ESC_OperForBurnOil']+bB['EI Pellet Drying']*C_Pellet/bB['Carbon Content Wood']

    # Pelletizing
    vo['E_CO2e_ESC_OperForBurnOil']=vo['E_CO2e_ESC_OperForBurnOil']+bB['EI Pellet Pelletizing']*C_Pellet/bB['Carbon Content Wood']

    # Pellet sieving
    vo['E_CO2e_ESC_OperForBurnOil']=vo['E_CO2e_ESC_OperForBurnOil']+bB['EI Pellet Seiving']*C_Pellet/bB['Carbon Content Wood']

    #--------------------------------------------------------------------------
    # Transport Mill -> Distribution Hub (Vancouver)
    # *** Exclude raw log exports, as they originate at mills on the ocean ***
    #--------------------------------------------------------------------------

    Mass=(vo['C_ToLumber']+vo['C_ToPlywood']+vo['C_ToOSB']+vo['C_ToMDF']+vo['C_ToPaper']+vo['C_ToPelletExport']+vo['C_ToPelletDomRNG'])/bB['Carbon Content Wood']/bB['Moisture Content Lumber']

    vo['E_CO2e_ET_OperForBurnOil']=vo['E_CO2e_ET_OperForBurnOil']+bB['Distance Mill To Distribution Hub']*bB['Emission Intensity Transport Rail']*Mass

    #--------------------------------------------------------------------------
    # Transport Hub -> Market
    #--------------------------------------------------------------------------

    # Solid wood products

    Mass=(vo['C_ToLumber']+vo['C_ToPlywood']+vo['C_ToOSB']+vo['C_ToMDF']+vo['C_ToPaper'])/bB['Carbon Content Wood']/bB['Moisture Content Lumber']

    E_HubToMarket_Water=bB['Fraction Solid Wood Product Water Dest 1']*bB['Distance Solid Wood Product Water Dest 1']*bB['Emission Intensity Transport Water (Bulk)']*Mass+ \
        bB['Fraction Solid Wood Product Water Dest 1']*bB['Distance Solid Wood Product Water Dest 1']*bB['Emission Intensity Transport Water (Bulk)']*Mass+ \
        bB['Fraction Solid Wood Product Water Dest 1']*bB['Distance Solid Wood Product Water Dest 1']*bB['Emission Intensity Transport Water (Bulk)']*Mass

    E_HubToMarket_Rail=bB['Fraction Solid Wood Product Rail Dest 1']*bB['Distance Solid Wood Product Rail Dest 1']*bB['Emission Intensity Transport Rail']*Mass+ \
        bB['Fraction Solid Wood Product Rail Dest 1']*bB['Distance Solid Wood Product Rail Dest 1']*bB['Emission Intensity Transport Rail']*Mass+ \
        bB['Fraction Solid Wood Product Rail Dest 1']*bB['Distance Solid Wood Product Rail Dest 1']*bB['Emission Intensity Transport Rail']*Mass

    E_HubToMarket_Truck=bB['Fraction Solid Wood Product Truck Dest 1']*bB['Distance Solid Wood Product Truck Dest 1']*bB['Emission Intensity Transport Truck']*Mass+ \
        bB['Fraction Solid Wood Product Truck Dest 1']*bB['Distance Solid Wood Product Truck Dest 1']*bB['Emission Intensity Transport Truck']*Mass+ \
        bB['Fraction Solid Wood Product Truck Dest 1']*bB['Distance Solid Wood Product Truck Dest 1']*bB['Emission Intensity Transport Truck']*Mass

    vo['E_CO2e_ET_OperForBurnOil']=vo['E_CO2e_ET_OperForBurnOil']+E_HubToMarket_Water+E_HubToMarket_Rail+E_HubToMarket_Truck

    # Log exports

    Mass=(vo['C_ToLogExport'])/bB['Carbon Content Wood']/bB['Moisture Content Wood']

    vo['E_CO2e_ET_OperForBurnOil']=vo['E_CO2e_ET_OperForBurnOil']+bB['Distance LogExport Water Dest 1']*bB['Emission Intensity Transport Water (Bulk)']*Mass

    # Pellets

    Mass=vo['C_ToPelletExport']/bB['Carbon Content Wood']

    E_HubToMarket_Water=bB['Fraction PelletExport Water Dest 1']*bB['Distance PelletExport Water Dest 1']*bB['Emission Intensity Transport Water (Bulk)']*Mass+ \
        bB['Fraction PelletExport Water Dest 1']*bB['Distance PelletExport Water Dest 1']*bB['Emission Intensity Transport Water (Bulk)']*Mass+ \
        bB['Fraction PelletExport Water Dest 1']*bB['Distance PelletExport Water Dest 1']*bB['Emission Intensity Transport Water (Bulk)']*Mass

    E_HubToMarket_Rail=bB['Fraction PelletExport Rail Dest 1']*bB['Distance PelletExport Rail Dest 1']*bB['Emission Intensity Transport Rail']*Mass+ \
        bB['Fraction PelletExport Rail Dest 1']*bB['Distance PelletExport Rail Dest 1']*bB['Emission Intensity Transport Rail']*Mass+ \
        bB['Fraction PelletExport Rail Dest 1']*bB['Distance PelletExport Rail Dest 1']*bB['Emission Intensity Transport Rail']*Mass

    # E_HubToMarket_Truck=bB['Fraction PelletExport Truck Dest 1']*bB['Distance PelletExport Truck Dest 1']*bB['Emission Intensity Transport Truck']*Mass+ \
    #     bB['Fraction PelletExport Truck Dest 1']*bB['Distance PelletExport Truck Dest 1']*bB['Emission Intensity Transport Truck']*Mass+ \
    #     bB['Fraction PelletExport Truck Dest 1']*bB['Distance PelletExport Truck Dest 1']*bB['Emission Intensity Transport Truck']*Mass

    vo['E_CO2e_ET_OperForBurnOil']=vo['E_CO2e_ET_OperForBurnOil']+E_HubToMarket_Water+E_HubToMarket_Rail

    # *** OLD ***
    #EI_Harvest=meta['Param']['BEV']['Biophysical']['Emission Intensity of Harvesting']
    #E=EI_Harvest*(vo['V_ToMillMerchTotal'][iT,iHarvest]+vo['V_ToMillNonMerch'][iT,iHarvest])
    #vo['E_CO2e_ESC_OperationBurnOil'][iT,iHarvest]=0.5*E
    #vo['E_CO2e_ET_OperationBurnOil'][iT,iHarvest]=0.5*E

    return vo

#%% Calculate substitution effects

def SubstitutionEffects(meta,pNam,vi,vo):

    # Extract parameters
    bB=meta['Param']['BEV']['Biophysical']
    bS=meta['Param']['BEV']['Substitution']

    # Electricity conversion efficiency ratio
    Ratio_EC=bB['Electrical Conversion Efficiency of Pellet Electricity Plant (>25MW)']/bB['Electrical Conversion Efficiency of Coal Electricity Plant']

    #----------------------------------------------------------------------
    # Domestic facility power generation (MgC/ha) to (green tonne/ha)
    #----------------------------------------------------------------------

    Yield_PowerFacilityDom=vo['C_ToPowerFacilityDom']/bB['Density Wood']/bB['Moisture Content Wood']

    # Yield to energy (GJ/ha)
    GJ_PowerFacilityDom=bB['Energy Content Wood (0% moisture)']*Yield_PowerFacilityDom

    GJ_PowerFacilityDom=GJ_PowerFacilityDom*Ratio_EC

    E_Sub_CoalForBioenergy_PowerFacilityDom=bS['PowerFacilityDomFracDisplacingCoal']*bB['Emission Intensity Coal']/1000*GJ_PowerFacilityDom
    E_Sub_DieselForBioenergy_PowerFacilityDom=bS['PowerFacilityDomFracDisplacingDiesel']*bB['Emission Intensity Diesel']/1000*GJ_PowerFacilityDom
    E_Sub_GasForBioenergy_PowerFacilityDom=bS['PowerFacilityDomFracDisplacingNaturalGas']*bB['Emission Intensity Natural Gas']/1000*GJ_PowerFacilityDom
    E_Sub_OilForBioenergy_PowerFacilityDom=bS['PowerFacilityDomFracDisplacingOil']*bB['Emission Intensity Oil']/1000*GJ_PowerFacilityDom

    #----------------------------------------------------------------------
    # Foreign facility power generation (MgC/ha) to (green tonne/ha)
    #----------------------------------------------------------------------

    # Yield (green tonnes/ha)
    Yield_PowerFacilityExport=vo['C_ToPowerFacilityExport']/bB['Density Wood']/bB['Moisture Content Wood']

    # Yield to energy (GJ/ha)
    GJ_PowerFacilityExport=bB['Energy Content Wood (0% moisture)']*Yield_PowerFacilityExport

    GJ_PowerFacilityExport=GJ_PowerFacilityExport*Ratio_EC

    E_Sub_CoalForBioenergy_PowerFacilityExport=bS['PowerFacilityExportFracDisplacingCoal']*bB['Emission Intensity Coal']/1000*GJ_PowerFacilityExport
    E_Sub_DieselForBioenergy_PowerFacilityExport=bS['PowerFacilityExportFracDisplacingDiesel']*bB['Emission Intensity Diesel']/1000*GJ_PowerFacilityExport
    E_Sub_GasForBioenergy_PowerFacilityExport=bS['PowerFacilityExportFracDisplacingNaturalGas']*bB['Emission Intensity Natural Gas']/1000*GJ_PowerFacilityExport
    E_Sub_OilForBioenergy_PowerFacilityExport=bS['PowerFacilityExportFracDisplacingOil']*bB['Emission Intensity Oil']/1000*GJ_PowerFacilityExport

    #----------------------------------------------------------------------
    # Independent power producers (MgC/ha) to (green tonne/ha)
    #----------------------------------------------------------------------

    # Yield (green tonnes/ha)
    Yield_PowerGrid=vo['C_ToPowerGrid']/bB['Density Wood']/bB['Moisture Content Wood']

    # Yield to energy (GJ/ha)
    GJ_PowerGrid=bB['Energy Content Wood (0% moisture)']*Yield_PowerGrid

    GJ_PowerGrid=GJ_PowerGrid*Ratio_EC

    E_Sub_CoalForBioenergy_PowerGrid=bS['PowerGridFracDisplacingCoal']*bB['Emission Intensity Coal']/1000*GJ_PowerGrid
    E_Sub_DieselForBioenergy_PowerGrid=bS['PowerGridFracDisplacingDiesel']*bB['Emission Intensity Diesel']/1000*GJ_PowerGrid
    E_Sub_GasForBioenergy_PowerGrid=bS['PowerGridFracDisplacingNaturalGas']*bB['Emission Intensity Natural Gas']/1000*GJ_PowerGrid
    E_Sub_OilForBioenergy_PowerGrid=bS['PowerGridFracDisplacingOil']*bB['Emission Intensity Oil']/1000*GJ_PowerGrid

    #----------------------------------------------------------------------
    # Pellet exports (MgC/ha) to (kiln dried tonne/ha)
    #----------------------------------------------------------------------

    # Yield (kiln dired tonnes/ha)
    Yield_PelletExport=vo['C_ToPelletExport']/bB['Density Wood']

    # Yield to energy (GJ/ha)
    GJ_PelletExport=bB['Energy Content Wood (Kiln-dried)']*Yield_PelletExport

    GJ_PelletExport=GJ_PelletExport*Ratio_EC

    E_Sub_CoalForBioenergy_PelletExport=bS['PelletExportFracDisplacingCoal']*bB['Emission Intensity Coal']/1000*GJ_PelletExport
    E_Sub_DieselForBioenergy_PelletExport=bS['PelletExportFracDisplacingDiesel']*bB['Emission Intensity Diesel']/1000*GJ_PelletExport
    E_Sub_GasForBioenergy_PelletExport=bS['PelletExportFracDisplacingNaturalGas']*bB['Emission Intensity Natural Gas']/1000*GJ_PelletExport
    E_Sub_OilForBioenergy_PelletExport=bS['PelletExportFracDisplacingOil']*bB['Emission Intensity Oil']/1000*GJ_PelletExport

    #----------------------------------------------------------------------
    # Pellet domestic grid (MgC/ha) to (kiln dried tonne/ha)
    #----------------------------------------------------------------------

    # Yield (kiln dired tonnes/ha)
    Yield_PelletDomGrid=vo['C_ToPelletDomGrid']/bB['Density Wood']

    # Yield to energy (GJ/ha)
    GJ_PelletDomGrid=bB['Energy Content Wood (Kiln-dried)']*Yield_PelletDomGrid

    GJ_PelletDomGrid=GJ_PelletDomGrid*Ratio_EC

    E_Sub_CoalForBioenergy_PelletDomGrid=bS['PelletDomGridFracDisplacingCoal']*bB['Emission Intensity Coal']/1000*GJ_PelletDomGrid
    E_Sub_DieselForBioenergy_PelletDomGrid=bS['PelletDomGridFracDisplacingDiesel']*bB['Emission Intensity Diesel']/1000*GJ_PelletDomGrid
    E_Sub_GasForBioenergy_PelletDomGrid=bS['PelletDomGridFracDisplacingNaturalGas']*bB['Emission Intensity Natural Gas']/1000*GJ_PelletDomGrid
    E_Sub_OilForBioenergy_PelletDomGrid=bS['PelletDomGridFracDisplacingOil']*bB['Emission Intensity Oil']/1000*GJ_PelletDomGrid

    #----------------------------------------------------------------------
    # Pellet domestic RNG (MgC/ha) to (kiln dried tonne/ha)
    #----------------------------------------------------------------------

    # Yield (kiln dired tonnes/ha)
    Yield_PelletDomRNG=vo['C_ToPelletDomRNG']/bB['Density Wood']

    # Yield to energy (GJ/ha)
    GJ_PelletDomRNG=bB['Energy Content Wood (Kiln-dried)']*Yield_PelletDomRNG

    E_Sub_CoalForBioenergy_PelletDomRNG=bS['PelletDomRNGFracDisplacingCoal']*bB['Emission Intensity Coal']/1000*GJ_PelletDomRNG
    E_Sub_DieselForBioenergy_PelletDomRNG=bS['PelletDomRNGFracDisplacingDiesel']*bB['Emission Intensity Diesel']/1000*GJ_PelletDomRNG
    E_Sub_GasForBioenergy_PelletDomRNG=bS['PelletDomRNGFracDisplacingNaturalGas']*bB['Emission Intensity Natural Gas']/1000*GJ_PelletDomRNG
    E_Sub_OilForBioenergy_PelletDomRNG=bS['PelletDomRNGFracDisplacingOil']*bB['Emission Intensity Oil']/1000*GJ_PelletDomRNG

    #----------------------------------------------------------------------
    # Domestic firewood (MgC/ha) to (air dried tonne/ha)
    #----------------------------------------------------------------------

    # Yield (air dried tonnes/ha)
    Yield_FirewoodDom=vo['C_ToFirewoodDom']/bB['Density Wood']

    # Yield to energy (GJ/ha)
    GJ_FirewoodDom=bB['Energy Content Wood (0% moisture)']*Yield_FirewoodDom

    E_Sub_CoalForBioenergy_FirewoodDom=bS['FirewoodDomFracDisplacingCoal']*bB['Emission Intensity Coal']/1000*GJ_FirewoodDom
    E_Sub_DieselForBioenergy_FirewoodDom=bS['FirewoodDomFracDisplacingDiesel']*bB['Emission Intensity Diesel']/1000*GJ_FirewoodDom
    E_Sub_GasForBioenergy_FirewoodDom=bS['FirewoodDomFracDisplacingNaturalGas']*bB['Emission Intensity Natural Gas']/1000*GJ_FirewoodDom
    E_Sub_OilForBioenergy_FirewoodDom=bS['FirewoodDomFracDisplacingOil']*bB['Emission Intensity Oil']/1000*GJ_FirewoodDom

    #----------------------------------------------------------------------
    # Foreign firewood (MgC/ha) to (air dried tonne/ha)
    #----------------------------------------------------------------------

    # Yield (air dried tonnes/ha)
    Yield_FirewoodExport=vo['C_ToFirewoodExport']/bB['Density Wood']

    # Yield to energy (GJ/ha)
    GJ_FirewoodExport=bB['Energy Content Wood (0% moisture)']*Yield_FirewoodExport

    E_Sub_CoalForBioenergy_FirewoodExport=bS['FirewoodExportFracDisplacingCoal']*bB['Emission Intensity Coal']/1000*GJ_FirewoodExport
    E_Sub_DieselForBioenergy_FirewoodExport=bS['FirewoodExportFracDisplacingDiesel']*bB['Emission Intensity Diesel']/1000*GJ_FirewoodExport
    E_Sub_GasForBioenergy_FirewoodExport=bS['FirewoodExportFracDisplacingNaturalGas']*bB['Emission Intensity Natural Gas']/1000*GJ_FirewoodExport
    E_Sub_OilForBioenergy_FirewoodExport=bS['FirewoodExportFracDisplacingOil']*bB['Emission Intensity Oil']/1000*GJ_FirewoodExport

    #--------------------------------------------------------------------------
    # Substitution of fossil fuels for bioenergy
    # *** Save as positive and then change sign in post-processing ***
    #--------------------------------------------------------------------------

    vo['E_CO2e_SUB_CoalForBioenergy']= \
        (E_Sub_CoalForBioenergy_PowerFacilityDom + \
         E_Sub_CoalForBioenergy_PowerFacilityExport + \
         E_Sub_CoalForBioenergy_PowerGrid + \
         E_Sub_CoalForBioenergy_PelletExport + \
         E_Sub_CoalForBioenergy_PelletDomGrid + \
         E_Sub_CoalForBioenergy_PelletDomRNG + \
         E_Sub_CoalForBioenergy_FirewoodDom + \
         E_Sub_CoalForBioenergy_FirewoodExport)

    vo['E_CO2e_SUB_OilForBioenergy']= \
        (E_Sub_OilForBioenergy_PowerFacilityDom + \
         E_Sub_OilForBioenergy_PowerFacilityExport + \
         E_Sub_OilForBioenergy_PowerGrid + \
         E_Sub_OilForBioenergy_PelletExport + \
         E_Sub_OilForBioenergy_FirewoodDom + \
         E_Sub_OilForBioenergy_FirewoodExport + \
         E_Sub_DieselForBioenergy_PowerFacilityDom + \
         E_Sub_DieselForBioenergy_PowerFacilityExport + \
         E_Sub_DieselForBioenergy_PowerGrid + \
         E_Sub_DieselForBioenergy_PelletExport + \
         E_Sub_DieselForBioenergy_PelletDomGrid + \
         E_Sub_DieselForBioenergy_PelletDomRNG + \
         E_Sub_DieselForBioenergy_FirewoodDom + \
         E_Sub_DieselForBioenergy_FirewoodExport)

    vo['E_CO2e_SUB_GasForBioenergy']= \
        (E_Sub_GasForBioenergy_PowerFacilityDom + \
         E_Sub_GasForBioenergy_PowerFacilityExport + \
         E_Sub_GasForBioenergy_PowerGrid + \
         E_Sub_GasForBioenergy_PelletExport + \
         E_Sub_GasForBioenergy_PelletDomGrid + \
         E_Sub_GasForBioenergy_PelletDomRNG + \
         E_Sub_GasForBioenergy_FirewoodDom + \
         E_Sub_GasForBioenergy_FirewoodExport)

    vo['E_CO2e_SUB_PowerFacilityDom']= \
        (E_Sub_CoalForBioenergy_PowerFacilityDom + \
         E_Sub_OilForBioenergy_PowerFacilityDom + \
         E_Sub_GasForBioenergy_PowerFacilityDom)

    vo['E_CO2e_SUB_PowerFacilityExport']= \
        (E_Sub_CoalForBioenergy_PowerFacilityExport + \
         E_Sub_OilForBioenergy_PowerFacilityExport + \
         E_Sub_GasForBioenergy_PowerFacilityExport)

    vo['E_CO2e_SUB_PowerGrid']= \
        (E_Sub_CoalForBioenergy_PowerGrid + \
         E_Sub_OilForBioenergy_PowerGrid + \
         E_Sub_GasForBioenergy_PowerGrid)

    vo['E_CO2e_SUB_PelletExport']= \
        (E_Sub_CoalForBioenergy_PelletExport + \
         E_Sub_OilForBioenergy_PelletExport + \
         E_Sub_GasForBioenergy_PelletExport)

    vo['E_CO2e_SUB_PelletDomGrid']= \
        (E_Sub_CoalForBioenergy_PelletDomGrid + \
         E_Sub_OilForBioenergy_PelletDomGrid + \
         E_Sub_GasForBioenergy_PelletDomGrid)

    vo['E_CO2e_SUB_PelletDomRNG']= \
        (E_Sub_CoalForBioenergy_PelletDomRNG + \
         E_Sub_OilForBioenergy_PelletDomRNG + \
         E_Sub_GasForBioenergy_PelletDomRNG)

    vo['E_CO2e_SUB_FirewoodDom']= \
        (E_Sub_CoalForBioenergy_FirewoodDom + \
         E_Sub_OilForBioenergy_FirewoodDom + \
         E_Sub_GasForBioenergy_FirewoodDom)

    vo['E_CO2e_SUB_FirewoodExport']= \
        (E_Sub_CoalForBioenergy_FirewoodExport + \
         E_Sub_OilForBioenergy_FirewoodExport + \
         E_Sub_GasForBioenergy_FirewoodExport)

    #----------------------------------------------------------------------
    # Substitution for structural wood produts (tCO2e/ha)
    #----------------------------------------------------------------------

    # Sawnwood (t DM)
    ODT_Sawnwood=(1/bB['Carbon Content Wood'])*vo['C_ToLumber']

    # Panel (t DM)
    ODT_Panel=(1/bB['Carbon Content Wood'])*(vo['C_ToPlywood']+vo['C_ToOSB']+vo['C_ToMDF'])

    # Residuals (tonnes)
    #Residuals=0

    # Substitution of concrete for structural wood
    fS_Concrete=bS['SawnwoodFracDisplacingConcrete']*bS['DisplacementRatio_ConcreteForSawnwood']*ODT_Sawnwood
    fP_Concrete=bS['PanelFracDisplacingConcrete']*bS['DisplacementRatio_ConcreteForPanel']*ODT_Panel
    #fR=0#bS['ResidualsFracDisplacingConcrete']*bS['DisplacementRatio_ConcreteForResiduals']*Residuals
    vo['ODT Concrete']=fS_Concrete+fP_Concrete#-fR

    # Substitution of steel for structural wood
    fS_Steel=bS['SawnwoodFracDisplacingSteel']*bS['DisplacementRatio_SteelForSawnwood']*ODT_Sawnwood
    fP_Steel=bS['PanelFracDisplacingSteel']*bS['DisplacementRatio_SteelForPanel']*ODT_Panel
    #fR=0#bS['ResidualsFracDisplacingSteel']*bS['DisplacementRatio_SteelForResiduals']*Residuals
    vo['ODT Steel']=fS_Steel+fP_Steel#-fR

    # Substitution of aluminum for structural wood
    fS_Aluminum=bS['SawnwoodFracDisplacingAluminum']*bS['DisplacementRatio_AluminumForSawnwood']*ODT_Sawnwood
    fP_Aluminum=bS['PanelFracDisplacingAluminum']*bS['DisplacementRatio_AluminumForPanel']*ODT_Panel
    #fR=0#bS['ResidualsFracDisplacingAluminum']*bS['DisplacementRatio_AluminumForResiduals']*Residuals
    vo['ODT Aluminum']=fS_Aluminum+fP_Aluminum#-fR

    # Substitution of plastics for structural wood
    fS_Plastic=bS['SawnwoodFracDisplacingPlastic']*bS['DisplacementRatio_PlasticForSawnwood']*ODT_Sawnwood
    fP_Plastic=bS['PanelFracDisplacingPlastic']*bS['DisplacementRatio_PlasticForPanel']*ODT_Panel
    #fR=0#bS['ResidualsFracDisplacingPlastic']*bS['DisplacementRatio_PlasticForResiduals']*Residuals
    vo['ODT Plastic']=fS_Plastic+fP_Plastic#-fR

    # Substitution of textiles for structural wood
    fS_Textile=bS['SawnwoodFracDisplacingTextile']*bS['DisplacementRatio_TextileForSawnwood']*ODT_Sawnwood
    fP_Textile=bS['PanelFracDisplacingTextile']*bS['DisplacementRatio_TextileForPanel']*ODT_Panel
    #fR=0#bS['ResidualsFracDisplacingTextile']*bS['DisplacementRatio_TextileForResiduals']*Residuals
    vo['ODT Textile']=fS_Textile+fP_Textile#-fR

    # Emissions from productoin of structural materials
    vo['E_CO2e_SUB_Concrete']=bB['Emission Intensity Concrete']*vo['ODT Concrete']
    vo['E_CO2e_SUB_Steel']=bB['Emission Intensity Steel']*vo['ODT Steel']
    vo['E_CO2e_SUB_Aluminum']=bB['Emission Intensity Aluminum']*vo['ODT Aluminum']
    vo['E_CO2e_SUB_Plastic']=bB['Emission Intensity Plastic']*vo['ODT Plastic']
    vo['E_CO2e_SUB_Textile']=bB['Emission Intensity Textile']*vo['ODT Textile']

    # Emissions from sawnwood and Panel
    vo['E_CO2e_SUB_Sawnwood']=bB['Emission Intensity Concrete']*fS_Concrete+ \
        bB['Emission Intensity Steel']*fS_Steel + \
        bB['Emission Intensity Aluminum']*fS_Aluminum + \
        bB['Emission Intensity Plastic']*fS_Plastic + \
        bB['Emission Intensity Textile']*fS_Textile

    vo['E_CO2e_SUB_Panel']=bB['Emission Intensity Concrete']*fP_Concrete+ \
        bB['Emission Intensity Steel']*fP_Steel + \
        bB['Emission Intensity Aluminum']*fP_Aluminum + \
        bB['Emission Intensity Plastic']*fP_Plastic + \
        bB['Emission Intensity Textile']*fP_Textile

    # Emissions from structural matierials, tallied by feedstock (and calcination)
    vo['E_CO2e_SUB_CoalForWood']= \
        bS['FracConcreteEmissionsFromCoal']*vo['E_CO2e_SUB_Concrete']+ \
        bS['FracSteelEmissionsFromCoal']*vo['E_CO2e_SUB_Steel']+ \
        bS['FracAluminumEmissionsFromCoal']*vo['E_CO2e_SUB_Aluminum']+ \
        bS['FracPlasticEmissionsFromCoal']*vo['E_CO2e_SUB_Plastic']+ \
        bS['FracTextileEmissionsFromCoal']*vo['E_CO2e_SUB_Textile']

    vo['E_CO2e_SUB_OilForWood']= \
        bS['FracConcreteEmissionsFromOil']*vo['E_CO2e_SUB_Concrete']+ \
        bS['FracSteelEmissionsFromOil']*vo['E_CO2e_SUB_Steel']+ \
        bS['FracAluminumEmissionsFromOil']*vo['E_CO2e_SUB_Aluminum']+ \
        bS['FracPlasticEmissionsFromOil']*vo['E_CO2e_SUB_Plastic']+ \
        bS['FracTextileEmissionsFromOil']*vo['E_CO2e_SUB_Textile']

    vo['E_CO2e_SUB_GasForWood']= \
        bS['FracConcreteEmissionsFromGas']*vo['E_CO2e_SUB_Concrete']+ \
        bS['FracSteelEmissionsFromGas']*vo['E_CO2e_SUB_Steel']+ \
        bS['FracAluminumEmissionsFromGas']*vo['E_CO2e_SUB_Aluminum']+ \
        bS['FracPlasticEmissionsFromGas']*vo['E_CO2e_SUB_Plastic']+ \
        bS['FracTextileEmissionsFromGas']*vo['E_CO2e_SUB_Textile']

    vo['E_CO2e_SUB_Calcination']=bS['FracConcreteEmissionsFromCalcination']*vo['E_CO2e_SUB_Concrete']

    #--------------------------------------------------------------------------
    # Back-calculate production of fossil fuel consumption from operational use
    # and substitution effects (GJ)
    # *** Moved to post-processing ***
    #--------------------------------------------------------------------------

    return vo