'''

BCHWP Batch - Version 1 

A Python version of the British Columbia Harvested Wood Product Calculator 
(Dymond 2012).

Notes:

'''

def BCHWP_Batch(meta,vi):
    
    # Initialize output variables
    vo={}
    vo['C_Pro_Pools']=np.zeros( (meta['N Time'],meta['N Space'],meta['N Pools']) )
    
    # Loop through calendar years
    for iT in range(vi['Time'].size):
    
        #--------------------------------------------------------------------------
        # Ecosystems --> Mills
        # Note: Mills act as transient reservoirs that are not tracked over time
        #--------------------------------------------------------------------------
    
        Mill_FuelFromChips=psl['HWP']['RemovedMerchToFuel']*vo['C_RemovedMerch'][iT,:] + \
                 psl['HWP']['RemovedNonMerchToFuel']*vo['C_RemovedNonMerch'][iT,:] + \
                 psl['HWP']['RemovedSnagStemToFuel']*vo['C_RemovedSnagStem'][iT,:]    
        Mill_PulpFromChips=psl['HWP']['RemovedMerchToPulp']*vo['C_RemovedMerch'][iT,:] + \
                 psl['HWP']['RemovedNonMerchToPulp']*vo['C_RemovedNonMerch'][iT,:] + \
                 psl['HWP']['RemovedSnagStemToPulp']*vo['C_RemovedSnagStem'][iT,:]    
        Mill_Lumber=psl['HWP']['RemovedMerchToLumber']*vo['C_RemovedMerch'][iT,:] + \
                 psl['HWP']['RemovedNonMerchToLumber']*vo['C_RemovedNonMerch'][iT,:] + \
                 psl['HWP']['RemovedSnagStemToLumber']*vo['C_RemovedSnagStem'][iT,:]
        Mill_Plywood=psl['HWP']['RemovedMerchToPlywood']*vo['C_RemovedMerch'][iT,:] + \
                 psl['HWP']['RemovedNonMerchToPlywood']*vo['C_RemovedNonMerch'][iT,:] + \
                 psl['HWP']['RemovedSnagStemToPlywood']*vo['C_RemovedSnagStem'][iT,:]
        Mill_OSB=psl['HWP']['RemovedMerchToOSB']*vo['C_RemovedMerch'][iT,:] + \
                 psl['HWP']['RemovedNonMerchToOSB']*vo['C_RemovedNonMerch'][iT,:] + \
                 psl['HWP']['RemovedSnagStemToOSB']*vo['C_RemovedSnagStem'][iT,:]
        Mill_MDF=psl['HWP']['RemovedMerchToMDF']*vo['C_RemovedMerch'][iT,:] + \
                 psl['HWP']['RemovedNonMerchToMDF']*vo['C_RemovedNonMerch'][iT,:] + \
                 psl['HWP']['RemovedSnagStemToMDF']*vo['C_RemovedSnagStem'][iT,:]
        Mill_Firewood=psl['HWP']['RemovedMerchToFirewood']*vo['C_RemovedMerch'][iT,:] + \
                 psl['HWP']['RemovedNonMerchToFirewood']*vo['C_RemovedNonMerch'][iT,:] + \
                 psl['HWP']['RemovedSnagStemToFirewood']*vo['C_RemovedSnagStem'][iT,:]
        Mill_Cants=psl['HWP']['RemovedMerchToCants']*vo['C_RemovedMerch'][iT,:] + \
                 psl['HWP']['RemovedNonMerchToCants']*vo['C_RemovedNonMerch'][iT,:] + \
                 psl['HWP']['RemovedSnagStemToCants']*vo['C_RemovedSnagStem'][iT,:]
     
        #--------------------------------------------------------------------------
        # Mills --> In-use pools or other mills
        #--------------------------------------------------------------------------
        
        # Transfer mill fibre to single-family homes
        ip=meta['iPP']['SFH']
        vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT-1,:,ip] + \
            psl['HWP']['LumberToSFH']*Mill_Lumber + \
            psl['HWP']['PlywoodToSFH']*Mill_Plywood + \
            psl['HWP']['OSBToSFH']*Mill_OSB + \
            psl['HWP']['MDFToSFH']*Mill_MDF
    
        # Transfer mill fibre to multi-family homes
        ip=meta['iPP']['MFH']
        vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT-1,:,ip] + \
            psl['HWP']['LumberToMFH']*Mill_Lumber + \
            psl['HWP']['PlywoodToMFH']*Mill_Plywood + \
            psl['HWP']['OSBToMFH']*Mill_OSB + \
            psl['HWP']['MDFToMFH']*Mill_MDF
    
        # Transfer mill fibre to commercial
        ip=meta['iPP']['Comm']
        vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT-1,:,ip] + \
            psl['HWP']['LumberToCom']*Mill_Lumber + \
            psl['HWP']['PlywoodToCom']*Mill_Plywood + \
            psl['HWP']['OSBToCom']*Mill_OSB + \
            psl['HWP']['MDFToCom']*Mill_MDF
    
        # Transfer mill fibre to furniture
        ip=meta['iPP']['Furn']
        vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT-1,:,ip] + \
            psl['HWP']['LumberToFurn']*Mill_Lumber + \
            psl['HWP']['PlywoodToFurn']*Mill_Plywood + \
            psl['HWP']['OSBToFurn']*Mill_OSB + \
            psl['HWP']['MDFToFurn']*Mill_MDF
    
        # Transfer mill fibre to shipping
        ip=meta['iPP']['Ship']
        vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT-1,:,ip] + \
            psl['HWP']['LumberToShip']*Mill_Lumber + \
            psl['HWP']['PlywoodToShip']*Mill_Plywood + \
            psl['HWP']['OSBToShip']*Mill_OSB + \
            psl['HWP']['MDFToShip']*Mill_MDF
    
        # Transfer mill fibre to repairs
        ip=meta['iPP']['Repairs']
        vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT-1,:,ip] + \
            psl['HWP']['LumberToRepairs']*Mill_Lumber + \
            psl['HWP']['PlywoodToRepairs']*Mill_Plywood + \
            psl['HWP']['OSBToRepairs']*Mill_OSB + \
            psl['HWP']['MDFToRepairs']*Mill_MDF
    
        # Transfer mill fibre to other
        ip=meta['iPP']['Other']
        vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT-1,:,ip] + \
            psl['HWP']['LumberToOther']*Mill_Lumber + \
            psl['HWP']['PlywoodToOther']*Mill_Plywood + \
            psl['HWP']['OSBToOther']*Mill_OSB + \
            psl['HWP']['MDFToOther']*Mill_MDF
    
        # Transfer mill fibre to pulp mill
        Mill_Pulp=Mill_PulpFromChips + \
            psl['HWP']['LumberToPulp']*Mill_Lumber + \
            psl['HWP']['PlywoodToPulp']*Mill_Plywood
                     
        # Transfer pulp mill fibre to paper
        ip=meta['iPP']['Paper']
        vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT-1,:,ip] + \
            psl['HWP']['PulpToPaper']*Mill_Pulp
    
        # Transfer mill fibre to fuel 
        ip=meta['iPP']['Fuel']
        vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT-1,:,ip] + \
            Mill_FuelFromChips + \
            psl['HWP']['LumberToFuel']*Mill_Lumber + \
            psl['HWP']['PlywoodToFuel']*Mill_Plywood + \
            psl['HWP']['OSBToFuel']*Mill_OSB + \
            psl['HWP']['MDFToFuel']*Mill_MDF
    
        # Transfer firewood to firewood pool
        ip=meta['iPP']['Firewood']
        vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT-1,:,ip] + Mill_Firewood
    
        # Transfer mill fibre to cants
        ip=meta['iPP']['Cants']
        vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT-1,:,ip] + Mill_Cants
    
        # Transfer pulp mill carbon to pulp-mill effluent
        ip=meta['iPP']['EffluentPulp']
        vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT-1,:,ip] + \
             psl['HWP']['PulpToEffluent']*Mill_Pulp   
    
        #--------------------------------------------------------------------------
        # Update dump and landfill reservoirs
        #--------------------------------------------------------------------------

        vo['C_Pro_Pools'][iT,:,meta['iPP']['DumpWood']]=vo['C_Pro_Pools'][iT-1,:,meta['iPP']['DumpWood']]
        vo['C_Pro_Pools'][iT,:,meta['iPP']['DumpPaper']]=vo['C_Pro_Pools'][iT-1,:,meta['iPP']['DumpPaper']]
        vo['C_Pro_Pools'][iT,:,meta['iPP']['LandfillWoodDegradable']]=vo['C_Pro_Pools'][iT-1,:,meta['iPP']['LandfillWoodDegradable']]
        vo['C_Pro_Pools'][iT,:,meta['iPP']['LandfillWoodNonDegradable']]=vo['C_Pro_Pools'][iT-1,:,meta['iPP']['LandfillWoodNonDegradable']]    
        vo['C_Pro_Pools'][iT,:,meta['iPP']['LandfillPaperDegradable']]=vo['C_Pro_Pools'][iT-1,:,meta['iPP']['LandfillPaperDegradable']]
        vo['C_Pro_Pools'][iT,:,meta['iPP']['LandfillPaperNonDegradable']]=vo['C_Pro_Pools'][iT-1,:,meta['iPP']['LandfillPaperNonDegradable']]    
        
        #--------------------------------------------------------------------------
        # Single-family homes --> dump and landfill
        #--------------------------------------------------------------------------
           
        # Turnover
        ip=meta['iPP']['SFH']
        C_retired=psl['HWP']['SFH_tr']*vo['C_Pro_Pools'][iT-1,:,ip]
    
        # Remove carbon
        vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] - C_retired
    
        # Transfer carbon to dump wood
        ip=meta['iPP']['DumpWood']
        vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + psl['HWP']['SFHToDumpWood']*C_retired
    
        # Transfer carbon to landfill (degradable)
        ip=meta['iPP']['LandfillWoodDegradable']
        vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + psl['HWP']['SFHToLandfillWood']*psl['HWP']['ToLandfillWoodDegradableFrac']*C_retired
    
        # Transfer carbon to landfill (non-degradable)
        ip=meta['iPP']['LandfillWoodNonDegradable']
        vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + psl['HWP']['SFHToLandfillWood']*(1-psl['HWP']['ToLandfillWoodDegradableFrac'])*C_retired
    
        #--------------------------------------------------------------------------
        # Multi-family homes --> dump and landfill
        #--------------------------------------------------------------------------
    
        # Turnover
        ip=meta['iPP']['MFH']
        C_retired=psl['HWP']['MFH_tr']*vo['C_Pro_Pools'][iT,:,ip]
    
        # Remove carbon
        vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] - C_retired
    
        # Transfer carbon to dump wood
        ip=meta['iPP']['DumpWood']
        vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + psl['HWP']['MFHToDumpWood']*C_retired
    
        # Transfer carbon to landfill (degradable)
        ip=meta['iPP']['LandfillWoodDegradable']
        vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + psl['HWP']['MFHToLandfillWood']*psl['HWP']['ToLandfillWoodDegradableFrac']*C_retired
    
        # Transfer carbon to landfill (non-degradable)
        ip=meta['iPP']['LandfillWoodNonDegradable']
        vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + psl['HWP']['MFHToLandfillWood']*(1-psl['HWP']['ToLandfillWoodDegradableFrac'])*C_retired
    
        #--------------------------------------------------------------------------
        # Commercial building --> dump and landfill
        #--------------------------------------------------------------------------
    
        # Turnover
        ip=meta['iPP']['Comm']
        C_retired=psl['HWP']['Comm_tr']*vo['C_Pro_Pools'][iT-1,:,ip]
    
        # Remove carbon
        vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] - C_retired
    
        # Transfer carbon to dump wood
        ip=meta['iPP']['DumpWood']
        vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + psl['HWP']['CommToDumpWood']*C_retired
    
        # Transfer carbon to landfill (degradable)
        ip=meta['iPP']['LandfillWoodDegradable']
        vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + psl['HWP']['CommToLandfillWood']*psl['HWP']['ToLandfillWoodDegradableFrac']*C_retired
    
        # Transfer carbon to landfill (non-degradable)
        ip=meta['iPP']['LandfillWoodNonDegradable']
        vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + psl['HWP']['CommToLandfillWood']*(1-psl['HWP']['ToLandfillWoodDegradableFrac'])*C_retired
    
        #--------------------------------------------------------------------------
        # Furniture --> dump and landfill
        #--------------------------------------------------------------------------
    
        # Turnover
        ip=meta['iPP']['Furn']
        C_retired=psl['HWP']['Furn_tr']*vo['C_Pro_Pools'][iT-1,:,ip]
    
        # Remove carbon
        vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] - C_retired
    
        # Transfer carbon to dump wood
        ip=meta['iPP']['DumpWood']
        vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + psl['HWP']['FurnToDumpWood']*C_retired
        
        # Transfer carbon to landfill (degradable)
        ip=meta['iPP']['LandfillWoodDegradable']
        vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + psl['HWP']['FurnToLandfillWood']*psl['HWP']['ToLandfillWoodDegradableFrac']*C_retired
    
        # Transfer carbon to landfill (non-degradable)
        ip=meta['iPP']['LandfillWoodNonDegradable']
        vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + psl['HWP']['FurnToLandfillWood']*(1-psl['HWP']['ToLandfillWoodDegradableFrac'])*C_retired
    
    
        #--------------------------------------------------------------------------
        # Shipping --> dump and landfill
        #--------------------------------------------------------------------------
    
        # Turnover
        ip=meta['iPP']['Ship']
        C_retired=psl['HWP']['Ship_tr']*vo['C_Pro_Pools'][iT-1,:,ip]
    
        # Remove carbon
        vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] - C_retired
        
        # Transfer carbon to dump wood
        ip=meta['iPP']['DumpWood']
        vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + psl['HWP']['ShipToDumpWood']*C_retired
    
        # Transfer carbon to landfill (degradable)
        ip=meta['iPP']['LandfillWoodDegradable']
        vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + psl['HWP']['ShipToLandfillWood']*psl['HWP']['ToLandfillWoodDegradableFrac']*C_retired
    
        # Transfer carbon to landfill (non-degradable)
        ip=meta['iPP']['LandfillWoodNonDegradable']
        vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + psl['HWP']['ShipToLandfillWood']*(1-psl['HWP']['ToLandfillWoodDegradableFrac'])*C_retired
    
        #--------------------------------------------------------------------------
        # Repairs --> dump and landfill
        #--------------------------------------------------------------------------
    
        # Turnover
        ip=meta['iPP']['Repairs']
        C_retired=psl['HWP']['Repairs_tr']*vo['C_Pro_Pools'][iT-1,:,ip]
    
        # Remove carbon
        vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] - C_retired
    
        # Transfer carbon to dump wood
        ip=meta['iPP']['DumpWood']
        vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + psl['HWP']['RepairsToDumpWood']*C_retired
    
        # Transfer carbon to landfill (degradble)
        ip=meta['iPP']['LandfillWoodDegradable']
        vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + psl['HWP']['RepairsToLandfillWood']*psl['HWP']['ToLandfillWoodDegradableFrac']*C_retired
    
        # Transfer carbon to landfill (non-degradable)
        ip=meta['iPP']['LandfillWoodNonDegradable']
        vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + psl['HWP']['RepairsToLandfillWood']*(1-psl['HWP']['ToLandfillWoodDegradableFrac'])*C_retired
    
        #--------------------------------------------------------------------------
        # Other --> dump and landfill
        #--------------------------------------------------------------------------
    
        # Turnover
        ip=meta['iPP']['Other']
        C_retired=psl['HWP']['Other_tr']*vo['C_Pro_Pools'][iT-1,:,ip]
    
        # Remove carbon
        vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] - C_retired
    
        # Transfer carbon to dump wood
        ip=meta['iPP']['DumpWood']
        vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + psl['HWP']['OtherToDumpWood']*C_retired
    
        # Transfer carbon to landfill (degradble)
        ip=meta['iPP']['LandfillWoodDegradable']
        vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + psl['HWP']['OtherToLandfillWood']*psl['HWP']['ToLandfillWoodDegradableFrac']*C_retired
    
        # Transfer carbon to landfill (non-degradable)
        ip=meta['iPP']['LandfillWoodNonDegradable']
        vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + psl['HWP']['OtherToLandfillWood']*(1-psl['HWP']['ToLandfillWoodDegradableFrac'])*C_retired
    
        #--------------------------------------------------------------------------
        # Cants --> dump and landfill
        #--------------------------------------------------------------------------
    
        # Turnover
        ip=meta['iPP']['Cants']
        C_retired=psl['HWP']['Cants_tr']*vo['C_Pro_Pools'][iT-1,:,ip]
    
        # Remove carbon
        vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] - C_retired
    
        # Transfer carbon to dump wood
        ip=meta['iPP']['DumpWood']
        vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + psl['HWP']['CantsToDumpWood']*C_retired
    
        # Transfer carbon to landfill (degradble)
        ip=meta['iPP']['LandfillWoodDegradable']
        vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + psl['HWP']['CantsToLandfillWood']*psl['HWP']['ToLandfillWoodDegradableFrac']*C_retired
    
        # Transfer carbon to landfill (non-degradable)
        ip=meta['iPP']['LandfillWoodNonDegradable']
        vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + psl['HWP']['CantsToLandfillWood']*(1-psl['HWP']['ToLandfillWoodDegradableFrac'])*C_retired
    
        #--------------------------------------------------------------------------
        # Paper --> dump and landfill
        #--------------------------------------------------------------------------
    
        # Turnover (with adjustment for recycling)
        ip=meta['iPP']['Paper']
        C_retired=(1-psl['HWP']['PaperRecycleRate'])*psl['HWP']['Paper_tr']*vo['C_Pro_Pools'][iT,:,ip]
    
        # Remove carbon
        vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] - C_retired
    
        # Transfer to dump
        ip=meta['iPP']['DumpPaper']
        vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + psl['HWP']['PaperToDumpPaper']*C_retired
    
        # Transfer to landfill (degradable)
        ip=meta['iPP']['LandfillWoodDegradable']
        vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + psl['HWP']['PaperToLandfillPaper']*psl['HWP']['ToLandfillPaperDegradableFrac']*C_retired
    
        # Transfer to landfill (non-degradable)
        ip=meta['iPP']['LandfillWoodNonDegradable']
        vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + psl['HWP']['PaperToLandfillPaper']*(1-psl['HWP']['ToLandfillPaperDegradableFrac'])*C_retired
        
        #--------------------------------------------------------------------------
        # Emissions from fuel combustion
        #--------------------------------------------------------------------------
    
        # Turnover
        ip=meta['iPP']['Fuel']
        C_retired=psl['HWP']['Fuel_tr']*vo['C_Pro_Pools'][iT,:,ip]
    
        # Remove carbon
        vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] - C_retired
    
        # Emissions of CO2 from fuel use
        ip=meta['iPP']['E_CO2']
        vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + psl['HWP']['FuelCombustionFracEmitCO2']*C_retired
    
        # Emissions of CH4 from fuel use
        ip=meta['iPP']['E_CH4']
        vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + (1-psl['HWP']['FuelCombustionFracEmitCO2'])*C_retired
    
        #--------------------------------------------------------------------------
        # Emissions from firewood combustion
        #--------------------------------------------------------------------------
    
        # Turnover
        ip=meta['iPP']['Firewood']
        C_retired=psl['HWP']['Firewood_tr']*vo['C_Pro_Pools'][iT,:,ip]
    
        # Remove carbon
        vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] - C_retired
        
        # Emissions of CO2
        ip=meta['iPP']['E_CO2']
        vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip] + C_retired
    
        # Emissions of CH4
        #ip=meta['iPP']['E_CH4']
        #vo['C_Pro_Pools[iT,:,ip]=vo['C_Pro_Pools[iT,:,ip] + (1-psl['FuelCombustionFracEmitCO2'])*C_retired
       
        #--------------------------------------------------------------------------
        # Emissions from pulp effluent
        #--------------------------------------------------------------------------
         
        # Emissions from pulp effluent (CO2 from aerobic decomposition)
        ip=meta['iPP']['EffluentPulp']
        c_emitted=psl['HWP']['EffluentPulp_tr']*vo['C_Pro_Pools'][iT,:,ip]
    
        # Remove emitted carbon
        vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip]-c_emitted
    
        # Add emitted carbon to CO2 emission "pool"
        ip=meta['iPP']['E_CO2']
        vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip]+c_emitted    
    
        #--------------------------------------------------------------------------
        # Emissions from dump wood
        #--------------------------------------------------------------------------
    
        # Turnover
        ip=meta['iPP']['DumpWood']
        c_emitted=psl['HWP']['DumpWood_tr']*vo['C_Pro_Pools'][iT,:,ip]
    
        # Removal
        vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip]-c_emitted
    
        # Add to emissions (CO2 emission from aerobic decomposition)
        ip=meta['iPP']['E_CO2']
        vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip]+c_emitted
    
        #--------------------------------------------------------------------------
        # Emissions from dump paper
        #--------------------------------------------------------------------------
    
        # Turnover
        ip=meta['iPP']['DumpPaper']
        c_emitted=psl['HWP']['DumpPaper_tr']*vo['C_Pro_Pools'][iT,:,ip]
    
        # Removal
        vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip]-c_emitted
    
        # Add to emissions (CO2 emission from aerobic decomposition)
        ip=meta['iPP']['E_CO2']
        vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip]+c_emitted
    
        #--------------------------------------------------------------------------
        # Emissions from landfill degradable wood
        #--------------------------------------------------------------------------
               
        # Turnover
        ip=meta['iPP']['LandfillWoodDegradable']
        c_emitted=psl['HWP']['LandfillWoodDegradable_tr']*vo['C_Pro_Pools'][iT,:,ip]
    
        # Removal
        vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip]-c_emitted
    
        # Add to emissions (50% CO2 emissions during anaerobic decomposition)
        ip=meta['iPP']['E_CO2']
        vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip]+psl['HWP']['LandfillDegradableFracEmitCO2']*c_emitted
    
        # Add to emissions (50% "potential" CH4 emissions during anaerobic decomposition)
        #E_ch4_pot=(1-psl['HWP']['LandfillDegradableFracEmitCO2)*c_emitted
    
        # Adjustment for proportion of degradable landfills with gas collection systems, 
        # efficiency of system, and methane oxided to CO2 from the landfill cover
        ch4_emitted=c_emitted*((1-psl['HWP']['LandfillMethaneEmit_GasColSysProp'])-psl['HWP']['LandfillMethaneOxidizedToCO2']*(1-psl['HWP']['LandfillMethaneEmit_GasColSysProp'])) + \
            c_emitted*psl['HWP']['LandfillMethaneEmit_GasColSysProp']*((1-psl['HWP']['LandfillMethaneEmit_GasColSysEffic'])-psl['HWP']['LandfillMethaneOxidizedToCO2']*(1-psl['HWP']['LandfillMethaneEmit_GasColSysProp']))
        
        ip=meta['iPP']['E_CH4']
        vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip]+ch4_emitted
    
        #--------------------------------------------------------------------------
        # Emissions from landfill degradable paper
        #--------------------------------------------------------------------------
    
        # Turnover
        ip=meta['iPP']['LandfillPaperDegradable']
        c_emitted=psl['HWP']['LandfillWoodDegradable_tr']*vo['C_Pro_Pools'][iT,:,ip]
    
        # Removal
        vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip]-c_emitted
        
        # Add to emissions (50% CO2 emissions during anaerobic decomposition)
        ip=meta['iPP']['E_CO2']
        vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip]+psl['HWP']['LandfillDegradableFracEmitCO2']*c_emitted
    
        # Add to emissions (50% "potential" CH4 emissions during anaerobic decomposition)
        #E_ch4_pot=(1-psl['HWP']['LandfillDegradableFracEmitCO2)*c_emitted
    
        # Adjustment for proportion of degradable landfills with gas collection systems, 
        # efficiency of system, and methane oxided to CO2 from the landfill cover
        ch4_emitted=c_emitted*((1-psl['HWP']['LandfillMethaneEmit_GasColSysProp'])-psl['HWP']['LandfillMethaneOxidizedToCO2']*(1-psl['HWP']['LandfillMethaneEmit_GasColSysProp'])) + \
            c_emitted*psl['HWP']['LandfillMethaneEmit_GasColSysProp']*((1-psl['HWP']['LandfillMethaneEmit_GasColSysEffic'])-psl['HWP']['LandfillMethaneOxidizedToCO2']*(1-psl['HWP']['LandfillMethaneEmit_GasColSysProp']))
        
        ip=meta['iPP']['E_CH4']
        vo['C_Pro_Pools'][iT,:,ip]=vo['C_Pro_Pools'][iT,:,ip]+ch4_emitted
    
    return meta,vo
