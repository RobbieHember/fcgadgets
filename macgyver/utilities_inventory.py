
'''

FCGADGETS - FOREST INVENTORY UTILITIES

===============================================================================
Instructions for how to update forest inventory data:
===============================================================================

 1) Manually download geodatabases from BCGW in ArcGIS (60 min)
    - if export not working, use copy and paste 
    - store files at the paths indicated in "DefineInventoryLayersAndVariables"

 2) Define inventory layers and variables (<1min))
    - Run: LayerInfo=invu.DefineInventoryLayersAndVariables()
    - make sure folder release dates are updated

 3) Build and save LUTs (20 min)
    - Run: invu.BuildForestInventoryLUTs(LayerInfo)

 4) Extract openings and forest cover polygons 
    - Run: RecoverMissingATUGeometries()
    - that can be used when activity spatial is missing
    
'''

#%% Import modules

import numpy as np
import pandas as pd
import geopandas as gpd
import fiona
import gdal
import time
import gc as garc
from fcgadgets.macgyver import utilities_gis as gis
from fcgadgets.macgyver import utilities_general as gu
from fcgadgets.cbrunner import cbrun_utilities as cbu

#%% Look at contents of geodatabase

# fiona.listlayers(r'C:\Users\rhember\Documents\Data\Basemaps\Basemaps.gdb')
# fiona.listlayers(r'C:\Users\rhember\Documents\Data\ForestInventory\LandUse\20210401\LandUse.gdb')
# fiona.listlayers(r'C:\Users\rhember\Documents\Data\ForestInventory\Results\20210401\Results.gdb')
# fiona.listlayers(r'C:\Users\rhember\Documents\Data\ForestInventory\Results\20210401\ForestCover.gdb')
# fiona.listlayers(r'C:\Users\rhember\Documents\Data\ForestInventory\VRI\20210401\VRI.gdb')
# fiona.listlayers(r'C:\Users\rhember\Documents\Data\ForestInventory\Disturbances\20210401\Disturbances.gdb')

#%% Define inventory layers and varialbes
# The "Field List" variable contains touples containing the variable name and 
# a flag indicating whether it is string (1) or numberic (0)

def DefineInventoryLayersAndVariables():

    # Define paths to geodatabase files
    PathInResultsFull=r'C:\Users\rhember\Documents\Data\ForestInventory\Results\20210930'
    PathInDisturbancesFull=r'C:\Users\rhember\Documents\Data\ForestInventory\Disturbances\20210930'
    PathInVRIFull=r'C:\Users\rhember\Documents\Data\ForestInventory\VRI\20210930'
    PathInLUPFull=r'C:\Users\rhember\Documents\Data\ForestInventory\LandUse\20210930'

    # Initialize a list of layer information
    LayerInfo=[]

    # Add layer and define required variables
    d={}
    d['Layer Name']='RSLT_ACTIVITY_TREATMENT_SVW'
    d['Path']=PathInResultsFull
    d['File Name']='Results.gdb'
    d['Field List']=[('ACTIVITY_TREATMENT_UNIT_ID',0,'float32'),
     ('OPENING_ID',0,'float32'), \
     ('SILV_BASE_CODE',1,'int16'), \
     ('SILV_TECHNIQUE_CODE',1,'int16'), \
     ('SILV_METHOD_CODE',1,'int16'), \
     ('SILV_OBJECTIVE_CODE_1',1,'int16'), \
     ('SILV_FUND_SOURCE_CODE',1,'int16'), \
     ('ATU_COMPLETION_DATE',0,'<U20'), \
     ('ACTUAL_TREATMENT_AREA',0,'float32'), \
     ('ACTUAL_TREATMENT_COST',0,'float32'), \
     ('ACTUAL_PLANTED_NUMBER',0,'float32'), \
     ('DISTURBANCE_CODE',1,'int16'), \
     ('SILV_SYSTEM_CODE',1,'int16'), \
     ('SILV_SYSTEM_VARIANT_CODE',1,'int16'), \
     ('FIA_PROJECT_ID',1,'int32'), \
     ('RESULTS_IND',0,'<U20')]
    d['LUT']={}
    for i in d['Field List']: d['LUT'][i[0]]=[] 
    LayerInfo.append(d)
    
    # Add layer and define required variables
    d={}
    d['Layer Name']='RSLT_OPENING_SVW'
    d['Path']=PathInResultsFull
    d['File Name']='Results.gdb'
    d['Field List']=[('OPENING_ID',0,'int32'), \
     ('REGION_NAME',1,'int16'), \
     ('DISTRICT_NAME',1,'int16'), \
     ('TIMBER_MARK',1,'int32'), \
     ('DISTURBANCE_START_DATE',0,'<U20'), \
     ('DISTURBANCE_END_DATE',0,'<U20'), \
     ('DENUDATION_1_DISTURBANCE_CODE',1,'int16'), \
     ('DENUDATION_1_SILV_SYSTEM_CODE',1,'int16'), \
     ('DENUDATION_1_SILV_VARIANT_CODE',1,'int16'), \
     ('DENUDATION_1_CUT_PHASE_CODE',1,'int16'), \
     ('DENUDATION_1_COMPLETION_DATE',0,'<U20'), \
     ('DENUDATION_2_DISTURBANCE_CODE',1,'int16'), \
     ('DENUDATION_2_SILV_SYSTEM_CODE',1,'int16'), \
     ('DENUDATION_2_SILV_VARIANT_CODE',1,'int16'), \
     ('DENUDATION_2_CUT_PHASE_CODE',1,'int16'), \
     ('DENUDATION_2_COMPLETION_DATE',0,'<U20'), \
     ('FEATURE_AREA',0,'float32'), \
     ('OPENING_GROSS_AREA',0,'float32'), \
     ('CUTTING_PERMIT_ID',1,'int16'), \
     ('TIMBER_MARK',1,'int16'), \
     ('PREV_TREE_SPECIES1_CODE',1,'int16'), \
     ('PREV_TREE_SPECIES2_CODE',1,'int16'), \
     ('PREV_STOCKING_STATUS_CODE',1,'int16'), \
     ('PREV_AGE_CLASS_CODE',1,'int16'), \
     ('PREV_HEIGHT_CLASS_CODE',1,'int16'), \
     ('PREV_SITE_INDEX',0,'float32'), \
     ('PREV_SITE_INDEX_SOURCE_CODE',1,'int16'), \
     ('MAX_ALLOW_PERMNT_ACCESS_PCT',0,'float32')]
    d['LUT']={}
    for i in d['Field List']: d['LUT'][i[0]]=[]
    LayerInfo.append(d)

    # Add layer and define required variables
    # *** Notes:OBSERVATION_DATE is empty ***
    d={}
    d['Layer Name']='RSLT_FOREST_COVER_INV_SVW'
    d['Path']=PathInResultsFull
    d['File Name']='Results.gdb'
    d['Field List']=[('OPENING_ID',0,'float32'), \
     ('FOREST_COVER_ID',0,'float32'), \
     ('I_FOREST_COVER_LAYER_ID',0,'float32'), \
     ('SILV_POLYGON_NUMBER',1,'int32'), \
     ('SILV_POLYGON_AREA',0,'float32'), \
     ('STOCKING_STANDARD_UNIT_ID',0,'float32'), \
     ('BGC_ZONE_CODE',1,'int16'), \
     ('BGC_SUBZONE_CODE',1,'int16'), \
     ('BGC_VARIANT',1,'int16'), \
     ('BGC_PHASE',1,'int16'), \
     ('BEC_SITE_SERIES',1,'int16'), \
     ('SITE_INDEX',0,'float32'), \
     ('I_SPECIES_AGE_1',0,'float32'), \
     ('I_SPECIES_CODE_1',1,'int16'), \
     ('I_SPECIES_CODE_2',1,'int16'), \
     ('I_SPECIES_CODE_3',1,'int16'), \
     ('I_SPECIES_CODE_4',1,'int16'), \
     ('I_SPECIES_CODE_5',1,'int16'), \
     ('I_SPECIES_PERCENT_1',0,'int16'), \
     ('I_SPECIES_PERCENT_2',0,'int16'), \
     ('I_SPECIES_PERCENT_3',0,'int16'), \
     ('I_SPECIES_PERCENT_4',0,'int16'), \
     ('I_SPECIES_PERCENT_5',0,'int16'), \
     ('I_TOTAL_STEMS_PER_HA',0,'float32'), \
     ('I_TOTAL_WELL_SPACED_STEMS_HA',0,'float32'), \
     ('I_WELL_SPACED_STEMS_PER_HA',0,'float32'), \
     ('I_FREE_GROWING_STEMS_PER_HA',0,'float32'), \
     ('I_CROWN_CLOSURE_PERCENT',0,'float32'), \
     ('REFERENCE_YEAR',0,'float32'), \
     ('STOCKING_STATUS_CODE',1,'int16'), \
     ('STOCKING_TYPE_CODE',1,'int16'), \
     ('FOREST_COVER_WHO_CREATED',1,'int16'), \
     ('FOREST_COVER_WHO_UPDATED',1,'int16'), \
     ('FOREST_COVER_WHEN_CREATED',0,'<U20'), \
     ('FOREST_COVER_WHEN_UPDATED',0,'<U20')]
    d['LUT']={}
    for field in d['Field List']: d['LUT'][field[0]]=[] 
    LayerInfo.append(d)

    # Add layer and define required variables
    # *** Notes:OBSERVATION_DATE is empty ***
    d={}
    d['Layer Name']='RSLT_FOREST_COVER_SILV_SVW'
    d['Path']=PathInResultsFull
    d['File Name']='Results.gdb'
    d['Field List']=[('OPENING_ID',0,'float32'), \
     ('FOREST_COVER_ID',0,'float32'), \
     ('SILV_POLYGON_NUMBER',1,'int32'), \
     ('SILV_POLYGON_AREA',0,'float32'), \
     ('STOCKING_STANDARD_UNIT_ID',0,'float32'), \
     ('SITE_INDEX',0,'float32'), \
     ('S_SPECIES_AGE_1',0,'float32'), \
     ('S_SPECIES_CODE_1',1,'int16'), \
     ('S_SPECIES_CODE_2',1,'int16'), \
     ('S_SPECIES_CODE_3',1,'int16'), \
     ('S_SPECIES_CODE_4',1,'int16'), \
     ('S_SPECIES_CODE_5',1,'int16'), \
     ('S_SPECIES_PERCENT_1',0,'int16'), \
     ('S_SPECIES_PERCENT_2',0,'int16'), \
     ('S_SPECIES_PERCENT_3',0,'int16'), \
     ('S_SPECIES_PERCENT_4',0,'int16'), \
     ('S_SPECIES_PERCENT_5',0,'int16'), \
     ('S_TOTAL_STEMS_PER_HA',0,'float32'), \
     ('S_TOTAL_WELL_SPACED_STEMS_HA',0,'float32'), \
     ('S_WELL_SPACED_STEMS_PER_HA',0,'float32'), \
     ('S_FREE_GROWING_STEMS_PER_HA',0,'float32'), \
     ('S_CROWN_CLOSURE_PERCENT',0,'float32'), \
     ('REFERENCE_YEAR',0,'float32'), \
     ('STOCKING_STATUS_CODE',1,'int16'), \
     ('STOCKING_TYPE_CODE',1,'int16'), \
     ('FOREST_COVER_WHO_CREATED',1,'int16'), \
     ('FOREST_COVER_WHO_UPDATED',1,'int16'), \
     ('FOREST_COVER_WHEN_CREATED',0,'<U20'), \
     ('FOREST_COVER_WHEN_UPDATED',0,'<U20')]
    d['LUT']={}
    for field in d['Field List']: d['LUT'][field[0]]=[] 
    LayerInfo.append(d)

    # Add layer and define required variables
    # *** Notes:OBSERVATION_DATE is empty ***
    d={}
    d['Layer Name']='RSLT_FOREST_COVER_RESERVE_SVW'
    d['Path']=PathInResultsFull
    d['File Name']='Results.gdb'
    d['Field List']=[('OPENING_ID',0,'float32'), \
     ('SILV_RESERVE_CODE',1,'int16'), \
     ('SILV_RESERVE_OBJECTIVE_CODE',1,'int16')]
    d['LUT']={}
    for field in d['Field List']: d['LUT'][field[0]]=[] 
    LayerInfo.append(d)

    # Add layer and define required variables
    d={}
    d['Layer Name']='RSLT_PLANTING_SVW'
    d['Path']=PathInResultsFull
    d['File Name']='Results.gdb'
    d['Field List']=[('OPENING_ID',0,'float32'), \
     ('ACTIVITY_TREATMENT_UNIT_ID',0,'float32'), \
     ('NUMBER_PLANTED',0,'float32'), \
     ('FIA_PROJECT_ID',1,'int32'), \
     ('SILV_TREE_SPECIES_CODE',1,'int16'), \
     ('SEEDLOT_NUMBER',0,'float32'), \
     ('VEG_LOT_ID',1,'int32'), \
     ('ACTUAL_TREATMENT_AREA',0,'float32'), \
     ('ATU_COMPLETION_DATE',0,'<U20')]
    d['LUT']={}
    for field in d['Field List']: d['LUT'][field[0]]=[]  
    LayerInfo.append(d)

    # Add layer and define required variables
    d={}
    d['Layer Name']='PEST_INFESTATION_POLY'
    d['Path']=PathInDisturbancesFull
    d['File Name']='Disturbances.gdb'
    d['Field List']=[('PEST_SEVERITY_CODE',1,'int16'), \
     ('PEST_SPECIES_CODE',1,'int16'), \
     ('CAPTURE_YEAR',0,'int16')]
    d['LUT']={}
    for field in d['Field List']: d['LUT'][field[0]]=[] 
    LayerInfo.append(d)

    # Add layer and define required variables
    d={}
    d['Layer Name']='PROT_HISTORICAL_FIRE_POLYS_SP'
    d['Path']=PathInDisturbancesFull
    d['File Name']='Disturbances.gdb'
    d['Field List']=[
     ('FIRE_YEAR',0,'int16'), \
     ('FIRE_CAUSE',1,'int16'), \
     ('FIRE_DATE',0,'<U20'), \
     ('FIRE_NUMBER',0,'<U20')]
    d['LUT']={}
    for field in d['Field List']: d['LUT'][field[0]]=[] 
    LayerInfo.append(d)

    # Add layer and define required variables
    d={}
    d['Layer Name']='VEG_COMP_LYR_R1_POLY'
    d['Path']=PathInVRIFull
    d['File Name']='VRI.gdb'
    d['Field List']=[('BCLCS_LEVEL_1',1,'int16'), \
     ('BCLCS_LEVEL_2',1,'int16'), \
     ('BCLCS_LEVEL_3',1,'int16'), \
     ('BCLCS_LEVEL_4',1,'int16'), \
     ('BCLCS_LEVEL_5',1,'int16'), \
     ('REFERENCE_YEAR',0,'int16'), \
     ('PROJECTED_DATE',0,'<U20'), \
     ('LAND_COVER_CLASS_CD_1',1,'int16'), \
     ('EST_COVERAGE_PCT_1',0,'float32'), \
     ('OPENING_ID',0,'float32'), \
     ('BEC_ZONE_CODE',1,'int16'), \
     ('BEC_SUBZONE',1,'int16'), \
     ('BEC_VARIANT',1,'int16'), \
     ('EARLIEST_NONLOGGING_DIST_TYPE',1,'int32'), \
     ('EARLIEST_NONLOGGING_DIST_DATE',1,'int32'), \
     ('STAND_PERCENTAGE_DEAD',0,'float32'), \
     ('FREE_TO_GROW_IND',1,'int16'), \
     ('EST_SITE_INDEX',0,'float32'), \
     ('SITE_INDEX',0,'float32'), \
     ('VRI_LIVE_STEMS_PER_HA',0,'float32'), \
     ('VRI_DEAD_STEMS_PER_HA',0,'float32'), \
     ('SPECIES_CD_1',1,'int16'),('SPECIES_PCT_1',0,'int16'), \
     ('SPECIES_CD_2',1,'int16'),('SPECIES_PCT_2',0,'int16'), \
     ('SPECIES_CD_3',1,'int16'),('SPECIES_PCT_3',0,'int16'), \
     ('SPECIES_CD_4',1,'int16'),('SPECIES_PCT_4',0,'int16'), \
     ('SPECIES_CD_5',1,'int16'),('SPECIES_PCT_5',0,'int16'), \
     ('SPECIES_CD_6',1,'int16'),('SPECIES_PCT_6',0,'int16'), \
     ('PROJ_AGE_1',0,'float32'), \
     ('PROJ_AGE_2',0,'float32'), \
     ('PROJ_HEIGHT_1',0,'float32'), \
     ('PROJ_HEIGHT_2',0,'float32'), \
     ('LINE_6_SITE_PREP_HISTORY',1,'int32'), \
     ('LINE_7B_DISTURBANCE_HISTORY',1,'int32'), \
     ('LINE_8_PLANTING_HISTORY',1,'int32'), \
     ('LIVE_VOL_PER_HA_SPP1_125',0,'float32'), \
     ('LIVE_VOL_PER_HA_SPP2_125',0,'float32'), \
     ('LIVE_VOL_PER_HA_SPP3_125',0,'float32'), \
     ('LIVE_VOL_PER_HA_SPP4_125',0,'float32'), \
     ('LIVE_VOL_PER_HA_SPP5_125',0,'float32'), \
     ('LIVE_VOL_PER_HA_SPP6_125',0,'float32'), \
     ('LIVE_STAND_VOLUME_125',0,'float32'), \
     ('DEAD_STAND_VOLUME_125',0,'float32'), \
     ('WHOLE_STEM_BIOMASS_PER_HA',0,'float32'), \
     ('BRANCH_BIOMASS_PER_HA',0,'float32'), \
     ('FOLIAGE_BIOMASS_PER_HA',0,'float32'), \
     ('BARK_BIOMASS_PER_HA',0,'float32')]
    d['LUT']={}
    for field in d['Field List']: d['LUT'][field[0]]=[] 
    LayerInfo.append(d)

    # Add layer and define required variables
    d={}
    d['Layer Name']='VEG_CONSOLIDATED_CUT_BLOCKS_SP'
    d['Path']=PathInDisturbancesFull
    d['File Name']='Disturbances.gdb'
    d['Field List']=[('OPENING_ID',0,'float32'), \
     ('HARVEST_YEAR',0,'int16'), \
     ('AREA_HA',0,'float32'), \
     ('DISTURBANCE_START_DATE',0,'<U20'), \
     ('DISTURBANCE_END_DATE',0,'<U20')]
    d['LUT']={}
    for field in d['Field List']: d['LUT'][field[0]]=[] 
    LayerInfo.append(d)

    # Add layer and define required variables
    d={}
    d['Layer Name']='VEG_BURN_SEVERITY_SP'
    d['Path']=PathInDisturbancesFull
    d['File Name']='Disturbances.gdb'
    d['Field List']=[('FIRE_YEAR',0,'int16'), \
     ('BURN_SEVERITY_RATING',1,'int16')]
    d['LUT']={}
    for field in d['Field List']: d['LUT'][field[0]]=[] 
    LayerInfo.append(d)
    
    # Add layer and define required variables
    d={}
    d['Layer Name']='RMP_PLAN_LEGAL_POLY_SVW'
    d['Path']=PathInLUPFull
    d['File Name']='LandUse.gdb'
    d['Field List']=[('STRGC_LAND_RSRCE_PLAN_NAME',1,'int16'), \
             ('LEGAL_FEAT_OBJECTIVE',1,'int16')]
    d['LUT']={}
    for field in d['Field List']: d['LUT'][field[0]]=[] 
    LayerInfo.append(d)
    
    # Add layer and define required variables
    d={}
    d['Layer Name']='RMP_PLAN_NON_LEGAL_POLY_SVW'
    d['Path']=PathInLUPFull
    d['File Name']='LandUse.gdb'
    d['Field List']=[('STRGC_LAND_RSRCE_PLAN_NAME',1,'int16'), \
             ('NON_LEGAL_FEAT_OBJECTIVE',1,'int16')]
    d['LUT']={}
    for field in d['Field List']: d['LUT'][field[0]]=[] 
    LayerInfo.append(d)
    
    # Add layer and define required variables
    d={}
    d['Layer Name']='RMP_OGMA_LEGAL_CURRENT_SVW'
    d['Path']=PathInLUPFull
    d['File Name']='LandUse.gdb'
    d['Field List']=[('LEGAL_OGMA_PROVID',1,'int16')]
    d['LUT']={}
    for field in d['Field List']: d['LUT'][field[0]]=[] 
    LayerInfo.append(d)
    
    # Add layer and define required variables
    d={}
    d['Layer Name']='WCP_UNGULATE_WINTER_RANGE_SP'
    d['Path']=PathInLUPFull
    d['File Name']='LandUse.gdb'
    d['Field List']=[('UWR_NUMBER',1,'int16')]
    d['LUT']={}
    for field in d['Field List']: d['LUT'][field[0]]=[] 
    LayerInfo.append(d)
    
    # Add layer and define required variables
    d={}
    d['Layer Name']='TA_PARK_ECORES_PA_SVW'
    d['Path']=PathInLUPFull
    d['File Name']='LandUse.gdb'
    d['Field List']=[('PROTECTED_LANDS_NAME',1,'int32'), \
         ('PROTECTED_LANDS_CODE',1,'int32'), \
         ('PROTECTED_LANDS_DESIGNATION',1,'int32'), \
         ('PARK_CLASS',1,'int32')]
    d['LUT']={}
    for field in d['Field List']: d['LUT'][field[0]]=[] 
    LayerInfo.append(d)
    
    # Add layer and define required variables
    #d={}
    #d['Layer Name']='WCL_CONSERVATION_LANDS_SP'
    #d['Path']=PathInLUPFull
    #d['File Name']='LandUse.gdb'
    #d['Field List']=[]
    #d['LUT']={}
    #for field in d['Field List']: d['LUT'][field[0]]=[] 
    #LayerInfo.append(d)
    
    #fiona.listlayers(PathInLUPFull + '\\Landuse.gdb')
    #lyr=fiona.open(PathInLUPFull + '\\Landuse.gdb',layer='TA_PARK_ECORES_PA_SVW')
    #lyr.schema
    
    return LayerInfo

#%% BUILD LOOK-UP TABLES FOR THE CODES IN INVENTORY LAYERS

def BuildForestInventoryLUTs(LayerInfo):

    t0=time.time()
    for iLyr in range(len(LayerInfo)):

        # Start counting time
        t_start=time.time()
        
        cnt=0
        
        # Loop through features in layer
        with fiona.open(LayerInfo[iLyr]['Path'] + '\\' + LayerInfo[iLyr]['File Name'],layer=LayerInfo[iLyr]['Layer Name']) as source:    
            
            for feat in source:
   
                for fnam,flag,dtype in LayerInfo[iLyr]['Field List']: 
                
                    # Only continue if it is a string variable
                    if flag==0: 
                        continue
                
                    # Only continue if the variable is populated
                    if feat['properties'][fnam]==None: 
                        continue
                    
                    # Add to code list
                    LayerInfo[iLyr]['LUT'][fnam].append(feat['properties'][fnam])
                
        # Check time
        t_ela=time.time()-t_start
        print(t_ela)
            
        # Remove duplicates
        for fnam,flag,dtype in LayerInfo[iLyr]['Field List']:
            if flag!=0:
                LayerInfo[iLyr]['LUT'][fnam]=np.unique(LayerInfo[iLyr]['LUT'][fnam]).tolist()
    
        # Build look-up table
        lut={}
        for fnam,flag,dtype in LayerInfo[iLyr]['Field List']:
            L=len(LayerInfo[iLyr]['LUT'][fnam])
            if L!=0:
                lut[fnam]={}
                for i in range(L):
                    lut[fnam][LayerInfo[iLyr]['LUT'][fnam][i]]=np.array(i+1,dtype=dtype,ndmin=1)
        
        # Save
        gu.opickle(LayerInfo[iLyr]['Path'] + '\\LUTs_' + LayerInfo[iLyr]['Layer Name'] +'.pkl',lut)

    #--------------------------------------------------------------------------
    # Standardize codes for species
    #--------------------------------------------------------------------------
    
    # Import data
    for iLyr in range(len(LayerInfo)):        
        if LayerInfo[iLyr]['Layer Name']=='VEG_COMP_LYR_R1_POLY':
            lut_vri=gu.ipickle(LayerInfo[iLyr]['Path'] + '\\LUTs_VEG_COMP_LYR_R1_POLY.pkl')
        if LayerInfo[iLyr]['Layer Name']=='RSLT_FOREST_COVER_INV_SVW':
            lut_fci=gu.ipickle(LayerInfo[iLyr]['Path'] + '\\LUTs_RSLT_FOREST_COVER_INV_SVW.pkl')    
        if LayerInfo[iLyr]['Layer Name']=='RSLT_FOREST_COVER_SILV_SVW':
            lut_fcs=gu.ipickle(LayerInfo[iLyr]['Path'] + '\\LUTs_RSLT_FOREST_COVER_SILV_SVW.pkl')        
        if LayerInfo[iLyr]['Layer Name']=='RSLT_PLANTING_SVW':
            lut_pl=gu.ipickle(LayerInfo[iLyr]['Path'] + '\\LUTs_RSLT_PLANTING_SVW.pkl')  
            
    # Add everything to a list
    cd=list(lut_vri['SPECIES_CD_1'].keys())
    cd=cd+list(lut_vri['SPECIES_CD_2'].keys())
    cd=cd+list(lut_vri['SPECIES_CD_3'].keys())
    cd=cd+list(lut_vri['SPECIES_CD_4'].keys())
    cd=cd+list(lut_vri['SPECIES_CD_5'].keys())
    cd=cd+list(lut_vri['SPECIES_CD_6'].keys())
    cd=cd+list(lut_fci['I_SPECIES_CODE_1'].keys())
    cd=cd+list(lut_fci['I_SPECIES_CODE_2'].keys())
    cd=cd+list(lut_fci['I_SPECIES_CODE_3'].keys())
    cd=cd+list(lut_fci['I_SPECIES_CODE_4'].keys())
    cd=cd+list(lut_fci['I_SPECIES_CODE_5'].keys())            
    cd=cd+list(lut_fcs['S_SPECIES_CODE_1'].keys())
    cd=cd+list(lut_fcs['S_SPECIES_CODE_2'].keys())
    cd=cd+list(lut_fcs['S_SPECIES_CODE_3'].keys())
    cd=cd+list(lut_fcs['S_SPECIES_CODE_4'].keys())
    cd=cd+list(lut_fcs['S_SPECIES_CODE_5'].keys())            
    cd=cd+list(lut_pl['SILV_TREE_SPECIES_CODE'].keys())
            
    # Get unique list
    uCode=np.unique(cd)
    
    # Save new standardized list
    for iLyr in range(len(LayerInfo)):        
        
        if LayerInfo[iLyr]['Layer Name']=='VEG_COMP_LYR_R1_POLY':
            lut_vri['SPECIES_CD_1']={}
            lut_vri['SPECIES_CD_2']={}
            lut_vri['SPECIES_CD_3']={}
            lut_vri['SPECIES_CD_4']={}
            lut_vri['SPECIES_CD_5']={}
            lut_vri['SPECIES_CD_6']={}
            for i in range(len(uCode)):
                lut_vri['SPECIES_CD_1'][uCode[i]]=i+1
                lut_vri['SPECIES_CD_2'][uCode[i]]=i+1
                lut_vri['SPECIES_CD_3'][uCode[i]]=i+1
                lut_vri['SPECIES_CD_4'][uCode[i]]=i+1
                lut_vri['SPECIES_CD_5'][uCode[i]]=i+1
                lut_vri['SPECIES_CD_6'][uCode[i]]=i+1
            gu.opickle(LayerInfo[iLyr]['Path'] + '\\LUTs_VEG_COMP_LYR_R1_POLY.pkl',lut_vri)
        
        if LayerInfo[iLyr]['Layer Name']=='RSLT_FOREST_COVER_INV_SVW':
            lut_fci['I_SPECIES_CODE_1']={}
            lut_fci['I_SPECIES_CODE_2']={}
            lut_fci['I_SPECIES_CODE_3']={}
            lut_fci['I_SPECIES_CODE_4']={}
            lut_fci['I_SPECIES_CODE_5']={}            
            for i in range(len(uCode)):
                lut_fci['I_SPECIES_CODE_1'][uCode[i]]=i+1
                lut_fci['I_SPECIES_CODE_2'][uCode[i]]=i+1
                lut_fci['I_SPECIES_CODE_3'][uCode[i]]=i+1
                lut_fci['I_SPECIES_CODE_4'][uCode[i]]=i+1
                lut_fci['I_SPECIES_CODE_5'][uCode[i]]=i+1
            gu.opickle(LayerInfo[iLyr]['Path'] + '\\LUTs_RSLT_FOREST_COVER_INV_SVW.pkl',lut_fci)
        
        if LayerInfo[iLyr]['Layer Name']=='RSLT_FOREST_COVER_SILV_SVW':
            lut_fcs['S_SPECIES_CODE_1']={}
            lut_fcs['S_SPECIES_CODE_2']={}
            lut_fcs['S_SPECIES_CODE_3']={}
            lut_fcs['S_SPECIES_CODE_4']={}
            lut_fcs['S_SPECIES_CODE_5']={}            
            for i in range(len(uCode)):
                lut_fcs['S_SPECIES_CODE_1'][uCode[i]]=i+1
                lut_fcs['S_SPECIES_CODE_2'][uCode[i]]=i+1
                lut_fcs['S_SPECIES_CODE_3'][uCode[i]]=i+1
                lut_fcs['S_SPECIES_CODE_4'][uCode[i]]=i+1
                lut_fcs['S_SPECIES_CODE_5'][uCode[i]]=i+1
            gu.opickle(LayerInfo[iLyr]['Path'] + '\\LUTs_RSLT_FOREST_COVER_SILV_SVW.pkl',lut_fcs)
        
        if LayerInfo[iLyr]['Layer Name']=='RSLT_PLANTING_SVW':
            lut_pl['SILV_TREE_SPECIES_CODE']={}
            for i in range(len(uCode)):
                lut_pl['SILV_TREE_SPECIES_CODE'][uCode[i]]=i+1
            gu.opickle(LayerInfo[iLyr]['Path'] + '\\LUTs_RSLT_PLANTING_SVW.pkl',lut_pl)
    print((time.time()-t0)/60)

#%% EXTRACT NUMERIC TIME VECTORS FROM DATE STRINGS IN RESULTS

def ExtractDateStringsFromRESULTS(lyr_nam,data):
    
    if (lyr_nam=='RSLT_ACTIVITY_TREATMENT_SVW') | (lyr_nam=='RSLT_PLANTING_SVW'):
        
        ds=data['ATU_COMPLETION_DATE']
        L=ds.size
        data['Year']=-999*np.ones(L,dtype='int16')
        data['Month']=-999*np.ones(L,dtype='int16')
        data['Day']=-999*np.ones(L,dtype='int16')
        for i in range(L):
            if ds[i]!='':
                data['Year'][i]=int(ds[i][0:4])
                data['Month'][i]=int(ds[i][5:7])
                data['Day'][i]=int(ds[i][8:10])
        
        # Remove original date string
        del data['ATU_COMPLETION_DATE']
            
    elif (lyr_nam=='RSLT_OPENING_SVW'):
        
        L=data['OPENING_ID'].size
        ds=data['DISTURBANCE_START_DATE']                    
        data['Year_DistStart']=-999*np.ones(L,dtype='int16')
        data['Month_DistStart']=-999*np.ones(L,dtype='int16')
        data['Day_DistStart']=-999*np.ones(L,dtype='int16')
        for i in range(L):            
            if ds[i]!='':
                data['Year_DistStart'][i]=int(ds[i][0:4])
                data['Month_DistStart'][i]=int(ds[i][5:7])
                data['Day_DistStart'][i]=int(ds[i][8:10])
        
        ds=data['DISTURBANCE_END_DATE']
        data['Year_DistEnd']=-999*np.ones(L,dtype='int16')
        data['Month_DistEnd']=-999*np.ones(L,dtype='int16')
        data['Day_DistEnd']=-999*np.ones(L,dtype='int16')
        for i in range(L):            
            if ds[i]!='':
                data['Year_DistEnd'][i]=int(ds[i][0:4])
                data['Month_DistEnd'][i]=int(ds[i][5:7])
                data['Day_DistEnd'][i]=int(ds[i][8:10])
        
        ds=data['DENUDATION_1_COMPLETION_DATE']
        data['Year_Denu1_Comp']=-999*np.ones(L,dtype='int16')
        data['Month_Denu1_Comp']=-999*np.ones(L,dtype='int16')
        data['Day_Denu1_Comp']=-999*np.ones(L,dtype='int16')
        for i in range(L):            
            if ds[i]!='':
                data['Year_Denu1_Comp'][i]=int(ds[i][0:4])
                data['Month_Denu1_Comp'][i]=int(ds[i][5:7])
                data['Day_Denu1_Comp'][i]=int(ds[i][8:10])
        
        ds=data['DENUDATION_2_COMPLETION_DATE']
        data['Year_Denu2_Comp']=-999*np.ones(L,dtype='int16')
        data['Month_Denu2_Comp']=-999*np.ones(L,dtype='int16')
        data['Day_Denu2_Comp']=-999*np.ones(L,dtype='int16')
        for i in range(L):            
            if ds[i]!='':
                data['Year_Denu2_Comp'][i]=int(ds[i][0:4])
                data['Month_Denu2_Comp'][i]=int(ds[i][5:7])
                data['Day_Denu2_Comp'][i]=int(ds[i][8:10])
            
        # Remove original date string
        del data['DISTURBANCE_START_DATE']       
        del data['DISTURBANCE_END_DATE']
        del data['DENUDATION_1_COMPLETION_DATE']
        del data['DENUDATION_2_COMPLETION_DATE']
            
    elif (lyr_nam=='RSLT_FOREST_COVER_INV_SVW') | (lyr_nam=='RSLT_FOREST_COVER_SILV_SVW'):
        
        L=data['OPENING_ID'].size        
        ds=data['FOREST_COVER_WHEN_CREATED']
        data['Year_Created']=-999*np.ones(L,dtype='int16')
        data['Month_Created']=-999*np.ones(L,dtype='int16')
        data['Day_Created']=-999*np.ones(L,dtype='int16')
        for i in range(L):     
            if ds[i]!='':
                data['Year_Created'][i]=int(ds[i][0:4])
                data['Month_Created'][i]=int(ds[i][5:7])
                data['Day_Created'][i]=int(ds[i][8:10])
            
        ds=data['FOREST_COVER_WHEN_UPDATED']
        data['Year_Updated']=-999*np.ones(L,dtype='int16')
        data['Month_Updated']=-999*np.ones(L,dtype='int16')
        data['Day_Updated']=-999*np.ones(L,dtype='int16')
        for i in range(L):            
            if ds[i]!='':
                data['Year_Updated'][i]=int(ds[i][0:4])
                data['Month_Updated'][i]=int(ds[i][5:7])
                data['Day_Updated'][i]=int(ds[i][8:10])
            
        # Remove original date string
        del data['FOREST_COVER_WHEN_CREATED']     
        del data['FOREST_COVER_WHEN_UPDATED']
        
    elif (lyr_nam=='PROT_HISTORICAL_FIRE_POLYS_SP'):
        
        L=data['FIRE_YEAR'].size
        ds=data['FIRE_DATE']
        data['Year']=-999*np.ones(L,dtype='int16')
        data['Month']=-999*np.ones(L,dtype='int16')
        data['Day']=-999*np.ones(L,dtype='int16')
        for i in range(L):
            if ds[i]!='':
                data['Year'][i]=int(ds[i][0:4])
                data['Month'][i]=int(ds[i][5:7])
                data['Day'][i]=int(ds[i][8:10])
        
        # Remove original date string
        del data['FIRE_DATE']
    
    elif (lyr_nam=='VEG_CONSOLIDATED_CUT_BLOCKS_SP'):
        
        L=data['HARVEST_YEAR'].size
        ds=data['DISTURBANCE_START_DATE']
        
        data['Year_Start']=-999*np.ones(L,dtype='int16')
        data['Month_Start']=-999*np.ones(L,dtype='int16')
        data['Day_Start']=-999*np.ones(L,dtype='int16')        
        for i in range(L):            
            if ds[i]!='':
                data['Year_Start'][i]=int(ds[i][0:4])
                data['Month_Start'][i]=int(ds[i][5:7])
                data['Day_Start'][i]=int(ds[i][8:10])
        
        ds=data['DISTURBANCE_END_DATE']
        data['Year_End']=-999*np.ones(L,dtype='int16')
        data['Month_End']=-999*np.ones(L,dtype='int16')
        data['Day_End']=-999*np.ones(L,dtype='int16')        
        for i in range(L):
            if ds[i]!='':
                data['Year_End'][i]=int(ds[i][0:4])
                data['Month_End'][i]=int(ds[i][5:7])
                data['Day_End'][i]=int(ds[i][8:10])    
        
        # Remove original date string
        del data['DISTURBANCE_START_DATE']
        del data['DISTURBANCE_END_DATE']
    
    elif (lyr_nam=='VEG_COMP_LYR_R1_POLY'):
        
        L=data['PROJECTED_DATE'].size
        ds=data['PROJECTED_DATE']
        data['Year']=-999*np.ones(L,dtype='int16')
        data['Month']=-999*np.ones(L,dtype='int16')
        data['Day']=-999*np.ones(L,dtype='int16')
        for i in range(L):            
            if ds[i]!='':
                data['Year'][i]=int(ds[i][0:4])
                data['Month'][i]=int(ds[i][5:7])
                data['Day'][i]=int(ds[i][8:10])
        
        # Remove original date string
        del data['PROJECTED_DATE']
            
    return data


#%% Recover missing ATU layer geometries

def RecoverMissingATUGeometries(meta):
    
    #--------------------------------------------------------------------------
    # Query AT layer for features with missing spatial information
    # Only consider planting, direct seeding, site prep and pest control
    # 2 min
    #--------------------------------------------------------------------------

    t0=time.time()
    
    lyr=fiona.open(PathInResultsFull + '\\Results.gdb',layer='RSLT_ACTIVITY_TREATMENT_SVW')
    L=len(lyr)
    atu_full={}
    atu_full['OPENING_ID']=np.zeros(L)
    atu_full['SpatialMissing']=np.zeros(L)
    atu_full['Flag Included']=np.zeros(L)
    atu_full['Year']=np.zeros(L)
    atu_full['ACTUAL_TREATMENT_AREA']=np.zeros(L)
    atu_full['ACTIVITY_TREATMENT_UNIT_ID']=np.zeros(L)    
    atu_full['SBC']=np.array(['empty' for _ in range(L)],dtype=object)
    atu_full['SMC']=np.array(['empty' for _ in range(L)],dtype=object)
    atu_full['STC']=np.array(['empty' for _ in range(L)],dtype=object)
    cnt=0
    with fiona.open(PathInResultsFull + '\\Results.gdb',layer='RSLT_ACTIVITY_TREATMENT_SVW') as source:
        for feat in source:
            
            prp=feat['properties']
            
            # Opening ID
            atu_full['OPENING_ID'][cnt]=prp['OPENING_ID']
            atu_full['ACTIVITY_TREATMENT_UNIT_ID'][cnt]=prp['ACTIVITY_TREATMENT_UNIT_ID']
            
            if (prp['SILV_BASE_CODE']=='PL') & (prp['SILV_TECHNIQUE_CODE']!='SE') & (prp['SILV_TECHNIQUE_CODE']!='CG') & (prp['SILV_METHOD_CODE']!='LAYOT') & (prp['RESULTS_IND']=='Y') | \
                (prp['SILV_BASE_CODE']=='DS') & (prp['RESULTS_IND']=='Y') | \
                (prp['SILV_BASE_CODE']=='SP') & (prp['RESULTS_IND']=='Y') | \
                (prp['SILV_BASE_CODE']=='PC') & (prp['RESULTS_IND']=='Y') | \
                (prp['SILV_BASE_CODE']=='FE') & (prp['SILV_TECHNIQUE_CODE']=='CA') & (prp['SILV_METHOD_CODE']=='HELI') & (prp['RESULTS_IND']=='Y') & (prp['ACTUAL_TREATMENT_AREA']!=None) & (prp['ATU_COMPLETION_DATE']!=None) & (prp['SILV_FUND_SOURCE_CODE']!=None):
                atu_full['Flag Included'][cnt]=1
            
            # Is spatial missing
            if feat['geometry']==None:
                atu_full['SpatialMissing'][cnt]=1      
            
            # SBC
            atu_full['SBC'][cnt]=prp['SILV_BASE_CODE']
            atu_full['STC'][cnt]=prp['SILV_TECHNIQUE_CODE']
            atu_full['SMC'][cnt]=prp['SILV_METHOD_CODE']
            
            # Area
            if prp['ACTUAL_TREATMENT_AREA']!=None:
                atu_full['ACTUAL_TREATMENT_AREA'][cnt]=prp['ACTUAL_TREATMENT_AREA']
            
            # Add year
            Year=0
            if prp['ATU_COMPLETION_DATE']!=None:
                Year=int(prp['ATU_COMPLETION_DATE'][0:4])
            atu_full['Year'][cnt]=Year
            
            # Update counter
            cnt=cnt+1
    
    #--------------------------------------------------------------------------
    # Extract missing AT geometries
    #--------------------------------------------------------------------------
    
    # Index to missing and included entries
    iMisAT=np.where( (atu_full['SpatialMissing']==1) & (atu_full['Flag Included']==1) )[0]
    
    atu_mis=atu_full.copy()
    for k in atu_mis.keys():
        atu_mis[k]=atu_mis[k][iMisAT]
    
    # Save
    gu.opickle(PathInResultsFull + '\\missing_geo_atu_list.pkl',atu_mis)
    
    print((time.time()-t0)/60)
    
    #--------------------------------------------------------------------------
    # Get missing AT geometries from opening layer
    # 3 min
    #--------------------------------------------------------------------------
    
    t0=time.time()
    
    u=np.unique(atu_mis['OPENING_ID'])    
    op_mis={}
    for i in range(u.size):
        op_mis[int(u[i])]=[]
        
    with fiona.open(PathInResultsFull + '\\Results.gdb',layer='RSLT_OPENING_SVW') as source:
        for feat in source:
            
            if (feat['geometry']==None):
                continue
            
            id=int(feat['properties']['OPENING_ID'])
            
            # Index to AT entries that correspond to this feature and have missing geometries
            indMis=np.where( (atu_mis['OPENING_ID']==id) )[0]
            if indMis.size==0:
                continue
            
            if id in op_mis:
                op_mis[id].append(feat['geometry'])
    
    gu.opickle(PathInResultsFull + '\\missing_geo_op_geos.pkl',op_mis)
    
    print((time.time()-t0)/60)
    
    #--------------------------------------------------------------------------
    # Forest cover by OPENING_ID (only saving geometries for artificial status)
    # Takes 5 min
    #--------------------------------------------------------------------------
    
    t0=time.time()
    
    u=np.unique(atu_mis['OPENING_ID'])    
    fc_mis={}
    for i in range(u.size):
        fc_mis[u[i]]=[]
        
    with fiona.open(PathInResultsFull + '\\Results.gdb',layer='RSLT_FOREST_COVER_INV_SVW') as source:
        for feat in source:
            
            if (feat['geometry']==None):
                continue
            
            if (feat['properties']['STOCKING_TYPE_CODE']!='ART'):
                continue
            
            id=feat['properties']['OPENING_ID']
            
            # Index to AT entries that correspond to this feature and have missing geometries
            indMis=np.where( (atu_mis['OPENING_ID']==id) )[0]
            if indMis.size==0:
                continue
            
            if id in fc_mis:
                fc_mis[id].append(feat['geometry'])
    
    gu.opickle(PathInResultsFull + '\\missing_geo_fc_geos.pkl',fc_mis)
    
    print((time.time()-t0)/60)    

    
    #--------------------------------------------------------------------------
    # Find matching FC geometries for each AT entry with missing spatial
    # Takes 5 min
    #--------------------------------------------------------------------------
    
    atu_mis=gu.ipickle(PathInResultsFull + '\\missing_geo_atu_list.pkl')
    fc_mis=gu.ipickle(PathInResultsFull + '\\missing_geo_fc_geos.pkl')
    
    t0=time.time()
    
    atu_mis['IdxToFC']=[None]*atu_mis['OPENING_ID'].size
    FlgAgreement=np.zeros(atu_mis['OPENING_ID'].size)
    IsDone=np.zeros(atu_mis['OPENING_ID'].size)
    
    for i in range(atu_mis['OPENING_ID'].size):
    
        if IsDone[i]==1:
            continue
        if len(fc_mis[atu_mis['OPENING_ID'][i]])==0:
            continue
    
        ind_atu=np.where( (atu_mis['OPENING_ID']==atu_mis['OPENING_ID'][i]) )[0]
            
        A_at=atu_mis['ACTUAL_TREATMENT_AREA'][ind_atu]  
    
        L=len(fc_mis[atu_mis['OPENING_ID'][i]])
        A_fc=np.zeros(L)
        for j in range(L):
            feat_fc={}; feat_fc['properties']={}; feat_fc['geometry']=fc_mis[atu_mis['OPENING_ID'][i]][j]
            gdf_fc=gpd.GeoDataFrame.from_features([feat_fc])
            A_fc[j]=np.round(gdf_fc.area.values[0]/10000)
        
        # Find best matches
        dist=np.abs(A_at[:, np.newaxis]-A_fc)
        Closest=dist.argmin(axis=1)
        A_fc_BestMatch=A_fc[Closest]
        E=A_at-A_fc_BestMatch

        iMatch=np.where( (np.abs(E)<=2) )[0]
        if iMatch.size>0:
            for j in range(iMatch.size):
                atu_mis['IdxToFC'][ind_atu[iMatch[j]]]=[ Closest[iMatch[j]] ]
            FlgAgreement[ind_atu[iMatch]]=1
        
        iNoMatch=np.where( (np.abs(E)>2) )[0]
        if iNoMatch.size>0:    
            idx=np.arange(0,L)
            idx=idx[~np.isin(idx,Closest[iMatch])]
            for j in range(iNoMatch.size):            
                atu_mis['IdxToFC'][ind_atu[iNoMatch[j]]]=list(idx)

        # Update progress tracker
        IsDone[ind_atu]=1
    
        #print(i)

    print((time.time()-t0)/60)    

    # Success rate
    ind=np.where(FlgAgreement==1)[0]
    ind.size/atu_mis['OPENING_ID'].size*100

    # Eliminate none values
    for i in range(len(atu_mis['IdxToFC'])):
        if atu_mis['IdxToFC'][i]==None:
            atu_mis['IdxToFC'][i]=[]
            
    gu.opickle(PathInResultsFull + '\\missing_geo_atu_list.pkl',atu_mis)

#%% Recover missing FC layer geometries

def RecoverMissingFCGeometries(meta):
    
    #--------------------------------------------------------------------------
    # Query AT layer for features with missing spatial information
    # Only consider planting, direct seeding, site prep and pest control
    # 14 min
    #--------------------------------------------------------------------------

    t0=time.time()
    
    lyr=fiona.open(PathInResultsFull + '\\Results.gdb',layer='RSLT_FOREST_COVER_INV_SVW')
    L=len(lyr)
    fcinv_full={}
    fcinv_full['OPENING_ID']=np.zeros(L)
    fcinv_full['SpatialMissing']=np.zeros(L)
    fcinv_full['SILV_POLYGON_NUMBER']=np.array(['' for _ in range(L)],dtype=object)
    fcinv_full['SILV_POLYGON_NET_AREA']=np.zeros(L)
    cnt=0
    with fiona.open(PathInResultsFull + '\\Results.gdb',layer='RSLT_FOREST_COVER_INV_SVW') as source:
        for feat in source:            
            prp=feat['properties']            
            fcinv_full['OPENING_ID'][cnt]=prp['OPENING_ID']
            fcinv_full['SILV_POLYGON_NUMBER'][cnt]=prp['SILV_POLYGON_NUMBER']
            fcinv_full['SILV_POLYGON_NET_AREA'][cnt]=prp['SILV_POLYGON_NET_AREA']
            if feat['geometry']==None:
                fcinv_full['SpatialMissing'][cnt]=1
            cnt=cnt+1
    
    #--------------------------------------------------------------------------
    # Get missing FC geometries from VRI layer
    #--------------------------------------------------------------------------
      
    # Isolate missing entries
    iMis=np.where(fcinv_full['SpatialMissing']==1)[0]
    
    # Unique list of openings with missing geometry in the FC layer
    dMisFC={}
    dMisFC['Unique Openings with Missing FC Geom']=np.unique(fcinv_full['OPENING_ID'][iMis]).astype('int32')
    
    # Initialize list that will store geometries from VRI corresponding with each
    # unique Opening ID with missing FC layer geometry
    dMisFC['Geom from VRI']=[None]*dMisFC['Unique Openings with Missing FC Geom'].size
    
    # Loop through VRI and retrieve missing geometries
    with fiona.open(PathInVRIFull + '\\VRI.gdb',layer='VEG_COMP_LYR_R1_POLY') as source:
        for feat in source:
            
            geom=feat['geometry']
            
            if (geom==None):
                continue
            
            prp=feat['properties']
            
            if prp['OPENING_ID']==None:
                continue
            
            id=prp['OPENING_ID']
            
            # Index to FC entries that correspond to this feature and have missing geometries
            iMisFC=np.where( (dMisFC['Unique Openings with Missing FC Geom']==prp['OPENING_ID']) )[0]
            
            if iMisFC.size==0:
                continue
        
            # Add area to geom dictionary
            geom['Hectares']=prp['GEOMETRY_Area']/10000
        
            if dMisFC['Geom from VRI'][iMisFC[0]]==None:
                dMisFC['Geom from VRI'][iMisFC[0]]=[geom]
            else:
                dMisFC['Geom from VRI'][iMisFC[0]].append(geom)

    # Save
    gu.opickle(PathInResultsFull + '\\missing_geo_fc_list.pkl',dMisFC)

    print((time.time()-t0)/60)

    #--------------------------------------------------------------------------
    # Add 
    #--------------------------------------------------------------------------
    
#    with fiona.open(PathInResultsFull + '\\Results.gdb',layer='RSLT_FOREST_COVER_INV_SVW') as source:
#        for feat in source:            
#            prp=feat['properties']            
#
#            if feat['geometry']!=None:
#                continue
#
#            iMis_fc=np.where( (dMisFC['Unique Openings with Missing FC Geom']==prp['OPENING_ID']) )[0]
#            iMis_fc=iMis_fc[0]
#            
#            if len(dMisFC['Geom from VRI'][iMis_fc])>1:
#                break
#            
#            D=np.zeros(len(dMisFC['Geom from VRI'][iMis_fc]))
#            for iV in range(len(dMisFC['Geom from VRI'][iMis_fc])):
#                D[iV]=np.abs(prp['SILV_POLYGON_AREA']-dMisFC['Geom from VRI'][iMis_fc][iV]['Hectares'])
#            iMinD=np.where(D==np.min(D))[0]
#            iMinD=iMinD[0]
#            
#            geom=dMisFC['Geom from VRI'][iMis_fc][iMinD]
#            
#            geom1=GetPolygonsFromFionaFeature(geom)
#            
#            gdf=gpd.GeoDataFrame(geom1,crs=gdf_bm.crs)
            
#%% Explore contents of fiona feature

def GetPolygonsFromFionaFeature(geom):
    
    # Polygons stored as numpy arrays
    dPoly=[]
    dPoly_inner=[]
    
    # Stored as shapely geometry object within dictionary, within list
    #gdf=gpd.GeoDataFrame(ListOfPolygons,crs=gdf_bm.crs)
    #gdf=gdf.set_geometry('geometry')
    #gdf_atu_polygons.to_file(filename=meta['Paths']['Project'] + '\\Geospatial\\atu_polygons.geojson',driver='GeoJSON')

    geom_out=[]
    
    cnt=0
    
    coords0=geom['coordinates']
    
    for iPoly in range(len(coords0)):
        
        coords1=coords0[iPoly]
            
        x=[]; 
        y=[];
        for k in range(len(coords1[0])):
            x.append(coords1[0][k][0])
            y.append(coords1[0][k][1])
            
        l=[]
        for k in range(len(x)):
            l.append([x[k],y[k]])
        lp=Polygon(l)             
                    
        A=gis.PolyArea(x,y)/10000
                
        dct={}
        dct['Hectares']=A.copy()
        dct['x']=np.array(x).copy()
        dct['y']=np.array(y).copy()
        dct['x_Centroid']=lp.centroid.x
        dct['y_Centroid']=lp.centroid.y
        dPoly.append(dct.copy())
            
        # Add to geodataframe
        dp0={}
        dp0['ID']=0
            
        gdf_outer=gpd.GeoDataFrame(crs=gdf_bm.crs)
        gdf_outer.loc[0,'geometry']=Polygon(coords1[0])
        
        if len(coords1)>1:
            # Adjust polygons if they have an inner ring
            for j in range(1,len(coords1)):
                gdf_inner=gpd.GeoDataFrame(crs=gdf_bm.crs)
                gdf_inner.loc[0,'geometry']=Polygon(coords1[j])                
                try:
                    # This very rarely receives "'NoneType' object has no attribute 'intersection'"
                    # Not sure why
                    gdf_outer=gpd.overlay(gdf_outer,gdf_inner,how='difference')
                    
                    x=[]; 
                    y=[];
                    for k in range(len(coords1[j])):
                        x.append(coords1[j][k][0])
                        y.append(coords1[j][k][1])
                    dct={}
                    dct['x']=np.array(x).copy()
                    dct['y']=np.array(y).copy()
                    dPoly_inner.append(dct.copy())
                except:
                    pass
            
        dp0['geometry']=gdf_outer.loc[0,'geometry']        
        geom_out.append(dp0)

    #gdf=gpd.GeoDataFrame(ListOfPolygons,crs=gdf_bm.crs)
    #gdf=gdf.set_geometry('geometry')

    return geom_out

#%% ADD PLANTING INFO TO DMEC

# Create function to avoid duplication of bulky code
def AddPlantingWithNoData(d_nd):
    for i in range(1,6):
        d_nd['PL_SPECIES_CD' + str(i)]=np.append(d_nd['PL_SPECIES_CD' + str(i)],-999)
        d_nd['PL_SPECIES_PCT' + str(i)]=np.append(d_nd['PL_SPECIES_PCT' + str(i)],-999)
        d_nd['PL_SPECIES_GW' + str(i)]=np.append(d_nd['PL_SPECIES_GW' + str(i)],-999)
    return d_nd

#%% PREPARE DISTURBANCE MANAGEMENT ENVENT HISTORY

def PrepDMEC(idx,meta,atu,pl,op,fcinv,vri,cut,fire,burnsev,pest):
    
    # Initiate disturbance-management event history
    dmec=[None]*meta['Project']['N Stand']
    
    # Specify flag indicating whether subsetting occurs
    if np.isin('iKeep',list(meta['Project'].keys()))==True:
        # Tile project, only keeping a subset           
        flag_subset=1
    else:
        # Run all entries
        flag_subset=0
    
    for iStand0 in range(meta['Project']['N Stand']):
        
        # Index to stand depends on subsetting
        if flag_subset==1:
            iStand=meta['Project']['iKeep'][iStand0]
        else:
            iStand=iStand0
        
        # Indices to stand
        indS_atu=idx['atu'][iStand]
        indS_pl=idx['pl'][iStand]
        indS_op=idx['op'][iStand]
        indS_burnsev=idx['burnsev'][iStand]
        indS_fire=idx['fire'][iStand]
        indS_cut=idx['cut'][iStand]
        indS_pest=idx['pest'][iStand]
    
        # Initialize dictionary for management history list        
        dmec0={}
        dmec0['Year']=-999*np.ones(1,dtype='float')
        dmec0['Month']=-999*np.ones(1,dtype='int16')
        dmec0['Day']=-999*np.ones(1,dtype='int16')
        dmec0['ID_Type']=-999*np.ones(1,dtype='int16')
        dmec0['MortalityFactor']=-999*np.ones(1,dtype='int16')
        dmec0['GrowthFactor']=-999*np.ones(1,dtype='int16')
        dmec0['SILV_FUND_SOURCE_CODE']=-999*np.ones(1,dtype='int16')
        dmec0['FIA_PROJECT_ID']=-999*np.ones(1,dtype='int16')
        dmec0['OPENING_ID']=-999*np.ones(1,dtype='int16')
        dmec0['ACTUAL_TREATMENT_AREA']=-999*np.ones(1,dtype='int32')
        dmec0['ACTUAL_PLANTED_NUMBER']=-999*np.ones(1,dtype='int32')
        dmec0['PL_SPECIES_CD1']=-999*np.ones(1,dtype='int16')
        dmec0['PL_SPECIES_CD2']=-999*np.ones(1,dtype='int16')
        dmec0['PL_SPECIES_CD3']=-999*np.ones(1,dtype='int16')
        dmec0['PL_SPECIES_CD4']=-999*np.ones(1,dtype='int16')
        dmec0['PL_SPECIES_CD5']=-999*np.ones(1,dtype='int16')
        dmec0['PL_SPECIES_PCT1']=-999*np.ones(1,dtype='int16')
        dmec0['PL_SPECIES_PCT2']=-999*np.ones(1,dtype='int16')
        dmec0['PL_SPECIES_PCT3']=-999*np.ones(1,dtype='int16')
        dmec0['PL_SPECIES_PCT4']=-999*np.ones(1,dtype='int16')
        dmec0['PL_SPECIES_PCT5']=-999*np.ones(1,dtype='int16')
        dmec0['PL_SPECIES_GW1']=-999*np.ones(1,dtype='int16')
        dmec0['PL_SPECIES_GW2']=-999*np.ones(1,dtype='int16')
        dmec0['PL_SPECIES_GW3']=-999*np.ones(1,dtype='int16')
        dmec0['PL_SPECIES_GW4']=-999*np.ones(1,dtype='int16')
        dmec0['PL_SPECIES_GW5']=-999*np.ones(1,dtype='int16')
    
        #--------------------------------------------------------------------------
        # ATU layer
        #--------------------------------------------------------------------------
    
        if indS_atu!=None:
    
            for i in range(indS_atu['Index'].size):
            
                iA=indS_atu['Index'][i]
            
                #ind=np.where(atu['IdxToGrid']==iA)[0]
                
                # Convert RESULTS activity to fcgadgets event type
                dNames={}
                dNames['SILV_BASE_CODE']=cbu.lut_n2s(meta['LUT']['ATU']['SILV_BASE_CODE'],atu['SILV_BASE_CODE'][iA])
                dNames['SILV_TECHNIQUE_CODE']=cbu.lut_n2s(meta['LUT']['ATU']['SILV_TECHNIQUE_CODE'],atu['SILV_TECHNIQUE_CODE'][iA])
                dNames['SILV_METHOD_CODE']=cbu.lut_n2s(meta['LUT']['ATU']['SILV_METHOD_CODE'],atu['SILV_METHOD_CODE'][iA])        
                dNames['SILV_OBJECTIVE_CODE_1']=cbu.lut_n2s(meta['LUT']['ATU']['SILV_OBJECTIVE_CODE_1'],atu['SILV_OBJECTIVE_CODE_1'][iA])
                Name_Type=cbu.QueryResultsActivity(dNames)[0]
                
                try:
                    ID_Type=meta['LUT']['Dist'][Name_Type]
                except:
                    ID_Type=-999
                 
                # Skip surveys (added to avoid surveys being included in tile runs)
                if Name_Type=='Surveys':
                    continue
                
                # Populate DMEC dictionary
                dmec0['Year']=np.append(dmec0['Year'],np.round(atu['Year'][iA]+atu['Month'][iA]/12,decimals=2))
                dmec0['Month']=np.append(dmec0['Month'],atu['Month'][iA])
                dmec0['Day']=np.append(dmec0['Day'],atu['Day'][iA])
                dmec0['ID_Type']=np.append(dmec0['ID_Type'],ID_Type)                
                dmec0['MortalityFactor']=np.append(dmec0['MortalityFactor'],np.array(0,dtype='int16'))
                dmec0['GrowthFactor']=np.append(dmec0['GrowthFactor'],np.array(0,dtype='int16'))
                dmec0['SILV_FUND_SOURCE_CODE']=np.append(dmec0['SILV_FUND_SOURCE_CODE'],atu['SILV_FUND_SOURCE_CODE'][iA])
                dmec0['OPENING_ID']=np.append(dmec0['OPENING_ID'],atu['OPENING_ID'][iA])
                dmec0['FIA_PROJECT_ID']=np.append(dmec0['FIA_PROJECT_ID'],atu['FIA_PROJECT_ID'][iA])
                dmec0['ACTUAL_TREATMENT_AREA']=np.append(dmec0['ACTUAL_TREATMENT_AREA'],np.maximum(1,atu['ACTUAL_TREATMENT_AREA'][iA]))
                dmec0['ACTUAL_PLANTED_NUMBER']=np.append(dmec0['ACTUAL_PLANTED_NUMBER'],atu['ACTUAL_PLANTED_NUMBER'][iA])
                dmec0=AddPlantingWithNoData(dmec0)
    
        #--------------------------------------------------------------------------
        # Planting layer
        # This section isn't adding new entries, just populating empty planting entries.
        #--------------------------------------------------------------------------
    
        if indS_pl!=None:
        
            uYear=np.unique(pl['Year'][indS_pl['Index']])
            for iYear in range(uYear.size):
            
                indYear=np.where(pl['Year'][indS_pl['Index']]==uYear[iYear])[0]
                indYear=indS_pl['Index'][indYear]
    
                # Species code
                cd0=pl['SILV_TREE_SPECIES_CODE'][indYear]
        
                # Number planted      
                n0=pl['NUMBER_PLANTED'][indYear]
            
                # Percent
                pct0=n0/np.sum(n0)*100
            
                # Genetic worth           
                gw0=np.zeros(n0.size)
                for j in range(n0.size):
                    ind=np.where(meta['Param']['BE']['Genetic Worth']['SEEDLOT_NUMBER']==pl['SEEDLOT_NUMBER'][indYear[j]])[0]
                    if ind.size!=0:
                        gw0[j]=meta['Param']['BE']['Genetic Worth']['GENETIC_WORTH_RTNG'][ind[0]]
    
                # Dissolve to unique species combinations and calculate weighted average genetic worth
                cd1=np.unique(cd0)
                pct1=np.zeros(cd1.size)
                gw1=np.zeros(cd1.size)    
                for iCD in range(cd1.size):
                    ind=np.where(cd0==cd1[iCD])[0]
                    pct1[iCD]=np.sum(pct0[ind])            
                    gw1[iCD]=np.sum(pct0[ind]*gw0[ind])/np.sum(pct0[ind])
        
                    # Sort
                    indSorted=np.argsort(-pct1)
                    cd1=cd1[indSorted]
                    pct1=pct1[indSorted]
                    gw1=gw1[indSorted]
        
                # Convert percent to integer
                pct1=np.array(pct1,dtype=int)
        
                ind_dmec=np.where( (np.floor(dmec0['Year'])==uYear[iYear]) & 
                          (dmec0['ID_Type']==meta['LUT']['Dist']['Planting']) )[0]
        
                # Populate structure
                if cd1.size>0:
                    dmec0['PL_SPECIES_CD1'][ind_dmec]=cd1[0]
                    dmec0['PL_SPECIES_PCT1'][ind_dmec]=pct1[0]
                    dmec0['PL_SPECIES_GW1'][ind_dmec]=gw1[0]
                if cd1.size>1:
                    dmec0['PL_SPECIES_CD2'][ind_dmec]=cd1[1]
                    dmec0['PL_SPECIES_PCT2'][ind_dmec]=pct1[1]
                    dmec0['PL_SPECIES_GW2'][ind_dmec]=gw1[1]
                if cd1.size>2:
                    dmec0['PL_SPECIES_CD3'][ind_dmec]=cd1[2]
                    dmec0['PL_SPECIES_PCT3'][ind_dmec]=pct1[2]
                    dmec0['PL_SPECIES_GW3'][ind_dmec]=gw1[2]
                if cd1.size>3:
                    dmec0['PL_SPECIES_CD4'][ind_dmec]=cd1[3]
                    dmec0['PL_SPECIES_PCT4'][ind_dmec]=pct1[3]
                    dmec0['PL_SPECIES_GW4'][ind_dmec]=gw1[3]
                if cd1.size>4:
                    dmec0['PL_SPECIES_CD5'][ind_dmec]=cd1[4]
                    dmec0['PL_SPECIES_PCT5'][ind_dmec]=pct1[4]
                    dmec0['PL_SPECIES_GW5'][ind_dmec]=gw1[4]
    
        #--------------------------------------------------------------------------
        # Wildfire
        # Notes:
        # Burn severity does not come with month. In FCI contect, 2017 fires often 
        # occurred after planting. Assume month is August for now.
        #--------------------------------------------------------------------------
        
        if indS_burnsev!=None:        
        
            indS=indS_burnsev['Index']
            
            for i in range(indS.size):
            
                # Only continue if nothing yet added through burnsev layer
                ind=np.where( ( np.floor(dmec0['Year'])==np.floor(burnsev['FIRE_YEAR'][indS[i]]) ) & (dmec0['ID_Type']==meta['LUT']['Dist']['Wildfire']) )[0]
                if ind.size>0: 
                    continue
                
                bsr=burnsev['BURN_SEVERITY_RATING'][indS[i]]
                if bsr==meta['LUT']['BS']['BURN_SEVERITY_RATING']['Low']: 
                    Severity=50
                elif bsr==meta['LUT']['BS']['BURN_SEVERITY_RATING']['Medium']:
                    Severity=90
                elif bsr==meta['LUT']['BS']['BURN_SEVERITY_RATING']['High']: 
                    Severity=100
                else:
                    Severity=5
            
                Month=8
                dmec0['Year']=np.append(dmec0['Year'],burnsev['FIRE_YEAR'][indS[i]]+Month/12)
                dmec0['Month']=np.append(dmec0['Month'],Month)
                dmec0['Day']=np.append(dmec0['Day'],-999)
                dmec0['ID_Type']=np.append(dmec0['ID_Type'],meta['LUT']['Dist']['Wildfire'])
                dmec0['MortalityFactor']=np.append(dmec0['MortalityFactor'],Severity)
                dmec0['GrowthFactor']=np.append(dmec0['GrowthFactor'],0)
                dmec0['SILV_FUND_SOURCE_CODE']=np.append(dmec0['SILV_FUND_SOURCE_CODE'],0)
                dmec0['OPENING_ID']=np.append(dmec0['OPENING_ID'],-999)
                dmec0['FIA_PROJECT_ID']=np.append(dmec0['FIA_PROJECT_ID'],-999)   
                dmec0['ACTUAL_TREATMENT_AREA']=np.append(dmec0['ACTUAL_TREATMENT_AREA'],-999)
                dmec0['ACTUAL_PLANTED_NUMBER']=np.append(dmec0['ACTUAL_PLANTED_NUMBER'],-999)
                dmec0=AddPlantingWithNoData(dmec0)
    
        if indS_fire!=None:
            
            indS=indS_fire['Index']
        
            for i in range(indS.size):
        
                # Only continue if nothing yet added through burnsev layer
                ind=np.where( ( np.floor(dmec0['Year'])==np.floor(fire['FIRE_YEAR'][indS[i]]) ) & (dmec0['ID_Type']==meta['LUT']['Dist']['Wildfire']) )[0]
                if ind.size>0: 
                    continue
            
                dmec0['Year']=np.append(dmec0['Year'],fire['FIRE_YEAR'][indS[i]]+fire['Month'][indS[i]]/12)
                dmec0['Month']=np.append(dmec0['Month'],fire['Month'][indS[i]])
                dmec0['Day']=np.append(dmec0['Day'],fire['Day'][indS[i]])
                dmec0['ID_Type']=np.append(dmec0['ID_Type'],meta['LUT']['Dist']['Wildfire'])            
                dmec0['MortalityFactor']=np.append(dmec0['MortalityFactor'],50)
                dmec0['GrowthFactor']=np.append(dmec0['GrowthFactor'],0)
                dmec0['SILV_FUND_SOURCE_CODE']=np.append(dmec0['SILV_FUND_SOURCE_CODE'],0)
                dmec0['OPENING_ID']=np.append(dmec0['OPENING_ID'],-999)            
                dmec0['FIA_PROJECT_ID']=np.append(dmec0['FIA_PROJECT_ID'],-999)   
                dmec0['ACTUAL_TREATMENT_AREA']=np.append(dmec0['ACTUAL_TREATMENT_AREA'],-999)
                dmec0['ACTUAL_PLANTED_NUMBER']=np.append(dmec0['ACTUAL_PLANTED_NUMBER'],999)
                dmec0=AddPlantingWithNoData(dmec0)
        
        #--------------------------------------------------------------------------
        # Consolidated cutblocks
        #------------------------------------------------------------------------
    
        if indS_cut!=None:
        
            indS=indS_cut['Index']
        
            for i in range(indS.size):
                
                # Add the harvest
                yr_harv=cut['HARVEST_YEAR'][indS[i]]+1/12
                dmec0['Year']=np.append(dmec0['Year'],yr_harv)
                dmec0['Month']=np.append(dmec0['Month'],-999)
                dmec0['Day']=np.append(dmec0['Day'],-999)
                dmec0['ID_Type']=np.append(dmec0['ID_Type'],meta['LUT']['Dist']['Harvest'])
                dmec0['MortalityFactor']=np.append(dmec0['MortalityFactor'],100)
                dmec0['GrowthFactor']=np.append(dmec0['GrowthFactor'],0)
                dmec0['SILV_FUND_SOURCE_CODE']=np.append(dmec0['SILV_FUND_SOURCE_CODE'],0)
                dmec0['OPENING_ID']=np.append(dmec0['OPENING_ID'],-999)
                dmec0['FIA_PROJECT_ID']=np.append(dmec0['FIA_PROJECT_ID'],-999)   
                dmec0['ACTUAL_TREATMENT_AREA']=np.append(dmec0['ACTUAL_TREATMENT_AREA'],-999)
                dmec0['ACTUAL_PLANTED_NUMBER']=np.append(dmec0['ACTUAL_PLANTED_NUMBER'],-999)
                dmec0=AddPlantingWithNoData(dmec0)
                
                # Add the slashpile burn
                yr_burn=yr_harv+0.1
                dmec0['Year']=np.append(dmec0['Year'],yr_burn)
                dmec0['Month']=np.append(dmec0['Month'],-999)
                dmec0['Day']=np.append(dmec0['Day'],-999)
                dmec0['ID_Type']=np.append(dmec0['ID_Type'],meta['LUT']['Dist']['Slashpile Burn'])
                dmec0['MortalityFactor']=np.append(dmec0['MortalityFactor'],100)
                dmec0['GrowthFactor']=np.append(dmec0['GrowthFactor'],0)
                dmec0['SILV_FUND_SOURCE_CODE']=np.append(dmec0['SILV_FUND_SOURCE_CODE'],0)
                dmec0['OPENING_ID']=np.append(dmec0['OPENING_ID'],-999)
                dmec0['FIA_PROJECT_ID']=np.append(dmec0['FIA_PROJECT_ID'],-999)   
                dmec0['ACTUAL_TREATMENT_AREA']=np.append(dmec0['ACTUAL_TREATMENT_AREA'],-999)
                dmec0['ACTUAL_PLANTED_NUMBER']=np.append(dmec0['ACTUAL_PLANTED_NUMBER'],-999)
                dmec0=AddPlantingWithNoData(dmec0)

        #--------------------------------------------------------------------------
        # Harvesting from RESULTS opening layer
        #--------------------------------------------------------------------------
    
        if indS_op!=None: 
        
            iExisting=np.where( (dmec0['ID_Type']==meta['LUT']['Dist']['Harvest']) )[0]
            YearExisting=np.floor(dmec0['Year'][iExisting])
        
            indS=indS_op['Index']
        
            for i in range(indS.size):
        
                Year=np.floor(op['Year_Denu1_Comp'][indS[i]])
                Month=op['Month_Denu1_Comp'][indS[i]]
                cd=cbu.lut_n2s(meta['LUT']['OP']['DENUDATION_1_DISTURBANCE_CODE'],op['DENUDATION_1_DISTURBANCE_CODE'][indS[i]])
            
                # Determine whether it should be added, or whether it was already in the
                # consolidated cutblocks DB
                flg_add=0
                if (cd=='L') | (cd=='S'):
                    if iExisting.size==0:
                        flg_add=1
                    else:
                        if np.isin(Year,YearExisting)==False:
                            flg_add=1
            
                if flg_add==1:
                    yr_harv=Year+Month/12
                    dmec0['Year']=np.append(dmec0['Year'],yr_harv)
                    dmec0['Month']=np.append(dmec0['Month'],Month)
                    dmec0['Day']=np.append(dmec0['Day'],-999)
                    dmec0['ID_Type']=np.append(dmec0['ID_Type'],meta['LUT']['Dist']['Harvest'])
                    dmec0['MortalityFactor']=np.append(dmec0['MortalityFactor'],100)
                    dmec0['GrowthFactor']=np.append(dmec0['GrowthFactor'],0)
                    dmec0['SILV_FUND_SOURCE_CODE']=np.append(dmec0['SILV_FUND_SOURCE_CODE'],0)
                    dmec0['OPENING_ID']=np.append(dmec0['OPENING_ID'],-999)
                    dmec0['FIA_PROJECT_ID']=np.append(dmec0['FIA_PROJECT_ID'],-999)   
                    dmec0['ACTUAL_TREATMENT_AREA']=np.append(dmec0['ACTUAL_TREATMENT_AREA'],-999)
                    dmec0['ACTUAL_PLANTED_NUMBER']=np.append(dmec0['ACTUAL_PLANTED_NUMBER'],-999)
                    dmec0=AddPlantingWithNoData(dmec0)
                    
                    # Add the slashpile burn
                    yr_burn=yr_harv+0.1
                    dmec0['Year']=np.append(dmec0['Year'],yr_burn)
                    dmec0['Month']=np.append(dmec0['Month'],-999)
                    dmec0['Day']=np.append(dmec0['Day'],-999)
                    dmec0['ID_Type']=np.append(dmec0['ID_Type'],meta['LUT']['Dist']['Slashpile Burn'])
                    dmec0['MortalityFactor']=np.append(dmec0['MortalityFactor'],100)
                    dmec0['GrowthFactor']=np.append(dmec0['GrowthFactor'],0)
                    dmec0['SILV_FUND_SOURCE_CODE']=np.append(dmec0['SILV_FUND_SOURCE_CODE'],0)
                    dmec0['OPENING_ID']=np.append(dmec0['OPENING_ID'],-999)
                    dmec0['FIA_PROJECT_ID']=np.append(dmec0['FIA_PROJECT_ID'],-999)   
                    dmec0['ACTUAL_TREATMENT_AREA']=np.append(dmec0['ACTUAL_TREATMENT_AREA'],-999)
                    dmec0['ACTUAL_PLANTED_NUMBER']=np.append(dmec0['ACTUAL_PLANTED_NUMBER'],-999)
                    dmec0=AddPlantingWithNoData(dmec0)
        
                Year=np.floor(op['Year_Denu2_Comp'][indS[i]])
                Month=op['Month_Denu2_Comp'][indS[i]]
                cd=cbu.lut_n2s(meta['LUT']['OP']['DENUDATION_2_DISTURBANCE_CODE'],op['DENUDATION_2_DISTURBANCE_CODE'][indS[i]])
            
                # Determine whether it should be added, or whether it was already in the
                # consolidated cutblocks DB
                flg_add=0
                if (cd=='L') | (cd=='S'):
                    if iExisting.size==0:
                        flg_add=1
                    else:
                        if np.isin(Year,YearExisting)==False:
                            flg_add=1
            
                if flg_add==1:
                    yr_harv=Year+Month/12
                    dmec0['Year']=np.append(dmec0['Year'],yr_harv)
                    dmec0['Month']=np.append(dmec0['Month'],Month)
                    dmec0['Day']=np.append(dmec0['Day'],-999)
                    dmec0['ID_Type']=np.append(dmec0['ID_Type'],meta['LUT']['Dist']['Harvest'])
                    dmec0['MortalityFactor']=np.append(dmec0['MortalityFactor'],100)
                    dmec0['GrowthFactor']=np.append(dmec0['GrowthFactor'],0)
                    dmec0['SILV_FUND_SOURCE_CODE']=np.append(dmec0['SILV_FUND_SOURCE_CODE'],0)
                    dmec0['OPENING_ID']=np.append(dmec0['OPENING_ID'],-999)
                    dmec0['FIA_PROJECT_ID']=np.append(dmec0['FIA_PROJECT_ID'],-999)   
                    dmec0['ACTUAL_TREATMENT_AREA']=np.append(dmec0['ACTUAL_TREATMENT_AREA'],-999)
                    dmec0['ACTUAL_PLANTED_NUMBER']=np.append(dmec0['ACTUAL_PLANTED_NUMBER'],-999)
                    dmec0=AddPlantingWithNoData(dmec0)
                    
                    # Add the slashpile burn
                    yr_burn=yr_harv+0.1
                    dmec0['Year']=np.append(dmec0['Year'],yr_burn)
                    dmec0['Month']=np.append(dmec0['Month'],-999)
                    dmec0['Day']=np.append(dmec0['Day'],-999)
                    dmec0['ID_Type']=np.append(dmec0['ID_Type'],meta['LUT']['Dist']['Slashpile Burn'])
                    dmec0['MortalityFactor']=np.append(dmec0['MortalityFactor'],100)
                    dmec0['GrowthFactor']=np.append(dmec0['GrowthFactor'],0)
                    dmec0['SILV_FUND_SOURCE_CODE']=np.append(dmec0['SILV_FUND_SOURCE_CODE'],0)
                    dmec0['OPENING_ID']=np.append(dmec0['OPENING_ID'],-999)
                    dmec0['FIA_PROJECT_ID']=np.append(dmec0['FIA_PROJECT_ID'],-999)   
                    dmec0['ACTUAL_TREATMENT_AREA']=np.append(dmec0['ACTUAL_TREATMENT_AREA'],-999)
                    dmec0['ACTUAL_PLANTED_NUMBER']=np.append(dmec0['ACTUAL_PLANTED_NUMBER'],-999)
                    dmec0=AddPlantingWithNoData(dmec0)
    
        #--------------------------------------------------------------------------
        # Knockdown from RESULTS opening layer
        #--------------------------------------------------------------------------
    
        if indS_op!=None: 
        
            iExisting=np.where( (dmec0['ID_Type']==meta['LUT']['Dist']['Knockdown']) )[0]
            YearExisting=np.floor(dmec0['Year'][iExisting])
        
            indS=indS_op['Index']
        
            for i in range(indS.size):
        
                Year=np.floor(op['Year_Denu1_Comp'][indS[i]])
                Month=op['Month_Denu1_Comp'][indS[i]]
                cd=cbu.lut_n2s(meta['LUT']['OP']['DENUDATION_1_DISTURBANCE_CODE'],op['DENUDATION_1_DISTURBANCE_CODE'][indS[i]])
            
                # Determine whether it should be added, or whether it was already in the
                # consolidated cutblocks DB
                flg_add=0
                if (cd=='R'):
                    if iExisting.size==0:
                        flg_add=1
                    else:
                        if np.isin(Year,YearExisting)==False:
                            flg_add=1
            
                if flg_add==1:
                    yr_rehab=Year+Month/12
                    dmec0['Year']=np.append(dmec0['Year'],yr_rehab)
                    dmec0['Month']=np.append(dmec0['Month'],Month)
                    dmec0['Day']=np.append(dmec0['Day'],-999)
                    dmec0['ID_Type']=np.append(dmec0['ID_Type'],meta['LUT']['Dist']['Knockdown'])
                    dmec0['MortalityFactor']=np.append(dmec0['MortalityFactor'],100)
                    dmec0['GrowthFactor']=np.append(dmec0['GrowthFactor'],0)
                    dmec0['SILV_FUND_SOURCE_CODE']=np.append(dmec0['SILV_FUND_SOURCE_CODE'],0)
                    dmec0['OPENING_ID']=np.append(dmec0['OPENING_ID'],-999)
                    dmec0['FIA_PROJECT_ID']=np.append(dmec0['FIA_PROJECT_ID'],-999)   
                    dmec0['ACTUAL_TREATMENT_AREA']=np.append(dmec0['ACTUAL_TREATMENT_AREA'],-999)
                    dmec0['ACTUAL_PLANTED_NUMBER']=np.append(dmec0['ACTUAL_PLANTED_NUMBER'],-999)
                    dmec0=AddPlantingWithNoData(dmec0)
                    
                    # Add the slashpile burn
                    yr_burn=yr_rehab+0.1
                    dmec0['Year']=np.append(dmec0['Year'],yr_burn)
                    dmec0['Month']=np.append(dmec0['Month'],-999)
                    dmec0['Day']=np.append(dmec0['Day'],-999)
                    dmec0['ID_Type']=np.append(dmec0['ID_Type'],meta['LUT']['Dist']['Slashpile Burn'])
                    dmec0['MortalityFactor']=np.append(dmec0['MortalityFactor'],100)
                    dmec0['GrowthFactor']=np.append(dmec0['GrowthFactor'],0)
                    dmec0['SILV_FUND_SOURCE_CODE']=np.append(dmec0['SILV_FUND_SOURCE_CODE'],0)
                    dmec0['OPENING_ID']=np.append(dmec0['OPENING_ID'],-999)
                    dmec0['FIA_PROJECT_ID']=np.append(dmec0['FIA_PROJECT_ID'],-999)   
                    dmec0['ACTUAL_TREATMENT_AREA']=np.append(dmec0['ACTUAL_TREATMENT_AREA'],-999)
                    dmec0['ACTUAL_PLANTED_NUMBER']=np.append(dmec0['ACTUAL_PLANTED_NUMBER'],-999)
                    dmec0=AddPlantingWithNoData(dmec0)
        
                Year=np.floor(op['Year_Denu2_Comp'][indS[i]])
                Month=op['Month_Denu2_Comp'][indS[i]]
                cd=cbu.lut_n2s(meta['LUT']['OP']['DENUDATION_2_DISTURBANCE_CODE'],op['DENUDATION_2_DISTURBANCE_CODE'][indS[i]])
            
                # Determine whether it should be added, or whether it was already in the
                # consolidated cutblocks DB
                flg_add=0
                if (cd=='R'):
                    if iExisting.size==0:
                        flg_add=1
                    else:
                        if np.isin(Year,YearExisting)==False:
                            flg_add=1
            
                if flg_add==1:
                    yr_rehab=Year+Month/12
                    dmec0['Year']=np.append(dmec0['Year'],yr_rehab)
                    dmec0['Month']=np.append(dmec0['Month'],Month)
                    dmec0['Day']=np.append(dmec0['Day'],-999)
                    dmec0['ID_Type']=np.append(dmec0['ID_Type'],meta['LUT']['Dist']['Knockdown'])
                    dmec0['MortalityFactor']=np.append(dmec0['MortalityFactor'],100)
                    dmec0['GrowthFactor']=np.append(dmec0['GrowthFactor'],0)
                    dmec0['SILV_FUND_SOURCE_CODE']=np.append(dmec0['SILV_FUND_SOURCE_CODE'],0)
                    dmec0['OPENING_ID']=np.append(dmec0['OPENING_ID'],-999)
                    dmec0['FIA_PROJECT_ID']=np.append(dmec0['FIA_PROJECT_ID'],-999)   
                    dmec0['ACTUAL_TREATMENT_AREA']=np.append(dmec0['ACTUAL_TREATMENT_AREA'],-999)
                    dmec0['ACTUAL_PLANTED_NUMBER']=np.append(dmec0['ACTUAL_PLANTED_NUMBER'],-999)
                    dmec0=AddPlantingWithNoData(dmec0)
                    
                    # Add the slashpile burn
                    yr_burn=yr_rehab+0.1
                    dmec0['Year']=np.append(dmec0['Year'],yr_burn)
                    dmec0['Month']=np.append(dmec0['Month'],-999)
                    dmec0['Day']=np.append(dmec0['Day'],-999)
                    dmec0['ID_Type']=np.append(dmec0['ID_Type'],meta['LUT']['Dist']['Slashpile Burn'])
                    dmec0['MortalityFactor']=np.append(dmec0['MortalityFactor'],100)
                    dmec0['GrowthFactor']=np.append(dmec0['GrowthFactor'],0)
                    dmec0['SILV_FUND_SOURCE_CODE']=np.append(dmec0['SILV_FUND_SOURCE_CODE'],0)
                    dmec0['OPENING_ID']=np.append(dmec0['OPENING_ID'],-999)
                    dmec0['FIA_PROJECT_ID']=np.append(dmec0['FIA_PROJECT_ID'],-999)   
                    dmec0['ACTUAL_TREATMENT_AREA']=np.append(dmec0['ACTUAL_TREATMENT_AREA'],-999)
                    dmec0['ACTUAL_PLANTED_NUMBER']=np.append(dmec0['ACTUAL_PLANTED_NUMBER'],-999)
                    dmec0=AddPlantingWithNoData(dmec0)
    
        #--------------------------------------------------------------------------    
        # Forest health
        # Adjust severity later based on host species
        #--------------------------------------------------------------------------
    
        # Turn on or off for QA purposes
        flg=1
    
        if (indS_pest!=None) & (flg==1):
        
            indS=indS_pest['Index']
        
            # See above section for a breakdown of the most abundant pests
            for i in range(indS.size):
            
                iYr=indS[i]
            
                if (pest['CAPTURE_YEAR'][iYr]==0):
                    continue                    
            
                # Don't include trace
                if (pest['PEST_SEVERITY_CODE'][iYr]==meta['LUT']['Pest']['PEST_SEVERITY_CODE']['T'][0]):
                    continue
            
                dmec0['Year']=np.append(dmec0['Year'],pest['CAPTURE_YEAR'][iYr])
                dmec0['Month']=np.append(dmec0['Month'],-999)
                dmec0['Day']=np.append(dmec0['Day'],-999)
                dmec0['SILV_FUND_SOURCE_CODE']=np.append(dmec0['SILV_FUND_SOURCE_CODE'],0)
                dmec0['OPENING_ID']=np.append(dmec0['OPENING_ID'],-999)
                dmec0['FIA_PROJECT_ID']=np.append(dmec0['FIA_PROJECT_ID'],-999)   
                dmec0['ACTUAL_TREATMENT_AREA']=np.append(dmec0['ACTUAL_TREATMENT_AREA'],-999)
                dmec0['ACTUAL_PLANTED_NUMBER']=np.append(dmec0['ACTUAL_PLANTED_NUMBER'],-999)
                dmec0=AddPlantingWithNoData(dmec0)
            
                # Shorten for repeated use
                psp=pest['PEST_SPECIES_CODE'][iYr]
                sev=pest['PEST_SEVERITY_CODE'][iYr]
                sev_s=cbu.lut_n2s(meta['LUT']['Pest']['PEST_SEVERITY_CODE'],sev)[0]
            
                # Mountain pine beetle - Adjust severity based on fractin of pine later using: AdjustSpeciesSpecificMortality         
                if (psp==meta['LUT']['Pest']['PEST_SPECIES_CODE']['IBM']):
                    dmec0['ID_Type']=np.append(dmec0['ID_Type'],meta['LUT']['Dist']['IBM'])
                    ind=np.where( (meta['Param']['BE']['DistBySC']['Name']=='IBM') & (meta['Param']['BE']['DistBySC']['SeverityCD']==sev_s) )[0]
                    dmec0['MortalityFactor']=np.append(dmec0['MortalityFactor'],meta['Param']['BE']['DistBySC']['MortalityFactor'][ind])
                    dmec0['GrowthFactor']=np.append(dmec0['GrowthFactor'],meta['Param']['BE']['DistBySC']['GrowthFactor'][ind])
                
                # Western Balsam Bark Beetle
                # (https://www2.gov.bc.ca/gov/content/industry/forestry/managing-our-forest-resources/forest-health/forest-pests/bark-beetles/western-balsam-bark-beetle)
                elif (psp==meta['LUT']['Pest']['PEST_SPECIES_CODE']['IBB']):
                    dmec0['ID_Type']=np.append(dmec0['ID_Type'],meta['LUT']['Dist']['IBB'])
                    ind=np.where( (meta['Param']['BE']['DistBySC']['Name']=='IBB') & (meta['Param']['BE']['DistBySC']['SeverityCD']==sev_s) )[0]
                    dmec0['MortalityFactor']=np.append(dmec0['MortalityFactor'],meta['Param']['BE']['DistBySC']['MortalityFactor'][ind])
                    dmec0['GrowthFactor']=np.append(dmec0['GrowthFactor'],meta['Param']['BE']['DistBySC']['GrowthFactor'][ind])
            
                # Douglas-fir beetle
                # (https://www2.gov.bc.ca/gov/content/industry/forestry/managing-our-forest-resources/forest-health/forest-pests/bark-beetles/western-balsam-bark-beetle)
                elif (psp==meta['LUT']['Pest']['PEST_SPECIES_CODE']['IBD']):
                    dmec0['ID_Type']=np.append(dmec0['ID_Type'],meta['LUT']['Dist']['IBD'])
                    ind=np.where( (meta['Param']['BE']['DistBySC']['Name']=='IBD') & (meta['Param']['BE']['DistBySC']['SeverityCD']==sev_s) )[0]
                    dmec0['MortalityFactor']=np.append(dmec0['MortalityFactor'],meta['Param']['BE']['DistBySC']['MortalityFactor'][ind])
                    dmec0['GrowthFactor']=np.append(dmec0['GrowthFactor'],meta['Param']['BE']['DistBySC']['GrowthFactor'][ind])
                
                # Spruce beetle
                # (https://www2.gov.bc.ca/gov/content/industry/forestry/managing-our-forest-resources/forest-health/forest-pests/bark-beetles/western-balsam-bark-beetle)
                elif (psp==meta['LUT']['Pest']['PEST_SPECIES_CODE']['IBS']):
                    dmec0['ID_Type']=np.append(dmec0['ID_Type'],meta['LUT']['Dist']['IBS'])
                    ind=np.where( (meta['Param']['BE']['DistBySC']['Name']=='IBS') & (meta['Param']['BE']['DistBySC']['SeverityCD']==sev_s) )[0]
                    dmec0['MortalityFactor']=np.append(dmec0['MortalityFactor'],meta['Param']['BE']['DistBySC']['MortalityFactor'][ind])
                    dmec0['GrowthFactor']=np.append(dmec0['GrowthFactor'],meta['Param']['BE']['DistBySC']['GrowthFactor'][ind])
                
                # Western spruce budworm 
                # Populate severity with the ID for severity - it will be revised below
                # to model mortality that occurs from repeated infestation
                elif (psp==meta['LUT']['Pest']['PEST_SPECIES_CODE']['IDW']):
                    dmec0['ID_Type']=np.append(dmec0['ID_Type'],meta['LUT']['Dist']['IDW'])
                    ind=np.where( (meta['Param']['BE']['DistBySC']['Name']=='IDW') & (meta['Param']['BE']['DistBySC']['SeverityCD']==sev_s) )[0]
                    dmec0['MortalityFactor']=np.append(dmec0['MortalityFactor'],sev)
                    dmec0['GrowthFactor']=np.append(dmec0['GrowthFactor'],meta['Param']['BE']['DistBySC']['GrowthFactor'][ind])
                
                # Other
                else:
                    dmec0['ID_Type']=np.append(dmec0['ID_Type'],meta['LUT']['Dist']['Beetles'])
                    dmec0['MortalityFactor']=np.append(dmec0['MortalityFactor'],5)
                    dmec0['GrowthFactor']=np.append(dmec0['GrowthFactor'],0)
    
        #--------------------------------------------------------------------------
        # Append to management history list
        #--------------------------------------------------------------------------
    
        dmec[iStand0]=dmec0
    
    return dmec

#%% Is FCI funded in DMEC 
    
def IsFCIFunding(meta,dmec,uPP):
    
    # Convert list of unique FCI PP number codes to IDs - it could be done with codes
    # but this is much faster
    id_uPP=np.zeros(uPP.size)
    for i in range(uPP.size):
        try:
            id_uPP[i]=meta['LUT']['ATU']['FIA_PROJECT_ID'][uPP[i]][0]
        except:
            pass
    
    for iS in range(len(dmec)):
        dmec[iS]['FCI Funded']=np.zeros(dmec[iS]['Year'].size)
        for iY in range(dmec[iS]['Year'].size):
            if dmec[iS]['SILV_FUND_SOURCE_CODE'][iY]==meta['LUT']['ATU']['SILV_FUND_SOURCE_CODE']['FCE']:
                dmec[iS]['FCI Funded'][iY]=1
                continue
            if dmec[iS]['SILV_FUND_SOURCE_CODE'][iY]==meta['LUT']['ATU']['SILV_FUND_SOURCE_CODE']['FCM']:
                dmec[iS]['FCI Funded'][iY]=1            
                continue
            if np.isin(dmec[iS]['FIA_PROJECT_ID'][iY],id_uPP)==True:
                dmec[iS]['FCI Funded'][iY]=1
    
    return dmec

#%% LOAD Look-up-tables

def Load_LUTs(meta):
    
    # Initialize LUTs dictionary
    if 'LUT' not in meta:
        meta['LUT']={}
    
    # Import distubance type
    p=gu.ReadExcel(meta['Paths']['Model Code'] + '\\Parameters\\' + 'Parameters_Disturbances.xlsx')
     
    meta['LUT']['Dist']={}
    for i in range(p['Name'].size):
        meta['LUT']['Dist'][p['Name'][i]]=p['ID'][i]
    
    # BGC zone     
    #LUT_BGC_Zone={}
    #for i in range(len(par['BGC_ZONE']['CODE_BGC_ZONE'])):
    #    LUT_BGC_Zone[par['BGC_ZONE']['CODE_BGC_ZONE'][i]]=par['BGC_ZONE']['ID_BGC_ZONE'][i]
    
    # Added this to accommodate jupyter notebook demos - will need updating periodically
    if 'Results' not in meta['Paths']:
        meta['Paths']['Results']=r'C:\Users\rhember\Documents\Data\ForestInventory\Results\20210930'
    if 'VRI' not in meta['Paths']:
        meta['Paths']['VRI']=r'C:\Users\rhember\Documents\Data\ForestInventory\VRI\20210930'
    if 'Disturbances' not in meta['Paths']:
        meta['Paths']['Disturbances']=r'C:\Users\rhember\Documents\Data\ForestInventory\Disturbances\20210930'
    if 'LandUse' not in meta['Paths']:
        meta['Paths']['LandUse']=r'C:\Users\rhember\Documents\Data\ForestInventory\LandUse\20210930'
    
    meta['LUT']['ATU']=gu.ipickle(meta['Paths']['Results'] + '\\LUTs_RSLT_ACTIVITY_TREATMENT_SVW.pkl')
    meta['LUT']['OP']=gu.ipickle(meta['Paths']['Results'] + '\\LUTs_RSLT_OPENING_SVW.pkl')
    meta['LUT']['PL']=gu.ipickle(meta['Paths']['Results'] + '\\LUTs_RSLT_PLANTING_SVW.pkl')
    meta['LUT']['FC_I']=gu.ipickle(meta['Paths']['Results'] + '\\LUTs_RSLT_FOREST_COVER_INV_SVW.pkl')
    meta['LUT']['FC_S']=gu.ipickle(meta['Paths']['Results'] + '\\LUTs_RSLT_FOREST_COVER_SILV_SVW.pkl')
    meta['LUT']['VRI']=gu.ipickle(meta['Paths']['VRI'] + '\\LUTs_VEG_COMP_LYR_R1_POLY.pkl')
    meta['LUT']['BS']=gu.ipickle(meta['Paths']['Disturbances'] + '\\LUTs_VEG_BURN_SEVERITY_SP.pkl')
    meta['LUT']['Pest']=gu.ipickle(meta['Paths']['Disturbances'] + '\\LUTs_PEST_INFESTATION_POLY.pkl')
    try:
        meta['LUT']['FC_R']=gu.ipickle(meta['Paths']['Results'] + '\\LUTs_RSLT_FOREST_COVER_RESERVE_SVW.pkl')
        meta['LUT']['LU NL']=gu.ipickle(meta['Paths']['LandUse'] + '\\LUTs_RMP_PLAN_NON_LEGAL_POLY_SVW.pkl')
        meta['LUT']['LU L']=gu.ipickle(meta['Paths']['LandUse'] + '\\LUTs_RMP_PLAN_LEGAL_POLY_SVW.pkl')
        meta['LUT']['PARK']=gu.ipickle(meta['Paths']['LandUse'] + '\\LUTs_TA_PARK_ECORES_PA_SVW.pkl')
        meta['LUT']['OGMA']=gu.ipickle(meta['Paths']['LandUse'] + '\\LUTs_RMP_OGMA_LEGAL_CURRENT_SVW.pkl')
        meta['LUT']['UWR']=gu.ipickle(meta['Paths']['LandUse'] + '\\LUTs_WCP_UNGULATE_WINTER_RANGE_SP.pkl')   
    except:
        pass
    
    
    meta['LUT']['TIPSY']={}
    meta['LUT']['TIPSY']['FIZ']={'C':np.array(1,dtype=int),'I':np.array(2,dtype=int)}
    meta['LUT']['TIPSY']['regeneration_method']={'C':np.array(1,dtype=int),'N':np.array(2,dtype=int),'P':np.array(3,dtype=int)}
    
    # Species (for Sawtooth)
    #meta['LUT']['SRS']={}
    #for i in range(len(par['SRS']['SRS_CD'])):
    #    meta['LUT']['SRS'][par['SRS']['SRS_CD'][i]]=par['SRS']['SRS_ID'][i]
    
    return meta

#%% LOAD PARAMTERS

def Load_Params(meta):
    
    par={}
    
    par['Dist']=gu.ReadExcel(meta['Paths']['Model Code'] + '\\Parameters\\Parameters_Disturbances.xlsx')
    
    par['DistBySC']=gu.ReadExcel(meta['Paths']['Model Code'] + '\\Parameters\\Parameters_DisturbanceBySeverityClass.xlsx')
    
    par['GW']=gu.ReadExcel(meta['Paths']['Model Code'] + '\\Parameters\\Parameters_Seedlot_GW.xlsx')
    
    par['WF']={}
    wf=gu.ReadExcel(meta['Paths']['Model Code'] + '\\Parameters\\Parameters_WildfireStatsMod.xlsx')
    for i in range(wf['Name'].size):
        try:
            par['WF'][wf['Name'][i]]=wf['Value'][i].astype(float)
        except:
            par['WF'][wf['Name'][i]]=wf['Value'][i]
            
    par['IBM']={}
    ibm=gu.ReadExcel(meta['Paths']['Model Code'] + '\\Parameters\\Parameters_IBMStatsMod.xlsx')
    for i in range(ibm['Name'].size):
        try:
            par['IBM'][ibm['Name'][i]]=ibm['Value'][i].astype(float)
        except:
            par['IBM'][ibm['Name'][i]]=ibm['Value'][i]        
         
    # Import historical harvest reconstruction (really basic)
    par['Ph_Simp']=gu.ReadExcel(meta['Paths']['Taz Datasets'] + '\\Harvest Stats and Scenarios\\HarvestHistoricalProbabilitySimple.xlsx')
            
    return par

#%% EXCLUDE DUPLICATE EVENTS

def Exclude_Duplicate_Events(meta,dmec):
    for iStand in range(meta['Project']['N Stand']):
        if dmec[iStand]==None:
            continue
        for key in meta['LUT']['Dist'].keys():
            indType=np.where(dmec[iStand]['ID_Type']==meta['LUT']['Dist'][key])[0]
            if indType.size==0:
                continue            
            uYear=np.unique(np.floor(dmec[iStand]['Year'][indType]))
            for iYear in range(uYear.size):
                ind=np.where( (dmec[iStand]['ID_Type']==meta['LUT']['Dist'][key]) & (np.floor(dmec[iStand]['Year'])==uYear[iYear]) )[0]
                dmec[iStand]['ID_Type'][ind[1:]]=-999
    return dmec

#%% EXCLUDE UNIDENTIFIED EVENTS

def Exclude_Unidentified_Events(meta,dmec):

    for iStand in range(meta['Project']['N Stand']):
        if dmec[iStand]==None:
            continue
        ind=np.where(dmec[iStand]['ID_Type']!=-999)[0]
        for key in dmec[iStand]:
            dmec[iStand][key]=dmec[iStand][key][ind]
    return dmec

#%% REMOVE SLASHPILE BURNS IN SELECT BGC ZONES

def Remove_SlashpileBurns_From_Select_Zones(meta,dmec,ba):

    for iStand in range(meta['Project']['N Stand']):
        
        if dmec[iStand]==None:
            continue
        
        if (ba['BEC_ZONE_CODE'][iStand]==meta['LUT']['VRI']['BEC_ZONE_CODE']['CWH']) | (ba['BEC_ZONE_CODE'][iStand]==meta['LUT']['VRI']['BEC_ZONE_CODE']['ICH']):
            ind=np.where(dmec[iStand]['ID_Type']!=meta['LUT']['Dist']['Slashpile Burn'])[0]
            if ind.size>0:
                for key in dmec[iStand]:
                    if key=='ScnAffected':
                        continue
                    dmec[iStand][key]=dmec[iStand][key][ind]
    
    return dmec

#%% ENSURE EVERY STAND HAS A MODERN DISTURBANCE

def Ensure_Every_Stand_Has_Modern_Disturbance(meta,dmec,name_dist,severity,StringsToFill):
    
    for iStand in range(meta['Project']['N Stand']):
        
        if dmec[iStand]==None:
            continue
        
        if dmec[iStand]['Year'].size==0:
            #print(iStand)
            #break
            r=np.random.randint(1700,2000)
            dmec[iStand]['Year']=np.append(dmec[iStand]['Year'],r)
            dmec[iStand]['ID_Type']=np.append(dmec[iStand]['ID_Type'],meta['LUT']['Dist'][name_dist])
            dmec[iStand]['MortalityFactor']=np.append(dmec[iStand]['MortalityFactor'],np.array(severity,dtype='int16'))
            dmec[iStand]['GrowthFactor']=np.append(dmec[iStand]['GrowthFactor'],np.array(0,dtype='int16'))
            for v in StringsToFill:
                dmec[iStand][v]=np.append(dmec[iStand][v],-999)
                
    return dmec

#%% ENSURE DISTURBANCE PRECEDES AERIAL FERTILIZATION
# So that age at fert is specified.

def Ensure_Fert_Preceded_By_Disturbance(meta,dmec,ba,StringsToFill):

    ListOfTestedDist=[meta['LUT']['Dist']['Wildfire'],meta['LUT']['Dist']['Harvest'],
            meta['LUT']['Dist']['Knockdown'],meta['LUT']['Dist']['Harvest Salvage'],
            meta['LUT']['Dist']['Beetles'],meta['LUT']['Dist']['IBM'],meta['LUT']['Dist']['IBB'],
            meta['LUT']['Dist']['IBD'],meta['LUT']['Dist']['IBS']]
    
    for iStand in range(meta['Project']['N Stand']):
        
        if dmec[iStand]==None:
            continue
    
        iFert=np.where( (dmec[iStand]['ID_Type']==meta['LUT']['Dist']['Fertilization Aerial']) )[0]
        
        if iFert.size==0: 
            continue
        
        iFert=iFert[0]
        
        # Index to events prior to first fertilization with 100% mortality
        ind=np.where( (dmec[iStand]['Year']<=dmec[iStand]['Year'][iFert]) & (dmec[iStand]['MortalityFactor']==100) & np.isin(dmec[iStand]['ID_Type'],ListOfTestedDist) )[0]
        
        if (ind.size==0):
            
            # Gapfill with VRI
            Year=ba['Year'][iStand]-ba['PROJ_AGE_1'][iStand]
            
            if Year<=0:
                # Assume mean of 36 + random variation
                r=36+np.random.randint(-6,high=6)
                Year=dmec[iStand]['Year'][iFert]-r
            
            # Add harvest
            dmec[iStand]['Year']=np.append(dmec[iStand]['Year'],Year)
            dmec[iStand]['ID_Type']=np.append(dmec[iStand]['ID_Type'],meta['LUT']['Dist']['Harvest'])
            dmec[iStand]['MortalityFactor']=np.append(dmec[iStand]['MortalityFactor'],100)
            dmec[iStand]['GrowthFactor']=np.append(dmec[iStand]['GrowthFactor'],np.array(0,dtype='int16'))
            for v in StringsToFill:
                dmec[iStand][v]=np.append(dmec[iStand][v],-999)  
            
            # Add slashpile burn
            dmec[iStand]['Year']=np.append(dmec[iStand]['Year'],Year+0.1)
            dmec[iStand]['ID_Type']=np.append(dmec[iStand]['ID_Type'],meta['LUT']['Dist']['Slashpile Burn'])
            dmec[iStand]['MortalityFactor']=np.append(dmec[iStand]['MortalityFactor'],100)
            dmec[iStand]['GrowthFactor']=np.append(dmec[iStand]['GrowthFactor'],np.array(0,dtype='int16'))
            for v in StringsToFill:
                dmec[iStand][v]=np.append(dmec[iStand][v],-999)  
            
            # Add planting
            dmec[iStand]['Year']=np.append(dmec[iStand]['Year'],Year+0.2)
            dmec[iStand]['ID_Type']=np.append(dmec[iStand]['ID_Type'],meta['LUT']['Dist']['Planting'])
            dmec[iStand]['MortalityFactor']=np.append(dmec[iStand]['MortalityFactor'],0)
            dmec[iStand]['GrowthFactor']=np.append(dmec[iStand]['GrowthFactor'],0)
            for v in StringsToFill:
                dmec[iStand][v]=np.append(dmec[iStand][v],-999)
    
    return dmec

#def Ensure_Fert_Preceded_By_Disturbance_old(meta,dmec,th_sev_last_dist,AgeAtFert,StringsToFill):
#
#    ListOfTestedDist=[meta['LUT']['Dist']['Wildfire'],meta['LUT']['Dist']['Harvest'],
#            meta['LUT']['Dist']['Knockdown'],meta['LUT']['Dist']['Harvest Salvage'],
#            meta['LUT']['Dist']['Beetles'],meta['LUT']['Dist']['IBM'],meta['LUT']['Dist']['IBB'],
#            meta['LUT']['Dist']['IBD'],meta['LUT']['Dist']['IBS']]
#    
#    for iStand in range(meta['Project']['N Stand']):
#        
#        if dmec[iStand]==None:
#            continue
#    
#        iA=np.where( (dmec[iStand]['ID_Type']==meta['LUT']['Dist']['Fertilization Aerial']) )[0]
#        if iA.size==0: 
#            continue
#        
#        # Index to events prior to first fertilization with 100% mortality
#        ind=np.where( (dmec[iStand]['Year']<=dmec[iStand]['Year'][iA[0]]) & (dmec[iStand]['MortalityFactor']==100) & np.isin(dmec[iStand]['ID_Type'],ListOfTestedDist) )[0]
#        
#        # Add random component to age at fert
#        AgeAtFertPlusRandom=AgeAtFert+np.random.randint(-6,high=6)
#        
#        #dYear=dmec[iStand]['Year'][iA[0]]-dmec[iStand]['Year'][ind[-1]]
#        if (ind.size==0):
#            
#            Year=dmec[iStand]['Year'][iA[0]]-AgeAtFertPlusRandom
#            
#            # Add harvest
#            dmec[iStand]['Year']=np.append(dmec[iStand]['Year'],Year)
#            dmec[iStand]['ID_Type']=np.append(dmec[iStand]['ID_Type'],meta['LUT']['Dist']['Harvest'])
#            dmec[iStand]['MortalityFactor']=np.append(dmec[iStand]['MortalityFactor'],100)
#            dmec[iStand]['GrowthFactor']=np.append(dmec[iStand]['GrowthFactor'],np.array(0,dtype='int16'))
#            for v in StringsToFill:
#                dmec[iStand][v]=np.append(dmec[iStand][v],-999)  
#            
#            # Add slashpile burn
#            dmec[iStand]['Year']=np.append(dmec[iStand]['Year'],Year+0.1)
#            dmec[iStand]['ID_Type']=np.append(dmec[iStand]['ID_Type'],meta['LUT']['Dist']['Slashpile Burn'])
#            dmec[iStand]['MortalityFactor']=np.append(dmec[iStand]['MortalityFactor'],100)
#            dmec[iStand]['GrowthFactor']=np.append(dmec[iStand]['GrowthFactor'],np.array(0,dtype='int16'))
#            for v in StringsToFill:
#                dmec[iStand][v]=np.append(dmec[iStand][v],-999)  
#            
#            # Add planting
#            dmec[iStand]['Year']=np.append(dmec[iStand]['Year'],Year+0.2)
#            dmec[iStand]['ID_Type']=np.append(dmec[iStand]['ID_Type'],meta['LUT']['Dist']['Planting'])
#            dmec[iStand]['MortalityFactor']=np.append(dmec[iStand]['MortalityFactor'],0)
#            dmec[iStand]['GrowthFactor']=np.append(dmec[iStand]['GrowthFactor'],0)
#            for v in StringsToFill:
#                dmec[iStand][v]=np.append(dmec[iStand][v],-999)            
#    
#    return dmec

#%% IDW FIX SEVERITY
# The dmec was populated with the numeric severity ID. Mortality only occurs 
# following repeated outrbreak years. 

def IDW_Fix_Severity(meta,dmec):

    for iStand in range(meta['Project']['N Stand']):
        
        if dmec[iStand]==None:
            continue
        
        # Save a frozen version of the initial severity
        Severity_Frozen=dmec[iStand]['MortalityFactor'].copy()
    
        for iA in range(dmec[iStand]['Year'].size):
            
            if dmec[iStand]['ID_Type'][iA]=='IDW':
            
                # By default, apply mortality rate in the year of defoliation
                sev1=dmec[iStand]['MortalityFactor'][iA]
                sev_s1=cbu.lut_n2s(meta['LUT']['Pest']['PEST_SEVERITY_CODE'],sev1)[0]
                ind=np.where( (meta['Param']['BE']['DistBySC']['Name']=='IDW') & (meta['Param']['BE']['DistBySC']['SeverityCD']==sev_s1) )[0]
                Mortality1=meta['Param']['BE']['DistBySV']['MortalityFactor'][ind]

                if iA>0:
                    if dmec[iStand]['ID_Type'][iA-1]=='IDW':
                    
                        # If it is back to back infestation, adjust mortality accordingly
                        sev0=Severity_Frozen[iA-1]
                        sev_s0=cbu.lut_n2s(meta['LUT']['Pest']['PEST_SEVERITY_CODE'],sev0)[0]
                        if (sev_s0=='M') & (sev_s1=='M'):
                            ind=np.where( (meta['Param']['BE']['DistBySC']['Name']=='IDW') & (meta['Param']['BE']['DistBySC']['SeverityCD']=='MM') )[0]
                            Mortality1=meta['Param']['BE']['DistBySV']['MortalityFactor'][ind]
                        elif (sev_s0=='M') & (sev_s1=='S'):
                            ind=np.where( (meta['Param']['BE']['DistBySC']['Name']=='IDW') & (meta['Param']['BE']['DistBySC']['SeverityCD']=='MS') )[0]
                            Mortality1=meta['Param']['BE']['DistBySV']['MortalityFactor'][ind]
                        elif (sev_s0=='M') & (sev_s1=='V'):
                            ind=np.where( (meta['Param']['BE']['DistBySC']['Name']=='IDW') & (meta['Param']['BE']['DistBySC']['SeverityCD']=='MV') )[0]
                            Mortality1=meta['Param']['BE']['DistBySV']['MortalityFactor'][ind]
                        elif (sev_s0=='S') & (sev_s1=='S'):
                            ind=np.where( (meta['Param']['BE']['DistBySC']['Name']=='IDW') & (meta['Param']['BE']['DistBySC']['SeverityCD']=='SS') )[0]
                            Mortality1=meta['Param']['BE']['DistBySV']['MortalityFactor'][ind]
                        elif (sev_s0=='S') & (sev_s1=='V'):
                            ind=np.where( (meta['Param']['BE']['DistBySC']['Name']=='IDW') & (meta['Param']['BE']['DistBySC']['SeverityCD']=='SV') )[0]
                            Mortality1=meta['Param']['BE']['DistBySV']['MortalityFactor'][ind]
                        elif (sev_s0=='V') & (sev_s1=='V'):
                            ind=np.where( (meta['Param']['BE']['DistBySC']['Name']=='IDW') & (meta['Param']['BE']['DistBySC']['SeverityCD']=='VV') )[0]
                            Mortality1=meta['Param']['BE']['DistBySV']['MortalityFactor'][ind]
                        else:
                            pass                    
                dmec[iStand]['MortalityFactor'][iA]=np.array(Mortality1,dtype='int16')
    return dmec

#%% ADD OLDEST KNOWN DISTURBANCE FROM VRI
# This really doesn't work well. RETIRED!!

def Add_Oldest_Disturbance_From_VRI(meta,dmec,idx,vri):
    
    for iStand in range(meta['Project']['N Stand']):
        if idx['vri'][iStand]==None:
            continue
        ind=idx['vri'][iStand]['Index'][0]
        DOE=vri['Year'][ind]-vri['PROJ_AGE_1'][ind]
        flg=0
        if dmec[iStand]['Year'].size==0:
            flg==1
        else:
            if (DOE>0) & (DOE<np.min(dmec[iStand]['Year'])):
                flg==1
        if flg==1:
            dmec[iStand]['Year']=np.append(dmec[iStand]['Year'],DOE)
            dmec[iStand]['ID_Type']=np.append(dmec[iStand]['ID_Type'],meta['LUT']['Dist']['Wildfire'])
            dmec[iStand]['MortalityFactor']=np.append(dmec[iStand]['MortalityFactor'],np.array(100,dtype='int16'))
            dmec[iStand]['GrowthFactor']=np.append(dmec[iStand]['GrowthFactor'],np.array(0,dtype='int16'))
            for v in StringsToFill:
                dmec[iStand][v]=np.append(dmec[iStand][v],-999)

    # Put events in order of calendar date
    for iStand in range(meta['Project']['N Stand']):
        d=dmec[iStand].copy()
        ord=np.argsort(d['Year'])
        for key in d.keys():
            d[key]=d[key][ord]
        dmec[iStand]=d.copy()
        
    return dmec

#%% CLEAN SPECIES COMPOSITION

def Clean_Species_Composition(meta,dmec,vri,fcinv):
    
    # List of code pairs (code to change, code to change to)
    ListS=[('SXW','SX')]
    
    # Fix dmec
    for iStand in range(meta['Project']['N Stand']):   
        for iSpc in range(len(ListS)):       
            # Disturbance/management inventory
            n0=meta['LUT']['PL']['SILV_TREE_SPECIES_CODE'][ListS[iSpc][0]]
            n1=meta['LUT']['PL']['SILV_TREE_SPECIES_CODE'][ListS[iSpc][1]]
            for k in range(dmec[iStand]['PL_SPECIES_CD1'].size):
                if dmec[iStand]['PL_SPECIES_CD1'][k]==n0: dmec[iStand]['PL_SPECIES_CD1'][k]=n1
                if dmec[iStand]['PL_SPECIES_CD2'][k]==n0: dmec[iStand]['PL_SPECIES_CD2'][k]=n1
                if dmec[iStand]['PL_SPECIES_CD3'][k]==n0: dmec[iStand]['PL_SPECIES_CD3'][k]=n1
                if dmec[iStand]['PL_SPECIES_CD4'][k]==n0: dmec[iStand]['PL_SPECIES_CD4'][k]=n1
                if dmec[iStand]['PL_SPECIES_CD5'][k]==n0: dmec[iStand]['PL_SPECIES_CD5'][k]=n1
 
    # Fix VRI
    for iStand in range(vri['SPECIES_CD_1'].size):       
        for iSpc in range(len(ListS)):    
            for k in range(6):
                n0=meta['LUT']['VRI']['SPECIES_CD_' + str(k+1)][ListS[iSpc][0]];
                n1=meta['LUT']['VRI']['SPECIES_CD_' + str(k+1)][ListS[iSpc][1]]        
                if vri['SPECIES_CD_' + str(k+1)][iStand]==n0: 
                    vri['SPECIES_CD_' + str(k+1)][iStand]=n1
        
    # Forest cover
    for iStand in range(fcinv['I_SPECIES_CODE_1'].size):       
        for iSpc in range(len(ListS)):         
            for k in range(5):
                n0=meta['LUT']['FC_I']['I_SPECIES_CODE_' + str(k+1)][ListS[iSpc][0]]
                n1=meta['LUT']['FC_I']['I_SPECIES_CODE_' + str(k+1)][ListS[iSpc][1]]        
                if fcinv['I_SPECIES_CODE_' + str(k+1)][iStand]==n0: 
                    fcinv['I_SPECIES_CODE_' + str(k+1)][iStand]=n1
                
                #n0=meta['LUT']['FC_S']['S_SPECIES_CODE_' + str(k+1)][ListS[j][0]]
                #n1=meta['LUT']['FC_S']['S_SPECIES_CODE_' + str(k+1)][ListS[j][1]]        
                #if fcSd['S_SPECIES_CODE_' + str(k+1)][i]==n0: 
                #    fcSd['S_SPECIES_CODE_' + str(k+1)][i]=n1
    
    return meta,dmec,vri,fcinv

#%% CREATE BEST AVAILABLE INVENTORY

def CreateBestAvailableInventory(meta,vri,fcinv,flag_projects,idx,sxy):

    #--------------------------------------------------------------------------
    # Initialize best-available (gap-filled) inventory
    #--------------------------------------------------------------------------
    
    ba={} 
    ba['FIZ']=meta['LUT']['TIPSY']['FIZ']['I']*np.ones(meta['Project']['N Stand'])
    ba['LAND_COVER_CLASS_CD_1']=-999*np.ones(meta['Project']['N Stand'])
    ba['BEC_ZONE_CODE']=meta['LUT']['VRI']['BEC_ZONE_CODE']['SBS']*np.ones(meta['Project']['N Stand'])
    ba['Spc_CD1']=-999*np.ones(meta['Project']['N Stand'])
    ba['Spc_CD2']=-999*np.ones(meta['Project']['N Stand'])
    ba['Spc_CD3']=-999*np.ones(meta['Project']['N Stand'])
    ba['Spc_CD4']=-999*np.ones(meta['Project']['N Stand'])
    ba['Spc_CD5']=-999*np.ones(meta['Project']['N Stand'])
    ba['Spc_Pct1']=-999*np.ones(meta['Project']['N Stand'])
    ba['Spc_Pct2']=-999*np.ones(meta['Project']['N Stand'])
    ba['Spc_Pct3']=-999*np.ones(meta['Project']['N Stand'])
    ba['Spc_Pct4']=-999*np.ones(meta['Project']['N Stand'])
    ba['Spc_Pct5']=-999*np.ones(meta['Project']['N Stand'])
    ba['SI']=-999*np.ones(meta['Project']['N Stand'])
    ba['VRI_LIVE_STEMS_PER_HA']=-999*np.ones(meta['Project']['N Stand'])
    ba['PROJ_AGE_1']=-999*np.ones(meta['Project']['N Stand'])
    ba['Year']=-999*np.ones(meta['Project']['N Stand'])
    ba['LIVE_STAND_VOLUME_125']=-999*np.ones(meta['Project']['N Stand'])
    
    # Also keep track of the data source percentages
    basp={}
    basp['BEC_ZONE_CODE']={}
    basp['FIZ']={}
    basp['Spc_CD1']={}
    basp['SI']={}

    #--------------------------------------------------------------------------
    # Specify flag indicating whether subsetting occurs
    #--------------------------------------------------------------------------
    
    if np.isin('iKeep',list(meta['Project'].keys()))==True:
        # Tile project, only keeping a subset           
        flag_subset=1
    else:
        # Run all entries
        flag_subset=0
    
    #--------------------------------------------------------------------------
    # Best-available BEC zone and FIZ 
    #--------------------------------------------------------------------------
    
    # Fill with VRI
    N_tot=0
    for iStand0 in range(meta['Project']['N Stand']):
    
        if flag_subset==1:
            iStand=meta['Project']['iKeep'][iStand0]
        else:
            iStand=iStand0
        
        if idx['vri'][iStand]==None:
            continue
    
        ind0=idx['vri'][iStand]['Index'][0]
        N_tot=N_tot+ind0.size
        
        ba['LAND_COVER_CLASS_CD_1'][iStand0]=vri['LAND_COVER_CLASS_CD_1'][ind0] 
        
        ba['BEC_ZONE_CODE'][iStand0]=vri['BEC_ZONE_CODE'][ind0]        
        
        if (vri['BEC_ZONE_CODE'][ind0]==meta['LUT']['VRI']['BEC_ZONE_CODE']['CWH']) | \
            (vri['BEC_ZONE_CODE'][ind0]==meta['LUT']['VRI']['BEC_ZONE_CODE']['CDF']) | \
            (vri['BEC_ZONE_CODE'][ind0]==meta['LUT']['VRI']['BEC_ZONE_CODE']['MH']):
            ba['FIZ'][iStand0]=meta['LUT']['TIPSY']['FIZ']['C']
        else:
            ba['FIZ'][iStand0]=meta['LUT']['TIPSY']['FIZ']['I']
    
        ba['VRI_LIVE_STEMS_PER_HA'][iStand0]=vri['VRI_LIVE_STEMS_PER_HA'][ind0] 
        ba['PROJ_AGE_1'][iStand0]=vri['PROJ_AGE_1'][ind0] 
        ba['Year'][iStand0]=vri['Year'][ind0] 
        ba['LIVE_STAND_VOLUME_125'][iStand0]=vri['LIVE_STAND_VOLUME_125'][ind0]
    
    basp['BEC_ZONE_CODE']['From VRI']=N_tot/ba['SI'].size*100
    basp['FIZ']['From VRI']=N_tot/ba['SI'].size*100
    
    # Fill with global assumption
    ind=np.where(ba['BEC_ZONE_CODE']<=0)[0]
    if ind.size>0:
        ba['BEC_ZONE_CODE'][ind]=meta['LUT']['VRI']['BEC_ZONE_CODE']['SBS']
        ba['FIZ'][ind]=meta['LUT']['TIPSY']['FIZ']['I']
    
        basp['BEC_ZONE_CODE']['From global gap filling assumption']=ind.size/ba['SI'].size*100
        basp['FIZ']['From global gap filling assumption']=ind.size/ba['SI'].size*100
    
    #--------------------------------------------------------------------------
    # Best-available species composition (not planted)
    #--------------------------------------------------------------------------

    # Fill with forest cover
    N_tot=0
    for iStand0 in range(meta['Project']['N Stand']):
        
        if flag_subset==1:
            iStand=meta['Project']['iKeep'][iStand0]
        else:
            iStand=iStand0
        
        if idx['fcinv'][iStand]==None:
            continue
        
        ind0=idx['fcinv'][iStand]['Index'][0]
        
        for iSpc in range(5):
            id=fcinv['I_SPECIES_CODE_' + str(iSpc+1)][ind0]
            if id>=0:
                ba['Spc_CD' + str(iSpc+1)][iStand0]=id
                ba['Spc_Pct' + str(iSpc+1)][iStand0]=fcinv['I_SPECIES_PERCENT_' + str(iSpc+1)][ind0]
        
        N_tot=N_tot+1
    
    basp['Spc_CD1']['From FC inventory']=N_tot/ba['SI'].size*100
    
    # Fill with VRI
    N_tot=0
    for iStand0 in range(meta['Project']['N Stand']):  
        
        if flag_subset==1:
            iStand=meta['Project']['iKeep'][iStand0]
        else:
            iStand=iStand0
        
        if ba['Spc_CD1'][iStand0]>0:
            # Already populated with FC inventory layer
            continue
        
        if idx['vri'][iStand]==None:
            continue
        
        ind0=idx['vri'][iStand]['Index'][0]
        for iSpc in range(5):
            id=vri['SPECIES_CD_' + str(iSpc+1)][ind0]
            if id>=0:
                ba['Spc_CD' + str(iSpc+1)][iStand0]=id
                ba['Spc_Pct' + str(iSpc+1)][iStand0]=vri['SPECIES_PCT_' + str(iSpc+1)][ind0]
        N_tot=N_tot+1
    
    basp['Spc_CD1']['From VRI']=N_tot/ba['SI'].size*100

    # Fill with regional assumptions  
    N_tot=0
    for iStand0 in range(meta['Project']['N Stand']):
        
        if flag_subset==1:
            iStand=meta['Project']['iKeep'][iStand0]
        else:
            iStand=iStand0
        
        if ba['Spc_CD1'][iStand0]>0:
            # Already populated with FC inventory layer
            continue
        
        if ba['FIZ'][iStand0]==meta['LUT']['TIPSY']['FIZ']['C']:
            ba['Spc_CD1'][iStand0]=meta['LUT']['VRI']['SPECIES_CD_1']['FD']
            ba['Spc_Pct1'][iStand0]=100
        else:
            ba['Spc_CD1'][iStand0]=meta['LUT']['VRI']['SPECIES_CD_1']['PL']
            ba['Spc_Pct1'][iStand0]=100 
        
        N_tot=N_tot+1
    
    basp['Spc_CD1']['From regional assumptions']=N_tot/ba['SI'].size*100

    #--------------------------------------------------------------------------
    # Import site productivity layer    
    #--------------------------------------------------------------------------
    
    spl={}
    spl['SI_SPL']=-999*np.ones(meta['Project']['N Stand'])
    
    if 'xlim' in sxy:
        
        # Tiled project, we can just clip rasters to tile boundaries -> goes super fast
        
        z=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\SiteProductivityLayer\Site_Prod_Fd.tif')
        z=gis.ClipRaster(z,sxy['xlim'],sxy['ylim'])
        Site_Prod_Fd=z['Data'].flatten()
        del z
        garc.collect()
    
        z=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\SiteProductivityLayer\Site_Prod_Pl.tif')
        z=gis.ClipRaster(z,sxy['xlim'],sxy['ylim'])    
        Site_Prod_Pl=z['Data'].flatten()
        del z
        garc.collect()
    
        z=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\SiteProductivityLayer\Site_Prod_Sx.tif')
        z=gis.ClipRaster(z,sxy['xlim'],sxy['ylim'])    
        Site_Prod_Sx=z['Data'].flatten()
        del z
        garc.collect()
    
        # Populate dictionary with nearest estimate
        for iStand0 in range(meta['Project']['N Stand']):
        
            if flag_subset==1:
                iStand=meta['Project']['iKeep'][iStand0]
            else:
                iStand=iStand0
        
            if (ba['Spc_CD1'][iStand0]==meta['LUT']['VRI']['SPECIES_CD_1']['FD']) | \
                (ba['Spc_CD1'][iStand0]==meta['LUT']['VRI']['SPECIES_CD_1']['FDI']) | \
                (ba['Spc_CD1'][iStand0]==meta['LUT']['VRI']['SPECIES_CD_1']['FDC']):
                if Site_Prod_Fd[iStand]>0: 
                    spl['SI_SPL'][iStand0]=Site_Prod_Fd[iStand]
            elif (ba['Spc_CD1'][iStand0]==meta['LUT']['VRI']['SPECIES_CD_1']['PL']) | \
                (ba['Spc_CD1'][iStand0]==meta['LUT']['VRI']['SPECIES_CD_1']['PLI']) | \
                (ba['Spc_CD1'][iStand0]==meta['LUT']['VRI']['SPECIES_CD_1']['PLC']):
                if Site_Prod_Pl[iStand]>0: 
                    spl['SI_SPL'][iStand0]=Site_Prod_Pl[iStand]
            elif (ba['Spc_CD1'][iStand0]==meta['LUT']['VRI']['SPECIES_CD_1']['SX']):
                if Site_Prod_Sx[iStand]>0: 
                    spl['SI_SPL'][iStand0]=Site_Prod_Sx[iStand]
    
        del Site_Prod_Fd,Site_Prod_Pl,Site_Prod_Sx
        
    else:
        
        # Not a tiled project
        
        z=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\SiteProductivityLayer\Site_Prod_Fd.tif')
        x=z['X'][0,:]
        y=z['Y'][:,0]
        Site_Prod_Fd=np.squeeze(z['Data'])
        del z
        garc.collect()
    
        z=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\SiteProductivityLayer\Site_Prod_Pl.tif')
        Site_Prod_Pl=np.squeeze(z['Data'])
        del z
        garc.collect()
    
        z=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\SiteProductivityLayer\Site_Prod_Sx.tif')
        Site_Prod_Sx=np.squeeze(z['Data'])
        del z
        garc.collect()
    
        # Populate dictionary with nearest estimate
        for iStand0 in range(meta['Project']['N Stand']):
        
            if flag_subset==1:
                iStand=meta['Project']['iKeep'][iStand0]
            else:
                iStand=iStand0
        
            adx=np.abs(sxy['x'][iStand]-x)
            ady=np.abs(sxy['y'][iStand]-y)
            ix=np.where(adx==np.min(adx))[0]
            iy=np.where(ady==np.min(ady))[0]
            ind=np.ix_(iy,ix)

            if (ba['Spc_CD1'][iStand0]==meta['LUT']['VRI']['SPECIES_CD_1']['FD']) | \
                (ba['Spc_CD1'][iStand0]==meta['LUT']['VRI']['SPECIES_CD_1']['FDI']) | \
                (ba['Spc_CD1'][iStand0]==meta['LUT']['VRI']['SPECIES_CD_1']['FDC']):
                if Site_Prod_Fd[ind][0][0]>0: 
                    spl['SI_SPL'][iStand0]=Site_Prod_Fd[ind][0][0]
            elif (ba['Spc_CD1'][iStand0]==meta['LUT']['VRI']['SPECIES_CD_1']['PL']) | \
                (ba['Spc_CD1'][iStand0]==meta['LUT']['VRI']['SPECIES_CD_1']['PLI']) | \
                (ba['Spc_CD1'][iStand0]==meta['LUT']['VRI']['SPECIES_CD_1']['PLC']):
                if Site_Prod_Pl[ind][0][0]>0: 
                    spl['SI_SPL'][iStand0]=np.maximum(0,Site_Prod_Pl[ind][0][0])
            elif (ba['Spc_CD1'][iStand0]==meta['LUT']['VRI']['SPECIES_CD_1']['SX']):
                if Site_Prod_Sx[ind][0][0]>0: 
                    spl['SI_SPL'][iStand0]=np.maximum(0,Site_Prod_Sx[ind][0][0])
    
        del x,y,Site_Prod_Fd,Site_Prod_Pl,Site_Prod_Sx
    
    garc.collect()

    #--------------------------------------------------------------------------
    # Best-available site index
    #--------------------------------------------------------------------------
    
    # Populate with site productivity layer
    ind=np.where( (spl['SI_SPL']>5) )[0]
    ba['SI'][ind]=spl['SI_SPL'][ind]
    basp['SI']['From site productivity layer']=ind.size/ba['SI'].size*100

    # Populate with SI from forest cover inventory
    N_tot=0
    for iStand0 in range(meta['Project']['N Stand']):
        
        if flag_subset==1:
            iStand=meta['Project']['iKeep'][iStand0]
        else:
            iStand=iStand0
        
        if ba['SI'][iStand0]>0:
            # Already populated with Site Productivity Layer
            continue
        
        if idx['fcinv'][iStand]==None:
            # FC layer info unavailable
            continue    
        
        ind0=idx['fcinv'][iStand]['Index'][0]
        if fcinv['SITE_INDEX'][ind0]>0:
            # Populate with FC layer info
            ba['SI'][iStand0]=fcinv['SITE_INDEX'][ind0]
            N_tot=N_tot+ind0.size
    
    basp['SI']['From FC Inventory']=N_tot/ba['SI'].size*100

    # Where Productivity Layer and Forest Cover are missing, try populating with
    # the average value from each project
    if flag_projects==1:
        ind=np.where(ba['SI']<5)[0]
        u=np.unique(sxy['ID_atu_polygons'][ind])
        N_tot=0
        for i in range(u.size):
            iGood=np.where((sxy['ID_atu_polygons']==u[i]) & (ba['SI']>5))[0]
            iBad=np.where((sxy['ID_atu_polygons']==u[i]) & (ba['SI']<=5))[0]
            if (iGood.size!=0) & (iBad.size!=0):
                mu=np.mean(ba['SI'][iGood])
                #print(mu)
                ba['SI'][iBad]=mu
                N_tot=N_tot+iGood.size
        basp['SI']['From treatment area average']=N_tot/ba['SI'].size*100

    # Where there is no FC estimate and no SPL layer, polulate with VRI estimate
    N_tot=0
    for iStand0 in range(meta['Project']['N Stand']):
        
        if flag_subset==1:
            iStand=meta['Project']['iKeep'][iStand0]
        else:
            iStand=iStand0
        
        if ba['SI'][iStand0]>0:
            continue
        
        if idx['vri'][iStand]==None:
            continue
        
        ind0=idx['vri'][iStand]['Index'][0]
        if vri['SITE_INDEX'][ind0]>0:
            ba['SI'][iStand0]=vri['SITE_INDEX'][ind0]
            N_tot=N_tot+ind0.size
    
    basp['SI']['From VRI']=N_tot/ba['SI'].size*100

    # Where there is nothing, populate with regional averages
    
    iToFill1=np.where( (ba['SI']<0) & (ba['FIZ']==meta['LUT']['TIPSY']['FIZ']['I']) )[0]
    iGood=np.where( (ba['SI']>0) & (ba['FIZ']==meta['LUT']['TIPSY']['FIZ']['I']) )[0]
    ba['SI'][iToFill1]=np.mean(ba['SI'][iGood])
    
    iToFill2=np.where( (ba['SI']<0) & (ba['FIZ']==meta['LUT']['TIPSY']['FIZ']['C']) )[0]
    iGood=np.where( (ba['SI']>0) & (ba['FIZ']==meta['LUT']['TIPSY']['FIZ']['C']) )[0]
    ba['SI'][iToFill2]=np.mean(ba['SI'][iGood])
    
    basp['SI']['From regional averages']=(iToFill1.size+iToFill2.size)/ba['SI'].size*100

    # Where there are still missing values assume mean of whole dataset
    ind=np.where(ba['SI']<0)[0]
    if ind.size>0:
        iGood=np.where( (ba['SI']>0) )[0]
        ba['SI'][ind]=np.mean(ba['SI'][iGood])
        basp['SI']['From global average']=ind.size/ba['SI'].size*100

    # Convert to integer (for BatchTIPSY)
    ba['SI']=np.round(ba['SI']).astype(int)

    return ba,basp

#%% Get unique growth curves
    
def ExtractUniqueGrowthCurves(meta,gc,GC_Variable_List):

    ugc={}
    ugc['GC_Variable_List']=np.array(GC_Variable_List)[3:]

    # Calculate unique stand types
    ugc['Full']=np.zeros((int(4e6),len(GC_Variable_List)))

    cnt=0
    for iStand in range(meta['Project']['N Stand']):
        for iScn in range(meta['Project']['N Scenario']):
            gc0=gc[iStand][iScn]
            for iGC in range(gc0['ID_GC'].size):
                for k in range(len(GC_Variable_List)):
                    key=GC_Variable_List[k]
                    ugc['Full'][cnt,k]=gc0[key][iGC]
                cnt=cnt+1
    ugc['Full']=ugc['Full'][0:cnt,:]

    # Unique growth curves
    # The 'Inverse' variable acts as the crosswalk between the full and unique gc arrays
    ugc['Unique'],ugc['Index'],ugc['Inverse']=np.unique(ugc['Full'][:,3:],return_index=True,return_inverse=True,axis=0)

    return ugc

#%% Adjust species-specific mortality

def AdjustSpeciesSpecificMortality(meta,dmec,gc,iB):

    # Species affected sets
    Pest_List=['IBM','IBB','IBS','IBD','IDW']
    
    SA_List=[None]*len(Pest_List)
    for iPest in range(len(Pest_List)):
        ind=np.where( (meta['Param']['BE']['DistBySC']['Name']==Pest_List[iPest]) )[0][0]
        SA_List[iPest]=np.array([meta['Param']['BE']['DistBySC']['SpcCD1'][ind],meta['Param']['BE']['DistBySC']['SpcCD2'][ind],
                 meta['Param']['BE']['DistBySC']['SpcCD3'][ind],meta['Param']['BE']['DistBySC']['SpcCD4'][ind],
                 meta['Param']['BE']['DistBySC']['SpcCD5'][ind],meta['Param']['BE']['DistBySC']['SpcCD6'][ind]])

    for iStand in range(meta['Project']['N Stand']):      
    
        for iYr in range(dmec[iStand]['Year'].size):
        
            for iPest in range(len(Pest_List)):
        
                if dmec[iStand]['ID_Type'][iYr]==meta['LUT']['Dist'][Pest_List[iPest]]:
            
                    ind_GC=int(dmec[iStand]['ID_GC'][iB][iYr]-1)
            
                    scd=[None]*4
                    scd[0]=cbu.lut_n2s(meta['LUT']['VRI']['SPECIES_CD_1'],gc[iStand][iB]['s1'][ind_GC])
                    scd[1]=cbu.lut_n2s(meta['LUT']['VRI']['SPECIES_CD_1'],gc[iStand][iB]['s2'][ind_GC])
                    scd[2]=cbu.lut_n2s(meta['LUT']['VRI']['SPECIES_CD_1'],gc[iStand][iB]['s3'][ind_GC])
                    scd[3]=cbu.lut_n2s(meta['LUT']['VRI']['SPECIES_CD_1'],gc[iStand][iB]['s4'][ind_GC])
            
                    spct=[None]*4
                    spct[0]=gc[iStand][iB]['p1'][ind_GC]
                    spct[1]=gc[iStand][iB]['p2'][ind_GC]
                    spct[2]=gc[iStand][iB]['p3'][ind_GC]
                    spct[3]=gc[iStand][iB]['p4'][ind_GC]            
            
                    PercentAffected=0
                    for i in range(4):
                        if np.isin(scd[i],SA_List[iPest])==True:
                            PercentAffected=PercentAffected+spct[i]
                    dmec[iStand]['MortalityFactor'][iYr]=(PercentAffected/100)*dmec[iStand]['MortalityFactor'][iYr]        

    return dmec

#%% Put DMEC events in order

def PutEventsInOrder(dmec,meta):    
    for iStand in range(meta['Project']['N Stand']):
        d=dmec[iStand].copy()
        ord=np.argsort(d['Year'])
        for key in d.keys():
            d[key]=d[key][ord]
        dmec[iStand]=d.copy()
    return dmec

#%% Export AT Layer data to spreadhseet
## OLD    
#def ExportSummaryByGridCell(meta,atu_multipolygons,dAdmin,sxy,atu,fcinv,vri,pl,op,include_planting,project_name):
#    
#    #--------------------------------------------------------------------------
#    # ATU
#    #--------------------------------------------------------------------------
#    
#    d={}
#    d['IdxToSXY']=atu['IdxToSXY'].copy()
#    d['ID_Multipolygon']=atu['IdxToSXY'].copy()
#    d['Year']=atu['Year'].copy()
#    d['Month']=atu['Month'].copy()
#    d['OPENING_ID']=atu['OPENING_ID'].copy()
#    d['Activity_Type']=np.array(['empty' for _ in range(atu['IdxToSXY'].size)],dtype=object)  
#    
#    d['AEF_ATU']=np.zeros(atu['IdxToSXY'].size)    
#    for i in range(len(atu_multipolygons)):
#        ind1=np.where(sxy['ID_atu_multipolygons']==i)[0]
#        nxy=ind1.size
#        A=atu_multipolygons[i]['ACTUAL_TREATMENT_AREA']
#        if A==None:
#            A=0.00001
#        for j in range(ind1.size):
#            ind2=np.where(atu['IdxToSXY']==ind1[j])[0]
#            d['AEF_ATU'][ind2]=np.round(A/nxy,3)
#            d['ID_Multipolygon'][ind2]=i    
#        
#    d['FIA_PROJECT_ID']=np.array(['empty' for _ in range(atu['IdxToSXY'].size)],dtype=object)    
#    u=np.unique(atu['FIA_PROJECT_ID'])
#    for i in range(u.size):
#        ind=np.where(atu['FIA_PROJECT_ID']==u[i])[0]
#        d['FIA_PROJECT_ID'][ind]=cbu.lut_n2s(meta['LUT']['ATU']['FIA_PROJECT_ID'],u[i])
#        
#    d['FSC']=np.array(['empty' for _ in range(atu['IdxToSXY'].size)],dtype=object)    
#    for i in range(atu['Year'].size):
#        d['FSC'][i]=cbu.lut_n2s(meta['LUT']['ATU']['SILV_FUND_SOURCE_CODE'],atu['SILV_FUND_SOURCE_CODE'][i])[0]    
#    d['DistCD']=np.array(['empty' for _ in range(atu['IdxToSXY'].size)],dtype=object)    
#    for i in range(atu['Year'].size):
#        d['DistCD'][i]=cbu.lut_n2s(meta['LUT']['ATU']['DISTURBANCE_CODE'],atu['DISTURBANCE_CODE'][i])[0]    
#    d['SBC']=np.array(['empty' for _ in range(atu['IdxToSXY'].size)],dtype=object)
#    for i in range(atu['Year'].size):
#        d['SBC'][i]=cbu.lut_n2s(meta['LUT']['ATU']['SILV_BASE_CODE'],atu['SILV_BASE_CODE'][i])[0]
#    d['SMC']=np.array(['empty' for _ in range(atu['IdxToSXY'].size)],dtype=object)
#    for i in range(atu['Year'].size):
#        d['SMC'][i]=cbu.lut_n2s(meta['LUT']['ATU']['SILV_METHOD_CODE'],atu['SILV_METHOD_CODE'][i])[0]
#    d['STC']=np.array(['empty' for _ in range(atu['IdxToSXY'].size)],dtype=object)
#    for i in range(atu['Year'].size):
#        d['STC'][i]=cbu.lut_n2s(meta['LUT']['ATU']['SILV_TECHNIQUE_CODE'],atu['SILV_TECHNIQUE_CODE'][i])[0]
#    d['SilvObjectiveCode1']=np.array(['empty' for _ in range(atu['IdxToSXY'].size)],dtype=object)
#    for i in range(atu['Year'].size):
#        d['SilvObjectiveCode1'][i]=cbu.lut_n2s(meta['LUT']['ATU']['SILV_OBJECTIVE_CODE_1'],atu['SILV_OBJECTIVE_CODE_1'][i])[0]    
#    
#    # Add VRI
#    d['BGCz']=np.array(['empty' for _ in range(atu['IdxToSXY'].size)],dtype=object)
#    d['BGCsz']=np.array(['empty' for _ in range(atu['IdxToSXY'].size)],dtype=object)
#    d['BGCv']=np.array(['empty' for _ in range(atu['IdxToSXY'].size)],dtype=object)
#    for i in range(d['IdxToSXY'].size):
#        ind=np.where(vri['IdxToSXY']==d['IdxToSXY'][i])[0]
#        if ind.size==0:
#            continue
#        d['BGCz'][i]=cbu.lut_n2s(meta['LUT']['VRI']['BEC_ZONE_CODE'],vri['BEC_ZONE_CODE'][ind[0]])[0]
#        d['BGCsz'][i]=cbu.lut_n2s(meta['LUT']['VRI']['BEC_SUBZONE'],vri['BEC_SUBZONE'][ind[0]])[0]
#        d['BGCv'][i]=cbu.lut_n2s(meta['LUT']['VRI']['BEC_VARIANT'],vri['BEC_VARIANT'][ind[0]])[0]    
#    
#    # Add opening variables
#    d['District']=np.array(['empty' for _ in range(atu['IdxToSXY'].size)],dtype=object)
#    for i in range(d['IdxToSXY'].size):
#        ind=np.where(op['IdxToSXY']==d['IdxToSXY'][i])[0]
#        if ind.size==0:
#            continue
#        d['District'][i]=cbu.lut_n2s(meta['LUT']['OP']['DISTRICT_NAME'],op['DISTRICT_NAME'][ind[0]])[0]
#    
#    # Add activity type
#    if project_name=='FCI':
#        
#        u=np.unique(d['FIA_PROJECT_ID'])
#        for i in range(u.size):
#            ind1=np.where(d['FIA_PROJECT_ID']==u[i])[0]
#            if ind1.size==0:
#                continue
#            ind2=np.where( (dAdmin['PP Number']==u[i]) & (dAdmin['Removed']!='Removed') )[0]
#            if ind2.size==0:
#                continue
#            u2=np.unique(d['IdxToSXY'][ind1])
#            for j in range(u2.size):
#                ind3=np.where(d['IdxToSXY']==u2[j])[0]
#                for k in range(ind3.size):
#                    d['Activity_Type'][ind3[k]]=dAdmin['Activity Type'][ind2[0]]
#                    
#        u=np.unique(d['OPENING_ID'])
#        for i in range(u.size):
#            ind1=np.where(d['OPENING_ID']==u[i])[0]
#            if ind1.size==0:
#                continue
#            ind2=np.where( (dAdmin['OPENING_ID']==u[i]) & (dAdmin['Removed']!='Removed') )[0]
#            if ind2.size==0:
#                continue
#            u2=np.unique(d['IdxToSXY'][ind1])
#            for j in range(u2.size):
#                ind3=np.where(d['IdxToSXY']==u2[j])[0]
#                for k in range(ind3.size):
#                    d['Activity_Type'][ind3[k]]=dAdmin['Activity Type'][ind2[0]]            
#    
#    elif project_name=='ReforestationNonOb':
#    
#        nam=['No Planting','SL','KD','UNDER','NSR Backlog','Unclassified']
#        for i in range(d['IdxToSXY'].size):
#            d['Activity_Type'][i]=nam[int(meta['ProjectType'][int(d['IdxToSXY'][i])])]
#    
#    # Add Planting
#    
#    if include_planting=='On':
#    
#        d['Pl_SPH']=np.round(atu['ACTUAL_PLANTED_NUMBER']/atu['ACTUAL_TREATMENT_AREA'])
#        ind=np.where( (d['SBC']!='PL') & (d['STC']!='PL') )[0]
#        d['Pl_SPH'][ind]=0
#    
#        # Add planting info
#        num_of_spc=10
#        for i in range(num_of_spc):
#            d['Pl_Spc' + str(i+1) + '_CD']=np.array([' ' for _ in range(atu['IdxToSXY'].size)],dtype=object)
#            d['Pl_Spc' + str(i+1) + '_Pct']=np.array([' ' for _ in range(atu['IdxToSXY'].size)],dtype=object)
#            #d['Pl_Spc' + str(i+1) + '_NumTree']=np.array([' ' for _ in range(atu['IdxToSXY'].size)],dtype=object)
#            d['Pl_Spc' + str(i+1) + '_SeedLot']=np.array([' ' for _ in range(atu['IdxToSXY'].size)],dtype=object)
#    
#        for i in range(atu['IdxToSXY'].size):
#            if (d['SBC'][i]!='PL') & (d['STC'][i]!='PL'):
#                continue        
#            ind=np.where( (pl['IdxToSXY']==d['IdxToSXY'][i]) & (pl['OPENING_ID']==d['OPENING_ID'][i]) & (pl['Year']==d['Year'][i]) )[0]
#            tot_pl=np.sum(pl['NUMBER_PLANTED'][ind])
#            if ind.size>0:
#                Ord=np.flip(np.argsort(pl['NUMBER_PLANTED'][ind]))
#                for j in range(ind.size):
#                    if j>num_of_spc-1:
#                        continue
#                    ind0=ind[Ord[j]]
#                    d['Pl_Spc' + str(j+1) + '_CD'][i]=cbu.lut_n2s(meta['LUT']['PL']['SILV_TREE_SPECIES_CODE'],pl['SILV_TREE_SPECIES_CODE'][ind0])[0]
#                    d['Pl_Spc' + str(j+1) + '_Pct'][i]=np.round(pl['NUMBER_PLANTED'][ind0]/tot_pl*100)
#                    #d['Pl_Spc' + str(j+1) + '_NumTree'][i]=pl['NUMBER_PLANTED'][ind0]
#                    d['Pl_Spc' + str(j+1) + '_SeedLot'][i]=pl['SEEDLOT_NUMBER'][ind0]
#    
#    # Convert to dataframe
#    df_atu=pd.DataFrame.from_dict(d)
#    
#    #--------------------------------------------------------------------------
#    # FC inventory
#    #--------------------------------------------------------------------------
#    
#    d={}
#    d['IdxToSXY']=fcinv['IdxToSXY'].copy()
#    d['ID_Multipolygon']=np.zeros(fcinv['IdxToSXY'].size)    
#    for i in range(len(atu_multipolygons)):
#        ind1=np.where(sxy['ID_atu_multipolygons']==i)[0]
#        for j in range(ind1.size):
#            ind2=np.where(fcinv['IdxToSXY']==ind1[j])[0]
#            d['ID_Multipolygon'][ind2]=i 
#    d['Year']=fcinv['REFERENCE_YEAR'].copy()
#    d['Activity_Type']=np.array(['empty' for _ in range(fcinv['IdxToSXY'].size)],dtype=object)  
#    d['SI']=fcinv['SITE_INDEX'].copy()
#    d['I_SPH']=fcinv['I_TOTAL_STEMS_PER_HA'].copy()
#    d['I_Spc1_CD']=np.array(['empty' for _ in range(fcinv['IdxToSXY'].size)],dtype=object)
#    d['I_Spc2_CD']=np.array(['empty' for _ in range(fcinv['IdxToSXY'].size)],dtype=object)
#    for i in range(fcinv['IdxToSXY'].size):
#        d['I_Spc1_CD'][i]=cbu.lut_n2s(meta['LUT']['VRI']['SPECIES_CD_1'],fcinv['I_SPECIES_CODE_1'][i])[0]
#        try:
#            d['I_Spc2_CD'][i]=cbu.lut_n2s(meta['LUT']['VRI']['SPECIES_CD_1'],fcinv['I_SPECIES_CODE_2'][i])[0]
#        except:
#            pass
#
#    # Add VRI
#    d['BGCz']=np.array(['empty' for _ in range(fcinv['IdxToSXY'].size)],dtype=object)
#    d['BGCsz']=np.array(['empty' for _ in range(fcinv['IdxToSXY'].size)],dtype=object)
#    d['BGCv']=np.array(['empty' for _ in range(fcinv['IdxToSXY'].size)],dtype=object)
#    for i in range(d['IdxToSXY'].size):
#        ind=np.where(vri['IdxToSXY']==d['IdxToSXY'][i])[0]
#        if ind.size==0:
#            continue
#        d['BGCz'][i]=cbu.lut_n2s(meta['LUT']['VRI']['BEC_ZONE_CODE'],vri['BEC_ZONE_CODE'][ind[0]])[0]
#        d['BGCsz'][i]=cbu.lut_n2s(meta['LUT']['VRI']['BEC_SUBZONE'],vri['BEC_SUBZONE'][ind[0]])[0]
#        d['BGCv'][i]=cbu.lut_n2s(meta['LUT']['VRI']['BEC_VARIANT'],vri['BEC_VARIANT'][ind[0]])[0]  
#
#    # Add opening variables
#    d['District']=np.array(['empty' for _ in range(fcinv['IdxToSXY'].size)],dtype=object)
#    for i in range(d['IdxToSXY'].size):
#        ind=np.where(op['IdxToSXY']==d['IdxToSXY'][i])[0]
#        if ind.size==0:
#            continue
#        #Natural Resource District
#        d['District'][i]=cbu.lut_n2s(meta['LUT']['OP']['DISTRICT_NAME'],op['DISTRICT_NAME'][ind[0]])[0]
#
#    df_fcinv=pd.DataFrame.from_dict(d)
#    
#    #--------------------------------------------------------------------------
#    # Merge ATU and FCI
#    #--------------------------------------------------------------------------
#    
#    df=df_atu.merge(df_fcinv,how='outer',on=('IdxToSXY','Year','ID_Multipolygon','Activity_Type','BGCz','BGCsz','BGCv','District'))
#    
#    #--------------------------------------------------------------------------
#    # Save
#    #--------------------------------------------------------------------------
#    
#    df=df.sort_values(by=['IdxToSXY','Year','Month'])
#    
#    df.to_excel(meta['Paths']['Project'] + '\\Inputs\\SummaryAttributesBySXY.xlsx',index=False)   
#
#    return

#%% Timber harvesting land base 

def DefineTHLB(meta,ba,dmec,fcres,lul,ogmal,park):
    
    meta['Project']['THLB']={}
    
    #------------------------------------------------------------------------------
    # Fix for tiled
    #------------------------------------------------------------------------------
    
    if 'IdxToSXY' not in park:
        fcres['IdxToSXY']=fcres['IdxToGrid']
        lul['IdxToSXY']=lul['IdxToGrid']
        ogmal['IdxToSXY']=ogmal['IdxToGrid']
        park['IdxToSXY']=park['IdxToGrid']
    
    #------------------------------------------------------------------------------
    # Initialize THLB flags (THLB=1,Non-THLB=0)
    #------------------------------------------------------------------------------

    # Initially assume everything is in the THLB
    meta['Project']['THLB']['Actual']=np.ones((meta['Project']['N Time'],meta['Project']['N Stand']))

    # Define second THLB with areas that have been removed due to value diversification
    meta['Project']['THLB']['Baseline']=meta['Project']['THLB']['Actual'].copy()

    # Index to stands that are uneconomic
    iUneconomic=np.where(ba['SI']<=5)[0]

    # Remove uneconomic stands from THLB
    meta['Project']['THLB']['Actual'][:,iUneconomic]=0
    meta['Project']['THLB']['Baseline'][:,iUneconomic]=0

    # Idenify stands that have been harvested
    has_been_harvested=np.zeros(meta['Project']['N Stand'])
    for i in range(len(dmec)):
        ind=np.where(dmec[i]['ID_Type']==meta['LUT']['Dist']['Harvest'])[0]
        if ind.size>0:
            has_been_harvested[i]=1

    # Index to stands that have not been harvested
    iNoHarv=np.where( (has_been_harvested==0) & (ba['SI']>5) )[0]

    # Use the ratio of THLB to non-THLB as an indicator of what will be harvested
    # among remaining primary forest
    ratio_thlb=22/55 # ratio of THLB to total forest (SOF)

    corr=iUneconomic.size/meta['Project']['N Stand']

    # Probability of evading harvest
    if iNoHarv.size>0:
        p_evade=(1-ratio_thlb-corr)*(meta['Project']['N Stand']/iNoHarv.size)
    else:
        p_evade=(1-ratio_thlb-corr)

    # Random prediction of whether it will evade harvesting
    iRem=np.where(np.random.random(iNoHarv.size)<p_evade)[0]
    meta['Project']['THLB']['Actual'][:,iNoHarv[iRem]]=0
    meta['Project']['THLB']['Baseline'][:,iNoHarv[iRem]]=0

    # np.sum(meta['Project']['THLB']['Actual'][0,:])/meta['Project']['N Stand']

    #------------------------------------------------------------------------------
    # Define the year of transition from THLB to non-THLB for specific LU types
    #------------------------------------------------------------------------------

    # Look at abundance of each LUL type
    d=meta['LUT']['LU L']['LEGAL_FEAT_OBJECTIVE'].copy()
    for k in meta['LUT']['LU L']['LEGAL_FEAT_OBJECTIVE'].keys():
        ind=np.where( (lul['LEGAL_FEAT_OBJECTIVE']==meta['LUT']['LU L']['LEGAL_FEAT_OBJECTIVE'][k]) )[0]
        d[k]=ind.size
    #ds={k: v for k,v in sorted(d.items(), key=lambda item: item[1])}

    # List of legal land use plan types that transition to conservation (this needs to be vetted by experts)
    ListCon=['Visually Sensitive Areas','Connectivity Corridors','Scenic Areas','Caribou Winter Habitat Zones',
      'Tourism Areas','Visual Quality','High Value Grizzly Bear Habitat','No Timber Harvesting Areas','Landscape Corridors',
      'Critical Deer Winter Range','Sensitive Watershed','Water Management Units','High Value Wetlands for Moose','Telkwa Caribou Recovery Area',
      'Caribou Migration Corridor','High Biodiversity Emphasis Areas','Scenic Corridors']

    # Calculate area that will transition to non-THLB
    A=0
    
    # Add regional land use plans    
    for i in range(len(ListCon)):
        ind=np.where( (lul['LEGAL_FEAT_OBJECTIVE']==meta['LUT']['LU L']['LEGAL_FEAT_OBJECTIVE'][ListCon[i]]) )[0]
        A=A+ind.size
    
    # Add parks
    A=A+park['IdxToSXY'].size
    
    # Add legal OGMAs
    A=A+ogmal['IdxToSXY'].size
    
    #print(A)
    print(str(A/meta['Project']['N Stand']*100) + '% of the sample transitioned to non-THLB.')

    # Look at reserves
    N={}
    for k in meta['LUT']['FC_R']['SILV_RESERVE_CODE'].keys():
        ind=np.where(fcres['SILV_RESERVE_CODE']==meta['LUT']['FC_R']['SILV_RESERVE_CODE'][k])[0]
        N[k]=ind.size

    # Initialize year of transition
    thlb_YearTransitionOut=np.zeros(meta['Project']['N Stand'])

    # Define year of transition from legal land use plan objectives
    for i in range(len(ListCon)):
        ind=np.where( (lul['LEGAL_FEAT_OBJECTIVE']==meta['LUT']['LU L']['LEGAL_FEAT_OBJECTIVE'][ListCon[i]]) )[0]
        thlb_YearTransitionOut[lul['IdxToSXY'][ind]]=2010

    # Define year of transition from parks
    thlb_YearTransitionOut[park['IdxToSXY']]=1995

    # Define year of transition from legal old-growth OGMAs
    thlb_YearTransitionOut[ogmal['IdxToSXY']]=2010

    # Apply transition to actual THLB
    for j in range(thlb_YearTransitionOut.size):
        if thlb_YearTransitionOut[j]>0:
            it=np.where( (meta['Year']>=thlb_YearTransitionOut[j]) )[0]
            meta['Project']['THLB']['Actual'][it,j]=0

    # Adjust the baseline so that simulated harvesting between 1995 and 2020 only
    # occurs in areas where the THLB was affected by value diversification
    for year in range(1990,2021,1):
        iT=np.where(meta['Year']==year)[0]
        iS=np.where( (meta['Project']['THLB']['Baseline'][iT,:]==1) & (meta['Project']['THLB']['Actual'][iT,:]==1) )[1]
        meta['Project']['THLB']['Baseline'][iT,iS]=0

    flg=0
    if flg==1:
        plt.figure(2)
        plt.plot(np.sum(meta['Project']['THLB']['Actual'],axis=1)/meta['Project']['N Stand'])
        plt.plot(np.sum(meta['Project']['THLB']['Baseline'],axis=1)/meta['Project']['N Stand'],'--')
        
    return meta

#%% Load sparse geospatiatial inputs

def LoadSparseGeospatialInputs(meta):
    
    try:        
        sxy=gu.ipickle(meta['Paths']['Geospatial'] + '\\sxy.pkl')
    except:
        # Make this work for tiled projects as well
        sxy=[]    
    atu=gu.ipickle(meta['Paths']['Geospatial'] + '\\RSLT_ACTIVITY_TREATMENT_SVW.pkl')
    op=gu.ipickle(meta['Paths']['Geospatial'] + '\\RSLT_OPENING_SVW.pkl')
    burnsev=gu.ipickle(meta['Paths']['Geospatial'] + '\\VEG_BURN_SEVERITY_SP.pkl')
    vri=gu.ipickle(meta['Paths']['Geospatial'] + '\\VEG_COMP_LYR_R1_POLY.pkl')
    fcinv=gu.ipickle(meta['Paths']['Geospatial'] + '\\RSLT_FOREST_COVER_INV_SVW.pkl')
    fcsilv=gu.ipickle(meta['Paths']['Geospatial'] + '\\RSLT_FOREST_COVER_SILV_SVW.pkl')
    fcres=gu.ipickle(meta['Paths']['Geospatial'] + '\\RSLT_FOREST_COVER_RESERVE_SVW.pkl')
    pl=gu.ipickle(meta['Paths']['Geospatial'] + '\\RSLT_PLANTING_SVW.pkl')
    fire=gu.ipickle(meta['Paths']['Geospatial'] + '\\PROT_HISTORICAL_FIRE_POLYS_SP.pkl')
    pest=gu.ipickle(meta['Paths']['Geospatial'] + '\\PEST_INFESTATION_POLY.pkl')
    cut=gu.ipickle(meta['Paths']['Geospatial'] + '\\VEG_CONSOLIDATED_CUT_BLOCKS_SP.pkl')
    lul=gu.ipickle(meta['Paths']['Geospatial'] + '\\RMP_PLAN_LEGAL_POLY_SVW.pkl')
    park=gu.ipickle(meta['Paths']['Geospatial'] + '\\TA_PARK_ECORES_PA_SVW.pkl')
    ogmal=gu.ipickle(meta['Paths']['Geospatial'] + '\\RMP_OGMA_LEGAL_CURRENT_SVW.pkl')

    return sxy,atu,op,burnsev,vri,fcinv,fcsilv,fcres,pl,fire,pest,cut,lul,park,ogmal

#%% Load index to sparse grid
    
def LoadSparseGridIndex(meta):
    
    idx={}
    idx['atu']=gu.ipickle(meta['Paths']['Geospatial'] + '\\RSLT_ACTIVITY_TREATMENT_SVW_IdxToInv.pkl')
    idx['pl']=gu.ipickle(meta['Paths']['Geospatial'] + '\\RSLT_PLANTING_SVW_IdxToInv.pkl')
    idx['op']=gu.ipickle(meta['Paths']['Geospatial'] + '\\RSLT_OPENING_SVW_IdxToInv.pkl')
    idx['burnsev']=gu.ipickle(meta['Paths']['Geospatial'] + '\\VEG_BURN_SEVERITY_SP_IdxToInv.pkl')
    idx['vri']=gu.ipickle(meta['Paths']['Geospatial'] + '\\VEG_COMP_LYR_R1_POLY_IdxToInv.pkl')
    idx['fcinv']=gu.ipickle(meta['Paths']['Geospatial'] + '\\RSLT_FOREST_COVER_INV_SVW_IdxToInv.pkl')
    idx['fire']=gu.ipickle(meta['Paths']['Geospatial'] + '\\PROT_HISTORICAL_FIRE_POLYS_SP_IdxToInv.pkl')
    idx['pest']=gu.ipickle(meta['Paths']['Geospatial'] + '\\PEST_INFESTATION_POLY_IdxToInv.pkl')
    idx['cut']=gu.ipickle(meta['Paths']['Geospatial'] + '\\VEG_CONSOLIDATED_CUT_BLOCKS_SP_IdxToInv.pkl')
    idx['lul']=gu.ipickle(meta['Paths']['Geospatial'] + '\\RMP_PLAN_NON_LEGAL_POLY_SVW_IdxToInv.pkl')
    idx['park']=gu.ipickle(meta['Paths']['Geospatial'] + '\\TA_PARK_ECORES_PA_SVW_IdxToInv.pkl')
    idx['ogmal']=gu.ipickle(meta['Paths']['Geospatial'] + '\\RMP_OGMA_LEGAL_CURRENT_SVW_IdxToInv.pkl')
    
    return idx

#%% Export summary by multipolygon

def ExportSummaryActivities_ByMP(meta,par,atu_multipolygons,uMP,sxy,atu,op,burnsev,vri,fcinv,fcsilv,fcres,pl,fire,pest,cut,lul,park,ogmal,dmec,ba,clm):

    # Get land cover class names
    lcc=gu.ReadExcel(r'C:\Users\rhember\Documents\Code_Python\fcgadgets\cbrunner\Parameters\Parameters_VRI.xlsx','LCC')
    
    # Initialize structure
    dMP={}

    # Add ID_Multipolygon
    dMP['ID_Multipolygon']=np.arange(0,len(atu_multipolygons))

    dMP['See Overlapping MPs']=np.zeros(dMP['ID_Multipolygon'].size)
    for i in range(len(atu_multipolygons)):
        ind=np.where(sxy['ID_atu_multipolygons']==i)[0]
        if ind.size==0:
            dMP['See Overlapping MPs'][i]=1

    # Initialize first set of variables from atu_multipolygon
    vr=['FIA_PROJECT_ID','OPENING_ID','Year','ACTUAL_TREATMENT_AREA','ACTUAL_PLANTED_NUMBER','SILV_BASE_CODE','SILV_TECHNIQUE_CODE',
        'SILV_METHOD_CODE', 'SILV_OBJECTIVE_CODE_1','GeomFromOpLyr', 'GeomFromFcLyr','ACTUAL_PLANTED_NUMBER']

    for k in vr:
        for i in range(1000):
            if atu_multipolygons[i][k]!=None:
                break
        if (type(atu_multipolygons[i][k])==np.int32) | (type(atu_multipolygons[i][k])==int) | (type(atu_multipolygons[i][k])==np.float64):
            dMP[k]=-999*np.ones(dMP['ID_Multipolygon'].size)
        else:
            dMP[k]=np.array(['' for _ in range(dMP['ID_Multipolygon'].size)],dtype=object)

    # Populate with contents of atu_multipolygon
    for i in range(len(atu_multipolygons)):
        for k in vr:
            if atu_multipolygons[i][k]==None:
                if (type(dMP[k][i])==np.int32) | (type(dMP[k][i])==int) | (type(dMP[k][i])==np.float64):
                    dMP[k][i]=np.nan
                else:
                    dMP[k][i]=''
            else:
                dMP[k][i]=atu_multipolygons[i][k]

    # Convert actual planted number to planting density
    dMP['ACTUAL_PLANTED_NUMBER']=dMP['ACTUAL_PLANTED_NUMBER']/dMP['ACTUAL_TREATMENT_AREA']
    for i in range(dMP['ACTUAL_PLANTED_NUMBER'].size):
        dMP['ACTUAL_PLANTED_NUMBER'][i]=np.round(dMP['ACTUAL_PLANTED_NUMBER'][i])
    #ind=np.where(dMP['ACTUAL_PLANTED_NUMBER']>0)[0]
    #dMP['ACTUAL_PLANTED_NUMBER'][ind]=int(dMP['ACTUAL_PLANTED_NUMBER'][ind])

    # Initiatlize Project Type
    dMP['Project Type']=np.array(['' for _ in range(dMP['ID_Multipolygon'].size)],dtype=object)
    for i in range(len(atu_multipolygons)):
        dMP['Project Type'][i]=meta['Project']['ProjectTypeByMP'][i]

    # Add land cover class
    dMP['LCC 1']=np.array(['' for _ in range(dMP['OPENING_ID'].size)],dtype=object)
    dMP['LCC 2']=np.array(['' for _ in range(dMP['OPENING_ID'].size)],dtype=object)
    for i in range(len(atu_multipolygons)):
        ind=np.where(sxy['ID_atu_multipolygons']==i)[0]
        if ind.size==0:
            continue
        cd=[]
        for j in range(ind.size):
            cd.append(cbu.lut_n2s(meta['LUT']['VRI']['LAND_COVER_CLASS_CD_1'],ba['LAND_COVER_CLASS_CD_1'][ind[j]])[0])
    
        nam=['']*len(cd)
        for j in range(len(cd)):
            ind=np.where(lcc['Code']==cd[j])[0]
            if ind.size>0:
                nam[j]=lcc['Name'][ind[0]]
    
        u,c=np.unique(nam,return_counts=True)
        ord=np.flip(np.argsort(c))
        for j in range(np.minimum(2,u.size)):
            dMP['LCC ' + str(j+1)][i]=nam[ord[j]] 
    
    # Add BGC zone
    dMP['BGC Zone 1']=np.array(['' for _ in range(dMP['OPENING_ID'].size)],dtype=object)
    dMP['BGC Zone 2']=np.array(['' for _ in range(dMP['OPENING_ID'].size)],dtype=object)
    dMP['BGC Zone 3']=np.array(['' for _ in range(dMP['OPENING_ID'].size)],dtype=object)
    for i in range(len(atu_multipolygons)):
        ind=np.where(sxy['ID_atu_multipolygons']==i)[0]
        if ind.size==0:
            continue
        cd=[]
        for j in range(ind.size):
            cd.append(cbu.lut_n2s(meta['LUT']['VRI']['BEC_ZONE_CODE'],ba['BEC_ZONE_CODE'][ind[j]])[0])
    
        u,c=np.unique(cd,return_counts=True)
        ord=np.flip(np.argsort(c))
        for j in range(np.minimum(3,u.size)):
            dMP['BGC Zone ' + str(j+1)][i]=cd[ord[j]]

    # Add SI from best-available
    dMP['SI_ba_mean']=-999*np.ones(dMP['ID_Multipolygon'].size)
    dMP['SI_ba_min']=-999*np.ones(dMP['ID_Multipolygon'].size)
    dMP['SI_ba_max']=-999*np.ones(dMP['ID_Multipolygon'].size)
    for i in range(len(atu_multipolygons)):
        ind=np.where(sxy['ID_atu_multipolygons']==i)[0]
        if ind.size==0:
            continue
        dMP['SI_ba_mean'][i]=np.round(np.mean(ba['SI'][ind]),decimals=1)
        dMP['SI_ba_min'][i]=np.min(ba['SI'][ind])
        dMP['SI_ba_max'][i]=np.max(ba['SI'][ind])

    # Add forest cover inventory update year
    for i in range(5):
        dMP['FC' + str(i+1) + ' Year']=-999*np.ones(dMP['ID_Multipolygon'].size)
        dMP['FC' + str(i+1) + ' Status']=np.array(['' for _ in range(dMP['OPENING_ID'].size)],dtype=object)
        dMP['FC' + str(i+1) + ' Type']=np.array(['' for _ in range(dMP['OPENING_ID'].size)],dtype=object)
        dMP['FC' + str(i+1) + ' CC']=-999*np.ones(dMP['ID_Multipolygon'].size)
        dMP['FC' + str(i+1) + ' SPH Tot']=-999*np.ones(dMP['ID_Multipolygon'].size)
        dMP['FC' + str(i+1) + ' SPH WS']=-999*np.ones(dMP['ID_Multipolygon'].size)
        dMP['FC' + str(i+1) + ' Sp1 I']=np.array(['' for _ in range(dMP['OPENING_ID'].size)],dtype=object)
        dMP['FC' + str(i+1) + ' Sp2 I']=np.array(['' for _ in range(dMP['OPENING_ID'].size)],dtype=object)

    for i in range(len(atu_multipolygons)):
    
        ind=np.where(sxy['ID_atu_multipolygons']==i)[0]
    
        # Inventory label variables
        ind_fc=np.array([],dtype=int)
        for j in range(ind.size):
            ind0=np.where(fcinv['IdxToSXY']==ind[j])[0]
            ind_fc=np.append(ind_fc,ind0)
            
        uYear,uIndex=np.unique(fcinv['Year_Updated'][ind_fc]+fcinv['Month_Updated'][ind_fc]/13,return_index=True)
        ord=np.flip(np.argsort(uYear))
    
        for j in range(np.minimum(5,uYear.size)):
            ind_fc0=ind_fc[uIndex[ord[j]]]
            dMP['FC' + str(j+1) + ' Year'][i]=uYear[ord[j]]
            
            dMP['FC' + str(j+1) + ' Status'][i]=cbu.lut_n2s(meta['LUT']['FC_I']['STOCKING_STATUS_CODE'],fcinv['STOCKING_STATUS_CODE'][ind_fc0])[0]
            try:
                dMP['FC' + str(j+1) + ' Type'][i]=cbu.lut_n2s(meta['LUT']['FC_I']['STOCKING_TYPE_CODE'],fcinv['STOCKING_TYPE_CODE'][ind_fc0])[0]
            except:
                pass
            
            dMP['FC' + str(j+1) + ' CC'][i]=fcinv['I_CROWN_CLOSURE_PERCENT'][ind_fc0]
            dMP['FC' + str(j+1) + ' SPH Tot'][i]=fcinv['I_TOTAL_STEMS_PER_HA'][ind_fc0]
            #dMP['FC' + str(j+1) + ' SPH WS'][i]=fcinv['I_TOTAL_WELL_SPACED_STEMS_HA'][ind_fc0]
            if fcinv['I_SPECIES_CODE_1'][ind_fc0]>0:
                dMP['FC' + str(j+1) + ' Sp1 I'][i]=cbu.lut_n2s(meta['LUT']['VRI']['SPECIES_CD_1'],fcinv['I_SPECIES_CODE_1'][ind_fc0])[0]
                if fcinv['I_SPECIES_CODE_2'][ind_fc0]>0:
                    dMP['FC' + str(j+1) + ' Sp2 I'][i]=cbu.lut_n2s(meta['LUT']['VRI']['SPECIES_CD_1'],fcinv['I_SPECIES_CODE_2'][ind_fc0])[0]

        # Silviculture label variables    
        ind_fc=np.array([],dtype=int)
        for j in range(ind.size):
            ind0=np.where(fcsilv['IdxToSXY']==ind[j])[0]
            ind_fc=np.append(ind_fc,ind0)
    
        uYear,uIndex=np.unique(fcsilv['REFERENCE_YEAR'][ind_fc],return_index=True)
        ord=np.flip(np.argsort(uYear))
    
        for j in range(np.minimum(5,uYear.size)):
            ind_fc0=ind_fc[uIndex[ord[j]]]
            dMP['FC' + str(j+1) + ' SPH WS'][i]=fcsilv['S_TOTAL_WELL_SPACED_STEMS_HA'][ind_fc0]
            #if fcsilv['S_SPECIES_CODE_1'][ind_fc0]>0:
        #    dMP['FC' + str(j+1) + ' Sp1 S'][i]=cbu.lut_n2s(meta['LUT']['VRI']['SPECIES_CD_1'],fcsilv['S_SPECIES_CODE_1'][ind_fc0])[0]
        #    dMP['FC' + str(j+1) + ' Sp2 S'][i]=cbu.lut_n2s(meta['LUT']['VRI']['SPECIES_CD_1'],fcsilv['S_SPECIES_CODE_2'][ind_fc0])[0]

    # Add planting info
    num_of_spc=5
    for i in range(num_of_spc):
        dMP['PL_S' + str(i+1) + '_CD']=np.array([' ' for _ in range(dMP['OPENING_ID'].size)],dtype=object)
        dMP['PL_S' + str(i+1) + '_PCT']=np.array([' ' for _ in range(dMP['OPENING_ID'].size)],dtype=object)
        dMP['PL_S' + str(i+1) + '_GW']=np.array([' ' for _ in range(dMP['OPENING_ID'].size)],dtype=object)

    for iMP in range(len(atu_multipolygons)):
    
        indSXY=np.where(sxy['ID_atu_multipolygons']==iMP)[0]
        
        ind_pl=np.array([],dtype=int)
        for j in range(indSXY.size):
            ind0=np.where( (pl['IdxToSXY']==indSXY[j]) & (pl['Year']==dMP['Year'][iMP]) )[0]
            ind_pl=np.append(ind_pl,ind0)
    
        if ind_pl.size==0:
            continue

        #tot_pl=np.sum(pl['NUMBER_PLANTED'][ind_pl])
        #Ord=np.flip(np.argsort(pl['NUMBER_PLANTED'][ind_pl]))
    
        # Genetic worth from seedlot
        gw0=np.zeros(ind_pl.size)
        for j in range(ind_pl.size):
            ind=np.where(par['GW']['SEEDLOT_NUMBER']==pl['SEEDLOT_NUMBER'][ind_pl[j]])[0]
            if ind.size!=0:
                gw0[j]=par['GW']['GENETIC_WORTH_RTNG'][ind[0]]
    
        # Unique species
        u=np.unique(pl['SILV_TREE_SPECIES_CODE'][ind_pl])
    
        # Number for each unique code
        nt=np.zeros(u.size)
        gw=np.zeros(u.size)
        for j in range(u.size):
            ind=np.where(pl['SILV_TREE_SPECIES_CODE'][ind_pl]==u[j])[0]
            nt[j]=np.sum(pl['NUMBER_PLANTED'][ind_pl[ind]])
            gw[j]=np.mean(gw0[ind])
        pct=nt/np.sum(nt)*100
        ord=np.flip(np.argsort(nt))
    
        for j in range(np.minimum(u.size,num_of_spc)):
            dMP['PL_S' + str(j+1) + '_CD'][iMP]=cbu.lut_n2s(meta['LUT']['PL']['SILV_TREE_SPECIES_CODE'],u[ord[j]])[0]
            dMP['PL_S' + str(j+1) + '_PCT'][iMP]=np.round(pct[ord[j]])
            dMP['PL_S' + str(j+1) + '_GW'][iMP]=np.round(gw[ord[j]])

    # Add age at time of fertilization
    dMP['Age @ Fert mean']=-999*np.ones(dMP['ID_Multipolygon'].size)
    try:
        for i in range(len(atu_multipolygons)):
            ind=np.where(sxy['ID_atu_multipolygons']==i)[0]
            if ind.size==0:
                continue
            dMP['Age @ Fert mean'][i]=int(np.mean(meta['AgeAtFert'][ind]))
    except:
        pass

    # Add district
    u=np.unique(dMP['OPENING_ID'])
    dMP['District']=np.array(['' for _ in range(dMP['OPENING_ID'].size)],dtype=object)
    for iU in range(u.size):
        ind1=np.where(dMP['OPENING_ID']==u[iU])[0]
        ind2=np.where(op['OPENING_ID']==u[iU])[0]
        if ind2.size>0:
            dMP['District'][ind1]=cbu.lut_n2s(meta['LUT']['OP']['DISTRICT_NAME'],op['DISTRICT_NAME'][ind2[0]])[0][0:-26]

    # Add GHG benefit
    #it=np.where( (tv>=2018) & (tv<=2050) )[0]
    #dMP['GHG Benefit 2050']=0*np.ones(dMP['ID_Multipolygon'].size)
    #for iMP in range(uMP.size):
    #    ghgb=MosByMP[1]['v2']['Mean']['Sec_NGHGB']['Ensemble Mean'][:,iMP]-MosByMP[0]['v2']['Mean']['Sec_NGHGB']['Ensemble Mean'][:,iMP]
    #    dMP['GHG Benefit 2050'][uMP[iMP]]=int(np.sum(ghgb[it]))

    # Delete any unwanted variables
    #del d['GeomFromCutLyr']

    # Not sure why layot it making it through query
    ind=np.where(dMP['SILV_METHOD_CODE']!='LAYOT')[0]
    for k in dMP.keys():
        dMP[k]=dMP[k][ind]

    # Add disturbance history
    nD=int(40)
    for iD in range(nD):
        dMP['Event' + str(iD+1) + '_Year']=np.array(['' for _ in range(dMP['ID_Multipolygon'].size)],dtype=object)
        dMP['Event' + str(iD+1) + '_Type']=np.array(['' for _ in range(dMP['ID_Multipolygon'].size)],dtype=object)
    for iMP in range(uMP.size):
        indMP=np.where(dMP['ID_Multipolygon']==uMP[iMP])[0]
        indS=np.where(sxy['ID_atu_multipolygons']==uMP[iMP])[0]
        Year0=np.array([])
        Type0=np.array([])
        Mort0=np.array([])
        BaseAffected0=np.array([])
        ProjAffected0=np.array([])
        for iS in range(indS.size):
            Year0=np.append(Year0,dmec[indS[iS]]['Year'])
            Mort0=np.append(Mort0,dmec[indS[iS]]['MortalityFactor'])
            Type0=np.append(Type0,dmec[indS[iS]]['ID_Type'])
            BaseAffected0=np.append(BaseAffected0,dmec[indS[iS]]['ScnAffected'][0])
            ProjAffected0=np.append(ProjAffected0,dmec[indS[iS]]['ScnAffected'][1])
    
        u=np.flip(np.unique(np.column_stack([Year0,Type0,Mort0,BaseAffected0,ProjAffected0]),axis=0),axis=0)    
        for iE in range(np.minimum(nD,u.shape[0])):
            dMP['Event' + str(iE+1) + '_Year'][indMP]=str(np.round(u[iE,0],decimals=2))
            dMP['Event' + str(iE+1) + '_Type'][indMP]=str(cbu.lut_n2s(meta['LUT']['Dist'],u[iE,1])[0]) + ' ' + str(int(u[iE,2])) + ' (' + str(int(u[iE,3])) + '-' + str(int(u[iE,4])) + ')'

    # Add climate
    dMP['tmin_ann']=np.zeros(dMP['ID_Multipolygon'].size)
    dMP['tmean_gs']=np.zeros(dMP['ID_Multipolygon'].size)
    dMP['ws_gs']=np.zeros(dMP['ID_Multipolygon'].size)
    try:
        for iMP in range(uMP.size): 
            indMP=np.where(dMP['ID_Multipolygon']==uMP[iMP])[0]
            indSXY=np.where(sxy['ID_atu_multipolygons']==iMP)[0]
            dMP['tmin_ann'][indMP]=np.mean(clm['tmin_ann'][indSXY])
            dMP['tmean_gs'][indMP]=np.mean(clm['tmean_gs'][indSXY])
            dMP['ws_gs'][indMP]=np.mean(clm['ws_gs'][indSXY])
    except:
        pass
        
    # Add dictonary to dataframe
    df=pd.DataFrame(dMP)

    # Remove -999
    df[df==-999]=np.nan
    df[df==-9999]=np.nan

    # Save
    path=meta['Paths']['Project'] + '\\Inputs\\SummaryActivities_ByMP.xlsx'
    df.to_excel(path,index=False)  

    return dMP

#%% Export AT Layer data to spreadhseet

def ExportSummaryActivities_BySXY(meta,par,atu_multipolygons,sxy,atu,op,burnsev,vri,fcinv,fcsilv,fcres,pl,fire,pest,cut,lul,park,ogmal):
    
    # Get land cover class names
    lcc=gu.ReadExcel(r'C:\Users\rhember\Documents\Code_Python\fcgadgets\cbrunner\Parameters\Parameters_VRI.xlsx','LCC')   
    
    #--------------------------------------------------------------------------
    # Reverse crosswalk for FC silv polygon (too slow using lut_num2)
    #--------------------------------------------------------------------------
    
    spn_id_I=np.zeros(len(meta['LUT']['FC_I']['SILV_POLYGON_NUMBER']))
    spn_cd_I=np.array(['' for _ in range(spn_id_I.size)],dtype=object)
    cnt=0
    for k in meta['LUT']['FC_I']['SILV_POLYGON_NUMBER'].keys():
        spn_id_I[cnt]=meta['LUT']['FC_I']['SILV_POLYGON_NUMBER'][k]
        spn_cd_I[cnt]=str(k)
        cnt=cnt+1
    lut_spn_I={}  
    for i in range(spn_id_I.size):
        lut_spn_I[int(spn_id_I[i])]=spn_cd_I[i]
        
    spn_id_S=np.zeros(len(meta['LUT']['FC_S']['SILV_POLYGON_NUMBER']))
    spn_cd_S=np.array(['' for _ in range(spn_id_S.size)],dtype=object)
    cnt=0
    for k in meta['LUT']['FC_I']['SILV_POLYGON_NUMBER'].keys():
        spn_id_S[cnt]=meta['LUT']['FC_I']['SILV_POLYGON_NUMBER'][k]
        spn_cd_S[cnt]=str(k)
        cnt=cnt+1
    lut_spn_S={}  
    for i in range(spn_id_S.size):
        lut_spn_S[int(spn_id_S[i])]=spn_cd_S[i]
     
    #--------------------------------------------------------------------------
    # Reverse crosswalk for FC speces (too slow using lut_num2)
    #--------------------------------------------------------------------------    
    #   if fcinv['I_SPECIES_CODE_1'][ind_fc0]>0:
    #    dMP['FC' + str(j+1) + ' Sp1 I'][i]=cbu.lut_n2s(meta['LUT']['VRI']['SPECIES_CD_1'],fcinv['I_SPECIES_CODE_1'][ind_fc0])[0]
       
    spc_id=np.zeros(len(meta['LUT']['VRI']['SPECIES_CD_1']))
    spc_cd=np.array(['' for _ in range(spc_id.size)],dtype=object)
    cnt=0
    for k in meta['LUT']['VRI']['SPECIES_CD_1'].keys():
        spc_id[cnt]=meta['LUT']['VRI']['SPECIES_CD_1'][k]
        spc_cd[cnt]=str(k)
        cnt=cnt+1
    lut_spc={}  
    for i in range(spc_id.size):
        lut_spc[int(spc_id[i])]=spc_cd[i]
    
    #--------------------------------------------------------------------------
    # Reverse crosswalk for FC damage agents (too slow using lut_num2)
    #--------------------------------------------------------------------------
    
    da_id=np.zeros(len(meta['LUT']['Pest']['PEST_SPECIES_CODE']))
    da_cd=np.array(['' for _ in range(da_id.size)],dtype=object)
    cnt=0
    for k in meta['LUT']['Pest']['PEST_SPECIES_CODE'].keys():
        da_id[cnt]=meta['LUT']['Pest']['PEST_SPECIES_CODE'][k]
        da_cd[cnt]=str(k)
        cnt=cnt+1
    lut_da={}  
    for i in range(da_id.size):
        lut_da[int(da_id[i])]=da_cd[i]    
        
    #--------------------------------------------------------------------------
    # Initialize dictionary
    #--------------------------------------------------------------------------
    
    d={}
    n_init=int(2e6)
    cnt=0
    
    d['ID_SXY']=-999*np.ones(n_init)
    d['ID_Multipolygon']=-999*np.ones(n_init)
    d['OPENING_ID']=-999*np.ones(n_init)
    d['SILV_POLYGON_NUMBER']=np.array(['' for _ in range(n_init)],dtype=object)
    d['Geom Source']=np.array(['' for _ in range(n_init)],dtype=object)
    d['Area']=-999*np.ones(n_init)
    d['AEF_ATU']=np.zeros(n_init)
    d['LCC']=np.array(['' for _ in range(n_init)],dtype=object)
    d['BGCz']=np.array(['' for _ in range(n_init)],dtype=object)
    d['BGCsz']=np.array(['' for _ in range(n_init)],dtype=object)
    d['BGCv']=np.array(['' for _ in range(n_init)],dtype=object)
    d['Year']=-999*np.ones(n_init)
    d['Dist Type']=np.array(['' for _ in range(n_init)],dtype=object)
    #d['Dist Mort']=np.array(['' for _ in range(n_init)],dtype=object)
    d['SBC']=np.array(['' for _ in range(n_init)],dtype=object)
    d['STC']=np.array(['' for _ in range(n_init)],dtype=object)
    d['SMC']=np.array(['' for _ in range(n_init)],dtype=object)
    d['SOC']=np.array(['' for _ in range(n_init)],dtype=object)
    d['FSC']=np.array(['' for _ in range(n_init)],dtype=object)
    #d['Area AT']=-999*np.ones(n_init)
    
    d['FC ID']=-999*np.ones(n_init)
    d['FC Status']=np.array(['' for _ in range(n_init)],dtype=object)    
    d['FC Type']=np.array(['' for _ in range(n_init)],dtype=object)    
    d['FC CC']=-999*np.ones(n_init)
    d['FC SPH TOT']=-999*np.ones(n_init)
    d['FC SPH TWS']=-999*np.ones(n_init)
    d['FC SPH WS']=-999*np.ones(n_init)
    d['FC SPH FG']=-999*np.ones(n_init)
    
    d['FC I Spc1 CD']=np.array([' ' for _ in range(n_init)],dtype=object)
    d['FC I Spc2 CD']=np.array([' ' for _ in range(n_init)],dtype=object)
    d['FC S Spc1 CD']=np.array([' ' for _ in range(n_init)],dtype=object)
    d['FC S Spc2 CD']=np.array([' ' for _ in range(n_init)],dtype=object)
    
    d['FC I DA1 CD']=np.array([' ' for _ in range(n_init)],dtype=object)
    d['FC I DA1 PCT']=-999*np.ones(n_init)
    d['FC I DA2 CD']=np.array([' ' for _ in range(n_init)],dtype=object)
    d['FC I DA2 PCT']=-999*np.ones(n_init)
    d['FC I DA3 CD']=np.array([' ' for _ in range(n_init)],dtype=object)
    d['FC I DA3 PCT']=-999*np.ones(n_init)
    d['FC S DA1 CD']=np.array([' ' for _ in range(n_init)],dtype=object)
    d['FC S DA1 PCT']=-999*np.ones(n_init)
    d['FC S DA2 CD']=np.array([' ' for _ in range(n_init)],dtype=object)
    d['FC S DA2 PCT']=-999*np.ones(n_init)
    d['FC S DA3 CD']=np.array([' ' for _ in range(n_init)],dtype=object)
    d['FC S DA3 PCT']=-999*np.ones(n_init)
    
    d['PL SPH']=-999*np.ones(n_init)    
    num_of_spc=5
    for i in range(num_of_spc):
        d['PL_S' + str(i+1) + '_CD']=np.array([' ' for _ in range(n_init)],dtype=object)
        d['PL_S' + str(i+1) + '_PCT']=np.array([' ' for _ in range(n_init)],dtype=object)
        d['PL_S' + str(i+1) + '_GW']=np.array([' ' for _ in range(n_init)],dtype=object)
    
    d['VRI RefYear']=-999*np.ones(n_init)
    d['VRI FREE_TO_GROW_IND']=np.array([' ' for _ in range(n_init)],dtype=object)
    d['VRI EST_COVERAGE_PCT_1']=-999*np.ones(n_init)
    d['VRI_LIVE_STEMS_PER_HA']=-999*np.ones(n_init)
    d['VRI_PROJ_AGE_1']=-999*np.ones(n_init)
    d['VRI LIVE_STAND_VOLUME_125']=-999*np.ones(n_init)
    
    #--------------------------------------------------------------------------
    # Loop through and populate
    #--------------------------------------------------------------------------
    
    for iSXY in range(sxy['x'].size):
        
        #----------------------------------------------------------------------
        # ATU
        #----------------------------------------------------------------------
        indAT=np.where(atu['IdxToSXY']==iSXY)[0]
        for iAT in range(indAT.size):
            ind=indAT[iAT]
            d['ID_SXY'][cnt]=iSXY
            d['ID_Multipolygon'][cnt]=sxy['ID_atu_multipolygons'][iSXY]
            d['OPENING_ID'][cnt]=sxy['OPENING_ID'][iSXY]
            d['Year'][cnt]=atu['Year'][ind]+atu['Month'][ind]/13
            d['Dist Type'][cnt]=cbu.lut_n2s(meta['LUT']['ATU']['DISTURBANCE_CODE'],atu['DISTURBANCE_CODE'][ind])[0]
            d['SBC'][cnt]=cbu.lut_n2s(meta['LUT']['ATU']['SILV_BASE_CODE'],atu['SILV_BASE_CODE'][ind])[0]
            d['STC'][cnt]=cbu.lut_n2s(meta['LUT']['ATU']['SILV_TECHNIQUE_CODE'],atu['SILV_TECHNIQUE_CODE'][ind])[0]
            d['SMC'][cnt]=cbu.lut_n2s(meta['LUT']['ATU']['SILV_METHOD_CODE'],atu['SILV_METHOD_CODE'][ind])[0]            
            d['SOC'][cnt]=cbu.lut_n2s(meta['LUT']['ATU']['SILV_OBJECTIVE_CODE_1'],atu['SILV_OBJECTIVE_CODE_1'][ind])[0]   
            d['FSC'][cnt]=cbu.lut_n2s(meta['LUT']['ATU']['SILV_FUND_SOURCE_CODE'],atu['SILV_FUND_SOURCE_CODE'][ind])[0]
            d['Area'][cnt]=np.round(atu['ACTUAL_TREATMENT_AREA'][ind],decimals=1)
            if (atu['ACTUAL_PLANTED_NUMBER'][ind])>0 & (atu['ACTUAL_TREATMENT_AREA'][ind]>0):
                d['PL SPH'][cnt]=atu['ACTUAL_PLANTED_NUMBER'][ind]/atu['ACTUAL_TREATMENT_AREA'][ind]
            cnt=cnt+1
        
        #----------------------------------------------------------------------
        # Forest cover inventory layer
        #----------------------------------------------------------------------
        
        indFC=np.where(fcinv['IdxToSXY']==iSXY)[0]
        for iFC in range(indFC.size):
            ind=indFC[iFC]
            #if np.isin(fcinv['FOREST_COVER_ID'][ind],d['FOREST_COVER_ID'])==True:
            #    continue
            d['ID_SXY'][cnt]=iSXY
            d['ID_Multipolygon'][cnt]=sxy['ID_atu_multipolygons'][iSXY]
            d['OPENING_ID'][cnt]=sxy['OPENING_ID'][iSXY]
            d['Area'][cnt]=fcinv['SILV_POLYGON_AREA'][ind]
            if fcinv['SILV_POLYGON_NUMBER'][ind]!=-9999:
                d['SILV_POLYGON_NUMBER'][cnt]=lut_spn_I[fcinv['SILV_POLYGON_NUMBER'][ind]]
            d['Year'][cnt]=fcinv['Year_Updated'][ind]+fcinv['Month_Updated'][ind]/13
            d['FC ID'][cnt]=fcinv['FOREST_COVER_ID'][ind]
            
            try:
                d['FC Status'][cnt]=cbu.lut_n2s(meta['LUT']['FC_I']['STOCKING_STATUS_CODE'],fcinv['STOCKING_STATUS_CODE'][ind])[0]      
            except:
                pass
            
            try:
                d['FC Type'][cnt]=cbu.lut_n2s(meta['LUT']['FC_I']['STOCKING_TYPE_CODE'],fcinv['STOCKING_TYPE_CODE'][ind])[0]
            except:
                pass
            
            d['FC CC'][cnt]=fcinv['I_CROWN_CLOSURE_PERCENT'][ind]
            d['FC SPH TOT'][cnt]=fcinv['I_TOTAL_STEMS_PER_HA'][ind]
            d['FC SPH TWS'][cnt]=fcinv['I_TOTAL_WELL_SPACED_STEMS_HA'][ind]
            #d['FC SPH WS'][cnt]=fcinv['I_WELL_SPACED_STEMS_HA'][ind]
            d['FC SPH FG'][cnt]=fcinv['I_FREE_GROWING_STEMS_PER_HA'][ind]
            
            if fcinv['I_SPECIES_CODE_1'][ind]>-999:
                d['FC I Spc1 CD'][cnt]=lut_spc[fcinv['I_SPECIES_CODE_1'][ind]]
            if fcinv['I_SPECIES_CODE_2'][ind]>-999:
                d['FC I Spc2 CD'][cnt]=lut_spc[fcinv['I_SPECIES_CODE_2'][ind]]   
            
            if fcinv['I_DA1 CD'][ind]!=-9999:
                d['FC I DA1 CD'][cnt]=lut_da[fcinv['I_DA1 CD'][ind]]
                d['FC I DA1 PCT'][cnt]=fcinv['I_DA1 PCT'][ind]
            if fcinv['I_DA2 CD'][ind]!=-9999:
                d['FC I DA2 CD'][cnt]=lut_da[fcinv['I_DA2 CD'][ind]]
                d['FC I DA2 PCT'][cnt]=fcinv['I_DA2 PCT'][ind]
            if fcinv['I_DA3 CD'][ind]!=-9999:
                d['FC I DA3 CD'][cnt]=lut_da[fcinv['I_DA3 CD'][ind]]
                d['FC I DA3 PCT'][cnt]=fcinv['I_DA3 PCT'][ind]
            if fcinv['S_DA1 CD'][ind]!=-9999:
                d['FC S DA1 CD'][cnt]=lut_da[fcinv['S_DA1 CD'][ind]]
                d['FC S DA1 PCT'][cnt]=fcinv['S_DA1 PCT'][ind]
            if fcinv['S_DA2 CD'][ind]!=-9999:
                d['FC S DA2 CD'][cnt]=lut_da[fcinv['S_DA2 CD'][ind]]
                d['FC S DA2 PCT'][cnt]=fcinv['S_DA2 PCT'][ind]
            if fcinv['S_DA3 CD'][ind]!=-9999:
                d['FC S DA3 CD'][cnt]=lut_da[fcinv['S_DA3 CD'][ind]]
                d['FC S DA3 PCT'][cnt]=fcinv['S_DA3 PCT'][ind] 

            cnt=cnt+1
        #----------------------------------------------------------------------
        # Add harvest evnets from consolidated cutblocks DB
        #----------------------------------------------------------------------
        
        indC=np.where(cut['IdxToSXY']==iSXY)[0]
        for iC in range(indC.size):
            ind=indC[iC]
            d['ID_SXY'][cnt]=iSXY
            d['ID_Multipolygon'][cnt]=sxy['ID_atu_multipolygons'][iSXY]
            d['OPENING_ID'][cnt]=sxy['OPENING_ID'][iSXY]
            d['Year'][cnt]=cut['HARVEST_YEAR'][ind]
            try:
                d['Area'][cnt]=cut['AREA_HA'][ind]
            except:
                pass        
            d['Dist Type'][cnt]='Harvest'
            cnt=cnt+1
        
        #----------------------------------------------------------------------
        # Add wildfire
        #----------------------------------------------------------------------
        
        indC=np.where(fire['IdxToSXY']==iSXY)[0]
        for iC in range(indC.size):
            ind=indC[iC]
            d['ID_SXY'][cnt]=iSXY
            d['ID_Multipolygon'][cnt]=sxy['ID_atu_multipolygons'][iSXY]
            d['OPENING_ID'][cnt]=sxy['OPENING_ID'][iSXY]
            d['Year'][cnt]=fire['FIRE_YEAR'][ind]   
            d['Dist Type'][cnt]='Wildfire'
            cnt=cnt+1
    
    #--------------------------------------------------------------------------
    # Truncate dictionary  
    #--------------------------------------------------------------------------
    
    for k in d.keys():
        d[k]=d[k][0:cnt]
    
    #--------------------------------------------------------------------------
    # Add forest cover silv layers
    #--------------------------------------------------------------------------
    
    for iD in range(d['ID_SXY'].size):
        
        if d['FC ID'][iD]<0:
            continue
        
        if d['FC Status'][iD]=='':
            continue
        
        indFC=np.where(fcsilv['FOREST_COVER_ID']==d['FC ID'][iD])[0]
        if indFC.size==0:
            continue        
        
        for iFC in range(indFC.size):
            ind=indFC[iFC]
            if fcsilv['S_TOTAL_WELL_SPACED_STEMS_HA'][ind]>0:
                d['FC SPH TWS'][iD]=fcsilv['S_TOTAL_WELL_SPACED_STEMS_HA'][ind]
            if fcsilv['S_FREE_GROWING_STEMS_PER_HA'][ind]>0:
                d['FC SPH FG'][iD]=fcsilv['S_FREE_GROWING_STEMS_PER_HA'][ind]
            
            if fcsilv['S_SPECIES_CODE_1'][ind]>-999:
                d['FC S Spc1 CD'][iD]=lut_spc[fcsilv['S_SPECIES_CODE_1'][ind]]
            if fcsilv['S_SPECIES_CODE_2'][ind]>-999:
                d['FC S Spc2 CD'][iD]=lut_spc[fcsilv['S_SPECIES_CODE_2'][ind]]
    
    #--------------------------------------------------------------------------
    # Add info from atu_multipolygon 
    #--------------------------------------------------------------------------
    
    uMP=np.unique(d['ID_Multipolygon']).astype(int)
    for iMP in range(uMP.size):
        
        indD=np.where(d['ID_Multipolygon']==uMP[iMP])[0]
        
        if atu_multipolygons[uMP[iMP]]['GeomFromOpLyr']==1:
            d['Geom Source'][indD]='OP'
        elif atu_multipolygons[uMP[iMP]]['GeomFromFcLyr']==1:
            d['Geom Source'][indD]='FC ART'
        else:
            d['Geom Source'][indD]='AT'
        
        A=atu_multipolygons[uMP[iMP]]['ACTUAL_TREATMENT_AREA']
        if A!=None:
            indS=np.where(sxy['ID_atu_multipolygons']==uMP[iMP])[0]
            d['AEF_ATU'][indD]=np.round(A/indS.size,3)
    
    #--------------------------------------------------------------------------
    # Add VRI
    #--------------------------------------------------------------------------
    
    for iD in range(d['ID_SXY'].size):
        indVRI=np.where(vri['IdxToSXY']==d['ID_SXY'][iD])[0]
        for iVRI in range(indVRI.size):
            ind=indVRI[iVRI]
            cd=cbu.lut_n2s(meta['LUT']['VRI']['LAND_COVER_CLASS_CD_1'],vri['LAND_COVER_CLASS_CD_1'][ind])[0]
            ind1=np.where(lcc['Code']==cd)[0]
            if ind1.size>0:
                d['LCC'][iD]=lcc['Name'][ind1[0]]
            d['BGCz'][iD]=cbu.lut_n2s(meta['LUT']['VRI']['BEC_ZONE_CODE'],vri['BEC_ZONE_CODE'][ind])[0]
            d['BGCsz'][iD]=cbu.lut_n2s(meta['LUT']['VRI']['BEC_SUBZONE'],vri['BEC_SUBZONE'][ind])[0]
            d['BGCv'][iD]=str(cbu.lut_n2s(meta['LUT']['VRI']['BEC_VARIANT'],vri['BEC_VARIANT'][ind])[0])
            
            d['VRI RefYear'][iD]=vri['REFERENCE_YEAR'][ind]
            d['VRI FREE_TO_GROW_IND'][iD]=cbu.lut_n2s(meta['LUT']['VRI']['FREE_TO_GROW_IND'],vri['FREE_TO_GROW_IND'][ind])[0]
            d['VRI EST_COVERAGE_PCT_1'][iD]=vri['EST_COVERAGE_PCT_1'][ind]
            d['VRI_LIVE_STEMS_PER_HA'][iD]=vri['VRI_LIVE_STEMS_PER_HA'][ind]
            d['VRI_PROJ_AGE_1'][iD]=vri['PROJ_AGE_1'][ind]
            d['VRI LIVE_STAND_VOLUME_125'][iD]=vri['LIVE_STAND_VOLUME_125'][ind]
    
    #--------------------------------------------------------------------------
    # Add planting info
    #--------------------------------------------------------------------------
    
    for iD in range(d['ID_SXY'].size):
        
        if d['SBC'][iD]!='PL':
            continue
        
        ind_pl=np.where( (pl['IdxToSXY']==d['ID_SXY'][iD]) & (pl['Year']==np.floor(d['Year'][iD])) )[0]
    
        if ind_pl.size==0:
            continue
    
        #tot_pl=np.sum(pl['NUMBER_PLANTED'][ind_pl])
        #Ord=np.flip(np.argsort(pl['NUMBER_PLANTED'][ind_pl]))
    
        # Genetic worth from seedlot
        gw0=np.zeros(ind_pl.size)
        for j in range(ind_pl.size):
            ind=np.where(par['GW']['SEEDLOT_NUMBER']==pl['SEEDLOT_NUMBER'][ind_pl[j]])[0]
            if ind.size!=0:
                gw0[j]=par['GW']['GENETIC_WORTH_RTNG'][ind[0]]
    
        # Unique species
        u=np.unique(pl['SILV_TREE_SPECIES_CODE'][ind_pl])
    
        # Number for each unique code
        nt=np.zeros(u.size)
        gw=np.zeros(u.size)
        for j in range(u.size):
            ind=np.where(pl['SILV_TREE_SPECIES_CODE'][ind_pl]==u[j])[0]
            nt[j]=np.sum(pl['NUMBER_PLANTED'][ind_pl[ind]])
            gw[j]=np.mean(gw0[ind])
        pct=nt/np.sum(nt)*100
        ord=np.flip(np.argsort(nt))
    
        for j in range(np.minimum(u.size,num_of_spc)):
            d['PL_S' + str(j+1) + '_CD'][iD]=cbu.lut_n2s(meta['LUT']['PL']['SILV_TREE_SPECIES_CODE'],u[ord[j]])[0]
            d['PL_S' + str(j+1) + '_PCT'][iD]=np.round(pct[ord[j]])
            d['PL_S' + str(j+1) + '_GW'][iD]=np.round(gw[ord[j]])
    
    #--------------------------------------------------------------------------
    # Rounding
    #--------------------------------------------------------------------------
    
    d['Year']=np.round(d['Year'],decimals=2)
    d['PL SPH']=d['PL SPH'].astype(int)
    
    #--------------------------------------------------------------------------
    # Add dictonary to dataframe
    #--------------------------------------------------------------------------
    
    df=pd.DataFrame(d)

    # Sort
    df=df.sort_values(by=['ID_SXY','Year'])

    #--------------------------------------------------------------------------
    # Remove -999 (setting to nan will appear as blank)
    #--------------------------------------------------------------------------
    
    df[df==-999]=np.nan
    df[df==-9999]=np.nan
    df[df=='Unidentified']=np.nan

    #--------------------------------------------------------------------------
    # Save
    #--------------------------------------------------------------------------
    
    path=meta['Paths']['Project'] + '\\Inputs\\SummaryActivities_BySXY.xlsx'
    df.to_excel(path,index=False)
        
    return d


#%% Export summary attributes by PP number

def ExportSummaryActivities_ByPP(meta,dAdmin,dMP,dFCI_Summary):

    MaterialMilestoneTypes=['Aerial Spray','Browse Protectors','Direct Seeding',
            'Disc Trenching','Dwarf Mistletoe Control','Fertilization Aerial',
            'Fertilization Hand','Fertilization Teabag','Incremental Haul','Knockdown','Planting','Ripping','Thinning']

    dPP={}
    dPP['FIA_PROJECT_ID']=np.unique(dMP['FIA_PROJECT_ID'])
    dPP['Num Openings']=np.nan*np.ones(dPP['FIA_PROJECT_ID'].size)
    dPP['Year Start']=np.nan*np.ones(dPP['FIA_PROJECT_ID'].size)
    dPP['Year Last']=np.nan*np.ones(dPP['FIA_PROJECT_ID'].size)
    dPP['Project Type']=np.array(['None' for _ in range(dPP['FIA_PROJECT_ID'].size)],dtype=object)
    dPP['Milestone Type']=np.array(['None' for _ in range(dPP['FIA_PROJECT_ID'].size)],dtype=object)
    dPP['BGC Zone 1']=np.array(['None' for _ in range(dPP['FIA_PROJECT_ID'].size)],dtype=object)
    dPP['BGC Zone 2']=np.array(['None' for _ in range(dPP['FIA_PROJECT_ID'].size)],dtype=object)
    dPP['ACTUAL_TREATMENT_AREA']=np.nan*np.ones(dPP['FIA_PROJECT_ID'].size)
    dPP['ACTUAL_PLANTED_NUMBER']=np.nan*np.ones(dPP['FIA_PROJECT_ID'].size)
    dPP['Age @ Fert']=np.nan*np.ones(dPP['FIA_PROJECT_ID'].size)
    dPP['SI BA AreaWeighted']=np.nan*np.ones(dPP['FIA_PROJECT_ID'].size)
    dPP['Year FC Update']=np.nan*np.ones(dPP['FIA_PROJECT_ID'].size)
    dPP['SPH FC Inv']=np.nan*np.ones(dPP['FIA_PROJECT_ID'].size)
    dPP['SPH Planted']=np.nan*np.ones(dPP['FIA_PROJECT_ID'].size)
    dPP['Spc1 FC Inv']=np.array(['None' for _ in range(dPP['FIA_PROJECT_ID'].size)],dtype=object)
    dPP['Spc2 FC Inv']=np.array(['None' for _ in range(dPP['FIA_PROJECT_ID'].size)],dtype=object)
    dPP['Spc3 FC Inv']=np.array(['None' for _ in range(dPP['FIA_PROJECT_ID'].size)],dtype=object)
    dPP['Spc1 Planted']=np.array(['None' for _ in range(dPP['FIA_PROJECT_ID'].size)],dtype=object)
    dPP['Spc2 Planted']=np.array(['None' for _ in range(dPP['FIA_PROJECT_ID'].size)],dtype=object)
    dPP['Spc3 Planted']=np.array(['None' for _ in range(dPP['FIA_PROJECT_ID'].size)],dtype=object)
    dPP['GHGB50 (tCO2e/ha)']=np.nan*np.ones(dPP['FIA_PROJECT_ID'].size)
    dPP['GeomFromOpLyr']=np.nan*np.ones(dPP['FIA_PROJECT_ID'].size)
    dPP['GeomFromFcLyr']=np.nan*np.ones(dPP['FIA_PROJECT_ID'].size)
    dPP['District']=np.array(['None' for _ in range(dPP['FIA_PROJECT_ID'].size)],dtype=object)
    dPP['Fire Num']=np.array(['None' for _ in range(dPP['FIA_PROJECT_ID'].size)],dtype=object)
    
    for i in range(dPP['FIA_PROJECT_ID'].size):
        
        # Indices
        iMP=np.where(dMP['FIA_PROJECT_ID']==dPP['FIA_PROJECT_ID'][i])[0]    
        iAdmin=np.where( (dAdmin['Estimation Method']=='Completed') & (dAdmin['PP Number']==dPP['FIA_PROJECT_ID'][i]) & (dAdmin['Removed']!='Removed') & (np.isin(dAdmin['Milestone Type'],MaterialMilestoneTypes)==True) )[0]
        iD2=np.where( (dFCI_Summary['PP Number']==dPP['FIA_PROJECT_ID'][i]) & (dFCI_Summary['Estimation Method']=='Completed') )[0]
        
        dPP['Num Openings'][i]=np.unique(dMP['OPENING_ID'][iMP]).size
    
        dPP['Year Start'][i]=np.min(dMP['Year'][iMP])
        dPP['Year Last'][i]=np.max(dMP['Year'][iMP])

        uBGCZ,cBGCZ=np.unique(dMP['BGC Zone 1'][iMP],return_counts=True)
        ord=np.flip(np.argsort(cBGCZ))
        dPP['BGC Zone 1'][i]=uBGCZ[ord[0]]
        if ord.size>1:
            dPP['BGC Zone 2'][i]=uBGCZ[ord[1]]        
    
        dPP['ACTUAL_TREATMENT_AREA'][i]=np.nansum(dMP['ACTUAL_TREATMENT_AREA'][iMP])
        
        dPP['ACTUAL_PLANTED_NUMBER'][i]=np.nansum(dMP['ACTUAL_PLANTED_NUMBER'][iMP])
        
        dPP['SI BA AreaWeighted'][i]=np.sum(dMP['SI_ba_mean'][iMP]*dMP['ACTUAL_TREATMENT_AREA'][iMP])/np.sum(dMP['ACTUAL_TREATMENT_AREA'][iMP])
    
        dPP['Age @ Fert'][i]=np.sum(dMP['Age @ Fert mean'][iMP]*dMP['ACTUAL_TREATMENT_AREA'][iMP])/np.sum(dMP['ACTUAL_TREATMENT_AREA'][iMP])
    
        if iAdmin.size>0:
            dPP['Project Type'][i]=dAdmin['Project Type'][iAdmin[0]]
            dPP['Milestone Type'][i]=dAdmin['Milestone Type'][iAdmin[0]]
         
        dPP['District'][i]=dMP['District'][iMP[0]]
    
        dPP['GeomFromOpLyr'][i]=np.nansum(dMP['GeomFromOpLyr'][iMP])
        dPP['GeomFromFcLyr'][i]=np.nansum(dMP['GeomFromFcLyr'][iMP])
    
        # Add GHG benefit from FCI DB   
        if dFCI_Summary['Area (ha)'][iD2]>0:
            dPP['GHGB50 (tCO2e/ha)'][i]=dFCI_Summary['GHG Benefit Cumu 2050 (tCO2e)'][iD2]/dFCI_Summary['Area (ha)'][iD2]

    # Convert to dataframe
    df=pd.DataFrame.from_dict(dPP)
    df.to_excel(meta['Paths']['Project'] + '\\Inputs\\SummaryActivities_ByPP.xlsx',index=False)  
    
    return dPP

#%% QA plot time series by Multipolygon
    
def QA_Plot_ByMultiPolygon(meta,uMP,ivlMP,iScnForArea,ivlT,tv,it,MosByMP,iB,iP):
    
    for iMP in range(0,uMP.size,ivlMP):
    
        # Get area affected for multipolygon    
        A=cbu.AreaAffectedInSingleMultipolygon(meta,iScnForArea,ivlT,tv,MosByMP,iMP)
        
        #atu_multipolygons[uMP[iMP]]
    
        lw1=1; cle1=[0,0,1]; cle2=[1,0,0]; ms=3; aw=0.28; ah=0.22;
        xlim=[tv[it[0]],tv[it[-1]]]; xticks=np.arange(1500,2200,20);
    
        plt.close('all'); fig,ax=plt.subplots(4,3,figsize=gu.cm2inch(28,20))
   
        pl_d=[None]*len(A['Nat Dist']); nams_d=[None]*len(A['Nat Dist']);
        for i in range(len(A['Nat Dist'])):
            bottom=0; 
            if i!=0:
                for j in range(i):
                    bottom=bottom+A['Nat Dist'][j]['Data']
            pl_d[i]=ax[0,0].bar(A['tv'],A['Nat Dist'][i]['Data'],ivlT,color=A['Nat Dist'][i]['Color'],bottom=bottom)
            nams_d[i]=A['Nat Dist'][i]['Name']
        ax[0,0].legend(pl_d,nams_d,loc='upper left',bbox_to_anchor=(0.05,0.99),labelspacing=0.12,facecolor=[1,1,1],frameon=False)
        ax[0,0].set(position=[0.04,0.77,aw,ah],xlim=xlim,xticks=xticks,ylabel='Area affected (ha)');

        pl_m=[None]*len(A['Management']); nams_m=[None]*len(A['Management']);
        for i in range(len(A['Management'])):
            bottom=0; 
            if i!=0:
                for j in range(i):
                    bottom=bottom+A['Management'][j]['Data']
            pl_m[i]=ax[0,1].bar(A['tv'],A['Management'][i]['Data'],ivlT,color=A['Management'][i]['Color'],bottom=bottom)
            nams_m[i]=A['Management'][i]['Name']
        ax[0,1].legend(pl_m,nams_m,loc='upper left',bbox_to_anchor=(0.05,0.99),labelspacing=0.12,facecolor=[1,1,1],frameon=False)
        ax[0,1].set(position=[0.37,0.77,aw,ah],xlim=xlim,xticks=xticks,xlabel='',ylabel='Area affected (ha)');

        ax[0,2].plot(tv,MosByMP[iB]['v1']['Mean']['A']['Ensemble Mean'][:,iMP],'-',color=cle1,linewidth=lw1)
        ax[0,2].plot(tv,MosByMP[iP]['v1']['Mean']['A']['Ensemble Mean'][:,iMP],'--',color=cle2,linewidth=lw1)
        ax[0,2].set(position=[0.71,0.77,aw,ah],xlim=xlim,xticks=xticks,xlabel='',ylabel='Age, years')
    
        ax[1,0].plot(tv,MosByMP[iB]['v1']['Mean']['C_Biomass_Tot']['Ensemble Mean'][:,iMP],'-',color=cle1,linewidth=lw1)
        ax[1,0].plot(tv,MosByMP[iP]['v1']['Mean']['C_Biomass_Tot']['Ensemble Mean'][:,iMP],'--',color=cle2,linewidth=lw1)
        ax[1,0].set(position=[0.04,0.53,aw,ah],xlim=xlim,xlabel='',ylabel='Biomass (MgC/ha)')
    
        ax[1,1].plot(tv,MosByMP[iB]['v1']['Mean']['C_DeadWood_Tot']['Ensemble Mean'][:,iMP],'-',color=cle1,linewidth=lw1)
        ax[1,1].plot(tv,MosByMP[iP]['v1']['Mean']['C_DeadWood_Tot']['Ensemble Mean'][:,iMP],'--',color=cle2,linewidth=lw1)
        ax[1,1].set(position=[0.37,0.53,aw,ah],xlim=xlim,xticks=xticks,xlabel='',ylabel='Dead wood (MgC/ha)')
        
        ax[1,2].plot(tv,MosByMP[iB]['v1']['Mean']['C_G_Gross_Tot']['Ensemble Mean'][:,iMP],'-',color=cle1,linewidth=lw1)
        ax[1,2].plot(tv,MosByMP[iP]['v1']['Mean']['C_G_Gross_Tot']['Ensemble Mean'][:,iMP],'--',color=cle2,linewidth=lw1)
        ax[1,2].set(position=[0.71,0.53,aw,ah],xlim=xlim,xticks=xticks,xlabel='',ylabel='NPP (MgC/ha/yr)')
    
        ax[2,0].plot(tv,MosByMP[iB]['v1']['Mean']['C_RH_Tot']['Ensemble Mean'][:,iMP],'-',color=cle1,linewidth=lw1,markersize=ms+1)
        ax[2,0].plot(tv,MosByMP[iP]['v1']['Mean']['C_RH_Tot']['Ensemble Mean'][:,iMP],'--',color=cle2,linewidth=lw1,markersize=ms)
        ax[2,0].set(position=[0.04,0.285,aw,ah],xlim=xlim,xticks=xticks,xlabel='',ylabel='RH (MgC/ha/yr)')
    
        ax[2,1].plot(tv,MosByMP[iB]['v1']['Mean']['E_CO2e_LULUCF_Wildfire']['Ensemble Mean'][:,iMP],'o-',color=cle1,linewidth=lw1,markersize=ms+1)
        ax[2,1].plot(tv,MosByMP[iP]['v1']['Mean']['E_CO2e_LULUCF_Wildfire']['Ensemble Mean'][:,iMP],'s--',color=cle2,linewidth=lw1,markersize=ms)
        ax[2,1].set(position=[0.37,0.285,aw,ah],xlim=xlim,xticks=xticks,xlabel='',ylabel='Wildfire emissions (MgC/ha/yr)')
        
        ax[2,2].plot(tv,MosByMP[iB]['v1']['Mean']['E_CO2e_LULUCF_OpenBurning']['Ensemble Mean'][:,iMP],'o-',color=cle1,linewidth=lw1,markersize=ms+1)
        ax[2,2].plot(tv,MosByMP[iP]['v1']['Mean']['E_CO2e_LULUCF_OpenBurning']['Ensemble Mean'][:,iMP],'s--',color=cle2,linewidth=lw1,markersize=ms)
        ax[2,2].set(position=[0.71,0.285,aw,ah],xlim=xlim,xticks=xticks,xlabel='',ylabel='Open burning emissions (MgC/ha/yr)')
    
        ax[3,0].plot(tv,MosByMP[iB]['v1']['Mean']['C_ToMill']['Ensemble Mean'][:,iMP],'o-',color=cle1,linewidth=lw1,markersize=ms+1)
        ax[3,0].plot(tv,MosByMP[iP]['v1']['Mean']['C_ToMill']['Ensemble Mean'][:,iMP],'s--',color=cle2,linewidth=lw1,markersize=ms)
        ax[3,0].set(position=[0.04,0.04,aw,ah],xlim=xlim,xticks=xticks,xlabel='',ylabel='Removals (MgC/ha)')
        
        ax[3,1].plot(tv,MosByMP[iB]['v1']['Mean']['E_CO2e_LULUCF_HWP']['Ensemble Mean'][:,iMP],'-',color=cle1,linewidth=lw1)
        ax[3,1].plot(tv,MosByMP[iP]['v1']['Mean']['E_CO2e_LULUCF_HWP']['Ensemble Mean'][:,iMP],'--',color=cle2,linewidth=lw1)
        ax[3,1].set(position=[0.37,0.04,aw,ah],xlim=xlim,xticks=xticks,xlabel='',ylabel='Product emissions (MgC/ha)')
    
        ax[3,2].plot(tv,MosByMP[iB]['v1']['Mean']['E_CO2e_AGHGB_WOSub']['Ensemble Mean'][:,iMP],'-',color=cle1,linewidth=lw1)
        ax[3,2].plot(tv,MosByMP[iP]['v1']['Mean']['E_CO2e_AGHGB_WOSub']['Ensemble Mean'][:,iMP],'--',color=cle2,linewidth=lw1)
        ax[3,2].set(position=[0.71,0.04,aw,ah],xlim=xlim,xticks=xticks,xlabel='',ylabel='Sector GHG balance (MgC/ha)')
    
        #gu.axletters(ax,plt,0.01,0.91)
        
        try:
            pt=meta['LUT']['ProjectType'][meta['ProjectTypeByMP'][uMP[iMP]]]
            gu.PrintFig(meta['Paths']['Figures'] + '\\BySparseGridSample\\MP' + str(uMP[iMP]) + '_' + pt,'png',200)
        except:
            gu.PrintFig(meta['Paths']['Figures'] + '\\BySparseGridSample\\MP' + str(uMP[iMP]),'png',200)
    
    return

#%% QA - Plot time series by PP number

def QA_PlotTimeSeries_ByPP(meta,dMP,MosByMP,tv,iB,iP):
    
    List=[]
    for iPr in range(dMP['FIA_PROJECT_ID'].size):
        
        if (np.isin(dMP['FIA_PROJECT_ID'][iPr],List)==True):
            continue
        
        # Find multipolygons for tis project
        if (dMP['FIA_PROJECT_ID'][iPr]!=''):
            indP=np.where(uP_ByMP==dMP['FIA_PROJECT_ID'][iPr])[0]
            nam=dMP['FIA_PROJECT_ID'][iPr]
            List.append(dMP['FIA_PROJECT_ID'][iPr])
        else:
            indP=np.where(uO_ByMP==dMP['OPENING_ID'][iPr])[0]
            nam=int(dMP['OPENING_ID'][iPr])
        
        lw1=1; cle1=[0,0,1]; cle2=[1,0,0]; ms=3; aw=0.28; ah=0.22;
        xlim=[tv[it[0]],tv[it[-1]]]; xticks=np.arange(1500,2200,20);
    
        plt.close('all'); fig,ax=plt.subplots(3,2,figsize=gu.cm2inch(28,20))  

        ax[0,0].plot(tv,np.mean(MosByMP[iB]['v1']['Mean']['A']['Ensemble Mean'][:,indP],axis=1),'-',color=cle1,linewidth=lw1)
        ax[0,0].plot(tv,np.mean(MosByMP[iP]['v1']['Mean']['A']['Ensemble Mean'][:,indP],axis=1),'--',color=cle2,linewidth=lw1)
        ax[0,0].set(xlim=xlim,xticks=xticks,xlabel='',ylabel='Age, years')
    
        ax[0,1].plot(tv,np.mean(MosByMP[iB]['v1']['Mean']['C_Biomass_Tot']['Ensemble Mean'][:,indP],axis=1),'-',color=cle1,linewidth=lw1)
        ax[0,1].plot(tv,np.mean(MosByMP[iP]['v1']['Mean']['C_Biomass_Tot']['Ensemble Mean'][:,indP],axis=1),'--',color=cle2,linewidth=lw1)
        ax[0,1].set(xlim=xlim,xlabel='',ylabel='Biomass (MgC/ha)')
    
        ax[1,0].plot(tv,np.mean(MosByMP[iB]['v1']['Mean']['E_CO2e_LULUCF_Wildfire']['Ensemble Mean'][:,indP],axis=1),'-',color=cle1,linewidth=lw1)
        ax[1,0].plot(tv,np.mean(MosByMP[iP]['v1']['Mean']['E_CO2e_LULUCF_Wildfire']['Ensemble Mean'][:,indP],axis=1),'--',color=cle2,linewidth=lw1)
        ax[1,0].set(xlim=xlim,xlabel='',ylabel='Wildfire emissions (tCO2e/ha/yr)')
        
        ax[1,1].plot(tv,np.mean(MosByMP[iB]['v1']['Mean']['C_ToMill']['Ensemble Mean'][:,indP],axis=1),'-',color=cle1,linewidth=lw1)
        ax[1,1].plot(tv,np.mean(MosByMP[iP]['v1']['Mean']['C_ToMill']['Ensemble Mean'][:,indP],axis=1),'--',color=cle2,linewidth=lw1)
        ax[1,1].set(xlim=xlim,xlabel='',ylabel='Removals (tCO2e/ha/yr)')
    
        ax[2,0].plot(tv,np.mean(MosByMP[iB]['v1']['Mean']['E_CO2e_AGHGB_WOSub']['Ensemble Mean'][:,indP],axis=1),'-',color=cle1,linewidth=lw1)
        ax[2,0].plot(tv,np.mean(MosByMP[iP]['v1']['Mean']['E_CO2e_AGHGB_WOSub']['Ensemble Mean'][:,indP],axis=1),'--',color=cle2,linewidth=lw1)
        ax[2,0].set(xlim=xlim,xlabel='',ylabel='AGHGB (tCO2e/ha/yr)')
        
        d=np.mean(MosByMP[iP]['v1']['Mean']['E_CO2e_AGHGB_WOSub']['Ensemble Mean'][:,indP],axis=1)-np.mean(MosByMP[iB]['v1']['Mean']['E_CO2e_AGHGB_WOSub']['Ensemble Mean'][:,indP],axis=1)
        ax[2,1].plot(tv,d,'g-',linewidth=lw1)
        ax[2,1].set(xlim=xlim,xlabel='',ylabel='$\Delta$ AGHGB (tCO2e/ha/yr)')
        plt.tight_layout()
        gu.PrintFig(meta['Paths']['Figures'] + '\\ByPPNumber\\Project_' + str(nam),'png',200)

    return

#%% Export event summary by MP

def ExportSummaryEvents_ByMP():

    n=np.array(1e5,dtype=int)
    dE={}
    dE['ID Multipolygon']=np.zeros(n)
    dE['Year']=np.zeros(n)
    dE['Baseline']=np.array(['' for _ in range(n)],dtype=object)
    dE['Project']=np.array(['' for _ in range(n)],dtype=object)
    cnt=0
    for iMP in range(uMP.size):
        indS=np.where(sxy['ID_atu_multipolygons']==uMP[iMP])[0]
        Year0=np.array([])
        Type0=np.array([])
        Mort0=np.array([])
        BaseAffected0=np.array([])
        ProjAffected0=np.array([])
        for iS in range(indS.size):
            Year0=np.append(Year0,dmec[indS[iS]]['Year'])
            Mort0=np.append(Mort0,dmec[indS[iS]]['MortalityFactor'])
            Type0=np.append(Type0,dmec[indS[iS]]['ID_Type'])
            BaseAffected0=np.append(BaseAffected0,dmec[indS[iS]]['ScnAffected'][0])
            ProjAffected0=np.append(ProjAffected0,dmec[indS[iS]]['ScnAffected'][1])
    
        u=np.unique(np.column_stack([Year0,Type0,Mort0,BaseAffected0,ProjAffected0]),axis=0)
    
        for iE in range(u.shape[0]):
            dE['ID Multipolygon'][cnt]=uMP[iMP]
            dE['Year'][cnt]=u[iE,0]
            if u[iE,3]==1:
                dE['Baseline'][cnt]=cbu.lut_n2s(meta['LUT']['Dist'],u[iE,1])[0] + '-' + str(int(u[iE,2]))
            if u[iE,4]==1:
                dE['Project'][cnt]=cbu.lut_n2s(meta['LUT']['Dist'],u[iE,1])[0] + '-' + str(int(u[iE,2]))
            cnt=cnt+1

    for k in dE.keys():
        dE[k]=dE[k][0:cnt-1]

    df=pd.DataFrame(dE)
    df.to_excel(meta['Paths']['Project'] + '\\Inputs\\SummaryEvents_ByMP.xlsx',index=False)
    
    return

#%% Export summary events by SXY

def ExportSummaryEvents_BySXY():

    n=np.array(2e6,dtype=int)
    dE={}
    dE['ID SXY']=np.zeros(n)
    dE['ID Multipolygon']=np.zeros(n)
    dE['OPENING_ID']=np.zeros(n)
    dE['Year']=np.zeros(n)
    dE['Baseline']=np.array(['' for _ in range(n)],dtype=object)
    dE['Project']=np.array(['' for _ in range(n)],dtype=object)
    cnt=0
    for iMP in range(uMP.size):
        indS=np.where(sxy['ID_atu_multipolygons']==uMP[iMP])[0]
        for iS in range(indS.size):
            Year0=dmec[indS[iS]]['Year']
            Mort0=dmec[indS[iS]]['MortalityFactor']
            Type0=dmec[indS[iS]]['ID_Type']
            BaseAffected0=dmec[indS[iS]]['ScnAffected'][0]
            ProjAffected0=dmec[indS[iS]]['ScnAffected'][1]
            u=np.unique(np.column_stack([Year0,Type0,Mort0,BaseAffected0,ProjAffected0]),axis=0)
    
            for iE in range(u.shape[0]):
                dE['ID SXY'][cnt]=indS[iS]
                dE['ID Multipolygon'][cnt]=uMP[iMP]
                dE['OPENING_ID'][cnt]=atu_multipolygons[uMP[iMP]]['OPENING_ID']
                dE['Year'][cnt]=u[iE,0]
                if u[iE,3]==1:
                    dE['Baseline'][cnt]=cbu.lut_n2s(meta['LUT']['Dist'],u[iE,1])[0] + '-' + str(int(u[iE,2]))
                if u[iE,4]==1:
                    dE['Project'][cnt]=cbu.lut_n2s(meta['LUT']['Dist'],u[iE,1])[0] + '-' + str(int(u[iE,2]))
                cnt=cnt+1

    for k in dE.keys():
        dE[k]=dE[k][0:cnt-1]    

    df=pd.DataFrame(dE)
    df.to_excel(meta['Paths']['Project'] + '\\Inputs\\SummaryEvents_BySXY.xlsx',index=False)  
    
    return

#%% Export area-weighted average conditions by Project Type

def ExportSummaryAttributes_AreaWeighted_ByActivityType(meta,dMP):

    uAT=np.unique(dMP['Project Type'])

    dS={}
    dS['Project Type']=uAT
    dS['SI coast']=np.nan*np.ones(uAT.size)
    dS['SPH coast']=np.nan*np.ones(uAT.size)
    dS['SI interior']=np.nan*np.ones(uAT.size)
    dS['SPH interior']=np.nan*np.ones(uAT.size)

    for i in range(uAT.size):
        
        ind=np.where( (dMP['Project Type']==uAT[i]) & (dMP['BGC Zone 1']=='CWH') & (dMP['SI_ba_mean']>0) | (dMP['Project Type']==uAT[i]) & (dMP['BGC Zone 1']=='ICH') & (dMP['SI_ba_mean']>0) )[0]
        if ind.size>0:
            dS['SI coast'][i]=np.sum(dMP['SI_ba_mean'][ind]*dMP['ACTUAL_TREATMENT_AREA'][ind])/np.sum(dMP['ACTUAL_TREATMENT_AREA'][ind])
        
        ind=np.where( (dMP['Project Type']==uAT[i]) & (dMP['BGC Zone 1']=='CWH') & (dMP['ACTUAL_PLANTED_NUMBER']>0) | (dMP['Project Type']==uAT[i]) & (dMP['BGC Zone 1']=='ICH') & (dMP['ACTUAL_PLANTED_NUMBER']>0) )[0]
        if ind.size>0:
            dS['SPH coast'][i]=np.sum(dMP['ACTUAL_PLANTED_NUMBER'][ind]*dMP['ACTUAL_TREATMENT_AREA'][ind])/np.sum(dMP['ACTUAL_TREATMENT_AREA'][ind])
        
        ind=np.where( (dMP['Project Type']==uAT[i]) & (dMP['BGC Zone 1']!='CWH') & (dMP['BGC Zone 1']!='ICH') & (dMP['SI_ba_mean']>0) )[0]
        if ind.size>0:
            dS['SI interior'][i]=np.sum(dMP['SI_ba_mean'][ind]*dMP['ACTUAL_TREATMENT_AREA'][ind])/np.sum(dMP['ACTUAL_TREATMENT_AREA'][ind])
        
        ind=np.where( (dMP['Project Type']==uAT[i]) & (dMP['BGC Zone 1']!='CWH') & (dMP['BGC Zone 1']!='ICH') & (dMP['ACTUAL_PLANTED_NUMBER']>0) )[0]
        if ind.size>0:            
            dS['SPH interior'][i]=np.sum(dMP['ACTUAL_PLANTED_NUMBER'][ind]*dMP['ACTUAL_TREATMENT_AREA'][ind])/np.sum(dMP['ACTUAL_TREATMENT_AREA'][ind])

    df=pd.DataFrame(dS)
    df.to_excel(meta['Paths']['Project'] + '\\Inputs\\SummaryAttributes_AreaWeighted_ByActivityType.xlsx',index=False)
    
    return

#%% Summary by BGC
# 68% of completed LCELF funds (excluding FFT) are CWH

def Summary_ByBGC(meta,dMP):
    
    uAT=np.unique(dMP['Project Type'])
    
    for iAT in range(uAT.size):
        
        indAT=np.where(dMP['Project Type']==uAT[iAT])[0]
        
        # Unique BGC classes
        uC=np.unique(np.array([dMP['BGCz'][indAT],dMP['BGCsz'][indAT],dMP['BGCv'][indAT]]).T,axis=0)
        
        for iC in range(uC.size):
            ind=np.where( (dMP['Project Type']==uAT[iAT]) & (dMP['BGCz'][indAT],dMP['BGCsz'][indAT],dMP['BGCv'][indAT]) )[0]
    
    ind=np.where( (dMP['SILV_BASE_CODE']=='FE') )[0]
    uBGC=np.unique(dMP['BGC Zone Mode'][ind])

    A=np.zeros(uBGC.size)
    for i in range(uBGC.size):
        ind=np.where( (dMP['BGC Zone Mode']==uBGC[i]) & (dMP['SILV_BASE_CODE']=='FE') )[0]
        A[i]=np.sum(dMP['ACTUAL_TREATMENT_AREA'][ind])

    ind=np.where( (uBGC=='CWH') )[0]
    A_cwh=np.sum(A[ind])
    A_tot=np.sum(A)
    A_cwh/A_tot*100

    # Area-weighted site index = 26 (May 3, 2021)
    ind=np.where( (dMP['SILV_BASE_CODE']=='FE') )[0]
    SI_w=np.sum(dMP['ACTUAL_TREATMENT_AREA'][ind]*dMP['SI_ba_mean'][ind])/np.sum(dMP['ACTUAL_TREATMENT_AREA'][ind])
    SI_w
    
    return

#%% Summary attributes by SXY (only one row per grid cell)

def ExportSummaryAttributes_BySXY(meta,sxy,atu_multipolygons,vri):
    
    # Get land cover class names
    lcc=gu.ReadExcel(r'C:\Users\rhember\Documents\Code_Python\fcgadgets\cbrunner\Parameters\Parameters_VRI.xlsx','LCC')
    
    d={}
    n_init=int(5e5)
    cnt=0
    
    d['ID_SXY']=-999*np.ones(n_init)
    d['ID_Multipolygon']=-999*np.ones(n_init)
    d['OPENING_ID']=-999*np.ones(n_init)
    d['SILV_POLYGON_NUMBER']=np.array(['' for _ in range(n_init)],dtype=object)
    d['Geom Source']=np.array(['' for _ in range(n_init)],dtype=object)
    d['AEF_ATU']=np.zeros(n_init)
    d['LCC']=np.array(['' for _ in range(n_init)],dtype=object)
    d['BGCz']=np.array(['' for _ in range(n_init)],dtype=object)
    d['BGCsz']=np.array(['' for _ in range(n_init)],dtype=object)
    d['BGCv']=np.array(['' for _ in range(n_init)],dtype=object)
    d['SBC']=np.array(['' for _ in range(n_init)],dtype=object)
    d['STC']=np.array(['' for _ in range(n_init)],dtype=object)
    d['SMC']=np.array(['' for _ in range(n_init)],dtype=object)
    
    for iSXY in range(sxy['x'].size):        
        d['ID_SXY'][cnt]=iSXY
        d['ID_Multipolygon'][cnt]=sxy['ID_atu_multipolygons'][iSXY]
        d['OPENING_ID'][cnt]=sxy['OPENING_ID'][iSXY]
        cnt=cnt+1
    
    # Truncate dictionary  
    for k in d.keys():
        d[k]=d[k][0:cnt]
    
    # Add info from atu_multipolygon 
    uMP=np.unique(d['ID_Multipolygon']).astype(int)
    for iMP in range(uMP.size):        
        indD=np.where(d['ID_Multipolygon']==uMP[iMP])[0]        
        if atu_multipolygons[uMP[iMP]]['GeomFromOpLyr']==1:
            d['Geom Source'][indD]='OP'
        elif atu_multipolygons[uMP[iMP]]['GeomFromFcLyr']==1:
            d['Geom Source'][indD]='FC ART'
        else:
            d['Geom Source'][indD]='AT'        
        A=atu_multipolygons[uMP[iMP]]['ACTUAL_TREATMENT_AREA']
        if A!=None:
            #indS=np.where(sxy['ID_atu_multipolygons']==uMP[iMP])[0]
            d['AEF_ATU'][indD]=np.round(A/indD.size,3)
        d['SBC'][indD]=atu_multipolygons[uMP[iMP]]['SILV_BASE_CODE']
        d['STC'][indD]=atu_multipolygons[uMP[iMP]]['SILV_TECHNIQUE_CODE']
        d['SMC'][indD]=atu_multipolygons[uMP[iMP]]['SILV_METHOD_CODE']

    for iD in range(d['ID_SXY'].size):
        indVRI=np.where(vri['IdxToSXY']==d['ID_SXY'][iD])[0]
        for iVRI in range(indVRI.size):
            ind=indVRI[iVRI]
            cd=cbu.lut_n2s(meta['LUT']['VRI']['LAND_COVER_CLASS_CD_1'],vri['LAND_COVER_CLASS_CD_1'][ind])[0]
            ind1=np.where(lcc['Code']==cd)[0]
            if ind1.size>0:
                d['LCC'][iD]=lcc['Name'][ind1[0]]
            d['BGCz'][iD]=cbu.lut_n2s(meta['LUT']['VRI']['BEC_ZONE_CODE'],vri['BEC_ZONE_CODE'][ind])[0]
            d['BGCsz'][iD]=cbu.lut_n2s(meta['LUT']['VRI']['BEC_SUBZONE'],vri['BEC_SUBZONE'][ind])[0]
            d['BGCv'][iD]=str(cbu.lut_n2s(meta['LUT']['VRI']['BEC_VARIANT'],vri['BEC_VARIANT'][ind])[0])
    
    # Add dictonary to dataframe
    df=pd.DataFrame(d)

    # Remove -999 (setting to nan will appear as blank)
    df[df==-999]=np.nan
    df[df==-9999]=np.nan
    df[df=='Unidentified']=np.nan

    # Save
    path=meta['Paths']['Project'] + '\\Inputs\\SummaryAttributes_BySXY.xlsx'
    df.to_excel(path,index=False)
        
    return d

#%% Export attribute summary by BGC

def ExportAttributes_ByBGC(meta):
      
    # Import attribute summary by sparse xy coords
    path=meta['Paths']['Project'] + '\\Inputs\\SummaryAttributes_BySXY.xlsx'
    dS=gu.ReadExcel(path)
    
    # Unique BGC classes
    u=np.unique(np.array([dS['BGCz'],dS['BGCsz'],dS['BGCv']]).T,axis=0)

    # Tabulate areas
    A_PL=np.zeros(u.shape[0])
    A_FE=np.zeros(u.shape[0])
    A_SP=np.zeros(u.shape[0])
    A_PC=np.zeros(u.shape[0])
    for i in range(u.shape[0]):
        ind=np.where( (dS['SBC']=='PL') & (dS['BGCz']==u[i,0]) & (dS['BGCsz']==u[i,1]) & (dS['BGCv']==np.array(u[i,2],dtype=float)) )[0]
        A_PL[i]=np.round(np.sum(dS['AEF_ATU'][ind]),2)
        ind=np.where( (dS['SBC']=='FE') & (dS['BGCz']==u[i,0]) & (dS['BGCsz']==u[i,1]) & (dS['BGCv']==np.array(u[i,2],dtype=float)) )[0]    
        A_FE[i]=np.round(np.sum(dS['AEF_ATU'][ind]),2)
        ind=np.where( (dS['SBC']=='SP') & (dS['BGCz']==u[i,0]) & (dS['BGCsz']==u[i,1]) & (dS['BGCv']==np.array(u[i,2],dtype=float)) )[0]    
        A_SP[i]=np.round(np.sum(dS['AEF_ATU'][ind]),2)
        ind=np.where( (dS['SBC']=='PC') & (dS['BGCz']==u[i,0]) & (dS['BGCsz']==u[i,1]) & (dS['BGCv']==np.array(u[i,2],dtype=float)) )[0]    
        A_PC[i]=np.round(np.sum(dS['AEF_ATU'][ind]),2)
    
    sts=np.column_stack([u,A_PL,A_SP,A_PC,A_FE])
    
    df=pd.DataFrame(data=sts,columns=['Zone','Subzone','Variant','Area PL','Area SP','Area PC','Area FE'])
    
    path=meta['Paths']['Project'] + '\\Inputs\\SummaryActivities_ByBGC.xlsx'
    df.to_excel(path,index=False)
    
    return

#%% Add Archive FC to FC Inventory dictionary

# Instructions to ensure all rows are retained when downloading large DBF files in ArcGIS:
# Properties -> Query Builder -> "FOREST_COVER_ID is not Null" -> Verify -> Get unique values -> Ok 
# -> Export to file gdb

# Requirements:
# RSLT_FOREST_COVER
# RSLT_FOREST_COVER_ARCHIVE
# RSLT_FOREST_COVER_LAYER
# RSLT_FOREST_COVER_LYER_ARCHIVE
# RSLT_FOREST_COVER_LYER_SPC_ARC
# RSLT_FORHEALTH_RSLT
# RSLT_FORHEALTH_RSLT_ARCH

def ForestCover_AddArchive(meta,fcinv,fcsilv):

    # Import look-up 
    meta=Load_LUTs(meta)
    
    #--------------------------------------------------------------------------
    # Did the download include all data??
    # 2021-12-20: Seems to be working!
    #--------------------------------------------------------------------------
    flg=0
    if flg==1:
        path=r'C:\Users\rhember\Documents\Data\ForestInventory\Results\20210930\ForestCover.gdb' 
        fiona.listlayers(path)
        
        id=-663760000
        with fiona.open(path,layer='RSLT_FOREST_COVER') as source:
            for feat in source:
                #break
                if feat['properties']['OPENING_ID']==id:
                    #print('Found it!')
                    #print(feat['properties']['ENTRY_TIMESTAMP'])
                    print(feat['properties']['UPDATE_TIMESTAMP'])
                    #print(feat['properties']['REFERENCE_YEAR'])
          
        with fiona.open(path,layer='RSLT_FOREST_COVER_ARCHIVE') as source:
            for feat in source:
                #break
                if feat['properties']['OPENING_ID']==id:
                    #print('Found it!')    
                    #print(feat['properties']['ENTRY_TIMESTAMP'])
                    print(feat['properties']['UPDATE_TIMESTAMP'])
                    print(feat['properties']['FOREST_COVER_ID'])
                    #print(feat['properties']['REFERENCE_YEAR'])
        
        #id=-663760005
        #id=147986
        id=1774790
        with fiona.open(path,layer='RSLT_FOREST_COVER_LYER_ARCHIVE') as source:
            for feat in source:
                #break
                if feat['properties']['FOREST_COVER_ID']==id:
                    #print('Found it!')    
                    #print(feat['properties']['ENTRY_TIMESTAMP'])
                    print(feat['properties']['UPDATE_TIMESTAMP'])
                    #print(feat['properties']['REFERENCE_YEAR'])

    
    path=meta['Paths']['Results'] + '\\ForestCover.gdb'    
    #cnt2=0
    #with fiona.open(path,layer='RSLT_FOREST_COVER_LYER_ARCHIVE') as source:
    #    for feat in source:
    #        cnt2=cnt2+1
    
    #print(cnt1)
    #print(cnt2)
    
    #fiona.listlayers(path)

    #--------------------------------------------------------------------------
    # Get FC info
    #--------------------------------------------------------------------------
    
    n=int(5e6)
    dFC={}
    dFC['FOREST_COVER_ID']=np.zeros(n)
    dFC['OPENING_ID']=np.zeros(n)
    dFC['SILV_POLYGON_NO']=np.array(['' for _ in range(n)],dtype=object)
    dFC['SILV_POLYGON_AREA']=np.zeros(n)
    dFC['REFERENCE_YEAR']=np.zeros(n)
    dFC['STOCKING_TYPE_CODE']=np.array(['' for _ in range(n)],dtype=object)
    dFC['STOCKING_STATUS_CODE']=np.array(['' for _ in range(n)],dtype=object)
    dFC['ENTRY_YEAR']=np.zeros(n)
    dFC['ENTRY_MONTH']=np.zeros(n)
    dFC['UPDATE_YEAR']=np.zeros(n)
    dFC['UPDATE_MONTH']=np.zeros(n)
    cnt=0
    with fiona.open(path,layer='RSLT_FOREST_COVER') as source:
        for feat in source:
            prp=feat['properties']
            dFC['FOREST_COVER_ID'][cnt]=prp['FOREST_COVER_ID']
            dFC['OPENING_ID'][cnt]=prp['OPENING_ID']
            dFC['SILV_POLYGON_NO'][cnt]=prp['SILV_POLYGON_NO']
            dFC['SILV_POLYGON_AREA'][cnt]=prp['SILV_POLYGON_AREA']
            dFC['REFERENCE_YEAR'][cnt]=prp['REFERENCE_YEAR']
            dFC['STOCKING_TYPE_CODE'][cnt]=prp['STOCKING_TYPE_CODE']
            dFC['STOCKING_STATUS_CODE'][cnt]=prp['STOCKING_STATUS_CODE']
            dFC['ENTRY_YEAR'][cnt]=int(prp['ENTRY_TIMESTAMP'][0:4])
            dFC['ENTRY_MONTH'][cnt]=int(prp['ENTRY_TIMESTAMP'][5:7])
            dFC['UPDATE_YEAR'][cnt]=int(prp['UPDATE_TIMESTAMP'][0:4])
            dFC['UPDATE_MONTH'][cnt]=int(prp['UPDATE_TIMESTAMP'][5:7])
            cnt=cnt+1       
    print('Forest cover length = ' + str(cnt))
    
    #--------------------------------------------------------------------------
    # Get FC archive info
    #--------------------------------------------------------------------------

    with fiona.open(path,layer='RSLT_FOREST_COVER_ARCHIVE') as source:
        for feat in source:
            prp=feat['properties']
            dFC['FOREST_COVER_ID'][cnt]=prp['FOREST_COVER_ID']
            dFC['OPENING_ID'][cnt]=prp['OPENING_ID']
            dFC['SILV_POLYGON_NO'][cnt]=prp['SILV_POLYGON_NO']
            dFC['SILV_POLYGON_AREA'][cnt]=prp['SILV_POLYGON_AREA']
            dFC['REFERENCE_YEAR'][cnt]=prp['REFERENCE_YEAR']
            dFC['STOCKING_TYPE_CODE'][cnt]=prp['STOCKING_TYPE_CODE']
            dFC['STOCKING_STATUS_CODE'][cnt]=prp['STOCKING_STATUS_CODE']
            dFC['ENTRY_YEAR'][cnt]=int(prp['ENTRY_TIMESTAMP'][0:4])
            dFC['ENTRY_MONTH'][cnt]=int(prp['ENTRY_TIMESTAMP'][5:7])
            dFC['UPDATE_YEAR'][cnt]=int(prp['UPDATE_TIMESTAMP'][0:4])
            dFC['UPDATE_MONTH'][cnt]=int(prp['UPDATE_TIMESTAMP'][5:7])
            cnt=cnt+1
    
    # Truncate
    for k in dFC.keys():
        dFC[k]=dFC[k][0:cnt]
    
    print('Forest cover archive length = ' + str(cnt))
    
    #--------------------------------------------------------------------------
    # Add year and month because REFERENCE YEAR HAS NO MONTH 
    #--------------------------------------------------------------------------
       
    for i in range(fcinv['FOREST_COVER_ID'].size):
        if fcinv['Year_Created'][i]>0:
            continue
        ind=np.where(dFC['FOREST_COVER_ID']==fcinv['FOREST_COVER_ID'][i])[0]
        if ind.size>0:
            fcinv['Year_Created'][i]=dFC['ENTRY_YEAR'][ind]
            fcinv['Month_Created'][i]=dFC['ENTRY_MONTH'][ind]
            fcinv['Year_Updated'][i]=dFC['UPDATE_YEAR'][ind]
            fcinv['Month_Updated'][i]=dFC['UPDATE_MONTH'][ind]
       
    for i in range(fcsilv['FOREST_COVER_ID'].size):
        if fcsilv['Year_Created'][i]>0:
            continue
        ind=np.where(dFC['FOREST_COVER_ID']==fcsilv['FOREST_COVER_ID'][i])[0]
        if ind.size>0:
            fcsilv['Year_Created'][i]=dFC['ENTRY_YEAR'][ind]
            fcsilv['Month_Created'][i]=dFC['ENTRY_MONTH'][ind]
            fcsilv['Year_Updated'][i]=dFC['UPDATE_YEAR'][ind]
            fcsilv['Month_Updated'][i]=dFC['UPDATE_MONTH'][ind]
        
    #--------------------------------------------------------------------------
    # Get FC layer info
    #--------------------------------------------------------------------------
    
    n=int(5e6)
    dFCL={}
    dFCL['FOREST_COVER_ID']=np.zeros(n)
    dFCL['FOREST_COVER_LAYER_ID']=np.zeros(n)
    dFCL['CROWN_CLOSURE_PCT']=np.zeros(n)
    dFCL['FCLA_TOTAL_STEMS_PER_HA']=np.zeros(n)    
    dFCL['TOTAL_WELL_SPACED_STEMS_P']=np.zeros(n)
    dFCL['WELL_SPACED_STEMS_PER_HA']=np.zeros(n)
    dFCL['FREE_GROWING_STEMS_PER_HA']=np.zeros(n)
    dFCL['FOREST_COVER_LAYER_CODE']=np.array(['' for _ in range(n)],dtype=object)  
    dFCL['STOCKING_STATUS_CODE']=np.array(['' for _ in range(n)],dtype=object)  
    dFCL['STOCKING_TYPE_CODE']=np.array(['' for _ in range(n)],dtype=object)  
    #dFCL['ARCHIVE_YEAR']=np.zeros(n)
    dFCL['ENTRY_TIME']=np.zeros(n)
    dFCL['UPDATE_TIME']=np.zeros(n)
    cnt=0
    with fiona.open(path,layer='RSLT_FOREST_COVER_LAYER') as source:
        for feat in source:
            prp=feat['properties']
            dFCL['FOREST_COVER_ID'][cnt]=prp['FOREST_COVER_ID']
            dFCL['FOREST_COVER_LAYER_ID'][cnt]=prp['FOREST_COVER_LAYER_ID']
            #dFCL['ARCHIVE_YEAR'][cnt]=int(prp['ARCHIVE_DA'][0:4])
            dFCL['CROWN_CLOSURE_PCT'][cnt]=prp['CROWN_CLOSURE_PCT']
            dFCL['FCLA_TOTAL_STEMS_PER_HA'][cnt]=prp['TOTAL_STEMS_PER_HA']                        
            dFCL['TOTAL_WELL_SPACED_STEMS_P'][cnt]=prp['FCLR_TOTAL_WELL_SPACED_STEMS_P']
            dFCL['WELL_SPACED_STEMS_PER_HA'][cnt]=prp['WELL_SPACED_STEMS_PER_HA']
            dFCL['FREE_GROWING_STEMS_PER_HA'][cnt]=prp['FREE_GROWING_STEMS_PER_HA']
            dFCL['FOREST_COVER_LAYER_CODE'][cnt]=prp['FOREST_COVER_LAYER_CODE']        
            dFCL['STOCKING_STATUS_CODE'][cnt]=prp['STOCKING_STATUS_CODE']
            dFCL['STOCKING_TYPE_CODE'][cnt]=prp['STOCKING_TYPE_CODE']
            dFCL['ENTRY_TIME'][cnt]=int(prp['ENTRY_TIMESTAMP'][0:4])
            dFCL['UPDATE_TIME'][cnt]=int(prp['UPDATE_TIMESTAMP'][0:4])
            cnt=cnt+1
    
    with fiona.open(path,layer='RSLT_FOREST_COVER_LYER_ARCHIVE') as source:
        for feat in source:
            prp=feat['properties']
            dFCL['FOREST_COVER_ID'][cnt]=prp['FOREST_COVER_ID']
            dFCL['FOREST_COVER_LAYER_ID'][cnt]=prp['FOREST_COVER_LAYER_ID']
            #dFCL['ARCHIVE_YEAR'][cnt]=int(prp['ARCHIVE_DA'][0:4])
            dFCL['CROWN_CLOSURE_PCT'][cnt]=prp['CROWN_CLOSURE_PCT']
            dFCL['FCLA_TOTAL_STEMS_PER_HA'][cnt]=prp['FCLA_TOTAL_STEMS_PER_HA']
            dFCL['TOTAL_WELL_SPACED_STEMS_P'][cnt]=prp['TOTAL_WELL_SPACED_STEMS_P']
            dFCL['WELL_SPACED_STEMS_PER_HA'][cnt]=prp['WELL_SPACED_STEMS_PER_HA']
            dFCL['FREE_GROWING_STEMS_PER_HA'][cnt]=prp['FREE_GROWING_STEMS_PER_HA']
            dFCL['FOREST_COVER_LAYER_CODE'][cnt]=prp['FOREST_COVER_LAYER_CODE']        
            dFCL['STOCKING_STATUS_CODE'][cnt]=prp['STOCKING_STATUS_CODE']
            dFCL['STOCKING_TYPE_CODE'][cnt]=prp['STOCKING_TYPE_CODE']
            dFCL['ENTRY_TIME'][cnt]=int(prp['ENTRY_TIMESTAMP'][0:4])
            dFCL['UPDATE_TIME'][cnt]=int(prp['UPDATE_TIMESTAMP'][0:4])
            cnt=cnt+1
    
    for k in dFCL.keys():
        dFCL[k]=dFCL[k][0:cnt]
    
    print('Forest cover layer archive length = ' + str(cnt))
    
    #--------------------------------------------------------------------------
    # Get FC layer archive species info
    #--------------------------------------------------------------------------
    
    n=int(5e6)
    dFC_Spc={}
    dFC_Spc['FOREST_COVER_ID']=np.zeros(n)
    dFC_Spc['FOREST_COVER_LAYER_ID']=np.zeros(n)
    #dFC_Spc['ARCHIVE_YEAR']=np.zeros(n)
    dFC_Spc['TREE_SPECIES_CODE']=np.array(['' for _ in range(n)],dtype=object)
    dFC_Spc['TREE_SPECIES_PCT']=np.zeros(n)
    #dFC_Spc['SPECIES_ORDER']=np.zeros(n)
    cnt=0
    with fiona.open(path,layer='RSLT_FOREST_COVER_LYER_SPC_ARC') as source:
        for feat in source:
            prp=feat['properties']
            dFC_Spc['FOREST_COVER_ID'][cnt]=prp['FOREST_COVER_ID']
            dFC_Spc['FOREST_COVER_LAYER_ID'][cnt]=prp['FOREST_COVER_LAYER_ID']
            #dFC_Spc['ARCHIVE_YEAR'][cnt]=int(prp['ARCHIVE_DATE'][0:4])
            dFC_Spc['TREE_SPECIES_CODE'][cnt]=prp['TREE_SPECIES_CODE']
            dFC_Spc['TREE_SPECIES_PCT'][cnt]=prp['TREE_SPECIES_PCT']
            #dFC_Spc['SPECIES_ORDER'][cnt]=prp['SPECIES_OR']
            cnt=cnt+1
    
    for k in dFC_Spc.keys():
        dFC_Spc[k]=dFC_Spc[k][0:cnt]
    
    #--------------------------------------------------------------------------
    # Find forest cover for each opening (Inventory)
    #--------------------------------------------------------------------------
    
    uO=np.unique(np.column_stack([fcinv['OPENING_ID'],fcinv['SILV_POLYGON_NUMBER']]),axis=0)
    
    # QA:
    #ind0=np.where(dFC['FOREST_COVER_ID']==767135)[0]
   
    #ind0=np.where(dFCL['FOREST_COVER_ID']==767135)[0]
    
    #pn=meta['LUT']['FC_I']['SILV_POLYGON_NUMBER'][dFC['SILV_POLYGON_NO'][ind0][0]]
    #ind2=np.where( (uO[:,0]==dFC['OPENING_ID'][ind0]) & (uO[:,1]==pn) )[0]
    
    # ind=np.where(atu['OPENING_ID']==np.float32(15241) )[0]
    # ind=np.where(fcinv['OPENING_ID']==np.float32(15241) )[0]
    #fcinv['REFERENCE_YEAR'][ind]
    
    #ind0=np.where(fcinv['FOREST_COVER_ID']==767135)[0]

    n=int(1e6)
    dToAdd={}
    dToAdd['IdxToSXY']=np.zeros(n,dtype=np.int32)
    dToAdd['REFERENCE_YEAR']=np.zeros(n,dtype=np.float32)
    dToAdd['FOREST_COVER_ID']=np.zeros(n,dtype=np.float32)
    dToAdd['SILV_POLYGON_NUMBER']=np.zeros(n,dtype=np.float32)
    dToAdd['SILV_POLYGON_AREA']=np.zeros(n,dtype=np.float32)
    dToAdd['STOCKING_STATUS_CODE']=-9999*np.ones(n)
    dToAdd['STOCKING_TYPE_CODE']=-9999*np.ones(n)
    dToAdd['I_TOTAL_STEMS_PER_HA']=-9999*np.ones(n)
    dToAdd['I_TOTAL_WELL_SPACED_STEMS_HA']=-9999*np.ones(n)
    dToAdd['I_FREE_GROWING_STEMS_PER_HA']=-9999*np.ones(n)
    dToAdd['I_CROWN_CLOSURE_PERCENT']=-9999*np.ones(n)
    dToAdd['I_SPECIES_CODE_1']=-9999*np.ones(n)
    dToAdd['I_SPECIES_CODE_2']=-9999*np.ones(n)
    dToAdd['I_SPECIES_CODE_3']=-9999*np.ones(n)
    dToAdd['I_SPECIES_CODE_4']=-9999*np.ones(n)
    dToAdd['I_SPECIES_CODE_5']=-9999*np.ones(n)
    dToAdd['I_SPECIES_PERCENT_1']=-9999*np.ones(n)
    dToAdd['I_SPECIES_PERCENT_2']=-9999*np.ones(n)
    dToAdd['I_SPECIES_PERCENT_3']=-9999*np.ones(n)
    dToAdd['I_SPECIES_PERCENT_4']=-9999*np.ones(n)
    dToAdd['I_SPECIES_PERCENT_5']=-9999*np.ones(n)
    dToAdd['Year_Created']=-9999*np.ones(n)
    dToAdd['Month_Created']=-9999*np.ones(n)
    dToAdd['Year_Updated']=-9999*np.ones(n)
    dToAdd['Month_Updated']=-9999*np.ones(n)
    cnt=0
    for iO in range(uO.shape[0]):
    
        if (uO[iO,0]==-9999):
            continue
        
        # Index to SXY from existing dictionaries
        IdxToSXY=np.array([])
        ind=np.where( (fcinv['OPENING_ID']==uO[iO,0]) & (fcinv['SILV_POLYGON_NUMBER']==uO[iO,1]) )[0]
        IdxToSXY=np.append(IdxToSXY,fcinv['IdxToSXY'][ind])
        #ind=np.where(atu['OPENING_ID']==uO[iO])[0]
        #IdxToSXY=np.append(IdxToSXY,atu['IdxToSXY'][ind])
        #IdxToSXY=np.unique(IdxToSXY)
    
        if IdxToSXY.size==0:
            continue
    
        #iO=np.where(uO[:,0]==15241)[0]
        # iO=iO[0]
        #cbu.lut_n2s(meta['LUT']['FC_I']['SILV_POLYGON_NUMBER'],uO[iO[0],1])[0]
    
        # Silv polygon number
        cd=cbu.lut_n2s(meta['LUT']['FC_I']['SILV_POLYGON_NUMBER'],uO[iO,1])[0]
    
        # Index to forest cover archive        
        indFC=np.where( (dFC['OPENING_ID']==uO[iO,0]) & (dFC['SILV_POLYGON_NO']==cd) )[0]
        # dFC['REFERENCE_YEAR'][indFC]
        # dFC['ENTRY_TIME'][indFC]
        # dFC['UPDATE_TIME'][indFC]
        #dFC['SILV_POLYGON_NO'][indFC]
        
        if indFC.size==0:
            continue

        for iFCA in range(indFC.size):
        
            iFCA0=indFC[iFCA]
            
            # Index to forest cover layer archive
            indFCL=np.where( (dFCL['FOREST_COVER_ID']==dFC['FOREST_COVER_ID'][iFCA0]) & (dFCL['FOREST_COVER_LAYER_CODE']=='I') )[0]
        
            # Index to forest cover species archive
            indSp=np.where( (dFC_Spc['FOREST_COVER_LAYER_ID']==dFCL['FOREST_COVER_LAYER_ID'][indFCL]) )[0]
            
            for iSXY in range(IdxToSXY.size):
                   
                dToAdd['IdxToSXY'][cnt]=IdxToSXY[iSXY]
                dToAdd['REFERENCE_YEAR'][cnt]=dFC['REFERENCE_YEAR'][iFCA0]
                dToAdd['Year_Created'][cnt]=dFC['ENTRY_YEAR'][iFCA0]
                dToAdd['Month_Created'][cnt]=dFC['ENTRY_MONTH'][iFCA0]
                dToAdd['Year_Updated'][cnt]=dFC['UPDATE_YEAR'][iFCA0]
                dToAdd['Month_Updated'][cnt]=dFC['UPDATE_MONTH'][iFCA0]
                dToAdd['FOREST_COVER_ID'][cnt]=dFC['FOREST_COVER_ID'][iFCA0]
                dToAdd['SILV_POLYGON_NUMBER'][cnt]=uO[iO,1]
                dToAdd['SILV_POLYGON_AREA'][cnt]=dFC['SILV_POLYGON_AREA'][iFCA0]
                try:
                    dToAdd['STOCKING_STATUS_CODE'][cnt]=meta['LUT']['FC_I']['STOCKING_STATUS_CODE'][dFC['STOCKING_STATUS_CODE'][iFCA0]][0]
                except:
                    pass
                
                try:
                    dToAdd['STOCKING_TYPE_CODE'][cnt]=meta['LUT']['FC_I']['STOCKING_TYPE_CODE'][dFC['STOCKING_TYPE_CODE'][iFCA0]][0]
                except:
                    pass
                    
                if indFCL.size>0:                    
                    dToAdd['I_TOTAL_STEMS_PER_HA'][cnt]=np.nanmean(dFCL['FCLA_TOTAL_STEMS_PER_HA'][indFCL])
                    dToAdd['I_TOTAL_WELL_SPACED_STEMS_HA'][cnt]=np.nanmean(dFCL['TOTAL_WELL_SPACED_STEMS_P'][indFCL])
                    dToAdd['I_FREE_GROWING_STEMS_PER_HA'][cnt]=np.nanmean(dFCL['FREE_GROWING_STEMS_PER_HA'][indFCL])
                    dToAdd['I_CROWN_CLOSURE_PERCENT'][cnt]=np.nanmean(dFCL['CROWN_CLOSURE_PCT'][indFCL])
                    
                    # Species
                    ord=np.flip(np.argsort(dFC_Spc['TREE_SPECIES_PCT'][indSp]))
                    for iSp in range(np.minimum(4,indSp.size)):
                        cd=dFC_Spc['TREE_SPECIES_CODE'][ indSp[ord[iSp]] ]
                        try:
                            dToAdd['I_SPECIES_CODE_' + str(iSp+1)][cnt]=meta['LUT']['VRI']['SPECIES_CD_1'][cd]
                        except:
                            pass   
                        dToAdd['I_SPECIES_PERCENT_' + str(iSp+1)][cnt]=dFC_Spc['TREE_SPECIES_PCT'][ indSp[ord[iSp]] ]
                    
                cnt=cnt+1
    
    # Truncate    
    for k in dToAdd.keys():
        dToAdd[k]=dToAdd[k][0:cnt]
    
    # Add to fcinv dictionary
    for k in fcinv.keys():
        if k in dToAdd:
            fcinv[k]=np.append(fcinv[k],dToAdd[k])
        else:
            fcinv[k]=np.append(fcinv[k],-9999*np.ones(dToAdd['IdxToSXY'].size))
    
    # Remove duplicates
    ToKeep=np.zeros(fcinv['FOREST_COVER_ID'].size)
    for i in range(fcinv['FOREST_COVER_ID'].size):
        ind=np.where(fcinv['FOREST_COVER_ID']==fcinv['FOREST_COVER_ID'][i])[0]
        if ind.size>1:
            ind2=np.where(fcinv['SILV_POLYGON_NUMBER'][ind]>0)[0]
            if ind2.size>0:
                ToKeep[ind[ind2]]=1
            else:
                ToKeep[i]=1
    ikp=np.where(ToKeep==1)[0]
    for k in fcinv.keys():
        fcinv[k]=fcinv[k][ikp]
    
    print('Done inv')
    
    #--------------------------------------------------------------------------
    # Find forest cover for each opening (Silv label)
    #--------------------------------------------------------------------------
    
    uO=np.unique(np.column_stack([fcsilv['OPENING_ID'],fcsilv['SILV_POLYGON_NUMBER']]),axis=0)

    n=int(1e6)
    dToAdd={}
    dToAdd['IdxToSXY']=np.zeros(n,dtype=np.int32)
    dToAdd['FOREST_COVER_ID']=np.zeros(n,dtype=np.float32)
    dToAdd['SILV_POLYGON_NUMBER']=np.zeros(n,dtype=np.float32)
    dToAdd['SILV_POLYGON_AREA']=np.zeros(n,dtype=np.float32)  
    dToAdd['REFERENCE_YEAR']=np.zeros(n,dtype=np.float32)
    dToAdd['STOCKING_STATUS_CODE']=-9999*np.ones(n)
    dToAdd['STOCKING_TYPE_CODE']=-9999*np.ones(n)
    dToAdd['S_TOTAL_STEMS_PER_HA']=-9999*np.ones(n)
    dToAdd['S_TOTAL_WELL_SPACED_STEMS_HA']=-9999*np.ones(n)
    dToAdd['S_FREE_GROWING_STEMS_PER_HA']=-9999*np.ones(n)
    dToAdd['S_CROWN_CLOSURE_PERCENT']=-9999*np.ones(n)
    dToAdd['S_SPECIES_CODE_1']=-9999*np.ones(n)
    dToAdd['S_SPECIES_CODE_2']=-9999*np.ones(n)
    dToAdd['S_SPECIES_CODE_3']=-9999*np.ones(n)
    dToAdd['S_SPECIES_CODE_4']=-9999*np.ones(n)
    dToAdd['S_SPECIES_CODE_5']=-9999*np.ones(n)
    dToAdd['S_SPECIES_PERCENT_1']=-9999*np.ones(n)
    dToAdd['S_SPECIES_PERCENT_2']=-9999*np.ones(n)
    dToAdd['S_SPECIES_PERCENT_3']=-9999*np.ones(n)
    dToAdd['S_SPECIES_PERCENT_4']=-9999*np.ones(n)
    dToAdd['S_SPECIES_PERCENT_5']=-9999*np.ones(n)
    dToAdd['Year_Created']=-9999*np.ones(n)
    dToAdd['Month_Created']=-9999*np.ones(n)
    dToAdd['Year_Updated']=-9999*np.ones(n)
    dToAdd['Month_Updated']=-9999*np.ones(n)
    cnt=0
    for iO in range(uO.shape[0]):
    
        if (uO[iO,0]==-9999):
            continue
        
        # Index to SXY from existing dictionaries
        IdxToSXY=np.array([])
        ind=np.where( (fcsilv['OPENING_ID']==uO[iO,0]) & (fcsilv['SILV_POLYGON_NUMBER']==uO[iO,1]) )[0]
        IdxToSXY=np.append(IdxToSXY,fcsilv['IdxToSXY'][ind])
    
        if IdxToSXY.size==0:
            continue
    
        # Index to forest cover archive
        #indFC=np.where( dFC['OPENING_ID']==uO[iO] )[0]
        cd=cbu.lut_n2s(meta['LUT']['FC_S']['SILV_POLYGON_NUMBER'],uO[iO,1])[0]        
        indFC=np.where( (dFC['OPENING_ID']==uO[iO,0]) & (dFC['SILV_POLYGON_NO']==cd) )[0]
    
        if indFC.size==0:
            continue

        for iFCA in range(indFC.size):
        
            iFCA0=indFC[iFCA]
            
            # Index to forest cover layer archive
            indFCL=np.where( (dFCL['FOREST_COVER_ID']==dFC['FOREST_COVER_ID'][iFCA0]) & (dFCL['FOREST_COVER_LAYER_CODE']=='S') )[0]
        
            # Index to forest cover species archive
            indSp=np.where( (dFC_Spc['FOREST_COVER_LAYER_ID']==dFCL['FOREST_COVER_LAYER_ID'][indFCL]) )[0]
            
            for iSXY in range(IdxToSXY.size):
                   
                dToAdd['IdxToSXY'][cnt]=IdxToSXY[iSXY]
                dToAdd['REFERENCE_YEAR'][cnt]=dFC['REFERENCE_YEAR'][iFCA0]
                dToAdd['Year_Created'][cnt]=dFC['ENTRY_YEAR'][iFCA0]
                dToAdd['Month_Created'][cnt]=dFC['ENTRY_MONTH'][iFCA0]
                dToAdd['Year_Updated'][cnt]=dFC['UPDATE_YEAR'][iFCA0]
                dToAdd['Month_Updated'][cnt]=dFC['UPDATE_MONTH'][iFCA0]
                dToAdd['FOREST_COVER_ID'][cnt]=dFC['FOREST_COVER_ID'][iFCA0]
                dToAdd['SILV_POLYGON_NUMBER'][cnt]=uO[iO,1]
                dToAdd['SILV_POLYGON_AREA'][cnt]=dFC['SILV_POLYGON_AREA'][iFCA0]
                try:
                    dToAdd['STOCKING_STATUS_CODE'][cnt]=meta['LUT']['FC_S']['STOCKING_STATUS_CODE'][dFC['STOCKING_STATUS_CODE'][iFCA0]][0]
                except:
                    pass
                
                try:
                    dToAdd['STOCKING_TYPE_CODE'][cnt]=meta['LUT']['FC_S']['STOCKING_TYPE_CODE'][dFC['STOCKING_TYPE_CODE'][iFCA0]][0]
                except:
                    pass
                    
                if indFCL.size>0:                    
                    dToAdd['S_TOTAL_STEMS_PER_HA'][cnt]=np.nanmean(dFCL['FCLA_TOTAL_STEMS_PER_HA'][indFCL])
                    dToAdd['S_TOTAL_WELL_SPACED_STEMS_HA'][cnt]=np.nanmean(dFCL['TOTAL_WELL_SPACED_STEMS_P'][indFCL])
                    dToAdd['S_FREE_GROWING_STEMS_PER_HA'][cnt]=np.nanmean(dFCL['FREE_GROWING_STEMS_PER_HA'][indFCL])
                    dToAdd['S_CROWN_CLOSURE_PERCENT'][cnt]=np.nanmean(dFCL['CROWN_CLOSURE_PCT'][indFCL])
                    
                    # Species
                    ord=np.flip(np.argsort(dFC_Spc['TREE_SPECIES_PCT'][indSp]))
                    for iSp in range(np.minimum(4,indSp.size)):
                        cd=dFC_Spc['TREE_SPECIES_CODE'][ indSp[ord[iSp]] ]
                        try:
                            dToAdd['S_SPECIES_CODE_' + str(iSp+1)][cnt]=meta['LUT']['VRI']['SPECIES_CD_1'][cd]
                        except:
                            pass   
                        dToAdd['S_SPECIES_PERCENT_' + str(iSp+1)][cnt]=dFC_Spc['TREE_SPECIES_PCT'][ indSp[ord[iSp]] ]
                    
                cnt=cnt+1
    
    # Truncate    
    for k in dToAdd.keys():
        dToAdd[k]=dToAdd[k][0:cnt]

    # Add to fcsilv dictionary
    for k in fcsilv.keys():
        if k in dToAdd:
            fcsilv[k]=np.append(fcsilv[k],dToAdd[k])
        else:
            fcsilv[k]=np.append(fcsilv[k],-9999*np.ones(dToAdd['IdxToSXY'].size))

    # Remove duplicates    
    ToKeep=np.zeros(fcsilv['FOREST_COVER_ID'].size)
    for i in range(fcsilv['FOREST_COVER_ID'].size):
        ind=np.where(fcsilv['FOREST_COVER_ID']==fcsilv['FOREST_COVER_ID'][i])[0]
        if ind.size>1:
            ind2=np.where(fcsilv['SILV_POLYGON_NUMBER'][ind]>0)[0]
            if ind2.size>0:
                ToKeep[ind[ind2]]=1
            else:
                ToKeep[i]=1
    ikp=np.where(ToKeep==1)[0]
    for k in fcsilv.keys():
        fcsilv[k]=fcsilv[k][ikp]

    return fcinv,fcsilv

#%% Add forest health info to FC layer (run after adding forest cover archive)

def ForestCover_AddForestHealth(meta,fcinv):
    
    # Import look-up 
    meta=Load_LUTs(meta)
    
    # Path to forest cover archive geodatabase
    path=meta['Paths']['Results'] + '\\ForestCover.gdb' 
    #fiona.listlayers(path)

    #--------------------------------------------------------------------------
    # Initialize forest health variables
    #--------------------------------------------------------------------------
    
    fcinv['I_DA1 CD']=-9999*np.ones(fcinv['FOREST_COVER_ID'].size)
    fcinv['I_DA1 PCT']=np.zeros(fcinv['FOREST_COVER_ID'].size)
    fcinv['I_DA2 CD']=-9999*np.ones(fcinv['FOREST_COVER_ID'].size)
    fcinv['I_DA2 PCT']=np.zeros(fcinv['FOREST_COVER_ID'].size)
    fcinv['I_DA3 CD']=-9999*np.ones(fcinv['FOREST_COVER_ID'].size)
    fcinv['I_DA3 PCT']=np.zeros(fcinv['FOREST_COVER_ID'].size)
    
    fcinv['S_DA1 CD']=-9999*np.ones(fcinv['FOREST_COVER_ID'].size)
    fcinv['S_DA1 PCT']=np.zeros(fcinv['FOREST_COVER_ID'].size)
    fcinv['S_DA2 CD']=-9999*np.ones(fcinv['FOREST_COVER_ID'].size)
    fcinv['S_DA2 PCT']=np.zeros(fcinv['FOREST_COVER_ID'].size)
    fcinv['S_DA3 CD']=-9999*np.ones(fcinv['FOREST_COVER_ID'].size)
    fcinv['S_DA3 PCT']=np.zeros(fcinv['FOREST_COVER_ID'].size)

    #--------------------------------------------------------------------------
    # Add current forest health
    #--------------------------------------------------------------------------

    with fiona.open(path,layer='RSLT_FORHEALTH_RSLT') as source:
        for feat in source:
            prp=feat['properties']

            ind=np.where(fcinv['FOREST_COVER_ID']==prp['FOREST_COVER_ID'])[0]
            if ind.size==0:
                continue
          
            try:
                id=meta['LUT']['Pest']['PEST_SPECIES_CODE'][prp['SILV_DAMAGE_AGENT_CODE']][0]
            except:
                continue
            
            indL=np.where(fcinv['I_FOREST_COVER_LAYER_ID']==prp['FOREST_COVER_LAYER_ID'])[0]
            if indL.size>0:
                # Inventory label
                if fcinv['I_DA1 CD'][ind[0]]==-9999:
                    fcinv['I_DA1 CD'][ind]=id
                    fcinv['I_DA1 PCT'][ind]=prp['INCIDENCE_PCT']
                elif fcinv['I_DA2 CD'][ind[0]]==-9999:
                    fcinv['I_DA2 CD'][ind]=id
                    fcinv['I_DA2 PCT'][ind]=prp['INCIDENCE_PCT']
                elif fcinv['I_DA3 CD'][ind[0]]==-9999:
                    fcinv['I_DA3 CD'][ind]=id
                    fcinv['I_DA3 PCT'][ind]=prp['INCIDENCE_PCT']
            else:
                # Silv label
                if fcinv['S_DA1 CD'][ind[0]]==-9999:
                    fcinv['S_DA1 CD'][ind]=id
                    fcinv['S_DA1 PCT'][ind]=prp['INCIDENCE_PCT']
                elif fcinv['S_DA2 CD'][ind[0]]==-9999:
                    fcinv['S_DA2 CD'][ind]=id
                    fcinv['S_DA2 PCT'][ind]=prp['INCIDENCE_PCT']
                elif fcinv['S_DA3 CD'][ind[0]]==-9999:
                    fcinv['S_DA3 CD'][ind]=id
                    fcinv['S_DA3 PCT'][ind]=prp['INCIDENCE_PCT']  
        
    #--------------------------------------------------------------------------
    # Add archive forest health
    #--------------------------------------------------------------------------

    with fiona.open(path,layer='RSLT_FORHEALTH_RSLT_ARCHIVE') as source:
        for feat in source:
            prp=feat['properties']
            ind=np.where(fcinv['FOREST_COVER_ID']==prp['FOREST_COVER_ID'])[0]
            if ind.size==0:
                continue
            
            try:
                id=meta['LUT']['Pest']['PEST_SPECIES_CODE'][prp['SILV_DAMAGE_AGENT_CODE']][0]
            except:
                continue
            
            indL=np.where(fcinv['I_FOREST_COVER_LAYER_ID']==prp['FOREST_COVER_LAYER_ID'])[0]
            if indL.size>0:
                # Inventory label
                if fcinv['I_DA1 CD'][ind[0]]==-9999:
                    fcinv['I_DA1 CD'][ind]=id
                    fcinv['I_DA1 PCT'][ind]=prp['INCIDENCE_PCT']
                elif fcinv['I_DA2 CD'][ind[0]]==-9999:
                    fcinv['I_DA2 CD'][ind]=id
                    fcinv['I_DA2 PCT'][ind]=prp['INCIDENCE_PCT']
                elif fcinv['I_DA3 CD'][ind[0]]==-9999:
                    fcinv['I_DA3 CD'][ind]=id
                    fcinv['I_DA3 PCT'][ind]=prp['INCIDENCE_PCT']
            else:
                # Silv label
                if fcinv['S_DA1 CD'][ind[0]]==-9999:
                    fcinv['S_DA1 CD'][ind]=id
                    fcinv['S_DA1 PCT'][ind]=prp['INCIDENCE_PCT']
                elif fcinv['S_DA2 CD'][ind[0]]==-9999:
                    fcinv['S_DA2 CD'][ind]=id
                    fcinv['S_DA2 PCT'][ind]=prp['INCIDENCE_PCT']
                elif fcinv['S_DA3 CD'][ind[0]]==-9999:
                    fcinv['S_DA3 CD'][ind]=id
                    fcinv['S_DA3 PCT'][ind]=prp['INCIDENCE_PCT']  
    
    return fcinv