
'''

FCGADGETS - FOREST INVENTORY UTILITIES

'''

#%% Import modules

import numpy as np
import pandas as pd
import geopandas as gpd
import fiona
import gdal
import time
import gc as garc
from fcgadgets.utilities import utilities_gis as gis
from fcgadgets.utilities import utilities_general as gu
from fcgadgets.cbrunner import cbrun_utilities as cbu

#%% Define strings that frequently need to be populated with zeros

StringsToFill=['Month','Day','FCI_Funded','FIA_PROJECT_ID','OPENING_ID','ACTUAL_TREATMENT_AREA','ACTUAL_PLANTED_NUMBER', \
               'PL_SPECIES_CD1','PL_SPECIES_PCT1','PL_SPECIES_GW1','PL_SPECIES_CD2','PL_SPECIES_PCT2','PL_SPECIES_GW2', \
               'PL_SPECIES_CD3','PL_SPECIES_PCT3','PL_SPECIES_GW3','PL_SPECIES_CD4','PL_SPECIES_PCT4','PL_SPECIES_GW4', \
               'PL_SPECIES_CD5','PL_SPECIES_PCT5','PL_SPECIES_GW5']

'''============================================================================
UPDATE FOREST INVENTORY DATA
    1) Download geodatabases from BC Data Catelogue (60 min)
        - store files at the paths indicated in "DefineInventoryLayersAndVariables"
    2) Define inventory layers and variable
    3) Build and save LUTs (20 min)
============================================================================'''

# 1) Manually download geodatabases in ArcGIS (if export is not working, try copy and paste)

# 2) Run function to define inventory layers and variables (make sure folder 
#    release dates are updated)
# LayerInfo=invu.DefineInventoryLayersAndVariables()

# 3) Run function to build and save LUTs
# invu.BuildForestInventoryLUTs(LayerInfo)


'''============================================================================
DEFINE INVENTORY LAYERS AND VARIABLES
The "Field List" variable contains touples containing the variable name and 
a flag indicating whether it is string (1) or numberic (0)
============================================================================'''

def DefineInventoryLayersAndVariables():

    # Define paths to geodatabase files
    PathInResultsFull=r'C:\Users\rhember\Documents\Data\ForestInventory\Results\20200430'
    PathInDisturbancesFull=r'C:\Users\rhember\Documents\Data\ForestInventory\Disturbances\20200430'
    PathInVRIFull=r'C:\Users\rhember\Documents\Data\ForestInventory\VRI\20200430'
    PathInLUPFull=r'C:\Users\rhember\Documents\Data\ForestInventory\LandUse\20200706'

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
     ('DENUDATION_1_COMPLETION_DATE',0,'<U20'), \
     ('DENUDATION_2_DISTURBANCE_CODE',1,'int16'), \
     ('DENUDATION_2_SILV_SYSTEM_CODE',1,'int16'), \
     ('DENUDATION_2_SILV_VARIANT_CODE',1,'int16'), \
     ('DENUDATION_2_COMPLETION_DATE',0,'<U20')]
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
     ('SITE_INDEX',0,'float32'), \
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
     ('I_FREE_GROWING_STEMS_PER_HA',0,'float32'), \
     ('I_CROWN_CLOSURE_PERCENT',0,'float32'), \
     ('REFERENCE_YEAR',0,'float32'), \
     ('STOCKING_STATUS_CODE',1,'int16'), \
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
     ('SITE_INDEX',0,'float32'), \
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
     ('S_FREE_GROWING_STEMS_PER_HA',0,'float32'), \
     ('S_CROWN_CLOSURE_PERCENT',0,'float32'), \
     ('REFERENCE_YEAR',0,'float32'), \
     ('STOCKING_STATUS_CODE',1,'int16'), \
     ('FOREST_COVER_WHO_CREATED',1,'int16'), \
     ('FOREST_COVER_WHO_UPDATED',1,'int16'), \
     ('FOREST_COVER_WHEN_CREATED',0,'<U20'), \
     ('FOREST_COVER_WHEN_UPDATED',0,'<U20')]
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
    
    #fiona.listlayers(PathInLUPFull + '\\Landuse.gdb')
    #lyr=fiona.open(PathInLUPFull + '\\Landuse.gdb',layer='TA_PARK_ECORES_PA_SVW')
    #lyr.schema
    
    return LayerInfo

'''============================================================================
BUILD LOOK-UP TABLES FOR THE CODES IN INVENTORY LAYERS
============================================================================'''

def BuildForestInventoryLUTs(LayerInfo):

    #for iLyr in range(6,7):
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


'''============================================================================
EXTRACT NUMERIC TIME VECTORS FROM DATE STRINGS IN RESULTS
============================================================================'''

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

# Old version (may be required for FCI from inventory project)

def ExtractDateStringsFromRESULTS_old(lyr_nam,data):
    
    if (lyr_nam=='RSLT_ACTIVITY_TREATMENT_SVW') | (lyr_nam=='RSLT_PLANTING_SVW'):
        
        for i in range(len(data)):
            if data[i]==None:
                continue
            L=data[i]['OPENING_ID'].size
            data[i]['Year']=-999*np.ones(L,dtype='int16')
            data[i]['Month']=-999*np.ones(L,dtype='int16')
            data[i]['Day']=-999*np.ones(L,dtype='int16')
            for j in range(L):
                ds=data[i]['ATU_COMPLETION_DATE'][j]
                if ds!='-999':
                    data[i]['Year'][j]=int(ds[0:4])
                    data[i]['Month'][j]=int(ds[5:7])
                    data[i]['Day'][j]=int(ds[8:10])
            # Remove original date string
            del data[i]['ATU_COMPLETION_DATE']
            
    elif (lyr_nam=='RSLT_OPENING_SVW'):
        
        for i in range(len(data)):
            if data[i]==None:
                continue
            L=data[i]['OPENING_ID'].size
            
            data[i]['Year_DistStart']=-999*np.ones(L,dtype='int16')
            data[i]['Month_DistStart']=-999*np.ones(L,dtype='int16')
            data[i]['Day_DistStart']=-999*np.ones(L,dtype='int16')
            for j in range(L):
                ds=data[i]['DISTURBANCE_START_DATE'][j]
                if ds!='-999':
                    data[i]['Year_DistStart'][j]=int(ds[0:4])
                    data[i]['Month_DistStart'][j]=int(ds[5:7])
                    data[i]['Day_DistStart'][j]=int(ds[8:10])
            
            data[i]['Year_DistEnd']=-999*np.ones(L,dtype='int16')
            data[i]['Month_DistEnd']=-999*np.ones(L,dtype='int16')
            data[i]['Day_DistEnd']=-999*np.ones(L,dtype='int16')
            for j in range(L):
                ds=data[i]['DISTURBANCE_END_DATE'][j]
                if ds!='-999':
                    data[i]['Year_DistEnd'][j]=int(ds[0:4])
                    data[i]['Month_DistEnd'][j]=int(ds[5:7])
                    data[i]['Day_DistEnd'][j]=int(ds[8:10])
                    
            data[i]['Year_Denu1_Comp']=-999*np.ones(L,dtype='int16')
            data[i]['Month_Denu1_Comp']=-999*np.ones(L,dtype='int16')
            data[i]['Day_Denu1_Comp']=-999*np.ones(L,dtype='int16')
            for j in range(L):
                ds=data[i]['DENUDATION_1_COMPLETION_DATE'][j]
                if ds!='-999':
                    data[i]['Year_Denu1_Comp'][j]=int(ds[0:4])
                    data[i]['Month_Denu1_Comp'][j]=int(ds[5:7])
                    data[i]['Day_Denu1_Comp'][j]=int(ds[8:10])
                    
            data[i]['Year_Denu2_Comp']=-999*np.ones(L,dtype='int16')
            data[i]['Month_Denu2_Comp']=-999*np.ones(L,dtype='int16')
            data[i]['Day_Denu2_Comp']=-999*np.ones(L,dtype='int16')
            for j in range(L):
                ds=data[i]['DENUDATION_2_COMPLETION_DATE'][j]
                if ds!='-999':
                    data[i]['Year_Denu2_Comp'][j]=int(ds[0:4])
                    data[i]['Month_Denu2_Comp'][j]=int(ds[5:7])
                    data[i]['Day_Denu2_Comp'][j]=int(ds[8:10])        
            
            # Remove original date string
            del data[i]['DISTURBANCE_START_DATE']       
            del data[i]['DISTURBANCE_END_DATE']
            del data[i]['DENUDATION_1_COMPLETION_DATE']
            del data[i]['DENUDATION_2_COMPLETION_DATE']
            
    elif (lyr_nam=='RSLT_FOREST_COVER_INV_SVW') | (lyr_nam=='RSLT_FOREST_COVER_SILV_SVW'):
        
        for i in range(len(data)):
            if data[i]==None:
                continue
            L=data[i]['OPENING_ID'].size
            
            data[i]['Year_Created']=-999*np.ones(L,dtype='int16')
            data[i]['Month_Created']=-999*np.ones(L,dtype='int16')
            data[i]['Day_Created']=-999*np.ones(L,dtype='int16')
            for j in range(L):
                ds=data[i]['FOREST_COVER_WHEN_CREATED'][j]
                if ds!='-999':
                    data[i]['Year_Created'][j]=int(ds[0:4])
                    data[i]['Month_Created'][j]=int(ds[5:7])
                    data[i]['Day_Created'][j]=int(ds[8:10])
            
            data[i]['Year_Updated']=-999*np.ones(L,dtype='int16')
            data[i]['Month_Updated']=-999*np.ones(L,dtype='int16')
            data[i]['Day_Updated']=-999*np.ones(L,dtype='int16')
            for j in range(L):
                ds=data[i]['FOREST_COVER_WHEN_UPDATED'][j]
                if ds!='-999':
                    data[i]['Year_Updated'][j]=int(ds[0:4])
                    data[i]['Month_Updated'][j]=int(ds[5:7])
                    data[i]['Day_Updated'][j]=int(ds[8:10])
            
            # Remove original date string
            del data[i]['FOREST_COVER_WHEN_CREATED']       
            del data[i]['FOREST_COVER_WHEN_UPDATED']
            
    elif (lyr_nam=='PROT_HISTORICAL_FIRE_POLYS_SP'):
        
        for i in range(len(data)):
            if data[i]==None:
                continue
            L=data[i]['FIRE_YEAR'].size
            data[i]['Year']=-999*np.ones(L,dtype='int16')
            data[i]['Month']=-999*np.ones(L,dtype='int16')
            data[i]['Day']=-999*np.ones(L,dtype='int16')
            for j in range(L):
                ds=data[i]['FIRE_DATE'][j]
                if ds!='-999':
                    data[i]['Year'][j]=int(ds[0:4])
                    data[i]['Month'][j]=int(ds[5:7])
                    data[i]['Day'][j]=int(ds[8:10])
            # Remove original date string
            del data[i]['FIRE_DATE']
    
    elif (lyr_nam=='VEG_CONSOLIDATED_CUT_BLOCKS_SP'):
        
        for i in range(len(data)):
            if data[i]==None:
                continue
            L=data[i]['HARVEST_YEAR'].size
            data[i]['Year_Start']=-999*np.ones(L,dtype='int16')
            data[i]['Month_Start']=-999*np.ones(L,dtype='int16')
            data[i]['Day_Start']=-999*np.ones(L,dtype='int16')
            data[i]['Year_End']=-999*np.ones(L,dtype='int16')
            data[i]['Month_End']=-999*np.ones(L,dtype='int16')
            data[i]['Day_End']=-999*np.ones(L,dtype='int16')
            for j in range(L):
                ds=data[i]['DISTURBANCE_START_DATE'][j]
                if ds!='-999':
                    data[i]['Year_Start'][j]=int(ds[0:4])
                    data[i]['Month_Start'][j]=int(ds[5:7])
                    data[i]['Day_Start'][j]=int(ds[8:10])
                ds=data[i]['DISTURBANCE_END_DATE'][j]
                if ds!='-999':
                    data[i]['Year_End'][j]=int(ds[0:4])
                    data[i]['Month_End'][j]=int(ds[5:7])
                    data[i]['Day_End'][j]=int(ds[8:10])    
            # Remove original date string
            del data[i]['DISTURBANCE_START_DATE']
            del data[i]['DISTURBANCE_END_DATE']
    
    elif (lyr_nam=='VEG_COMP_LYR_R1_POLY'):
        
        for i in range(len(data)):
            if data[i]==None:
                continue
            L=data[i]['PROJECTED_DATE'].size
            data[i]['Year']=-999*np.ones(L,dtype='int16')
            data[i]['Month']=-999*np.ones(L,dtype='int16')
            data[i]['Day']=-999*np.ones(L,dtype='int16')
            for j in range(L):
                ds=data[i]['PROJECTED_DATE'][j]
                if ds!='-999':
                    data[i]['Year'][j]=int(ds[0:4])
                    data[i]['Month'][j]=int(ds[5:7])
                    data[i]['Day'][j]=int(ds[8:10])
            # Remove original date string
            del data[i]['PROJECTED_DATE']
            
    return data


'''============================================================================
GET MISSING ATU LAYER GEOMETRIES FROM OPENING LAYER
============================================================================'''

def GetMissingATUGeometries(meta):

    Paths={}
    meta['Paths']['Project']=r'C:\Users\rhember\Documents\Data\FCI_Projects\FCI_SparseGrid'
    meta['Paths']['Results']=r'C:\Users\rhember\Documents\Data\ForestInventory\Results\20200430'
    meta['Paths']['Results File']=meta['Paths']['Results'] + '\\Results.gdb'
    
    # Define name of layer
    fiona.listlayers(meta['Paths']['Results File'])
    nam_lyr='RSLT_ACTIVITY_TREATMENT_SVW'
    lyr=fiona.open(meta['Paths']['Results File'],layer=nam_lyr)
    L=len(lyr)
    scma=lyr.schema

    # Create a dictionary to store openings with no AT spatial info
    atu_mis={}
    atu_mis['OPENING_ID']=np.zeros(L)
    atu_mis['SpatialMissing']=np.zeros(L)
    cnt=0
    with fiona.open(meta['Paths']['Results File'],layer=nam_lyr) as source:
        for feat in source:
            atu_mis['OPENING_ID'][cnt]=feat['properties']['OPENING_ID']        
            if feat['geometry']==None:
                atu_mis['SpatialMissing'][cnt]=1            
            cnt=cnt+1

    # Import a complete list of openings from the OPENING layer
    nam_lyr='RSLT_OPENING_SVW'
    lyr=fiona.open(meta['Paths']['Results File'],layer=nam_lyr)
    L=len(lyr)
    op_mis={}
    op_mis['OPENING_ID']=np.zeros(L)
    cnt=0
    with fiona.open(meta['Paths']['Results File'],layer=nam_lyr) as source:
        for feat in source:
            op_mis['OPENING_ID'][cnt]=feat['properties']['OPENING_ID']
            cnt=cnt+1

    # Define indices to missing AT spatial openings and the same openings
    # in the OPENING layer
    # *** There are a small number of AT openings that are missing from the OP 
    # layer ***
    iMisAT=np.where(atu_mis['SpatialMissing']==1)[0]
    iMisOp=np.zeros(iMisAT.size)
    mis=0
    for i in range(iMisAT.size):
        ind0=np.where(op_mis['OPENING_ID']==atu_mis['OPENING_ID'][iMisAT[i]])[0]
        if ind0.size==0:
            mis=mis+1
            continue
        iMisOp[i]=ind0

    # Create a list that will contain the geometries for the missing AT rows
    at_geo_from_op=[None]*atu_mis['OPENING_ID'].size
    cnt=0
    with fiona.open(meta['Paths']['Results File'],layer=nam_lyr) as source:
        for feat in source:
            ind=np.where(iMisOp==cnt)[0]
            if ind.size!=0:
                at_geo_from_op[iMisAT[ind[0]]]=feat['geometry']
            cnt=cnt+1

    gu.opickle(meta['Paths']['Results'] + '\\atu_mis.pkl',atu_mis)
    gu.opickle(meta['Paths']['Results'] + '\\at_geo_from_op.pkl',at_geo_from_op)

'''============================================================================
ADD PLANTING INFO TO DMEH
============================================================================'''

# Create function to avoid duplication of bulky code
def AddPlantingWithNoData(d_nd):
    for i in range(1,6):
        d_nd['PL_SPECIES_CD' + str(i)]=np.append(d_nd['PL_SPECIES_CD' + str(i)],-999)
        d_nd['PL_SPECIES_PCT' + str(i)]=np.append(d_nd['PL_SPECIES_PCT' + str(i)],-999)
        d_nd['PL_SPECIES_GW' + str(i)]=np.append(d_nd['PL_SPECIES_GW' + str(i)],-999)
    return d_nd

'''============================================================================
PREPARE DISTURBANCE MANAGEMENT ENVENT HISTORY
============================================================================'''

def PrepDMEH(idx,meta,par,atu,pl,op,fcinv,vri,cut,fire,burnsev,pest):
    
    # Initiate disturbance-management event history
    dmeh=[None]*meta['N Stand']
    
    # Specify flag indicating whether subsetting occurs
    if np.isin('iKeep',list(meta.keys()))==True:
        # Tile project, only keeping a subset           
        flag_subset=1
    else:
        # Run all entries
        flag_subset=0
    
    for iStand0 in range(meta['N Stand']):
        
        if flag_subset==1:
            iStand=meta['iKeep'][iStand0]
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
        dmeh0={}
        dmeh0['Year']=-999*np.ones(1,dtype='float')
        dmeh0['Month']=-999*np.ones(1,dtype='int16')
        dmeh0['Day']=-999*np.ones(1,dtype='int16')
        dmeh0['ID_Type']=-999*np.ones(1,dtype='int16')
        dmeh0['MortalityFactor']=-999*np.ones(1,dtype='int16')
        dmeh0['GrowthFactor']=-999*np.ones(1,dtype='int16')
        dmeh0['FCI_Funded']=-999*np.ones(1,dtype='int16')
        dmeh0['FIA_PROJECT_ID']=-999*np.ones(1,dtype='int16')
        dmeh0['OPENING_ID']=-999*np.ones(1,dtype='int16')
        dmeh0['ACTUAL_TREATMENT_AREA']=-999*np.ones(1,dtype='int32')
        dmeh0['ACTUAL_PLANTED_NUMBER']=-999*np.ones(1,dtype='int32')
        dmeh0['PL_SPECIES_CD1']=-999*np.ones(1,dtype='int16')
        dmeh0['PL_SPECIES_CD2']=-999*np.ones(1,dtype='int16')
        dmeh0['PL_SPECIES_CD3']=-999*np.ones(1,dtype='int16')
        dmeh0['PL_SPECIES_CD4']=-999*np.ones(1,dtype='int16')
        dmeh0['PL_SPECIES_CD5']=-999*np.ones(1,dtype='int16')
        dmeh0['PL_SPECIES_PCT1']=-999*np.ones(1,dtype='int16')
        dmeh0['PL_SPECIES_PCT2']=-999*np.ones(1,dtype='int16')
        dmeh0['PL_SPECIES_PCT3']=-999*np.ones(1,dtype='int16')
        dmeh0['PL_SPECIES_PCT4']=-999*np.ones(1,dtype='int16')
        dmeh0['PL_SPECIES_PCT5']=-999*np.ones(1,dtype='int16')
        dmeh0['PL_SPECIES_GW1']=-999*np.ones(1,dtype='int16')
        dmeh0['PL_SPECIES_GW2']=-999*np.ones(1,dtype='int16')
        dmeh0['PL_SPECIES_GW3']=-999*np.ones(1,dtype='int16')
        dmeh0['PL_SPECIES_GW4']=-999*np.ones(1,dtype='int16')
        dmeh0['PL_SPECIES_GW5']=-999*np.ones(1,dtype='int16')
    
        #--------------------------------------------------------------------------
        # ATU layer
        #--------------------------------------------------------------------------
    
        if indS_atu!=None:
    
            for i in range(indS_atu['Index'].size):
            
                iA=indS_atu['Index'][i]
            
                #ind=np.where(atu['IdxToGrid']==iA)[0]
                
                # Convert RESULTS activity to fcgadgets event type
                dNames={}
                dNames['SILV_BASE_CODE']=cbu.lut_n2s(meta['LUT ATU']['SILV_BASE_CODE'],atu['SILV_BASE_CODE'][iA])
                dNames['SILV_TECHNIQUE_CODE']=cbu.lut_n2s(meta['LUT ATU']['SILV_TECHNIQUE_CODE'],atu['SILV_TECHNIQUE_CODE'][iA])
                dNames['SILV_METHOD_CODE']=cbu.lut_n2s(meta['LUT ATU']['SILV_METHOD_CODE'],atu['SILV_METHOD_CODE'][iA])        
                dNames['SILV_OBJECTIVE_CODE_1']=cbu.lut_n2s(meta['LUT ATU']['SILV_OBJECTIVE_CODE_1'],atu['SILV_OBJECTIVE_CODE_1'][iA])
                Name_Type=cbu.QueryResultsActivity(dNames)[0]
                try:
                    ID_Type=meta['LUT Dist'][Name_Type]
                except:
                    ID_Type=-999
                 
                # Skip surveys (added to avoid surveys being included in tile runs)
                if Name_Type=='Surveys':
                    continue
                    
                # Define FCI funded activities
                FCI_Funded=0
                if (atu['SILV_FUND_SOURCE_CODE'][iA]==meta['LUT ATU']['SILV_FUND_SOURCE_CODE']['FCE']) | \
                    (atu['SILV_FUND_SOURCE_CODE'][iA]==meta['LUT ATU']['SILV_FUND_SOURCE_CODE']['FCM']):
                    FCI_Funded=1    
                    # Criterion taken out - find another way - modifiy the FSC code instead prior to running this script
                    #(np.isin(cbu.lut_n2s(meta['LUT ATU']['FIA_PROJECT_ID'],atu['FIA_PROJECT_ID'][iA]),meta['uPP'])==True):
                    
                dmeh0['Year']=np.append(dmeh0['Year'],np.round(atu['Year'][iA]+atu['Month'][iA]/12,decimals=2))
                dmeh0['Month']=np.append(dmeh0['Month'],atu['Month'][iA])
                dmeh0['Day']=np.append(dmeh0['Day'],atu['Day'][iA])
                dmeh0['ID_Type']=np.append(dmeh0['ID_Type'],ID_Type)                
                dmeh0['MortalityFactor']=np.append(dmeh0['MortalityFactor'],np.array(0,dtype='int16'))
                dmeh0['GrowthFactor']=np.append(dmeh0['GrowthFactor'],np.array(0,dtype='int16'))
                dmeh0['FCI_Funded']=np.append(dmeh0['FCI_Funded'],FCI_Funded)
                dmeh0['OPENING_ID']=np.append(dmeh0['OPENING_ID'],atu['OPENING_ID'][iA])
                dmeh0['FIA_PROJECT_ID']=np.append(dmeh0['FIA_PROJECT_ID'],atu['FIA_PROJECT_ID'][iA])
                dmeh0['ACTUAL_TREATMENT_AREA']=np.append(dmeh0['ACTUAL_TREATMENT_AREA'],np.maximum(1,atu['ACTUAL_TREATMENT_AREA'][iA]))
                dmeh0['ACTUAL_PLANTED_NUMBER']=np.append(dmeh0['ACTUAL_PLANTED_NUMBER'],atu['ACTUAL_PLANTED_NUMBER'][iA])
                dmeh0=AddPlantingWithNoData(dmeh0)
    
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
                    ind=np.where(par['GW']['SEEDLOT_NUMBER']==pl['SEEDLOT_NUMBER'][indYear[j]])[0]
                    if ind.size!=0:
                        gw0[j]=par['GW']['GENETIC_WORTH_RTNG'][ind[0]]
    
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
        
                ind_dmeh=np.where( (np.floor(dmeh0['Year'])==uYear[iYear]) & 
                          (dmeh0['ID_Type']==meta['LUT Dist']['Planting']) )[0]
        
                # Populate q1 structure
                if cd1.size>0:
                    dmeh0['PL_SPECIES_CD1'][ind_dmeh]=cd1[0]
                    dmeh0['PL_SPECIES_PCT1'][ind_dmeh]=pct1[0]
                    dmeh0['PL_SPECIES_GW1'][ind_dmeh]=gw1[0]
                if cd1.size>1:
                    dmeh0['PL_SPECIES_CD2'][ind_dmeh]=cd1[1]
                    dmeh0['PL_SPECIES_PCT2'][ind_dmeh]=pct1[1]
                    dmeh0['PL_SPECIES_GW2'][ind_dmeh]=gw1[1]
                if cd1.size>2:
                    dmeh0['PL_SPECIES_CD3'][ind_dmeh]=cd1[2]
                    dmeh0['PL_SPECIES_PCT3'][ind_dmeh]=pct1[2]
                    dmeh0['PL_SPECIES_GW3'][ind_dmeh]=gw1[2]
                if cd1.size>3:
                    dmeh0['PL_SPECIES_CD4'][ind_dmeh]=cd1[3]
                    dmeh0['PL_SPECIES_PCT4'][ind_dmeh]=pct1[3]
                    dmeh0['PL_SPECIES_GW4'][ind_dmeh]=gw1[3]
                if cd1.size>4:
                    dmeh0['PL_SPECIES_CD5'][ind_dmeh]=cd1[4]
                    dmeh0['PL_SPECIES_PCT5'][ind_dmeh]=pct1[4]
                    dmeh0['PL_SPECIES_GW5'][ind_dmeh]=gw1[4]
    
        #--------------------------------------------------------------------------
        # Wildfire
        # Notes:
        # Burn severity does not come with month. In FCI contect, 2017 fires often 
        # occurred after planting. Assume month is August for now.
        #--------------------------------------------------------------------------
        
        if indS_burnsev!=None:        
        
            indS=indS_burnsev['Index']
            
            for i in range(indS.size):
            
                bsr=burnsev['BURN_SEVERITY_RATING'][indS[i]]
                if bsr==meta['LUT BS']['BURN_SEVERITY_RATING']['Low']: 
                    Severity=50
                elif bsr==meta['LUT BS']['BURN_SEVERITY_RATING']['Medium']:
                    Severity=90
                elif bsr==meta['LUT BS']['BURN_SEVERITY_RATING']['High']: 
                    Severity=100
                else:
                    Severity=5
            
                Month=8
                dmeh0['Year']=np.append(dmeh0['Year'],burnsev['FIRE_YEAR'][indS[i]]+Month/12)
                dmeh0['Month']=np.append(dmeh0['Month'],Month)
                dmeh0['Day']=np.append(dmeh0['Day'],-999)
                dmeh0['ID_Type']=np.append(dmeh0['ID_Type'],meta['LUT Dist']['Wildfire'])
                dmeh0['MortalityFactor']=np.append(dmeh0['MortalityFactor'],Severity)
                dmeh0['GrowthFactor']=np.append(dmeh0['GrowthFactor'],0)
                dmeh0['FCI_Funded']=np.append(dmeh0['FCI_Funded'],0)
                dmeh0['OPENING_ID']=np.append(dmeh0['OPENING_ID'],-999)
                dmeh0['FIA_PROJECT_ID']=np.append(dmeh0['FIA_PROJECT_ID'],-999)   
                dmeh0['ACTUAL_TREATMENT_AREA']=np.append(dmeh0['ACTUAL_TREATMENT_AREA'],-999)
                dmeh0['ACTUAL_PLANTED_NUMBER']=np.append(dmeh0['ACTUAL_PLANTED_NUMBER'],-999)
                dmeh0=AddPlantingWithNoData(dmeh0)
    
        if indS_fire!=None:
            
            indS=indS_fire['Index']
        
            for i in range(indS.size):
        
                # Only continue if nothing yet added through burnsev layer
                ind=np.where( (np.floor(dmeh0['Year'])==np.floor(fire['FIRE_YEAR'][indS[i]])) & 
                             (dmeh0['ID_Type']==meta['LUT Dist']['Wildfire']) )[0]
                if ind.size>0: 
                    continue
            
                dmeh0['Year']=np.append(dmeh0['Year'],fire['FIRE_YEAR'][indS[i]]+fire['Month'][indS[i]]/12)
                dmeh0['Month']=np.append(dmeh0['Month'],fire['Month'][indS[i]])
                dmeh0['Day']=np.append(dmeh0['Day'],fire['Day'][indS[i]])
                dmeh0['ID_Type']=np.append(dmeh0['ID_Type'],meta['LUT Dist']['Wildfire'])            
                dmeh0['MortalityFactor']=np.append(dmeh0['MortalityFactor'],50)
                dmeh0['GrowthFactor']=np.append(dmeh0['GrowthFactor'],0)
                dmeh0['FCI_Funded']=np.append(dmeh0['FCI_Funded'],0)
                dmeh0['OPENING_ID']=np.append(dmeh0['OPENING_ID'],-999)            
                dmeh0['FIA_PROJECT_ID']=np.append(dmeh0['FIA_PROJECT_ID'],-999)   
                dmeh0['ACTUAL_TREATMENT_AREA']=np.append(dmeh0['ACTUAL_TREATMENT_AREA'],-999)
                dmeh0['ACTUAL_PLANTED_NUMBER']=np.append(dmeh0['ACTUAL_PLANTED_NUMBER'],999)
                dmeh0=AddPlantingWithNoData(dmeh0)
        
        #--------------------------------------------------------------------------
        # Consolidated cutblocks
        #------------------------------------------------------------------------
    
        if indS_cut!=None:
        
            indS=indS_cut['Index']
        
            for i in range(indS.size):
                
                # Add the harvest
                yr_harv=cut['HARVEST_YEAR'][indS[i]]+1/12
                dmeh0['Year']=np.append(dmeh0['Year'],yr_harv)
                dmeh0['Month']=np.append(dmeh0['Month'],-999)
                dmeh0['Day']=np.append(dmeh0['Day'],-999)
                dmeh0['ID_Type']=np.append(dmeh0['ID_Type'],meta['LUT Dist']['Harvest'])
                dmeh0['MortalityFactor']=np.append(dmeh0['MortalityFactor'],100)
                dmeh0['GrowthFactor']=np.append(dmeh0['GrowthFactor'],0)
                dmeh0['FCI_Funded']=np.append(dmeh0['FCI_Funded'],0)
                dmeh0['OPENING_ID']=np.append(dmeh0['OPENING_ID'],-999)
                dmeh0['FIA_PROJECT_ID']=np.append(dmeh0['FIA_PROJECT_ID'],-999)   
                dmeh0['ACTUAL_TREATMENT_AREA']=np.append(dmeh0['ACTUAL_TREATMENT_AREA'],-999)
                dmeh0['ACTUAL_PLANTED_NUMBER']=np.append(dmeh0['ACTUAL_PLANTED_NUMBER'],-999)
                dmeh0=AddPlantingWithNoData(dmeh0)
                
                # Add the slashpile burn
                yr_burn=yr_harv+0.1
                dmeh0['Year']=np.append(dmeh0['Year'],yr_burn)
                dmeh0['Month']=np.append(dmeh0['Month'],-999)
                dmeh0['Day']=np.append(dmeh0['Day'],-999)
                dmeh0['ID_Type']=np.append(dmeh0['ID_Type'],meta['LUT Dist']['Slashpile Burn'])
                dmeh0['MortalityFactor']=np.append(dmeh0['MortalityFactor'],100)
                dmeh0['GrowthFactor']=np.append(dmeh0['GrowthFactor'],0)
                dmeh0['FCI_Funded']=np.append(dmeh0['FCI_Funded'],0)
                dmeh0['OPENING_ID']=np.append(dmeh0['OPENING_ID'],-999)
                dmeh0['FIA_PROJECT_ID']=np.append(dmeh0['FIA_PROJECT_ID'],-999)   
                dmeh0['ACTUAL_TREATMENT_AREA']=np.append(dmeh0['ACTUAL_TREATMENT_AREA'],-999)
                dmeh0['ACTUAL_PLANTED_NUMBER']=np.append(dmeh0['ACTUAL_PLANTED_NUMBER'],-999)
                dmeh0=AddPlantingWithNoData(dmeh0)

        #--------------------------------------------------------------------------
        # Harvesting from RESULTS opening layer
        #--------------------------------------------------------------------------
    
        if indS_op!=None: 
        
            iExisting=np.where( (dmeh0['ID_Type']==meta['LUT Dist']['Harvest']) )[0]
            YearExisting=np.floor(dmeh0['Year'][iExisting])
        
            indS=indS_op['Index']
        
            for i in range(indS.size):
        
                Year=np.floor(op['Year_Denu1_Comp'][indS[i]])
                Month=op['Month_Denu1_Comp'][indS[i]]
                cd=cbu.lut_n2s(meta['LUT OP']['DENUDATION_1_DISTURBANCE_CODE'],op['DENUDATION_1_DISTURBANCE_CODE'][indS[i]])
            
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
                    dmeh0['Year']=np.append(dmeh0['Year'],yr_harv)
                    dmeh0['Month']=np.append(dmeh0['Month'],Month)
                    dmeh0['Day']=np.append(dmeh0['Day'],-999)
                    dmeh0['ID_Type']=np.append(dmeh0['ID_Type'],meta['LUT Dist']['Harvest'])
                    dmeh0['MortalityFactor']=np.append(dmeh0['MortalityFactor'],100)
                    dmeh0['GrowthFactor']=np.append(dmeh0['GrowthFactor'],0)
                    dmeh0['FCI_Funded']=np.append(dmeh0['FCI_Funded'],0)
                    dmeh0['OPENING_ID']=np.append(dmeh0['OPENING_ID'],-999)
                    dmeh0['FIA_PROJECT_ID']=np.append(dmeh0['FIA_PROJECT_ID'],-999)   
                    dmeh0['ACTUAL_TREATMENT_AREA']=np.append(dmeh0['ACTUAL_TREATMENT_AREA'],-999)
                    dmeh0['ACTUAL_PLANTED_NUMBER']=np.append(dmeh0['ACTUAL_PLANTED_NUMBER'],-999)
                    dmeh0=AddPlantingWithNoData(dmeh0)
                    
                    # Add the slashpile burn
                    yr_burn=yr_harv+0.1
                    dmeh0['Year']=np.append(dmeh0['Year'],yr_burn)
                    dmeh0['Month']=np.append(dmeh0['Month'],-999)
                    dmeh0['Day']=np.append(dmeh0['Day'],-999)
                    dmeh0['ID_Type']=np.append(dmeh0['ID_Type'],meta['LUT Dist']['Slashpile Burn'])
                    dmeh0['MortalityFactor']=np.append(dmeh0['MortalityFactor'],100)
                    dmeh0['GrowthFactor']=np.append(dmeh0['GrowthFactor'],0)
                    dmeh0['FCI_Funded']=np.append(dmeh0['FCI_Funded'],0)
                    dmeh0['OPENING_ID']=np.append(dmeh0['OPENING_ID'],-999)
                    dmeh0['FIA_PROJECT_ID']=np.append(dmeh0['FIA_PROJECT_ID'],-999)   
                    dmeh0['ACTUAL_TREATMENT_AREA']=np.append(dmeh0['ACTUAL_TREATMENT_AREA'],-999)
                    dmeh0['ACTUAL_PLANTED_NUMBER']=np.append(dmeh0['ACTUAL_PLANTED_NUMBER'],-999)
                    dmeh0=AddPlantingWithNoData(dmeh0)
        
                Year=np.floor(op['Year_Denu2_Comp'][indS[i]])
                Month=op['Month_Denu2_Comp'][indS[i]]
                cd=cbu.lut_n2s(meta['LUT OP']['DENUDATION_2_DISTURBANCE_CODE'],op['DENUDATION_2_DISTURBANCE_CODE'][indS[i]])
            
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
                    dmeh0['Year']=np.append(dmeh0['Year'],yr_harv)
                    dmeh0['Month']=np.append(dmeh0['Month'],Month)
                    dmeh0['Day']=np.append(dmeh0['Day'],-999)
                    dmeh0['ID_Type']=np.append(dmeh0['ID_Type'],meta['LUT Dist']['Harvest'])
                    dmeh0['MortalityFactor']=np.append(dmeh0['MortalityFactor'],100)
                    dmeh0['GrowthFactor']=np.append(dmeh0['GrowthFactor'],0)
                    dmeh0['FCI_Funded']=np.append(dmeh0['FCI_Funded'],0)
                    dmeh0['OPENING_ID']=np.append(dmeh0['OPENING_ID'],-999)
                    dmeh0['FIA_PROJECT_ID']=np.append(dmeh0['FIA_PROJECT_ID'],-999)   
                    dmeh0['ACTUAL_TREATMENT_AREA']=np.append(dmeh0['ACTUAL_TREATMENT_AREA'],-999)
                    dmeh0['ACTUAL_PLANTED_NUMBER']=np.append(dmeh0['ACTUAL_PLANTED_NUMBER'],-999)
                    dmeh0=AddPlantingWithNoData(dmeh0)
                    
                    # Add the slashpile burn
                    yr_burn=yr_harv+0.1
                    dmeh0['Year']=np.append(dmeh0['Year'],yr_burn)
                    dmeh0['Month']=np.append(dmeh0['Month'],-999)
                    dmeh0['Day']=np.append(dmeh0['Day'],-999)
                    dmeh0['ID_Type']=np.append(dmeh0['ID_Type'],meta['LUT Dist']['Slashpile Burn'])
                    dmeh0['MortalityFactor']=np.append(dmeh0['MortalityFactor'],100)
                    dmeh0['GrowthFactor']=np.append(dmeh0['GrowthFactor'],0)
                    dmeh0['FCI_Funded']=np.append(dmeh0['FCI_Funded'],0)
                    dmeh0['OPENING_ID']=np.append(dmeh0['OPENING_ID'],-999)
                    dmeh0['FIA_PROJECT_ID']=np.append(dmeh0['FIA_PROJECT_ID'],-999)   
                    dmeh0['ACTUAL_TREATMENT_AREA']=np.append(dmeh0['ACTUAL_TREATMENT_AREA'],-999)
                    dmeh0['ACTUAL_PLANTED_NUMBER']=np.append(dmeh0['ACTUAL_PLANTED_NUMBER'],-999)
                    dmeh0=AddPlantingWithNoData(dmeh0)
    
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
            
                dmeh0['Year']=np.append(dmeh0['Year'],pest['CAPTURE_YEAR'][iYr])
                dmeh0['Month']=np.append(dmeh0['Month'],-999)
                dmeh0['Day']=np.append(dmeh0['Day'],-999)
                dmeh0['FCI_Funded']=np.append(dmeh0['FCI_Funded'],0)
                dmeh0['OPENING_ID']=np.append(dmeh0['OPENING_ID'],-999)
                dmeh0['FIA_PROJECT_ID']=np.append(dmeh0['FIA_PROJECT_ID'],-999)   
                dmeh0['ACTUAL_TREATMENT_AREA']=np.append(dmeh0['ACTUAL_TREATMENT_AREA'],-999)
                dmeh0['ACTUAL_PLANTED_NUMBER']=np.append(dmeh0['ACTUAL_PLANTED_NUMBER'],-999)
                dmeh0=AddPlantingWithNoData(dmeh0)
            
                # Shorten for repeated use
                psp=pest['PEST_SPECIES_CODE'][iYr]
                sev=pest['PEST_SEVERITY_CODE'][iYr]
                sev_s=cbu.lut_n2s(meta['LUT Pest']['PEST_SEVERITY_CODE'],sev)[0]
            
                # Mountain pine beetle - Adjust severity based on fractin of pine later using: AdjustSpeciesSpecificMortality         
                if (psp==meta['LUT Pest']['PEST_SPECIES_CODE']['IBM']):
                    dmeh0['ID_Type']=np.append(dmeh0['ID_Type'],meta['LUT Dist']['IBM'])
                    ind=np.where( (par['DistBySC']['Name']=='IBM') & (par['DistBySC']['SeverityCD']==sev_s) )[0]
                    dmeh0['MortalityFactor']=np.append(dmeh0['MortalityFactor'],par['DistBySC']['MortalityFactor'][ind])
                    dmeh0['GrowthFactor']=np.append(dmeh0['GrowthFactor'],par['DistBySC']['GrowthFactor'][ind])
                
                # Western Balsam Bark Beetle
                # (https://www2.gov.bc.ca/gov/content/industry/forestry/managing-our-forest-resources/forest-health/forest-pests/bark-beetles/western-balsam-bark-beetle)
                elif (psp==meta['LUT Pest']['PEST_SPECIES_CODE']['IBB']):
                    dmeh0['ID_Type']=np.append(dmeh0['ID_Type'],meta['LUT Dist']['IBB'])
                    ind=np.where( (par['DistBySC']['Name']=='IBB') & (par['DistBySC']['SeverityCD']==sev_s) )[0]
                    dmeh0['MortalityFactor']=np.append(dmeh0['MortalityFactor'],par['DistBySC']['MortalityFactor'][ind])
                    dmeh0['GrowthFactor']=np.append(dmeh0['GrowthFactor'],par['DistBySC']['GrowthFactor'][ind])
            
                # Douglas-fir beetle
                # (https://www2.gov.bc.ca/gov/content/industry/forestry/managing-our-forest-resources/forest-health/forest-pests/bark-beetles/western-balsam-bark-beetle)
                elif (psp==meta['LUT Pest']['PEST_SPECIES_CODE']['IBD']):
                    dmeh0['ID_Type']=np.append(dmeh0['ID_Type'],meta['LUT Dist']['IBD'])
                    ind=np.where( (par['DistBySC']['Name']=='IBD') & (par['DistBySC']['SeverityCD']==sev_s) )[0]
                    dmeh0['MortalityFactor']=np.append(dmeh0['MortalityFactor'],par['DistBySC']['MortalityFactor'][ind])
                    dmeh0['GrowthFactor']=np.append(dmeh0['GrowthFactor'],par['DistBySC']['GrowthFactor'][ind])
                
                # Spruce beetle
                # (https://www2.gov.bc.ca/gov/content/industry/forestry/managing-our-forest-resources/forest-health/forest-pests/bark-beetles/western-balsam-bark-beetle)
                elif (psp==meta['LUT Pest']['PEST_SPECIES_CODE']['IBS']):
                    dmeh0['ID_Type']=np.append(dmeh0['ID_Type'],meta['LUT Dist']['IBS'])
                    ind=np.where( (par['DistBySC']['Name']=='IBS') & (par['DistBySC']['SeverityCD']==sev_s) )[0]
                    dmeh0['MortalityFactor']=np.append(dmeh0['MortalityFactor'],par['DistBySC']['MortalityFactor'][ind])
                    dmeh0['GrowthFactor']=np.append(dmeh0['GrowthFactor'],par['DistBySC']['GrowthFactor'][ind])
                
                # Western spruce budworm 
                # Populate severity with the ID for severity - it will be revised below
                # to model mortality that occurs from repeated infestation
                elif (psp==meta['LUT Pest']['PEST_SPECIES_CODE']['IDW']):
                    dmeh0['ID_Type']=np.append(dmeh0['ID_Type'],meta['LUT Dist']['IDW'])
                    ind=np.where( (par['DistBySC']['Name']=='IDW') & (par['DistBySC']['SeverityCD']==sev_s) )[0]
                    dmeh0['MortalityFactor']=np.append(dmeh0['MortalityFactor'],sev)
                    dmeh0['GrowthFactor']=np.append(dmeh0['GrowthFactor'],par['DistBySC']['GrowthFactor'][ind])
                
                # Other
                else:
                    dmeh0['ID_Type']=np.append(dmeh0['ID_Type'],meta['LUT Dist']['Beetles'])
                    dmeh0['MortalityFactor']=np.append(dmeh0['MortalityFactor'],5)
                    dmeh0['GrowthFactor']=np.append(dmeh0['GrowthFactor'],0)
    
        #--------------------------------------------------------------------------
        # Append to management history list
        #--------------------------------------------------------------------------
    
        dmeh[iStand0]=dmeh0
    
    return dmeh

'''============================================================================
LOAD Look-up-tables
============================================================================'''

def Load_LUTs(meta):
    
    meta['LUT ATU']=gu.ipickle(meta['Paths']['Results'] + '\\LUTs_RSLT_ACTIVITY_TREATMENT_SVW.pkl')
    meta['LUT OP']=gu.ipickle(meta['Paths']['Results'] + '\\LUTs_RSLT_OPENING_SVW.pkl')
    meta['LUT PL']=gu.ipickle(meta['Paths']['Results'] + '\\LUTs_RSLT_PLANTING_SVW.pkl')
    meta['LUT FC_I']=gu.ipickle(meta['Paths']['Results'] + '\\LUTs_RSLT_FOREST_COVER_INV_SVW.pkl')
    meta['LUT FC_S']=gu.ipickle(meta['Paths']['Results'] + '\\LUTs_RSLT_FOREST_COVER_SILV_SVW.pkl')
    meta['LUT VRI']=gu.ipickle(meta['Paths']['VRI'] + '\\LUTs_VEG_COMP_LYR_R1_POLY.pkl')
    meta['LUT BS']=gu.ipickle(meta['Paths']['Disturbances'] + '\\LUTs_VEG_BURN_SEVERITY_SP.pkl')
    meta['LUT Pest']=gu.ipickle(meta['Paths']['Disturbances'] + '\\LUTs_PEST_INFESTATION_POLY.pkl')
    try:
        meta['LUT LU NL']=gu.ipickle(meta['Paths']['LandUse'] + '\\LUTs_RMP_PLAN_NON_LEGAL_POLY_SVW.pkl')
        meta['LUT LU L']=gu.ipickle(meta['Paths']['LandUse'] + '\\LUTs_RMP_PLAN_LEGAL_POLY_SVW.pkl')
        meta['LUT PARK']=gu.ipickle(meta['Paths']['LandUse'] + '\\LUTs_TA_PARK_ECORES_PA_SVW.pkl')
        meta['LUT OGMA']=gu.ipickle(meta['Paths']['LandUse'] + '\\LUTs_RMP_OGMA_LEGAL_CURRENT_SVW.pkl')
        meta['LUT UWR']=gu.ipickle(meta['Paths']['LandUse'] + '\\LUTs_WCP_UNGULATE_WINTER_RANGE_SP.pkl')   
    except:
        pass
    meta['LUT TIPSY']={}
    meta['LUT TIPSY']['FIZ']={'C':np.array(1,dtype=int),'I':np.array(2,dtype=int)}
    meta['LUT TIPSY']['regeneration_method']={'C':np.array(1,dtype=int),'N':np.array(2,dtype=int),'P':np.array(3,dtype=int)}
    
    return meta

'''============================================================================
LOAD PARAMTERS
============================================================================'''

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
    par['Ph_Simp']=gu.ReadExcel(r'G:\My Drive\Data\Harvest\Historical BC Harvesting\\HarvestHistoricalProbabilitySimple.xlsx')        
            
    return par

'''============================================================================
EXCLUDE DUPLICATE EVENTS
============================================================================'''

def Exclude_Duplicate_Events(meta,dmeh):
    for iStand in range(meta['N Stand Full']):
        if dmeh[iStand]==None:
            continue
        for key in meta['LUT Dist'].keys():
            ind=np.where(dmeh[iStand]['ID_Type']==meta['LUT Dist'][key])[0]
            if ind.size==0:
                continue
            uYear=np.unique(np.floor(dmeh[iStand]['Year'][ind]))
            for iYear in range(uYear.size):
                ind=np.where( (dmeh[iStand]['ID_Type']==meta['LUT Dist'][key]) & (np.floor(dmeh[iStand]['Year'])==uYear[iYear]) )[0]
                dmeh[iStand]['ID_Type'][ind[1:]]=-999
    return dmeh

'''============================================================================
EXCLUDE UNIDENTIFIED EVENTS
============================================================================'''

def Exclude_Unidentified_Events(meta,dmeh):

    for iStand in range(meta['N Stand Full']):
        if dmeh[iStand]==None:
            continue
        ind=np.where(dmeh[iStand]['ID_Type']!=-999)[0]
        for key in dmeh[iStand]:
            dmeh[iStand][key]=dmeh[iStand][key][ind]
    return dmeh

'''============================================================================
REMOVE SLASHPILE BURNS IN SELECT BGC ZONES
============================================================================'''

def Remove_SlashpileBurns_From_Select_Zones(meta,dmeh,ba):

    for iStand in range(meta['N Stand Full']):
        
        if dmeh[iStand]==None:
            continue
        
        if (ba['BEC_ZONE_CODE'][iStand]==meta['LUT VRI']['BEC_ZONE_CODE']['CWH']) | (ba['BEC_ZONE_CODE'][iStand]==meta['LUT VRI']['BEC_ZONE_CODE']['ICH']):
            ind=np.where(dmeh[iStand]['ID_Type']!=meta['LUT Dist']['Slashpile Burn'])[0]
            if ind.size>0:
                for key in dmeh[iStand]:
                    if key=='ScnAffected':
                        continue
                    dmeh[iStand][key]=dmeh[iStand][key][ind]
    
    return dmeh


'''============================================================================
ENSURE EVERY STAND HAS A MODERN DISTURBANCE
============================================================================'''

def Ensure_Every_Stand_Has_Modern_Disturbance(meta,dmeh,name_dist,severity):
    for iStand in range(meta['N Stand Full']):
        if dmeh[iStand]==None:
            continue
        if dmeh[iStand]['Year'].size==0:
            #print(iStand)
            #break
            r=np.random.randint(1700,2000)
            dmeh[iStand]['Year']=np.append(dmeh[iStand]['Year'],r)
            dmeh[iStand]['ID_Type']=np.append(dmeh[iStand]['ID_Type'],meta['LUT Dist'][name_dist])
            dmeh[iStand]['MortalityFactor']=np.append(dmeh[iStand]['MortalityFactor'],np.array(severity,dtype='int16'))
            dmeh[iStand]['GrowthFactor']=np.append(dmeh[iStand]['GrowthFactor'],np.array(0,dtype='int16'))
            for v in StringsToFill:
                dmeh[iStand][v]=np.append(dmeh[iStand][v],-999)
    return dmeh

'''============================================================================
ENSURE DISTURBANCE PRECEDES AERIAL FERTILIZATION
So that age at fert is specified.
============================================================================'''

def Ensure_Fert_Preceded_By_Disturbance(meta,dmeh,th_sev_last_dist):

    # Assumptions
    AgeAtFert=35
    #ID_Type_Prev=meta['LUT Dist']['Harvest']
    #sev_Prev=100

    Ntot=np.array([])
    Nfix=np.array([])
    for iStand in range(meta['N Stand Full']):
        
        if dmeh[iStand]==None:
            continue
    
        iA=np.where( (dmeh[iStand]['ID_Type']==meta['LUT Dist']['Fertilization Aerial']) )[0]
        if iA.size==0: 
            continue
    
        # QA: Look at what events occur before fert
        #print(cbu.lut_n2s(meta['LUT Dist'],dmeh[iStand]['ID_Type'][iA[0]-1]))
        #print(dmeh[iStand]['MortalityFactor'][iA[0]-1])    

        Ntot=np.append(Ntot,1)
        
        if (dmeh[iStand]['Year'][iA[0]]==np.min(dmeh[iStand]['Year'])) | (np.max(dmeh[iStand]['MortalityFactor'][0:iA[0]+1])<th_sev_last_dist):        
        
            Nfix=np.append(Nfix,1)    
            
            # Add harvest
            Year=dmeh[iStand]['Year'][iA[0]]-AgeAtFert
            dmeh[iStand]['Year']=np.append(dmeh[iStand]['Year'],Year)
            dmeh[iStand]['ID_Type']=np.append(dmeh[iStand]['ID_Type'],meta['LUT Dist']['Harvest'])
            dmeh[iStand]['MortalityFactor']=np.append(dmeh[iStand]['MortalityFactor'],100)
            dmeh[iStand]['GrowthFactor']=np.append(dmeh[iStand]['GrowthFactor'],np.array(0,dtype='int16'))
            for v in StringsToFill:
                dmeh[iStand][v]=np.append(dmeh[iStand][v],-999)  
            
            # Add slashpile burn
            dmeh[iStand]['Year']=np.append(dmeh[iStand]['Year'],Year+0.1)
            dmeh[iStand]['ID_Type']=np.append(dmeh[iStand]['ID_Type'],meta['LUT Dist']['Slashpile Burn'])
            dmeh[iStand]['MortalityFactor']=np.append(dmeh[iStand]['MortalityFactor'],100)
            dmeh[iStand]['GrowthFactor']=np.append(dmeh[iStand]['GrowthFactor'],np.array(0,dtype='int16'))
            for v in StringsToFill:
                dmeh[iStand][v]=np.append(dmeh[iStand][v],-999)  

    # What fraction of fertilization events had no record of a stand-replacing disturbance
    print(str(np.sum(Nfix)/np.sum(Ntot)) + ' of study area required generating a pre-fert disturbance.')

    #Age=vri['PROJ_AGE_1'][iStand]+dmeh[iStand]['Year'][0]-vri['Year'][iStand]

    # Put events in order of calendar date
    for iStand in range(meta['N Stand Full']):
        d=dmeh[iStand].copy()
        ord=np.argsort(d['Year'])
        for key in d.keys():
            d[key]=d[key][ord]
        dmeh[iStand]=d.copy()
    
    return meta,dmeh

'''============================================================================
IDW FIX SEVERITY
# The DMEH was populated with the numeric severity ID. Mortality only occurs 
# following repeated outrbreak years. 
============================================================================'''

def IDW_Fix_Severity(meta,dmeh,par):

    for iStand in range(meta['N Stand Full']):
        
        if dmeh[iStand]==None:
            continue
        
        # Save a frozen version of the initial severity
        Severity_Frozen=dmeh[iStand]['MortalityFactor'].copy()
    
        for iA in range(dmeh[iStand]['Year'].size):
            
            if dmeh[iStand]['ID_Type'][iA]=='IDW':
            
                # By default, apply mortality rate in the year of defoliation
                sev1=dmeh[iStand]['MortalityFactor'][iA]
                sev_s1=cbu.lut_n2s(meta['LUT Pest']['PEST_SEVERITY_CODE'],sev1)[0]
                ind=np.where( (par['DistBySC']['Name']=='IDW') & (par['DistBySC']['SeverityCD']==sev_s1) )[0]
                Mortality1=par['DistBySV']['MortalityFactor'][ind]

                if iA>0:
                    if dmeh[iStand]['ID_Type'][iA-1]=='IDW':
                    
                        # If it is back to back infestation, adjust mortality accordingly
                        sev0=Severity_Frozen[iA-1]
                        sev_s0=cbu.lut_n2s(meta['LUT Pest']['PEST_SEVERITY_CODE'],sev0)[0]
                        if (sev_s0=='M') & (sev_s1=='M'):
                            ind=np.where( (par['DistBySC']['Name']=='IDW') & (par['DistBySC']['SeverityCD']=='MM') )[0]
                            Mortality1=par['DistBySV']['MortalityFactor'][ind]
                        elif (sev_s0=='M') & (sev_s1=='S'):
                            ind=np.where( (par['DistBySC']['Name']=='IDW') & (par['DistBySC']['SeverityCD']=='MS') )[0]
                            Mortality1=par['DistBySV']['MortalityFactor'][ind]
                        elif (sev_s0=='M') & (sev_s1=='V'):
                            ind=np.where( (par['DistBySC']['Name']=='IDW') & (par['DistBySC']['SeverityCD']=='MV') )[0]
                            Mortality1=par['DistBySV']['MortalityFactor'][ind]
                        elif (sev_s0=='S') & (sev_s1=='S'):
                            ind=np.where( (par['DistBySC']['Name']=='IDW') & (par['DistBySC']['SeverityCD']=='SS') )[0]
                            Mortality1=par['DistBySV']['MortalityFactor'][ind]
                        elif (sev_s0=='S') & (sev_s1=='V'):
                            ind=np.where( (par['DistBySC']['Name']=='IDW') & (par['DistBySC']['SeverityCD']=='SV') )[0]
                            Mortality1=par['DistBySV']['MortalityFactor'][ind]
                        elif (sev_s0=='V') & (sev_s1=='V'):
                            ind=np.where( (par['DistBySC']['Name']=='IDW') & (par['DistBySC']['SeverityCD']=='VV') )[0]
                            Mortality1=par['DistBySV']['MortalityFactor'][ind]
                        else:
                            pass                    
                dmeh[iStand]['MortalityFactor'][iA]=np.array(Mortality1,dtype='int16')
    return dmeh

'''============================================================================
ADD OLDEST KNOWN DISTURBANCE FROM VRI
This really doesn't work well.
============================================================================'''

def Add_Oldest_Disturbance_From_VRI(meta,dmeh,idx,vri):
    
    for iStand in range(meta['N Stand Full']):
        if idx['vri'][iStand]==None:
            continue
        ind=idx['vri'][iStand]['Index'][0]
        DOE=vri['Year'][ind]-vri['PROJ_AGE_1'][ind]
        flg=0
        if dmeh[iStand]['Year'].size==0:
            flg==1
        else:
            if (DOE>0) & (DOE<np.min(dmeh[iStand]['Year'])):
                flg==1
        if flg==1:
            dmeh[iStand]['Year']=np.append(dmeh[iStand]['Year'],DOE)
            dmeh[iStand]['ID_Type']=np.append(dmeh[iStand]['ID_Type'],meta['LUT Dist']['Wildfire'])
            dmeh[iStand]['MortalityFactor']=np.append(dmeh[iStand]['MortalityFactor'],np.array(100,dtype='int16'))
            dmeh[iStand]['GrowthFactor']=np.append(dmeh[iStand]['GrowthFactor'],np.array(0,dtype='int16'))
            for v in StringsToFill:
                dmeh[iStand][v]=np.append(dmeh[iStand][v],-999)

    # Put events in order of calendar date
    for iStand in range(meta['N Stand Full']):
        d=dmeh[iStand].copy()
        ord=np.argsort(d['Year'])
        for key in d.keys():
            d[key]=d[key][ord]
        dmeh[iStand]=d.copy()
        
    return dmeh

'''============================================================================
CLEAN SPECIES COMPOSITION
============================================================================'''

def Clean_Species_Composition(meta,dmeh,vri,fcinv):
    
    # List of code pairs (code to change, code to change to)
    ListS=[('SXW','SX')]
    
    # Fix DMEH
    for iStand in range(meta['N Stand Full']):   
        for iSpc in range(len(ListS)):       
            # Disturbance/management inventory
            n0=meta['LUT PL']['SILV_TREE_SPECIES_CODE'][ListS[iSpc][0]]
            n1=meta['LUT PL']['SILV_TREE_SPECIES_CODE'][ListS[iSpc][1]]
            for k in range(dmeh[iStand]['PL_SPECIES_CD1'].size):
                if dmeh[iStand]['PL_SPECIES_CD1'][k]==n0: dmeh[iStand]['PL_SPECIES_CD1'][k]=n1
                if dmeh[iStand]['PL_SPECIES_CD2'][k]==n0: dmeh[iStand]['PL_SPECIES_CD2'][k]=n1
                if dmeh[iStand]['PL_SPECIES_CD3'][k]==n0: dmeh[iStand]['PL_SPECIES_CD3'][k]=n1
                if dmeh[iStand]['PL_SPECIES_CD4'][k]==n0: dmeh[iStand]['PL_SPECIES_CD4'][k]=n1
                if dmeh[iStand]['PL_SPECIES_CD5'][k]==n0: dmeh[iStand]['PL_SPECIES_CD5'][k]=n1
 
    # Fix VRI
    for iStand in range(vri['SPECIES_CD_1'].size):       
        for iSpc in range(len(ListS)):    
            for k in range(6):
                n0=meta['LUT VRI']['SPECIES_CD_' + str(k+1)][ListS[iSpc][0]];
                n1=meta['LUT VRI']['SPECIES_CD_' + str(k+1)][ListS[iSpc][1]]        
                if vri['SPECIES_CD_' + str(k+1)][iStand]==n0: 
                    vri['SPECIES_CD_' + str(k+1)][iStand]=n1
        
    # Forest cover
    for iStand in range(fcinv['I_SPECIES_CODE_1'].size):       
        for iSpc in range(len(ListS)):         
            for k in range(5):
                n0=meta['LUT FC_I']['I_SPECIES_CODE_' + str(k+1)][ListS[iSpc][0]]
                n1=meta['LUT FC_I']['I_SPECIES_CODE_' + str(k+1)][ListS[iSpc][1]]        
                if fcinv['I_SPECIES_CODE_' + str(k+1)][iStand]==n0: 
                    fcinv['I_SPECIES_CODE_' + str(k+1)][iStand]=n1
                
                #n0=meta['LUT FC_S']['S_SPECIES_CODE_' + str(k+1)][ListS[j][0]]
                #n1=meta['LUT FC_S']['S_SPECIES_CODE_' + str(k+1)][ListS[j][1]]        
                #if fcSd['S_SPECIES_CODE_' + str(k+1)][i]==n0: 
                #    fcSd['S_SPECIES_CODE_' + str(k+1)][i]=n1
    
    return meta,dmeh,vri,fcinv

'''============================================================================
CREATE BEST AVAILABLE INVENTORY
============================================================================'''

def CreateBestAvailableInventory(meta,vri,fcinv,flag_projects,idx,sxy):

    #--------------------------------------------------------------------------
    # Initialize best-available (gap-filled) inventory
    #--------------------------------------------------------------------------
    
    ba={} 
    ba['FIZ']=meta['LUT TIPSY']['FIZ']['I']*np.ones(meta['N Stand'])
    ba['BEC_ZONE_CODE']=meta['LUT VRI']['BEC_ZONE_CODE']['SBS']*np.ones(meta['N Stand'])
    ba['Spc_CD1']=-999*np.ones(meta['N Stand'])
    ba['Spc_CD2']=-999*np.ones(meta['N Stand'])
    ba['Spc_CD3']=-999*np.ones(meta['N Stand'])
    ba['Spc_CD4']=-999*np.ones(meta['N Stand'])
    ba['Spc_CD5']=-999*np.ones(meta['N Stand'])
    ba['Spc_Pct1']=-999*np.ones(meta['N Stand'])
    ba['Spc_Pct2']=-999*np.ones(meta['N Stand'])
    ba['Spc_Pct3']=-999*np.ones(meta['N Stand'])
    ba['Spc_Pct4']=-999*np.ones(meta['N Stand'])
    ba['Spc_Pct5']=-999*np.ones(meta['N Stand'])
    ba['SI']=-999*np.ones(meta['N Stand'])

    # Also keep track of the data source percentages
    basp={}
    basp['BEC_ZONE_CODE']={}
    basp['FIZ']={}
    basp['Spc_CD1']={}
    basp['SI']={}

    #--------------------------------------------------------------------------
    # Specify flag indicating whether subsetting occurs
    #--------------------------------------------------------------------------
    
    if np.isin('iKeep',list(meta.keys()))==True:
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
    for iStand0 in range(meta['N Stand']):
    
        if flag_subset==1:
            iStand=meta['iKeep'][iStand0]
        else:
            iStand=iStand0
        
        if idx['vri'][iStand]==None:
            continue
    
        ind0=idx['vri'][iStand]['Index'][0]
        N_tot=N_tot+ind0.size
        
        ba['BEC_ZONE_CODE'][iStand0]=vri['BEC_ZONE_CODE'][ind0]        
        
        if (vri['BEC_ZONE_CODE'][ind0]==meta['LUT VRI']['BEC_ZONE_CODE']['CWH']) | \
            (vri['BEC_ZONE_CODE'][ind0]==meta['LUT VRI']['BEC_ZONE_CODE']['CDF']) | \
            (vri['BEC_ZONE_CODE'][ind0]==meta['LUT VRI']['BEC_ZONE_CODE']['MH']):
            ba['FIZ'][iStand0]=meta['LUT TIPSY']['FIZ']['C']
        else:
            ba['FIZ'][iStand0]=meta['LUT TIPSY']['FIZ']['I']
    
    basp['BEC_ZONE_CODE']['From VRI']=N_tot/ba['SI'].size*100
    basp['FIZ']['From VRI']=N_tot/ba['SI'].size*100
    
    # Fill with global assumption
    ind=np.where(ba['BEC_ZONE_CODE']<=0)[0]
    if ind.size>0:
        ba['BEC_ZONE_CODE'][ind]=meta['LUT VRI']['BEC_ZONE_CODE']['SBS']
        ba['FIZ'][ind]=meta['LUT TIPSY']['FIZ']['I']
    
        basp['BEC_ZONE_CODE']['From global gap filling assumption']=ind.size/ba['SI'].size*100
        basp['FIZ']['From global gap filling assumption']=ind.size/ba['SI'].size*100
    
    #--------------------------------------------------------------------------
    # Best-available species composition (not planted)
    #--------------------------------------------------------------------------

    # Fill with forest cover
    N_tot=0
    for iStand0 in range(meta['N Stand']):
        
        if flag_subset==1:
            iStand=meta['iKeep'][iStand0]
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
    for iStand0 in range(meta['N Stand']):  
        
        if flag_subset==1:
            iStand=meta['iKeep'][iStand0]
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
    for iStand0 in range(meta['N Stand']):
        
        if flag_subset==1:
            iStand=meta['iKeep'][iStand0]
        else:
            iStand=iStand0
        
        if ba['Spc_CD1'][iStand0]>0:
            # Already populated with FC inventory layer
            continue
        
        if ba['FIZ'][iStand0]==meta['LUT TIPSY']['FIZ']['C']:
            ba['Spc_CD1'][iStand0]=meta['LUT VRI']['SPECIES_CD_1']['FD']
            ba['Spc_Pct1'][iStand0]=100
        else:
            ba['Spc_CD1'][iStand0]=meta['LUT VRI']['SPECIES_CD_1']['PL']
            ba['Spc_Pct1'][iStand0]=100 
        
        N_tot=N_tot+1
    
    basp['Spc_CD1']['From regional assumptions']=N_tot/ba['SI'].size*100

    #--------------------------------------------------------------------------
    # Import site productivity layer    
    #--------------------------------------------------------------------------
    
    spl={}
    spl['SI_SPL']=-999*np.ones(meta['N Stand'])
    
    z=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\SiteProductivityLayer\Site_Prod_Fd.tif')
    x=z['X'][0,:]
    y=z['Y'][:,0]
    Site_Prod_Fd=z['Data']
    del z
    garc.collect()
    
    z=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\SiteProductivityLayer\Site_Prod_Pl.tif')
    Site_Prod_Pl=z['Data']
    del z
    garc.collect()
    
    z=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\SiteProductivityLayer\Site_Prod_Sx.tif')
    Site_Prod_Sx=z['Data']
    del z
    garc.collect()
    
    # Populate dictionary with nearest estimate
    for iStand0 in range(meta['N Stand']):
        
        if flag_subset==1:
            iStand=meta['iKeep'][iStand0]
        else:
            iStand=iStand0
        
        adx=np.abs(sxy['x'][iStand]-x)
        ady=np.abs(sxy['y'][iStand]-y)
        ix=np.where(adx==np.min(adx))[0]
        iy=np.where(ady==np.min(ady))[0]
        ind=np.ix_(iy,ix)
        if (ba['Spc_CD1'][iStand0]==meta['LUT VRI']['SPECIES_CD_1']['FD']) | \
            (ba['Spc_CD1'][iStand0]==meta['LUT VRI']['SPECIES_CD_1']['FDI']) | \
            (ba['Spc_CD1'][iStand0]==meta['LUT VRI']['SPECIES_CD_1']['FDC']):
            if Site_Prod_Fd[ind][0][0]>0: 
                spl['SI_SPL'][iStand0]=Site_Prod_Fd[ind][0][0]
        elif (ba['Spc_CD1'][iStand0]==meta['LUT VRI']['SPECIES_CD_1']['PL']) | \
            (ba['Spc_CD1'][iStand0]==meta['LUT VRI']['SPECIES_CD_1']['PLI']) | \
            (ba['Spc_CD1'][iStand0]==meta['LUT VRI']['SPECIES_CD_1']['PLC']):
            if Site_Prod_Pl[ind][0][0]>0: 
                spl['SI_SPL'][iStand0]=np.maximum(0,Site_Prod_Pl[ind][0][0])
        elif (ba['Spc_CD1'][iStand0]==meta['LUT VRI']['SPECIES_CD_1']['SX']):
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
    for iStand0 in range(meta['N Stand']):
        
        if flag_subset==1:
            iStand=meta['iKeep'][iStand0]
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
    for iStand0 in range(meta['N Stand']):
        
        if flag_subset==1:
            iStand=meta['iKeep'][iStand0]
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
    
    iToFill1=np.where( (ba['SI']<0) & (ba['FIZ']==meta['LUT TIPSY']['FIZ']['I']) )[0]
    iGood=np.where( (ba['SI']>0) & (ba['FIZ']==meta['LUT TIPSY']['FIZ']['I']) )[0]
    ba['SI'][iToFill1]=np.mean(ba['SI'][iGood])
    
    iToFill2=np.where( (ba['SI']<0) & (ba['FIZ']==meta['LUT TIPSY']['FIZ']['C']) )[0]
    iGood=np.where( (ba['SI']>0) & (ba['FIZ']==meta['LUT TIPSY']['FIZ']['C']) )[0]
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

#%% ADJUST SPECIES-SPECIFIC MORTALITY

def AdjustSpeciesSpecificMortality(meta,dmeh,par,gc,iB):

    # Species affected sets
    Pest_List=['IBM','IBB','IBS','IBD','IDW']
    
    SA_List=[None]*len(Pest_List)
    for iPest in range(len(Pest_List)):
        ind=np.where( (par['DistBySC']['Name']==Pest_List[iPest]) )[0][0]
        SA_List[iPest]=np.array([par['DistBySC']['SpcCD1'][ind],par['DistBySC']['SpcCD2'][ind],
                 par['DistBySC']['SpcCD3'][ind],par['DistBySC']['SpcCD4'][ind],
                 par['DistBySC']['SpcCD5'][ind],par['DistBySC']['SpcCD6'][ind]])

    for iStand in range(meta['N Stand Full']):      
    
        for iYr in range(dmeh[iStand]['Year'].size):
        
            for iPest in range(len(Pest_List)):
        
                if dmeh[iStand]['ID_Type'][iYr]==meta['LUT Dist'][Pest_List[iPest]]:
            
                    ind_GC=int(dmeh[iStand]['ID_GC'][iB][iYr]-1)
            
                    scd=[None]*4
                    scd[0]=cbu.lut_n2s(meta['LUT VRI']['SPECIES_CD_1'],gc[iStand][iB]['s1'][ind_GC])
                    scd[1]=cbu.lut_n2s(meta['LUT VRI']['SPECIES_CD_1'],gc[iStand][iB]['s2'][ind_GC])
                    scd[2]=cbu.lut_n2s(meta['LUT VRI']['SPECIES_CD_1'],gc[iStand][iB]['s3'][ind_GC])
                    scd[3]=cbu.lut_n2s(meta['LUT VRI']['SPECIES_CD_1'],gc[iStand][iB]['s4'][ind_GC])
            
                    spct=[None]*4
                    spct[0]=gc[iStand][iB]['p1'][ind_GC]
                    spct[1]=gc[iStand][iB]['p2'][ind_GC]
                    spct[2]=gc[iStand][iB]['p3'][ind_GC]
                    spct[3]=gc[iStand][iB]['p4'][ind_GC]            
            
                    PercentAffected=0
                    for i in range(4):
                        if np.isin(scd[i],SA_List[iPest])==True:
                            PercentAffected=PercentAffected+spct[i]
                    dmeh[iStand]['MortalityFactor'][iYr]=(PercentAffected/100)*dmeh[iStand]['MortalityFactor'][iYr]        

    return dmeh
