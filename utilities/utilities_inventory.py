
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

#%% Look at contents of geodatabase

# fiona.listlayers(r'C:\Users\rhember\Documents\Data\Basemaps\Basemaps.gdb')
# fiona.listlayers(r'C:\Users\rhember\Documents\Data\ForestInventory\Conservation.gdb')
# fiona.listlayers(r'C:\Users\rhember\Documents\Data\ForestInventory\LandUse\20200706\LandUse.gdb')
# fiona.listlayers(r'C:\Users\rhember\Documents\Data\ForestInventory\Results\20210208\Results.gdb')

#%% Define strings that frequently need to be populated with zeros

#StringsToFill=['Month','Day','SILV_FUND_SOURCE_CODE','FIA_PROJECT_ID','OPENING_ID','ACTUAL_TREATMENT_AREA','ACTUAL_PLANTED_NUMBER', \
#               'PL_SPECIES_CD1','PL_SPECIES_PCT1','PL_SPECIES_GW1','PL_SPECIES_CD2','PL_SPECIES_PCT2','PL_SPECIES_GW2', \
#               'PL_SPECIES_CD3','PL_SPECIES_PCT3','PL_SPECIES_GW3','PL_SPECIES_CD4','PL_SPECIES_PCT4','PL_SPECIES_GW4', \
#               'PL_SPECIES_CD5','PL_SPECIES_PCT5','PL_SPECIES_GW5']


''' 
Update forest inventory data

 1) Manually download geodatabases from BCGW in ArcGIS (60 min)
    - if export not working, use copy and paste 
    - store files at the paths indicated in "DefineInventoryLayersAndVariables"

 2) Define inventory layers and variables (<1min))
    - Run: LayerInfo=invu.DefineInventoryLayersAndVariables()
    - make sure folder release dates are updated

 3) Build and save LUTs (20 min)
    - Run: invu.BuildForestInventoryLUTs(LayerInfo)
'''

#%% Define inventory layers and varialbes
# The "Field List" variable contains touples containing the variable name and 
# a flag indicating whether it is string (1) or numberic (0)

def DefineInventoryLayersAndVariables():

    # Define paths to geodatabase files
    PathInResultsFull=r'C:\Users\rhember\Documents\Data\ForestInventory\Results\20210208'
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
    #for iLyr in range(2,4):
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
    
    # Prepare meta['Paths']
    
    # fiona.listlayers(meta['Paths']['Results'] + '\\Results.gdb')
    
    #--------------------------------------------------------------------------
    # Query AT layer for features with missing spatial information
    # Only consider planting, direct seeding, site prep and pest control
    # 2 min
    #--------------------------------------------------------------------------

    t0=time.time()
    
    lyr=fiona.open(meta['Paths']['Results'] + '\\Results.gdb',layer='RSLT_ACTIVITY_TREATMENT_SVW')
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
    with fiona.open(meta['Paths']['Results'] + '\\Results.gdb',layer='RSLT_ACTIVITY_TREATMENT_SVW') as source:
        for feat in source:
            
            prp=feat['properties']
            
            # Opening ID
            atu_full['OPENING_ID'][cnt]=prp['OPENING_ID']
            atu_full['ACTIVITY_TREATMENT_UNIT_ID'][cnt]=prp['ACTIVITY_TREATMENT_UNIT_ID']
            
            if (prp['SILV_BASE_CODE']=='PL') & (prp['SILV_TECHNIQUE_CODE']!='SE') & (prp['SILV_TECHNIQUE_CODE']!='CG') & (prp['SILV_METHOD_CODE']!='LAYOT') & (prp['RESULTS_IND']=='Y') | \
                (prp['SILV_BASE_CODE']=='DS') & (prp['RESULTS_IND']=='Y') | \
                (prp['SILV_BASE_CODE']=='SP') & (prp['RESULTS_IND']=='Y') | \
                (prp['SILV_BASE_CODE']=='PC') & (prp['RESULTS_IND']=='Y'):
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
    gu.opickle(meta['Paths']['Results'] + '\\missing_geo_atu_list.pkl',atu_mis)
    
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
        
    with fiona.open(meta['Paths']['Results'] + '\\Results.gdb',layer='RSLT_OPENING_SVW') as source:
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
    
    gu.opickle(meta['Paths']['Results'] + '\\missing_geo_op_geos.pkl',op_mis)
    
    print((time.time()-t0)/60)
    
#    t0=time.time()
#    at_geo_from_op=[None]*atu_mis['OPENING_ID'].size
#    with fiona.open(meta['Paths']['Results'] + '\\Results.gdb',layer='RSLT_OPENING_SVW') as source:
#        for feat in source:
#            
#            # Index to AT entries with missing spatial that matches the spatial of this feature from the the opening layer
#            ind0=np.where( (atu_mis['OPENING_ID']==feat['properties']['OPENING_ID']) )[0]
#            
#            if ind0.size==0:
#                continue
#            
#            # There are not multiple instances of openings so you do not need to
#            # accumulate hits in a list
#            for i in range(ind0.size):
#                at_geo_from_op[ind0[i]]=feat['geometry']
#    
#    # Save
#    gu.opickle(meta['Paths']['Results'] + '\\at_geo_from_op.pkl',at_geo_from_op)
#    
#    print((time.time()-t0)/60)
        
    #--------------------------------------------------------------------------
    # Forest cover (only saving geometries for artificial status)
    # Too big if you save all geometries
    # Takes 5 min
    #--------------------------------------------------------------------------
    
#    t0=time.time()
#    
#    lyr=fiona.open(meta['Paths']['Results'] + '\\Results.gdb',layer='RSLT_FOREST_COVER_INV_SVW')
#    L=len(lyr)    
#    
#    at_geo_from_fc=[None]*atu_mis['OPENING_ID'].size
#    ref_year_fc=[None]*atu_mis['OPENING_ID'].size   
#    year_updated_fc=[None]*atu_mis['OPENING_ID'].size   
#
#    with fiona.open(meta['Paths']['Results'] + '\\Results.gdb',layer='RSLT_FOREST_COVER_INV_SVW') as source:
#        for feat in source:
#            
#            if (feat['geometry']==None):
#                continue
#            
#            if (feat['properties']['STOCKING_TYPE_CODE']!='ART'):
#                continue
#            
#            # Index to AT entries that correspond to this feature and have missing geometries
#            ind0=np.where( (atu_mis['OPENING_ID']==feat['properties']['OPENING_ID']) )[0]
#            
#            if ind0.size==0:
#                continue
#            
#            YearUpdated=0
#            if feat['properties']['FOREST_COVER_WHEN_UPDATED']!=None:
#                YearUpdated=int(feat['properties']['FOREST_COVER_WHEN_UPDATED'][0:4]) 
#            
#            for i in range(ind0.size):
#                
#                # There can be multiple instances so accumulate hits in a list
#                if at_geo_from_fc[ind0[i]]==None:                    
#                    #at_geo_from_fc[ind0[i]]=feat['geometry']
#                    at_geo_from_fc[ind0[i]]=[feat['geometry']]
#                    ref_year_fc[ind0[i]]=[feat['properties']['REFERENCE_YEAR']]
#                    year_updated_fc[ind0[i]]=[YearUpdated]
#                else:
#                    at_geo_from_fc[ind0[i]].append(feat['geometry'])
#                    ref_year_fc[ind0[i]].append(feat['properties']['REFERENCE_YEAR'])
#                    year_updated_fc[ind0[i]].append(YearUpdated)
#    
#    # Save
#    gu.opickle(meta['Paths']['Results'] + '\\at_geo_from_fcinv.pkl',at_geo_from_fc)
#    gu.opickle(meta['Paths']['Results'] + '\\ref_year_fc.pkl',ref_year_fc)
#    gu.opickle(meta['Paths']['Results'] + '\\year_updated_fc.pkl',year_updated_fc)
#            
#    print((time.time()-t0)/60)    
    
    #--------------------------------------------------------------------------
    # Forest cover by OPENING_ID (only saving geometries for artificial status)
    # Takes 5 min
    #--------------------------------------------------------------------------
    
    t0=time.time()
    
    u=np.unique(atu_mis['OPENING_ID'])    
    fc_mis={}
    for i in range(u.size):
        fc_mis[u[i]]=[]
        
    with fiona.open(meta['Paths']['Results'] + '\\Results.gdb',layer='RSLT_FOREST_COVER_INV_SVW') as source:
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
    
    gu.opickle(meta['Paths']['Results'] + '\\missing_geo_fc_geos.pkl',fc_mis)
    
    print((time.time()-t0)/60)    

    
    #--------------------------------------------------------------------------
    # Find matching FC geometries for each AT entry with missing spatial
    # Takes 5 min
    #--------------------------------------------------------------------------
    
    atu_mis=gu.ipickle(meta['Paths']['Results'] + '\\missing_geo_atu_list.pkl')
    fc_mis=gu.ipickle(meta['Paths']['Results'] + '\\missing_geo_fc_geos.pkl')
    
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
            
    gu.opickle(meta['Paths']['Results'] + '\\missing_geo_atu_list.pkl',atu_mis)

##%% Scraps from toubleshooting missing geometries:
#    
#    #%% Analyze missing planting spatial
#
#c=[0,1,0,1]; #c=[0,1,2,3,0,1,2,3,0,1,2,3]
#r=[0,0,1,1]; #r=[0,0,0,0,1,1,1,1,2,2,2,2]
#
## Import ATU data
#atu_mis=gu.ipickle(meta['Paths']['Results'] + '\\atu_mis.pkl')
#
#at_geo_from_op=gu.ipickle(meta['Paths']['Results'] + '\\at_geo_from_op.pkl')
#
#garc.collect()
#
#
## 5 percent where opening spatial cannot be found for AT entries
#N=0 
#for i in range(len(at_geo_from_op)):
#    if at_geo_from_op[i]==None:
#        N=N+1
#N/atu_mis['OPENING_ID'].size
#
#
## i=150 # Easy
## i=2300 # Easy
#
#for i in range(atu_mis['OPENING_ID'].size):
#
#    if at_geo_from_op[i]==None:
#        continue
#    
#    # i=2300
#
#    #print(atu_mis['OPENING_ID'][i])
#    
#    feat_op={}; feat_op['properties']=prp; feat_op['geometry']=at_geo_from_op[i]
#    gdf_op=gpd.GeoDataFrame.from_features([feat_op])
#    A_op=np.round(gdf_op.area.values[0]/1000)
#    
#    ind_atu=np.where( (atu_mis['OPENING_ID']==atu_mis['OPENING_ID'][i]) & (atu_mis['SBC']=='PL') )[0]
#    atu_mis['SBC'][ind_atu]
#    atu_mis['Year'][ind_atu]
#    A_at=atu_mis['ACTUAL_TREATMENT_AREA'][ind_atu]        
#
#    if at_geo_from_fc[i]==None:
#        continue
#
#    A_fc=np.zeros(len(at_geo_from_fc[i]))
#    Year_fc=np.zeros(len(at_geo_from_fc[i]))
#    for j in range(len(at_geo_from_fc[i])):
#        feat_fc={}; feat_fc['properties']=prp; feat_fc['geometry']=at_geo_from_fc[i][j]
#        gdf_fc=gpd.GeoDataFrame.from_features([feat_fc])
#        A_fc[j]=np.round(gdf_fc.area.values[0]/10000)
#        Year_fc[j]=ref_year_fc[i][j]
#
#    if ind_atu.size==A_fc.size:
#        dist=np.abs(A_at[:, np.newaxis]-A_fc)
#        potentialClosest=dist.argmin(axis=1)
#        sD[i]=np.sum(A_at-A_fc[potentialClosest])
#
#
#    idx=-999*np.ones(ind_atu.size)
#    for j in range(ind_atu.size):
#        if (atu_mis['SBC'][ind_atu[j]]=='PL'):
#            ind1=np.where(A_fc==atu_mis['ACTUAL_TREATMENT_AREA'][ind_atu[j]])[0]
#            if ind1.size>0:
#                idx[j]=ind1
#
#
#    plt.close('all')
#    fig,ax=plt.subplots(1,figsize=gu.cm2inch(12,9))
#    gdf_op.plot(ax=ax,facecolor='None',linewidth=3,edgecolor='k',linestyle='-')
#    for j in range(len(at_geo_from_fc[i])):
#        gdf_op.plot(ax=ax,facecolor='None',linewidth=3,edgecolor='k',linestyle='-') 
#        feat_fc={}; feat_fc['properties']=prp; feat_fc['geometry']=at_geo_from_fc[i][j]
#        gdf_fc=gpd.GeoDataFrame.from_features([feat_fc])
#        gdf_fc.plot(ax=ax,facecolor='g',edgecolor='y',linewidth=1,linestyle='-') 
#          
#
#    plt.close('all')
#    fig,ax=plt.subplots(2,2,figsize=gu.cm2inch(12,9))
#    for j in range(4):
#        gdf_op.plot(ax=ax[r[j],c[j]],facecolor='None',linewidth=3,edgecolor='k',linestyle='-')
#        ax[r[j],c[j]].set_yticklabels(''); ax[r[j],c[j]].set_xticklabels('')
#    if at_geo_from_fc[i]!=None:
#        Atot_fc=0
#        for j in range(4):
#            if j+1>len(at_geo_from_fc[i]):
#                continue
#            gdf_op.plot(ax=ax[r[j],c[j]],facecolor='None',linewidth=3,edgecolor='k',linestyle='-') 
#            feat_fc={}; feat_fc['properties']=prp; feat_fc['geometry']=at_geo_from_fc[i][j]
#            gdf_fc=gpd.GeoDataFrame.from_features([feat_fc])
#            gdf_fc.plot(ax=ax[r[j],c[j]],facecolor='None',edgecolor='c',linewidth=2,linestyle='-') 
#            ax[r[j],c[j]].set_title('Year: ' + str(ref_year_fc[i][j]) + ', Area: ' + str(np.round(gdf_fc.area.values[0]/10000)) + ' ha',fontsize=7)
#            ax[r[j],c[j]].set_xticks([]); ax[r[j],c[j]].set_yticks([])
#            Atot_fc=Atot_fc+gdf_fc.area.values[0]/10000
#    plt.tight_layout()
#    gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Planting\ReconMissingGeom_Eg1','png',900)

    
#%% ADD PLANTING INFO TO DMEC

# Create function to avoid duplication of bulky code
def AddPlantingWithNoData(d_nd):
    for i in range(1,6):
        d_nd['PL_SPECIES_CD' + str(i)]=np.append(d_nd['PL_SPECIES_CD' + str(i)],-999)
        d_nd['PL_SPECIES_PCT' + str(i)]=np.append(d_nd['PL_SPECIES_PCT' + str(i)],-999)
        d_nd['PL_SPECIES_GW' + str(i)]=np.append(d_nd['PL_SPECIES_GW' + str(i)],-999)
    return d_nd

#%% PREPARE DISTURBANCE MANAGEMENT ENVENT HISTORY

def PrepDMEC(idx,meta,par,atu,pl,op,fcinv,vri,cut,fire,burnsev,pest):
    
    # Initiate disturbance-management event history
    dmec=[None]*meta['N Stand']
    
    # Specify flag indicating whether subsetting occurs
    if np.isin('iKeep',list(meta.keys()))==True:
        # Tile project, only keeping a subset           
        flag_subset=1
    else:
        # Run all entries
        flag_subset=0
    
    for iStand0 in range(meta['N Stand']):
        
        # Index to stand depends on subsetting
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
                ind=np.where( (np.floor(dmec0['Year'])==np.floor(fire['FIRE_YEAR'][indS[i]])) & 
                             (dmec0['ID_Type']==meta['LUT']['Dist']['Wildfire']) )[0]
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
                    ind=np.where( (par['DistBySC']['Name']=='IBM') & (par['DistBySC']['SeverityCD']==sev_s) )[0]
                    dmec0['MortalityFactor']=np.append(dmec0['MortalityFactor'],par['DistBySC']['MortalityFactor'][ind])
                    dmec0['GrowthFactor']=np.append(dmec0['GrowthFactor'],par['DistBySC']['GrowthFactor'][ind])
                
                # Western Balsam Bark Beetle
                # (https://www2.gov.bc.ca/gov/content/industry/forestry/managing-our-forest-resources/forest-health/forest-pests/bark-beetles/western-balsam-bark-beetle)
                elif (psp==meta['LUT']['Pest']['PEST_SPECIES_CODE']['IBB']):
                    dmec0['ID_Type']=np.append(dmec0['ID_Type'],meta['LUT']['Dist']['IBB'])
                    ind=np.where( (par['DistBySC']['Name']=='IBB') & (par['DistBySC']['SeverityCD']==sev_s) )[0]
                    dmec0['MortalityFactor']=np.append(dmec0['MortalityFactor'],par['DistBySC']['MortalityFactor'][ind])
                    dmec0['GrowthFactor']=np.append(dmec0['GrowthFactor'],par['DistBySC']['GrowthFactor'][ind])
            
                # Douglas-fir beetle
                # (https://www2.gov.bc.ca/gov/content/industry/forestry/managing-our-forest-resources/forest-health/forest-pests/bark-beetles/western-balsam-bark-beetle)
                elif (psp==meta['LUT']['Pest']['PEST_SPECIES_CODE']['IBD']):
                    dmec0['ID_Type']=np.append(dmec0['ID_Type'],meta['LUT']['Dist']['IBD'])
                    ind=np.where( (par['DistBySC']['Name']=='IBD') & (par['DistBySC']['SeverityCD']==sev_s) )[0]
                    dmec0['MortalityFactor']=np.append(dmec0['MortalityFactor'],par['DistBySC']['MortalityFactor'][ind])
                    dmec0['GrowthFactor']=np.append(dmec0['GrowthFactor'],par['DistBySC']['GrowthFactor'][ind])
                
                # Spruce beetle
                # (https://www2.gov.bc.ca/gov/content/industry/forestry/managing-our-forest-resources/forest-health/forest-pests/bark-beetles/western-balsam-bark-beetle)
                elif (psp==meta['LUT']['Pest']['PEST_SPECIES_CODE']['IBS']):
                    dmec0['ID_Type']=np.append(dmec0['ID_Type'],meta['LUT']['Dist']['IBS'])
                    ind=np.where( (par['DistBySC']['Name']=='IBS') & (par['DistBySC']['SeverityCD']==sev_s) )[0]
                    dmec0['MortalityFactor']=np.append(dmec0['MortalityFactor'],par['DistBySC']['MortalityFactor'][ind])
                    dmec0['GrowthFactor']=np.append(dmec0['GrowthFactor'],par['DistBySC']['GrowthFactor'][ind])
                
                # Western spruce budworm 
                # Populate severity with the ID for severity - it will be revised below
                # to model mortality that occurs from repeated infestation
                elif (psp==meta['LUT']['Pest']['PEST_SPECIES_CODE']['IDW']):
                    dmec0['ID_Type']=np.append(dmec0['ID_Type'],meta['LUT']['Dist']['IDW'])
                    ind=np.where( (par['DistBySC']['Name']=='IDW') & (par['DistBySC']['SeverityCD']==sev_s) )[0]
                    dmec0['MortalityFactor']=np.append(dmec0['MortalityFactor'],sev)
                    dmec0['GrowthFactor']=np.append(dmec0['GrowthFactor'],par['DistBySC']['GrowthFactor'][ind])
                
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
    
    # Open connection to parameter database    
    par=gu.ipickle(meta['Paths']['Model Code'] + '\\Parameters\\Parameters.pkl')
    
    # Import distubance type        
    meta['LUT']['Dist']={}
    for i in range(len(par['Disturbances']['Name'])):
        meta['LUT']['Dist'][par['Disturbances']['Name'][i]]=par['Disturbances']['ID'][i]
    
    # BGC zone     
    #LUT_BGC_Zone={}
    #for i in range(len(par['BGC_ZONE']['CODE_BGC_ZONE'])):
    #    LUT_BGC_Zone[par['BGC_ZONE']['CODE_BGC_ZONE'][i]]=par['BGC_ZONE']['ID_BGC_ZONE'][i]
    
    # Added this to accommodate jupyter notebook demos - will need updating periodically
    if 'Results' not in meta['Paths']:
        meta['Paths']['Results']=r'C:\Users\rhember\Documents\Data\ForestInventory\Results\20210208'
    if 'VRI' not in meta['Paths']:
        meta['Paths']['VRI']=r'C:\Users\rhember\Documents\Data\ForestInventory\VRI\20200430'    
    if 'Disturbances' not in meta['Paths']:
        meta['Paths']['Disturbances']=r'C:\Users\rhember\Documents\Data\ForestInventory\Disturbances\20200430'
    if 'LandUse' not in meta['Paths']:
        meta['Paths']['LandUse']=r'C:\Users\rhember\Documents\Data\ForestInventory\LandUse\20200706'    
    
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
    meta['LUT']['SRS']={}
    for i in range(len(par['SRS']['SRS_CD'])):
        meta['LUT']['SRS'][par['SRS']['SRS_CD'][i]]=par['SRS']['SRS_ID'][i]
    
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
    for iStand in range(meta['N Stand Full']):
        if dmec[iStand]==None:
            continue
        for key in meta['LUT']['Dist'].keys():
            ind=np.where(dmec[iStand]['ID_Type']==meta['LUT']['Dist'][key])[0]
            if ind.size==0:
                continue
            uYear=np.unique(np.floor(dmec[iStand]['Year'][ind]))
            for iYear in range(uYear.size):
                ind=np.where( (dmec[iStand]['ID_Type']==meta['LUT']['Dist'][key]) & (np.floor(dmec[iStand]['Year'])==uYear[iYear]) )[0]
                dmec[iStand]['ID_Type'][ind[1:]]=-999
    return dmec

#%% EXCLUDE UNIDENTIFIED EVENTS

def Exclude_Unidentified_Events(meta,dmec):

    for iStand in range(meta['N Stand Full']):
        if dmec[iStand]==None:
            continue
        ind=np.where(dmec[iStand]['ID_Type']!=-999)[0]
        for key in dmec[iStand]:
            dmec[iStand][key]=dmec[iStand][key][ind]
    return dmec

#%% REMOVE SLASHPILE BURNS IN SELECT BGC ZONES

def Remove_SlashpileBurns_From_Select_Zones(meta,dmec,ba):

    for iStand in range(meta['N Stand Full']):
        
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
    for iStand in range(meta['N Stand Full']):
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

def Ensure_Fert_Preceded_By_Disturbance(meta,dmec,th_sev_last_dist,AgeAtFert,StringsToFill):

    ListOfTestedDist=[meta['LUT']['Dist']['Wildfire'],meta['LUT']['Dist']['Harvest'],
            meta['LUT']['Dist']['Knockdown'],meta['LUT']['Dist']['Salvage Logging'],
            meta['LUT']['Dist']['Beetles'],meta['LUT']['Dist']['IBM'],meta['LUT']['Dist']['IBB'],
            meta['LUT']['Dist']['IBD'],meta['LUT']['Dist']['IBS']]
    
    for iStand in range(meta['N Stand Full']):
        
        if dmec[iStand]==None:
            continue
    
        iA=np.where( (dmec[iStand]['ID_Type']==meta['LUT']['Dist']['Fertilization Aerial']) )[0]
        if iA.size==0: 
            continue
        
        # Index to events prior to first fertilization with 100% mortality
        ind=np.where( (dmec[iStand]['Year']<=dmec[iStand]['Year'][iA[0]]) & (dmec[iStand]['MortalityFactor']==100) & np.isin(dmec[iStand]['ID_Type'],ListOfTestedDist) )[0]
        
        #dYear=dmec[iStand]['Year'][iA[0]]-dmec[iStand]['Year'][ind[-1]]
        if (ind.size==0):
            
            # Add harvest
            Year=dmec[iStand]['Year'][iA[0]]-AgeAtFert
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

#%% IDW FIX SEVERITY
# The dmec was populated with the numeric severity ID. Mortality only occurs 
# following repeated outrbreak years. 

def IDW_Fix_Severity(meta,dmec,par):

    for iStand in range(meta['N Stand Full']):
        
        if dmec[iStand]==None:
            continue
        
        # Save a frozen version of the initial severity
        Severity_Frozen=dmec[iStand]['MortalityFactor'].copy()
    
        for iA in range(dmec[iStand]['Year'].size):
            
            if dmec[iStand]['ID_Type'][iA]=='IDW':
            
                # By default, apply mortality rate in the year of defoliation
                sev1=dmec[iStand]['MortalityFactor'][iA]
                sev_s1=cbu.lut_n2s(meta['LUT']['Pest']['PEST_SEVERITY_CODE'],sev1)[0]
                ind=np.where( (par['DistBySC']['Name']=='IDW') & (par['DistBySC']['SeverityCD']==sev_s1) )[0]
                Mortality1=par['DistBySV']['MortalityFactor'][ind]

                if iA>0:
                    if dmec[iStand]['ID_Type'][iA-1]=='IDW':
                    
                        # If it is back to back infestation, adjust mortality accordingly
                        sev0=Severity_Frozen[iA-1]
                        sev_s0=cbu.lut_n2s(meta['LUT']['Pest']['PEST_SEVERITY_CODE'],sev0)[0]
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
                dmec[iStand]['MortalityFactor'][iA]=np.array(Mortality1,dtype='int16')
    return dmec

#%% ADD OLDEST KNOWN DISTURBANCE FROM VRI
# This really doesn't work well. RETIRED!!

def Add_Oldest_Disturbance_From_VRI(meta,dmec,idx,vri):
    
    for iStand in range(meta['N Stand Full']):
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
    for iStand in range(meta['N Stand Full']):
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
    for iStand in range(meta['N Stand Full']):   
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
    ba['FIZ']=meta['LUT']['TIPSY']['FIZ']['I']*np.ones(meta['N Stand'])
    ba['BEC_ZONE_CODE']=meta['LUT']['VRI']['BEC_ZONE_CODE']['SBS']*np.ones(meta['N Stand'])
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
        
        if (vri['BEC_ZONE_CODE'][ind0]==meta['LUT']['VRI']['BEC_ZONE_CODE']['CWH']) | \
            (vri['BEC_ZONE_CODE'][ind0]==meta['LUT']['VRI']['BEC_ZONE_CODE']['CDF']) | \
            (vri['BEC_ZONE_CODE'][ind0]==meta['LUT']['VRI']['BEC_ZONE_CODE']['MH']):
            ba['FIZ'][iStand0]=meta['LUT']['TIPSY']['FIZ']['C']
        else:
            ba['FIZ'][iStand0]=meta['LUT']['TIPSY']['FIZ']['I']
    
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
    spl['SI_SPL']=-999*np.ones(meta['N Stand'])
    
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
        for iStand0 in range(meta['N Stand']):
        
            if flag_subset==1:
                iStand=meta['iKeep'][iStand0]
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
    for iStand in range(meta['N Stand']):
        for iScn in range(meta['N Scenario']):
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

def AdjustSpeciesSpecificMortality(meta,dmec,par,gc,iB):

    # Species affected sets
    Pest_List=['IBM','IBB','IBS','IBD','IDW']
    
    SA_List=[None]*len(Pest_List)
    for iPest in range(len(Pest_List)):
        ind=np.where( (par['DistBySC']['Name']==Pest_List[iPest]) )[0][0]
        SA_List[iPest]=np.array([par['DistBySC']['SpcCD1'][ind],par['DistBySC']['SpcCD2'][ind],
                 par['DistBySC']['SpcCD3'][ind],par['DistBySC']['SpcCD4'][ind],
                 par['DistBySC']['SpcCD5'][ind],par['DistBySC']['SpcCD6'][ind]])

    for iStand in range(meta['N Stand Full']):      
    
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
    for iStand in range(meta['N Stand Full']):
        d=dmec[iStand].copy()
        ord=np.argsort(d['Year'])
        for key in d.keys():
            d[key]=d[key][ord]
        dmec[iStand]=d.copy()
    return dmec

#%% Export AT Layer data to spreadhseet
    
def ExportSummaryByGridCell(meta,atu):
    
    n=atu['Year'].size+50000
    fl=np.zeros(50000)
    
    d={}
    d['IdxToSXY']=np.append(atu['IdxToSXY'],fl)
    d['ID_Multipolygon']=np.append(atu['IdxToSXY'],fl)
    d['Year']=np.append(atu['Year'],fl)
    d['Month']=np.append(atu['Month'],fl)
    d['OPENING_ID']=np.append(atu['OPENING_ID'],fl)
    d['Activity_Type']=np.array(['empty' for _ in range(n)],dtype=object)  
    
    d['AEF_ATU']=np.zeros(n)    
    for i in range(len(atu_multipolygons)):
        ind1=np.where(sxy['ID_atu_multipolygons']==i)[0]
        nxy=ind1.size
        A=atu_multipolygons[i]['ACTUAL_TREATMENT_AREA']
        if A==None:
            A=0.00001
        for j in range(ind1.size):
            ind2=np.where(atu['IdxToSXY']==ind1[j])[0]
            d['AEF_ATU'][ind2]=np.round(A/nxy,3)
            d['ID_Multipolygon'][ind2]=i    
        
    d['FIA_PROJECT_ID']=np.array(['empty' for _ in range(n)],dtype=object)    
    u=np.unique(atu['FIA_PROJECT_ID'])
    for i in range(u.size):
        ind=np.where(atu['FIA_PROJECT_ID']==u[i])[0]
        d['FIA_PROJECT_ID'][ind]=cbu.lut_n2s(meta['LUT']['ATU']['FIA_PROJECT_ID'],u[i])    
    #d['FIA_PROJECT_ID']=np.array(['empty' for _ in range(n)],dtype=object)    
    #for i in range(atu['Year'].size):
    #    d['FIA_PROJECT_ID'][i]=cbu.lut_n2s(meta['LUT']['ATU']['FIA_PROJECT_ID'],atu['FIA_PROJECT_ID'][i])[0]   
        
    d['FSC']=np.array(['empty' for _ in range(n)],dtype=object)    
    for i in range(atu['Year'].size):
        d['FSC'][i]=cbu.lut_n2s(meta['LUT']['ATU']['SILV_FUND_SOURCE_CODE'],atu['SILV_FUND_SOURCE_CODE'][i])[0]    
    d['DistCD']=np.array(['empty' for _ in range(n)],dtype=object)    
    for i in range(atu['Year'].size):
        d['DistCD'][i]=cbu.lut_n2s(meta['LUT']['ATU']['DISTURBANCE_CODE'],atu['DISTURBANCE_CODE'][i])[0]    
    d['SBC']=np.array(['empty' for _ in range(n)],dtype=object)
    for i in range(atu['Year'].size):
        d['SBC'][i]=cbu.lut_n2s(meta['LUT']['ATU']['SILV_BASE_CODE'],atu['SILV_BASE_CODE'][i])[0]
    d['SMC']=np.array(['empty' for _ in range(n)],dtype=object)
    for i in range(atu['Year'].size):
        d['SMC'][i]=cbu.lut_n2s(meta['LUT']['ATU']['SILV_METHOD_CODE'],atu['SILV_METHOD_CODE'][i])[0]
    d['STC']=np.array(['empty' for _ in range(n)],dtype=object)
    for i in range(atu['Year'].size):
        d['STC'][i]=cbu.lut_n2s(meta['LUT']['ATU']['SILV_TECHNIQUE_CODE'],atu['SILV_TECHNIQUE_CODE'][i])[0]
    
    d['BGCz']=np.array(['empty' for _ in range(n)],dtype=object)
    d['BGCsz']=np.array(['empty' for _ in range(n)],dtype=object)
    d['BGCv']=np.array(['empty' for _ in range(n)],dtype=object)
    
    d['SI_FCinv']=np.zeros(n)
    d['I_SPH_FCinv']=np.zeros(n)
    d['I_Spc1_CD']=np.array(['empty' for _ in range(n)],dtype=object)
    d['I_Spc2_CD']=np.array(['empty' for _ in range(n)],dtype=object)
    
    d['Pl_SPH']=np.append(np.round(atu['ACTUAL_PLANTED_NUMBER']/atu['ACTUAL_TREATMENT_AREA']),fl)
    ind=np.where( (d['SBC']!='PL') & (d['STC']!='PL') )[0]
    d['Pl_SPH'][ind]=0
    
    cnt=atu['Year'].size
    for i in range(fcinv['IdxToSXY'].size):
        ind=np.where( (d['IdxToSXY']==fcinv['IdxToSXY'][i]) )[0]
        if ind.size==0:
            # Don't add forest cover where there is no ATU info for that grid cell 
            # *** I have no idea why this happens!! ***
            continue
        
        ind=np.where( (d['IdxToSXY']==fcinv['IdxToSXY'][i]) & (d['OPENING_ID']==fcinv['OPENING_ID'][i]) & (d['Year']==fcinv['REFERENCE_YEAR'][i]) )[0]
        if ind.size>0:
            for j in range(ind.size):
                d['I_SPH_FCinv'][ind[j]]=fcinv['I_TOTAL_STEMS_PER_HA'][i]
                d['SI_FCinv'][ind[j]]=fcinv['SITE_INDEX'][i]
                d['I_Spc1_CD'][ind[j]]=cbu.lut_n2s(meta['LUT']['VRI']['SPECIES_CD_1'],fcinv['I_SPECIES_CODE_1'][i])
                d['I_Spc2_CD'][ind[j]]=cbu.lut_n2s(meta['LUT']['VRI']['SPECIES_CD_1'],fcinv['I_SPECIES_CODE_2'][i])
        else:
            d['IdxToSXY'][cnt]=fcinv['IdxToSXY'][i]
            d['Year'][cnt]=fcinv['REFERENCE_YEAR'][i]
            d['OPENING_ID'][cnt]=fcinv['OPENING_ID'][i]
            d['I_SPH_FCinv'][cnt]=fcinv['I_TOTAL_STEMS_PER_HA'][i]
            d['SI_FCinv'][cnt]=fcinv['SITE_INDEX'][i]
            d['I_Spc1_CD'][cnt]=cbu.lut_n2s(meta['LUT']['VRI']['SPECIES_CD_1'],fcinv['I_SPECIES_CODE_1'][i])
            d['I_Spc2_CD'][cnt]=cbu.lut_n2s(meta['LUT']['VRI']['SPECIES_CD_1'],fcinv['I_SPECIES_CODE_2'][i])
            cnt=cnt+1
    
    # Get rid of empty values
    for k in d.keys():
        d[k]=d[k][0:cnt]
      
    # Add VRI
    for i in range(d['IdxToSXY'].size):
        ind=np.where(vri['IdxToSXY']==d['IdxToSXY'][i])[0]
        if ind.size==0:
            continue
        d['BGCz'][i]=cbu.lut_n2s(meta['LUT']['VRI']['BEC_ZONE_CODE'],vri['BEC_ZONE_CODE'][ind[0]])[0]
        d['BGCsz'][i]=cbu.lut_n2s(meta['LUT']['VRI']['BEC_SUBZONE'],vri['BEC_SUBZONE'][ind[0]])[0]
        d['BGCv'][i]=cbu.lut_n2s(meta['LUT']['VRI']['BEC_VARIANT'],vri['BEC_VARIANT'][ind[0]])[0]    
        
    # Add planting info    
    for i in range(20):
        d['Pl_Spc' + str(i+1) + '_CD']=np.array([' ' for _ in range(d['IdxToSXY'].size)],dtype=object)
        d['Pl_Spc' + str(i+1) + '_Pct']=np.array([' ' for _ in range(d['IdxToSXY'].size)],dtype=object)
        d['Pl_Spc' + str(i+1) + '_NumTree']=np.array([' ' for _ in range(d['IdxToSXY'].size)],dtype=object)
        d['Pl_Spc' + str(i+1) + '_SeedLot']=np.array([' ' for _ in range(d['IdxToSXY'].size)],dtype=object)
    
    for i in range(d['IdxToSXY'].size):
        if (d['SBC'][i]!='PL') & (d['STC'][i]!='PL'):
            continue        
        ind=np.where( (pl['IdxToSXY']==d['IdxToSXY'][i]) & (pl['OPENING_ID']==d['OPENING_ID'][i]) & (pl['Year']==d['Year'][i]) )[0]
        tot_pl=np.sum(pl['NUMBER_PLANTED'][ind])
        if ind.size>0:
            Ord=np.flip(np.argsort(pl['NUMBER_PLANTED'][ind]))
            for j in range(ind.size):
                if j>19:
                    continue
                ind0=ind[Ord[j]]
                d['Pl_Spc' + str(j+1) + '_CD'][i]=cbu.lut_n2s(meta['LUT']['PL']['SILV_TREE_SPECIES_CODE'],pl['SILV_TREE_SPECIES_CODE'][ind0])[0]
                d['Pl_Spc' + str(j+1) + '_Pct'][i]=np.round(pl['NUMBER_PLANTED'][ind0]/tot_pl*100)
                d['Pl_Spc' + str(j+1) + '_NumTree'][i]=pl['NUMBER_PLANTED'][ind0]
                d['Pl_Spc' + str(j+1) + '_SeedLot'][i]=pl['SEEDLOT_NUMBER'][ind0]
    
    # Activity type
    #flg_source='FCI'
    flg_source='ReforestationNonOb'
    
    if flg_source=='FCI':
        
        u=np.unique(d['FIA_PROJECT_ID'])
        for i in range(u.size):
            ind1=np.where(d['FIA_PROJECT_ID']==u[i])[0]
            if ind1.size==0:
                continue
            ind2=np.where(dAdmin['PP Number']==u[i])[0]
            if ind2.size==0:
                continue
            u2=np.unique(d['IdxToSXY'][ind1])
            for j in range(u2.size):
                ind3=np.where(d['IdxToSXY']==u2[j])[0]
                for k in range(ind3.size):
                    d['Activity_Type'][ind3[k]]=dAdmin['Project Type'][ind2[0]]
    
    elif flg_source=='ReforestationNonOb':
        
        nam=['No Planting','SL','KD','UNDER','Unidentified']
        for i in range(d['IdxToSXY'].size):
            d['Activity_Type'][i]=nam[int(meta['ProjectType'][int(d['IdxToSXY'][i])])]
       
    # Save
    df=pd.DataFrame.from_dict(d)    
    df=df.sort_values(by=['IdxToSXY','Year','Month'])
    df.to_excel(meta['Paths']['Project'] + '\\Inputs\\SummarySiteAndEventsByGridCell.xlsx',index=False)
    
    return

#%% Timber harvesting land base 

def DefineTHLB(meta,ba,dmec,fcres,lul,ogmal,park):
    
    #------------------------------------------------------------------------------
    # Initialize THLB flags (THLB=1,Non-THLB=0)
    #------------------------------------------------------------------------------

    # Initially assume everything is in the THLB
    thlb_flag_Actual=np.ones((meta['N Time'],meta['N Stand Full']))

    # Define second THLB with areas that have been removed due to value diversification
    thlb_flag_Baseline=thlb_flag_Actual.copy()

    # Index to stands that are uneconomic
    iUneconomic=np.where(ba['SI']<=5)[0]

    # Remove uneconomic stands from THLB
    thlb_flag_Actual[:,iUneconomic]=0
    thlb_flag_Baseline[:,iUneconomic]=0

    # Idenify stands that have been harvested
    has_been_harvested=np.zeros(meta['N Stand Full'])
    for i in range(len(dmec)):
        ind=np.where(dmec[i]['ID_Type']==meta['LUT']['Dist']['Harvest'])[0]
        if ind.size>0:
            has_been_harvested[i]=1

    # Index to stands that have not been harvested
    iNoHarv=np.where( (has_been_harvested==0) & (ba['SI']>5) )[0]

    # Use the ratio of THLB to non-THLB as an indicator of what will be harvested
    # among remaining primary forest
    ratio_thlb=22/55 # ratio of THLB to total forest (SOF)

    corr=iUneconomic.size/meta['N Stand Full']

    # Probability of evading harvest
    if iNoHarv.size>0:
        p_evade=(1-ratio_thlb-corr)*(meta['N Stand Full']/iNoHarv.size)
    else:
        p_evade=(1-ratio_thlb-corr)

    # Random prediction of whether it will evade harvesting
    iRem=np.where(np.random.random(iNoHarv.size)<p_evade)[0]
    thlb_flag_Actual[:,iNoHarv[iRem]]=0
    thlb_flag_Baseline[:,iNoHarv[iRem]]=0

    # np.sum(thlb_flag_Actual[0,:])/meta['N Stand Full']

    #------------------------------------------------------------------------------
    # Define the year of transition from THLB to non-THLB for specific LU types
    #------------------------------------------------------------------------------

    # Look at abundance of each LUL type
    d=meta['LUT']['LU L']['LEGAL_FEAT_OBJECTIVE'].copy()
    for k in meta['LUT']['LU L']['LEGAL_FEAT_OBJECTIVE'].keys():
        ind=np.where( (lul['LEGAL_FEAT_OBJECTIVE']==meta['LUT']['LU L']['LEGAL_FEAT_OBJECTIVE'][k]) )[0]
        d[k]=ind.size
    ds={k: v for k,v in sorted(d.items(), key=lambda item: item[1])}

    # List of legal land use plan types that transition to conservatoin (this needs to be vetted by experts)
    ListCon=['Visually Sensitive Areas','Connectivity Corridors','Scenic Areas','Caribou Winter Habitat Zones',
      'Tourism Areas','Visual Quality','High Value Grizzly Bear Habitat','No Timber Harvesting Areas','Landscape Corridors',
      'Critical Deer Winter Range','Sensitive Watershed','Water Management Units','High Value Wetlands for Moose','Telkwa Caribou Recovery Area',
      'SRMZ3:Caribou Migration Corridor Sub-Zone','Caribou Migration Corridor','High Biodiversity Emphasis Areas','Scenic Corridors']

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
    print(str(A/meta['N Stand Full']*100) + '% of the sample transitioned to non-THLB.')

    # Look at reserves
    N={}
    for k in meta['LUT']['FC_R']['SILV_RESERVE_CODE'].keys():
        ind=np.where(fcres['SILV_RESERVE_CODE']==meta['LUT']['FC_R']['SILV_RESERVE_CODE'][k])[0]
        N[k]=ind.size

    # Initialize year of transition
    thlb_YearTransitionOut=np.zeros(meta['N Stand Full'])

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
            thlb_flag_Actual[it,j]=0

    # Adjust the baseline so that simulated harvesting between 1995 and 2020 only
    # occurs in areas where the THLB was affected by value diversification
    for year in range(1990,2021,1):
        iT=np.where(meta['Year']==year)[0]
        iS=np.where( (thlb_flag_Baseline[iT,:]==1) & (thlb_flag_Actual[iT,:]==1) )[1]
        thlb_flag_Baseline[iT,iS]=0

    flg=0
    if flg==1:
        plt.figure(2)
        plt.plot(np.sum(thlb_flag_Actual,axis=1)/meta['N Stand Full'])
        plt.plot(np.sum(thlb_flag_Baseline,axis=1)/meta['N Stand Full'],'--')
        
    return thlb_flag_Actual,thlb_flag_Baseline

#%% Load sparse geospatiatial inputs

def LoadSparseGeospatialInputs(meta):
    
    sxy=gu.ipickle(meta['Paths']['Geospatial'] + '\\sxy.pkl')
    atu=gu.ipickle(meta['Paths']['Geospatial'] + '\\RSLT_ACTIVITY_TREATMENT_SVW.pkl')
    op=gu.ipickle(meta['Paths']['Geospatial'] + '\\RSLT_OPENING_SVW.pkl')
    burnsev=gu.ipickle(meta['Paths']['Geospatial'] + '\\VEG_BURN_SEVERITY_SP.pkl')
    vri=gu.ipickle(meta['Paths']['Geospatial'] + '\\VEG_COMP_LYR_R1_POLY.pkl')
    fcinv=gu.ipickle(meta['Paths']['Geospatial'] + '\\RSLT_FOREST_COVER_INV_SVW.pkl')
    fcres=gu.ipickle(meta['Paths']['Geospatial'] + '\\RSLT_FOREST_COVER_RESERVE_SVW.pkl')
    pl=gu.ipickle(meta['Paths']['Geospatial'] + '\\RSLT_PLANTING_SVW.pkl')
    fire=gu.ipickle(meta['Paths']['Geospatial'] + '\\PROT_HISTORICAL_FIRE_POLYS_SP.pkl')
    pest=gu.ipickle(meta['Paths']['Geospatial'] + '\\PEST_INFESTATION_POLY.pkl')
    cut=gu.ipickle(meta['Paths']['Geospatial'] + '\\VEG_CONSOLIDATED_CUT_BLOCKS_SP.pkl')
    lul=gu.ipickle(meta['Paths']['Geospatial'] + '\\RMP_PLAN_LEGAL_POLY_SVW.pkl')
    park=gu.ipickle(meta['Paths']['Geospatial'] + '\\TA_PARK_ECORES_PA_SVW.pkl')
    ogmal=gu.ipickle(meta['Paths']['Geospatial'] + '\\RMP_OGMA_LEGAL_CURRENT_SVW.pkl')

    return sxy,atu,op,burnsev,vri,fcinv,fcres,pl,fire,pest,cut,lul,park,ogmal