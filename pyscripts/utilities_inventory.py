
'''

FOREST INVENTORY UTILITIES

'''

import numpy as np
import pandas as pd
import geopandas as gpd
import fiona
import matplotlib.pyplot as plt
import fiona
import time
from fcgadgets.pyscripts import utilities_gis as gis
from fcgadgets.pyscripts import utilities_general as gu
from fcgadgets.pyscripts import utilities_inventory as invu

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
     ('FIA_PROJECT_ID',1,'int32')]
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
     ('FIRE_DATE',0,'<U20')]
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
    
    return LayerInfo


'''============================================================================
BUILD LOOK-UP TABLES FOR THE CODES IN INVENTORY LAYERS
============================================================================'''

def BuildForestInventoryLUTs(LayerInfo):

    #for iLyr in range(0,4):
    for iLyr in range(len(LayerInfo)):    

        # I didn't include this in the last pull of RESULTS because I dont think we use it?
        #if LayerInfo[iLyr]['Layer Name']=='RSLT_FOREST_COVER_SILV_SVW': 
        #    continue
    
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