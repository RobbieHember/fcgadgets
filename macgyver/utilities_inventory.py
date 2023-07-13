
'''
INVENTORY UTILITIES
'''

#%% Import modules

import numpy as np
import pandas as pd
import geopandas as gpd
import fiona
import time
import gc as garc
import matplotlib.pyplot as plt
import scipy.io as spio
from fcgadgets.macgyver import utilities_gis as gis
from fcgadgets.macgyver import utilities_general as gu
from fcgadgets.cbrunner import cbrun_utilities as cbu
from fcgadgets.taz import aspatial_stat_models as asm

#%% Import variables

def ImportVariables(meta,pNam):

    #--------------------------------------------------------------------------
    # Define best-available (gap-filled) inventory
    #--------------------------------------------------------------------------

    print('Creating best-available variables from inventory')
    t0=time.time()

    inv={}

    # BGC zone (no longer using VRI because there are errors and gaps)
    z=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\BEC_BIOGEOCLIMATIC_POLY\\ZONE.tif')
    inv['ID_BGCZ']=gis.UpdateGridCellsize(z,meta['Geos']['RGSF'])['Data'][ meta['Geos']['iMask'] ]

    # Region code
    inv['Region Code']=meta['LUT']['Region']['Interior']*np.ones(meta[pNam]['Project']['N Stand'],dtype='int8')
    cd=[meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE']['CDF'],
        meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE']['CWH'],
        meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE']['MH']]
    ind=np.where(np.isin(inv['ID_BGCZ'],cd)==True)
    inv['Region Code'][ind]=meta['LUT']['Region']['Coast']

    # Age
    if 'Custom Age Source' in meta[pNam]['Project']:
        z=gis.OpenGeoTiff(meta[pNam]['Project']['Custom Age Source'])
        inv['Age']=gis.UpdateGridCellsize(z,meta['Geos']['RGSF'])['Data'][ meta['Geos']['iMask'] ]
    else:
        z=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\PROJ_AGE_1.tif')
        inv['Age']=gis.UpdateGridCellsize(z,meta['Geos']['RGSF'])['Data'][ meta['Geos']['iMask'] ]

    # Site index
    if 'Custom SI Source' in meta[pNam]['Project']:
        z=gis.OpenGeoTiff(meta[pNam]['Project']['Custom SI Source'])
        inv['SI']=gis.UpdateGridCellsize(z,meta['Geos']['RGSF'])['Data'][ meta['Geos']['iMask'] ]
    else:
        # From BGC Zone table
        inv['SI']=18*np.ones(meta[pNam]['Project']['N Stand'])
        u=np.unique(inv['ID_BGCZ'])
        for iU in range(u.size):
            cd=cbu.lut_n2s(meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE'],u[iU])
            ind1=np.where(inv['ID_BGCZ']==u[iU])[0]
            ind2=np.where(meta['Param']['BE']['BGC Zone Averages']['Name']==cd)[0]
            inv['SI'][ind1]=meta['Param']['BE']['BGC Zone Averages']['SI SME'][ind2]

    # Species
    if 'Custom Species Source' in meta[pNam]['Project']:
        # Use custom species
        z=gis.OpenGeoTiff(meta[pNam]['Project']['Custom Species Source'])
        inv['Spc1_ID']=gis.UpdateGridCellsize(z,meta['Geos']['RGSF'])['Data'][ meta['Geos']['iMask'] ]
        inv['Spc2_ID']=-999*np.ones(inv['Spc1_ID'].size,dtype='int16')
        inv['Spc3_ID']=-999*np.ones(inv['Spc1_ID'].size,dtype='int16')
        inv['Spc4_ID']=-999*np.ones(inv['Spc1_ID'].size,dtype='int16')
        inv['Spc5_ID']=-999*np.ones(inv['Spc1_ID'].size,dtype='int16')
        inv['Spc1_P']=100*np.ones(inv['Spc1_ID'].size,dtype='int16')
        inv['Spc2_P']=-999*np.ones(inv['Spc1_ID'].size,dtype='int16')
        inv['Spc3_P']=-999*np.ones(inv['Spc1_ID'].size,dtype='int16')
        inv['Spc4_P']=-999*np.ones(inv['Spc1_ID'].size,dtype='int16')
        inv['Spc5_P']=-999*np.ones(inv['Spc1_ID'].size,dtype='int16')
    else:
        # Use species from VRI
        for s in range(5):
            ss=str(s+1)
            z=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\SPECIES_CD_' + ss + '.tif')
            inv['Spc' + ss + '_ID']=gis.UpdateGridCellsize(z,meta['Geos']['RGSF'])['Data'][ meta['Geos']['iMask'] ]
            z=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\SPECIES_PCT_' + ss + '.tif')
            inv['Spc' + ss + '_P']=gis.UpdateGridCellsize(z,meta['Geos']['RGSF'])['Data'][ meta['Geos']['iMask'] ]

        # Clean species
        ind=np.where( (inv['Spc1_ID']==0) )[0]
        inv['Spc1_ID'][ind]=meta['LUT']['VEG_COMP_LYR_R1_POLY']['SPECIES_CD_1']['PL']
        inv['Spc1_P'][ind]=100
        for s in range(5):
            ss=str(s+1)
            ind=np.where( (inv['Spc' + ss + '_ID']==0) )[0]
            inv['Spc' + ss + '_ID'][ind]=-999
            inv['Spc' + ss + '_P'][ind]=-999

    inv['LAND_COVER_CLASS_CD_1']=-999*np.ones(meta[pNam]['Project']['N Stand'])

    inv['Year']=-999*np.ones(meta[pNam]['Project']['N Stand'])

    z=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Disturbances\\HarvestProbability.tif')
    inv['P Harvest Weight']=gis.UpdateGridCellsize(z,meta['Geos']['RGSF'])['Data'][ meta['Geos']['iMask'] ]

    P_Harvest_Max=0.486 # (%/yr)
    inv['P Harvest Weight']=inv['P Harvest Weight']/P_Harvest_Max

    z=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\REARs.tif')
    inv['THLB Layer']=gis.UpdateGridCellsize(z,meta['Geos']['RGSF'])['Data'][ meta['Geos']['iMask'] ]

    #z=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\TreeDensityClass.tif')
    #inv['TreeDensityClass']=gis.UpdateGridCellsize(z,meta['Geos']['RGSF'])['Data'][ meta['Geos']['iMask'] ]

    #z=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VEG_CONSOLIDATED_CUT_BLOCKS_SP\\HARVEST_YEAR_MaskAll.tif')
    #inv['HarvestMask']=gis.UpdateGridCellsize(z,meta['Geos']['RGSF'])['Data'][ meta['Geos']['iMask'] ]

    print(str(np.round((time.time()-t0)/60,decimals=1)) + ' min')

    #--------------------------------------------------------------------------
    # Import disturbance and management event chronology
    #--------------------------------------------------------------------------

    print('Preparing Disturbance/management event chronology')
    t0=time.time()

    # Initiate disturbance-management event history
    dmec=[None]*meta[pNam]['Project']['N Stand']

    # Initialize a counter for the number of events recorded for each stand
    cnt_e=np.zeros(meta[pNam]['Project']['N Stand'],'int16')

    # Initialize each stand with 100 events
    N_init=100

    for iStand in range(meta[pNam]['Project']['N Stand']):
        dmec[iStand]={}
        dmec[iStand]['Year']=-999*np.ones(N_init,dtype='float')
        dmec[iStand]['Month']=-999*np.ones(N_init,dtype='int16')
        dmec[iStand]['Day']=-999*np.ones(N_init,dtype='int16')
        dmec[iStand]['ID Event Type']=-999*np.ones(N_init,dtype='int16')
        dmec[iStand]['Mortality Factor']=-999*np.ones(N_init,dtype='int16')
        dmec[iStand]['Growth Factor']=-999*np.ones(N_init,dtype='int16')
        dmec[iStand]['SILV_FUND_SOURCE_CODE']=-999*np.ones(N_init,dtype='int16')
        dmec[iStand]['Planted SPH']=-999*np.ones(N_init,dtype='int16')
        for iSpc in range(6):
            dmec[iStand]['PL_SPECIES_CD' + str(iSpc+1)]=-999*np.ones(N_init,dtype='int16')
            dmec[iStand]['PL_SPECIES_PCT' + str(iSpc+1)]=-999*np.ones(N_init,dtype='int16')
            dmec[iStand]['PL_SPECIES_GW' + str(iSpc+1)]=-999*np.ones(N_init,dtype='int16')

    if meta[pNam]['Project']['DMEC Method']=='From Events':

        #----------------------------------------------------------------------
        # Events are taken from historical record (with or w/o gap-filling of age)
        #----------------------------------------------------------------------

        # Add wildfire observations
        for iY in range(6):
            Year=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\PROT_HISTORICAL_FIRE_POLYS_SP\\FIRE_YEAR_' + str(iY+1) + '_Year.tif')['Data'][ meta['Geos']['iMask'] ]
            ind=np.where(Year>0)[0]
            for iStand in ind:
                dmec[iStand]['Year'][cnt_e[iStand]]=Year[iStand]
                dmec[iStand]['ID Event Type'][cnt_e[iStand]]=meta['LUT']['Event']['Wildfire']
                dmec[iStand]['Mortality Factor'][cnt_e[iStand]]=100
                cnt_e[iStand]=cnt_e[iStand]+1

        # Add IBM observations
        for iY in range(10):
            Year=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\PEST_INFESTATION_POLY\\PEST_SEVERITY_CODE_IBM_' + str(iY+1) + '_Year.tif')['Data'][ meta['Geos']['iMask'] ]
            ind=np.where(Year>0)[0]
            for iStand in ind:
                dmec[iStand]['Year'][cnt_e[iStand]]=Year[iStand]
                dmec[iStand]['ID Event Type'][cnt_e[iStand]]=meta['LUT']['Event']['IBM']
                dmec[iStand]['Mortality Factor'][cnt_e[iStand]]=100
                cnt_e[iStand]=cnt_e[iStand]+1

        # Add harvest observations
        for iY in range(6):
            Year=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VEG_CONSOLIDATED_CUT_BLOCKS_SP\\HARVEST_YEAR_' + str(iY+1) + '_Year.tif')['Data'][ meta['Geos']['iMask'] ]
            ind=np.where(Year>0)[0]
            for iStand in ind:
                dmec[iStand]['Year'][cnt_e[iStand]]=Year[iStand]
                dmec[iStand]['ID Event Type'][cnt_e[iStand]]=meta['LUT']['Event']['Harvest']
                dmec[iStand]['Mortality Factor'][cnt_e[iStand]]=100
                print(dmec[iStand]['Year'][cnt_e[iStand]])
                print(dmec[iStand]['Mortality Factor'][cnt_e[iStand]])
                cnt_e[iStand]=cnt_e[iStand]+1

        # Add planting observations
        for iY in range(6):
            Year=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_FirstOnly_' + str(iY+1) + '_Year.tif')['Data'][ meta['Geos']['iMask'] ]
            FSC=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_FirstOnly_' + str(iY+1) + '_SILV_FUND_SOURCE_CODE.tif')['Data'][ meta['Geos']['iMask'] ]
            SPH_Planted=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_FirstOnly_' + str(iY+1) + '_SPH_Planted.tif')['Data'][ meta['Geos']['iMask'] ]
            cd={}; pct={}; gw={}
            for iSpc in range(6):
                cd[iSpc+1]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_FirstOnly_' + str(iY+1) + '_PL_SPECIES_CD' + str(iSpc+1) + '.tif')['Data'][ meta['Geos']['iMask'] ]
                pct[iSpc+1]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_FirstOnly_' + str(iY+1) + '_PL_SPECIES_PCT' + str(iSpc+1) + '.tif')['Data'][ meta['Geos']['iMask'] ]
                gw[iSpc+1]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_FirstOnly_' + str(iY+1) + '_PL_SPECIES_GW' + str(iSpc+1) + '.tif')['Data'][ meta['Geos']['iMask'] ]

            ind=np.where(Year>0)[0]
            for iStand in ind:
                dmec[iStand]['Year'][cnt_e[iStand]]=Year[iStand]
                dmec[iStand]['ID Event Type'][cnt_e[iStand]]=meta['LUT']['Event']['Planting']
                dmec[iStand]['Mortality Factor'][cnt_e[iStand]]=0
                dmec[iStand]['SILV_FUND_SOURCE_CODE'][cnt_e[iStand]]=FSC[iStand]
                dmec[iStand]['Planted SPH'][cnt_e[iStand]]=SPH_Planted[iStand]
                for iS in range(6):
                    dmec[iStand]['PL_SPECIES_CD' + str(iSpc+1)][cnt_e[iStand]]=cd[iSpc+1][iStand]
                    dmec[iStand]['PL_SPECIES_PCT' + str(iSpc+1)][cnt_e[iStand]]=pct[iSpc+1][iStand]
                    dmec[iStand]['PL_SPECIES_GW' + str(iSpc+1)][cnt_e[iStand]]=gw[iSpc+1][iStand]
                cnt_e[iStand]=cnt_e[iStand]+1

        # Add aerial fertilization observations
        for iY in range(6):
            Year=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\FECA_' + str(iY+1) + '_Year.tif')['Data'][ meta['Geos']['iMask'] ]
            ind=np.where(Year>0)[0]
            for iStand in ind:
                dmec[iStand]['Year'][cnt_e[iStand]]=Year[iStand]
                dmec[iStand]['ID Event Type'][cnt_e[iStand]]=meta['LUT']['Event']['Fertilization Aerial']
                dmec[iStand]['Mortality Factor'][cnt_e[iStand]]=0
                cnt_e[iStand]=cnt_e[iStand]+1

    elif meta[pNam]['Project']['DMEC Method']=='From Inventory Age':

        #----------------------------------------------------------------------
        # Carbon stocks in the project year are being constrained by observations
        # through input of specified age and productivity.
        #----------------------------------------------------------------------

        dType=meta['LUT']['Event']['Harvest']*np.ones(inv['Age'].size)
        dType[np.where(np.random.random(inv['Age'].size)<0.1)]=meta['LUT']['Event']['Wildfire']

        for iStand in range(meta[pNam]['Project']['N Stand']):
            dmec[iStand]['Year'][cnt_e[iStand]]=meta[pNam]['Project']['Year Project']-inv['Age'][iStand]
            dmec[iStand]['ID Event Type'][cnt_e[iStand]]=dType[iStand]
            dmec[iStand]['Mortality Factor'][cnt_e[iStand]]=100
            # if 'OAF1 Modified' in inv:
            #     if inv['OAF1 Modified'][iStand]!=inv['OAF1 Default'][iStand]:
            #         dmec[iStand]['Growth Factor'][cnt_e[iStand]]=1
            #         dmec[iStand]['ID Event Type'][cnt_e[iStand]]=meta['LUT']['Event']['LUC']
            cnt_e[iStand]=cnt_e[iStand]+1

    # Remove excess size from each stand
    for iStand in range(meta[pNam]['Project']['N Stand']):
        ind=np.where(dmec[iStand]['Year']!=-999)
        for k in dmec[iStand].keys():
            dmec[iStand][k]=dmec[iStand][k][ind]

    print(str(np.round((time.time()-t0)/60,decimals=1)) + ' min')

    #--------------------------------------------------------------------------
    # Exclude duplicate events
    #--------------------------------------------------------------------------

    if meta[pNam]['Project']['Exclude duplicate events']=='On':
        print('Removing duplicate events from DMEC')
        t0=time.time()
        dmec=Exclude_Duplicate_Events(meta,pNam,dmec)
        print(str(np.round((time.time()-t0)/60,decimals=1)) + ' min')

    #--------------------------------------------------------------------------
    # Put events in order
    # *** Must be done before several other processing steps ***
    #--------------------------------------------------------------------------

    print('Putting events in order')
    t0=time.time()
    dmec=PutEventsInOrder(meta,pNam,dmec)
    print(str(np.round((time.time()-t0)/60,decimals=1)) + ' min')

    #--------------------------------------------------------------------------
    # Ensure that a stand-replacing disturbance precedes fertilization so that age
    #--------------------------------------------------------------------------

    # Assume previous disturbance must be missing. Assume fertilization occurs at age
    # 35. Assume previous disturbance was harvest.
    # *** Must be in order of calendar date first ***
    # *** This will not work if the previous harvest or fire had a severity < 100. ***
    # Only applies to cbrunner when fertilization is simulated from TIPSY

    if meta[pNam]['Project']['Ensure aerial fert is preceded by disturbance']=='On':
        print('Ensure stand-replacing disturbance precedes fertilization')
        t0=time.time()
        dmec=Ensure_Fert_Preceded_By_Disturbance(meta,pNam,dmec,inv)
        dmec=PutEventsInOrder(meta,pNam,dmec)
        print(str(np.round((time.time()-t0)/60,decimals=1)) + ' min')

    #--------------------------------------------------------------------------
    # Ensure knockdown events are followed by slashpile events
    # *** Not written. ***
    #--------------------------------------------------------------------------

    if meta[pNam]['Project']['Revise some knockdown to include slashpile burn']=='On':
        print('Add slashpile burns after knockdown')
        #dmec=Revise_Knockdown_to_Include_Slashpile_Burning(meta,dmec)

    #--------------------------------------------------------------------------
    # Ensure every stand has a modern disturbance
    # So that there is at least one event
    #--------------------------------------------------------------------------

    if meta[pNam]['Project']['Ensure every stand has a modern disturbance']=='On':
        name_dist='Wildfire'
        severity=100
        print('Ensure every stand has a disturbance in the modern era')
        t0=time.time()
        dmec=Ensure_Every_Stand_Has_Modern_Disturbance(meta,pNam,dmec,name_dist,severity)
        print(str(np.round((time.time()-t0)/60,decimals=1)) + ' min')

    #--------------------------------------------------------------------------
    # Gap-fill stands with no event history based on VRI age
    # So that there is at least one event
    #--------------------------------------------------------------------------

    if meta[pNam]['Project']['Gap-fill DMEC with inventory age']=='On':
        print('Gap-fill DMEC with inventory age')
        t0=time.time()
        dmec=GapFill_DMEC_WithAge(meta,pNam,dmec,inv)
        print(str(np.round((time.time()-t0)/60,decimals=1)) + ' min')

    #--------------------------------------------------------------------------
    # IDW - Western Spruce Budworm - Fix severity
    #--------------------------------------------------------------------------
    # The dmec was populated with the numeric severity ID. Mortality only occurs
    # following repeated outrbreak years.

    # if meta[pNam]['Project']['Fix severity of western spruce budworm']=='On':
    #     print('Adjust severity of IDW')
    #     t0=time.time()
    #     dmec=IDW_Fix_Severity(meta,dmec)
    #     print(str(np.round((time.time()-t0)/60,decimals=1)) + ' min')

    #--------------------------------------------------------------------------
    # Reduce number of growth curves by adjusting site index
    #--------------------------------------------------------------------------

    if meta[pNam]['Project']['Revise SI to reduce num of growth curves']=='On':
        print('Reducing the number of growth curves by lowering the precision of site index')
        t0=time.time()
        #ba=ReduceVariationInSiteIndex(meta,ba)
        print(str(np.round((time.time()-t0)/60,decimals=1)) + ' min')

    #--------------------------------------------------------------------------
    # Clean species composition - TIPSY will not recognize certain codes
    #--------------------------------------------------------------------------

    # print('Cleaning species composition')
    # t0=time.time()
    # #meta,dmec,vri,fcinv=Clean_Species_Composition(meta,dmec,vri,fcinv)
    # print(str(np.round((time.time()-t0)/60,decimals=1)) + ' min')

    #--------------------------------------------------------------------------
    # Remove slashpile burns from areas where they don't often burn
    #--------------------------------------------------------------------------

    # if meta[pNam]['Project']['Remove slashpile burns from select zones']=='On':
    #     print('Removing slashpile burns from select BGC zones')
    #     t0=time.time()
    #     dmec=Remove_SlashpileBurns_From_Select_Zones(meta,dmec,ba)
    #     print(str(np.round((time.time()-t0)/60,decimals=1)) + ' min')

    #--------------------------------------------------------------------------
    # Define project type (salvage, knockdown, underplanting, NSR backlog)
    #--------------------------------------------------------------------------

    # if meta[pNam]['Project']['Special Attribution Method']=='NOSE':
    #     print('Defining types of stand establishment')
    #     t0=time.time()
    #     dmec=DefineTypeOfStandEstablishment(meta,dmec)
    #     dmec=PutEventsInOrder(meta,pNam,dmec)
    #     print(str(np.round((time.time()-t0)/60,decimals=1)) + ' min')

    #--------------------------------------------------------------------------
    # Expand DMEC to each scenario
    #--------------------------------------------------------------------------

    print('Expanding DMEC for each scenario')
    t0=time.time()
    dmec=ExpandDMEC(meta,pNam,dmec)
    print(str(np.round((time.time()-t0)/60,decimals=1)) + ' min')

    #--------------------------------------------------------------------------
    # Import land surface classification
    #--------------------------------------------------------------------------

    print('Importing land surface classification')
    t0=time.time()
    if meta[pNam]['Project']['Land Surface Class Dependent']=='Yes':
        lsc=gu.ipickle(meta['Paths'][pNam]['Data'] + '\\Inputs\\lsc.pkl')
    else:
        lsc={}
    print(str(np.round((time.time()-t0)/60,decimals=1)) + ' min')

    print('Adding land surface change to DMEC')
    t0=time.time()
    dmec=AddLandSurfaceChangesToDMEC(meta,pNam,dmec,lsc)
    print(str(np.round((time.time()-t0)/60,decimals=1)) + ' min')

    print('Putting events in order')
    t0=time.time()
    dmec=PutEventsInOrder(meta,pNam,dmec)
    print(str(np.round((time.time()-t0)/60,decimals=1)) + ' min')

    return meta,dmec,inv,lsc

#%% Exclude duplicate events from DMEC

def Exclude_Duplicate_Events(meta,pNam,dmec):
    for iStand in range(meta[pNam]['Project']['N Stand']):
        if dmec[iStand]==None:
            continue
        for key in meta['LUT']['Event'].keys():
            indType=np.where( (dmec[iStand]['ID Event Type']==meta['LUT']['Event'][key]) )[0]
            if indType.size==0:
                continue
            uYear=np.unique(np.floor(dmec[iStand]['Year'][indType]))
            for iYear in range(uYear.size):
                ind=np.where( (dmec[iStand]['ID Event Type']==meta['LUT']['Event'][key]) & (np.floor(dmec[iStand]['Year'])==uYear[iYear]) )[0]
                dmec[iStand]['ID Event Type'][ind[1:]]=-999
    return dmec

#%% Put DMEC events in order

def PutEventsInOrder(meta,pNam,dmec):

    def count_lists(l):
        return sum(1 + count_lists(i) for i in l if isinstance(i,list))

    if count_lists(dmec)==0:

        for iStand in range(meta[pNam]['Project']['N Stand']):
            #iStand=4584
            d=dmec[iStand].copy()

            ord=np.argsort(d['Year'])

            # Fix index to inciting events
            if 'IndIncitingEvent' in d.keys():
                indOld=np.where(d['IndIncitingEvent']>=0)[0]
                indNew=np.where(ord==d['IndIncitingEvent'][indOld])[0]
                d['IndIncitingEvent'][indOld]=indNew

            for key in d.keys():
                d[key]=d[key][ord]

            dmec[iStand]=d.copy()

    else:

        for iScn in range(meta[pNam]['Project']['N Scenario']):

            for iStand in range(meta[pNam]['Project']['N Stand']):

                d=dmec[iScn][iStand].copy()

                ord=np.argsort(d['Year'])

                # Fix index to inciting events
                if 'IndIncitingEvent' in d.keys():
                    indOld=np.where(d['IndIncitingEvent']>=0)[0]
                    indNew=np.where(ord==d['IndIncitingEvent'][indOld])[0]
                    d['IndIncitingEvent'][indOld]=indNew

                for key in d.keys():
                    d[key]=d[key][ord]

                dmec[iScn][iStand]=d.copy()

    return dmec

#%% Ensure disturbance prededes fertilization so age is specified
# So that age at fert is specified.

def Ensure_Fert_Preceded_By_Disturbance(meta,pNam,dmec,ba):

    ListOfTestedDist=[meta['LUT']['Event']['Wildfire'],
                      meta['LUT']['Event']['Harvest'],
                      meta['LUT']['Event']['Knockdown'],
                      meta['LUT']['Event']['Harvest Salvage'],
                      meta['LUT']['Event']['Beetles'],
                      meta['LUT']['Event']['IBM'],
                      meta['LUT']['Event']['IBB'],
                      meta['LUT']['Event']['IBD'],
                      meta['LUT']['Event']['IBS']]

    for iStand in range(meta[pNam]['Project']['N Stand']):

        if dmec[iStand]==None:
            continue

        iFert=np.where( (dmec[iStand]['ID Event Type']==meta['LUT']['Event']['Fertilization Aerial']) )[0]

        if iFert.size==0:
            continue

        iFert=iFert[0]

        # Index to events prior to first fertilization with 100% mortality
        ind=np.where( (dmec[iStand]['Year']<dmec[iStand]['Year'][iFert]) & (dmec[iStand]['Mortality Factor']==100) & np.isin(dmec[iStand]['ID Event Type'],ListOfTestedDist) )[0]

        if ind.size>0:
            continue

        print('Adding a harvest before nutrient application')

        # Assume mean of 38 + random variation (planting is 2 years after harvest)
        r=38+np.random.randint(-6,high=6)
        Year=dmec[iStand]['Year'][iFert]-r

        # Add harvest
        for k in dmec[iStand].keys():
            dmec[iStand][k]=np.append(dmec[iStand][k],-999)
        dmec[iStand]['Year'][-1]=Year
        dmec[iStand]['ID Event Type'][-1]=meta['LUT']['Event']['Harvest']
        dmec[iStand]['Mortality Factor'][-1]=100
        dmec[iStand]['Growth Factor'][-1]=-999

        # Add slashpile burn
        for k in dmec[iStand].keys():
            dmec[iStand][k]=np.append(dmec[iStand][k],-999)
        dmec[iStand]['Year'][-1]=Year+1
        dmec[iStand]['ID Event Type'][-1]=meta['LUT']['Event']['Slashpile Burn']
        dmec[iStand]['Mortality Factor'][-1]=100
        dmec[iStand]['Growth Factor'][-1]=-999

        # Add planting
        for k in dmec[iStand].keys():
            dmec[iStand][k]=np.append(dmec[iStand][k],-999)
        dmec[iStand]['Year'][-1]=Year+2
        dmec[iStand]['ID Event Type'][-1]=meta['LUT']['Event']['Planting']
        dmec[iStand]['Mortality Factor'][-1]=0
        dmec[iStand]['Growth Factor'][-1]=-999

    return dmec

#%% Expand DMEC to each scenario

def ExpandDMEC(meta,pNam,dmec_in):
    dmec_out=[None]*meta[pNam]['Project']['N Scenario']
    for iScn in range(meta[pNam]['Project']['N Scenario']):
        dmec_out[iScn]=dmec_in.copy()
    return dmec_out

#%% Reduce number of growth curves by adjusting site index

def ReduceVariationInSiteIndex(meta,pNam,inv):
    trig=0
    for i in range(1,55,2):
        ind=np.where(inv['SI']==i)[0]
        if trig==0:
            inv['SI'][ind]=inv['SI'][ind]+1
            trig=1
        else:
            inv['SI'][ind]=inv['SI'][ind]-1
            trig=0
    return inv

#%% Ensure every stand has a modern disturbance event

def Ensure_Every_Stand_Has_Modern_Disturbance(meta,pNam,dmec,name_dist,severity):

    for iStand in range(meta[pNam]['Project']['N Stand']):

        if dmec[iStand]==None:
            continue

        if dmec[iStand]['Year'].size==0:
            #print(iStand)
            #break
            r=np.random.randint(1700,2000)
            dmec[iStand]['Year']=np.append(dmec[iStand]['Year'],r)
            dmec[iStand]['ID Event Type']=np.append(dmec[iStand]['ID Event Type'],meta['LUT']['Event'][name_dist])
            dmec[iStand]['Mortality Factor']=np.append(dmec[iStand]['Mortality Factor'],np.array(severity,dtype='int16'))
            dmec[iStand]['Growth Factor']=np.append(dmec[iStand]['Growth Factor'],np.array(0,dtype='int16'))
            if 'FCI Funded' in dmec[iStand]:
                dmec[iStand]['FCI Funded']=np.append(dmec[iStand]['FCI Funded'],np.array(0,dtype='int16'))
            for v in meta['Core']['StringsToFill']:
                dmec[iStand][v]=np.append(dmec[iStand][v],-999)

    return dmec

#%% Generate DMEC from estimate of stand age

def GapFill_DMEC_WithAge(meta,pNam,dmec,inv):

    for iStand in range(meta[pNam]['Project']['N Stand']):

        if dmec[iStand]==None:
            continue

        if dmec[iStand]['Mortality Factor'].size>0:
            if np.max(dmec[iStand]['Mortality Factor']==100):
                # If there is a stand-replacing disturbance on record, no need to proceed
                continue

        if inv['Age'][iStand]>=0:

            r=np.random.random(1)
            if r<0.33:
                type=meta['LUT']['Event']['Wildfire']
            else:
                type=meta['LUT']['Event']['IBM']

            for k in dmec[iStand].keys():
                dmec[iStand][k]=np.append(dmec[iStand][k],-999)

            dmec[iStand]['Year'][-1]=meta[pNam]['Project']['Year Project']-inv['Age'][iStand]
            dmec[iStand]['ID Event Type'][-1]=type
            dmec[iStand]['Mortality Factor'][-1]=np.array(100,dtype='int16')
            dmec[iStand]['Growth Factor'][-1]=-999

    return dmec

#%% Remove slashpile burns in select BGC zones

def Remove_SlashpileBurns_From_Select_Zones(meta,pNam,dmec,inv):

    for iStand in range(meta[pNam]['Project']['N Stand']):

        if dmec[iStand]==None:
            continue

        if (inv['ID_BGCZ'][iStand]==meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE']['CWH']) | (inv['ID_BGCZ'][iStand]==meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE']['ICH']):
            ind=np.where(dmec[iStand]['ID Event Type']!=meta['LUT']['Event']['Slashpile Burn'])[0]
            if ind.size>0:
                for key in dmec[iStand]:
                    if key=='Scenario Affected':
                        continue
                    dmec[iStand][key]=dmec[iStand][key][ind]

    return dmec

#%% Define type of stnad establishment (salvage, knockdown, underplanting, NSR backlog)

# Create index to the inciting disturbance:
# It is critical that the above steps (ensuring there is a previous
# stand-replacing disturbance) work before this will work properly.

# If it is a salvage logging project, change the inciting disturbance (harvest)
# to be "Harvest Salvage" so that it removes a higher percentage of snags. I
# checked that a combo of "Harvest Salvage" + "Slashpile Burn" is consistent with
# the custom harvest (with slashpile burn) used in the salvage demo (March 2021).

# The mortality corrections can be overridden by the adjustment of species-specific
# mortality (Adjust species-specific mortality='Off' to avoid this)

def DefineTypeOfStandEstablishment(meta,pNam,dmec):

    # Threshold search period for disturbances prior to stand establishment
    LeadUpPeriod=20

    # Look up table for project type
    meta['LUT']['SE Type']={'Not simulated':0,'SL':1,'KD':2,'SP':3,'NSR Backlog':4,'Unclassified':5}

    # Initialize project type
    meta[pNam]['Project']['SE Type']=np.zeros(meta[pNam]['Project']['N Stand'],dtype=int)

    for iStand in range(meta[pNam]['Project']['N Stand']):

        L=dmec[iStand]['Year'].size

        # Initialize variable
        dmec[iStand]['IndIncitingEvent']=-999*np.ones(L,dtype=int)

        # Status
        StatNO=np.zeros(L,dtype=int)
        for iA in range(L):
            if QueryNonObStandEstablishment(meta,pNam,dmec[iStand]['SILV_FUND_SOURCE_CODE'][iA])==True:
                StatNO[iA]=1

        # Index to stand establishment events
        #iEstab=np.where( (dmec[iStand]['ID Event Type']==meta['LUT']['Event']['Planting']) & (StatNO==True) | \
        #                (dmec[iStand]['ID Event Type']==meta['LUT']['Event']['Direct Seeding']) & (StatNO==True) )[0]

        # Exclude direct seeding
        iEstab=np.where( (dmec[iStand]['ID Event Type']==meta['LUT']['Event']['Planting']) & (StatNO==True) )[0]

        if iEstab.size==0:
            continue

        # If multiple planting events, focus on...
        iEstab=iEstab[0]

        # Define a lead-up period
        iLeadUp=np.where( (dmec[iStand]['Year']>=dmec[iStand]['Year'][iEstab]-LeadUpPeriod) & (dmec[iStand]['Year']<=dmec[iStand]['Year'][iEstab]) )[0]
        Year_LeadUp=dmec[iStand]['Year'][iLeadUp]
        id=dmec[iStand]['ID Event Type'][iLeadUp]
        mort=dmec[iStand]['Mortality Factor'][iLeadUp]

        # Index to inciting event types in the lead up period
        iH=np.where( (id==meta['LUT']['Event']['Harvest']) |  (id==meta['LUT']['Event']['Harvest Salvage']) )[0]
        iH_All=iH.copy()

        # There can be multiple harvests - try choosing the last one
        if iH.size>1:
            iH=iH[-1]

        # Index to knockdown
        iKD=np.where(id==meta['LUT']['Event']['Knockdown'])[0]

        # Index to wildfire
        iWF=np.where(id==meta['LUT']['Event']['Wildfire'])[0]

        # Index to insects
        iI=np.where( (id==meta['LUT']['Event']['Beetles']) | (id==meta['LUT']['Event']['IBM']) | (id==meta['LUT']['Event']['IBD']) | (id==meta['LUT']['Event']['IBB']) | (id==meta['LUT']['Event']['IBS']) )[0]

        # Index to wildfire and insects
        iWF_and_I=np.where( (id==meta['LUT']['Event']['Wildfire']) | (id==meta['LUT']['Event']['Beetles']) | (id==meta['LUT']['Event']['IBM']) | (id==meta['LUT']['Event']['IBD']) | (id==meta['LUT']['Event']['IBB']) | (id==meta['LUT']['Event']['IBS']) )[0]

        # Add a negative one so that it is never empty
        iHa=np.append(-1,iH)
        iKDa=np.append(-1,iKD)
        iWFa=np.append(-1,iWF)
        #iIa=np.append(-1,iI)

        if (iH.size>0) & (Year_LeadUp[iH]>=1987) & (np.max(iHa)>=np.max(iKDa)) & (np.max(iHa)>=np.max(iWFa)):

            #----------------------------------------------------------------------
            # Salvage logging
            #----------------------------------------------------------------------

            meta[pNam]['Project']['SE Type'][iStand]=meta['LUT']['SE Type']['SL']

            # Find the event that incited the funded stand establishment
            # If the event has more than 75% mortality, use it. If it is less than
            # 75%, assume disturbance databases have underestimated mortality and
            # change mortality to 80%.

            if (iWF.size>0) & (iI.size==0):

                # Wildfire only
                iMaxMort=np.where(mort[iWF]==np.max(mort[iWF]))[0]
                if iMaxMort.size>1:
                    iMaxMort=iMaxMort[0]

                dmec[iStand]['IndIncitingEvent'][iEstab]=iLeadUp[iWF[iMaxMort][0]]#iLeadUp[iWF[iMaxMort]]

                if mort[iMaxMort]<75:
                    # Unrealistically low, fix
                    dmec[iStand]['Mortality Factor'][iLeadUp[iWF[iMaxMort]]]=80

            elif (iWF.size==0) & (iI.size>0):

                # Insects only
                iMaxMort=np.where(mort[iI]==np.max(mort[iI]))[0]
                if iMaxMort.size>1:
                    iMaxMort=iMaxMort[0]

                dmec[iStand]['IndIncitingEvent'][iEstab]=iLeadUp[iI[iMaxMort]]

                if mort[iI[iMaxMort]]<75:
                    # Unrealistically low, fix
                    dmec[iStand]['Mortality Factor'][iLeadUp[iI[iMaxMort]]]=80

            elif (iWF_and_I.size>0):

                # Wildfire and insects
                iMaxMort=np.where(mort[iWF_and_I]==np.max(mort[iWF_and_I]))[0]
                if iMaxMort.size>1:
                    iMaxMort=iMaxMort[0]

                dmec[iStand]['IndIncitingEvent'][iEstab]=iLeadUp[iWF_and_I[iMaxMort]]

                if mort[iWF_and_I[iMaxMort]]<75:
                    # Unrealistically low, fix
                    dmec[iStand]['Mortality Factor'][iLeadUp[iWF_and_I[iMaxMort]]]=80

            elif (iWF.size==0) & (iI.size==0):

                # Disturbance databases presumably misses the inciting event (eg fireguards), create
                # an event five years before the harvest

                ID_GapFill=meta['LUT']['Event']['Wildfire']
                Year_GapFill=dmec[iStand]['Year'][iLeadUp[iH]]-5
                Mort_GapFill=85

                dmec[iStand]['Year']=np.append(dmec[iStand]['Year'],Year_GapFill)
                dmec[iStand]['ID Event Type']=np.append(dmec[iStand]['ID Event Type'],ID_GapFill)
                dmec[iStand]['Mortality Factor']=np.append(dmec[iStand]['Mortality Factor'],Mort_GapFill)
                dmec[iStand]['Growth Factor']=np.append(dmec[iStand]['Growth Factor'],-999)
                for v in meta['Core']['StringsToFill']:
                    dmec[iStand][v]=np.append(dmec[iStand][v],-999)
                dmec[iStand]['IndIncitingEvent']=np.append(dmec[iStand]['IndIncitingEvent'],-999)

                iIncite=np.where( (dmec[iStand]['ID Event Type']==ID_GapFill) & (dmec[iStand]['Year']==Year_GapFill) )[0]
                dmec[iStand]['IndIncitingEvent'][iEstab]=iIncite

            # Change harvest to salvage so that more dead trees are taken
            dmec[iStand]['ID Event Type'][iLeadUp[iH_All]]=meta['LUT']['Event']['Harvest Salvage']

        elif (iH.size>0) & (Year_LeadUp[iH]<1987) & (np.max(iHa)>=np.max(iKDa)) & (np.max(iHa)>=np.max(iWFa)):

            #----------------------------------------------------------------------
            # NSR backlog
            #----------------------------------------------------------------------

            meta[pNam]['Project']['SE Type'][iStand]=meta['LUT']['SE Type']['NSR Backlog']

            # There is no inciting event - it is regen failure or poor effort that leads
            # to the decision to plant by gov programs. Add dummy event just before the
            # non-ob stand establishment event. Assume trace beetles

            ID_GapFill=meta['LUT']['Event']['Regen Failure']
            #Year_GapFill=dmec[iStand]['Year'][iEstab]-0.1
            Year_GapFill=dmec[iStand]['Year'][iLeadUp[iH]]+0.1
            Mort_GapFill=100

            dmec[iStand]['Year']=np.append(dmec[iStand]['Year'],Year_GapFill)
            dmec[iStand]['ID Event Type']=np.append(dmec[iStand]['ID Event Type'],ID_GapFill)
            dmec[iStand]['Mortality Factor']=np.append(dmec[iStand]['Mortality Factor'],Mort_GapFill)
            dmec[iStand]['Growth Factor']=np.append(dmec[iStand]['Growth Factor'],-999)
            dmec[iStand]['IndIncitingEvent']=np.append(dmec[iStand]['IndIncitingEvent'],-999)
            for v in meta['Core']['StringsToFill']:
                dmec[iStand][v]=np.append(dmec[iStand][v],-999)

            ind=np.where( (dmec[iStand]['ID Event Type']==ID_GapFill) & (dmec[iStand]['Year']==Year_GapFill) )[0]
            dmec[iStand]['IndIncitingEvent'][iEstab]=ind

            #        # As initial plantings presumably fail, add a regen failure event.
            #
            #        ID_GapFill=meta['LUT']['Event']['Regen Failure'].copy()
            #        Year_GapFill=dmec[iStand]['Year'][iH]+1
            #        Mort_GapFill=100
            #
            #        dmec[iStand]['Year']=np.append(dmec[iStand]['Year'],Year_GapFill)
            #        dmec[iStand]['ID Event Type']=np.append(dmec[iStand]['ID Event Type'],ID_GapFill)
            #        dmec[iStand]['Mortality Factor']=np.append(dmec[iStand]['Mortality Factor'],Mort_GapFill)
            #        dmec[iStand]['Growth Factor']=np.append(dmec[iStand]['Growth Factor'],-999)
            #        for v in meta['Core']['StringsToFill']:
            #            dmec[iStand][v]=np.append(dmec[iStand][v],-999)
            #        dmec[iStand]['IndIncitingEvent']=np.append(dmec[iStand]['IndIncitingEvent'],-999)
            #
            #        # Put everything back in order
            #        dTmp=dmec[iStand].copy()
            #        ord=np.argsort(dTmp['Year'])
            #        for key in dTmp.keys():
            #            dTmp[key]=dTmp[key][ord]
            #            dmec[iStand]=dTmp.copy()

        elif (iKD.size>0) & (np.max(iKDa)>=np.max(iHa)) & (np.max(iKDa)>=np.max(iWFa)):

            #----------------------------------------------------------------------
            # Knockdown
            #----------------------------------------------------------------------

            meta[pNam]['Project']['SE Type'][iStand]=meta['LUT']['SE Type']['KD']

            # Find the event that incited the funded stand establishment
            # If the event has more than 75% mortality, use it. If it is less than
            # 75%, assume disturbance databases have underestimated mortality and
            # change mortality to 90%.
            if (iWF.size>0) & (iI.size==0):

                # Wildfire only
                iMaxMort=np.where(mort[iWF]==np.max(mort[iWF]))[0]
                if iMaxMort.size>0:
                    iMaxMort=iMaxMort[0]

                dmec[iStand]['IndIncitingEvent'][iEstab]=iLeadUp[iWF[iMaxMort]]

                if mort[iWF[iMaxMort]]<75:
                    # Unrealistically low, fix
                    dmec[iStand]['Mortality Factor'][iLeadUp[iWF[iMaxMort]]]=90

            elif (iWF.size==0) & (iI.size>0):

                # Insects only
                iMaxMort=np.where(mort[iI]==np.max(mort[iI]))[0]
                if iMaxMort.size>0:
                    iMaxMort=iMaxMort[0]

                if mort[iI[iMaxMort]]>=75:
                    dmec[iStand]['IndIncitingEvent'][iEstab]=iLeadUp[iI[iMaxMort]]
                else:
                    dmec[iStand]['Mortality Factor'][iLeadUp[iI[iMaxMort]]]=90
                    dmec[iStand]['IndIncitingEvent'][iEstab]=iLeadUp[iI[iMaxMort]]

            elif (iWF_and_I.size>0):

                # Wildfire and insects
                iMaxMort=np.where(mort[iWF_and_I]==np.max(mort[iWF_and_I]))[0]
                if iMaxMort.size>0:
                    iMaxMort=iMaxMort[0]

                dmec[iStand]['IndIncitingEvent'][iEstab]=iLeadUp[iWF_and_I[iMaxMort]]
                if mort[iWF_and_I[iMaxMort]]<75:
                    # Unrealistically low, fix
                    dmec[iStand]['Mortality Factor'][iLeadUp[iWF_and_I[iMaxMort]]]=90

        elif (iWF.size>0) & (np.max(iWFa)>=np.max(iHa)) & (np.max(iWFa)>=np.max(iKDa)):

            #----------------------------------------------------------------------
            # Straight planting
            #----------------------------------------------------------------------

            meta[pNam]['Project']['SE Type'][iStand]=meta['LUT']['SE Type']['SP']

            # Find the event that incited the funded stand establishment
            # If the event has more than 75% mortality, use it. If it is less than
            # 90%, assume disturbance databases have underestimated mortality and
            # change mortality to 100%.

            iMaxMort=np.where(mort[iWF]==np.max(mort[iWF]))[0]
            if iMaxMort.size>1:
                iMaxMort=iMaxMort[0]

            dmec[iStand]['IndIncitingEvent'][iEstab]=iLeadUp[iWF[iMaxMort]]
            if mort[iWF[iMaxMort]]<90:
                # Unrealistically low, fix
                dmec[iStand]['Mortality Factor'][iLeadUp[iWF[iMaxMort]]]=95

        else:
            pass
#            #----------------------------------------------------------------------
#            # Unclassified -> straight planting or salvage
#            # No inciting event was found - assume a mix of salvage and straight planting
#            #----------------------------------------------------------------------
#
#            rn=np.random.random(1)[0]
#
#            if rn<0.5:
#
#                # Straight planting
#                meta[pNam]['Project']['SE Type'][iStand]=meta['LUT']['SE Type']['SP']
#                ID_GapFill=meta['LUT']['Event']['Wildfire']
#                Year_GapFill=dmec[iStand]['Year'][iEstab]-5
#                Mort_GapFill=100
#                dmec[iStand]['Year']=np.append(dmec[iStand]['Year'],Year_GapFill)
#                dmec[iStand]['ID Event Type']=np.append(dmec[iStand]['ID Event Type'],ID_GapFill)
#                dmec[iStand]['Mortality Factor']=np.append(dmec[iStand]['Mortality Factor'],Mort_GapFill)
#                dmec[iStand]['Growth Factor']=np.append(dmec[iStand]['Growth Factor'],-999)
#                dmec[iStand]['IndIncitingEvent']=np.append(dmec[iStand]['IndIncitingEvent'],-999)
#                for v in meta['Core']['StringsToFill']:
#                    dmec[iStand][v]=np.append(dmec[iStand][v],-999)
#                ind=np.where( (dmec[iStand]['ID Event Type']==ID_GapFill) & (dmec[iStand]['Year']==Year_GapFill) )[0]
#                dmec[iStand]['IndIncitingEvent'][iEstab]=ind
#
#            else:
#
#                # Salvage logging
#                meta[pNam]['Project']['SE Type'][iStand]=meta['LUT']['SE Type']['SL']
#                ID_GapFill=meta['LUT']['Event']['IBM']
#                Year_GapFill=dmec[iStand]['Year'][iEstab]-7
#                Mort_GapFill=75
#                dmec[iStand]['Year']=np.append(dmec[iStand]['Year'],Year_GapFill)
#                dmec[iStand]['ID Event Type']=np.append(dmec[iStand]['ID Event Type'],ID_GapFill)
#                dmec[iStand]['Mortality Factor']=np.append(dmec[iStand]['Mortality Factor'],Mort_GapFill)
#                dmec[iStand]['Growth Factor']=np.append(dmec[iStand]['Growth Factor'],-999)
#                dmec[iStand]['IndIncitingEvent']=np.append(dmec[iStand]['IndIncitingEvent'],-999)
#                for v in meta['Core']['StringsToFill']:
#                    dmec[iStand][v]=np.append(dmec[iStand][v],-999)
#                ind=np.where( (dmec[iStand]['ID Event Type']==ID_GapFill) & (dmec[iStand]['Year']==Year_GapFill) )[0]
#                dmec[iStand]['IndIncitingEvent'][iEstab]=ind
#
#                # Now add harvest
#                ID_GapFill=meta['LUT']['Event']['Harvest Salvage']
#                Year_GapFill=dmec[iStand]['Year'][iEstab]-2
#                Mort_GapFill=100
#                dmec[iStand]['Year']=np.append(dmec[iStand]['Year'],Year_GapFill)
#                dmec[iStand]['ID Event Type']=np.append(dmec[iStand]['ID Event Type'],ID_GapFill)
#                dmec[iStand]['Mortality Factor']=np.append(dmec[iStand]['Mortality Factor'],Mort_GapFill)
#                dmec[iStand]['Growth Factor']=np.append(dmec[iStand]['Growth Factor'],-999)
#                dmec[iStand]['IndIncitingEvent']=np.append(dmec[iStand]['IndIncitingEvent'],-999)
#                for v in meta['Core']['StringsToFill']:
#                    dmec[iStand][v]=np.append(dmec[iStand][v],-999)

    #--------------------------------------------------------------------------
    # Look at summary by type
    #--------------------------------------------------------------------------
    flg=0
    if flg==1:

        N=np.zeros(len(meta['LUT']['SE Type']))
        for i in range(len(meta['LUT']['SE Type'])):
            N[i]=np.where(meta[pNam]['Project']['SE Type']==i)[0].size
        plt.close('all')
        fig,ax=plt.subplots(1,figsize=gu.cm2inch(16,7))
        ax.bar(np.arange(0,N.size),N/np.sum(N)*100)
        ax.set(xticks=np.arange(0,N.size),xticklabels=np.array(['None','SL','KD','SP','NSR Backlog','Unclassified']))

    return dmec

#%% Add changes in land surface classfication to DMEC

def AddLandSurfaceChangesToDMEC(meta,pNam,dmec,lsc):

    for iScn in range(meta[pNam]['Project']['N Scenario']):

        # Name of LS scenario for each scenario
        nam=meta[pNam]['Scenario'][iScn]['Land Surface Scenario']

        # Only continue if using events from change in land surface class
        if nam=='None':
            continue

        # Index to LSC scenario
        for i in range(len(lsc['Scenarios'])):
            if lsc['Scenarios'][i]['Name']==nam:
                idx=i
                break

        if nam!='None':

            #Cover=np.reshape(lsc['Scenarios'][idx]['Cover'].copy(),(lsc['tv'].size,meta[pNam]['Project']['N Stand']))
            Use=np.reshape(lsc['Scenarios'][idx]['Use'].copy(),(lsc['tv'].size,meta[pNam]['Project']['N Stand']))

            #----------------------------------------------------------------------
            # Fuel breaks
            #----------------------------------------------------------------------

            nam='Fuel Break'
            indS=np.unique(np.where( Use==meta['LUT']['LSC']['Use'][nam] )[1])
            if indS.size>0:
                for i in range(indS.size):

                    iS=indS[i]
                    iT=np.where(Use[:,iS]==meta['LUT']['LSC']['Use'][nam])[0][0]

                    # Add harvest
                    dmec[iScn][iS]['Year']=np.append(dmec[iScn][iS]['Year'],lsc['tv'][iT])
                    dmec[iScn][iS]['ID Event Type']=np.append(dmec[iScn][iS]['ID Event Type'],meta['LUT']['Event']['Harvest'])
                    dmec[iScn][iS]['Mortality Factor']=np.append(dmec[iScn][iS]['Mortality Factor'],100)
                    dmec[iScn][iS]['Growth Factor']=np.append(dmec[iScn][iS]['Growth Factor'],-999)
                    if 'IndIncitingEvent' in dmec[iScn][iS]:
                        dmec[iScn][iS]['IndIncitingEvent']=np.append(dmec[iScn][iS]['IndIncitingEvent'],-999)
                    for v in meta['Core']['StringsToFill']:
                        dmec[iScn][iS][v]=np.append(dmec[iScn][iS][v],-999)

                    # Add slashpile burn
                    dmec[iScn][iS]['Year']=np.append(dmec[iScn][iS]['Year'],lsc['tv'][iT]+1)
                    dmec[iScn][iS]['ID Event Type']=np.append(dmec[iScn][iS]['ID Event Type'],meta['LUT']['Event']['Slashpile Burn'])
                    dmec[iScn][iS]['Mortality Factor']=np.append(dmec[iScn][iS]['Mortality Factor'],100)
                    dmec[iScn][iS]['Growth Factor']=np.append(dmec[iScn][iS]['Growth Factor'],-999)
                    if 'IndIncitingEvent' in dmec[iScn][iS]:
                        dmec[iScn][iS]['IndIncitingEvent']=np.append(dmec[iScn][iS]['IndIncitingEvent'],-999)
                    for v in meta['Core']['StringsToFill']:
                        dmec[iScn][iS][v]=np.append(dmec[iScn][iS][v],-999)

                    # Add planting
                    dmec[iScn][iS]['Year']=np.append(dmec[iScn][iS]['Year'],lsc['tv'][iT]+2)
                    dmec[iScn][iS]['ID Event Type']=np.append(dmec[iScn][iS]['ID Event Type'],meta['LUT']['Event']['Planting'])
                    dmec[iScn][iS]['Mortality Factor']=np.append(dmec[iScn][iS]['Mortality Factor'],0)
                    dmec[iScn][iS]['Growth Factor']=np.append(dmec[iScn][iS]['Growth Factor'],-999)
                    if 'IndIncitingEvent' in dmec[iScn][iS]:
                        dmec[iScn][iS]['IndIncitingEvent']=np.append(dmec[iScn][iS]['IndIncitingEvent'],-999)
                    for v in meta['Core']['StringsToFill']:
                        dmec[iScn][iS][v]=np.append(dmec[iScn][iS][v],-999)
                    dmec[iScn][iS]['PL_SPECIES_CD1'][-1]=meta['LUT']['VEG_COMP_LYR_R1_POLY']['Spc1_ID']['AT']
                    dmec[iScn][iS]['PL_SPECIES_PCT1'][-1]=100

                    # Add harvest
                    RotationLength=14
                    for iR in range(1,10):
                        dmec[iScn][iS]['Year']=np.append(dmec[iScn][iS]['Year'],lsc['tv'][iT]+2+iR*RotationLength)
                        dmec[iScn][iS]['ID Event Type']=np.append(dmec[iScn][iS]['ID Event Type'],meta['LUT']['Event']['Harvest'])
                        dmec[iScn][iS]['Mortality Factor']=np.append(dmec[iScn][iS]['Mortality Factor'],100)
                        dmec[iScn][iS]['Growth Factor']=np.append(dmec[iScn][iS]['Growth Factor'],-999)
                        if 'IndIncitingEvent' in dmec[iScn][iS]:
                            dmec[iScn][iS]['IndIncitingEvent']=np.append(dmec[iScn][iS]['IndIncitingEvent'],-999)
                        for v in meta['Core']['StringsToFill']:
                            dmec[iScn][iS][v]=np.append(dmec[iScn][iS][v],-999)

            #----------------------------------------------------------------------
            # Energy Production
            #----------------------------------------------------------------------

            nam='Energy Production'
            indS=np.unique(np.where( Use==meta['LUT']['LSC']['Use'][nam] )[1])
            if indS.size>0:
                for i in range(indS.size):

                    iS=indS[i]
                    iT=np.where(Use[:,iS]==meta['LUT']['LSC']['Use'][nam])[0][0]

                    # Add harvest
                    dmec[iScn][iS]['Year']=np.append(dmec[iScn][iS]['Year'],lsc['tv'][iT])
                    dmec[iScn][iS]['ID Event Type']=np.append(dmec[iScn][iS]['ID Event Type'],meta['LUT']['Event']['Harvest'])
                    dmec[iScn][iS]['Mortality Factor']=np.append(dmec[iScn][iS]['Mortality Factor'],100)
                    dmec[iScn][iS]['Growth Factor']=np.append(dmec[iScn][iS]['Growth Factor'],-999)
                    if 'IndIncitingEvent' in dmec[iScn][iS]:
                        dmec[iScn][iS]['IndIncitingEvent']=np.append(dmec[iScn][iS]['IndIncitingEvent'],-999)
                    for v in meta['Core']['StringsToFill']:
                        dmec[iScn][iS][v]=np.append(dmec[iScn][iS][v],-999)

                    # Add slashpile burn
                    dmec[iScn][iS]['Year']=np.append(dmec[iScn][iS]['Year'],lsc['tv'][iT]+1)
                    dmec[iScn][iS]['ID Event Type']=np.append(dmec[iScn][iS]['ID Event Type'],meta['LUT']['Event']['Slashpile Burn'])
                    dmec[iScn][iS]['Mortality Factor']=np.append(dmec[iScn][iS]['Mortality Factor'],100)
                    dmec[iScn][iS]['Growth Factor']=np.append(dmec[iScn][iS]['Growth Factor'],-999)
                    if 'IndIncitingEvent' in dmec[iScn][iS]:
                        dmec[iScn][iS]['IndIncitingEvent']=np.append(dmec[iScn][iS]['IndIncitingEvent'],-999)
                    for v in meta['Core']['StringsToFill']:
                        dmec[iScn][iS][v]=np.append(dmec[iScn][iS][v],-999)

                    # Add planting
                    dmec[iScn][iS]['Year']=np.append(dmec[iScn][iS]['Year'],lsc['tv'][iT]+2)
                    dmec[iScn][iS]['ID Event Type']=np.append(dmec[iScn][iS]['ID Event Type'],meta['LUT']['Event']['Planting'])
                    dmec[iScn][iS]['Mortality Factor']=np.append(dmec[iScn][iS]['Mortality Factor'],0)
                    dmec[iScn][iS]['Growth Factor']=np.append(dmec[iScn][iS]['Growth Factor'],-999)
                    if 'IndIncitingEvent' in dmec[iScn][iS]:
                        dmec[iScn][iS]['IndIncitingEvent']=np.append(dmec[iScn][iS]['IndIncitingEvent'],-999)
                    for v in meta['Core']['StringsToFill']:
                        dmec[iScn][iS][v]=np.append(dmec[iScn][iS][v],-999)
                    dmec[iScn][iS]['PL_SPECIES_CD1'][-1]=meta['LUT']['VEG_COMP_LYR_R1_POLY']['Spc1_ID']['AT']
                    dmec[iScn][iS]['PL_SPECIES_PCT1'][-1]=100

                    # Add harvest
                    RotationLength=14
                    for iR in range(1,10):
                        dmec[iScn][iS]['Year']=np.append(dmec[iScn][iS]['Year'],lsc['tv'][iT]+2+iR*RotationLength)
                        dmec[iScn][iS]['ID Event Type']=np.append(dmec[iScn][iS]['ID Event Type'],meta['LUT']['Event']['Harvest'])
                        dmec[iScn][iS]['Mortality Factor']=np.append(dmec[iScn][iS]['Mortality Factor'],100)
                        dmec[iScn][iS]['Growth Factor']=np.append(dmec[iScn][iS]['Growth Factor'],-999)
                        if 'IndIncitingEvent' in dmec[iScn][iS]:
                            dmec[iScn][iS]['IndIncitingEvent']=np.append(dmec[iScn][iS]['IndIncitingEvent'],-999)
                        for v in meta['Core']['StringsToFill']:
                            dmec[iScn][iS][v]=np.append(dmec[iScn][iS][v],-999)

    return dmec

#%% Query for non-obligation events

def QueryNonObStandEstablishment(meta,pNam,ID_FSC):

    ListOfNonObFSC=['FTL','FTM','RBM','RBL','FR','VG','FIL','FID','FIM','S', \
                    'FRP','XXX','O','GFS','IR','FES','FCE','FCM']

    if ID_FSC==0:
        String_FSC=['Unidentified']
    else:
        String_FSC=cbu.lut_n2s(meta['LUT']['RSLT_ACTIVITY_TREATMENT_SVW']['SILV_FUND_SOURCE_CODE'],ID_FSC)

    return np.isin(String_FSC,ListOfNonObFSC)

#%% Adjust species-specific mortality

# Make sure iScn_Actual is the index to the actual inventory. If you run it with
# a counterfactual scenario, it may encounter stands with only one GC. Yet, this
# method requires a previous GC.

def AdjustSpeciesSpecificMortality(meta,pNam,dmec,gc,iScn_Actual):

    # Species affected sets
    Pest_List=['IBM','IBB','IBS','IBD','IDW']

    SA_List=[None]*len(Pest_List)
    for iPest in range(len(Pest_List)):
        ind=np.where( (meta['Param']['BE']['DistBySC']['Name']==Pest_List[iPest]) )[0][0]
        SA_List[iPest]=np.array([meta['Param']['BE']['DistBySC']['SpcCD1'][ind],meta['Param']['BE']['DistBySC']['SpcCD2'][ind],
                 meta['Param']['BE']['DistBySC']['SpcCD3'][ind],meta['Param']['BE']['DistBySC']['SpcCD4'][ind],
                 meta['Param']['BE']['DistBySC']['SpcCD5'][ind],meta['Param']['BE']['DistBySC']['SpcCD6'][ind]])

    for iScn in range(meta[pNam]['Project']['N Scenario']):

        for iStand in range(meta[pNam]['Project']['N Stand']):

            for iYr in range(dmec[iScn][iStand]['Year'].size):

                for iPest in range(len(Pest_List)):

                    if dmec[iScn][iStand]['ID Event Type'][iYr]==meta['LUT']['Event'][Pest_List[iPest]]:

                        ind_GC=int(dmec[iScn_Actual][iStand]['ID_GC'][iYr]-1)

                        scd=[None]*4
                        try:
                            scd[0]=cbu.lut_n2s(meta['LUT']['VEG_COMP_LYR_R1_POLY']['Spc1_ID'],gc[iScn_Actual][iStand]['s1'][ind_GC])
                        except:
                            print(gc[iScn_Actual][iStand]['s1'])
                            print(iYr)
                            print(dmec[iScn_Actual][iStand]['ID_GC'])
                            print(ind_GC)

                        scd[1]=cbu.lut_n2s(meta['LUT']['VEG_COMP_LYR_R1_POLY']['Spc1_ID'],gc[iScn_Actual][iStand]['s2'][ind_GC])
                        scd[2]=cbu.lut_n2s(meta['LUT']['VEG_COMP_LYR_R1_POLY']['Spc1_ID'],gc[iScn_Actual][iStand]['s3'][ind_GC])
                        scd[3]=cbu.lut_n2s(meta['LUT']['VEG_COMP_LYR_R1_POLY']['Spc1_ID'],gc[iScn_Actual][iStand]['s4'][ind_GC])

                        spct=[None]*4
                        spct[0]=gc[iScn_Actual][iStand]['p1'][ind_GC]
                        spct[1]=gc[iScn_Actual][iStand]['p2'][ind_GC]
                        spct[2]=gc[iScn_Actual][iStand]['p3'][ind_GC]
                        spct[3]=gc[iScn_Actual][iStand]['p4'][ind_GC]

                        PercentAffected=0
                        for i in range(4):
                            if np.isin(scd[i],SA_List[iPest])==True:
                                PercentAffected=PercentAffected+spct[i]
                        dmec[iScn][iStand]['Mortality Factor'][iYr]=(PercentAffected/100)*dmec[iScn][iStand]['Mortality Factor'][iYr]

    return dmec

#%% Get unique growth curves

def ExtractUniqueGrowthCurves(meta,pNam,gc):

    ugc={}
    ugc['GC_Variable_List']=np.array(meta['Modules']['GYM']['GC_Variable_List'])[3:]

    # Calculate unique stand types
    ugc['Full']=np.zeros((int(4e6),len(meta['Modules']['GYM']['GC_Variable_List'])))

    cnt=0
    for iScn in range(meta[pNam]['Project']['N Scenario']):
        for iStand in range(meta[pNam]['Project']['N Stand']):
            gc0=gc[iScn][iStand]
            for iGC in range(gc0['ID_GC'].size):
                for k in range(len(meta['Modules']['GYM']['GC_Variable_List'])):
                    key=meta['Modules']['GYM']['GC_Variable_List'][k]
                    ugc['Full'][cnt,k]=gc0[key][iGC]
                cnt=cnt+1
    ugc['Full']=ugc['Full'][0:cnt,:]

    # Unique growth curves
    # The 'Inverse' variable acts as the crosswalk between the full and unique gc arrays
    ugc['Unique'],ugc['Index'],ugc['Inverse']=np.unique(ugc['Full'][:,3:],return_index=True,return_inverse=True,axis=0)

    return ugc

#%% Process project inputs 2

def ProcessProjectInputs2(meta,pNam,inv,dmec):

    #--------------------------------------------------------------------------
    # Indicate which scenario is affected by events
    # *** Make sure scenarios have been defined as "baseline" or "actual" ***
    #--------------------------------------------------------------------------

    for iScn in range(meta[pNam]['Project']['N Scenario']):

        for iStand in range(meta[pNam]['Project']['N Stand']):

            # Initialize indicator for each scenario
            dmec[iScn][iStand]['Scenario Affected']=np.zeros(dmec[iScn][iStand]['Year'].size,dtype='int16')

            if meta[pNam]['Project']['Special Attribution Method']=='None':

                # All events occur in all scenarios
                for iT in range(dmec[iScn][iStand]['Year'].size):
                    dmec[iScn][iStand]['Scenario Affected'][iT]=1

            elif meta[pNam]['Project']['Special Attribution Method']=='Actual':

                meta[pNam]['Project']['Activities To Exclude From Baseline']=np.array([meta['LUT']['Event']['Direct Seeding'],
                   meta['LUT']['Event']['Dwarf Mistletoe Control'],
                   meta['LUT']['Event']['Knockdown'],
                   meta['LUT']['Event']['Ripping'],
                   meta['LUT']['Event']['Disc Trenching'],
                   meta['LUT']['Event']['Harvest'],
                   meta['LUT']['Event']['Harvest Salvage'],
                   meta['LUT']['Event']['Thinning'],
                   meta['LUT']['Event']['Aerial Spray'],
                   meta['LUT']['Event']['Planting'],
                   meta['LUT']['Event']['Fertilization Aerial'],
                   meta['LUT']['Event']['Fertilization Hand'],
                   meta['LUT']['Event']['Fertilization Teabag'],
                   meta['LUT']['Event']['Slashpile Burn'],
                   meta['LUT']['Event']['Prescribed Burn']])

                for iT in range(dmec[iScn][iStand]['Year'].size):

                    if (np.isin(iScn,meta[pNam]['Project']['Baseline Indices'])==False):

                        # All events impact
                        dmec[iScn][iStand]['Scenario Affected'][iT]=1

                    elif (np.isin(iScn,meta[pNam]['Project']['Baseline Indices'])==True) & (np.isin(dmec[iScn][iStand]['ID Event Type'][iT],meta[pNam]['Project']['Activities To Exclude From Baseline'])==False):

                        # Events not in the list are added
                        dmec[iScn][iStand]['Scenario Affected'][iT]=1

                    else:
                        pass

            elif meta[pNam]['Project']['Special Attribution Method']=='BAU':

                meta[pNam]['Project']['Activities To Exclude From Baseline']=np.array([
                   meta['LUT']['Event']['Direct Seeding'],
                   meta['LUT']['Event']['Dwarf Mistletoe Control'],
                   meta['LUT']['Event']['Knockdown'],
                   meta['LUT']['Event']['Ripping'],
                   meta['LUT']['Event']['Disc Trenching'],
                   meta['LUT']['Event']['Harvest'],
                   meta['LUT']['Event']['Harvest Salvage'],
                   meta['LUT']['Event']['Thinning'],
                   meta['LUT']['Event']['Aerial Spray'],
                   meta['LUT']['Event']['Planting'],
                   meta['LUT']['Event']['Fertilization Aerial'],
                   meta['LUT']['Event']['Fertilization Hand'],
                   meta['LUT']['Event']['Fertilization Teabag'],
                   meta['LUT']['Event']['Slashpile Burn'],
                   meta['LUT']['Event']['Prescribed Burn']])

                meta[pNam]['Project']['Activities To Exclude From Actual']=np.array([
                   meta['LUT']['Event']['Dwarf Mistletoe Control'],
                   meta['LUT']['Event']['Ripping'],
                   meta['LUT']['Event']['Disc Trenching'],
                   meta['LUT']['Event']['Thinning'],
                   meta['LUT']['Event']['Aerial Spray'],
                   meta['LUT']['Event']['Fertilization Aerial'],
                   meta['LUT']['Event']['Fertilization Hand'],
                   meta['LUT']['Event']['Fertilization Teabag'],
                   meta['LUT']['Event']['Prescribed Burn']])

                # Obligation status for events of stand iStand
                StatusNO=np.zeros(dmec[iScn][iStand]['Year'].size,dtype=int)
                for iT in range(dmec[iScn][iStand]['Year'].size):
                    if QueryNonObStandEstablishment(meta,dmec[iScn][iStand]['SILV_FUND_SOURCE_CODE'][iT])==True:
                        StatusNO[iT]=1

                for iT in range(dmec[iScn][iStand]['Year'].size):

                    if (np.isin(iScn,meta[pNam]['Project']['Actual Indices'])==True) & (np.isin(dmec[iScn][iStand]['ID Event Type'][iT],meta[pNam]['Project']['Activities To Exclude From Actual'])==False) & (StatusNO[iT]==0):

                        # All but excluded events
                        dmec[iScn][iStand]['Scenario Affected'][iT]=1

                    elif (np.isin(iScn,meta[pNam]['Project']['Baseline Indices'])==True) & (np.isin(dmec[iScn][iStand]['ID Event Type'][iT],meta[pNam]['Project']['Activities To Exclude From Baseline'])==False):

                        # Events not in the list are added
                        dmec[iScn][iStand]['Scenario Affected'][iT]=1

                    else:
                        pass

            elif meta[pNam]['Project']['Special Attribution Method']=='NOSE':

                # Non-obligation stand establishment

                # List of activities that will be excluded from reforestation baseline
                meta[pNam]['Project']['Activities To Exclude From Baseline']=np.array([
                    meta['LUT']['Event']['Planting'],
                    meta['LUT']['Event']['Direct Seeding'],
                    meta['LUT']['Event']['Harvest'],
                    meta['LUT']['Event']['Harvest Salvage'],
                    meta['LUT']['Event']['Knockdown'],
                    meta['LUT']['Event']['Slashpile Burn'],
                    meta['LUT']['Event']['Disc Trenching'],
                    meta['LUT']['Event']['Ripping'],
                    meta['LUT']['Event']['Dwarf Mistletoe Control'] ])

                # Obligation status for events of stand iStand
                StatusNO=np.zeros(dmec[iScn][iStand]['Year'].size,dtype=int)
                for iA in range(dmec[iScn][iStand]['Year'].size):
                    if QueryNonObStandEstablishment(meta,dmec[iScn][iStand]['SILV_FUND_SOURCE_CODE'][iA])==True:
                        StatusNO[iA]=1

                # Index to stand establishment events
                iEstabNO=np.where( (dmec[iScn][iStand]['ID Event Type']==meta['LUT']['Event']['Planting']) & (StatusNO==True) | \
                                 (dmec[iScn][iStand]['ID Event Type']==meta['LUT']['Event']['Direct Seeding']) & (StatusNO==True) )[0]

                if iEstabNO.size>0:

                    # If multiple planting events, focus on the first instance
                    iEstabNO=iEstabNO[0]

                    # Index to inciting event
                    iIncitingEvent=dmec[iScn][iStand]['IndIncitingEvent'][iEstabNO]

                for iT in range(dmec[iScn][iStand]['Year'].size):

                    if (np.isin(iScn,meta[pNam]['Project']['Baseline Indices'])==True) & (iT>=iIncitingEvent) & (np.isin(dmec[iScn][iStand]['ID Event Type'][iT],meta[pNam]['Project']['Activities To Exclude From Baseline'])==True):
                        # Events only impact the project scenarios
                        dmec[iScn][iStand]['Scenario Affected'][iT]=1
                    else:
                        # Otherwise, everything occurs in all scenarios
                        dmec[iScn][iStand]['Scenario Affected'][iT]=1

            elif meta[pNam]['Project']['Special Attribution Method']=='Nutrient Management':

                for iT in range(dmec[iScn][iStand]['Year'].size):

                    if (np.isin(iScn,meta[pNam]['Project']['Actual Indices'])==True):
                        dmec[iScn][iStand]['Scenario Affected'][iT]=1
                    else:
                        if dmec[iScn][iStand]['ID Event Type'][iT]!=meta['LUT']['Event']['Fertilization Aerial']:
                            dmec[iScn][iStand]['Scenario Affected'][iT]=1

    #--------------------------------------------------------------------------
    # Parameterize growth curves
    #--------------------------------------------------------------------------

    print('Preparing growth curves')
    t0=time.time()

    # Prepare parameters for natural establishment (by BGC)
    pByBGC=gu.ReadExcel(meta['Paths']['Model']['Code'] + '\\Parameters\\Parameters_ByBGC.xlsx')
    SPH_Init_Nat=1500*np.ones(meta[pNam]['Project']['N Stand'],dtype='int16')
    RegenDelay_Nat=1500*np.ones(meta[pNam]['Project']['N Stand'],dtype='int16')
    u=np.unique(inv['ID_BGCZ'])
    for iU in range(u.size):
        ind0=np.where(inv['ID_BGCZ']==u[iU])[0]
        ind1=np.where(pByBGC['Name']==cbu.lut_n2s(meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE'],u[iU])[0])[0]
        SPH_Init_Nat[ind0]=pByBGC['Natural Initial Tree Density'][ind1[0]]
        RegenDelay_Nat[ind0]=pByBGC['Natural Regeneration Delay'][ind1[0]]

    # # Operational Adjustment Factor 1

    # flg=0
    # if flg==1:
    #     oaf1=0.9*np.ones(meta[pNam]['Project']['N Stand'])
    #     ind=np.where(inv['TreeDensityClass']==1)[0]
    #     oaf1[ind]=0.60
    #     ind=np.where(inv['TreeDensityClass']==2)[0]
    #     oaf1[ind]=0.75
    #     ind=np.where(inv['TreeDensityClass']==3)[0]
    #     oaf1[ind]=0.9

    # flg=0
    # if flg==1:
    #     oaf1=1.0*np.ones(meta[pNam]['Project']['N Stand'])
    #     ind=np.where(inv['HarvestMask']==0)[0]
    #     oaf1[ind]=0.5
    #     ind=np.where(inv['HarvestMask']>0)[0]
    #     oaf1[ind]=1.0
    #     inv['SI'][ind]=inv['SI'][ind]+2

    # Initialize list
    gc=[None]*meta[pNam]['Project']['N Scenario']

    for iScn in range(meta[pNam]['Project']['N Scenario']):

        gc[iScn]=[None]*meta[pNam]['Project']['N Stand']

        for iStand in range(meta[pNam]['Project']['N Stand']):

            # Index to growth curve
            cnt_gc=0

            # Initialize growth curve identifiers in DMEC
            dmec[iScn][iStand]['ID_GC']=1*np.ones(dmec[iScn][iStand]['Year'].size)

            # Initialize growth curve info
            gc[iScn][iStand]={}
            for key in meta['Modules']['GYM']['GC_Variable_List']:
                gc[iScn][iStand][key]=-999*np.ones(12)

            if inv['Region Code'][iStand]==meta['LUT']['Region']['Coast']:
                fiz=meta['LUT']['TIPSY']['FIZ']['C']
            else:
                fiz=meta['LUT']['TIPSY']['FIZ']['I']

            #--------------------------------------------------------------------------
            # Add pre-contact growth curve
            #--------------------------------------------------------------------------

            gc[iScn][iStand]['ID_Stand'][cnt_gc]=iStand
            gc[iScn][iStand]['ID_Scn'][cnt_gc]=iScn
            gc[iScn][iStand]['ID_GC'][cnt_gc]=1
            gc[iScn][iStand]['regeneration_method'][cnt_gc]=meta['LUT']['TIPSY']['regeneration_method']['N']
            gc[iScn][iStand]['s1'][cnt_gc]=inv['Spc1_ID'][iStand]
            gc[iScn][iStand]['p1'][cnt_gc]=inv['Spc1_P'][iStand]
            gc[iScn][iStand]['i1'][cnt_gc]=inv['SI'][iStand]
            gc[iScn][iStand]['s2'][cnt_gc]=inv['Spc2_ID'][iStand]
            gc[iScn][iStand]['p2'][cnt_gc]=inv['Spc2_P'][iStand]
            gc[iScn][iStand]['s3'][cnt_gc]=inv['Spc3_ID'][iStand]
            gc[iScn][iStand]['p3'][cnt_gc]=inv['Spc3_P'][iStand]
            gc[iScn][iStand]['s4'][cnt_gc]=inv['Spc4_ID'][iStand]
            gc[iScn][iStand]['p4'][cnt_gc]=inv['Spc4_P'][iStand]
            gc[iScn][iStand]['s5'][cnt_gc]=inv['Spc5_ID'][iStand]
            gc[iScn][iStand]['p5'][cnt_gc]=inv['Spc5_P'][iStand]
            gc[iScn][iStand]['init_density'][cnt_gc]=SPH_Init_Nat[iStand]
            gc[iScn][iStand]['regen_delay'][cnt_gc]=RegenDelay_Nat[iStand]
            gc[iScn][iStand]['oaf1'][cnt_gc]=meta['Modules']['GYM']['OAF1 Default']
            gc[iScn][iStand]['oaf2'][cnt_gc]=meta['Modules']['GYM']['OAF2 Default']
            gc[iScn][iStand]['bec_zone'][cnt_gc]=inv['ID_BGCZ'][iStand]
            gc[iScn][iStand]['FIZ'][cnt_gc]=fiz
            cnt_gc=cnt_gc+1

            #--------------------------------------------------------------------------
            # Add events from disturbance/management event history
            #--------------------------------------------------------------------------

            for iYr in range(dmec[iScn][iStand]['Year'].size):

                # Calculate planting density
                PlantingDensity=int(dmec[iScn][iStand]['Planted SPH'][iYr])
                if inv['Region Code'][iStand]==meta['LUT']['Region']['Interior']:
                    PlantingDensity=np.minimum(2400,np.maximum(1400,PlantingDensity))
                else:
                    PlantingDensity=np.minimum(2400,np.maximum(900,PlantingDensity))
                PlantingDensity=int(PlantingDensity)

                # Create a flag that indicates whether there are back-to-back planting
                # Back to back planting I think occurs in some cases because they go back
                # and add a bit. You can tell by looking at the treatment area - if the
                # second planting treatment area is tiny compared to the first, you could
                # ignore it I guess. Certainly not ideal, but I don't see a work around.
                # We are going to ignore the second planting for now.
                Flag_PlantingBackToBack=0
                if iYr>0:
                    if (dmec[iScn][iStand]['ID Event Type'][iYr]==meta['LUT']['Event']['Planting']):
                        Yr=dmec[iScn][iStand]['Year'][iYr]
                        indPrevPL=np.where( (dmec[iScn][iStand]['Year']==Yr-1) & (dmec[iScn][iStand]['ID Event Type']==meta['LUT']['Event']['Planting']) )[0]
                        if indPrevPL.size>0:
                            Flag_PlantingBackToBack=1

                # Index to previous disturbance for fertilization
                #IndPrevDistForFert=int(dmec[iStand]['IndPrevDistForFert'][iYr])

                # Non-obligation status
                StatusNO=QueryNonObStandEstablishment(meta,pNam,dmec[iScn][iStand]['SILV_FUND_SOURCE_CODE'][iYr])

                if (meta[pNam]['Project']['Special Attribution Method']=='None') | (meta[pNam]['Project']['Special Attribution Method']=='Actual') | (meta[pNam]['Project']['Special Attribution Method']=='Nutrient Management'):

                    #--------------------------------------------------------------
                    # None or Actual or Nutrient Management
                    #--------------------------------------------------------------

                    if (dmec[iScn][iStand]['Scenario Affected'][iYr]==1) & (dmec[iScn][iStand]['ID Event Type'][iYr]==meta['LUT']['Event']['Planting']):

                        dmec[iScn][iStand]['ID_GC'][iYr:]=dmec[iScn][iStand]['ID_GC'][iYr]+1

                        gc[iScn][iStand]['ID_Stand'][cnt_gc]=iStand
                        gc[iScn][iStand]['ID_Scn'][cnt_gc]=iScn
                        gc[iScn][iStand]['ID_GC'][cnt_gc]=dmec[iScn][iStand]['ID_GC'][iYr]
                        gc[iScn][iStand]['regeneration_method'][cnt_gc]=meta['LUT']['TIPSY']['regeneration_method']['P']
                        gc[iScn][iStand]['init_density'][cnt_gc]=int(PlantingDensity)
                        gc[iScn][iStand]['regen_delay'][cnt_gc]=0
                        gc[iScn][iStand]['i1'][cnt_gc]=inv['SI'][iStand]

                        if (dmec[iScn][iStand]['PL_SPECIES_CD1'][iYr]>0) & (dmec[iScn][iStand]['PL_SPECIES_PCT1'][iYr]>0):

                            # *** Adjust site index if it is energy production ***
                            if (dmec[iScn][iStand]['PL_SPECIES_CD1'][iYr]==meta['LUT']['VEG_COMP_LYR_R1_POLY']['Spc1_ID']['AT']):
                                gc[iScn][iStand]['i1'][cnt_gc]=30
                                gc[iScn][iStand]['init_density'][cnt_gc]=2000

                            # Using planting info if it exists
                            gc[iScn][iStand]['s1'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_CD1'][iYr]
                            gc[iScn][iStand]['p1'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_PCT1'][iYr]
                            gc[iScn][iStand]['gain1'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_GW1'][iYr]
                            gc[iScn][iStand]['selage1'][cnt_gc]=10

                            gc[iScn][iStand]['s2'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_CD2'][iYr]
                            gc[iScn][iStand]['p2'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_PCT2'][iYr]
                            gc[iScn][iStand]['gain2'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_GW2'][iYr]
                            gc[iScn][iStand]['selage2'][cnt_gc]=10

                            gc[iScn][iStand]['s3'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_CD3'][iYr]
                            gc[iScn][iStand]['p3'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_PCT3'][iYr]
                            gc[iScn][iStand]['gain3'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_GW3'][iYr]
                            gc[iScn][iStand]['selage3'][cnt_gc]=10

                            gc[iScn][iStand]['s4'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_CD4'][iYr]
                            gc[iScn][iStand]['p4'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_PCT4'][iYr]
                            gc[iScn][iStand]['gain4'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_GW4'][iYr]
                            gc[iScn][iStand]['selage4'][cnt_gc]=10

                            gc[iScn][iStand]['s5'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_CD5'][iYr]
                            gc[iScn][iStand]['p5'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_PCT5'][iYr]
                            gc[iScn][iStand]['gain5'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_GW5'][iYr]
                            gc[iScn][iStand]['selage5'][cnt_gc]=10

                        else:

                            # Otherwise assume best-available inventory spc. comp.
                            gc[iScn][iStand]['s1'][cnt_gc]=inv['Spc_CD1'][iStand]
                            gc[iScn][iStand]['p1'][cnt_gc]=inv['Spc_Pct1'][iStand]

                            gc[iScn][iStand]['s2'][cnt_gc]=inv['Spc_CD2'][iStand]
                            gc[iScn][iStand]['p2'][cnt_gc]=inv['Spc_Pct2'][iStand]

                            gc[iScn][iStand]['s3'][cnt_gc]=inv['Spc_CD3'][iStand]
                            gc[iScn][iStand]['p3'][cnt_gc]=inv['Spc_Pct3'][iStand]

                            gc[iScn][iStand]['s4'][cnt_gc]=inv['Spc_CD4'][iStand]
                            gc[iScn][iStand]['p4'][cnt_gc]=inv['Spc_Pct4'][iStand]

                            gc[iScn][iStand]['s5'][cnt_gc]=inv['Spc_CD5'][iStand]
                            gc[iScn][iStand]['p5'][cnt_gc]=inv['Spc_Pct5'][iStand]

                        gc[iScn][iStand]['oaf1'][cnt_gc]=meta['Modules']['GYM']['OAF1 Default']
                        gc[iScn][iStand]['oaf2'][cnt_gc]=meta['Modules']['GYM']['OAF2 Default']
                        gc[iScn][iStand]['bec_zone'][cnt_gc]=inv['ID_BGCZ'][iStand]
                        gc[iScn][iStand]['FIZ'][cnt_gc]=fiz

                        # Update counter
                        cnt_gc=cnt_gc+1

                elif (meta[pNam]['Project']['Special Attribution Method']=='BAU'):

                    #--------------------------------------------------------------
                    # Harvesting
                    # Exclude incremental improvements in silviculture, including:
                    #  - genetic gains
                    #--------------------------------------------------------------

                    if (dmec[iScn][iStand]['Scenario Affected'][iYr]==1) & (dmec[iScn][iStand]['ID Event Type'][iYr]==meta['LUT']['Event']['Planting']):

                        dmec[iScn][iStand]['ID_GC'][iYr:]=dmec[iScn][iStand]['ID_GC'][iYr]+1

                        gc[iScn][iStand]['ID_Stand'][cnt_gc]=iStand
                        gc[iScn][iStand]['ID_Scn'][cnt_gc]=iScn
                        gc[iScn][iStand]['ID_GC'][cnt_gc]=dmec[iScn][iStand]['ID_GC'][iYr]
                        gc[iScn][iStand]['regeneration_method'][cnt_gc]=meta['LUT']['TIPSY']['regeneration_method']['P']
                        gc[iScn][iStand]['init_density'][cnt_gc]=int(PlantingDensity)
                        gc[iScn][iStand]['regen_delay'][cnt_gc]=0
                        gc[iScn][iStand]['i1'][cnt_gc]=inv['SI'][iStand]

                        if (dmec[iScn][iStand]['PL_SPECIES_CD1'][iYr]>0) & (dmec[iScn][iStand]['PL_SPECIES_PCT1'][iYr]>0):

                            # Using planting info if it exists
                            gc[iScn][iStand]['s1'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_CD1'][iYr]
                            gc[iScn][iStand]['p1'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_PCT1'][iYr]

                            gc[iScn][iStand]['s2'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_CD2'][iYr]
                            gc[iScn][iStand]['p2'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_PCT2'][iYr]

                            gc[iScn][iStand]['s3'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_CD3'][iYr]
                            gc[iScn][iStand]['p3'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_PCT3'][iYr]

                            gc[iScn][iStand]['s4'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_CD4'][iYr]
                            gc[iScn][iStand]['p4'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_PCT4'][iYr]

                            gc[iScn][iStand]['s5'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_CD5'][iYr]
                            gc[iScn][iStand]['p5'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_PCT5'][iYr]

                        else:

                            # Otherwise assume best-available inventory spc. comp.
                            gc[iScn][iStand]['s1'][cnt_gc]=inv['Spc_CD1'][iStand]
                            gc[iScn][iStand]['p1'][cnt_gc]=inv['Spc_Pct1'][iStand]

                            gc[iScn][iStand]['s2'][cnt_gc]=inv['Spc_CD2'][iStand]
                            gc[iScn][iStand]['p2'][cnt_gc]=inv['Spc_Pct2'][iStand]

                            gc[iScn][iStand]['s3'][cnt_gc]=inv['Spc_CD3'][iStand]
                            gc[iScn][iStand]['p3'][cnt_gc]=inv['Spc_Pct3'][iStand]

                            gc[iScn][iStand]['s4'][cnt_gc]=inv['Spc_CD4'][iStand]
                            gc[iScn][iStand]['p4'][cnt_gc]=inv['Spc_Pct4'][iStand]

                            gc[iScn][iStand]['s5'][cnt_gc]=inv['Spc_CD5'][iStand]
                            gc[iScn][iStand]['p5'][cnt_gc]=inv['Spc_Pct5'][iStand]

                        gc[iScn][iStand]['oaf1'][cnt_gc]=meta['Modules']['GYM']['OAF1 Default']
                        gc[iScn][iStand]['oaf2'][cnt_gc]=meta['Modules']['GYM']['OAF2 Default']
                        gc[iScn][iStand]['bec_zone'][cnt_gc]=inv['ID_BGCZ'][iStand]
                        gc[iScn][iStand]['FIZ'][cnt_gc]=fiz

                        # Update counter
                        cnt_gc=cnt_gc+1

                elif (meta[pNam]['Project']['Special Attribution Method']=='NOSE'):

                    #--------------------------------------------------------------
                    # Non-obligation stand establishment
                    #--------------------------------------------------------------

                    # Index to event that incited NO stand establishment
                    iIncitingNOSE=int(dmec[iScn][iStand]['IndIncitingEvent'][iYr])

                    #----------------------------------------------------------------------
                    # Planting (non-obligation)
                    #----------------------------------------------------------------------

                    if (dmec[iScn][iStand]['Scenario Affected'][iYr]==1) & (dmec[iScn][iStand]['ID Event Type'][iYr]==meta['LUT']['Event']['Planting']) & (Flag_PlantingBackToBack==0) & (StatusNO==True):

                        if (np.isin(iScn,meta[pNam]['Project']['Baseline Indices'])==True):

                            # Baseline scenarios:

                            # Growth curve update
                            if iIncitingNOSE>=0:
                                # Baseline scenarios, transition at time of inciting evnet
                                dmec[iScn][iStand]['ID_GC'][iIncitingNOSE:]=dmec[iScn][iStand]['ID_GC'][iIncitingNOSE]+1
                            else:
                                print('Problem with NOSE project - cant find an inciting event (unitentified SE Type)!')
                                dmec[iScn][iStand]['ID_GC'][iYr:]=dmec[iScn][iStand]['ID_GC'][iYr]+1

                            gc[iScn][iStand]['regeneration_method'][cnt_gc]=meta['LUT']['TIPSY']['regeneration_method']['N']
                            if meta[pNam]['Project']['SE Type'][iStand]==3:
                                gc[iScn][iStand]['init_density'][cnt_gc]=200
                                gc[iScn][iStand]['regen_delay'][cnt_gc]=5
                            else:
                                gc[iScn][iStand]['init_density'][cnt_gc]=1800
                                gc[iScn][iStand]['regen_delay'][cnt_gc]=1

                        else:

                            # Project scenarios:

                            # Growth curve update: Project scenarios with planting at iYr
                            dmec[iScn][iStand]['ID_GC'][iYr:]=dmec[iScn][iStand]['ID_GC'][iYr]+1

                            gc[iScn][iStand]['regeneration_method'][cnt_gc]=meta['LUT']['TIPSY']['regeneration_method']['P']
                            gc[iScn][iStand]['init_density'][cnt_gc]=int(PlantingDensity)
                            gc[iScn][iStand]['regen_delay'][cnt_gc]=0

                        gc[iScn][iStand]['ID_Stand'][cnt_gc]=iStand
                        gc[iScn][iStand]['ID_Scn'][cnt_gc]=iScn
                        gc[iScn][iStand]['ID_GC'][cnt_gc]=dmec[iScn][iStand]['ID_GC'][iYr]
                        gc[iScn][iStand]['i1'][cnt_gc]=inv['SI'][iStand]
                        gc[iScn][iStand]['oaf1'][cnt_gc]=meta['Modules']['GYM']['OAF1 Default']
                        gc[iScn][iStand]['oaf2'][cnt_gc]=meta['Modules']['GYM']['OAF2 Default']
                        gc[iScn][iStand]['bec_zone'][cnt_gc]=inv['ID_BGCZ'][iStand]
                        gc[iScn][iStand]['FIZ'][cnt_gc]=fiz

                        if (np.isin(iScn,meta[pNam]['Project']['Baseline Indices'])==False) & (dmec[iScn][iStand]['PL_SPECIES_CD1'][iYr]>0) & (dmec[iScn][iStand]['PL_SPECIES_PCT1'][iYr]>0):

                            # Include genetic gain
                            gc[iScn][iStand]['s1'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_CD1'][iYr]
                            gc[iScn][iStand]['p1'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_PCT1'][iYr]
                            gc[iScn][iStand]['gain1'][cnt_gc]=int(dmec[iScn][iStand]['PL_SPECIES_GW1'][iYr])
                            gc[iScn][iStand]['selage1'][cnt_gc]=10

                            gc[iScn][iStand]['s2'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_CD2'][iYr]
                            gc[iScn][iStand]['p2'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_PCT2'][iYr]
                            gc[iScn][iStand]['gain2'][cnt_gc]=int(dmec[iScn][iStand]['PL_SPECIES_GW2'][iYr])
                            gc[iScn][iStand]['selage2'][cnt_gc]=10

                            gc[iScn][iStand]['s3'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_CD3'][iYr]
                            gc[iScn][iStand]['p3'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_PCT3'][iYr]
                            gc[iScn][iStand]['gain3'][cnt_gc]=int(dmec[iScn][iStand]['PL_SPECIES_GW3'][iYr])
                            gc[iScn][iStand]['selage3'][cnt_gc]=10

                            gc[iScn][iStand]['s4'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_CD4'][iYr]
                            gc[iScn][iStand]['p4'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_PCT4'][iYr]
                            gc[iScn][iStand]['gain4'][cnt_gc]=int(dmec[iScn][iStand]['PL_SPECIES_GW4'][iYr])
                            gc[iScn][iStand]['selage4'][cnt_gc]=10

                            gc[iScn][iStand]['s5'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_CD5'][iYr]
                            gc[iScn][iStand]['p5'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_PCT5'][iYr]
                            gc[iScn][iStand]['gain5'][cnt_gc]=int(dmec[iScn][iStand]['PL_SPECIES_GW5'][iYr])
                            gc[iScn][iStand]['selage5'][cnt_gc]=10

                        else:

                            # Baseline
                            gc[iScn][iStand]['s1'][cnt_gc]=inv['Spc_CD1'][iStand]
                            gc[iScn][iStand]['p1'][cnt_gc]=inv['Spc_Pct1'][iStand]

                            gc[iScn][iStand]['s2'][cnt_gc]=inv['Spc_CD2'][iStand]
                            gc[iScn][iStand]['p2'][cnt_gc]=inv['Spc_Pct2'][iStand]

                            gc[iScn][iStand]['s3'][cnt_gc]=inv['Spc_CD3'][iStand]
                            gc[iScn][iStand]['p3'][cnt_gc]=inv['Spc_Pct3'][iStand]

                            gc[iScn][iStand]['s4'][cnt_gc]=inv['Spc_CD4'][iStand]
                            gc[iScn][iStand]['p4'][cnt_gc]=inv['Spc_Pct4'][iStand]

                            gc[iScn][iStand]['s5'][cnt_gc]=inv['Spc_CD5'][iStand]
                            gc[iScn][iStand]['p5'][cnt_gc]=inv['Spc_Pct5'][iStand]

                        # Update counter
                        cnt_gc=cnt_gc+1

                    #----------------------------------------------------------------------
                    # Planting (obligation)
                    #----------------------------------------------------------------------

                    if (dmec[iScn][iStand]['Scenario Affected'][iYr]==1) & (dmec[iScn][iStand]['ID Event Type'][iYr]==meta['LUT']['Event']['Planting']) & (Flag_PlantingBackToBack==0) & (StatusNO==False):

                        dmec[iScn][iStand]['ID_GC'][iYr:]=dmec[iScn][iStand]['ID_GC'][iYr]+1

                        gc[iScn][iStand]['ID_Stand'][cnt_gc]=iStand
                        gc[iScn][iStand]['ID_Scn'][cnt_gc]=iScn
                        gc[iScn][iStand]['ID_GC'][cnt_gc]=dmec[iScn][iStand]['ID_GC'][iYr]
                        gc[iScn][iStand]['regeneration_method'][cnt_gc]=meta['LUT']['TIPSY']['regeneration_method']['P']
                        gc[iScn][iStand]['i1'][cnt_gc]=inv['SI'][iStand]
                        gc[iScn][iStand]['init_density'][cnt_gc]=int(PlantingDensity)
                        gc[iScn][iStand]['regen_delay'][cnt_gc]=0
                        gc[iScn][iStand]['oaf1'][cnt_gc]=meta['Modules']['GYM']['OAF1 Default']
                        gc[iScn][iStand]['oaf2'][cnt_gc]=meta['Modules']['GYM']['OAF2 Default']
                        gc[iScn][iStand]['bec_zone'][cnt_gc]=inv['ID_BGCZ'][iStand]
                        gc[iScn][iStand]['FIZ'][cnt_gc]=fiz

                        if (dmec[iScn][iStand]['PL_SPECIES_CD1'][iYr]>0) & (dmec[iScn][iStand]['PL_SPECIES_PCT1'][iYr]>0):
                            gc[iScn][iStand]['s1'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_CD1'][iYr]
                            gc[iScn][iStand]['p1'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_PCT1'][iYr]
                            gc[iScn][iStand]['gain1'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_GW1'][iYr]
                            gc[iScn][iStand]['selage1'][cnt_gc]=10

                            gc[iScn][iStand]['s2'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_CD2'][iYr]
                            gc[iScn][iStand]['p2'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_PCT2'][iYr]
                            gc[iScn][iStand]['gain2'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_GW2'][iYr]
                            gc[iScn][iStand]['selage2'][cnt_gc]=10

                            gc[iScn][iStand]['s3'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_CD3'][iYr]
                            gc[iScn][iStand]['p3'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_PCT3'][iYr]
                            gc[iScn][iStand]['gain3'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_GW3'][iYr]
                            gc[iScn][iStand]['selage3'][cnt_gc]=10

                            gc[iScn][iStand]['s4'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_CD4'][iYr]
                            gc[iScn][iStand]['p4'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_PCT4'][iYr]
                            gc[iScn][iStand]['gain4'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_GW4'][iYr]
                            gc[iScn][iStand]['selage4'][cnt_gc]=10

                            gc[iScn][iStand]['s5'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_CD5'][iYr]
                            gc[iScn][iStand]['p5'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_PCT5'][iYr]
                            gc[iScn][iStand]['gain5'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_GW5'][iYr]
                            gc[iScn][iStand]['selage5'][cnt_gc]=10
                        else:
                            gc[iScn][iStand]['s1'][cnt_gc]=inv['Spc_CD1'][iStand]
                            gc[iScn][iStand]['p1'][cnt_gc]=inv['Spc_Pct1'][iStand]

                            gc[iScn][iStand]['s2'][cnt_gc]=inv['Spc_CD2'][iStand]
                            gc[iScn][iStand]['p2'][cnt_gc]=inv['Spc_Pct2'][iStand]

                            gc[iScn][iStand]['s3'][cnt_gc]=inv['Spc_CD3'][iStand]
                            gc[iScn][iStand]['p3'][cnt_gc]=inv['Spc_Pct3'][iStand]

                            gc[iScn][iStand]['s4'][cnt_gc]=inv['Spc_CD4'][iStand]
                            gc[iScn][iStand]['p4'][cnt_gc]=inv['Spc_Pct4'][iStand]

                            gc[iScn][iStand]['s5'][cnt_gc]=inv['Spc_CD5'][iStand]
                            gc[iScn][iStand]['p5'][cnt_gc]=inv['Spc_Pct5'][iStand]

                        # Update counter
                        cnt_gc=cnt_gc+1

            # Add a planted stand if missing:
            # Some stands may have no recorded historical management, but then they
            # may be harvested on-the-fly. Add a second growth curve for the on-the-fly
            # planting.
            if cnt_gc==1:

                gc[iScn][iStand]['ID_Stand'][cnt_gc]=iStand
                gc[iScn][iStand]['ID_Scn'][cnt_gc]=iScn
                gc[iScn][iStand]['ID_GC'][cnt_gc]=np.max(gc[iScn][iStand]['ID_GC'])+1
                gc[iScn][iStand]['regeneration_method'][cnt_gc]=meta['LUT']['TIPSY']['regeneration_method']['P']
                gc[iScn][iStand]['i1'][cnt_gc]=inv['SI'][iStand]
                gc[iScn][iStand]['init_density'][cnt_gc]=int(1500)
                gc[iScn][iStand]['regen_delay'][cnt_gc]=0
                gc[iScn][iStand]['oaf1'][cnt_gc]=meta['Modules']['GYM']['OAF1 Default']
                gc[iScn][iStand]['oaf2'][cnt_gc]=meta['Modules']['GYM']['OAF2 Default']
                gc[iScn][iStand]['bec_zone'][cnt_gc]=inv['ID_BGCZ'][iStand]
                gc[iScn][iStand]['FIZ'][cnt_gc]=fiz

                # The natural species composition may not be realistic. Use regional
                # default planting composition

                gain=15

                if (inv['ID_BGCZ'][iStand]==meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE']['CWH']) | (inv['ID_BGCZ'][iStand]==meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE']['CDF']):

                    # Coastal
                    gc[iScn][iStand]['s1'][cnt_gc]=meta['LUT']['VEG_COMP_LYR_R1_POLY']['SPECIES_CD_1']['FDC']
                    gc[iScn][iStand]['p1'][cnt_gc]=70
                    gc[iScn][iStand]['gain1'][cnt_gc]=gain
                    gc[iScn][iStand]['selage1'][cnt_gc]=10
                    gc[iScn][iStand]['s2'][cnt_gc]=meta['LUT']['VEG_COMP_LYR_R1_POLY']['SPECIES_CD_1']['CW']
                    gc[iScn][iStand]['p2'][cnt_gc]=15
                    gc[iScn][iStand]['gain2'][cnt_gc]=gain
                    gc[iScn][iStand]['selage2'][cnt_gc]=10
                    gc[iScn][iStand]['s3'][cnt_gc]=meta['LUT']['VEG_COMP_LYR_R1_POLY']['SPECIES_CD_1']['BA']
                    gc[iScn][iStand]['p3'][cnt_gc]=15
                    gc[iScn][iStand]['gain3'][cnt_gc]=gain
                    gc[iScn][iStand]['selage3'][cnt_gc]=10

                elif (inv['ID_BGCZ'][iStand]==meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE']['ICH']):

                    # Interior wetbelt
                    gc[iScn][iStand]['s1'][cnt_gc]=meta['LUT']['VEG_COMP_LYR_R1_POLY']['SPECIES_CD_1']['CW']
                    gc[iScn][iStand]['p1'][cnt_gc]=70
                    gc[iScn][iStand]['gain1'][cnt_gc]=gain
                    gc[iScn][iStand]['selage1'][cnt_gc]=10
                    gc[iScn][iStand]['s2'][cnt_gc]=meta['LUT']['VEG_COMP_LYR_R1_POLY']['SPECIES_CD_1']['HW']
                    gc[iScn][iStand]['p2'][cnt_gc]=15
                    gc[iScn][iStand]['gain2'][cnt_gc]=gain
                    gc[iScn][iStand]['selage2'][cnt_gc]=10

                else:

                    # Interior
                    gc[iScn][iStand]['s1'][cnt_gc]=meta['LUT']['VEG_COMP_LYR_R1_POLY']['SPECIES_CD_1']['PL']
                    gc[iScn][iStand]['p1'][cnt_gc]=60
                    gc[iScn][iStand]['gain1'][cnt_gc]=gain
                    gc[iScn][iStand]['selage1'][cnt_gc]=10
                    gc[iScn][iStand]['s2'][cnt_gc]=meta['LUT']['VEG_COMP_LYR_R1_POLY']['SPECIES_CD_1']['SW']
                    gc[iScn][iStand]['p2'][cnt_gc]=40
                    gc[iScn][iStand]['gain2'][cnt_gc]=gain
                    gc[iScn][iStand]['selage2'][cnt_gc]=10

    # Get rid of rows with no info
    for iScn in range(meta[pNam]['Project']['N Scenario']):
        for iStand in range(meta[pNam]['Project']['N Stand']):
            ind=np.where(gc[iScn][iStand]['ID_Stand']!=-999)[0]
            for key in meta['Modules']['GYM']['GC_Variable_List']:
                gc[iScn][iStand][key]=gc[iScn][iStand][key][ind]

    print(str(np.round((time.time()-t0)/60,decimals=1)) + ' min')

    #--------------------------------------------------------------------------
    # Adjust mortality factors that only affect specific tree species
    #--------------------------------------------------------------------------

    if meta[pNam]['Project']['Adjust species-specific mortality']=='On':
        print('Adjusting mortality based on species-specific pests')
        t0=time.time()
        dmec=AdjustSpeciesSpecificMortality(meta,pNam,dmec,gc,meta[pNam]['Project']['Actual Indices'][0])
        print(str(np.round((time.time()-t0)/60,decimals=1)) + ' min')

    #--------------------------------------------------------------------------
    # Extract a set of unique growth curves
    # Decompose the full set of stands into a subset of unique stand types.
    # Exclude the first three columns, as they are all different.
    #--------------------------------------------------------------------------

    print('Extracting unique growth curves')
    t0=time.time()
    ugc=ExtractUniqueGrowthCurves(meta,pNam,gc)
    print(str(np.round((time.time()-t0)/60,decimals=1)) + ' min')

    #--------------------------------------------------------------------------
    # Export to BatchTIPSY parameter spreadsheet
    #--------------------------------------------------------------------------

    print('Exporting BatchTIPSY parameters to spreadsheet')
    t0=time.time()
    cbu.Write_BatchTIPSY_Input_Spreadsheet(meta,pNam,ugc)
    print(str(np.round((time.time()-t0)/60,decimals=1)) + ' min')

    #--------------------------------------------------------------------------
    # Populate BatchTIPSY.exe input variable (.dat) file
    #--------------------------------------------------------------------------

    print('Creating BatchTIPSY.exe input varialbe (.dat) file')
    t0=time.time()
    cbu.Write_BatchTIPSY_Input_File(meta,pNam)
    print(str(np.round((time.time()-t0)/60,decimals=1)) + ' min')

    return gc,ugc,dmec,inv

#%% Process project inputs 3

def ProcessProjectInputs3(meta,pNam,inv,dmec,lsc,gc,ugc):

    #--------------------------------------------------------------------------
    # Timber harvesting land base
    #--------------------------------------------------------------------------

    print('Defining timber harvesting landbase')
    t0=time.time()
    thlb=DefineTHLB(meta,pNam,inv,dmec,lsc)

    print(str(np.round((time.time()-t0)/60,decimals=1)) + ' min')

    # QA scenarios
    def Plot_THLB_Scenarios_For_QA():
        iScn=0
        plt.close('all')
        plt.plot(meta[pNam]['Year'],np.sum(thlb[iScn]['Actual'],axis=1),'r--')
        plt.plot(meta[pNam]['Year'],np.sum(thlb[iScn]['Baseline'],axis=1),'b-')

        plt.close('all')
        plt.plot(meta[pNam]['Year'],np.sum(thlb[iScn]['Actual'],axis=1),'r--')
        plt.plot(meta[pNam]['Year'],np.sum(thlb[iScn]['Actual WithDef'],axis=1),'c-')
        return

    #--------------------------------------------------------------------------
    # Prepare inventory
    #--------------------------------------------------------------------------

    print('Preparing inventory input files')
    t0=time.time()
    for iScn in range(meta[pNam]['Project']['N Scenario']):

        # Loop through batches, saving inventory to file
        for iBat in range(meta[pNam]['Project']['N Batch']):

            # Index to batch
            indBat=cbu.IndexToBatch(meta[pNam],iBat)
            N_StandsInBatch=len(indBat)

            inv1={}
            for k in inv.keys():
                inv1[k]=np.zeros((1,N_StandsInBatch),dtype='int16')
                inv1[k][0,:]=inv[k][indBat]

            # Land surface classification
            nam=meta[pNam]['Scenario'][iScn]['Land Surface Scenario']
            idx=LSC_Scenario_Crosswalk(lsc,nam)

            inv1['LSC']={}
            if nam!='None':
                inv1['LSC']['tv']=lsc['tv']
                inv1['LSC']['Cover']=np.reshape(lsc['Scenarios'][idx]['Cover'].copy(),(lsc['tv'].size,meta[pNam]['Project']['N Stand']))[:,indBat]
                inv1['LSC']['Use']=np.reshape(lsc['Scenarios'][idx]['Use'].copy(),(lsc['tv'].size,meta[pNam]['Project']['N Stand']))[:,indBat]

            # Timber harvesting landbase (1=yes, 0=no)
            inv1['THLB']=thlb[iScn][ meta[pNam]['Scenario'][iScn]['THLB Scenario'] ][:,indBat]

            # Temperature will be updated automatically
            inv1['MAT']=4*np.ones((1,N_StandsInBatch))

            # Sawtooth species-region samples
            if meta[pNam]['Project']['Biomass Module']=='Sawtooth':
                inv1['Srs1_ID']=meta['LUT']['Spc'][meta[pNam]['Scenario'][iScn]['SRS1_CD']]*np.ones((1,N_StandsInBatch),dtype='int16')
            else:
                inv1['Srs1_ID']=9999*np.ones((1,N_StandsInBatch),dtype='int16')

            inv1['Spc1_ID']=0*np.ones((1,N_StandsInBatch),dtype='int16')
            inv1['Spc1_Pct']=0*np.ones((1,N_StandsInBatch),dtype='int16')
            inv1['Spc2_ID']=0*np.ones((1,N_StandsInBatch),dtype='int16')
            inv1['Spc2_Pct']=0*np.ones((1,N_StandsInBatch),dtype='int16')
            inv1['Spc3_ID']=0*np.ones((1,N_StandsInBatch),dtype='int16')
            inv1['Spc3_Pct']=0*np.ones((1,N_StandsInBatch),dtype='int16')
            inv1['Spc3_ID']=0*np.ones((1,N_StandsInBatch),dtype='int16')
            inv1['Spc3_Pct']=0*np.ones((1,N_StandsInBatch),dtype='int16')

            # Probability of harvest (%/yr) from spatial map
            inv1['P Harvest Weight']=inv['P Harvest Weight'][indBat]

            # Save
            gu.opickle(meta['Paths'][pNam]['Input Scenario'][iScn] + '\\Inventory_Bat' + cbu.FixFileNum(iBat) + '.pkl',inv1)

    print(str(np.round((time.time()-t0)/60,decimals=1)) + ' min')

    #--------------------------------------------------------------------------
    # Simulate wildfires
    #--------------------------------------------------------------------------

    if (meta[pNam]['Scenario'][iScn]['Wildfire Status Pre-modern']=='On') | (meta[pNam]['Scenario'][iScn]['Wildfire Status Modern']=='On') | (meta[pNam]['Scenario'][iScn]['Wildfire Status Future']=='On'):
        if meta[pNam]['Project']['Frozen Ensembles Status']=='Off':
            print('Generating wildfire information')
            t0=time.time()
            asm.SimulateWildfireFromAAO(meta,pNam,inv)
            print(str(np.round((time.time()-t0)/60,decimals=1)) + ' min')

    #--------------------------------------------------------------------------
    # Simulate IBM
    #--------------------------------------------------------------------------

    if (meta[pNam]['Scenario'][iScn]['MPB Status Pre-modern']=='On') | (meta[pNam]['Scenario'][iScn]['MPB Status Modern']=='On') | (meta[pNam]['Scenario'][iScn]['MPB Status Future']=='On'):
        if meta[pNam]['Project']['Frozen Ensembles Status']=='Off':
            print('Generating MPB information')
            t0=time.time()
            asm.SimulateIBMFromAAO(meta,pNam,inv)
            print(str(np.round((time.time()-t0)/60,decimals=1)) + ' min')

    #--------------------------------------------------------------------------
    # Prepare disturbance/management event chronology
    #--------------------------------------------------------------------------

    print('Preparing DMEC input files')
    t0=time.time()
    for iEns in range(meta[pNam]['Project']['N Ensemble']):
        for iScn in range(meta[pNam]['Project']['N Scenario']):

            if (meta[pNam]['Scenario'][iScn]['Wildfire Status Pre-modern']=='On') | (meta[pNam]['Scenario'][iScn]['Wildfire Status Modern']=='On') | (meta[pNam]['Scenario'][iScn]['Wildfire Status Future']=='On'):

                # Import wildfire from aspatial stats model
                if meta[pNam]['Project']['Frozen Ensembles Status']=='Off':
                    wf_sim=gu.ipickle(meta['Paths'][pNam]['Data'] + '\\Inputs\\Ensembles\\wf_sim_Scn' + cbu.FixFileNum(iScn) + '_Ens' + cbu.FixFileNum(iEns) + '.pkl')
                else:
                    wf_sim=gu.ipickle(meta[pNam]['Project']['Frozen Ensembles Path'] + '\\wf_sim_Scn' + cbu.FixFileNum(iScn) + '_Ens' + cbu.FixFileNum(iEns) + '.pkl')
                if 'idx' in wf_sim:
                    idx=wf_sim['idx']
                    tmp=wf_sim.copy()
                    for v in ['Occurrence','Mortality']:
                        wf_sim[v]=np.zeros((meta[pNam]['Project']['N Time'],meta[pNam]['Project']['N Stand']),dtype='int8')
                        wf_sim[v][idx[0],idx[1]]=tmp[v]
                    del tmp

                # Import wildfire from onset-spread model
                if 'Use Wildfire from OSM' in meta[pNam]['Project']:

                    idx=gu.ipickle(meta['Paths'][pNam]['Data'] + '\\Inputs\\Ensembles\\wf_sim_osm_Scn' + cbu.FixFileNum(iScn) + '_Ens' + cbu.FixFileNum(iEns) + '.pkl')
                    wf_osm={}
                    wf_osm['Raw']=np.zeros( (lsc['Scenarios'][0]['Cover'].shape) ,dtype='int8')
                    wf_osm['Raw'][idx]=1
                    wf_osm['Raw']=np.reshape(wf_osm['Raw'],(lsc['tv'].size,meta[pNam]['Project']['N Stand']))

                    wf_osm['Occurrence']=np.zeros((meta[pNam]['Project']['N Time'],meta[pNam]['Project']['N Stand']),dtype='int8')
                    wf_osm['Mortality']=np.zeros((meta[pNam]['Project']['N Time'],meta[pNam]['Project']['N Stand']),dtype='int8')
                    for iT in range(lsc['tv'].size):
                        if lsc['tv'][iT]<2022:
                            continue
                        indT=np.where(meta[pNam]['Year']==lsc['tv'][iT])[0]
                        indS=np.where(wf_osm['Raw'][iT,:]==1)[0]
                        wf_osm['Occurrence'][indT,indS]=1
                        wf_osm['Mortality'][indT,indS]=100

            if (meta[pNam]['Scenario'][iScn]['MPB Status Pre-modern']=='On') | (meta[pNam]['Scenario'][iScn]['MPB Status Modern']=='On') | (meta[pNam]['Scenario'][iScn]['MPB Status Future']=='On'):

                # Import simulated mountain pine beetle
                if meta[pNam]['Project']['Frozen Ensembles Status']=='Off':
                    ibm_sim=gu.ipickle(meta['Paths'][pNam]['Data'] + '\\Inputs\\Ensembles\\ibm_sim_Scn' + cbu.FixFileNum(iScn) + '_Ens' + cbu.FixFileNum(iEns) + '.pkl')
                else:
                    ibm_sim=gu.ipickle(meta[pNam]['Project']['Frozen Ensembles Path'] + '\\ibm_sim_Scn' + cbu.FixFileNum(iScn) + '_Ens' + cbu.FixFileNum(iEns) + '.pkl')
                if 'idx' in ibm_sim:
                    idx=ibm_sim['idx']
                    tmp=ibm_sim.copy()
                    for v in ['Occurrence','Mortality']:
                        ibm_sim[v]=np.zeros((meta[pNam]['Project']['N Time'],meta[pNam]['Project']['N Stand']),dtype='int8')
                        ibm_sim[v][idx[0],idx[1]]=tmp[v]
                    del tmp

            for iBat in range(meta[pNam]['Project']['N Batch']):

                # Index to batch
                indBat=cbu.IndexToBatch(meta[pNam],iBat)

                # Initialize dictionary
                ec={}
                ec['ID Event Type']=np.zeros((meta[pNam]['Year'].size,indBat.size,meta['Core']['Max Events Per Year']),dtype='int16')
                ec['Mortality Factor']=np.zeros((meta[pNam]['Year'].size,indBat.size,meta['Core']['Max Events Per Year']),dtype='int16')
                ec['Growth Factor']=np.zeros((meta[pNam]['Year'].size,indBat.size,meta['Core']['Max Events Per Year']),dtype='int16')
                ec['ID Growth Curve']=np.zeros((meta[pNam]['Year'].size,indBat.size,meta['Core']['Max Events Per Year']),dtype='int16')
                tv=np.arange(meta[pNam]['Project']['Year Start'],meta[pNam]['Project']['Year End']+1,1)

                #----------------------------------------------------------
                # Spinup with constant return interval
                #----------------------------------------------------------

                if (meta[pNam]['Project']['Spinup Status']=='On'):

                    inv0=gu.ipickle(meta['Paths'][pNam]['Input Scenario'][iScn] + '\\Inventory_Bat' + cbu.FixFileNum(iBat) + '.pkl')

                    for iS in range(indBat.size):

                        # Index to stand
                        iStandFull=indBat[iS]

                        # Pre-industrial disturbance interval

                        # Old: regional
                        #if inv0['Region Code'][0,iS]==meta['LUT']['Region']['Coast']:
                        #    ivl_pi=300
                        #else:
                        #    ivl_pi=125

                        # New BGC zone-specific
                        cd=cbu.lut_n2s(meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE'],inv0['ID_BGCZ'][0,iS])[0]
                        ind=np.where(meta['Param']['BE']['SpinupRI']['Name']==cd)[0]
                        ivl_pi=meta['Param']['BE']['SpinupRI']['Value'][ind]

                        # Timing of transition between pre-industrial and modern periods
                        YearRef=dmec[iScn][iStandFull]['Year'][0]
                        #YearRef=np.random.randint(1775,high=1920,size=1,dtype=int)
                        AgeRef=ivl_pi

                        YrRegCyc=np.arange(YearRef-AgeRef-100*ivl_pi,YearRef-AgeRef+ivl_pi,ivl_pi)
                        Year=YrRegCyc[np.where(YrRegCyc>=meta[pNam]['Year'][0])[0]]
                        ID_Type=meta['LUT']['Event']['Wildfire']*np.ones(Year.size)
                        MortF=100*np.ones(Year.size)
                        GrowthF=0*np.ones(Year.size)
                        ID_GrowthCurve=1*np.ones(Year.size)
                        ec=cbu.CompileEvents(ec,tv,iS,ID_Type,Year,MortF,GrowthF,ID_GrowthCurve)

                #------------------------------------------------------------------
                # Add modern era events
                #------------------------------------------------------------------

                for iS in range(indBat.size):

                    # Index to stand
                    iStandFull=indBat[iS]

                    ind=np.where(dmec[iScn][iStandFull]['Scenario Affected']==1)[0]
                    if ind.size>0:
                        ID_Type=dmec[iScn][iStandFull]['ID Event Type'][ind]
                        Year=dmec[iScn][iStandFull]['Year'][ind]
                        MortF=dmec[iScn][iStandFull]['Mortality Factor'][ind]
                        GrowthF=dmec[iScn][iStandFull]['Growth Factor'][ind]
                        ID_GrowthCurve=dmec[iScn][iStandFull]['ID_GC'][ind]
                        ec=cbu.CompileEvents(ec,tv,iS,ID_Type,Year,MortF,GrowthF,ID_GrowthCurve)

                #------------------------------------------------------------------
                # Add simulated wildfire
                #------------------------------------------------------------------

                if (meta[pNam]['Scenario'][iScn]['Wildfire Status Pre-modern']=='On') | (meta[pNam]['Scenario'][iScn]['Wildfire Status Modern']=='On') | (meta[pNam]['Scenario'][iScn]['Wildfire Status Future']=='On'):

                    Occ=wf_sim['Occurrence'][:,indBat].copy()
                    Mort=wf_sim['Mortality'][:,indBat]
                    for iEPY in range(meta['Core']['Max Events Per Year']):

                        # Index to available spots with simulated wildfire mortality
                        ind=np.where( (ec['ID Event Type'][:,:,iEPY]==0) & (Occ==1) )

                        # Populate
                        ec['ID Event Type'][ind[0],ind[1],iEPY]=meta['LUT']['Event']['Wildfire']
                        ec['Mortality Factor'][ind[0],ind[1],iEPY]=Mort[ind[0],ind[1]]
                        ec['Growth Factor'][ind[0],ind[1],iEPY]=-999
                        ec['ID Growth Curve'][ind[0],ind[1],iEPY]=1

                        # Eliminate occurrence so that it is not populated again as loop
                        # through events per year continues
                        Occ[ind[0],ind[1]]=0

                    # Add simulated wildfire from onset-spread model
                    if 'Use Wildfire from OSM' in meta[pNam]['Project']:
                        Occ=wf_osm['Occurrence'][:,indBat].copy()
                        Mort=wf_osm['Mortality'][:,indBat]
                        for iEPY in range(meta['Core']['Max Events Per Year']):

                            # Index to available spots with simulated wildfire mortality
                            ind=np.where( (ec['ID Event Type'][:,:,iEPY]==0) & (Occ==1) )

                            # Populate
                            ec['ID Event Type'][ind[0],ind[1],iEPY]=meta['LUT']['Event']['Wildfire']
                            ec['Mortality Factor'][ind[0],ind[1],iEPY]=Mort[ind[0],ind[1]]
                            ec['Growth Factor'][ind[0],ind[1],iEPY]=-999
                            ec['ID Growth Curve'][ind[0],ind[1],iEPY]=1

                            # Eliminate occurrence so that it is not populated again as loop
                            # through events per year continues
                            Occ[ind[0],ind[1]]=0

                #------------------------------------------------------------------
                # Add simulated MPB
                #------------------------------------------------------------------

                if (meta[pNam]['Scenario'][iScn]['MPB Status Pre-modern']=='On') | (meta[pNam]['Scenario'][iScn]['MPB Status Modern']=='On') | (meta[pNam]['Scenario'][iScn]['MPB Status Future']=='On'):

                    Occ=ibm_sim['Occurrence'][:,indBat].copy()
                    Mort=ibm_sim['Mortality'][:,indBat]
                    for iEPY in range(meta['Core']['Max Events Per Year']):

                        # Index to available spots with simulated wildfire mortality
                        ind=np.where( (ec['ID Event Type'][:,:,iEPY]==0) & (Occ==1) )

                        # Populate
                        ec['ID Event Type'][ind[0],ind[1],iEPY]=meta['LUT']['Event']['IBM']
                        ec['Mortality Factor'][ind[0],ind[1],iEPY]=Mort[ind[0],ind[1]]
                        ec['Growth Factor'][ind[0],ind[1],iEPY]=-999
                        ec['ID Growth Curve'][ind[0],ind[1],iEPY]=1

                        # Eliminate occurrence so that it is not populated again as loop
                        # through events per year continues
                        Occ[ind[0],ind[1]]=0

                #------------------------------------------------------------------
                # Add future scheduled NOSE
                #------------------------------------------------------------------

                if 'NOSE Future' in meta[pNam]['Project']:

                    c,ia,ib=np.intersect1d(indBat,meta[pNam]['Project']['NOSE Future']['Stand Index'],return_indices=True)

                    if ia.size>0:

                        for iP in range(ia.size):

                            YearIncite=meta[pNam]['Project']['NOSE Future']['Year'][ib[iP]]
                            iT=np.where(tv==YearIncite)[0]
                            iAvailable=np.where(ec['ID Event Type'][iT,ia[iP],:]==0)[0]
                            ec['ID Event Type'][iT,ia[iP],iAvailable]=meta['LUT']['Event']['GC Switch']
                            ec['Mortality Factor'][iT,ia[iP],iAvailable]=0
                            ec['Growth Factor'][iT,ia[iP],iAvailable]=-10
                            ec['ID Growth Curve'][iT,ia[iP],iAvailable]=1

                            if np.isin(iScn,meta[pNam]['Project']['Actual Indices'])==True:
                                TimeBetweenFireAndPlant=2
                                iT=np.where(tv==YearIncite+TimeBetweenFireAndPlant)[0]
                                iAvailable=np.where(ec['ID Event Type'][iT,ia[iP],:]==0)[0]

                                ID_GC=np.max(gc[iScn][indBat[ia[iP]]]['ID_GC'])

                                ec['ID Event Type'][iT,ia[iP],iAvailable]=meta['LUT']['Event']['Planting']
                                ec['Mortality Factor'][iT,ia[iP],iAvailable]=0
                                ec['Growth Factor'][iT,ia[iP],iAvailable]=10
                                ec['ID Growth Curve'][iT,ia[iP],iAvailable]=ID_GC

                #------------------------------------------------------------------
                # Compress by indexing into the elements with information
                #------------------------------------------------------------------

                ec['idx']=np.where(ec['ID Event Type']>0)
                ec['ID Event Type']=ec['ID Event Type'][ec['idx']]
                ec['Mortality Factor']=ec['Mortality Factor'][ec['idx']]
                ec['Growth Factor']=ec['Growth Factor'][ec['idx']]
                ec['ID Growth Curve']=ec['ID Growth Curve'][ec['idx']]

                gu.opickle(meta['Paths'][pNam]['Input Scenario'][iScn] + '\\Events_Ens' + cbu.FixFileNum(iEns) + '_Bat' + cbu.FixFileNum(iBat) + '.pkl',ec)

    print(str(np.round((time.time()-t0)/60,decimals=1)) + ' min')

    #--------------------------------------------------------------------------
    # Prepare growth curves
    #--------------------------------------------------------------------------

    print('Preparing growth curve input files')
    t0=time.time()
    cbu.PrepGrowthCurvesUniqueForCBR(meta,pNam,ugc)
    print(str(np.round((time.time()-t0)/60,decimals=1)) + ' min')

    #--------------------------------------------------------------------------
    # Save data
    #--------------------------------------------------------------------------

    print('Saving input files')
    t0=time.time()
    gu.opickle(meta['Paths'][pNam]['Data'] + '\\Inputs\\Metadata.pkl',meta)
    gu.opickle(meta['Paths'][pNam]['Data'] + '\\Inputs\\Metadata_backup.pkl',meta)
    gu.opickle(meta['Paths'][pNam]['Data'] + '\\Inputs\\dmec.pkl',dmec)
    gu.opickle(meta['Paths'][pNam]['Data'] + '\\Inputs\\inv.pkl',inv)
    gu.opickle(meta['Paths'][pNam]['Data'] + '\\Inputs\\gc.pkl',gc)
    gu.opickle(meta['Paths'][pNam]['Data'] + '\\Inputs\\ugc.pkl',ugc)
    gu.opickle(meta['Paths'][pNam]['Data'] + '\\Inputs\\thlb.pkl',thlb)
    print(str(np.round((time.time()-t0)/60,decimals=1)) + ' min')

    #--------------------------------------------------------------------------
    # Delete all output files
    #--------------------------------------------------------------------------

    print('Deleting any output files')
    t0=time.time()
    cbu.DeleteAllOutputFiles(meta,pNam)
    print(str(np.round((time.time()-t0)/60,decimals=1)) + ' min')

    return meta,dmec,inv,thlb

#%% Timber harvesting land base

def DefineTHLB(meta,pNam,inv,dmec,lsc):

    thlb=[{}]*meta[pNam]['Project']['N Scenario']

    for iScn in range(meta[pNam]['Project']['N Scenario']):

        #------------------------------------------------------------------------------
        # Initialize THLB flags (THLB=1,Non-THLB=0)
        #------------------------------------------------------------------------------

        # Initially assume everything is in the THLB
        thlb[iScn]['Actual']=np.ones((meta[pNam]['Project']['N Time'],meta[pNam]['Project']['N Stand']),dtype='int8')
        thlb[iScn]['Baseline']=thlb[iScn]['Actual'].copy()

        thlb[iScn]['Scn1 Actual']=thlb[iScn]['Actual'].copy()
        thlb[iScn]['Scn1 Baseline']=thlb[iScn]['Actual'].copy()

        # Index to stands that are uneconomic
        iUneconomic=np.where(inv['SI']<=5)[0]
        if iUneconomic.size>0:
            # Remove uneconomic stands from THLB
            thlb[iScn]['Actual'][:,iUneconomic]=0
            thlb[iScn]['Baseline'][:,iUneconomic]=0
            thlb[iScn]['Scn1 Actual'][:,iUneconomic]=0
            thlb[iScn]['Scn1 Baseline'][:,iUneconomic]=0

        # Idenify stands that have been harvested
        has_been_harvested=np.zeros(meta[pNam]['Project']['N Stand'])
        for iStand in range(len(dmec)):
            ind=np.where( (dmec[iScn][iStand]['ID Event Type']==meta['LUT']['Event']['Harvest']) | (dmec[iScn][iStand]['ID Event Type']==meta['LUT']['Event']['Harvest Salvage']) )[0]
            if ind.size>0:
                has_been_harvested[iStand]=1

        # Index to stands that have not been harvested
        iNoHarv=np.where( (has_been_harvested==0) & (inv['SI']>5) )[0]

        # Use the ratio of THLB to non-THLB as an indicator of what will be harvested
        # among remaining primary forest
        ratio_thlb=22/60 # ratio of THLB to total forest (SOF)

        # Ratio of uneconomic to total (needed to adjust probability)
        corr=iUneconomic.size/meta[pNam]['Project']['N Stand']

        # Probability of evading harvest
        if iNoHarv.size>0:
            p_evade=(1-ratio_thlb-corr)*(meta[pNam]['Project']['N Stand']/iNoHarv.size)
        else:
            p_evade=(1-ratio_thlb-corr)

        # Random prediction of whether it will evade harvesting
        iRem=np.where(np.random.random(iNoHarv.size)<p_evade)[0]
        thlb[iScn]['Actual'][:,iNoHarv[iRem]]=0
        thlb[iScn]['Baseline'][:,iNoHarv[iRem]]=0

        thlb[iScn]['Scn1 Actual'][:,iNoHarv[iRem]]=0
        thlb[iScn]['Scn1 Baseline'][:,iNoHarv[iRem]]=0

        # np.sum(thlb['Actual'][0,:])/meta[pNam]['Project']['N Stand']

        #------------------------------------------------------------------------------
        # Actual (New based on REAR layer)
        #------------------------------------------------------------------------------

        # Initialize year of transition
        thlb_YearTransitionOut=np.zeros(meta[pNam]['Project']['N Stand'])

        ind=np.where(inv['THLB Layer']==1)[0]
        thlb_YearTransitionOut[ind]=1990

        # Conservation from land surface classification
        name=meta[pNam]['Scenario'][iScn]['Land Surface Scenario']
        if name!='None':
            idx=LSC_Scenario_Crosswalk(lsc,name)
            Use=np.reshape(lsc['Scenarios'][idx]['Use'].copy(),(lsc['tv'].size,meta[pNam]['Project']['N Stand']))
            ind=np.where( Use==meta['LUT']['LSC']['Use']['Conservation Consistent'] )
            if ind[0].size>0:
                for i in range(ind[0].size):
                    thlb_YearTransitionOut[ind[1][i]]=lsc['tv'][ind[0][i]]

        # Apply transition to actual THLB
        for j in range(thlb_YearTransitionOut.size):
            if thlb_YearTransitionOut[j]>0:
                it=np.where( (meta[pNam]['Year']>=thlb_YearTransitionOut[j]) )[0]
                thlb[iScn]['Actual'][it,j]=0
                thlb[iScn]['Scn1 Actual'][it,j]=0

        #------------------------------------------------------------------------------
        # Scn1 Actual (with deferrals + random areas to achieve 30 by 30)
        #------------------------------------------------------------------------------

        thlb_YearTransitionOut=np.zeros(meta[pNam]['Project']['N Stand'])

        ind=np.where(inv['THLB Layer']==2)[0]
        thlb_YearTransitionOut[ind]=2023

        # Apply transition to actual THLB
        for j in range(thlb_YearTransitionOut.size):
            if thlb_YearTransitionOut[j]>0:
                it=np.where( (meta[pNam]['Year']>=thlb_YearTransitionOut[j]) )[0]
                thlb[iScn]['Scn1 Actual'][it,j]=0

        #------------------------------------------------------------------------------
        # Baselines
        #------------------------------------------------------------------------------

        # Adjust the baseline so that simulated harvesting between 1995 and 2022 only
        # occurs in areas where the THLB was affected by value diversification
        for year in range(1990,2023,1):

            iT=np.where(meta[pNam]['Year']==year)[0]

            iS=np.where( (thlb[iScn]['Baseline'][iT,:]==1) & (thlb[iScn]['Actual'][iT,:]==1) )[1]
            thlb[iScn]['Baseline'][iT,iS]=0

            iS=np.where( (thlb[iScn]['Scn1 Baseline'][iT,:]==1) & (thlb[iScn]['Scn1 Actual'][iT,:]==1) )[1]
            thlb[iScn]['Scn1 Baseline'][iT,iS]=0

    return thlb

#%% Get index to a scenario in the land surface class list

def LSC_Scenario_Crosswalk(lsc,name):
    if 'Scenarios' not in lsc:
        return
    for i in range(len(lsc['Scenarios'])):
        if lsc['Scenarios'][i]['Name']==name:
            return i