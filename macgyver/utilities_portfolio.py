'''
Portfolio Utilities

'''

#%% Import modules

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import shutil
import openpyxl
from fcgadgets.macgyver import utilities_general as gu
from fcgadgets.cbrunner import cbrun_utilities as cbu
from fcgadgets.silviculture import economics as econo

#%% Import portfolio inputs from spreadsheet

def ImportPortfolio(meta):

    meta['Project']={}
    meta['Project']['N Scenario']=2
    
    #--------------------------------------------------------------------------
    # Import activity types
    #--------------------------------------------------------------------------
    
    d=pd.read_excel(meta['Paths']['Project'] + '\\Inputs\\ProjectConfig.xlsx',sheet_name='Activity Types',skiprows=0).to_dict('split')
    
    meta['Project']['Activities']={}
    
    # Index to columns to keep
    a=np.array(d['data'][1])
    ind=np.where( (a!='Activity description') & (a!='nan') )[0]
    for i in range(len(d['data'])):
        if d['data'][i][0][-1]!=':':
            tmp=np.array([])
            for j in range(len(ind)):
                tmp=np.append(tmp,d['data'][i][ind[j]])
            meta['Project']['Activities'][d['data'][i][0]]=tmp
    
    meta['Project']['Activities']['Activity ID']=meta['Project']['Activities']['Activity ID'].astype(int)
    
    #--------------------------------------------------------------------------
    # Import implementation levels
    #--------------------------------------------------------------------------
    
    d=pd.read_excel(meta['Paths']['Project'] + '\\Inputs\\ProjectConfig.xlsx',sheet_name='Implementation Levels',skiprows=1).to_dict('split')
    meta['Project']['AIL']={}
    
    meta['Project']['AIL']['Year']=np.arange(d['data'][1][0],d['data'][-1][0]+1,1,dtype=int)
    
    meta['Project']['AIL']['BAU']=np.zeros((meta['Project']['AIL']['Year'].size,12))
    meta['Project']['AIL']['CAP']=np.zeros((meta['Project']['AIL']['Year'].size,12))
    for i in range(0,meta['Project']['AIL']['Year'].size):
        cnt=0
        for j in range(1,12):
            meta['Project']['AIL']['BAU'][i,cnt]=d['data'][i+1][j]
            cnt=cnt+1
        cnt=0
        for j in range(12,22):    
            meta['Project']['AIL']['CAP'][i,cnt]=d['data'][i+1][j]
            cnt=cnt+1
    
    meta['Project']['AIL']['BAU']=np.nan_to_num(meta['Project']['AIL']['BAU'])
    meta['Project']['AIL']['CAP']=np.nan_to_num(meta['Project']['AIL']['CAP'])
    meta['Project']['AIL']['BAU']=meta['Project']['AIL']['BAU'].astype(int)
    meta['Project']['AIL']['CAP']=meta['Project']['AIL']['CAP'].astype(int)
    
    indBAU=np.where(np.sum(meta['Project']['AIL']['BAU'],axis=0)>0)[0]
    indCAP=np.where(np.sum(meta['Project']['AIL']['CAP'],axis=0)>0)[0]
    aid=np.array(d['data'][0][1:25])
    meta['Project']['AIL']['Activity ID']=np.append(aid[indBAU],aid[indCAP+11])
    meta['Project']['AIL']['Activity ID']=meta['Project']['AIL']['Activity ID'].astype(int)
    
    #--------------------------------------------------------------------------
    # Index to activity ID and Year for each unique stand
    #--------------------------------------------------------------------------
    
    # Number of unique activities
    meta['Project']['N AT BAU']=np.sum(np.sum(meta['Project']['AIL']['BAU'],axis=0)>0)
    meta['Project']['N AT CAP']=np.sum(np.sum(meta['Project']['AIL']['CAP'],axis=0)>0)
    meta['Project']['N AT']=meta['Project']['N AT BAU']+meta['Project']['N AT CAP']
    
    # Number of years with activities
    meta['Project']['N AT Years']=meta['Project']['AIL']['Year'].size
    
    return meta

#%% Prepare inputs for BatchTIPSY.exe

def PrepareInputsForBatchTIPSY(meta):
    
    # Create a function that will return the column corresponding to a variable name
    fin=meta['Paths']['Model Code'] + '\\Parameters\\GrowthCurvesTIPSY_Parameters_Template.xlsx'
    df_frmt=pd.read_excel(fin,sheet_name='Sheet1')
    gy_labels=df_frmt.loc[5,:].values

    def GetColumn(lab):
        ind=np.where(gy_labels==lab)[0]
        return int(ind+1)

    def isnan(x):        
        if (x.dtype!='float') | (x.dtype!='float'):
            y=False
        else:
            y=np.isnan(x)
        return y

    # Generate new spreadsheet from template
    fout=meta['Paths']['Project'] + '\\Inputs\\GrowthCurvesTIPSY_Parameters.xlsx'
    shutil.copy(fin,fout)
    
    # Load file, delete everything - start clean!
    xfile=openpyxl.load_workbook(fout)
    sheet=xfile['Sheet1']
    N_headers=7

    # Initiate counter
    cnt=1

    # meta['Project']['N Stand']=5

    # Define growth curve 1 (pre-project curves)
    for iStand in range(meta['Project']['N Stand']): 
        
        # Index to activity
        iActivity=meta['Project']['Idx']['AT'][iStand]
        
        for iScn in range(meta['Project']['N Scenario']):
        
            sheet.cell(row=cnt+N_headers,column=1).value=cnt  
            sheet.cell(row=cnt+N_headers,column=2).value=iStand+1 # Stand ID
            sheet.cell(row=cnt+N_headers,column=3).value=iScn+1 # Scenario ID
            sheet.cell(row=cnt+N_headers,column=4).value=1 # Growth curve
            
            vnam='regeneration_method'; vc=GetColumn(vnam)
            sheet.cell(row=cnt+N_headers,column=vc).value='N' # Regeneration type (N, C, P)                             
            
            vnam='s1'; vc=GetColumn(vnam)
            sheet.cell(row=cnt+N_headers,column=vc).value=meta['Project']['Activities']['Previous_Spc1_CD'][iActivity] # Species 1 code            
            
            vnam='p1'; vc=GetColumn(vnam)
            sheet.cell(row=cnt+N_headers,column=vc).value=meta['Project']['Activities']['Previous_Spc1_PCT'][iActivity] # Species 1 percent           
            
            vnam='i1'; vc=GetColumn(vnam)
            sheet.cell(row=cnt+N_headers,column=vc).value=meta['Project']['Activities']['Site index (m)'][iActivity] # Species 1 site index            
            
            if isnan(meta['Project']['Activities']['Previous_Spc2_CD'][iActivity])==False:
                vnam='s2'; vc=GetColumn(vnam)
                sheet.cell(row=cnt+N_headers,column=vc).value=meta['Project']['Activities']['Previous_Spc2_CD'][iActivity] # Species 2 code             
            
                vnam='p2'; vc=GetColumn(vnam)
                sheet.cell(row=cnt+N_headers,column=vc).value=meta['Project']['Activities']['Previous_Spc2_PCT'][iActivity] # Species 2 percent            
            
            if isnan(meta['Project']['Activities']['Previous_Spc3_CD'][iActivity])==False:
                vnam='s3'; vc=GetColumn(vnam)
                sheet.cell(row=cnt+N_headers,column=vc).value=meta['Project']['Activities']['Previous_Spc3_CD'][iActivity] # Species 3 code             
                
                vnam='p3'; vc=GetColumn(vnam)
                sheet.cell(row=cnt+N_headers,column=vc).value=meta['Project']['Activities']['Previous_Spc3_PCT'][iActivity] # Species 3 percent
            
            if isnan(meta['Project']['Activities']['Previous_Spc4_CD'][iActivity])==False:
                vnam='s4'; vc=GetColumn(vnam)
                sheet.cell(row=cnt+N_headers,column=vc).value=meta['Project']['Activities']['Previous_Spc4_CD'][iActivity] # Species 3 code             
                
                vnam='p4'; vc=GetColumn(vnam)
                sheet.cell(row=cnt+N_headers,column=vc).value=meta['Project']['Activities']['Previous_Spc4_PCT'][iActivity] # Species 3 percent
            
            vnam='init_density'; vc=GetColumn(vnam)            
            sheet.cell(row=cnt+N_headers,column=vc).value=int(meta['Project']['Activities']['Previous initial density (SPH)'][iActivity]) # Planting density    
            
            vnam='regen_delay'; vc=GetColumn(vnam)
            sheet.cell(row=cnt+N_headers,column=vc).value=0 # Regeneration delay            
            
            vnam='oaf1'; vc=GetColumn(vnam)
            sheet.cell(row=cnt+N_headers,column=vc).value=meta['Project']['Activities']['OAF1'][iActivity] # OAF1                        
            
            vnam='oaf2'; vc=GetColumn(vnam)
            sheet.cell(row=cnt+N_headers,column=vc).value=meta['Project']['Activities']['OAF2'][iActivity] # OAF2
            
            vnam='bec_zone'; vc=GetColumn(vnam)
            sheet.cell(row=cnt+N_headers,column=vc).value=meta['Project']['Activities']['BGC Zone'][iActivity] # BEC zone
            
            vnam='FIZ'; vc=GetColumn(vnam)
            if (meta['Project']['Activities']['BGC Zone'][iActivity]=='CWH') | (meta['Project']['Activities']['BGC Zone'][iActivity]=='CDF'):
                fiz='C'
            else:
                fiz='I'
            sheet.cell(row=cnt+N_headers,column=vc).value=fiz # FIZ
        
            # Update counter
            cnt=cnt+1 

    # Define growth curve 2 (baseline or project scenario)
    for iStand in range(meta['Project']['N Stand']):
        
        # Index to activity
        iActivity=meta['Project']['Idx']['AT'][iStand]
        
        for iScn in range(meta['Project']['N Scenario']):
    
            sheet.cell(row=cnt+N_headers,column=1).value=cnt
            sheet.cell(row=cnt+N_headers,column=2).value=iStand+1 # Stand ID
            sheet.cell(row=cnt+N_headers,column=3).value=iScn+1 # Scenario ID
            sheet.cell(row=cnt+N_headers,column=4).value=2 # Growth curve
        
            if meta['Project']['Activities']['New regeneration method'][iActivity]=='NAT':
                rm='N'
            else:
                rm='P'        
            vnam='regeneration_method'; vc=GetColumn(vnam)
            sheet.cell(row=cnt+N_headers,column=vc).value=rm # Regeneration type (N, C, P)                             
            
            vnam='s1'; vc=GetColumn(vnam)
            sheet.cell(row=cnt+N_headers,column=vc).value=meta['Project']['Activities']['New_Spc1_CD'][iActivity] # Species 1 code            
            
            vnam='p1'; vc=GetColumn(vnam)
            sheet.cell(row=cnt+N_headers,column=vc).value=meta['Project']['Activities']['New_Spc1_PCT'][iActivity] # Species 1 percent           
            
            vnam='i1'; vc=GetColumn(vnam)
            sheet.cell(row=cnt+N_headers,column=vc).value=meta['Project']['Activities']['Site index (m)'][iActivity] # Species 1 site index            
            
            if isnan(meta['Project']['Activities']['New_Spc2_CD'][iActivity])==False:
                vnam='s2'; vc=GetColumn(vnam)
                sheet.cell(row=cnt+N_headers,column=vc).value=meta['Project']['Activities']['New_Spc2_CD'][iActivity] # Species 2 code             
            
                vnam='p2'; vc=GetColumn(vnam)
                sheet.cell(row=cnt+N_headers,column=vc).value=meta['Project']['Activities']['New_Spc2_PCT'][iActivity] # Species 2 percent            
            
            if isnan(meta['Project']['Activities']['New_Spc3_CD'][iActivity])==False:
                vnam='s3'; vc=GetColumn(vnam)
                sheet.cell(row=cnt+N_headers,column=vc).value=meta['Project']['Activities']['New_Spc3_CD'][iActivity] # Species 3 code             
                
                vnam='p3'; vc=GetColumn(vnam)
                sheet.cell(row=cnt+N_headers,column=vc).value=meta['Project']['Activities']['New_Spc3_PCT'][iActivity] # Species 3 percent
            
            if isnan(meta['Project']['Activities']['New_Spc4_CD'][iActivity])==False:
                vnam='s4'; vc=GetColumn(vnam)
                sheet.cell(row=cnt+N_headers,column=vc).value=meta['Project']['Activities']['New_Spc4_CD'][iActivity] # Species 3 code             
                
                vnam='p4'; vc=GetColumn(vnam)
                sheet.cell(row=cnt+N_headers,column=vc).value=meta['Project']['Activities']['New_Spc4_PCT'][iActivity] # Species 3 percent
            
            vnam='init_density'; vc=GetColumn(vnam)            
            sheet.cell(row=cnt+N_headers,column=vc).value=int(meta['Project']['Activities']['New initial density (SPH)'][iActivity]) # Planting density    
            
            vnam='regen_delay'; vc=GetColumn(vnam)
            sheet.cell(row=cnt+N_headers,column=vc).value=meta['Project']['Activities']['New regeneration delay (years)'][iActivity] # Regeneration delay            
            
            vnam='oaf1'; vc=GetColumn(vnam)
            sheet.cell(row=cnt+N_headers,column=vc).value=meta['Project']['Activities']['OAF1'][iActivity] # OAF1                        
            
            vnam='oaf2'; vc=GetColumn(vnam)
            sheet.cell(row=cnt+N_headers,column=vc).value=meta['Project']['Activities']['OAF2'][iActivity] # OAF2
            
            vnam='bec_zone'; vc=GetColumn(vnam)
            sheet.cell(row=cnt+N_headers,column=vc).value=meta['Project']['Activities']['BGC Zone'][iActivity] # BEC zone
            
            vnam='FIZ'; vc=GetColumn(vnam)
            if (meta['Project']['Activities']['BGC Zone'][iActivity]=='CWH') | (meta['Project']['Activities']['BGC Zone'][iActivity]=='CDF'):
                fiz='C'
            else:
                fiz='I'
            sheet.cell(row=cnt+N_headers,column=vc).value=fiz # FIZ
    
        # Update counter
        cnt=cnt+1 

    # Save to spreadsheet
    xfile.save(fout)

    # Populate BatchTIPSY.exe input variable (.dat) file 
    cbu.Write_BatchTIPSY_Input_File(meta)
    
    return

#%% Prepare inventory

def PrepareInventory(meta):
    
    for iScn in range(meta['Project']['N Scenario']):
    
        # Loop through batches, saving inventory to file
        for iBat in range(meta['Project']['N Batch']):
      
            inv={}
    
            # Index to batch
            indBat=cbu.IndexToBatch(meta,iBat)    
            N_StandsInBatch=len(indBat)
    
            # Initialize inventory variables
            inv['Lat']=np.zeros((1,N_StandsInBatch))
            inv['Lon']=np.zeros((1,N_StandsInBatch))
            inv['X']=inv['Lat']
            inv['Y']=inv['Lon']
        
            # BEC zone
            inv['ID_BECZ']=np.zeros((1,N_StandsInBatch),dtype=np.int)
            for i in range(N_StandsInBatch):
                cd=meta['Project']['Activities']['BGC Zone'][meta['Project']['Idx']['AT'][indBat[i]]]
                inv['ID_BECZ'][0,indBat[i]]=meta['LUT']['VRI']['BEC_ZONE_CODE'][cd]
    
            # Timber harvesting landbase (1=yes, 0=no)
            inv['THLB']=1*np.ones((meta['Year'].size,N_StandsInBatch))
        
            # Temperature will be updated automatically
            inv['MAT']=4*np.ones((1,N_StandsInBatch))
            
            if meta['Project']['Biomass Module']=='Sawtooth':
                inv['Srs1_ID']=meta['LUT']['Spc'][meta['Scenario'][iScn]['SRS1_CD']]*np.ones((1,N_StandsInBatch),dtype=np.int)
            else:
                inv['Srs1_ID']=9999*np.ones((1,N_StandsInBatch),dtype=np.int)
                inv['Srs1_Pct']=100*np.ones((1,N_StandsInBatch),dtype=np.int)
                inv['Srs2_ID']=0*np.ones((1,N_StandsInBatch),dtype=np.int)
                inv['Srs2_Pct']=0*np.ones((1,N_StandsInBatch),dtype=np.int)
                inv['Srs3_ID']=0*np.ones((1,N_StandsInBatch),dtype=np.int)
                inv['Srs3_Pct']=0*np.ones((1,N_StandsInBatch),dtype=np.int)
                inv['Spc1_ID']=0*np.ones((1,N_StandsInBatch),dtype=np.int)
                inv['Spc1_Pct']=0*np.ones((1,N_StandsInBatch),dtype=np.int)
                inv['Spc2_ID']=0*np.ones((1,N_StandsInBatch),dtype=np.int)
                inv['Spc2_Pct']=0*np.ones((1,N_StandsInBatch),dtype=np.int)
                inv['Spc3_ID']=0*np.ones((1,N_StandsInBatch),dtype=np.int)
                inv['Spc3_Pct']=0*np.ones((1,N_StandsInBatch),dtype=np.int)

            # Save
            gu.opickle(meta['Paths']['Input Scenario'][iScn] + '\\Inventory_Bat' + cbu.FixFileNum(iBat) + '.pkl',inv)
    
    return

#%% Prepare event chronology

def PrepareEventChronology(meta):
    
    #--------------------------------------------------------------------------
    # Import custom harvest inputs
    #--------------------------------------------------------------------------
    
    meta['Harvest Custom']={}
    for i in range(meta['Project']['N AT']):
        
        iA=meta['Project']['Idx']['AT'][i]
        
        d={}
        d['BiomassMerch_Affected']=meta['Project']['Activities']['% of overstory that is felled'][iA]
        d['BiomassMerch_Removed']=meta['Project']['Activities']['% removed (merch biomass)'][iA]
        d['BiomassMerch_Burned']=meta['Project']['Activities']['% burned on site (merch biomass)'][iA]
        d['BiomassMerch_LeftOnSite']=meta['Project']['Activities']['% left to decay on site (merch biomass)'][iA]
        d['BiomassNonMerch_Affected']=meta['Project']['Activities']['% of overstory that is felled'][iA]
        d['BiomassNonMerch_Removed']=meta['Project']['Activities']['% removed (non merch biomass)'][iA]
        d['BiomassNonMerch_Burned']=meta['Project']['Activities']['% burned on site (non merch biomass)'][iA]
        d['BiomassNonMerch_LeftOnSite']=meta['Project']['Activities']['% left to decay on site (non merch biomass)'][iA]
        d['Snags_Affected']=meta['Project']['Activities']['% of overstory that is felled'][iA]
        d['Snags_Removed']=meta['Project']['Activities']['% removed (snags)'][iA]
        d['Snags_Burned']=meta['Project']['Activities']['% burned on site (snags)'][iA]
        d['Snags_LeftOnSite']=meta['Project']['Activities']['% left to decay on site (snags)'][iA]
        
        d['RemovedMerchToPulp']=5
        d['RemovedMerchToFuel']=5
        d['RemovedMerchToLumber']=5
        d['RemovedMerchToPlywood']=5
        d['RemovedMerchToOSB']=5
        d['RemovedMerchToMDF']=5
        d['RemovedMerchToCants']=5
        d['RemovedMerchToFirewood']=5
        d['RemovedNonMerchToFuel']=5
        d['RemovedNonMerchToLumber']=5
        d['RemovedNonMerchToPlywood']=5
        d['RemovedNonMerchToOSB']=5
        d['RemovedNonMerchToMDF']=5
        d['RemovedNonMerchToPulp']=5
        d['RemovedNonMerchToCants']=5
        d['RemovedNonMerchToFirewood']=5
        d['RemovedSnagStemToFuel']=5
        d['RemovedSnagStemToLumber']=5
        d['RemovedSnagStemToPlywood']=5
        d['RemovedSnagStemToOSB']=5
        d['RemovedSnagStemToMDF']=5
        d['RemovedSnagStemToPulp']=5
        d['RemovedSnagStemToCants']=5
        d['RemovedSnagStemToFirewood']=5
    
        # Populate custom harvest event
        meta['Harvest Custom'][int(i+1)]=d
    
    #--------------------------------------------------------------------------
    # Generate event chronology
    #--------------------------------------------------------------------------    
    
    for iScn in range(meta['Project']['N Scenario']):
        
        for iEns in range(meta['Project']['N Ensemble']):        
            
            for iBat in range(meta['Project']['N Batch']):
    
                # Index to batch
                indBat=cbu.IndexToBatch(meta,iBat)
                N_StandsInBatch=len(indBat)
    
                tv=np.arange(meta['Project']['Year Start'],meta['Project']['Year End']+1,1)
    
                # Initialize dictionary
                ec={}
                ec['ID_Type']=np.zeros((meta['Year'].size,indBat.size,meta['Core']['Max Events Per Year']),dtype='int16')
                ec['MortalityFactor']=np.zeros((meta['Year'].size,indBat.size,meta['Core']['Max Events Per Year']),dtype='int16')
                ec['GrowthFactor']=np.zeros((meta['Year'].size,indBat.size,meta['Core']['Max Events Per Year']),dtype='int16')
                ec['ID_GrowthCurve']=np.zeros((meta['Year'].size,indBat.size,meta['Core']['Max Events Per Year']),dtype='int16')
            
                for iS in range(N_StandsInBatch): 
                    
                    # Index to portfolio tables
                    iY=meta['Project']['Idx']['Year'][indBat[iS]]
                    iA=meta['Project']['Idx']['AT'][indBat[iS]]
                    
                    #----------------------------------------------------------
                    # Add spinup events
                    #----------------------------------------------------------
                    
                    ivl_spin=meta['Project']['Activities']['Spinup Disturbance Return Inverval'][iA]  
                    type_spin=meta['Project']['Activities']['Spinup Disturbance Type'][iA]
                    gcid_spin=meta['Project']['Activities']['Spinup Growth Curve ID'][iA]
                    
                    if np.isnan(meta['Project']['Activities']['Year of preceding disturbance'][iA])==0:                        
                        YearRef=meta['Project']['Activities']['Year of preceding disturbance'][iA]                        
                    else:
                        YearRef=meta['Project']['AIL']['Year'][iY]-meta['Project']['Activities']['Years before felling and/or planting'][iA]
                    
                    AgeRef=meta['Project']['Activities']['Stand age at time of preceding disturbance (years)'][iA]
                    
                    if AgeRef>=0:
                        Year=np.arange(YearRef-AgeRef-100*ivl_spin,YearRef-AgeRef+ivl_spin,ivl_spin)        
                    else:
                        Year1=meta['Project']['Year Start']+ivl_spin
                        Year2=meta['Spinup Year End']
                        Year=np.arange(Year1,Year2+1,ivl_spin)
                    
                    for iYr in range(Year.size):
                        iT=np.where(tv==Year[iYr])[0]
                        ec['ID_Type'][iT,iS,0]=meta['LUT']['Dist'][type_spin]
                        ec['MortalityFactor'][iT,iS,0]=100
                        ec['GrowthFactor'][iT,iS,0]=0
                        ec['ID_GrowthCurve'][iT,iS,0]=gcid_spin
                    
                    #----------------------------------------------------------
                    # Add preceding event
                    #----------------------------------------------------------
                    
                    type_pre=meta['Project']['Activities']['Preceding disturbance type'][iA]
                    mort_pre=meta['Project']['Activities']['% of biomass affected duirng preceding disturbance'][iA]
                    
                    iT=np.where(tv==YearRef)[0]
                    ec['ID_Type'][iT,iS,0]=meta['LUT']['Dist'][type_pre]
                    ec['MortalityFactor'][iT,iS,0]=mort_pre
                    ec['GrowthFactor'][iT,iS,0]=0
                    ec['ID_GrowthCurve'][iT,iS,0]=1
                    
                    #----------------------------------------------------------
                    # Add felling
                    #----------------------------------------------------------
                    
                    if iScn>0:
                        Year=meta['Project']['AIL']['Year'][iY]
                        iT=np.where(tv==Year)[0]
                        ec['ID_Type'][iT,iS,0]=meta['LUT']['Dist']['Harvest Custom ' + str(iA+1)]
                        ec['MortalityFactor'][iT,iS,0]=meta['Project']['Activities']['% of overstory that is felled'][iA]
                        ec['GrowthFactor'][iT,iS,0]=0
                        ec['ID_GrowthCurve'][iT,iS,0]=1
                        
                    #----------------------------------------------------------
                    # Add planting
                    #----------------------------------------------------------
                    
                    if iScn>0:
                        Year=meta['Project']['AIL']['Year'][iY]+meta['Project']['Activities']['Average time between felling and planting (years)'][iA]
                        iT=np.where(tv==Year)[0]
                        ec['ID_Type'][iT,iS,0]=meta['LUT']['Dist']['Planting']
                        ec['MortalityFactor'][iT,iS,0]=100
                        ec['GrowthFactor'][iT,iS,0]=0
                        ec['ID_GrowthCurve'][iT,iS,0]=2
            
                #--------------------------------------------------------------
                # Save to file            
                #--------------------------------------------------------------
                
                gu.opickle(meta['Paths']['Input Scenario'][iScn] + '\\Events_Ens' + cbu.FixFileNum(iEns) + '_Bat' + cbu.FixFileNum(iBat) + '.pkl',ec)
    
    return

#%% Import results
    
def ImportResults(meta):
    
    # Define time vector
    tv=np.arange(meta['Project']['Year Start Saving'],meta['Project']['Year End']+1,1)
    
    # Import cbrunner outputs
    v1=cbu.LoadScenarioResults(meta,[0,1])

    # Calculate GHG balance
    v2,meta=cbu.CalculateGHGBalance(v1,meta)        
    
    # Extract annual implementation level (ha/year)
    iBAU=np.where(np.sum(meta['Project']['AIL']['BAU'],axis=0)>0)[0]
    iCAP=np.where(np.sum(meta['Project']['AIL']['CAP'],axis=0)>0)[0]
    ail=np.column_stack([meta['Project']['AIL']['BAU'][:,iBAU],meta['Project']['AIL']['CAP'][:,iCAP]])

    # Initialize GHG balance for each activity type
    vAT=[None]*meta['Project']['N Scenario']
    
    for iScn in range(meta['Project']['N Scenario']): 
        
        # Calculate cashflow
        iEns=0 
        iBat=0
        inv=gu.ipickle(meta['Paths']['Input Scenario'][iScn] + '\\Inventory_Bat' + cbu.FixFileNum(iBat) + '.pkl')
        ec=gu.ipickle(meta['Paths']['Input Scenario'][iScn] + '\\Events_Ens' + cbu.FixFileNum(iEns) + '_Bat' + cbu.FixFileNum(iBat) + '.pkl')            
        ec=cbu.EventChronologyDecompress(meta,ec,iScn,iEns,iBat)
        econ=econo.CalculateNetRevenue(meta,iScn,iEns,iBat,inv,ec,v1[iScn])
        
        # Initialize dictionary to store results for each unique AT
        vAT[iScn]={}
        vAT[iScn]['Sum']={}
        vAT[iScn]['Per Ha']={}
        
        List=['V_StemMerch','C_RemovedMerch','C_RemovedNonMerch','C_RemovedSnagStem']
        for k in List:         
            vAT[iScn]['Per Ha'][k]=np.zeros( (v2[iScn]['Year'].size,meta['Project']['N AT Years'],meta['Project']['N AT']) )
            vAT[iScn]['Sum'][k]=np.zeros( (v2[iScn]['Year'].size,meta['Project']['N AT Years'],meta['Project']['N AT']) )            
            for iYr in range(meta['Project']['AIL']['Year'].size):                
                for iAT in range(meta['Project']['N AT']):                
                    indATY=np.where( (meta['Project']['Idx']['Year']==iYr) & (meta['Project']['Idx']['AT Unique']==iAT) )[0]                
                    if ail[iYr,iAT]>0:
                        vAT[iScn]['Per Ha'][k][:,iYr,iAT]=np.nanmean(v1[iScn][k][:,indATY],axis=1)                    
                    vAT[iScn]['Sum'][k][:,iYr,iAT]=ail[iYr,iAT]*np.nanmean(v1[iScn][k][:,indATY],axis=1)            
            vAT[iScn]['Per Ha'][k]=np.nanmean(vAT[iScn]['Per Ha'][k],axis=1)
            vAT[iScn]['Sum'][k]=np.nansum(vAT[iScn]['Sum'][k],axis=1)
        
        for k in v2[0].keys():
            if k=='Year':
                continue            
            vAT[iScn]['Per Ha'][k]=np.zeros( (v2[iScn]['Year'].size,meta['Project']['N AT Years'],meta['Project']['N AT']) )
            vAT[iScn]['Sum'][k]=np.zeros( (v2[iScn]['Year'].size,meta['Project']['N AT Years'],meta['Project']['N AT']) )            
            for iYr in range(meta['Project']['AIL']['Year'].size):                
                for iAT in range(meta['Project']['N AT']):                
                    indATY=np.where( (meta['Project']['Idx']['Year']==iYr) & (meta['Project']['Idx']['AT Unique']==iAT) )[0]                
                    if ail[iYr,iAT]>0:
                        vAT[iScn]['Per Ha'][k][:,iYr,iAT]=np.nanmean(v2[iScn][k][:,indATY],axis=1)                    
                    vAT[iScn]['Sum'][k][:,iYr,iAT]=ail[iYr,iAT]*np.nanmean(v2[iScn][k][:,indATY],axis=1)            
            vAT[iScn]['Per Ha'][k]=np.nanmean(vAT[iScn]['Per Ha'][k],axis=1)
            vAT[iScn]['Sum'][k]=np.nansum(vAT[iScn]['Sum'][k],axis=1)
        
        for k in econ.keys():            
            if k=='Year':
                continue            
            vAT[iScn]['Per Ha'][k]=np.zeros( (v2[iScn]['Year'].size,meta['Project']['N AT Years'],meta['Project']['N AT']) )
            vAT[iScn]['Sum'][k]=np.zeros( (v2[iScn]['Year'].size,meta['Project']['N AT Years'],meta['Project']['N AT']) )            
            for iYr in range(meta['Project']['AIL']['Year'].size):                
                for iAT in range(meta['Project']['N AT']):                
                    indATY=np.where( (meta['Project']['Idx']['Year']==iYr) & (meta['Project']['Idx']['AT Unique']==iAT) )[0]                
                    if ail[iYr,iAT]>0:
                        vAT[iScn]['Per Ha'][k][:,iYr,iAT]=np.nanmean(econ[k][:,indATY],axis=1)                    
                    vAT[iScn]['Sum'][k][:,iYr,iAT]=ail[iYr,iAT]*np.nanmean(econ[k][:,indATY],axis=1)            
            vAT[iScn]['Per Ha'][k]=np.nanmean(vAT[iScn]['Per Ha'][k],axis=1)
            vAT[iScn]['Sum'][k]=np.nansum(vAT[iScn]['Sum'][k],axis=1)

    # Results by portfolio
    iBAU=meta['Project']['Idx']['BAU AT']
    iCAP=meta['Project']['Idx']['CAP AT']
    
    vBAU=[None]*meta['Project']['N Scenario']
    vCAP=[None]*meta['Project']['N Scenario']
    for iScn in range(meta['Project']['N Scenario']):
        vBAU[iScn]={}
        vBAU[iScn]['Sum']={}
        vBAU[iScn]['W Ave']={}        
        vCAP[iScn]={}
        vCAP[iScn]['Sum']={}
        vCAP[iScn]['W Ave']={}
        for k in vAT[0]['Sum'].keys():
            if k=='Year':
                continue
            vBAU[iScn]['Sum'][k]=np.sum(vAT[iScn]['Sum'][k][:,iBAU],axis=1)
            vBAU[iScn]['W Ave'][k]=0
            vCAP[iScn]['Sum'][k]=np.sum(vAT[iScn]['Sum'][k][:,iCAP],axis=1)
            vCAP[iScn]['W Ave'][k]=0 

    return tv,vAT,vBAU,vCAP,meta

#%% Export tabular results
    
def ExportResultsTabular(meta,vAT,tv,t0,t1):
    
    #vList=['Year Calendar','Sec_NGHGB','Revenue Net','Revenue Gross','Cost Total']
    sList=['BAU','CAP']

    iT=np.where( (tv>=t0) & (tv<=t1) )[0]

    for iScn in range(meta['Project']['N Scenario']):

        for iAT in range(meta['Project']['N AT']):
            d={}
            d['Year Calendar']=tv[iT]
            for k in vAT[iScn]['Sum'].keys():
                d[k]=vAT[iScn]['Sum'][k][iT,iAT]
            df=pd.DataFrame.from_dict(d)
            df.to_excel(meta['Paths']['Project'] + '\\Outputs\\Activity' + str(iAT+1) + '_Scenario' + sList[iScn] + '_Sum.xlsx')
            
            d={}
            d['Year Calendar']=tv[iT]
            for k in vAT[iScn]['Per Ha'].keys():
                d[k]=vAT[iScn]['Per Ha'][k][iT,iAT]
            df=pd.DataFrame.from_dict(d)
            df.to_excel(meta['Paths']['Project'] + '\\Outputs\\Activity' + str(iAT+1) + '_Scenario' + sList[iScn] + '_PerHa.xlsx')
    
    return

#%% Plot results

def PlotResults(meta,vAT,vBAU,vCAP,tv,it,iB,iP):
    
    t_disc=np.maximum(0,tv-2021)

    fig,ax=plt.subplots(7,3,figsize=gu.cm2inch(18,20))
    
    cl_b=[0,0,0.8]; cl_p=[1,0,0]; cl_d=[0.4,0.8,0]
    
    #------------------------------------------------------------------------
    vr='Eco_Biomass'
    ax[0,0].plot(tv[it],vBAU[iB]['Sum'][vr][it]/1e6,'-',color=cl_b,label='Baseline')
    ax[0,0].plot(tv[it],vBAU[iP]['Sum'][vr][it]/1e6,'--',color=cl_p,label='Project')
    ax[0,0].set(ylabel='Biomass')
    ax[0,0].set_title('BAU Scenario')
    ax[0,0].legend(loc='lower right',frameon=False,facecolor=None)
    
    vr='Eco_Biomass'
    ax[0,1].plot(tv[it],vCAP[iB]['Sum'][vr][it]/1e6,'-',color=cl_b)
    ax[0,1].plot(tv[it],vCAP[iP]['Sum'][vr][it]/1e6,'--',color=cl_p)
    ax[0,1].set(ylabel='Biomass')
    ax[0,1].set_title('CAP Scenario')
    
    vr='Eco_Biomass'
    #ax[0,2].plot(tv[it],(vCAP[iB]['Sum'][vr][it]-vBAU[iB]['Sum'][vr][it])/1e6,'-',color=cl_d)
    ax[0,2].set(ylabel='$\Delta$ Biomass')
    ax[0,2].set_title('Scenario Difference')
    
    #------------------------------------------------------------------------
    vr='Sec_NGHGB'
    ax[1,0].plot(tv[it],vBAU[iB]['Sum'][vr][it]/1e6,'-',color=cl_b)
    ax[1,0].plot(tv[it],vBAU[iP]['Sum'][vr][it]/1e6,'--',color=cl_p)
    ax[1,0].set(ylabel='Annual GHG balance')
    
    vr='Sec_NGHGB'
    ax[1,1].plot(tv[it],vCAP[iB]['Sum'][vr][it]/1e6,'-',color=cl_b)
    ax[1,1].plot(tv[it],vCAP[iP]['Sum'][vr][it]/1e6,'--',color=cl_p)
    ax[1,1].set(ylabel='Annual GHG balance')
    
    #------------------------------------------------------------------------
    vr='Sec_NGHGB'
    ax[2,0].plot(tv[it],np.cumsum(vBAU[iB]['Sum'][vr][it]/1e6),'-',color=cl_b)
    ax[2,0].plot(tv[it],np.cumsum(vBAU[iP]['Sum'][vr][it]/1e6),'--',color=cl_p)
    ax[2,0].set(ylabel='Cumu. GHG balance')
    
    ax[2,1].plot(tv[it],np.cumsum(vCAP[iB]['Sum'][vr][it]/1e6),'-',color=cl_b)
    ax[2,1].plot(tv[it],np.cumsum(vCAP[iP]['Sum'][vr][it]/1e6),'--',color=cl_p)
    ax[2,1].set(ylabel='Cumu. GHG balance')
    
    #------------------------------------------------------------------------
    vr='Sec_NGHGB'
    a=np.cumsum( (vBAU[iP]['Sum'][vr][it]-vBAU[iB]['Sum'][vr][it])/1e6 )
    ax[3,0].plot(tv[it],a,'-',color=cl_d)
    ax[3,0].set(ylabel='$\Delta$ Cumu. GHG benefit')
    
    b=np.cumsum( (vCAP[iP]['Sum'][vr][it]-vCAP[iB]['Sum'][vr][it])/1e6 )
    ax[3,1].plot(tv[it],b,'-',color=cl_d)
    ax[3,1].set(ylabel='$\Delta$ Cumu. GHG benefit')
    
    ax[3,2].plot(tv[it],b-a,'-',color=cl_d)
    ax[3,2].set(ylabel='$\Delta$ Cumu. GHG benefit')
    
    #------------------------------------------------------------------------
    
    cost_b_bau=(vBAU[iB]['Sum']['Cost Planting'][it]+vBAU[iB]['Sum']['Cost Knockdown'][it]+vBAU[iB]['Sum']['Cost Ripping'][it])/1e6
    cost_p_bau=(vBAU[iP]['Sum']['Cost Planting'][it]+vBAU[iP]['Sum']['Cost Knockdown'][it]+vBAU[iP]['Sum']['Cost Ripping'][it])/1e6
    cost_b_bau=(cost_b_bau)/(1+0.03**t_disc[it])
    cost_p_bau=(cost_p_bau)/(1+0.03**t_disc[it])
    ax[4,0].plot(tv[it],cost_b_bau,'-',color=cl_b)
    ax[4,0].plot(tv[it],cost_p_bau,'--',color=cl_p)
    ax[4,0].set(ylabel='Cost (CDN$M/year)')
    
    cost_b_cap=(vCAP[iB]['Sum']['Cost Planting'][it]+vCAP[iB]['Sum']['Cost Knockdown'][it]+vCAP[iB]['Sum']['Cost Ripping'][it])/1e6
    cost_p_cap=(vCAP[iP]['Sum']['Cost Planting'][it]+vCAP[iP]['Sum']['Cost Knockdown'][it]+vCAP[iP]['Sum']['Cost Ripping'][it])/1e6
    cost_b_cap=(cost_b_cap)/(1+0.03**t_disc[it])
    cost_p_cap=(cost_p_cap)/(1+0.03**t_disc[it])
    ax[4,1].plot(tv[it],cost_b_cap,'-',color=cl_b)
    ax[4,1].plot(tv[it],cost_p_cap,'--',color=cl_p)
    ax[4,1].set(ylabel='Cost (CDN$M/year)')
    
    ax[4,2].plot(tv[it],cost_p_cap-cost_p_bau,'-',color=cl_d)
    ax[4,2].set(ylabel='$\Delta$ cost (Million CAD/year)')
    
    nr_b_bau=np.cumsum( (vBAU[iB]['Sum']['Revenue Net'][it]/1e6)/(1+0.03**t_disc[it]) )
    nr_p_bau=np.cumsum( (vBAU[iP]['Sum']['Revenue Net'][it]/1e6)/(1+0.03**t_disc[it]) )
    nr_b_cap=np.cumsum( (vCAP[iB]['Sum']['Revenue Net'][it]/1e6)/(1+0.03**t_disc[it]) )
    nr_p_cap=np.cumsum( (vCAP[iP]['Sum']['Revenue Net'][it]/1e6)/(1+0.03**t_disc[it]) )
    ax[5,0].plot(tv[it],nr_b_bau,'-',color=cl_b)
    ax[5,0].plot(tv[it],nr_p_bau,'--',color=cl_p)
    ax[5,0].set(ylabel='Cumu. net revenue')
    
    ax[5,1].plot(tv[it],nr_b_cap,'-',color=cl_b)
    ax[5,1].plot(tv[it],nr_p_cap,'--',color=cl_p)
    ax[5,1].set(ylabel='Cumu. net revenue')
    
    ax[5,2].plot(tv[it],nr_p_cap-nr_p_bau,'-',color=cl_d)
    ax[5,2].set(ylabel='Cumu. net revenue')
    
    #--------------------------------------------------------------------------
    # Harvest volume
    
    v_b_bau=vBAU[iB]['Sum']['C_RemovedMerch'][it]+vBAU[iB]['Sum']['C_RemovedNonMerch'][it]+vBAU[iB]['Sum']['C_RemovedSnagStem'][it]
    v_p_bau=vBAU[iP]['Sum']['C_RemovedMerch'][it]+vBAU[iP]['Sum']['C_RemovedNonMerch'][it]+vBAU[iP]['Sum']['C_RemovedSnagStem'][it]
    
    v_b_cap=vCAP[iB]['Sum']['C_RemovedMerch'][it]+vCAP[iB]['Sum']['C_RemovedNonMerch'][it]+vCAP[iB]['Sum']['C_RemovedSnagStem'][it]
    v_p_cap=vCAP[iP]['Sum']['C_RemovedMerch'][it]+vCAP[iP]['Sum']['C_RemovedNonMerch'][it]+vCAP[iP]['Sum']['C_RemovedSnagStem'][it]
    
    ax[6,0].plot(tv[it],v_b_bau/1e6,'-',color=cl_b)
    ax[6,0].plot(tv[it],v_p_bau/1e6,'--',color=cl_p)
    ax[6,0].set(ylabel='Harvest volume (m3/yr)')
    
    ax[6,1].plot(tv[it],v_b_bau/1e6,'-',color=cl_b)
    ax[6,1].plot(tv[it],v_p_bau/1e6,'--',color=cl_p)
    ax[6,1].set(ylabel='Harvest volume (m3/yr)')
    
    ax[6,2].plot(tv[it],(v_p_cap-v_p_bau)/1e6,'-',color=cl_d)
    ax[6,2].set(ylabel='$\Delta$ harvest volume (m3/yr)')    
    
    plt.tight_layout()
    
    return
