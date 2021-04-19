'''
Portfolio Tool

'''

#%% Import modules

import numpy as np
import pandas as pd
import shutil
import openpyxl
from fcgadgets.macgyver import utilities_general as gu
from fcgadgets.macgyver import utilities_gis as gis
from fcgadgets.macgyver import utilities_inventory as invu
from fcgadgets.cbrunner import cbrun_utilities as cbu

meta={}
meta['Paths']={}
meta['Paths']['Project']=r'C:\Users\rhember\Documents\Data\FCI_Projects\FCI_Portfolio'
meta['Paths']['Model Code']=r'C:\Users\rhember\Documents\Code_Python\fcgadgets\cbrunner'

#%% Import portfolio inputs from spreadsheet

def ImportPortfolio(meta):

    meta['N Scenario']=2
    
    meta['Portfolio']={}
    
    #--------------------------------------------------------------------------
    # Import activity types
    #--------------------------------------------------------------------------
    
    d=pd.read_excel(meta['Paths']['Project'] + '\\Inputs\\Portfolio Inputs.xlsx',sheet_name='Activity Types',skiprows=1).to_dict('split')
    meta['Portfolio']['Activities']={}
    
    meta['Portfolio']['Activities']['Portfolio Code']=np.array(d['columns'][1:])
    for i in range(len(meta['Portfolio']['Activities']['Portfolio Code'])):
        meta['Portfolio']['Activities']['Portfolio Code'][i]=np.array(meta['Portfolio']['Activities']['Portfolio Code'][i][0:3])
    
    for i in range(len(d['data'])):
        meta['Portfolio']['Activities'][d['data'][i][0]]=np.array(d['data'][i][1:])
    
    #--------------------------------------------------------------------------
    # Import implementation levels
    #--------------------------------------------------------------------------
    
    d=pd.read_excel(meta['Paths']['Project'] + '\\Inputs\\Portfolio Inputs.xlsx',sheet_name='Annual Implementation Level',skiprows=2).to_dict('split')
    meta['Portfolio']['AIL']={}
    meta['Portfolio']['AIL']['Year']=np.arange(d['data'][0][0],d['data'][-1][0]+1,1,dtype=int)
    meta['Portfolio']['AIL']['BAU']=np.zeros((meta['Portfolio']['AIL']['Year'].size,10))
    meta['Portfolio']['AIL']['CAP']=np.zeros((meta['Portfolio']['AIL']['Year'].size,10))
    for i in range(meta['Portfolio']['AIL']['Year'].size):
        cnt=0
        for j in range(1,11):
            meta['Portfolio']['AIL']['BAU'][i,cnt]=d['data'][i][j]
            cnt=cnt+1
        cnt=0
        for j in range(11,21):    
            meta['Portfolio']['AIL']['CAP'][i,cnt]=d['data'][i][j]
            cnt=cnt+1
    meta['Portfolio']['AIL']['BAU']=np.nan_to_num(meta['Portfolio']['AIL']['BAU'])
    meta['Portfolio']['AIL']['CAP']=np.nan_to_num(meta['Portfolio']['AIL']['CAP'])
    
    #--------------------------------------------------------------------------
    # Index to activity ID and Year for each unique stand
    #--------------------------------------------------------------------------
    
    meta['Portfolio']['N AT']=meta['Portfolio']['Activities']['Activity ID:'].size
    meta['Portfolio']['N Year']=meta['Portfolio']['AIL']['Year'].size
    
    meta['N Stand']=meta['Portfolio']['N AT']*meta['Portfolio']['N Year']
    
    meta['Portfolio']['Indices']={}
    meta['Portfolio']['Indices']['Activity']=np.zeros(meta['N Stand'],dtype=int)
    meta['Portfolio']['Indices']['Year']=np.zeros(meta['N Stand'],dtype=int)
    cnt=0
    for iA in range(meta['Portfolio']['N AT']):
        for iY in range(meta['Portfolio']['N Year']):
            meta['Portfolio']['Indices']['Activity'][cnt]=iA
            meta['Portfolio']['Indices']['Year'][cnt]=iY
            cnt=cnt+1
    
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

    # meta['N Stand Full']=5

    # Define growth curve 1 (pre-project curves)
    for iStand in range(meta['N Stand Full']): 
        
        # Index to activity
        iActivity=meta['Portfolio']['Indices']['Activity'][iStand]
        
        for iScn in range(meta['N Scenario']):
        
            sheet.cell(row=cnt+N_headers,column=1).value=cnt  
            sheet.cell(row=cnt+N_headers,column=2).value=iStand+1 # Stand ID
            sheet.cell(row=cnt+N_headers,column=3).value=iScn+1 # Scenario ID
            sheet.cell(row=cnt+N_headers,column=4).value=1 # Growth curve
            
            vnam='regeneration_method'; vc=GetColumn(vnam)
            sheet.cell(row=cnt+N_headers,column=vc).value='N' # Regeneration type (N, C, P)                             
            
            vnam='s1'; vc=GetColumn(vnam)
            sheet.cell(row=cnt+N_headers,column=vc).value=meta['Portfolio']['Activities']['Previous_Spc1_CD'][iActivity] # Species 1 code            
            
            vnam='p1'; vc=GetColumn(vnam)
            sheet.cell(row=cnt+N_headers,column=vc).value=meta['Portfolio']['Activities']['Previous_Spc1_PCT'][iActivity] # Species 1 percent           
            
            vnam='i1'; vc=GetColumn(vnam)
            sheet.cell(row=cnt+N_headers,column=vc).value=meta['Portfolio']['Activities']['Site index (m)'][iActivity] # Species 1 site index            
            
            if isnan(meta['Portfolio']['Activities']['Previous_Spc2_CD'][iActivity])==False:
                vnam='s2'; vc=GetColumn(vnam)
                sheet.cell(row=cnt+N_headers,column=vc).value=meta['Portfolio']['Activities']['Previous_Spc2_CD'][iActivity] # Species 2 code             
            
                vnam='p2'; vc=GetColumn(vnam)
                sheet.cell(row=cnt+N_headers,column=vc).value=meta['Portfolio']['Activities']['Previous_Spc2_PCT'][iActivity] # Species 2 percent            
            
            if isnan(meta['Portfolio']['Activities']['Previous_Spc3_CD'][iActivity])==False:
                vnam='s3'; vc=GetColumn(vnam)
                sheet.cell(row=cnt+N_headers,column=vc).value=meta['Portfolio']['Activities']['Previous_Spc3_CD'][iActivity] # Species 3 code             
                
                vnam='p3'; vc=GetColumn(vnam)
                sheet.cell(row=cnt+N_headers,column=vc).value=meta['Portfolio']['Activities']['Previous_Spc3_PCT'][iActivity] # Species 3 percent
            
            if isnan(meta['Portfolio']['Activities']['Previous_Spc4_CD'][iActivity])==False:
                vnam='s4'; vc=GetColumn(vnam)
                sheet.cell(row=cnt+N_headers,column=vc).value=meta['Portfolio']['Activities']['Previous_Spc4_CD'][iActivity] # Species 3 code             
                
                vnam='p4'; vc=GetColumn(vnam)
                sheet.cell(row=cnt+N_headers,column=vc).value=meta['Portfolio']['Activities']['Previous_Spc4_PCT'][iActivity] # Species 3 percent
            
            vnam='init_density'; vc=GetColumn(vnam)            
            sheet.cell(row=cnt+N_headers,column=vc).value=int(meta['Portfolio']['Activities']['Previous initial density (SPH)'][iActivity]) # Planting density    
            
            vnam='regen_delay'; vc=GetColumn(vnam)
            sheet.cell(row=cnt+N_headers,column=vc).value=0 # Regeneration delay            
            
            vnam='oaf1'; vc=GetColumn(vnam)
            sheet.cell(row=cnt+N_headers,column=vc).value=meta['Portfolio']['Activities']['OAF1'][iActivity] # OAF1                        
            
            vnam='oaf2'; vc=GetColumn(vnam)
            sheet.cell(row=cnt+N_headers,column=vc).value=meta['Portfolio']['Activities']['OAF2'][iActivity] # OAF2
            
            vnam='bec_zone'; vc=GetColumn(vnam)
            sheet.cell(row=cnt+N_headers,column=vc).value=meta['Portfolio']['Activities']['BGC Zone'][iActivity] # BEC zone
            
            vnam='FIZ'; vc=GetColumn(vnam)
            if (meta['Portfolio']['Activities']['BGC Zone'][iActivity]=='CWH') | (meta['Portfolio']['Activities']['BGC Zone'][iActivity]=='CDF'):
                fiz='C'
            else:
                fiz='I'
            sheet.cell(row=cnt+N_headers,column=vc).value=fiz # FIZ
        
            # Update counter
            cnt=cnt+1 

    # Define growth curve 2 (baseline or project scenario)
    for iStand in range(meta['N Stand Full']):
        
        # Index to activity
        iActivity=meta['Portfolio']['Indices']['Activity'][iStand]
        
        for iScn in range(meta['N Scenario']):
    
            sheet.cell(row=cnt+N_headers,column=1).value=cnt
            sheet.cell(row=cnt+N_headers,column=2).value=iStand+1 # Stand ID
            sheet.cell(row=cnt+N_headers,column=3).value=iScn+1 # Scenario ID
            sheet.cell(row=cnt+N_headers,column=4).value=2 # Growth curve
        
            if meta['Portfolio']['Activities']['New regeneration method'][iActivity]=='NAT':
                rm='N'
            else:
                rm='P'        
            vnam='regeneration_method'; vc=GetColumn(vnam)
            sheet.cell(row=cnt+N_headers,column=vc).value=rm # Regeneration type (N, C, P)                             
            
            vnam='s1'; vc=GetColumn(vnam)
            sheet.cell(row=cnt+N_headers,column=vc).value=meta['Portfolio']['Activities']['New_Spc1_CD'][iActivity] # Species 1 code            
            
            vnam='p1'; vc=GetColumn(vnam)
            sheet.cell(row=cnt+N_headers,column=vc).value=meta['Portfolio']['Activities']['New_Spc1_PCT'][iActivity] # Species 1 percent           
            
            vnam='i1'; vc=GetColumn(vnam)
            sheet.cell(row=cnt+N_headers,column=vc).value=meta['Portfolio']['Activities']['Site index (m)'][iActivity] # Species 1 site index            
            
            if isnan(meta['Portfolio']['Activities']['New_Spc2_CD'][iActivity])==False:
                vnam='s2'; vc=GetColumn(vnam)
                sheet.cell(row=cnt+N_headers,column=vc).value=meta['Portfolio']['Activities']['New_Spc2_CD'][iActivity] # Species 2 code             
            
                vnam='p2'; vc=GetColumn(vnam)
                sheet.cell(row=cnt+N_headers,column=vc).value=meta['Portfolio']['Activities']['New_Spc2_PCT'][iActivity] # Species 2 percent            
            
            if isnan(meta['Portfolio']['Activities']['New_Spc3_CD'][iActivity])==False:
                vnam='s3'; vc=GetColumn(vnam)
                sheet.cell(row=cnt+N_headers,column=vc).value=meta['Portfolio']['Activities']['New_Spc3_CD'][iActivity] # Species 3 code             
                
                vnam='p3'; vc=GetColumn(vnam)
                sheet.cell(row=cnt+N_headers,column=vc).value=meta['Portfolio']['Activities']['New_Spc3_PCT'][iActivity] # Species 3 percent
            
            if isnan(meta['Portfolio']['Activities']['New_Spc4_CD'][iActivity])==False:
                vnam='s4'; vc=GetColumn(vnam)
                sheet.cell(row=cnt+N_headers,column=vc).value=meta['Portfolio']['Activities']['New_Spc4_CD'][iActivity] # Species 3 code             
                
                vnam='p4'; vc=GetColumn(vnam)
                sheet.cell(row=cnt+N_headers,column=vc).value=meta['Portfolio']['Activities']['New_Spc4_PCT'][iActivity] # Species 3 percent
            
            vnam='init_density'; vc=GetColumn(vnam)            
            sheet.cell(row=cnt+N_headers,column=vc).value=int(meta['Portfolio']['Activities']['New initial density (SPH)'][iActivity]) # Planting density    
            
            vnam='regen_delay'; vc=GetColumn(vnam)
            sheet.cell(row=cnt+N_headers,column=vc).value=meta['Portfolio']['Activities']['New regeneration delay (years)'][iActivity] # Regeneration delay            
            
            vnam='oaf1'; vc=GetColumn(vnam)
            sheet.cell(row=cnt+N_headers,column=vc).value=meta['Portfolio']['Activities']['OAF1'][iActivity] # OAF1                        
            
            vnam='oaf2'; vc=GetColumn(vnam)
            sheet.cell(row=cnt+N_headers,column=vc).value=meta['Portfolio']['Activities']['OAF2'][iActivity] # OAF2
            
            vnam='bec_zone'; vc=GetColumn(vnam)
            sheet.cell(row=cnt+N_headers,column=vc).value=meta['Portfolio']['Activities']['BGC Zone'][iActivity] # BEC zone
            
            vnam='FIZ'; vc=GetColumn(vnam)
            if (meta['Portfolio']['Activities']['BGC Zone'][iActivity]=='CWH') | (meta['Portfolio']['Activities']['BGC Zone'][iActivity]=='CDF'):
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
    
    for iScn in range(meta['N Scenario']):
    
        # Loop through batches, saving inventory to file
        for iBat in range(meta['N Batch']):
      
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
                cd=meta['Portfolio']['Activities']['BGC Zone'][meta['Portfolio']['Indices']['Activity'][indBat[i]]]
                inv['ID_BECZ'][0,indBat[i]]=meta['LUT']['VRI']['BEC_ZONE_CODE'][cd]
    
            # Timber harvesting landbase (1=yes, 0=no)
            inv['THLB']=1*np.ones((meta['Year'].size,N_StandsInBatch))
        
            # Temperature will be updated automatically
            inv['MAT']=4*np.ones((1,N_StandsInBatch))
            
            if meta['Biomass Module']=='Sawtooth':
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
    
    for iScn in range(meta['N Scenario']):
        
        for iEns in range(meta['N Ensemble']):        
            
            for iBat in range(meta['N Batch']):
    
                # Index to batch
                indBat=cbu.IndexToBatch(meta,iBat)
                N_StandsInBatch=len(indBat)
    
                tv=np.arange(meta['Year Start'],meta['Year End']+1,1)
    
                # Initialize dictionary
                ec={}
                ec['ID_Type']=np.zeros((meta['Year'].size,indBat.size,meta['Max Events Per Year']),dtype='int16')
                ec['MortalityFactor']=np.zeros((meta['Year'].size,indBat.size,meta['Max Events Per Year']),dtype='int16')
                ec['GrowthFactor']=np.zeros((meta['Year'].size,indBat.size,meta['Max Events Per Year']),dtype='int16')
                ec['ID_GrowthCurve']=np.zeros((meta['Year'].size,indBat.size,meta['Max Events Per Year']),dtype='int16')
            
                #--------------------------------------------------------------
                # Simulate wildfire occurrence and severity from Taz
                #--------------------------------------------------------------
                    
#                if meta['Scenario'][iScn]['AAO Wildfire Status']=='On':
#                        
#                    # Import inventory to get BGC zone
#                    inv=gu.ipickle(meta['Paths']['Input Scenario'][iScn] + '\\Inventory_Bat' + FixFileNum(iBat) + '.pkl')
#                        
#                    # Prepare required parameters dictionary
#                    # Import fcgadgets parameters
#                    par=invu.Load_Params(meta)
#                    par['WF']['Scenario ID']=meta['Scenario'][iScn]['AAO Wildfire Scenario ID']
#                    par['WF']['Exclude simulations during modern period']='On'
#                    par['WF']['Exclude simulations during historical period']='Off'
#                    par['WF']['Exclude simulations during future period']='Off'
#                    
#                    # Think about moving this into the par dictionary
#                    method_occ='DirectFromParetoDraw'
#                    
#                    # Normally, this would be run here, but this approach assumes
#                    # that stands have been swapped with ensembles (when running
#                    # from spreadsheet). Instead, it is being re-run for each stand
#                    wf_sim=asm.GenerateWildfireEnsembleFromAAO(meta,par,inv['ID_BECZ'],method_occ)  
            
                for iS in range(N_StandsInBatch): 
                    
                    # Index to portfolio tables
                    iY=meta['Portfolio']['Indices']['Year'][indBat[iS]]
                    iA=meta['Portfolio']['Indices']['Activity'][indBat[iS]]
                    
                    #----------------------------------------------------------
                    # Add spinup events
                    #----------------------------------------------------------
                    
                    ivl_spin=meta['Spinup Disturbance Return Inverval']                
                    YearRef=meta['Portfolio']['Year'][iY]
                    AgeRef=meta['Portfolio']['Activity']['Mean stand age at time of inciting disturbance (years)'][iA]
                    if AgeRef>=0:
                        Year=np.arange(YearRef-AgeRef-100*ivl_spin,YearRef-AgeRef+ivl_spin,ivl_spin)        
                    else:
                        Year1=meta['Year Start']+ivl_spin
                        Year2=meta['Spinup Year End']
                        Year=np.arange(Year1,Year2+1,meta['Spinup Disturbance Return Inverval'])
                    
                    for iYr in range(Year.size):
                        iT=np.where(tv==Year[iYr])[0]
                        ec['ID_Type'][iT,iS,0]=meta['LUT']['Dist'][meta['Spinup Disturbance Type']]
                        ec['MortalityFactor'][iT,iS,0]=100
                        ec['GrowthFactor'][iT,iS,0]=0
                        ec['ID_GrowthCurve'][iT,iS,0]=meta['Spinup Growth Curve ID']
                    
                    #----------------------------------------------------------
                    # Add simulated constant disturbances
                    #----------------------------------------------------------
               
                    # Historical disturbance from simulation 1
                    ri=meta['Scenario'][iScn]['ReturnInterval1_Hist_DisFromSim']
                    if (ri!=0) & (np.isnan(ri)==False):
                    
                        p_Dist=1/ri
                        p_Rand=np.random.uniform(0,1,size=(meta['Year'].size))        
                        it=np.where((p_Rand<p_Dist) & (meta['Year']>meta['Spinup Year End']) & (meta['Year']<meta['Year Project']))[0]
                        Year=meta['Year'][it]                    
                        for iYr in range(Year.size):
                            iT=np.where(tv==Year[iYr])[0]
                            ec['ID_Type'][iT,iS,0]=meta['LUT']['Dist'][meta['Scenario'][iScn]['Type1_Hist_DisFromSim']]
                            ec['MortalityFactor'][iT,iS,0]=100
                            ec['GrowthFactor'][iT,iS,0]=0
                            ec['ID_GrowthCurve'][iT,iS,0]=meta['Spinup Growth Curve ID']
                    
                    # Historical disturbance from simulation 2
                    ri=meta['Scenario'][iScn]['ReturnInterval2_Hist_DisFromSim']
                    if (ri!=0) & (np.isnan(ri)==False):
                    
                        p_Dist=1/ri
                        p_Rand=np.random.uniform(0,1,size=(meta['Year'].size))        
                        it=np.where((p_Rand<p_Dist) & (meta['Year']>meta['Spinup Year End']) & (meta['Year']<meta['Year Project']))[0]
                        Year=meta['Year'][it]                    
                        for iYr in range(Year.size):
                            iT=np.where(tv==Year[iYr])[0]
                            ec['ID_Type'][iT,iS,0]=meta['LUT']['Dist'][meta['Scenario'][iScn]['Type2_Hist_DisFromSim']]
                            ec['MortalityFactor'][iT,iS,0]=100
                            ec['GrowthFactor'][iT,iS,0]=0
                            ec['ID_GrowthCurve'][iT,iS,0]=meta['Spinup Growth Curve ID']
      
                    #----------------------------------------------------------
                    # Add events from inventory
                    #----------------------------------------------------------
                    
                    for iYr in range(1,7):
                        
                        if np.isnan(meta['Scenario'][iScn]['Year' + str(iYr) + '_DisFromInv'])==True:
                            continue
            
                        # If IDW, convert IDW class to growth and mortality factor
                        sc=np.array(['IDW-T','IDW-L','IDW-M','IDW-S','IDW-V','IDW-MM','IDW-MS','IDW-MV','IDW-SS','IDW-SV','IDW-VV'])
                        flg_i=0
                        indSc=np.where(sc==meta['Scenario'][iScn]['Type' + str(iYr) + '_DisFromInv'])[0]
                        if indSc.size!=0:
                            if flg_i==0:
                                dfParDistBySC=pd.read_excel(meta['Paths']['Model Code'] + '\\Parameters\\Parameters_DisturbanceBySeverityClass.xlsx')
                                flg_i=1
                            indPar=np.where( (dfParDistBySC['Name']=='IDW') & (dfParDistBySC['SeverityCD']==sc[indSc[0]][4:]) )[0]
                            ID_TypeN=meta['LUT']['Dist']['IDW']
                            MF=dfParDistBySC.loc[indPar,'MortalityFactor']
                            GF=dfParDistBySC.loc[indPar,'GrowthFactor']
                        else:
                            ID_TypeS=meta['Scenario'][iScn]['Type' + str(iYr) + '_DisFromInv']
                            ID_TypeN=meta['LUT']['Dist'][ID_TypeS]
                            MF=meta['Scenario'][iScn]['Severity' + str(iYr) + '_DisFromInv']
                            GF=0
            
                        Year=meta['Scenario'][iScn]['Year' + str(iYr) + '_DisFromInv']
                        iT=np.where(tv==Year)[0]
                        
                        iE=np.where(ec['ID_Type'][iT,iS,:]==0)[1]

                        ec['ID_Type'][iT,iS,iE[0]]=ID_TypeN
                        ec['MortalityFactor'][iT,iS,iE[0]]=MF
                        ec['GrowthFactor'][iT,iS,iE[0]]=GF
                        ec['ID_GrowthCurve'][iT,iS,iE[0]]=meta['Scenario'][iScn]['GrowthCurve' + str(iYr) + '_DisFromInv']

                    #----------------------------------------------------------
                    # Add simulated constant future disturbances
                    #----------------------------------------------------------
               
                    # Future disturbance from simulation 1
                    ri=meta['Scenario'][iScn]['ReturnInterval1_Fut_DisFromSim']
                    if (ri!=0) & (np.isnan(ri)==False):
                    
                        p_Dist=1/ri
                        p_Rand=np.random.uniform(0,1,size=(meta['Year'].size))        
                        it=np.where((p_Rand<p_Dist) & (meta['Year']>meta['Year Project']))[0]
                        Year=meta['Year'][it]                    
                        for iYr in range(Year.size):
                            iT=np.where(tv==Year[iYr])[0]
                            ec['ID_Type'][iT,iS,0]=meta['LUT']['Dist'][meta['Scenario'][iScn]['Type1_Fut_DisFromSim']]
                            ec['MortalityFactor'][iT,iS,0]=100
                            ec['GrowthFactor'][iT,iS,0]=0
                            ec['ID_GrowthCurve'][iT,iS,0]=meta['Spinup Growth Curve ID']
                    
                    # Future disturbance from simulation 2
                    ri=meta['Scenario'][iScn]['ReturnInterval2_Fut_DisFromSim']
                    if (ri!=0) & (np.isnan(ri)==False):
                    
                        p_Dist=1/ri
                        p_Rand=np.random.uniform(0,1,size=(meta['Year'].size))        
                        it=np.where((p_Rand<p_Dist) & (meta['Year']>meta['Year Project']))[0]
                        Year=meta['Year'][it]                    
                        for iYr in range(Year.size):
                            iT=np.where(tv==Year[iYr])[0]
                            ec['ID_Type'][iT,iS,0]=meta['LUT']['Dist'][meta['Scenario'][iScn]['Type2_Fut_DisFromSim']]
                            ec['MortalityFactor'][iT,iS,0]=100
                            ec['GrowthFactor'][iT,iS,0]=0
                            ec['ID_GrowthCurve'][iT,iS,0]=meta['Spinup Growth Curve ID']    
            
                    #----------------------------------------------------------
                    # Add simulated wildfire from Taz
                    #----------------------------------------------------------
                    
                    if meta['Scenario'][iScn]['AAO Wildfire Status']=='On':
            
                        ind=np.where(wf_sim['Occurrence'][:,iS]==1)[0]
                        if ind.size==0:
                            continue
                        
                        ID_Type=meta['LUT']['Dist']['Wildfire']*np.ones(ind.size)
                        Year=tv[ind]
                        MortF=wf_sim['Mortality'][ind,iS]
                        GrowthF=0*np.ones(ind.size)
                        ID_GrowthCurve=1*np.ones(ind.size)
                            
                        for iYr in range(Year.size):
                            iT=np.where(tv==Year[iYr])[0]
                            ec['ID_Type'][iT,iS,0]=ID_Type[iYr]
                            ec['MortalityFactor'][iT,iS,0]=MortF[iYr]
                            ec['GrowthFactor'][iT,iS,0]=GrowthF[iYr]
                            ec['ID_GrowthCurve'][iT,iS,0]=ID_GrowthCurve[iYr]
    
                    #----------------------------------------------------------
                    # Add simulated MPB from Taz
                    #----------------------------------------------------------

            
                #--------------------------------------------------------------
                # Save to file            
                #--------------------------------------------------------------
                
                gu.opickle(meta['Paths']['Input Scenario'][iScn] + '\\Events_Ens' + cbu.FixFileNum(iEns) + '_Bat' + cbu.FixFileNum(iBat) + '.pkl',ec)
    
    return