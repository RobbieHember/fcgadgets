
import os
import pickle
import gc
import time
import datetime
import numpy as np
import matplotlib.pyplot as plt

from fcgadgets.pyscripts import utilities_general as gu
from fcgadgets.cbrunner.cbrun_utilities import *
from fcgadgets.cbrunner.cbrun_annproc import *

'''============================================================================
RUN PROJECT
============================================================================'''

def RunProject(meta):

    # Loop through scenarios
    for iScn in range(0,meta['N Scenario']):
        
        # Loop through ensembles
        for iEns in range(0,meta['N Ensemble']):
        
            # Loop through batches
            for iBat in range(0,meta['N Batch']):
                
                t0=time.time()
                #iScn=0
                #iEns=0
                #iBat=0
                
                # Path to output data file
                pthAC=meta['Path Output Scenario'][iScn] + '\\Data_Scn' + FixFileNum(iScn) + '_Ens' + FixFileNum(iEns) + '_Bat' + FixFileNum(iBat) + '.pkl'
                
                # Path to temporary "working on" file -> tells other instances to back off!
                pthWO=meta['Path Output Scenario'][iScn] + '\\WorkingOn_Scn' + FixFileNum(iScn) + '_Ens' + FixFileNum(iEns) + '_Bat' + FixFileNum(iBat) + '.pkl'
                
                # Only proceed if the file does not exist running multiple instances
                if meta['Skip Completed Runs']=='On':                    
                    
                    if os.path.exists(pthAC):
                        continue                    
                    if os.path.exists(pthWO):
                        continue        
                    
                    # Create a WorkingOn file to let other instances know that this batch is
                    # in progress
                    gu.opickle(pthWO,[])                    
            
                print('Running Scenario ' + FixFileNum(iScn) + ', Ensemble ' + FixFileNum(iEns) + ', Batch ' + FixFileNum(iBat))
                
                # Initialize stands
                meta,vi,vo=InitializeStands(meta,iScn,iEns,iBat)
                
                # Import parameters
                vi,meta=ImportParameters(meta,vi)
                
                # Extract stand-level parameters from metadata structure
                psl=Bunch(meta['psl'])
              
                # Indices to ecosystem pools
                iEP=meta['iEP']
                
                # Biomass dynamics from Sawtooth
                if meta['Biomass Module']=='Sawtooth':                    
                    for iS in range(meta['N Stand']):                        
                        vo=BiomassFromSawtooth(iScn,iS,vi,vo,meta,iEP)
                                
                # Loop through time intervals (start in second time step)
                for iT in range(1,meta['N Time']):            
                    
                    # Biomass dynamics
                    if meta['Biomass Module']=='TIPSY':
                        
                        # FCI standard approach:
                        vo=BiomassFromBatchTIPSY(iScn,iT,vi,vo,psl,meta,iEP)
                    
                    elif meta['Biomass Module']=='TIPSY_SpecialAdjustments_EP703':
                        
                        # *** Special version of module writen for EP703 study ***
                        vo=BiomassFromBatchTIPSY_SpecialAdjustments_EP703(iScn,iT,vi,vo,psl,meta,iEP)
                    
                    # Calculate annual dead organic matter dynamics
                    vo=DOMFromCBM08(iT,vi,vo,psl,iEP)
                    
                    # Calculate prescribed disturbances
                    vo=Taz(iT,vi,vo,psl,meta,iEP)
                    
                    # Calculate products sector                    
                    if meta['Time'][iT]>=1800:
                        
                        # No need to run this before a certain date
                        vo=HWPFromDymond12(iT,vi,vo,psl,meta)
                
                # Export simulation results to file
                ExportSimulation(meta,vi,vo,iScn,iEns,iBat,psl,iEP)
                
                # Delete 'working on' file
                if meta['Skip Completed Runs']=='On':   
                    os.remove(pthWO)
                    
                # Delete variables
                del vi,vo
                gc.collect()
                
                # Track simulation time
                t1=time.time()
                print(t1-t0)
                
'''============================================================================
INITIALIZE STANDS
============================================================================'''

#------------------------------------------------------------------------------
# List of pools
#------------------------------------------------------------------------------
#Index:  0           1              2         3        4      5            6          7                 8                    9               10             11           12         13        14        15        16         17           18       19      20      21        22        23       24        25          26             27
#Name:  'StemMerch','StemNonMerch','Foliage','Branch','Bark','RootCoarse','RootFine','FelledStemMerch','FelledStemNonMerch','FelledFoliage','FelledBranch','FelledBark','LitterVF','LitterF','LitterM','LitterS','SnagStem','SnagBranch','SoilVF','SoilF','SoilS','ECO2asC','ECH4asC','ECOasC','EN2OasC','MillMerch','MillNonMerch','MillStemSnag'

def InitializeStands(meta,iScn,iEns,iBat):
    
    #--------------------------------------------------------------------------
    # Import input variables
    #--------------------------------------------------------------------------

    # Input variables dictionary
    vi={}
        
    vi['tv']=meta['Time']
    tv=meta['Time']
    
    # Import inventory
    vi['Inv']=gu.ipickle(meta['Path Input Scenario'][iScn] + '\\Inventory_Bat' + FixFileNum(iBat) + '.pkl')
            
    # Import disturbance history
    vi['DH']=gu.ipickle(meta['Path Input Scenario'][iScn] + '\\Disturbance_Ens' + FixFileNum(iEns) + '_Bat' + FixFileNum(iBat) + '.pkl')
    
    # Import growth curves
    if meta['Biomass Module']=='TIPSY':        
        
        # Import growth curve 1
        vi['GC1']=gu.ipickle(meta['Path Input Scenario'][iScn] + '\\GrowthCurve1_Bat' + FixFileNum(iBat) + '.pkl')
        
        # Set active growth curve to growth curve 1
        vi['GCA']=vi['GC1'].copy()
        
        # Initialize an indicator of the active growth curve ID
        vi['ID_GCA']=np.ones(meta['N Stand'])
        
        # Import growth curve 2
        vi['GC2']=gu.ipickle(meta['Path Input Scenario'][iScn] + '\\GrowthCurve2_Bat' + FixFileNum(iBat) + '.pkl')
        
        # Import growth curve 3
        vi['GC3']=gu.ipickle(meta['Path Input Scenario'][iScn] + '\\GrowthCurve3_Bat' + FixFileNum(iBat) + '.pkl')
        
    else:        
        
        vi['GCA']=0
        vi['GC1']=0
        vi['GC2']=0
        vi['GC3']=0
    
    #--------------------------------------------------------------------------
    # Ensure simultaneous occurences of harvest and slashpile burning are input
    # in a realistic order
    #--------------------------------------------------------------------------

    for iDH in range(len(vi['DH'])):
        u=np.unique(vi['DH'][iDH]['Year'])
        for iu in range(u.size):
            iUY=np.where(vi['DH'][iDH]['Year']==u[iu])[0]
            if iUY.size>1:
                iH=np.where(vi['DH'][iDH]['ID_Type'][iUY]==meta['LUT Dist']['Harvest'])[0]
                iB=np.where(vi['DH'][iDH]['ID_Type'][iUY]==meta['LUT Dist']['Slashpile Burn'])[0]
                if (iH.size>0) & (iB.size>0):
                    if iB<iH:
                        # Unrealistic sequence of slashpile burn and harvest, fix order                  
                        ID_Type_h=vi['DH'][iDH]['ID_Type'][iUY][iH]
                        Year_h=vi['DH'][iDH]['Year'][iUY][iH]
                        Severity_h=vi['DH'][iDH]['Severity'][iUY][iH]
                        ID_GrowthCurve_h=vi['DH'][iDH]['ID_GrowthCurve'][iUY][iH]
                        
                        ID_Type_b=vi['DH'][iDH]['ID_Type'][iUY][iB]
                        Year_b=vi['DH'][iDH]['Year'][iUY][iB]
                        Severity_b=vi['DH'][iDH]['Severity'][iUY][iB]
                        ID_GrowthCurve_b=vi['DH'][iDH]['ID_GrowthCurve'][iUY][iB]
                        
                        vi['DH'][iDH]['ID_Type'][iUY][iH]=ID_Type_b
                        vi['DH'][iDH]['Year'][iUY][iH]=Year_b
                        vi['DH'][iDH]['Severity'][iUY][iH]=Severity_b
                        vi['DH'][iDH]['ID_GrowthCurve'][iUY][iH]=ID_GrowthCurve_b
                        
                        vi['DH'][iDH]['ID_Type'][iUY][iB]=ID_Type_h
                        vi['DH'][iDH]['Year'][iUY][iB]=Year_h
                        vi['DH'][iDH]['Severity'][iUY][iB]=Severity_h
                        vi['DH'][iDH]['ID_GrowthCurve'][iUY][iB]=ID_GrowthCurve_h    

    #--------------------------------------------------------------------------
    # Disturbance year is stored as a float, but it needs to be an integer once
    # in the model
    #--------------------------------------------------------------------------

    for iDH in range(len(vi['DH'])):
        vi['DH'][iDH]['Year']=vi['DH'][iDH]['Year'].astype(int)

    #--------------------------------------------------------------------------
    # Update project configuration
    #--------------------------------------------------------------------------
    
    meta['N Time']=len(tv)
    meta['N Stand']=vi['Inv']['Lat'].shape[1]       
    
    #--------------------------------------------------------------------------
    # Create a monotonic increase in different growth curves with each disturbance 
    # event
    #--------------------------------------------------------------------------
    
    for j in range(meta['N Stand']):        
        try:
            n=vi['DH'][j]['ID_GrowthCurve'].shape[0]
        except:
            print(iBat)
            print(j)
        vi['DH'][j]['ID_GrowthCurveM']=1*np.ones((n,))        
        d=np.diff(vi['DH'][j]['ID_GrowthCurve'])
        cnt=1
        for k in range(0,d.shape[0]):
            if d[k]!=0:
                cnt=cnt+1
                vi['DH'][j]['ID_GrowthCurveM'][k+1]=cnt 
            else:
                vi['DH'][j]['ID_GrowthCurveM'][k+1]=vi['DH'][j]['ID_GrowthCurveM'][k]
    
    #--------------------------------------------------------------------------
    # Identify the start and end of a fixed spin up disturbance interval
    # This is used to fast-track through the spinup period.
    #--------------------------------------------------------------------------
    
    if meta['Biomass Module']=='Sawtooth':        
        d=np.append(0,np.diff(vi['DH'][0]['Year']))
        ind=np.where(d==meta['Spinup Disturbance Return Inverval'])[0]
        meta['SpinupSpanFastTrack']=[np.min(vi['DH'][0]['Year'][ind])+meta['Spinup Disturbance Return Inverval'],np.max(vi['DH'][0]['Year'][ind])]
    
    #--------------------------------------------------------------------------
    # Initialize flag for fixing negative net growth. When TIPSY yields negative
    # net growth, the fluxes of gross growth and mortality need adjustment. 
    # This flag helps achieve that.
    #--------------------------------------------------------------------------
    
    meta['FlagNegNetGrowth']=np.zeros(meta['N Stand'])
    meta['G_Net_PriorToBreakup']=np.zeros((meta['N Stand'],7))
    
    #--------------------------------------------------------------------------
    # Initialize output variables
    #--------------------------------------------------------------------------
        
    # Output variables dictionary   
    vo={}
    
    # Get dimensions    
    m=meta['N Time']
    n=meta['N Stand']
    o=meta['N Pools Eco']
    
    # Stand age (i.e. time since stand-replacing disturbance)    
    vo['A']=np.zeros((m,n))
        
    # Carbon density of ecosystem (Mg C ha-1)
    # -> 3-D matrix: Time x Stand x Carbon pool
    vo['C_Eco_Pools']=np.zeros((m,n,o))
        
    # Carbon flux densities (Mg C ha-1 yr-1)
    vo['C_NPP']=np.zeros((m,n,o))
    vo['C_G_Net']=np.zeros((m,n,o))
    vo['C_G_Gross']=np.zeros((m,n,o))
    vo['C_M_Reg']=np.zeros((m,n,o))
    vo['C_M_Inv_Fir']=np.zeros((m,n))
    vo['C_M_Inv_Ins']=np.zeros((m,n))
    vo['C_M_Inv_Pat']=np.zeros((m,n))
    vo['C_M_Inv_Har']=np.zeros((m,n))
    vo['C_M_Inv_Win']=np.zeros((m,n))
    vo['C_M_Sim_Reg']=np.zeros((m,n))
    vo['C_M_Sim_Fir']=np.zeros((m,n))
    vo['C_M_Sim_Ins']=np.zeros((m,n))
    vo['C_M_Sim_Pat']=np.zeros((m,n))
    vo['C_M_Sim_Har']=np.zeros((m,n))
    vo['C_M_Sim_Win']=np.zeros((m,n))     
    vo['C_LF']=np.zeros((m,n,o))    
    vo['C_RH']=np.zeros((m,n,o))
    vo['C_E_FireAsCO2']=np.zeros((m,n))
    vo['C_E_FireAsCH4']=np.zeros((m,n))
    vo['C_E_FireAsCO']=np.zeros((m,n))
    vo['C_E_FireAsN2O']=np.zeros((m,n))
    vo['C_E_OperationsAsCO2']=np.zeros((m,n))
    vo['C_RemovedMerch']=np.zeros((m,n))
    vo['C_RemovedNonMerch']=np.zeros((m,n))
    vo['C_RemovedSnagStem']=np.zeros((m,n))
    
    # Carbon density of products sector (Mg C ha-1)
    # -> 3-D matrix: Time x Stand x Carbon pool
    vo['C_Pro_Pools']=np.zeros((m,n,meta['N Pools Pro']))  
    
    # Summary of aboveground biomass
    vo['C_BiomassAG']=np.zeros((m,n))    
    
    # Stand density
    vo['N']=np.zeros((m,n))
    
    # Change in stand density (stems ha-1 yr-1)
    vo['N_R']=np.zeros((m,n))
    vo['N_M_Tot']=np.zeros((m,n))
    vo['N_M_Reg']=np.zeros((m,n))
    vo['N_M_Inv_Fir']=np.zeros((m,n))
    vo['N_M_Inv_Ins']=np.zeros((m,n))
    vo['N_M_Inv_Pat']=np.zeros((m,n))
    vo['N_M_Inv_Har']=np.zeros((m,n))
    vo['N_M_Inv_Win']=np.zeros((m,n))
    vo['N_M_Sim_Reg']=np.zeros((m,n))
    vo['N_M_Sim_Fir']=np.zeros((m,n))
    vo['N_M_Sim_Ins']=np.zeros((m,n))
    vo['N_M_Sim_Pat']=np.zeros((m,n))
    vo['N_M_Sim_Har']=np.zeros((m,n))
    vo['N_M_Sim_Win']=np.zeros((m,n))
    
    # Mean of tree attributes
    vo['TreeMean_A']=np.zeros((m,n))
    vo['TreeMean_H']=np.zeros((m,n))
    vo['TreeMean_D']=np.zeros((m,n))
    vo['TreeMean_Csw']=np.zeros((m,n))
    vo['TreeMean_Csw_G']=np.zeros((m,n))
    
    vo['V_StemMerch']=np.zeros((m,n))
    
    #--------------------------------------------------------------------------
    # Initialize a log book that will record various diagnostics, warnings, 
    # and error flags
    #--------------------------------------------------------------------------
    
    meta['Logbook']=list()
    
    return meta,vi,vo

'''============================================================================
IMPORT PARAMETERS
============================================================================'''

def ImportParameters(meta,vi):
      
    #--------------------------------------------------------------------------
    # Open parameter database (updated manually by running the UpdateParameters
    # module)
    #--------------------------------------------------------------------------
    
    par=gu.ipickle(meta['Path Model Code'] + '\\Parameters\\Parameters.pkl')
    
    #--------------------------------------------------------------------------
    # Initialize parameter structures
    #--------------------------------------------------------------------------
    
    psl={}
    ptl={}
    #ptl=Parameter
         
    #--------------------------------------------------------------------------
    # Biomass allometry paramaters (stand level)
    #--------------------------------------------------------------------------
    
    # Find index to species/BGC combination
    d=par['BiomassAllomSL']
    for i in range(len(d['Code_Spc1_PSP'].values())):
        spc=d['Code_Spc1_PSP'][i]
        bgc=d['Code_BGC_PSP'][i]
        if (spc=='All') & (bgc=='SBS'):
            break    
    psl['bASL_StemToF1']=d['B_F1'][i]    
    psl['bASL_StemToF2']=d['B_F2'][i]
    psl['bASL_MerchBarkFrac']=d['MerchBarkFrac'][i]
    
    #--------------------------------------------------------------------------
    # Biomass turnover paramaters
    #--------------------------------------------------------------------------
    
    # Define biomass pool names
    str_Pool=meta['Name Pools Eco'][0:7]
    
    # Initialize attributes       
    for i in range(len(str_Pool)):        
        Name='bTR_' + str_Pool[i]
        Value=np.zeros((1,meta['N Stand']))
        psl[Name]=Value
    
    # Populate parameter structure from database (ID_Species,ID_Pool,Value)
    for i in range(len(str_Pool)):
        Name='bTR_' + str_Pool[i]
        Value=np.zeros((1,meta['N Stand']))
        for j in range(meta['N Stand']):             
            ind=np.where((par['BiomassTurnover']['ID_Species']==vi['Inv']['Spc1_ID'][0,j]) & \
                              (par['BiomassTurnover']['ID_Pool']==i))[0]
            Value[0,j]=par['BiomassTurnover']['Value'][ind]
        psl[Name]=Value
    
    #--------------------------------------------------------------------------
    # Biophysical paramaters
    #--------------------------------------------------------------------------
    
    for i in range(len(par['Biophysical']['Handle'])):
        psl['b' + par['Biophysical']['Handle'][i]]=par['Biophysical']['Value'][i]
    
    #--------------------------------------------------------------------------
    # Inter-pool transfer parameters
    #--------------------------------------------------------------------------
    
    # Populate parameter structure
    for i in range(len(par['InterPoolFluxes']['Name'])):        
        Name='bIPF_' + par['InterPoolFluxes']['Name'][i]
        Value=par['InterPoolFluxes']['Value'][i]
        psl[Name]=Value    
    
    #--------------------------------------------------------------------------
    # Decomposition parameters
    #--------------------------------------------------------------------------
    
    for i in range(len(par['Decomposition']['Name_Pool'])):
        Name='bDec_' + par['Decomposition']['Name_Pool'][i] + '_PhysTransRate'
        Value=par['Decomposition']['PhysTransRate'][i]               
        psl[Name]=Value
        Name='bDec_' + par['Decomposition']['Name_Pool'][i] + '_Rten'
        Value=par['Decomposition']['Rten'][i]        
        psl[Name]=Value
        Name='bDec_' + par['Decomposition']['Name_Pool'][i] + '_Qten'
        Value=par['Decomposition']['Qten'][i]
        psl[Name]=Value        
        
    #--------------------------------------------------------------------------
    # Disturbance parameters
    #--------------------------------------------------------------------------
    
    for i in par['Disturbances'].keys():        
        psl['bDist_' + i]=np.array(list(par['Disturbances'][i].values()))
     
    #--------------------------------------------------------------------------
    # Harvest wood products parameters
    #--------------------------------------------------------------------------
        
    for i in range(len(par['HWP']['Name'])):
        
        Name='HWP_' + par['HWP']['Name'][i]
        Value=par['HWP']['Value'][i]
        
        # Convert half life to turnover rate
        if Name[-2:]=='hl':            
            Name=Name[0:-2] + 'tr'
            Value=1-np.exp(-np.log(2)/Value)        
        psl[Name]=Value
    
    if meta['Biomass Module']=='Sawtooth':
    
        #--------------------------------------------------------------------------
        # Sawtooth key between provincial species codes adn species-region sample 
        # codes
        #--------------------------------------------------------------------------
    
        cn=list(par['SRS'].keys())
        KeySrsToSpc=list()
        for i in range(par['SRS']['SRS_CD'].size):
            dic={}
            for j in range(len(cn)):
                Name=cn[j]
                Value=par['SRS'][cn[j]][i]
                dic.update({Name:Value})
            KeySrsToSpc.append(dic)    
        ptl.KeySrsToSpc=KeySrsToSpc    
    
        #--------------------------------------------------------------------------
        # Populate species based on speices-region sample code
        #--------------------------------------------------------------------------
    
        for i in range(vi['Inv']['Srs1_ID'].shape[1]):
            # Species 1
            cd=ptl.KeySrsToSpc[vi['Inv']['Srs1_ID'][0,i]-1]['Species_CD_BC']
            ind=np.where(par['SpeciesVRI']['SPECIES_CODE']==cd)[0]
            vi['Inv']['Spc1_ID'][0,i]=par['SpeciesVRI']['ID_SPECIES'][ind]
        
            # Species 2
            cd=ptl.KeySrsToSpc[vi['Inv']['Srs2_ID'][0,i]-1]['Species_CD_BC']
            if np.isnan(cd)==False:
                ind=np.where(par['SpeciesVRI']['SPECIES_CODE']==cd)[0]
                vi['Inv']['Spc2_ID'][0,i]=par['SpeciesVRI']['ID_SPECIES'][ind]
        
            # Species 3
            cd=ptl.KeySrsToSpc[vi['Inv']['Srs3_ID'][0,i]-1]['Species_CD_BC']
            if np.isnan(cd)==False:
                ind=np.where(par['SpeciesVRI']['SPECIES_CODE']==cd)[0]
                vi['Inv']['Spc3_ID'][0,i]=par['SpeciesVRI']['ID_SPECIES'][ind]             
    
        #----------------------------------------------------------------------
        # Sawtooth parameters
        #----------------------------------------------------------------------
        
        if meta['Biomass Module']=='Sawtooth':
            
            tab='TreeAllometry'
            cn=list(par[tab].keys())
            Allom=list()
            for i in range(par[tab]['SRS_CD'].size):
                dic={}
                for j in range(len(cn)):
                    Name=cn[j]
                    Value=par[tab][cn[j]][i]
                    dic.update({Name:Value})
                Allom.append(dic)
        
            tab='R_Def1'
            cn=list(par[tab].keys())
            R_Def1=list()
            for i in range(par[tab]['SRS_CD'].size):
                dic={}
                for j in range(len(cn)):
                    Name=cn[j]
                    Value=par[tab][cn[j]][i]
                    dic.update({Name:Value})
                R_Def1.append(dic) 
        
            tab='M_Def1'
            cn=list(par[tab].keys())
            M_Def1=list()
            for i in range(par[tab]['SRS_CD'].size):
                dic={}
                for j in range(len(cn)):
                    Name=cn[j]
                    Value=par[tab][cn[j]][i]
                    dic.update({Name:Value})
                M_Def1.append(dic) 
    
            tab='G_Def1'
            cn=list(par[tab].keys())
            G_Def1=list()
            for i in range(par[tab]['SRS_CD'].size):
                dic={}
                for j in range(len(cn)):
                    Name=cn[j]
                    Value=par[tab][cn[j]][i]
                    dic.update({Name:Value})
                G_Def1.append(dic)     
    
            # Add to tree-level parameter structure
            ptl['Allom']=Allom
            ptl['R_Type']='Def1'
            ptl['R_Coef']=R_Def1
            ptl['M_Type']='Def1'
            ptl['M_Coef']=M_Def1
            ptl['G_Type']='Def1'
            ptl['G_Coef']=G_Def1
            #ptl.Allom=Allom
            #ptl.R_Type=[]
            #ptl.R_Coef=[]
            #ptl.R_Type.append('Def1')
            #ptl.R_Coef.append(R_Def1)
            #ptl.M_Type=[]
            #ptl.M_Coef=[]
            #ptl.M_Type.append('Def1')
            #ptl.M_Coef.append(M_Def1)
            #ptl.G_Type=[]
            #ptl.G_Coef=[]
            #ptl.G_Type.append('Def1')
            #ptl.G_Coef.append(G_Def1)
    
    #--------------------------------------------------------------------------
    # Add to meta dictionary
    #--------------------------------------------------------------------------
    
    # Add stand=level parameter dictionary to meta dicationary
    meta['psl']=psl
    meta['ptl']=ptl
    
    #return vi,meta,ptl
    return vi,meta

'''============================================================================
EXPORT SIMULATION RESULTS
============================================================================'''

def ExportSimulation(meta,vi,vo,iScn,iEns,iBat,psl,iEP):
      
    #--------------------------------------------------------------------------
    # Save project metadata
    #--------------------------------------------------------------------------
    
    if (iBat==0) & (iScn==0):
        gu.opickle(meta['Path Output Scenario'][iScn] + '\\Metadata.pkl',meta)
        
    #--------------------------------------------------------------------------
    # Save project data
    #--------------------------------------------------------------------------
    
    # Define time period to save to file
    it=np.where(vi['tv']>=meta['Year Start Saving'])[0]
       
    # Define scale factor
    ScaleFactor=meta['Scale Factor Export']
    
    for i in vo.keys():
        vo[i]=ScaleFactor*vo[i]
        vo[i]=vo[i].astype(int)
    Year=ScaleFactor*vi['tv'][it]
    Year=Year.astype(int)
    
    # Add data to dictionary
    if meta['Save Biomass Pools']=='Yes':
        
        dat={'ScaleFactor':ScaleFactor,
             'Year':Year,
             'A':vo['A'][it,:],
             'C_Eco_Pools':vo['C_Eco_Pools'][it,:,:],
             'C_Pro_Pools':vo['C_Pro_Pools'][it,:,:],
             'C_G_Gross':vo['C_G_Gross'][it,:,:],
             'C_G_Net':vo['C_G_Net'][it,:,:],         
             'C_M_Reg':vo['C_M_Reg'][it,:,:],
             'C_LF':vo['C_LF'][it,:,:],
             'C_NPP':vo['C_NPP'][it,:,:],
             'C_RH':vo['C_RH'][it,:,:],
             'C_M_Inv_Fir':vo['C_M_Inv_Fir'][it,:],
             'C_M_Inv_Har':vo['C_M_Inv_Har'][it,:],
             'C_M_Inv_Ins':vo['C_M_Inv_Ins'][it,:],
             'C_M_Inv_Pat':vo['C_M_Inv_Pat'][it,:],
             'C_M_Inv_Win':vo['C_M_Inv_Win'][it,:],
             'C_RemovedMerch':vo['C_RemovedMerch'][it,:],
             'C_RemovedNonMerch':vo['C_RemovedNonMerch'][it,:],
             'C_RemovedSnagStem':vo['C_RemovedSnagStem'][it,:],
             'C_E_FireAsCH4':vo['C_E_FireAsCH4'][it,:],
             'C_E_FireAsCO':vo['C_E_FireAsCO'][it,:],
             'C_E_FireAsCO2':vo['C_E_FireAsCO2'][it,:],
             'C_E_FireAsN2O':vo['C_E_FireAsN2O'][it,:],
             'C_E_OperationsAsCO2':vo['C_E_OperationsAsCO2'][it,:],
             'V_StemMerch':vo['V_StemMerch'][it,:]}
    
    elif meta['Save Biomass Pools']=='No':
          
        dat={'ScaleFactor':ScaleFactor,
             'Year':Year,
             'A':vo['A'][it,:],
             'C_Eco_Pools':np.sum(vo['C_Eco_Pools'][it,:,:],axis=2),
             'C_Pro_Pools':np.sum(vo['C_Pro_Pools'][it,:,:],axis=2),
             'C_G_Gross':np.sum(vo['C_G_Gross'][it,:,:],axis=2),
             'C_G_Net':np.sum(vo['C_G_Net'][it,:,:],axis=2),
             'C_M_Reg':np.sum(vo['C_M_Reg'][it,:,:],axis=2),
             'C_LF':np.sum(vo['C_LF'][it,:,:],axis=2),
             'C_NPP':np.sum(vo['C_NPP'][it,:,:],axis=2),
             'C_RH':np.sum(vo['C_RH'][it,:,:],axis=2),
             'C_M_Inv_Fir':vo['C_M_Inv_Fir'][it,:],
             'C_M_Inv_Har':vo['C_M_Inv_Har'][it,:],
             'C_M_Inv_Ins':vo['C_M_Inv_Ins'][it,:],
             'C_M_Inv_Pat':vo['C_M_Inv_Pat'][it,:],
             'C_M_Inv_Win':vo['C_M_Inv_Win'][it,:],
             'C_RemovedMerch':vo['C_RemovedMerch'][it,:],
             'C_RemovedNonMerch':vo['C_RemovedNonMerch'][it,:],
             'C_RemovedSnagStem':vo['C_RemovedSnagStem'][it,:],
             'C_E_FireAsCH4':vo['C_E_FireAsCH4'][it,:],
             'C_E_FireAsCO':vo['C_E_FireAsCO'][it,:],
             'C_E_FireAsCO2':vo['C_E_FireAsCO2'][it,:],
             'C_E_FireAsN2O':vo['C_E_FireAsN2O'][it,:],
             'C_E_OperationsAsCO2':vo['C_E_OperationsAsCO2'][it,:],
             'V_StemMerch':vo['V_StemMerch'][it,:]} 
        
        # Add pool categories
        dat['Eco_Biomass']=np.sum(vo['C_Eco_Pools'][:,:,iEP['BiomassTotal']],axis=2)[it,:]
        dat['Eco_BiomassAG']=np.sum(vo['C_Eco_Pools'][:,:,iEP['BiomassAboveground']],axis=2)[it,:]
        dat['Eco_Litter']=np.sum(vo['C_Eco_Pools'][:,:,iEP['Litter']],axis=2)[it,:]
        dat['Eco_DeadWood']=np.sum(vo['C_Eco_Pools'][:,:,iEP['DeadWood']],axis=2)[it,:]
        dat['Eco_Soil']=np.sum(vo['C_Eco_Pools'][:,:,iEP['Soil']],axis=2)[it,:]
        dat['Eco_Total']=dat['Eco_Biomass']+dat['Eco_Litter']+dat['Eco_DeadWood']+dat['Eco_Soil']
        dat['Pro_InUse']=np.sum(vo['C_Pro_Pools'][:,:,0:10],axis=2)[it,:]
        dat['Pro_DumpLandfill']=np.sum(vo['C_Pro_Pools'][:,:,10:17],axis=2)[it,:]

        # Combine emissions from product sector, express as tCO2e/ha/yr
        co2_to_c=meta['psl']['bRatio_CO2_to_C']
        gwp_co2=1
        gwp_ch4=meta['psl']['bGWP_CH4_AR5']
        f_co2=vo['C_Pro_Pools'][it,:,17]
        f_ch4=vo['C_Pro_Pools'][it,:,18]
        dat['Pro_Emissions_co2e']=gwp_co2*co2_to_c*f_co2+gwp_ch4*co2_to_c*f_ch4
        
    # Sawtooth variables
    if meta['Biomass Module']=='Sawtooth':
        dat['N']=vo['N'][it,:]
        dat['N_M']=vo['N_M_Tot'][it,:]
        dat['N_R']=vo['N_R'][it,:]
        dat['C_M_Sim_Reg']=vo['C_M_Sim_Reg'][it,:]
        dat['C_M_Sim_Fir']=vo['C_M_Sim_Fir'][it,:]
        dat['C_M_Sim_Har']=vo['C_M_Sim_Har'][it,:]
        dat['C_M_Sim_Ins']=vo['C_M_Sim_Ins'][it,:]
        dat['C_M_Sim_Pat']=vo['C_M_Sim_Pat'][it,:]
        dat['C_M_Sim_Win']=vo['C_M_Sim_Win'][it,:]
        dat['TreeMean_A']=vo['TreeMean_A'][it,:]
        dat['TreeMean_Csw']=vo['TreeMean_Csw'][it,:]
        dat['TreeMean_Csw_G']=vo['TreeMean_Csw_G'][it,:]
        dat['TreeMean_D']=vo['TreeMean_D'][it,:]
        dat['TreeMean_H']=vo['TreeMean_H'][it,:]
    
    # Save
    gu.opickle(meta['Path Output Scenario'][iScn] + '\\Data_Scn' + FixFileNum(iScn) + '_Ens' + FixFileNum(iEns) + '_Bat' + FixFileNum(iBat) + '.pkl',dat)
    
    #--------------------------------------------------------------------------
    # Save diagnostics
    #--------------------------------------------------------------------------
    
    if (iBat==0) & (iScn==0):
    
        cnams=['Value','Variable']
        df=pd.DataFrame(columns=cnams)
    
        # Model version
        df.loc[0]=[meta['Path Model Code'],'Version']    
    
        # Date
        now=datetime.datetime.now()
        df.loc[1]=[now,'Run date']
    
        # Simulation time
        #meta['t_Sim']=time.time()-meta['t_SimStart']
        #df.loc[2]=[meta['t_Sim']/60,'Simulation time (min)']
    
        # Check conservation of mass (stock change = NEBP)
        # *** only done for the first stand ***
        NPP=np.sum(vo['C_NPP'][:,0,0:7],axis=1)
        RH=np.sum(vo['C_RH'][:,0,7:16],axis=1)
        E=vo['C_E_FireAsCO2'][:,0]+vo['C_E_FireAsCO'][:,0]+vo['C_E_FireAsCH4'][:,0]+vo['C_E_FireAsN2O'][:,0]
        R=vo['C_RemovedMerch'][:,0]+vo['C_RemovedNonMerch'][:,0]
        x=np.sum(vo['C_Eco_Pools'][:,0,0:16],axis=1)
        y=np.sum(vo['C_Eco_Pools'][0,0,0:16],axis=0)+np.cumsum(NPP-RH-R-E)
        D_abs=np.mean(np.abs(y-x))
        D_rel=np.mean(np.abs(y-x)/np.maximum(0.000001,x)*100)
        df.loc[3]=[D_abs,'Mean absolute difference between stock change and NECB (MgC/ha)']
        df.loc[4]=[D_rel,'Mean relative difference between stock change and NECB (%)']    
        #plt.close(6)
        #plt.figure(6)
        #plt.plot(x)
        #plt.plot(y)
    
        # Year that slow soil carbon pool reaches dynamic equilibrium        
        def runningMean(x, N):
            y=np.zeros((len(x),))
            for ctr in range(len(x)):
                y[ctr]=np.sum(x[ctr:(ctr+N)])
            return y/N
    
        #plt.plot(runningMean(vo['C_Eco_Pools[:,0,15],150))
        #x=np.diff(vo['C_Eco_Pools'][:,0,15])/vo['C_Eco_Pools'][0:-1,0,15]*100
        #rmx=runningMean(x,150)
        #plt.figure(33)
        #plt.plot(rmx)
        #iDE=np.min(np.where(rmx<0.001))
        #yrDE=vi['tv'][iDE]
        #tDE=meta.N_t-iDE
        #df.loc[5]=[yrDE,'Equilibrium in SoilS reached at (year)']
        #df.loc[6]=[tDE,'Equlibrium in Soil S reached at (time before end of simulation)']
        
        # Save
        pthoutD=meta['Path Output Scenario'][iScn] + '\\Diagnostics.xlsx'
        writer=pd.ExcelWriter(pthoutD)
        df.to_excel(writer,'Sheet1')
        writer.save()
    
    return

