# Import Python modules
import os
import pickle
import gc

from cbrun_annproc import *
from cbrun_export import *
from cbrun_utilities import *

def RunProject(meta):

    # Loop through scenarios
    for iScn in range(0,meta['N Scenario']): 
        
        # Loop through ensembles
        for iEns in range(0,meta['N Ensemble']):
        
            # Loop through batches
            for iBat in range(0,meta['N Batch']):
        
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
                    fout=open(pthWO,'wb')
                    pickle.dump([],fout); fout.close();
            
                print('Running Scenario ' + FixFileNum(iScn) + ', Ensemble ' + FixFileNum(iEns) + ', Batch ' + FixFileNum(iBat))
                
                # Initialize stands
                meta,vi,vo=InitializeStands(meta,iScn,iEns,iBat)
                
                # Import parameters
                vi,meta,ptl=ImportParameters(meta,vi)
                
                # Add stand-level parameters to metadata structure
                psl=Bunch(meta['psl'])
              
                # Biomass dynamics from Sawtooth
                if meta['Biomass Module']=='Sawtooth':
                    for iS in range(meta['N Stand']):
                        vo=BiomassFromSawtooth(iScn,iS,vi,vo,ptl,meta)
                                
                # Loop through time intervals (start in second time step)
                for iT in range(1,meta['N Time']):
            
                    if meta['Biomass Module']=='TIPSY':
                        vo=BiomassFromTIPSY(iScn,iT,vi,vo,psl,meta)
                    
                    # Calculate annual dead organic matter dynamics
                    vo=DeadOrganicMatterDynamics(iT,vi,vo,psl)
    
                    # Calculate prescribed disturbances
                    vo=Disturbances(iT,vi,vo,psl,meta)
        
                    # Calculate products sector
                    vo=ProductSector(iT,vi,vo,psl,meta)
        
                # Export simulation
                ExportVariables(meta,vi,vo,iScn,iEns,iBat,psl)
                
                # Delete working on file
                if meta['Skip Completed Runs']=='On':   
                    os.remove(pthWO)
                    
                # Delete variables
                del vi,vo
                gc.collect()

# List of pools
# Number 1           2              3         4        5      6            7          8          9         10        11        12         13           14      15       16      17           18     19     20    21   22    23          24             25
#Index:  0           1              2         3        4      5            6          7          8         9         10        11         12           13      14       15      16           17     18     19    20   21    22          23             24
#Name:  'StemMerch','StemNonMerch','Foliage','Branch','Bark','RootCoarse','RootFine','LitterVF','LitterF','LitterM','LitterS','SnagStem','SnagBranch','SoilVF','SoilF','SoilS','BlackCarbon','Peat','CO2','CH4','CO','N2O','MillMerch','MillNonMerch','MillStemSnag'

'''============================================================================
INITIALIZE STANDS
============================================================================'''

def InitializeStands(meta,iScn,iEns,iBat):
    
    #--------------------------------------------------------------------------
    # IMPORT INVENTORY, GROWTH CURVES, AND DISTURBANCE HISTORY
    #--------------------------------------------------------------------------
        
    # Import inventory
    fin=open(meta['Path Input Scenario'][iScn] + '\\Inventory_Bat' + FixFileNum(iBat) + '.pkl','rb')
    inv=pickle.load(fin); fin.close()            
            
    # Import disturbance history
    fin=open(meta['Path Input Scenario'][iScn] + '\\Disturbance_Ens' + FixFileNum(iEns) + '_Bat' + FixFileNum(iBat) + '.pkl','rb')
    dh=pickle.load(fin); fin.close()   
    
    # Import growth curves
    if meta['Biomass Module']=='TIPSY':
        
        # Import growth curve 1
        fin=open(meta['Path Input Scenario'][iScn] + '\\GrowthCurve1_Bat' + FixFileNum(iBat) + '.pkl','rb')
        gc1=pickle.load(fin); fin.close()
        gcA=gc1.copy()
    
        # Import growth curve 2
        fin=open(meta['Path Input Scenario'][iScn] + '\\GrowthCurve2_Bat' + FixFileNum(iBat) + '.pkl','rb')
        gc2=pickle.load(fin); fin.close()
        
        # Import growth curve 3
        fin=open(meta['Path Input Scenario'][iScn] + '\\GrowthCurve3_Bat' + FixFileNum(iBat) + '.pkl','rb')
        gc3=pickle.load(fin); fin.close()
        
    else:        
        gcA=0
        gc1=0
        gc2=0
        gc3=0
    
    # Calculate time vector
    tv=meta['Time']
          
    # Add inventory and growth curves to input variable structure
    vi=InputVariable(tv,inv,dh,gcA,gc1,gc2,gc3)
    
    # Update project configuration
    meta['N Time']=len(tv)
    meta['N Stand']=vi.inv.Lat.shape[1]       
    
    # Create a monotonic increase in different growth curves with each disturbance event    
    for j in range(meta['N Stand']):        
        n=vi.dh[j].ID_GrowthCurve.shape[0]
        vi.dh[j].ID_GrowthCurveM=1*np.ones((n,))        
        d=np.diff(vi.dh[j].ID_GrowthCurve)
        cnt=1
        for k in range(0,d.shape[0]):
            if d[k]!=0:
                cnt=cnt+1
                vi.dh[j].ID_GrowthCurveM[k+1]=cnt 
            else:
                vi.dh[j].ID_GrowthCurveM[k+1]=vi.dh[j].ID_GrowthCurveM[k]
    
    #--------------------------------------------------------------------------
    # INITIALIZE OUTPUT VARIABLES
    #--------------------------------------------------------------------------
        
    # Output variables   
    vo=OutputVariable
    
    # Stand age (i.e. time since stand-replacing disturbance)    
    vo.A=np.zeros((meta['N Time'],meta['N Stand']))
    
    # Stand density
    vo.N=np.zeros((meta['N Time'],meta['N Stand']))
    
    # Change in stand density (stems ha-1 yr-1)
    vo.N_R=np.zeros((meta['N Time'],meta['N Stand']))
    vo.N_M=np.zeros((meta['N Time'],meta['N Stand']))
    vo.N_M_Inv_Fir=np.zeros((meta['N Time'],meta['N Stand']))
    vo.N_M_Inv_Ins=np.zeros((meta['N Time'],meta['N Stand']))
    vo.N_M_Inv_Pat=np.zeros((meta['N Time'],meta['N Stand']))
    vo.N_M_Inv_Har=np.zeros((meta['N Time'],meta['N Stand']))
    vo.N_M_Inv_Win=np.zeros((meta['N Time'],meta['N Stand']))
    vo.N_M_Sim_Dir=np.zeros((meta['N Time'],meta['N Stand']))
    vo.N_M_Sim_Fir=np.zeros((meta['N Time'],meta['N Stand']))
    vo.N_M_Sim_Ins=np.zeros((meta['N Time'],meta['N Stand']))
    vo.N_M_Sim_Pat=np.zeros((meta['N Time'],meta['N Stand']))
    vo.N_M_Sim_Har=np.zeros((meta['N Time'],meta['N Stand']))
    vo.N_M_Sim_Win=np.zeros((meta['N Time'],meta['N Stand']))
    
    # Carbon density of ecosystem (Mg C ha-1)
    # -> 3-D matrix: Time x Stand x Carbon pool
    vo.C_Eco_Pools=np.zeros((meta['N Time'],meta['N Stand'],meta['N Pools Eco']))
        
    # Carbon flux densities (Mg C ha-1 yr-1)
    vo.C_NPP=np.zeros((meta['N Time'],meta['N Stand'],meta['N Pools Eco']))
    vo.C_G_Net=np.zeros((meta['N Time'],meta['N Stand'],meta['N Pools Eco']))
    vo.C_G_Gross=np.zeros((meta['N Time'],meta['N Stand'],meta['N Pools Eco']))
    vo.C_M=np.zeros((meta['N Time'],meta['N Stand'],meta['N Pools Eco']))
    vo.C_M_Inv_Fir=np.zeros((meta['N Time'],meta['N Stand']))
    vo.C_M_Inv_Ins=np.zeros((meta['N Time'],meta['N Stand']))
    vo.C_M_Inv_Pat=np.zeros((meta['N Time'],meta['N Stand']))
    vo.C_M_Inv_Har=np.zeros((meta['N Time'],meta['N Stand']))
    vo.C_M_Inv_Win=np.zeros((meta['N Time'],meta['N Stand']))
    vo.C_M_Sim_Dir=np.zeros((meta['N Time'],meta['N Stand']))
    vo.C_M_Sim_Fir=np.zeros((meta['N Time'],meta['N Stand']))
    vo.C_M_Sim_Ins=np.zeros((meta['N Time'],meta['N Stand']))
    vo.C_M_Sim_Pat=np.zeros((meta['N Time'],meta['N Stand']))
    vo.C_M_Sim_Har=np.zeros((meta['N Time'],meta['N Stand']))
    vo.C_M_Sim_Win=np.zeros((meta['N Time'],meta['N Stand']))     
    vo.C_LF=np.zeros((meta['N Time'],meta['N Stand'],meta['N Pools Eco']))    
    vo.C_RH=np.zeros((meta['N Time'],meta['N Stand'],meta['N Pools Eco']))
    vo.C_CombustCO2=np.zeros((meta['N Time'],meta['N Stand']))
    vo.C_CombustCH4=np.zeros((meta['N Time'],meta['N Stand']))
    vo.C_CombustCO=np.zeros((meta['N Time'],meta['N Stand']))
    vo.C_CombustN2O=np.zeros((meta['N Time'],meta['N Stand']))
    vo.C_E_Operations=np.zeros((meta['N Time'],meta['N Stand']))
    vo.C_RemovedMerch=np.zeros((meta['N Time'],meta['N Stand']))
    vo.C_RemovedNonMerch=np.zeros((meta['N Time'],meta['N Stand']))    
    vo.C_RemovedSnagStem=np.zeros((meta['N Time'],meta['N Stand']))    
    
    # Carbon density of products sector (Mg C ha-1)
    # -> 3-D matrix: Time x Stand x Carbon pool
    vo.C_Pro_Pools=np.zeros((meta['N Time'],meta['N Stand'],meta['N Pools Pro']))  
    
    # Summary of aboveground biomass
    vo.C_BiomassAG=np.zeros((meta['N Time'],meta['N Stand']))    
    
    # Mean of tree attributes
    vo.TreeMean_A=np.zeros((meta['N Time'],meta['N Stand']))
    vo.TreeMean_H=np.zeros((meta['N Time'],meta['N Stand']))
    vo.TreeMean_D=np.zeros((meta['N Time'],meta['N Stand']))
    vo.TreeMean_Csw=np.zeros((meta['N Time'],meta['N Stand']))
    vo.TreeMean_Csw_G=np.zeros((meta['N Time'],meta['N Stand']))
    
    vo.V_StemMerch=np.zeros((meta['N Time'],meta['N Stand']))
    
    return meta,vi,vo

'''============================================================================
IMPORT PARAMETERS
============================================================================'''

def ImportParameters(meta,vi):
      
    # Open connection to parameter database    
    fin=open(meta['Path Model Code'] + '\\Parameters\\Parameters.pkl','rb')
    par=pickle.load(fin)
    fin.close()
    
    # Initialize parameter structure
    #psl=Parameter
    psl={}
    
    # Initialize parameter structure
    ptl=Parameter
         
    #--------------------------------------------------------------------------
    # BIOMASS Allometry (STAND LEVEL)
    #--------------------------------------------------------------------------
    
    # Find index to species/BGC combination
    d=par['BiomassAllomSL']
    for i in range(len(d['Code_Spc1_PSP'].values())):
        spc=d['Code_Spc1_PSP'][i]
        bgc=d['Code_BGC_PSP'][i]
        if (spc=='All') & (bgc=='SBS'):
            break    
    psl['bAllo_StemToF1']=d['B_F1'][i]    
    psl['bAllo_StemToF2']=d['B_F2'][i]
    psl['bAllo_StemToBr1']=d['B_Br1'][i]
    psl['bAllo_StemToBr2']=d['B_Br2'][i]
    psl['bAllo_StemToBk1']=d['B_Bk1'][i]
    psl['bAllo_StemToBk2']=d['B_Bk2'][i]
    
    #--------------------------------------------------------------------------
    # BIOMASS TURNOVER
    #--------------------------------------------------------------------------
    
    # Define biomass pool names
    str_Pool=meta['Name Pools Eco'][0:7]
    
    # Initialize attributes       
    for i in range(len(str_Pool)):        
        Name='bTR_' + str_Pool[i]
        Value=np.zeros((1,meta['N Stand']))
        #setattr(psl,Name,Value) 
        psl[Name]=Value
    
    # Populate parameter structure from database (ID_Species,ID_Pool,Value)
    for i in range(len(str_Pool)):
        Name='bTR_' + str_Pool[i]
        Value=np.zeros((1,meta['N Stand']))
        for j in range(meta['N Stand']):             
            ind=np.where((par['BiomassTurnover']['ID_Species']==vi.inv.Spc1_ID[0,j]) & \
                              (par['BiomassTurnover']['ID_Pool']==i))[0]
            Value[0,j]=par['BiomassTurnover']['Value'][ind]
        #setattr(psl,Name,Value)
        psl[Name]=Value
    
    #--------------------------------------------------------------------------
    # BIOPHYSICAL
    #--------------------------------------------------------------------------
    
    for i in range(len(par['Biophysical']['Handle'])):
        #setattr(psl,'b' + par['Biophysical']['Handle'][i],par['Biophysical']['Value'][i])
        psl['b' + par['Biophysical']['Handle'][i]]=par['Biophysical']['Value'][i]
    
    #--------------------------------------------------------------------------
    # INTER-POOL TRANSFER FRACTIONS
    #--------------------------------------------------------------------------
    
    # Populate parameter structure
    for i in range(len(par['InterPoolFluxes']['Name'])):
        Name='bIPF_' + par['InterPoolFluxes']['Name'][i]
        Value=par['InterPoolFluxes']['Value'][i]
        #setattr(psl,Name,Value) 
        psl[Name]=Value
    
    #--------------------------------------------------------------------------
    # DECOMPOSITION
    #--------------------------------------------------------------------------
    
    for i in range(len(par['Decomposition']['Name_Pool'])):
        Name='bDec_' + par['Decomposition']['Name_Pool'][i] + '_Rten'
        Value=par['Decomposition']['Rten'][i]
        #setattr(psl,Name,Value) 
        psl[Name]=Value

        Name='bDec_' + par['Decomposition']['Name_Pool'][i] + '_Qten'
        Value=par['Decomposition']['Qten'][i]
        #setattr(psl,Name,Value) 
        psl[Name]=Value
        
        Name='bDec_' + par['Decomposition']['Name_Pool'][i] + '_PhysTransRate'
        Value=par['Decomposition']['PhysTransRate'][i]               
        #setattr(psl,Name,Value) 
        psl[Name]=Value
        
    #--------------------------------------------------------------------------
    # DISTURBANCE
    #--------------------------------------------------------------------------
    
    for i in par['Disturbances'].keys():
        #setattr(psl,'bDist_' + i,np.array(list(par['Disturbances'][i].values())))
        psl['bDist_' + i]=np.array(list(par['Disturbances'][i].values()))
     
    #--------------------------------------------------------------------------
    # HARVESTED WOOD PRODUCTS
    #--------------------------------------------------------------------------
        
    for i in range(len(par['HWP']['Name'])):
        
        Name='HWP_' + par['HWP']['Name'][i]
        Value=par['HWP']['Value'][i]
        
        # Convert half life to turnover rate
        if Name[-2:]=='hl':            
            Name=Name[0:-2] + 'tr'
            Value=1-np.exp(-np.log(2)/Value)        
        #setattr(psl,Name,Value)
        psl[Name]=Value
    
    if meta['Biomass Module']=='Sawtooth':
    
        #--------------------------------------------------------------------------
        # SAWTOOTH KEY BETWEEN PROVINCIAL SPECIES AND SPECIES-REGION-SAMPLES
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
        # POPULATE SPECIES BASED ON SPECIES_REGION_SAMPLES
        #--------------------------------------------------------------------------
    
        for i in range(vi.inv.Srs1_ID.shape[1]):
            # Species 1
            cd=ptl.KeySrsToSpc[vi.inv.Srs1_ID[0,i]-1]['Species_CD_BC']
            ind=np.where(par['SpeciesVRI']['SPECIES_CODE']==cd)[0]
            vi.inv.Spc1_ID[0,i]=par['SpeciesVRI']['ID_SPECIES'][ind]
        
            # Species 2
            cd=ptl.KeySrsToSpc[vi.inv.Srs2_ID[0,i]-1]['Species_CD_BC']
            if np.isnan(cd)==False:
                ind=np.where(par['SpeciesVRI']['SPECIES_CODE']==cd)[0]
                vi.inv.Spc2_ID[0,i]=par['SpeciesVRI']['ID_SPECIES'][ind]
        
            # Species 3
            cd=ptl.KeySrsToSpc[vi.inv.Srs3_ID[0,i]-1]['Species_CD_BC']
            if np.isnan(cd)==False:
                ind=np.where(par['SpeciesVRI']['SPECIES_CODE']==cd)[0]
                vi.inv.Spc3_ID[0,i]=par['SpeciesVRI']['ID_SPECIES'][ind]             
    
        #--------------------------------------------------------------------------
        # SAWTOOTH PARAMETERS
        #--------------------------------------------------------------------------
        
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
            ptl.Allom=Allom
            ptl.R_Type=[]
            ptl.R_Coef=[]
            ptl.R_Type.append('Def1')
            ptl.R_Coef.append(R_Def1)
            ptl.M_Type=[]
            ptl.M_Coef=[]
            ptl.M_Type.append('Def1')
            ptl.M_Coef.append(M_Def1)
            ptl.G_Type=[]
            ptl.G_Coef=[]
            ptl.G_Type.append('Def1')
            ptl.G_Coef.append(G_Def1)
    
    # Add stand=level parameter dictionary to meta dicationary
    meta['psl']=psl
    
    return vi,meta,ptl

'''============================================================================
DEFINE CLASSES
============================================================================'''

# Define inventory data structure
class InventoryData:
    def __init__(self,Lat,Lon,X,Y,Srs1_ID,Srs1_Pct,Srs2_ID,Srs2_Pct,Srs3_ID,Srs3_Pct, \
                 Spc1_ID,Spc1_Pct,Spc2_ID,Spc2_Pct,Spc3_ID,Spc3_Pct,ID_BECZ,MAT):
        self.Lat=Lat
        self.Lon=Lon
        self.X=X
        self.Y=Y
        self.Srs1_ID=Srs1_ID
        self.Srs1_Pct=Srs1_Pct
        self.Srs2_ID=Srs2_ID
        self.Srs2_Pct=Srs2_Pct
        self.Srs3_ID=Srs3_ID
        self.Srs3_Pct=Srs3_Pct
        self.Spc1_ID=Spc1_ID
        self.Spc1_Pct=Spc1_Pct
        self.Spc2_ID=Spc2_ID
        self.Spc2_Pct=Spc2_Pct
        self.Spc3_ID=Spc3_ID
        self.Spc3_Pct=Spc3_Pct
        self.ID_BECZ=ID_BECZ
        self.MAT=MAT  

# Define disturbance events structure
class Events:
    def __init__(self,Year,ID_Type,Severity,ID_GrowthCurve,Biomass_Affected_Pct, \
                 Biomass_Merch_Removed_Pct,Biomass_Merch_Burned_Pct,Biomass_Merch_LeftOnSite_Pct, \
                 Biomass_NonMerch_Removed_Pct,Biomass_NonMerch_Burned_Pct,Biomass_NonMerch_LeftOnSite_Pct, \
                 Snags_Affected_Pct,Snags_Removed_Pct,Snags_Burned_Pct,Snags_LeftOnSite_Pct, \
                 RemovedMerchToFuel_Pct,RemovedMerchToLumber_Pct,RemovedMerchToPlywood_Pct, \
                 RemovedMerchToOSB_Pct,RemovedMerchToMDF_Pct,RemovedMerchToPulp_Pct, \
                 RemovedMerchToFirewood_Pct,RemovedMerchToCants_Pct,RemovedNonMerchToFuel_Pct,RemovedNonMerchToLumber_Pct, \
                 RemovedNonMerchToPlywood_Pct,RemovedNonMerchToOSB_Pct,RemovedNonMerchToMDF_Pct, \
                 RemovedNonMerchToPulp_Pct,RemovedNonMerchToFirewood_Pct,RemovedNonMerchToCants_Pct,RemovedSnagStemToFuel_Pct, \
                 RemovedSnagStemToLumber_Pct,RemovedSnagStemToPlywood_Pct,RemovedSnagStemToOSB_Pct, \
                 RemovedSnagStemToMDF_Pct,RemovedSnagStemToPulp_Pct,RemovedSnagStemToFirewood_Pct,RemovedSnagStemToCants_Pct):
        self.Year=Year
        self.ID_Type=ID_Type
        self.Severity=Severity
        self.ID_GrowthCurve=ID_GrowthCurve
        self.Biomass_Affected_Pct=Biomass_Affected_Pct
        self.Biomass_Merch_Removed_Pct=Biomass_Merch_Removed_Pct
        self.Biomass_Merch_Burned_Pct=Biomass_Merch_Burned_Pct
        self.Biomass_Merch_LeftOnSite_Pct=Biomass_Merch_LeftOnSite_Pct
        self.Biomass_NonMerch_Removed_Pct=Biomass_NonMerch_Removed_Pct
        self.Biomass_NonMerch_Burned_Pct=Biomass_NonMerch_Burned_Pct
        self.Biomass_NonMerch_LeftOnSite_Pct=Biomass_NonMerch_LeftOnSite_Pct
        self.Snags_Affected_Pct=Snags_Affected_Pct
        self.Snags_Removed_Pct=Snags_Removed_Pct
        self.Snags_Burned_Pct=Snags_Burned_Pct
        self.Snags_LeftOnSite_Pct=Snags_LeftOnSite_Pct
        self.RemovedMerchToFuel_Pct=RemovedMerchToFuel_Pct
        self.RemovedMerchToLumber_Pct=RemovedMerchToLumber_Pct
        self.RemovedMerchToPlywood_Pct=RemovedMerchToPlywood_Pct
        self.RemovedMerchToOSB_Pct=RemovedMerchToOSB_Pct
        self.RemovedMerchToMDF_Pct=RemovedMerchToMDF_Pct
        self.RemovedMerchToPulp_Pct=RemovedMerchToPulp_Pct
        self.RemovedMerchToFirewood_Pct=RemovedMerchToFirewood_Pct
        self.RemovedMerchToCants_Pct=RemovedMerchToCants_Pct
        self.RemovedNonMerchToFuel_Pct=RemovedNonMerchToFuel_Pct
        self.RemovedNonMerchToLumber_Pct=RemovedNonMerchToLumber_Pct
        self.RemovedNonMerchToPlywood_Pct=RemovedNonMerchToPlywood_Pct
        self.RemovedNonMerchToOSB_Pct=RemovedNonMerchToOSB_Pct
        self.RemovedNonMerchToMDF_Pct=RemovedNonMerchToMDF_Pct
        self.RemovedNonMerchToPulp_Pct=RemovedNonMerchToPulp_Pct
        self.RemovedNonMerchToFirewood_Pct=RemovedNonMerchToFirewood_Pct
        self.RemovedNonMerchToCants_Pct=RemovedNonMerchToCants_Pct
        self.RemovedSnagStemToFuel_Pct=RemovedSnagStemToFuel_Pct
        self.RemovedSnagStemToLumber_Pct=RemovedSnagStemToLumber_Pct
        self.RemovedSnagStemToPlywood_Pct=RemovedSnagStemToPlywood_Pct
        self.RemovedSnagStemToOSB_Pct=RemovedSnagStemToOSB_Pct
        self.RemovedSnagStemToMDF_Pct=RemovedSnagStemToMDF_Pct
        self.RemovedSnagStemToPulp_Pct=RemovedSnagStemToPulp_Pct
        self.RemovedSnagStemToFirewood_Pct=RemovedSnagStemToFirewood_Pct
        self.RemovedSnagStemToCants_Pct=RemovedSnagStemToCants_Pct

# Define input variable structure
class InputVariable:
    def __init__(self,tv,inv,dh,GrowthCurveA,GrowthCurve1,GrowthCurve2,GrowthCurve3):
        self.tv=tv
        self.inv=inv        
        self.dh=dh
        self.GrowthCurveA=GrowthCurveA
        self.GrowthCurve1=GrowthCurve1
        self.GrowthCurve2=GrowthCurve2
        self.GrowthCurve3=GrowthCurve3

# Define output variable sturcture
class OutputVariable:
    def __init__(self, \
    A,N,N_R,N_M,N_M_Inv_Fir,N_M_Inv_Ins,N_M_Inv_Pat,N_M_Inv_Har,N_M_Inv_Win, \
    N_M_Sim_Dir,N_M_Sim_Fir,N_M_Sim_Ins,N_M_Sim_Pat,N_M_Sim_Har,N_M_Sim_Win, \
    C_NPP,C_G_Gross,C_G_Net,C_M_Inv_Fir,C_M_Inv_Ins,C_M_Inv_Pat,C_M_Inv_Har, \
    C_M_Inv_Win,C_M_Sim_Dir,C_M_Sim_Fir,C_M_Sim_Ins,C_M_Sim_Pat,C_M_Sim_Har, \
    C_M_Sim_Win,C_LF,C_RH,C_CombustCO2,C_CombustCH4,C_CombustCO,C_CombustN2O, \
    C_E_Operations, \
    C_RemovedMerch,C_RemovedNonMerch,C_RemovedSnagStem,C_BiomassAG,C_Eco_Pools, \
    C_Pro_Pools,TreeMean_A,TreeMean_H,TreeMean_D,TreeMean_Csw,TreeMean_Gsw,V_StemMerch):
        self.A=A        
        self.N=N
        self.N_R=N_R
        self.N_M_Inv_Fir=N_M_Inv_Fir
        self.N_M_Inv_Ins=N_M_Inv_Ins
        self.N_M_Inv_Pat=N_M_Inv_Pat
        self.N_M_Inv_Har=N_M_Inv_Har
        self.N_M_Inv_Win=N_M_Inv_Win
        self.N_M_Sim_Dir=N_M_Sim_Dir
        self.N_M_Sim_Fir=N_M_Sim_Fir
        self.N_M_Sim_Ins=N_M_Sim_Ins
        self.N_M_Sim_Pat=N_M_Sim_Pat
        self.N_M_Sim_Har=N_M_Sim_Har
        self.N_M_Sim_Win=N_M_Sim_Win     
        self.C_NPP=C_NPP
        self.C_G_Gross=C_G_Gross
        self.C_G_Net=C_G_Net
        self.C_M_Inv_Fir=C_M_Inv_Fir
        self.C_M_Inv_Ins=C_M_Inv_Ins
        self.C_M_Inv_Pat=C_M_Inv_Pat
        self.C_M_Inv_Har=C_M_Inv_Har
        self.C_M_Inv_Win=C_M_Inv_Win
        self.C_M_Sim_Dir=C_M_Sim_Dir
        self.C_M_Sim_Fir=C_M_Sim_Fir
        self.C_M_Sim_Ins=C_M_Sim_Ins
        self.C_M_Sim_Pat=C_M_Sim_Pat
        self.C_M_Sim_Har=C_M_Sim_Har
        self.C_M_Sim_Win=C_M_Sim_Win        
        self.C_LF=C_LF                
        self.C_RH=C_RH        
        self.C_CombustCH4=C_CombustCH4
        self.C_CombustCO=C_CombustCO
        self.C_CombustN2O=C_CombustN2O
        self.C_E_Operations=C_E_Operations
        self.C_RemovedMerch=C_RemovedMerch
        self.C_RemovedNonMerch=C_RemovedNonMerch
        self.C_RemovedSnagStem=C_RemovedSnagStem
        self.C_BiomassAG=C_BiomassAG        
        self.C_Eco_Pools=C_Eco_Pools
        self.C_Pro_Pools=C_Pro_Pools
        self.TreeMean_A=TreeMean_A
        self.TreeMean_H=TreeMean_H
        self.TreeMean_D=TreeMean_D
        self.TreeMean_Csw=TreeMean_Csw
        self.TreeMean_Csw_G=TreeMean_Csw_G
        self.V_StemMerch=V_StemMerch

# Define tree level structure
class TreeVariable:
    def __init__(self, \
    p,ID_Srs,ID_Decid,A,N_R,N_M_Sim_Dir,N_M_Sim_Fir,N_M_Sim_Ins,N_M_Sim_Pat, \
    N_M_Sim_Har,N_M_Sim_Win,N_M_Inv_Dir,N_M_Inv_Fir,N_M_Inv_Ins,N_M_Inv_Pat, \
    N_M_Inv_Har,N_M_Inv_Win,H,D,Csw,Csw_Larger,Csw_G):
        self.p=p
        self.ID_Species=ID_Species
        self.ID_Decid=ID_Decid
        self.A=A
        self.N_R=N_R
        self.N_M_Sim_Dir=N_M_Sim_Dir
        self.N_M_Sim_Fir=N_M_Sim_Fir
        self.N_M_Sim_Ins=N_M_Sim_Ins
        self.N_M_Sim_Pat=N_M_Sim_Pat
        self.N_M_Sim_Har=N_M_Sim_Har
        self.N_M_Sim_Win=N_M_Sim_Win
        self.N_M_Inv_Dir=N_M_Inv_Dir
        self.N_M_Inv_Fir=N_M_Inv_Fir
        self.N_M_Inv_Ins=N_M_Inv_Ins
        self.N_M_Inv_Pat=N_M_Inv_Pat
        self.N_M_Inv_Har=N_M_Inv_Har
        self.N_M_Inv_Win=N_M_Inv_Win
        self.H=H
        self.D=D
        self.Csw=Csw
        self.Csw_Larger=Csw_Larger
        self.Csw_G=Csw_G

# Define parameter class
class Parameter(object):
    pass

