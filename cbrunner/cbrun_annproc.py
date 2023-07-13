
#%% Import python modules

import numpy as np
from fcgadgets.macgyver import utilities_general as gu
from fcgadgets.cbrunner import cbrun_utilities as cbu
from fcgadgets.hardhat import nutrient_application as napp
from fcgadgets.taz import aspatial_stat_models as asm
import warnings
import time
from scipy import stats

#%% Biomass dynamics

def Biomass_FromGYM(meta,pNam,iScn,iBat,iT,vi,vo,iEP):

    # Update stand age
    vo['A'][iT,:]=vo['A'][iT-1,:]+1

    #--------------------------------------------------------------------------
    # Net growth of aboveground biomass
    #--------------------------------------------------------------------------

    # Index to growth curves at current age of stand
    iAge=np.minimum(vo['A'][iT,:],meta['Modules']['GYM']['BatchTIPSY Maximum Age'])-1

    # Convert to integer
    iAge=iAge.astype(int)

    # Extract net growth
    # Notes: I tried to make this faster using unravel but it was way slower
    NetGrowth=np.zeros( (meta[pNam]['Project']['Batch Size'][iBat],6) )
    for iS in range(meta[pNam]['Project']['Batch Size'][iBat]):
        NetGrowth[iS,:]=vi['GC']['Active'][iAge[iS],iS,:].copy()

    # Growth enhancement
    if 'Growth Enhancement Status' in meta[pNam]['Scenario'][iScn]:
        if meta[pNam]['Scenario'][iScn]['Growth Enhancement Status']=='On':

            # if 'Growth Enhancement' not in meta[pNam]['Project']:
            #     #tv=np.arange(1800,2151,1)
            #     tv=meta[pNam]['Year'].copy()
            #     fG_PreI=0.5
            #     fG_PosI=1.0
            #     fG=fG_PreI*np.ones(tv.size)
            #     it1=np.where(tv==1860)[0]
            #     it2=np.where(tv==2010)[0]
            #     ind=np.arange(it1,it2+1,dtype='int16')

            #     fG[ind]=np.linspace(fG_PreI,fG_PosI,ind.size)
            #     fG[it2[0]:]=fG_PosI
            #     #plt.plot(tv,fG,'b-')
            #     meta[pNam]['Project']['Growth Enhancement']=fG

            NetGrowth=np.tile(meta[pNam]['Project']['Growth Enhancement'][iT,:],(NetGrowth.shape[1],1)).T*NetGrowth

    # # *** Special EvalAtPlots With GaiaSL Modifier***
    # t0=1950
    # t1=2020
    # m0=0.5
    # m1_young=2.0
    # m1_old=0.6
    # Age0=75
    # Age1=25
    # m1=gu.Clamp(m0+(Age1-Age0)/(Age1-Age0)*(vo['A'][iT,:]-Age0),m1_old,m1_young)
    # fG=gu.Clamp(m0+(m1-m0)/(t1-t0)*(meta[pNam]['Year'][iT]-t0),m0,m1)
    # NetGrowth=np.tile(fG,(NetGrowth.shape[1],1)).T*NetGrowth
    # # *** Special ***

    # # *** Special EvalAtPlots With Age Modifier***
    # b=-0.1
    # ymn=0.45
    # ymx=1.7
    # A_inf=50
    # #xhat=np.arange(1,301,1)
    # fG=ymn+(ymx-ymn)*1/(1+np.exp(-b*(vo['A'][iT,:]-A_inf)))
    # NetGrowth=np.tile(fG,(NetGrowth.shape[1],1)).T*NetGrowth
    # # *** Special ***

    # Net growth of total stemwood
    Gnet_Stem=NetGrowth[:,iEP['StemMerch']]+NetGrowth[:,iEP['StemNonMerch']]

    # Net growth of foliage
    NetGrowth[:,iEP['Foliage']]=Gnet_Stem*(meta['Param']['BEV']['Biomass Allometry']['Gf1']+(meta['Param']['BEV']['Biomass Allometry']['Gf2']-meta['Param']['BEV']['Biomass Allometry']['Gf1'])*np.exp(-meta['Param']['BEV']['Biomass Allometry']['Gf3']*vo['A'][iT,:]))

    # Net growth of branches
    NetGrowth[:,iEP['Branch']]=Gnet_Stem*(meta['Param']['BEV']['Biomass Allometry']['Gbr1']+(meta['Param']['BEV']['Biomass Allometry']['Gbr2']-meta['Param']['BEV']['Biomass Allometry']['Gbr1'])*np.exp(-meta['Param']['BEV']['Biomass Allometry']['Gbr3']*vo['A'][iT,:]))

    # Net growth of bark
    NetGrowth[:,iEP['Bark']]=Gnet_Stem*(meta['Param']['BEV']['Biomass Allometry']['Gbk1']+(meta['Param']['BEV']['Biomass Allometry']['Gbk2']-meta['Param']['BEV']['Biomass Allometry']['Gbk1'])*np.exp(-meta['Param']['BEV']['Biomass Allometry']['Gbk3']*vo['A'][iT,:]))

    # Add net growth to output variable structure
    # Oddly, using meta['iEP']['BiomassAboveground'] will invert the dimensions
    # of C_G_Net - don't change it.
    vo['C_G_Net'][iT,:,0:5]=NetGrowth[:,0:5]

    # Total net growth of root biomass (Li et al. 2003, Eq. 4)
    G_Net_Root_Total=0.22*np.sum(vo['C_G_Net'][iT,:,0:5],axis=1)

    # Fine root fraction should decline with total biomass, but estimating it based on
    # size causes a lot of problems. Think about creating new equation.
    vo['C_G_Net'][iT,:,iEP['RootCoarse']]=(1-0.072)*G_Net_Root_Total
    vo['C_G_Net'][iT,:,iEP['RootFine']]=0.072*G_Net_Root_Total

    # Update live stemwood merchantable volume
    vo['V_MerchLive'][iT,:]=vo['V_MerchLive'][iT-1,:]+NetGrowth[:,5]#.astype(float)*meta['Modules']['GYM']['Scale Factor']

    # Update dead stemwood merchantable volume
    vo['V_MerchDead'][iT,:]=vo['V_MerchDead'][iT-1,:]

    # Update total (live+dead) stemwood merchantable volume
    vo['V_MerchTotal'][iT,:]=vo['V_MerchLive'][iT,:]+vo['V_MerchDead'][iT,:]

    #--------------------------------------------------------------------------
    # Nutrient application effects to net growth
    #--------------------------------------------------------------------------

    # Index to stands that are stimulated by nutrient application
    meta['Modules']['Nutrient Management']['iApplication']=np.where( meta['Modules']['Nutrient Management']['ResponseCounter']>0 )[0]

    # Adjust N application response counter
    if meta['Modules']['Nutrient Management']['iApplication'].size>0:
        vi,vo,meta=napp.UpdateStatus(vi,vo,iT,meta,'UpdateCounter')

    # Adjust root net growth
    if meta['Modules']['Nutrient Management']['iApplication'].size>0:
        vi,vo,meta=napp.UpdateStatus(vi,vo,iT,meta,'BelowgroundNetGrowth')

    #--------------------------------------------------------------------------
    # Add net growth to biomass pools
    #--------------------------------------------------------------------------

    # Stemwood
    vo['C_Eco_Pools'][iT,:,iEP['StemMerch']]=np.maximum(0,vo['C_Eco_Pools'][iT-1,:,iEP['StemMerch']]+vo['C_G_Net'][iT,:,iEP['StemMerch']])
    vo['C_Eco_Pools'][iT,:,iEP['StemNonMerch']]=np.maximum(0,vo['C_Eco_Pools'][iT-1,:,iEP['StemNonMerch']]+vo['C_G_Net'][iT,:,iEP['StemNonMerch']])

    # Foliage
    # *** Notes 1: TIPSY foliage appears to be really low, use alternative equations
    # instead. As net growth of foliage has been redefined, we must recalculate
    # net growth. Adjustment is applied in utilities script. ***
    vo['C_Eco_Pools'][iT,:,iEP['Foliage']]=np.maximum(0,vo['C_Eco_Pools'][iT-1,:,iEP['Foliage']]+vo['C_G_Net'][iT,:,iEP['Foliage']])

    # Branches
    vo['C_Eco_Pools'][iT,:,iEP['Branch']]=np.maximum(0,vo['C_Eco_Pools'][iT-1,:,iEP['Branch']]+vo['C_G_Net'][iT,:,iEP['Branch']])

    # Bark
    vo['C_Eco_Pools'][iT,:,iEP['Bark']]=np.maximum(0,vo['C_Eco_Pools'][iT-1,:,iEP['Bark']]+vo['C_G_Net'][iT,:,iEP['Bark']])

    # Coarse roots
    vo['C_Eco_Pools'][iT,:,iEP['RootCoarse']]=np.maximum(0,vo['C_Eco_Pools'][iT-1,:,iEP['RootCoarse']]+vo['C_G_Net'][iT,:,iEP['RootCoarse']])

    # Fine roots
    vo['C_Eco_Pools'][iT,:,iEP['RootFine']]=np.maximum(0,vo['C_Eco_Pools'][iT-1,:,iEP['RootFine']]+vo['C_G_Net'][iT,:,iEP['RootFine']])

    # # Old root process based on size instead of growth - caused many problems
    # flg=0
    # if flg==1:

    #     # Calculate aboveground biomass
    #     vo['C_BiomassAG'][iT,:]=np.sum(vo['C_Eco_Pools'][iT,:,0:5],axis=1)

    #     # Total root biomass (Li et al. 2003, Eq. 4)
    #     BiomassRootTotal=0.222*vo['C_BiomassAG'][iT,:]

    #     # Fine root biomass (Li et al. 2003, Eq. 6)
    #     Pf=0.072+0.354*np.exp(-0.06*(2*BiomassRootTotal))
    #     vo['C_Eco_Pools'][iT,:,iEP['RootFine']]=Pf*BiomassRootTotal

    #     # Coarse root biomass
    #     Pc=1-Pf
    #     vo['C_Eco_Pools'][iT,:,iEP['RootCoarse']]=np.maximum(0,Pc*BiomassRootTotal)

    #     # Calculate net growth of roots from change in pools
    #     #Gnet_RC=vo['C_Eco_Pools'][iT,:,iEP['RootCoarse']]-vo['C_Eco_Pools'][iT-1,:,iEP['RootCoarse']]
    #     #Gnet_RF=vo['C_Eco_Pools'][iT,:,iEP['RootFine']]-vo['C_Eco_Pools'][iT-1,:,iEP['RootFine']]

    #--------------------------------------------------------------------------
    # Biomass turnover
    #--------------------------------------------------------------------------

    # Calculate biomass loss due to mortality
    vo['C_M_Reg'][iT,:,0]=meta['Param']['BEV']['Biomass Turnover']['StemMerch']*vo['C_Eco_Pools'][iT,:,0]
    vo['C_M_Reg'][iT,:,1]=meta['Param']['BEV']['Biomass Turnover']['StemMerch']*vo['C_Eco_Pools'][iT,:,1]
    vo['C_M_Reg'][iT,:,2]=meta['Param']['BEV']['Biomass Turnover']['StemMerch']*vo['C_Eco_Pools'][iT,:,2]
    vo['C_M_Reg'][iT,:,3]=meta['Param']['BEV']['Biomass Turnover']['StemMerch']*vo['C_Eco_Pools'][iT,:,3]
    vo['C_M_Reg'][iT,:,4]=meta['Param']['BEV']['Biomass Turnover']['StemMerch']*vo['C_Eco_Pools'][iT,:,4]
    vo['C_M_Reg'][iT,:,5]=meta['Param']['BEV']['Biomass Turnover']['StemMerch']*vo['C_Eco_Pools'][iT,:,5]
    vo['C_M_Reg'][iT,:,6]=meta['Param']['BEV']['Biomass Turnover']['StemMerch']*vo['C_Eco_Pools'][iT,:,6]

    # Adjust mortality to account for N application response
    if meta['Modules']['Nutrient Management']['iApplication'].size>0:
        vi,vo,meta=napp.UpdateStatus(vi,vo,iT,meta,'Mortality')

    # Apply optional factor to mortality
    #if meta['Scenario Switch']['Mortality Factor Status'][iScn]=='On':
    #    mortality_factor=meta[pNam]['Scenario'][iScn]['Mortality Factor'][iT]
    #    if mortality_factor!=0:
    #        # Normal operation
    #        vo['C_M_Reg'][iT,:,:]=net_growth_factor*vo['C_M_Reg'][iT,:,:]

    #--------------------------------------------------------------------------
    # If TIPSY indicates negative net growth (e.g., alder break-up), revise so
    # that net growth during the period of negative values equals the value in the
    # previous timestep.
    #--------------------------------------------------------------------------

    # Keep track of what regular mortality would be before it is affected by
    # correction for catastrophic mortality
    C_M_Reg=vo['C_M_Reg'][iT,:,:].copy()

    flg=0
    if flg==1:

        # Define a threshold level of negative growth so it isn't triggered by noise
        NegGrowthThreshold=-0.25
        #NegGrowthThreshold=-1e6

        # Find stands with negative net growth
        iNegNetG=np.where(vo['C_G_Net'][iT,:,0]<NegGrowthThreshold)[0]

        if iNegNetG.size>0:

            # If it is the first instance of negative net growth: 91) record net growth
            # of the preceeding timestep and (2) set flag = 1.
            iSwitchFlag=np.where(meta[pNam]['Project']['FlagNegNetGrowth'][iNegNetG]==0)[0]
            meta[pNam]['Project']['FlagNegNetGrowth'][iNegNetG[iSwitchFlag]]=1
            meta[pNam]['Project']['G_Net_PriorToBreakup'][iNegNetG[iSwitchFlag],0:7]=vo['C_G_Net'][iT-1,iNegNetG[iSwitchFlag],0:7]

            d=vo['C_G_Net'][iT,iNegNetG,0:7]-meta[pNam]['Project']['G_Net_PriorToBreakup'][iNegNetG,:]

            CToTransfer=np.zeros((meta[pNam]['Project']['Batch Size'][iBat],7))
            CToTransfer[iNegNetG,:]=-1*d

            vo['C_G_Net'][iT,:,0:7]=vo['C_G_Net'][iT,:,0:7]+CToTransfer
            vo['C_M_Reg'][iT,:,0:7]=vo['C_M_Reg'][iT,:,0:7]+CToTransfer

            # # Logbook entry
            # for i in range(iNegNetG.size):
            #     txt='Scenario:' + str(iScn) + ', Stand:' + str(iNegNetG[i]) + ', Time:' + str(iT) + ', Negative net growth received from GY model is being adjusted.'
            #     meta['Logbook'].append(txt)

    #--------------------------------------------------------------------------
    # Litterfall
    #--------------------------------------------------------------------------

    # Setting turnover as a function of age will decouple NPP from net growth.
    Aref=meta['Param']['BEV']['Biomass Turnover']['BiomassTurnoverAgeRef']
    fA=-meta['Param']['BEV']['Biomass Turnover']['BiomassTurnoverAgeDependence']

    # Calculate foliage biomass turnover due to litterfall
    tr=np.maximum(0.001,fA*(vo['A'][iT,:]-Aref)+meta['Param']['BEV']['Biomass Turnover']['Foliage'])
    vo['C_LF'][iT,:,iEP['Foliage']]=tr*vo['C_Eco_Pools'][iT,:,iEP['Foliage']]

    # Calculate branch biomass turnover due to litterfall
    tr=np.maximum(0.001,fA*(vo['A'][iT,:]-Aref)+meta['Param']['BEV']['Biomass Turnover']['Branch'])
    vo['C_LF'][iT,:,iEP['Branch']]=tr*vo['C_Eco_Pools'][iT,:,iEP['Branch']]

    # Calculate bark biomass turnover due to litterfall
    tr=np.maximum(0.001,fA*(vo['A'][iT,:]-Aref)+meta['Param']['BEV']['Biomass Turnover']['Bark'])
    vo['C_LF'][iT,:,iEP['Bark']]=meta['Param']['BEV']['Biomass Turnover']['Bark']*vo['C_Eco_Pools'][iT,:,iEP['Bark']]

    # Calculate coarse root biomass turnover due to litterfall
    tr=np.maximum(0.001,fA*(vo['A'][iT,:]-Aref)+meta['Param']['BEV']['Biomass Turnover']['RootCoarse'])
    vo['C_LF'][iT,:,iEP['RootCoarse']]=tr*vo['C_Eco_Pools'][iT,:,iEP['RootCoarse']]

    # Calculate fine root biomass turnover due to litterfall
    tr=np.maximum(0.001,fA*(vo['A'][iT,:]-Aref)+meta['Param']['BEV']['Biomass Turnover']['RootFine'])
    vo['C_LF'][iT,:,iEP['RootFine']]=tr*vo['C_Eco_Pools'][iT,:,iEP['RootFine']]

    # Adjust litterfall to account for N application response
    if meta['Modules']['Nutrient Management']['iApplication'].size>0:
        vi,vo,meta=napp.UpdateStatus(vi,vo,iT,meta,'Litterfall')

    #--------------------------------------------------------------------------
    # Update summary variables
    #--------------------------------------------------------------------------

    # Update gross growth
    vo['C_G_Gross'][iT,:,:]=vo['C_G_Net'][iT,:,:]+C_M_Reg

    # Update NPP
    #vo['C_NPP'][iT,:,:]=vo['C_G_Net'][iT,:,:]+C_M_Reg+vo['C_LF'][iT,:,:]

    return vo

#%% Annual tree biomass dynamics from Sawtooth

def BiomassFromSawtooth(meta,pNam,iScn,iS,vi,vo,iEP):

    t0=time.time()

    # Turn warnings off to avoid divide by zero warning
    warnings.filterwarnings('ignore')

    # *************************************************************************
    # What if we only track the current and previous time step - memory useage would
    # be low enough to vectorize stands
    # Requires constant population of stand level summary variables at the end
    # of each time step.
    #**************************************************************************

    # Define demensions
    N_time=meta[pNam]['Project']['N Time']

    # Maximum number of trees per stand
    N_tree=meta['Param']['BE']['Sawtooth']['Core']['Max SPH'].astype(int)

    # Multiplier to convert kg C to Mg C
    cm=1000

    #--------------------------------------------------------------------------
    # Initialize tree-level variables
    #--------------------------------------------------------------------------

    # List of variables
    vL=['ID_SRS','ID_Decid','A','H','D','N_R','Csw','Csw_Larger','Csw_G',
        'N_M_Reg','N_M_Fir','N_M_Ins','N_M_Pat','N_M_Har','N_M_Win',
        'Csw_M_Reg','Csw_M_Fir','Csw_M_Ins','Csw_M_Pat','Csw_M_Har','Csw_M_Win']

    # Populate
    tl={}
    for v in vL:
        if v[0:2]=='ID':
            tl[v]=np.zeros((N_time,N_tree),dtype=np.int)
        else:
            tl[v]=np.zeros((N_time,N_tree),dtype=np.float)

    #--------------------------------------------------------------------------
    # Populate species ID
    #--------------------------------------------------------------------------

    # Generate random number vector
    rp=np.random.permutation(N_tree).astype(int)

    # Species 1
    n1=np.ceil(vi['Inv']['SRS1_PCT'][0,iS]/100*N_tree).astype(int)
    tl['ID_SRS'][:,rp[0:n1]]=vi['Inv']['SRS1_ID'][0,iS]

    # Species 2
    if vi['Inv']['SRS2_PCT'][0,iS]>0:
        n2=np.ceil(vi['Inv']['SRS2_Pct'][0,iS]/100*N_tree).astype(int)
        tl['ID_SRS'][:,rp[n1:n2]]=vi['Inv']['SRS2_ID'][0,iS]

    # Species 3
    if vi['Inv']['SRS3_PCT'][0,iS]>0:
        tl['ID_SRS'][:,rp[n2:]]=vi['Inv']['SRS3_ID'][0,iS]

    #--------------------------------------------------------------------------
    # Initialize parameter vectors
    #--------------------------------------------------------------------------

    bA={}
    for k in meta['Param']['BE']['Sawtooth']['Allom'].keys():
        bA[k]=np.zeros(N_tree)

    bR={}
    for k in meta['Param']['BE']['Sawtooth']['Eq R'][ meta[pNam]['Scenario'][iScn]['Eq R CD'] ].keys():
        bR[k]=np.zeros(N_tree)

    bM={}
    for k in meta['Param']['BE']['Sawtooth']['Eq M'][ meta[pNam]['Scenario'][iScn]['Eq M CD'] ].keys():
        bM[k]=np.zeros(N_tree)

    bG={}
    for k in meta['Param']['BE']['Sawtooth']['Eq G'][ meta[pNam]['Scenario'][iScn]['Eq G CD'] ].keys():
        bG[k]=np.zeros(N_tree)

    #--------------------------------------------------------------------------
    # Populate parameter vectors based on species to each tree based on fractions
    # from inventory.
    #--------------------------------------------------------------------------

    uS=np.unique(tl['ID_SRS'])
    for iU in range(uS.size):

        ind0=np.where(tl['ID_SRS'][0,:]==uS[iU])[0]

        ind1=np.where(meta['Param']['BE']['Sawtooth']['Allom']['SRS_ID']==uS[iU])[0]
        for k in bA.keys():
            if (k=='SRS_ID') | (k=='SRS_CD'):
                continue
            bA[k][ind0]=meta['Param']['BE']['Sawtooth']['Allom'][k][ind1]

        ind1=np.where(meta['Param']['BE']['Sawtooth']['Eq R'][ meta[pNam]['Scenario'][iScn]['Eq R CD'] ]['SRS_ID']==uS[iU])[0]
        for k in bR.keys():
            if (k=='SRS_ID') | (k=='SRS_CD'):
                continue
            bR[k][ind0]=meta['Param']['BE']['Sawtooth']['Eq R'][ meta[pNam]['Scenario'][iScn]['Eq R CD'] ][k][ind1]

        ind1=np.where(meta['Param']['BE']['Sawtooth']['Eq M'][ meta[pNam]['Scenario'][iScn]['Eq M CD'] ]['SRS_ID']==uS[iU])[0]
        for k in bM.keys():
            if (k=='SRS_ID') | (k=='SRS_CD'):
                continue
            bM[k][ind0]=meta['Param']['BE']['Sawtooth']['Eq M'][ meta[pNam]['Scenario'][iScn]['Eq M CD'] ][k][ind1]

        ind1=np.where(meta['Param']['BE']['Sawtooth']['Eq G'][ meta[pNam]['Scenario'][iScn]['Eq G CD'] ]['SRS_ID']==uS[iU])[0]
        for k in bG.keys():
            if (k=='SRS_ID') | (k=='SRS_CD'):
                continue
            bG[k][ind0]=meta['Param']['BE']['Sawtooth']['Eq G'][ meta[pNam]['Scenario'][iScn]['Eq G CD'] ][k][ind1]

    #--------------------------------------------------------------------------
    # Germination
    #--------------------------------------------------------------------------

    # In the first time step, assume 1...seed_n seedlings germinate with initial
    # biomass drawn from a normal distribution, with a very small mean
    # biomass and a coefficient of variation of 10 percent. Force values to be positive.
    # Once a few live trees have been Invscribed, recruitment rate should
    # rise quickly because stand biomass (the primary driver of recruitment)
    # is low. This should lead to repopulation of the stand assuming mortality
    # does not exceed recruitment, and growth is viable.

    # Planted:
    # seed_n=1100;
    # seed_mu=0.06; % (kg C tree-1)

    # Natural:
    seed_n=meta['Param']['BE']['Sawtooth']['Core']['Natural Germinant SPH'].astype(int)
    seed_mass_mu=meta['Param']['BE']['Sawtooth']['Core']['Natural Germinant Mass Mean'] # (kg C tree-1)
    seed_mass_sd=meta['Param']['BE']['Sawtooth']['Core']['Natural Germinant Mass SD'] # (kg C tree-1)
    seed_min=meta['Param']['BE']['Sawtooth']['Core']['Natural Germinant Mass Min'] # (kg C tree-1)
    Csw_seed=np.random.normal(seed_mass_mu,seed_mass_sd,(1,seed_n))
    Csw_seed=np.maximum(seed_min,Csw_seed)

    # Populate
    tl['Csw'][0,0:seed_n]=Csw_seed

    # *** This needs updating ***
    tl['H'][0,0:seed_n]=bA['Cag2H1'][0:seed_n]*((1-np.exp(-bA['Cag2H2'][0:seed_n]*tl['Csw'][0,0:seed_n])))**(1/(1-bA['Cag2H3'][0:seed_n]))

    #--------------------------------------------------------------------------
    # Reshape allometric parameters
    #--------------------------------------------------------------------------

    # Foliage biomass
    Csw2Cf1=np.tile( np.reshape(bA['Csw2Cf1'],(1,-1)),(N_time,1) )
    Csw2Cf2=np.tile( np.reshape(bA['Csw2Cf2'],(1,-1)),(N_time,1) )

    # Branch biomass
    Csw2Cbk1=np.tile( np.reshape(bA['Csw2Cbk1'],(1,-1)),(N_time,1) )
    Csw2Cbk2=np.tile( np.reshape(bA['Csw2Cbk2'],(1,-1)),(N_time,1) )

    # Branch biomass
    Csw2Cbr1=np.tile( np.reshape(bA['Csw2Cbr1'],(1,-1)),(N_time,1) )
    Csw2Cbr2=np.tile( np.reshape(bA['Csw2Cbr2'],(1,-1)),(N_time,1) )

    #--------------------------------------------------------------------------
    # Function for calculating total biomass from stemwood biomass
    #--------------------------------------------------------------------------

    def ExpandBiomass(Csw_tot,Csw_merch):

        # Merchantable stemwood
        SWm=Csw_merch*meta['Param']['BE']['Sawtooth']['Core']['Merch Stemwood Fraction']

        # Non-merchantable stemwood
        SWnm=Csw_tot-SWm

        # Foliage biomass
        F=Csw_tot*Csw2Cf1*Csw_tot**Csw2Cf2

        # Bark biomass
        Bk=Csw_tot*Csw2Cbk1*Csw_tot**Csw2Cbk2

        # Branch biomass
        Br=Csw_tot*Csw2Cbr1*Csw_tot**Csw2Cbr2

        # Total aboveground biomass
        AG=Csw_tot+F+Bk+Br

        # Total root biomass (Li et al. 2003)

        # Conifer
        R=0.222*AG

        # Deciduous (equation from Li et al. XXXX)
        ind=np.where(tl['ID_Decid']==1)[0]
        if ind.size!=0:
            R[ind]=1.576*AG[ind]**0.615

        # Fine root biomass (equation from Li et al. XXXX)
        Rf=R*(0.072+0.354*np.exp(-0.06*(2*R)))

        # Coarse root biomass
        Rc=R-Rf

        # Total biomass
        Tot=AG+R

        return SWm,SWnm,F,Bk,Br,Rc,Rf,Tot

    #--------------------------------------------------------------------------
    # Loop through time intervals (start in second time step)
    #--------------------------------------------------------------------------

    for iT in range(1,N_time):

        #----------------------------------------------------------------------
        # Apply first disturbance cycle to subsequent disturbance cycles during
        # the spinup period to save time.
        # *** Not working properly - just needs maintenance ***
        #----------------------------------------------------------------------

#        if (vi['tv'][iT]>meta[pNam]['Project']['SpinupSpanFastTrack'][iScn]['Start']+1) & (vi['tv'][iT]<=meta[pNam]['Project']['SpinupSpanFastTrack'][iScn]['End']+1):
#
#            iT0=iT-meta[pNam]['Project']['Spinup Disturbance Return Inverval']
#
#            tl['A'][iT,:]=tl['A'][iT0,:]
#            tl['H'][iT,:]=tl['H'][iT0,:]
#            tl['D'][iT,:]=tl['D'][iT0,:]
#            tl['N_R'][iT,:]=tl['N_R'][iT0,:]
#            tl['N_M_Reg'][iT,:]=tl['N_M_Reg'][iT0,:]
#            tl['Csw'][iT,:]=tl['Csw'][iT0,:]
#            tl['Csw_Larger'][iT,:]=tl['Csw_Larger'][iT0,:]
#            tl['Csw_G'][iT,:]=tl['Csw_G'][iT0,:]
#            tl['Csw_M_Reg'][iT,:]=tl['Csw_M_Reg'][iT0,:]
#            vo['A'][iT,iS]=vo['A'][iT0,iS]
#
#            continue

        #----------------------------------------------------------------------

        # Index to live trees
        iLive=np.array(np.where(~np.isnan(tl['Csw'][iT-1,:]))).flatten()
        nLive=iLive.size

        # Index to dead trees
        #iDead=np.array(np.where(np.isnan(tl['Csw'][iT-1,:])))
        iDead=np.where(np.isnan(tl['Csw'][iT-1,:])==True)[0]
        nDead=iDead.size

        # Create random number vectors that will be used for calculation of
        # annual probabilities of recruitment and mortality
        rLive=np.random.uniform(0,1,nLive)
        rDead=np.random.uniform(0,1,nDead)

        # Update stand age (i.e., time since stand-replacing disturbance)
        vo['A'][iT,iS]=vo['A'][iT-1,iS]+1

        # Update tree age
        tl['A'][iT,iLive]=tl['A'][iT-1,iLive]+1

        tmp=tl['D'][iT-1,iLive]
        tmp=tmp[tmp>=0]
        try:
            kde=stats.gaussian_kde(tmp)
        except:
            print(tmp.shape)
        try:
            vo['DBH_Class'][iT,iS,:]=kde(meta['Core']['Sawtooth']['DBH Classes'])
        except:
            print(vo['DBH_Class'][iT,iS,:].shape)
            print(meta['Core']['Sawtooth']['DBH Classes'].shape)

        #data=np.random.random(100)
        #kde=stats.gaussian_kde(data)
        #a=kde(np.linspace(-2,2,100))

        #----------------------------------------------------------------------
        # Calculate predictor variables used in the equations of recruitment,
        # mortality and growth
        #----------------------------------------------------------------------

        # Tree age
        A=tl['A'][iT-1,:]

        # Aboveground biomass of individual trees from t-1 (kg C tree-1)
        Csw=tl['Csw'][iT-1,:]

        # Tree height (m)
        H=tl['H'][iT-1,:]

        # Stand age
        SA=np.nanmean(tl['A'][iT-1,:])

        # Stand-level biomass from t-1 (Mg C ha-1)
        SCsw=np.nansum(tl['Csw'][iT-1,:])/cm

        # Stand density from t-1 (stems ha-1)
        SN=nLive

        # Biomass of larger trees (Mg C ha-1)
        ListToSort=np.zeros((N_tree,2))
        ListToSort[:,0]=np.arange(0,N_tree,1)
        ListToSort[:,1]=Csw/cm
        SortedList=ListToSort[ListToSort[:,1].argsort()]
        SortedList=np.flip(SortedList,0)
        tmp=np.reshape(np.nancumsum(SortedList[:,1])-SortedList[:,1]/1000,(N_tree,1))
        SortedList2=np.append(SortedList,tmp,axis=1)

        Csw_Larger=np.zeros(N_tree)
        Csw_Larger[SortedList2[:,0].astype(int)]=SortedList2[:,2]

        tl['Csw_Larger'][iT,:]=Csw_Larger

        #----------------------------------------------------------------------
        # Probability of recruitment
        #----------------------------------------------------------------------

        # Only do this if there are living trees
        if iDead.size!=0:

            if meta[pNam]['Scenario'][iScn]['Eq R CD']=='Def1':
                SCsw_z=(SCsw-bR['SB_mu'])/bR['SB_sig']
                lgit=bR['Int']+bR['SB']*SCsw_z
                Pr=(np.exp(lgit)/(1+np.exp(lgit)))

            # Establish trees based on annual probability of recruitment.
            # Initial values of biomass and height are set low, arbitrary until
            # they can be set according to actual observations
            # at end of first year.

            iRec=np.where(Pr[iDead]>=rDead)[0]

            tl['A'][iT:,iDead[iRec]]=1
            tl['H'][iT:,iDead[iRec]]=0.1
            tl['Csw'][iT:,iDead[iRec]]=0.05
            tl['N_R'][iT,iDead[iRec]]=1

        #----------------------------------------------------------------------
        # Growth of stemwood biomass (kg C tree-1 yr-1)
        #----------------------------------------------------------------------

        if meta[pNam]['Scenario'][iScn]['Eq G CD']=='Def1':

            # Standardization
            LnCsw_z=(np.log(Csw)-bG['LnB_mu'])/bG['LnB_sig']
            Csw_z=(Csw-bG['B_mu'])/bG['B_sig']
            SA_z=(A-bG['SA_mu'])/bG['SA_sig']
            SCswLT_z=(Csw_Larger-bG['SBLT_mu'])/bG['SBLT_sig']
            SCsw_z=(SCsw-bG['SB_mu'])/bG['SB_sig']

            # Add all effects to intercept
            yhat=bG['Int'] + bG['LnB']*LnCsw_z + bG['B']*Csw_z + bG['SA']*SA_z + bG['SBLT']*SCswLT_z + bG['SB']*SCsw_z

            # Back-transform and apply log correction
            yhat=np.exp(yhat)

            # Cap unrealistic growth
            #yhat[yhat>G_max]=G_max

            # Populate tree level structure with growth predictions
            tl['Csw_G'][iT,iLive]=yhat[iLive]

        #----------------------------------------------------------------------
        # Update state variables
        #----------------------------------------------------------------------

        # Update stemwood biomass (kg C tree-1)
        tl['Csw'][iT,iLive]=tl['Csw'][iT-1,iLive]+tl['Csw_G'][iT,iLive]

        # Update tree height (m)
        # *** This needs updating ***
        tl['H'][iT,:]=bA['Cag2H1']*((1-np.exp(-bA['Cag2H2']*tl['Csw'][iT,:]))**(1/(1-bA['Cag2H3'])))

        # Update diameter (cm)
        tl['D'][iT,:]=(tl['Csw'][iT,:]/0.5/meta['Param']['BE']['Sawtooth']['Core']['D_to_Bsw_b0'])**(1/meta['Param']['BE']['Sawtooth']['Core']['D_to_Bsw_b1'])

        #----------------------------------------------------------------------
        # Probability of tree mortality (regular)
        #----------------------------------------------------------------------

        if meta[pNam]['Scenario'][iScn]['Eq M CD']=='Def1':

            Csw_z=(Csw-bM['B_mu'])/bM['B_sig']
            Csw2_z=(Csw**2-bM['B2_mu'])/bM['B2_sig']
            SA_z=(A-bM['SA_mu'])/bM['SA_sig']
            SCswLT_z=(Csw_Larger-bM['SBLT_mu'])/bM['SBLT_sig']
            SCsw_z=(SCsw-bM['SB_mu'])/bM['SB_sig']

            lgit=bM['Int'] + bM['B']*Csw_z + bM['B2']*Csw2_z + bM['SA']*SA_z + bM['SBLT']*SCswLT_z + bM['SB']*SCsw_z

            Pm_Sim_Reg=1*(np.exp(lgit)/(1+np.exp(lgit)))
            Pm_Sim_Ins=np.zeros(N_tree)
            Pm_Sim_Pat=np.zeros(N_tree)

        #----------------------------------------------------------------------
        # Remove biomass of trees that died
        #----------------------------------------------------------------------

        # Remove biomass of trees that died directly from competition and
        # environmental conditions
        iKill=np.where(Pm_Sim_Reg[iLive]>=rLive)[0]
        Csw_M=tl['Csw'][iT,iLive[iKill]].copy()
        tl['Csw_M_Reg'][iT,iLive[iKill]]=Csw_M
        tl['N_M_Reg'][iT,iLive[iKill]]=1
        tl['A'][iT:,iLive[iKill]]=np.nan
        tl['H'][iT:,iLive[iKill]]=np.nan
        tl['D'][iT+1:,iLive[iKill]]=np.nan
        tl['Csw'][iT:,iLive[iKill]]=np.nan
        rLive[iKill]=10; # Update rLive to avoid double-counting mortality

        # Remove biomass of trees that die due to insect attack
        iKill=np.where(Pm_Sim_Ins[iLive]>=rLive)[0]
        Csw_M=tl['Csw'][iT,iLive[iKill]].copy()
        tl['Csw_M_Ins'][iT,iLive[iKill]]=Csw_M
        tl['N_M_Ins'][iT,iLive[iKill]]=1
        tl['A'][iT:,iLive[iKill]]=np.nan
        tl['H'][iT:,iLive[iKill]]=np.nan
        tl['D'][iT+1:,iLive[iKill]]=np.nan
        tl['Csw'][iT:,iLive[iKill]]=np.nan
        rLive[iKill]=10

        # Remove biomass of trees that die due to pathogen infection
        iKill=np.where(Pm_Sim_Pat[iLive]>=rLive)[0]
        Csw_M=tl['Csw'][iT,iLive[iKill]].copy()
        tl['Csw_M_Pat'][iT,iLive[iKill]]=Csw_M
        tl['N_M_Pat'][iT,iLive[iKill]]=1
        tl['A'][iT:,iLive[iKill]]=np.nan
        tl['H'][iT:,iLive[iKill]]=np.nan
        tl['D'][iT+1:,iLive[iKill]]=np.nan
        tl['Csw'][iT:,iLive[iKill]]=np.nan
        rLive[iKill]=10

        #----------------------------------------------------------------------
        # Mortality from inventory data sources
        #----------------------------------------------------------------------

        # Loop through events in year
        for iE in range(meta['Core']['Max Events Per Year']):

            # Event type IDs for the iE'th event of the year
            ID_Type=vi['EC']['ID Event Type'][iT,iS,iE].copy()

            # Total affected biomass carbon
            MF=vi['EC']['Mortality Factor'][iT,iS,iE].copy()

            if MF>0:

                iLive=np.where(np.isnan(tl['Csw'][iT,:])==False)[0]

                nKill=int(MF*iLive.size)

                if ID_Type==meta['LUT']['Event']['Wildfire']:
                    # Wildfire-------------------------------------------------

                    iKill=np.arange(0,nKill)
                    tl['Csw_M_Fir'][iT,iLive[iKill]]=tl['Csw'][iT,iLive[iKill]].copy()
                    tl['N_M_Fir'][iT,iLive[iKill]]=1
                    vo['A'][iT,iS]=0

                elif (ID_Type==meta['LUT']['Event']['Harvest']) | (ID_Type==meta['LUT']['Event']['Harvest Salvage']):
                    # Harvest--------------------------------------------------

                    iKill=np.arange(0,nKill)
                    tl['Csw_M_Har'][iT,iLive[iKill]]=tl['Csw'][iT,iLive[iKill]].copy()
                    tl['N_M_Har'][iT,iLive[iKill]]=1
                    vo['A'][iT,iS]=0

                elif (ID_Type==meta['LUT']['Event']['Sawtooth Commercial Thinning']):
                    # Commercial thinning--------------------------------------

                    idx_sorted=np.flip(np.argsort(tl['Csw'][iT,iLive]))
                    iKill=idx_sorted[0:nKill]
                    tl['Csw_M_Har'][iT,iLive[iKill]]=tl['Csw'][iT,iLive[iKill]].copy()
                    tl['N_M_Har'][iT,iLive[iKill]]=1

                elif (ID_Type==meta['LUT']['Event']['Beetles']) | (ID_Type==meta['LUT']['Event']['IBM']) | (ID_Type==meta['LUT']['Event']['Defoliators']):
                    # Insects--------------------------------------------------

                    iKill=np.arange(0,nKill)
                    tl['Csw_M_Ins'][iT,iLive[iKill]]=tl['Csw'][iT,iLive[iKill]].copy()
                    tl['N_M_Ins'][iT,iLive[iKill]]=1

                elif (ID_Type==meta['LUT']['Event']['Sawtooth IDW']):
                    # Sawtooth IDW---------------------------------------------

                    Csw_interm=np.percentile(tl['Csw'][iT,iLive],33)
                    idx_sorted=np.argsort( np.abs(tl['Csw'][iT,iLive]-Csw_interm) )
                    iKill=idx_sorted[0:nKill]
                    tl['Csw_M_Ins'][iT,iLive[iKill]]=tl['Csw'][iT,iLive[iKill]].copy()
                    tl['N_M_Ins'][iT,iLive[iKill]]=1

                else:
                    pass

                # Update tree variables
                tl['A'][iT:,iLive[iKill]]=0
                tl['H'][iT:,iLive[iKill]]=np.nan
                tl['D'][iT+1:,iLive[iKill]]=np.nan
                tl['Csw'][iT:,iLive[iKill]]=np.nan

    #--------------------------------------------------------------------------
    # Calculate biomass of other tissues from stemwood biomass
    #--------------------------------------------------------------------------

    # Calculate merch stemwood biomass
    iMerch=np.where(tl['D']>12.5)

    Csw_merch=np.zeros(tl['Csw'].shape); Csw_merch[iMerch]=tl['Csw'][iMerch]
    Csw_merch_M_Reg=np.zeros(tl['Csw'].shape); Csw_merch_M_Reg[iMerch]=tl['Csw_M_Reg'][iMerch]
    Csw_merch_M_Fir=np.zeros(tl['Csw'].shape); Csw_merch_M_Fir[iMerch]=tl['Csw_M_Fir'][iMerch]
    Csw_merch_M_Ins=np.zeros(tl['Csw'].shape); Csw_merch_M_Ins[iMerch]=tl['Csw_M_Ins'][iMerch]
    Csw_merch_M_Pat=np.zeros(tl['Csw'].shape); Csw_merch_M_Pat[iMerch]=tl['Csw_M_Pat'][iMerch]
    Csw_merch_M_Har=np.zeros(tl['Csw'].shape); Csw_merch_M_Har[iMerch]=tl['Csw_M_Har'][iMerch]

    # Standing biomass
    tl['Cswm'],tl['Cswnm'],tl['Cf'],tl['Cbk'],tl['Cbr'],tl['Crc'],tl['Crf'],tl['Ctot']=ExpandBiomass(tl['Csw'],Csw_merch)

    # Biomass loss due to regular mortality
    tl['Cswm_M_Reg'],tl['Cswnm_M_Reg'],tl['Cf_M_Reg'],tl['Cbk_M_Reg'],tl['Cbr_M_Reg'],tl['Crc_M_Reg'],tl['Crf_M_Reg'],tl['Ctot_M_Reg']=ExpandBiomass(tl['Csw_M_Reg'],Csw_merch_M_Reg)

    # Biomass loss due to wildfire
    tl['Cswm_M_Fir'],tl['Cswnm_M_Fir'],tl['Cf_M_Fir'],tl['Cbk_M_Fir'],tl['Cbr_M_Fir'],tl['Crc_M_Fir'],tl['Crf_M_Fir'],tl['Ctot_M_Fir']=ExpandBiomass(tl['Csw_M_Fir'],Csw_merch_M_Fir)

    # Biomass loss due to insects
    tl['Cswm_M_Ins'],tl['Cswnm_M_Ins'],tl['Cf_M_Ins'],tl['Cbk_M_Ins'],tl['Cbr_M_Ins'],tl['Crc_M_Ins'],tl['Crf_M_Ins'],tl['Ctot_M_Ins']=ExpandBiomass(tl['Csw_M_Ins'],Csw_merch_M_Ins)

    # Biomass loss due to pathogens
    tl['Cswm_M_Pat'],tl['Cswnm_M_Pat'],tl['Cf_M_Pat'],tl['Cbk_M_Pat'],tl['Cbr_M_Pat'],tl['Crc_M_Pat'],tl['Crf_M_Pat'],tl['Ctot_M_Pat']=ExpandBiomass(tl['Csw_M_Pat'],Csw_merch_M_Pat)

    # Biomass loss due to harvest
    tl['Cswm_M_Har'],tl['Cswnm_M_Har'],tl['Cf_M_Har'],tl['Cbk_M_Har'],tl['Cbr_M_Har'],tl['Crc_M_Har'],tl['Crf_M_Har'],tl['Ctot_M_Har']=ExpandBiomass(tl['Csw_M_Har'],Csw_merch_M_Har)

    # Biomass loss due to wind
    #tl['Cswm_M_Win'],tl['Cswnm_M_Win'],tl['Cf_M_Win'],tl['Cbk_M_Win'],tl['Cbr_M_Win'],tl['Crc_M_Win'],tl['Crf_M_Win'],tl['Ctot_M_Win']=ExpandBiomass(tl['Csw_M_Win'])

    #--------------------------------------------------------------------------
    # Biomass turnover (litterfall)
    #--------------------------------------------------------------------------

    # Foliage turnover
    tl['Cf_LF']=0.100*tl['Cf']

    # Branch turnover
    tl['Cbr_LF']=0.035*tl['Cbr']

    # Bark turnover
    tl['Cbk_LF']=0.035*tl['Cbk']

    # Fine root turnover
    tl['Crf_LF']=0.641*tl['Crf']

    # Coarse root turnover
    tl['Crc_LF']=0.02*tl['Crc']

    #--------------------------------------------------------------------------
    # Populate stand variables (individual-tree statistics)
    #--------------------------------------------------------------------------

    # Mean tree age
    vo['TreeMean_A'][:,iS]=np.nanmean(tl['A'],axis=1)

    # Mean tree diamter at breast height
    vo['TreeMean_D'][:,iS]=np.nanmean(tl['D'],axis=1)

    # Mean tree height (m)
    vo['TreeMean_H'][:,iS]=np.nanmean(tl['H'],axis=1)

    # Mean tree stemwood biomass (kg C)
    vo['TreeMean_Csw'][:,iS]=np.maximum(0,np.nanmean(tl['Csw'],axis=1))

    # Mean stemwood biomass growth (kg C yr-1)
    vo['TreeMean_Csw_G'][:,iS]=np.maximum(0,np.nanmean(tl['Csw_G'],axis=1))

    #--------------------------------------------------------------------------
    # Populate stand-level demographics
    #--------------------------------------------------------------------------

    # Stand density (stems ha-1)
    vo['N'][:,iS]=np.sum(~np.isnan(tl['Csw']),axis=1)

    # Stand density at beginning of time step
    N_time0=np.append(vo['N'][0,iS],vo['N'][0:-1,iS])

    # Demographic recruitment rate (% yr-1)
    vo['N_R'][:,iS]=np.minimum(100,np.sum(tl['N_R'],axis=1)/N_time0*100)

    # Demographic mortality rate (% yr-1)

    # Need to update N0 or else stand-replacing disturbances will have lower
    # relative mortality rate
    # *** The other solution is to move the prescribed disturbances above the
    # background mortality ***
    #N0=N0-Sum_N_M_Sim_Reg-Sum_N_M_Sim_Ins-Sum_N_M_Sim_Pat;

    N_M_Tot=tl['N_M_Reg']+tl['N_M_Fir']+tl['N_M_Ins']+tl['N_M_Pat']+tl['N_M_Har']+tl['N_M_Win']

    vo['N_M_Tot'][:,iS]=np.maximum(0,np.minimum(100,np.sum(N_M_Tot,axis=1)/N_time0*100))

    #    vo['N_M_Reg'][:,iS]=np.maximum(0,np.sum(tl['N_M_Reg'],axis=1)/N_time0*100)
    #    vo['N_M_Fir'][:,iS]=np.maximum(0,np.sum(tl['N_M_Fir'],axis=1)/N_time0*100)
    #    vo['N_M_Ins'][:,iS]=np.maximum(0,np.sum(tl['N_M_Ins'],axis=1)/N_time0*100)
    #    vo['N_M_Pat'][:,iS]=np.maximum(0,np.sum(tl['N_M_Pat'],axis=1)/N_time0*100)
    #    vo['N_M_Win'][:,iS]=np.maximum(0,np.sum(tl['N_M_Win'],axis=1)/N_time0*100)
    #    vo['N_M_Har'][:,iS]=np.maximum(0,np.sum(tl['N_M_Har'],axis=1)/N_time0*100)

    #--------------------------------------------------------------------------
    # Populate stand-level biomass
    #--------------------------------------------------------------------------

    # Gross growth of biomass (Mg C ha-1 yr-1)
    vo['C_G_Gross'][1:,iS,0]=np.nansum(np.maximum(0,np.diff(tl['Cswm'],axis=0)),axis=1)/cm
    vo['C_G_Gross'][1:,iS,1]=np.nansum(np.maximum(0,np.diff(tl['Cswnm'],axis=0)),axis=1)/cm
    vo['C_G_Gross'][1:,iS,2]=np.nansum(np.maximum(0,np.diff(tl['Cf'],axis=0)),axis=1)/cm
    vo['C_G_Gross'][1:,iS,3]=np.nansum(np.maximum(0,np.diff(tl['Cbr'],axis=0)),axis=1)/cm
    vo['C_G_Gross'][1:,iS,4]=np.nansum(np.maximum(0,np.diff(tl['Cbk'],axis=0)),axis=1)/cm
    vo['C_G_Gross'][1:,iS,5]=np.nansum(np.maximum(0,np.diff(tl['Crc'],axis=0)),axis=1)/cm
    vo['C_G_Gross'][1:,iS,6]=np.nansum(np.maximum(0,np.diff(tl['Crf'],axis=0)),axis=1)/cm

    # Biomass turnover (Mg C ha-1 yr-1)
    vo['C_LF'][:,iS,2]=np.nansum(tl['Cf_LF'],axis=1)/cm
    vo['C_LF'][:,iS,3]=np.nansum(tl['Cbr_LF'],axis=1)/cm
    vo['C_LF'][:,iS,4]=np.nansum(tl['Cbk_LF'],axis=1)/cm
    vo['C_LF'][:,iS,5]=np.nansum(tl['Crc_LF'],axis=1)/cm
    vo['C_LF'][:,iS,6]=np.nansum(tl['Crf_LF'],axis=1)/cm

    # Biomass loss due to regular mortality (Mg C ha-1 yr-1)
    vo['C_M_Reg'][:,iS,0]=np.nansum(tl['Cswm_M_Reg'],axis=1)/cm
    vo['C_M_Reg'][:,iS,1]=np.nansum(tl['Cswnm_M_Reg'],axis=1)/cm
    vo['C_M_Reg'][:,iS,2]=np.nansum(tl['Cf_M_Reg'],axis=1)/cm
    vo['C_M_Reg'][:,iS,3]=np.nansum(tl['Cbr_M_Reg'],axis=1)/cm
    vo['C_M_Reg'][:,iS,4]=np.nansum(tl['Cbk_M_Reg'],axis=1)/cm
    vo['C_M_Reg'][:,iS,5]=np.nansum(tl['Crc_M_Reg'],axis=1)/cm
    vo['C_M_Reg'][:,iS,6]=np.nansum(tl['Crf_M_Reg'],axis=1)/cm

    # Biomass loss due to all mortality (Mg C ha-1 yr-1)
    vo['C_M_Tot'][:,iS,0]=np.nansum(tl['Cswm_M_Reg']+tl['Cswm_M_Fir']+tl['Cswm_M_Ins']+tl['Cswm_M_Pat']+tl['Cswm_M_Har'],axis=1)/cm
    vo['C_M_Tot'][:,iS,1]=np.nansum(tl['Cswnm_M_Reg']+tl['Cswnm_M_Fir']+tl['Cswnm_M_Ins']+tl['Cswnm_M_Pat']+tl['Cswnm_M_Har'],axis=1)/cm
    vo['C_M_Tot'][:,iS,2]=np.nansum(tl['Cf_M_Reg']+tl['Cf_M_Fir']+tl['Cswnm_M_Ins']+tl['Cswnm_M_Pat']+tl['Cswnm_M_Har'],axis=1)/cm
    vo['C_M_Tot'][:,iS,3]=np.nansum(tl['Cbr_M_Reg']+tl['Cbr_M_Fir']+tl['Cbr_M_Ins']+tl['Cbr_M_Pat']+tl['Cbr_M_Har'],axis=1)/cm
    vo['C_M_Tot'][:,iS,4]=np.nansum(tl['Cbk_M_Reg']+tl['Cbk_M_Fir']+tl['Cbk_M_Ins']+tl['Cbk_M_Pat']+tl['Cbk_M_Har'],axis=1)/cm
    vo['C_M_Tot'][:,iS,5]=np.nansum(tl['Crc_M_Reg']+tl['Crc_M_Fir']+tl['Crc_M_Ins']+tl['Crc_M_Pat']+tl['Crc_M_Har'],axis=1)/cm
    vo['C_M_Tot'][:,iS,6]=np.nansum(tl['Crf_M_Reg']+tl['Crf_M_Fir']+tl['Crf_M_Ins']+tl['Crf_M_Pat']+tl['Crf_M_Har'],axis=1)/cm

    #    vo['C_M_Fir'][:,iS]=np.nansum(tl['Ctot_M_Fir'],axis=1)/cm
    #    vo['C_M_Ins'][:,iS]=np.nansum(tl['Ctot_M_Ins'],axis=1)/cm
    #    vo['C_M_Pat'][:,iS]=np.nansum(tl['Ctot_M_Pat'],axis=1)/cm
    #    vo['C_M_Win'][:,iS]=np.nansum(tl['Ctot_M_Win'],axis=1)/cm
    #    vo['C_M_Har'][:,iS]=np.nansum(tl['Ctot_M_Har'],axis=1)/cm

    # Net growth by tissue (Mg C ha-1 yr-1)
    vo['C_G_Net'][:,iS,0:7]=vo['C_G_Gross'][:,iS,0:7]-vo['C_M_Reg'][:,iS,0:7]

    # Ecosystem carbon pools by tissue (Mg C ha-1 yr-1)
    vo['C_Eco_Pools'][:,iS,0]=np.maximum(0,np.nansum(tl['Cswm'],axis=1))/cm
    vo['C_Eco_Pools'][:,iS,1]=np.maximum(0,np.nansum(tl['Cswnm'],axis=1))/cm
    vo['C_Eco_Pools'][:,iS,2]=np.maximum(0,np.nansum(tl['Cf'],axis=1))/cm
    vo['C_Eco_Pools'][:,iS,3]=np.maximum(0,np.nansum(tl['Cbr'],axis=1))/cm
    vo['C_Eco_Pools'][:,iS,4]=np.maximum(0,np.nansum(tl['Cbk'],axis=1))/cm
    vo['C_Eco_Pools'][:,iS,5]=np.maximum(0,np.nansum(tl['Crc'],axis=1))/cm
    vo['C_Eco_Pools'][:,iS,6]=np.maximum(0,np.nansum(tl['Crf'],axis=1))/cm

    #--------------------------------------------------------------------------
    # Populate stand-level merch stemwood volume
    #--------------------------------------------------------------------------

    # Merchantable stemwood volume of live trees
    vo['V_MerchLive'][:,iS]=np.maximum(0,np.nansum(tl['Cswm'],axis=1))/cm/meta['Param']['BEV']['Biophysical']['Ratio_Wood_C_to_DM']/0.45

    #t1=time.time()
    #print((t1-t0))

    return vo

#%% STAND DEAD ORGANIC MATTER DYNAMICS

def DOM_LikeKurzetal2009(meta,pNam,iT,iBat,vi,vo,iEP):

    # Extract parameters
    bIPF=meta['Param']['BEV']['Inter Pool Fluxes']
    bDec=meta['Param']['BEV']['Decomp']

    #--------------------------------------------------------------------------
    # Flux of carbon between biomass pools and dead organic matter pools
    # Current DOM pools = previous DOM pools + biomass turnover
    #--------------------------------------------------------------------------

    # Transfer biomass turnover to very fast litter pool
    vo['C_Eco_Pools'][iT,:,iEP['LitterVF']]=vo['C_Eco_Pools'][iT-1,:,iEP['LitterVF']]+ \
        bIPF['FoliageLitToLitterVF']*vo['C_LF'][iT,:,iEP['Foliage']]+ \
        bIPF['RootFineLitToLitterVF']*vo['C_LF'][iT,:,iEP['RootFine']]+ \
        bIPF['FoliageMorToLitterVF']*vo['C_M_Reg'][iT,:,iEP['Foliage']]

    # Transfer biomass turnover to fast litter pool
    vo['C_Eco_Pools'][iT,:,iEP['LitterF']]=vo['C_Eco_Pools'][iT-1,:,iEP['LitterF']]+ \
        bIPF['BranchLitToLitterF']*vo['C_LF'][iT,:,iEP['Branch']]+ \
        bIPF['BarkLitToLitterF']*vo['C_LF'][iT,:,iEP['Bark']]+ \
        bIPF['RootCoarseLitToLitterF']*vo['C_LF'][iT,:,iEP['RootCoarse']]+ \
        bIPF['BranchMorToLitterF']*vo['C_M_Reg'][iT,:,iEP['Branch']]+ \
        bIPF['StemNonMerchLitToLitterF']*vo['C_LF'][iT,:,iEP['StemNonMerch']]

    # Transfer biomass turnover to medium litter pool
    vo['C_Eco_Pools'][iT,:,iEP['LitterM']]=vo['C_Eco_Pools'][iT-1,:,iEP['LitterM']]+ \
        bIPF['BarkLitToLitterM']*vo['C_LF'][iT,:,iEP['Bark']]+ \
        bIPF['BarkMorToLitterM']*vo['C_M_Reg'][iT,:,iEP['Bark']]

    # Transfer biomass turnover to slow litter pool
    vo['C_Eco_Pools'][iT,:,iEP['LitterS']]=vo['C_Eco_Pools'][iT-1,:,iEP['LitterS']]+ \
        bIPF['BarkLitToLitterS']*vo['C_LF'][iT,:,iEP['Bark']]+ \
        bIPF['BarkMorToLitterS']*vo['C_M_Reg'][iT,:,iEP['Bark']]

    # Transfer biomass turnover to snag stemwood pool
    vo['C_Eco_Pools'][iT,:,iEP['SnagStem']]=vo['C_Eco_Pools'][iT-1,:,iEP['SnagStem']]+ \
        bIPF['StemMerchMorToSnagStem']*vo['C_M_Reg'][iT,:,iEP['StemMerch']]+ \
        bIPF['StemNonMerchMorToSnagStem']*vo['C_M_Reg'][iT,:,iEP['StemNonMerch']]

    # Transfer biomass turnover to snag branch pool
    vo['C_Eco_Pools'][iT,:,iEP['SnagBranch']]=vo['C_Eco_Pools'][iT-1,:,iEP['SnagBranch']]+ \
        bIPF['BranchLitToSnagBranch']*vo['C_LF'][iT,:,iEP['Branch']]+ \
        bIPF['BranchMorToSnagBranch']*vo['C_M_Reg'][iT,:,iEP['Branch']]+ \
        bIPF['BarkMorToSnagBranch']*vo['C_M_Reg'][iT,:,iEP['Bark']]+ \
        bIPF['StemNonMerchLitToSnagBranch']*vo['C_LF'][iT,:,iEP['Bark']]

    # Transfer biomass turnover to very fast soil pool
    vo['C_Eco_Pools'][iT,:,iEP['SoilVF']]=vo['C_Eco_Pools'][iT-1,:,iEP['SoilVF']]+ \
        bIPF['RootFineLitToSoilVF']*vo['C_LF'][iT,:,iEP['RootFine']] +\
        bIPF['RootCoarseMorToSoilVF']*vo['C_M_Reg'][iT,:,iEP['RootCoarse']]+ \
        bIPF['RootFineMorToSoilVF']*vo['C_M_Reg'][iT,:,iEP['RootFine']]

    # Transfer biomass turnover to fast soil pool
    vo['C_Eco_Pools'][iT,:,iEP['SoilF']]=vo['C_Eco_Pools'][iT-1,:,iEP['SoilF']]+ \
        bIPF['RootCoarseLitToSoilF']*vo['C_LF'][iT,:,iEP['RootCoarse']]+ \
        bIPF['RootCoarseMorToSoilF']*vo['C_M_Reg'][iT,:,iEP['RootCoarse']]

    # Transfer biomass turnover to slow soil pool
    vo['C_Eco_Pools'][iT,:,iEP['SoilS']]=vo['C_Eco_Pools'][iT-1,:,iEP['SoilS']]

    #--------------------------------------------------------------------------
    # Update piles
    # The above section updates DOM pools, but not piled pools
    #--------------------------------------------------------------------------

    vo['C_Eco_Pools'][iT,:,iEP['PiledStemMerch']]=vo['C_Eco_Pools'][iT-1,:,iEP['PiledStemMerch']]
    vo['C_Eco_Pools'][iT,:,iEP['PiledStemNonMerch']]=vo['C_Eco_Pools'][iT-1,:,iEP['PiledStemNonMerch']]
    vo['C_Eco_Pools'][iT,:,iEP['PiledBranch']]=vo['C_Eco_Pools'][iT-1,:,iEP['PiledBranch']]
    vo['C_Eco_Pools'][iT,:,iEP['PiledBark']]=vo['C_Eco_Pools'][iT-1,:,iEP['PiledBark']]
    vo['C_Eco_Pools'][iT,:,iEP['PiledSnagStem']]=vo['C_Eco_Pools'][iT-1,:,iEP['PiledSnagStem']]
    vo['C_Eco_Pools'][iT,:,iEP['PiledSnagBranch']]=vo['C_Eco_Pools'][iT-1,:,iEP['PiledSnagBranch']]

    #--------------------------------------------------------------------------
    # Decomposition
    # Vectorizing the parts of this that can be vectorized was attempted, but
    # runtime actually increased a bit.
    #--------------------------------------------------------------------------

    # Prepare air temperature for respiration calculation
    Tref=10
    fT=(vi['Inv']['MAT']-Tref)/10

    # Respiration rate - Note that these terms do not equal the
    # atmosphere-bound efflux from heterotrophic respiration as a fraction is
    # emitted to the atmosphere and the remaining fraction is reorganized
    # within the ecosystem.
    meta[pNam]['Project']['R_LitterVF']=bDec['LitterVF_R10']*vo['C_Eco_Pools'][iT,:,iEP['LitterVF']].flatten()*bDec['LitterVF_Q10']**fT
    meta[pNam]['Project']['R_LitterF']=bDec['LitterF_R10']*vo['C_Eco_Pools'][iT,:,iEP['LitterF']].flatten()*bDec['LitterF_Q10']**fT
    meta[pNam]['Project']['R_LitterM']=bDec['LitterM_R10']*vo['C_Eco_Pools'][iT,:,iEP['LitterM']].flatten()*bDec['LitterM_Q10']**fT
    meta[pNam]['Project']['R_LitterS']=bDec['LitterS_R10']*vo['C_Eco_Pools'][iT,:,iEP['LitterS']].flatten()*bDec['LitterS_Q10']**fT
    meta[pNam]['Project']['R_SnagStem']=bDec['SnagStem_R10']*vo['C_Eco_Pools'][iT,:,iEP['SnagStem']].flatten()*bDec['SnagStem_Q10']**fT
    meta[pNam]['Project']['R_SnagBranch']=bDec['SnagBranch_R10']*vo['C_Eco_Pools'][iT,:,iEP['SnagBranch']].flatten()*bDec['SnagBranch_Q10']**fT
    meta[pNam]['Project']['R_SoilVF']=bDec['SoilVF_R10']*vo['C_Eco_Pools'][iT,:,iEP['SoilVF']].flatten()*bDec['SoilVF_Q10']**fT
    meta[pNam]['Project']['R_SoilF']=bDec['SoilF_R10']*vo['C_Eco_Pools'][iT,:,iEP['SoilF']].flatten()*bDec['SoilF_Q10']**fT
    meta[pNam]['Project']['R_SoilS']=bDec['SoilS_R10']*vo['C_Eco_Pools'][iT,:,iEP['SoilS']].flatten()*bDec['SoilS_Q10']**fT

    meta[pNam]['Project']['R_PiledStemMerch']=bDec['Piled_R10']*vo['C_Eco_Pools'][iT,:,iEP['PiledStemMerch']].flatten()*bDec['Piled_Q10']**fT
    meta[pNam]['Project']['R_PiledStemNonMerch']=bDec['Piled_R10']*vo['C_Eco_Pools'][iT,:,iEP['PiledStemNonMerch']].flatten()*bDec['Piled_Q10']**fT
    meta[pNam]['Project']['R_PiledBranch']=bDec['Piled_R10']*vo['C_Eco_Pools'][iT,:,iEP['PiledBranch']].flatten()*bDec['Piled_Q10']**fT
    meta[pNam]['Project']['R_PiledBark']=bDec['Piled_R10']*vo['C_Eco_Pools'][iT,:,iEP['PiledBark']].flatten()*bDec['Piled_Q10']**fT
    meta[pNam]['Project']['R_PiledSnagStem']=bDec['Piled_R10']*vo['C_Eco_Pools'][iT,:,iEP['PiledSnagStem']].flatten()*bDec['Piled_Q10']**fT
    meta[pNam]['Project']['R_PiledSnagBranch']=bDec['Piled_R10']*vo['C_Eco_Pools'][iT,:,iEP['PiledSnagBranch']].flatten()*bDec['Piled_Q10']**fT

    # Adjust decomposition to account for N application response
    if meta['Modules']['Nutrient Management']['iApplication'].size>0:
        vi,vo,meta=napp.UpdateStatus(vi,vo,iT,meta,'HeterotrophicRespiration')

    # Remove respired carbon from source DOM pools
    vo['C_Eco_Pools'][iT,:,iEP['LitterVF']]=vo['C_Eco_Pools'][iT,:,iEP['LitterVF']]-meta[pNam]['Project']['R_LitterVF']
    vo['C_Eco_Pools'][iT,:,iEP['LitterF']]=vo['C_Eco_Pools'][iT,:,iEP['LitterF']]-meta[pNam]['Project']['R_LitterF']
    vo['C_Eco_Pools'][iT,:,iEP['LitterM']]=vo['C_Eco_Pools'][iT,:,iEP['LitterM']]-meta[pNam]['Project']['R_LitterM']
    vo['C_Eco_Pools'][iT,:,iEP['LitterS']]=vo['C_Eco_Pools'][iT,:,iEP['LitterS']]-meta[pNam]['Project']['R_LitterS']
    vo['C_Eco_Pools'][iT,:,iEP['SnagStem']]=vo['C_Eco_Pools'][iT,:,iEP['SnagStem']]-meta[pNam]['Project']['R_SnagStem']
    vo['C_Eco_Pools'][iT,:,iEP['SnagBranch']]=vo['C_Eco_Pools'][iT,:,iEP['SnagBranch']]-meta[pNam]['Project']['R_SnagBranch']
    vo['C_Eco_Pools'][iT,:,iEP['SoilVF']]=vo['C_Eco_Pools'][iT,:,iEP['SoilVF']]-meta[pNam]['Project']['R_SoilVF']
    vo['C_Eco_Pools'][iT,:,iEP['SoilF']]=vo['C_Eco_Pools'][iT,:,iEP['SoilF']]-meta[pNam]['Project']['R_SoilF']
    vo['C_Eco_Pools'][iT,:,iEP['SoilS']]=vo['C_Eco_Pools'][iT,:,iEP['SoilS']]-meta[pNam]['Project']['R_SoilS']

    vo['C_Eco_Pools'][iT,:,iEP['PiledStemMerch']]=vo['C_Eco_Pools'][iT,:,iEP['PiledStemMerch']]-meta[pNam]['Project']['R_PiledStemMerch']
    vo['C_Eco_Pools'][iT,:,iEP['PiledStemNonMerch']]=vo['C_Eco_Pools'][iT,:,iEP['PiledStemNonMerch']]-meta[pNam]['Project']['R_PiledStemNonMerch']
    vo['C_Eco_Pools'][iT,:,iEP['PiledBranch']]=vo['C_Eco_Pools'][iT,:,iEP['PiledBranch']]-meta[pNam]['Project']['R_PiledBranch']
    vo['C_Eco_Pools'][iT,:,iEP['PiledBark']]=vo['C_Eco_Pools'][iT,:,iEP['PiledBark']]-meta[pNam]['Project']['R_PiledBark']
    vo['C_Eco_Pools'][iT,:,iEP['PiledSnagStem']]=vo['C_Eco_Pools'][iT,:,iEP['PiledSnagStem']]-meta[pNam]['Project']['R_PiledSnagStem']
    vo['C_Eco_Pools'][iT,:,iEP['PiledSnagBranch']]=vo['C_Eco_Pools'][iT,:,iEP['PiledSnagBranch']]-meta[pNam]['Project']['R_PiledSnagBranch']

    # Re-define decayed fast litter
    vo['C_Eco_Pools'][iT,:,iEP['LitterF']]=vo['C_Eco_Pools'][iT,:,iEP['LitterF']]+bIPF['SnagBranchToLitterF']*meta[pNam]['Project']['R_SnagBranch']

    # Re-define decayed medium litter
    vo['C_Eco_Pools'][iT,:,iEP['LitterM']]=vo['C_Eco_Pools'][iT,:,iEP['LitterM']]+bIPF['SnagStemToLitterM']*meta[pNam]['Project']['R_SnagStem']

    # Re-define decayed slow litter
    vo['C_Eco_Pools'][iT,:,iEP['LitterS']]=vo['C_Eco_Pools'][iT,:,iEP['LitterS']]+bIPF['LitterVFToLitterS']*meta[pNam]['Project']['R_LitterVF']
    vo['C_Eco_Pools'][iT,:,iEP['LitterS']]=vo['C_Eco_Pools'][iT,:,iEP['LitterS']]+bIPF['LitterFToLitterS']*meta[pNam]['Project']['R_LitterF']
    vo['C_Eco_Pools'][iT,:,iEP['LitterS']]=vo['C_Eco_Pools'][iT,:,iEP['LitterS']]+bIPF['LitterMToLitterS']*meta[pNam]['Project']['R_LitterM']

    # Re-define decayed slow soil
    vo['C_Eco_Pools'][iT,:,iEP['SoilS']]=vo['C_Eco_Pools'][iT,:,iEP['SoilS']]+bIPF['SoilVFToSoilS']*meta[pNam]['Project']['R_SoilVF']
    vo['C_Eco_Pools'][iT,:,iEP['SoilS']]=vo['C_Eco_Pools'][iT,:,iEP['SoilS']]+bIPF['SoilFToSoilS']*meta[pNam]['Project']['R_SoilF']

    # Heterotrophic respiration
    vo['C_RH'][iT,:,iEP['LitterVF']]=bIPF['LitterVFToCO2']*meta[pNam]['Project']['R_LitterVF']
    vo['C_RH'][iT,:,iEP['LitterF']]=bIPF['LitterFToCO2']*meta[pNam]['Project']['R_LitterF']
    vo['C_RH'][iT,:,iEP['LitterM']]=bIPF['LitterMToCO2']*meta[pNam]['Project']['R_LitterM']
    vo['C_RH'][iT,:,iEP['LitterS']]=bIPF['LitterSToCO2']*meta[pNam]['Project']['R_LitterS']
    vo['C_RH'][iT,:,iEP['SnagStem']]=bIPF['SnagStemToCO2']*meta[pNam]['Project']['R_SnagStem']
    vo['C_RH'][iT,:,iEP['SnagBranch']]=bIPF['SnagBranchToCO2']*meta[pNam]['Project']['R_SnagBranch']
    vo['C_RH'][iT,:,iEP['SoilVF']]=bIPF['SoilVFToCO2']*meta[pNam]['Project']['R_SoilVF']
    vo['C_RH'][iT,:,iEP['SoilF']]=bIPF['SoilFToCO2']*meta[pNam]['Project']['R_SoilF']
    vo['C_RH'][iT,:,iEP['SoilS']]=bIPF['SoilSToCO2']*meta[pNam]['Project']['R_SoilS']

    vo['C_RH'][iT,:,iEP['PiledStemMerch']]=bIPF['PiledToCO2']*meta[pNam]['Project']['R_PiledStemMerch']
    vo['C_RH'][iT,:,iEP['PiledStemNonMerch']]=bIPF['PiledToCO2']*meta[pNam]['Project']['R_PiledStemNonMerch']
    vo['C_RH'][iT,:,iEP['PiledBranch']]=bIPF['PiledToCO2']*meta[pNam]['Project']['R_PiledBranch']
    vo['C_RH'][iT,:,iEP['PiledBark']]=bIPF['PiledToCO2']*meta[pNam]['Project']['R_PiledBark']
    vo['C_RH'][iT,:,iEP['PiledSnagStem']]=bIPF['PiledToCO2']*meta[pNam]['Project']['R_PiledSnagStem']
    vo['C_RH'][iT,:,iEP['PiledSnagBranch']]=bIPF['PiledToCO2']*meta[pNam]['Project']['R_PiledSnagBranch']

    #--------------------------------------------------------------------------
    # Physical transfer
    #--------------------------------------------------------------------------

    PT_LitterS=bDec['LitterS_PhysTransRate']*vo['C_Eco_Pools'][iT,:,iEP['LitterS']]
    PT_SnagStem=bDec['SnagStem_PhysTransRate']*vo['C_Eco_Pools'][iT,:,iEP['SnagStem']]
    PT_SnagBranch=bDec['SnagBranch_PhysTransRate']*vo['C_Eco_Pools'][iT,:,iEP['SnagBranch']]

    # Remove carbon that is physically transferred
    vo['C_Eco_Pools'][iT,:,iEP['LitterS']]=vo['C_Eco_Pools'][iT,:,iEP['LitterS']]-PT_LitterS
    vo['C_Eco_Pools'][iT,:,iEP['SnagStem']]=vo['C_Eco_Pools'][iT,:,iEP['SnagStem']]-PT_SnagStem
    vo['C_Eco_Pools'][iT,:,iEP['SnagBranch']]=vo['C_Eco_Pools'][iT,:,iEP['SnagBranch']]-PT_SnagBranch

    # Add decomposed carbon to more decomposed DOM pools
    vo['C_Eco_Pools'][iT,:,iEP['SoilS']]=vo['C_Eco_Pools'][iT,:,iEP['SoilS']]+PT_LitterS
    vo['C_Eco_Pools'][iT,:,iEP['LitterF']]=vo['C_Eco_Pools'][iT,:,iEP['LitterF']]+PT_SnagBranch
    vo['C_Eco_Pools'][iT,:,iEP['LitterM']]=vo['C_Eco_Pools'][iT,:,iEP['LitterM']]+PT_SnagStem

    #--------------------------------------------------------------------------
    # litter decomposition
    # *** Added for fertilization study ***
    #--------------------------------------------------------------------------

    #vo['C_Eco_Pools'][iT,:,28]=meta[pNam]['Project']['R_LitterVF']+meta[pNam]['Project']['R_LitterF']+meta[pNam]['Project']['R_LitterM']+meta[pNam]['Project']['R_LitterS']

    return vo

#%% Disturbance and management events (from TAZ)

def Events_FromTaz(meta,pNam,iT,iScn,iEns,iBat,vi,vo,iEP):

    #t0=time.time()
    # Predict stand breakup (on the fly)
    if meta[pNam]['Scenario'][iScn]['Breakup Status']=='On':
        vi=asm.PredictStandBreakup_OnTheFly(meta,pNam,vi,iT,iEns,vo['A'][iT,:])

    # Predict historical harvesting (on the fly)
    if meta[pNam]['Scenario'][iScn]['Harvest Status Historical']=='On':
        if vi['tv'][iT]<meta[pNam]['Scenario'][iScn]['Harvest Year Transition']:
            Period='Historical'
            #Volume=vo['V_MerchLive'][iT,:]+2*(1/0.45)*vo['C_Eco_Pools'][iT,:,iEP['SnagStem']]
            Volume=vo['V_MerchTotal'][iT,:]
            vi=asm.PredictHarvesting_OnTheFly(meta,pNam,vi,iT,iScn,iEns,Volume,Period)

    # Predict future harvesting (on the fly)
    if meta[pNam]['Scenario'][iScn]['Harvest Status Future']=='On':
        if vi['tv'][iT]>=meta[pNam]['Scenario'][iScn]['Harvest Year Transition']:
            Period='Future'
            #Volume=vo['V_MerchLive'][iT,:]+2*(1/0.45)*vo['C_Eco_Pools'][iT,:,iEP['SnagStem']]
            Volume=vo['V_MerchTotal'][iT,:]
            vi=asm.PredictHarvesting_OnTheFly(meta,pNam,vi,iT,iScn,iEns,Volume,Period)

    # Predict future nutrient application (on the fly)
    if meta[pNam]['Scenario'][iScn]['Nutrient Application Status']=='On':
        if vi['tv'][iT]>=meta[pNam]['Project']['Year Project']:
            vi=napp.ScheduleNutrientApplication(meta,pNam,vi,vo,iT,iScn,iEns,iBat)

    # Initialize indicator of aerial nutrient application
    flag_nutrient_application=np.zeros(meta[pNam]['Project']['Batch Size'][iBat])

    #t1=time.time()
    #meta[pNam]['Project']['Run Time Summary']['Test3']=meta[pNam]['Project']['Run Time Summary']['Test3']+t1-t0

    # Check to see how many events occur in this time step (don't do more than necessary)
    NumEventsInTimeStep=np.sum(np.sum(vi['EC']['ID Event Type'][iT,:,:]>0,axis=0)>0)

    # Keep track of total tree biomass at the start of each annual time step
    C_Biomass_t0=np.sum(vo['C_Eco_Pools'][iT,:,iEP['BiomassTotal']],axis=0)

    # Loop through events in year
    for iE in range(NumEventsInTimeStep): # meta['Core']['Max Events Per Year']

        # Event type IDs for the iE'th event of the year
        ID_Type=vi['EC']['ID Event Type'][iT,:,iE].copy()

        # Total affected biomass carbon
        MortalityFactor=vi['EC']['Mortality Factor'][iT,:,iE].copy()

        # Record stands with aerial nutrient application
        #t0=time.time()
        iApp=np.where( (ID_Type==meta['LUT']['Event']['Fertilization Aerial']) | (ID_Type==meta['LUT']['Event']['Fertilization Hand']) )[0]
        flag_nutrient_application[iApp]=1
        #t1=time.time()
        #meta[pNam]['Project']['Run Time Summary']['Test9']=meta[pNam]['Project']['Run Time Summary']['Test9']+t1-t0

        #----------------------------------------------------------------------
        # Get event-specific parameters
        #----------------------------------------------------------------------
        #t0=time.time()
        u,idx,inver=np.unique(ID_Type,return_index=True,return_inverse=True)
        b={}
        for k in meta['Param']['BEV']['Event'][1].keys():
            bU=np.zeros(u.size)
            for iU in range(u.size):
                if u[iU]==0:
                    continue
                bU[iU]=meta['Param']['BEV']['Event'][u[iU]][k]
            b[k]=bU[inver]
        #t1=time.time()
        #meta[pNam]['Project']['Run Time Summary']['Test1']=meta[pNam]['Project']['Run Time Summary']['Test1']+t1-t0

        #----------------------------------------------------------------------
        # Adjust event-specific parameters to reflect time- and region-specific
        # fate of felled material and fate of removed fibre
        #----------------------------------------------------------------------
        #t0=time.time()
        # Index to harvesting
        iHarvest=np.where( (ID_Type==meta['LUT']['Event']['Harvest']) | (ID_Type==meta['LUT']['Event']['Harvest Salvage']) )[0]

        # Adjust fate of felled material parameters
        if iHarvest.size>0:

            # Index to time-dependent fate of felled materials
            iT_P=np.where(meta['Param']['BE']['Felled Fate']['Year']==meta[pNam]['Year'][iT])[0]

            # Simulations may exceed the timeframe of the felled fate parameters
            # If so, set to the last year
            if iT_P.size==0:
                iT_P=-1

            for k in meta['Param']['BEV']['Felled Fate'].keys():
                b[k][iHarvest]=meta['Param']['BEV']['Felled Fate'][k][iT_P,iHarvest]
        #t1=time.time()
        #meta[pNam]['Project']['Run Time Summary']['Test2']=meta[pNam]['Project']['Run Time Summary']['Test2']+t1-t0

        #----------------------------------------------------------------------
        # Define the amount of each pool that is affected by the event
        #----------------------------------------------------------------------

        # Fractions
        FracBiomassMerchAffected=b['BiomassMerch_Affected']*MortalityFactor
        FracBiomassNonMerchAffected=b['BiomassNonMerch_Affected']*MortalityFactor
        FracSnagsAffected=b['Snags_Affected']*MortalityFactor

        # Affected biomass carbon
        if meta[pNam]['Project']['Biomass Module']=='Sawtooth':
            Affected_StemMerch=vo['C_M_Tot'][iT,:,iEP['StemMerch']]
            Affected_StemNonMerch=vo['C_M_Tot'][iT,:,iEP['StemNonMerch']]
            Affected_Foliage=vo['C_M_Tot'][iT,:,iEP['Foliage']]
            Affected_Branch=vo['C_M_Tot'][iT,:,iEP['Branch']]
            Affected_Bark=vo['C_M_Tot'][iT,:,iEP['Bark']]
            Affected_RootCoarse=vo['C_M_Tot'][iT,:,iEP['RootCoarse']]
            Affected_RootFine=vo['C_M_Tot'][iT,:,iEP['RootFine']]
        else:
            Affected_StemMerch=FracBiomassMerchAffected*vo['C_Eco_Pools'][iT,:,iEP['StemMerch']]
            Affected_StemNonMerch=FracBiomassNonMerchAffected*vo['C_Eco_Pools'][iT,:,iEP['StemNonMerch']]
            Affected_Foliage=FracBiomassNonMerchAffected*vo['C_Eco_Pools'][iT,:,iEP['Foliage']]
            Affected_Branch=FracBiomassNonMerchAffected*vo['C_Eco_Pools'][iT,:,iEP['Branch']]
            Affected_Bark=FracBiomassNonMerchAffected*vo['C_Eco_Pools'][iT,:,iEP['Bark']]
            Affected_RootCoarse=FracBiomassNonMerchAffected*vo['C_Eco_Pools'][iT,:,iEP['RootCoarse']]
            Affected_RootFine=FracBiomassNonMerchAffected*vo['C_Eco_Pools'][iT,:,iEP['RootFine']]

        # Partition bark into merch and non-merch components
        Affected_BarkMerch=0.85*Affected_Bark
        Affected_BarkNonMerch=(1.00-0.85)*Affected_Bark

        # Sum up total affected non-merchantable biomass
        Affected_TotNonMerch=Affected_StemNonMerch+Affected_Branch+Affected_BarkNonMerch

        # Snags
        Affected_SnagStem=FracSnagsAffected*vo['C_Eco_Pools'][iT,:,iEP['SnagStem']]
        Affected_SnagBranch=FracSnagsAffected*vo['C_Eco_Pools'][iT,:,iEP['SnagBranch']]

        # All carbon
        Affected_All=Affected_StemMerch+Affected_StemNonMerch+Affected_Foliage+Affected_Branch+Affected_Bark+Affected_RootCoarse+Affected_RootFine

        # Live merch. stemwood volume
        Affected_VolumeStemMerchLive=FracBiomassMerchAffected*vo['V_MerchLive'][iT,:]

        # Dead merch. stemwood volume
        Affected_VolumeStemMerchDead=FracSnagsAffected*vo['V_MerchDead'][iT,:]

        #----------------------------------------------------------------------
        # Calculate mortality
        #----------------------------------------------------------------------
        #t0=time.time()

        # Total mortality
        vo['C_M_Dist'][iT,:]=vo['C_M_Dist'][iT,:]+Affected_All

        # Mortality by category (lumping stands together)
        uType=np.unique(ID_Type)
        for iType in range(len(uType)):
            if uType[iType]==0:
                continue
            ind=np.where(ID_Type==uType[iType])[0]
            Type=cbu.lut_n2s(meta['LUT']['Event'],uType[iType])[0]
            vo['C_M_ByAgent'][Type][iT,ind]=vo['C_M_ByAgent'][Type][iT,ind]+Affected_All[ind]

        #t1=time.time()
        #meta[pNam]['Project']['Run Time Summary']['Test4']=meta[pNam]['Project']['Run Time Summary']['Test4']+t1-t0

        #----------------------------------------------------------------------
        # Remove affected amount from each pool
        #----------------------------------------------------------------------
        #t0=time.time()
        if meta[pNam]['Project']['Biomass Module']!='Sawtooth':
            # Remove carbon from affected biomass pools
            vo['C_Eco_Pools'][iT,:,iEP['StemMerch']]=vo['C_Eco_Pools'][iT,:,iEP['StemMerch']]-Affected_StemMerch
            vo['C_Eco_Pools'][iT,:,iEP['StemNonMerch']]=vo['C_Eco_Pools'][iT,:,iEP['StemNonMerch']]-Affected_StemNonMerch
            vo['C_Eco_Pools'][iT,:,iEP['Foliage']]=vo['C_Eco_Pools'][iT,:,iEP['Foliage']]-Affected_Foliage
            vo['C_Eco_Pools'][iT,:,iEP['Branch']]=vo['C_Eco_Pools'][iT,:,iEP['Branch']]-Affected_Branch
            vo['C_Eco_Pools'][iT,:,iEP['Bark']]=vo['C_Eco_Pools'][iT,:,iEP['Bark']]-Affected_Bark
            vo['C_Eco_Pools'][iT,:,iEP['RootCoarse']]=vo['C_Eco_Pools'][iT,:,iEP['RootCoarse']]-Affected_RootCoarse
            vo['C_Eco_Pools'][iT,:,iEP['RootFine']]=vo['C_Eco_Pools'][iT,:,iEP['RootFine']]-Affected_RootFine

            # Remove carbon from snag pools
            vo['C_Eco_Pools'][iT,:,iEP['SnagStem']]=vo['C_Eco_Pools'][iT,:,iEP['SnagStem']]-Affected_SnagStem
            vo['C_Eco_Pools'][iT,:,iEP['SnagBranch']]=vo['C_Eco_Pools'][iT,:,iEP['SnagBranch']]-Affected_SnagBranch

            # Remove stemwood merch volume
            vo['V_MerchLive'][iT,:]=np.maximum(0,vo['V_MerchLive'][iT,:]-Affected_VolumeStemMerchLive)

        #----------------------------------------------------------------------
        # We have not explicity tracked a partition of snags into merch and
        # non merch, but these categories may be removed in different degrees
        # Move a proportion of snag stemwood to non-merch biomass
        #----------------------------------------------------------------------

        Transfer_SnagStemToNonMerch=meta['Param']['BEV']['Inter Pool Fluxes']['FractionSnagStemThatIsNonMerch']*Affected_SnagStem
        Affected_SnagStem=Affected_SnagStem-Transfer_SnagStemToNonMerch
        Affected_TotNonMerch=Affected_TotNonMerch+Transfer_SnagStemToNonMerch

        #----------------------------------------------------------------------
        # Carbon that is removed (ie sent to mills)
        #----------------------------------------------------------------------

        # Merch biomass to mill - of the total amount of biomass affected,
        # Add bark
        Biomass_Merch_Removed=b['BiomassMerch_Removed']*(Affected_StemMerch+Affected_BarkMerch)
        vo['C_ToMillMerch'][iT,:]=vo['C_ToMillMerch'][iT,:]+Biomass_Merch_Removed

        # NonMerch biomass to mill - of the total amount of biomass affected,
        # what fraction of non-merch biomass was sent to the mill?
        # - NonMerch = NonMerchStem + Foliage + Branch + Bark
        Biomass_NonMerch_Removed=b['BiomassNonMerch_Removed']*Affected_TotNonMerch
        vo['C_ToMillNonMerch'][iT,:]=vo['C_ToMillNonMerch'][iT,:]+Biomass_NonMerch_Removed

        # Snag stemwood to mill
        SnagStem_Removed=b['Snags_Removed']*Affected_SnagStem
        vo['C_ToMillSnagStem'][iT,:]=vo['C_ToMillSnagStem'][iT,:]+SnagStem_Removed

        # Snag branches to mill
        SnagBranch_Removed=b['Snags_Removed']*Affected_SnagBranch
        vo['C_ToMillNonMerch'][iT,:]=vo['C_ToMillNonMerch'][iT,:]+SnagBranch_Removed

        #----------------------------------------------------------------------
        # Volume that is removed (ie sent to mills)
        #----------------------------------------------------------------------

        # Conversion factor
        cf=(1/meta['Param']['BEV']['Biophysical']['Density Wood'])*(1/meta['Param']['BEV']['Biophysical']['Carbon Content Wood'])

        if meta[pNam]['Project']['Biomass Module']=='Sawtooth':
            vo['V_ToMillMerchLive'][iT,:]=cf*vo['C_ToMillMerch'][iT,:]
            vo['V_ToMillMerchDead'][iT,:]=cf*vo['C_ToMillSnagStem'][iT,:]
        else:
            vo['V_ToMillMerchLive'][iT,:]=vo['V_ToMillMerchLive'][iT,:]+b['BiomassMerch_Removed']*Affected_VolumeStemMerchLive
            vo['V_ToMillMerchDead'][iT,:]=vo['V_ToMillMerchDead'][iT,:]+b['Snags_Removed']*Affected_VolumeStemMerchDead

        # Total merch stemwood volume removed
        vo['V_ToMillMerchTotal'][iT,:]=vo['V_ToMillMerchLive'][iT,:]+vo['V_ToMillMerchDead'][iT,:]

        # Non-mech volume removed
        vo['V_ToMillNonMerch'][iT,:]=cf*vo['C_ToMillNonMerch'][iT,:]

        #----------------------------------------------------------------------
        # Carbon that is left dispersed on site (after felling or wind storms)
        #----------------------------------------------------------------------

        # Foliage transferred directly to very fast litter
        vo['C_Eco_Pools'][iT,:,iEP['LitterVF']]=vo['C_Eco_Pools'][iT,:,iEP['LitterVF']]+Affected_Foliage

        # Stem, branch and bark carbon transfered to medium and fast litter pools
        vo['C_Eco_Pools'][iT,:,iEP['LitterM']]=vo['C_Eco_Pools'][iT,:,iEP['LitterM']]+b['BiomassMerch_LeftOnSite']*Affected_StemMerch
        vo['C_Eco_Pools'][iT,:,iEP['LitterM']]=vo['C_Eco_Pools'][iT,:,iEP['LitterM']]+b['BiomassNonMerch_LeftOnSite']*Affected_StemNonMerch
        vo['C_Eco_Pools'][iT,:,iEP['LitterF']]=vo['C_Eco_Pools'][iT,:,iEP['LitterF']]+b['BiomassMerch_LeftOnSite']*Affected_Branch
        vo['C_Eco_Pools'][iT,:,iEP['LitterF']]=vo['C_Eco_Pools'][iT,:,iEP['LitterF']]+b['BiomassMerch_LeftOnSite']*Affected_Bark

        vo['C_Eco_Pools'][iT,:,iEP['LitterM']]=vo['C_Eco_Pools'][iT,:,iEP['LitterM']]+b['Snags_LeftOnSite']*Affected_SnagStem
        vo['C_Eco_Pools'][iT,:,iEP['LitterF']]=vo['C_Eco_Pools'][iT,:,iEP['LitterF']]+b['Snags_LeftOnSite']*Affected_SnagBranch

        # Roots transferred directly to DOM
        vo['C_Eco_Pools'][iT,:,iEP['LitterF']]=vo['C_Eco_Pools'][iT,:,iEP['LitterF']]+0.5*Affected_RootCoarse
        vo['C_Eco_Pools'][iT,:,iEP['SoilF']]=vo['C_Eco_Pools'][iT,:,iEP['SoilF']]+0.5*Affected_RootCoarse
        vo['C_Eco_Pools'][iT,:,iEP['LitterVF']]=vo['C_Eco_Pools'][iT,:,iEP['LitterVF']]+0.5*Affected_RootFine
        vo['C_Eco_Pools'][iT,:,iEP['SoilVF']]=vo['C_Eco_Pools'][iT,:,iEP['SoilVF']]+0.5*Affected_RootFine

        #----------------------------------------------------------------------
        # Carbon that is piled
        # (a small fraction of piled wood is collected for firewood)
        #----------------------------------------------------------------------

        PiledStemMerch=b['BiomassMerch_Piled']*Affected_StemMerch
        PiledSnagStem=b['Snags_Piled']*Affected_SnagStem

        vo['C_Eco_Pools'][iT,:,iEP['PiledStemMerch']]=vo['C_Eco_Pools'][iT,:,iEP['PiledStemMerch']]+(1.0-meta['Param']['BEV']['HWP']['PiledStemwoodToFirewoodDom'])*PiledStemMerch
        vo['C_Eco_Pools'][iT,:,iEP['PiledStemNonMerch']]=vo['C_Eco_Pools'][iT,:,iEP['PiledStemNonMerch']]+b['BiomassNonMerch_Piled']*Affected_StemNonMerch
        vo['C_Eco_Pools'][iT,:,iEP['PiledBranch']]=vo['C_Eco_Pools'][iT,:,iEP['PiledBranch']]+b['BiomassNonMerch_Piled']*Affected_Branch
        vo['C_Eco_Pools'][iT,:,iEP['PiledBark']]=vo['C_Eco_Pools'][iT,:,iEP['PiledBark']]+b['BiomassNonMerch_Piled']*Affected_Bark
        vo['C_Eco_Pools'][iT,:,iEP['PiledSnagStem']]=vo['C_Eco_Pools'][iT,:,iEP['PiledSnagStem']]+(1.0-meta['Param']['BEV']['HWP']['PiledStemwoodToFirewoodDom'])*PiledSnagStem
        vo['C_Eco_Pools'][iT,:,iEP['PiledSnagBranch']]=vo['C_Eco_Pools'][iT,:,iEP['PiledSnagBranch']]+b['Snags_Piled']*Affected_SnagBranch
        vo['C_ToFirewoodDom'][iT,:]=vo['C_ToFirewoodDom'][iT,:]+meta['Param']['BEV']['HWP']['PiledStemwoodToFirewoodDom']*(PiledStemMerch+PiledSnagStem)

        # Affected carbon
        StemMerch=b['PiledStemMerch_Burned']*vo['C_Eco_Pools'][iT,:,iEP['PiledStemMerch']]
        StemNonMerch=b['PiledStemNonMerch_Burned']*vo['C_Eco_Pools'][iT,:,iEP['PiledStemNonMerch']]
        Branch=b['PiledBranch_Burned']*vo['C_Eco_Pools'][iT,:,iEP['PiledBranch']]
        Bark=b['PiledBark_Burned']*vo['C_Eco_Pools'][iT,:,iEP['PiledBark']]
        SnagStem=b['PiledSnagStem_Burned']*vo['C_Eco_Pools'][iT,:,iEP['PiledSnagStem']]
        SnagBranch=b['PiledSnagBranch_Burned']*vo['C_Eco_Pools'][iT,:,iEP['PiledSnagBranch']]
        Total_Burned=StemMerch+StemNonMerch+Branch+Bark+SnagStem+SnagBranch
        NonMerch_Burned=StemNonMerch+Branch+Bark+SnagBranch

        # Remove affected carbon
        vo['C_Eco_Pools'][iT,:,iEP['PiledStemMerch']]=vo['C_Eco_Pools'][iT,:,iEP['PiledStemMerch']]-StemMerch
        vo['C_Eco_Pools'][iT,:,iEP['PiledStemNonMerch']]=vo['C_Eco_Pools'][iT,:,iEP['PiledStemNonMerch']]-StemNonMerch
        vo['C_Eco_Pools'][iT,:,iEP['PiledBranch']]=vo['C_Eco_Pools'][iT,:,iEP['PiledBranch']]-Branch
        vo['C_Eco_Pools'][iT,:,iEP['PiledBark']]=vo['C_Eco_Pools'][iT,:,iEP['PiledBark']]-Bark
        vo['C_Eco_Pools'][iT,:,iEP['PiledSnagStem']]=vo['C_Eco_Pools'][iT,:,iEP['PiledSnagStem']]-SnagStem
        vo['C_Eco_Pools'][iT,:,iEP['PiledSnagBranch']]=vo['C_Eco_Pools'][iT,:,iEP['PiledSnagBranch']]-SnagBranch

        # Add to fire emissions
        vo['C_E_OpenBurningAsCO2'][iT,:]=vo['C_E_OpenBurningAsCO2'][iT,:]+meta['Param']['BEV']['Biophysical']['CombFrac_CO2']*Total_Burned
        vo['C_E_OpenBurningAsCH4'][iT,:]=vo['C_E_OpenBurningAsCH4'][iT,:]+meta['Param']['BEV']['Biophysical']['CombFrac_CH4']*Total_Burned
        vo['C_E_OpenBurningAsCO'][iT,:]=vo['C_E_OpenBurningAsCO'][iT,:]+meta['Param']['BEV']['Biophysical']['CombFrac_CO']*Total_Burned

        # If it is slashpile burning, track it accordingly
        vo['C_ToSlashpileBurnTot'][iT,:]=vo['C_ToSlashpileBurnTot'][iT,:]+Total_Burned
        vo['C_ToSlashpileBurnNonMerch'][iT,:]=vo['C_ToSlashpileBurnNonMerch'][iT,:]+NonMerch_Burned

        #----------------------------------------------------------------------
        # Carbon that is moved from biomass to snags
        #----------------------------------------------------------------------

        # Merch biomass that is killed by wildfire
        vo['C_Eco_Pools'][iT,:,iEP['SnagStem']]=vo['C_Eco_Pools'][iT,:,iEP['SnagStem']]+b['BiomassMerch_ToSnag']*(Affected_StemMerch+Affected_BarkMerch)

        # Non-merch biomass that is killed by wildfire
        vo['C_Eco_Pools'][iT,:,iEP['SnagBranch']]=vo['C_Eco_Pools'][iT,:,iEP['SnagBranch']]+b['BiomassNonMerch_ToSnag']*(Affected_StemNonMerch+Affected_BarkNonMerch)

        #----------------------------------------------------------------------
        # Biomass and DOM burned in wildfire
        #----------------------------------------------------------------------

        Biomass_Merch_Burned=b['BiomassMerch_Burned']*(Affected_StemMerch+Affected_BarkMerch)
        Biomass_NonMerch_Burned=b['BiomassNonMerch_Burned']*Affected_TotNonMerch
        SnagStem_Burned=b['Snags_Burned']*Affected_SnagStem
        SnagBranch_Burned=b['Snags_Burned']*Affected_SnagBranch
        Total_Burned=Biomass_Merch_Burned+Biomass_NonMerch_Burned+SnagStem_Burned+SnagBranch_Burned

        vo['C_E_WildfireAsCO2'][iT,:]=vo['C_E_WildfireAsCO2'][iT,:]+meta['Param']['BEV']['Biophysical']['CombFrac_CO2']*Total_Burned
        vo['C_E_WildfireAsCH4'][iT,:]=vo['C_E_WildfireAsCH4'][iT,:]+meta['Param']['BEV']['Biophysical']['CombFrac_CH4']*Total_Burned
        vo['C_E_WildfireAsCO'][iT,:]=vo['C_E_WildfireAsCO'][iT,:]+meta['Param']['BEV']['Biophysical']['CombFrac_CO']*Total_Burned

        #----------------------------------------------------------------------
        # Update total stemwood volume
        #----------------------------------------------------------------------

        vo['V_MerchDead'][iT,:]=(1/meta['Param']['BEV']['Biophysical']['Density Wood'])*(1/meta['Param']['BEV']['Biophysical']['Carbon Content Wood'])*vo['C_Eco_Pools'][iT,:,iEP['SnagStem']]
        vo['V_MerchTotal'][iT,:]=vo['V_MerchLive'][iT,:]+vo['V_MerchDead'][iT,:]
        #t1=time.time()
        #meta[pNam]['Project']['Run Time Summary']['Test8']=meta[pNam]['Project']['Run Time Summary']['Test8']+t1-t0
        #----------------------------------------------------------------------
        # Update stand age
        #----------------------------------------------------------------------
        #t0=time.time()
        if (meta[pNam]['Project']['Biomass Module']!='Sawtooth') & (meta[pNam]['Project']['Partial Mortality Affects Age']=='On'):

            # List of exception event types (that will remain at the same age)
            Exceptions_to_Partial_Mortality=[  ]

            # Index to types where age will change
            ind=np.where( (np.isin(ID_Type,Exceptions_to_Partial_Mortality)==False) )[0]

            # Assume oldest trees were most affected, reduce age in prportion
            # with mortality rate
            # Not always realistic, but see how net growth is affected.
            vo['A'][iT,ind]=vo['A'][iT,ind]*(1-MortalityFactor)

        # Ensure stand-replacing events reset stand age to 0
        vo['A'][iT,(MortalityFactor==1)]=0

        # Ensure planting resets to age 0
        vo['A'][iT,(ID_Type==meta['LUT']['Event']['Planting']) | (ID_Type==meta['LUT']['Event']['Direct Seeding'])]=0
        #t1=time.time()
        #meta[pNam]['Project']['Run Time Summary']['Test7']=meta[pNam]['Project']['Run Time Summary']['Test7']+t1-t0

        #----------------------------------------------------------------------
        # Transition to new growth curve
        #----------------------------------------------------------------------

        # Initialize a flag that indicates whether the growth curve changes
        #flg_gc_change=0
        #if vi['EH'][iS]['ID_GrowthCurveM'][indDist[iDist]]!=vi['ID_GCA'][iS]:
        #    flg_gc_change=1
        #t0=time.time()

        # Only applies to BatchTIPSY
        if meta[pNam]['Project']['Biomass Module']=='BatchTIPSY':

            for iGC in range(meta['Modules']['GYM']['N Growth Curves']):

                # Don't alter growth curve for fertilization - index to non-fertilization events
                ind=np.where( (vi['EC']['ID Growth Curve'][iT,:,iE]==meta['Modules']['GYM']['ID GC Unique'][iGC]) & (ID_Type!=meta['LUT']['Event']['Fertilization Aerial']) | \
                              (vi['EC']['ID Growth Curve'][iT,:,iE]==meta['Modules']['GYM']['ID GC Unique'][iGC]) & (ID_Type!=meta['LUT']['Event']['Fertilization Hand']) )[0]

                if ind.size>0:
                    # Notes: This crashes often when the GY model has been accidently run with the wrong project info
                    # - because the project path in BatchTipsy is wrong.
                    vi['GC']['Active'][:,ind,:]=vi['GC'][ meta['Modules']['GYM']['ID GC Unique'][iGC] ][:,ind,:].astype(float)*meta['Modules']['GYM']['Scale Factor']
                    vi['GC']['ID_GCA'][ind]=int(meta['Modules']['GYM']['ID GC Unique'][iGC])
        #t1=time.time()
        #meta[pNam]['Project']['Run Time Summary']['Test5']=meta[pNam]['Project']['Run Time Summary']['Test5']+t1-t0

        #----------------------------------------------------------------------
        # Impose regen failure
        #----------------------------------------------------------------------

        # Only applies to BatchTIPSY
        if meta[pNam]['Project']['Biomass Module']=='BatchTIPSY':
            iFailure=np.where(ID_Type==meta['LUT']['Event']['Regen Failure'])[0]
            if iFailure.size>0:
                vi['GC']['Active'][:,iFailure,:]=0
                vi['GC'][1][:,iFailure,:]=0

        #------------------------------------------------------------------
        # Update net growth (in response to lethal events)
        #
        # Partial disturbances likely have a lasting impact on gross growth
        # by removing growing stock. Without ad hoc
        # adjustment, this will not be reflected in the growth curve from TIPSY.
        #
        # This will apply to instances where Mortality=100% only if the
        # user does not specify a change in growth curve.
        # With this functionality, the user has the options of:
        #  1) Changing the growth curve
        #  2) Assuming a gradual recovery to the previous curve, as defined
        #     by the default recovery half life for that disturbance type.
        #------------------------------------------------------------------

        # *** This is not up to date - will crash if attempted - needs update. ***
        flg=0

        if (flg==1) & (meta[pNam]['Project']['Biomass Module']!='Sawtooth'):

            # Half life
            hl=meta['Param']['BEV']['bDist_GrowthRecovery_HL']

            # Only proceed if:
            #   1) the disturbance has a lasting growth impact/recovery
            #   2) the user is not specifying a change in growth curve
            if (hl>0) & (flg_gc_change==0):

                # Extract net growth for active growth curve
                NetGrowth=vi['GCA'][:,iS,:].copy().astype(float)*meta['Modules']['GYM']['Scale Factor']

                # Growth pre-event
                G_pre=NetGrowth.copy()

                # Growth post-event
                G_post=(1-Biomass_Affected_Frac)*NetGrowth.copy()

                # Difference in growth
                dG=G_pre-G_post

                # Age vector
                A=np.arange(1,meta['Modules']['GYM']['BatchTIPSY Maximum Age']+1,1)

                # Time since disturbance
                TSD=np.tile(A-vo['A'][iT,iS],(G_pre.shape[1],1)).T

                # Relative effect
                fTSD=1/(1+np.exp(-0.5*(TSD-hl)))
                fTSD[0:int(vo['A'][iT,iS]),:]=0

                # Growth recovery
                G_Recovery=fTSD*dG

                # New estimate of net growth
                NetGrowthNew=G_post+G_Recovery

                # Add back to dictionary
                NetGrowthNew=NetGrowthNew*meta['Modules']['GYM']['Scale Factor']
                NetGrowthNew=NetGrowthNew.astype(np.int16)

                vi['GCA'][:,iS,:]=NetGrowthNew

        #----------------------------------------------------------------------
        # Apply growth factors (in response to non-lethal events)
        # Growth factor should comes in DMEC as percent effect
        # 10 = 10% increase, -10 = 10% decrease
        #----------------------------------------------------------------------

        #t0=time.time()
        flg=1
        if (flg==1) & (meta[pNam]['Project']['Biomass Module']!='Sawtooth'):

            GrowthFactor=vi['EC']['Growth Factor'][iT,:,iE].copy()

            indAdj=np.where( (GrowthFactor!=-999) & (GrowthFactor!=0) )[0]

            if (indAdj.size>0):

                Age=np.arange(0,meta['Modules']['GYM']['BatchTIPSY Maximum Age']+1,1)

                # Convert growth factor to respoonse ratio
                GrowthResponseRatio=GrowthFactor.astype(float)/100+1.0

                NetGrowth=vi['GC']['Active'].copy()

                for iAdj in range(indAdj.size):

                    if ID_Type[iAdj]==meta['LUT']['Event']['IDW']:

                        # Only a temporary change in growth (see severity class table)
                        # *** I have not checked to see if this is working properly ***
                        ResponsePeriod=6

                        iResponse=np.where( (Age>=vo['A'][iT,iAdj]) & (Age<=vo['A'][iT,iAdj]+ResponsePeriod) )[0]

                        NetGrowth[iResponse,iAdj,:]=GrowthResponseRatio[iAdj]*NetGrowth[iResponse,iAdj,:]

                    else:

                        # A permanent change in growth
                        NetGrowth[:,iAdj,:]=GrowthResponseRatio[iAdj]*NetGrowth[:,iAdj,:]

                vi['GC']['Active']=NetGrowth

        #t1=time.time()
        #meta[pNam]['Project']['Run Time Summary']['Test6']=meta[pNam]['Project']['Run Time Summary']['Test6']+t1-t0

    #--------------------------------------------------------------------------
    # Relative gravimetric mortality (%)
    #--------------------------------------------------------------------------

    for k in meta['LUT']['Event'].keys():
        vo['C_M_Pct_ByAgent'][k][iT,:]=vo['C_M_ByAgent'][k][iT,:]/C_Biomass_t0*100

    #--------------------------------------------------------------------------
    # Aerial nutrient application events
    #--------------------------------------------------------------------------

    # Generate index to stands that were fertilized
    meta['Modules']['Nutrient Management']['iApplication']=np.where(flag_nutrient_application==1)[0]

    if meta['Modules']['Nutrient Management']['iApplication'].size>0:

        # Adjust net growth of aboveground biomass
        vi,vo,meta=napp.UpdateStatus(vi,vo,iT,meta,'AbovegroundNetGrowth')

        # Adjust emissions
        vi,vo,meta=napp.UpdateStatus(vi,vo,iT,meta,'Emissions')

    return vo,vi

#%% Harvested wood products sector

def HWP_Update21(meta,pNam,iT,iBat,vi,vo):

    #--------------------------------------------------------------------------
    # Index to time-dependent parameters
    #--------------------------------------------------------------------------

    # This will be the same between end use and removed fate because the time vectors
    # are the same.
    # *** Don't change one time vector without changing the other. ***
    iT_P=np.where(meta['Param']['BE']['HWP End Use']['Year']==meta[pNam]['Year'][iT])[0]

    # Removed fate
    bRF={}
    for k in meta['Param']['BEV']['Removed Fate'].keys():
        if iT_P.size!=0:
            bRF[k]=meta['Param']['BEV']['Removed Fate'][k][iT_P,:].flatten()
        else:
            # Simulations may extend beyond the 2100 limit of the time-dependent HWP variables
            bRF[k]=meta['Param']['BEV']['Removed Fate'][k][-1,:].flatten()

    # End uses
    bEU={}
    for k in meta['Param']['BEV']['HWP End Use'].keys():
        if iT_P.size!=0:
            bEU[k]=meta['Param']['BEV']['HWP End Use'][k][iT_P,:].flatten()
        else:
            # Simulations may extend beyond the 2100 limit of the time-dependent HWP variables
            bEU[k]=meta['Param']['BEV']['HWP End Use'][k][-1,:].flatten()

    #--------------------------------------------------------------------------
    # Carbon transferred from forest to mills or direct to end-uses
    #--------------------------------------------------------------------------

    C_ToPulpMill=bRF['RemovedMerchToPulpMill']*vo['C_ToMillMerch'][iT,:] + \
                 bRF['RemovedNonMerchToPulpMill']*vo['C_ToMillNonMerch'][iT,:] + \
                 bRF['RemovedSnagStemToPulpMill']*vo['C_ToMillSnagStem'][iT,:]

    C_ToPelletMill=bRF['RemovedMerchToPelletMill']*vo['C_ToMillMerch'][iT,:] + \
                 bRF['RemovedNonMerchToPelletMill']*vo['C_ToMillNonMerch'][iT,:] + \
                 bRF['RemovedSnagStemToPelletMill']*vo['C_ToMillSnagStem'][iT,:]

    C_ToLumberMill=bRF['RemovedMerchToLumberMill']*vo['C_ToMillMerch'][iT,:] + \
                 bRF['RemovedNonMerchToLumberMill']*vo['C_ToMillNonMerch'][iT,:] + \
                 bRF['RemovedSnagStemToLumberMill']*vo['C_ToMillSnagStem'][iT,:]

    C_ToPlywoodMill=bRF['RemovedMerchToPlywoodMill']*vo['C_ToMillMerch'][iT,:] + \
                 bRF['RemovedNonMerchToPlywoodMill']*vo['C_ToMillNonMerch'][iT,:] + \
                 bRF['RemovedSnagStemToPlywoodMill']*vo['C_ToMillSnagStem'][iT,:]

    C_ToOSBMill=bRF['RemovedMerchToOSBMill']*vo['C_ToMillMerch'][iT,:] + \
                 bRF['RemovedNonMerchToOSBMill']*vo['C_ToMillNonMerch'][iT,:] + \
                 bRF['RemovedSnagStemToOSBMill']*vo['C_ToMillSnagStem'][iT,:]

    C_ToMDFMill=bRF['RemovedMerchToMDFMill']*vo['C_ToMillMerch'][iT,:] + \
                 bRF['RemovedNonMerchToMDFMill']*vo['C_ToMillNonMerch'][iT,:] + \
                 bRF['RemovedSnagStemToMDFMill']*vo['C_ToMillSnagStem'][iT,:]

    C_ToLogExport=bRF['RemovedMerchToLogExport']*vo['C_ToMillMerch'][iT,:] + \
                 bRF['RemovedNonMerchToLogExport']*vo['C_ToMillNonMerch'][iT,:] + \
                 bRF['RemovedSnagStemToLogExport']*vo['C_ToMillSnagStem'][iT,:]

    C_ToPowerGrid=bRF['RemovedMerchToIPP']*vo['C_ToMillMerch'][iT,:] + \
                 bRF['RemovedNonMerchToIPP']*vo['C_ToMillNonMerch'][iT,:] + \
                 bRF['RemovedSnagStemToIPP']*vo['C_ToMillSnagStem'][iT,:]

    C_ToFirewoodCollection=bRF['RemovedMerchToFirewood']*vo['C_ToMillMerch'][iT,:] + \
                 bRF['RemovedNonMerchToFirewood']*vo['C_ToMillNonMerch'][iT,:] + \
                 bRF['RemovedSnagStemToFirewood']*vo['C_ToMillSnagStem'][iT,:]

    #C_ChipperMill=meta['Param']['BEV']['Removed Fate']['RemovedMerchToChipperMill']*vo['C_ToMillMerch'][iT,:] + \
    #             meta['Param']['BEV']['Removed Fate']['RemovedNonMerchToChipperMill']*vo['C_ToMillNonMerch'][iT,:] + \
    #             meta['Param']['BEV']['Removed Fate']['RemovedSnagStemToChipperMill']*vo['C_ToMillSnagStem'][iT,:]

    #C_PolePostMill=bRF['RemovedMerchToPolePostMill']*vo['C_ToMillMerch'][iT,:] + \
    #             bRF['RemovedNonMerchToPolePostMill']*vo['C_ToMillNonMerch'][iT,:] + \
    #             bRF['RemovedSnagStemToPolePostMill']*vo['C_ToMillSnagStem'][iT,:]

    #C_ShakeShingleMill=bRF['RemovedMerchToShakeShingleMill']*vo['C_ToMillMerch'][iT,:] + \
    #             bRF['RemovedNonMerchToShakeShingleMill']*vo['C_ToMillNonMerch'][iT,:] + \
    #             bRF['RemovedSnagStemToShakeShingleMill']*vo['C_ToMillSnagStem'][iT,:]

    #--------------------------------------------------------------------------
    # Carbon transferred from mill to mill
    #--------------------------------------------------------------------------

    C_ToPulpMill=C_ToPulpMill+C_ToLumberMill*meta['Param']['BEV']['HWP']['LumberMillToPulpMill']

    C_ToMDFMill=C_ToMDFMill+C_ToLumberMill*meta['Param']['BEV']['HWP']['LumberMillToMDFMill']

    C_ToPelletMill=C_ToPelletMill+C_ToLumberMill*meta['Param']['BEV']['HWP']['LumberMillToPelletMill']

    #--------------------------------------------------------------------------
    # Carbon transferred from LumberMill logs to log Export
    #--------------------------------------------------------------------------

    C_ToLogExport=C_ToLogExport+bEU['LumberMillToLogExport']*C_ToLumberMill

    #--------------------------------------------------------------------------
    # Carbon transferred from log Export to firewood
    #--------------------------------------------------------------------------

    C_ToFirewoodExport=meta['Param']['BEV']['HWP']['LogExportToFirewood']*C_ToLogExport

    #--------------------------------------------------------------------------
    # Carbon transferred to paper
    #--------------------------------------------------------------------------

    C_ToPaper=meta['Param']['BEV']['HWP']['PulpMillToPaper']*C_ToPulpMill

    #--------------------------------------------------------------------------
    # Carbon transferred to pulp effluent
    #--------------------------------------------------------------------------

    C_ToPulpEffluent=meta['Param']['BEV']['HWP']['PulpMillToEffluent']*C_ToPulpMill

    #--------------------------------------------------------------------------
    # Carbon transferred to pellets
    #--------------------------------------------------------------------------

    C_ToPelletExport=meta['Param']['BEV']['HWP']['PelletMillToPelletExport']*C_ToPelletMill
    C_ToPelletDomGrid=meta['Param']['BEV']['HWP']['PelletMillToDomGrid']*C_ToPelletMill
    C_ToPelletDomRNG=meta['Param']['BEV']['HWP']['PelletMillToDomRNG']*C_ToPelletMill

    #--------------------------------------------------------------------------
    # Carbon transferred to facility power
    #--------------------------------------------------------------------------

    C_ToPowerFacilityDom=bEU['LumberMillToPowerFacility']*C_ToLumberMill + \
        meta['Param']['BEV']['HWP']['PulpMillToPowerFacility']*C_ToPulpMill + \
        bEU['PlywoodMillToPowerFacility']*C_ToPlywoodMill + \
        bEU['OSBMillToPowerFacility']*C_ToOSBMill

    C_ToPowerFacilityExport=meta['Param']['BEV']['HWP']['LogExportToPowerFacility']*C_ToLogExport

    #--------------------------------------------------------------------------
    # Carbon transferred to grid by independent power producers
    #--------------------------------------------------------------------------

    C_ToPowerGrid=C_ToPowerGrid+bEU['LumberMillToIPP']*C_ToLumberMill + \
        meta['Param']['BEV']['HWP']['PulpMillToIPP']*C_ToPulpMill + \
        bEU['PlywoodMillToIPP']*C_ToPlywoodMill + \
        bEU['OSBMillToIPP']*C_ToOSBMill

    #--------------------------------------------------------------------------
    # Log-size effect
    #--------------------------------------------------------------------------

    #    bLogSizeEffect=0.02
    #
    #    iLSE=np.where(vo['LogSizeEnhancement'][iT,:]>0)[0]
    #
    #    if iLSE.size>0:
    #
    #        C_Transfer=bLogSizeEffect*C_Paper[iLSE]
    #        C_Paper[iLSE]=C_Paper[iLSE]-C_Transfer
    #        C_LumberMill[iLSE]=C_LumberMill[iLSE]+C_Transfer
    #
    #        C_Transfer=bLogSizeEffect*C_Pellet[iLSE]
    #        C_Pellet[iLSE]=C_Pellet[iLSE]-C_Transfer
    #        C_LumberMill[iLSE]=C_LumberMill[iLSE]+C_Transfer
    #
    #        C_Transfer=bLogSizeEffect*C_PowerFacilityDom[iLSE]
    #        C_PowerFacilityDom[iLSE]=C_PowerFacilityDom[iLSE]-C_Transfer
    #        C_LumberMill[iLSE]=C_LumberMill[iLSE]+C_Transfer
    #
    #        C_Transfer=bLogSizeEffect*C_PowerFacilityFor[iLSE]
    #        C_PowerFacilityFor[iLSE]=C_PowerFacilityFor[iLSE]-C_Transfer
    #        C_LumberMill[iLSE]=C_LumberMill[iLSE]+C_Transfer
    #
    #        C_Transfer=bLogSizeEffect*C_PowerGrid[iLSE]
    #        C_PowerGrid[iLSE]=C_PowerGrid[iLSE]-C_Transfer
    #        C_LumberMill[iLSE]=C_LumberMill[iLSE]+C_Transfer

    #--------------------------------------------------------------------------
    # Production of single-family homes
    #--------------------------------------------------------------------------

    C_LumberMillToSFH=C_ToLumberMill*bEU['LumberMillToSFH']
    C_PlywoodMillToSFH=C_ToPlywoodMill*bEU['PlywoodMillToSFH']
    C_OSBMillToSFH=C_ToOSBMill*bEU['OSBMillToSFH']
    C_MDFMillToSFH=C_ToMDFMill*bEU['MDFMillToSFH']
    C_LogExportToSFH=C_ToLogExport*meta['Param']['BEV']['HWP']['LogExportToSFH']
    #C_ShakeShingleToSFH=C_ShakeShingleMill*meta['Param']['BEV']['HWP']['ShakeShingleMillToSFH']

    #--------------------------------------------------------------------------
    # Production of multi-family homes
    #--------------------------------------------------------------------------

    C_LumberMillToMFH=C_ToLumberMill*bEU['LumberMillToMFH']
    C_PlywoodMillToMFH=C_ToPlywoodMill*bEU['PlywoodMillToMFH']
    C_OSBMillToMFH=C_ToOSBMill*bEU['OSBMillToMFH']
    C_MDFMillToMFH=C_ToMDFMill*bEU['MDFMillToMFH']
    C_LogExportToMFH=C_ToLogExport*meta['Param']['BEV']['HWP']['LogExportToMFH']
    #C_ShakeShingleToMFH=C_ShakeShingleMill*meta['Param']['BEV']['HWP']['ShakeShingleMillToMFH']

    #--------------------------------------------------------------------------
    # Production of commercial buildings
    #--------------------------------------------------------------------------

    C_LumberMillToCom=C_ToLumberMill*bEU['LumberMillToCom']
    C_PlywoodMillToCom=C_ToPlywoodMill*bEU['PlywoodMillToCom']
    C_OSBMillToCom=C_ToOSBMill*bEU['OSBMillToCom']
    C_MDFMillToCom=C_ToMDFMill*bEU['MDFMillToCom']
    C_LogExportToCom=C_ToLogExport*meta['Param']['BEV']['HWP']['LogExportToCom']
    #C_ShakeShingleToCom=C_ShakeShingleMill*meta['Param']['BEV']['HWP']['ShakeShingleMillToCom']

    #--------------------------------------------------------------------------
    # Production of furniture
    #--------------------------------------------------------------------------

    C_LumberMillToFurn=C_ToLumberMill*bEU['LumberMillToFurn']
    C_PlywoodMillToFurn=C_ToPlywoodMill*bEU['PlywoodMillToFurn']
    C_OSBMillToFurn=C_ToOSBMill*bEU['OSBMillToFurn']
    C_MDFMillToFurn=C_ToMDFMill*bEU['MDFMillToFurn']
    C_LogExportToFurn=C_ToLogExport*meta['Param']['BEV']['HWP']['LogExportToFurn']

    #--------------------------------------------------------------------------
    # Production of shipping containers
    #--------------------------------------------------------------------------

    C_LumberMillToShip=C_ToLumberMill*bEU['LumberMillToShip']
    C_PlywoodMillToShip=C_ToPlywoodMill*bEU['PlywoodMillToShip']
    C_OSBMillToShip=C_ToOSBMill*bEU['OSBMillToShip']
    C_MDFMillToShip=C_ToMDFMill*bEU['MDFMillToShip']
    C_LogExportToShip=C_ToLogExport*meta['Param']['BEV']['HWP']['LogExportToShip']

    #--------------------------------------------------------------------------
    # Production of repairs
    #--------------------------------------------------------------------------

    C_LumberMillToRepairs=C_ToLumberMill*bEU['LumberMillToRepairs']
    C_PlywoodMillToRepairs=C_ToPlywoodMill*bEU['PlywoodMillToRepairs']
    C_OSBMillToRepairs=C_ToOSBMill*bEU['OSBMillToRepairs']
    C_MDFMillToRepairs=C_ToMDFMill*bEU['MDFMillToRepairs']
    C_LogExportToRepairs=C_ToLogExport*meta['Param']['BEV']['HWP']['LogExportToRepairs']

    #--------------------------------------------------------------------------
    # Production of other
    #--------------------------------------------------------------------------

    C_LumberMillToOther=C_ToLumberMill*bEU['LumberMillToOther']
    C_PlywoodMillToOther=C_ToPlywoodMill*bEU['PlywoodMillToOther']
    C_OSBMillToOther=C_ToOSBMill*bEU['OSBMillToOther']
    C_MDFMillToOther=C_ToMDFMill*bEU['MDFMillToOther']
    C_LogExportToOther=C_ToLogExport*meta['Param']['BEV']['HWP']['LogExportToOther']
    #C_PolePostMillToOther=PolePostMill*meta['Param']['BEV']['HWP']['PolePostMillToOther']

    #--------------------------------------------------------------------------
    # Track sales (for economic modelling)
    # *** Conservation of mass is not being retained in the below set of variables
    # - log Export retain the full amount of log Export, but some of it was
    # also added to foreign power facility and foreign firewood.
    # Economics only consider log Export, substitution effects only consider
    # foreign power facility ad foreign firewood. ***
    #--------------------------------------------------------------------------

    vo['C_ToLumber'][iT,:]=C_LumberMillToSFH+C_LumberMillToMFH+C_LumberMillToCom+C_LumberMillToFurn+C_LumberMillToShip+C_LumberMillToRepairs+C_LumberMillToOther

    vo['C_ToPlywood'][iT,:]=C_PlywoodMillToSFH+C_PlywoodMillToMFH+C_PlywoodMillToCom+C_PlywoodMillToFurn+C_PlywoodMillToShip+C_PlywoodMillToRepairs+C_PlywoodMillToOther

    vo['C_ToOSB'][iT,:]=C_OSBMillToSFH+C_OSBMillToMFH+C_OSBMillToCom+C_OSBMillToFurn+C_OSBMillToShip+C_OSBMillToRepairs+C_OSBMillToOther

    vo['C_ToMDF'][iT,:]=C_MDFMillToSFH+C_MDFMillToMFH+C_MDFMillToCom+C_MDFMillToFurn+C_MDFMillToShip+C_MDFMillToRepairs+C_MDFMillToOther

    vo['C_ToPaper'][iT,:]=C_ToPaper

    vo['C_ToPowerFacilityDom'][iT,:]=C_ToPowerFacilityDom

    vo['C_ToPowerFacilityExport'][iT,:]=C_ToPowerFacilityExport

    vo['C_ToPowerGrid'][iT,:]=C_ToPowerGrid

    vo['C_ToPelletExport'][iT,:]=C_ToPelletExport
    vo['C_ToPelletDomGrid'][iT,:]=C_ToPelletDomGrid
    vo['C_ToPelletDomRNG'][iT,:]=C_ToPelletDomRNG

    # *** There is firewood taken directly from forest ecosystems (see events function) ***
    vo['C_ToFirewoodDom'][iT,:]=vo['C_ToFirewoodDom'][iT,:]+C_ToFirewoodCollection

    vo['C_ToFirewoodExport'][iT,:]=C_ToFirewoodExport

    vo['C_ToLogExport'][iT,:]=C_ToLogExport

    #--------------------------------------------------------------------------
    # Carbon transferred from mills to in-use products
    #--------------------------------------------------------------------------

    # Transfer mill fibre to single-family homes
    iP=meta['Core']['iPP']['SFH']
    vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT-1,:,iP] + \
        C_LumberMillToSFH + \
        C_PlywoodMillToSFH + \
        C_OSBMillToSFH + \
        C_MDFMillToSFH + \
        C_LogExportToSFH

    # Transfer mill fibre to multi-family homes
    iP=meta['Core']['iPP']['MFH']
    vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT-1,:,iP] + \
        C_LumberMillToMFH + \
        C_PlywoodMillToMFH + \
        C_OSBMillToMFH + \
        C_MDFMillToMFH + \
        C_LogExportToMFH

    # Transfer mill fibre to commercial
    iP=meta['Core']['iPP']['Comm']
    vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT-1,:,iP] + \
        C_LumberMillToCom + \
        C_PlywoodMillToCom + \
        C_OSBMillToCom + \
        C_MDFMillToCom + \
        C_LogExportToCom

    # Transfer mill fibre to furniture
    iP=meta['Core']['iPP']['Furn']
    vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT-1,:,iP] + \
        C_LumberMillToFurn + \
        C_PlywoodMillToFurn + \
        C_OSBMillToFurn + \
        C_MDFMillToFurn + \
        C_LogExportToFurn

    # Transfer mill fibre to shipping
    iP=meta['Core']['iPP']['Ship']
    vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT-1,:,iP] + \
        C_LumberMillToShip + \
        C_PlywoodMillToShip + \
        C_OSBMillToShip + \
        C_MDFMillToShip + \
        C_LogExportToShip

    # Transfer mill fibre to repairs
    iP=meta['Core']['iPP']['Repairs']
    vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT-1,:,iP] + \
        C_LumberMillToRepairs + \
        C_PlywoodMillToRepairs + \
        C_OSBMillToRepairs + \
        C_MDFMillToRepairs + \
        C_LogExportToRepairs

    # Transfer mill fibre to other
    iP=meta['Core']['iPP']['Other']
    vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT-1,:,iP] + \
        C_LumberMillToOther + \
        C_PlywoodMillToOther + \
        C_OSBMillToOther + \
        C_MDFMillToOther + \
        C_LogExportToOther

    # Transfer pulp mill fibre to paper
    iP=meta['Core']['iPP']['Paper']
    vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT-1,:,iP]+C_ToPaper

    # Transfer mill fibre to power facitity, domestic
    iP=meta['Core']['iPP']['PowerFacilityDom']
    vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT-1,:,iP]+C_ToPowerFacilityDom

    # Transfer mill fibre to power facitity, foreign
    iP=meta['Core']['iPP']['PowerFacilityExport']
    vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT-1,:,iP]+C_ToPowerFacilityExport

    # Transfer mill fibre to pellet Export
    iP=meta['Core']['iPP']['PelletExport']
    vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT-1,:,iP]+C_ToPelletExport

    # Transfer mill fibre to pellet domestic grid
    iP=meta['Core']['iPP']['PelletDomGrid']
    vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT-1,:,iP]+C_ToPelletDomGrid

    # Transfer mill fibre to pellet domestic RNG
    iP=meta['Core']['iPP']['PelletDomRNG']
    vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT-1,:,iP]+C_ToPelletDomRNG

    # Transfer domestic firewood to domestic firewood pool
    iP=meta['Core']['iPP']['FirewoodDom']
    vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT-1,:,iP]+C_ToFirewoodCollection

    # Transfer foreign firewood to foreign firewood pool
    iP=meta['Core']['iPP']['FirewoodExport']
    vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT-1,:,iP]+C_ToFirewoodExport

    # Transfer pulp mill carbon to pulp-mill effluent
    iP=meta['Core']['iPP']['EffluentPulp']
    vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT-1,:,iP]+C_ToPulpEffluent

    #--------------------------------------------------------------------------
    # Update dump and landfill reservoirs
    #--------------------------------------------------------------------------

    vo['C_Pro_Pools'][iT,:,meta['Core']['iPP']['DumpWood']]=vo['C_Pro_Pools'][iT-1,:,meta['Core']['iPP']['DumpWood']]
    vo['C_Pro_Pools'][iT,:,meta['Core']['iPP']['DumpPaper']]=vo['C_Pro_Pools'][iT-1,:,meta['Core']['iPP']['DumpPaper']]
    vo['C_Pro_Pools'][iT,:,meta['Core']['iPP']['LandfillWoodDegradable']]=vo['C_Pro_Pools'][iT-1,:,meta['Core']['iPP']['LandfillWoodDegradable']]
    vo['C_Pro_Pools'][iT,:,meta['Core']['iPP']['LandfillWoodNonDegradable']]=vo['C_Pro_Pools'][iT-1,:,meta['Core']['iPP']['LandfillWoodNonDegradable']]
    vo['C_Pro_Pools'][iT,:,meta['Core']['iPP']['LandfillPaperDegradable']]=vo['C_Pro_Pools'][iT-1,:,meta['Core']['iPP']['LandfillPaperDegradable']]
    vo['C_Pro_Pools'][iT,:,meta['Core']['iPP']['LandfillPaperNonDegradable']]=vo['C_Pro_Pools'][iT-1,:,meta['Core']['iPP']['LandfillPaperNonDegradable']]

    #--------------------------------------------------------------------------
    # Single-family homes --> dump and landfill
    #--------------------------------------------------------------------------

    # Turnover
    iP=meta['Core']['iPP']['SFH']
    C_retired=meta['Param']['BEV']['HWP']['SFH_tr']*vo['C_Pro_Pools'][iT-1,:,iP]

    # Remove carbon
    vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT,:,iP] - C_retired

    # Transfer carbon to dump wood
    iP=meta['Core']['iPP']['DumpWood']
    vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT,:,iP] + meta['Param']['BEV']['HWP']['SFHToDumpWood']*C_retired

    # Transfer carbon to landfill (degradable)
    iP=meta['Core']['iPP']['LandfillWoodDegradable']
    vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT,:,iP] + meta['Param']['BEV']['HWP']['SFHToLandfillWood']*meta['Param']['BEV']['HWP']['ToLandfillWoodDegradableFrac']*C_retired

    # Transfer carbon to landfill (non-degradable)
    iP=meta['Core']['iPP']['LandfillWoodNonDegradable']
    vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT,:,iP] + meta['Param']['BEV']['HWP']['SFHToLandfillWood']*(1-meta['Param']['BEV']['HWP']['ToLandfillWoodDegradableFrac'])*C_retired

    #--------------------------------------------------------------------------
    # Multi-family homes --> dump and landfill
    #--------------------------------------------------------------------------

    # Turnover
    iP=meta['Core']['iPP']['MFH']
    C_retired=meta['Param']['BEV']['HWP']['MFH_tr']*vo['C_Pro_Pools'][iT,:,iP]

    # Remove carbon
    vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT,:,iP] - C_retired

    # Transfer carbon to dump wood
    iP=meta['Core']['iPP']['DumpWood']
    vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT,:,iP] + meta['Param']['BEV']['HWP']['MFHToDumpWood']*C_retired

    # Transfer carbon to landfill (degradable)
    iP=meta['Core']['iPP']['LandfillWoodDegradable']
    vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT,:,iP] + meta['Param']['BEV']['HWP']['MFHToLandfillWood']*meta['Param']['BEV']['HWP']['ToLandfillWoodDegradableFrac']*C_retired

    # Transfer carbon to landfill (non-degradable)
    iP=meta['Core']['iPP']['LandfillWoodNonDegradable']
    vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT,:,iP] + meta['Param']['BEV']['HWP']['MFHToLandfillWood']*(1-meta['Param']['BEV']['HWP']['ToLandfillWoodDegradableFrac'])*C_retired

    #--------------------------------------------------------------------------
    # Commercial building --> dump and landfill
    #--------------------------------------------------------------------------

    # Turnover
    iP=meta['Core']['iPP']['Comm']
    C_retired=meta['Param']['BEV']['HWP']['Comm_tr']*vo['C_Pro_Pools'][iT-1,:,iP]

    # Remove carbon
    vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT,:,iP] - C_retired

    # Transfer carbon to dump wood
    iP=meta['Core']['iPP']['DumpWood']
    vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT,:,iP] + meta['Param']['BEV']['HWP']['CommToDumpWood']*C_retired

    # Transfer carbon to landfill (degradable)
    iP=meta['Core']['iPP']['LandfillWoodDegradable']
    vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT,:,iP] + meta['Param']['BEV']['HWP']['CommToLandfillWood']*meta['Param']['BEV']['HWP']['ToLandfillWoodDegradableFrac']*C_retired

    # Transfer carbon to landfill (non-degradable)
    iP=meta['Core']['iPP']['LandfillWoodNonDegradable']
    vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT,:,iP] + meta['Param']['BEV']['HWP']['CommToLandfillWood']*(1-meta['Param']['BEV']['HWP']['ToLandfillWoodDegradableFrac'])*C_retired

    #--------------------------------------------------------------------------
    # Furniture --> dump and landfill
    #--------------------------------------------------------------------------

    # Turnover
    iP=meta['Core']['iPP']['Furn']
    C_retired=meta['Param']['BEV']['HWP']['Furn_tr']*vo['C_Pro_Pools'][iT-1,:,iP]

    # Remove carbon
    vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT,:,iP] - C_retired

    # Transfer carbon to dump wood
    iP=meta['Core']['iPP']['DumpWood']
    vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT,:,iP] + meta['Param']['BEV']['HWP']['FurnToDumpWood']*C_retired

    # Transfer carbon to landfill (degradable)
    iP=meta['Core']['iPP']['LandfillWoodDegradable']
    vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT,:,iP] + meta['Param']['BEV']['HWP']['FurnToLandfillWood']*meta['Param']['BEV']['HWP']['ToLandfillWoodDegradableFrac']*C_retired

    # Transfer carbon to landfill (non-degradable)
    iP=meta['Core']['iPP']['LandfillWoodNonDegradable']
    vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT,:,iP] + meta['Param']['BEV']['HWP']['FurnToLandfillWood']*(1-meta['Param']['BEV']['HWP']['ToLandfillWoodDegradableFrac'])*C_retired

    #--------------------------------------------------------------------------
    # Shipping --> dump and landfill
    #--------------------------------------------------------------------------

    # Turnover
    iP=meta['Core']['iPP']['Ship']
    C_retired=meta['Param']['BEV']['HWP']['Ship_tr']*vo['C_Pro_Pools'][iT-1,:,iP]

    # Remove carbon
    vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT,:,iP] - C_retired

    # Transfer carbon to dump wood
    iP=meta['Core']['iPP']['DumpWood']
    vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT,:,iP] + meta['Param']['BEV']['HWP']['ShipToDumpWood']*C_retired

    # Transfer carbon to landfill (degradable)
    iP=meta['Core']['iPP']['LandfillWoodDegradable']
    vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT,:,iP] + meta['Param']['BEV']['HWP']['ShipToLandfillWood']*meta['Param']['BEV']['HWP']['ToLandfillWoodDegradableFrac']*C_retired

    # Transfer carbon to landfill (non-degradable)
    iP=meta['Core']['iPP']['LandfillWoodNonDegradable']
    vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT,:,iP] + meta['Param']['BEV']['HWP']['ShipToLandfillWood']*(1-meta['Param']['BEV']['HWP']['ToLandfillWoodDegradableFrac'])*C_retired

    #--------------------------------------------------------------------------
    # Repairs --> dump and landfill
    #--------------------------------------------------------------------------

    # Turnover
    iP=meta['Core']['iPP']['Repairs']
    C_retired=meta['Param']['BEV']['HWP']['Repairs_tr']*vo['C_Pro_Pools'][iT-1,:,iP]

    # Remove carbon
    vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT,:,iP] - C_retired

    # Transfer carbon to dump wood
    iP=meta['Core']['iPP']['DumpWood']
    vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT,:,iP] + meta['Param']['BEV']['HWP']['RepairsToDumpWood']*C_retired

    # Transfer carbon to landfill (degradble)
    iP=meta['Core']['iPP']['LandfillWoodDegradable']
    vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT,:,iP] + meta['Param']['BEV']['HWP']['RepairsToLandfillWood']*meta['Param']['BEV']['HWP']['ToLandfillWoodDegradableFrac']*C_retired

    # Transfer carbon to landfill (non-degradable)
    iP=meta['Core']['iPP']['LandfillWoodNonDegradable']
    vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT,:,iP] + meta['Param']['BEV']['HWP']['RepairsToLandfillWood']*(1-meta['Param']['BEV']['HWP']['ToLandfillWoodDegradableFrac'])*C_retired

    #--------------------------------------------------------------------------
    # Other --> dump and landfill
    #--------------------------------------------------------------------------

    # Turnover
    iP=meta['Core']['iPP']['Other']
    C_retired=meta['Param']['BEV']['HWP']['Other_tr']*vo['C_Pro_Pools'][iT-1,:,iP]

    # Remove carbon
    vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT,:,iP] - C_retired

    # Transfer carbon to dump wood
    iP=meta['Core']['iPP']['DumpWood']
    vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT,:,iP] + meta['Param']['BEV']['HWP']['OtherToDumpWood']*C_retired

    # Transfer carbon to landfill (degradble)
    iP=meta['Core']['iPP']['LandfillWoodDegradable']
    vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT,:,iP] + meta['Param']['BEV']['HWP']['OtherToLandfillWood']*meta['Param']['BEV']['HWP']['ToLandfillWoodDegradableFrac']*C_retired

    # Transfer carbon to landfill (non-degradable)
    iP=meta['Core']['iPP']['LandfillWoodNonDegradable']
    vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT,:,iP] + meta['Param']['BEV']['HWP']['OtherToLandfillWood']*(1-meta['Param']['BEV']['HWP']['ToLandfillWoodDegradableFrac'])*C_retired

    #--------------------------------------------------------------------------
    # Paper --> dump and landfill
    #--------------------------------------------------------------------------

    # Turnover (with adjustment for recycling)
    iP=meta['Core']['iPP']['Paper']
    C_retired=(1-meta['Param']['BEV']['HWP']['PaperRecycleRate'])*meta['Param']['BEV']['HWP']['Paper_tr']*vo['C_Pro_Pools'][iT,:,iP]

    # Remove carbon
    vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT,:,iP] - C_retired

    # Transfer to dump
    iP=meta['Core']['iPP']['DumpPaper']
    vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT,:,iP] + meta['Param']['BEV']['HWP']['PaperToDumpPaper']*C_retired

    # Transfer to landfill (degradable)
    iP=meta['Core']['iPP']['LandfillWoodDegradable']
    vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT,:,iP] + meta['Param']['BEV']['HWP']['PaperToLandfillPaper']*meta['Param']['BEV']['HWP']['ToLandfillPaperDegradableFrac']*C_retired

    # Transfer to landfill (non-degradable)
    iP=meta['Core']['iPP']['LandfillWoodNonDegradable']
    vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT,:,iP] + meta['Param']['BEV']['HWP']['PaperToLandfillPaper']*(1-meta['Param']['BEV']['HWP']['ToLandfillPaperDegradableFrac'])*C_retired

    #--------------------------------------------------------------------------
    # Emissions from combustion during domestic power generation
    #--------------------------------------------------------------------------

    # Turnover
    iP=meta['Core']['iPP']['PowerFacilityDom']
    C_emitted=meta['Param']['BEV']['HWP']['Energy_tr']*vo['C_Pro_Pools'][iT,:,iP]

    # Remove carbon
    vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT,:,iP] - C_emitted

    # Emissions
    E_CO2e_AsCO2=meta['Param']['BEV']['Biophysical']['Ratio_CO2_to_C']*meta['Param']['BEV']['HWP']['EnergyCombustionFracEmitCO2']*C_emitted
    E_CH4=meta['Param']['BEV']['Biophysical']['Ratio_CH4_to_C']*(1-meta['Param']['BEV']['HWP']['EnergyCombustionFracEmitCO2'])*C_emitted
    E_CO2e_AsCH4=meta['Param']['BE']['Biophysical']['GWP_CH4_AR5']*E_CH4
    vo['E_CO2e_ESC_Bioenergy'][iT,:]=vo['E_CO2e_ESC_Bioenergy'][iT,:] + (E_CO2e_AsCO2+E_CO2e_AsCH4)
    vo['E_CO2e_ESC_BioenergyPowerFacilityDom'][iT,:]=vo['E_CO2e_ESC_BioenergyPowerFacilityDom'][iT,:] + (E_CO2e_AsCO2+E_CO2e_AsCH4)

    vo['Atm_CO2_In'][iT,:]=vo['Atm_CO2_In'][iT,:]+E_CO2e_AsCO2
    vo['Atm_CH4_In'][iT,:]=vo['Atm_CH4_In'][iT,:]+E_CH4

    #--------------------------------------------------------------------------
    # Emissions from combustion during foreign power generation
    #--------------------------------------------------------------------------

    # Turnover
    iP=meta['Core']['iPP']['PowerFacilityExport']
    C_emitted=meta['Param']['BEV']['HWP']['Energy_tr']*vo['C_Pro_Pools'][iT,:,iP]

    # Remove carbon
    vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT,:,iP] - C_emitted

    # Emissions
    E_CO2e_AsCO2=meta['Param']['BEV']['Biophysical']['Ratio_CO2_to_C']*meta['Param']['BEV']['HWP']['EnergyCombustionFracEmitCO2']*C_emitted
    E_CH4=meta['Param']['BEV']['Biophysical']['Ratio_CH4_to_C']*(1-meta['Param']['BEV']['HWP']['EnergyCombustionFracEmitCO2'])*C_emitted
    E_CO2e_AsCH4=meta['Param']['BE']['Biophysical']['GWP_CH4_AR5']*E_CH4
    vo['E_CO2e_ESC_Bioenergy'][iT,:]=vo['E_CO2e_ESC_Bioenergy'][iT,:] + (E_CO2e_AsCO2+E_CO2e_AsCH4)
    vo['E_CO2e_ESC_BioenergyPowerFacilityExport'][iT,:]=vo['E_CO2e_ESC_BioenergyPowerFacilityExport'][iT,:] + (E_CO2e_AsCO2+E_CO2e_AsCH4)

    vo['Atm_CO2_In'][iT,:]=vo['Atm_CO2_In'][iT,:]+E_CO2e_AsCO2
    vo['Atm_CH4_In'][iT,:]=vo['Atm_CH4_In'][iT,:]+E_CH4

    #--------------------------------------------------------------------------
    # Emissions from combustion of pellet Export
    #--------------------------------------------------------------------------

    # Turnover
    iP=meta['Core']['iPP']['PelletExport']
    C_emitted=meta['Param']['BEV']['HWP']['Energy_tr']*vo['C_Pro_Pools'][iT,:,iP]

    # Remove carbon
    vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT,:,iP] - C_emitted

    # Emissions
    E_CO2e_AsCO2=meta['Param']['BEV']['Biophysical']['Ratio_CO2_to_C']*meta['Param']['BEV']['HWP']['EnergyCombustionFracEmitCO2']*C_emitted
    E_CH4=meta['Param']['BEV']['Biophysical']['Ratio_CH4_to_C']*(1-meta['Param']['BEV']['HWP']['EnergyCombustionFracEmitCO2'])*C_emitted
    E_CO2e_AsCH4=meta['Param']['BE']['Biophysical']['GWP_CH4_AR5']*E_CH4
    vo['E_CO2e_ESC_Bioenergy'][iT,:]=vo['E_CO2e_ESC_Bioenergy'][iT,:] + (E_CO2e_AsCO2+E_CO2e_AsCH4)
    vo['E_CO2e_ESC_BioenergyPelletExport'][iT,:]=vo['E_CO2e_ESC_BioenergyPelletExport'][iT,:] + (E_CO2e_AsCO2+E_CO2e_AsCH4)

    vo['Atm_CO2_In'][iT,:]=vo['Atm_CO2_In'][iT,:]+E_CO2e_AsCO2
    vo['Atm_CH4_In'][iT,:]=vo['Atm_CH4_In'][iT,:]+E_CH4

    #--------------------------------------------------------------------------
    # Emissions from combustion of pellet domestic grid
    #--------------------------------------------------------------------------

    # Turnover
    iP=meta['Core']['iPP']['PelletDomGrid']
    C_emitted=meta['Param']['BEV']['HWP']['Energy_tr']*vo['C_Pro_Pools'][iT,:,iP]

    # Remove carbon
    vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT,:,iP] - C_emitted

    # Emissions
    E_CO2e_AsCO2=meta['Param']['BEV']['Biophysical']['Ratio_CO2_to_C']*meta['Param']['BEV']['HWP']['EnergyCombustionFracEmitCO2']*C_emitted
    E_CH4=meta['Param']['BEV']['Biophysical']['Ratio_CH4_to_C']*(1-meta['Param']['BEV']['HWP']['EnergyCombustionFracEmitCO2'])*C_emitted
    E_CO2e_AsCH4=meta['Param']['BE']['Biophysical']['GWP_CH4_AR5']*E_CH4
    vo['E_CO2e_ESC_Bioenergy'][iT,:]=vo['E_CO2e_ESC_Bioenergy'][iT,:] + (E_CO2e_AsCO2+E_CO2e_AsCH4)
    vo['E_CO2e_ESC_BioenergyPelletDomGrid'][iT,:]=vo['E_CO2e_ESC_BioenergyPelletDomGrid'][iT,:] + (E_CO2e_AsCO2+E_CO2e_AsCH4)

    vo['Atm_CO2_In'][iT,:]=vo['Atm_CO2_In'][iT,:]+E_CO2e_AsCO2
    vo['Atm_CH4_In'][iT,:]=vo['Atm_CH4_In'][iT,:]+E_CH4

    #--------------------------------------------------------------------------
    # Emissions from combustion of pellet domestic RNG
    #--------------------------------------------------------------------------

    # Turnover
    iP=meta['Core']['iPP']['PelletDomRNG']
    C_emitted=meta['Param']['BEV']['HWP']['Energy_tr']*vo['C_Pro_Pools'][iT,:,iP]

    # Remove carbon
    vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT,:,iP] - C_emitted

    # Emissions
    E_CO2e_AsCO2=meta['Param']['BEV']['Biophysical']['Ratio_CO2_to_C']*meta['Param']['BEV']['HWP']['EnergyCombustionFracEmitCO2']*C_emitted
    E_CH4=meta['Param']['BEV']['Biophysical']['Ratio_CH4_to_C']*(1-meta['Param']['BEV']['HWP']['EnergyCombustionFracEmitCO2'])*C_emitted
    E_CO2e_AsCH4=meta['Param']['BE']['Biophysical']['GWP_CH4_AR5']*E_CH4
    vo['E_CO2e_ESC_Bioenergy'][iT,:]=vo['E_CO2e_ESC_Bioenergy'][iT,:] + (E_CO2e_AsCO2+E_CO2e_AsCH4)
    vo['E_CO2e_ESC_BioenergyPelletDomRNG'][iT,:]=vo['E_CO2e_ESC_BioenergyPelletDomRNG'][iT,:] + (E_CO2e_AsCO2+E_CO2e_AsCH4)

    vo['Atm_CO2_In'][iT,:]=vo['Atm_CO2_In'][iT,:]+E_CO2e_AsCO2
    vo['Atm_CH4_In'][iT,:]=vo['Atm_CH4_In'][iT,:]+E_CH4

    #--------------------------------------------------------------------------
    # Emissions from combustion of domestic firewood
    #--------------------------------------------------------------------------

    # Turnover
    iP=meta['Core']['iPP']['FirewoodDom']
    C_emitted=meta['Param']['BEV']['HWP']['Firewood_tr']*vo['C_Pro_Pools'][iT,:,iP]

    # Remove carbon
    vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT,:,iP] - C_emitted

    # Emissions
    E_CO2e_AsCO2=meta['Param']['BEV']['Biophysical']['Ratio_CO2_to_C']*meta['Param']['BEV']['HWP']['EnergyCombustionFracEmitCO2']*C_emitted
    E_CH4=meta['Param']['BEV']['Biophysical']['Ratio_CH4_to_C']*(1-meta['Param']['BEV']['HWP']['EnergyCombustionFracEmitCO2'])*C_emitted
    E_CO2e_AsCH4=meta['Param']['BE']['Biophysical']['GWP_CH4_AR5']*E_CH4
    vo['E_CO2e_ESC_Bioenergy'][iT,:]=vo['E_CO2e_ESC_Bioenergy'][iT,:] + (E_CO2e_AsCO2+E_CO2e_AsCH4)
    vo['E_CO2e_ESC_BioenergyFirewoodDom'][iT,:]=vo['E_CO2e_ESC_BioenergyFirewoodDom'][iT,:] + (E_CO2e_AsCO2+E_CO2e_AsCH4)

    vo['Atm_CO2_In'][iT,:]=vo['Atm_CO2_In'][iT,:]+E_CO2e_AsCO2
    vo['Atm_CH4_In'][iT,:]=vo['Atm_CH4_In'][iT,:]+E_CH4

    #--------------------------------------------------------------------------
    # Emissions from combustion of foreign firewood
    #--------------------------------------------------------------------------

    # Turnover
    iP=meta['Core']['iPP']['FirewoodExport']
    C_emitted=meta['Param']['BEV']['HWP']['Firewood_tr']*vo['C_Pro_Pools'][iT,:,iP]

    # Remove carbon
    vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT,:,iP] - C_emitted

    # Emissions
    E_CO2e_AsCO2=meta['Param']['BEV']['Biophysical']['Ratio_CO2_to_C']*meta['Param']['BEV']['HWP']['EnergyCombustionFracEmitCO2']*C_emitted
    E_CH4=meta['Param']['BEV']['Biophysical']['Ratio_CH4_to_C']*(1-meta['Param']['BEV']['HWP']['EnergyCombustionFracEmitCO2'])*C_emitted
    E_CO2e_AsCH4=meta['Param']['BE']['Biophysical']['GWP_CH4_AR5']*E_CH4
    vo['E_CO2e_ESC_Bioenergy'][iT,:]=vo['E_CO2e_ESC_Bioenergy'][iT,:] + (E_CO2e_AsCO2+E_CO2e_AsCH4)
    vo['E_CO2e_ESC_BioenergyFirewoodExport'][iT,:]=vo['E_CO2e_ESC_BioenergyFirewoodExport'][iT,:] + (E_CO2e_AsCO2+E_CO2e_AsCH4)

    vo['Atm_CO2_In'][iT,:]=vo['Atm_CO2_In'][iT,:]+E_CO2e_AsCO2
    vo['Atm_CH4_In'][iT,:]=vo['Atm_CH4_In'][iT,:]+E_CH4

    #--------------------------------------------------------------------------
    # Emissions from pulp effluent
    #--------------------------------------------------------------------------

    # Emissions from pulp effluent (CO2 from aerobic decomposition)
    iP=meta['Core']['iPP']['EffluentPulp']
    C_emitted=meta['Param']['BEV']['HWP']['EffluentPulp_tr']*vo['C_Pro_Pools'][iT,:,iP]

    # Remove emitted carbon
    vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT,:,iP]-C_emitted

    # Add emitted carbon to CO2 emission "pool"
    # Emissions
    E_CO2e_AsCO2=meta['Param']['BEV']['Biophysical']['Ratio_CO2_to_C']*meta['Param']['BEV']['HWP']['EnergyCombustionFracEmitCO2']*C_emitted
    vo['E_CO2e_LULUCF_HWP'][iT,:]=vo['E_CO2e_LULUCF_HWP'][iT,:] + E_CO2e_AsCO2

    vo['Atm_CO2_In'][iT,:]=vo['Atm_CO2_In'][iT,:]+E_CO2e_AsCO2
    #vo['Atm_CH4_In'][iT,:]=vo['Atm_CH4_In'][iT,:]+E_CH4

    #--------------------------------------------------------------------------
    # Decomposition of dump wood
    #--------------------------------------------------------------------------

    # Turnover
    iP=meta['Core']['iPP']['DumpWood']
    C_emitted=meta['Param']['BEV']['HWP']['DumpWood_tr']*vo['C_Pro_Pools'][iT,:,iP]

    # Removal
    vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT,:,iP]-C_emitted

    E_CO2e_AsCO2=meta['Param']['BEV']['Biophysical']['Ratio_CO2_to_C']*C_emitted

    # Add to emissions (CO2 emission from aerobic decomposition)
    vo['E_CO2e_LULUCF_HWP'][iT,:]=vo['E_CO2e_LULUCF_HWP'][iT,:] + E_CO2e_AsCO2

    vo['Atm_CO2_In'][iT,:]=vo['Atm_CO2_In'][iT,:]+E_CO2e_AsCO2
    #vo['Atm_CH4_In'][iT,:]=vo['Atm_CH4_In'][iT,:]+E_CH4

    #--------------------------------------------------------------------------
    # Decomposition of dump paper
    #--------------------------------------------------------------------------

    # Turnover
    iP=meta['Core']['iPP']['DumpPaper']
    C_emitted=meta['Param']['BEV']['HWP']['DumpPaper_tr']*vo['C_Pro_Pools'][iT,:,iP]

    # Removal
    vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT,:,iP]-C_emitted

    E_CO2e_AsCO2=meta['Param']['BEV']['Biophysical']['Ratio_CO2_to_C']*C_emitted

    # Add to emissions (CO2 emission from aerobic decomposition)
    vo['E_CO2e_LULUCF_HWP'][iT,:]=vo['E_CO2e_LULUCF_HWP'][iT,:] + E_CO2e_AsCO2

    vo['Atm_CO2_In'][iT,:]=vo['Atm_CO2_In'][iT,:]+E_CO2e_AsCO2
    #vo['Atm_CH4_In'][iT,:]=vo['Atm_CH4_In'][iT,:]+E_CH4

    #--------------------------------------------------------------------------
    # Decomposition of landfill degradable wood
    #--------------------------------------------------------------------------

    # Turnover
    iP=meta['Core']['iPP']['LandfillWoodDegradable']
    C_emitted=meta['Param']['BEV']['HWP']['LandfillWoodDegradable_tr']*vo['C_Pro_Pools'][iT,:,iP]

    # Removal
    vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT,:,iP]-C_emitted

    # Add to emissions
    E_CO2e_AsCO2=meta['Param']['BEV']['Biophysical']['Ratio_CO2_to_C']*meta['Param']['BEV']['HWP']['LandfillDegradableFracEmitCO2']*C_emitted

    vo['E_CO2e_LULUCF_HWP'][iT,:]=vo['E_CO2e_LULUCF_HWP'][iT,:]+E_CO2e_AsCO2

    # Adjustment for proportion of degradable landfills with gas collection systems,
    # efficiency of system, and methane oxided to CO2 from the landfill cover
    gcsp=meta['Param']['BEV']['HWP']['LandfillMethaneEmit_GasColSysProp']
    c1=1-gcsp
    c2=1-meta['Param']['BEV']['HWP']['LandfillMethaneEmit_GasColSysEffic']
    E_C_AsCH4=C_emitted*(c1-meta['Param']['BEV']['HWP']['LandfillMethaneOxidizedToCO2']*c1)+C_emitted*gcsp*(c2-meta['Param']['BEV']['HWP']['LandfillMethaneOxidizedToCO2']*c1)

    E_CH4=meta['Param']['BEV']['Biophysical']['Ratio_CH4_to_C']*E_C_AsCH4

    E_CO2e_AsCH4=meta['Param']['BEV']['Biophysical']['GWP_CH4_AR5']*E_CH4

    vo['E_CO2e_LULUCF_HWP'][iT,:]=vo['E_CO2e_LULUCF_HWP'][iT,:]+E_CO2e_AsCH4

    vo['Atm_CO2_In'][iT,:]=vo['Atm_CO2_In'][iT,:]+E_CO2e_AsCO2
    vo['Atm_CH4_In'][iT,:]=vo['Atm_CH4_In'][iT,:]+E_CH4

    #--------------------------------------------------------------------------
    # Decomposition of landfill degradable paper
    #--------------------------------------------------------------------------

    # Turnover
    iP=meta['Core']['iPP']['LandfillPaperDegradable']
    c_emitted=meta['Param']['BEV']['HWP']['LandfillWoodDegradable_tr']*vo['C_Pro_Pools'][iT,:,iP]

    # Removal
    vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT,:,iP]-c_emitted

    # Add to emissions
    E_CO2e_AsCO2=meta['Param']['BEV']['Biophysical']['Ratio_CO2_to_C']*meta['Param']['BEV']['HWP']['LandfillDegradableFracEmitCO2']*C_emitted

    vo['E_CO2e_LULUCF_HWP'][iT,:]=vo['E_CO2e_LULUCF_HWP'][iT,:]+E_CO2e_AsCO2

    # Adjustment for proportion of degradable landfills with gas collection systems,
    # efficiency of system, and methane oxided to CO2 from the landfill cover
    E_C_AsCH4=C_emitted*(c1-meta['Param']['BEV']['HWP']['LandfillMethaneOxidizedToCO2']*c1)+C_emitted*gcsp*(c2-meta['Param']['BEV']['HWP']['LandfillMethaneOxidizedToCO2']*c1)
    E_CO2e_AsCH4=meta['Param']['BEV']['Biophysical']['GWP_CH4_AR5']*E_C_AsCH4

    E_CH4=meta['Param']['BEV']['Biophysical']['Ratio_CH4_to_C']*E_C_AsCH4

    E_CO2e_AsCH4=meta['Param']['BEV']['Biophysical']['GWP_CH4_AR5']*E_CH4

    vo['E_CO2e_LULUCF_HWP'][iT,:]=vo['E_CO2e_LULUCF_HWP'][iT,:]+E_CO2e_AsCH4

    vo['Atm_CO2_In'][iT,:]=vo['Atm_CO2_In'][iT,:]+E_CO2e_AsCO2
    vo['Atm_CH4_In'][iT,:]=vo['Atm_CH4_In'][iT,:]+E_CH4

    return vo

#%% Grassland module

def Biomass_FromGrasses(meta,pNam,iScn,iBat,iT,vi,vo,iEP):

    # Aboveground biomass (foliage)
    vo['C_G_Gross'][iT,:,iEP['Foliage']]=vo['C_G_Gross'][iT,:,iEP['Foliage']]+1.5
    vo['C_LF'][iT,:,iEP['Foliage']]=vo['C_LF'][iT,:,iEP['Foliage']]+1.5

    # Belowground biomass
    if vo['C_Eco_Pools'][iT,:,iEP['RootCoarse']]>15:
        vo['C_G_Gross'][iT,:,iEP['RootCoarse']]=+1.5
    else:
        vo['C_G_Gross'][iT,:,iEP['RootCoarse']]=+1.75

    vo['C_LF'][iT,:,iEP['RootCoarse']]=+1.5

    # Net growth

    G_net=vo['C_G_Gross'][iT,:,iEP['RootCoarse']]-vo['C_LF'][iT,:,iEP['RootCoarse']]

    #vo['C_G_Net'][iT,:,iEP['RootCoarse']]=vo['C_G_Net'][iT,:,iEP['RootCoarse']]+G_net

    vo['C_Eco_Pools'][iT,:,iEP['RootCoarse']]=vo['C_Eco_Pools'][iT,:,iEP['RootCoarse']]+G_net

    return vo


