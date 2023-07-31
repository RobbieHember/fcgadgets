'''
Aspatial Statistical Models of Disturbance Events

'''

#%% Import modules

import numpy as np
import scipy.stats as stats
from fcgadgets.macgyver import util_general as gu
from fcgadgets.cbrunner import cbrun_util as cbu
import warnings
warnings.filterwarnings("ignore")

#%% Simulate wildfire occurrence and severity from Annual area of occurrence models

def SimulateWildfireFromAAO(meta,pNam,inv):

    # Ensure BGC zone has the right key name
    if 'ID_BECZ' in inv:
        inv['ID_BGCZ']=inv['ID_BECZ']

    # Import wildfire stats (by BGC zone)
    wfss=gu.ipickle(meta['Paths']['Model']['Taz Datasets'] + '\\Wildfire Stats and Scenarios\\Wildfire_Stats_Scenarios_By_BGCZ.pkl')
    tv_wfss=np.arange(-2000,2201,1)

    # Prepare mortality probability coefficients
    beta_mort_pi=np.cumsum([meta['Param']['BE']['Taz']['WF']['p_Unburned_pi'],meta['Param']['BE']['Taz']['WF']['p_Low_pi'],meta['Param']['BE']['Taz']['WF']['p_Medium_pi'],meta['Param']['BE']['Taz']['WF']['p_High_pi']])
    beta_mort_obs=np.cumsum([meta['Param']['BE']['Taz']['WF']['p_Unburned_obs'],meta['Param']['BE']['Taz']['WF']['p_Low_obs'],meta['Param']['BE']['Taz']['WF']['p_Medium_obs'],meta['Param']['BE']['Taz']['WF']['p_High_obs']])
    beta_mort_fut=np.cumsum([meta['Param']['BE']['Taz']['WF']['p_Unburned_fut'],meta['Param']['BE']['Taz']['WF']['p_Low_fut'],meta['Param']['BE']['Taz']['WF']['p_Medium_fut'],meta['Param']['BE']['Taz']['WF']['p_High_fut']])

    # Extract wildfire scenario IDs
    id_wf=np.zeros(meta[pNam]['Project']['N Scenario'],dtype=int)
    for iScn in range(meta[pNam]['Project']['N Scenario']):
        id_wf[iScn]=meta[pNam]['Scenario'][iScn]['Wildfire Scenario ID']

    # Unique wildfire scenario IDs
    u_id_wf=np.unique(id_wf)

    for iEns in range(meta[pNam]['Project']['N Ensemble']):

        #----------------------------------------------------------------------
        # Generate random numbers (the same among scenarios, different by ensemble)
        #----------------------------------------------------------------------

        rn_oc=np.random.random((meta[pNam]['Project']['N Time'],meta[pNam]['Project']['N Stand']))
        rn_sev=np.random.random((meta[pNam]['Project']['N Time'],meta[pNam]['Project']['N Stand']))

        #----------------------------------------------------------------------
        # Occurrence
        # If all project scenarios run the same wildfire scenario, only calculate
        # Pocc once outside of scenario loop to save time. If it varies, then
        # it needs to be done within the scenario loop
        #----------------------------------------------------------------------

        # All WF scenarios are the same, so calculate Pocc once outside of loop
        P_oc=[None]*u_id_wf.size
        for iWF in range(u_id_wf.size):

            P_oc[iWF]=np.zeros((meta[pNam]['Project']['N Time'],meta[pNam]['Project']['N Stand']))

            # Get unique BGC zone
            uZone=np.unique(inv['ID_BGCZ'])

            for iZone in range(uZone.size):

                namZone=cbu.lut_n2s(meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE'],uZone[iZone])[0]
                indZone=np.where(inv['ID_BGCZ']==uZone[iZone])[0]

                # Use the specified scenario
                Po_Det=wfss[namZone]['Po_Det_WF_Scn' + str(u_id_wf[iWF])]

                for iT in range(meta[pNam]['Year'].size):

                    # Adjust shape parameter to match specified annual probability of
                    # occurrence from the deterministic component
                    ind_scn=np.where(tv_wfss==meta[pNam]['Year'][iT])[0]
                    beta=wfss[namZone]['Beta_Pareto_Cal'].copy()
                    Scale=wfss[namZone]['Pareto_scale_to_match_Po_mu'][1]*Po_Det[ind_scn]+wfss[namZone]['Pareto_scale_to_match_Po_mu'][0]
                    beta[1]=-np.abs(Scale)
                    beta[2]=np.abs(Scale)

                    # Draw of annual area burned from Pareto distribution
                    N_t=1
                    P_oc[iWF][iT,indZone]=stats.pareto.rvs(beta[0],loc=beta[1],scale=beta[2],size=N_t)

        #----------------------------------------------------------------------
        # Populate for scenarios
        #----------------------------------------------------------------------

        for iScn in range(meta[pNam]['Project']['N Scenario']):

            # Index to unique wildfire scenario ID
            iU_wf=np.where(u_id_wf==id_wf[iScn])[0][0]

            # Initialize annual probability of occurrence (final with deterministic and
            # random components)
            wf_sim={}
            wf_sim['Occurrence']=np.zeros((meta[pNam]['Year'].size,meta[pNam]['Project']['N Stand']),dtype='int8')
            wf_sim['Mortality']=np.zeros((meta[pNam]['Year'].size,meta[pNam]['Project']['N Stand']),dtype='int8')

            # Populate occurrence
            iOc=np.where(rn_oc<P_oc[iU_wf])
            wf_sim['Occurrence'][iOc]=1

            # Exclude modern period
            if meta[pNam]['Scenario'][iScn]['Wildfire Status Pre-modern']=='Off':
                ind=np.where( (meta[pNam]['Year']<1921) )[0]
                wf_sim['Occurrence'][ind,:]=0

            if meta[pNam]['Scenario'][iScn]['Wildfire Status Modern']=='Off':
                ind=np.where( (meta[pNam]['Year']>=1921) & (meta[pNam]['Year']<meta[pNam]['Project']['Year Project']) )[0]
                wf_sim['Occurrence'][ind,:]=0

            if meta[pNam]['Scenario'][iScn]['Wildfire Status Future']=='Off':
                ind=np.where( (meta[pNam]['Year']>=meta[pNam]['Project']['Year Project']) )[0]
                wf_sim['Occurrence'][ind,:]=0

            #------------------------------------------------------------------
            # Mortality
            #------------------------------------------------------------------

            for iT in range(meta[pNam]['Year'].size):

                indOc=np.where( (wf_sim['Occurrence'][iT,:]>0) )[0]
                if indOc.size==0:
                    continue

                if (meta[pNam]['Year'][iT]<1920):
                    beta=beta_mort_pi
                elif  (meta[pNam]['Scenario'][iScn]['Wildfire Scenario ID']==1):
                    beta=beta_mort_pi
                elif (meta[pNam]['Year'][iT]>1920) & (meta[pNam]['Year'][iT]<=meta[pNam]['Project']['Year Project']) & (meta[pNam]['Scenario'][iScn]['Wildfire Scenario ID']!=1):
                    beta=beta_mort_obs
                elif (meta[pNam]['Year'][iT]>meta[pNam]['Project']['Year Project']) & (meta[pNam]['Scenario'][iScn]['Wildfire Scenario ID']==2):
                    beta=beta_mort_obs
                else:
                    beta=beta_mort_fut

                # Unburned
                ind=np.where(rn_sev[iT,indOc]<beta[0])[0]
                if ind.size>0:
                    wf_sim['Mortality'][iT,indOc[ind]]=5

                # Low severity
                ind=np.where( (rn_sev[iT,indOc]>=beta[0]) & (rn_sev[iT,indOc]<beta[1]) )[0]
                if ind.size>0:
                    wf_sim['Mortality'][iT,indOc[ind]]=50

                # Medium severity
                ind=np.where( (rn_sev[iT,indOc]>=beta[1]) & (rn_sev[iT,indOc]<beta[2]) )[0]
                if ind.size>0:
                    wf_sim['Mortality'][iT,indOc[ind]]=90

                # High severity
                ind=np.where( (rn_sev[iT,indOc]<beta[2]) )[0]
                if ind.size>0:
                    wf_sim['Mortality'][iT,indOc[ind]]=100

            # Compress to sparse vectors
            wf_sim_sparse={}
            wf_sim_sparse['idx']=np.where(wf_sim['Occurrence']>0)
            wf_sim_sparse['Occurrence']=wf_sim['Occurrence'][wf_sim_sparse['idx']]
            wf_sim_sparse['Mortality']=wf_sim['Mortality'][wf_sim_sparse['idx']]

            # Save
            fout=meta['Paths'][pNam]['Data'] + '\\Inputs\\Ensembles\\wf_sim_Scn' + cbu.FixFileNum(iScn) + '_Ens' + cbu.FixFileNum(iEns) + '.pkl'
            gu.opickle(fout,wf_sim_sparse)

#            #------------------------------------------------------------------
#            # If the project is scheduling future NOSE, schedule here for both
#            # baseline and actual scenario
#            #------------------------------------------------------------------
#
#            if meta[pNam]['Scenario'][iScn]['NOSE Future Status']=='On':
#
#                Year=np.tile( np.reshape(meta[pNam]['Year'],(meta[pNam]['Year'].size,-1)) ,(1,meta[pNam]['Project']['N Stand']))
#
#                rn=np.random.random(wf_sim['Occurrence'].shape)
#
#                ind=np.where( (Year>=meta[pNam]['Project']['Year Project']) & (wf_sim['Occurrence']==1) & (wf_sim['Mortality'][iT,:]>=80) & (rn<=meta[pNam]['Scenario'][iScn]['NOSE Future Prob']) )
#
#                meta[pNam]['Project']['NOSE Future']={}
#                meta[pNam]['Project']['NOSE Future']['Year']=Year[ind]
#                meta[pNam]['Project']['NOSE Future']['Stand Index']=ind[1]

    return

#%% Generate wildfire occurrence and severity from Annual area of occurrence models

def SimulateWildfireFromAAO_StandsActAsEnsembles(meta,pNam,inv,iScn):

    # Import wildfire stats (by BGC zone)
    wfss=gu.ipickle(meta['Paths']['Model']['Taz Datasets'] + '\\Wildfire Stats and Scenarios\\Wildfire_Stats_Scenarios_By_BGCZ.pkl')
    tv_wfss=np.arange(-2000,2201,1)

    # Prepare mortality probability coefficients
    beta_pi=np.cumsum([meta['Param']['BE']['Taz']['WF']['p_Unburned_pi'],meta['Param']['BE']['Taz']['WF']['p_Low_pi'],meta['Param']['BE']['Taz']['WF']['p_Medium_pi'],meta['Param']['BE']['Taz']['WF']['p_High_pi']])
    beta_obs=np.cumsum([meta['Param']['BE']['Taz']['WF']['p_Unburned_obs'],meta['Param']['BE']['Taz']['WF']['p_Low_obs'],meta['Param']['BE']['Taz']['WF']['p_Medium_obs'],meta['Param']['BE']['Taz']['WF']['p_High_obs']])

    # Get deterministic component of Po (specific to BGC zone)
    namZone=cbu.lut_n2s(meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE'],inv['ID_BECZ'][0,0])[0]
    P_oc_Det=wfss[namZone]['Po_Det_WF_Scn' + str(int(meta['Param']['BE']['Taz']['WF']['Scenario ID']))]

    # Initialize annual probability of occurrence (final with deterministic and
    # random components)
    wf_sim={}
    wf_sim['Occurrence']=np.zeros((meta[pNam]['Year'].size,meta[pNam]['Project']['N Stand']),dtype='int8')
    wf_sim['Mortality']=np.zeros((meta[pNam]['Year'].size,meta[pNam]['Project']['N Stand']),dtype='int8')

    # Occurrence (by BGC zone) (the same among scenarios, different by ensemble)
    P_oc=np.zeros((meta[pNam]['Year'].size,meta[pNam]['Project']['N Stand']))
    for iT in range(meta[pNam]['Year'].size):

        if (meta[pNam]['Scenario'][iScn]['Wildfire Status Pre-modern']=='Off') & (meta[pNam]['Year'][iT]<=1920):
            continue

        if (meta[pNam]['Scenario'][iScn]['Wildfire Status Modern']=='Off') & (meta[pNam]['Year'][iT]>=1920) & (meta[pNam]['Year'][iT]>meta[pNam]['Project']['Year Project']):
            continue

        if (meta[pNam]['Scenario'][iScn]['Wildfire Status Future']=='Off') & (meta[pNam]['Year'][iT]>=meta[pNam]['Project']['Year Project']):
            continue

        # Adjust shape parameter to match specified annual probability of
        # occurrence from the deterministic component
        ind_scn=np.where(tv_wfss==meta[pNam]['Year'][iT])[0]
        beta=wfss[namZone]['Beta_Pareto_Cal'].copy()
        Scale=wfss[namZone]['Pareto_scale_to_match_Po_mu'][1]*P_oc_Det[ind_scn]+wfss[namZone]['Pareto_scale_to_match_Po_mu'][0]
        beta[1]=-Scale
        beta[2]=Scale

        # Draw of annual area burned from Pareto distribution
        P_oc[iT,:]=stats.pareto.rvs(beta[0],loc=beta[1],scale=beta[2],size=meta[pNam]['Project']['N Stand'])

    for iS in range(meta[pNam]['Project']['N Stand']):

        # Generate random numbers (the same among scenarios, different by ensemble)
        rn_oc=np.random.random(meta[pNam]['Project']['N Time'])
        rn_sev=np.random.random(meta[pNam]['Project']['N Time'])

        for iT in range(meta[pNam]['Year'].size):

            if rn_oc[iT]<P_oc[iT,iS]:

                # Populate occurrence
                wf_sim['Occurrence'][iT,iS]=1

                # Severity
                if meta[pNam]['Year'][iT]<1920:
                    beta=beta_pi
                else:
                    beta=beta_obs

                # Unburned
                if (rn_sev[iT]<beta[0]):
                    wf_sim['Mortality'][iT,iS]=5

                # Low severity
                if (rn_sev[iT]>=beta[0]) & (rn_sev[iT]<beta[1]):
                    wf_sim['Mortality'][iT,iS]=50

                # Medium severity
                if (rn_sev[iT]>=beta[1]) & (rn_sev[iT]<beta[2]):
                    wf_sim['Mortality'][iT,iS]=90

                # High severity
                if (rn_sev[iT]<beta[2]):
                    wf_sim['Mortality'][iT,iS]=100

    return wf_sim

#%% Generate Mountain Pine Beetle impacts from annual area of occurrence model

def SimulateIBMFromAAO(meta,pNam,inv):

    # Import IBM stats
    ibmss=gu.ipickle(meta['Paths']['Model']['Taz Datasets'] + '\\Beetle Stats and Scenarios\\IBM_Stats_Scenarios_By_BGCZ.pkl')

    # Prepare mortality probability coefficients
    beta_obs=np.cumsum([meta['Param']['BE']['Taz']['IBM']['p_Trace_obs'],
                       meta['Param']['BE']['Taz']['IBM']['p_Low_obs'],
                       meta['Param']['BE']['Taz']['IBM']['p_Medium_obs'],
                       meta['Param']['BE']['Taz']['IBM']['p_Severe_obs'],
                       meta['Param']['BE']['Taz']['IBM']['p_VerySevere_obs']])

    for iEns in range(meta[pNam]['Project']['N Ensemble']):

        #----------------------------------------------------------------------
        # Generate random numbers (the same among scenarios, different by ensemble)
        #----------------------------------------------------------------------

        rn_oc=np.random.random((meta[pNam]['Project']['N Time'],meta[pNam]['Project']['N Stand']))
        rn_sev=np.random.random((meta[pNam]['Project']['N Time'],meta[pNam]['Project']['N Stand']))

        #----------------------------------------------------------------------
        # Generate probability of occurrence (the same among scenarios, different by ensemble)
        #----------------------------------------------------------------------

        P_oc=np.zeros((meta[pNam]['Project']['N Time'],meta[pNam]['Project']['N Stand']))

        uZone=np.unique(inv['ID_BGCZ'])

        for iZone in range(uZone.size):

            indZone=np.where(inv['ID_BGCZ']==uZone[iZone])[0]

            namZone=cbu.lut_n2s(meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE'],uZone[iZone])[0]

            # Alternative model
            beta=ibmss[namZone]['Beta_Pareto_Alt'].copy()

            for iT in range(meta[pNam]['Year'].size):
                N_t=1
                P_oc[iT,indZone]=stats.pareto.rvs(beta[0],loc=beta[1],scale=beta[2],size=N_t)

        # # *** Add historical Mountain Pine Beetle outbreaks ***
        # yrL=[1835,1865,1875,1905,1915,1945,1955]
        # for iY in yrL:
        #     iT=np.where(meta[pNam]['Year']==iY)[0]
        #     P_oc[iT,:]=0.1

        #--------------------------------------------------------------------------
        # Loop through scenarios
        #--------------------------------------------------------------------------

        for iScn in range(meta[pNam]['Project']['N Scenario']):

            # Initialize
            ibm_sim={}
            ibm_sim['Occurrence']=np.zeros((meta[pNam]['Year'].size,meta[pNam]['Project']['N Stand']),dtype='int16')
            ibm_sim['Mortality']=np.zeros((meta[pNam]['Year'].size,meta[pNam]['Project']['N Stand']),dtype='int16')

            # Populate occurrence
            iOc=np.where(rn_oc<P_oc)
            ibm_sim['Occurrence'][iOc]=1

            # Exclude modern period
            if meta[pNam]['Scenario'][iScn]['MPB Status Pre-modern']=='Off':
                ind=np.where( (meta[pNam]['Year']<1951) )[0]
                ibm_sim['Occurrence'][ind,:]=0

            if meta[pNam]['Scenario'][iScn]['MPB Status Modern']=='Off':
                ind=np.where( (meta[pNam]['Year']>=1951) & (meta[pNam]['Year']<meta[pNam]['Project']['Year Project']) )[0]
                ibm_sim['Occurrence'][ind,:]=0

            if meta[pNam]['Scenario'][iScn]['MPB Status Future']=='Off':
                ind=np.where( (meta[pNam]['Year']>=meta[pNam]['Project']['Year Project']) )[0]
                ibm_sim['Occurrence'][ind,:]=0

            # Severity / mortality

            # Index to occurrence
            iOcc=np.where( (ibm_sim['Occurrence']>0) )

            # Add pre-inventory severity
            Mort=np.zeros((meta[pNam]['Year'].size,meta[pNam]['Project']['N Stand']),dtype='int16')
            Mort[iOcc[0],iOcc[1]]=GetMortalityFromIBMSeverity(rn_sev[iOcc],beta_obs)
            it=np.where(meta[pNam]['Year']<1951)[0]
            ibm_sim['Mortality'][it,:]=Mort[it,:]

            # Add post-inventory severity
            Mort=np.zeros((meta[pNam]['Year'].size,meta[pNam]['Project']['N Stand']),dtype='int16')
            Mort[iOcc[0],iOcc[1]]=GetMortalityFromIBMSeverity(rn_sev[iOcc],beta_obs)
            it=np.where(meta[pNam]['Year']>meta[pNam]['Project']['Year Project'])[0]
            ibm_sim['Mortality'][it,:]=Mort[it,:]

            # Compress to sparse vectors
            ibm_sim_sparse={}
            ibm_sim_sparse['idx']=np.where(ibm_sim['Occurrence']>0)
            ibm_sim_sparse['Occurrence']=ibm_sim['Occurrence'][ibm_sim_sparse['idx']]
            ibm_sim_sparse['Mortality']=ibm_sim['Mortality'][ibm_sim_sparse['idx']]

            # Save
            fout=meta['Paths'][pNam]['Data'] + '\\Inputs\\Ensembles\\ibm_sim_Scn' + cbu.FixFileNum(iScn) + '_Ens' + cbu.FixFileNum(iEns) + '.pkl'
            gu.opickle(fout,ibm_sim_sparse)

#%% Mortality from burn severity rating
# Mortality numbers come from DisturbanceBySeverityClass spreadsheet

def GetMortalityFromIBMSeverity(rn,beta):

    y=np.zeros(rn.size)

    for i in range(rn.size):

        if rn[i]<beta[0]:
            # Trace
            y[i]=1
        elif (rn[i]>=beta[0]) & (rn[i]<beta[1]):
            # Low
            y[i]=5
        elif (rn[i]>=beta[1]) & (rn[i]<beta[2]):
            # Medium
            y[i]=25
        elif (rn[i]>=beta[2]) & (rn[i]<beta[3]):
            # Severe
            y[i]=50
        elif (rn[i]>=beta[3]):
            # Severe
            y[i]=75

    return y

#%% Simulate probability of stand breakup based on age

def PredictStandBreakup_OnTheFly(meta,pNam,vi,iT,iEns,Age):

    # Plot function:
    flg=0
    if flg==1:
        beta=[-0.025,250]
        Age=np.arange(1,500)
        Po=1/(1+np.exp(beta[0]*(Age-beta[1])))
        fig,ax=plt.subplots(1,figsize=gu.cm2inch(7.8,7))
        ax.plot(Age,Po,'k-',linewidth=0.75,label='Default model')
        ax.set(position=[0.11,0.11,0.88,0.88],xlim=[0,500],xticks=np.arange(0,550,50),xlabel='Age, years',ylabel='Annual probability of breakup')
        ax.legend(loc='upper left',bbox_to_anchor=(0.06,0.92),frameon=False,facecolor='w')
        ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both')

    Po=1/(1+np.exp(meta['Param']['BE']['On The Fly']['Pa Breakup Saturation']*(Age-meta['Param']['BE']['On The Fly']['Pa Breakup Inflection'])))
    rn=meta[pNam]['Project']['On the Fly']['Random Numbers']['Breakup'][iT,:]

    indAffected=np.where(rn<Po)[0]
    if indAffected.size>0:
        for iA in indAffected:
            iAvailable=np.where(vi['EC']['ID Event Type'][iT,iA,:]==0)[0]
            if iAvailable.size>0:
                iE=iAvailable[0]+0
                if rn[iA]<0.6*Po[iA]:
                    Type=meta['LUT']['Event']['Disease Root']
                else:
                    Type=meta['LUT']['Event']['Wind']
                vi['EC']['ID Event Type'][iT,iA,iE]=Type
                vi['EC']['Mortality Factor'][iT,iA,iE]=1.0
                vi['EC']['ID Growth Curve'][iT,iA,iE]=1
            else:
                print('No space left in event chronology for on-the-fly breakup!')

    return vi

#%% Simulate probability of harvesting on the fly

def PredictHarvesting_OnTheFly(meta,pNam,vi,iT,iScn,iEns,V_Merch,Period):

    # Indicator of THLB (THLB=1, Non-THLB=0)
    flag_thlb=vi['Inv']['THLB'][iT,:]

    # Indicator of energy production (don't harvest on the fly, it is pre-scheduled)
    if meta[pNam]['Project']['Land Surface Class Dependent']!='No':
        iT_lsc=np.where(vi['Inv']['LSC']['tv']==meta[pNam]['Year'][iT])[0]
        if iT_lsc.size>0:
            flag_ep=1*( (vi['Inv']['LSC']['Use'][iT_lsc,:]==meta['LUT']['LSC']['Use']['Fuel Break']) | (vi['Inv']['LSC']['Use'][iT_lsc,:]==meta['LUT']['LSC']['Use']['Energy Production']) )
        else:
            flag_ep=np.ones(flag_thlb.size,dtype=int)
    else:
        flag_ep=np.ones(flag_thlb.size,dtype=int)

    if Period=='Historical':

        # Historical period

        bH=[0.00085,4.85,0.32,1968]
        f1=bH[0]*np.maximum(0,(meta[pNam]['Year'][iT]-1800)/100)**bH[1]
        f2=(1/(1+np.exp(bH[2]*(meta[pNam]['Year'][iT]-bH[3]))))
        Pa_H_Sat=f1*f2

        # Plot
        flg=0
        if flg==1:
            t=np.arange(1700,2001,1)
            f1=bH[0]*np.maximum(0,(t-1800)/100)**bH[1]
            f2=(1/(1+np.exp(bH[2]*(t-bH[3]))))
            Pa_H_Sat=f1*f2
            plt.plot(t,Pa_H_Sat,'g--',lw=1.5)

    else:

        # Future period
        if 'Pa Harvest Sat' in meta[pNam]['Scenario'][iScn]:
            Pa_H_Sat=meta[pNam]['Scenario'][iScn]['Pa Harvest Sat']*vi['Inv']['Harvest Index']
        else:
            Pa_H_Sat=meta['Param']['BE']['On The Fly']['Pa_Harvest_Sat']*vi['Inv']['Harvest Index']

    # Inflection point
    if 'Pa Harvest Inf' in meta[pNam]['Scenario'][iScn]:

        # Default has been overriden with project-specific value
        Pa_H_Inf=meta[pNam]['Scenario'][iScn]['Pa Harvest Inf']

    else:

        # Use default
        Pa_H_Inf=meta['Param']['BE']['On The Fly']['Pa_Harvest_Inflection']

    # Shape parameter
    Pa_H_Shape=meta['Param']['BE']['On The Fly']['Pa_Harvest_Shape']

    # Plot function:
    flg=0
    if flg==1:

        beta=[0.03,-0.04,400]
        V_Merch=np.arange(1,1200)
        Po=beta[0]*(1/(1+np.exp(beta[1]*(V_Merch-beta[2]))))

        plt.close('all')
        fig,ax=plt.subplots(1,figsize=gu.cm2inch(7.8,7))
        ax.plot(V_Merch,Po*100,'k-',linewidth=0.75,label='Harvest on-the-fly model 1')
        ax.set(position=[0.1,0.12,0.87,0.86],xlim=[0,800],xticks=np.arange(0,1300,100),xlabel='Merchantable volume (m$^3$ ha$^-$$^1$)', \
               ylim=[0,5],ylabel='Annual probability of harvest (%)')
        ax.legend(loc='upper left',bbox_to_anchor=(0.06,0.92),frameon=False,facecolor='w')
        ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both')
        gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Harvest\taz_ann_prob_harvest','png',500)

    # Annual probability of occurrence
    Po=Pa_H_Sat*(1/(1+np.exp(Pa_H_Shape*(V_Merch-Pa_H_Inf))))

    # Random number
    #rn=np.random.random(V_Merch.size)
    rn=meta[pNam]['Project']['On the Fly']['Random Numbers']['Harvest'][iT,:]

    # Occurrence
    Oc=flag_ep*flag_thlb*np.floor(np.minimum(1,Po/rn))

    # Index to occurrence
    indS=np.where(Oc==1)[0]

    if indS.size>0:
        for i in range(indS.size):
            iAvailable=np.where(vi['EC']['ID Event Type'][iT,indS[i],:]==0)[0]
            if iAvailable.size>0:
                iE=iAvailable[0]+0
                vi['EC']['ID Event Type'][iT,indS[i],iE]=meta['LUT']['Event']['Harvest']
                vi['EC']['Mortality Factor'][iT,indS[i],iE]=1.0
                vi['EC']['ID Growth Curve'][iT,indS[i],iE]=1
                iE=iAvailable[0]+1 # changing this to zero will cause the harvest to be overwritten
                vi['EC']['ID Event Type'][iT,indS[i],iE]=meta['LUT']['Event']['Slashpile Burn']
                vi['EC']['Mortality Factor'][iT,indS[i],iE]=1.0
                vi['EC']['ID Growth Curve'][iT,indS[i],iE]=2
                iE=iAvailable[0]+2
                vi['EC']['ID Event Type'][iT,indS[i],iE]=meta['LUT']['Event']['Planting']
                vi['EC']['Mortality Factor'][iT,indS[i],iE]=1.0
                vi['EC']['ID Growth Curve'][iT,indS[i],iE]=2

    return vi