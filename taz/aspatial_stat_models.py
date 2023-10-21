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

#%%
def PredictWildfire_OnTheFly(meta,pNam,vi,iT,iScn,iEns):
   
    sH='H' + str(meta[pNam]['Scenario'][iScn]['Wildfire Scn Pre-obs'])
    sF='F' + str(meta[pNam]['Scenario'][iScn]['Wildfire Scn Future'])
    
    if (meta[pNam]['Scenario'][iScn]['Wildfire Scn Pre-obs']==-9999) & (meta[pNam]['Scenario'][iScn]['Wildfire Scn Future']!=-9999):
        # The project is trying to include the future but it crashes because no history was set. 
        # Just hardwire the historical scenario - it will not actually be used.
        sH='H2'
    
    iT_wf=np.where(meta['Modules']['Taz']['wfss']['tv']==meta[pNam]['Year'][iT])[0]
    #Wildfire Scn Pre-obs
    Po=np.zeros(meta[pNam]['Project']['indBat'].size)
    for zone in vi['lsat']['bgcz idx'].keys():        
        Po_Det=meta['Modules']['Taz']['wfss']['By BGC Zone'][zone]['Po Det Scenarios'][sH][sF][iT_wf]        
        beta=meta['Modules']['Taz']['wfss']['By BGC Zone'][zone]['Beta_Pareto_Cal'].copy()
        b0=meta['Modules']['Taz']['wfss']['By BGC Zone'][zone]['Pareto_scale_to_match_Po_mu'][0]
        b1=meta['Modules']['Taz']['wfss']['By BGC Zone'][zone]['Pareto_scale_to_match_Po_mu'][1]
        Scale=b1*Po_Det+b0
        beta[1]=-np.abs(Scale)
        beta[2]=np.abs(Scale)
        N_t=1
        Po[vi['lsat']['bgcz idx'][zone]]=stats.pareto.rvs(beta[0],loc=beta[1],scale=beta[2],size=N_t)
    rn=np.random.random(meta[pNam]['Project']['indBat'].size)
    
    # Index to occurrence
    indS=np.where(rn<Po)[0]
    if indS.size>0:
        for i in range(indS.size):
            iAvailable=np.where(vi['EC']['ID Event Type'][iT,indS[i],:]==0)[0]
            if iAvailable.size>0:
                iE=iAvailable[0]+0
                vi['EC']['ID Event Type'][iT,indS[i],iE]=meta['LUT']['Event']['Wildfire']
                vi['EC']['Mortality Factor'][iT,indS[i],iE]=0.59
                vi['EC']['ID Growth Curve'][iT,indS[i],iE]=1    
    return vi

#%%
def PredictBeetles_OnTheFly(meta,pNam,vi,iT,iScn,iEns):   
    Po=meta['Modules']['Taz']['IBM Pareto']['Po'][iT,iEns]
    rn=np.random.random(meta[pNam]['Project']['indBat'].size)
    indS=np.where(rn<Po)[0]
    if indS.size>0:
        for i in range(indS.size):
            iAvailable=np.where(vi['EC']['ID Event Type'][iT,indS[i],:]==0)[0]
            if iAvailable.size>0:
                iE=iAvailable[0]+0
                vi['EC']['ID Event Type'][iT,indS[i],iE]=meta['LUT']['Event']['Mountain Pine Beetle']
                vi['EC']['Mortality Factor'][iT,indS[i],iE]=meta['Modules']['Taz']['IBM Pareto']['Mortality Fraction']
                vi['EC']['ID Growth Curve'][iT,indS[i],iE]=1    
    return vi

#%% Simulate probability of wind disturbance as a function of age
def PredictWind_OnTheFly(meta,pNam,vi,iT,iEns,Age):    
    flg=0
    if flg==1:
        beta=[-0.04,400]
        Age=np.arange(1,500)
        Po=1/(1+np.exp(beta[0]*(Age-beta[1])))
        fig,ax=plt.subplots(1,figsize=gu.cm2inch(7.8,7))
        ax.plot(Age,Po,'k-',linewidth=0.75,label='Default model')
        ax.set(position=[0.11,0.11,0.88,0.88],xlim=[0,500],xticks=np.arange(0,550,50),xlabel='Age, years',ylabel='Annual probability of breakup')
        ax.legend(loc='upper left',bbox_to_anchor=(0.06,0.92),frameon=False,facecolor='w')
        ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both')

    Po=meta['Param']['BE']['On The Fly']['Wind Saturation']*(1/(1+np.exp(meta['Param']['BE']['On The Fly']['Wind Shape']*(Age-meta['Param']['BE']['On The Fly']['Wind Inflection']))))
    rn=np.random.random(Age.size)
    #rn=meta[pNam]['Project']['On the Fly']['Random Numbers']['Breakup'][iT,:]    

    # Control disturbances before harvest
    if meta[pNam]['Project']['Scenario Source']!='Spreadsheet':
        yr=meta[pNam]['Year'][iT]
        yrH=vi['lsat']['Year Harvest First']    
        dY=yr-yrH
        ind=np.where( (dY>-300) & (dY<=0) )[0]
        Po[ind]=0

    indAffected=np.where(rn<Po)[0]
    if indAffected.size>0:
        for iA in indAffected:
            iAvailable=np.where(vi['EC']['ID Event Type'][iT,iA,:]==0)[0]
            if iAvailable.size>0:
                iE=iAvailable[0]+0
                vi['EC']['ID Event Type'][iT,iA,iE]=meta['LUT']['Event']['Wind']
                vi['EC']['Mortality Factor'][iT,iA,iE]=1.0
                vi['EC']['ID Growth Curve'][iT,iA,iE]=1
            else:
                print('No space left in event chronology for on-the-fly event!')

    return vi

#%% Simulate probability of wind disturbance as a function of age
def PredictDisease_OnTheFly(meta,pNam,vi,iT,iEns,Age):    
    flg=0
    if flg==1:
        beta=[-0.04,400]
        Age=np.arange(1,500)
        Po=1/(1+np.exp(beta[0]*(Age-beta[1])))
        fig,ax=plt.subplots(1,figsize=gu.cm2inch(7.8,7))
        ax.plot(Age,Po,'k-',linewidth=0.75,label='Default model')
        ax.set(position=[0.11,0.11,0.88,0.88],xlim=[0,500],xticks=np.arange(0,550,50),xlabel='Age, years',ylabel='Annual probability of disease')
        ax.legend(loc='upper left',bbox_to_anchor=(0.06,0.92),frameon=False,facecolor='w')
        ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both')

    Po=meta['Param']['BE']['On The Fly']['Disease Saturation']*(1/(1+np.exp(meta['Param']['BE']['On The Fly']['Disease Shape']*(Age-meta['Param']['BE']['On The Fly']['Disease Inflection']))))
    rn=np.random.random(Age.size)

    # Control disturbances before harvest
    if meta[pNam]['Project']['Scenario Source']!='Spreadsheet':
        yr=meta[pNam]['Year'][iT]
        yrH=vi['lsat']['Year Harvest First']    
        dY=yr-yrH
        ind=np.where( (dY>-300) & (dY<=0) )[0]
        Po[ind]=0

    indAffected=np.where(rn<Po)[0]
    if indAffected.size>0:
        for iA in indAffected:
            iAvailable=np.where(vi['EC']['ID Event Type'][iT,iA,:]==0)[0]
            if iAvailable.size>0:
                iE=iAvailable[0]+0
                vi['EC']['ID Event Type'][iT,iA,iE]=meta['LUT']['Event']['Disease Root']
                vi['EC']['Mortality Factor'][iT,iA,iE]=1.0
                vi['EC']['ID Growth Curve'][iT,iA,iE]=1
            else:
                print('No space left in event chronology for on-the-fly event!')

    return vi

#%%
def PredictFrost_OnTheFly(meta,pNam,vi,iT,iScn,iEns):   
    Po=meta['Modules']['Taz']['Frost Pareto']['Po'][iT,iEns]
    rn=np.random.random(meta[pNam]['Project']['indBat'].size)
    indS=np.where(rn<Po)[0]
    if indS.size>0:
        for i in range(indS.size):
            iAvailable=np.where(vi['EC']['ID Event Type'][iT,indS[i],:]==0)[0]
            if iAvailable.size>0:
                iE=iAvailable[0]+0
                vi['EC']['ID Event Type'][iT,indS[i],iE]=meta['LUT']['Event']['Frost Snow Ice Hail']
                vi['EC']['Mortality Factor'][iT,indS[i],iE]=meta['Modules']['Taz']['Frost Pareto']['Mortality Fraction']
                vi['EC']['ID Growth Curve'][iT,indS[i],iE]=1    
    return vi

#%% Simulate probability of harvesting on the fly
def PredictHarvesting_OnTheFly(meta,pNam,vi,iT,iScn,iEns,V_Merch,Period):        

    # Indicator of THLB (THLB=1, Non-THLB=0)
    flag_thlb=vi['lsat']['THLB'][iT,:]

    # # Indicator of energy production (don't harvest on the fly, it is pre-scheduled)
    # if meta[pNam]['Project']['Land Surface Class Dependent']!='No':
    #     iT_lsc=np.where(vi['lsat']['LSC']['tv']==meta[pNam]['Year'][iT])[0]
    #     if iT_lsc.size>0:
    #         flag_ep=1*( (vi['lsat']['LSC']['Use'][iT_lsc,:]==meta['LUT']['LSC']['Use']['Fuel Break']) | (vi['lsat']['LSC']['Use'][iT_lsc,:]==meta['LUT']['LSC']['Use']['Energy Production']) )
    #     else:
    #         flag_ep=np.ones(flag_thlb.size,dtype=int)
    # else:
    #     flag_ep=np.ones(flag_thlb.size,dtype=int)

    # Saturation parameter
    if Period=='Historical':

        # Historical period

        #bH=[0.0007,5.05,0.32,1975]
        bH=[0.00085,5.15,0.32,1975]
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
            Pa_H_Sat=meta[pNam]['Scenario'][iScn]['Pa Harvest Sat']*vi['lsat']['Harvest Index']
        else:
            Pa_H_Sat=meta['Param']['BE']['On The Fly']['Pa_Harvest_Sat']*vi['lsat']['Harvest Index']

    # Inflection point
    if Period=='Historical':
        Pa_H_Inf=500
    else:
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
    rn=np.random.random(V_Merch.size)
    #rn=meta[pNam]['Project']['On the Fly']['Random Numbers']['Harvest'][iT,:]

    # if Period=='Historical':
    #     # This will make sure the historical harvest simulations don't occur on
    #     # the footprint of recorded historical harvesting
    #     flag_HistHarv=-1*(np.minimum(1,vi['lsat']['Year Harvest First'])-1)        
    # else:
    #     flag_HistHarv=np.ones(rn.size)

    # Occurrence
    #Oc=flag_ep*flag_thlb*flag_HistHarv*np.floor(np.minimum(1,Po/rn))
    #Oc=flag_ep*flag_thlb*np.floor(np.minimum(1,Po/rn))
    Oc=flag_thlb*np.floor(np.minimum(1,Po/rn))

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
                try:
                    # This can crash in the last 2 time steps so try it
                    iE=iAvailable[0]+1 # changing this to zero will cause the harvest to be overwritten
                    vi['EC']['ID Event Type'][iT,indS[i],iE]=meta['LUT']['Event']['Slashpile Burn']
                    vi['EC']['Mortality Factor'][iT,indS[i],iE]=1.0
                    vi['EC']['ID Growth Curve'][iT,indS[i],iE]=2
                    iE=iAvailable[0]+2
                    vi['EC']['ID Event Type'][iT,indS[i],iE]=meta['LUT']['Event']['Planting']
                    vi['EC']['Mortality Factor'][iT,indS[i],iE]=1.0
                    vi['EC']['ID Growth Curve'][iT,indS[i],iE]=2
                except:
                    pass

    return vi

