
import numpy as np
import pandas as pd
import pickle
import matplotlib.pyplot as plt
from fcgadgets.pyscripts import utilities_general as gu
from fcgadgets.cbrunner.cbrun_utilities import *


'''============================================================================
PLOT DISTURBANCE RATE TIME SERIES
============================================================================'''

def PlotDisturbanceRates(meta,iScn,tStart):

    for iBat in range(meta['N Batch']):
        
        iEns=0
        
        # Import disturbance history
        pth=meta['Path Input Scenario'][iScn] + '\\Disturbance_Ens' + FixFileNum(iEns) + '_Bat' + FixFileNum(iBat) + '.pkl'
        fin=open(pth,'rb')
        dh=pickle.load(fin); fin.close()
    
        DistTypes=list(meta['LUT Dist'].keys())
        DistID=np.array(list(meta['LUT Dist'].values()))        
        
        # Initialize area disturbed by disturbance type (hectares)
        A=np.zeros((meta['Time'].size,len(DistTypes)))
    
        # Loop through time intervals
        for i in range(len(dh)):        
            for j in range(dh[i]['Year'].size):            
                iT=np.where(meta['Time']==dh[i]['Year'][j])[0]            
                iD=np.where(DistID==dh[i]['ID_Type'][j])[0]
                A[iT,iD]=A[iT,iD]+1                
        
        AR=A/meta['N Stand']*100 # Convert to percent
    
    AR_ma=gu.movingave(AR,10,'Historical')
    
    # Plot
    it=np.where(meta['Time']>=tStart)[0]
                
    plt.close('all')
    fig,ax=plt.subplots(3,1,figsize=gu.cm2inch(15,12))
    # Wildfire
    ax[0].bar(meta['Time'][it],AR[it,0],width=1,color=[0.8,0,0])
    ax[0].plot(meta['Time'][it],AR_ma[it,0],linewidth=1.25,color=[0.5,0,0])
    ax[0].set(position=[0.07,0.7,0.9,0.27],xlim=[tStart,meta['Time'][-1]+1],ylabel='Area disturbed (%)')
    # Insects
    ax[1].bar(meta['Time'][it],AR[it,4]+AR[it,5],width=1,color=[0.8,0,0])
    ax[1].plot(meta['Time'][it],AR_ma[it,4]+AR_ma[it,5],linewidth=1.25,color=[0.5,0,0])
    ax[1].set(position=[0.07,0.38,0.9,0.27],xlim=[tStart,meta['Time'][-1]+1],ylabel='Area disturbed (%)')
    # Harvest
    ax[2].bar(meta['Time'][it],AR[it,1]+AR[it,2]+AR[it,7],width=1,color=[0.8,0,0])
    ax[2].plot(meta['Time'][it],AR_ma[it,1]+AR_ma[it,2]+AR_ma[it,7],linewidth=1.25,color=[0.5,0,0])
    ax[2].set(position=[0.07,0.07,0.9,0.27],xlim=[tStart,meta['Time'][-1]+1],ylabel='Area disturbed (%)',xlabel='Time, years')
    gu.axletters(ax,plt,0.02,0.87,['Wildfire','Insects','Harvesting'],0.03)
    fig.savefig(meta['Path Project'] + '\\Outputs\\Figures\\ts_mean_Disturbance1.png',format='png',dpi=600)
    
    # Ten-year moving average
    plt.close('all')
    fig,ax=plt.subplots(1,figsize=gu.cm2inch(15,7))
    plt.plot(meta['Time'][it],AR_ma[it,0],color=[0.75,0,0],linewidth=2,label='Wildfire')
    plt.plot(meta['Time'][it],AR_ma[it,1]+AR_ma[it,2],color=[0.27,0.49,0.77],linewidth=2,label='Harvest')
    plt.plot(meta['Time'][it],AR_ma[it,4],color=[0,0.8,0],linewidth=2,label='Beetles')
    plt.plot(meta['Time'][it],AR_ma[it,5],color=[0.5,1,0],linewidth=2,label='Defoliators')
    ax.set(xlim=[tStart,meta['Time'][-1]+1],ylabel='Area disturbed (%)',xlabel='Time, years')
    plt.legend(loc='upper right')
    plt.grid()
    fig.savefig(meta['Path Project'] + '\\Outputs\\Figures\\ts_mean_Disturbance2.png',format='png',dpi=600)

    return

