'''
HWP END-USE SCENARIOS

'''

#%% Import modules

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from fcgadgets.macgyver import utilities_general as gu

#%% Graphics parameters

fs=7
params={'font.sans-serif':'Arial','font.size':fs,'axes.edgecolor':'black','axes.labelsize':fs,'axes.labelcolor':'black','axes.titlesize':fs,'axes.linewidth':0.5,'lines.linewidth':0.5,
        'text.color':'black','xtick.color':'black','xtick.labelsize':fs,'xtick.major.width':0.5,'xtick.major.size':3,'xtick.direction':'in','ytick.color':'black','ytick.labelsize':fs,
        'ytick.major.width':0.5,'ytick.major.size':3,'ytick.direction':'in','legend.fontsize':fs,'savefig.dpi':300,'savefig.transparent':True}
plt.rcParams.update(params)

#%% Import BaseCase HWP parameters

dDef0=gu.ReadExcel(r'C:\Users\rhember\Documents\Code_Python\fcgadgets\cbrunner\Parameters\Parameters_HWP.xlsx','Default')

# Variables to keep
vrL=['SawMillToPulpMill','SawMillToMDFMill','SawMillToPelletMill','SawMillToPowerFacility','SawMillToIPP','SawMillToLogExport','SawMillToSFH','SawMillToMFH','SawMillToCom',
  'SawMillToFurn','SawMillToShip','SawMillToRepairs','SawMillToOther','PlywoodMillToPowerFacility','PlywoodMillToIPP','PlywoodMillToSFH','PlywoodMillToMFH','PlywoodMillToCom',
  'PlywoodMillToFurn','PlywoodMillToShip','PlywoodMillToRepairs','PlywoodMillToOther','OSBMillToPowerFacility','OSBMillToIPP','OSBMillToSFH','OSBMillToMFH','OSBMillToCom',
  'OSBMillToFurn','OSBMillToShip','OSBMillToRepairs','OSBMillToOther','MDFMillToSFH','MDFMillToMFH','MDFMillToCom','MDFMillToFurn','MDFMillToShip','MDFMillToRepairs','MDFMillToOther']

# *** If you change the time vector, change it consistently for other variable scenarios - it needs to be consistent ***
tv=np.arange(1850,2101,1)

#%% Initialize dictionary 

regL=['Coast','Interior','GFS22']
scnL=['BaseCase','S1','S2','S3']

d={}
d['Year']=tv
for scn in scnL:
    d[scn]={}
    for reg in regL:
        d[scn][reg]={}
        for i in range(dDef0['Name'].size):
            if np.isin(dDef0['Name'][i],vrL)==True:
                d[scn][reg][ dDef0['Name'][i] ]=np.ones(tv.size)*dDef0[reg][i]

#%% Scenario 1

for reg in regL:
    
    #------------------------------------------------------------------------------
    # Define change in fraction of building materials
    #------------------------------------------------------------------------------
            
    # Assumed change in building output
    tStart=2022
    tEnd=2030
    
    # Change in building materials
    DeltaBuildingsTot=1.1
    DeltaBuildings=np.ones(tv.size)
    
    iT=np.where( (tv>=tStart) & (tv<=tEnd) )[0]
    DeltaBuildings[iT]=np.linspace(1,DeltaBuildingsTot,iT.size)
    iT=np.where( (tv>tEnd) )[0]
    DeltaBuildings[iT]=DeltaBuildingsTot
    
    # Plot assumptions
    #plt.plot(tv,DeltaBuildings,'b-')
    
    #------------------------------------------------------------------------------
    # Sawmills
    #------------------------------------------------------------------------------
    
    nm='SawMillTo'
    
    # Pool groups 
    BuildingL=['SFH','MFH','Com']
    NonBuildingL=['PulpMill','MDFMill','PelletMill','PowerFacility','IPP','LogExport','Furn','Ship','Repairs','Other']
    
    # Adjust building pools
    TotalChangeInBuildingTrans=np.zeros(tv.size)
    for ip in BuildingL:
        x=d['BaseCase'][reg][nm + ip]
        d['S1'][reg][nm + ip]=DeltaBuildings*x
        TotalChangeInBuildingTrans=TotalChangeInBuildingTrans+DeltaBuildings*x-DeltaBuildings[0]*x
    #plt.plot(tv,TotalChangeInBuildingTrans,'-')
    
    # Adjust non-building pools to compensate
    for ip in NonBuildingL:
        x=d['BaseCase'][reg][nm + ip]
        d['S1'][reg][nm + ip]=x-TotalChangeInBuildingTrans/len(NonBuildingL)
        
    # QA (check conservation of mass)
    y=np.zeros(tv.size)
    for ip in BuildingL:
        y=y+d['S1'][reg][nm + ip]
    for ip in NonBuildingL:
        y=y+d['S1'][reg][nm + ip]
    
    flg=0
    if flg==1:
        plt.close('all')
        plt.plot(tv,y,'-')
    
    #------------------------------------------------------------------------------
    # Plywood mills
    #------------------------------------------------------------------------------
        
    nm='PlywoodMillTo'
    
    # Pool groups 
    BuildingL=['SFH','MFH','Com']
    NonBuildingL=['PowerFacility','IPP','Furn','Ship','Repairs','Other']
    
    # Adjust building pools
    TotalChangeInBuildingTrans=np.zeros(tv.size)
    for ip in BuildingL:
        x=d['BaseCase'][reg][nm + ip]
        d['S1'][reg][nm + ip]=DeltaBuildings*x
        TotalChangeInBuildingTrans=TotalChangeInBuildingTrans+DeltaBuildings*x-DeltaBuildings[0]*x
    #plt.plot(tv,TotalChangeInBuildingTrans,'-')
    
    # Adjust non-building pools to compensate
    for ip in NonBuildingL:
        x=d['BaseCase'][reg][nm + ip]
        d['S1'][reg][nm + ip]=x-TotalChangeInBuildingTrans/len(NonBuildingL)
        
    # QA (check conservation of mass)
    y=np.zeros(tv.size)
    for ip in BuildingL:
        y=y+d['S1'][reg][nm + ip]
    for ip in NonBuildingL:
        y=y+d['S1'][reg][nm + ip]
    
    flg=0
    if flg==1:
        plt.close('all')
        plt.plot(tv,y,'-')
    
    #------------------------------------------------------------------------------
    # OSB mills
    #------------------------------------------------------------------------------
        
    nm='OSBMillTo'
    
    # Pool groups 
    BuildingL=['SFH','MFH','Com']
    NonBuildingL=['PowerFacility','IPP','Furn','Ship','Repairs','Other']
    
    # Adjust building pools
    TotalChangeInBuildingTrans=np.zeros(tv.size)
    for ip in BuildingL:
        x=d['BaseCase'][reg][nm + ip]
        d['S1'][reg][nm + ip]=DeltaBuildings*x
        TotalChangeInBuildingTrans=TotalChangeInBuildingTrans+DeltaBuildings*x-DeltaBuildings[0]*x
    #plt.plot(tv,TotalChangeInBuildingTrans,'-')
    
    # Adjust non-building pools to compensate
    for ip in NonBuildingL:
        x=d['BaseCase'][reg][nm + ip]
        d['S1'][reg][nm + ip]=x-TotalChangeInBuildingTrans/len(NonBuildingL)
        
    # QA (check conservation of mass)
    y=np.zeros(tv.size)
    for ip in BuildingL:
        y=y+d['S1'][reg][nm + ip]
    for ip in NonBuildingL:
        y=y+d['S1'][reg][nm + ip]
    
    flg=0
    if flg==1:
        plt.close('all')
        plt.plot(tv,y,'-')
    
    #------------------------------------------------------------------------------
    # MDF mills
    #------------------------------------------------------------------------------
        
    nm='MDFMillTo'
    
    # Pool groups 
    BuildingL=['SFH','MFH','Com']
    NonBuildingL=['Furn','Ship','Repairs','Other']
    
    # Adjust building pools
    TotalChangeInBuildingTrans=np.zeros(tv.size)
    for ip in BuildingL:
        x=d['BaseCase'][reg][nm + ip]
        d['S1'][reg][nm + ip]=DeltaBuildings*x
        TotalChangeInBuildingTrans=TotalChangeInBuildingTrans+DeltaBuildings*x-DeltaBuildings[0]*x
    #plt.plot(tv,TotalChangeInBuildingTrans,'-')
    
    # Adjust non-building pools to compensate
    for ip in NonBuildingL:
        x=d['BaseCase'][reg][nm + ip]
        d['S1'][reg][nm + ip]=x-TotalChangeInBuildingTrans/len(NonBuildingL)
        
    # QA (check conservation of mass)
    y=np.zeros(tv.size)
    for ip in BuildingL:
        y=y+d['S1'][reg][nm + ip]
    for ip in NonBuildingL:
        y=y+d['S1'][reg][nm + ip]
    
    flg=0
    if flg==1:
        plt.close('all')
        plt.plot(tv,y,'-')

#%% Scenario 2

for reg in regL:

    #------------------------------------------------------------------------------
    # Define change in fraction of building materials
    #------------------------------------------------------------------------------
            
    # Assumed change in building output
    tStart=2022
    tEnd=2030
    
    # Change in building materials
    DeltaBuildingsTot=2.0
    DeltaBuildings=np.ones(tv.size)
    
    iT=np.where( (tv>=tStart) & (tv<=tEnd) )[0]
    DeltaBuildings[iT]=np.linspace(1,DeltaBuildingsTot,iT.size)
    iT=np.where( (tv>tEnd) )[0]
    DeltaBuildings[iT]=DeltaBuildingsTot
    
    # Plot assumptions
    #plt.plot(tv,DeltaBuildings,'b-')
    
    #------------------------------------------------------------------------------
    # Sawmills
    #------------------------------------------------------------------------------
    
    nm='SawMillTo'
    
    # Pool groups 
    BuildingL=['SFH','MFH','Com']
    NonBuildingL=['PulpMill','MDFMill','PelletMill','PowerFacility','IPP','LogExport','Furn','Ship','Repairs','Other']
    
    # Adjust building pools
    TotalChangeInBuildingTrans=np.zeros(tv.size)
    for ip in BuildingL:
        x=d['BaseCase'][reg][nm + ip]
        d['S2'][reg][nm + ip]=DeltaBuildings*x
        TotalChangeInBuildingTrans=TotalChangeInBuildingTrans+DeltaBuildings*x-DeltaBuildings[0]*x
    #plt.plot(tv,TotalChangeInBuildingTrans,'-')
    
    # Adjust non-building pools to compensate
    for ip in NonBuildingL:
        x=d['BaseCase'][reg][nm + ip]
        d['S2'][reg][nm + ip]=x-TotalChangeInBuildingTrans/len(NonBuildingL)
        
    # QA (check conservation of mass)
    y=np.zeros(tv.size)
    for ip in BuildingL:
        y=y+d['S2'][reg][nm + ip]
    for ip in NonBuildingL:
        y=y+d['S2'][reg][nm + ip]
    
    flg=0
    if flg==1:
        plt.close('all')
        plt.plot(tv,y,'-')
    
    #------------------------------------------------------------------------------
    # Plywood mills
    #------------------------------------------------------------------------------
        
    nm='PlywoodMillTo'
    
    # Pool groups 
    BuildingL=['SFH','MFH','Com']
    NonBuildingL=['PowerFacility','IPP','Furn','Ship','Repairs','Other']
    
    # Adjust building pools
    TotalChangeInBuildingTrans=np.zeros(tv.size)
    for ip in BuildingL:
        x=d['BaseCase'][reg][nm + ip]
        d['S2'][reg][nm + ip]=DeltaBuildings*x
        TotalChangeInBuildingTrans=TotalChangeInBuildingTrans+DeltaBuildings*x-DeltaBuildings[0]*x
    #plt.plot(tv,TotalChangeInBuildingTrans,'-')
    
    # Adjust non-building pools to compensate
    for ip in NonBuildingL:
        x=d['BaseCase'][reg][nm + ip]
        d['S2'][reg][nm + ip]=x-TotalChangeInBuildingTrans/len(NonBuildingL)
        
    # QA (check conservation of mass)
    y=np.zeros(tv.size)
    for ip in BuildingL:
        y=y+d['S2'][reg][nm + ip]
    for ip in NonBuildingL:
        y=y+d['S2'][reg][nm + ip]
    
    flg=0
    if flg==1:
        plt.close('all')
        plt.plot(tv,y,'-')
    
    #------------------------------------------------------------------------------
    # OSB mills
    #------------------------------------------------------------------------------
        
    nm='OSBMillTo'
    
    # Pool groups 
    BuildingL=['SFH','MFH','Com']
    NonBuildingL=['PowerFacility','IPP','Furn','Ship','Repairs','Other']
    
    # Adjust building pools
    TotalChangeInBuildingTrans=np.zeros(tv.size)
    for ip in BuildingL:
        x=d['BaseCase'][reg][nm + ip]
        d['S2'][reg][nm + ip]=DeltaBuildings*x
        TotalChangeInBuildingTrans=TotalChangeInBuildingTrans+DeltaBuildings*x-DeltaBuildings[0]*x
    #plt.plot(tv,TotalChangeInBuildingTrans,'-')
    
    # Adjust non-building pools to compensate
    for ip in NonBuildingL:
        x=d['BaseCase'][reg][nm + ip]
        d['S2'][reg][nm + ip]=x-TotalChangeInBuildingTrans/len(NonBuildingL)
        
    # QA (check conservation of mass)
    y=np.zeros(tv.size)
    for ip in BuildingL:
        y=y+d['S2'][reg][nm + ip]
    for ip in NonBuildingL:
        y=y+d['S2'][reg][nm + ip]
    
    flg=0
    if flg==1:
        plt.close('all')
        plt.plot(tv,y,'-')
    
    #------------------------------------------------------------------------------
    # MDF mills
    #------------------------------------------------------------------------------
        
    nm='MDFMillTo'
    
    # Pool groups 
    BuildingL=['SFH','MFH','Com']
    NonBuildingL=['Furn','Ship','Repairs','Other']
    
    # Adjust building pools
    TotalChangeInBuildingTrans=np.zeros(tv.size)
    for ip in BuildingL:
        x=d['BaseCase'][reg][nm + ip]
        d['S2'][reg][nm + ip]=DeltaBuildings*x
        TotalChangeInBuildingTrans=TotalChangeInBuildingTrans+DeltaBuildings*x-DeltaBuildings[0]*x
    #plt.plot(tv,TotalChangeInBuildingTrans,'-')
    
    # Adjust non-building pools to compensate
    for ip in NonBuildingL:
        x=d['BaseCase'][reg][nm + ip]
        d['S2'][reg][nm + ip]=x-TotalChangeInBuildingTrans/len(NonBuildingL)
        
    # QA (check conservation of mass)
    y=np.zeros(tv.size)
    for ip in BuildingL:
        y=y+d['S2'][reg][nm + ip]
    for ip in NonBuildingL:
        y=y+d['S2'][reg][nm + ip]
    
    flg=0
    if flg==1:
        plt.close('all')
        plt.plot(tv,y,'-')

#%% Scenario 3

for reg in regL:
    
    #------------------------------------------------------------------------------
    # Define change in fraction of building materials
    #------------------------------------------------------------------------------
            
    # Assumed change in building output
    tStart=2022
    tEnd=2030
    
    # Change in building materials
    DeltaBuildingsTot=3.0
    DeltaBuildings=np.ones(tv.size)
    
    iT=np.where( (tv>=tStart) & (tv<=tEnd) )[0]
    DeltaBuildings[iT]=np.linspace(1,DeltaBuildingsTot,iT.size)
    iT=np.where( (tv>tEnd) )[0]
    DeltaBuildings[iT]=DeltaBuildingsTot
    
    # Plot assumptions
    #plt.plot(tv,DeltaBuildings,'b-')
    
    #------------------------------------------------------------------------------
    # Sawmills
    #------------------------------------------------------------------------------
    
    nm='SawMillTo'
    
    # Pool groups 
    BuildingL=['SFH','MFH','Com']
    NonBuildingL=['PulpMill','MDFMill','PelletMill','PowerFacility','IPP','LogExport','Furn','Ship','Repairs','Other']
    
    # Adjust building pools
    TotalChangeInBuildingTrans=np.zeros(tv.size)
    for ip in BuildingL:
        x=d['BaseCase'][reg][nm + ip]
        d['S3'][reg][nm + ip]=DeltaBuildings*x
        TotalChangeInBuildingTrans=TotalChangeInBuildingTrans+DeltaBuildings*x-DeltaBuildings[0]*x
    #plt.plot(tv,TotalChangeInBuildingTrans,'-')
    
    # Adjust non-building pools to compensate
    for ip in NonBuildingL:
        x=d['BaseCase'][reg][nm + ip]
        d['S3'][reg][nm + ip]=x-TotalChangeInBuildingTrans/len(NonBuildingL)
        
    # QA (check conservation of mass)
    y=np.zeros(tv.size)
    for ip in BuildingL:
        y=y+d['S3'][reg][nm + ip]
    for ip in NonBuildingL:
        y=y+d['S3'][reg][nm + ip]
    
    flg=0
    if flg==1:
        plt.close('all')
        plt.plot(tv,y,'-')
    
    #------------------------------------------------------------------------------
    # Plywood mills
    #------------------------------------------------------------------------------
        
    nm='PlywoodMillTo'
    
    # Pool groups 
    BuildingL=['SFH','MFH','Com']
    NonBuildingL=['PowerFacility','IPP','Furn','Ship','Repairs','Other']
    
    # Adjust building pools
    TotalChangeInBuildingTrans=np.zeros(tv.size)
    for ip in BuildingL:
        x=d['BaseCase'][reg][nm + ip]
        d['S3'][reg][nm + ip]=DeltaBuildings*x
        TotalChangeInBuildingTrans=TotalChangeInBuildingTrans+DeltaBuildings*x-DeltaBuildings[0]*x
    #plt.plot(tv,TotalChangeInBuildingTrans,'-')
    
    # Adjust non-building pools to compensate
    for ip in NonBuildingL:
        x=d['BaseCase'][reg][nm + ip]
        d['S3'][reg][nm + ip]=x-TotalChangeInBuildingTrans/len(NonBuildingL)
        
    # QA (check conservation of mass)
    y=np.zeros(tv.size)
    for ip in BuildingL:
        y=y+d['S3'][reg][nm + ip]
    for ip in NonBuildingL:
        y=y+d['S3'][reg][nm + ip]
    
    flg=0
    if flg==1:
        plt.close('all')
        plt.plot(tv,y,'-')
    
    #------------------------------------------------------------------------------
    # OSB mills
    #------------------------------------------------------------------------------
        
    nm='OSBMillTo'
    
    # Pool groups 
    BuildingL=['SFH','MFH','Com']
    NonBuildingL=['PowerFacility','IPP','Furn','Ship','Repairs','Other']
    
    # Adjust building pools
    TotalChangeInBuildingTrans=np.zeros(tv.size)
    for ip in BuildingL:
        x=d['BaseCase'][reg][nm + ip]
        d['S3'][reg][nm + ip]=DeltaBuildings*x
        TotalChangeInBuildingTrans=TotalChangeInBuildingTrans+DeltaBuildings*x-DeltaBuildings[0]*x
    #plt.plot(tv,TotalChangeInBuildingTrans,'-')
    
    # Adjust non-building pools to compensate
    for ip in NonBuildingL:
        x=d['BaseCase'][reg][nm + ip]
        d['S3'][reg][nm + ip]=x-TotalChangeInBuildingTrans/len(NonBuildingL)
        
    # QA (check conservation of mass)
    y=np.zeros(tv.size)
    for ip in BuildingL:
        y=y+d['S3'][reg][nm + ip]
    for ip in NonBuildingL:
        y=y+d['S3'][reg][nm + ip]
    
    flg=0
    if flg==1:
        plt.close('all')
        plt.plot(tv,y,'-')
    
    #------------------------------------------------------------------------------
    # MDF mills
    #------------------------------------------------------------------------------
        
    nm='MDFMillTo'
    
    # Pool groups 
    BuildingL=['SFH','MFH','Com']
    NonBuildingL=['Furn','Ship','Repairs','Other']
    
    # Adjust building pools
    TotalChangeInBuildingTrans=np.zeros(tv.size)
    for ip in BuildingL:
        x=d['BaseCase'][reg][nm + ip]
        d['S3'][reg][nm + ip]=DeltaBuildings*x
        TotalChangeInBuildingTrans=TotalChangeInBuildingTrans+DeltaBuildings*x-DeltaBuildings[0]*x
    #plt.plot(tv,TotalChangeInBuildingTrans,'-')
    
    # Adjust non-building pools to compensate
    for ip in NonBuildingL:
        x=d['BaseCase'][reg][nm + ip]
        d['S3'][reg][nm + ip]=x-TotalChangeInBuildingTrans/len(NonBuildingL)
        
    # QA (check conservation of mass)
    y=np.zeros(tv.size)
    for ip in BuildingL:
        y=y+d['S3'][reg][nm + ip]
    for ip in NonBuildingL:
        y=y+d['S3'][reg][nm + ip]
    
    flg=0
    if flg==1:
        plt.close('all')
        plt.plot(tv,y,'-')

#%% Save

gu.opickle(r'C:\Users\rhember\Documents\Code_Python\fcgadgets\cbrunner\Parameters\Variables_HWP_EndUses.pkl',d)


#%% Plot

reg='GFS22'

yBD=np.zeros(tv.size)
for ip in BuildingL:
    yBD=yBD+d['BaseCase'][reg][nm + ip]
yBS1=np.zeros(tv.size)
for ip in BuildingL:
    yBS1=yBS1+d['S1'][reg][nm + ip]
yBS2=np.zeros(tv.size)
for ip in BuildingL:
    yBS2=yBS2+d['S2'][reg][nm + ip]    

yNBD=np.zeros(tv.size)
for ip in NonBuildingL:
    yNBD=yNBD+d['BaseCase'][reg][nm + ip]
yNBS1=np.zeros(tv.size)
for ip in NonBuildingL:
    yNBS1=yNBS1+d['S1'][reg][nm + ip]
yNBS2=np.zeros(tv.size)
for ip in NonBuildingL:
    yNBS2=yNBS2+d['S2'][reg][nm + ip]    
 
plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(16,6.5)); ms=3; lw=1.25
ax[0].plot(tv,yBD,'b-',label='BaseCase',lw=lw)
ax[0].plot(tv,yBS1,'g--',label='Scenario 1',lw=lw)
ax[0].plot(tv,yBS2,'r-.',label='Scenario 2',lw=lw)
ax[0].set(ylim=[0,1],yticks=np.arange(0,1.2,0.2),ylabel='Building material fraction',xlim=[1900,2100],xticks=np.arange(1500,2200,25),xlabel='Time, years')
ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both')
ax[0].legend(loc='lower left',facecolor=[1,1,1],frameon=False)
ax[1].plot(tv,yNBD,'b-',label='BaseCase',lw=lw)
ax[1].plot(tv,yNBS1,'g--',label='Scenario 1',lw=lw)
ax[1].plot(tv,yNBS2,'r-.',label='Scenario 2',lw=lw)
ax[1].set(ylim=[0,1],yticks=np.arange(0,1.2,0.2),ylabel='Non-building material fraction',xlim=[1900,2100],xticks=np.arange(1500,2200,25),xlabel='Time, years')
ax[1].yaxis.set_ticks_position('both'); ax[1].xaxis.set_ticks_position('both')

gu.axletters(ax,plt,0.028,0.9,LetterStyle='Caps',FontWeight='Bold') # 

fig.tight_layout()
gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Harvest\End Use Scenarios\End Use Scenarios ' + reg,'png',900)
