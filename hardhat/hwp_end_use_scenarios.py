'''
HWP END-USE SCENARIOS

'''

#%% Import modules

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from fcgadgets.macgyver import utilities_general as gu

#%% Graphics parameters

gp=gu.SetGraphics('Manuscript')

#%% Import BaseCase HWP parameters

dDef0=gu.ReadExcel(r'C:\Users\rhember\Documents\Code_Python\fcgadgets\cbrunner\Parameters\Parameters_HWP.xlsx','Main')

# Variables to keep
vrL=['SawMillToPulpMill','SawMillToMDFMill','SawMillToPelletMill','SawMillToPowerFacility','SawMillToIPP','SawMillToLogExport','SawMillToSFH','SawMillToMFH','SawMillToCom','SawMillToFurn','SawMillToShip','SawMillToRepairs','SawMillToOther',
  'PlywoodMillToPowerFacility','PlywoodMillToIPP','PlywoodMillToSFH','PlywoodMillToMFH','PlywoodMillToCom','PlywoodMillToFurn','PlywoodMillToShip','PlywoodMillToRepairs','PlywoodMillToOther',
  'OSBMillToPowerFacility','OSBMillToIPP','OSBMillToSFH','OSBMillToMFH','OSBMillToCom','OSBMillToFurn','OSBMillToShip','OSBMillToRepairs','OSBMillToOther',
  'MDFMillToSFH','MDFMillToMFH','MDFMillToCom','MDFMillToFurn','MDFMillToShip','MDFMillToRepairs','MDFMillToOther']

# *** If you change the time vector, change it consistently for other variable scenarios - it needs to be consistent ***
tv=np.arange(1850,2101,1)

#%% Initialize dictionary

# Regions
regL=['Coast',
      'Interior',
      'GFS22',
      'Burn Uneconomic',
      'Leave Uneconomic Standing',
      'Sawlogs and Pellets',
      'Energy Production Pellets']

# Scenarios
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
    tStart=2021
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
    tStart=2021
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

gu.opickle(r'C:\Users\rhember\Documents\Code_Python\fcgadgets\cbrunner\Parameters\Variables_HWP_EndUse.pkl',d)

