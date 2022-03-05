'''
HWP END-USE SCENARIOS

'''

#%% Import modules

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from fcgadgets.macgyver import utilities_general as gu

#%% Import default HWP parameters

dDef0=gu.ReadExcel(r'C:\Users\rhember\Documents\Code_Python\fcgadgets\cbrunner\Parameters\Parameters_HWP.xlsx','Default')

# Variables to keep
vrL=['SawMillToPulpMill','SawMillToMDFMill','SawMillToPelletMill','SawMillToPowerFacility','SawMillToIPP','SawMillToLogExport','SawMillToSFH','SawMillToMFH','SawMillToCom',
  'SawMillToFurn','SawMillToShip','SawMillToRepairs','SawMillToOther','PlywoodMillToPowerFacility','PlywoodMillToIPP','PlywoodMillToSFH','PlywoodMillToMFH','PlywoodMillToCom',
  'PlywoodMillToFurn','PlywoodMillToShip','PlywoodMillToRepairs','PlywoodMillToOther','OSBMillToPowerFacility','OSBMillToIPP','OSBMillToSFH','OSBMillToMFH','OSBMillToCom',
  'OSBMillToFurn','OSBMillToShip','OSBMillToRepairs','OSBMillToOther','MDFMillToSFH','MDFMillToMFH','MDFMillToCom','MDFMillToFurn','MDFMillToShip','MDFMillToRepairs','MDFMillToOther']

tv=np.arange(1850,2101,1)

dDef={}
for i in range(dDef0['Name'].size):
    if np.isin(dDef0['Name'][i],vrL)==True:
        dDef[ dDef0['Name'][i] ]=np.ones(tv.size)*dDef0['Best Estimate'][i]
del dDef0  

#%% Initialize database

dS1={}
for k in dDef.keys():
    dS1[k]=np.ones(tv.size)*dDef[k]

# Assumed change in building output
DeltaBuildingsTot=3.0
tStart=2022
tEnd=2030
DeltaBuildings=np.ones(tv.size)

iT=np.where( (tv>=tStart) & (tv<=tEnd) )[0]
DeltaBuildings[iT]=np.linspace(1,DeltaBuildingsTot,iT.size)
iT=np.where( (tv>tEnd) )[0]
DeltaBuildings[iT]=DeltaBuildingsTot

# Plot assumptions
#plt.plot(tv,DeltaBuildings,'b-')

plt.plot(tv,1/DeltaBuildings,'b-')

#%% Sawmills

nm='SawMillTo'

# Pool groups 
BuildingL=['SFH','MFH','Com']
NonBuildingL=['PulpMill','MDFMill','PelletMill','PowerFacility','IPP','LogExport','Furn','Ship','Repairs','Other']

# Adjust building pools
TotalChangeInBuildingTrans=np.zeros(tv.size)
for ip in BuildingL:
    x=dDef[nm + ip]
    dS1[nm + ip]=DeltaBuildings*x
    TotalChangeInBuildingTrans=TotalChangeInBuildingTrans+DeltaBuildings*x-DeltaBuildings[0]*x
#plt.plot(tv,TotalChangeInBuildingTrans,'-')

# Adjust non-building pools to compensate
for ip in NonBuildingL:
    x=dDef[nm + ip]
    dS1[nm + ip]=x-TotalChangeInBuildingTrans/len(NonBuildingL)
    
# QA (check conservation of mass)
y=np.zeros(tv.size)
for ip in BuildingL:
    y=y+dS1[nm + ip]
for ip in NonBuildingL:
    y=y+dS1[nm + ip]

flg=0
if flg==1:
    plt.close('all')
    plt.plot(tv,y,'-')

#%% Plywood mills

nm='PlywoodMillTo'

# Pool groups 
BuildingL=['SFH','MFH','Com']
NonBuildingL=['PowerFacility','IPP','Furn','Ship','Repairs','Other']

# Adjust building pools
TotalChangeInBuildingTrans=np.zeros(tv.size)
for ip in BuildingL:
    x=dDef[nm + ip]
    dS1[nm + ip]=DeltaBuildings*x
    TotalChangeInBuildingTrans=TotalChangeInBuildingTrans+DeltaBuildings*x-DeltaBuildings[0]*x
#plt.plot(tv,TotalChangeInBuildingTrans,'-')

# Adjust non-building pools to compensate
for ip in NonBuildingL:
    x=dDef[nm + ip]
    dS1[nm + ip]=x-TotalChangeInBuildingTrans/len(NonBuildingL)
    
# QA (check conservation of mass)
y=np.zeros(tv.size)
for ip in BuildingL:
    y=y+dS1[nm + ip]
for ip in NonBuildingL:
    y=y+dS1[nm + ip]

flg=0
if flg==1:
    plt.close('all')
    plt.plot(tv,y,'-')

#%% OSB mills

nm='OSBMillTo'

# Pool groups 
BuildingL=['SFH','MFH','Com']
NonBuildingL=['PowerFacility','IPP','Furn','Ship','Repairs','Other']

# Adjust building pools
TotalChangeInBuildingTrans=np.zeros(tv.size)
for ip in BuildingL:
    x=dDef[nm + ip]
    dS1[nm + ip]=DeltaBuildings*x
    TotalChangeInBuildingTrans=TotalChangeInBuildingTrans+DeltaBuildings*x-DeltaBuildings[0]*x
#plt.plot(tv,TotalChangeInBuildingTrans,'-')

# Adjust non-building pools to compensate
for ip in NonBuildingL:
    x=dDef[nm + ip]
    dS1[nm + ip]=x-TotalChangeInBuildingTrans/len(NonBuildingL)
    
# QA (check conservation of mass)
y=np.zeros(tv.size)
for ip in BuildingL:
    y=y+dS1[nm + ip]
for ip in NonBuildingL:
    y=y+dS1[nm + ip]

flg=0
if flg==1:
    plt.close('all')
    plt.plot(tv,y,'-')

#%% MDF mills

nm='MDFMillTo'

# Pool groups 
BuildingL=['SFH','MFH','Com']
NonBuildingL=['Furn','Ship','Repairs','Other']

# Adjust building pools
TotalChangeInBuildingTrans=np.zeros(tv.size)
for ip in BuildingL:
    x=dDef[nm + ip]
    dS1[nm + ip]=DeltaBuildings*x
    TotalChangeInBuildingTrans=TotalChangeInBuildingTrans+DeltaBuildings*x-DeltaBuildings[0]*x
#plt.plot(tv,TotalChangeInBuildingTrans,'-')

# Adjust non-building pools to compensate
for ip in NonBuildingL:
    x=dDef[nm + ip]
    dS1[nm + ip]=x-TotalChangeInBuildingTrans/len(NonBuildingL)
    
# QA (check conservation of mass)
y=np.zeros(tv.size)
for ip in BuildingL:
    y=y+dS1[nm + ip]
for ip in NonBuildingL:
    y=y+dS1[nm + ip]

flg=0
if flg==1:
    plt.close('all')
    plt.plot(tv,y,'-')

#%% Save

d={}
d['Year']=tv
d['Default']=dDef
d['SSP 1']=dS1

gu.opickle(r'C:\Users\rhember\Documents\Code_Python\fcgadgets\cbrunner\Parameters\Variables_HWP_EndUses.pkl',d)



