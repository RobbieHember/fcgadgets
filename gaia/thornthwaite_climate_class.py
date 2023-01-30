'''
Thornthwaite climate classification
'''

#%% Import modules

import numpy as np
import matplotlib.pyplot as plt
import fcgadgets.macgyver.utilities_general as gu
import fcgadgets.macgyver.utilities_gis as gis

#%% Calculate variables

zRef=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Admin\BC_Land_Mask.tif')

pet=np.zeros(zRef['Data'].shape,dtype='float32')
prcp=np.zeros(zRef['Data'].shape,dtype='float32')
defi=np.zeros(zRef['Data'].shape,dtype='float32')
sur=np.zeros(zRef['Data'].shape,dtype='float32')

for mo in range(12):

  # Version 2:
  pet0=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Climate\Monthly\BC1ha_etp_tmw_norm_1971to2000_comp_hist_v1\BC1ha_etp_mon_norm_1971to2000_comp_hist_v1_' + str(mo+1) + '.tif')
  pet0=pet0['Data'].astype('float32')/100*30
  pet=pet+pet0

  prcp0=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Climate\Monthly\BC1ha_prcp_mon_norm_1971to2000_si_hist_v1\BC1ha_prcp_mon_norm_1971to2000_si_hist_v1_' + str(mo+1) + '.tif')
  prcp0=prcp0['Data'].astype('float32')
  prcp=prcp+prcp0

  ind=np.where(pet0>prcp0)
  defi[ind]=defi[ind]+pet0[ind]-prcp0[ind]

  ind=np.where(pet0<prcp0)
  sur[ind]=sur[ind]+prcp0[ind]-pet0[ind]

#%% Classification

z0=((sur-defi)/pet)*100

lab=['Hyper humid','Humid','Subhumid wet','Subhumid dry','Semi arid','Arid','Outside Boundary']

z1=7*np.ones(z0.shape,dtype='int8')
ind=np.where(z0>=100); z1[ind]=1;
ind=np.where( (z0>=20) & (z0<100) ); z1[ind]=2;
ind=np.where( (z0>=0) & (z0<20) ); z1[ind]=3;
ind=np.where( (z0>=-20) & (z0<0) ); z1[ind]=4;
ind=np.where( (z0>=-40) & (z0<-20) ); z1[ind]=5;
ind=np.where(z0<-40); z1[ind]=6;
ind=np.where(zRef['Data']==0); z1[ind]=7;

# plt.matshow(z1,clim=[1,7])

# Save
z2=zRef.copy()
z2['Data']=z1
gis.SaveGeoTiff(z2,r'C:\Users\rhember\Documents\Data\BC1ha\Climate\BC1ha_ThornthwaiteClimateClass_norm_1971to2000.tif')

#%% Classification (condensed)

z0=((sur-defi)/pet)*100

lab=['Hyper humid','Humid','Subhumid','Semi arid','Arid','Outside Boundary']

z1=6*np.ones(z0.shape,dtype='int8')
ind=np.where(z0>=100); z1[ind]=1;
ind=np.where( (z0>=20) & (z0<100) ); z1[ind]=2;
ind=np.where( (z0>=-20) & (z0<20) ); z1[ind]=3;
ind=np.where( (z0>=-40) & (z0<-20) ); z1[ind]=4;
ind=np.where(z0<-40); z1[ind]=5;
ind=np.where(zRef['Data']==0); z1[ind]=6;

# plt.matshow(z1,clim=[1,7])

# Save
z2=zRef.copy()
z2['Data']=z1
gis.SaveGeoTiff(z2,r'C:\Users\rhember\Documents\Data\BC1ha\Climate\BC1ha_ThornthwaiteClimateClassCondensed_norm_1971to2000.tif')