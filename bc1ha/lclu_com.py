#%% Import modules
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import cv2
import fcgadgets.macgyver.util_general as gu
import fcgadgets.macgyver.util_gis as gis
import fcgadgets.bc1ha.bc1ha_utils as u1ha
import fcgadgets.gaia.lclu_util as ulclu
#import fcgadgets.gaia.lclu_plot as plclu
gp=gu.SetGraphics('Manuscript')

#%% Import paths and look-up-tables
meta=u1ha.Init()
zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])

#%% Import CEC map 2020
ulclu.PrepareLandUseCEC(meta)

#%% Import NTEMS (age, land cover, harvest year)
ulclu.PrepareNTEMS(meta)

#%% Reclassify NTEMS and VRI according to LC Comp 1
ulclu.ReclassifyNTEMS_LandCover(meta)
ulclu.ReclassifyVRI_LandCover(meta)

#%% Derive land cover compilation 1 (year 2019)
ulclu.DeriveLandCoverComp1(meta)

#%% Derive land use compilation 1 (year 2019)
ulclu.DeriveLandUseComp1(meta)

#%% Derive land cover compilation 1 (year 1800)
ulclu.DeriveLandCoverComp1_1800(meta)

#%% Export summary of LC/LU areas (ha)
ulclu.ExportSummaryLCLU(meta)

#%% Derive land use change year for 1800-2019
ulclu.DeriveLandUseChangeYear_1800to2019(meta)

#%% Derive land cover / use compilation 1 (year 2029 scenarios)
ulclu.DeriveLandCoverLandUseComp1_2020to2049_Scenarios(meta)

#%% Derive deforestation mask
ulclu.DeriveLandCoverLandUseComp1_DeforstationMask(meta)

#%% Rasterize early land use change year
ulclu.RasterizeEaryLandUseChangeYear(meta)

#%% Global Forest Change Loss Year (from GEE)
ulclu.Download_GFC(meta)

#%% Filter Global Forest Change Loss Year (to remove known disturbances)
ulclu.FilterGFC_LossYear(meta)

#%% Compare forest area estimates
ulclu.CompareForestArea(meta)

#%% Global forest change analysis

#%%
def AnalyzeGlobalForestChange(meta):

	meta=u1ha.Init()
	zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])

	lut_comp1=meta['LUT']['Derived']['lc_comp1']

	# Import data
	z=u1ha.Import_Raster(meta,[],['lc_comp1_2019','gfc_tc2000','gfc_ly','fire_yl','ibm_yr','harv_yr_comp1','road_ften_s'],'Extract Grid')

	z['gfc_tc2000']=z['gfc_tc2000']*zRef['Data']
	z['gfc_ly']=z['gfc_ly']*zRef['Data']

	kernel=np.ones((1,1))
	z['harv_yr_comp1']=cv2.dilate(z['harv_yr_comp1'].astype(np.uint8),kernel,iterations=1)
	z['road_ften_s']=cv2.dilate(z['road_ften_s'].astype(np.uint8),kernel,iterations=1)

	tv=np.arange(2001,2023,1)
	A={}
	A['Total']=np.zeros(tv.size)
	A['Wildfire']=np.zeros(tv.size)
	A['Harvest']=np.zeros(tv.size)
	A['Road']=np.zeros(tv.size)
	A['IBM']=np.zeros(tv.size)
	for i,yr in enumerate(tv):
		print(yr)
		ind=np.where( (z['gfc_ly']==yr-2000) ); A['Total'][i]=ind[0].size/1e6
		ind=np.where( (z['gfc_ly']==yr-2000) & (z['fire_yl']>=yr-2) & (z['fire_yl']<=yr) ); A['Wildfire'][i]=ind[0].size/1e6
		ind=np.where( (z['gfc_ly']==yr-2000) & (np.abs(z['harv_yr_comp1']-yr)<=5) ); A['Harvest'][i]=ind[0].size/1e6
		ind=np.where( (z['gfc_ly']==yr-2000) & (z['road_ften_s']>0) ); A['Road'][i]=ind[0].size/1e6
		ind=np.where( (z['gfc_ly']==yr-2000) & (z['ibm_yr']>=2000) & (z['fire_yl']<2000) & (z['harv_yr_comp1']<2000) & (z['road_ften_s']==0) ); A['IBM'][i]=ind[0].size/1e6

	A['Undesignated']=A['Total']-A['Wildfire']-A['Harvest']-A['Road']-A['IBM']

	plt.close('all')
	plt.plot(tv,A['Total'],'-ko')
	plt.plot(tv,A['Wildfire'],'-rs')
	plt.plot(tv,A['Harvest'],'-b^')
	plt.plot(tv,A['Road'],'-cs')
	plt.plot(tv,A['IBM'],'-gd')
	plt.plot(tv,A['Undesignated'],'-yo')

	Am={}
	for k in A.keys():
		Am[k]=np.mean(A[k][(tv>2000) & (tv<2023)])
	Am

	yr=2012
	y=np.zeros(zRef['Data'].shape,dtype='int8')
	ind=np.where( (z['gfc_ly']==yr-2000) ); y[ind]=2
	ind=np.where( (y==2) & (z['fire_yl']==yr) ); y[ind]=1
	ind=np.where( (y==2) & (np.abs(z['harv_yr_comp1']-yr)<=5) ); y[ind]=1
	ind=np.where( (y==2) & (z['road_ften_s']>0) ); y[ind]=1
	ind=np.where( (y==2) & (z['ibm_yr']>=2000) ); y[ind]=1

	plt.matshow(y)

	return

#%% OLD ANALYSIS BELOW

#%% Import data

zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])

vList=['rd','lcc_cec_2010','lcc_cec_2020','harv_yr_con1','fire_yr','ibm_yr'] #'lcc1_c','gfcly','gfcly_filt',
z=u1ha.Import_Raster(meta,[],vList)

#%% Land use change analysis - BC wide

id=np.array(list(meta['LUT']['Derived']['lcc_cec_c'].values()))
# Areas at start and end
v=np.zeros( (2,id.size) )
for i in range(id.size):
    ind=np.where( (z['lcc_cec_2010']['Data']==id[i]) & (z['harv_yr_con1']['Data']==0) & (z['fire_yr']['Data']==0) & (z['ibm_yr']['Data']==0) )
    v[0,i]=ind[0].size
    ind=np.where( (z['lcc_cec_2020']['Data']==id[i]) & (z['harv_yr_con1']['Data']==0) & (z['fire_yr']['Data']==0) & (z['ibm_yr']['Data']==0) )
    v[1,i]=ind[0].size
df=pd.DataFrame(v)
df.to_excel(r'C:\Users\rhember\OneDrive - Government of BC\Manuscripts\Land Use Change Assessment\CEC Change Analysis\BC Wide\cec_c0_nodist.xlsx')

# Transitions
v=np.zeros( (id.size,id.size) )
for i in range(id.size):
    for j in range(id.size):
        ind=np.where( (z['lcc_cec_2010']['Data']==id[i]) & (z['lcc_cec_2020']['Data']==id[j]) & (z['harv_yr_con1']['Data']==0) & (z['fire_yr']['Data']==0) & (z['ibm_yr']['Data']==0) )
        v[i,j]=ind[0].size
df=pd.DataFrame(v)
df.to_excel(r'C:\Users\rhember\OneDrive - Government of BC\Manuscripts\Land Use Change Assessment\CEC Change Analysis\BC Wide\cec_c1_nodist.xlsx')

#%% Filter by regional district

ind=np.where(z['rd']['Data']==meta['LUT']['TA_REGIONAL_DISTRICTS_SVW']['REGIONAL_DISTRICT_NAME']['CAPITAL'])
for k in z.keys():
    z[k]['Data']=z[k]['Data'][ind]

#%% Land use change analysis - by regional district (based on compressed categories)

id=np.array(list(meta['LUT']['Derived']['lcc_cec_c'].values()))

rs={}
for k in meta['LUT']['TA_REGIONAL_DISTRICTS_SVW']['REGIONAL_DISTRICT_NAME'].keys():
    # Areas at start and end
    vI=np.zeros( (2,id.size) )
    for i in range(id.size):
        ind=np.where( (z['lcc_cec_2010']['Data']==id[i]) & (z['harv_yr_con1']['Data']==0) & (z['fire_yr']['Data']==0) & (z['ibm_yr']['Data']==0) & (z['rd']['Data']==meta['LUT']['TA_REGIONAL_DISTRICTS_SVW']['REGIONAL_DISTRICT_NAME'][k]) )
        vI[0,i]=ind[0].size
        ind=np.where( (z['lcc_cec_2020']['Data']==id[i]) & (z['harv_yr_con1']['Data']==0) & (z['fire_yr']['Data']==0) & (z['ibm_yr']['Data']==0) & (z['rd']['Data']==meta['LUT']['TA_REGIONAL_DISTRICTS_SVW']['REGIONAL_DISTRICT_NAME'][k]) )
        vI[1,i]=ind[0].size
    df=pd.DataFrame(vI)
    df.to_excel(r'C:\Users\rhember\OneDrive - Government of BC\Manuscripts\Land Use Change Assessment\CEC Change Analysis\cec_c0_nodist_' + k + '.xlsx')

    # Transitions
    vT=np.zeros( (id.size,id.size) )
    for i in range(id.size):
        for j in range(id.size):
            ind=np.where( (z['lcc_cec_2010']['Data']==id[i]) & (z['lcc_cec_2020']['Data']==id[j]) & (z['harv_yr_con1']['Data']==0) & (z['fire_yr']['Data']==0) & (z['ibm_yr']['Data']==0) & (z['rd']['Data']==meta['LUT']['TA_REGIONAL_DISTRICTS_SVW']['REGIONAL_DISTRICT_NAME'][k]) )
            vT[i,j]=ind[0].size
    df=pd.DataFrame(vT)
    df.to_excel(r'C:\Users\rhember\OneDrive - Government of BC\Manuscripts\Land Use Change Assessment\CEC Change Analysis\cec_c1_nodist_' + k + '.xlsx')

    # Net deforestation summary
    ind=np.array([meta['LUT']['Derived']['lcc_cec_c']['Cropland'],meta['LUT']['Derived']['lcc_cec_c']['Barren Ground'],meta['LUT']['Derived']['lcc_cec_c']['Urban']],dtype='int8')
    Ad=np.sum(vT[0,ind-1])
    Aa=np.sum(vT[ind-1,0])
    Pd=np.sum(vT[0,ind-1])/vI[0,0]*100
    Pa=np.sum(vT[ind-1,0])/vI[0,0]*100
    Pn=Pa-Pd

#%% Map deforestation
# Notes:
#   1) Forest to Grassland and Forest to Shrubland excluded due to suspicion of error
#   2) Areas with harvest excluded (harvest from NTEM)

z1=zRef.copy()
z1['Data']=np.zeros(zRef['Data'].shape,dtype='int8')
ind=np.where( (z['lcc_cec_2010']['Data']==meta['LUT']['Derived']['lcc_cec_c']['Forest']) & (z['lcc_cec_2020']['Data']==meta['LUT']['Derived']['lcc_cec_c']['Urban']) & (z['harv_yr_con1']['Data']==0) )
z1['Data'][ind]=1
ind=np.where( (z['lcc_cec_2010']['Data']==meta['LUT']['Derived']['lcc_cec_c']['Forest']) & (z['lcc_cec_2020']['Data']==meta['LUT']['Derived']['lcc_cec_c']['Cropland']) & (z['harv_yr_con1']['Data']==0) )
z1['Data'][ind]=2
ind=np.where( (z['lcc_cec_2010']['Data']==meta['LUT']['Derived']['lcc_cec_c']['Forest']) & (z['lcc_cec_2020']['Data']==meta['LUT']['Derived']['lcc_cec_c']['Barren Ground']) & (z['harv_yr_con1']['Data']==0) )
z1['Data'][ind]=3
#plt.matshow(z1['Data'])
gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\Disturbances\\LandUseChange_Deforestation_10to20_CEC.tif')

#%% Map Afforestation (from cropland, urban, or barren ground)

z1=zRef.copy()
z1['Data']=np.zeros(zRef['Data'].shape,dtype='int8')
ind=np.where( (z['lcc_cec_2010']['Data']==meta['LUT']['Derived']['lcc_cec_c']['Urban']) & (z['lcc_cec_2020']['Data']==meta['LUT']['Derived']['lcc_cec_c']['Forest']) & (z['harv_yr_con1']['Data']==0) )
z1['Data'][ind]=1
ind=np.where( (z['lcc_cec_2010']['Data']==meta['LUT']['Derived']['lcc_cec_c']['Cropland']) & (z['lcc_cec_2020']['Data']==meta['LUT']['Derived']['lcc_cec_c']['Forest']) & (z['harv_yr_con1']['Data']==0) )
z1['Data'][ind]=2
ind=np.where( (z['lcc_cec_2010']['Data']==meta['LUT']['Derived']['lcc_cec_c']['Barren Ground']) & (z['lcc_cec_2020']['Data']==meta['LUT']['Derived']['lcc_cec_c']['Forest']) & (z['harv_yr_con1']['Data']==0) )
z1['Data'][ind]=3
gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\Disturbances\\LandUseChange_Afforestation_10to20_CEC.tif')

# ind=np.where( (z['harv_yr_con1']['Data']>0) )
# z['Data'][ind]=2
# plt.matshow(z['Data'])

# ind=np.where(z['Data']==1)
# ind[1].size/1000000

#%% Map Forest to Grassland and Forest to Shrubland
# *** Way too much to be real! ***

z=zRef.copy()
z['Data']=np.zeros(zRef['Data'].shape,dtype='int8')
ind=np.where( (z['lcc_cec_2010']['Data']==meta['LUT']['Forest']) & (z['lcc_cec_2020']['Data']==meta['LUT']['Grassland']) )
z['Data'][ind]=1
ind=np.where( (z['lcc_cec_2010']['Data']==meta['LUT']['Forest']) & (z['lcc_cec_2020']['Data']==meta['LUT']['Shrubland']) )
z['Data'][ind]=1
ind=np.where( (z['harv_yr_con1']['Data']>0) )
z['Data'][ind]=2
plt.matshow(z['Data'])







ind=np.where(z['harv_yr_con1']['Data']>0)
ind[1].size/1000000

#%%

tv=np.arange(2001,2022,1)
#List=[meta['LUT']['Derived']['lcc_cec_c']['Urban'],meta['LUT']['Derived']['lcc_cec_c']['Barren Ground'],meta['LUT']['Derived']['lcc_cec_c']['Cropland'],meta['LUT']['Derived']['lcc_cec_c']['Grassland']]
List=[meta['LUT']['Derived']['lcc_cec_c']['Urban'],meta['LUT']['Derived']['lcc_cec_c']['Cropland']]
A1=np.zeros(tv.size)
A2=np.zeros(tv.size)
A3=np.zeros(tv.size)
for iT in range(tv.size):
    print(tv[iT])
    ind=np.where( (z['gfcly']['Data']==tv[iT]) )[0]
    A1[iT]=ind.size
    ind=np.where( (z['gfcly_filt']['Data']==tv[iT]) )[0]
    A2[iT]=ind.size
    #ind=np.where( (z['gfcly']['Data']==tv[iT]) & (np.abs(z['gfcly']['Data']-z['harv_yr_con1']['Data'])>5) & (z['lcc1_c']['Data']==meta['LUT']['Derived']['lcc1']['Forest']) )[0]
    ind=np.where( (z['gfcly']['Data']==tv[iT]) & (np.isin(z['lcc_cec_2020']['Data'],List)==True) )[0]
    A3[iT]=ind.size

plt.close('all')
plt.plot(tv,A1,'-bo')
plt.plot(tv,A2,'-gs')
plt.plot(tv,A3,'-cd')


