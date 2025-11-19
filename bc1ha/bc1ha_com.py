#%% Import modules
import os
import numpy as np
import time
import matplotlib.pyplot as plt
import pandas as pd
import fiona
import fcgadgets.macgyver.util_general as gu
import fcgadgets.macgyver.util_gis as gis
import fcgadgets.bc1ha.bc1ha_utils as u1ha
import fcgadgets.bc1ha.bc1ha_update as up1ha
import fcgadgets.bc1ha.bc1ha_climate as c1ha
gp=gu.SetGraphics('Manuscript')

#%% Import paths and look-up-tables
meta=u1ha.Init()
zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])

#%%
#d=gu.ReadExcel(r'G:\My Drive\Code_Python\fcgadgets\cbrunner\Parameters\Workbook_BCForestSector.xlsx')
#plt.plot(d['Year'],d['Production Electricity Total BC (PJ/yr)'],'-bo')

#%% Build look up tables (only do this once a year, takes 1.5 hours)
up1ha.BuildLUTsFromSourceGDBs(meta)

#%% Major rivers
u1ha.GetRiversMajor(meta)

#%% Simplify geometry of provincial geodatabases (faster useage)
u1ha.SimplifyProvincialGDBs(meta)

#%% Create land mask for BC
u1ha.GenerateLandMaskBC(meta)

#%% Digitize the boundary of TSAs
# The original is organized at the sub-TSA level so need to aggregate and then digitize
up1ha.DigitizeTSABoundaries(meta)

#%% Rasterize VRI (manually rasterize vri_feaid.tif in ArcGIS first)
up1ha.RasterizeVRI(meta,2024)

#%% Rasterize variables from source
# *** VRI and Forest Cover Inventory are too big - crash ***
u1ha.RasterizeFromSource(meta,zRef,'FADM_TSA','TSA_NUMBER_DESCRIPTION')
u1ha.RasterizeFromSource(meta,zRef,'FADM_TFL','FOREST_FILE_ID')
u1ha.RasterizeFromSource(meta,zRef,'BEC_BIOGEOCLIMATIC_POLY','ZONE')
u1ha.RasterizeFromSource(meta,zRef,'BEC_BIOGEOCLIMATIC_POLY','SUBZONE')
u1ha.RasterizeFromSource(meta,zRef,'BEC_NATURAL_DISTURBANCE_SV','NATURAL_DISTURBANCE_TYPE_CODE')
u1ha.RasterizeFromSource(meta,zRef,'BC_MAJOR_CITIES_POINTS_500M','NAME')
u1ha.RasterizeFromSource(meta,zRef,'F_OWN','OWNERSHIP_DESCRIPTION')
u1ha.RasterizeFromSource(meta,zRef,'FTEN_ROAD_SEGMENT_LINES_SVW','LIFE_CYCLE_STATUS_CODE')
u1ha.RasterizeFromSource(meta,zRef,'FTEN_ROAD_SEGMENT_LINES_SVW','RETIREMENT_DATE')
u1ha.RasterizeFromSource(meta,zRef,'FTEN_CUT_BLOCK_POLY_SVW','TIMBER_MARK')
u1ha.RasterizeFromSource(meta,zRef,'FTEN_CUT_BLOCK_POLY_SVW','CUT_BLOCK_ID')
u1ha.RasterizeFromSource(meta,zRef,'FWA_WETLANDS_POLY','WATERBODY_TYPE')
u1ha.RasterizeFromSource(meta,zRef,'FWA_MANMADE_WATERBODIES_POLY','WATERBODY_TYPE')
u1ha.RasterizeFromSource(meta,zRef,'DRA_DGTL_ROAD_ATLAS_MPAR_SP','DIGITAL_ROAD_ATLAS_LINE_ID')
u1ha.RasterizeFromSource(meta,zRef,'GBA_TRANSMISSION_LINES_SP','TRANSMISSION_LINE_ID')
u1ha.RasterizeFromSource(meta,zRef,'HSP_MJR_MINES_PERMTTD_AREAS_SP','STATUS_TYPE')
u1ha.RasterizeFromSource(meta,zRef,'NRC_POPULATED_PLACES_1M_SP','NAME')
u1ha.RasterizeFromSource(meta,zRef,'OATS_ALR_POLYS','ALR_POLY_ID')
u1ha.RasterizeFromSource(meta,zRef,'OGSR_TAP_PRIORITY_DEF_AREA_SP','OGSR_TPDA_SYSID')
u1ha.RasterizeFromSource(meta,zRef,'OGSR_TAP_PRIORITY_DEF_AREA_SP','PRIORITY_DEFERRAL_ID')
u1ha.RasterizeFromSource(meta,zRef,'RMP_OGMA_LEGAL_ALL_SVW','LEGAL_OGMA_PROVID')
u1ha.RasterizeFromSource(meta,zRef,'RSLT_OPENING_SVW','OPENING_ID')
u1ha.RasterizeFromSource(meta,zRef,'RSLT_ACTIVITY_TREATMENT_SVW','OPENING_ID')
u1ha.RasterizeFromSource(meta,zRef,'RSLT_ACTIVITY_TREATMENT_SVW','SILV_BASE_CODE')
u1ha.RasterizeFromSource(meta,zRef,'RSLT_ACTIVITY_TREATMENT_SVW','SILV_TECHNIQUE_CODE')
u1ha.RasterizeFromSource(meta,zRef,'RSLT_FOREST_COVER_RESERVE_SVW','SILV_RESERVE_CODE')
u1ha.RasterizeFromSource(meta,zRef,'TA_MUNICIPALITIES_SVW','MUNICIPALITY_NAME')
u1ha.RasterizeFromSource(meta,zRef,'TA_REGIONAL_DISTRICTS_SVW','REGIONAL_DISTRICT_NAME')
u1ha.RasterizeFromSource(meta,zRef,'TA_PROTECTED_LANDS_SV','PROTECTED_LANDS_DESIGNATION')
u1ha.RasterizeFromSource(meta,zRef,'VEG_BURN_SEVERITY_SP','FIRE_YEAR')
u1ha.RasterizeFromSource(meta,zRef,'VEG_BURN_SEVERITY_SP','BURN_SEVERITY_RATING')

#%% Generate sparse inputs
# This speeds import by an order of magnitude
u1ha.GenerateSparseInputs(meta,200,'Province')
u1ha.GenerateSparseInputs(meta,100,'Province')
u1ha.GenerateSparseInputs(meta,50,'Province')
u1ha.GenerateSparseInputs(meta,25,'Province')
u1ha.GenerateSparseInputs(meta,20,'Province')
u1ha.GenerateSparseInputs(meta,50,'BCFCS_LUC')
u1ha.GenerateSparseInputs(meta,10,'BCFCS_NOSE')
u1ha.GenerateSparseInputs(meta,20,'BCFCS_NOSE')
u1ha.GenerateSparseInputs(meta,10,'BCFCS_NMC')
u1ha.GenerateSparseInputs(meta,20,'BCFCS_NMC')
u1ha.GenerateSparseInputs(meta,1,'TSA_DawsonCreek')
u1ha.GenerateSparseInputs(meta,1,'Landscape_NicolaRiverWatershed')
u1ha.GenerateSparseInputs(meta,50,'BCFCS_CWH')
u1ha.GenerateSparseInputs(meta,20,'BCFCS_CWH')
u1ha.GenerateSparseInputs(meta,10,'BCFCS_CWH')
u1ha.GenerateSparseInputs(meta,50,'BCFCS_SBS')
u1ha.GenerateSparseInputs(meta,1,'BCFCS_Eval')
u1ha.GenerateSparseInputs(meta,1,'BCFCS_EvalAtPlots')
u1ha.GenerateSparseInputs(meta,1,'BCFCS_EvalAtCN')
u1ha.GenerateSparseInputs(meta,1,'BCFCS_EvalCoast')

#%% Gap-fill BGC Zone
u1ha.GapFillBGCZ(meta)

#%% Ecozones of Canada
u1ha.RasterizeEcozonesOfCanada(meta)

#%% Prepare CEC Land Use Map 2020
u1ha.PrepareLandUseCEC(meta)

#%% Prepare data from NTEMS (age, land cover, harvest year)
up1ha.PrepareNTEMS(meta)

#%% Derive land cover compilation 1 (year 2019)
u1ha.DeriveLandCoverComp1(meta)

#%% Derive land use compilation 1 (year 2019)
u1ha.DeriveLandUseComp1(meta)

#%% Derive land cover compilation 1 (year 1800)
u1ha.DeriveLandCoverComp1_1800(meta)

#%% Export summary of LC/LU areas (ha)
u1ha.ExportSummaryLCLU(meta)

#%% Derive land use change year for 1800-2019
u1ha.DeriveLandUseChangeYear_1800to2019(meta)

#%% Derive land cover / use compilation 1 (year 2029 scenarios)
u1ha.DeriveLandCoverLandUseComp1_2020to2049_Scenarios(meta)

#%% Derive deforestation mask
u1ha.DeriveLandCoverLandUseComp1_DeforstationMask(meta)

#%% Reclassify NTEMS and VRI according to LC Comp 1
up1ha.ReclassifyNTEMS_LandCover(meta)
u1ha.ReclassifyVRI_LandCover(meta)

#%% Upland-wetland Forest
u1ha.DeriveUplandWetlandForstMask(meta)
#z1=u1ha.Import_Raster(meta,[],['upwetf_ntems'])#gu.CountByCategories(z['upwetf_ntems']['Data'],'Percent')

#%% Rasterize OPENING ID FROM OPENING LAYER (1 hour)
up1ha.RasterizeOpeningID2(meta)

#%% Rasterize forest cover attributes
# Use FC ID created in ArcGIS to connect with attributes
up1ha.RasterizeForestCoverInventory(meta)

#%% Rasterize denudations in the ATU layer
up1ha.RasterizeDenudationsFromATU(meta)

#%% Wildfire perimiters
YearCurrent=2025
YearLast=2024
up1ha.RasterizeWildfirePerimitersHistorical(meta,zRef,YearLast) # Rasterize hisorical wildfire
#up1ha.RasterizeWildfirePerimitersCurrentYear(meta,YearCurrent) # Current year + preveous year (sometimes missing) - Need to add this to historical and then re-create the condensed layers
up1ha.DeriveWildfireComposite(meta,zRef,YearCurrent) # Combine information sources

#%% Burn severity rating
# Historical is rasterized from source
up1ha.RasterizeBurnSeverityCurrentYear(meta,2024) # Rasterize current-year burn severity rating
#up1ha.RasterizeBurnSeverityCanada(meta) # Reproject National Burn Severity data
u1ha.DeriveBurnSeverityCompilation(meta) # Derive burn severity to pair with wildfire perimeters

#%% Rasterize insect compilation from AOS (5 hours)
YearLast=2024
up1ha.RasterizeInsects(meta,zRef,YearLast) # Rasterize insects from AOS (Use InsectComp1 where possible, but this is still needed)
up1ha.RasterizeInsectComp1(meta,YearLast)
u1ha.CalcInsectComp1_TimeSeries(meta)
u1ha.SummarizeInsectComp1_ByTime(meta)

#%%
zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
zG=u1ha.Import_Raster(meta,[],['bgcz','spc1_vri23','fire_yl'],'Extract Grid')

ind1=np.where( (zG['bgcz']==meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE']['IDF']) & (np.isin(zG['spc1_vri23'],[meta['LUT']['VEG_COMP_LYR_R1_POLY']['SPECIES_CD_1']['FD'],meta['LUT']['VEG_COMP_LYR_R1_POLY']['SPECIES_CD_1']['FDI']])==True) )[0]
N1=ind1.size/1e6

ind2=np.where( (zG['bgcz']==meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE']['IDF']) & (np.isin(zG['spc1_vri23'],[meta['LUT']['VEG_COMP_LYR_R1_POLY']['SPECIES_CD_1']['FD'],meta['LUT']['VEG_COMP_LYR_R1_POLY']['SPECIES_CD_1']['FDI']])==True) & (zG['fire_yl']<1920) )[0]
N2=ind2.size/1e6

ind3=np.where( (zG['bgcz']==meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE']['IDF']) & (np.isin(zG['spc1_vri23'],[meta['LUT']['VEG_COMP_LYR_R1_POLY']['SPECIES_CD_1']['FD'],meta['LUT']['VEG_COMP_LYR_R1_POLY']['SPECIES_CD_1']['FDI']])==True) & (zG['fire_yl']>1920) )[0]
N3=ind3.size/1e6
print([N1,N2,N3])

N_yr=2024-1920
y1=0.84/2.33
y2=y1/N_yr

#%%

def SummarizeInsectComp1_ByTime(meta):

	zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
	zG=u1ha.Import_Raster(meta,[],['bgcz','spc1_vri02','spc2_vri02','spc3_vri02','spc1_pct_vri02','spc2_pct_vri02','spc3_pct_vri02'],'Extract Grid')

	pine_frac=np.zeros(zRef['Data'].shape)
	pineL=[meta['LUT']['VEG_COMP_LYR_R1_POLY']['SPECIES_CD_1']['PL'],meta['LUT']['VEG_COMP_LYR_R1_POLY']['SPECIES_CD_1']['PLI'],meta['LUT']['VEG_COMP_LYR_R1_POLY']['SPECIES_CD_1']['PY']]
	for i in range(1,4):
		ind=np.where( (np.isin(zG['spc' + str(i) + '_vri02'],pineL)==True) & (zG['spc' + str(i) + '_pct_vri02']>0) )
		pine_frac[ind]=zG['spc' + str(i) + '_pct_vri02'][ind].astype('float')/100

	idP=meta['Param']['Raw']['InsectComp1']['ID']
	bugP=meta['Param']['Raw']['InsectComp1']['Insect Name']
	sevP=meta['Param']['Raw']['InsectComp1']['Severity Code']

	sevL=np.array(['L','M','S','V'])
	mL=np.array([1,5,20,40,75])

	m_cumu=np.zeros(zG['bgcz'].shape)

	bugL=np.unique(bugP)

	# Initialize dictionary
	YearLast=2023
	d={}
	d['Year']=np.arange(1950,YearLast+1,1)
	d['Data']={}
	for bug in bugL:
		d['Data'][bug]={}
		for sev in sevL:
			d['Data'][bug][sev]={}
			d['Data'][bug][sev]['Province']=np.zeros(d['Year'].size)
			d['Data'][bug][sev]['ByBGC']={}
			for zone in meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE'].keys():
				d['Data'][bug][sev]['ByBGC'][zone]=np.zeros(d['Year'].size)

	N_Year=10
	for iY in range(N_Year):
		print(iY)
		yr0=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\PEST_INFESTATION_POLY\\InsectComp1_' + str(iY+1) + '_Year.tif')['Data']
		bug0=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\PEST_INFESTATION_POLY\\InsectComp1_' + str(iY+1) + '_Type.tif')['Data']

		ind0=np.where(yr0>0)
		yr0=yr0[ind0]
		bug0=bug0[ind0]
		zone0=zG['bgcz'][ind0]
		pine_frac0=pine_frac[ind0]

		uID=np.unique(bug0)
		for uBug in uID:
			iID=np.where(idP==uBug)[0]
			bugName=bugP[iID[0]]
			iSev=np.where( (sevL==sevP[iID[0]]) )[0]
			sevName=sevL[iSev[0]]

			if bugName=='Mountain Pine Beetle':
				ind1=np.where( (bug0==uBug) & (pine_frac0>0.25) )[0]
			else:
				ind1=np.where( (bug0==uBug) )[0]
			yr1=yr0[ind1]
			idx=gu.IndicesFromUniqueArrayValues(yr1)
			for i in idx.keys():
				iT=np.where(d['Year']==i)[0]
				d['Data'][ bugName ][ sevName ]['Province'][iT]=d['Data'][ bugName ][ sevName ]['Province'][iT]+idx[i].size

			for zone in meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE'].keys():
				if bugName=='Mountain Pine Beetle':
					ind1=np.where( (bug0==uBug) & (pine_frac0>0.25) & (zone0==meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE'][zone]) )[0]
				else:
					ind1=np.where( (bug0==uBug) & (zone0==meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE'][zone]) )[0]
				if ind1.size==0:
					continue
				yr1=yr0[ind1]
				idx=gu.IndicesFromUniqueArrayValues(yr1)
				for i in idx.keys():
					iT=np.where(d['Year']==i)[0]
					d['Data'][ bugName ][ sevName ]['ByBGC'][zone][iT]=d['Data'][ bugName ][ sevName ]['ByBGC'][zone][iT]+idx[i].size

	# Check that it workded
	flg=0
	if flg==1:
		#plt.close('all'); plt.plot(d['Year'],d['Data']['Mountain Pine Beetle']['S']['Province'],'-ko')
		plt.close('all'); plt.plot(d['Year'],d['Data']['Mountain Pine Beetle']['S']['ByBGC']['SBS'],'-bo')
		plt.close('all'); plt.plot(d['Year'],d['Data']['Mountain Pine Beetle']['S']['ByBGC']['IDF'],'-bo')
		#plt.close('all'); plt.plot(d['Year'],d['Data']['Douglas-fir Beetle']['S']['ByBGC']['IDF'],'-go')

	# Save
	gu.opickle(meta['Paths']['DB']['Beetles'] + '\\InsectComp1_SummaryByTime.pkl',d)

	return

#%% Rasterize harvesting
u1ha.RasterizeHarvest_CC(meta,zRef) # Rasterize consolidated cutblocks database
u1ha.DeriveHarvestEarlyYears(meta) # Harvest early reconstruction
u1ha.DeriveHarvestCompilation(meta) # Harvest consolidated (with and w/o early reconstruction)
u1ha.DeriveHarvest_SILV_SYSTEM_CODE(meta,zRef) # Harvest silviculture system code

#%% Tabular summary of activities - full query of RESULTS (list of everything), exports to spreadsheet (takes 10-15 min)
YearLast=2024
up1ha.Query_RESULTS(meta,YearLast)

#%%
pthin=meta['Paths']['GDB']['Results']
fiona.listlayers(pthin)
with fiona.open(pthin,layer='RSLT_OPENING_SVW') as source:
	for feat in source:
		#geom=dict(feat['geometry'].ite:ms())
		prop=dict(feat['properties'].items())
		if prop['OPENING_ID']==1707913:
			print(prop['OPENING_ID'])
			print(prop)
			break

#%% Extract geojson files for each activity (non-obligation only, no surveys, no audits, no denudations)
up1ha.VectorMapsByActivity_NO(meta)

#%% Rasterize knockdown
YearLast=2024
u1ha.RasterizeSilviculture(meta,np.array(['SP']),np.array([]),np.array(['CABLE','GUARD','HARV','MDOWN','PUSH']),np.array([]),'SP-KD',YearLast)
u1ha.Calc_AIL_Silv_ByFSC(meta,'SP-KD',YearLast)

#%% Rasterize road deactivation
# *** The spatial info for this appears to be nonsense ***
u1ha.RasterizeSilviculture(meta,np.array(['RD']),np.array(['DE']),np.array([]),np.array([]),'RD-DE',YearLast)
u1ha.Calc_AIL_Silv_ByFSC(meta,'RD-DE',YearLast)
u1ha.Calc_AIL_Silv_ByBGC(meta,'RD-DE',YearLast)

# *** The spatial info for this appears to be nonsense ***
# d=gu.ipickle(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\RD-DE_SummaryByTimeAndFSC.pkl')
#zRD=np.zeros(zRef['Data'].shape,dtype='int16')
#for i in range(4):
#	zYr=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\RD-DE_' + str(i+1) + '_Year.tif')['Data']
#	ind=np.where(zYr>0)
#	zRD[ind]=zYr[ind]
#plt.matshow(zRD)

#%% Rasterize direct seeding
vNam='DS'
YearLast=2024
u1ha.RasterizeSilviculture(meta,np.array(['DS']),np.array([]),np.array([]),np.array([]),vNam,YearLast)
u1ha.Calc_AIL_Silv_ByFSC(meta,vNam,YearLast)
u1ha.Calc_AIL_Silv_ByBGC(meta,vNam,YearLast)
#d=gu.ipickle(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + vNam + '_SummaryByTimeAndBGCZone.pkl')
d=gu.ipickle(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + vNam + '_SummaryByTimeAndFSC.pkl')

#%% Rasterize planting (year, funding source, composition, density, genetic worth)
vNam='PL_All'
YearLast=2024
up1ha.RasterizePlanting(meta,vNam,YearLast) # Compile location, SBC, STC from various sources
up1ha.RasterizePlantingLayer(meta,zRef) # Get attributes from planting layer (to get species and genetic worth)
up1ha.DerivePlantingYearLast(meta,zRef)
u1ha.DerivePlantingMaskNonOb(meta,zRef) # Derive non-obligation planting mask
u1ha.Calc_AIL_Silv_ByFSC(meta,vNam,YearLast) # AIL by funding source
u1ha.DerivePlantingStatsByTime(meta) # AIL summary
u1ha.Calc_AIL_ByFSC_And_ASET(meta,YearLast) # Planting by ASET and year for FIP

#%% Derive artificial stand establishment type (ASET)
u1ha.DeriveASETComp1(meta,YearLast)
u1ha.RasterizeASET(meta,YearLast) # Derive last instance of Artificial Stand Establishment Type (ASET)

#%% Rasterize nutrient management (fertilization)
vNam='FE-CA'
YearLast=2023
u1ha.RasterizeSilviculture(meta,np.array(['FE']),np.array(['CA']),np.array([]),np.array([]),vNam,YearLast)
u1ha.Calc_AIL_Silv_ByFSC(meta,vNam,YearLast)
u1ha.Calc_AIL_Silv_ByBGC(meta,vNam,YearLast)
#d=gu.ipickle(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\FE-CA_SummaryByTimeAndBGCZone.pkl')

#%% Rasterize mechanical site prep (ripping)
vNam='SP-RIP'
YearLast=2024
u1ha.RasterizeSilviculture(meta,np.array(['SP']),np.array(['ME']),np.array(['DISC','EXCAV','MOUND','LRIP','RRIP','WING']),np.array([]),vNam,YearLast)
u1ha.Calc_AIL_Silv_ByFSC(meta,vNam,YearLast)
u1ha.Calc_AIL_Silv_ByBGC(meta,vNam,YearLast)
d=gu.ipickle(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + vNam + '_SummaryByTimeAndBGCZone.pkl')

#%% Rasterize mechanical site prep (WING and DISC)
YearLast=2024
vNam='SP-WING'
u1ha.RasterizeSilviculture(meta,np.array(['SP']),np.array(['ME']),np.array(['WING']),np.array([]),vNam,YearLast)
u1ha.Calc_AIL_Silv_ByFSC(meta,vNam,YearLast)
u1ha.Calc_AIL_Silv_ByBGC(meta,vNam,YearLast)
d=gu.ipickle(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + vNam + '_SummaryByTimeAndBGCZone.pkl')

vNam='SP-DISC'
u1ha.RasterizeSilviculture(meta,np.array(['SP']),np.array(['ME']),np.array(['DISC']),np.array([]),vNam,YearLast)
u1ha.Calc_AIL_Silv_ByFSC(meta,vNam,YearLast)
u1ha.Calc_AIL_Silv_ByBGC(meta,vNam,YearLast)
d=gu.ipickle(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + vNam + '_SummaryByTimeAndBGCZone.pkl')

vNam='SP-MOUND'
u1ha.RasterizeSilviculture(meta,np.array(['SP']),np.array(['ME']),np.array(['MOUND']),np.array([]),vNam,YearLast)
u1ha.Calc_AIL_Silv_ByFSC(meta,vNam,YearLast)
u1ha.Calc_AIL_Silv_ByBGC(meta,vNam,YearLast)
d=gu.ipickle(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + vNam + '_SummaryByTimeAndBGCZone.pkl')

#%% Rasterize incremental haul (fibre utilization)
# Crashed due to lack of geometry for some reason ***
vNam='SP-CHAUL'
YearLast=2024
u1ha.RasterizeSilviculture(meta,np.array(['SP']),np.array(['ME']),np.array(['CHUAL']),np.array([]),vNam,YearLast)
u1ha.Calc_AIL_Silv_ByFSC(meta,vNam,YearLast)
u1ha.Calc_AIL_Silv_ByBGC(meta,vNam,YearLast)
d=gu.ipickle(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + vNam + '_SummaryByTimeAndBGCZone.pkl')

#%% Rasterize chip and leave on site
vNam='SP-TOP'
YearLast=2023
u1ha.RasterizeSilviculture(meta,np.array(['SP']),np.array(['ME']),np.array(['BURY','CSCAT','MULCH']),np.array([]),vNam,YearLast)
u1ha.Calc_AIL_Silv_ByFSC(meta,vNam,YearLast)
u1ha.Calc_AIL_Silv_ByBGC(meta,vNam,YearLast)
d=gu.ipickle(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + vNam + '_SummaryByTimeAndBGCZone.pkl')

#%% Rasterize pile burning
vNam='SP-PBURN'
YearLast=2024
up1ha.RasterizePileBurn(meta,vNam,YearLast)
u1ha.Calc_AIL_Silv_ByFSC(meta,vNam,YearLast)
u1ha.Calc_AIL_Silv_ByBGC(meta,vNam,YearLast)
d=gu.ipickle(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + vNam + '_SummaryByTimeAndBGCZone.pkl')
d=gu.ipickle(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + vNam + '_SummaryByTimeAndFSC.pkl')
# ufcs.Plot_AIL_PileBurn(meta,pNam)

#%% Rasterize piling
YearLast=2024
vNam='SP-PILE'
up1ha.RasterizePiling(meta,vNam,YearLast)
u1ha.Calc_AIL_Silv_ByFSC(meta,vNam,YearLast)
u1ha.Calc_AIL_Silv_ByBGC(meta,vNam,YearLast)

#%% Rasterize burn landindings
vNam='SP-BU-LAND'
u1ha.RasterizeSilviculture(meta,np.array(['SP']),np.array(['BU']),np.array(['LAND']),np.array([]),vNam,YearLast)
u1ha.Calc_AIL_Silv_ByFSC(meta,vNam,YearLast)
u1ha.Calc_AIL_Silv_ByBGC(meta,vNam,YearLast)
d=gu.ipickle(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + vNam + '_SummaryByTimeAndBGCZone.pkl')
d=gu.ipickle(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + vNam + '_SummaryByTimeAndFSC.pkl')
# ufcs.Plot_AIL_PileBurn(meta,pNam)

#%% Rasterize broadcast burning (237 min)
vNam='SP-BU-BROAD'
u1ha.RasterizeSilviculture(meta,np.array(['SP']),np.array(['BU']),np.array(['BROAD']),np.array([]),vNam,YearLast)

#%% Rasterize plot surveys
vNam='SU-PLOT'
u1ha.RasterizeSilviculture(meta,np.array(['SU']),np.array(['RA','RE','RG','SR','SU']),np.array(['PLOT']),np.array([]),vNam,YearLast)
u1ha.Calc_AIL_Silv_ByFSC(meta,vNam,YearLast)
u1ha.Calc_AIL_Silv_ByBGC(meta,vNam,YearLast)

#%% Rasterize aerial spray treatment
# Wait until historical treatments entered into RESULTS
up1ha.RasterizeBTKApplication(meta)
YearLast=2023
u1ha.Calc_AIL_Silv_ByFSC(meta,'PC-CA-IDW',YearLast)
#import fcgadgets.macgyver.util_fcs_graphs as ufcs
#ufcs.Plot_AIL_Silv(meta,'PCCA-IDW')

#%% Percent dead from cruises
u1ha.RasterizeCruisePercentDead(meta)

#%% Salvage harvest mask
u1ha.HarvestSalvageMaskFromCruise(meta)

#%% Derive range tenure consolidated
u1ha.DeriveRangeTenureMask(meta)

#%% Derive crown land mask
u1ha.DeriveCrownLandMask(meta)

#%% Derive BGC Zone / NDT Combination
u1ha.DeriveBGCZoneNDTCombo(meta)

#%% Derive Tree Density Class Compilation
# *** Ensure harvested areas are listed as dense ***
u1ha.DeriveTreeDensityClass(meta)
df=u1ha.CalculateStatsByTDC(meta) # Create a table of descriptive stats

#%% Rasterize early land use change year
u1ha.RasterizeEaryLandUseChangeYear(meta)

#%% Derive distance from roads
u1ha.DeriveDistanceFromRoads(meta)

#%% Derive distance from timber facilities
u1ha.DeriveDistanceFromFacility(meta)

#%% Derive Forest Cover Reserve compilation 1
u1ha.DeriveForestCoverReserveComp1(meta)

#%% Global Forest Change Loss Year (mosaic and reproject)
u1ha.RasterizeGFC_LossYear(meta)

#%% Filter Global Forest Change Loss Year (to remove known disturbances)
u1ha.FilterGFC_LossYear(meta)

#%% Climate (1ha)
c1ha.ImportNormalsFromClimateNA(meta)
c1ha.CalcSaturationVapourPressureNormal(meta)
c1ha.CalcActualVapourPressureNormalFromTemps(meta)
c1ha.CalcActualVapourPressureNormalBiasCorrected(meta)
c1ha.CalcVapourPressureDeficitNormal(meta)
c1ha.CalcSurfaceWaterBalanceNormals(meta)
c1ha.ClimateStatsByBGCZone(meta) # Extract mean climate data by BGC zone
c1ha.HoldridgeLifeZones(meta)

#%% Deciduous fraction
u1ha.DeriveBroadleafDeciduousFraction(meta)

#%% Derive access zones
u1ha.DeriveAccessZones(meta)

#%% Species groups
z0=u1ha.Import_Raster(meta,[],['refg','lc_comp1_2019','fire_yr','bsr_sc','spc1_vri23','geomorph'],'Extract Grid')

#%% Peatlands from national 250m map
fin=r'C:\Data\BC1ha\LandCoverUse\peat3r.tif'
fout=r'C:\Data\BC1ha\LandCoverUse\peat.tif'
gis.ReprojectRasterAndClipToRaster(fin,fout,meta['Paths']['bc1ha Ref Grid'],meta['Geos']['crs'])
z=gis.OpenGeoTiff(r'C:\Data\BC1ha\LandCoverUse\peat.tif')

#%%
z=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Terrain\\FromBurtWilliam\\geomorphons_search50.tif')
plt.matshow(z['Data'],clim=[0,11])

#%% QA
flg=0
if flg==1:
	z=u1ha.Import_Raster(meta,[],['refg','lcc1_c','PROJ_AGE_1','harv_yr_con1','fire_yr'])
	yr=2023-z['PROJ_AGE_1']['Data']
	plt.close('all')
	plt.matshow(yr,clim=[1800,2000])

#%% QA
flg=0
if flg==1:
	z=u1ha.Import_Raster(meta,[],['lc_comp1_2019','lu_comp1_2049s1','lu_comp1_2049s2','lu_comp1_2019','lc_comp1_2049s1','lc_comp1_2049s2'],'Extract Grid')

ind=np.where( (z['lc_comp1_2019']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (z['lu_comp1_2019']==meta['LUT']['Derived']['lu_comp1']['Conservation Natural']) |
			  (z['lc_comp1_2019']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (z['lu_comp1_2019']==meta['LUT']['Derived']['lu_comp1']['Conservation Consistent']) )
print(ind[0].size/1e6)
ind=np.where( (z['lc_comp1_2019']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (z['lc_comp1_2049s1']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (z['lu_comp1_2049s1']==meta['LUT']['Derived']['lu_comp1']['Conservation Natural']) |
			  (z['lc_comp1_2019']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (z['lc_comp1_2049s1']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (z['lu_comp1_2049s1']==meta['LUT']['Derived']['lu_comp1']['Conservation Consistent']) )
print(ind[0].size/1e6)
ind=np.where( (z['lc_comp1_2019']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (z['lc_comp1_2049s2']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (z['lu_comp1_2049s2']==meta['LUT']['Derived']['lu_comp1']['Conservation Natural']) |
			  (z['lc_comp1_2019']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (z['lc_comp1_2049s2']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (z['lu_comp1_2049s2']==meta['LUT']['Derived']['lu_comp1']['Conservation Consistent']) )
print(ind[0].size/1e6)

#%%
p=np.zeros(10)
for i in range(10):
	ind1=np.where( (z0['geomorph']==i+1) & (z0['lc_comp1_2019']==1) )
	ind2=np.where( (z0['geomorph']==i+1) & (z0['lc_comp1_2019']==1) & (z0['fire_yr']>0) & (z0['bsr_sc']>=3) )
	p[i]=ind2[0].size/ind1[0].size*100
plt.close('all')
fig,ax=plt.subplots(1,figsize=gu.cm2inch(14,8));
plt.bar(np.arange(10),p,0.8)

#%%
p=np.zeros((4,10))
for i in range(4):
	ind1=np.where( (z0['bsr_sc']==i+1) & (z0['fire_yr']>0) )
	tmp=z0['geomorph'][ind1]
	for j in range(10):
		ind2=np.where( (tmp==j+1) )
		p[i,j]=ind2[0].size/ind1[0].size
plt.close('all')
fig,ax=plt.subplots(1,figsize=gu.cm2inch(14,8));
plt.bar(np.arange(10),p[3,:]/p[0,:],0.8)

fig,ax=plt.subplots(1,figsize=gu.cm2inch(14,8));
plt.bar(np.arange(10),p[0,:],0.8)

#%%
p=np.zeros(4)
for i in range(4):
	ind1=np.where( (z0['bsr_sc']==i+1) & (z0['fire_yr']>0) )
	ind2=np.where( (z0['bsr_sc']==i+1) & (z0['fire_yr']>0) & (z0['spc1_vri23']==meta['LUT']['VEG_COMP_LYR_R1_POLY']['SPECIES_CD_1']['DG']) )
	p[i]=ind2[0].size/ind1[0].size
plt.bar(np.arange(4),p,0.8)

#%%
p=np.zeros(10)
for i in range(1,11):
	ind1=np.where( (z0['geomorph']==i) & (z0['lc_comp1_2019']==1) )
	ind2=np.where( (z0['geomorph']==i) & (z0['lc_comp1_2019']==1) & (z0['spc1_vri23']==meta['LUT']['VEG_COMP_LYR_R1_POLY']['SPECIES_CD_1']['PLI']) )
	p[i-1]=ind2[0].size/ind1[0].size*100
plt.close('all')
fig,ax=plt.subplots(1,figsize=gu.cm2inch(14,8));
plt.bar(np.arange(1,11),p,0.8)

#%% Balsam-Deciduous Compilation
vList=['refg','lc_comp1_2019','spc1_vri23','wshed','wshedgrp']
z0=u1ha.Import_Raster(meta,[],vList,'Extract Grid')

#c,N=gu.CountByCategories(z0['geomorph'][(z0['geomorph']>0)])
#plt.bar(c,N)
plt.matshow(z0['wshedgrp'])

#%%
z1=np.zeros(zRef['Data'].shape)
ind=np.where( (z0['geomorph']==4) )
z1[ind]=1
plt.matshow(z1)

#%%

vList=['refg','refyear','lc_comp1_2019','fire_2023']
z=u1ha.Import_Raster(meta,[],vList,'Extract Grid')

ind=np.where( (z['fire_2023']>0) & (z['refyear']>0) )
ind1=np.where( (z['fire_2023']>0) & (z['refyear']==0) )
ind2=np.where( (z['fire_2023']>0) )
ind3=np.where( (z['fire_2023']>0) & (z['refyear']>2019) )
ind3[0].size/ind2[0].size*100

u,N=gu.CountByCategories(z['lc_comp1_2019'][ind1],'Percent')

2023-np.mean(z['refyear'][ind])
plt.hist(z['refyear'][ind])

#%%

vList=['feca_yr','fire_yl','spc1_vri23']
z=u1ha.Import_Raster(meta,[],vList,'Extract Grid')

z1=zRef['Data']
ind=np.where( (z['feca_yr']>0) & (z['feca_yr']<2025-10) & (z['feca_yr']>z['fire_yl']) & (z['spc1_vri23']==meta['LUT']['VEG_COMP_LYR_R1_POLY']['SPECIES_CD_1']['SX']) |
			 (z['feca_yr']>0) & (z['feca_yr']<2025-10) & (z['feca_yr']>z['fire_yl']) & (z['spc1_vri23']==meta['LUT']['VEG_COMP_LYR_R1_POLY']['SPECIES_CD_1']['SW']) )
z1[ind]=2

ind=np.where( (z['feca_yr']==2024) & (z['feca_yr']>z['fire_yl']) & (z['spc1_vri23']==meta['LUT']['VEG_COMP_LYR_R1_POLY']['SPECIES_CD_1']['SX']) |
			 (z['feca_yr']==2024) & (z['feca_yr']>z['fire_yl']) & (z['spc1_vri23']==meta['LUT']['VEG_COMP_LYR_R1_POLY']['SPECIES_CD_1']['SW']) )
z1[ind]=3

plt.close('all'); plt.matshow(z1)