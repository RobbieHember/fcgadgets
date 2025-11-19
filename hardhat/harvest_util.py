'''
HARVEST ATTRIBUTES BY TIMBER MARK - UTILITIES
	- Tenures
	- HBS
	- Waste and REsidues
	- Timber Cruise
	- RESULTS
'''
#%% Import modules
from os import listdir
import numpy as np
from scipy import stats
import pandas as pd
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf
import fiona
import fcgadgets.macgyver.util_general as gu
import fcgadgets.macgyver.util_gis as gis
import fcgadgets.bc1ha.bc1ha_utils as u1ha

#%% Import HBS
def CompileHarvest_ByTMY(meta,YearLast):
	# https://www2.gov.bc.ca/gov/content/industry/forestry/competitive-forest-industry/timber-pricing/harvest-billing-system
	# The system appears to be down sometimes
	# -> By Date of Invoice,
	# -> Can only do 12 months at a time
	# -> select "Billing processed" only, select all billing types
	
	# Definition of volume types (from Ross Pavan)
	#   1.	Normal:  the volume to the sawmill (the volume calculation is based on weight to volume ratios)
	#   2.	Cruise Based:  the volume is derived from the cruise compilation (the billed volume never exceeds the cruise volume for a cutting permit, for a timber sale license the billed volume will never exceed the cruise volume plus any additional volume added after the fact (example: external R/W volume)
	#   3.	Waste: the volume is derived from the waste assessment of a scale based cutting authority
	#   4.	Beachcomb: volume scaled because of the beachcombing efforts of an individual or salvage logger

	gradeL=list(meta['LUT']['Derived']['LogGrade'].keys())
	binY=np.arange(2007,YearLast+1,1)

	# Initialize
	n=500000
	dTMY={}
	dTMY['TM']=np.array(['' for _ in range(n)],dtype=object)
	dTMY['District']=np.array(['' for _ in range(n)],dtype=object)
	dTMY['Tenure Type']=np.array(['' for _ in range(n)],dtype=object)
	dTMY['Tenure Holder']=np.array(['' for _ in range(n)],dtype=object)
	dTMY['FTEN_Status']=np.zeros(n)
	dTMY['N Openings']=np.zeros(n)
	dTMY['DISTURBANCE_GROSS_AREA']=np.zeros(n)
	dTMY['PLANNED_NET_BLOCK_AREA']=np.zeros(n)
	dTMY['Year']=np.zeros(n)
	dTMY['V Logs m3']=np.zeros(n)
	dTMY['V Logs Abs m3']=np.zeros(n)
	dTMY['V Logs m3/ha']=np.zeros(n)
	dTMY['V Logs Abs m3/ha']=np.zeros(n)
	dTMY['V Logs Waste m3']=np.zeros(n)
	dTMY['V Logs Waste m3/ha']=np.zeros(n)
	dTMY['V NonLog m3']=np.zeros(n)
	dTMY['V NonLog m3/ha']=np.zeros(n)
	dTMY['V NonLog Abs m3/ha']=np.zeros(n)
	dTMY['V NonLog Abs m3']=np.zeros(n)
	dTMY['V NonLog Waste m3/ha']=np.zeros(n)
	dTMY['V Logs Cruise m3/ha']=np.zeros(n)
	dTMY['Stump Logs $']=np.zeros(n)
	dTMY['Stump Logs Abs $']=np.zeros(n)
	dTMY['Stump Logs $/ha']=np.zeros(n)
	dTMY['Stump Logs Abs $/ha']=np.zeros(n)
	for Grade in gradeL:
		dTMY['V Logs Grade ' + Grade + ' m3']=np.zeros(n)
		dTMY['V Logs Grade ' + Grade + ' Abs m3']=np.zeros(n)
		dTMY['V Logs Grade ' + Grade + ' m3/ha']=np.zeros(n)
		dTMY['V Logs Grade ' + Grade + ' Abs m3/ha']=np.zeros(n)
		dTMY['Stump Logs Grade ' + Grade + ' $']=np.zeros(n)
		dTMY['Stump Logs Grade ' + Grade + ' Abs $']=np.zeros(n)
		dTMY['Stump Logs Grade ' + Grade + ' $/ha']=np.zeros(n)
		dTMY['Stump Logs Abs Grade ' + Grade + ' $/ha']=np.zeros(n)

	# Populate from annual summary files
	cnt=0
	for iY in range(binY.size):
		print(str(binY[iY]) + ' '  + str(cnt))
		df=pd.read_csv(meta['Paths']['DB']['HBS'] + '\\HBS ' + str(binY[iY]) + '.csv',header=None)
		dHB={}
		for i in range(len(df.columns)):
			dHB[i]=df.iloc[:,i].to_numpy()
	
		# Create indices for each TM
		dHB[0]=dHB[0].astype('U')
		iTM=gu.IndicesFromUniqueArrayValues(dHB[0])
		for k in iTM.keys():

			# Extract attributes for TM and year
			ind=iTM[k]
			TM=dHB[0][ind]
			Mat=dHB[6][ind]
			Type=dHB[9][ind]
			V=dHB[11][ind]
			V_abs=np.abs(dHB[11][ind])
			Stump=dHB[12][ind]
			Stump_abs=np.abs(dHB[12][ind])
			Grade=dHB[7][ind]
			Ten=dHB[19][ind]
			Holder=dHB[26][ind]
			Dist=dHB[15][ind]

			dTMY['Year'][cnt]=binY[iY]
			dTMY['TM'][cnt]=TM[0]
			dTMY['District'][cnt]=Dist[0]
			dTMY['Tenure Type'][cnt]=Ten[0]
			dTMY['Tenure Holder'][cnt]=Holder[0]

			# logs
			ind2=np.where( (Mat=='Logs') )[0]
			if ind2.size>0:
				dTMY['V Logs m3'][cnt]=np.round(np.sum(V[ind2]),decimals=0)
				dTMY['V Logs Abs m3'][cnt]=np.round(np.sum(V_abs[ind2]),decimals=0)
				dTMY['Stump Logs $'][cnt]=np.round(np.sum(Stump[ind2]),decimals=0)
				dTMY['Stump Logs Abs $'][cnt]=np.round(np.sum(Stump_abs[ind2]),decimals=0)
	
			# Logs by grade
			for iGrade in range(len(gradeL)):
				ind2=np.where( (Mat=='Logs') & (Grade==gradeL[iGrade])  )[0]
				if ind2.size>0:
					dTMY['V Logs Grade ' + gradeL[iGrade] + ' m3'][cnt]=np.round(np.sum(V_abs[ind2]))
					dTMY['V Logs Grade ' + gradeL[iGrade] + ' Abs m3'][cnt]=np.round(np.sum(V_abs[ind2]))
					dTMY['Stump Logs Grade ' + gradeL[iGrade] + ' $'][cnt]=np.round(np.sum(Stump[ind2]))
					dTMY['Stump Logs Grade ' + gradeL[iGrade] + ' Abs $'][cnt]=np.round(np.sum(Stump_abs[ind2]))

			# Logs Waste
			ind2=np.where( (Type=='WA') & (Mat=='Logs') | (Type=='WU') & (Mat=='Logs') )[0]
			if ind2.size>0:
				dTMY['V Logs Waste m3'][cnt]=np.round(np.sum(V[ind2]),decimals=0)
	
			# Non-logs
			ind2=np.where( (Mat!='Logs') )[0]
			if ind2.size>0:
				dTMY['V NonLog m3'][cnt]=np.round(np.sum(V[ind2]),decimals=0)
				dTMY['V NonLog Abs m3'][cnt]=np.round(np.sum(V_abs[ind2]),decimals=0)

			# Update counter
			cnt=cnt+1

	# Remove remaining empty space
	ikp=np.where(dTMY['TM']!='')[0]
	for k in dTMY.keys():
		dTMY[k]=dTMY[k][ikp]

	return dTMY

#%%
def GetHarvestAreaFromFTEN(meta,dTMY):
	iTM=gu.IndicesFromUniqueArrayValues(dTMY['TM'])
	n=500000
	dOP={}
	dOP['TM']=np.array(['' for _ in range(n)],dtype=object)
	dOP['OPENING_ID']=np.array(['' for _ in range(n)],dtype=object)
	dOP['DISTURBANCE_GROSS_AREA']=np.zeros(n)
	dOP['PLANNED_NET_BLOCK_AREA']=np.zeros(n)
	cntOP=0
	with fiona.open(meta['Paths']['GDB']['Disturbance'],layer='FTEN_CUT_BLOCK_OPEN_ADMIN') as source:
		for feat in source:
			p=feat['properties']
			if p['TIMBER_MARK']==None:
				continue
			dOP['TM'][cntOP]=p['TIMBER_MARK']
			if p['OPENING_ID']!=None:
				dOP['OPENING_ID'][cntOP]=p['OPENING_ID']
			dOP['DISTURBANCE_GROSS_AREA'][cntOP]=p['DISTURBANCE_GROSS_AREA']
			dOP['PLANNED_NET_BLOCK_AREA'][cntOP]=p['PLANNED_NET_BLOCK_AREA']
			cntOP=cntOP+1
	# Remove remaining empty space
	ikp=np.where(dOP['TM']!='')[0]
	for k in dOP.keys():
		dOP[k]=dOP[k][ikp]

	tm=dOP['TM'].copy().astype('U')
	for k in iTM.keys():
		ind1=iTM[k]
		ind2=np.where(tm==k)[0]
		if ind2.size==0:
			print('Missing')
			continue
		dTMY['DISTURBANCE_GROSS_AREA'][ind1]=np.nansum(dOP['DISTURBANCE_GROSS_AREA'][ind2])/ind1.size
		dTMY['PLANNED_NET_BLOCK_AREA'][ind1]=np.nansum(dOP['PLANNED_NET_BLOCK_AREA'][ind2])/ind1.size
		dTMY['N Openings'][ind1]=ind2.size

	return dTMY,dOP

#%%
def CalcHarvest_ByTM_FromTMY(meta,dTMY):
	iTMY=gu.IndicesFromUniqueArrayValues(dTMY['TM'])
	n=len(iTMY)
	dTM={}
	for k in dTMY.keys():
		if dTMY[k].dtype=='O':
			dTM[k]=np.tile(dTMY[k][0],n).astype('O')
		else:
			dTM[k]=np.tile(dTMY[k][0],n)
	dTM['TM']=np.array(list(iTMY.keys()),dtype='O')
	for tm in iTMY.keys():
		ind=np.where(dTM['TM']==tm)[0]
		for k in dTMY.keys():
			if dTMY[k].dtype=='O':
				dTM[k][ind]=np.array(dTMY[k][iTMY[tm]][0]).astype('O')
			else:
				dTM[k][ind]=np.sum(dTMY[k][iTMY[tm]])

	# Add derived variables
	del dTM['Year']
	dTM['Year Min']=np.zeros(dTM['TM'].size)
	dTM['Year Max']=np.zeros(dTM['TM'].size)
	for tm in iTMY.keys():
		ind=np.where(dTM['TM']==tm)[0]
		dTM['Year Min'][ind]=np.min(dTMY['Year'][iTMY[tm]])
		dTM['Year Max'][ind]=np.max(dTMY['Year'][iTMY[tm]])
	dTM['V Logs m3/ha']=dTM['V Logs m3']/dTM['DISTURBANCE_GROSS_AREA']
	dTM['V Logs Abs m3/ha']=dTM['V Logs Abs m3']/dTM['DISTURBANCE_GROSS_AREA']
	dTM['V Logs Waste m3/ha']=dTM['V Logs Waste m3']/dTM['DISTURBANCE_GROSS_AREA']
	dTM['V NonLog m3/ha']=dTM['V NonLog m3']/dTM['DISTURBANCE_GROSS_AREA']
	dTM['V NonLog Abs m3/ha']=dTM['V NonLog Abs m3']/dTM['DISTURBANCE_GROSS_AREA']
	gradeL=list(meta['LUT']['Derived']['LogGrade'].keys())
	for Grade in gradeL:
		dTM['V Logs Grade ' + Grade + ' m3/ha']=dTM['V Logs Grade ' + Grade + ' m3']/dTM['DISTURBANCE_GROSS_AREA']
		dTM['V Logs Grade ' + Grade + ' Abs m3/ha']=dTM['V Logs Grade ' + Grade + ' Abs m3']/dTM['DISTURBANCE_GROSS_AREA']
		dTM['Stump Logs Grade ' + Grade + ' $/ha']=dTM['Stump Logs Grade ' + Grade + ' $']/dTM['DISTURBANCE_GROSS_AREA']
		dTM['Stump Logs Grade ' + Grade + ' Abs $/ha']=dTM['Stump Logs Grade ' + Grade + ' Abs $']/dTM['DISTURBANCE_GROSS_AREA']

	return dTM

#%%
def ImportWaste(meta,dTMY,dTM):
	# Just inlcude status = "Billing Issued"
	# Don't download more than 3 years at a time
	# Don't download the regions - tempting but too big
	
	fL=listdir(meta['Paths']['DB']['Waste'])
	
	N=2000000
	dW={}
	dW['TM']=np.array(['' for _ in range(N)],dtype=object)
	dW['Year']=np.zeros(N)
	dW['Net Area Ha']=np.zeros(N)
	dW['Total A. Sawlog  Volume m3']=np.zeros(N)
	dW['A. Grade Y/4 Volume m3']=np.zeros(N)
	dW['Av. + Unav. All Grades Volume m3']=np.zeros(N)
	dW['Disp_A']=np.zeros(N)
	dW['Disp_V']=np.zeros(N)
	dW['Accu_A']=np.zeros(N)
	dW['Accu_V']=np.zeros(N)
	dW['Stan_A']=np.zeros(N)
	dW['Stan_V']=np.zeros(N)
	dW['Waste Billing']=np.zeros(N)
	cnt=0
	for iF in range(len(fL)):
		df=pd.read_excel(meta['Paths']['DB']['Waste'] + '\\' + fL[iF],skiprows=list(range(0,14)))
		#d=gu.ReadExcel(meta['Paths']['DB']['Waste'] + '\\' + fL[iF],skiprows=14 )
		dW['TM'][cnt:cnt+len(df)]=df['TM'].to_numpy()
		dW['Year'][cnt:cnt+len(df)]=np.array(fL[iF][4:8],dtype='float')
		dW['Net Area Ha'][cnt:cnt+len(df)]=df['Net Area Ha'].to_numpy()
		dW['Total A. Sawlog  Volume m3'][cnt:cnt+len(df)]=df['Total A. Sawlog  Volume m3'].to_numpy()
		dW['A. Grade Y/4 Volume m3'][cnt:cnt+len(df)]=df['A. Grade Y/4 Volume m3'].to_numpy()
		dW['Av. + Unav. All Grades Volume m3'][cnt:cnt+len(df)]=df['Av. + Unav. All Grades Volume m3'].to_numpy()
		dW['Disp_A'][cnt:cnt+len(df)]=df['Area\nHa'].to_numpy()
		dW['Disp_V'][cnt:cnt+len(df)]=df[' Volume m3'].to_numpy()
		dW['Accu_A'][cnt:cnt+len(df)]=df['Area\nHa.1'].to_numpy()
		dW['Accu_V'][cnt:cnt+len(df)]=df[' Volume m3.1'].to_numpy()
		dW['Stan_A'][cnt:cnt+len(df)]=df['Area\nHa.2'].to_numpy()
		dW['Stan_V'][cnt:cnt+len(df)]=df[' Volume m3.2'].to_numpy()
		dW['Waste Billing'][cnt:cnt+len(df)]=df['Waste $Billing'].to_numpy()
		cnt=cnt+len(df)

	# Truncate
	for k in dW.keys():
		dW[k]=dW[k][0:cnt]

	# It crashes later if this is not done
	tmU=dW['TM'].copy().astype('U')

	#==========================================================================
	# Add waste data to TMY database
	dTMY['Waste N Entries']=np.zeros(dTMY['TM'].size)
	dTMY['Waste Net Area Tot']=np.zeros(dTMY['TM'].size)
	dTMY['Waste Total m3']=np.zeros(dTMY['TM'].size)
	dTMY['Waste Total m3/ha']=np.zeros(dTMY['TM'].size)
	#dTMY['Waste Sawlog m3/ha']=np.zeros(dTMY['TM'].size)
	#dTMY['Waste GradeY/4 m3/ha']=np.zeros(dTMY['TM'].size)
	dTMY['Waste Dispersed m3/ha']=np.zeros(dTMY['TM'].size)
	dTMY['Waste Accumulation m3/ha']=np.zeros(dTMY['TM'].size)
	dTMY['Waste Standing m3/ha']=np.zeros(dTMY['TM'].size)
	dTMY['Waste Bill $']=np.zeros(dTMY['TM'].size)
	for i in range(dTMY['TM'].size):
		indW=np.where( (dW['TM']==dTMY['TM'][i]) & (dW['Year']==dTMY['Year'][i]) )[0]
		dTMY['Waste N Entries'][i]=indW.size
		if indW.size>0:
			A_tot=np.sum(dW['Net Area Ha'][indW])
			dTMY['Waste Net Area Tot'][i]=A_tot
			dTMY['Waste Total m3'][i]=np.round(np.sum(dW['Av. + Unav. All Grades Volume m3'][indW]),decimals=0)
			dTMY['Waste Total m3/ha'][i]=np.round(np.sum(dW['Av. + Unav. All Grades Volume m3'][indW])/A_tot,decimals=0)
			#dTMY['Waste Sawlog m3/ha'][i]=np.round(np.sum(dW['Total A. Sawlog  Volume m3'][indW])/A_tot,decimals=0)
			#dTMY['Waste GradeY/4 m3/ha'][i]=np.round(np.sum(dW['A. Grade Y/4 Volume m3'][indW])/A_tot,decimals=0)
			dTMY['Waste Dispersed m3/ha'][i]=np.round(np.sum(dW['Disp_V'][indW])/A_tot,decimals=0)
			dTMY['Waste Accumulation m3/ha'][i]=np.round(np.sum(dW['Accu_V'][indW])/A_tot,decimals=0)
			dTMY['Waste Standing m3/ha'][i]=np.round(np.sum(dW['Stan_V'][indW])/A_tot,decimals=0)
			dTMY['Waste Bill $'][i]=np.sum(dW['Waste Billing'][indW])
	
	# Percent fate of waste
	tot=dTMY['Waste Dispersed m3/ha']+dTMY['Waste Accumulation m3/ha']+dTMY['Waste Standing m3/ha']
	dTMY['Waste Dispersed %']=np.round(dTMY['Waste Dispersed m3/ha']/tot*100)
	dTMY['Waste Accumulation %']=np.round(dTMY['Waste Accumulation m3/ha']/tot*100)
	dTMY['Waste Standing %']=np.round(dTMY['Waste Standing m3/ha']/tot*100)

	#==========================================================================
	# Add waste data to TM database
	dTM['Waste N Entries']=np.zeros(dTM['TM'].size)
	dTM['Waste Net Area Tot']=np.zeros(dTM['TM'].size)
	dTM['Waste Total m3']=np.zeros(dTM['TM'].size)
	dTM['Waste Total m3/ha']=np.zeros(dTM['TM'].size)
	#dTM['Waste Sawlog m3/ha']=np.zeros(dTM['TM'].size)
	#dTM['Waste GradeY/4 m3/ha']=np.zeros(dTM['TM'].size)
	dTM['Waste Dispersed m3/ha']=np.zeros(dTM['TM'].size)
	dTM['Waste Accumulation m3/ha']=np.zeros(dTM['TM'].size)
	dTM['Waste Standing m3/ha']=np.zeros(dTM['TM'].size)
	dTM['Waste Bill $']=np.zeros(dTM['TM'].size)

	# Create indices for each TM
	iTM=gu.IndicesFromUniqueArrayValues(dTM['TM'])
	iW=gu.IndicesFromUniqueArrayValues(tmU)
	for tm in iTM.keys():
		if tm not in iW.keys():
			# TM not active this year, continue
			continue
		indW=iW[tm]
		dTM['Waste N Entries'][iTM[tm]]=indW.size
		if indW.size>0:
			A_tot=np.sum(dW['Net Area Ha'][indW])
			dTM['Waste Net Area Tot'][iTM[tm]]=A_tot
			dTM['Waste Total m3'][iTM[tm]]=np.round(np.sum(dW['Av. + Unav. All Grades Volume m3'][indW]),decimals=0)
			dTM['Waste Total m3/ha'][iTM[tm]]=np.round(np.sum(dW['Av. + Unav. All Grades Volume m3'][indW])/A_tot,decimals=0)
			#dTM['Waste Sawlog m3/ha'][iTM[tm]]=np.round(np.sum(dW['Total A. Sawlog  Volume m3'][indW])/A_tot,decimals=0)
			#dTM['Waste GradeY/4 m3/ha'][iTM[tm]]=np.round(np.sum(dW['A. Grade Y/4 Volume m3'][indW])/A_tot,decimals=0)
			dTM['Waste Dispersed m3/ha'][iTM[tm]]=np.round(np.sum(dW['Disp_V'][indW])/A_tot,decimals=0)
			dTM['Waste Accumulation m3/ha'][iTM[tm]]=np.round(np.sum(dW['Accu_V'][indW])/A_tot,decimals=0)
			dTM['Waste Standing m3/ha'][iTM[tm]]=np.round(np.sum(dW['Stan_V'][indW])/A_tot,decimals=0)
			dTM['Waste Bill $'][iTM[tm]]=np.sum(dW['Waste Billing'][indW])
	
	# Percent fate of waste
	tot=dTM['Waste Dispersed m3/ha']+dTM['Waste Accumulation m3/ha']+dTM['Waste Standing m3/ha']
	dTM['Waste Dispersed %']=np.round(dTM['Waste Dispersed m3/ha']/tot*100)
	dTM['Waste Accumulation %']=np.round(dTM['Waste Accumulation m3/ha']/tot*100)
	dTM['Waste Standing %']=np.round(dTM['Waste Standing m3/ha']/tot*100)

	return dTMY,dTM

#%%
def ImportTimberCruise(meta,dTMY,dTM):
	dTC=gu.ipickle(r'C:\Data\ECAS\Received 2023-04-04\rh_COMP_DATA_CleanCombined.pkl')
	vL=['NET_AREA','Year','Age','V Gross (m3/ha)','V Net (m3/ha)','Pct Dead Net','PCT Net','PCT_DIST','PCT_DECAY','PCT_WASTE','PCT_WASTE_BILL','PCT_BREAK','PCT_DWB','PCT_NET_2NDGROWTH','PCT_INSECT_VOL','PCT_NO_INSECT_M3','PCT_GREEN_INSECT_M3', \
		'PCT_RED_INSECT_M3','PCT_GREY_INSECT_M3','PCT_OTHER_INSECT_M3','X_DEFOLIATOR_LIVE_CAMB_PCT','Y_DEFOLIATOR_DEAD_CAMB_PCT']
	for v in vL:
		dTMY['Cruise_' + v]=np.nan*np.ones(dTMY['TM'].size)
		dTM['Cruise_' + v]=np.nan*np.ones(dTM['TM'].size)
	uTM=np.unique(dTC['PRIMARY_MARK'])
	for iTM in range(uTM.size):
		ind1=np.where(dTC['PRIMARY_MARK']==uTM[iTM])[0]
		ind2=np.where(dTMY['TM']==uTM[iTM])[0]
		for v in vL:
			dTMY['Cruise_' + v][ind2]=dTC[v][ind1]
		ind2=np.where(dTM['TM']==uTM[iTM])[0]
		for v in vL:
			dTM['Cruise_' + v][ind2]=dTC[v][ind1]
	return dTMY,dTM

#%%
def ImportRESULTS(meta,dOP,dTM):
	# fiona.listlayers(Paths['Results'])
	dOP['SILV_SYSTEM_CODE']=np.zeros(dOP['TM'].size,dtype='int8')
	with fiona.open(meta['Paths']['GDB']['Results'],layer='RSLT_OPENING_SVW') as source:
		for feat in source:
			p=feat['properties']
			ind=np.where(dOP['OPENING_ID']==p['OPENING_ID'])[0]
			if ind.size==0:
				continue
			#dOp['PREV_AGE_CLASS_CODE'][ind]=p['PREV_AGE_CLASS_CODE']
			if p['DENUDATION_1_SILV_SYSTEM_CODE']!=None:
				dOP['SILV_SYSTEM_CODE'][ind]=meta['LUT']['RSLT_OPENING_SVW']['DENUDATION_1_SILV_SYSTEM_CODE'][p['DENUDATION_1_SILV_SYSTEM_CODE']]
			if p['DENUDATION_2_SILV_SYSTEM_CODE']!=None:
				dOP['SILV_SYSTEM_CODE'][ind]=meta['LUT']['RSLT_OPENING_SVW']['DENUDATION_1_SILV_SYSTEM_CODE'][p['DENUDATION_2_SILV_SYSTEM_CODE']]

# 	# Convert age class to age (years)
# 	AgeClassAge=[10,30,50,70,90,110,130,195,250]
# 	dOp['Age RSLTS']=np.nan*np.ones(dOp['OPENING_ID'].size)
# 	for i in range(len(AgeClassAge)):
# 		ind=np.where(dOp['PREV_AGE_CLASS_CODE']==np.array(i-1))[0]
# 		dOp['Age RSLTS'][ind]=AgeClassAge[i]
	
	# Add to TM database
	dTM['SILV_SYSTEM_CODE ID 1']=np.zeros(dTM['TM'].size)
	dTM['SILV_SYSTEM_CODE % 1']=np.zeros(dTM['TM'].size)
	dTM['SILV_SYSTEM_CODE ID 2']=np.zeros(dTM['TM'].size)
	dTM['SILV_SYSTEM_CODE % 2']=np.zeros(dTM['TM'].size)
	dTM['SILV_SYSTEM_CODE ID 3']=np.zeros(dTM['TM'].size)
	dTM['SILV_SYSTEM_CODE % 3']=np.zeros(dTM['TM'].size)
	dTM['SILV_SYSTEM_CODE ID 4']=np.zeros(dTM['TM'].size)
	dTM['SILV_SYSTEM_CODE % 4']=np.zeros(dTM['TM'].size)
	dTM['SILV_SYSTEM_CODE ID 5']=np.zeros(dTM['TM'].size)
	dTM['SILV_SYSTEM_CODE % 5']=np.zeros(dTM['TM'].size)
	dTM['SILV_SYSTEM_CODE ID 6']=np.zeros(dTM['TM'].size)
	dTM['SILV_SYSTEM_CODE % 6']=np.zeros(dTM['TM'].size)
	for iTM in range(dTM['TM'].size):
		ind=np.where( (dOP['TM']==dTM['TM'][iTM]) )[0]
		if ind.size==0:
			continue
		Atot=np.sum(dOP['PLANNED_NET_BLOCK_AREA'][ind])
		A_pct=dOP['PLANNED_NET_BLOCK_AREA'][ind]/Atot*100
		ord=np.argsort(A_pct)
		A_pct[ord]
		ssc=dOP['SILV_SYSTEM_CODE'][ind][ord]
		u=np.unique(ssc)
		for iU in range(u.size):
			if iU+1>6:
				continue
			ind1=np.where(ssc==u[iU])[0]
			if ind1.size==1:
				dTM['SILV_SYSTEM_CODE ID ' + str(iU+1)][iTM]=ssc[ind1]
				dTM['SILV_SYSTEM_CODE % ' + str(iU+1)][iTM]=A_pct[ind1]
			else:
				dTM['SILV_SYSTEM_CODE ID ' + str(iU+1)][iTM]=ssc[ind1[0]]
				dTM['SILV_SYSTEM_CODE % ' + str(iU+1)][iTM]=np.sum(A_pct[ind1])

	# Add SSC names
	dTM['SSC']=np.array(['' for _ in range(dTM['TM'].size)],dtype=object)
	u=np.unique(dTM['SILV_SYSTEM_CODE ID 1'])
	for iU in range(u.size):
		if u[iU]==0:
			continue
		ind=np.where(dTM['SILV_SYSTEM_CODE ID 1']==u[iU])[0]
		dTM['SSC'][ind]=u1ha.lut_n2s(meta['LUT']['RSLT_OPENING_SVW']['DENUDATION_1_SILV_SYSTEM_CODE'],u[iU])
	return dTM

#%%
def CreateDerivedVariables(meta,dTMY,dTM):

	# By TMY

	# Region
	dTMY['Region']=np.array(['Interior' for _ in range(dTMY['TM'].size)],dtype=object)
	ind=np.where( (dTMY['District']=='DCR') | (dTMY['District']=='DNI') | (dTMY['District']=='DSI') | \
			(dTMY['District']=='DSC') | (dTMY['District']=='DSQ') | (dTMY['District']=='DQC') | \
			(dTMY['District']=='DCK') )[0]
	dTMY['Region'][ind]='Coast'

	# Volume felled (Volume + Waste)
	dTMY['V Felled m3']=dTMY['V Logs Abs m3']+dTMY['Waste Total m3']
	dTMY['V Felled m3/ha']=dTMY['V Logs Abs m3/ha']+dTMY['Waste Total m3/ha']

	# Waste Fraction
	dTMY['Waste Fraction']=dTMY['Waste Total m3/ha']/dTMY['V Felled m3/ha']
	dTMY['Waste Fraction From Tot']=dTMY['Waste Total m3']/dTMY['V Felled m3']

	# Sawlog recovery ratio (based on total felled, not total logs)
	dTMY['Sawlog Fraction']=np.zeros(dTMY['TM'].size)
	ind=np.where( (dTMY['Region']=='Interior') )[0]
	dTMY['Sawlog Fraction'][ind]=(dTMY['V Logs Grade 1 Abs m3/ha'][ind]+dTMY['V Logs Grade 2 Abs m3/ha'][ind])/dTMY['V Felled m3/ha'][ind]
	ind=np.where( (dTMY['Region']=='Coast') )[0]
	dTMY['Sawlog Fraction'][ind]=(dTMY['V Logs Grade B Abs m3/ha'][ind]+
						dTMY['V Logs Grade C Abs m3/ha'][ind]+
						dTMY['V Logs Grade D Abs m3/ha'][ind]+
						dTMY['V Logs Grade E Abs m3/ha'][ind]+
						dTMY['V Logs Grade F Abs m3/ha'][ind]+
						dTMY['V Logs Grade G Abs m3/ha'][ind]+
						dTMY['V Logs Grade H Abs m3/ha'][ind]+
						dTMY['V Logs Grade I Abs m3/ha'][ind]+
						dTMY['V Logs Grade J Abs m3/ha'][ind]+
						dTMY['V Logs Grade K Abs m3/ha'][ind]+
						dTMY['V Logs Grade L Abs m3/ha'][ind]+
						dTMY['V Logs Grade M Abs m3/ha'][ind])/dTMY['V Felled m3/ha'][ind]

	dTMY['Secondary Log Fraction']=1.0-dTMY['Waste Fraction']-dTMY['Sawlog Fraction']

	# By TM

	# Region
	dTM['Region']=np.array(['Interior' for _ in range(dTM['TM'].size)],dtype=object)
	ind=np.where( (dTM['District']=='DCR') | (dTM['District']=='DNI') | (dTM['District']=='DSI') | \
			(dTM['District']=='DSC') | (dTM['District']=='DSQ') | (dTM['District']=='DQC') | \
			(dTM['District']=='DCK') )[0]
	dTM['Region'][ind]='Coast'

	# Volume felled (Volume + Waste)
	#dTM['V Felled m3/ha']=dTM['V Logs Abs m3/ha']+dTM['V NonLog Abs m3/ha']+dTM['Waste Total m3/ha']
	dTM['V Felled m3']=dTM['V Logs Abs m3/ha']+dTM['Waste Total m3']
	dTM['V Felled m3/ha']=dTM['V Logs Abs m3/ha']+dTM['Waste Total m3/ha']

	# Waste rate
	dTM['Waste Fraction']=dTM['Waste Total m3/ha']/dTM['V Felled m3/ha']
	dTM['Waste Fraction From Tot']=dTM['Waste Total m3']/dTM['V Felled m3']

	# Sawlog recovery ratio
	dTM['Sawlog Fraction']=np.zeros(dTM['TM'].size)
	ind=np.where( (dTM['Region']=='Interior') )[0]
	dTM['Sawlog Fraction'][ind]=(dTM['V Logs Grade 1 Abs m3/ha'][ind]+dTM['V Logs Grade 2 Abs m3/ha'][ind])/dTM['V Felled m3/ha'][ind]
	ind=np.where( (dTM['Region']=='Coast') )[0]
	dTM['Sawlog Fraction'][ind]=(dTM['V Logs Grade B Abs m3/ha'][ind]+
						dTM['V Logs Grade C Abs m3/ha'][ind]+
						dTM['V Logs Grade D Abs m3/ha'][ind]+
						dTM['V Logs Grade E Abs m3/ha'][ind]+
						dTM['V Logs Grade F Abs m3/ha'][ind]+
						dTM['V Logs Grade G Abs m3/ha'][ind]+
						dTM['V Logs Grade H Abs m3/ha'][ind]+
						dTM['V Logs Grade I Abs m3/ha'][ind]+
						dTM['V Logs Grade J Abs m3/ha'][ind]+
						dTM['V Logs Grade K Abs m3/ha'][ind]+
						dTM['V Logs Grade L Abs m3/ha'][ind]+
						dTM['V Logs Grade M Abs m3/ha'][ind])/dTM['V Felled m3/ha'][ind]

	dTM['Secondary Log Fraction']=1.0-dTM['Waste Fraction']-dTM['Sawlog Fraction']

	# Discrepency in area values between HBS and WS
	dTM['Delta A %']=(dTM['DISTURBANCE_GROSS_AREA']-dTM['Waste Net Area Tot'])/dTM['DISTURBANCE_GROSS_AREA']*100

	return dTMY,dTM

#%% Add raster variables
def GetVariablesFromRasterDB(meta,dTMY,dTM):
	z=u1ha.Import_Raster(meta,[],['tm','bgcz','age_vri02'],'Extract Grid')
	idxTM=gu.IndicesFromUniqueArrayValues(z['tm'].flatten())
	iTM=gu.IndicesFromUniqueArrayValues(dTM['TM'])
	iTMY=gu.IndicesFromUniqueArrayValues(dTMY['TM'])
	bgcz=z['bgcz'].flatten()
	age=z['age_vri02'].flatten()
	dTMY['ID_BGCZ']=np.zeros(dTMY['TM'].size,dtype='int16')
	dTM['ID_BGCZ']=np.zeros(dTM['TM'].size,dtype='int16')
	dTM['age_vri02']=np.zeros(dTM['TM'].size,dtype='int16')
	N=0
	for k in idxTM.keys():
		try:
			tm=u1ha.lut_n2s(meta['LUT']['FTEN_CUT_BLOCK_POLY_SVW']['TIMBER_MARK'],k)[0]
		except:
			# This just crashes once times on '0'
			continue

		if tm not in iTM.keys():
			N=N+1
			continue
		Mode=stats.mode(bgcz[idxTM[k]])[0]
		dTM['ID_BGCZ'][iTM[tm]]=Mode
		dTM['age_vri02'][iTM[tm]]=np.mean(age[idxTM[k]])

		dTMY['ID_BGCZ'][iTMY[tm]]=Mode

	# Correct age
	Adj=dTM['Year Max']-2002
	ind=np.where( (dTM['age_vri02']>0) & (dTM['age_vri02']<300) & (dTM['Year Max']>0) )[0]
	dTM['AgeC']=np.nan*np.ones(dTM['TM'].size)
	dTM['AgeC'][ind]=dTM['age_vri02'][ind]+Adj[ind]

	# Add BGCZ names
	dTM['BGCZ']=np.array(['' for _ in range(dTM['TM'].size)],dtype=object)
	u=np.unique(dTM['ID_BGCZ'])
	for iU in range(u.size):
		if u[iU]==0:
			continue
		ind=np.where(dTM['ID_BGCZ']==u[iU])[0]
		dTM['BGCZ'][ind]=u1ha.lut_n2s(meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE'],u[iU])

	# Add BGCZ names
	dTMY['BGCZ']=np.array(['' for _ in range(dTMY['TM'].size)],dtype=object)
	u=np.unique(dTMY['ID_BGCZ'])
	for iU in range(u.size):
		if u[iU]==0:
			continue
		ind=np.where(dTMY['ID_BGCZ']==u[iU])[0]
		dTMY['BGCZ'][ind]=u1ha.lut_n2s(meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE'],u[iU])

	return dTMY,dTM

#%%
def CalcHarvestTimeSeries(meta,YearLast):
	# Area in Kha/year
	# Volume in Mm3/year

	# Initialize dictionary
	d={}
	d['Year']=np.arange(1850,YearLast+1,1)
	d['Province']={'Area':{},'Volume':{}}
	d['ByBGC']={'Area':{},'Volume':{}}

	d['ByBGC']['Area']['CC']={}
	d['ByBGC']['Area']['CC+EH']={}
	d['ByBGC']['Volume']['HBS']={}
	for zone in meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE'].keys():
		d['ByBGC']['Area']['CC'][zone]=np.nan*np.ones(d['Year'].size)
		d['ByBGC']['Area']['CC+EH'][zone]=np.nan*np.ones(d['Year'].size)
		d['ByBGC']['Volume']['HBS'][zone]=np.nan*np.ones(d['Year'].size)

	# Import rasters
	vList=['bgcz','harv_yr_early']
	zG=u1ha.Import_Raster(meta,[],vList,'Extract Grid')

	# Harvest area from consolidated cutblock DB
	lNam='VEG_CONSOLIDATED_CUT_BLOCKS_SP'
	vNam='HARVEST_YEAR'
	nPack=6
	tv,N=u1ha.TimeSeriesOccurrenceFromPackedEvents(meta,lNam,vNam,nPack)
	d['Province']['Area']['CC']=np.nan*np.ones(d['Year'].size)
	d['Province']['Area']['CC+EH']=np.nan*np.ones(d['Year'].size)
	for iT in range(d['Year'].size):
		ind=np.where(tv==d['Year'][iT])[0]
		if ind.size==0:
			continue
		if N[ind]>0:
			d['Province']['Area']['CC'][iT]=N[ind]/1000
			d['Province']['Area']['CC+EH'][iT]=N[ind]/1000

		# Add early harvesting
		ind=np.where(zG['harv_yr_early']==d['Year'][iT])[0]
		d['Province']['Area']['CC+EH'][iT]=d['Province']['Area']['CC+EH'][iT]+ind.size/1000

	# Harvest area from consolidated cutblock DB (by BGC zone)
	for i in range(nPack):
		zY=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\' + lNam + '\\' + vNam + '_' + str(i+1) + '_Year.tif')['Data']
		uY=np.unique(zY[zY>0])
		for zone in meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE'].keys():
			iZone=np.where( (zG['bgcz']==meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE'][zone]) )
			zYz=zY[iZone]
			zHE=zG['harv_yr_early'][iZone]
			for yr in uY:
				iT=np.where(d['Year']==yr)[0]
				iS=np.where( (zYz==yr) )
				if np.isnan(d['ByBGC']['Area']['CC'][zone][iT])==True:
					d['ByBGC']['Area']['CC'][zone][iT]=iS[0].size/1000
					d['ByBGC']['Area']['CC+EH'][zone][iT]=iS[0].size/1000
				else:
					d['ByBGC']['Area']['CC'][zone][iT]=d['ByBGC']['Area']['CC'][zone][iT]+iS[0].size/1000
					d['ByBGC']['Area']['CC+EH'][zone][iT]=d['ByBGC']['Area']['CC'][zone][iT]+iS[0].size/1000

				ind=np.where(zHE==d['Year'][iT])[0]
				d['ByBGC']['Area']['CC+EH'][zone][iT]=d['ByBGC']['Area']['CC+EH'][zone][iT]+ind.size/1000

	# Harvest area from National Forest Database
	d0=gu.ReadExcel(r'C:\Data\Harvest\Harvest Area\NFD - Area Harvested by ownership and harvesting method - EN FR.xlsx')
	d['Province']['Area']['NFD']=np.nan*np.ones(d['Year'].size)
	for i,yr in enumerate(d['Year']):
		ind=np.where( (d0['Jurisdiction']=='British Columbia') & (d0['Year']==yr) )[0]
		if ind.size==0:
			continue
		d['Province']['Area']['NFD'][i]=np.nansum(d0['Area (hectares)'][ind])/1000

	# Harvest area from cruise comp
	dCru=gu.ipickle(r'C:\Data\ECAS\Received 2023-04-04\rh_COMP_DATA_CleanCombined.pkl')
	d['Province']['Area']['Cruise']=np.nan*np.ones(d['Year'].size)
	for i,yr in enumerate(d['Year']):
		ind=np.where( (dCru['Year']==yr) )[0]
		if ind.size==0:
			continue
		d['Province']['Area']['Cruise'][i]=np.nansum(dCru['NET_AREA'][ind])/1000

	# Harvest area from NTEMS
	zH=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Disturbances\\Harvest_NTEM_Year.tif')['Data']
	zH=zH[zH>0]
	d['Province']['Area']['NTEMS']=np.nan*np.ones(d['Year'].size)
	for i,yr in enumerate(d['Year']):
		ind=np.where(zH==yr)[0]
		if ind.size==0:
			continue
		d['Province']['Area']['NTEMS'][i]=ind.size/1000

	# Harvest volume from FLNR website
	dV_FLNR=gu.ReadExcel(r'C:\Data\Harvest\SummaryDataChangeInTimberHarvest\bctimberharvest.xlsx')
	d['Province']['Volume']['FLNR']=np.nan*np.ones(d['Year'].size)
	for i,yr in enumerate(dV_FLNR['Year']):
		ind=np.where(d['Year']==yr)[0]
		if ind.size==0:
			continue
		d['Province']['Volume']['FLNR'][ind]=dV_FLNR['Total_harvest_millions_m3'][i]

	# Harvest volume from HBS
	dTMY=gu.ipickle(meta['Paths']['DB']['Harvest'] + '\\HavestSummary_ByTMAndYear.pkl')
	dT=CalcStatsByTime(meta,dTMY)
	d['Province']['Volume']['HBS']=np.nan*np.ones(d['Year'].size)
	for i,yr in enumerate(dT['Time']):
		ind=np.where(d['Year']==yr)[0]
		if ind.size==0:
			continue
		d['Province']['Volume']['HBS'][ind]=(dT['Sum']['V Logs Abs m3'][i]+dT['Sum']['V NonLog Abs m3'][i])/1e6

	for zone in meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE'].keys():
		dT=CalcStatsByTime(meta,dTMY,zone=zone)
		for i,yr in enumerate(dT['Time']):
			ind=np.where(d['Year']==yr)[0]
			if ind.size==0:
				continue
			if np.isnan(d['ByBGC']['Volume']['HBS'][zone][ind])==True:
				d['ByBGC']['Volume']['HBS'][zone][ind]=(dT['Sum']['V Logs Abs m3'][i]+dT['Sum']['V NonLog Abs m3'][i])/1e6
			else:
				d['ByBGC']['Volume']['HBS'][zone][ind]=d['ByBGC']['Volume']['HBS'][zone][ind]+(dT['Sum']['V Logs Abs m3'][i]+dT['Sum']['V NonLog Abs m3'][i])/1e6

	# QA that ByBGC values sum up to province
	# Notes: ByBGC less than province by 8%! Must be because of missing BGC zone data
	flg=0
	if flg==1:
		iT=np.where(d['Year']==2015)[0]
		V=0
		for zone in meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE'].keys():
			V=V+d['ByBGC']['Volume']['HBS'][zone][iT]
		print(V)
		print(d['Province']['Volume']['HBS'][iT])

	# Save time series of harvest area
	gu.opickle(meta['Paths']['DB']['Harvest'] + '\\HarvestTimeSeries.pkl',d)

	return

#%%
def CalcPilingRateByBGCZone(meta,dTM,**kwargs):
	Tot=dTM['Waste Accumulation %']+dTM['Waste Dispersed %']+dTM['Waste Standing %']
	#plt.hist(Tot)
	ikp=np.where( (dTM['Waste Fraction']>0) & (dTM['Waste Fraction']<=1) & \
				 (dTM['Waste N Entries']>0) & \
				(Tot==100) & \
				(dTM['BGCZ']!='') & \
				(dTM['Waste Accumulation %']>=0) & (dTM['Waste Accumulation %']<=100) )[0]
	sts1=gu.StatsByCategories(dTM['BGCZ'][ikp],dTM['Waste Accumulation %'][ikp])
	sts2=gu.StatsByCategories(dTM['BGCZ'][ikp],dTM['Waste Standing %'][ikp])
	sts3=gu.StatsByCategories(dTM['BGCZ'][ikp],dTM['Waste Dispersed %'][ikp])

	# Save to file
	if kwargs['save']=='On':
		d={}
		for i in range(sts1['u'].size):
			d[sts1['u'][i]]=[sts1['mu'][i]/100]
		df=pd.DataFrame.from_dict(d)
		df.to_excel(meta['Paths']['Model']['Parameters'] + '\\Parameters_PilingRate.xlsx',index=False)

	# Plot
	if kwargs['plot']=='On':
		ord=np.flip(np.argsort(sts1['mu']))
		for k in sts1.keys():
			sts1[k]=sts1[k][ord]
			sts2[k]=sts2[k][ord]
			sts3[k]=sts3[k][ord]
		plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(22,4.5));
		ax.bar(np.arange(sts1['u'].size),sts1['mu'],fc=[0.6,0.75,1],label='Piled')
		ax.bar(np.arange(sts1['u'].size),sts2['mu'],bottom=sts1['mu'],fc=[0.4,0.55,0.8],label='Left standing')
		ax.bar(np.arange(sts1['u'].size),sts3['mu'],bottom=sts1['mu']+sts2['mu'],fc=[0.2,0.35,0.6],label='Dispersed')
		ax.set(xticks=np.arange(sts1['u'].size),xticklabels=sts1['u'],ylabel='Frequency (%)',xlabel='BGC Zone',xlim=[-0.5,sts1['u'].size-0.5])
		ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
		ax.legend(loc='lower left',facecolor=[1,1,1],frameon=True);
		if meta['Graphics']['Print Figures']=='On':
			gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\PilingRatesByBGCZone','png',900)
	return

#%%
def WasteFractionRegression(meta,d,**kwargs):

	# Impute % dead
	ind=np.where( (np.isnan(d['Cruise_Pct Dead Net'])==True) )[0];
	d['Cruise_Pct Dead Net'][ind]=np.nanmean(d['Cruise_Pct Dead Net'])

	# Filter
	ikp=np.where( (d['Waste Fraction']>=0) & (d['Waste Fraction']<=1) & (d['Waste N Entries']>0) & \
				(d['Waste Accumulation %']>=0) & (d['Waste Accumulation %']<=100) & \
				(d['Year Max']>1950) & (d['Year Max']<=2022) & \
				(d['BGCZ']!='') & (d['Tenure Type']!='') & (d['SSC']!='') & \
				(np.isnan(d['Cruise_Pct Dead Net'])==False) )[0]
	#print(ikp.size)
	df=pd.DataFrame({'WR':100*d['Waste Fraction'][ikp],
					 'BGCZ':d['BGCZ'][ikp],
					 'Time':d['Year Max'][ikp],
					 'Age':d['AgeC'][ikp],
					 'SSC':d['SSC'][ikp],
					 'Ten':d['Tenure Type'][ikp],
					 'PctDead':d['Cruise_Pct Dead Net'][ikp]})
	#print(df.describe())
	md=smf.ols("WR~Time+Age+PctDead+C(BGCZ)+C(SSC)+C(Ten)",data=df)
	rs=md.fit(maxiter=100)
	b=rs.params

	# Show table
	if kwargs['print_table']=='On':
		print(rs.summary())

	# Save to file
	if kwargs['save']=='On':
		df=pd.DataFrame.from_dict(b)
		df.to_excel(meta['Paths']['Model']['Parameters'] + '\\Parameters_WasteFractionRegression.xlsx',index=False)
	
	# Plot behaviour of WR model
	if kwargs['plot']=='On':
		Age=np.arange(1,300)
		pdD=10; yrD=2020; bgcD='C(BGCZ)[T.SBS]'; sscD='C(SSC)[T.CLEAR]'; tenD='C(Ten)[T.Forest Licence]'
		plt.close('all'); fig,ax=plt.subplots(2,3,figsize=gu.cm2inch(22,8));
		PctDead=0; Year=yrD; bgc=bgcD; ssc=sscD; ten=tenD
		y=np.maximum(0,b['Intercept']+b['Age']*Age+b['Time']*Year+b['PctDead']*PctDead+b[ssc]+b[ten]+b[bgc])
		ax[0,0].plot(Age,y,'k-',label='% dead = 0')
		PctDead=100; Year=yrD; bgc=bgcD; ssc=sscD; ten=tenD
		y=np.maximum(0,b['Intercept']+b['Age']*Age+b['Time']*Year+b['PctDead']*PctDead+b[ssc]+b[ten]+b[bgc])
		ax[0,0].plot(Age,y,'k--',label='% dead = 100')
		ax[0,0].legend(loc='center left',facecolor=[1,1,1],frameon=False);
		ax[0,0].set(xlabel='Age, years',ylabel='Waste fraction (%)')
		ax[0,0].yaxis.set_ticks_position('both'); ax[0,0].xaxis.set_ticks_position('both'); ax[0,0].tick_params(length=meta['Graphics']['gp']['tickl'])
		
		PctDead=pdD; Year=2010; bgc=bgcD; ssc=sscD; ten=tenD
		y=np.maximum(0,b['Intercept']+b['Age']*Age+b['Time']*Year+b['PctDead']*PctDead+b[ssc]+b[ten]+b[bgc])
		ax[0,1].plot(Age,y,'k-',label='2010')
		PctDead=pdD; Year=2024; bgc=bgcD; ssc=sscD; ten=tenD
		y=np.maximum(0,b['Intercept']+b['Age']*Age+b['Time']*Year+b['PctDead']*PctDead+b[ssc]+b[ten]+b[bgc])
		ax[0,1].plot(Age,y,'k--',label='2024')
		ax[0,1].legend(loc='center left',facecolor=[1,1,1],frameon=False);
		ax[0,1].set(xlabel='Age, years',ylabel='Waste fraction (%)')
		ax[0,1].yaxis.set_ticks_position('both'); ax[0,1].xaxis.set_ticks_position('both'); ax[0,1].tick_params(length=meta['Graphics']['gp']['tickl'])
		
		PctDead=pdD; Year=yrD; bgc='C(BGCZ)[T.MH]'; ssc=sscD; ten=tenD
		y=np.maximum(0,b['Intercept']+b['Age']*Age+b['Time']*Year+b['PctDead']*PctDead+b[ssc]+b[ten]+b[bgc])
		ax[0,2].plot(Age,y,'k-',label='MH')
		PctDead=pdD; Year=yrD; bgc='C(BGCZ)[T.CWH]'; ssc=sscD; ten=tenD
		y=np.maximum(0,b['Intercept']+b['Age']*Age+b['Time']*Year+b['PctDead']*PctDead+b[ssc]+b[ten]+b[bgc])
		ax[0,2].plot(Age,y,'k--',label='CWH')
		PctDead=pdD; Year=yrD; bgc='C(BGCZ)[T.IDF]'; ssc=sscD; ten=tenD
		y=np.maximum(0,b['Intercept']+b['Age']*Age+b['Time']*Year+b['PctDead']*PctDead+b[ssc]+b[ten]+b[bgc])
		ax[0,2].plot(Age,y,'k:',label='IDF')
		PctDead=pdD; Year=yrD; bgc='C(BGCZ)[T.MS]'; ssc=sscD; ten=tenD
		y=np.maximum(0,b['Intercept']+b['Age']*Age+b['Time']*Year+b['PctDead']*PctDead+b[ssc]+b[ten]+b[bgc])
		ax[0,2].plot(Age,y,'k-.',label='MS')
		ax[0,2].legend(loc='center left',facecolor=[1,1,1],frameon=False);
		ax[0,2].set(xlabel='Age, years',ylabel='Waste fraction (%)')
		ax[0,2].yaxis.set_ticks_position('both'); ax[0,2].xaxis.set_ticks_position('both'); ax[0,2].tick_params(length=meta['Graphics']['gp']['tickl'])
		
		PctDead=pdD; Year=yrD; bgc=bgcD; ssc='C(SSC)[T.CLEAR]'; ten=tenD
		y=np.maximum(0,b['Intercept']+b['Age']*Age+b['Time']*Year+b['PctDead']*PctDead+b[ssc]+b[ten]+b[bgc])
		ax[1,0].plot(Age,y,'k-',label='CLEAR')
		PctDead=pdD; Year=yrD; bgc=bgcD; ssc='C(SSC)[T.IMCUT]'; ten=tenD
		y=np.maximum(0,b['Intercept']+b['Age']*Age+b['Time']*Year+b['PctDead']*PctDead+b[ssc]+b[ten]+b[bgc])
		ax[1,0].plot(Age,y,'k--',label='IMCUT')
		PctDead=pdD; Year=yrD; bgc=bgcD; ssc='C(SSC)[T.RETEN]'; ten=tenD
		y=np.maximum(0,b['Intercept']+b['Age']*Age+b['Time']*Year+b['PctDead']*PctDead+b[ssc]+b[ten]+b[bgc])
		ax[1,0].plot(Age,y,'k:',label='RETEN')
		PctDead=pdD; Year=yrD; bgc=bgcD; ssc='C(SSC)[T.SEEDT]'; ten=tenD
		y=np.maximum(0,b['Intercept']+b['Age']*Age+b['Time']*Year+b['PctDead']*PctDead+b[ssc]+b[ten]+b[bgc])
		ax[1,0].plot(Age,y,'k-.',label='SEEDT')
		PctDead=pdD; Year=yrD; bgc=bgcD; ssc='C(SSC)[T.SELEC]'; ten=tenD
		y=np.maximum(0,b['Intercept']+b['Age']*Age+b['Time']*Year+b['PctDead']*PctDead+b[ssc]+b[ten]+b[bgc])
		ax[1,0].plot(Age,y,'k-',lw=1.5,label='SELEC')
		ax[1,0].legend(loc='center left',facecolor=[1,1,1],frameon=False);
		ax[1,0].set(xlabel='Age, years',ylabel='Waste fraction (%)')
		ax[1,0].yaxis.set_ticks_position('both'); ax[1,0].xaxis.set_ticks_position('both'); ax[1,0].tick_params(length=meta['Graphics']['gp']['tickl'])
		
		PctDead=pdD; Year=yrD; bgc=bgcD; ssc=sscD; ten='C(Ten)[T.Forest Licence]'
		y=np.maximum(0,b['Intercept']+b['Age']*Age+b['Time']*Year+b['PctDead']*PctDead+b[ssc]+b[ten]+b[bgc])
		ax[1,1].plot(Age,y,'k-',label='Forest Licence')
		PctDead=pdD; Year=yrD; bgc=bgcD; ssc=sscD; ten='C(Ten)[T.Forestry Licence to Cut]'
		y=np.maximum(0,b['Intercept']+b['Age']*Age+b['Time']*Year+b['PctDead']*PctDead+b[ssc]+b[ten]+b[bgc])
		ax[1,1].plot(Age,y,'k--',label='Forestry Licence to Cut')
		PctDead=pdD; Year=yrD; bgc=bgcD; ssc=sscD; ten='C(Ten)[T.SB TSL S20 single mark]'
		y=np.maximum(0,b['Intercept']+b['Age']*Age+b['Time']*Year+b['PctDead']*PctDead+b[ssc]+b[ten]+b[bgc])
		ax[1,1].plot(Age,y,'k:',label='SB SSL S20 Single Mark')
		PctDead=pdD; Year=yrD; bgc=bgcD; ssc=sscD; ten='C(Ten)[T.Tree Farm Licence]'
		y=np.maximum(0,b['Intercept']+b['Age']*Age+b['Time']*Year+b['PctDead']*PctDead+b[ssc]+b[ten]+b[bgc])
		ax[1,1].plot(Age,y,'k-.',label='Tree Farm Licence')
		PctDead=pdD; Year=yrD; bgc=bgcD; ssc=sscD; ten='C(Ten)[T.Woodlot Licence]'
		y=np.maximum(0,b['Intercept']+b['Age']*Age+b['Time']*Year+b['PctDead']*PctDead+b[ssc]+b[ten]+b[bgc])
		ax[1,1].plot(Age,y,'k-',lw=1.5,label='Woodlot Licence')
		ax[1,1].set(xlabel='Age, years',ylabel='Waste fraction (%)')
		ax[1,1].legend(loc='center left',facecolor=[1,1,1],frameon=False);
		ax[1,1].yaxis.set_ticks_position('both'); ax[1,1].xaxis.set_ticks_position('both'); ax[1,1].tick_params(length=meta['Graphics']['gp']['tickl'])
		
		PctDead=np.arange(0,100,1); Year=yrD; bgc=bgcD; ssc=sscD; ten=tenD
		y=np.maximum(0,b['Intercept']+b['Age']*60+b['Time']*Year+b['PctDead']*PctDead+b[ssc]+b[ten]+b[bgc])
		ax[1,2].plot(PctDead,y,'k-',label='Age = 60')
		y=np.maximum(0,b['Intercept']+b['Age']*200+b['Time']*Year+b['PctDead']*PctDead+b[ssc]+b[ten]+b[bgc])
		ax[1,2].plot(PctDead,y,'k--',label='Age = 200')
		ax[1,2].set(xlabel='Proportion of dead trees (%)',ylabel='Waste fraction (%)')
		ax[1,2].legend(loc='center left',facecolor=[1,1,1],frameon=False);
		ax[1,2].yaxis.set_ticks_position('both'); ax[1,2].xaxis.set_ticks_position('both'); ax[1,2].tick_params(length=meta['Graphics']['gp']['tickl'])
		
		gu.axletters(ax,plt,0.03,0.88,FontColor=meta['Graphics']['gp']['cla'],LetterStyle=meta['Graphics']['Modelling']['AxesLetterStyle'],FontWeight=meta['Graphics']['Modelling']['AxesFontWeight'])
		plt.tight_layout()
	return

#%% Regression of Sawlog Fraction
def SawlogFractionRegression(meta,d,**kwargs):
	# Impute % dead
	ind=np.where( (np.isnan(d['Cruise_Pct Dead Net'])==True) )[0];
	d['Cruise_Pct Dead Net'][ind]=np.nanmean(d['Cruise_Pct Dead Net'])
	
	# Filter
	ikp=np.where( (d['Sawlog Fraction']>=0) & (d['Sawlog Fraction']<=1) & \
				(d['Year Max']>1950) & (d['Year Max']<=2022) & \
				(d['BGCZ']!='') & (d['Tenure Type']!='') & (d['SSC']!='') & \
				(np.isnan(d['Cruise_Pct Dead Net'])==False) )[0]
	#print(ikp.size)
	df=pd.DataFrame({'SR':100*d['Sawlog Fraction'][ikp],
					 'BGCZ':d['BGCZ'][ikp],
					 'Time':d['Year Max'][ikp],
					 'Age':d['AgeC'][ikp],
					 'SSC':d['SSC'][ikp],
					 'Ten':d['Tenure Type'][ikp],
					 'PctDead':d['Cruise_Pct Dead Net'][ikp]})

	md=smf.ols("SR~Time+Age+PctDead+C(BGCZ)+C(Ten)+C(SSC)",data=df)
	rs=md.fit(maxiter=100)
	b=rs.params

	# Show table
	if kwargs['print_table']=='On':
		print(rs.summary())

	# Save to file
	if kwargs['save']=='On':
		df=pd.DataFrame.from_dict(b)
		df.to_excel(meta['Paths']['Model']['Parameters'] + '\\Parameters_SawlogFractionRegression.xlsx',index=False)
	
	# Plot behaviour of model
	if kwargs['plot']=='On':
		Age=np.arange(1,300)
		pdD=10; yrD=2020; bgcD='C(BGCZ)[T.SBS]'; sscD='C(SSC)[T.CLEAR]'; tenD='C(Ten)[T.Forest Licence]'
		plt.close('all'); fig,ax=plt.subplots(2,3,figsize=gu.cm2inch(22,8));
		PctDead=0; Year=yrD; bgc=bgcD; ssc=sscD; ten=tenD
		y=np.maximum(0,b['Intercept']+b['Age']*Age+b['Time']*Year+b['PctDead']*PctDead+b[ssc]+b[ten]+b[bgc])
		ax[0,0].plot(Age,y,'k-',label='% dead = 0')
		PctDead=100; Year=yrD; bgc=bgcD; ssc=sscD; ten=tenD
		y=np.maximum(0,b['Intercept']+b['Age']*Age+b['Time']*Year+b['PctDead']*PctDead+b[ssc]+b[ten]+b[bgc])
		ax[0,0].plot(Age,y,'k--',label='% dead = 100')
		ax[0,0].legend(loc='center left',facecolor=[1,1,1],frameon=False);
		ax[0,0].set(xlabel='Age, years',ylabel='Sawlog fraction (%)')
		ax[0,0].yaxis.set_ticks_position('both'); ax[0,0].xaxis.set_ticks_position('both'); ax[0,0].tick_params(length=meta['Graphics']['gp']['tickl'])
		
		PctDead=pdD; Year=2010; bgc=bgcD; ssc=sscD; ten=tenD
		y=np.maximum(0,b['Intercept']+b['Age']*Age+b['Time']*Year+b['PctDead']*PctDead+b[ssc]+b[ten]+b[bgc])
		ax[0,1].plot(Age,y,'k-',label='2010')
		PctDead=pdD; Year=2024; bgc=bgcD; ssc=sscD; ten=tenD
		y=np.maximum(0,b['Intercept']+b['Age']*Age+b['Time']*Year+b['PctDead']*PctDead+b[ssc]+b[ten]+b[bgc])
		ax[0,1].plot(Age,y,'k--',label='2024')
		ax[0,1].legend(loc='center left',facecolor=[1,1,1],frameon=False);
		ax[0,1].set(xlabel='Age, years',ylabel='Sawlog fraction (%)')
		ax[0,1].yaxis.set_ticks_position('both'); ax[0,1].xaxis.set_ticks_position('both'); ax[0,1].tick_params(length=meta['Graphics']['gp']['tickl'])
		
		PctDead=pdD; Year=yrD; bgc='C(BGCZ)[T.MH]'; ssc=sscD; ten=tenD
		y=np.maximum(0,b['Intercept']+b['Age']*Age+b['Time']*Year+b['PctDead']*PctDead+b[ssc]+b[ten]+b[bgc])
		ax[0,2].plot(Age,y,'k-',label='MH')
		PctDead=pdD; Year=yrD; bgc='C(BGCZ)[T.CWH]'; ssc=sscD; ten=tenD
		y=np.maximum(0,b['Intercept']+b['Age']*Age+b['Time']*Year+b['PctDead']*PctDead+b[ssc]+b[ten]+b[bgc])
		ax[0,2].plot(Age,y,'k--',label='CWH')
		PctDead=pdD; Year=yrD; bgc='C(BGCZ)[T.IDF]'; ssc=sscD; ten=tenD
		y=np.maximum(0,b['Intercept']+b['Age']*Age+b['Time']*Year+b['PctDead']*PctDead+b[ssc]+b[ten]+b[bgc])
		ax[0,2].plot(Age,y,'k:',label='IDF')
		PctDead=pdD; Year=yrD; bgc='C(BGCZ)[T.MS]'; ssc=sscD; ten=tenD
		y=np.maximum(0,b['Intercept']+b['Age']*Age+b['Time']*Year+b['PctDead']*PctDead+b[ssc]+b[ten]+b[bgc])
		ax[0,2].plot(Age,y,'k-.',label='MS')
		ax[0,2].legend(loc='center left',facecolor=[1,1,1],frameon=False);
		ax[0,2].set(xlabel='Age, years',ylabel='Sawlog fraction (%)')
		ax[0,2].yaxis.set_ticks_position('both'); ax[0,2].xaxis.set_ticks_position('both'); ax[0,2].tick_params(length=meta['Graphics']['gp']['tickl'])
		
		PctDead=pdD; Year=yrD; bgc=bgcD; ssc='C(SSC)[T.CLEAR]'; ten=tenD
		y=np.maximum(0,b['Intercept']+b['Age']*Age+b['Time']*Year+b['PctDead']*PctDead+b[ssc]+b[ten]+b[bgc])
		ax[1,0].plot(Age,y,'k-',label='CLEAR')
		PctDead=pdD; Year=yrD; bgc=bgcD; ssc='C(SSC)[T.IMCUT]'; ten=tenD
		y=np.maximum(0,b['Intercept']+b['Age']*Age+b['Time']*Year+b['PctDead']*PctDead+b[ssc]+b[ten]+b[bgc])
		ax[1,0].plot(Age,y,'k--',label='IMCUT')
		PctDead=pdD; Year=yrD; bgc=bgcD; ssc='C(SSC)[T.RETEN]'; ten=tenD
		y=np.maximum(0,b['Intercept']+b['Age']*Age+b['Time']*Year+b['PctDead']*PctDead+b[ssc]+b[ten]+b[bgc])
		ax[1,0].plot(Age,y,'k:',label='RETEN')
		PctDead=pdD; Year=yrD; bgc=bgcD; ssc='C(SSC)[T.SEEDT]'; ten=tenD
		y=np.maximum(0,b['Intercept']+b['Age']*Age+b['Time']*Year+b['PctDead']*PctDead+b[ssc]+b[ten]+b[bgc])
		ax[1,0].plot(Age,y,'k-.',label='SEEDT')
		PctDead=pdD; Year=yrD; bgc=bgcD; ssc='C(SSC)[T.SELEC]'; ten=tenD
		y=np.maximum(0,b['Intercept']+b['Age']*Age+b['Time']*Year+b['PctDead']*PctDead+b[ssc]+b[ten]+b[bgc])
		ax[1,0].plot(Age,y,'k-',lw=1.5,label='SELEC')
		ax[1,0].legend(loc='center left',facecolor=[1,1,1],frameon=False);
		ax[1,0].set(xlabel='Age, years',ylabel='Sawlog fraction (%)')
		ax[1,0].yaxis.set_ticks_position('both'); ax[1,0].xaxis.set_ticks_position('both'); ax[1,0].tick_params(length=meta['Graphics']['gp']['tickl'])
		
		PctDead=pdD; Year=yrD; bgc=bgcD; ssc=sscD; ten='C(Ten)[T.Forest Licence]'
		y=np.maximum(0,b['Intercept']+b['Age']*Age+b['Time']*Year+b['PctDead']*PctDead+b[ssc]+b[ten]+b[bgc])
		ax[1,1].plot(Age,y,'k-',label='Forest Licence')
		PctDead=pdD; Year=yrD; bgc=bgcD; ssc=sscD; ten='C(Ten)[T.Forestry Licence to Cut]'
		y=np.maximum(0,b['Intercept']+b['Age']*Age+b['Time']*Year+b['PctDead']*PctDead+b[ssc]+b[ten]+b[bgc])
		ax[1,1].plot(Age,y,'k--',label='Forestry Licence to Cut')
		PctDead=pdD; Year=yrD; bgc=bgcD; ssc=sscD; ten='C(Ten)[T.SB TSL S20 single mark]'
		y=np.maximum(0,b['Intercept']+b['Age']*Age+b['Time']*Year+b['PctDead']*PctDead+b[ssc]+b[ten]+b[bgc])
		ax[1,1].plot(Age,y,'k:',label='SB SSL S20 Single Mark')
		PctDead=pdD; Year=yrD; bgc=bgcD; ssc=sscD; ten='C(Ten)[T.Tree Farm Licence]'
		y=np.maximum(0,b['Intercept']+b['Age']*Age+b['Time']*Year+b['PctDead']*PctDead+b[ssc]+b[ten]+b[bgc])
		ax[1,1].plot(Age,y,'k-.',label='Tree Farm Licence')
		PctDead=pdD; Year=yrD; bgc=bgcD; ssc=sscD; ten='C(Ten)[T.Woodlot Licence]'
		y=np.maximum(0,b['Intercept']+b['Age']*Age+b['Time']*Year+b['PctDead']*PctDead+b[ssc]+b[ten]+b[bgc])
		ax[1,1].plot(Age,y,'k-',lw=1.5,label='Woodlot Licence')
		ax[1,1].set(xlabel='Age, years',ylabel='Sawlog fraction (%)')
		ax[1,1].legend(loc='center left',facecolor=[1,1,1],frameon=False);
		ax[1,1].yaxis.set_ticks_position('both'); ax[1,1].xaxis.set_ticks_position('both'); ax[1,1].tick_params(length=meta['Graphics']['gp']['tickl'])
		
		PctDead=np.arange(0,100,1); Year=yrD; bgc=bgcD; ssc=sscD; ten=tenD
		y=np.maximum(0,b['Intercept']+b['Age']*60+b['Time']*Year+b['PctDead']*PctDead+b[ssc]+b[ten]+b[bgc])
		ax[1,2].plot(PctDead,y,'k-',label='Age = 60')
		y=np.maximum(0,b['Intercept']+b['Age']*200+b['Time']*Year+b['PctDead']*PctDead+b[ssc]+b[ten]+b[bgc])
		ax[1,2].plot(PctDead,y,'k--',label='Age = 200')
		ax[1,2].set(xlabel='Proportion of dead trees (%)',ylabel='Sawlog fraction (%)')
		ax[1,2].legend(loc='center left',facecolor=[1,1,1],frameon=False);
		ax[1,2].yaxis.set_ticks_position('both'); ax[1,2].xaxis.set_ticks_position('both'); ax[1,2].tick_params(length=meta['Graphics']['gp']['tickl'])
		
		gu.axletters(ax,plt,0.03,0.88,FontColor=meta['Graphics']['gp']['cla'],LetterStyle=meta['Graphics']['Modelling']['AxesLetterStyle'],FontWeight=meta['Graphics']['Modelling']['AxesFontWeight'])
		plt.tight_layout()
	return

#%%
def CalcStatsByDistrict(meta,dTM):
	iD=gu.IndicesFromUniqueArrayValues(dTM['District'])
	del iD['']
	distL=list(iD.keys())
	d={}
	for i in iD.keys():
		d[i]={}
		for v in dTM.keys():
			if dTM[v].dtype=='O': continue
			ikp=np.where(dTM[v][iD[i]]<2000)[0]
			d[i][v]=np.nanmean(dTM[v][iD[i]][ikp])
	df=pd.DataFrame(d)
	df.loc['Waste Fraction']
	df.loc['Sawlog Fraction']
	df.loc['V Felled m3/ha']
	return d

#%%
def CalcStatsByBGCZone(meta,dTM):
	iD=gu.IndicesFromUniqueArrayValues(dTM['BGCZ'])
	#del iD['']
	distL=list(iD.keys())
	d={}
	for i in iD.keys():
		d[i]={}
		for v in dTM.keys():
			if dTM[v].dtype=='O': continue
			d[i][v]=np.nanmean(dTM[v][iD[i]])
	df=pd.DataFrame(d)
	df.loc['Waste Fraction']
	return d

#%%
def CalcStatsBySILV_SYSTEM_CODE_And_Region(meta,dTM,reg):
	d={}
	for k in meta['LUT']['RSLT_OPENING_SVW']['DENUDATION_1_SILV_SYSTEM_CODE'].keys():
		#(dTM['ID TSA']==meta['LUT']['FADM_TSA']['TSA_NUMBER_DESCRIPTION']['Boundary TSA']) & \
		ind=np.where( (dTM['PLANNED_NET_BLOCK_AREA']>0) & \
			(dTM['Region']==reg) & \
			(dTM['Sawlog Fraction']>=0) & (dTM['Sawlog Fraction']<=1) & \
			(dTM['Waste Fraction']>=0) & (dTM['Waste Fraction']<=1) & \
			(dTM['V Logs m3/ha']>0) & (dTM['V Logs m3/ha']<3000) & \
			(dTM['SILV_SYSTEM_CODE % 1']>90) & \
			(dTM['SILV_SYSTEM_CODE ID 1']==meta['LUT']['RSLT_OPENING_SVW']['DENUDATION_1_SILV_SYSTEM_CODE'][k]) )[0]
		if ind.size>3:
			d[k]=[ind.size,
			  np.round(np.nanmean(dTM['V Felled m3/ha'][ind]),decimals=1),
			  np.round(np.nanmean(dTM['Sawlog Fraction'][ind]+dTM['Secondary Log Fraction'][ind])*100,decimals=1),
			  np.round(np.nanmean(dTM['Sawlog Fraction'][ind])*100,decimals=1),
			  np.round(np.nanmean(dTM['Secondary Log Fraction'][ind])*100,decimals=1),
			  np.round(np.nanmean(dTM['Waste Fraction'][ind])*100,decimals=1)]
			#d[k]=np.mean(dTM['Waste Fraction'][ind])
			#d[k]=np.mean(dTM['Sawlog Fraction'][ind])
		else:
			d[k]=[0,0,0,0,0,0]
	df=pd.DataFrame(d).T
	df.columns=['Sample size (# TMs)','Total removals (m3/ha)','Total removals (%)','Sawlogs (%)','Secondary logs (%)','Waste wood (%)']
	#df=df.astype(int)
	df=df.reset_index()
	df=df.rename(columns={'index':'SSC'})
	return df

#%%
def CalcStatsBySILV_SYSTEM_CODE_ForBGCZone(meta,dTM,zone):
	d={}
	for k in meta['LUT']['RSLT_OPENING_SVW']['DENUDATION_1_SILV_SYSTEM_CODE'].keys():
		ind=np.where( (dTM['BGCZ']==zone) & (dTM['V Logs m3/ha']>0) & (dTM['V Logs m3/ha']<3000) & (dTM['SILV_SYSTEM_CODE % 1']>90) & \
			   (dTM['Sawlog Fraction']>=0) & (dTM['Sawlog Fraction']<=1) & \
			   (dTM['SILV_SYSTEM_CODE ID 1']==meta['LUT']['RSLT_OPENING_SVW']['DENUDATION_1_SILV_SYSTEM_CODE'][k]) )[0]
		if ind.size>3:
			d[k]=[int(ind.size),
			 np.round(np.nanmean(dTM['V Felled m3/ha'][ind]),decimals=1),
			  np.round(np.nanmean(dTM['Sawlog Fraction'][ind]*100),decimals=0),
			  np.round(np.nanmean(1-dTM['Sawlog Fraction'][ind]-dTM['Waste Fraction'][ind])*100,decimals=0),
			  np.round(np.nanmean(dTM['Waste Fraction'][ind]*100),decimals=0)]
			#d[k]=np.mean(dTM['Waste Fraction'][ind])
			#d[k]=np.mean(dTM['Sawlog Fraction'][ind])
		else:
			d[k]=[np.nan,np.nan,np.nan,np.nan]
	#pd.DataFrame(d).T
	df=pd.DataFrame(d).T
	df.columns=['Sample size (# TMs)','Yield (m3/ha)','Waste wood (%)','Sawlog Fraction (%)']
	return df

#%%
def CalcStatsByTime(meta,dTMY,**kwargs):
	#ind=np.where( (dTM['ID_BGCZ']==6) & (dTM['V Logs Abs m3/ha']<1500) )[0]
	tv=np.arange(2007,2025,1)
	dT={'Time':tv,'Mean':{},'Sum':{}}
	for k in dTMY.keys():
		dT['Mean'][k]=np.zeros(tv.size)
		dT['Sum'][k]=np.zeros(tv.size)
	for iT in range(tv.size):
		if 'zone' in kwargs.keys():
			ikp=np.where( (dTMY['Year']==tv[iT]) & (dTMY['BGCZ']==kwargs['zone']) )[0]
		elif 'zones' in kwargs.keys():
			ikp=np.where( (dTMY['Year']==tv[iT]) & (np.isin(dTMY['BGCZ'],kwargs['zones'])==True) )[0]
		elif 'District' in kwargs.keys():
			ikp=np.where( (dTMY['Year']==tv[iT]) & (dTMY['District']==kwargs['District']) )[0]
		else:
			ikp=np.where( (dTMY['Year']==tv[iT]) )[0]
		for k in dTMY.keys():
			tmp=dTMY[k].copy()
			if (k=='Waste Fraction') | (k=='Sawlog Fraction'):
				ind=np.where((tmp<0) | (tmp>1.0))[0]
				tmp[ind]=np.nan
			if (dTMY[k].dtype=='O'):
				continue
			dT['Mean'][k][iT]=np.nanmean(dTMY[k][ikp])
			dT['Sum'][k][iT]=np.nansum(dTMY[k][ikp])
	return dT

#%%
def Plot_NetVolumeComparison(meta,dTM):
	list(dTM.keys())
	ind=np.where( (dTM['Cruise_V Net (m3/ha)']>0) & (dTM['V Logs Abs m3/ha']>0) & (dTM['V Logs Abs m3/ha']<1500) )[0]
	x=dTM['Cruise_V Net (m3/ha)'][ind]
	y=dTM['V Logs Abs m3/ha'][ind]
	rs,txt=gu.GetRegStats(x,y,'No Intercept')
	
	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(8,8))
	ax.plot(x,y,'b.',mew=0.25,ms=3,mfc='b',mec='w')
	ax.plot([0,2000],[0,2000],'k-',lw=2,color=[0.8,0.8,0.8])
	ax.text(1300,1300,'1:1',fontsize=10,ha='center',va='center')
	ax.plot(rs['xhat Line'],rs['yhat Line'],'r--',lw=1.5,color=[0,0,0])
	ax.text(1400,200,txt,fontsize=8,ha='right',va='center')
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=2)
	ax.set(xlabel='Net vol from cruise (m3/ha)',ylabel='Net volume (m3/ha)',xlim=[0,1500],ylim=[0,1500])
	plt.tight_layout()
	#gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Harvest\Timber Cruise Info\TimeSeries','png',900)
	return

#%%
def Plot_HarvestAreaTimeSeries(meta,tlim):

	d=gu.ipickle(meta['Paths']['DB']['Harvest'] + '\\HarvestTimeSeries.pkl')
	iT=np.where(d['Year']>=tlim[0])[0]

	rat=d['Province']['Area']['Cruise']/d['Province']['Area']['CC']
	ind=np.where( (rat>0) & (rat<2) )[0]
	#print(np.mean(rat[ind])) 0.71

	ms=2.5;
	cl=np.array([[0.1,0.3,0.6],[0.5,0.85,0],[1,0.5,0],[0.65,0.35,1],[0.6,1,0],[0,0,0]])

	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(22,8));
	ax.plot(d['Year'][iT],d['Province']['Area']['CC'][iT],'-ko',
		 mfc=meta['Graphics']['Colours']['rgb']['Blue Dark'],
		 mec=meta['Graphics']['Colours']['rgb']['Blue Dark'],
		 color=meta['Graphics']['Colours']['rgb']['Blue Dark'],
		 lw=meta['Graphics']['gp']['lw1'],
		 ms=ms,
		 label='Consolidated cutblocks database')
	#ax.plot(d['Year'][iT],d['Province']['Area']['NFD'][iT],'-ko',mfc=cl[0,:],mec=cl[0,:],color=cl[0,:],lw=meta['Graphics']['gp']['lw1'],ms=ms,label='National forestry database')
	#ax.plot(d['Year'][iT],d['Province']['Area']['Cruise'][iT],'--kd',mfc=cl[2,:],mec=cl[2,:],color=cl[2,:],lw=meta['Graphics']['gp']['lw1'],ms=ms,label='Cruise compilation')
	#ax.plot(d['Year'][iT],d['Province']['Area']['NTEMS'][iT],'-.k^',mfc=cl[3,:],mec=cl[3,:],color=cl[3,:],lw=meta['Graphics']['gp']['lw1'],ms=ms,label='NTEMS 2020')
	#ax.plot(d['Year'][iT],d['Area Planted RESULTS'],'-ks',mfc=cl[3,:],mec=cl[3,:],color=cl[3,:],lw=meta['Graphics']['gp']['lw1'],ms=3,label='Area planted (RESULTS)')
	ax.set(xticks=np.arange(1800,2120,5),yticks=np.arange(0,400,25),ylabel='Area harvested (Kha yr$^{-1}$)',xlabel='Time, years',ylim=[0,275],xlim=[tlim[0]-0.5,tlim[1]+0.5])
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=1.5)
	ax.legend(loc='upper left',facecolor=[1,1,1],frameon=False)
	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics'][ 'Print Figure Path'] + '\\AreaHarvestedBC','png',900)

	return

#%%
def Plot_HarvestVolumeTimeSeries(meta,dTMY):

	d=gu.ipickle(meta['Paths']['DB']['Harvest'] + '\\HarvestComp1_ByTime.pkl')

	dT=CalcStatsByTime(meta,dTMY)

	xlim=[2006,np.max(d['Year'])]
	ikp=np.where( (d['Year']>=xlim[0]) & (d['Year']<=xlim[1]) )[0]
	ikp2=np.where( (dT['Time']>=xlim[0]) & (dT['Time']<=xlim[1]) )[0]

	plt.close('all'); fig,ax=plt.subplots(2,1,figsize=gu.cm2inch(22,7.75))
	ax[0].plot(dT['Time'][ikp2],(dT['Sum']['V Logs Abs m3'][ikp2])/1e6,'-bo',label='HBS')
	#ax[0].plot(dT['Time'][ikp2],(dT['Sum']['V NonLog m3'][ikp2])/1e6,'-md',label='')
	ax[0].plot(d['Year'][ikp],d['Province']['Volume']['FLNR'][ikp],'r-',label='Tot. merch. vol. (FLNR)')
	#ax[0,0].plot(H_FLNR['Year'],H_FLNR['Total_harvest_millions_m3'],'r-',label='Tot. merch. vol. (FLNR)')
	ax[0].set(xticks=np.arange(1950,2025,2),ylabel='Merchantable volume (Mm$^3$ yr$^{-1}$)',xlabel='Time, calendar year',xlim=[xlim[0]-0.5,xlim[1]+0.5],ylim=[0,80])

	ax[1].plot(dT['Time'][ikp2],(dT['Sum']['V Logs Waste m3'][ikp2])/1e6,'-bo',label='HBS')
	ax[1].set(xticks=np.arange(1950,2025,2),ylabel='Waste wood volume (Mm$^3$ yr$^{-1}$)',xlabel='Time, calendar year',xlim=[xlim[0]-0.5,xlim[1]+0.5])

	#ax[0].legend(loc='lower left',frameon=False,facecolor=None,edgecolor='w')
	plt.tight_layout()
	gu.axletters(ax,plt,0.025,0.87,FontColor=meta['Graphics']['gp']['cla'],LetterStyle=meta['Graphics']['Modelling']['AxesLetterStyle'],FontWeight=meta['Graphics']['Modelling']['AxesFontWeight'])
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics'][ 'Print Figure Path'] + '\\VolumeHarvest_TS','png',900)

	# Per hectare
	plt.close('all'); fig,ax=plt.subplots(2,1,figsize=gu.cm2inch(22,7.75))
	ax[0].plot(dT['Time'][ikp2],(dT['Mean']['V Logs Abs m3/ha'][ikp2])/1e6,'-bo',label='HBS')
	ax[0].set(xticks=np.arange(1950,2025,2),ylabel='Merchantable volume (Mm$^3$ yr$^{-1}$)',xlabel='Time, calendar year',xlim=[xlim[0]-0.5,xlim[1]+0.5])

	ax[1].plot(dT['Time'][ikp2],(dT['Mean']['V Logs Waste m3/ha'][ikp2])/1e6,'-bo',label='HBS')
	ax[1].set(xticks=np.arange(1950,2025,2),ylabel='Waste wood volume (Mm$^3$ yr$^{-1}$)',xlabel='Time, calendar year',xlim=[xlim[0]-0.5,xlim[1]+0.5])

	#ax[0].legend(loc='lower left',frameon=False,facecolor=None,edgecolor='w')
	plt.tight_layout()
	gu.axletters(ax,plt,0.025,0.87,FontColor=meta['Graphics']['gp']['cla'],LetterStyle=meta['Graphics']['Modelling']['AxesLetterStyle'],FontWeight=meta['Graphics']['Modelling']['AxesFontWeight'])
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics'][ 'Print Figure Path'] + '\\VolumeHarvest_TS','png',900)
	return

#%%
def CalcHarvestVolumeByBGCZ(meta,dTMY,dTM):
	u=np.unique(dTM['BGCZ'])
	d={'m3/ha':{},'Mm3/yr':{},'Time':{},'Mm3/yr Mean':{}}
	for iU in range(u.size):
		if u[iU]=='':
			continue
		ikp=np.where( (dTM['BGCZ']==u[iU]) & (dTM['V Logs m3/ha']>0) & (dTM['V Logs m3/ha']<2000) )[0]
		d['m3/ha'][u[iU]]=np.nanmean(dTM['V Logs m3/ha'][ikp])
		dt=StatsByTime(meta,dTMY,zone=u[iU])
		iT=np.where( (dt['Time']>=2007) & (dt['Time']<=2022) )[0]
		d['Time'][u[iU]]=dt['Time'][iT]
		d['Mm3/yr'][u[iU]]=np.round(dt['Sum']['V Logs m3'][iT]/1e6,decimals=2)
		d['Mm3/yr Mean'][u[iU]]=np.round(np.mean(dt['Sum']['V Logs m3'][iT]/1e6),decimals=2)
	return d

#%%
def Plot_YieldByGrade(meta,dTM):
	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(7.8,6));
	reg='Interior'
	ikp=np.where( (dTM['PLANNED_NET_BLOCK_AREA']>0) & (dTM['Region']==reg) & (dTM['V Logs Abs m3/ha']<3000) )[0]
	labs=['1','2','3','4','5','6','7','8','Z']
	y=np.zeros(9)
	c=0; y[c]=np.nanmean(dTM['V Logs Grade 1 Abs m3/ha'][ikp])
	c=1; y[c]=np.nanmean(dTM['V Logs Grade 2 Abs m3/ha'][ikp])
	c=2; y[c]=np.nanmean(dTM['V Logs Grade 3 Abs m3/ha'][ikp])
	c=3; y[c]=np.nanmean(dTM['V Logs Grade 4 Abs m3/ha'][ikp])
	c=4; y[c]=np.nanmean(dTM['V Logs Grade 5 Abs m3/ha'][ikp])
	c=5; y[c]=np.nanmean(dTM['V Logs Grade 6 Abs m3/ha'][ikp])
	c=6; y[c]=np.nanmean(dTM['V Logs Grade 7 Abs m3/ha'][ikp])
	c=7; y[c]=np.nanmean(dTM['V Logs Grade 8 Abs m3/ha'][ikp])
	c=8; y[c]=np.nanmean(dTM['V Logs Grade Z Abs m3/ha'][ikp])
	ax.barh(np.arange(1,10),y,0.8)
	ax.set(yticks=np.arange(1,10),yticklabels=labs,xlabel='Average yield (m3/ha)',ylabel='Grade')
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
	plt.tight_layout()
	
	#plt.close('all');
	fig,ax=plt.subplots(1,figsize=gu.cm2inch(7.8,6));
	reg='Coast'
	ikp=np.where( (dTM['PLANNED_NET_BLOCK_AREA']>0) & (dTM['Region']==reg) & (dTM['V Logs Abs m3/ha']<3000) )[0]
	labs=['B','C','D','E','F','G','H','I','J','K','L','M','U','W','X','Y','Z']
	y=np.zeros(17)
	c=0; y[c]=np.nanmean(dTM['V Logs Grade B Abs m3/ha'][ikp])
	c=1; y[c]=np.nanmean(dTM['V Logs Grade C Abs m3/ha'][ikp])
	c=2; y[c]=np.nanmean(dTM['V Logs Grade D Abs m3/ha'][ikp])
	c=3; y[c]=np.nanmean(dTM['V Logs Grade E Abs m3/ha'][ikp])
	c=4; y[c]=np.nanmean(dTM['V Logs Grade F Abs m3/ha'][ikp])
	c=5; y[c]=np.nanmean(dTM['V Logs Grade G Abs m3/ha'][ikp])
	c=6; y[c]=np.nanmean(dTM['V Logs Grade H Abs m3/ha'][ikp])
	c=7; y[c]=np.nanmean(dTM['V Logs Grade I Abs m3/ha'][ikp])
	c=8; y[c]=np.nanmean(dTM['V Logs Grade J Abs m3/ha'][ikp])
	c=9; y[c]=np.nanmean(dTM['V Logs Grade K Abs m3/ha'][ikp])
	c=10; y[c]=np.nanmean(dTM['V Logs Grade L Abs m3/ha'][ikp])
	c=11; y[c]=np.nanmean(dTM['V Logs Grade M Abs m3/ha'][ikp])
	c=12; y[c]=np.nanmean(dTM['V Logs Grade U Abs m3/ha'][ikp])
	c=13; y[c]=np.nanmean(dTM['V Logs Grade W Abs m3/ha'][ikp])
	c=14; y[c]=np.nanmean(dTM['V Logs Grade X Abs m3/ha'][ikp])
	c=15; y[c]=np.nanmean(dTM['V Logs Grade Y Abs m3/ha'][ikp])
	c=16; y[c]=np.nanmean(dTM['V Logs Grade Z Abs m3/ha'][ikp])
	ax.barh(np.arange(1,18),y,0.8)
	ax.set(yticks=np.arange(1,18),yticklabels=labs,xlabel='Average yield (m3/ha)',ylabel='Grade')
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
	plt.tight_layout()
	return

#%%
def Calc_HarvestProbability(meta):
	# Import paths and look-up-tables
	meta=u1ha.Init()
	gp=gu.SetGraphics('Manuscript')

	# Import raster datasets

	# Land mask
	zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])

	zSlope=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Terrain\\slope.tif')
	#zElev=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Terrain\elevation.tif')
	zDTF=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Management\\DistanceFromForestryFacility.tif')
	zDTR=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Terrain\\DistanceFromRoads.tif')

	# Harvest mask
	zH=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Disturbances\\VEG_CONSOLIDATED_CUT_BLOCKS_SP_All.tif')

	# Filters and sampling

	ssi=50

	xMask=zRef['Data'][0::ssi,0::ssi].flatten()
	xSlope=zSlope['Data'][0::ssi,0::ssi].flatten().astype(float)
	#xElev=zElev['Data'][0::ssi,0::ssi].flatten().astype(float)
	xDTR=zDTR['Data'][0::ssi,0::ssi].flatten().astype(float)
	xDTF=zDTF['Data'][0::ssi,0::ssi].flatten().astype(float)
	xH=zH['Data'][0::ssi,0::ssi].flatten().astype(float)

	ikp=np.where(xMask>0)[0]

	# Prepare variables
	x=np.column_stack((xSlope[ikp],xDTR[ikp],xDTF[ikp]))
	y=xH[ikp]

	# z-score x variables
	xz,mu,sig=gu.zscore(x)

	# Regression analysis

	# Is it balanced
	print(np.sum(y)/y.size)

	x_train,x_test,y_train,y_test=train_test_split(xz,y,test_size=0.5,random_state=0)

	mod=LogisticRegression()#solver='liblinear',random_state=0)

	mod.fit(x_train,y_train)

	yhat=mod.predict(x_test)

	beta=mod.coef_[0]
	intercept=mod.intercept_[0]
	print(intercept)
	print(beta)

	score=mod.score(x_test,y_test)
	print(score)

	cm=metrics.confusion_matrix(y_test,yhat)
	print(cm)

	# Predict map of annual probability of harvest (%/yr)

	# Duration (Consolidated cutblocks goes from about 1955-2021)
	Ivl=2021-1954

	zSlope['Data']=np.maximum(0,zSlope['Data'])

	zHhat=zRef.copy()
	zHhat['Data']=zHhat['Data'].astype('float32')
	zHhat['Data']=intercept + \
	    (beta[0]*((zSlope['Data'].astype(float)-mu[0])/sig[0])) + \
	    (beta[1]*((zDTR['Data'].astype(float)-mu[1])/sig[1])) + \
	    (beta[2]*((zDTF['Data'].astype(float)-mu[2])/sig[2]))

	zHhat['Data']=(np.exp(zHhat['Data'])/(1+np.exp(zHhat['Data'])))/Ivl*100

	zHhat['Data'][zRef['Data']==0]=0

	# Plot histogram
	plt.close('all')
	plt.hist(zHhat['Data'][0::ssi,0::ssi].flatten())
	
	# Map annual probability of harvest (%/yr)
	ssi=5
	
	plt.close('all')
	plt.matshow(zHhat['Data'][0::ssi,0::ssi],clim=[0,0.5])
	plt.colorbar()
	
	print(np.percentile(100*zHhat['Data'][0::ssi,0::ssi],99.9))
	
	# Save Annal probability of harvest (%/yr)
	sf=1000
	z1=zRef.copy()
	z1['Data']=zHhat['Data']*sf
	z1['Data']=z1['Data'].astype('int32')
	gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\Disturbances\\HarvestProbability.tif')
	
	# THLB at various threshold probabilities
	zLC2=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\\VRI\lc2.tif')
	
	sf=1000
	zPh=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Disturbances\\HarvestProbability.tif')
	zPh['Data']=zPh['Data'].astype('float')/sf
	
	# Treed area
	ind0=np.where( (zLC2['Data']==4) )
	A_treed=ind0[0].size/1e6
	print(A_treed)
	
	# Profile of THLB vs. Ph threshold
	bin=np.arange(0.05,0.3,0.05)
	Abin=np.zeros(bin.size)
	for i in range(bin.size):
	    print(i)
	    ind1=np.where( (zLC2['Data']==4) & (zPh['Data']>bin[i]) )
	    Abin[i]=ind1[0].size/1e6
	
	plt.close('all')
	plt.plot(bin,Abin,'ob-')

	return

#%%
def Plot_TenureType_Frequency(meta,dTM):
	u,N=gu.CountByCategories(dTM['Tenure Type'],'Percent')
	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(14,8));
	ax.barh(np.arange(u.size),N)
	ax.set(yticks=np.arange(u.size),yticklabels=u)
	return

#%%
def Plot_SilvicultureSystemCode_Frequency(meta,dTM):
	u,N=gu.CountByCategories(dTM['SSC'],'Percent')
	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(14,14));
	ax.barh(np.arange(u.size),N)
	ax.set(yticks=np.arange(u.size),yticklabels=u)
	return

#%%
def Plot_RelationshipBtwnGrades7and8andTotalYield(dTM):
	u=np.unique(dTM['District'])
	for iU in range(u.size):
		#dT=hut.StatsByTime(meta,dTMY,District='DCS')
		dT=hut.StatsByTime(meta,dTMY,District=u[iU])
		#dT=hut.StatsByTime(meta,dTMY,District='DMH')
		
		x=(dT['Sum']['V Logs Grade 7 m3']+dT['Sum']['V Logs Grade 8 m3'])/1e6#/dT['Sum']['V Felled m3']
		y=dT['Sum']['V Felled m3']/1e6
		ikp=np.where( (dT['Time']>2010) & (x>0) & (y>0) )[0]
		if ikp.size==0:
			continue
		x=x[ikp]
		y=y[ikp]
		
		plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(7.8,7))
		ax.plot([-10000,10000],[-10000,10000],'k-',lw=2,color=[0.8,0.8,0.8])
		ax.plot(x,y,'o',ms=5,mec='w',mfc='k',mew=0.5)
		rs,txt=gu.GetRegStats(x,y)
		ax.plot(rs['xhat Line'],rs['yhat Line'],'r-')
		ax.text(1,1,rs['txt'],fontsize=7,ha='left')
		ax.text(20,20,'1:1',fontsize=7,ha='center')
		ax.set(xlabel='Grades 7 and 8',ylabel='Total yield',xticks=np.arange(-30,30,0.5),yticks=np.arange(-30,30,0.5),xlim=[0,5],ylim=[0,5])
		ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
		plt.tight_layout()
		gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Harvest\Salvage Logging vs Total Yield\plot_' + u[iU],'png',900)
	return

#%%
def Plot_GrossVsNetVolume(dTM):
	ind=np.where( (dTM['Region']=='Interior') & (dTM['Cruise_Year']>=2014) & (dTM['Cruise_V Net (m3/ha)']>0) & (dTM['Cruise_V Net (m3/ha)']<3000) & (dTM['Cruise_V Gross (m3/ha)']>0) & (dTM['Cruise_V Gross (m3/ha)']<3000) & (dTM['V Logs Abs m3/ha']<2500) )[0]
	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(8,8))
	ax.plot([0,1200],[0,1200],'k-',lw=2,color=[0.8,0.8,0.8])
	x=dTM['Cruise_V Gross (m3/ha)'][ind]
	y=dTM['Cruise_V Net (m3/ha)'][ind]
	ax.plot(x,y,'k.',ms=5,mfc='b',mec='w',mew=0.5)
	rs,txt=gu.GetRegStats(x,y,'No Intercept')
	ax.plot(rs['xhat Line'],rs['yhat Line'],'r--',lw=1.5,color=[0,0,0])
	ax.set(xlabel='Gross volume (m3/ha)',ylabel='Net volume (m3/ha)',xlim=[0,1200],ylim=[0,1200])
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=2)
	ax.legend(loc='lower left',facecolor=[1,1,1],frameon=False);
	plt.tight_layout()
	#gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Harvest\Timber Cruise Info\TimeSeries','png',900)
	return

#%%
def Plot_YieldVsPercentDead(meta,dTM):
	list(dTM.keys())
	
	ind=np.where( (dTM['Region']=='Interior') & (dTM['Cruise_Year']>=2014) & (dTM['Cruise_V Net (m3/ha)']>0) & (dTM['V Logs Abs m3/ha']<1500) )[0]
	x=dTM['Cruise_Pct Dead Net'][ind]
	bw=20; bin=np.arange(0,120,bw)
	
	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(8,8))
	y=dTM['Cruise_V Gross (m3/ha)'][ind]
	N,mu,med,sig,se=gu.discres(x,y,bw,bin)
	ax.plot(bin,mu,'go',ms=5,mfc='w',mec='g',label='Gross volume from cruise')
	y=dTM['Cruise_V Net (m3/ha)'][ind]
	N,mu,med,sig,se=gu.discres(x,y,bw,bin)
	ax.plot(bin,mu,'rs',ms=5,mfc='w',mec='r',label='Net volume from cruise')
	y=dTM['V Logs Abs m3/ha'][ind]
	N,mu,med,sig,se=gu.discres(x,y,bw,bin)
	ax.plot(bin,mu,'b^',ms=4,mfc='w',mec='b',label='Delivered volume from HBS')
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=2)
	ax.set(xlabel='Dead volume component (%)',ylabel='Net volume (m3/ha)')
	ax.legend(loc='lower left',facecolor=[1,1,1],frameon=False);
	plt.tight_layout()
	#gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Harvest\Timber Cruise Info\TimeSeries','png',900)
	return

#%%
def SummarizeAgeAtHarvest(dTM):
	meta=u1ha.Init()
	vList=['lc_comp1_2019','age_vri02','age_vri15','age_vri23','harv_yr_comp1','bgcz','rd','v_vri02']
	z0=u1ha.Import_Raster(meta,[],vList,'Extract Grid')
	tv=np.arange(2002,2025,1)
	d={}
	for dr in meta['LUT']['TA_REGIONAL_DISTRICTS_SVW']['REGIONAL_DISTRICT_NAME'].keys():
		d[dr]={}
		ind=np.where( (z0['lc_comp1_2019']==1) & (z0['harv_yr_comp1']>0) & (z0['rd']==meta['LUT']['TA_REGIONAL_DISTRICTS_SVW']['REGIONAL_DISTRICT_NAME'][dr]) )
		if ind[0].size==0:
			d[dr]=np.nan
			continue
		d0={}
		for k in z0.keys():
			d0[k]=z0[k][ind]
	
		A=np.zeros(tv.size)
		V=np.zeros(tv.size)
		N50=np.zeros(tv.size)
		N200=np.zeros(tv.size)
		for iT in range(tv.size):
			ind=np.where( (d0['lc_comp1_2019']==1) & (d0['harv_yr_comp1']==tv[iT]) & (d0['age_vri02']>0) & (d0['v_vri02']>0) )[0]
			if tv[iT]<2100:
				A[iT]=np.mean(d0['age_vri02'][ind])
				V[iT]=np.mean(d0['v_vri02'][ind])
				N50[iT]=np.sum(d0['age_vri02'][ind]<50)/ind.size*100
				N200[iT]=np.sum(d0['age_vri02'][ind]>200)/ind.size*100
			else:
				A[iT]=np.mean(d0['age_vri15'][ind])
				N50[iT]=np.sum(d0['age_vri15'][ind]<50)/ind.size*100
		rs,txt=gu.GetRegStats(tv,N50)
		d[dr]=rs['B'][1]
		plt.plot(tv,N50,'-o')
		#plt.close('all'); plt.plot(tv,A,'-bo')
		#plt.close('all'); plt.plot(tv,V,'-bo')
		#plt.close('all'); plt.plot(tv,N50,'-bo'); #plt.plot(tv,N200,'-gs')
		#plt.plot(rs['xhat Line'],rs['yhat Line'],'r-')
	return

#%%
def Plot_WasteRatioVsPercentDead(dTM):
	ind=np.where( (dTM['Region']=='Interior') & (dTM['Cruise_Year']>=2014) & (dTM['Cruise_V Net (m3/ha)']>0) & (dTM['V Logs Abs m3/ha']>0) & (dTM['V Logs Abs m3/ha']<2500) )[0]
	x=dTM['Cruise_Pct Dead Net'][ind]
	bw=20; bin=np.arange(0,120,bw)
	
	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(8,8))
	y=dTM['Waste Fraction'][ind]
	#y=dTM['Sawlog Ratio'][ind]
	N,mu,med,sig,se=gu.discres(x,y,bw,bin)
	ax.plot(bin,mu,'-bo',ms=5,mfc=[0.27,0.49,0.74],mec=[0.27,0.49,0.74],color=[0.27,0.49,0.74])
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=2)
	ax.set(xlabel='Dead volume fraction (%)',ylabel='Waste ratio (%)')
	plt.tight_layout()
	#gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Harvest\Timber Cruise Info\TimeSeries','png',900)
	return

#%%
def CalcStumpageByGrade(dTM):
	for gr in gradeL:
		#dTM['Stumpage Grade ' + gr]=dTM['Stump Logs Grade ' + gr + ' $']/(453*dTM['V Logs Grade ' + gr + ' Abs m3']/1000)
		dTM['Stumpage Grade ' + gr]=dTM['Stump Logs Grade ' + gr + ' $']/dTM['V Logs Grade ' + gr + ' m3']
	
	ind=np.where( (dTM['Stumpage Grade 1']>0) & (dTM['V Logs Grade 1 Abs m3']>0) )[0]
	plt.hist(dTM['Stumpage Grade 1'][ind])
	
	ind=np.where( (dTM['Stumpage Grade J']>0) & (dTM['V Logs Grade J Abs m3']>0) )[0]
	plt.hist(dTM['Stumpage Grade J'][ind])
	
	d={}
	for k in dTM.keys():
		try:
			ind=np.where( (dTM['Region']=='Coast') & (dTM[k]>0) & (dTM[k]<1000000) )[0]
			d[k]=np.nanmean(dTM[k][ind])
		except:
			pass
	
	d={}
	for k in dTM.keys():
		try:
			ind=np.where( (dTM['Region']=='Interior') & (dTM[k]>0) & (dTM[k]<1000000) )[0]
			d[k]=np.nanmean(dTM[k][ind])
		except:
			pass
	return d

#%% Add insect and fire damage
# Not done.


#%% Recover TSA (I did this for Shijing, it is really time consuming)
#meta=u1ha.Init()
# zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
# zTM=gis.OpenGeoTiff(r'C:\Data\BC1ha\FTEN_CUT_BLOCK_POLY_SVW\TIMBER_MARK.tif')
# zTSA=gis.OpenGeoTiff(r'C:\Data\BC1ha\FADM_TSA\TSA_NUMBER_DESCRIPTION.tif')
# from scipy import stats

# N_mis=0
# u=np.unique(dTM['TM'])
# dTM['ID TSA']=np.zeros(dTM['TM'].size,dtype='int16')
# for iU in range(u.size):
# 	try:
# 		id=meta['LUT']['FTEN_CUT_BLOCK_POLY_SVW']['TIMBER_MARK'][u[iU]]
# 	except:
# 		print('missing')
# 		N_mis=N_mis+1
# 		continue
# 	ind1=np.where(zTM['Data']==id)
# 	md=stats.mode(zTSA['Data'][ind1])[0]
# 	ind2=np.where(dTM['TM']==u[iU])[0]
# 	dTM['ID TSA'][ind2]=md
# gu.opickle(meta['Paths']['DB']['Harvest'] + '\\HavestSummary_ByTM_New.pkl',dTM)