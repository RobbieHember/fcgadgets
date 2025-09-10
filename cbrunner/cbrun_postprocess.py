#%% Import modules
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import geopandas as gpd
import os
import glob
import openpyxl
import copy
import gc as garc
import time
from shapely.geometry import Point
import scipy.stats as stats
import fcgadgets.cbrunner.cbrun_util as cbu
from fcgadgets.macgyver import util_general as gu
from fcgadgets.macgyver import util_gis as gis
from fcgadgets.macgyver import util_fcs_graphs as ufc
from fcgadgets.hardhat import economics as econ
from fcgadgets.taz import default_stat_event_sim as asm

#%% Calculate model output statistics for GHG variables (from points)

# Notes:

# You can't summarize the statistics in this script - the full set of ensemblers
# must be saved so that the statistics can be calculate for each scenario comparison.

# This also calculates age class distribution stats, but it does not include uncertainty
# (individual ensembles are not saved)

def Calc_MOS_GHG(meta,pNam,**kwargs):

	t0=time.time()
	print('Calculating model output statistics for GHG balance and age class distribution')

	if 'Scenarios' in kwargs.keys():
		scnL=kwargs['Scenarios']
	else:
		scnL=list(np.arange(0,meta[pNam]['Project']['N Scenario'],1).astype(int))

	# Time series of saved results
	tv_saving=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)

	# All non-economic values
	d1=cbu.LoadSingleOutputFile(meta,pNam,0,0,0)
	del d1['C_M_DistByAgent']
	del d1['C_M_DistByAgentPct']
	del d1['Year']

	if 'VariablesToKeep' not in kwargs.keys():
		v2include=list(d1.keys())
	else:
		v2include=kwargs['VariablesToKeep']

	# Add Sawtooth variables
	if meta[pNam]['Project']['Biomass Module']=='Sawtooth':
		v2include=v2include+['N','N_R','N_M_Tot','N_M_Reg','TreeMean_A','TreeMean_H','TreeMean_D','TreeMean_Csw','TreeMean_Csw_G']

	# Loop through scenarios
	for iScn0 in range(len(scnL)):

		iScn=scnL[iScn0]

		# Initialize data structures
		uPS_size=meta[pNam]['Project']['Strata']['Project Type']['Unique ID'].size
		uSS_size=meta[pNam]['Project']['Strata']['Spatial']['Unique ID'].size
		uYS_size=meta[pNam]['Project']['Strata']['Year']['Unique ID'].size
		uOS_size=meta[pNam]['Project']['Strata']['Other']['Unique ID'].size

		# Initialize GHGs
		ghg1={}
		if meta[pNam]['Project']['Scenario Source']=='Script':
			ghg1={}
			for k in v2include:
				ghg1[k]=-99*np.ones((tv_saving.size,meta[pNam]['Project']['N Ensemble'],uPS_size,uSS_size,uYS_size,uOS_size))
		else:
			# Using stands as ensembles for spreadsheet runs
			ghg1={}
			for k in v2include:
				ghg1[k]=-99*np.ones( (tv_saving.size,meta[pNam]['Project']['N Stand'],uPS_size,uSS_size,uYS_size,uOS_size) ,dtype='float64')

		# Age class distribution data
		acd={}
		acd['bwT']=10
		acd['binT']=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End'],acd['bwT'])
		acd['bwA']=5
		acd['binA']=np.arange(0,400,acd['bwA'])
		acd['Data']=np.zeros( (acd['binT'].size,acd['binA'].size,uPS_size,uSS_size,uYS_size,uOS_size) )

		# Land cover / land use
		lc1={}
		for k in meta['LUT']['Derived']['lc_comp1'].keys():
			lc1[k]=np.zeros(tv_saving.size,dtype='float64')
		lu1={}
		for k in meta['LUT']['Derived']['lu_comp1'].keys():
			lu1[k]=np.zeros(tv_saving.size,dtype='float64')

		# Loop through ensembles
		for iEns in range(meta[pNam]['Project']['N Ensemble']):

			# Initialize temporary data structure for full simulation
			ghg0={}
			for k in v2include:
				ghg0[k]=np.nan*np.empty((tv_saving.size,meta[pNam]['Project']['N Stand']))

			# Initialize age class distribution for this ensemble
			Age0=np.zeros( (acd['binT'].size,meta[pNam]['Project']['N Stand']) )

			for iBat in range(meta[pNam]['Project']['N Batch']):
				#print(str(iScn) + ' ' + str(iEns) + ' ' + str(iBat) )

				# Index to batch
				indBat=cbu.IndexToBatch(meta[pNam],iBat)

				# Include specific subset of stands
				if 'StandsToInclude' in kwargs.keys():
					iKeepStands=np.where(kwargs['StandsToInclude'][iEns,indBat]==1)[0]
				else:
					iKeepStands=np.arange(0,meta[pNam]['Project']['Batch Size'][iBat],1,dtype='int32')

				# Populate GHGs
				d1=cbu.LoadSingleOutputFile(meta,pNam,iScn,iEns,iBat)
				for k in v2include:
					if (d1[k].ndim>2):
						# Skip C pools with more than 2 dims
						continue
					try:
						ghg0[k][:,indBat[iKeepStands]]=d1[k][:,iKeepStands].copy()
					except:
						print(k)
						print(d1[k].shape)

				# Populate land cover / land use area
				if (meta[pNam]['Project']['Save List Type']!='Basic') & (meta[pNam]['Project']['Save List Type']!='Basic+Econ'):
					if iEns==0:
						for k in meta['LUT']['Derived']['lc_comp1'].keys():
							ind=np.where(d1['LandCover']==meta['LUT']['Derived']['lc_comp1'][k])
							a=np.zeros(d1['LandCover'].shape,dtype='int16')
							tv0=np.tile(tv_saving,(a.shape[1],1)).T
							a[ind]=tv0[ind]
							idx=gu.IndicesFromUniqueArrayValues(a.flatten())
							for yr in idx:
								if yr==0:
									continue
								iT=np.where(tv_saving==yr)[0]
								lc1[k][iT]=lc1[k][iT]+idx[yr].size
						for k in meta['LUT']['Derived']['lu_comp1'].keys():
							ind=np.where(d1['LandUse']==meta['LUT']['Derived']['lu_comp1'][k])
							a=np.zeros(d1['LandUse'].shape)
							tv0=np.tile(tv_saving,(a.shape[1],1)).T
							a[ind]=tv0[ind]
							idx=gu.IndicesFromUniqueArrayValues(a.flatten())
							for yr in idx:
								if yr==0:
									continue
								iT=np.where(tv_saving==yr)[0]
								lu1[k][iT]=lu1[k][iT]+idx[yr].size

				# Populate age class distribution
				for iT_acd in range(acd['binT'].size):
					indT_acd=np.where(tv_saving==acd['binT'][iT_acd])[0]
					Age0[iT_acd,indBat]=d1['A'][indT_acd,:]

			del d1
			garc.collect()

			# Summarize by project type, region, and time
			if 'StandsToInclude' in kwargs.keys():
				iKeepS=np.where(kwargs['StandsToInclude'][iEns,:]==1)[0]
			else:
				iKeepS=np.arange(0,meta[pNam]['Project']['N Stand'],1,dtype='int32')

			# Track number of stands in each stratum
			nStrat=np.zeros((uPS_size,uSS_size,uYS_size,uOS_size))

			sta=meta[pNam]['Project']['Strata']
			for iPS in range(uPS_size):
				for iSS in range(uSS_size):
					for iYS in range(uYS_size):
						for iOS in range(uOS_size):
							if (sta['Project Type']['Unique CD'][iPS]=='All') & (sta['Spatial']['Unique CD'][iSS]=='All') & (sta['Year']['Unique CD'][iYS]=='All') & (sta['Other']['Unique CD'][iOS]=='All'):
								# Index to all values for "All"
								indStrat=np.where( (iKeepS>-1) )[0]
							elif (sta['Project Type']['Unique CD'][iPS]!='All') & (sta['Spatial']['Unique CD'][iSS]=='All') & (sta['Year']['Unique CD'][iYS]=='All') & (sta['Other']['Unique CD'][iOS]=='All'):
								# Index to specific project type stratum
								indStrat=np.where( (sta['Project Type']['ID'][iKeepS]==sta['Project Type']['Unique ID'][iPS]) )[0]
							elif (sta['Project Type']['Unique CD'][iPS]=='All') & (sta['Spatial']['Unique CD'][iSS]!='All') & (sta['Year']['Unique CD'][iYS]=='All') & (sta['Other']['Unique CD'][iOS]=='All'):
								# Index to specific spatial stratum
								indStrat=np.where( (sta['Spatial']['ID'][iKeepS]==sta['Spatial']['Unique ID'][iSS]) )[0]
							elif (sta['Project Type']['Unique CD'][iPS]=='All') & (sta['Spatial']['Unique CD'][iSS]=='All') & (sta['Year']['Unique CD'][iYS]!='All') & (sta['Other']['Unique CD'][iOS]=='All'):
								# Index to specific year stratum
								indStrat=np.where( (sta['Year']['ID'][iKeepS]==sta['Year']['Unique ID'][iYS]) )[0]
							elif (sta['Project Type']['Unique CD'][iPS]=='All') & (sta['Spatial']['Unique CD'][iSS]=='All') & (sta['Year']['Unique CD'][iYS]=='All') & (sta['Other']['Unique CD'][iOS]!='All'):
								# Index to specific other stratum
								indStrat=np.where( (sta['Other']['ID'][iKeepS]==sta['Other']['Unique ID'][iOS]) )[0]
							elif (sta['Project Type']['Unique CD'][iPS]!='All') & (sta['Spatial']['Unique CD'][iSS]=='All') & (sta['Year']['Unique CD'][iYS]!='All') & (sta['Other']['Unique CD'][iOS]=='All'):
								# Index to specific combination of project type and year stratum
								indStrat=np.where( (sta['Project Type']['ID'][iKeepS]==sta['Project Type']['Unique ID'][iPS]) & (sta['Year']['ID'][iKeepS]==sta['Year']['Unique ID'][iYS]) )[0]
							elif (sta['Project Type']['Unique CD'][iPS]=='All') & (sta['Spatial']['Unique CD'][iSS]!='All') & (sta['Year']['Unique CD'][iYS]!='All') & (sta['Other']['Unique CD'][iOS]=='All'):
								# Index to specific combination of spatial and year stratum
								indStrat=np.where( (sta['Spatial']['ID'][iKeepS]==sta['Spatial']['Unique ID'][iSS]) & (sta['Year']['ID'][iKeepS]==sta['Year']['Unique ID'][iYS]) )[0]
							elif (sta['Project Type']['Unique CD'][iPS]=='All') & (sta['Spatial']['Unique CD'][iSS]=='All') & (sta['Year']['Unique CD'][iYS]!='All') & (sta['Other']['Unique CD'][iOS]!='All'):
								# Index to specific combination of year and other stratum
								indStrat=np.where( (sta['Other']['ID'][iKeepS]==sta['Other']['Unique ID'][iOS]) & (sta['Year']['ID'][iKeepS]==sta['Year']['Unique ID'][iYS]) )[0]
							elif (sta['Project Type']['Unique CD'][iPS]!='All') & (sta['Spatial']['Unique CD'][iSS]=='All') & (sta['Year']['Unique CD'][iYS]=='All') & (sta['Other']['Unique CD'][iOS]!='All'):
								# Index to specific combination of project type and other stratum
								indStrat=np.where( (sta['Project Type']['ID'][iKeepS]==sta['Project Type']['Unique ID'][iPS]) & (sta['Other']['ID'][iKeepS]==sta['Other']['Unique ID'][iOS]) )[0]
							elif (sta['Project Type']['Unique CD'][iPS]!='All') & (sta['Spatial']['Unique CD'][iSS]=='All') & (sta['Year']['Unique CD'][iYS]!='All') & (sta['Other']['Unique CD'][iOS]!='All'):
								# Index to specific combination of project type, year and other stratum
								indStrat=np.where( (sta['Project Type']['ID'][iKeepS]==sta['Project Type']['Unique ID'][iPS]) & (sta['Year']['ID'][iKeepS]==sta['Year']['Unique ID'][iYS]) & (sta['Other']['ID'][iKeepS]==sta['Other']['Unique ID'][iOS]) )[0]
							else:
								continue

							nStrat[iPS,iSS,iYS,iOS]=indStrat.size
	
							# Popuate GHG data
							if meta[pNam]['Project']['Scenario Source']=='Script':
								for k in v2include:
									if k=='A_Harvest':
										ind=np.where(ghg0[k]<=0)
										ghg0[k][ind]=np.nan
										ghg1[k][:,iEns,iPS,iSS,iYS,iOS]=np.nanmean(ghg0[k][:,iKeepS[indStrat]],axis=1)
									else:
										# Don't use nanmean and nansum - too slow
										ghg1[k][:,iEns,iPS,iSS,iYS,iOS]=np.mean(ghg0[k][:,iKeepS[indStrat]],axis=1)
	
							elif meta[pNam]['Project']['Scenario Source']=='Spreadsheet':
								# Treat stands as ensembles
								for iEnsAsStand in range(meta[pNam]['Project']['N Stand']):
									# Popuate GHG data
									for k in v2include:
										ghg1[k][:,iEnsAsStand,iPS,iSS,iYS,iOS]=ghg0[k][:,iEnsAsStand]
	
							else:
								print('Warning, scenario source not recognized!')
	
							# Add ensemble to age class distribution
							for iT_acd in range(acd['binT'].size):
								for iA_acd in range(acd['binA'].size):
									ind_acd=np.where( np.abs(Age0[iT_acd,iKeepS[indStrat]]-acd['binA'][iA_acd])<=acd['bwA']/2 )[0]
									acd['Data'][iT_acd,iA_acd,iPS,iSS,iYS,iOS]=acd['Data'][iT_acd,iA_acd,iPS,iSS,iYS,iOS]+ind_acd.size

		# Save stratum size
		gu.opickle(meta['Paths'][pNam]['Data'] + '\\Outputs\\MOS_ByStrata_NumStands_Scn' + str(iScn+1) + '.pkl',nStrat)

		# Convert GHG to sparse (no need to save a bunch of strata that are empty - combinations of strata can be added if needed)
		ikp=np.where(ghg1['A']!=-99)
		ghg2={}
		for k in ghg1.keys():
			ghg2[k]=ghg1[k][ikp]

		# Save sparse GHG and index
		gu.opickle(meta['Paths'][pNam]['Data'] + '\\Outputs\\MOS_ByStrata_GHGB_Scn' + str(iScn+1) + '.pkl',ghg2)
		gu.opickle(meta['Paths'][pNam]['Data'] + '\\Outputs\\MOS_ByStrata_GHGB_Scn' + str(iScn+1) + '_index.pkl',ikp)

		# Save land cover / land use areas (only for first ensemble because no stochasticity represented and not accommodating query functionality)
		lclu={'LC':lc1,'LU':lu1}
		gu.opickle(meta['Paths'][pNam]['Data'] + '\\Outputs\\MOS_ByStrata_LandCoverAndUse_Scn' + str(iScn+1) + '.pkl',lclu)

		# Save age class distribution
		acd['Data']=acd['Data']/meta[pNam]['Project']['N Ensemble']
		gu.opickle(meta['Paths'][pNam]['Data'] + '\\Outputs\\MOS_ByStrata_AgeClassDist_Scn' + str(iScn+1) + '.pkl',acd)

	t1=time.time()
	print(str((t1-t0)/60) + ' min')
	return

#%% Calculate model output statistics for econ variables (from points)
# Notes: You can't summarize the statistics in this script - the full set of ensemblers
# must be saved so that the statistics can be calculate for each scenario comparison.
def Calc_MOS_Econ(meta,pNam,**kwargs):

	t0=time.time()
	print('Calculating model output statistics for economic variables')

	# Keyword argumenst
	if 'ScenariosToInclude' in kwargs.keys():
		ScenariosToInclude=kwargs['ScenariosToInclude']
	else:
		# Default is to not save individual ensembles
		ScenariosToInclude=np.arange(0,meta[pNam]['Project']['N Scenario'])

	# Time series of saved results
	tv_saving=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)

	# All (excluding 'C_M_DistByAgent' and 'Year','C_M_Dist')
	v2include=['Cost Roads','Cost Harvest Overhead','Cost Harvest Felling and Piling','Cost Harvest Hauling','Cost Harvest Residuals',
			'Cost Milling','Cost Nutrient Management','Cost Planting','Cost Survey','Cost Knockdown','Cost Ripping','Cost PAS Deactivation',
			'Cost Pile Burn','Cost Total','Cost Silviculture Total','Revenue Lumber','Revenue Plywood',
			'Revenue OSB','Revenue MDF','Revenue Paper','Revenue PowerFacilityDom','Revenue PowerGrid','Revenue PelletExport','Revenue PelletDom','Revenue FirewoodDom',
			'Revenue LogExport','Revenue Gross','Revenue Net','Revenue Net Disc','Revenue Gross Disc','Cost Total Disc','Cost Nutrient Management Disc',
			'Cost Silviculture Total Disc','Cost Total_Cumulative','Cost Silviculture Total_Cumulative','Cost Nutrient Management_Cumulative','Cost Total Disc_Cumulative',
			'Cost Silviculture Total Disc_Cumulative','Cost Nutrient Management Disc_Cumulative','Revenue Gross_Cumulative','Revenue Gross Disc_Cumulative','Revenue Net_Cumulative',
			'Revenue Net Disc_Cumulative']

	# Loop through scenarios
	for iS in range(ScenariosToInclude.size):

		iScn=ScenariosToInclude[iS]

		# Initialize data structures
		uPS_size=meta[pNam]['Project']['Strata']['Project Type']['Unique ID'].size
		uSS_size=meta[pNam]['Project']['Strata']['Spatial']['Unique ID'].size
		uYS_size=meta[pNam]['Project']['Strata']['Year']['Unique ID'].size
		uOS_size=meta[pNam]['Project']['Strata']['Other']['Unique ID'].size

# 		Data1={}
# 		for k in v2include:
# 			Data1[k]=-99*np.ones( (tv_saving.size,meta[pNam]['Project']['N Ensemble'],uPS_size,uSS_size,uYS_size,uOS_size) ,dtype='float64')

		Data1={}
		if meta[pNam]['Project']['Scenario Source']=='Script':
			Data1={}
			for k in v2include:
				Data1[k]=-99*np.ones((tv_saving.size,meta[pNam]['Project']['N Ensemble'],uPS_size,uSS_size,uYS_size,uOS_size))
		else:
			# Using stands as ensembles for spreadsheet runs
			Data1={}
			for k in v2include:
				Data1[k]=-99*np.ones( (tv_saving.size,meta[pNam]['Project']['N Stand'],uPS_size,uSS_size,uYS_size,uOS_size) ,dtype='float64')

		# An option to skip economic calculations. The file will still be created, but all zeros
		if meta[pNam]['Project']['Skip Economics']=='Off':

			# Loop through ensembles
			for iEns in range(meta[pNam]['Project']['N Ensemble']):

				# Initialize temporary data structure for full simulation
				Data0={}
				for k in v2include:
					Data0[k]=np.nan*np.empty( (tv_saving.size,meta[pNam]['Project']['N Stand']) ,dtype='float64')

				# Import batches
				for iBat in range(meta[pNam]['Project']['N Batch']):
					# Index to batch
					indBat=cbu.IndexToBatch(meta[pNam],iBat)

					# Include specific subset of stands
					if 'StandsToInclude' in kwargs.keys():
						iKeepStands=np.where(kwargs['StandsToInclude'][iEns,indBat]==1)[0]
					else:
						iKeepStands=np.arange(0,meta[pNam]['Project']['Batch Size'][iBat],1,dtype=int)

					d1=cbu.LoadSingleOutputFile(meta,pNam,iScn,iEns,iBat)

					try:
						ec=gu.ipickle(meta['Paths'][pNam]['Input Scenario'][iScn] + '\\Modified_Events_Ens' + cbu.FixFileNum(iEns) + '_Bat' + cbu.FixFileNum(iBat) + '.pkl')
					except:
						ec=gu.ipickle(meta['Paths'][pNam]['Input Scenario'][iScn] + '\\Events_Ens' + cbu.FixFileNum(iEns) + '_Bat' + cbu.FixFileNum(iBat) + '.pkl')

					# Uncompress event chronology if it has been compressed
					ec=cbu.EventChronologyDecompress(meta,pNam,ec,iScn,iEns,iBat)

					# Land surface attributes
					lsat=gu.ipickle(meta['Paths'][pNam]['Input Scenario'][iScn] + '\\Inventory_Bat' + cbu.FixFileNum(iBat) + '.pkl')

					# Calculate annual cashflow
					econ1=econ.CashflowFromEventChronology(meta,pNam,iScn,iEns,iBat,lsat,ec,d1)

					for k in v2include:
						Data0[k][:,indBat[iKeepStands]]=econ1[k][:,iKeepStands].copy()

				del econ1,lsat,ec
				#garc.collect()

				# Summarize by project type, region, and time
				if 'StandsToInclude' in kwargs.keys():
					iKeepS=np.where(kwargs['StandsToInclude'][iEns,:]==1)[0]
				else:
					iKeepS=np.arange(0,meta[pNam]['Project']['N Stand'],1,dtype='int32')

				sta=meta[pNam]['Project']['Strata']
				for iPS in range(uPS_size):
					for iSS in range(uSS_size):
						for iYS in range(uYS_size):
							for iOS in range(uOS_size):
								if (sta['Project Type']['Unique CD'][iPS]=='All') & (sta['Spatial']['Unique CD'][iSS]=='All') & (sta['Year']['Unique CD'][iYS]=='All') & (sta['Other']['Unique CD'][iOS]=='All'):
									# Index to all values for "All"
									indStrat=np.where( (iKeepS>-1) )[0]
								elif (sta['Project Type']['Unique CD'][iPS]!='All') & (sta['Spatial']['Unique CD'][iSS]=='All') & (sta['Year']['Unique CD'][iYS]=='All') & (sta['Other']['Unique CD'][iOS]=='All'):
									# Index to specific project type stratum
									indStrat=np.where( (sta['Project Type']['ID'][iKeepS]==sta['Project Type']['Unique ID'][iPS]) )[0]
								elif (sta['Project Type']['Unique CD'][iPS]=='All') & (sta['Spatial']['Unique CD'][iSS]!='All') & (sta['Year']['Unique CD'][iYS]=='All') & (sta['Other']['Unique CD'][iOS]=='All'):
									# Index to specific spatial stratum
									indStrat=np.where( (sta['Spatial']['ID'][iKeepS]==sta['Spatial']['Unique ID'][iSS]) )[0]
								elif (sta['Project Type']['Unique CD'][iPS]=='All') & (sta['Spatial']['Unique CD'][iSS]=='All') & (sta['Year']['Unique CD'][iYS]!='All') & (sta['Other']['Unique CD'][iOS]=='All'):
									# Index to specific year stratum
									indStrat=np.where( (sta['Year']['ID'][iKeepS]==sta['Year']['Unique ID'][iYS]) )[0]
								elif (sta['Project Type']['Unique CD'][iPS]=='All') & (sta['Spatial']['Unique CD'][iSS]=='All') & (sta['Year']['Unique CD'][iYS]=='All') & (sta['Other']['Unique CD'][iOS]!='All'):
									# Index to specific other stratum
									indStrat=np.where( (sta['Other']['ID'][iKeepS]==sta['Other']['Unique ID'][iOS]) )[0]
								elif (sta['Project Type']['Unique CD'][iPS]!='All') & (sta['Spatial']['Unique CD'][iSS]=='All') & (sta['Year']['Unique CD'][iYS]!='All') & (sta['Other']['Unique CD'][iOS]=='All'):
									# Index to specific combination of project type and year stratum
									indStrat=np.where( (sta['Project Type']['ID'][iKeepS]==sta['Project Type']['Unique ID'][iPS]) & (sta['Year']['ID'][iKeepS]==sta['Year']['Unique ID'][iYS]) )[0]
								elif (sta['Project Type']['Unique CD'][iPS]=='All') & (sta['Spatial']['Unique CD'][iSS]!='All') & (sta['Year']['Unique CD'][iYS]!='All') & (sta['Other']['Unique CD'][iOS]=='All'):
									# Index to specific combination of spatial and year stratum
									indStrat=np.where( (sta['Spatial']['ID'][iKeepS]==sta['Spatial']['Unique ID'][iSS]) & (sta['Year']['ID'][iKeepS]==sta['Year']['Unique ID'][iYS]) )[0]
								elif (sta['Project Type']['Unique CD'][iPS]=='All') & (sta['Spatial']['Unique CD'][iSS]=='All') & (sta['Year']['Unique CD'][iYS]!='All') & (sta['Other']['Unique CD'][iOS]!='All'):
									# Index to specific combination of year and other stratum
									indStrat=np.where( (sta['Other']['ID'][iKeepS]==sta['Other']['Unique ID'][iOS]) & (sta['Year']['ID'][iKeepS]==sta['Year']['Unique ID'][iYS]) )[0]
								elif (sta['Project Type']['Unique CD'][iPS]!='All') & (sta['Spatial']['Unique CD'][iSS]=='All') & (sta['Year']['Unique CD'][iYS]=='All') & (sta['Other']['Unique CD'][iOS]!='All'):
									# Index to specific combination of project type and other stratum
									indStrat=np.where( (sta['Project Type']['ID'][iKeepS]==sta['Project Type']['Unique ID'][iPS]) & (sta['Other']['ID'][iKeepS]==sta['Other']['Unique ID'][iOS]) )[0]
								elif (sta['Project Type']['Unique CD'][iPS]!='All') & (sta['Spatial']['Unique CD'][iSS]=='All') & (sta['Year']['Unique CD'][iYS]!='All') & (sta['Other']['Unique CD'][iOS]!='All'):
									# Index to specific combination of project type, year and other stratum
									indStrat=np.where( (sta['Project Type']['ID'][iKeepS]==sta['Project Type']['Unique ID'][iPS]) & (sta['Year']['ID'][iKeepS]==sta['Year']['Unique ID'][iYS]) & (sta['Other']['ID'][iKeepS]==sta['Other']['Unique ID'][iOS]) )[0]
								else:
									continue

								#nStrat[iPS,iSS,iYS,iOS]=indStrat.size

								# Popuate data
								if meta[pNam]['Project']['Scenario Source']=='Script':
									for k in v2include:
										# Don't use nanmean and nansum - too slow
										Data1[k][:,iEns,iPS,iSS,iYS,iOS]=np.mean(Data0[k][:,iKeepS[indStrat]],axis=1)
		
								elif meta[pNam]['Project']['Scenario Source']=='Spreadsheet':
									# Treat stands as ensembles
									for iEnsAsStand in range(meta[pNam]['Project']['N Stand']):
										for k in v2include:
											Data1[k][:,iEnsAsStand,iPS,iSS,iYS,iOS]=Data0[k][:,iEnsAsStand]
		
								else:
									print('Warning, scenario source not recognized!')

# 				sta=meta[pNam]['Project']['Strata']
# 				for iPS in range(uPS_size):
# 					for iSS in range(uSS_size):
# 						for iYS in range(uYS_size):
# 							for iOS in range(uOS_size):
# 								if (sta['Project Type']['Unique CD'][iPS]=='All') & (sta['Spatial']['Unique CD'][iSS]=='All') & (sta['Year']['Unique CD'][iYS]=='All') & (sta['Other']['Unique CD'][iOS]=='All'):
# 									# Index to all values for "All"
# 									indStrat=np.where( (iKeepS>-1) )[0]
# 								elif (sta['Project Type']['Unique CD'][iPS]!='All') & (sta['Spatial']['Unique CD'][iSS]=='All') & (sta['Year']['Unique CD'][iYS]=='All') & (sta['Other']['Unique CD'][iOS]=='All'):
# 									# Index to specific project type stratum
# 									indStrat=np.where( (sta['Project Type']['ID'][iKeepS]==sta['Project Type']['Unique ID'][iPS]) )[0]
# 								elif (sta['Project Type']['Unique CD'][iPS]=='All') & (sta['Spatial']['Unique CD'][iSS]!='All') & (sta['Year']['Unique CD'][iYS]=='All') & (sta['Other']['Unique CD'][iOS]=='All'):
# 									# Index to specific spatial stratum
# 									indStrat=np.where( (sta['Spatial']['ID'][iKeepS]==sta['Spatial']['Unique ID'][iSS]) )[0]
# 								elif (sta['Project Type']['Unique CD'][iPS]=='All') & (sta['Spatial']['Unique CD'][iSS]=='All') & (sta['Year']['Unique CD'][iYS]!='All') & (sta['Other']['Unique CD'][iOS]=='All'):
# 									# Index to specific year stratum
# 									indStrat=np.where( (sta['Year']['ID'][iKeepS]==sta['Year']['Unique ID'][iYS]) )[0]
# 								elif (sta['Project Type']['Unique CD'][iPS]=='All') & (sta['Spatial']['Unique CD'][iSS]=='All') & (sta['Year']['Unique CD'][iYS]=='All') & (sta['Other']['Unique CD'][iOS]!='All'):
# 									# Index to specific other stratum
# 									indStrat=np.where( (sta['Other']['ID'][iKeepS]==sta['Other']['Unique ID'][iOS]) )[0]
# 								elif (sta['Project Type']['Unique CD'][iPS]!='All') & (sta['Spatial']['Unique CD'][iSS]=='All') & (sta['Year']['Unique CD'][iYS]!='All') & (sta['Other']['Unique CD'][iOS]=='All'):
# 									# Index to specific combination of project type and year stratum
# 									indStrat=np.where( (sta['Project Type']['ID'][iKeepS]==sta['Project Type']['Unique ID'][iPS]) & (sta['Year']['ID'][iKeepS]==sta['Year']['Unique ID'][iYS]) )[0]
# 								elif (sta['Project Type']['Unique CD'][iPS]=='All') & (sta['Spatial']['Unique CD'][iSS]!='All') & (sta['Year']['Unique CD'][iYS]!='All') & (sta['Other']['Unique CD'][iOS]=='All'):
# 									# Index to specific combination of spatial and year stratum
# 									indStrat=np.where( (sta['Spatial']['ID'][iKeepS]==sta['Spatial']['Unique ID'][iSS]) & (sta['Year']['ID'][iKeepS]==sta['Year']['Unique ID'][iYS]) )[0]
# 								elif (sta['Project Type']['Unique CD'][iPS]!='All') & (sta['Spatial']['Unique CD'][iSS]=='All') & (sta['Year']['Unique CD'][iYS]=='All') & (sta['Other']['Unique CD'][iOS]!='All'):
# 									# Index to specific combination of project type and other stratum
# 									indStrat=np.where( (sta['Project Type']['ID'][iKeepS]==sta['Project Type']['Unique ID'][iPS]) & (sta['Other']['ID'][iKeepS]==sta['Other']['Unique ID'][iOS]) )[0]
# 								elif (sta['Project Type']['Unique CD'][iPS]!='All') & (sta['Spatial']['Unique CD'][iSS]=='All') & (sta['Year']['Unique CD'][iYS]!='All') & (sta['Other']['Unique CD'][iOS]!='All'):
# 									# Index to specific combination of project type, year and other stratum
# 									indStrat=np.where( (sta['Project Type']['ID'][iKeepS]==sta['Project Type']['Unique ID'][iPS]) & (sta['Year']['ID'][iKeepS]==sta['Year']['Unique ID'][iYS]) & (sta['Other']['ID'][iKeepS]==sta['Other']['Unique ID'][iOS]) )[0]
# 								else:
# 									continue

# 								for k in v2include:
# 									# Don't use nanmean and nansum - too slow
# 									Data1[k][:,iEns,iPS,iSS,iYS,iOS]=np.mean(Data0[k][:,iKeepS[indStrat]],axis=1)

		# Convert to sparse (no need to save a bunch of strata that are empty - combinations of strata can be added if needed)
		ikp=np.where(Data1['Cost Total']!=-99)
		Data2={}
		for k in Data1.keys():
			Data2[k]=Data1[k][ikp]

		# Save
		gu.opickle(meta['Paths'][pNam]['Data'] + '\\Outputs\\MOS_ByStrata_Econ_Scn' + str(iScn+1) + '.pkl',Data2)
		gu.opickle(meta['Paths'][pNam]['Data'] + '\\Outputs\\MOS_ByStrata_Econ_Scn' + str(iScn+1) + '_index.pkl',ikp)

	t1=time.time()
	print(str((t1-t0)/60) + ' min')
	return

#%% Calculate model output statistics for area from points
# Notes: Unlike the GHG and Econ scripts, the area can calculate the statistics
# in this scripts (scenario comparisons don't need to include area)
def Calc_MOS_Area(meta,pNam,**kwargs):

	print('Calculating model output statistics for event areas')
	t0=time.time()

	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)

	# Time series of saved results
	tv_saving=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)

	# Full time series (for event chronology)
	tv_full=np.arange(meta[pNam]['Project']['Year Start'],meta[pNam]['Project']['Year End']+1,1)

	# Loop through scenarios
	for iScn in range(meta[pNam]['Project']['N Scenario']):

		# Initialize data structures
		uPS_size=meta[pNam]['Project']['Strata']['Project Type']['Unique ID'].size
		uSS_size=meta[pNam]['Project']['Strata']['Spatial']['Unique ID'].size
		uYS_size=meta[pNam]['Project']['Strata']['Year']['Unique ID'].size
		uOS_size=meta[pNam]['Project']['Strata']['Other']['Unique ID'].size

		Area1={}
		if meta[pNam]['Project']['Scenario Source']=='Script':
			for k in meta['LUT']['Event'].keys():
				Area1[k]=-99*np.ones( (tv_saving.size,meta[pNam]['Project']['N Ensemble'],uPS_size,uSS_size,uYS_size,uOS_size) ,dtype='float64')
		else:
			# Using stands as ensembles for spreadsheet runs
			Area1={}
			for k in meta['LUT']['Event'].keys():
				Area1[k]=-99*np.ones( (tv_saving.size,meta[pNam]['Project']['N Stand'],uPS_size,uSS_size,uYS_size,uOS_size) ,dtype='float64')

		# An option to skip economic calculations. The file will still be created, but all zeros
		#if meta[pNam]['Project']['Skip Economics']=='Off':

		# Loop through ensembles
		for iEns in range(meta[pNam]['Project']['N Ensemble']):

			# Initialize temporary data structure for full simulation
			Area0={}
			for k in meta['LUT']['Event'].keys():
				Area0[k]=np.nan*np.empty( (tv_saving.size,meta[pNam]['Project']['N Stand']) ,dtype='float64')

			Area0_Full={}
			for k in meta['LUT']['Event'].keys():
				Area0_Full[k]=np.zeros( (tv_saving.size,meta[pNam]['Project']['N Stand']) ,dtype='float64')

			# Import batches
			for iBat in range(meta[pNam]['Project']['N Batch']):

				# Index to batches
				indBat=cbu.IndexToBatch(meta[pNam],iBat)

				# Include specific subset of stands
				if 'StandsToInclude' in kwargs.keys():
					iKeepStands=np.where(kwargs['StandsToInclude'][iEns,indBat]==1)[0]
				else:
					iKeepStands=np.arange(0,meta[pNam]['Project']['Batch Size'][iBat],1,dtype=int)

				# Import event chronology
				if (meta[pNam]['Scenario'][iScn]['Harvest Sim Future Status']=='On') | (meta[pNam]['Scenario'][iScn]['Wildfire Sim Pre-obs Status']=='On') | (meta[pNam]['Scenario'][iScn]['Wildfire Sim Future Status']=='On') | (meta[pNam]['Scenario'][iScn]['IBM Sim Pre-obs Status']=='On') | (meta[pNam]['Scenario'][iScn]['IBM Sim Future Status']=='On') | (meta[pNam]['Scenario'][iScn]['Disease Sim Historical Status']=='On') | (meta[pNam]['Scenario'][iScn]['Disease Sim Future Status']=='On') |(meta[pNam]['Scenario'][iScn]['Wind Sim Historical Status']=='On') | (meta[pNam]['Scenario'][iScn]['Wind Sim Future Status']=='On'):
					ec=gu.ipickle(meta['Paths'][pNam]['Input Scenario'][iScn] + '\\Modified_Events_Ens' + cbu.FixFileNum(iEns) + '_Bat' + cbu.FixFileNum(iBat) + '.pkl')
				else:
					ec=gu.ipickle(meta['Paths'][pNam]['Input Scenario'][iScn] + '\\Events_Ens' + cbu.FixFileNum(iEns) + '_Bat' + cbu.FixFileNum(iBat) + '.pkl')

				# Uncompress event chronology if it has been compressed
				ec=cbu.EventChronologyDecompress(meta,pNam,ec,iScn,iEns,iBat)
				for iT in range(tv_full.size):
					indT=np.where(tv==tv_full[iT])[0]
					if indT.size==0:
						continue
					idDist0=np.squeeze(ec['ID Event Type'][iT,:,:])
					for iDistWithinYear in range(idDist0.shape[1]):
						idDist1=idDist0[:,iDistWithinYear]
						uidDist=np.unique(idDist1)
						if np.sum(uidDist)==0:
							continue
						for iU in range(uidDist.size):
							if uidDist[iU]==0:
								continue
							namDist=cbu.lut_n2s(meta['LUT']['Event'],uidDist[iU])[0]
							indDist=np.where(idDist1==uidDist[iU])[0]

							Area0_Full[namDist][indT,indBat[indDist]]=Area0_Full[namDist][indT,indBat[indDist]]+1

				# Pull results for subset of stands
				for k in meta['LUT']['Event'].keys():
					Area0[k][:,indBat[iKeepStands]]=Area0_Full[k][:,indBat[iKeepStands]]

			del ec
			garc.collect()

			# Summarize by project type, region, and time
			if 'StandsToInclude' in kwargs.keys():
				iKeepS=np.where(kwargs['StandsToInclude'][iEns,:]==1)[0]
			else:
				iKeepS=np.arange(0,meta[pNam]['Project']['N Stand'],1,dtype='int32')

			sta=meta[pNam]['Project']['Strata']
			for iPS in range(uPS_size):
				for iSS in range(uSS_size):
					for iYS in range(uYS_size):
						for iOS in range(uOS_size):
							if (sta['Project Type']['Unique CD'][iPS]=='All') & (sta['Spatial']['Unique CD'][iSS]=='All') & (sta['Year']['Unique CD'][iYS]=='All') & (sta['Other']['Unique CD'][iOS]=='All'):
								# Index to all values for "All"
								indStrat=np.where( (iKeepS>-1) )[0]
							elif (sta['Project Type']['Unique CD'][iPS]!='All') & (sta['Spatial']['Unique CD'][iSS]=='All') & (sta['Year']['Unique CD'][iYS]=='All') & (sta['Other']['Unique CD'][iOS]=='All'):
								# Index to specific project type stratum
								indStrat=np.where( (sta['Project Type']['ID'][iKeepS]==sta['Project Type']['Unique ID'][iPS]) )[0]
							elif (sta['Project Type']['Unique CD'][iPS]=='All') & (sta['Spatial']['Unique CD'][iSS]!='All') & (sta['Year']['Unique CD'][iYS]=='All') & (sta['Other']['Unique CD'][iOS]=='All'):
								# Index to specific spatial stratum
								indStrat=np.where( (sta['Spatial']['ID'][iKeepS]==sta['Spatial']['Unique ID'][iSS]) )[0]
							elif (sta['Project Type']['Unique CD'][iPS]=='All') & (sta['Spatial']['Unique CD'][iSS]=='All') & (sta['Year']['Unique CD'][iYS]!='All') & (sta['Other']['Unique CD'][iOS]=='All'):
								# Index to specific year stratum
								indStrat=np.where( (sta['Year']['ID'][iKeepS]==sta['Year']['Unique ID'][iYS]) )[0]
							elif (sta['Project Type']['Unique CD'][iPS]=='All') & (sta['Spatial']['Unique CD'][iSS]=='All') & (sta['Year']['Unique CD'][iYS]=='All') & (sta['Other']['Unique CD'][iOS]!='All'):
								# Index to specific other stratum
								indStrat=np.where( (sta['Other']['ID'][iKeepS]==sta['Other']['Unique ID'][iOS]) )[0]
							elif (sta['Project Type']['Unique CD'][iPS]!='All') & (sta['Spatial']['Unique CD'][iSS]=='All') & (sta['Year']['Unique CD'][iYS]!='All') & (sta['Other']['Unique CD'][iOS]=='All'):
								# Index to specific combination of project type and year stratum
								indStrat=np.where( (sta['Project Type']['ID'][iKeepS]==sta['Project Type']['Unique ID'][iPS]) & (sta['Year']['ID'][iKeepS]==sta['Year']['Unique ID'][iYS]) )[0]
							elif (sta['Project Type']['Unique CD'][iPS]=='All') & (sta['Spatial']['Unique CD'][iSS]!='All') & (sta['Year']['Unique CD'][iYS]!='All') & (sta['Other']['Unique CD'][iOS]=='All'):
								# Index to specific combination of spatial and year stratum
								indStrat=np.where( (sta['Spatial']['ID'][iKeepS]==sta['Spatial']['Unique ID'][iSS]) & (sta['Year']['ID'][iKeepS]==sta['Year']['Unique ID'][iYS]) )[0]
							elif (sta['Project Type']['Unique CD'][iPS]!='All') & (sta['Spatial']['Unique CD'][iSS]=='All') & (sta['Year']['Unique CD'][iYS]=='All') & (sta['Other']['Unique CD'][iOS]!='All'):
								# Index to specific combination of project type and other stratum
								indStrat=np.where( (sta['Project Type']['ID'][iKeepS]==sta['Project Type']['Unique ID'][iPS]) & (sta['Other']['ID'][iKeepS]==sta['Other']['Unique ID'][iOS]) )[0]
							elif (sta['Project Type']['Unique CD'][iPS]!='All') & (sta['Spatial']['Unique CD'][iSS]=='All') & (sta['Year']['Unique CD'][iYS]!='All') & (sta['Other']['Unique CD'][iOS]!='All'):
								# Index to specific combination of project type, year and other stratum
								indStrat=np.where( (sta['Project Type']['ID'][iKeepS]==sta['Project Type']['Unique ID'][iPS]) & (sta['Year']['ID'][iKeepS]==sta['Year']['Unique ID'][iYS]) & (sta['Other']['ID'][iKeepS]==sta['Other']['Unique ID'][iOS]) )[0]
							else:
								continue

							if meta[pNam]['Project']['Scenario Source']=='Script':
								for k in meta['LUT']['Event'].keys():
									Area1[k][:,iEns,iPS,iSS,iYS,iOS]=np.mean(Area0[k][:,iKeepS[indStrat]],axis=1)
							elif meta[pNam]['Project']['Scenario Source']=='Spreadsheet':
								# Treat stands as ensembles
								for k in meta['LUT']['Event'].keys():
									for iEnsAsStand in range(meta[pNam]['Project']['N Stand']):
										Area1[k][:,iEnsAsStand,iPS,iSS,iYS,iOS]=Area0[k][:,iEnsAsStand]
							else:
								pass

		# Convert to sparse (no need to save a bunch of strata that are empty - combinations of strata can be added if needed)
		ikp=np.where(Area1['Wildfire']!=-99)
		Area2={}
		for k in Area1.keys():
			Area2[k]=Area1[k][ikp]
		gu.opickle(meta['Paths'][pNam]['Data'] + '\\Outputs\\MOS_ByStrata_Area_Scn' + str(iScn+1) + '.pkl',Area2)
		gu.opickle(meta['Paths'][pNam]['Data'] + '\\Outputs\\MOS_ByStrata_Area_Scn' + str(iScn+1) + '_index.pkl',ikp)

	t1=time.time()
	print(str((t1-t0)/60) + ' min')
	return

#%%
def Calc_MOS_Area_OLD(meta,pNam,**kwargs):

	print('Calculating model output statistics for event areas')
	t0=time.time()

	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)

	# Time series of saved results
	tv_saving=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)

	# Full time series (for event chronology)
	tv_full=np.arange(meta[pNam]['Project']['Year Start'],meta[pNam]['Project']['Year End']+1,1)

	# Loop through scenarios
	for iScn in range(meta[pNam]['Project']['N Scenario']):

		# Initialize data structures
		uPS_size=meta[pNam]['Project']['Strata']['Project Type']['Unique ID'].size
		uSS_size=meta[pNam]['Project']['Strata']['Spatial']['Unique ID'].size
		uYS_size=meta[pNam]['Project']['Strata']['Year']['Unique ID'].size
		uOS_size=meta[pNam]['Project']['Strata']['Other']['Unique ID'].size

		Area1={}
		if meta[pNam]['Project']['Scenario Source']=='Script':
			for k in meta['LUT']['Event'].keys():
				Area1[k]=np.zeros( (tv_saving.size,meta[pNam]['Project']['N Ensemble'],uPS_size,uSS_size,uYS_size,uOS_size) ,dtype='float64')
		else:
			# Using stands as ensembles for spreadsheet runs
			Area1={}
			for k in meta['LUT']['Event'].keys():
				Area1[k]=np.zeros( (tv_saving.size,meta[pNam]['Project']['N Stand'],uPS_size,uSS_size,uYS_size,uOS_size) ,dtype='float64')

		# An option to skip economic calculations. The file will still be created, but all zeros
		#if meta[pNam]['Project']['Skip Economics']=='Off':

		# Loop through ensembles
		for iEns in range(meta[pNam]['Project']['N Ensemble']):

			# Initialize temporary data structure for full simulation
			Area0={}
			for k in meta['LUT']['Event'].keys():
				Area0[k]=np.nan*np.empty( (tv_saving.size,meta[pNam]['Project']['N Stand']) ,dtype='float64')

			Area0_Full={}
			for k in meta['LUT']['Event'].keys():
				Area0_Full[k]=np.zeros( (tv_saving.size,meta[pNam]['Project']['N Stand']) ,dtype='float64')

			# Import batches
			for iBat in range(meta[pNam]['Project']['N Batch']):

				# Index to batches
				indBat=cbu.IndexToBatch(meta[pNam],iBat)

				# Include specific subset of stands
				if 'StandsToInclude' in kwargs.keys():
					iKeepStands=np.where(kwargs['StandsToInclude'][iEns,indBat]==1)[0]
				else:
					iKeepStands=np.arange(0,meta[pNam]['Project']['Batch Size'][iBat],1,dtype=int)

				# Import event chronology
				if (meta[pNam]['Scenario'][iScn]['Harvest Sim Future Status']=='On') | (meta[pNam]['Scenario'][iScn]['Wildfire Sim Pre-obs Status']=='On') | (meta[pNam]['Scenario'][iScn]['Wildfire Sim Future Status']=='On') | (meta[pNam]['Scenario'][iScn]['IBM Sim Pre-obs Status']=='On') | (meta[pNam]['Scenario'][iScn]['IBM Sim Future Status']=='On') | (meta[pNam]['Scenario'][iScn]['Disease Sim Historical Status']=='On') | (meta[pNam]['Scenario'][iScn]['Disease Sim Future Status']=='On') |(meta[pNam]['Scenario'][iScn]['Wind Sim Historical Status']=='On') | (meta[pNam]['Scenario'][iScn]['Wind Sim Future Status']=='On'):
					ec=gu.ipickle(meta['Paths'][pNam]['Input Scenario'][iScn] + '\\Modified_Events_Ens' + cbu.FixFileNum(iEns) + '_Bat' + cbu.FixFileNum(iBat) + '.pkl')
				else:
					ec=gu.ipickle(meta['Paths'][pNam]['Input Scenario'][iScn] + '\\Events_Ens' + cbu.FixFileNum(iEns) + '_Bat' + cbu.FixFileNum(iBat) + '.pkl')

				# Uncompress event chronology if it has been compressed
				ec=cbu.EventChronologyDecompress(meta,pNam,ec,iScn,iEns,iBat)
				for iT in range(tv_full.size):
					indT=np.where(tv==tv_full[iT])[0]
					if indT.size==0:
						continue
					idDist0=np.squeeze(ec['ID Event Type'][iT,:,:])
					for iDistWithinYear in range(idDist0.shape[1]):
						idDist1=idDist0[:,iDistWithinYear]
						uidDist=np.unique(idDist1)
						if np.sum(uidDist)==0:
							continue
						for iU in range(uidDist.size):
							if uidDist[iU]==0:
								continue
							namDist=cbu.lut_n2s(meta['LUT']['Event'],uidDist[iU])[0]
							indDist=np.where(idDist1==uidDist[iU])[0]

							Area0_Full[namDist][indT,indBat[indDist]]=Area0_Full[namDist][indT,indBat[indDist]]+1

				# Pull results for subset of stands
				for k in meta['LUT']['Event'].keys():
					Area0[k][:,indBat[iKeepStands]]=Area0_Full[k][:,indBat[iKeepStands]]

			del ec
			garc.collect()

			# Summarize by project type, region, and time
			if 'StandsToInclude' in kwargs.keys():
				iKeepS=np.where(kwargs['StandsToInclude'][iEns,:]==1)[0]
			else:
				iKeepS=np.arange(0,meta[pNam]['Project']['N Stand'],1,dtype='int32')

			for iPS in range(uPS_size):
				for iSS in range(uSS_size):
					for iYS in range(uYS_size):
						for iOS in range(uOS_size):
							if (meta[pNam]['Project']['Strata']['Project Type']['Unique CD'][iPS]=='All') & (meta[pNam]['Project']['Strata']['Spatial']['Unique CD'][iSS]=='All') & (meta[pNam]['Project']['Strata']['Year']['Unique CD'][iYS]=='All') & (meta[pNam]['Project']['Strata']['Other']['Unique CD'][iOS]=='All'):
								# Index to all values for "All"
								indStrat=np.where(meta[pNam]['Project']['Strata']['Project Type']['ID'][iKeepS]!=77777)[0]
							elif (meta[pNam]['Project']['Strata']['Project Type']['Unique CD'][iPS]!='All') & (meta[pNam]['Project']['Strata']['Spatial']['Unique CD'][iSS]=='All') & (meta[pNam]['Project']['Strata']['Year']['Unique CD'][iYS]=='All') & (meta[pNam]['Project']['Strata']['Other']['Unique CD'][iOS]=='All'):
								# Index to specific project type stratum
								indStrat=np.where( (meta[pNam]['Project']['Strata']['Project Type']['ID'][iKeepS]==meta[pNam]['Project']['Strata']['Project Type']['Unique ID'][iPS]) )[0]
							elif (meta[pNam]['Project']['Strata']['Project Type']['Unique CD'][iPS]=='All') & (meta[pNam]['Project']['Strata']['Spatial']['Unique CD'][iSS]!='All') & (meta[pNam]['Project']['Strata']['Year']['Unique CD'][iYS]=='All') & (meta[pNam]['Project']['Strata']['Other']['Unique CD'][iOS]=='All'):
								# Index to specific spatial stratum
								indStrat=np.where( (meta[pNam]['Project']['Strata']['Spatial']['ID'][iKeepS]==meta[pNam]['Project']['Strata']['Spatial']['Unique ID'][iSS]) )[0]
							elif (meta[pNam]['Project']['Strata']['Project Type']['Unique CD'][iPS]=='All') & (meta[pNam]['Project']['Strata']['Spatial']['Unique CD'][iSS]=='All') & (meta[pNam]['Project']['Strata']['Year']['Unique CD'][iYS]!='All') & (meta[pNam]['Project']['Strata']['Other']['Unique CD'][iOS]=='All'):
								# Index to specific year stratum
								indStrat=np.where( (meta[pNam]['Project']['Strata']['Year']['ID'][iKeepS]==meta[pNam]['Project']['Strata']['Year']['Unique ID'][iSS]) )[0]
							elif (meta[pNam]['Project']['Strata']['Project Type']['Unique CD'][iPS]=='All') & (meta[pNam]['Project']['Strata']['Spatial']['Unique CD'][iSS]=='All') & (meta[pNam]['Project']['Strata']['Year']['Unique CD'][iYS]=='All') & (meta[pNam]['Project']['Strata']['Other']['Unique CD'][iOS]!='All'):
								# Index to specific other stratum
								indStrat=np.where( (meta[pNam]['Project']['Strata']['Other']['ID'][iKeepS]==meta[pNam]['Project']['Strata']['Other']['Unique ID'][iSS]) )[0]
							else:
								continue

							if meta[pNam]['Project']['Scenario Source']=='Script':
								for k in meta['LUT']['Event'].keys():
									Area1[k][:,iEns,iPS,iSS,iYS,iOS]=np.mean(Area0[k][:,iKeepS[indStrat]],axis=1)
							elif meta[pNam]['Project']['Scenario Source']=='Spreadsheet':
								# Treat stands as ensembles
								for k in meta['LUT']['Event'].keys():
									for iEnsAsStand in range(meta[pNam]['Project']['N Stand']):
										Area1[k][:,iEnsAsStand,iPS,iSS,iYS,iOS]=Area0[k][:,iEnsAsStand]
							else:
								pass

		# Calculate statistics
		Area2={}
		for k in Area1.keys():
			Area2[k]={}
			Area2[k]['Ensemble Mean']=np.mean(Area1[k],axis=1)
			Area2[k]['Ensemble SD']=np.std(Area1[k],axis=1)
			Area2[k]['Ensemble SE']=np.std(Area1[k],axis=1)/np.sqrt(meta[pNam]['Project']['N Ensemble'])
			Area2[k]['Ensemble P005']=np.percentile(Area1[k],0.5,axis=1)
			Area2[k]['Ensemble P025']=np.percentile(Area1[k],2.5,axis=1)
			Area2[k]['Ensemble P250']=np.percentile(Area1[k],25,axis=1)
			Area2[k]['Ensemble P750']=np.percentile(Area1[k],75,axis=1)
			Area2[k]['Ensemble P975']=np.percentile(Area1[k],97.5,axis=1)
			Area2[k]['Ensemble P995']=np.percentile(Area1[k],99.5,axis=1)

		# Save
		gu.opickle(meta['Paths'][pNam]['Data'] + '\\Outputs\\MOS_ByStrata_Area_Scn' + str(iScn+1) + '.pkl',Area2)

		del Area1,Area2
		garc.collect()

	t1=time.time()
	print(str((t1-t0)/60) + ' min')
	return

#%% Calculate mortality by agent (from points)
def Calc_MOS_MortByAgent(meta,pNam,**kwargs):

	t0=time.time()
	print('Calculating model output statistics for mortality by agent')

	# Key word arguments

	if 'StandsToInclude' in kwargs.keys():
		flag_stands_to_include=kwargs['StandsToInclude']
	else:
		flag_stands_to_include=[]

	if 'Scenarios' in kwargs.keys():
		scnL=kwargs['Scenarios']
	else:
		scnL=list(np.arange(0,meta[pNam]['Project']['N Scenario'],1).astype(int))

	# Time series of saved results
	tv_saving=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)

	# All non-economic values
	d0=cbu.LoadSingleOutputFile(meta,pNam,0,0,0)

	# Loop through scenarios
	for iScn0 in range(len(scnL)):

		iScn=scnL[iScn0]

		# Initialize data structures
		uPS_size=meta[pNam]['Project']['Strata']['Project Type']['Unique ID'].size
		uSS_size=meta[pNam]['Project']['Strata']['Spatial']['Unique ID'].size
		uYS_size=meta[pNam]['Project']['Strata']['Year']['Unique ID'].size
		uOS_size=meta[pNam]['Project']['Strata']['Other']['Unique ID'].size

		Data1={}
		for k in d0['C_M_DistByAgent'].keys():
			Data1[k]=np.zeros( (tv_saving.size,meta[pNam]['Project']['N Ensemble'],uPS_size,uSS_size,uYS_size,uOS_size) ,dtype='float32')
			Data1[k]=np.zeros( (tv_saving.size,meta[pNam]['Project']['N Ensemble'],uPS_size,uSS_size,uYS_size,uOS_size) ,dtype='float32')

		# Loop through ensembles
		for iEns in range(meta[pNam]['Project']['N Ensemble']):

			# Initialize temporary data structure for full simulation
			Data0={}
			for k in d0['C_M_DistByAgent'].keys():
				Data0[k]=np.nan*np.empty( (tv_saving.size,meta[pNam]['Project']['N Stand']) ,dtype='float32')

			for iBat in range(meta[pNam]['Project']['N Batch']):

				#print(str(iScn) + ' ' + str(iEns) + ' ' + str(iBat) )

				# Index to batch
				indBat=cbu.IndexToBatch(meta[pNam],iBat)

				# Include specific subset of stands
				if 'StandsToInclude' in kwargs.keys():
					iKeepStands=np.where(kwargs['StandsToInclude'][iEns,indBat]==1)[0]
				else:
					iKeepStands=np.arange(0,meta[pNam]['Project']['Batch Size'][iBat],1,dtype=int)

				d1=cbu.LoadSingleOutputFile(meta,pNam,iScn,iEns,iBat)

				for k in d1['C_M_DistByAgent'].keys():
					idx=d1['C_M_DistByAgent'][k]['idx']
					z=np.zeros( (tv_saving.size,indBat.size),dtype='float32')
					z[idx[0],idx[1]]=meta['Core']['Scale Factor C_M_DistByAgent']*d1['C_M_DistByAgent'][k]['M'].astype('float32')
					Data0[k][:,indBat[iKeepStands]]=z[:,iKeepStands].copy()
					#Data0[k][:,indBat[iKeepStands]]=d1[k][:,iKeepStands].copy()

			del d1
			garc.collect()

			# Summarize by project type, region, and time
			if 'StandsToInclude' in kwargs.keys():
				iKeepS=np.where(kwargs['StandsToInclude'][iEns,:]==1)[0]
			else:
				iKeepS=np.arange(0,meta[pNam]['Project']['N Stand'],1,dtype='int32')

			for iPS in range(uPS_size):
				for iSS in range(uSS_size):
					for iYS in range(uYS_size):
						for iOS in range(uOS_size):
							if (meta[pNam]['Project']['Strata']['Project Type']['Unique CD'][iPS]=='All') & (meta[pNam]['Project']['Strata']['Spatial']['Unique CD'][iSS]=='All') & (meta[pNam]['Project']['Strata']['Year']['Unique CD'][iYS]=='All') & (meta[pNam]['Project']['Strata']['Other']['Unique CD'][iOS]=='All'):
								# Index to all values for "All"
								indStrat=np.where(meta[pNam]['Project']['Strata']['Project Type']['ID'][iKeepS]!=77777)[0]
							elif (meta[pNam]['Project']['Strata']['Project Type']['Unique CD'][iPS]!='All') & (meta[pNam]['Project']['Strata']['Spatial']['Unique CD'][iSS]=='All') & (meta[pNam]['Project']['Strata']['Year']['Unique CD'][iYS]=='All') & (meta[pNam]['Project']['Strata']['Other']['Unique CD'][iOS]=='All'):
								# Index to specific project type stratum
								indStrat=np.where( (meta[pNam]['Project']['Strata']['Project Type']['ID'][iKeepS]==meta[pNam]['Project']['Strata']['Project Type']['Unique ID'][iPS]) )[0]
							elif (meta[pNam]['Project']['Strata']['Project Type']['Unique CD'][iPS]=='All') & (meta[pNam]['Project']['Strata']['Spatial']['Unique CD'][iSS]!='All') & (meta[pNam]['Project']['Strata']['Year']['Unique CD'][iYS]=='All') & (meta[pNam]['Project']['Strata']['Other']['Unique CD'][iOS]=='All'):
								# Index to specific spatial stratum
								indStrat=np.where( (meta[pNam]['Project']['Strata']['Spatial']['ID'][iKeepS]==meta[pNam]['Project']['Strata']['Spatial']['Unique ID'][iSS]) )[0]
							elif (meta[pNam]['Project']['Strata']['Project Type']['Unique CD'][iPS]=='All') & (meta[pNam]['Project']['Strata']['Spatial']['Unique CD'][iSS]=='All') & (meta[pNam]['Project']['Strata']['Year']['Unique CD'][iYS]!='All') & (meta[pNam]['Project']['Strata']['Other']['Unique CD'][iOS]=='All'):
								# Index to specific year stratum
								indStrat=np.where( (meta[pNam]['Project']['Strata']['Year']['ID'][iKeepS]==meta[pNam]['Project']['Strata']['Year']['Unique ID'][iSS]) )[0]
							elif (meta[pNam]['Project']['Strata']['Project Type']['Unique CD'][iPS]=='All') & (meta[pNam]['Project']['Strata']['Spatial']['Unique CD'][iSS]=='All') & (meta[pNam]['Project']['Strata']['Year']['Unique CD'][iYS]=='All') & (meta[pNam]['Project']['Strata']['Other']['Unique CD'][iOS]!='All'):
								# Index to specific year stratum
								indStrat=np.where( (meta[pNam]['Project']['Strata']['Year']['ID'][iKeepS]==meta[pNam]['Project']['Strata']['Year']['Unique ID'][iSS]) )[0]
							else:
								continue

							for k in Data1.keys():
								Data1[k][:,iEns,iPS,iSS,iYS,iOS]=np.mean(Data0[k][:,iKeepS[indStrat]],axis=1)

		# Save
		gu.opickle(meta['Paths'][pNam]['Data'] + '\\Outputs\\MOS_ByStrata_MortByAgent_Scn' + str(iScn+1) + '.pkl',Data1)

	t1=time.time()
	print(str((t1-t0)/60) + ' min')

	return

#%% Calculate MOS for mortality from points

# *** This is not by project type and region because cbrunner is not currently
# set up to save mortality by stand (too much data) ***

def Calc_MOS_Mortality_WhenStandsCombined(meta,pNam):

	print('Calculating model output statistics for mortality')

	# Time series of saved results
	tv_saving=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)

	# Initialize list
	mos=[None]*meta[pNam]['Project']['N Scenario']

	# Loop through scenarios
	for iScn in range(meta[pNam]['Project']['N Scenario']):

		#--------------------------------------------------------------------------
		# Initialize data structures
		#--------------------------------------------------------------------------

		d0=cbu.LoadSingleOutputFile(meta,pNam,0,0,0)

		mos[iScn]={}
		for k in d0['C_M_DistByAgent'].keys():
			mos[iScn][k]=np.zeros((tv_saving.size,meta[pNam]['Project']['N Ensemble']))

		# An option to skip economic calculations. The file will still be created, but all zeros
		#if meta[pNam]['Project']['Skip Economics']=='Off':

		# Loop through ensembles
		for iEns in range(meta[pNam]['Project']['N Ensemble']):

			# Initialize temporary data structure for full simulation
			Data={}
			for k in d0['C_M_DistByAgent'].keys():
				Data[k]=np.zeros(tv_saving.size)

			# Loop through batches and add to Data structure
			for iBat in range(meta[pNam]['Project']['N Batch']):
				d1=cbu.LoadSingleOutputFile(meta,pNam,iScn,iEns,iBat)
				for k in d0['C_M_DistByAgent'].keys():
					Data[k]=Data[k]+d1['C_M_DistByAgent'][k].flatten()

				del d1
				garc.collect()

			# Divide by N stand to get mean
			for k in d0['C_M_DistByAgent'].keys():
				Data[k]=Data[k]/meta[pNam]['Project']['N Stand']

			# Populate ensembles
			for k in d0['C_M_DistByAgent'].keys():
				mos[iScn][k][:,iEns]=Data[k].copy()

		# Average all ensembles
		for k in d0['C_M_DistByAgent'].keys():
			mos[iScn][k]=np.nanmean(mos[iScn][k],axis=1)

	# Save
	gu.opickle(meta['Paths'][pNam]['Data'] + '\\Outputs\\MOS_ByAgent_Mortality.pkl',mos)

	return

#%%
def Calc_AgeClassDistribution(meta,pNam,acd):

	# Time series of saved results
	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)

	acd['Data']=[None]*meta[pNam]['Project']['N Scenario']
	for iScn in range(meta[pNam]['Project']['N Scenario']):
		acd['Data'][iScn]=np.zeros((acd['binT'].size,acd['binA'].size))
		for iEns in range(meta[pNam]['Project']['N Ensemble']):
			A=np.zeros( (acd['binT'].size,meta[pNam]['Project']['N Stand']) )
			for iBat in range(meta[pNam]['Project']['N Batch']):
				indBat=cbu.IndexToBatch(meta[pNam],iBat)
				d1=cbu.LoadSingleOutputFile(meta,pNam,iScn,iEns,iBat)
				for iT in range(acd['binT'].size):
					indT=np.where(tv==acd['binT'][iT])[0]
					A[iT,indBat]=d1['A'][indT,:]
			del d1
			garc.collect()
			for iA in range(acd['binA'].size):
				for iT in range(acd['binT'].size):
					ind=np.where( np.abs(A[iT,:]-acd['binA'][iA])<=acd['bwA']/2 )[0]
					acd['Data'][iScn][iT,iA]=ind.size
		acd['Data'][iScn]=acd['Data'][iScn]/meta[pNam]['Project']['N Ensemble']

	gu.opickle(meta['Paths'][pNam]['Data'] + '\\Outputs\\MOS_Scn' + str(iScn+1) + '_AgeClassDistribution.pkl',acd)

	return acd

#%%
# def Calc_MOS_Map(meta,pNam,iScn,**kwargs):
# REPLACED BY "Calc_MOS_ByStand" (see below)
def Calc_MOS_ByStand(meta,mos,pNam,**kwargs):

	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)

	vL=['A','C_Biomass','C_Litter','C_Soil','C_Soil_OHorizon','C_DeadWood','C_G_Gross',
		'C_G_Net','C_M_Reg','C_M_Dist','C_LF','C_ToMillMerchGreen',
		'C_ToMillNonMerchGreen','C_ToMillMerchDead','C_ToMillNonMerchDead',
		'C_ToPileBurnTot','E_OpenBurning_ForestSector_Total',
		'E_NEB','E_NSB','E_NAB']

	operL=['Mean','Sum']

	d={}
	d['Scenarios']=[None]*len(kwargs['scnL'])
	for iScn in kwargs['scnL']:
		d['Scenarios'][iScn]={}
		for th in kwargs['time_horizon'].keys():
			d['Scenarios'][iScn][th]={}
			for oper in operL:
				d['Scenarios'][iScn][th][oper]={}
				for v in vL:
					d['Scenarios'][iScn][th][oper][v]=np.zeros(meta[pNam]['Project']['N Stand'])

	for iScn in kwargs['scnL']:
		for iEns in range(meta[pNam]['Project']['N Ensemble']):
			print(iEns)
			if iEns>2:
				continue
			for iBat in range(meta[pNam]['Project']['N Batch']):
				indBat=cbu.IndexToBatch(meta[pNam],iBat)
				d0=cbu.LoadSingleOutputFile(meta,pNam,iScn,iEns,iBat)
				for v in vL:
					if v not in d0.keys():
						continue
					if (type(d0[v])==dict) | (v=='Year'):
						continue
					for th in kwargs['time_horizon'].keys():
						iT=np.where( (tv>=kwargs['time_horizon'][th]['t0']) & (tv<=kwargs['time_horizon'][th]['t1']) )[0]
						d['Scenarios'][iScn][th]['Sum'][v][indBat]=d['Scenarios'][iScn][th]['Sum'][v][indBat]+np.sum(d0[v][iT,:],axis=0)
						d['Scenarios'][iScn][th]['Mean'][v][indBat]=d['Scenarios'][iScn][th]['Mean'][v][indBat]+np.mean(d0[v][iT,:],axis=0)

	for v in vL:
		d['Scenarios'][iScn][th]['Sum'][v]=d['Scenarios'][iScn][th]['Sum'][v]/meta[pNam]['Project']['N Ensemble']
		d['Scenarios'][iScn][th]['Mean'][v]=d['Scenarios'][iScn][th]['Mean'][v]/meta[pNam]['Project']['N Ensemble']

	# Calculate the deltas
	d['Delta']={}
	for cNam in mos[pNam]['Delta'].keys():
		d['Delta'][cNam]={}
		for th in kwargs['time_horizon'].keys():
			d['Delta'][cNam][th]={}
			for oper in operL:
				d['Delta'][cNam][th][oper]={}
				for v in vL:
					iB=mos[pNam]['Delta'][cNam]['iB']
					iP=mos[pNam]['Delta'][cNam]['iP']
					d['Delta'][cNam][th][oper][v]=d['Scenarios'][ iP ][th][oper][v]-d['Scenarios'][ iB ][th][oper][v]

	return d

#%% Import scenario data from points
# *** You can't use the nanmean and nanpercentile - way too slow ***
def Import_MOS_ByScnAndStrata_GHGEcon(meta,pNam):

	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)

	nS=[]
	for k in meta[pNam]['Project']['Strata'].keys():
		nS.append(meta[pNam]['Project']['Strata'][k]['Unique ID'].size)

	mos={}
	mos[pNam]={}
	mos[pNam]['Scenarios']=[None]*meta[pNam]['Project']['N Scenario']

	oper='Mean'
	Data=[None]*meta[pNam]['Project']['N Scenario']

	for iScn in range(meta[pNam]['Project']['N Scenario']):

		Data[iScn]={}
		Data[iScn][oper]={}
		mos[pNam]['Scenarios'][iScn]={}
		mos[pNam]['Scenarios'][iScn][oper]={}
		mos[pNam]['Scenarios'][iScn]['nStrat']=gu.ipickle(meta['Paths'][pNam]['Data'] + '\\Outputs\\MOS_ByStrata_NumStands_Scn' + str(iScn+1) + '.pkl')

		# Import GHG
		d0=gu.ipickle(meta['Paths'][pNam]['Data'] + '\\Outputs\\MOS_ByStrata_GHGB_Scn' + str(iScn+1) + '.pkl')
		idx=gu.ipickle(meta['Paths'][pNam]['Data'] + '\\Outputs\\MOS_ByStrata_GHGB_Scn' + str(iScn+1) + '_index.pkl')
		for k in d0.keys():
			if meta[pNam]['Project']['Scenario Source']=='Spreadsheet':
				N_Ens=meta[pNam]['Project']['N Stand']
			else:
				N_Ens=meta[pNam]['Project']['N Ensemble']

			d={}
			d[k]=-99*np.ones((tv.size,N_Ens,nS[0],nS[1],nS[2],nS[3]))
			d[k][idx]=d0[k]

			Data[iScn][oper][k]={}
			Data[iScn][oper][k]['Ensemble Mean']=np.mean(d[k],axis=1)
			#Data[iScn][oper][k]['Ensemble SD']=np.std(d[k],axis=1)
			Data[iScn][oper][k]['Ensemble SE']=np.std(d[k],axis=1)/np.sqrt(meta[pNam]['Project']['N Ensemble'])
			#Data[iScn][oper][k]['Ensemble P005']=np.percentile(d[k],0.5,axis=1)
			Data[iScn][oper][k]['Ensemble P025']=np.percentile(d[k],2.5,axis=1)
			Data[iScn][oper][k]['Ensemble P250']=np.percentile(d[k],25,axis=1)
			Data[iScn][oper][k]['Ensemble P750']=np.percentile(d[k],75,axis=1)
			Data[iScn][oper][k]['Ensemble P975']=np.percentile(d[k],97.5,axis=1)
			#Data[iScn][oper][k]['Ensemble P995']=np.percentile(d[k],99.5,axis=1)

		# Import Economics
		d0=gu.ipickle(meta['Paths'][pNam]['Data'] + '\\Outputs\\MOS_ByStrata_Econ_Scn' + str(iScn+1) + '.pkl')
		idx=gu.ipickle(meta['Paths'][pNam]['Data'] + '\\Outputs\\MOS_ByStrata_Econ_Scn' + str(iScn+1) + '_index.pkl')
		# Convert from sparse to full
		for k in d0.keys():
			d={}
			d[k]=-99*np.ones((tv.size,N_Ens,nS[0],nS[1],nS[2],nS[3]))
			d[k][idx]=d0[k]
			Data[iScn][oper][k]={}
			Data[iScn][oper][k]['Ensemble Mean']=np.mean(d[k],axis=1)
			#Data[iScn][oper][k]['Ensemble SD']=np.std(d[k],axis=1)
			Data[iScn][oper][k]['Ensemble SE']=np.std(d[k],axis=1)/np.sqrt(meta[pNam]['Project']['N Ensemble'])
			#Data[iScn][oper][k]['Ensemble P005']=np.percentile(d[k],0.5,axis=1)
			Data[iScn][oper][k]['Ensemble P025']=np.percentile(d[k],2.5,axis=1)
			Data[iScn][oper][k]['Ensemble P250']=np.percentile(d[k],25,axis=1)
			Data[iScn][oper][k]['Ensemble P750']=np.percentile(d[k],75,axis=1)
			Data[iScn][oper][k]['Ensemble P975']=np.percentile(d[k],97.5,axis=1)
			#Data[iScn][oper][k]['Ensemble P995']=np.percentile(d[k],99.5,axis=1)

		# Import area of events
		d0=gu.ipickle(meta['Paths'][pNam]['Data'] + '\\Outputs\\MOS_ByStrata_Area_Scn' + str(iScn+1) + '.pkl')
		idx=gu.ipickle(meta['Paths'][pNam]['Data'] + '\\Outputs\\MOS_ByStrata_Area_Scn' + str(iScn+1) + '_index.pkl')
		for k in d0.keys():
			d={}
			d[k]=-99*np.ones((tv.size,N_Ens,nS[0],nS[1],nS[2],nS[3]))
			d[k][idx]=d0[k]
			Data[iScn][oper]['Area_' + k]={}
			Data[iScn][oper]['Area_' + k]['Ensemble Mean']=np.mean(d[k],axis=1)

		# Mortality by Agent
		d=gu.ipickle(meta['Paths'][pNam]['Data'] + '\\Outputs\\MOS_ByStrata_MortbyAgent_Scn' + str(iScn+1) + '.pkl')
		for k in d.keys():
			nam=cbu.lut_n2s(meta['LUT']['Event'],k)[0]
			Data[iScn][oper]['C_M_' + nam]={}
			Data[iScn][oper]['C_M_' + nam]['Ensemble Mean']=np.mean(d[k],axis=1)
			#Data[iScn][oper]['C_M_' + nam]['Ensemble SD']=np.std(d[k],axis=1)
			Data[iScn][oper]['C_M_' + nam]['Ensemble SE']=np.std(d[k],axis=1)/np.sqrt(meta[pNam]['Project']['N Ensemble'])
			#Data[iScn][oper]['C_M_' + nam]['Ensemble P005']=np.percentile(d[k],0.5,axis=1)
			Data[iScn][oper]['C_M_' + nam]['Ensemble P025']=np.percentile(d[k],2.5,axis=1)
			Data[iScn][oper]['C_M_' + nam]['Ensemble P250']=np.percentile(d[k],25,axis=1)
			Data[iScn][oper]['C_M_' + nam]['Ensemble P750']=np.percentile(d[k],75,axis=1)
			Data[iScn][oper]['C_M_' + nam]['Ensemble P975']=np.percentile(d[k],97.5,axis=1)
			#Data[iScn][oper]['C_M_' + nam]['Ensemble P995']=np.percentile(d[k],99.5,axis=1)

		# Convert to sparse
		mos[pNam]['MOS Index']=np.where(Data[iScn][oper]['A']['Ensemble Mean']!=-99)
		for v in Data[iScn][oper].keys():
			mos[pNam]['Scenarios'][iScn][oper][v]={}
			for stat in Data[iScn][oper][v].keys():
				mos[pNam]['Scenarios'][iScn][oper][v][stat]=Data[iScn][oper][v][stat][ mos[pNam]['MOS Index'] ]

		# Sort variables
		myKeys=list(mos[pNam]['Scenarios'][iScn][oper].keys())
		myKeys.sort()
		mos[pNam]['Scenarios'][iScn][oper]={i:mos[pNam]['Scenarios'][iScn][oper][i] for i in myKeys}

		# Calculate sum and covert to sparse and sort
		mos[pNam]['Scenarios'][iScn]['Sum']={}
		for v in Data[iScn][oper].keys():
			mos[pNam]['Scenarios'][iScn]['Sum'][v]={}
			for stat in Data[iScn][oper][v].keys():
				y=Data[iScn][oper][v][stat]*mos[pNam]['Scenarios'][iScn]['nStrat']*meta[pNam]['Project']['AEF']
				mos[pNam]['Scenarios'][iScn]['Sum'][v][stat]=y[ mos[pNam]['MOS Index'] ]

		myKeys=list(mos[pNam]['Scenarios'][iScn]['Sum'].keys())
		myKeys.sort()
		mos[pNam]['Scenarios'][iScn]['Sum']={i:mos[pNam]['Scenarios'][iScn]['Sum'][i] for i in myKeys}

	return mos

#%% Import future scenario comparison (from points)
def Import_MOS_ByScnComparisonAndStrata(meta,pNam,mos):

	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)

	nS=[]
	for k in meta[pNam]['Project']['Strata'].keys():
		nS.append(meta[pNam]['Project']['Strata'][k]['Unique ID'].size)

	for cNam in mos[pNam]['Delta']:

		nStrat=mos[pNam]['Scenarios'][ mos[pNam]['Delta'][cNam]['iB'] ]['nStrat']
		nStrat=np.tile(nStrat,(tv.size,1,1,1,1))

		if meta[pNam]['Project']['Scenario Source']=='Spreadsheet':
			N_Ens=meta[pNam]['Project']['N Stand']
		else:
			N_Ens=meta[pNam]['Project']['N Ensemble']

		# Initialize
		mos[pNam]['Delta'][cNam]['Data']={}
		mos[pNam]['Delta'][cNam]['Data']['Mean']={}
		mos[pNam]['Delta'][cNam]['Data']['Sum']={}

		dC={}
		dC['Mean']={}

		# Add GHG emissions
		dB0=gu.ipickle(meta['Paths'][pNam]['Data'] + '\\Outputs\\MOS_ByStrata_GHGB_Scn' + str(mos[pNam]['Delta'][cNam]['iB']+1) + '.pkl')
		idx=gu.ipickle(meta['Paths'][pNam]['Data'] + '\\Outputs\\MOS_ByStrata_GHGB_Scn' + str(mos[pNam]['Delta'][cNam]['iB']+1) + '_index.pkl')

		dB={}
		for k in dB0.keys():
			dB[k]=np.zeros((tv.size,N_Ens,nS[0],nS[1],nS[2],nS[3]))
			dB[k][idx]=dB0[k]

		dP0=gu.ipickle(meta['Paths'][pNam]['Data'] + '\\Outputs\\MOS_ByStrata_GHGB_Scn' + str(mos[pNam]['Delta'][cNam]['iP']+1) + '.pkl')
		idx=gu.ipickle(meta['Paths'][pNam]['Data'] + '\\Outputs\\MOS_ByStrata_GHGB_Scn' + str(mos[pNam]['Delta'][cNam]['iP']+1) + '_index.pkl')
		dP={}
		for k in dP0.keys():
			dP[k]=np.zeros((tv.size,N_Ens,nS[0],nS[1],nS[2],nS[3]))
			dP[k][idx]=dP0[k]

		for k in dB.keys():
			dC['Mean'][k]={}
			dC['Mean'][k]['Ensemble Mean']=np.mean(dP[k]-dB[k],axis=1)
			#dC['Mean'][k]['Ensemble SD']=np.std(dP[k]-dB[k],axis=1)
			dC['Mean'][k]['Ensemble SE']=np.std(dP[k]-dB[k],axis=1)/np.sqrt(meta[pNam]['Project']['N Ensemble'])
			#dC['Mean'][k]['Ensemble P005']=np.percentile(dP[k]-dB[k],0.5,axis=1)
			dC['Mean'][k]['Ensemble P025']=np.percentile(dP[k]-dB[k],2.5,axis=1)
			dC['Mean'][k]['Ensemble P250']=np.percentile(dP[k]-dB[k],25,axis=1)
			dC['Mean'][k]['Ensemble P750']=np.percentile(dP[k]-dB[k],75,axis=1)
			dC['Mean'][k]['Ensemble P975']=np.percentile(dP[k]-dB[k],97.5,axis=1)
			#dC['Mean'][k]['Ensemble P995']=np.percentile(dP[k]-dB[k],99.5,axis=1)

		# Add Economics
		dB0=gu.ipickle(meta['Paths'][pNam]['Data'] + '\\Outputs\\MOS_ByStrata_Econ_Scn' + str(mos[pNam]['Delta'][cNam]['iB']+1) + '.pkl')
		idx=gu.ipickle(meta['Paths'][pNam]['Data'] + '\\Outputs\\MOS_ByStrata_Econ_Scn' + str(mos[pNam]['Delta'][cNam]['iB']+1) + '_index.pkl')
		dB={}
		for k in dB0.keys():
			dB[k]=np.zeros((tv.size,N_Ens,nS[0],nS[1],nS[2],nS[3]))
			dB[k][idx]=dB0[k]

		dP0=gu.ipickle(meta['Paths'][pNam]['Data'] + '\\Outputs\\MOS_ByStrata_Econ_Scn' + str(mos[pNam]['Delta'][cNam]['iP']+1) + '.pkl')
		idx=gu.ipickle(meta['Paths'][pNam]['Data'] + '\\Outputs\\MOS_ByStrata_Econ_Scn' + str(mos[pNam]['Delta'][cNam]['iP']+1) + '_index.pkl')
		dP={}
		for k in dP0.keys():
			dP[k]=np.zeros((tv.size,N_Ens,nS[0],nS[1],nS[2],nS[3]))
			dP[k][idx]=dP0[k]

		for k in dB.keys():
			dC['Mean'][k]={}
			dC['Mean'][k]['Ensemble Mean']=np.mean(dP[k]-dB[k],axis=1)
			#dC['Mean'][k]['Ensemble SD']=np.std(dP[k]-dB[k],axis=1)
			dC['Mean'][k]['Ensemble SE']=np.std(dP[k]-dB[k],axis=1)/np.sqrt(meta[pNam]['Project']['N Ensemble'])
			#dC['Mean'][k]['Ensemble P005']=np.percentile(dP[k]-dB[k],0.5,axis=1)
			dC['Mean'][k]['Ensemble P025']=np.percentile(dP[k]-dB[k],2.5,axis=1)
			dC['Mean'][k]['Ensemble P250']=np.percentile(dP[k]-dB[k],25,axis=1)
			dC['Mean'][k]['Ensemble P750']=np.percentile(dP[k]-dB[k],75,axis=1)
			dC['Mean'][k]['Ensemble P975']=np.percentile(dP[k]-dB[k],97.5,axis=1)
			#dC['Mean'][k]['Ensemble P995']=np.percentile(dP[k]-dB[k],99.5,axis=1)

		# Convert mean to sparse and sort
		for v in dC['Mean'].keys():
			mos[pNam]['Delta'][cNam]['Data']['Mean'][v]={}
			for stat in dC['Mean'][v].keys():
				mos[pNam]['Delta'][cNam]['Data']['Mean'][v][stat]=dC['Mean'][v][stat][ mos[pNam]['MOS Index'] ]

		myKeys=list(mos[pNam]['Delta'][cNam]['Data']['Mean'].keys())
		myKeys.sort()
		mos[pNam]['Delta'][cNam]['Data']['Mean']={i:mos[pNam]['Delta'][cNam]['Data']['Mean'][i] for i in myKeys}

		# Calculate sum and covert to sparse and sort
		for v in dC['Mean'].keys():
			mos[pNam]['Delta'][cNam]['Data']['Sum'][v]={}
			for stat in dC['Mean'][v].keys():
				y=dC['Mean'][v][stat]*nStrat*meta[pNam]['Project']['AEF']
				mos[pNam]['Delta'][cNam]['Data']['Sum'][v][stat]=y[ mos[pNam]['MOS Index'] ]

		myKeys=list(mos[pNam]['Delta'][cNam]['Data']['Sum'].keys())
		myKeys.sort()
		mos[pNam]['Delta'][cNam]['Data']['Sum']={i:mos[pNam]['Delta'][cNam]['Data']['Sum'][i] for i in myKeys}

	return mos

#%%
def GetMosScnVar(meta,pNam,mos,iScn,v,iPS,iSS,iYS,iOS):
	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	nS=[]
	for k in meta[pNam]['Project']['Strata'].keys():
		nS.append(meta[pNam]['Project']['Strata'][k]['Unique ID'].size)
	y={}
	y['Mean']={}
	y['Sum']={}
	for stat in mos[pNam]['Scenarios'][iScn]['Mean'][v].keys():
		y['Mean'][stat]=-99*np.ones((tv.size,nS[0],nS[1],nS[2],nS[3]))
		y['Mean'][stat][ mos[pNam]['MOS Index'] ]=mos[pNam]['Scenarios'][iScn]['Mean'][v][stat]
		y['Mean'][stat]=y['Mean'][stat][:,iPS,iSS,iYS,iOS]
		nStrat=mos[pNam]['Scenarios'][iScn]['nStrat']
		y['Sum'][stat]=y['Mean'][stat]*nStrat[iPS,iSS,iYS,iOS]*meta[pNam]['Project']['AEF']
	return y

#%%
def GetMosDeltaVar(meta,pNam,mos,cNam,v,iPS,iSS,iYS,iOS):
	# *** Sum has already applied AEF during import of the Deltas ***
	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	nS=[]
	for k in meta[pNam]['Project']['Strata'].keys():
		nS.append(meta[pNam]['Project']['Strata'][k]['Unique ID'].size)
	y={}
	for oper in ['Mean','Sum']:
		y[oper]={}
		for stat in mos[pNam]['Delta'][cNam]['Data'][oper][v].keys():
			y[oper][stat]=-99*np.ones((tv.size,nS[0],nS[1],nS[2],nS[3]))
			y[oper][stat][ mos[pNam]['MOS Index'] ]=mos[pNam]['Delta'][cNam]['Data'][oper][v][stat]
			y[oper][stat]=y[oper][stat][:,iPS,iSS,iYS,iOS]

	return y

#%% Import model output statistics for area (from points)
# Area isn't considered in scenario comparisons so the model output statistics
# have already been calculated and just need to be imported.
def Import_MOS_ByScnAndStrata_Area(meta,pNam,mos):
	for iScn in range(meta[pNam]['Project']['N Scenario']):
		d=gu.ipickle(meta['Paths'][pNam]['Data'] + '\\Outputs\\MOS_ByStrata_Area_Scn' + str(iScn+1) + '.pkl')
		for v in d.keys():
			mos[pNam]['Scenarios'][iScn]['Mean']['Area_' + v]=d[v]
	return mos

def Import_MOS_ByScnAndStrata_Area_OLD(meta,pNam,mos):
	for iScn in range(meta[pNam]['Project']['N Scenario']):
		d=gu.ipickle(meta['Paths'][pNam]['Data'] + '\\Outputs\\MOS_ByStrata_Area_Scn' + str(iScn+1) + '.pkl')
		for v in d.keys():
			mos[pNam]['Scenarios'][iScn]['Mean']['Area_' + v]=d[v]

		# Import number of stands per stratum
		nStrat=gu.ipickle(meta['Paths'][pNam]['Data'] + '\\Outputs\\MOS_ByStrata_NumStands_Scn' + str(iScn+1) + '.pkl')
		sPS=meta[pNam]['Project']['Strata']['Project Type']['Unique ID'].size
		sSS=meta[pNam]['Project']['Strata']['Spatial']['Unique ID'].size
		sYS=meta[pNam]['Project']['Strata']['Year']['Unique ID'].size

		for v in d.keys():
			mos[pNam]['Scenarios'][iScn]['Sum']['Area_' + v]={}
			for st in mos[pNam]['Scenarios'][iScn]['Mean']['Area_' + v].keys():
				mos[pNam]['Scenarios'][iScn]['Sum']['Area_' + v][st]=0*mos[pNam]['Scenarios'][iScn]['Mean']['Area_' + v][st]
				for iPS in range(sPS):
					for iSS in range(sSS):
						for iYS in range(sYS):
							mos[pNam]['Scenarios'][iScn]['Sum']['Area_' + v][st][:,iPS,iSS,iYS]=mos[pNam]['Scenarios'][iScn]['Mean']['Area_' + v][st][:,iPS,iSS,iYS]*nStrat[iPS,iSS,iYS]*meta[pNam]['Project']['AEF']

	return mos

#%%
def CalcMosByBGC(meta,pNam,lsat):

	tvSaved=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	vmL=['A','A_Harvest','C_Biomass','C_StemMerch','C_StemNonMerch','C_Foliage','C_Branch',
		'C_Bark','C_Root','C_DeadWood','C_G_Gross','C_G_Net_Reg','C_G_Net','C_LF',
		'C_M','C_M_Reg','C_M_Dist','C_M_Nat','C_M_Harv','C_Soil','C_Soil_OHorizon','C_NPP',
		'V_ToMill_MerchTotal']
	
	d=[None]*meta[pNam]['Project']['N Scenario']
	for iScn in range(meta[pNam]['Project']['N Scenario']):
		d[iScn]={}
		for zone in meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE'].keys():
			d[iScn][zone]={}
			for v in vmL:
				d[iScn][zone][v]=np.zeros( (tvSaved.size,) )
				
		for iEns in range(meta[pNam]['Project']['N Ensemble']):
			for iBat in range(meta[pNam]['Project']['N Batch']):
				indBat=cbu.IndexToBatch(meta[pNam],iBat)
				idxBat=gu.IndicesFromUniqueArrayValues(lsat['ID_BGCZ'][indBat])
				d1=cbu.LoadSingleOutputFile(meta,pNam,iScn,iEns,iBat)
				for i in idxBat.keys():
					zone=cbu.lut_n2s(meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE'],i)[0]
					for v in vmL:
						d[iScn][zone][v]=d[iScn][zone][v]+np.sum(d1[v][:,idxBat[i]],axis=1)

		# Convert from sum to average
		idx=gu.IndicesFromUniqueArrayValues(lsat['ID_BGCZ'])
		for v in vmL:
			if v=='V_ToMill_MerchTotal':
				continue
			for zone in meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE'].keys():
				id=meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE'][zone]
				if id in idx.keys():
					d[iScn][zone][v]=d[iScn][zone][v]/idx[id].size

		# Divide by number of ensembles
		for v in vmL:
			for zone in meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE'].keys():
				d[iScn][zone][v]=d[iScn][zone][v]/meta[pNam]['Project']['N Ensemble']
	gu.opickle(meta['Paths'][pNam]['Data'] + '\\Outputs\\MosByBGC.pkl',d)
	return d

#%%
def CalcMosByAgeAndRegion(meta,pNam,lsat):
	bw=25; bin=np.arange(bw,250+bw,bw)
	tvSaved=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)

	vmL=['C_Biomass','C_StemMerch','C_StemNonMerch','C_Foliage','C_Branch',
	'C_Bark','C_Root','C_DeadWood','C_G_Gross','C_G_Net_Reg','C_G_Net_Reg',
	'C_M','C_M_Reg','C_M_Dist','C_Soil','C_Soil_OHorizon','C_M_Harv','C_M_Nat','C_NPP']
	
	d=[None]*meta[pNam]['Project']['N Scenario']
	for iScn in range(meta[pNam]['Project']['N Scenario']):
		d[iScn]={}
		for reg in meta['LUT']['Region'].keys():
			d[iScn][reg]={}
			for v in vmL: 
				d[iScn][reg][v]={}
				d[iScn][reg][v]['N']=np.zeros( bin.size )
				d[iScn][reg][v]['mu']=np.zeros( bin.size )
				d[iScn][reg][v]['se']=np.zeros( bin.size )
	
		for iEns in range(meta[pNam]['Project']['N Ensemble']):
			for iBat in range(meta[pNam]['Project']['N Batch']):
				indBat=cbu.IndexToBatch(meta[pNam],iBat)
				idxBat=gu.IndicesFromUniqueArrayValues(lsat['Region Code'][indBat])
				d1=cbu.LoadSingleOutputFile(meta,pNam,iScn,iEns,iBat)
				for i in idxBat.keys():
					reg=cbu.lut_n2s(meta['LUT']['Region'],i)[0]
					for v in vmL:
						x=d1['A'][:,idxBat[i]].flatten()
						y=d1[v][:,idxBat[i]].flatten()
						N,mu,med,sig,se=gu.discres(x,y,bw,bin)
						d[iScn][reg][v]['N']=d[iScn][reg][v]['N']+N
						d[iScn][reg][v]['mu']=d[iScn][reg][v]['mu']+mu
						d[iScn][reg][v]['se']=d[iScn][reg][v]['se']+se

			# Divide by number of batches
			for reg in d[iScn].keys():
				for v in d[iScn][reg].keys():
					d[iScn][reg][v]['mu']=d[iScn][reg][v]['mu']/meta[pNam]['Project']['N Batch']
					d[iScn][reg][v]['se']=d[iScn][reg][v]['se']/meta[pNam]['Project']['N Batch']

		# Divide by number of ensembles
		for reg in d[iScn].keys():
			for v in d[iScn][reg].keys():
				d[iScn][reg][v]['N']=d[iScn][reg][v]['N']/meta[pNam]['Project']['N Ensemble']
				d[iScn][reg][v]['mu']=d[iScn][reg][v]['mu']/meta[pNam]['Project']['N Ensemble']
				d[iScn][reg][v]['se']=d[iScn][reg][v]['se']/meta[pNam]['Project']['N Ensemble']
	gu.opickle(meta['Paths'][pNam]['Data'] + '\\Outputs\\MosByAgeAndRegion.pkl',d)
	return d

#%% Combine completed and future
# projects must have the same time period saved
def CombineProjectMOSs(meta,mos,pNamNew,pNamL):
	mos[pNamNew]=copy.deepcopy(mos[pNamL[0]])
	for iScn in range(meta[pNamL[0]]['Project']['N Scenario']):
		for k1 in mos[pNamNew]['Scenarios'][iScn]['Sum'].keys():
			for k2 in mos[pNamNew]['Scenarios'][iScn]['Sum'][k1].keys():
				mos[pNamNew]['Scenarios'][iScn]['Sum'][k1][k2]=mos[pNamNew]['Scenarios'][iScn]['Sum'][k1][k2]*meta[pNamL[0]]['Project']['AEF']+mos[pNamL[1]]['Scenarios'][iScn]['Sum'][k1][k2]*meta[pNamL[1]]['Project']['AEF']
				mos[pNamNew]['Scenarios'][iScn]['Mean'][k1][k2]=mos[pNamNew]['Scenarios'][iScn]['Mean'][k1][k2]+mos[pNamL[1]]['Scenarios'][iScn]['Mean'][k1][k2]
	for cNam in mos[pNamNew]['Delta'].keys():
		for k1 in mos[pNamNew]['Delta'][cNam]['ByStrata']['Sum'].keys():
			for k2 in mos[pNamNew]['Delta'][cNam]['ByStrata']['Sum'][k1].keys():
				mos[pNamNew]['Delta'][cNam]['ByStrata']['Sum'][k1][k2]=mos[pNamNew]['Delta'][cNam]['ByStrata']['Sum'][k1][k2]*meta[pNamL[0]]['Project']['AEF']+mos[pNamL[1]]['Delta'][cNam]['ByStrata']['Sum'][k1][k2]*meta[pNamL[1]]['Project']['AEF']
				mos[pNamNew]['Delta'][cNam]['ByStrata']['Mean'][k1][k2]=mos[pNamNew]['Delta'][cNam]['ByStrata']['Mean'][k1][k2]+mos[pNamL[1]]['Delta'][cNam]['ByStrata']['Mean'][k1][k2]
	return mos

#%% Unpack ensemble stats from MOS

def UnpackEnsembleStatsFromMos(meta,pNam,mos):

	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)

	#--------------------------------------------------------------------------
	# Unpack contents for easy use
	#--------------------------------------------------------------------------

	mu_mos=[]
	for iScn in range(meta[pNam]['Project']['N Scenario']):
		d={}
		for k in mos[iScn]['v1']['Mean'].keys():
			d[k]=mos[iScn]['v1']['Mean'][k]['Ensemble Mean']
		for k in mos[iScn]['Cashflow']['Mean'].keys():
			d[k]=mos[iScn]['Cashflow']['Mean'][k]['Ensemble Mean']
		mu_mos.append(d)

	cil_mos=[]
	for iScn in range(meta[pNam]['Project']['N Scenario']):
		d={}
		for k in mos[iScn]['v1']['Mean'].keys():
			d[k]=mos[iScn]['v1']['Mean'][k]['Ensemble CIL']
		for k in mos[iScn]['Cashflow']['Mean'].keys():
			d[k]=mos[iScn]['Cashflow']['Mean'][k]['Ensemble CIL']
		cil_mos.append(d)

	cih_mos=[]
	for iScn in range(meta[pNam]['Project']['N Scenario']):
		d={}
		for k in mos[iScn]['v1']['Mean'].keys():
			d[k]=mos[iScn]['v1']['Mean'][k]['Ensemble CIH']
		for k in mos[iScn]['Cashflow']['Mean'].keys():
			d[k]=mos[iScn]['Cashflow']['Mean'][k]['Ensemble CIH']
		cih_mos.append(d)

	sdl_mos=[]
	for iScn in range(meta[pNam]['Project']['N Scenario']):
		d={}
		for k in mos[iScn]['v1']['Mean'].keys():
			d[k]=mos[iScn]['v1']['Mean'][k]['Ensemble Mean']-mos[iScn]['v1']['Mean'][k]['Ensemble SD']
		for k in mos[iScn]['Cashflow']['Mean'].keys():
			d[k]=mos[iScn]['Cashflow']['Mean'][k]['Ensemble Mean']-mos[iScn]['Cashflow']['Mean'][k]['Ensemble SD']
		sdl_mos.append(d)

	sdh_mos=[]
	for iScn in range(meta[pNam]['Project']['N Scenario']):
		d={}
		for k in mos[iScn]['v1']['Mean'].keys():
			d[k]=mos[iScn]['v1']['Mean'][k]['Ensemble Mean']+mos[iScn]['v1']['Mean'][k]['Ensemble SD']
		for k in mos[iScn]['Cashflow']['Mean'].keys():
			d[k]=mos[iScn]['Cashflow']['Mean'][k]['Ensemble Mean']+mos[iScn]['Cashflow']['Mean'][k]['Ensemble SD']
		sdh_mos.append(d)

#	p1_mos=[]
#	for iScn in range(meta[pNam]['Project']['N Scenario']):
#		d={}
#		for k in mos[iScn]['v1']['Sum'].keys():
#			d[k]=mos[iScn]['v1']['Sum'][k]['Ensemble P1']
#		for k in mos[iScn]['Cashflow']['Sum'].keys():
#			d[k]=mos[iScn]['Cashflow']['Sum'][k]['Ensemble P1']
#		p1_mos.append(d)
#
#	p10_mos=[]
#	for iScn in range(meta[pNam]['Project']['N Scenario']):
#		d={}
#		for k in mos[iScn]['v1']['Sum'].keys():
#			d[k]=mos[iScn]['v1']['Sum'][k]['Ensemble P10']
#		for k in mos[iScn]['Cashflow']['Sum'].keys():
#			d[k]=mos[iScn]['Cashflow']['Sum'][k]['Ensemble P10']
#		p10_mos.append(d)
#
#
#	p90_mos=[]
#	for iScn in range(meta[pNam]['Project']['N Scenario']):
#		d={}
#		for k in mos[iScn]['v1']['Sum'].keys():
#			d[k]=mos[iScn]['v1']['Sum'][k]['Ensemble P90']
#		for k in mos[iScn]['Cashflow']['Sum'].keys():
#			d[k]=mos[iScn]['Cashflow']['Sum'][k]['Ensemble P90']
#		p90_mos.append(d)
#
#	p99_mos=[]
#	for iScn in range(meta[pNam]['Project']['N Scenario']):
#		d={}
#		for k in mos[iScn]['v1']['Sum'].keys():
#			d[k]=mos[iScn]['v1']['Sum'][k]['Ensemble P99']
#		for k in mos[iScn]['Cashflow']['Sum'].keys():
#			d[k]=mos[iScn]['Cashflow']['Sum'][k]['Ensemble P99']
#		p99_mos.append(d)

	return tv,mu_mos,cil_mos,cih_mos,sdl_mos,sdh_mos