#%% Import modules
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import openpyxl
import gc as garc
import copy
import time
from fcgadgets.macgyver import util_general as gu
from fcgadgets.cbrunner import cbrun_util as cbu

#%% Calculate net revenue from events
def CashflowFromEventChronology(meta,pNam,iScn,iEns,iBat,inv,ec,v1):

	# Extract biophysical parameters
	b=meta['Param']['BEV']['Biophysical']

	# Time vector
	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Year'][-1]+1,1)
	tv_full=np.arange(meta[pNam]['Project']['Year Start'],meta[pNam]['Project']['Year End']+1,1)

	# Discounting variables
	r_disc=meta['Param']['BEV']['Economics']['Discount Rate']
	t_disc=np.maximum(0,tv-meta[pNam]['Project']['Year Project'])
	t_disc=np.tile(t_disc,(v1['A'].shape[1],1)).T

	# Import prices
	dCP=gu.ReadExcel(meta['Paths']['Model']['Code'] + '\\Parameters\\Parameters_CostsAndPrices.xlsx','Data')

	# Initialize dictionary of economic variables
	d={}
	
	#----------------------------------------------------------------------
	# Price
	#----------------------------------------------------------------------
	d['Price Lumber']=np.zeros(v1['A'].shape)
	d['Price Plywood']=np.zeros(v1['A'].shape)
	d['Price OSB']=np.zeros(v1['A'].shape)
	d['Price MDF']=np.zeros(v1['A'].shape)
	d['Price Newsprint']=np.zeros(v1['A'].shape)
	d['Price PowerFacilityDom']=np.zeros(v1['A'].shape)
	d['Price PowerGrid']=np.zeros(v1['A'].shape)
	d['Price PelletExport']=np.zeros(v1['A'].shape)
	d['Price PelletDom']=np.zeros(v1['A'].shape)
	d['Price LogExport']=np.zeros(v1['A'].shape)
	d['Price FirewoodDom']=np.zeros(v1['A'].shape)
	d['Exchange Rate US']=np.zeros(v1['A'].shape)
	d['Exchange Rate Euro']=np.zeros(v1['A'].shape)

	it0=np.where( (dCP['Year']>=tv[0]) & (dCP['Year']<=tv[-1]) )[0]
	it1=np.where( (tv>=dCP['Year'][0]) & (tv<=dCP['Year'][-1]) )[0]

	for iStand in range(v1['A'].shape[1]):
		d['Price Lumber'][it1,iStand]=dCP['Price Lumber SPF 2x4 (USD/mbf)'][it0]
		d['Price Plywood'][it1,iStand]=dCP['Price Plywood (CAD/000 sq ft)'][it0]
		d['Price OSB'][it1,iStand]=dCP['Price OSB (CAD/000 sq ft)'][it0]
		d['Price MDF'][it1,iStand]=dCP['Price MDF (CAD/000 sq ft)'][it0]
		d['Price Newsprint'][it1,iStand]=dCP['Price Newsprint (USD/tonne)'][it0]
		d['Price PowerFacilityDom'][it1,iStand]=dCP['Price Power Facility (CAD/MWh)'][it0]
		d['Price PowerGrid'][it1,iStand]=dCP['Price Power Grid (CAD/MWh)'][it0]
		d['Price PelletExport'][it1,iStand]=dCP['Price Pellet Export (Euro/MWh CIF)'][it0]
		d['Price PelletDom'][it1,iStand]=dCP['Price Pellet Domestic (CAD/MWh)'][it0]
		d['Price LogExport'][it1,iStand]=dCP['Price Log Export (CAD/m3)'][it0]
		d['Price FirewoodDom'][it1,iStand]=dCP['Price Firewood (CAD/ODT)'][it0]
		d['Exchange Rate US'][it1,iStand]=dCP['Exchange Rate (US to CAD)'][it0]
		d['Exchange Rate Euro'][it1,iStand]=dCP['Exchange Rate (Euro to CAD)'][it0]

	#----------------------------------------------------------------------
	# Cost
	#----------------------------------------------------------------------
	d['Cost Roads']=np.zeros(v1['A'].shape)
	d['Cost Harvest Overhead']=np.zeros(v1['A'].shape)
	d['Cost Harvest Felling and Piling']=np.zeros(v1['A'].shape)
	d['Cost Harvest Hauling']=np.zeros(v1['A'].shape)
	d['Cost Harvest Residuals']=np.zeros(v1['A'].shape)
	d['Cost Milling']=np.zeros(v1['A'].shape)
	d['Cost Nutrient Management']=np.zeros(v1['A'].shape)
	d['Cost Planting']=np.zeros(v1['A'].shape)
	d['Cost Survey']=np.zeros(v1['A'].shape)
	d['Cost Knockdown']=np.zeros(v1['A'].shape)
	d['Cost Ripping']=np.zeros(v1['A'].shape)
	d['Cost PAS Deactivation']=np.zeros(v1['A'].shape)
	d['Cost Slashpile Burn']=np.zeros(v1['A'].shape)
	d['Cost Aerial BTK Spray']=np.zeros(v1['A'].shape)
	#d['Harvest Vol Merch']=np.zeros(v1['A'].shape)
	#d['Harvest Vol Resid']=np.zeros(v1['A'].shape)

	for iStand in range(v1['A'].shape[1]):

		#----------------------------------------------------------------------
		# Nutrient management
		#----------------------------------------------------------------------
		for k in range(meta['Core']['Max Events Per Year']):

			try:
				ind=np.where( (ec['ID Event Type'][:,iStand,k]==meta['LUT']['Event']['Nutrient App Aerial']) )[0]
			except:
				print(list(ec.keys()))

			for i in ind:
				Year=tv_full[i]
				it0=np.where(dCP['Year']==Year)[0]
				it1=np.where(tv==Year)[0]
				if it1.size==0:
					continue
				d['Cost Nutrient Management'][it1,iStand]=dCP['Cost Nutrient Purchase (CAD/ha)'][it0]+ \
					dCP['Cost Nurtrient Application (CAD/ha)'][it0]+ \
					dCP['Cost Nutrient Overhead (CAD/ha)'][it0]

		#----------------------------------------------------------------------
		# Aerial BTK Spray (Btk)
		#----------------------------------------------------------------------
		for k in range(meta['Core']['Max Events Per Year']):

			ind=np.where( (ec['ID Event Type'][:,iStand,k]==meta['LUT']['Event']['Aerial BTK Spray']) )[0]

			if ind.size==0:
				continue

			Year=tv_full[ind]
			it0=np.where(dCP['Year']==Year)[0]
			it1=np.where(tv==Year)[0]

			if it1.size==0:
				continue

			d['Cost Aerial BTK Spray'][it1,iStand]=dCP['Cost Aerial BTK Spray (CAD/ha)'][it0]

		#----------------------------------------------------------------------
		# Planting
		#----------------------------------------------------------------------
		for k in range(meta['Core']['Max Events Per Year']):

			ind=np.where( (ec['ID Event Type'][:,iStand,k]==meta['LUT']['Event']['Planting']) )[0]

			if ind.size==0:
				continue

			for i_ind in range(ind.size):

				Year=tv_full[ind[i_ind]]
				it0=np.where(dCP['Year']==Year)[0]
				it1=np.where(tv==Year)[0]

				if it0.size==0:
					continue

				if it1.size==0:
					continue

				# Administrative costs
				try:
					dt=-2
					d['Cost Planting'][it1+dt,iStand]=d['Cost Planting'][it1+dt,iStand]+dCP['Cost Planting Admin (CAD/ha)'][it0+dt]
				except:
					pass

				# Seedling purchase
				try:
					dt=-1
					d['Cost Planting'][it1+dt,iStand]=d['Cost Planting'][it1+dt,iStand]+dCP['Cost Seedling Purchase (CAD/ha)'][it0+dt]
				except:
					pass

				# Planting
				dt=0
				d['Cost Planting'][it1+dt,iStand]=d['Cost Planting'][it1+dt,iStand]+dCP['Cost Planting (CAD/ha)'][it0+dt]

				# Survey 1
				dt=0
				d['Cost Survey'][it1+dt,iStand]=d['Cost Survey'][it1+dt,iStand]+dCP['Cost Survey (CAD/ha)'][it0+dt]

				# Survey 2 - this could creash if it extends beyond simulation period
				try:
					dt=+2
					d['Cost Survey'][it1+dt,iStand]=d['Cost Survey'][it1+dt,iStand]+dCP['Cost Survey (CAD/ha)'][it0+dt]
				except:
					pass

				# Survey 3 - this could creash if it extends beyond simulation period
				try:
					dt=+8
					d['Cost Survey'][it1+dt,iStand]=d['Cost Survey'][it1+dt,iStand]+dCP['Cost Survey (CAD/ha)'][it0+dt]
				except:
					pass

		#----------------------------------------------------------------------
		# Knockdown
		#----------------------------------------------------------------------
		for k in range(meta['Core']['Max Events Per Year']):

			ind=np.where( (ec['ID Event Type'][:,iStand,k]==meta['LUT']['Event']['Knockdown']) )[0]

			if ind.size==0:
				continue
			if ind.size>1:
				ind=ind[0]

			Year=tv_full[ind]
			it0=np.where(dCP['Year']==Year)[0]
			it1=np.where(tv==Year)[0]

			if it0.size==0:
				continue

			if it1.size==0:
				continue

			d['Cost Knockdown'][it1,iStand]=dCP['Cost Knockdown (CAD/ha)'][it0]

		#----------------------------------------------------------------------
		# Ripping
		#----------------------------------------------------------------------
		for k in range(meta['Core']['Max Events Per Year']):

			ind=np.where( (ec['ID Event Type'][:,iStand,k]==meta['LUT']['Event']['Mechanical Site Prep']) )[0]

			if ind.size==0:
				continue

			Year=tv_full[ind]
			it0=np.where(dCP['Year']==Year)[0]
			it1=np.where(tv==Year)[0]

			if it0.size==0:
				continue

			if it1.size==0:
				continue

			d['Cost Ripping'][it1,iStand]=dCP['Cost Ripping (CAD/ha)'][it0]

		#----------------------------------------------------------------------
		# Deactivation of permanent access structures
		#----------------------------------------------------------------------
		for k in range(meta['Core']['Max Events Per Year']):

			if 'PAS Deactivation' not in meta['LUT']['Event']:
				continue

			ind=np.where( (ec['ID Event Type'][:,iStand,k]==meta['LUT']['Event']['PAS Deactivation']) )[0]

			if ind.size==0:
				continue

			Year=tv_full[ind]
			it0=np.where(dCP['Year']==Year)[0]
			it1=np.where(tv==Year)[0]

			if it0.size==0:
				continue

			if it1.size==0:
				continue

			d['Cost PAS Deactivation'][it1,iStand]=dCP['Cost PAS Deactivation (CAD/ha)'][it0]

		#----------------------------------------------------------------------
		# Harvesting
		#----------------------------------------------------------------------
		for k in range(meta['Core']['Max Events Per Year']):

			ind=np.where( (ec['ID Event Type'][:,iStand,k]==meta['LUT']['Event']['Harvest']) )[0]

			if ind.size==0:
				continue

			for i_ind in range(ind.size):

				Year=tv_full[ind[i_ind]]
				it0=np.where(dCP['Year']==Year)[0]
				it1=np.where(tv==Year)[0]

				if (it1.size==0) | (it0.size==0):
					continue

				# Merchantable wood volume
				Removed_V=v1['V_ToMillMerchTotal'][it1,iStand]

				# Merchantable wood thousand board feet (mbf is thousand)
				Removed_mbf=b['Ratio bd ft Lumber per m3']*Removed_V/1000

				# Log-size effect
				if v1['LogSizeEnhancement'][it1,iStand]==0:
					lsef_skidding=1.0
					lsef_roads=1.0
				elif v1['LogSizeEnhancement'][it1,iStand]==1:
					lsef_skidding=0.98
					lsef_roads=0.97
				elif v1['LogSizeEnhancement'][it1,iStand]==2:
					lsef_skidding=0.96
					lsef_roads=0.94
				elif v1['LogSizeEnhancement'][it1,iStand]==3:
					lsef_skidding=0.94
					lsef_roads=0.91
				elif v1['LogSizeEnhancement'][it1,iStand]==4:
					lsef_skidding=0.92
					lsef_roads=0.88
				else:
					lsef_skidding=0.9
					lsef_roads=0.85

				#d['Harvest Vol Merch'][it1,iStand]=Removed_V

				if (inv['ID_BGCZ'][0,iStand]==meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE']['CWH']) | (inv['ID_BGCZ'][0,iStand]==meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE']['CDF']):

					d['Cost Harvest Overhead'][it1,iStand]=dCP['Cost Harvest Overhead Coast (CAD/ha)'][it0]

					d['Cost Harvest Hauling'][it1,iStand]=dCP['Cost Transportation Coast (CAD/m3)'][it0]*Removed_V

					d['Cost Roads'][it1,iStand]=dCP['Cost Roads Coast (CAD/mbf)'][it0]*lsef_roads*Removed_mbf

				else:

					d['Cost Harvest Overhead'][it1,iStand]=dCP['Cost Harvest Overhead Interior (CAD/ha)'][it0]

					d['Cost Harvest Hauling'][it1,iStand]=dCP['Cost Transportation Interior (CAD/m3)'][it0]*Removed_V

					d['Cost Roads'][it1,iStand]=dCP['Cost Roads Interior (CAD/mbf)'][it0]*lsef_roads*Removed_mbf

				d['Cost Harvest Felling and Piling'][it1,iStand]=dCP['Cost Harvest Felling and Piling (CAD/m3)'][it0]*lsef_skidding*Removed_V

				d['Cost Milling'][it1,iStand]=dCP['Cost Milling (CAD/mbf)'][it0]*Removed_mbf

				# Residual fibre
				d['Cost Harvest Residuals'][it1,iStand]=dCP['Cost Residual Haul and Grind (CAD/m3)'][it0]*v1['V_ToMillNonMerch'][it1,iStand]

		#----------------------------------------------------------------------
		# Slashpile burning
		#----------------------------------------------------------------------
		for k in range(meta['Core']['Max Events Per Year']):

			ind=np.where( (ec['ID Event Type'][:,iStand,k]==meta['LUT']['Event']['Slashpile Burn']) )[0]

			if ind.size==0:
				continue

			for i_ind in range(ind.size):

				Year=tv_full[ind[i_ind]]
				it0=np.where(dCP['Year']==Year)[0]
				it1=np.where(tv==Year)[0]

				if it0.size==0:
					continue

				if it1.size==0:
					continue

				Burned_C=v1['C_ToSlashpileBurnTot'][it1,iStand]
				Burned_V=Burned_C/b['Carbon Content Wood']/b['Density Wood']
				d['Cost Slashpile Burn'][it1,iStand]=(dCP['Cost Slashpile Burn (CAD/m3)'][it0])*Burned_V

	# Total cost
	d['Cost Total']=d['Cost Roads']+ \
		d['Cost Harvest Overhead']+ \
		d['Cost Harvest Felling and Piling']+ \
		d['Cost Harvest Hauling']+ \
		d['Cost Harvest Residuals']+ \
		d['Cost Milling']+ \
		d['Cost Nutrient Management']+ \
		d['Cost Planting']+ \
		d['Cost Survey']+ \
		d['Cost Ripping']+ \
		d['Cost PAS Deactivation']+ \
		d['Cost Slashpile Burn']+ \
		d['Cost Knockdown']+ \
		d['Cost Aerial BTK Spray']

	# Total silviculture
	d['Cost Silviculture Total']=d['Cost Harvest Residuals']+ \
		d['Cost Nutrient Management']+ \
		d['Cost Planting']+ \
		d['Cost Survey']+ \
		d['Cost Ripping']+ \
		d['Cost PAS Deactivation']+ \
		d['Cost Slashpile Burn']+ \
		d['Cost Knockdown']+ \
		d['Cost Aerial BTK Spray']

	# Thousand board feet lumber / ha
	mbf_Lumber=b['Ratio bd ft Lumber per m3']*v1['ODT Lumber']/b['Density Wood']/1000

	# Thousand sq ft plywood / ha
	Thousand_sq_ft_Plywood=b['Ratio sq ft Plywood per m3']*v1['ODT Plywood']/b['Density Wood']/1000

	# Thousand sq ft OSB / ha
	Thousand_sq_ft_OSB=b['Ratio sq ft OSB per m3']*v1['ODT OSB']/b['Density Wood']/1000

	# Thousand sq ft MDF / ha
	Thousand_sq_ft_MDF=b['Ratio sq ft MDF per m3']*v1['ODT MDF']/b['Density Wood']/1000

	# Gross revenue from sale of lumber: (CAD/ha) = (USD/mbf) * (CAD/USD) * (000 bd ft/ha)
	d['Revenue Lumber']=d['Price Lumber']*(1/d['Exchange Rate US'])*mbf_Lumber

	# Revenue from sale of plywood: (CAD/ha) = (CAD/000 sq ft) * (1000 sq ft/ha)
	d['Revenue Plywood']=d['Price Plywood']*Thousand_sq_ft_Plywood

	# Revenue from sale of OSB: (CAD/ha) = (CAD/000 sq ft) * (1000 sq ft/ha)
	d['Revenue OSB']=d['Price OSB']*Thousand_sq_ft_OSB

	# Revenue from sale of MDF: (CAD/ha) = (CAD/000 sq ft) * (sq ft/ha)
	d['Revenue MDF']=d['Price MDF']*Thousand_sq_ft_MDF

	# Revenue from sale of newsprint = (USD/ODT) * (CAD/USD) * (ODT/ha)
	d['Revenue Paper']=d['Price Newsprint']*(1/d['Exchange Rate US'])*v1['ODT Paper']

	# Revenue from domestic facility power
	d['Revenue PowerFacilityDom']=d['Price PowerFacilityDom']*b['Ratio MWh per GJ']*v1['GJ PowerFacilityDom']

	# Revenue from IPP sales
	d['Revenue PowerGrid']=d['Price PowerGrid']*b['Ratio MWh per GJ']*v1['GJ PowerGrid']

	# Revenue from sale of pellets (CAD/ha) = (Euro$/MWh) * (CAD/Euro$) * (MWh/ha)
	d['Revenue PelletExport']=d['Price PelletExport']*(1/d['Exchange Rate Euro'])*b['Ratio MWh per GJ']*v1['GJ PelletExport']

	d['Revenue PelletDom']=d['Price PelletDom']*b['Ratio MWh per GJ']*(v1['GJ PelletDomGrid']+v1['GJ PelletDomRNG'])

	# Revenue from sale of firewood (CAD/ha) = (CAD/ODT) * (ODT/ha)
	d['Revenue FirewoodDom']=d['Price FirewoodDom']*v1['ODT FirewoodDom']

	# Revenue from sale of log Export (CAD/ha) = (CAD/m3) * (m3/ha)
	d['Revenue LogExport']=d['Price LogExport']*v1['ODT LogExport']/b['Density Wood']

	# Remove all NaNs
	for k in d.keys():
		d[k]=np.nan_to_num(d[k])

	# Gross revenue ($)
	d['Revenue Gross']=d['Revenue Lumber']+d['Revenue Plywood']+d['Revenue OSB']+d['Revenue MDF']+ \
		d['Revenue Paper']+d['Revenue PowerFacilityDom']+d['Revenue PowerGrid']+d['Revenue PelletExport']+d['Revenue PelletDom']+d['Revenue FirewoodDom']+d['Revenue LogExport']

	# Net revenue
	d['Revenue Net']=d['Revenue Gross']-d['Cost Total']

	# Discounting
	d['Revenue Net Disc']=d['Revenue Net'].copy()/((1+r_disc)**t_disc)
	d['Revenue Gross Disc']=d['Revenue Gross'].copy()/((1+r_disc)**t_disc)
	d['Cost Total Disc']=d['Cost Total'].copy()/((1+r_disc)**t_disc)
	d['Cost Nutrient Management Disc']=d['Cost Nutrient Management'].copy()/((1+r_disc)**t_disc)
	d['Cost Silviculture Total Disc']=d['Cost Silviculture Total'].copy()/((1+r_disc)**t_disc)

	# Cumulative
	iT=np.where(tv>=meta[pNam]['Project']['Year Start Cumulative'])[0]
	List=['Cost Total','Cost Silviculture Total','Cost Nutrient Management', \
		  'Cost Total Disc','Cost Silviculture Total Disc','Cost Nutrient Management Disc', \
		  'Revenue Gross','Revenue Gross Disc', \
		  'Revenue Net','Revenue Net Disc']
	for vnam in List:
		d[vnam + '_cumu']=np.zeros(d[vnam].shape)
		d[vnam + '_cumu'][iT,:]=np.cumsum(d[vnam][iT,:],axis=0)

	return d

#%% Calculate economic results for each scenario
def CalculateEconomics(meta,v1):
	econ=[]
	for iScn in range(len(v1)):
		iEns=0;
		iBat=0;
		inv=gu.ipickle(meta['Paths']['Input Scenario'][iScn] + '\\Inventory_Bat' + cbu.FixFileNum(iBat) + '.pkl')
		ec=gu.ipickle(meta['Paths']['Input Scenario'][iScn] + '\\Events_Ens' + cbu.FixFileNum(iEns) + '_Bat' + cbu.FixFileNum(iBat) + '.pkl')
		ec=cbu.EventChronologyDecompress(meta,ec,iScn,iEns,iBat)
		econ.append(CalculateNetRevenue(meta,iScn,iEns,iBat,inv,ec,v1[iScn]))

	return econ
