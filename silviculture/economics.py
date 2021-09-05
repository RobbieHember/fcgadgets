
#%% Import modules

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import openpyxl
import gc as garc
import time
from fcgadgets.macgyver import utilities_general as gu
from fcgadgets.cbrunner import cbrun_utilities as cbu

#%% Calculate cashflow

def CalculateNetRevenue(meta,iScn,iEns,iBat,inv,ec,v1):
    
    # Extract economic parameters
    b=meta['Param']['BEV']['Econ']
    
    # Time vector
    tv=np.arange(meta['Project']['Year Start Saving'],meta['Year'][-1]+1,1)
    tv_full=np.arange(meta['Project']['Year Start'],meta['Project']['Year End']+1,1)
    
    # Import prices
    dPrice=gu.ReadExcel(meta['Paths']['Model Code'] + '\\Parameters\\Industrial_Product_Prices_BC.xlsx','Summary')
    
    # Import costs
    dCost=gu.ReadExcel(meta['Paths']['Model Code'] + '\\Parameters\\Forest_Operation_Costs_BC.xlsx','Summary')
        
    # Initialize dictionary of economic variables
    d={}
    
    #----------------------------------------------------------------------
    # Yield
    #----------------------------------------------------------------------
        
    # Convert lumber carbon (MgC/ha) to yield (000 bd ft)        
    d['Yield Lumber']=b['wood_C_to_DM']*b['wood_DM_to_m3']*b['wood_m3_to_bd_ft']*(1/1000)*v1['C_ToLumber']
        
    # Convert plywood carbon (MgC/ha) to yield (000 sq ft)
    d['Yield Plywood']=b['wood_C_to_DM']*b['wood_DM_to_m3']*b['sq_ft_plywood_per_m3']*(1/1000)*v1['C_ToPlywood']
        
    # Convert OSB carbon (MgC/ha) to yield (000 sq ft)
    d['Yield OSB']=b['wood_C_to_DM']*b['wood_DM_to_m3']*b['sq_ft_osb_per_m3']*(1/1000)*v1['C_ToOSB']
        
    # Convert MDF carbon (MgC/ha) to yield (000 sq ft)
    d['Yield MDF']=b['wood_C_to_DM']*b['wood_DM_to_m3']*b['sq_ft_mdf_per_m3']*(1/1000)*v1['C_ToMDF']
        
    # Convert paper carbon (MgC/ha) to yield (tonnes DM/ha)
    d['Yield Paper']=b['wood_C_to_DM']*v1['C_ToPaper']
        
    # Convert pellet carbon (MgC/ha) to yield (MWh/ha)
    d['Yield Pellets']=b['wood_C_to_DM']*v1['C_ToPellets']*b['GJ per ODT']*b['MWh per GJ']
    
    # Convert power carbon (MgC/ha) to yield (?)
    d['Yield Power']=b['wood_C_to_DM']*v1['C_ToPowerGeneration']
        
    #----------------------------------------------------------------------
    # Price
    #----------------------------------------------------------------------
        
    d['Price Lumber']=np.zeros(v1['A'].shape)
    d['Price Plywood']=np.zeros(v1['A'].shape)
    d['Price OSB']=np.zeros(v1['A'].shape)
    d['Price MDF']=np.zeros(v1['A'].shape)
    d['Price Newsprint']=np.zeros(v1['A'].shape)
    d['Price PowerGeneration']=np.zeros(v1['A'].shape)
    d['Price Pellets']=np.zeros(v1['A'].shape)
    d['Price Chips']=np.zeros(v1['A'].shape)
    d['Exchange Rate US']=np.zeros(v1['A'].shape)
    d['Exchange Rate Euro']=np.zeros(v1['A'].shape)
        
    it0=np.where( (dPrice['Year']>=tv[0]) & (dPrice['Year']<=tv[-1]) )[0]
    it1=np.where( (tv>=dPrice['Year'][0]) & (tv<=dPrice['Year'][-1]) )[0]
        
    for iStand in range(v1['A'].shape[1]):
        
        d['Price Lumber'][it1,iStand]=dPrice['Price Lumber SPF 2x4 (US$/mbf)'][it0]
        d['Price Plywood'][it1,iStand]=dPrice['Price Plywood (CDN$/000 sq ft)'][it0]
        d['Price OSB'][it1,iStand]=dPrice['Price OSB (CDN$/000 sq ft)'][it0]
        d['Price MDF'][it1,iStand]=dPrice['Price MDF (CDN$/000 sq ft)'][it0]
        d['Price Newsprint'][it1,iStand]=dPrice['Price Newsprint (US$/tonne)'][it0]
        #d['Price Chips'][it1,iStand]=dPrice['Price Chips (US$/ton)'][it0]
        d['Price PowerGeneration'][it1,iStand]=0.0 #dPrice['Price PowerGeneration (CDN$/tC)'][it0]
        d['Price Pellets'][it1,iStand]=dPrice['Price Pellets (Euro/MWh CIF)'][it0]
        d['Exchange Rate US'][it1,iStand]=dPrice['Exchange Rate (US to CDN)'][it0]
        d['Exchange Rate Euro'][it1,iStand]=dPrice['Exchange Rate (Euro to CDN)'][it0]

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
    d['Cost Slashpile Burn']=np.zeros(v1['A'].shape)
    d['Harvest Vol Merch']=np.zeros(v1['A'].shape)
    d['Harvest Vol Resid']=np.zeros(v1['A'].shape)
        
    for iStand in range(v1['A'].shape[1]):
        
        #----------------------------------------------------------------------
        # Nutrient management
        #----------------------------------------------------------------------
        
        for k in range(meta['Core']['Max Events Per Year']):
                
            ind=np.where( (ec['ID_Type'][:,iStand,k]==meta['LUT']['Dist']['Fertilization Aerial']) )[0]
                
            if ind.size==0:
                continue
                
            Year=tv_full[ind]
            it0=np.where(dPrice['Year']==Year)[0]
            it1=np.where(tv==Year)[0]
                
            if it1.size==0:
                continue
                
            d['Cost Nutrient Management'][it1,iStand]=dCost['Cost Nutrient Purchase (CDN$/ha)'][it0]+ \
                dCost['Cost Nurtrient Application (CDN$/ha)'][it0]+ \
                dCost['Cost Nutrient Overhead (CDN$/ha)'][it0]
        
        #----------------------------------------------------------------------
        # Planting
        #----------------------------------------------------------------------
        
        for k in range(meta['Core']['Max Events Per Year']):
                
            ind=np.where( (ec['ID_Type'][:,iStand,k]==meta['LUT']['Dist']['Planting']) )[0]
                
            if ind.size==0:
                continue
                
            Year=tv_full[ind]
            it0=np.where(dPrice['Year']==Year)[0]
            it1=np.where(tv==Year)[0]
              
            if it0.size==0:
                continue
            
            if it1.size==0:
                continue
            
            # Administrative costs
            try:
                dt=-2
                d['Cost Planting'][it1+dt,iStand]=d['Cost Planting'][it1+dt,iStand]+dCost['Cost Planting Admin (CDN$/ha)'][it0+dt]
            except:
                pass
            
            # Seedling purchase
            try:
                dt=-1
                d['Cost Planting'][it1+dt,iStand]=d['Cost Planting'][it1+dt,iStand]+dCost['Cost Seedling Purchase (CDN$/ha)'][it0+dt]
            except:
                pass
            
            # Planting
            dt=0
            d['Cost Planting'][it1+dt,iStand]=d['Cost Planting'][it1+dt,iStand]+dCost['Cost Planting (CDN$/ha)'][it0+dt]
            
            # Survey 1
            dt=0
            d['Cost Survey'][it1+dt,iStand]=d['Cost Planting'][it1+dt,iStand]+dCost['Cost Survey (CDN$/ha)'][it0+dt]
            
            # Survey 2 - this could creash if it extends beyond simulation period
            try:
                dt=+2
                d['Cost Survey'][it1+dt,iStand]=d['Cost Planting'][it1+dt,iStand]+dCost['Cost Survey (CDN$/ha)'][it0+dt]
            except:
                pass
            
            # Survey 3 - this could creash if it extends beyond simulation period
            try:
                dt=+8
                d['Cost Survey'][it1+dt,iStand]=d['Cost Planting'][it1+dt,iStand]+dCost['Cost Survey (CDN$/ha)'][it0+dt]
            except:
                pass
        
        #----------------------------------------------------------------------
        # Knockdown
        #----------------------------------------------------------------------
        
        for k in range(meta['Core']['Max Events Per Year']):
                
            ind=np.where( (ec['ID_Type'][:,iStand,k]==meta['LUT']['Dist']['Knockdown']) )[0]
                
            if ind.size==0:
                continue
                
            Year=tv_full[ind]
            it0=np.where(dPrice['Year']==Year)[0]
            it1=np.where(tv==Year)[0]
            
            if it0.size==0:
                continue
            
            if it1.size==0:
                continue
                
            d['Cost Knockdown'][it1,iStand]=dCost['Cost Knockdown (CDN$/ha)'][it0]
        
        #----------------------------------------------------------------------
        # Ripping
        #----------------------------------------------------------------------
        
        for k in range(meta['Core']['Max Events Per Year']):
                
            ind=np.where( (ec['ID_Type'][:,iStand,k]==meta['LUT']['Dist']['Ripping']) | (ec['ID_Type'][:,iStand,k]==meta['LUT']['Dist']['Disc Trenching']) )[0]
                
            if ind.size==0:
                continue
                
            Year=tv_full[ind]
            it0=np.where(dPrice['Year']==Year)[0]
            it1=np.where(tv==Year)[0]
              
            if it0.size==0:
                continue
            
            if it1.size==0:
                continue
                
            d['Cost Ripping'][it1,iStand]=dCost['Cost Ripping (CDN$/ha)'][it0] 
            
        #----------------------------------------------------------------------    
        # Harvesting
        #----------------------------------------------------------------------
        
        for k in range(meta['Core']['Max Events Per Year']):
                
            ind=np.where( (ec['ID_Type'][:,iStand,k]==meta['LUT']['Dist']['Harvest']) | \
                         (ec['ID_Type'][:,iStand,k]==meta['LUT']['Dist']['Harvest Salvage']) | \
                         (ec['ID_Type'][:,iStand,k]==meta['LUT']['Dist']['Harvest Custom 1']) | \
                         (ec['ID_Type'][:,iStand,k]==meta['LUT']['Dist']['Harvest Custom 2']) | \
                         (ec['ID_Type'][:,iStand,k]==meta['LUT']['Dist']['Harvest Custom 3']) | \
                         (ec['ID_Type'][:,iStand,k]==meta['LUT']['Dist']['Harvest Custom 4']) | \
                         (ec['ID_Type'][:,iStand,k]==meta['LUT']['Dist']['Harvest Custom 5']) | \
                         (ec['ID_Type'][:,iStand,k]==meta['LUT']['Dist']['Harvest Custom 6']) )[0]
                
            if ind.size==0:
                continue
                
            for i_ind in range(ind.size):                    
                    
                Year=tv_full[ind[i_ind]]
                it0=np.where(dCost['Year']==Year)[0]
                it1=np.where(tv==Year)[0]
                    
                if (it1.size==0) | (it0.size==0):
                    continue

                # Merchantable wood
                
                Removed_C=v1['C_ToMillMerch'][it1,iStand]+v1['C_ToMillSnagStem'][it1,iStand]                
                Removed_DM=b['wood_C_to_DM']*Removed_C                
                Removed_V=b['wood_DM_to_m3']*Removed_DM                
                Removed_mbf=b['wood_m3_to_bd_ft']*(1/1000)*Removed_V
                    
                d['Harvest Vol Merch'][it1,iStand]=Removed_V

                if (inv['ID_BECZ'][0,iStand]==meta['LUT']['VRI']['BEC_ZONE_CODE']['CWH']) | (inv['ID_BECZ'][0,iStand]==meta['LUT']['VRI']['BEC_ZONE_CODE']['CDF']):
                    d['Cost Harvest Overhead'][it1,iStand]=dCost['Cost Harvest Overhead Coast (CDN$/ha)'][it0]                            
                    d['Cost Harvest Hauling'][it1,iStand]=dCost['Cost Transportation Coast (CDN$/m3)'][it0]*Removed_V
                    d['Cost Roads'][it1,iStand]=dCost['Cost Roads Coast (CDN$/mbf)'][it0]*Removed_mbf
                        
                else:
                    d['Cost Harvest Overhead'][it1,iStand]=dCost['Cost Harvest Overhead Interior (CDN$/ha)'][it0]
                    d['Cost Harvest Hauling'][it1,iStand]=dCost['Cost Transportation Interior (CDN$/m3)'][it0]*Removed_V
                    d['Cost Roads'][it1,iStand]=dCost['Cost Roads Interior (CDN$/mbf)'][it0]*Removed_mbf
                    
                d['Cost Harvest Felling and Piling'][it1,iStand]=dCost['Cost Harvest Felling and Piling (CDN$/m3)'][it0]*Removed_V
                
                d['Cost Milling'][it1,iStand]=dCost['Cost Milling (CDN$/mbf)'][it0]*Removed_mbf
                
                # Residual fibre
                
                Removed_C=v1['C_ToMillNonMerch'][it1,iStand]
                Removed_DM=b['wood_C_to_DM']*Removed_C                
                Removed_V=b['wood_DM_to_m3']*Removed_DM                
                #Removed_mbf=b['wood_m3_to_bd_ft']*(1/1000)*Removed_V
                
                d['Harvest Vol Resid'][it1,iStand]=Removed_V
                
                d['Cost Harvest Residuals'][it1,iStand]=dCost['Cost Residual Haul and Grind (CDN$/m3)'][it0]*Removed_V
        
        #----------------------------------------------------------------------
        # Slashpile burning
        #----------------------------------------------------------------------
        
        for k in range(meta['Core']['Max Events Per Year']):
                
            ind=np.where( (ec['ID_Type'][:,iStand,k]==meta['LUT']['Dist']['Slashpile Burn']) )[0]
                
            if ind.size==0:
                continue
                
            Year=tv_full[ind]
            it0=np.where(dPrice['Year']==Year)[0]
            it1=np.where(tv==Year)[0]
              
            if it0.size==0:
                continue
            
            if it1.size==0:
                continue
                
            Burned_C=v1['C_ToSlashpileBurn'][it1,iStand]               
            Burned_V=b['wood_DM_to_m3']*b['wood_C_to_DM']*Burned_C            
            d['Cost Slashpile Burn'][it1,iStand]=dCost['Cost Ripping (CDN$/ha)'][it0]*Burned_V
        
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
        d['Cost Slashpile Burn']+ \
        d['Cost Knockdown']

    # Gross revenue from sale of lumber: (CDN$/ha) = (US$/mbf) * (CDN$/US$) * (mbf/ha)
    d['Revenue Lumber']=d['Price Lumber']*(1/d['Exchange Rate US'])*d['Yield Lumber']
        
    # Revenue from sale of plywood: (CDN$/ha) = (CDN$/000 sq ft) * (sq ft/ha)
    d['Revenue Plywood']=d['Price Plywood']*d['Yield Plywood']
        
    # Revenue from sale of OSB: (CDN$/ha) = (CDN$/000 sq ft) * (sq ft/ha)
    d['Revenue OSB']=d['Price OSB']*d['Yield OSB']
        
    # Revenue from sale of MDF: (CDN$/ha) = (CDN$/000 sq ft) * (sq ft/ha)
    d['Revenue MDF']=d['Price MDF']*d['Yield MDF']
        
    # Revenue from sale of newsprint = (US$/tonne) * (CDN$/US$) * (tonne/ha)
    d['Revenue Paper']=d['Price Newsprint']*(1/d['Exchange Rate US'])*d['Yield Paper']
        
    # Revenue from sale of power generation
    d['Revenue PowerGeneration']=0*d['Revenue Lumber']
    
    # Revenue from sale of pellets (CDN$/ha) = (Euro$/MWh) * (CDN$/Euro$) * (MWh/ha)
    d['Revenue Pellets']=d['Price Pellets']*(1/d['Exchange Rate Euro'])*d['Yield Pellets']
    
    # Revenue from sale of chips
    d['Revenue Chips']=0*d['Revenue Lumber']
        
    # Remove all NaNs
    for k in d.keys():
        d[k]=np.nan_to_num(d[k])
        
    # Gross revenue ($)
    d['Revenue Gross']=d['Revenue Lumber']+d['Revenue Plywood']+d['Revenue OSB']+d['Revenue MDF']+ \
        d['Revenue Paper']+d['Revenue PowerGeneration']+d['Revenue Pellets']+d['Revenue Chips']            
        
    # Net revenue
    d['Revenue Net']=d['Revenue Gross']-d['Cost Total']
        
    # Time
    #d['t']=np.tile(np.reshape(tv,(-1,1)),(1,v1['A'].shape[1]))
        
    # Present value
    #d['PV']=d['Net Revenue']/(1+b['Interest Rate'])**d['t']
    
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
