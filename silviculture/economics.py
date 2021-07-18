
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
    b=meta['Param']['Econ']
    
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
    d['Yield Lumber']=b['wood_C_to_DM']*b['wood_DM_to_m3']*b['wood_m3_to_bd_ft']*(1/1000)*v1['C_Lumber']
        
    # Convert plywood carbon (MgC/ha) to yield (000 sq ft)
    d['Yield Plywood']=b['wood_C_to_DM']*b['wood_DM_to_m3']*b['sq_ft_plywood_per_m3']*(1/1000)*v1['C_Plywood']
        
    # Convert OSB carbon (MgC/ha) to yield (000 sq ft)
    d['Yield OSB']=b['wood_C_to_DM']*b['wood_DM_to_m3']*b['sq_ft_osb_per_m3']*(1/1000)*v1['C_OSB']
        
    # Convert MDF carbon (MgC/ha) to yield (000 sq ft)
    d['Yield MDF']=b['wood_C_to_DM']*b['wood_DM_to_m3']*b['sq_ft_mdf_per_m3']*(1/1000)*v1['C_MDF']
        
    # Convert paper carbon (MgC/ha) to yield (tonnes DM/ha)
    d['Yield Newsprint']=b['wood_C_to_DM']*v1['C_Paper']
        
    # Convert bioenergy carbon (MgC/ha) to yield (?)
    d['Yield Pellets']=0*v1['C_Fuel']
        
    #----------------------------------------------------------------------
    # Price
    #----------------------------------------------------------------------
        
    d['Price Lumber']=np.zeros(v1['A'].shape)
    d['Price Plywood']=np.zeros(v1['A'].shape)
    d['Price OSB']=np.zeros(v1['A'].shape)
    d['Price MDF']=np.zeros(v1['A'].shape)
    d['Price Newsprint']=np.zeros(v1['A'].shape)
    d['Price Bioenergy']=np.zeros(v1['A'].shape)
    d['Exchange Rate']=np.zeros(v1['A'].shape)
        
    it0=np.where( (dPrice['Year']>=tv[0]) & (dPrice['Year']<=tv[-1]) )[0]
    it1=np.where( (tv>=dPrice['Year'][0]) & (tv<=dPrice['Year'][-1]) )[0]
        
    for iStand in range(v1['A'].shape[1]):
        
        d['Price Lumber'][it1,iStand]=dPrice['Price Lumber SPF 2x4 (US$/mbf)'][it0]
        d['Price Plywood'][it1,iStand]=dPrice['Price Plywood (CDN$/000 sq ft)'][it0]
        d['Price OSB'][it1,iStand]=dPrice['Price OSB (CDN$/000 sq ft)'][it0]
        d['Price MDF'][it1,iStand]=dPrice['Price MDF (CDN$/000 sq ft)'][it0]
        d['Price Newsprint'][it1,iStand]=dPrice['Price Newsprint (US$/tonne)'][it0]
        #d['Price Bioenergy'][it1,iStand]=dPrice['Price Bioenergy'][it0]
        d['Exchange Rate'][it1,iStand]=dPrice['Exchange Rate (US to CDN)'][it0]

    #----------------------------------------------------------------------
    # Cost 
    #----------------------------------------------------------------------
        
    d['Cost Nutrient Management']=np.zeros(v1['A'].shape)
    d['Cost Roads']=np.zeros(v1['A'].shape)
    d['Cost Harvest Overhead']=np.zeros(v1['A'].shape)
    d['Cost Harvest Felling and Piling']=np.zeros(v1['A'].shape)
    d['Cost Harvest Hauling']=np.zeros(v1['A'].shape)
    d['Cost Milling']=np.zeros(v1['A'].shape)        
    d['Harvest Volume']=np.zeros(v1['A'].shape)
        
    for iStand in range(v1['A'].shape[1]):
            
        # Nutrient management
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
            
        # Harvesting
        for k in range(meta['Core']['Max Events Per Year']):
                
            ind=np.where( (ec['ID_Type'][:,iStand,k]==meta['LUT']['Dist']['Harvest']) | \
                         (ec['ID_Type'][:,iStand,k]==meta['LUT']['Dist']['Harvest Salvage']) )[0]
                
            if ind.size==0:
                continue
                
            for i_ind in range(ind.size):                    
                    
                Year=tv_full[ind[i_ind]]
                it0=np.where(dCost['Year']==Year)[0]
                it1=np.where(tv==Year)[0]
                    
                if (it1.size==0) | (it0.size==0):
                    continue

                Removed_C=v1['C_RemovedMerch'][it1,iStand]+v1['C_RemovedNonMerch'][it1,iStand]+ \
                    v1['C_RemovedSnagStem'][it1,iStand]
                Removed_DM=b['wood_C_to_DM']*Removed_C
                Removed_V=b['wood_DM_to_m3']*Removed_DM
                Removed_mbf=b['wood_m3_to_bd_ft']*(1/1000)*Removed_V
                    
                d['Harvest Volume'][it1,iStand]=Removed_V

                if (inv['ID_BECZ'][0,iStand]==meta['LUT']['VRI']['BEC_ZONE_CODE']['CWH']) | (inv['ID_BECZ'][0,iStand]==meta['LUT']['VRI']['BEC_ZONE_CODE']['CDF']):
                    d['Cost Harvest Overhead'][it1,iStand]=dCost['Cost Harvest Overhead Coast (CDN$/ha)'][it0]                            
                    d['Cost Harvest Hauling'][it1,iStand]=dCost['Cost Harvest Hauling Coast (CDN$/m3)'][it0]*Removed_V
                    d['Cost Roads'][it1,iStand]=dCost['Cost Roads Coast (CDN$/mbf)'][it0]*Removed_mbf
                        
                else:
                    d['Cost Harvest Overhead'][it1,iStand]=dCost['Cost Harvest Overhead Interior (CDN$/ha)'][it0]
                    d['Cost Harvest Hauling'][it1,iStand]=dCost['Cost Harvest Hauling Interior (CDN$/m3)'][it0]*Removed_V
                    d['Cost Roads'][it1,iStand]=dCost['Cost Roads Interior (CDN$/mbf)'][it0]*Removed_mbf
                    
                d['Cost Harvest Felling and Piling'][it1,iStand]=dCost['Cost Harvest Felling and Piling (CDN$/m3)'][it0]*Removed_V             
                    
                d['Cost Milling'][it1,iStand]=dCost['Cost Milling (CDN$/mbf)'][it0]*Removed_mbf
        
    # Total cost
    d['Cost Total']=d['Cost Roads']+ \
        d['Cost Harvest Overhead']+ \
        d['Cost Harvest Felling and Piling']+ \
        d['Cost Harvest Hauling']+ \
        d['Cost Milling']+ \
        d['Cost Nutrient Management']

    # Gross revenue from sale of lumber: (CDN$/ha) = (US$/mbf) * (CDN$/US$) * (mbf/ha)
    d['Revenue Lumber']=d['Price Lumber']*(1/d['Exchange Rate'])*d['Yield Lumber']
        
    # Revenue from sale of plywood: (CDN$/ha) = (CDN$/000 sq ft) * (sq ft/ha)
    d['Revenue Plywood']=d['Price Plywood']*d['Yield Plywood']
        
    # Revenue from sale of OSB: (CDN$/ha) = (CDN$/000 sq ft) * (sq ft/ha)
    d['Revenue OSB']=d['Price OSB']*d['Yield OSB']
        
    # Revenue from sale of MDF: (CDN$/ha) = (CDN$/000 sq ft) * (sq ft/ha)
    d['Revenue MDF']=d['Price MDF']*d['Yield MDF']
        
    # Revenue from sale of newsprint: (US$/tonne) * (CDN$/US$) * (tonne/ha)
    d['Revenue Newsprint']=d['Price Newsprint']*(1/d['Exchange Rate'])*d['Yield Newsprint']
        
    # Revenue from sale of bioenergy
    d['Revenue Bioenergy']=0*d['Revenue Lumber']
        
    # Remove all NaNs
    for k in d.keys():
        d[k]=np.nan_to_num(d[k])
        
    # Gross revenue ($)
    d['Revenue Gross']=d['Revenue Lumber']+d['Revenue Plywood']+d['Revenue OSB']+d['Revenue MDF']+ \
        d['Revenue Newsprint']+d['Revenue Bioenergy']            
        
    # Net revenue
    d['Revenue Net']=d['Revenue Gross']-d['Cost Total']
        
    # Time
    #d['t']=np.tile(np.reshape(tv,(-1,1)),(1,v1['A'].shape[1]))
        
    # Present value
    #d['PV']=d['Net Revenue']/(1+b['Interest Rate'])**d['t']
    
    return d