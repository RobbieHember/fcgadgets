'''

ONSET AND SPREAD MODELS

'''

#%% Import modules

import numpy as np
import scipy.stats as stats
from fcgadgets.utilities import utilities_general as gu
from fcgadgets.cbrunner import cbrun_utilities as cbu
from fcgadgets.taz import aspatial_stat_models as gensm

#%% Spread modelling tools

def Import_OSM_Tools(metaOS):
    
    ost={}
    
    # Displacements from a cell to its eight nearest neighbours
    #neighbourhood=((-1,-1),(-1,0),(-1,1),(0,-1),(0, 1),(1,-1),(1,0),(1,1))
    ost['neighbourhood']=((-1,1),(0,1),(1,1),(1,0),(1, -1),(0,-1),(-1,-1),(-1,0))
    #direction=['NW','N','NE','E','SE','S','SW','W']
    ost['UNSUSEPTABLE']=0
    ost['SUSEPTABLE']=1    
    ost['ONSET']=2

    ost['non']=lambda s: s if s<0 else None
    ost['mom']=lambda s: max(0,s)
    
    metaOS['ost']=ost
    metaOS['N_t']=metaOS['Time'].size
    metaOS['N_Y']=metaOS['Grid_LCC'].shape[0]
    metaOS['N_X']=metaOS['Grid_LCC'].shape[1]
    metaOS['Shape']=(metaOS['N_Y'],metaOS['N_X'])
    metaOS['Size']=metaOS['N_Y']*metaOS['N_X']
    
    return metaOS

#%% Onset-spread model annual loop
    
def OSM_Annual(metaOS,z,DataAll):
    
    for iYear in range(metaOS['Time'].size):    
    
        # Run daily loop
        z=OSM_DailyLoop(z,metaOS,iYear)
    
        # Update the years since occurrence counter        
        ind=(metaOS['Grid_YearsSinceOcc']<metaOS['YearsToRegenerate'])
        metaOS['Grid_YearsSinceOcc'][ind]=metaOS['Grid_YearsSinceOcc'][ind]+1
    
        # Identify stands that have regenerated
        ind=(metaOS['Grid_YearsSinceOcc']==metaOS['YearsToRegenerate'])
        
        # Transfer recently-affected stands back to the suseptable class
        z['Data'][ind]=metaOS['ost']['SUSEPTABLE']
        
        # Adjust counter to be just above threshold regeneration period so that
        # they will not be found in queries
        metaOS['Grid_YearsSinceOcc'][ind]=metaOS['YearsToRegenerate']+1
        
        z['Asus'][iYear]=np.sum(z['Data']==1)
        
        # Update full matrix
        DataAll[iYear,:,:]=z['Data']    
        #print(metaOS['Time'][iYear],z['Aburn'][iYear])
    
    return metaOS,z,DataAll

#%% Onset-spread model daily loop

def OSM_DailyLoop(z,metaOS,iYear):
    
    if iYear==0:
        
        # Initialize conditions at start of year
        z['Data']=metaOS['Grid_LCC'].copy()
    
    else:
        
        # Probability of regeneration
        #rn=np.random.random(metaOS['Shape'])
        #z['Data'][(metaOS['Grid_LCC']==SUSEPTABLE) & (z['Data']==UNSUSEPTABLE) & (rn<metaOS['p_Regen'])]=SUSEPTABLE        
    
        # Transition ONSET to UNSUSEPTABLE        
        if metaOS['p_Continue']>0:
            rn=np.random.random(metaOS['Shape'])
            ind=(z['Data']==metaOS['ost']['ONSET']) & (rn>=metaOS['p_Continue'])
            z['Data'][ind]=metaOS['ost']['UNSUSEPTABLE']
        else:
            z['Data'][(z['Data']==metaOS['ost']['ONSET'])]=metaOS['ost']['UNSUSEPTABLE']
    
    neighbourhood=metaOS['ost']['neighbourhood']
    mom=metaOS['ost']['mom']
    non=metaOS['ost']['non']
    ONSET=metaOS['ost']['ONSET']
    SUSEPTABLE=metaOS['ost']['SUSEPTABLE']
    #UNSUSEPTABLE=metaOS['ost']['UNSUSEPTABLE']
    
    for iDay in range(metaOS['Season_length']):
                
        # Spread       
        rn=np.random.random(metaOS['Shape'])
        #wind=[0.9,1.1,1.25,1.15,1.15,0.75,0.8,0.9]
        #wind=[0.25,0.25,0.25,0.25,1,2,2,2]
        wind=np.ones(8)
        cnt_wind=0
        #np.sum(wind)
        for dx,dy in neighbourhood:
            rn1=rn*wind[cnt_wind]
            z_shift=z['Data'].copy()
            z_shift[mom(dy):non(dy),mom(dx):non(dx)]=z['Data'][mom(-dy):non(-dy),mom(-dx):non(-dx)]
            ind=(z['Data']==SUSEPTABLE) & (z_shift==ONSET) & (rn1<metaOS['Grid_p_Spread'])
            z['Data'][ind]=ONSET
            metaOS['Grid_YearsSinceOcc'][ind]=0
            cnt_wind=cnt_wind+1
            
        # Onset
        rn=np.random.random(metaOS['Shape'])
        ind=(metaOS['Grid_LCC']==1) & (rn<metaOS['p_Onset'])
        z['Data'][ind]=ONSET
        metaOS['Grid_YearsSinceOcc'][ind]=0
            
        # Summarize results    
        z['A_Occurrence'][iYear]=z['A_Occurrence'][iYear]+np.sum(z['Data']==metaOS['ONSET'])
    
    return z
