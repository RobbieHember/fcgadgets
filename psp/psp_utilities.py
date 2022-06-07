
#%% Import modules

import numpy as np
from scipy.io import loadmat

#%% Read tree level data by species

def Read_L3_TL_BySpc(spp):
    
    #spp=['SW']
    
    pthin=r'E:\Data\ForestInventory\PSP-NADB\Data\03_IntegratedDatabases\TreeLevel\BySpecies'
    
    for iS in range(len(spp)):
        
        z=loadmat(pthin + '\\FI_PSP_L3_TL_BySpc_' + spp[iS] + '.mat')
    
        d={}
        for iV in range(z['z']['Name'][0].size):
            d[z['z']['Name'][0][iV][0]]=z['z']['Data'][0][iV].flatten().astype(float)*z['z']['ScaleFactor'][0][iV][0][0]
        
    #--------------------------------------------------------------------------   
    # Fix stand age variable in BC
    #--------------------------------------------------------------------------
    
#    ind0=np.where(d['ID_DB']==1)[0]
#    
#    ID_Plot=d['ID_Plot'][ind0]
#    ID_Tree=d['ID_Tree'][ind0]
#    Age_t0_Stand=d['Age_t0_Stand'][ind0]
#    Age_t1_Stand=d['Age_t1_Stand'][ind0]
#    DT=d['DT'][ind0]
#    
#    u=np.unique( np.column_stack((ID_Plot,ID_Tree)),axis=1 )
#
#    for i in range(u.shape[0]):
#        ind=np.where( (ID_Plot==u[i,0]) & (ID_Tree==u[i,1]) )[0]
#        if ind.size==0:
#            continue
#        if ind.shape[0]==1:
#            continue
#        Age_t0_Stand[ind]=Age_t0_Stand[ind[0]]+np.append(0,np.cumsum(DT[ind[0:-1]]))
#        Age_t1_Stand[ind]=Age_t0_Stand[ind[0]]+np.cumsum(DT[ind])
#    
#    d['Age_t0_Stand'][ind0]=Age_t0_Stand
#    d['Age_t1_Stand'][ind0]=Age_t1_Stand        
        
    return d