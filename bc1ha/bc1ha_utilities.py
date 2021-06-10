
#%% Import modules

import numpy as np
import matplotlib.pyplot as plt
import fcgadgets.macgyver.utilities_general as gu
import fcgadgets.macgyver.utilities_gis as gis

#%% Get normal climate data

def GetClimateNormalFromGeoTiff(x,y,pthin,varl):

    d={}
    d['x']=x
    d['y']=y
    
    for iv in range(len(varl)):
    
        # Specify file name for each variable code and scale factor
        if varl[iv]=='tmean_gs':
            fin=pthin + '\\BC1ha_tmean_gs_norm_1971to2000_si_hist_v1.tif'
            sf=0.1
        elif varl[iv]=='tmin_ann':
            fin=pthin + '\\BC1ha_tmin_ann_norm_1971to2000_si_hist_v1.tif'
            sf=0.1
        elif varl[iv]=='etp_gs':
            fin=pthin + '\\BC1ha_etp_gs_norm_1971to2000_comp_hist_v1.tif'
            sf=0.1
        elif varl[iv]=='eta_gs':
            fin=pthin + '\\BC1ha_eta_gs_norm_1971to2000_comp_hist_v1.tif'
            sf=0.1            
        elif varl[iv]=='ws_gs':
            fin=pthin + '\\BC1ha_ws_gs_norm_1971to2000_comp_hist_v1.tif'
            sf=1
        
        # Load geotiff
        z=gis.OpenGeoTiff(fin)
    
        if iv==0:
            # Extract vector x and y
            xv=z['X'][0,:]
            yv=z['Y'][:,0]
        
            # Generate indices to x and y
            ix=np.zeros(x.size,dtype='int32')
            iy=np.zeros(x.size,dtype='int32')
            for i in range(x.size):
                dx=np.abs(xv-x[i])
                dy=np.abs(yv-y[i])
                ix[i]=np.where(dx==np.min(dx))[0][0]
                iy[i]=np.where(dy==np.min(dy))[0][0]
        
        # Extract nearest grid cell data
        d[varl[iv]]=np.nan*np.ones(x.size)
        for i in range(x.size):
            d[varl[iv]][i]=z['Data'][iy[i],ix[i]]

        # Covernt to float and apply scale factor
        d[varl[iv]]=d[varl[iv]].astype('float')*sf
    
    return d