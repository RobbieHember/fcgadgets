
#%% Import Modules

import bcdata
import os
import subprocess
import numpy as np
import fcgadgets.macgyver.utilities_general as gu

# # get a feature as geojson
# gdb=bcdata.get_data('WHSE_FOREST_VEGETATION.F_OWN',as_gdf=True)

# Check availability
flg=0
if flg==1:
    L=bcdata.list_tables()
    for iL in L:
        if np.char.find(iL,'FOREST_COVER')>=0:
            print(iL)

# Import list of layers to download
d=gu.ReadExcel(r'C:\Users\rhember\Documents\Code_Python\fcgadgets\cbrunner\Parameters\Data_Sources.xlsx')

#L=["WHSE_FOREST_TENURE.FTEN_ROAD_SEGMENT_LINES_SVW"]#,"WHSE_FOREST_TENURE.FTEN_ROAD_SEGMENT_LINES_SVW"]
namL=list(d['Layer Name'])

# Output location
pth_out='C:/Users/rhember/Documents/Data/ForestInventory/20230319'

for nam in namL:
    try:
        if os.path.exists(pth_out + '//' + nam + '.geojson')==True:
            print('Skipping ' + nam)
            continue
        cmd_str='bcdata dump ' + nam + ' > ' + pth_out + '/' + nam + '.geojson"'
        print(nam)
        subprocess.run(cmd_str,shell=True)
    except:
        print(nam + ' not working!')
