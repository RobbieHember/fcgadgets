'''

'''

#%% IMPORT MODULES

import sys
import numpy as np
from osgeo import gdal
from osgeo import osr
import matplotlib.pyplot as plt
import geopandas as gpd
from rasterio.transform import from_origin
import pandas as pd
import fiona
import rasterio
from rasterio import features
from shapely.geometry import Point, Polygon,box
from shapely import geometry

import fcgadgets.macgyver.utilities_general as gu
import fcgadgets.macgyver.utilities_gis as gis
from fcgadgets.macgyver import utilities_inventory as invu

#import contextily as ctx
from shapely.ops import cascaded_union
from geovoronoi.plotting import subplot_for_map, plot_voronoi_polys_with_points_in_area
from geovoronoi import voronoi_regions_from_coords, points_to_coords

#%% Import data

gdf_bc_bound=gpd.read_file(r'C:\Users\rhember\Documents\Data\Basemaps\Basemaps.gdb',layer='NRC_POLITICAL_BOUNDARIES_1M_SP')
gdf_tpf=gpd.read_file(r'C:\Users\rhember\Documents\Data\ForestInventory\Infrastructure.gdb',layer='GSR_TMBR_PRCSSING_FAC_SV')

zS=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Admin\FAIB_Standard.tif')

#%%

# Isolate lumber mills
#coord_mills=points_to_coords(gdf_tpf.geometry)

# Area
Area=np.zeros(len(gdf_bc_bound))
for i in range(len(gdf_bc_bound)):
    Area[i]=gdf_bc_bound.iloc[i].geometry.area
isort=np.flip(np.argsort(Area))

#%%

#box(W, S, E, N)
#bnd=box(zS['xmin'],zS['ymin'],zS['xmax'],zS['ymax'])
bnd=cascaded_union(gdf_bc_bound.iloc[isort[0]].geometry)

ikp0=np.where( (gdf_tpf['PRODUCT_CODE']=='LBR') )[0]
ikp=np.where( (gdf_tpf['PRODUCT_CODE']=='LBR') & (gdf_tpf.within(bnd)==True) )[0]
xy=points_to_coords(gdf_tpf.loc[ikp,'geometry'])

poly,pts=voronoi_regions_from_coords(xy,bnd)

gdf_catch=gpd.GeoDataFrame({'geometry':poly,'ID':1},crs=gdf_bc_bound.crs)

plt.close('all')
fig,ax=plt.subplots(1,figsize=gu.cm2inch(10,10))
gdf_bc_bound.plot(ax=ax)
gdf_catch.plot(ax=ax,edgecolor='k',facecolor='none',linewidth=0.5)
gdf_tpf.iloc[ikp0].plot(ax=ax,marker='o',markersize=5,facecolor='y')
