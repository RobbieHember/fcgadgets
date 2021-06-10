
#%% Import modules

import ipyleaflet as ipyl
import ipywidgets as ipyw
import geopandas as gpd
import fiona
import json
import pyproj
import numpy as np
from ipywidgets import Layout, HBox, VBox, FloatText
from shapely.geometry import Polygon,Point
from fcgadgets.macgyver import utilities_general as gu
#from fcgadgets.macgyver import utilities_gis as gis

#%%

def update_html(feature,**kwargs):
    html1.value='''<h4>OPENING_ID</h4>{}<h4>FIA_PROJECT_ID</h4>{}<h4>SILV_BASE_CODE</h4>{}
    '''.format(feature['properties']['OPENING_ID'],feature['properties']['FIA_PROJECT_ID'],feature['properties']['SILV_BASE_CODE'])

#%%
    
def Add_Project(m,id,gdf):
    
    op=[]; 
    op_json=[]; 
    op_geojson=[]
    for iOp in range(len(id)):
        op.append(gdf[gdf.FIA_PROJECT_ID==id[iOp]])
        op[iOp]=op[iOp].reset_index()
        op_json.append(op[iOp].geometry.__geo_interface__)
        for i in range(len(op_json[iOp]['features'])):
            for key in op[iOp].columns:
                if key=='geometry': 
                    continue
                op_json[iOp]['features'][i]['properties'][key]=op[iOp].loc[i,key]
        
        op_geojson.append(ipyl.GeoJSON(data=op_json[iOp],style={'color':'yellow','weight':5,'dashArray':'1','opacity':1,
                'fillColor':'yellow','fillOpacity':0.25},
                name='FIA_PROJECT_ID:' + str(id) ))
        
        m.add_layer(op_geojson[iOp])
        
        # Does not work inside external function
        #op_geojson[iOp].on_hover(update_html)
    
    return op,op_json,op_geojson

#%% Add opening spatial
    
def Add_Openings(m,id,gdf):
    op=[]; 
    op_json=[]; 
    op_geojson=[]
    for iOp in range(len(id)):
        op.append(gdf[gdf.OPENING_ID==id[iOp]])
        op[iOp]=op[iOp].reset_index()
        op_json.append(op[iOp].geometry.__geo_interface__)
        for i in range(len(op_json[iOp]['features'])):
            for key in op[iOp].columns:
                if key=='geometry': continue
                op_json[iOp]['features'][i]['properties'][key]=op[iOp].loc[i,key]
        
        op_geojson.append(ipyl.GeoJSON(data=op_json[iOp],
            style={'color':'yellow','weight':3,'dashArray':'1','opacity':1,
            'fillColor':'yellow','fillOpacity':0.04},name='Opening:' + str(id) ))
        
        m.add_layer(op_geojson[iOp])
    
    return op,op_json,op_geojson

#%% Add forect cover for opening
    
def Add_ForestCover_ForOpening(m,id,gdf):
    op=[]; 
    op_json=[]; 
    op_geojson=[]
    for iOp in range(len(id)):
        op.append(gdf[gdf.OPENING_ID==id[iOp]])
        op[iOp]=op[iOp].reset_index()
        op_json.append(op[iOp].geometry.__geo_interface__)
        for i in range(len(op_json[iOp]['features'])):
            for key in op[iOp].columns:
                if key=='geometry': continue
                op_json[iOp]['features'][i]['properties'][key]=op[iOp].loc[i,key]
        
        op_geojson.append(ipyl.GeoJSON(data=op_json[iOp],
            style={'color':'cyan','weight':1,'dashArray':'3','opacity':1,
            'fillColor':None,'fillOpacity':0},name='Opening:' + str(id) ))
        
        m.add_layer(op_geojson[iOp])
    
    return op,op_json,op_geojson

#%% Add forect cover for opening
    
def Add_SXY_ForOpening(m,id,gdf,gdf_sxy):
    
    for iOp in range(len(id)):
        gdf_sxy0=gdf_sxy[gdf_sxy.OPENING_ID==id[iOp]]    
        mc=ipyl.MarkerCluster(markers=[ipyl.Marker(location=geolocation.coords[0][::-1]) for geolocation in gdf_sxy0.geometry])
        m.add_layer(mc);
    
    return m

#%% 
        
def Add_Activity(m,at,gdf,ivl):   
    
    if at=='Fertilization Aerial':
        fa_gdf=gdf[(gdf.SILV_BASE_CODE=='FE') & (gdf.SILV_TECHNIQUE_CODE=='CA') & (gdf.SILV_METHOD_CODE=='HELI')]
    elif at=='Planting':
        fa_gdf=gdf[(gdf.SILV_BASE_CODE=='PL') & (gdf.SILV_METHOD_CODE!='LAYOT')]
        
    fa_gdf=fa_gdf.reset_index()
    fa_json=fa_gdf.geometry.__geo_interface__    
    for i in range(len(fa_json['features'])):
        for key in fa_gdf.columns:
            if key=='geometry': continue
            fa_json['features'][i]['properties'][key]=fa_gdf.loc[i,key]
    fa_geojson=ipyl.GeoJSON(data=fa_json,style={'color':'red','weight':2,'dashArray':'1','opacity':1,'fillColor':'red','fillOpacity':0.25},name='Fertilization Aerial')
    m.add_layer(fa_geojson)
    #fa_geojson.on_hover(update_html)
    
    return fa_geojson
    
#    # Exctract points that intersect the selected openings
#    id=fa_gdf['OPENING_ID'].unique()
#    pop=np.array([])
#    for i in range(0,len(id),ivl):
#        ind=np.where(fa_gdf['OPENING_ID']==id[i])[0]
#        pop=np.append(pop,fa_gdf.loc[ind,'ID_atu_pol'].values)
#    ind=np.array([])
#    for i in range(pop.size):
#        ind0=np.where(sxy['ID_atu_pol'].values==pop[i])[0]
#        ind=np.append(ind,ind0)
#    sxy0=sxy.iloc[ind]
#    sxy0=sxy0.reset_index()
#    
#    # Convert subset of intersecting grid cells to json and add to map
#    y_json=sxy0.geometry.__geo_interface__
#    for i in range(len(y_json['features'])):
#        for key in sxy0.columns:
#            if key=='geometry': continue
#            y_json['features'][i]['properties'][key]=sxy0.loc[i,key]
#
#    for i in range(len(y_json['features'])):
#        loc=(y_json['features'][i]['geometry']['coordinates'][1],y_json['features'][i]['geometry']['coordinates'][0])
#        circle_marker=ipyl.CircleMarker(location=loc,color="yellow",face_color=None,radius=2)
#        m.add_layer(circle_marker)
#    
#    fa_geojson.on_hover(update_html)
#    
#    return
    
