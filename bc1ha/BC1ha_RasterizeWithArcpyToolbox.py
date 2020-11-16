
#=============================================================================

# Run instructions:

# 1a) Move geodatabase files to Table of Contents, then export to shapefile and save in ArcGIS folder (not in Default geodatabase).
#       BEC.BIOGEOCLIMAITC_POLY --> bec.shp
#       VEG_COMP_POLY --> vegco.shp
#
# OR
#
# 1b) Read from previously prepared geodatabase (e.g., VRI.gdb) (appears to take way longer)


# 2) Open command prompt and run: E:\sw_nt\Python27\ArcGIS10.3\python.exe H:\GIS\RasterizeTo100haGrid.py

# Problems:
# None reported	

#=============================================================================

#-----------------------------------------------------------------------------
# Preliminary stuff
#-----------------------------------------------------------------------------

# Import modules
import arcpy as ap

# Set workspace
ap.env.workspace="C:/Users/rhember/Documents/ArcGIS"

# Cell size (m)
cc=100

PathOut="C:/Users/rhember/Documents/ArcGIS/BC1ha"

#-----------------------------------------------------------------------------
# Rasterize TSA layers
#-----------------------------------------------------------------------------

#fin="tsa.shp"

ln_in=""; ln_out="tsa"
#ap.PolygonToRaster_conversion(fin,ln_in,ln_out,"CELL_CENTER","None",cc)
#ap.RasterToOtherFormat_conversion(ln_out,PathOut,"TIFF")

#-----------------------------------------------------------------------------
# Rasterize BEC layers
#-----------------------------------------------------------------------------

fin="bec.shp"

ln_in="ZONE"; ln_out="becz"
#ap.PolygonToRaster_conversion(fin,ln_in,ln_out,"CELL_CENTER","None",cc)
#ap.RasterToOtherFormat_conversion(ln_out,PathOut,"TIFF")


ln_in="SUBZONE"; ln_out="becsz"
#ap.PolygonToRaster_conversion("bec.shp",ln_in,ln_out,"CELL_CENTER","None",cc)
#ap.RasterToOtherFormat_conversion(ln_out,PathOut,"TIFF")

ln_in="VARIANT"; ln_out="becv"
#ap.PolygonToRaster_conversion("bec.shp",ln_in,ln_out,"CELL_CENTER","None",cc)
#ap.RasterToOtherFormat_conversion(ln_out,PathOut,"TIFF")

#-----------------------------------------------------------------------------
# Rasterize VRI layers
# Source: VEG_COMP_LYR_R1_POLY
#-----------------------------------------------------------------------------

fin="C:/Users/rhember/Documents/ArcGIS/VRI.gdb/VEG_COMP_LYR_R1_POLY"

ln_in="REFERENCE_YEAR"; ln_out="refyr"
#ap.FeatureToRaster_conversion(fin,ln_in,ln_out,cc)
#ap.RasterToOtherFormat_conversion(ln_out,PathOut,"TIFF")

ln_in="LINE_7B_DISTURBANCE_HISTORY"; ln_out="histdist"
#ap.PolygonToRaster_conversion(fin,ln_in,ln_out,"CELL_CENTER","None",cc)
#ap.RasterToOtherFormat_conversion(ln_out,PathOut,"TIFF")

ln_in="LINE_8_PLANTING_HISTORY"; ln_out="histplant"
#ap.PolygonToRaster_conversion(fin,ln_in,ln_out,"CELL_CENTER","None",cc)
#ap.RasterToOtherFormat_conversion(ln_out,PathOut,"TIFF")

#ln_in="EARLIEST_NONLOGGING_DIST_TYPE"; ln_out="enldistt"
#ap.PolygonToRaster_conversion(fin,ln_in,ln_out,"CELL_CENTER","None",cc)
#ap.RasterToOtherFormat_conversion(ln_out,PathOut,"TIFF")

# *** This crashes in Arc for some reason ***
#ln_in="EARLIEST_NONLOGGING_DIST_DATE"; ln_out="enldiyr"
#ap.PolygonToRaster_conversion(fin,ln_in,ln_out,"CELL_CENTER","None",cc)
#ap.RasterToOtherFormat_conversion(ln_out,PathOut,"TIFF")

ln_in="STAND_PERCENTAGE_DEAD"; ln_out="pdead"
#ap.PolygonToRaster_conversion(fin,ln_in,ln_out,"CELL_CENTER","None",cc)
#ap.RasterToOtherFormat_conversion(ln_out,PathOut,"TIFF")

# Skipped
#ln_in="FREE_TO_GROW_IND"; ln_out="freetogro"
#ap.PolygonToRaster_conversion(fin,ln_in,ln_out,"CELL_CENTER","None",cc)
#ap.RasterToOtherFormat_conversion(ln_out,PathOut,"TIFF")

ln_in="SITE_INDEX"; ln_out="si"
#ap.PolygonToRaster_conversion(fin,ln_in,ln_out,"CELL_CENTER","None",cc)
#ap.RasterToOtherFormat_conversion(ln_out,PathOut,"TIFF")

ln_in="VRI_LIVE_STEMS_PER_HA"; ln_out="sphlive"
#ap.PolygonToRaster_conversion(fin,ln_in,ln_out,"CELL_CENTER","None",cc)
#ap.RasterToOtherFormat_conversion(ln_out,PathOut,"TIFF")

ln_in="VRI_DEAD_STEMS_PER_HA"; ln_out="sphdead"
#ap.PolygonToRaster_conversion(fin,ln_in,ln_out,"CELL_CENTER","None",cc)
#ap.RasterToOtherFormat_conversion(ln_out,PathOut,"TIFF")

ln_in="SPECIES_CD_1"; ln_out="spc_cd1"
#ap.PolygonToRaster_conversion(fin,ln_in,ln_out,"CELL_CENTER","None",cc)
#ap.RasterToOtherFormat_conversion(ln_out,PathOut,"TIFF")

ln_in="SPECIES_PCT_1"; ln_out="spc_p1"
#ap.PolygonToRaster_conversion(fin,ln_in,ln_out,"CELL_CENTER","None",cc)
#ap.RasterToOtherFormat_conversion(ln_out,PathOut,"TIFF")

ln_in="SPECIES_CD_2"; ln_out="spc_cd2"
#ap.PolygonToRaster_conversion(fin,ln_in,ln_out,"CELL_CENTER","None",cc)
#ap.RasterToOtherFormat_conversion(ln_out,PathOut,"TIFF")

ln_in="SPECIES_PCT_2"; ln_out="spc_p2"
#ap.PolygonToRaster_conversion(fin,ln_in,ln_out,"CELL_CENTER","None",cc)
#ap.RasterToOtherFormat_conversion(ln_out,PathOut,"TIFF")

ln_in="SPECIES_CD_3"; ln_out="spc_cd3"
#ap.PolygonToRaster_conversion(fin,ln_in,ln_out,"CELL_CENTER","None",cc)
#ap.RasterToOtherFormat_conversion(ln_out,PathOut,"TIFF")

ln_in="SPECIES_PCT_3"; ln_out="spc_p3"
#ap.PolygonToRaster_conversion(fin,ln_in,ln_out,"CELL_CENTER","None",cc)
#ap.RasterToOtherFormat_conversion(ln_out,PathOut,"TIFF")

ln_in="PROJ_AGE_1"; ln_out="age1"
#ap.PolygonToRaster_conversion(fin,ln_in,ln_out,"CELL_CENTER","None",cc)
#ap.RasterToOtherFormat_conversion(ln_out,PathOut,"TIFF")

ln_in="PROJ_HEIGHT_1"; ln_out="h1"
#ap.PolygonToRaster_conversion(fin,ln_in,ln_out,"CELL_CENTER","None",cc)
#ap.RasterToOtherFormat_conversion(ln_out,PathOut,"TIFF")

ln_in="LIVE_VOLUME_125cm"; ln_out="vlive125"
#ap.PolygonToRaster_conversion(fin,ln_in,ln_out,"CELL_CENTER","None",cc)
#ap.RasterToOtherFormat_conversion(ln_out,PathOut,"TIFF")


#-----------------------------------------------------------------------------
# Rasterize burn severity
# *** No longer doing it this way because it will not capture overlapping polygons ***
# Source: WHSE_FOREST_VEGETATION.VEG_BURN_SEVERITY_SP
#-----------------------------------------------------------------------------

#fin="C:/Users/rhember/Documents/ArcGIS/"
#ln_in=""; ln_out="burnsev_class"
#ap.FeatureToRaster_conversion(fin,ln_in,ln_out,cc)
#ap.RasterToOtherFormat_conversion(ln_out,PathOut,"TIFF")

#ln_in=""; ln_out="burnsev_yr"
#ap.FeatureToRaster_conversion(fin,ln_in,ln_out,cc)
#ap.RasterToOtherFormat_conversion(ln_out,PathOut,"TIFF")

#-----------------------------------------------------------------------------
# Rasterize wildfire year
# *** No longer doing it this way because it will not capture overlapping polygons ***
# Source: WHSE_FOREST_VEGETATION.VEG_BURN_SEVERITY_SP
#-----------------------------------------------------------------------------

#fin="C:/Users/rhember/Documents/ArcGIS/Disturbances.gdb/PROT_HISTORICAL_FIRE_POLYS_SP"

#ln_in="FIRE_YEAR"; ln_out="fire_yr"
#ap.FeatureToRaster_conversion(fin,ln_in,ln_out,cc)
#ap.RasterToOtherFormat_conversion(ln_out,PathOut,"TIFF")

#fin="C:/Users/rhember/Documents/ArcGIS/Disturbances.gdb/PROT_CURRENT_FIRE_POLYS_SP"
#ln_in="FIRE_YEAR"; ln_out="fire_yr_c"
#ap.FeatureToRaster_conversion(fin,ln_in,ln_out,cc)
#ap.RasterToOtherFormat_conversion(ln_out,PathOut,"TIFF")


#-----------------------------------------------------------------------------
# Template:
#-----------------------------------------------------------------------------

#ln_in=""; ln_out=""
#ap.PolygonToRaster_conversion(fin,ln_in,ln_out,"CELL_CENTER","None",cc)
#ap.RasterToOtherFormat_conversion(ln_out,PathOut,"TIFF")
