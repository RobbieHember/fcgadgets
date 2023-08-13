# bc1ha
## PURPOSE
The **bc1ha** repository supports geospatial analysis in British Columbiaâ€™s forest sector at 100m spatial resolution. The repository has two goals: 1) streamline an annual update cycle to download and 
rasterize information that is made available through BC Data Catalogue and 2) generate derived variables, where varying degrees of manipulation (e.g. gap-filling, consolidation of multiple data sources) 
have been applied to create the input variables used to drive simulation models. 
<br><br>
Derived variables are often a work in progress and should be used with caution.
<br><br>
bc1ha_grid.py: Command script that runs the annual update steps.
<br>
bc1ha_util.py: Supporting functions
<br>
bc1ha_map_roi: Graphics by region of interest

## ANNUAL UPDATE STEPS
- Download required layers and store as geodatabases for land use, land cover, disturbances, RESULTS, and VRI. The layers that are included should be listed in the Data Sources spreadsheet in the **fcgadgets.cbrunner.parameters** repository. (8+ hours)

- Update look-up-tables (LUTs) using **BuildLUTsFromSourceGDBs**. (7+ hours)

- Rasterize VRI. VRI is too big (32GB RAM) to be rasterized in Python. (10 hours)
    - In ArcGIS, rasterize the feature ID. Convert to TIFF. 
    - Export to local machine. Use **ClipToRaster_ByFile** to standardize grid extent.
    - Use custom function from bc1ha_util.py specifically designed to rasterize VRI, using the feature ID as input. 
    - The **CreateIdForCategoricalVariable** function generates numerical IDs for a categorical variable in GDBs using the LUTs. 

- Rasterize variables from the Forest Cover Inventory (FCI) layer. FCI is too big (32GB RAM) to be rasterized in Python. 
    - In ArcGIS, rasterize Forest Cover ID. 
    - Convert to TIFF. 
    - Export to local machine. 
    - Use **ClipToRaster_ByFile** to standardize grid extent. 
    - Use **RasterizeForestCoverInventory** to rasterize attributes using the raster feature ID as key. (1 hours)

- Rasterize other time-independent variables (that are not too big for Python) use **RasterizeFromSource** function to rasterize most variables.(0.5 hours)

- Rasterize opening ID from RESULTS OPENING layer. There is spatial overlap of openings going back in time so **RasterizeOpeningID** will generate two raster variables OPENING_ID_1 and OPENING_ID_2.

- Rasterize wildfire occurrence (PROT_HISTORICAL_FIRE_POLYS_SP) (0.1 hour)

- Rasterize insect occurrence from Aerial Overview Survey (PEST_INFESTATION_POLY) (0.5 hour)

- Rasterize harvest occurrence
    - VEG_CONSOLIDATED_CUT_BLOCKS_SP
    - NTEMS: In ArcGIS resample NTEMS variables to 100m and clip to BC. Then run **ReprojectDataFromNTEMS**.

- Rasterize planting using **RasterizePlanting** (1 hour)
    - Observations of area planted come from the ACTUAL_AREA_PLANTED variable (A_pl), where RESULTS_IND is “Y” and SILV_METHOD_CODE is not “LAYOT”. Over 1960-2022, there were two entries in the activity layer of results where A_pl was not reported. These are excluded from analysis.
    - Where possible, use geometry from RSLT_ACTIVITY_TREATMENT_SVW. For planting with no reported geometry, sequentially estimate it from the area within the opening with STOCKING_TYPE_CLASS  = “Artificial" from the Forest Cover Inventory layer (when the artificial area is within 2% of the planting area listed in the RESULTS activity layer (A_pl). When that fails, find the spatial geometry for the opening from the OPENING_SVW layer and randomly assign a proportion of the opening consistent with A_pl. This is first done with the first opening ID and then the second opening ID variable (OPENING_ID_2). When that fails, find the spatial geometry of the opening from VRI and randomly assign a proportion that is planted consistent with A_pl. 

- Rasterize planting layer using **RasterizePlantingLayer** (0.5 hours). This generates rasters for species codes, percents and genetic worths.

- Digitize Timber Supply Area (TSA) boundaries if they have changed. (< 0.1 hours)

- Run scripts that generate derived variables. (~1 hour)

## DERIVED VARIABLES
### Regeneration Type Compilation: Representing artificial stand establishment, including planting and direct seeding, is aided by the silviculture base code (SBC), silviculture technique code (STC), and silviculture method code (SMC) provided in the RESULTS activity layer. The Regeneration Type Compilation 1 (RTC1) re-classifies stand establishment events into regeneration types to facilitate use in models. (0.5 hours)
- Back-to-back planting: Sometimes planting projects span multiple years, yet the precise spatial location is sometimes only approximated (aspatially) within the opening (see section on rasterizing planting events). By include a class for back-to-back planting, models can exclude the excessive number of events, focusing on the first or last instance.
- 

###

## SPARSE SUBSAMPLES
Prepare sparse inputs at commonly-used spatial subsampling resolutions, including 1km, 5km, 10km and 20km. When models are run at these resolutions, **fcgadgets** will automatically draw on these sparse inputs, which will 
lower computing time during project development.

## License

    Copyright 2020 Province of British Columbia

    Licensed under the Apache License, Version 2.0 (the "License");
    you may not use this file except in compliance with the License.
    You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS,
    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and limitations under the License.

------------------------------------------------------------------------

