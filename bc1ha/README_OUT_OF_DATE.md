# bc1ha
## PURPOSE
The **bc1ha** repository supports geospatial analysis in British Columbia’s forest sector at 100m spatial resolution. The repository has two goals: 1) streamline an annual update cycle to download and 
rasterize information that is made available through BC Data Catalogue and 2) generate derived variables, where varying degrees of manipulation (e.g. gap-filling, consolidation of multiple data sources) 
have been applied to create the input variables used to drive simulation models. Derived variables are often a work in progress and should be used with caution.
<br>
- Command annual updates with **bc1ha_grid**.
- Perform individual processing steps with functions in **bc1ha_util**
- Visualize results with functions in **bc1ha_map_roi**

## ANNUAL UPDATE STEPS
- Download required layers and store as geodatabases for land use, land cover, disturbances, RESULTS, and VRI. The layers that are included should be listed in the Data Sources spreadsheet in the **fcgadgets.cbrunner.parameters** repository. (8+ hours)

- Update look-up-tables (LUTs) using **BuildLUTsFromSourceGDBs**. (7+ hours)

- Simplify some basemaps used for graphing with **SimplifyProvincialGDBs**. (<0.1 hour)

- Mask land area within political border of BC with **GenerateLandMaskBC**. (<0.1 hour)

- Digitize Timber Supply Area (TSA) boundaries if they have changed using **DigitizeTSABoundaries**. (< 0.1 hours)

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

- Rasterize wildfire occurrence (PROT_HISTORICAL_FIRE_POLYS_SP). If needed, rasterize current year fires with **RasterizeWildfireCurrentYear** (0.1 hour)

- Rasterize insect occurrence from Aerial Overview Survey (PEST_INFESTATION_POLY) (0.5 hour)

- Rasterize harvest occurrence:
    - Rasterize harvest year from VEG_CONSOLIDATED_CUT_BLOCKS_SP. (0.2 hours)
    - Rasterize harvest year from NTEMS: In ArcGIS resample NTEMS variables to 100m and clip to BC. Then run **ReprojectDataFromNTEMS**. (0.1 hours)
- Rasterize planting using **RasterizePlanting** (1 hour)
    - Observations of area planted come from the ACTUAL_AREA_PLANTED variable (A_pl), where RESULTS_IND is “Y” and SILV_METHOD_CODE is not “LAYOT”. Over 1960-2022, there were two entries in the activity layer of results where A_pl was not reported. These are excluded from analysis.
    - Where possible, use geometry from RSLT_ACTIVITY_TREATMENT_SVW. For planting with no reported geometry, sequentially estimate it from the area within the opening with STOCKING_TYPE_CLASS  = “Artificial" from the Forest Cover Inventory layer (when the artificial area is within 2% of the planting area listed in the RESULTS activity layer (A_pl). When that fails, find the spatial geometry for the opening from the OPENING_SVW layer and randomly assign a proportion of the opening consistent with A_pl. This is first done with the first opening ID and then the second opening ID variable (OPENING_ID_2). When that fails, find the spatial geometry of the opening from VRI and randomly assign a proportion that is planted consistent with A_pl. 
- Rasterize planting layer using **RasterizePlantingLayer** (0.5 hours). This generates rasters for species codes, percents and genetic worths.
- Rasterize variables from Timber Cruise Compilation:
    - Rasterize percent dead from timber cruise data using **RasterizeCruisePercentDead** (<0.1 hour)
- Run scripts that generate derived variables. (~1 hour)

## DERIVED VARIABLES
### Biogeoclimatic Zone (Gap-Filled)
Gap-fill biogeoclimatic zone for areas not classified as treed.
### Artificial Regeneration Type Compilation 1
Artificial stand establishment consists of mostly planting and a small amount of direct seeding. To understand how planting events affect the carbon balance, events are classified into regeneration types in the Artificial Regeneration Type Compilation 1 (ART-C1). The regeneration types help to understand the context of the planting. The strategy uses information from the silviculture base code (SBC), silviculture technique code (STC), and silviculture method code (SMC) provided in the RESULTS activity layer. It also relies on other variables that describe what preceded each planting event, whether it be a harvest, wildfire, beetles, or a previous planting event. (0.5 hours)
- Back-to-back planting: Sometimes planting projects span multiple years, yet the precise spatial location is sometimes only approximated (aspatially) within the opening (see section on rasterizing planting events). It is also possible that some Fill Planting or Replanting are mistakenly reported as planting, which would be captured in back-to-back planting. The back-to-back planting type is generally excluded from models, which will focus on the first instance.
- Salvage and Planting: Planting events preceded by harvest, which itself was preceded by a catastrophic natural disturbance. For non-obligation planting (as indicated by government-funded funding source codes), a planting event preceded by harvest is always salvage. For licensee planting, the regeneration type is only classified as Salvage and Planting if the percent dead volume (as indicated by Timber Cruise Compilation) exceeds a specified threshold. As we only have access to timber cruises for 2015 on, the ART-C1 likely underestimates the abundance of licensee salvage logging prior to 2015. 
- Straight-to-planting Post Wildfire: Planting events preceded by wildfire. This type is likely mostly “underplanting.”
- Straight-to-planting Post Beetles: Planting events preceded by beetle outbreak with more than light severity. This type may include a high level of mis-classification if beetles happen to be recorded between harvest and planting, or between wildfire and planting.
- Replanting: Replanting indicates that the previously-planted stand failed. These events are classified strictly based on STC = “RP”. 
- Fill Planting: Fill planting indicates efforts to fill-in any areas not meeting a specific regeneration standard. These events are classified strictly based on STC = “FP”. 
- Harvest and Planting NSR Backlog: Planting by government programs in areas harvested prior to 1987.
- Harvest and Planting: Planting in response to obligation to reforest following commercial harvesting.
- Road Rehabilitation: Planting with STC = “RR”.
- Direct Seeding: Direct seeding indicates that the stand was established by seeding. These events are classified strictly based on SBC = “DS”.

### Non-obligation Stand Establishment Mask
- Mask non-obligation stand establishment with **Planting_NonOb_Mask**. (<0.1 hours)

### Biogeoclimatic Zone (Gap-Filled)-Natural Disturbance Type Combination
- A combination of BGC zone (using the gap-filled version) and natural disturbance type. Based on **DeriveBGCZoneNDTCombo**. (<0.1 hours)

### Tree Density Class Consolidated 1
- Generation of a consolidated variable describing tree density class. This dataset is a revised version of the Level 5 BC Land Surface Classification Scheme from VRI. It is created using **ConsolidateTreeDensityClass**.

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

