# bc1ha
## PURPOSE
The **bc1ha** repository supports geospatial analysis in British Columbia’s forest sector at 100m spatial resolution. The repository has two goals: 1) streamline an annual update cycle to download and 
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
- Rasterize insect occurrence from Aerial Overview Survey (PEST_INFESTATION_POLY) (0.1 hour)
- Rasterize harvest occurrence
    - 



 mechanical site prep, knockdown, planting (1 hour)
- Digitize Timber Supply Area (TSA) boundaries if they have changed. (1 min)
- Run scripts that generate derived variables. (~1 hour)

## DERIVED VARIABLES
### Planting Compilation
While it is mandatory to report the spatial geometry of planting for governnment-funded planting post 2017, the spatial geometry for other planting is not necessarily tracked explicitly. To reconstruct planting spatial,
we sequentially compiled it from 1) spatial geometry when reported in the activity layer of RESULTS; 2) the area within the opening that is classified as "Artificial" in the Stocking Type Class" of the Forest Cover Inventory 
layer (when the area classified as Artificial is within 10% of the planting area listed in the RESULTS activity layer); 3) a random draw of areas from the spatial geometry for the opening from the OPENING_SVW layer.

Over 1960-2022, there were two entries in the activity layer of results where the area planted was not reported. These are excluded from analysis.

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

