# bc1ha
## PURPOSE
The **bc1ha** repository supports geospatial analysis in British Columbia’s forest sector at 100m spatial resolution. The repository has two goals: 1) streamline an annual update cycle to download and 
rasterize information that is made available through BC Data Catalogue and 2) generate derived variables, where varying degrees of manipulation (e.g. gap-filling, consolidation of multiple source data sources) 
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

- Update look-up-tables (LUTs). (7+ hours)

- Rasterize VRI. VRI is too big to be rasterized in Python (wtih 32 GB RAM). In ArcGIS, rasterize the feature ID or the Forest Cover ID. Convert to TIFF. Export to local machine. 
Use ClipToRaster_ByFile to standardize grid extent. Use custom function from bc1ha_util.py specifically designed to rasterize VRI, using the feature ID as input. (10 hours)

- Rasterize Forest Cover Inventory (FCI). FCI is too big to be rasterized in Python (wtih 32 GB RAM). In ArcGIS, rasterize the feature ID or the Forest Cover ID. Convert to TIFF. Export to local machine. 
Use ClipToRaster_ByFile to standardize grid extent. Use custom function from bc1ha_util.py specifically designed to rasterize FCI, using the feature ID as input. (1 hours)

- Other time-independent variables (that are not too big for Python), use RasterizeFromSource function to rasterize variables. (20 min)

- Rasterize time-dependent variables, including wildfire, beetles, defoliators, harvest, mechanical site prep, knockdown, planting (1 hour)

- Digitize Timber Supply Area (TSA) boundaries if they have changed. (1 min)

- Run scripts that generate derived variables. (~1 hour)

## LIST OF DERIVED VARIABLES
Additions forthcoming.

## SPARSE SUBSAMPLES
Prepare sparse inputs at common spatial subsampling resolutions, including 1km, 5km, 10km and 20km. When models are run at these resolutions, **fcgadgets** will automatically draw on these sparse inputs, which will 
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

