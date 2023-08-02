# bc1ha
## PURPOSE
The **bc1ha** repository supports geospatial analysis in British Columbia’s forest sector at 100m spatial resolution. The repository has two goals: 1) streamline an annual update cycle to download and rasterize information that is made available through BC Data Catalogue; 2) produce derived variables through cleaning, consolidating, smoothing, massaging, "moving forward, not backward; upward, not forward; and always twirling, twirling, twirling towards" preparation of high-quality model input variables. 
<br><br>
Derived variables are often a work in progress and should be used with caution.
<br><br>
bc1ha_grid.py: Command script that runs the annual update steps.
<br>
bc1ha_util.py: Supporting functions
<br>
bc1ha_map_roi: Graphics by region of interest

## ANNUAL UPDATE STEPS
Step 1: Download required layers and store as geodatabases for land use, land cover, disturbances, RESULTS, and VRI. The layers that are included should be listed in the Data Sources spreadsheet in the **fcgadgets.cbrunner.parameters** repository. (8+ hours)

Step 2: Update look-up-tables (LUTs). (7+ hours)

Step 3: VRI and Forest Cover variables are too big to be rasterized in Python (wtih 32 GB RAM). In ArcGIS, rasterize the feature ID or the Forest Cover ID. Convert to TIFF. Export to local machine. Use ClipToRaster_ByFile to standardize grid extent. Use custom functions from bc1ha_util.py, with the feature ID as input, to rasterize the necessary variables. (10 hours)

Step 4: Digitize Timber Supply Area (TSA) boundaries if they have changed. (1 min)

Step 5: Run scripts that generate derived variables. (~1 hour)

## DERIVED VARIABLES


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

