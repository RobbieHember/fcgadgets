# bc1ha
## INTRODUCTION
The **bc1ha** repository supports geospatial analysis in British Columbia’s forest sector at 100m spatial resolution.
<br>
<br>
## WORKFLOW
The scripts are focused on streamlining an annual update cycle. The bc1ha_grid.py is a script that runs the annual update steps. The bc1ha_util stores supporting functions.

Step 1: Download required layers and store as geodatabases for land use, land cover, disturbances, RESULTS, and VRI. (8+ hours).

Step 2: Update look-up-tables (LUTs). (7+ hours). 

Step 3: VRI and Forest Cover variables are too big to be rasterized in Python (wtih 32 GB RAM). In ArcGIS, rasterize the feature ID or the Forest Cover ID. Convert to TIFF. Export to local machine. Use ClipToRaster_ByFile to standardize grid extent.

Step 4: Digitize Timber Supply Area (TSA) boundaries if they have changed.

Step 5: Use the 

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

