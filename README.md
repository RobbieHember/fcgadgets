* <a id="devex-badge" rel="Exploration" href="https://github.com/BCDevExchange/assets/blob/master/README.md"><img alt="Being designed and built, but in the lab. May change, disappear, or be buggy." style="border-width:0" src="https://assets.bcdevexchange.org/images/badges/exploration.svg" title="Being designed and built, but in the lab. May change, disappear, or be buggy." /></a>

# fcgadgets
## Features
The fcgadgets package supports greenhouse gas (GHG) balance estimation, accounting, and reporting for climate change mitigation projects in British Columbia’s Land Use, Land Use Change and Forestry (LULUCF) sector. It is a community-based computer simulation model for estimating the GHG balance of treatment areas under various management scenarios. The fcgadgets.utilities module stores script templates and custom utilities that support application of cbrunner. 
cbrunner is written in the Python 3/Jupyter environment, benefiting from stable integrated libraries for simulation modelling, geographical information systems, data analytics, and application deployment (Downey, 2017). 

The repository was designed for intermediate and advanced analysts who value: 
* Transparent methods;
* Efficient workflow through custom links with BC’s spatial forest inventories and growth and yield models;
* Modular representation of biophysical processes;
* Equivalent principles and standards applied in Canada’s National GHG Inventory.

As examples, the demos subpackage walks through demonstrations of Forest Carbon Initiative (FCI) project types:
* Aerial fertilization
* Underplanting fire-impacted forests
* Salvage logging beetle-impacted forests under varying fibre utilization intensities

## Workflow
There are four ways to apply cbrunner depending on the nature of the desired project. Small projects – with fewer than 1,500 combinations of locations or scenarios – can be run from a Jupyter Notebook. The work simply involves populating two Excel spreadsheets with the input variables and parameters. Bigger projects are scripted in Python and can adopt existing templates for projects that focus on running simulations at point locations, or across scattered polygons, or across continuous regular grids.
![image info](./images/fcgadgets_runoptions.png)

## Core modules used in cbrunner
Biophysical processes in cbrunner are represented by a variety of Python methods within the fcgadgets.cbrunner.cbrun_annproc module. The cbrunner module can also call modules specifically developed to simulate natural disturbances, including fcgadgets.taz and forest management activities, including fcgadgets.activities.
Biomass dynamics:
* BiomassFromTASSorTIPSY: Integration with the TASS/TIPSY growth and yield software application (https://www2.gov.bc.ca/gov/content/industry/forestry/managing-our-forest-resources/forest-inventory/growth-and-yield-modelling)
* BiomassFromSawtooth:* Integration with an individual-tree, distance-independent model of biomass dynamics (Hember et al., 2019; Hember and Kurz, 2018).
Dead organic matter dynamics:
* DOM_From_CBM08: Representation of carbon cycling through snags, coarse woody debris, organic soil horizon, and mineral soil horizon following methods similar to the Carbon Budget Model of the Canadian Forest Sector (Kurz et al., 2009).
Natural disturbance and management events
* DisturbanceAndManagement: A custom Python method designed to represent wildfire, insects, disease, and management treatments.
Harvested wood products:
* HWP_From_Dymond12: Representation of GHG balance for fibre that is removed from forest ecosystems. This module aims to capture the dynamics described by the BC Harvested Wood Products model version 1 (Dymond, 2012).

## Additional modules
### Management activities
* nutrient_addition: Representation of GHG balance responses to aerial applications of Urea.
### Natural disturbances
* taz: Statistical models that represent natural disturbances.

## References
Downey, A.B., 2017. Modeling and Simulation in Python – Green Tea Press, 2.3. ed. Green Tea Press, Needham, Massaschusetts.

Dymond, C.C., 2012. Forest carbon in North America: annual storage and emissions from British Columbia’s harvest, 1965-2065. Carbon Balance and Management 7, (24 July 2012)-(24 July 2012).

Hember, R.A., Kurz, W.A., 2018. Low tree-growth elasticity of forest biomass indicated by an individual-based model. Forests 9, 21. https://doi.org/10.3390/f9010021

Hember, R.A., Kurz, W.A., Girardin, M.P., 2019. Tree Ring Reconstructions of Stemwood Biomass Indicate Increases in the Growth Rate of Black Spruce Trees Across Boreal Forests of Canada. Journal of Geophysical Research: Biogeosciences 124, 2460–2480. https://doi.org/10.1029/2018JG004573

Kurz, W.A., Dymond, C.C., White, T.M., Stinson, G., Shaw, C.H., Rampley, G.J., Smyth, C., Simpson, B.N., Neilson, E.T., Trofymow, J.A., Metsaranta, J., Apps, M.J., 2009. CBM-CFS3: A model of carbon-dynamics in forestry and land-use change implementing IPCC standards. Ecological Modelling 220, 480–504. https://doi.org/10.1016/j.ecolmodel.2008.10.018 


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

