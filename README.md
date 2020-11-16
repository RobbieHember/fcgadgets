# fcgadgets
## Features
The fcgadgets package supports greenhouse gas (GHG) balance estimation, accounting, and reporting in British Columbia’s forest sector. It features a computer simulation model, **cbrunner**, that draws on established and custom methods to simulate net forest sector GHG balance. The package is written in the Python 3/Jupyter environment, benefiting from stable integrated libraries for simulation modelling, geographical information systems, data analytics, and application deployment (Downey, 2017). 

The package features a versatile and streamlined workflow through integration with BC databases and models: 
* Vegetation Resource Inventory (VRI);
* Reporting Silviculture Updates and Land Status Tracking System (RESULTS);
* Growth and yield models;

The aim of the package is to maintain the principles and standards applied in Canada’s National GHG Inventory.

As examples, the demos subpackage walks through demonstrations of Forest Carbon Initiative (FCI) project types:
* Aerial fertilization
* Underplanting fire-impacted forests
* Salvage logging beetle-impacted forests under varying fibre utilization intensities

## Users
The fcgadgets 

## Dependencies
Much of the in-house background science and analysis in support of the fcgadgets package is organized in the **fcexplore** package. The spatial reference system for many projects relies on information and processing from the BC1ha package.

## List of modules & methods
### BiomassFromTASSorTIPSY: 
* Simulates tree biomass dynamics on an annual basis based on inputs of net biomass growth from the TASS/TIPSY growth and yield software application (https://www2.gov.bc.ca/gov/content/industry/forestry/managing-our-forest-resources/forest-inventory/growth-and-yield-modelling)
* Default settings assume inputs generated with BatchTIPSY.exe, but this can be overridden to input tables generated with TASS
### BiomassFromSawtooth:
* Integration with an individual-tree, distance-independent model of biomass dynamics (Hember et al., 2019; Hember and Kurz, 2018).
### DOM_From_CBM08: 
* Representation of carbon cycling through snags, coarse woody debris, organic soil horizon, and mineral soil horizon following methods similar to the Carbon Budget Model of the Canadian Forest Sector (Kurz et al., 2009).
### DisturbanceAndManagement: 
* A custom Python method designed to represent wildfire, insects, disease, and management treatments.
Harvested wood products:
### HWP_From_Dymond12: 
* Representation of GHG balance for fibre that is removed from forest ecosystems. This module aims to capture the dynamics described by the BC Harvested Wood Products model version 1 (Dymond, 2012).
### nutrient_addition:
* Representation of GHG balance responses to aerial applications of Urea.
### taz: 
* Statistical models that represent natural disturbances.

## Workflow
There are four ways to apply cbrunner depending on the nature of the desired project. Small projects – with fewer than 1,500 combinations of locations or scenarios – can be run from a Jupyter Notebook. The work simply involves populating two Excel spreadsheets with the input variables and parameters. Bigger projects are scripted in Python and can adopt existing templates for projects that focus on running simulations at point locations, or across scattered polygons, or across continuous regular grids.
![image info](./images/fcgadgets_runoptions.png)

### Small projects (with Jupyter Notebooks)
When projects consist of fewer than 1,500 unique combinations of stand and/or scenario (e.g., see FCI project demos), cbrunner can be applied with assumptions set directly in the ProjectConfig.xlsx spreadsheet. It is recommended that the project script for small projects make use of Jupyter Notebooks (.ipynb files).
Project workflow:
1. Define project-level parameters in ProjectConfig.xlsx.
2. Define scenario-level parameters in ProjectConfig.xlsx.
3. Parameterize input parameters for BatchTIPSY.exe in the file, GrowthCurvesTIPSY_Parameters.xlsx. This can be done manually for small projects.
4. Convert input parameters to a format readable by BatchTIPSY.exe (automated by running fcgadgets.cbrunner.cbrun_utilities.py.BuildTIPSYInputs).
5. Run BatchTIPSY.exe.
6. Prepare inventory. (automated)
7. Prepare disturbance and management event history. (automated)
8. Prepare growth curves. (automated)
9. Run the simulation and save the outputs by calling the method RunProject. 
10. Import output variables to analysis session by calling LoadScenarioResults. 
11. Calculate GHG balance variables, including net sector greenhouse gas balance by calling the method CalculateGHGBalance.

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

