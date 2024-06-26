# fcgadgets
## INTRODUCTION
The Forest Carbon Gadgets (<b>fcgadgets</b>) repository supports estimation, accounting, and reporting of greenhouse gas (GHG) emissions in British Columbia�s forest sector. The repository is a toolbox that offers flexible processing. It exports the variables required to meet international reporting standards for the forest sector, or complete life cycle assessments.
<br>
<br>
The repository was developed to: 
* Achieve transparency and reproduceability as commitments to open government and open 
science standards
* Share knowledge, methods and limitations
* Automate and streamline workflow
* Promote a diverse ecosystem of existing and new modelling approaches
* Support complex policy decisions in land resource management

The <b>fcgadgets</b> repository was written in the Python programming language, benefitting from integrated libraries for simulation modelling, geographical information systems, data analytics, and application deployment (Downey, 2017). The <b>fcgadgets</b> repository was designed for a community that conducts forest carbon modelling on a full-time time basis. Users must be fluent in the Python language. Trying to apply <b>fcgadgets</b> without assistance is not advised. That said, it is possible for novice users to set up small projects that demonstrate dynamics for a single site. 

## PLUG-AND-PLAY MODULARITY
The repository allows for comprehensive representation of processes and new science by connecting a constellation of supporting modules.
![image info](./images/fcgadgets_constellation.png)

## CBRUNNER
<b>cbrunner</b> is a computer simulation model that estimates the greenhouse gas (GHG) balance of the forest sector, including forest ecosystems and wood products. The annual net flux of GHGs between the forest sector and the atmosphere is estimated by simulating several biophysical processes each year, including the biomass dynamics of trees, the decay and physical transformation of dead organic matter, the impact of natural disturbances, harvest removals, silvicultural treatments, and 
nutrient applications. 
![image info](./images/fcgadgets_annual_processes.png)

The model achieves this with a set of plug-and-play functions found in **cbrun_annproc.py**:
### Tree Biomass Dynamics (from Growth and Yield models): 
* Simulates tree biomass dynamics on an annual basis based on inputs of net biomass growth from the [TASS/TIPSY growth and yield software application](https://www2.gov.bc.ca/gov/content/industry/forestry/managing-our-forest-resources/forest-inventory/growth-and-yield-modelling).
* Default settings assume inputs generated with BatchTIPSY.exe, but this can be overridden to input tables generated with TASS
* Total stemwood growth is frequently zero for as much as 25 years during early stand development. This leads to underestimation of early biomass production when using allometric relationships between stemwood and other biomass pools. To avoid this, initial inputs of stemwood growth for the first 30 years of stand development are replaced with exponential increase in total stemwood biomass from 0 to the prediction at age 30. The exponential coefficient is solved such that net growth over the 30-year period will match that originally predicted by the GY model.
### Tree Biomass Dynamics (from Sawtooth):
* Simulates biomass dynamics of individual trees (Hember et al., 2019; Hember and Kurz, 2018)
* Distance-independent representation of resource competition
* Driven by equations of annual aboveground biomass growth, annual probability of recruitment, and annual probability of mortality
* Equations are fitted against species/region samples
### Dead Wood, Litter and Soil Dynamics
* Simulates cycling of organic carbon through:
	* Dead wood (snags and coarse woody debris);
	* Litter (organic soil horizon); 
	* Soil (mineral soil horizon);
	* Felled & piled materials
* Based on methods described by Kurz et al. (2009) and Shaw et al. (2014)
* Includes additional representation of piles
### Disturbance and Management Events: 
* This method imposes changes caused by natural disturbances and management events
* All events are defined by an event ID, decimal year, mortality factor, growth factor, and the ID of the growth curve that represents the new stand
* It is driven by the event chronology, which has two potential sources:
	* Prescribed by the user as input variables in the Disturbance and Management Event Chronology (DMEC)
	* Optional on-the-fly simulation of natural disturbances or management activities (based on functions of age or merchantable volume at the beginning of the year)

### Product Dynamics
* Representation of the annual GHG fluxes that arise from fibre that is removed from forest ecosystems
* Fate of removed fibre
* Product types
* Scenarios describing change in the fate of removed fibre and product types

### Geological Dynamics
This method represents annual GHG fluxes associated with:
* Forest sector operations (e.g., use of fossil fuels during hauling)
* Substitution of fossil fuels and cement for wood products 

### Model Structure
The <b>cbrunner</b> model has a hierarchical structure of forest stands, batches, scenarios, and ensembles:

N<sub>Simulation</sub> = N<sub>Stands</sub> � N<sub>Batches</sub> � N<sub>Scenarios</sub> � N<sub>Ensembles</sub>

Forest stands are the primary modelling unit in GHG estimation methods, and define an area of homogeneous conditions at the time a project is established (i.e., treatment area). Each stand is described by an inventory record, disturbance and management event chronology (DMEC), and age response functions of forest growth (if using a GY model). 

Projects with N<sub>Stands</sub> > 1,500 are segmented internally into batches that are run in sequence in order to work within the memory limits of individual work machines. Batch size (e.g., 1,500) is adjustable, but the batch size that optimizes simulation runtime, tends to be ~1,500 stands per unique combination of scenario and ensemble. 

### Scenarios Comparisons
Projects that explore climate change impacts or mitigation activities invariably consider multiple hypothetical scenarios for each forest stand. The hierarchical structure and post-processing scripts are
therefore built around running and comparing multiple scenarios.

### Uncertainty and Ensemble Forecasting
The <b>cbrunner</b> model adopts a probabilistic framework to accommodate processes with both deterministic and random components, as well as uncertainty analysis. Multiple ensembles
occur when project configuration specifies a stochastic component to simulations. This generally only occurs if users incorporate simulations of the annual 
probability of tree mortality or annual probability of tree recruitment. 

### Working with Growth & Yield Models
The <b>cbrunner</b> model can be driven with output from [TASS/TIPSY growth and yield (GY) modelling applications](https://www2.gov.bc.ca/gov/content/industry/forestry/managing-our-forest-resources/forest-inventory/growth-and-yield-modelling) using automated functions in **cbrun_utilities.py**.
* Prepare input parameters that are required to run BatchTIPSY.exe in a spreadsheet with **Write_BatchTIPSY_Input_Spreadsheet**
* Convert the spreadsheet of input parameters to the format expected by BatchTIPSY.exe using **Write_BatchTIPSY_Input_File**
* Convert the output from BatchTIPSY.exe to pickle files that will be read by cbrunner using **PostProcessBatchTIPSY**
* Convert the output from TASS for use in cbrunner with **GetTASSCurves**
* Import GY curves into a work session with **Import_BatchTIPSY_Output**

### Model Output Statistics (MOS)
Once simulations are complete, use a series of functions in <b>cbrun_util.py</b> to summarize model output statistics (MOS).
* Import simulation output variables for a given scenario, ensemble, and batch using **LoadSingleOutputFile**
* Import simulation output variables for a given scenario using **LoadScenarioResults**
* Calculate the mean and variance of ensemble simulations using **ModelOutputStats**

## MACGYVER
The <b>macgyver</b> toolbox contains custom scripts that compile information sources and prepare projects that use <b>cbrunner</b>. If pre-processing steps are similar among a wide range of project types, the goal is to store the scripts here for shared useage. 
* Pre-processing script template to prepare <b>cbrunner</b> inputs for a:
	* Sample of points
	* Sample of polygons
	* Tile or multi-tile project
* Methods for processing spatial information from: 
	* Vegetation Resource Inventory (VRI)
	* Reporting Silviculture Updates and Land Status Tracking System (RESULTS)
	* Wildfire perimiter and burn severity databases
	* Aerial overview (forest insects and disease) survey
	* Strategic land and resource plans
	* ClimateNA base-period mean climate
	* Growth and yield models

### util_inventory.py
The general workflow of <b>cbrunner</b> projects rely on the use of look-up tables (LUTs) for each variable in the inventory layers within Results.gdb, VRI.gdb, Disturbance.gdb, and LandUse.gdb.

### util_general.py
This module contains general utilities for workflow in Python.

### util_gis.py
This module contains utilities for performing spatial analysis in Python.

## TAZ
Forest sector GHG balance simulations depend on realistic variation of natural disturbances over space and time. While inventory records provide much of the information needed 
to represent natural disturbances over the modern era, additional simulations are needed to represent disturbances over the pre-inventory and future periods. The <b>taz</b> subpackage was developed to improve representation of disturbances in carbon models. It consists of statistics and scenarios of disturbance that were developed using a combination of observed constraints and probabilistic models. Despite high prediction uncertainty, using the pre-defined scenarios ensures that  representation of natural disturbances is grounded by available observations and science-informed 
scenarios, consistent across project studies, and supported by documentation.
### aspatial_stat_models.py 
* Statistical models of breakup as a function of stand age
* Statistical models of harvest as a function of standing merchantable volume 
* Annual area of occurrence (AAO) models of wildfire and beetles 
### onset_spread_models.py
* Spatially explicit simulations of events based on annual probability of onset and spread

## HARDHAT
The <b>hardhat</b> repository contains resources for representing effects of forest management on forest sector GHG balance.
### Nutrient Applications
The <b>nutrient_application</b> module contains functions that update annual nutrient status called by the <b>cbrunner</b> model and a function that schedules hypothetical nutrient applications (during the future period of simulation).
* Representation of GHG balance responses to aerial applications of Urea
* Schedule aerial nutrient applications with specified stand selection criteria
### economics.py
* Calculate cashflow from implementation of forest management events

## BC1HA
* The <b>bc1ha</b> repository supports raster processing on a 1 hectare regular grid of British Columbia.<br>
* The <b>bc1ha</b> raster framework is controlled with scripts in the <b>bc1ha</b> repository. The ordered workflow is shown in bc1ha_command.py, which calls functions in bc1ha_util.py.
* On an annual basis, the required vector data are downloaded from the BC Data Catalogue and stored in local geodatabases. The layers that are downloaded from BC Data Catalogue are listed in: cbrunner/Parameters/Table_BCFCS_DataSources.xlsx.
* The variables that are rasterized directly from layers on BC Data Catalogue are listed in: cbrunner/Parameters/Table_BCFCS_BC1haRasterVariableList.xlsx.
* Upon download, the vector data are rasterized. Categorical variables must be converted into numerical data. These numerical variables record an �ID� for each category and are stored as 8-, 16- or 32-bit integers depending on the required precision. This is achieved using look-up-tables (LUTs). The LUTs for categorical variables sourced directly from geodatabases are generated automatically during the annual update process. The LUT for any derived categorical variables (e.g., land cover compilation 1) is manually created in a spreadsheet stored in the <b>cbrunner</b> Parameters folder. 
* The LUTs are imported into the meta dictionary at the onset of each project and can be accessed through the meta dictionary, for example:
meta[�LUT�][� BEC_BIOGEOCLIMATIC_POLY�][�ZONE�][�CWH�]
returns ID = 6
* The ID associated with a specific category can be found using the lut_n2s function:
cbu.lut_n2s(6,meta[�LUT�][� BEC_BIOGEOCLIMATIC_POLY�][�ZONE�])
returning category �CWH�

## PROJECT WORKFLOW
There are two ways to apply <b>cbrunner</b>. In small projects, the land surface attributes and disturbance and management event chronology are specified in the project configuration spreadsheet by the user. This approach is suitable for projects with less than 20 stands or less than 20 scenarios. In bigger georeferenced projects, where modelling is conducted over thousands of hectares, the model is run from a python script. Big projects can subsample from landscapes to reduce computation time.
![image info](./images/fcgadgets_project_types.png)

A filing system needs to accommodate three information sources:
![image info](./images/fcgadgets_filingsystem.png)

## RUNNING SMALL PROJECTS
* If users only want to run a few stands, or a few scenarios or a small number of combinations, projects can be controlled entirely by the ProjectConfig.xlsx spreadsheet. To do so, �Scenario Source� in the Project Parameter tab must be set to �Spreadsheet� (as opposed to �Script�). 
* In this configuration, users must manually specify inventory variables in the �Scenario Parameters� tab and in the BatchTIPSY Parameters spreadsheet. 
* Users must manually specify the disturbance and management event chronologies in the Scenario Parameters tab.

## RUNNING BIG (GEOSPATIAL) PROJECTS
* <b>fcgadgets</b> is designed to run geospatial projects within a raster framework. All existing projects are connected to a code repository (<b>bc1ha</b>) and accompanying database of 1-hectare raster geotiffs. Projects that wish to run at a higher resolution can copy the bc1ha and emulate its functionality at a different desired resolution.
* To run geospatial projects, �Scenario Source� in the Project Parameter tab must be set to �Script� (as opposed to �Spreadsheet�).
* The spatial boundaries of a project are defined by a region of interest (ROI). The ROI is defined by parameters in the Project Parameters tab of the project configuration spreadsheet. Example:
       * �Name Project�: User specified
       * �ROI Source�: Options include �TSA�, �BGC Zone�, �Regional District�, or Custom name
       * �ROI Elements�: This allows one to specify multiple TSAs, or multiple BGC Zones
* The ROI is defined by binary mask spanning the extent of a regular grid. Simulation will be confined to any irregular shapes defined by the binary mask. 
* While the native spatial resolution of the bc1ha raster framework is 1 hectare, users can specify the �regular grid sampling frequency (RGSF)� in the Project Parameters tab. For example, setting RGSF = 50 will set up a project that runs simulations at every 50th hectare grid cell across the ROI. The project automatically records the Area expansion factor (AEF), which is used in post-processing to scale carbon reservoirs and fluxes to the native 1-hectare resolution.
* Every time a user specifies a new RGSF (e.g., 50 ha), the user must run a script from the bc1ha repository that creates �sparse� samples by subsetting the native 1-hectare geotiffs and storing them as .pkl files in �BC1ha/Sparse� folder for subsequent use. This step takes time, but once it is done once, all subsequent projects working at that RGSF will run faster because they only need to open the smaller .pkl files instead of the native 1-hectare grids. (Remember to keep sparse files updated to match annual updates of bc1ha.)
* In the Project Parameters tab, the �Batch Interval� value specifies the number of spatial units that will be run in each batch. The default is 1500, which was found to optimize run time. 
* Without special measures taken during pre-processing to identify reporting strata, the model output statistics (MOS) that are calculated during post processing of a project will only report values representing the full ROI. User can specify �Reporting Strata� during pre-processing, telling fcgadgets to calculate MOS for specific stratifications of the data. Currently, there is flexibility to stratify MOS in three different ways:
       * Project Type: Use this stratum to generate MOS by various attributes (e.g. obligation vs. non-obligation silviculture).
       * Year: Use this stratum to generate MOS by calendar year (e.g. results specific to events that occurred in 2023).
       * Spatial: Use this stratum to generate MOS by subregion of the ROI (e.g., timber supply areas)

## REFERENCES
Downey, A.B., 2017. Modeling and Simulation in Python � Green Tea Press, 2.3. ed. Green Tea Press, Needham, Massaschusetts.

Dymond, C.C., 2012. Forest carbon in North America: annual storage and emissions from British Columbia�s harvest, 1965-2065. Carbon Balance and Management 7, (24 July 2012)-(24 July 2012).

Hember, R.A., Kurz, W.A., 2018. Low tree-growth elasticity of forest biomass indicated by an individual-based model. Forests 9, 21. https://doi.org/10.3390/f9010021

Hember, R.A., Kurz, W.A., Girardin, M.P., 2019. Tree Ring Reconstructions of Stemwood Biomass Indicate Increases in the Growth Rate of Black Spruce Trees Across Boreal Forests of Canada. Journal of Geophysical Research: Biogeosciences 124, 2460�2480. https://doi.org/10.1029/2018JG004573

Kurz, W.A., Dymond, C.C., White, T.M., Stinson, G., Shaw, C.H., Rampley, G.J., Smyth, C., Simpson, B.N., Neilson, E.T., Trofymow, J.A., Metsaranta, J., Apps, M.J., 2009. CBM-CFS3: A model of carbon-dynamics in forestry and land-use change implementing IPCC standards. Ecological Modelling 220, 480�504. https://doi.org/10.1016/j.ecolmodel.2008.10.018 

Shaw, C.H., Hilger, A.B., Metsaranta, J., Kurz, W.A., Russo, G., Eichel, F., Stinson, G., Smyth, C., Filiatrault, M., 2014. Evaluation of simulated estimates of forest ecosystem carbon stocks using ground plot data from Canada�s National Forest Inventory. Ecological Modelling 272, 323�347.

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

