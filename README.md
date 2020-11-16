<a id="devex-badge" rel="Exploration" href="https://github.com/BCDevExchange/assets/blob/master/README.md"><img alt="Being designed and built, but in the lab. May change, disappear, or be buggy." style="border-width:0" src="https://assets.bcdevexchange.org/images/badges/exploration.svg" title="Being designed and built, but in the lab. May change, disappear, or be buggy." /></a>

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

## Module list
Biophysical processes in cbrunner are represented by a variety of Python methods within the fcgadgets.cbrunner.cbrun_annproc module. The cbrunner module can also call modules specifically developed to simulate natural disturbances, including fcgadgets.taz and forest management activities, including fcgadgets.activities.

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

