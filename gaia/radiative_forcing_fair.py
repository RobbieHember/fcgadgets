'''
ESTIMATE RADIATIVE FORCING FROM GHG EMISSIONS
'''

#%% Import modules
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from fair import FAIR
from fair.interface import fill,initialise
from fair.io import read_properties

#%% Project managment
path=r'C:\Data\FAIR\Pulses'

#%%
f=FAIR(ch4_method='thornhill2021')
f.define_time(1750,2100,1) # Define time period

#%% Define scenarios
sL=['ssp370']#,'ssp245','ssp370','ssp434','ssp460','ssp534-over','ssp585']
f.define_scenarios(sL)

#%% Import data
dfConfig=pd.read_csv(path + '//4xCO2_cummins_ebm3.csv')
df_volcanic=pd.read_csv(path + '//volcanic_ERF_monthly_175001-201912.csv',index_col='year')

#%% Define configs
models=dfConfig['model'].unique()
configs=[]
for imodel, model in enumerate(models):
	for run in dfConfig.loc[dfConfig['model']==model,'run']:
		configs.append(f"{model}_{run}")
f.define_configs(configs)

#%% Define species
#species,properties=read_properties(path + '//species_configs_properties.csv')
species,properties=read_properties()
f.define_species(species,properties)

#%% Create arrays
f.allocate()

#%% Populate inputs
f.fill_species_configs()

#%% Populate emissions
f.fill_from_rcmip()
f.fill_from_csv(emissions_file=path + '//emissions.csv')

# Overwrite volcanic emissions
df_volcanic[1750:].head()
volcanic_forcing=np.zeros(351)
volcanic_forcing[:271]=df_volcanic[1749:].groupby(np.ceil(df_volcanic[1749:].index) // 1).mean().squeeze().values
fill(f.forcing,volcanic_forcing[:,None,None],specie="Volcanic") # sometimes need to expand the array

# The concentration in the first timestep will be set to the baseline concentration, which are the IPCC AR6 1750 values.
initialise(f.concentration,f.species_configs['baseline_concentration'])
initialise(f.forcing,0)
initialise(f.temperature,0)
initialise(f.cumulative_emissions,0)
initialise(f.airborne_emissions,0)

#%% Populate climate configurations
models=dfConfig['model'].unique()
seed=1355763
for config in configs:
	model, run = config.split('_')
	condition = (dfConfig['model']==model) & (dfConfig['run']==run)
	fill(f.climate_configs['ocean_heat_capacity'], dfConfig.loc[condition, 'C1':'C3'].values.squeeze(), config=config)
	fill(f.climate_configs['ocean_heat_transfer'], dfConfig.loc[condition, 'kappa1':'kappa3'].values.squeeze(), config=config)
	fill(f.climate_configs['deep_ocean_efficacy'], dfConfig.loc[condition, 'epsilon'].values[0], config=config)
	fill(f.climate_configs['gamma_autocorrelation'], dfConfig.loc[condition, 'gamma'].values[0], config=config)
	fill(f.climate_configs['sigma_eta'], dfConfig.loc[condition, 'sigma_eta'].values[0], config=config)
	fill(f.climate_configs['sigma_xi'], dfConfig.loc[condition, 'sigma_xi'].values[0], config=config)
	fill(f.climate_configs['stochastic_run'], True, config=config)
	fill(f.climate_configs['use_seed'], True, config=config)
	fill(f.climate_configs['seed'], seed, config=config)
	seed = seed + 399

#%% Run
f.run()

#%% Export results
d={}
d['Time']=f.timebounds
plt.close('all')
for scn in sL:
	d[scn]={}
	d[scn]['Emissions']={}
	for sp in species:
		d[scn]['Emissions'][sp]=np.mean(f.emissions.loc[dict(scenario=scn,specie=sp)].to_numpy(),axis=1)
	d[scn]['Concentrations']={}
	for sp in species:
		d[scn]['Concentrations'][sp]=np.mean(f.concentration.loc[dict(scenario=scn,specie=sp)].to_numpy(),axis=1)
	#d[scn]['Concentrations']['CO2']=np.mean(f.concentration.loc[dict(scenario=scn,specie='CO2')].to_numpy(),axis=1)
	#d[scn]['Concentrations']['CH4']=np.mean(f.concentration.loc[dict(scenario=scn,specie='CH4')].to_numpy(),axis=1)
	#d[scn]['Concentrations']['N2O']=np.mean(f.concentration.loc[dict(scenario=scn,specie='N2O')].to_numpy(),axis=1)
	d[scn]['Temp']=np.mean(f.temperature.loc[dict(scenario=scn,layer=0)].to_numpy(),axis=1)
	d[scn]['RF Total']=np.mean(f.forcing_sum.loc[dict(scenario=scn)].to_numpy(),axis=1)
	#plt.plot(d['Time'],d[scn]['RF Total'])

	#plt.plot(d['Time'][1:],d[scn]['Emissions']['CO2 FFI'])
	#plt.plot(d['Time'][1:],d[scn]['Emissions']['CO2 AFOLU'])
	#plt.plot(d['Time'][1:],d[scn]['Emissions']['CO'])

	#plt.plot(d['Time'],d[scn]['Concentrations']['CO2'],'k-')

#%%
#d['ssp126']['Emissions']['CO2']-d['ssp119']['Emissions']['CO2']

#plt.close('all');plt.plot(d['Time'],d['ssp126']['Concentrations']['CO2']-d['ssp119']['Concentrations']['CO2'],'b-')

#plt.close('all');plt.plot(d['Time'],d['ssp126']['RF Total']-d['ssp119']['RF Total'],'b-')

plt.close('all');plt.plot(d['Time'],d['ssp126']['Temp']-d['ssp119']['Temp'],'b-')

