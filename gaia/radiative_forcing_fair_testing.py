'''
EXPLORING FAIR
'''

#%% Import modules
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from fair import FAIR
from fair.interface import fill,initialise
from fair.io import read_properties
from fair.earth_params import seconds_per_year

'''
PULSE BY YEAR
'''

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

#%%





'''
CMIP6 RCPs
'''

#%% Project managment
path=r'C:\Data\FAIR\Anthropocene'

#%%
f=FAIR(ch4_method='thornhill2021')
f.define_time(1750,2100,1) # Define time period

#%% Define scenarios
sL=['ssp119','ssp126','ssp245','ssp370','ssp434','ssp460','ssp534-over','ssp585']
f.define_scenarios(sL)

#%% Define configs
df=pd.read_csv(path + '//4xCO2_cummins_ebm3.csv')
models=df['model'].unique()
configs=[]
for imodel, model in enumerate(models):
	for run in df.loc[df['model']==model,'run']:
		configs.append(f"{model}_{run}")
f.define_configs(configs)

#%% Define species
#species,properties=read_properties(path + '//species_configs_properties.csv')
#f.define_species(species,properties)
species,properties=read_properties()
#species = list(properties.keys())
f.define_species(species,properties)

#%% Create arrays
f.allocate()

#%% Populate inputs
f.fill_species_configs()
#fill(f.species_configs['unperturbed_lifetime'], 10.8537568, specie='CH4')
#fill(f.species_configs['baseline_emissions'], 19.01978312, specie='CH4')
#fill(f.species_configs['baseline_emissions'], 0.08602230754, specie='N2O')

#%% Populate emissions

# Import prepared CMPIP6 emissions from RCMIP
f.fill_from_rcmip()

# Overwrite volcanic emissions
df_volcanic=pd.read_csv(path + '//volcanic_ERF_monthly_175001-201912.csv',index_col='year')
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
df=pd.read_csv(path + '//4xCO2_cummins_ebm3.csv')
models=df['model'].unique()
seed=1355763
for config in configs:
	model, run = config.split('_')
	condition = (df['model']==model) & (df['run']==run)
	fill(f.climate_configs['ocean_heat_capacity'], df.loc[condition, 'C1':'C3'].values.squeeze(), config=config)
	fill(f.climate_configs['ocean_heat_transfer'], df.loc[condition, 'kappa1':'kappa3'].values.squeeze(), config=config)
	fill(f.climate_configs['deep_ocean_efficacy'], df.loc[condition, 'epsilon'].values[0], config=config)
	fill(f.climate_configs['gamma_autocorrelation'], df.loc[condition, 'gamma'].values[0], config=config)
	fill(f.climate_configs['sigma_eta'], df.loc[condition, 'sigma_eta'].values[0], config=config)
	fill(f.climate_configs['sigma_xi'], df.loc[condition, 'sigma_xi'].values[0], config=config)
	fill(f.climate_configs['stochastic_run'], True, config=config)
	fill(f.climate_configs['use_seed'], True, config=config)
	fill(f.climate_configs['seed'], seed, config=config)
	seed = seed + 399

#%% Run
f.run()

#%% Export results
d={}
d['Time']=f.timebounds
d['Scn']={}
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

	plt.plot(d['Time'],d[scn]['Concentrations']['CO2'],label=scn)

plt.xlabel('year')
#plt.ylabel('Temperature anomaly (K)')
plt.legend()

#%% Export emissions
for scn in sL:
	df=pd.DataFrame.from_dict(d[scn]['Emissions']).to_excel(r'C:\Data\FAIR\Pulses\cmip6_emissions_' + scn + '.xlsx')

#%%

'''
SIMPLT RUN EXAMPLE
'''

#%% Project managment
path=r'C:\Data\FAIR\Simple Run'

#%%
f=FAIR()
f.define_time(2000,2050,1) # Define time period

#%% Define scenarios
scnL=['abrupt','ramp']
f.define_scenarios(scnL)

#%% Define three scenarios
f.define_configs(["high","central","low"])
#f.configs

#%% Define species
species,properties=read_properties(path + '//species_configs_properties.csv')
f.define_species(species,properties)

#%% Create arrays
f.allocate()

#%% Populate inputs
f.fill_from_csv(
    emissions_file=path + '//emissions.csv',
    concentration_file=path + '//concentration.csv',
    forcing_file=path + '//forcing.csv')

#f.emissions.loc[(dict(specie="CO2 FFI", scenario="abrupt"))] = 38

#%% Define initial conditions
initialise(f.concentration, 278.3,specie='CO2')
initialise(f.forcing, 0)
initialise(f.temperature, 0)
initialise(f.cumulative_emissions, 0)
initialise(f.airborne_emissions, 0)
initialise(f.ocean_heat_content_change, 0)

#%% Populate climate parameters
f.fill_species_configs(path + '//species_configs_properties.csv')
f.override_defaults(path + '//configs_ensemble.csv')

#%% Run
f.run()

#%% Export results
d={}
d['Time']=f.timebounds
d['Scn']={}
for scn in scnL:
	d[scn]={}
	d[scn]['CO2']=f.concentration.loc[dict(scenario=scn,specie='CO2')].to_numpy()
	d[scn]['CH4']=f.concentration.loc[dict(scenario=scn,specie='CH4')].to_numpy()

#%%
plt.plot(d['Time'],d['CO2'],'ko')

#%%


#%%
plt.plot(f.timebounds, f.temperature.loc[dict(scenario='abrupt', layer=0)], label=f.configs)
plt.xlabel('year')
plt.ylabel('Temperature anomaly (K)')
plt.legend()

#%%
plt.plot(f.timebounds,f.concentration.loc[dict(scenario='ramp', specie='CO2')], label=f.configs)
plt.title('Ramp scenario: CO2')
plt.xlabel('year')
plt.ylabel('CO2 (ppm)')
plt.legend()

#%%
plt.plot(f.timebounds,f.forcing.loc[dict(scenario='ramp', specie='Aerosol-cloud interactions')], label=f.configs)
plt.title('Ramp scenario: forcing')
plt.xlabel('year')
plt.ylabel('ERF from aerosol-cloud interactions (W m$^{-2}$)')
plt.legend()

plt.plot(f.timebounds, f.forcing_sum.loc[dict(scenario='ramp')], label=f.configs)
plt.title('Ramp scenario: forcing')
plt.xlabel('year')
plt.ylabel('Total ERF (W m$^{-2}$)')
plt.legend()

plt.plot(f.timebounds, f.temperature.loc[dict(scenario='abrupt', layer=0)], label=f.configs)
plt.title('Abrupt scenario: temperature')
plt.xlabel('year')
plt.ylabel('Temperature anomaly (K)')
plt.legend()

plt.plot(f.timebounds, f.forcing_sum.loc[dict(scenario='abrupt')], label=f.configs)
plt.title('Abrupt scenario: forcing')
plt.xlabel('year')
plt.ylabel('Total ERF (W m$^{-2}$)')
plt.legend()

plt.plot(f.timebounds, f.concentration.loc[dict(scenario='abrupt', specie='CO2')], label=f.configs)
plt.title('Abrupt scenario: CO2')
plt.xlabel('year')
plt.ylabel('CO2 (ppm)')
plt.legend()
f.species_configs['g0'].loc[dict(specie='CO2')]
f.forcing[-1, :, 1, :]



#%%
#%%












#%% Define methods, time horizon
f=FAIR(ch4_method="Thornhill2021")
f.define_time(1750,2300,1)  # start, end, step

#%% Define scenarios
scenarios = [
    "high-extension", 
    "high-overshoot",
    "medium-overshoot", 
    "medium-extension", 
    "low", 
    "verylow",
    "verylow-overshoot"]
f.define_scenarios(scenarios)

#%% Define parameters
fair_params_1_4_1_file = 'data/calibrated_constrained_ensemble/calibrated_constrained_parameters_calibration1.4.1.csv'

df_configs = pd.read_csv(fair_params_1_4_1_file, index_col=0)
configs = df_configs.index  # this is used as a label for the "config" axis
f.define_configs(configs)

#%% Define species
fair_species_configs_1_4_1_file = 'data/calibrated_constrained_ensemble/species_configs_properties_calibration1.4.1.csv'
species, properties = read_properties(filename=fair_species_configs_1_4_1_file)
f.define_species(species, properties)

#%%
f.allocate()


#%%
f.fill_from_csv(
    emissions_file='data/calibrated_constrained_ensemble/extensions_1750-2500.csv',
    forcing_file='data/calibrated_constrained_ensemble/volcanic_solar.csv',
)

#%%
f.emissions

#%%


#%%


#%%