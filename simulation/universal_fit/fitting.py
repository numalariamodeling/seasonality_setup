"""
Run 30 years of simulations with varying temperature and habitat type.
"""

import os
from datetime import datetime as dt
from functools import partial

import emod_api.campaign as camp
import emod_api.config.default_from_schema_no_validation as dfs
import emod_api.demographics.Demographics as Demog
import emodpy_malaria.demographics.MalariaDemographics as Demographics
import emodpy_malaria.malaria_config as conf
import emodpy_malaria.malaria_config as malaria_config
from emod_api.demographics.DemographicsTemplates import CrudeRate
from emodpy.emod_task import EMODTask
from emodpy_malaria.interventions.treatment_seeking import add_treatment_seeking
import emodpy_malaria.interventions.usage_dependent_bednet as itn 
import emodpy_malaria.interventions.drug_campaign as dc    
from emodpy_malaria.reporters.builtin import (
    add_event_recorder,
    add_malaria_summary_report,
    add_report_node_demographics,
    add_report_vector_stats
)
from idmtools.builders import SimulationBuilder
from idmtools.core.platform_factory import Platform
from idmtools.entities import Suite
from idmtools.entities.experiment import Experiment
import pandas as pd

import manifest
from sweeping import CfgFn, ItvFn, sweep_functions


# batch_1: BK 60 years beginning in 1960 (TEST)
# rnd0 -id 0beba4e5-2b73-4662-8c7b-116a05f208dd
# rnd1 -id 

# Global variables
batch = "batch_1"
weatherset = 'tempshift_0'
rnd = 1

DS = "BK"
ds = DS.lower()
sim_years = 60
sim_start_year = 1960
num_seeds = 1

# Functions
def set_climate(config, shift):
    """
    Set climate to specific files, currently hardcoded.
    """
    config.parameters.Air_Temperature_Offset = shift
    config.parameters.Land_Temperature_Offset = shift

    return {"Temperature_Shift": shift}


def set_species_from_file(config, vector_file):
    """
    Set species and habitat level, assuming 1 species and 2 habitat type
    (1 input 1 constant)
    """
    # Multiple-Species Capability
    vdf = pd.read_csv(os.path.join(manifest.input_dir,DS,vector_file))                                              
    vdf = vdf[vdf['node_id']== vdf['node_id'].unique()[0]]                                                                 
    s = [species for species in vdf['species']]                                                                            
    conf.add_species(config, manifest, s) 
    for r in range(len(s)):          
        conf.set_species_param(config, 
                               species = vdf['species'][r],
                               parameter='Anthropophily',
                               value=vdf['anthropophily'][r],
                               overwrite=True)
        conf.set_species_param(config, 
                               species = vdf['species'][r],
                               parameter='Indoor_Feeding_Fraction',
                               value=vdf['indoor_feeding'][r],
                               overwrite=True)                                                                                     

    return

def set_param_fn(config):
    """
    This function is a callback that is passed to emod-api.config to set config
    parameters, including the malaria defaults.
    """
    config = conf.set_team_defaults(config, manifest)

    #set_species(config)
    set_species_from_file(config, vector_file="vectors.csv")

    config.parameters.Simulation_Duration = sim_years * 365 + 60
    config.parameters.Run_Number = 0
    config.parameters.x_Temporary_Larval_Habitat = 1
    config.parameters.Birth_Rate_Dependence = "FIXED_BIRTH_RATE"
    config.parameters.Age_Initialization_Distribution_Type = "DISTRIBUTION_COMPLEX"

    config.parameters.Air_Temperature_Filename = os.path.join(
        "climate", "air_temperature_daily.bin"
    )
    config.parameters.Land_Temperature_Filename = os.path.join(
        "climate", "air_temperature_daily.bin"
    )
    config.parameters.Rainfall_Filename = os.path.join(
        "climate", "rainfall_daily.bin"
    )
    config.parameters.Relative_Humidity_Filename = os.path.join(
        "climate", "relative_humidity_daily.bin"
    )
    config.parameters.Custom_Individual_Events = ["Received_Treatment","Received_SMC","Bednet_Using","Bednet_Discarded","Bednet_Got_New_One"]

    return config


def set_param(simulation, param, value):
    """
    Set specific parameter value
    Args:
        simulation: idmtools Simulation
        param: parameter
        value: new value
    Returns:
        dict
    """
    return simulation.task.set_parameter(param, value)


def build_camp():
    """
    This function builds a campaign input file for the DTK using emod_api.
    """
    camp.set_schema(manifest.schema_file)

    return camp


def build_demog():
    """
    This function builds a demographics input file for the DTK using emod_api.
    """

    new_nodes = [Demog.Node(lat=1.857481477,lon=6.444799744, pop=1000, name=DS, forced_id=1)]
    demog = Demographics.MalariaDemographics(
        nodes=new_nodes,
        idref=DS,
        init_prev=0.2,
        include_biting_heterogeneity=True,
    )
    demog.SetEquilibriumVitalDynamics(CrudeRate(38.0))
    demog.SetEquilibriumAgeDistFromBirthAndMortRates(CrudeRate(38.0), CrudeRate(38.0))
    demog.SetBirthRate(CrudeRate(38.0 * 1000))

    return demog


def add_outputs(task):
    """
    Requesting reports/outputs to the task.
    """
    for year in range(sim_years - 5, sim_years):
        start_day = 0 + 365 * year
        sim_year = sim_start_year + year
        add_malaria_summary_report(
            task,
            manifest,
            start_day=start_day,
            end_day=365 + year * 365,
            reporting_interval=30,
            age_bins=[0.25, 5, 115],
            max_number_reports=13,
            pretty_format=True,
            filename_suffix=f"Monthly{sim_year}",
        )

    # add_vector_habitat_report(task, manifest)
    add_event_recorder(
        task,
        event_list=["Received_Treatment","Bednet_Using","Bednet_Got_New_One","Bednet_Discarded", "Received_SMC"],
        start_day=(sim_years - 10) * 365,
        end_day=sim_years * 365,
        min_age_years=0,
        max_age_years=50,
    )

    add_report_node_demographics(task, manifest,
                                 age_bins = [0.25,5,50],
                                 stratify_by_gender = False)

    add_report_vector_stats(task, manifest,
                            species_list = ['arabiensis','funestus','gambiae'],
                            stratify_by_species = True,
                            include_death_state = False,
                            include_wolbachia = False,
                            include_gestation = False,
                            include_microsporidia = False)


    return

def set_habitat_from_file(config, hab_multipliers, vector_file, samp_id, hab_base=1e8, const_base=1e6):
    """
    Set mosquitoes so all species have all three habitat types, and set multipliers universally,
    Baseline: 1e8 for specific habitat, 1e6 for constant habitat.
    hab_multipliers always a list of three, corresponding to CONSTANT,
    TEMPORARY_RAINFALL and WATER_VEGETATION
    """
    vdf = pd.read_csv(os.path.join(manifest.input_dir,DS,vector_file))                                                                                                          
    s = [species for species in vdf['species']]   

    for r in range(len(s)):  
        habitat1 = dfs.schema_to_config_subnode(
            manifest.schema_file, ["idmTypes", "idmType:VectorHabitat"]
        )
        habitat1.parameters.Habitat_Type = "CONSTANT"
        malaria_config.set_species_param(
            config, vdf['species'][r], "Habitats", habitat1.parameters, overwrite=True
        )

        habitat2 = dfs.schema_to_config_subnode(
            manifest.schema_file, ["idmTypes", "idmType:VectorHabitat"]
        )
        habitat2.parameters.Habitat_Type = "TEMPORARY_RAINFALL"
        malaria_config.set_species_param(config, vdf['species'][r], "Habitats", habitat2.parameters)

        habitat3 = dfs.schema_to_config_subnode(
            manifest.schema_file, ["idmTypes", "idmType:VectorHabitat"]
        )
        habitat3.parameters.Habitat_Type = "WATER_VEGETATION"
        malaria_config.set_species_param(config, vdf['species'][r], "Habitats", habitat3.parameters)

        conf.set_max_larval_capacity(
            config, vdf['species'][r], "CONSTANT", const_base * (vdf['fraction'][r] * vdf['constant'][r]) * hab_multipliers[0]
        )
        conf.set_max_larval_capacity(
            config, vdf['species'][r], "TEMPORARY_RAINFALL", hab_base * (vdf['fraction'][r] * vdf['temp_rain'][r]) * hab_multipliers[1]
        )
        conf.set_max_larval_capacity(
            config, vdf['species'][r], "WATER_VEGETATION", hab_base * (vdf['fraction'][r] * vdf['water_veg'][r]) * hab_multipliers[2]
        )

    return {
        "CONST_Multiplier": hab_multipliers[0],
        "TEMPR_Multiplier": hab_multipliers[1],
        "WATEV_Multiplier": hab_multipliers[2],
        "Sample_ID": samp_id
    }
        

def add_cm_from_file(campaign, cm_file):
    ##### Case management ##################################################################################################                             
    if os.path.exists(os.path.join(manifest.input_dir,DS,"interventions",cm_file)):                                 #
        cm_df = pd.read_csv(os.path.join(manifest.input_dir,DS,"interventions",cm_file))                            # Read cm_file
        nodes = [n for n in cm_df['node_id'].unique()]                                                                     #
        for node in nodes:                                                                                                 # For each node...
            cm_df = cm_df[cm_df['node_id']==node]                                                                          #   
            for year in cm_df['year']:                                                                                     # For each year...
                  sub_df = cm_df[cm_df['year'] == year].reset_index()                                                      #
                  targets = [] 
                  #print(year)                                                                                             #
                  for r in range(len(sub_df)) :                                                                            # ... Build a set of targets from rows of cm_file ...
                      cm_coverage_by_age =  {'trigger': str(sub_df['trigger'][r]),                                         #
                                             'coverage': float(sub_df['coverage'][r]),                                     #
                                             'agemin': float(sub_df['age_min'][r]),                                        #
                                             'agemax': float(sub_df['age_max'][r]),                                        #                                        #
                                             'rate': float(sub_df['rate'][r])}                                             #
                      targets.append(cm_coverage_by_age)                                                                   #
                  #print(targets)
                  add_treatment_seeking(campaign, node_ids = [int(node)],                                                   # ... Add treatment seeking to simulation.
                                        start_day = int(sub_df['start_day'][0]),                                        #
                                        duration = int(sub_df['duration'][0]),                                          #
                                        drug=['Artemether','Lumefantrine'],                                             #   **Hard-Coded treatment with AL**
                                        targets=targets,                                                                #
                                        broadcast_event_name="Received_Treatment")                                      #
    return{"CM_file":cm_file}

def add_itn_from_file(campaign, itn_file, itn_age_file, itn_season_file):
    ##### Bednets ##################################################################################################                             
    if os.path.exists(os.path.join(manifest.input_dir,DS,"interventions",itn_file)):                                 #
        itn_df = pd.read_csv(os.path.join(manifest.input_dir,DS,"interventions", itn_file))                         # Read itn_file
        nodes = [n for n in itn_df['node_id'].unique()]
        for node in nodes:
            itn_df = itn_df[itn_df['node_id']==node]                                                          # Filter to burnin or pickup
            if len(itn_df) > 0:                                                                                            #
                itn_age = pd.read_csv(os.path.join(manifest.input_dir,DS,"interventions",itn_age_file))             # Read age dependence file
                itn_season = pd.read_csv(os.path.join(manifest.input_dir,DS,"interventions",itn_season_file))       # Read seasonal dependence file - **currently assumes same seasonal pattern for all years**
                itn_seasonal_usage = {"Times": list(itn_season['season_time']),                                            #
                                      "Values":list(itn_season['season_usage'])}                                           #
                for year in itn_df['year']:                                                                                # For each year ...
                    sub_df = itn_df[itn_df['year']==year].reset_index()                                                    #
                    itn_discard_config = {"Expiration_Period_Distribution": "WEIBULL_DISTRIBUTION",                        # ...Set discard distribution based on itn_file
                                          "Expiration_Period_Kappa": float(sub_df['discard_k'][0]),                        #
                                          "Expiration_Period_Lambda": float(sub_df['discard_l'][0])}                       #
                    itn_age_year = itn_age[itn_age['year']==year]                                                          # ...Set age-dependence
                    itn_age_bins = itn_age_year['age']                                                                     #
                    itn_age_usage = itn_age_year['age_usage']                                                              #
                    itn.add_scheduled_usage_dependent_bednet(campaign, node_ids = [int(node)],                                 # ...Add scheduled usage-dependent bednets with specified:
                                                             intervention_name = "UsageDependentBednet",                   # 
                                                             start_day = int(sub_df['start_day'][0]),                      # - distribution timing 
                                                             demographic_coverage = float(sub_df['coverage'][0]),          # - distribution coverage
                                                     			   killing_initial_effect = float(sub_df['kill_effect'][0]),     # - insecticide properties
                                                     			   killing_decay_time_constant = int(sub_df['kill_decay'][0]),   # 
                                                     			   blocking_initial_effect = float(sub_df['block_effect'][0]),   # 
                                                     			   blocking_decay_time_constant=int(sub_df['block_decay'][0]),   # 
                                                             age_dependence = {"Times": list(itn_age_bins),                # - Age effect on usage 
                                                                               "Values": list(itn_age_usage)},             #
                                                             seasonal_dependence = itn_seasonal_usage,                     # - Seasonal effect on usage
                                                             discard_config = itn_discard_config)                          # - Discard probability
    return{"ITN_file":itn_file, "ITN_age_file": itn_age_file, "ITN_season_file":itn_season_file}

    
def add_smc_from_file(campaign, smc_file, leak_age_max, leak_multiplier):
    if os.path.exists(os.path.join(manifest.input_dir,DS,"interventions",smc_file)):                                #
        smc_df = pd.read_csv(os.path.join(manifest.input_dir,DS,"interventions",smc_file), encoding='latin')        # Read smc_file
        nodes = [n for n in smc_df['node_id'].unique()]                                                                    #
        for node in nodes:                                                                                                 #
            smc_df = smc_df[smc_df['node_id']==node]                                                                       #                                                       # Filter to burnin or pickup
            if len(smc_df) > 0:                                                                                            #
                for r in range(len(smc_df)):                                                                               # For each row (round) in smc_file...
                    dc.add_drug_campaign(campaign, campaign_type="MDA", drug_code="SPA", node_ids =[int(node)],                     # ... Add SMC with drug SPA
                                         start_days=[int(smc_df['start_day'][r])],                                         # ... on this day
                                         repetitions=1,                                                                    # ... once
                                         coverage=float(smc_df['coverage'][r]),                                            # ... at this coverage
                                         target_group={'agemin': 0.25, 'agemax': 5},                                       # ... for children age 0.25 to 5
                                         receiving_drugs_event_name="Received_SMC")                                        #
                    dc.add_drug_campaign(campaign, campaign_type="MDA", drug_code="SPA", node_ids = [int(node)],                    # and ... Add SMC "leak" for children ages 5 to 6
                                         start_days=[int(smc_df['start_day'][r])],                                         #
                                         repetitions = 1,                                                                  #
                                         coverage=float(smc_df['coverage'][r]) * leak_multiplier,                          # Older children have half the coverage **(hard-coded)**     
                                         target_group={'agemin': 5, 'agemax': leak_age_max},                               # Ages 5-6
                                         receiving_drugs_event_name="Received_SMC")                                        #
    return{"SMC_file":smc_file, "SMC_leak_age_max":leak_age_max, "SMC_leak_multiplier":leak_multiplier}


def config_task(platform):
    # create EMODTask
    print("Creating EMODTask (from files)...")
    task = EMODTask.from_default2(
        config_path=None,
        eradication_path=manifest.eradication_path,
        campaign_builder=build_camp,
        schema_path=manifest.schema_file,
        param_custom_cb=set_param_fn,
        ep4_custom_cb=None,
        demog_builder=build_demog,
    )

    # set the singularity image to be used when running this experiment
    task.set_sif(manifest.SIF_PATH, platform)

    # add weather directory as an asset
    indir = os.path.join(f"{DS}", f"{DS}_10yr")
    task.config.parameters.Birth_Rate_Dependence = "FIXED_BIRTH_RATE"
    task.common_assets.add_directory(
        os.path.join(manifest.input_dir, indir, weatherset), relative_path="climate"
    )

    # add output
    add_outputs(task)

    return task


def config_sweep_builder(df):
    builder = SimulationBuilder()

    sweeps = [
        [
            CfgFn(
                set_habitat_from_file,
                hab_multipliers=row[
                    ["CONST_Multiplier", "TEMPR_Multiplier", "WATEV_Multiplier"]
                ],
                samp_id=row['samp_id'],
                vector_file="vectors.csv"
            ),
            CfgFn(set_climate, shift=row["Temperature_Shift"]),
            ItvFn(add_cm_from_file, cm_file ="case_management.csv"),
            ItvFn(add_itn_from_file, itn_file = "ITN.csv", itn_age_file="ITN_age.csv", itn_season_file="ITN_season.csv"),
            ItvFn(add_smc_from_file, smc_file="SMC.csv", leak_age_max=6, leak_multiplier=0.5),
            partial(set_param, param="Run_Number", value=s)
        ]
        for rid, row in df.iterrows()
        for s in range(100, 100 + num_seeds)
    ]

    builder.add_sweep_definition(sweep_functions, sweeps)
    print(builder.count)

    return [builder]


def general_sim(selected_platform):
    """
    This function is designed to be a parameterized version of the sequence of things we
    do every time we run an emod experiment.
    """

    # Set platform and associated values, such as the maximum number of jobs to run at
    # one time
    platform = Platform(
        selected_platform,
        job_directory=manifest.job_directory,
        partition="short",
        time="2:00:00",
        account="p30781",
        modules=["singularity"],
        max_running_jobs=1500,
    )

    # create experiment from builder
    homepath = os.path.expanduser("~")
    user = homepath.split("/")[2]
    suite_name = f"{user}_{DS}_univ_season"
    suite = Suite(name=suite_name)
    suite.uid = suite_name
    samp_df = pd.read_csv(f"simulation_output/{DS}/{batch}/rnd{rnd}.csv")
    #samp_df = samp_df.loc[0:4,:]

    task = config_task(platform)
    builder = config_sweep_builder(df=samp_df)

    now = dt.now()
    now_str = now.strftime("%Y_%m_%d_%H_%M_%S")
    exp_name = f"fit_rnd{rnd}"
    exp_id = f"{exp_name}_{now_str}"
    experiment = Experiment.from_builder(builder, task, name=exp_name)
    experiment.uid = exp_id
    suite.add_experiment(experiment)

    # The last step is to call run() on the ExperimentManager to run the simulations.
    experiment.run(wait_until_done=True, platform=platform)

    # Check result
    if not experiment.succeeded:
        print(f"Experiment {experiment.uid} failed.\n")
        exit()

    print(f"Experiment {experiment.uid} succeeded.")


if __name__ == "__main__":
    # If you don't have Eradication, uncomment out the following to download Eradication
    # import pathlib
    # import emod_malaria.bootstrap as dtk

    # dtk.setup(pathlib.Path(manifest.eradication_path).parent)
    selected_platform = "SLURM_LOCAL"
    general_sim(selected_platform)
