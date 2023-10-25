import argparse
import os

from idmtools.analysis.analyze_manager import AnalyzeManager
from idmtools.core import ItemType

from simulation.analyzer_emodpy.analyzer_collection import (
    EventReporterAnalyzer,
    MonthlyPfPRAnalyzer,
    InsetChartAnalyzer,
    EventReporterSummaryAnalyzer,
    NodeDemographicsAnalyzer,
    VectorStatsAnalyzer
)

sweep_variables = ['Run_Number', 'Sample_ID']
sim_years=60



def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-name", dest="expt_name", type=str, required=True)
    parser.add_argument("-id", dest="expt_id", type=str, required=True)

    return parser.parse_args()


def analyze_experiment(platform, expt_id, wdir):
    if not os.path.exists(wdir):
        os.makedirs(wdir)

    analyzers = []
    # custom analyzers
    analyzers.append(EventReporterSummaryAnalyzer(sweep_variables=sweep_variables,
                                                 working_dir=wdir, 
                                                 time_cutoff=(sim_years-10)*365,
                                                 #event_list=["Received_Treatment"],
                                                 event_list=["Received_Treatment", "Bednet_Using","Bednet_Got_New_One","Bednet_Discarded", "Received_SMC"],
                                                 output_filename="event_counts"))
    analyzers.append(NodeDemographicsAnalyzer(sweep_variables=sweep_variables,
                                              working_dir=wdir,
                                              output_filename="age_population",
                                              time_cutoff=(sim_years-10)*365))
    analyzers.append(InsetChartAnalyzer(sweep_variables=sweep_variables,
                                       working_dir=wdir,
                                       start_day=(sim_years-10)*365,
                                       channels=["PCR Parasite Prevalence", "Adult Vectors", "Daily Bites per Human", "Air Temperature", "Rainfall", "Daily EIR", "Infectious Vectors"]))

    analyzers.append(VectorStatsAnalyzer(sweep_variables=sweep_variables,
                                         working_dir=wdir,
                                         start_time=(sim_years-10)*365, 
                                         end_time=sim_years*365)
                                         )  
                                                                       
    # Don't change these - used for fitting #
    analyzers.append(MonthlyPfPRAnalyzer(sweep_variables=sweep_variables,
                                        working_dir=wdir,
                                        start_year=2015,
                                        end_year=2020)) 
    analyzers.append(EventReporterAnalyzer(sweep_variables=sweep_variables,
                                          working_dir=wdir, 
                                          time_cutoff=(sim_years-5)*365,
                                          event_list=["Received_Treatment"],
                                          output_filename="events"))
    
    manager = AnalyzeManager(platform=platform,
                             configuration={},
                             ids=[(expt_id, ItemType.EXPERIMENT)],
                             analyzers=analyzers,
                             partial_analyze_ok=True,
                             max_workers=16)
    manager.analyze()


if __name__ == "__main__":
    
    from idmtools.core.platform_factory import Platform
    import manifest

    args = parse_args()
    platform = Platform('SLURM_LOCAL', job_directory=manifest.job_directory)
    outdir = args.expt_name
    analyze_experiment(platform, 
                       args.expt_id,
                       os.path.join(manifest.output_dir, outdir))
