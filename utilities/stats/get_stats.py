import IMP
import IMP.rmf
import RMF
import IMP.atom
import os,sys,math
import numpy as np
import glob
import datetime
import subprocess
import stats_helper

##################
#### main
#################

expt = sys.argv[1]

bio_system = sys.argv[2]

required_representation =sys.argv[3]

expts_dir = sys.argv[4]

config_file = os.path.join(os.path.expanduser('~'),"optrep/input/"+bio_system,bio_system+".config."+expt)

configs_dict = stats_helper.parse_config_file(config_file)

out_file=open(os.path.join(os.path.expanduser('~'),"optrep/expts/stats",bio_system+".expt"+expt+".res"+required_representation+".stats.txt"),"w")

os.chdir(os.path.join(expts_dir,'expt'+expt,bio_system))

if required_representation == "max":
    optimalRepresentationDir = 'r'+str(max([int(cg.lstrip('r')) for cg in glob.glob('r*')]))

else:
    optimalRepresentationDir = 'r'+required_representation
    
os.chdir(optimalRepresentationDir)

# 1. Get sampling efficiency
machine_run_on = subprocess.check_output(["awk","$1=="+expt+" && $2==\""+bio_system+"\" {print $3}",os.path.join(os.path.expanduser('~'),"optrep/expts/info.expts")]).strip()

avg_sampling_time, std_err_sampling_time = stats_helper.get_sampling_time(machine_run_on,expt)

print >>out_file,"Average sampling time and std error in seconds on ",machine_run_on," : %.2f %.2f" %(float(avg_sampling_time),float(std_err_sampling_time))

# 2. Get representation precision 
all_beads_stats,data_beads_stats,non_data_beads_stats=stats_helper.get_representation_resolution_and_precision(configs_dict["PROTEINS_TO_OPTIMIZE_LIST"],configs_dict["DOMAINS_TO_OPTIMIZE_LIST"],xlink_file=os.path.join(configs_dict["INPUT_DIR"],configs_dict["XLINKS_FILE"]))

print >>out_file,"All bead average resolution, diameter, weighted diameter :%.2f %.2f %.2f" %(all_beads_stats["avg_num_residues"],all_beads_stats["avg_dia"],all_beads_stats["avg_weighted_dia"])                                                                                                                               
                                                                                                                                 
print >>out_file,"Data bead average resolution, diameter, weighted diameter :%.2f %.2f %.2f" %(data_beads_stats["avg_num_residues"],data_beads_stats["avg_dia"],data_beads_stats["avg_weighted_dia"])

print >>out_file,"Non-data bead average resolution, diameter, weighted diameter :%.2f %.2f %.2f" %(non_data_beads_stats["avg_num_residues"],non_data_beads_stats["avg_dia"],non_data_beads_stats["avg_weighted_dia"])

# 3. Get fit to data
data_type =configs_dict["GOOD_SCORING_MODEL_CRITERIA_LIST"]  # for each criterion/data set/data type

if data_type=="Crosslinks":
    avg_distance_between_xlink_beads,std_err_distance_between_xlink_beads = stats_helper.get_fit_to_xlinks(float(configs_dict["GOOD_SCORING_MODEL_MEMBER_UPPER_THRESHOLDS_LIST"]),configs_dict["GOOD_SCORING_MODEL_KEYWORD_LIST"])
    # TODO assuming only 1 xlink type that is related to 1 xlink file. There could be multiple in principle, though!
    
    print >>out_file,"Average distance and std error between xlink'ed beads:%.2f %.2f" %(float(avg_distance_between_xlink_beads),float(std_err_distance_between_xlink_beads))

os.chdir('good_scoring_models')

sampling_precision=stats_helper.get_sampling_precision('model_sample_ids.txt',configs_dict["PROTEINS_TO_OPTIMIZE_LIST"],float(configs_dict["GRID_SIZE"]))

print >>out_file,"Sampling precision %.2f" %(sampling_precision)

out_file.close()
