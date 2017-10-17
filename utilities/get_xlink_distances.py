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

os.chdir(os.path.join(expts_dir,'expt'+expt,bio_system))

if required_representation == "max":
    optimalRepresentationDir = 'r'+str(max([int(cg.lstrip('r')) for cg in glob.glob('r*')]))

else:
    optimalRepresentationDir = 'r'+required_representation
    
os.chdir(optimalRepresentationDir)
# 3. Get fit to data
data_type =configs_dict["GOOD_SCORING_MODEL_CRITERIA_LIST"]  # for each criterion/data set/data type

if data_type=="Crosslinks":
    avg_distance_between_xlink_beads,std_err_distance_between_xlink_beads = stats_helper.get_fit_to_xlinks(float(configs_dict["GOOD_SCORING_MODEL_MEMBER_UPPER_THRESHOLDS_LIST"]),configs_dict["GOOD_SCORING_MODEL_KEYWORD_LIST"])
