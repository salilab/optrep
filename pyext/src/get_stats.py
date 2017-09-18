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

class Violations(object):

    def __init__(self, xlinkfile,threshold,xlink_keyword):

        self.violation_threshold  = threshold 
        self.xlink_keyword=xlink_keyword
        self.violation_counts = {}  # dictionary with a key per restraint

        xf = open(xlinkfile)
        for ln in xf.readlines(): 
            if ln.startswith('PROTEIN1') or ln.startswith('res1'):
                continue
            (res1,pos1,res2,pos2) =ln.strip().split(',')[0:4]
        
            self.violation_counts[(res1,pos1,res2,pos2)]=0 
        xf.close()

    def get_xlink_distances(self, lndict,xlinkType="XLDSS"):
        # get the minimum distances for each xlink'ed pair in a single model

        minimum_distances = {}
        dist_over_threshold = {}
        mindist_copies = {}

        for k in self.violation_counts:
            minimum_distances[k] = 10000.0 # some big number
            dist_over_threshold[k] = 0.0 
            mindist_copies[k] = ""
        
        for rst in lndict:
            items = rst.split('|')
            if items[0]!=self.xlink_keyword or items[1]!=xlinkType:
                continue
            (res1,pos1,res2,pos2) = items[3:7]
                     
            if '.' in res1:
                resonly1 = res1.split('.')[0]  # separate the residue name from copy number
            else:
                resonly1 = res1
            if '.' in res2:
                resonly2 = res2.split('.')[0]   # separate the residue name from copy number
            else:
                resonly2=res2
        
            modeldistance = float(lndict[rst])
       
            if modeldistance<minimum_distances[(resonly1,pos1,resonly2,pos2)]:
                minimum_distances[(resonly1,pos1,resonly2,pos2)]=modeldistance
                mindist_copies[(resonly1,pos1,resonly2,pos2)]=res1+"|"+pos1+"|"+res2+"|"+pos2 #write in the appropriate way
            
        for k in minimum_distances:
            if minimum_distances[k]<=self.violation_threshold:
                dist_over_threshold[k] = 0.0                
        
            else:
                dist_over_threshold[k] = minimum_distances[k]-self.violation_threshold
	   	    
        return minimum_distances,mindist_copies,dist_over_threshold
   
def get_run_replica_model_numbers(ids_file):
    
    run_replica_model_ids=[]
    idf=open(ids_file,'r')
    for ln in idf.readlines():
        fields = ln.strip().split()
        run_replica_model_ids.append((fields[1],fields[2],fields[3]))
 
    idf.close()
    
    return run_replica_model_ids
    
def get_fit_to_xlinks(xlink_file,xlink_threshold,xlink_keyword):
    
    gsm_ids_list = get_run_replica_model_numbers(os.path.join('good_scoring_models','model_ids_scores.txt'))
    
    Analysis = Violations(xlink_file,xlink_threshold,xlink_keyword) 
      
    for count,(runid,replicaid,modelid) in enumerate(gsm_ids_list):
        
        stat_file_line=subprocess.check_output(["awk","NR=="+int(modelid+1),"{print}",os.path.join("run."+runid,"output","stat."+replicaid+".out")])
        
        (xlink_distances_model,minimum_distance_copy_model,xlink_distances_over_threshold_model)=Analysis.get_xlink_distances(eval(stat_file_line.strip('\n'))) # return a dictionary keyed by crosslinks, and values equal to values for a model
        
        if count==0:
            distance_between_xlink_beads={}
            distance_over_threshold_between_xlink_beads={}
            
            for k in xlink_distances_model:
                distance_between_xlink_beads[k]=[]
                distance_over_threshold_between_xlink_beads[k]=[]
            
        for k in xlink_distances_model:
                distance_between_xlink_beads[k].append(xlink_distances_model[k])
                distance_over_threshold_between_xlink_beads[k].append(xlink_distances_over_threshold_model[k])
     
    # average over each xlink. 
    avg_distances=[]
    avg_distances_over_threshold=[]
    
    for k in distance_between_xlink_beads:
        avg_distances.append(sum(distances_between_xlink_beads[k])/float(len(distances_between_xlink_beads[k])))
        avg_distances_over_threshold.append(sum(distances_between_xlink_beads_over_threshold[k])/float(len(distances_between_xlink_beads_over_threshold[k])))
        
        print avg_distances[-1],avg_distances_over_threshold[-1]
    # delete till here
     
    # get average over all xlinks and all models
    num_xlinks_in_all_models=0.0
    sum_xlink_distances_across_all_models=0.0
    sum_xlink_distances_over_thresholds_across_all_models=0.0
    
    for k in distances_between_xlink_beads:
        sum_xlink_distances_across_all_models += sum(distances_between_xlink_beads[k])
        sum_xlink_distances_over_thresholds_across_all_models += sum(distances_over_threshold_between_xlink_beads[k])
        
        num_xlinks_in_all_models+=len(distances_between_xlink_beads[k])
        
    avg_distances_across_all_xlinks=sum_xlink_distances_across_all_models/float(num_xlinks_in_all_models)
    
    avg_distances_over_threshold_across_all_xlinks=sum_xlink_distances_over_thresholds_across_all_models/float(num_xlinks_in_all_models)
    
    return(avg_distances_across_all_xlinks,avg_distances_over_threshold_across_all_xlinks)

##################
#### main
#################

expt = sys.argv[1]

bio_system = sys.argv[2]

required_representation =sys.argv[3]

config_file = os.path.join(os.path.expanduser('~'),"optrep/input/"+bio_system,bio_system+".config."+expt)

configs_dict = parse_config_file(config_file)

out_file=open(os.path.join(os.path.expanduser('~'),"optrep/expts/stats",bio_system+"."+expt+".stats.txt"),"w")

os.chdir(os.path.join(os.path.expanduser('~'),'optrep/expts/expt'+expt,bio_system))

if required_representation == "max":
    optimalRepresentationDir = 'r'+str(max([int(cg.lstrip('r')) for cg in glob.glob('r*')]))

else:
    optimalRepresentationDir = 'r'+required_representation
    
os.chdir(optimalRepresentationDir)

# 1. Get sampling efficiency
machine_run_on = subprocess.check_output(["awk","$1=="+expt+" && $2==\""+bio_system+"\" {print $3}",os.path.join(os.path.expanduser('~'),"optrep/expts/info.expts")]).strip()

avg_sampling_time = stats_helper.get_sampling_time(machine_run_on,expt)

print >>out_file,"Average sampling time in seconds on ",machine_run_on," : %.2f" %(float(avg_sampling_time))

# 2. Get representation precision 
all_beads_stats,data_beads_stats=stats_helper.get_representation_resolution_and_precision(configs_dict["PROTEINS_TO_OPTIMIZE_LIST"],configs_dict["DOMAINS_TO_OPTIMIZE_LIST"],xlink_file=os.path.join(configs_dict["INPUT_DIR"],configs_dict["XLINKS_FILE"]))

print >>out_file,"All bead average resolution, diameter, weighted diameter :%.2f %.2f %.2f" %(all_beads_stats["avg_num_residues"],all_beads_stats["avg_dia"],all_beads_stats["avg_weighted_dia"])                                                                                                                               
                                                                                                                                 
print >>out_file,"Data bead average resolution, diameter, weighted diameter :%.2f %.2f %.2f" %(data_beads_stats["avg_num_residues"],data_beads_stats["avg_dia"],data_beads_stats["avg_weighted_dia"])

# 3. Get fit to data
for i,data_type in enumerate(configs_dict["GOOD_SCORING_MODEL_CRITERIA_LIST"]): # for each criterion/data set/data type
    
    if data_type=="Crosslinks":
        stat_file_keyword=configs_dict["GOOD_SCORING_MODEL_KEYWORD_LIST"][i]
        
        violation_threshold=float(configs_dict["GOOD_SCORING_MODEL_MEMBER_UPPER_THRESHOLD"][i])
        
        (avg_distance_between_xlink_beads,avg_distance_over_threshold_between_xlink_beads) = get_fit_to_xlinks(os.path.join(configs_dict["INPUT_DIR"],configs_dict["XLINKS_FILE"]),violation_threshold,stat_file_keyword)
        # TODO assuming only 1 xlink type that is related to 1 xlink file. There could be multiple in principle, though!
        
        print >>out_file,"Average distance between xlink'ed beads:",avg_distance_between_xlink_beads
        print >>out_file,"Average distance over threshold between xlink'ed beads:",avg_distance_over_threshold_between_xlink_beads

out_file.close()
