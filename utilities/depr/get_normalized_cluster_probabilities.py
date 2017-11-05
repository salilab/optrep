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

def get_sampling_precision_from_stats_file(stats_file):
    
    sf = open(stats_file,'r')
    for ln in sf.readlines():
        if ln.startswith("Sampling precision"):
            return(float(ln.strip().split()[2]))
   
    sf.close()
    
    return 100000.00 # some arbitrary big number
    
def get_model_probabilities(score_type,models_to_average):
    
    all_probs_dict={} # store all probs in 2 formats: as a dict and as a list 
    
    all_probs_list = []
    
    # first parse the file of all model scores
    if score_type=="total_score":
        score_index=3
    else:
        score_index=4 # data posterior only!
        
    scores_file=open('all_model_scores.out','r')
    for ln in scores_file.readlines():
        fields = ln.strip().split()
        score=float(fields[score_index])
        
        prob=math.exp(-1.0*score)
        all_probs_list.append(prob)
        
        runid=int(fields[0])
        replicaid=int(fields[1])
        frameid=int(fields[2])
        
        all_probs_dict[(runid,replicaid,frameid)]=prob
        
        #print runid,replicaid,frameid,score,prob
    
    scores_file.close()
    
    # get probs for gsms only as a list, indexed by model number 
    gsm_probs_list=[]
    gsm_ids_file=open('model_ids_scores.txt','r')
    for ln in gsm_ids_file.readlines():
        (runid,replicaid,frameid)=(int(r) for r in ln.strip().split()[1:4])
        gsm_probs_list.append(all_probs_dict[(runid,replicaid,frameid)])
        #print runid,replicaid,frameid,gsm_probs_list[-1]
    gsm_ids_file.close()
    
    return(gsm_probs_list,all_probs_list) 
    
def get_aggregate(lst,averaging_type):
    if averaging_type=="direct_average":
        return sum(lst)
    else:
        return sum(lst)/float(len(lst))

##################
#### main
#################

expt = sys.argv[1]

bio_system = sys.argv[2]

required_representation =sys.argv[3]

expts_dir = sys.argv[4]

functn = sys.argv[5] # "cluster" or "probability"

score_type = sys.argv[6] # total_score or data_posterior 

averaging_type = sys.argv[7] # direct_average or average_of_average

models_to_average = sys.argv[8] # avg over retained_clusters or all_gsms or all_models 

config_file = os.path.join(os.path.expanduser('~'),"optrep/input/"+bio_system,bio_system+".config."+expt)

configs_dict = stats_helper.parse_config_file(config_file)

os.chdir(os.path.join(expts_dir,'expt'+expt,bio_system))

if required_representation == "max":
    optimalRepresentationDir = 'r'+str(max([int(cg.lstrip('r')) for cg in glob.glob('r*')]))

else:
    optimalRepresentationDir = 'r'+required_representation
    
os.chdir(os.path.join(optimalRepresentationDir,'good_scoring_models'))

if functn == "cluster":
    # Step 1. Get the significant clusters and members of the clusters
    sampling_precision=get_sampling_precision_from_stats_file(os.path.join(expts_dir,'stats',bio_system+".expt"+expt+".res"+required_representation+".stats.txt"))

    stats_helper.get_clusters_at_given_precision('model_sample_ids.txt',configs_dict["PROTEINS_TO_OPTIMIZE_LIST"],sampling_precision)

    exit(1)

else: 
     
    cluster_members = stats_helper.read_cluster_members()
     
    # Step 2. Get the model scores from file
    gsm_probs,all_model_probs = get_model_probabilities(score_type,models_to_average)

    # Step 3a. Get probability of each cluster 
    prob_models_clusters=[] # list of list of probabilities of models in each cluster

    for clus in cluster_members:
        prob_models_clusters.append([])
        
        for mem in clus: 
            # depends on how we do the averaging! 
            prob_models_clusters[-1].append(gsm_probs[mem])    
        
    # Step 3b. Normalize probabilities of each cluster
    normalized_cluster_probs=[]
     
    for clus_probs in prob_models_clusters:
        
         clus_numerator = get_aggregate(clus_probs,averaging_type)
            
         if models_to_average == "all_models":
             clus_denominator = get_aggregate(all_model_probs,averaging_type)
         elif models_to_average == "all_gsms":
             clus_denominator = get_aggregate(gsm_probs,averaging_type)
         elif models_to_average == "retained_clusters":
             retained_cluster_probs=[]
             for sublist in prob_models_clusters:
                 for itm in sublist:
                    retained_cluster_probs.append(itm)
 
             clus_denominator = get_aggregate(retained_cluster_probs,averaging_type)

         normalized_cluster_probs.append(clus_numerator/clus_denominator)
     
    # Step 4. Write to file
    output_file=open('cluster_probs_'+score_type+'_'+averaging_type+'_'+models_to_average+'.txt','w')
    for elem in normalized_cluster_probs:
         print >>output_file,elem
    output_file.close()
     
    exit(1)
