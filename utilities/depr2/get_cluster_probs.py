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
    
def get_model_probabilities():
    
    all_probs_dict={} 
   
    score_index=3 # total_score
        
    scores_file=open('all_model_scores.out','r')
    for ln in scores_file.readlines():
        fields = ln.strip().split()
        score=float(fields[score_index])
        
        prob=math.exp(-1.0*score)
        
        runid=int(fields[0])
        replicaid=int(fields[1])
        frameid=int(fields[2])
        
        all_probs_dict[(runid,replicaid,frameid)]=prob
    
    scores_file.close()
    
    # get probs for gsms only as a list, indexed by model number 
    gsm_probs_list=[]
    gsm_ids_file=open('model_ids_scores.txt','r')
    for ln in gsm_ids_file.readlines():
        (runid,replicaid,frameid)=(int(r) for r in ln.strip().split()[1:4])
        gsm_probs_list.append(all_probs_dict[(runid,replicaid,frameid)])
   
    gsm_ids_file.close()
    
    return(gsm_probs_list) 

##################
#### main
#################

expt = sys.argv[1]

bio_system = sys.argv[2]

required_representation =sys.argv[3]

expts_dir = os.path.join(os.path.expanduser('~'),"optrep/expts")

config_file = os.path.join(os.path.expanduser('~'),"optrep/input/"+bio_system,bio_system+".config."+expt)

configs_dict = stats_helper.parse_config_file(config_file)

os.chdir(os.path.join(expts_dir,'expt'+expt,bio_system))

if required_representation == "max":
    optimalRepresentationDir = 'r'+str(max([int(cg.lstrip('r')) for cg in glob.glob('r*')]))

else:
    optimalRepresentationDir = 'r'+required_representation
    
os.chdir(os.path.join(optimalRepresentationDir,'good_scoring_models'))

# Step 1. Get the significant clusters and members of the clusters
if required_representation=="30":
    required_representation="max"

sampling_precision=get_sampling_precision_from_stats_file(os.path.join(expts_dir,'stats',bio_system+".expt"+expt+".res"+required_representation+".stats.txt"))

stats_helper.get_clusters_at_given_precision('model_sample_ids.txt',configs_dict["PROTEINS_TO_OPTIMIZE_LIST"],sampling_precision)

cluster_members = stats_helper.read_cluster_members()
 
# Step 2. Get the model scores from file
gsm_probs = get_model_probabilities()

# Step 3a. Get probability of each cluster 
prob_models_clusters=[] # list of list of probabilities of models in each cluster

retained_cluster_probs=[]

for clus in cluster_members:
    prob_models_clusters.append([])
    
    for mem in clus: 
        # depends on how we do the averaging! 
        prob_models_clusters[-1].append(gsm_probs[mem])    
        retained_cluster_probs.append(gsm_probs[mem])

# Type 1 sum_prob_clus_by_sum_retained_clus_prob
# Step 3b. Normalize probabilities of each cluster
normalized_cluster_probs=[]

for clus_probs in prob_models_clusters:

    clus_numerator = sum(clus_probs)
    clus_denominator = sum(retained_cluster_probs)

    normalized_cluster_probs.append(clus_numerator/clus_denominator)

# Step 4. Write to file
output_file=open('cluster_probs_1.txt','w')
for elem in normalized_cluster_probs:
    print >>output_file,elem
output_file.close()

#Type 2 pop_clus_by_sum_pop_retained_clus 
# Step 3b. Normalize probabilities of each cluster
normalized_cluster_populations=[]

for clus_probs in prob_models_clusters:

    clus_numerator = float(len(clus_probs))
    clus_denominator = float(len(retained_cluster_probs))

    normalized_cluster_populations.append(clus_numerator/clus_denominator)

# Step 4. Write to file
output_file=open('cluster_probs_2.txt','w')
for elem in normalized_cluster_populations:
    print >>output_file,elem
output_file.close()
            
#Type 3 avg_prob_by_sum_avg_retained_clus_prob 
# Step 3b. Normalize probabilities of each cluster
avg_cluster_probs=[]

for clus_probs in prob_models_clusters:
    avg_cluster_probs.append(sum(clus_probs)/float(len(clus_probs)))
    
# divide each by sum of average
normalized_cluster_probs=[i/sum(avg_cluster_probs) for i in avg_cluster_probs]

# Step 4. Write to file
output_file=open('cluster_probs_3.txt','w')
for elem in normalized_cluster_probs:
    print >>output_file,elem
output_file.close()
