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
import IMP.optrep
import IMP.optrep.BeadMapBuilder
    
def get_gsm_loglike():
    
    all_loglike_dict={} 
   
    score_index=4 # log likelihood
        
    scores_file=open('all_model_scores.out','r')
    for ln in scores_file.readlines():
        fields = ln.strip().split()
        score=float(fields[score_index])
        
        runid=int(fields[0])
        replicaid=int(fields[1])
        frameid=int(fields[2])
        
        all_loglike_dict[(runid,replicaid,frameid)]=score
    
    scores_file.close()
    
    # get probs for gsms only as a list, indexed by model number 
    gsm_probs_list=[]
    gsm_ids_file=open('model_ids_scores.txt','r')
    for ln in gsm_ids_file.readlines():
        (runid,replicaid,frameid)=(int(r) for r in ln.strip().split()[1:4])
        gsm_probs_list.append(all_loglike_dict[(runid,replicaid,frameid)])
   
    gsm_ids_file.close()
    
    return(gsm_probs_list) 

def get_number_of_params(bead_precisions_file):
    
    bpf=open(bead_precisions_file,'r')
    
    num_params=0;
    for ln in bpf.readlines():
        num_params+=1
        
    bpf.close()
    return(num_params)
    
    
##################
#### main
##################

bio_system = sys.argv[1]

expt = sys.argv[2]

rep = sys.argv[3]

# 1: sum of probs of cluster by sum of probs over all retained clusters that are mapped on to an r1 cluster.
# 2: population of cluster by sum of populations of all retained clusters that are mapped on to an r1 cluster.
# 3: average of cluster by sum of averages of all retained clusters. NOT USED any more. 
# 4: clustering both model and target ensembles to get the same number of clusters?? E.g. 5 clusters?  NOT USED any more. 
# 5: comparing the prob of a target ensemble GSM with the prob of a GSM in the model ensemble. Target GSMs are normalized like usual. Model ensemble GSMs are normalized by sum over all GSMs that are mapped. 

models_dir = os.path.join(os.path.expanduser('~'),"optrep/expts/expt"+expt,bio_system,'r'+rep,"good_scoring_models")

os.chdir(models_dir)

num_params = get_number_of_params("../bead_precisions_"+rep+".txt")

gsm_log_like = get_gsm_loglike()

sample_size=len(gsm_log_like)
avg_log_like=sum(gsm_log_like)/float(len(gsm_log_like))

aic = 2.0*num_params - 2.0*avg_log_like

bic = num_params*math.log(sample_size) - 2.0*avg_log_like

if rep=="30":
    rep="non-uni_30"
    
print bio_system,rep,sample_size,"aic:",aic,"bic:",bic

    