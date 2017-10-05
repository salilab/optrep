import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os,sys,math,string
import numpy as np
import stats_helper

def get_precisions_per_residue(protein,domain,bead_map_file,bead_precisions_file):
    
    bead_starts_and_ends = []
    precisions_per_residue=[]
    
    bmf=open(bead_map_file,'r')
    for ln in bmf.readlines():
        fields=ln.strip().split()
        
        if fields[0]==protein and fields[1]==domain:
            bead_starts_and_ends.append((int(fields[2]),int(fields[3])))
        
    bmf.close()    
         
    bpf=open(bead_precisions_file,'r')
    
    for ln in bpf.readlines():
        fields=ln.strip().split()
        
        if fields[0]==protein and fields[1]==domain:
            bead_index = int(fields[2])
            
            for residue in range(bead_starts_and_ends[bead_index][0],bead_starts_and_ends[bead_index][1]+1):
                precisions_per_residue.append(residue,fields[3])
      
    bpf.close()
    
    return precisions_per_residue
    
    
    
##################
#### main
#################

expt = sys.argv[1]

bio_system = sys.argv[2]

config_file = os.path.join(os.path.expanduser('~'),"optrep/input/"+bio_system,bio_system+".config."+expt)

configs_dict = stats_helper.parse_config_file(config_file)

xlink_file=configs_dict["XLINKS_FILE"]

proteins_list = configs_dict["PROTEINS_TO_OPTIMIZE_LIST"]
domains_list = configs_dict["DOMAINS_TO_OPTIMIZE_LIST"]
resolutions_list=configs_dict["RESOLUTIONS_LIST"]

precisions_by_residue={}

data_dir = os.path.join(os.path.expanduser('~'),"optrep/expts/expt"+expt,bio_system)

topology_file=os.path.join(config_params["INPUT_DIR"],config_params["TOPOLOGY_FILE"])

for protein,domain in zip(proteins_list,domains_list):
   
    for resolution in resolutions_list:
        precisions_by_residue[protein+"_"+domain+"_"+resolution] = get_precisions_per_residue(protein,domain,os.path.join(data_dir,"r"+resolution,"bead_map_"+resolution+".txt"),os.path.join(data_dir,"r"+resolution,"bead_precisions_"+resolution+".txt"))
        
        print precisions_by_residue
        
        
        
        
