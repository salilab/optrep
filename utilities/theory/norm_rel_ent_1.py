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
        
        if prob==0.0:
            prob=sys.float_info.min
        
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

def get_ligand_coords(rmf_file,proteins_to_optimize,get_mass=False):
    # should return coords in the order that PMI created it, which should be in the order of the topology file.
    # which is why it would match the bead map! 
    
    coords=[]
    masses=[]
    
    m = IMP.Model()
    inf = RMF.open_rmf_file_read_only(rmf_file)
    h = IMP.rmf.create_hierarchies(inf, m)[0]
    IMP.rmf.load_frame(inf, 0)

    for state in h.get_children():
        for component in state.get_children():
            if stats_helper.included_protein(component.get_name(),proteins_to_optimize):
                for leaf in IMP.core.get_leaves(component):
                    p=IMP.core.XYZ(leaf.get_particle())
                    coords.append(p.get_coordinates())
                  
                    if get_mass:
                        masses.append(IMP.atom.Mass(leaf.get_particle()).get_mass())
                    
    if get_mass:
        return coords,masses  
    else:
        return coords
    
def read_cluster_coords(models_dir,proteins_to_optimize):
    
    ccf=open(os.path.join(models_dir,'cluster_centers.txt'),'r')
    cluster_centers = [ln.strip() for ln in ccf.readlines()]
    
    ccf.close()
    
    cluster_coords=[]
    
    for i,mdl in enumerate(cluster_centers):
                     
        if i==0:
            model_coords,model_masses = get_ligand_coords(os.path.join(models_dir,mdl+".rmf3"),proteins_to_optimize,get_mass=True)        
            cluster_masses=model_masses
        else:
            model_coords = get_ligand_coords(os.path.join(models_dir,mdl+".rmf3"),proteins_to_optimize)
            
        cluster_coords.append(model_coords)
    
    return cluster_coords,cluster_masses
    
def convert_fg_to_cg(coords,masses,fg_bmb,cg_bmb,proteins,domains):
    
    reduced_coords=[]
    
    fg_bead_index=0
        
    for prot,dom in zip(proteins,domains): # for each protein_domain that was CG'ed
    
        for cg_bead in cg_bmb.bead_maps[(prot,dom)]: # for each bead of the protein_domain
        
            centroid = IMP.algebra.Vector3D(0.0,0.0,0.0)
            mass_sum = 0.0
            
            # take center of mass of component beads
            for count in range(cg_bead[0],cg_bead[1]+1):
                #assuming FG beads have analogous order as CG beads
                #assuming CG'ing along backbone
                #assuming the comparison is with r1 
                
                centroid = centroid + masses[fg_bead_index]*coords[fg_bead_index] 
             
                mass_sum += masses[fg_bead_index]
                
                fg_bead_index+=1
                
            centroid=centroid/mass_sum 
            
            reduced_coords.append(centroid) 
    
    return(reduced_coords)
    
def coarse_grain_each_cluster_center(high_res_cluster_coords,high_res_masses,high_res_bead_map_file,cg_bead_map_file,proteins,domains):
    
    fg_bmb = IMP.optrep.BeadMapBuilder.BeadMapBuilder()
    fg_bmb.set_bead_map_from_beadmap_file(high_res_bead_map_file) 
    
    cg_bmb= IMP.optrep.BeadMapBuilder.BeadMapBuilder()
    cg_bmb.set_bead_map_from_beadmap_file(cg_bead_map_file)
    
    # CG according to the cg_bead_map
    reduced_cluster_coords=[]
    
    for cluster_coord in high_res_cluster_coords:
        # for each high res cluster
             
        reduced_cluster_coords.append(convert_fg_to_cg(cluster_coord,high_res_masses,fg_bmb,cg_bmb,proteins,domains))
        
    return reduced_cluster_coords
    
def get_cluster_mapping(fg_cluster_coords,cg_cluster_coords):
    
    fg_to_cg_mapping=[-1 for i in range(len(fg_cluster_coords))] # list with 1 elem per FG cluster. Each elem is the index of the CG cluster that it maps to.
    omega=[0 for i in range(len(cg_cluster_coords))] # multiplicity for each CG cluster (number of FG clusters that map to it)
    
    for ifg,fg_cluster_coord in enumerate(fg_cluster_coords):
        nearest_cg_cluster_index=1000 #arbitrary big number
        distance_to_nearest_cg_cluster=10000.0
        
        for icg,cg_cluster_coord in enumerate(cg_cluster_coords):
            
            curr_distance = IMP.algebra.get_rmsd(fg_cluster_coord,cg_cluster_coord) 
                        
            if curr_distance<distance_to_nearest_cg_cluster:
                distance_to_nearest_cg_cluster = curr_distance
                nearest_cg_cluster_index = icg
                
        fg_to_cg_mapping[ifg]= nearest_cg_cluster_index
        #print ifg,distance_to_nearest_cg_cluster
        omega[nearest_cg_cluster_index]+=1
      
            
    return fg_to_cg_mapping,omega

def read_cluster_probabilities(txt_file):
    
    clus_probs=[]
    tf=open(txt_file,'r')
    for ln in tf.readlines():
        clus_probs.append(float(ln.strip()))
    
    return clus_probs

##################
#### main
##################

bio_system = sys.argv[1]

expt_high_res = sys.argv[2]

rep_high_res = sys.argv[3]

expt_cg = sys.argv[4]

rep_cg = sys.argv[5]

prob_type="1_norm"

# 1: sum of probs of cluster by sum of probs over all retained clusters that are mapped on to an r1 cluster.
# 2: population of cluster by sum of populations of all retained clusters that are mapped on to an r1 cluster.
# 3: average of cluster by sum of averages of all retained clusters. NOT USED any more. 
# 4: clustering both model and target ensembles to get the same number of clusters?? E.g. 5 clusters?  
# 5: comparing the prob of a target ensemble GSM with the prob of a GSM in the model ensemble. Target GSMs are normalized like usual. Model ensemble GSMs are normalized by sum over all GSMs that are mapped. 

expts_dir = os.path.join(os.path.expanduser('~'),"optrep/expts")

cg_config_file = os.path.join(os.path.expanduser('~'),"optrep/input/"+bio_system,bio_system+".config."+expt_cg)

configs_dict = stats_helper.parse_config_file(cg_config_file)

proteins_to_optimize = configs_dict["PROTEINS_TO_OPTIMIZE_LIST"] 

domains_to_optimize = configs_dict["DOMAINS_TO_OPTIMIZE_LIST"] 

# Step 0. Do clustering, get the significant clusters and members of the clusters
# This is already done    
#sampling_precision=get_sampling_precision_from_stats_file(os.path.join(expts_dir,'stats',bio_system+".expt"+expt+".res"+name_representation+".stats.txt"))
#stats_helper.get_clusters_at_given_precision('model_sample_ids.txt',configs_dict["PROTEINS_TO_OPTIMIZE_LIST"],sampling_precision)

# Step 1. Load the cluster coordinates for both representations
high_res_models_dir = os.path.join(os.path.expanduser('~'),"optrep/expts/expt"+expt_high_res,bio_system,'r'+rep_high_res,"good_scoring_models")
cg_models_dir = os.path.join(os.path.expanduser('~'),"optrep/expts/expt"+expt_cg,bio_system,'r'+rep_cg,"good_scoring_models")

cg_bead_map_file = os.path.join(os.path.expanduser('~'),"optrep/expts/expt"+expt_cg,bio_system,'r'+rep_cg,'bead_map_'+rep_cg+'.txt')
high_res_bead_map_file = os.path.join(os.path.expanduser('~'),"optrep/expts/expt"+expt_high_res,bio_system,'r'+rep_high_res,'bead_map_'+rep_high_res+'.txt')

high_res_cluster_coords,high_res_masses = read_cluster_coords(high_res_models_dir,proteins_to_optimize)

cg_cluster_coords,cg_masses = read_cluster_coords(cg_models_dir,proteins_to_optimize)

# Step 2. Get the CG version of the high res coords and mapping between clusters of FG and CG coordinates
reduced_high_res_cluster_coords = coarse_grain_each_cluster_center(high_res_cluster_coords,high_res_masses,high_res_bead_map_file,cg_bead_map_file,proteins_to_optimize,domains_to_optimize)

fg_to_cg_mapping,omega_cg_clusters = get_cluster_mapping(reduced_high_res_cluster_coords,cg_cluster_coords)

prob_models_mapped_clusters={}

retained_mapped_cluster_probs={}

old_to_new_cg_cluster_index=[] # for each old cluster index what is the new cluster index? after removing clusters that are not mapped?

for representation,expt in zip([rep_high_res,rep_cg],[expt_high_res,expt_cg]):

    os.chdir(os.path.join(expts_dir,'expt'+expt,bio_system))

    os.chdir(os.path.join(os.path.expanduser('~'),"optrep/expts/expt"+expt,bio_system,'r'+representation,"good_scoring_models"))

    cluster_members = stats_helper.read_cluster_members()
 
    # Step 3. Get the model scores from file
    gsm_probs = get_model_probabilities()

    # Step 3a. Get probability of each cluster 
    prob_models_mapped_clusters[representation]=[] # list of list of probabilities of models in each cluster

    retained_mapped_cluster_probs[representation]=[]    

    for ic,clus in enumerate(cluster_members):
    
        if representation==rep_cg and omega_cg_clusters[ic]<1: # skip unmapped CG clusters 
            old_to_new_cg_cluster_index.append(-1) # mark this cluster as unindexed in the new set of mapped clusters
            continue
     
        prob_models_mapped_clusters[representation].append([])
    
        for mem in clus: 
            prob_models_mapped_clusters[representation][-1].append(gsm_probs[mem])    
            retained_mapped_cluster_probs[representation].append(gsm_probs[mem])
            
        if representation==rep_cg:
            old_to_new_cg_cluster_index.append(len(prob_models_mapped_clusters[representation])-1) # last index 

#print len(old_to_new_cg_cluster_index), len(prob_models_mapped_clusters[rep_cg])
            
# Step 3b. Normalize probabilities of each cluster
normalized_cluster_probs={}

for representation in prob_models_mapped_clusters:
    normalized_cluster_probs[representation]=[]
    
    for clus_probs in prob_models_mapped_clusters[representation]:

        clus_numerator = sum(clus_probs)
        clus_denominator = sum(retained_mapped_cluster_probs[representation])

        normalized_cluster_probs[representation].append(clus_numerator/clus_denominator)

# Step 4. Calculate relative entropy
relative_entropy=0.0
     
target_probs=normalized_cluster_probs[rep_high_res]

model_probs=normalized_cluster_probs[rep_cg]

print len(target_probs),len(model_probs)
    
for itarget in range(len(target_probs)):
    
    imodel_old = fg_to_cg_mapping[itarget]
    imodel_new = old_to_new_cg_cluster_index[imodel_old] # skipped some clusters so need new index
    
    if omega_cg_clusters[imodel_old]==0:
        exit(1)
    
    relative_entropy += target_probs[itarget]*math.log(target_probs[itarget]*float(omega_cg_clusters[imodel_old])/model_probs[imodel_new]) #not the same index for omega and model_probs
    
    #print itarget,imodel_old,imodel_new,target_probs[itarget],model_probs[imodel_new],relative_entropy

if rep_cg=="30":
    rep_cg="non-uni_30"
    
print bio_system,rep_high_res,rep_cg,"relative_ent_"+prob_type,relative_entropy #,relative_entropy_without_omega
    
    
    
