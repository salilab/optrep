import IMP
import IMP.rmf
import RMF
import IMP.atom
import os,sys,math
import numpy 
import glob
import datetime
import subprocess
import scipy.stats

def get_sampling_precision(model_sample_id_file,proteins_list,grid_size):
        # get the model index of run1 and run2 models separately
        run1_all_models,run2_all_models=get_run_identity(model_sample_id_file)
        all_models=run1_all_models+run2_all_models

        total_num_models = len(all_models)

        # get cutoffs from all vs all distances
        distmat_full=get_all_vs_all_rmsd(proteins_list)

        # parameters to set
        cutoffs_list=get_cutoffs_list(distmat_full,grid_size,total_num_models)

        sampling_precision=100000.0

        for c in cutoffs_list:
            cluster_centers,cluster_members=precision_cluster(distmat_full,total_num_models,c)

            ctable,retained_clusters=get_contingency_table(len(cluster_centers),cluster_members,all_models,run1_all_models,run2_all_models)
            
            (pval,cramersv)=test_sampling_convergence(ctable,total_num_models)

            percent_explained= percent_ensemble_explained(ctable,total_num_models)
            
            if percent_explained>80.0 and (pval>0.05 or cramersv<0.10):
                if sampling_precision > c:
                    sampling_precision=c
            else:
                sampling_precision = 100000.0

        return sampling_precision
    
def get_all_vs_all_rmsd(proteins_list):
    
    #load coordinates
    conform=[]
    numModels=0
    weight = []
    
    for fl in glob.glob("*.rmf3"):
        partcoord = []
        m = IMP.Model()
        inf = RMF.open_rmf_file_read_only(fl)
        h = IMP.rmf.create_hierarchies(inf, m)[0]
        IMP.rmf.load_frame(inf, 0)
      
        for state in h.get_children():
            for component in state.get_children():
                if included_protein(component.get_name(),proteins_list):
                    
                    for leaf in IMP.core.get_leaves(component):
                        
                        p=IMP.core.XYZ(leaf.get_particle())
                        partcoord.append(p.get_coordinates())
                        
                        if numModels==0:
                            if len(IMP.atom.Fragment(leaf).get_residue_indexes()) == 0:
                                weight.append(1.0)
                            else:
                                weight.append(float(len(IMP.atom.Fragment(leaf).get_residue_indexes())))

        conform.append(partcoord)
        numModels += 1
                   
    weight=numpy.array(weight)
  
    # Calculate the distance matrix
    distmat=numpy.zeros((numModels,numModels))

    for i in range(numModels-1):
        
        for j in range(i+1,numModels):
             dist=IMP.algebra.get_weighted_rmsd(conform[i],conform[j], weight)
             distmat[i][j]=dist
    
    return distmat

def get_run_identity(sample_id_file):
    # whether a run belongs to run1 or run2
        
    run1_models=[]
    run2_models=[]
    sif=open(sample_id_file,'r')

    for ln in sif.readlines():
        model_index=int(ln.strip().split()[0])
        sample_number=int(ln.strip().split()[1])
    
        if sample_number==1:
            run1_models.append(model_index)
        elif sample_number==2:
            run2_models.append(model_index)
    
    sif.close()
   
    return run1_models,run2_models

def get_cutoffs_list(distmat,gridSize,numModels):

	maxdist=0.0
	mindist=1000000.0
	for i in range(numModels-1):
    		for j in range(i+1,numModels):
        		if distmat[i][j]>maxdist:
            			maxdist=distmat[i][j]

        		if distmat[i][j]<mindist:
           			 mindist=distmat[i][j]

    
	cutoffs=numpy.arange(mindist,maxdist,gridSize) # or maxdist/2.0, 5.0 s

	return cutoffs

def precision_cluster(distmat,numModels,rmsd_cutoff):
    #STEP 2. Populate the neighbors ofa given model
    neighbors=[]
    for count in range(numModels):
        neighbors.append([count])  # model is a neighbor of itself
 
    for i in range(numModels-1):
        for j in range(i+1,numModels):
        
            if distmat[i][j]<=rmsd_cutoff: # accepted to be a neighbor
                #print i,j,distmat[i][j]

                neighbors[i].append(j)
                neighbors[j].append(i)
     
    #STEP 3. Get the weightiest cluster, and iterate
    unclustered=[]
    boolUnclustered=[]
    for i in range(numModels):
        unclustered.append(i)
        boolUnclustered.append(True)

    cluster_members=[] # list of lists : one list per cluster
    cluster_centers=[]

    while len(unclustered)>0:
        # get cluster with maximum weight
        max_neighbors=0
        currcenter=-1
        for eachu in unclustered:  # if multiple clusters have same maxweight this tie is broken arbitrarily! 
            if len(neighbors[eachu])>max_neighbors:
                max_neighbors=len(neighbors[eachu])
                currcenter=eachu   
   
        #form a new cluster with u and its neighbors
        cluster_centers.append(currcenter)
        cluster_members.append([n for n in neighbors[currcenter]]) 

        #update neighbors 
        for n in neighbors[currcenter]:
            #removes the neighbor from the pool
            unclustered.remove(n) #first occurence of n is removed. 
            boolUnclustered[n]=False # clustered

        for n in neighbors[currcenter]:
            for unn in neighbors[n]: #unclustered neighbor
                if not boolUnclustered[unn]:
                    continue
                neighbors[unn].remove(n)
   
    return cluster_centers,cluster_members

def get_contingency_table(num_clusters,cluster_members,models_subset,run1_models,run2_models):

	full_ctable=numpy.zeros((num_clusters,2))
		
	for ic,cluster in enumerate(cluster_members):
		for member in cluster:
			model_index=models_subset[member]

			if model_index in run1_models:
                                #print "run1",model_index
                                full_ctable[ic][0]+=1.0
			elif model_index in run2_models:
				#print "run2",model_index
                                full_ctable[ic][1]+=1.0

	## now normalize by number of models in each run
	numModelsRun1 = float(numpy.sum(full_ctable,axis=0)[0])
	numModelsRun2 = float(numpy.sum(full_ctable,axis=0)[1])

  	reduced_ctable=[]
	
	retained_clusters=[]

	for i in range(num_clusters):
		if full_ctable[i][0]<=10.0 or full_ctable[i][1]<=10.0:
			continue
		reduced_ctable.append([full_ctable[i][0],full_ctable[i][1]])
		retained_clusters.append(i)

	return numpy.array(reduced_ctable),retained_clusters

def test_sampling_convergence(contingency_table,total_num_models):

    if len(contingency_table)==0:
        return 0.0,1.0

    ct = numpy.transpose(contingency_table) 
    
    [chisquare,pvalue,dof,expected]=scipy.stats.chi2_contingency(ct)
    
    if dof==0.0:
            cramersv=0.0 #converged, one single cluster
    else:
            cramersv=math.sqrt(chisquare/float(total_num_models))
                
    return(pvalue,cramersv)

def percent_ensemble_explained(ctable,total_num_models):
       
        if len(ctable)==0:
            return 0.0
        percent_clustered=float(numpy.sum(ctable,axis=0)[0]+numpy.sum(ctable,axis=0)[1])*100.0/float(total_num_models)
        
        return percent_clustered
    
class Violations(object):

    def __init__(self, threshold,xlink_keyword,example_file):

        self.violation_threshold  = threshold 
        self.xlink_keyword=xlink_keyword
             
        self.required_keys={}
        
        ef=open(example_file,'r')
        header_dict = eval(ef.readlines()[0].strip())
        for k in header_dict:
            if header_dict[k].startswith(xlink_keyword):
                items = header_dict[k].split("|")
                (res1,pos1,res2,pos2) = items[3:7]
            
                self.required_keys[k]=(res1,pos1,res2,pos2)
        
        ef.close()
  
    def get_xlink_distances(self, lndict):
        # get the minimum distances for each xlink'ed pair in a single model

        distances = {}
       
        for k in self.required_keys:
            distances[self.required_keys[k]]=[]
      
        for rst in lndict:
            if not rst in self.required_keys:
                continue
                
            distances[self.required_keys[rst]] = float(lndict[rst])
       
        return distances
   
def get_run_replica_model_numbers(ids_file):
    
    run_replica_model_ids=[]
    idf=open(ids_file,'r')
    for ln in idf.readlines():
        fields = ln.strip().split()
        run_replica_model_ids.append((fields[1],fields[2],fields[3]))
 
    idf.close()
    
    return run_replica_model_ids
    
def get_fit_to_xlinks(xlink_threshold,xlink_keyword):
  
    gsm_ids_list = get_run_replica_model_numbers(os.path.join('good_scoring_models','model_ids_scores.txt'))
    
    Analysis = Violations(xlink_threshold,xlink_keyword,'run.1/output/stat.0.out') # all other files look the same
      
    distance_between_xlink_beads={Analysis.required_keys[k]:[] for k in Analysis.required_keys}
    
    for runid,replicaid,modelid in gsm_ids_list:
            
        stat_file_line=subprocess.check_output(["awk","NR=="+str(int(modelid)+1)+" {print} ",os.path.join("run."+runid,"output","stat."+replicaid+".out")]) # with awk and subprocess, always have awk in one string, the condition and action in one string and the filename in another string, that works! 
   
        xlink_distances_model=Analysis.get_xlink_distances(eval(stat_file_line.strip('\n'))) # return a dictionary keyed by crosslinks, and values equal to values for a model
        
        for k in xlink_distances_model:
                distance_between_xlink_beads[k].append(xlink_distances_model[k])
   
    # get average over all xlinks and all models
    num_xlinks_in_all_models=0.0
    sum_xlink_distances_across_all_models=0.0
  
    for k in distance_between_xlink_beads:
        sum_xlink_distances_across_all_models += sum(distance_between_xlink_beads[k])
       
        num_xlinks_in_all_models+=len(distance_between_xlink_beads[k])
        
    avg_distances_across_all_xlinks=sum_xlink_distances_across_all_models/float(num_xlinks_in_all_models)
    
    return(avg_distances_across_all_xlinks)

def get_sampling_time(machine_run_on,expt):
    
    curr_dir=os.getcwd().split('/')[-1]
  
    if machine_run_on == "bass": # on local machine things were returned directly from 
        sampling_time_file = open(os.path.join('average_sampling_time_res'+curr_dir.lstrip('r')+'_expt'+expt+'.txt'),'r')
        
        avg_sampling_time = sampling_time_file.readlines()[0].strip()
        
        sampling_time_file.close()
        
        return(avg_sampling_time)
    
    else:
        sampling_times = []

        date_format = "%m/%d/%Y_%H:%M:%S"

        for sampling_dir in glob.glob(os.path.join('run.*')):
        
            curr_job_times=[]
            # get the output file for the job, assumes there was only 1 output file.
            
            job_out_file = glob.glob(os.path.join(sampling_dir,'s*.o[0-9]*'))
            
            if len(job_out_file)>1:
                print "Note that there are two sampling output files. Manually delete one of them!"
                print sampling_dir,job_out_file
                exit(1)
           
            jof=open(job_out_file[0],'r')
         
            for ln in jof.readlines():
                if ln.startswith("Datetime:"):
                   curr_job_times.append(ln.strip().lstrip("Datetime:"))
                    
            jof.close()
            
            start_time=datetime.datetime.strptime(curr_job_times[0],date_format)
            end_time=datetime.datetime.strptime(curr_job_times[1],date_format)
            
            #print sampling_dir, job_out_file[0],start_time,end_time
            
            diff = end_time - start_time
            
            #sampling_times.append(diff.days*86400 + diff.seconds)
            
            sampling_times.append(diff.total_seconds()) # not available in python versions lower than 2.7
        
        avg_sampling_time = sum(sampling_times)/float(len(sampling_times))
        
        return(avg_sampling_time)

def parse_config_file(config_file):
    
    cf=open(config_file,'r')
    
    configs_dict ={}
    
    for ln in cf.readlines():
    
        if ln.startswith("#") or not ln.strip(): # ignore comment and blank lines
            continue
        
        fields=ln.strip().split("=")
        
        if fields[0]=="RESOLUTIONS_LIST":
            configs_dict[fields[0]]=fields[1].split() # key to list map
            
        else :
            if fields[1].startswith('~'): # location of a file or directory
     
                fields[1]=os.path.expanduser('~')+fields[1].lstrip('~') #os.path.join did not work here for some weird reason!
      
            configs_dict[fields[0]]=fields[1] # just a single key to string map

    cf.close()    
    
    return configs_dict

def included_protein(component_name,proteins_list):
    
    for pl in proteins_list:
        if pl in component_name:
            return pl
    
    return False

def get_beads_of_xlink_residues(xlink_file,proteins_list):
    
    xlink_beads={}
    
    for pl in proteins_list:
        xlink_beads[pl]=[]
        
    xf = open(xlink_file,'r')
    for ln in xf.readlines():
        fields=ln.strip().split(',')
    
        if fields[1] in xlink_beads and fields[0] not in xlink_beads[fields[1]]:
            xlink_beads[fields[1]].append(int(fields[0]))
            
        if fields[3] in xlink_beads and fields[2] not in xlink_beads[fields[3]]:
            xlink_beads[fields[3]].append(int(fields[2]))
   
    xf.close()
    
    return xlink_beads    

def is_xlink_bead(xlink_beads_protein,leaf):
    
    if len(IMP.atom.Fragment(leaf).get_residue_indexes()) == 0:
        residue_index = IMP.atom.Residue(leaf).get_index() 
      
        if residue_index in xlink_beads_protein:
            return True
        
    else:
        for residue_index in IMP.atom.Fragment(leaf).get_residue_indexes(): 
            if residue_index in xlink_beads_protein:
                return True
    
    return False
    
def get_representation_resolution_and_precision(proteins_list,domains_list,xlink_file=None):
    # actually we are ignoring the domains list for now, but this could be useful in future
    
    all_beads_stats={"residues_per_bead":[],"diameter_per_bead":[],"avg_num_residues":0.0,"avg_dia":0.0,"avg_weighted_dia":0.0}  
    data_beads_stats={"residues_per_bead":[],"diameter_per_bead":[],"avg_num_residues":0.0,"avg_dia":0.0,"avg_weighted_dia":0.0}  
    non_data_beads_stats={"residues_per_bead":[],"diameter_per_bead":[],"avg_num_residues":0.0,"avg_dia":0.0,"avg_weighted_dia":0.0}  
    
    if xlink_file:
        xlink_beads = get_beads_of_xlink_residues(xlink_file,proteins_list) # return a dictionary by protein
   
    rmf_file = os.path.join('good_scoring_models','0.rmf3')
    
    m = IMP.Model()
    
    fl = RMF.open_rmf_file_read_only(rmf_file)
    
    h = IMP.rmf.create_hierarchies(fl,m)[0]
    
    IMP.rmf.load_frame(fl,0)
    
    for state in h.get_children():
        for component in state.get_children():
            prot = included_protein(component.get_name(),proteins_list)
            if prot: 
                
                for leaf in IMP.core.get_leaves(component):
                  
                    if len(IMP.atom.Fragment(leaf).get_residue_indexes()) == 0: # bead with only one residue
                        curr_bead_resolution = 1.0
                       
                    else:
                        curr_bead_resolution = float(len(IMP.atom.Fragment(leaf).get_residue_indexes()))
                 
                    curr_bead_dia = IMP.core.XYZR(leaf).get_radius()*2.0
                    
                    all_beads_stats["residues_per_bead"].append(curr_bead_resolution)
                    all_beads_stats["diameter_per_bead"].append(curr_bead_dia)
                    
                    if is_xlink_bead(xlink_beads[prot],leaf):
                        data_beads_stats["residues_per_bead"].append(curr_bead_resolution)
                        data_beads_stats["diameter_per_bead"].append(curr_bead_dia)
                    else:
                        non_data_beads_stats["residues_per_bead"].append(curr_bead_resolution)
                        non_data_beads_stats["diameter_per_bead"].append(curr_bead_dia)
                           
    for dc in [all_beads_stats,data_beads_stats,non_data_beads_stats]:
        dc["avg_num_residues"]=sum(dc["residues_per_bead"])/float(len(dc["residues_per_bead"]))
        dc["avg_dia"] = sum(dc["diameter_per_bead"])/float(len(dc["diameter_per_bead"]))
        dc["avg_weighted_dia"] = float(sum([r*d for (r,d) in zip(dc["residues_per_bead"],dc["diameter_per_bead"])]))/float(sum(dc["residues_per_bead"]))

    return all_beads_stats,data_beads_stats,non_data_beads_stats
    
