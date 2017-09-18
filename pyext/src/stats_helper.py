import IMP
import IMP.rmf
import RMF
import IMP.atom
import os,sys,math
import numpy as np
import glob
import datetime
import subprocess

def get_sampling_time(machine_run_on,expt):
    
   
    if machine_run_on == "bass": # on local machine things were returned directly from 
        sampling_time_file = open(os.path.join('average_sampling_time_res'+optimalRepresentationDir.lstrip('r')+'_expt'+expt+'.txt'),'r')
        
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
                    
                           
    for dc in [all_beads_stats,data_beads_stats]:
        dc["avg_num_residues"]=sum(dc["residues_per_bead"])/float(len(dc["residues_per_bead"]))
        dc["avg_dia"] = sum(dc["diameter_per_bead"])/float(len(dc["diameter_per_bead"]))
        dc["avg_weighted_dia"] = float(sum([r*d for (r,d) in zip(dc["residues_per_bead"],dc["diameter_per_bead"])]))/float(sum(dc["residues_per_bead"]))

    return all_beads_stats,data_beads_stats
    
