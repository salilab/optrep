from __future__ import print_function
import IMP
import IMP.optrep
import IMP.optrep.SamplingPrecisionEstimator
import os,sys,string,math
import argparse
import IMP.parallel
import IMP.optrep.parallel_tasks
from operator import itemgetter

def running_on_cluster():
    import distutils.spawn
    return distutils.spawn.find_executable('qsub') is not None

def parse_args():
    
    parser = argparse.ArgumentParser(description="Get the sampling precision per bead and identify imprecise beads. Usage: estimate_sampling_precision.py -n <number of cores> -pl <protein list: e.g. A B> (Number of elements in protein and domain lists should match.)  -dl <domain list: e.g. A_1 B_1> -rd <run_directory_for_sampling> -tf <topology_file> -gs <grid_size> -xs <xscale> -ao <all beads precisions file> -io <imprecise bead precisions file>. Flag -h for more details.")
    
    
    parser.add_argument("-n","--num_cores",dest="num_cores",type=int,help="number of cores to run on")
    
    parser.add_argument("-pl","--protein_list",dest="proteins_to_update",nargs="+",type=str,help="proteins to calcualte precision on")
    
    parser.add_argument("-dl","--domain_list",dest="domains_to_update",nargs="+",type=str,help="domain to calculate precision on")

    parser.add_argument("-rd","--run_directory",dest="run_dir",help="directory in which good scoring models from sampling are stored") 
    
    parser.add_argument("-tf","--topology_file",dest="topology_file",help="file created before running sampling, that stores info by protein/domain")
    
    parser.add_argument("-gs","--grid_size",dest="grid_size",type=float,help="Size of grid for clustering and determining sampling precision") 
    
    parser.add_argument("-xs","--xscale",dest="xscale",type=float,help="Sampling precision can be utmost xscale times the representation precision plus 2.0") 
    
    parser.add_argument("-o","--beads_outfile",dest="beads_output_file",help="File where precision values of all beads is written")
    
    result = parser.parse_args()
 
    return result

def sort_beads_by_index(bead_precisions_dict):
    
    for (prot,dom) in bead_precisions_dict:
        bead_precisions_dict[(prot,dom)]=sorted(bead_precisions_dict[(prot,dom)],key=itemgetter(0))
    
    return bead_precisions_dict

def output_precisions_to_file(bead_precisions,all_beads_output_file):

    with open(all_beads_output_file, 'w') as abof:
        for (prot,dom) in bead_precisions:
            for bead_data in bead_precisions[(prot,dom)]:
                print(prot,dom,bead_data[0],".2f" %(bead_data[1]),bead_data[2],
                      file=abof)
  
    
def estimate_sampling_precision():
     # Authors: Shruthi Viswanath
     
    # process input
    arg=parse_args()


    ### list of (protein, domain) tuples corresponding to topology file
    ### could be changed to PMI selection later
    components_to_update=[]
    # initialize the dictionary of bead precisions which stores a list of bead precisions for every bead
    bead_precisions = {}
    for prot,dom in zip(arg.proteins_to_update,arg.domains_to_update):
        components_to_update.append((prot,dom))
        bead_precisions[(prot,dom)]=[]
    
    #Earlier very slow code!
    spe=IMP.optrep.SamplingPrecisionEstimator.SamplingPrecisionEstimator(arg.run_dir,components_to_update) # the components to calculate precision are the 
    ### same as the components to update
        
    spe.load_coordinates_and_bead_sizes_from_model_files()
        
    spe.estimate_perbead_sampling_precision(grid_size=arg.grid_size)
    
    spe.get_imprecise_beads(xscale=arg.xscale)

    spe.print_bead_precisions(outfile=arg.output_file)
    
 
if __name__ == "__main__" :
    estimate_sampling_precision()
