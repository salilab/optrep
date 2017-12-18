import sys
sys.dont_write_bytecode=True

import IMP
import IMP.optrep
#import IMP.optrep.SamplingPrecisionEstimator
import os,string,math
import argparse

def parse_args():
    
    parser = argparse.ArgumentParser(description="Get the sampling precision per bead and identify imprecise beads. Usage: estimate_sampling_precision.py -n <number of cores for the job> -cn <current core number (starting from 0)> -pl <protein list: e.g. A B> (Number of elements in protein and domain lists should match.)  -dl <domain list: e.g. A_1 B_1> -rd <run_directory_for_sampling> -tf <topology_file> -gs <grid_size> -xs <xscale> -ao <all beads precisions file> -io <imprecise bead precisions file>. Flag -h for more details.")
    
    parser.add_argument("-n","--num_cores",dest="num_cores",type=int,help="the total number of cores in the multicore job")
    
    parser.add_argument("-cn","--core_number",dest="core_number",type=int,help="the SGE TASK ID-1 of the job")
    
    parser.add_argument("-pl","--protein_list",dest="proteins_to_update",nargs="+",type=str,help="proteins to calcualte precision on")
    
    parser.add_argument("-dl","--domain_list",dest="domains_to_update",nargs="+",type=str,help="domain to calculate precision on")

    parser.add_argument("-rd","--run_directory",dest="run_dir",help="directory in which good scoring models from sampling are stored") 
    
    parser.add_argument("-tf","--topology_file",dest="topology_file",help="file created before running sampling, that stores info by protein/domain")
    
    parser.add_argument("-gs","--grid_size",dest="grid_size",type=float,help="Size of grid for clustering and determining sampling precision") 
    
    parser.add_argument("-xs","--xscale",dest="xscale",type=float,help="Sampling precision can be utmost xscale times the representation precision plus linear cutoff") 
    parser.add_argument("-lc","--linear_cutoff",dest="linear_cutoff",type=float,help="Sampling precision can be utmost xscale times the representation precision plus linear cutoff")
     
    
    parser.add_argument("-o","--outfile_prefix",dest="output_file_prefix",help="Prefix of file where precision values of all beads is written")
    
    result = parser.parse_args()
 
    return result


 
def get_number_of_global_beads(components_to_update,run_dir,topology_file):
    """ 
    Given a set of protein/domain components, find the total number of beads across proteins.
    """ 
    # outer () is for argument, [] is for vector, inner () is for string pair.
    
    components_to_update=IMP.optrep.ProteinDomainList(components_to_update)

    spe=IMP.optrep.SPE(topology_file,run_dir,components_to_update)
    
    spe.load_coordinates_and_bead_sizes_from_model_files(True) # loads only first model to get the number of beads in the system. 
    
    return(spe.number_of_global_beads_) 
    
def estimate_sampling_precision():
     # Authors: Shruthi Viswanath
     
    # process input
    arg=parse_args()

    ### list of (protein, domain) tuples corresponding to topology file
    ### could be changed to PMI selection later
    components_to_update=[]
   
    for prot,dom in zip(arg.proteins_to_update,arg.domains_to_update):
        components_to_update.append((prot,dom))
   
    num_global_beads = get_number_of_global_beads(components_to_update,arg.run_dir,arg.topology_file)

    # Add tasks with different input parameters
    num_beads_per_core=math.ceil(float(num_global_beads)/float(arg.num_cores))

    start_bead=0
    for i in range(arg.core_number):
        start_bead=start_bead+num_beads_per_core
    
    end_bead=int(min(start_bead+num_beads_per_core-1,num_global_beads-1))
    
    spe=IMP.optrep.SPE(arg.topology_file,arg.run_dir,components_to_update)
    
    spe.load_coordinates_and_bead_sizes_from_model_files()
        
    spe.print_to_file_precision_for_range_of_beads(start_bead, end_bead,
    arg.grid_size, arg.xscale,arg.linear_cutoff,arg.output_file_prefix+"."+str(arg.core_number))
       
 
if __name__ == "__main__" :
    estimate_sampling_precision()
