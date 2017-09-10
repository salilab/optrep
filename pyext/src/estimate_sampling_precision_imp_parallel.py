
import IMP
import IMP.optrep
#import IMP.optrep.SamplingPrecisionEstimator
import os,sys,string,math
import argparse
import IMP.parallel
import IMP.optrep.parallel_tasks
from operator import itemgetter

def running_on_cluster():
    import distutils.spawn
    return distutils.spawn.find_executable('qsub') is not None

def parse_args():
    
    parser = argparse.ArgumentParser(description="Get the sampling precision per bead and identify imprecise beads. Usage: estimate_sampling_precision.py -n <number of cores> -pl <protein list: e.g. A B> (Number of elements in protein and domain lists should match.)  -dl <domain list: e.g. A_1 B_1> -rd <run_directory_for_sampling> -tf <topology_file> -gs <grid_size> -xs <xscale> -o <bead precisions file>. Flag -h for more details.")
    
    
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

    abof = open(all_beads_output_file,'w')
  
    for (prot,dom) in bead_precisions:
        
        for bead_data in bead_precisions[(prot,dom)]:
            print >>abof,prot,dom,bead_data[0],bead_data[1],bead_data[2]
    abof.close()
  
    
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
    #spe=IMP.optrep.SamplingPrecisionEstimator.SamplingPrecisionEstimator(arg.run_dir,components_to_update) # the components to calculate precision are the 
    ### same as the components to update
        
    #spe.load_coordinates_and_bead_sizes_from_model_files()
        
    #spe.estimate_perbead_sampling_precision(grid_size=arg.grid_size)
    
    #spe.get_imprecise_beads(xscale=arg.xscale)

    #spe.print_bead_precisions(outfile=arg.output_file)
    
    # Set up a Manager to keep track of slaves and our tasks
    m = IMP.parallel.Manager()

    #m = IMP.parallel.Manager(python="module load sali-libraries; module load openmpi-1.6-nodlopen; ~/imp-clean/build/setup_environment.sh python",host="172.18.1.209")

    # Add slaves 
    if running_on_cluster():
        s = IMP.parallel.SGEQsubSlaveArray(arg.num_cores,'-l arch=linux-x64 -q lab.q -l hostname="i*"')
    
    else:
        for i in range(arg.num_cores):
            s = IMP.parallel.LocalSlave()
    
    m.add_slave(s)

    # Generate a context (an environment on each slave in which tasks will be
    # run). Provide a setup function for this context. 
    c = m.get_context(IMP.optrep.parallel_tasks.slave_setup)

    num_global_beads = IMP.optrep.parallel_tasks.master_setup(components_to_update,arg.run_dir,arg.topology_file)

    # Add tasks with different input parameters
    num_beads_per_core=math.ceil(float(num_global_beads)/float(arg.num_cores))

    start_bead=0
    for i in range(arg.num_cores):
        ##print start_bead,min(start_bead+num_beads_per_core-1,num_global_beads-1)
        
        c.add_task(IMP.optrep.parallel_tasks.SlaveTask(components_to_update,arg.topology_file,arg.run_dir,arg.grid_size,arg.xscale,start_bead,min(start_bead+num_beads_per_core-1,num_global_beads-1)))
        start_bead=start_bead+num_beads_per_core
    

    # Run all tasks, distributed between the slaves. Get the results in
    # the order they are returned (not necessarily the order they were created).
    for bead_precision_sublist in c.get_results_unordered():
        for bead_string in bead_precision_sublist:
            (prot,dom,bead_index,sampling_precision,is_imprecise)=bead_string.strip().split()  
            bead_precisions[(prot,dom)].append((int(bead_index),float(sampling_precision),int(is_imprecise)))

    # after all the output is obtained in random order, need to sort them by bead and print the new beadmap file.
    bead_precisions = sort_beads_by_index(bead_precisions)
    
    # outputs 2 files, one with all precisions and one with imprecise beads only, in a format that the BeadMapBuilder can understand
    output_precisions_to_file(bead_precisions,arg.beads_output_file) 
    
        
 
if __name__ == "__main__" :
    estimate_sampling_precision()
