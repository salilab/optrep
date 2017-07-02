
import IMP
import IMP.optrep
import IMP.optrep.SamplingPrecisionEstimator
import os,sys,string,math
import argparse

def parse_args():
    
    parser = argparse.ArgumentParser(description="Get the sampling precision per bead and identify imprecise beads. Usage: estimate_sampling_precision.py -pl <protein list: e.g. A B> (Number of elements in protein and domain lists should match.)  -dl <domain list: e.g. A_1 B_1> -rd <run_directory_for_sampling> -gs <grid_size> -xs <xscale>. Flag -h for more details.")
    parser.add_argument("-pl","--protein_list",dest="proteins_to_update",nargs="+",type=string,help="proteins to calcualte precision on")
    
    parser.add_argument("-dl","--domain_list",dest="domains_to_update",nrgs="+",type=string,help="domain to calculate precision on")

    parser.add_argument("-rd","--run_directory",dest="run_directory",help="directory in which sampling results are stored") 
    
    parser.add_argument("-gs","--grid_size",dest="grid_size",help="Size of grid for clustering and determining sampling precision") 
    
    parser.add_argument("-xs","--xscale",dest="xscale",help="Sampling precision can be utmost xscale times the representation precision plus 2.0") 
    
    parser.add_argument("-o","--outfile",dest="output_file",help="File where precision values and identity of imprecise beads is written")
                 
    result = parser.parse_args()
 
    return result
    
def estimate_sampling_precision():
     # Authors: Shruthi Viswanath
     
    # process input
    args=parse_args()
    
    ### list of (protein, domain) tuples corresponding to topology file
    ### could be changed to PMI selection later
    components_to_update=[]
    for prot,dom in zip(arg.proteins_to_update,arg,domains_to_update):
        components_to_update.append((prot,dom))
        
    spe=IMP.optrep.SamplingPrecisionEstimator.SamplingPrecisionEstimator(os.path.join(arg.run_dir,"good_scoring_models"),components_to_update) # the components to calculate precision are the 
    ## same as the components to update
        
    spe.load_coordinates_and_bead_sizes_from_model_files()
        
    spe.estimate_perbead_sampling_precision(grid_size=arg.grid_size)
    
    spe.get_imprecise_beads(xscale=arg.xscale)

    spe.print_bead_precisions(outfile=arg.output_file)

if __name__ == "__main__" :
    estimate_sampling_precision()
