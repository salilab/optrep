import IMP
import IMP.optrep
#import IMP.optrep.SamplingPrecisionEstimator
import os,sys,string,math
import argparse
from operator import itemgetter

def parse_args():
    
    parser = argparse.ArgumentParser(description="Collect bead precisions from outputs of parallel jobs, sort them and output in a single file. Usage: collate_sampling_precision.py -n <number_of_cores> -pl <protein list: e.g. A B> (Number of elements in protein and domain lists should match.)  -dl <domain list: e.g. A_1 B_1> -ip <input file prefix> -o <output file>. Flag -h for more details.")
    
    parser.add_argument("-n","--number_of_cores",dest="num_cores",type=int,help="the total number of cores for the parallel sampling precision estimation")
    
    parser.add_argument("-pl","--protein_list",dest="proteins_to_update",nargs="+",type=str,help="proteins to calcualte precision on")
    
    parser.add_argument("-dl","--domain_list",dest="domains_to_update",nargs="+",type=str,help="domain to calculate precision on")
  
    parser.add_argument("-ip","--input_file_prefix",dest="input_file_prefix",type=str,help="Prefix of files where precision values of bead (sublists) are written")
    
    parser.add_argument("-o","--output_file",dest="output_file",type=str,help="Name of single output file containing info of all beads in order")
    
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
  
def populate_bead_precisions(bead_precisions_dict,num_cores,input_file_prefix):
    
    for i in range(num_cores):
        
        sublist_file = open(input_file_prefix+"."+str(i),'r')
        
        for bead_string in sublist_file.readlines():
            (prot,dom,bead_index,sampling_precision,is_imprecise)=bead_string.strip().split()  
            bead_precisions_dict[(prot,dom)].append((int(bead_index),float(sampling_precision),int(is_imprecise)))
        
        sublist_file.close()    
    
    return bead_precisions_dict
          
def collate_sampling_precisions():
    ''' Collect bead precisions from outputs of parallel jobs, sort them and output in a single file.
    '''

    arg = parse_args()
    # initialize the dictionary of bead precisions which stores a list of bead precisions for every bead
    bead_precisions = {}
    for prot,dom in zip(arg.proteins_to_update,arg.domains_to_update):
        bead_precisions[(prot,dom)]=[]
 
    # collect and store bead precisions
    bead_precisions=populate_bead_precisions(bead_precisions,arg.num_cores,arg.input_file_prefix)
    
    # after all the output is obtained in random order, need to sort them by bead and print the new beadmap file.
    bead_precisions = sort_beads_by_index(bead_precisions)
    
    # outputs 2 files, one with all precisions and one with imprecise beads only, in a format that the BeadMapBuilder can understand
    output_precisions_to_file(bead_precisions,arg.output_file) 
        
 
if __name__ == "__main__" :
    collate_sampling_precisions()
