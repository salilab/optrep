from __future__ import print_function
import IMP
import IMP.optrep
import IMP.optrep.BeadMapBuilder
import IMP.optrep.GoodScoringModelSelector
import IMP.optrep.SamplingPrecisionEstimator
import os,sys,string,math
import argparse

def parse_args():
    
    parser = argparse.ArgumentParser(description="Get the next iteration's representation (bead map), based on results (precision) of the previous iteration. Example of usage: (1) When creating a new beadmap: set_next_representation.py -u create -r 1 -tf <topology_file> (2) When you have previously created a beadmap and need to update it based on the results of sampling:  set_next_representation.py -u update -r <next_resolution> -bmf <previous_beadmap_file> -pf <precions_file>.Flag -h for more details.")
    
    parser.add_argument("-u","--usage",dest="usage",required=True,type=str,help="Say create for creating a beadmap from scratch from the topology file and update if we are in a future iteration and are updating a beadmap")

    parser.add_argument("-r","--next_resolution",dest="next_resolution",required=True,type=int,help="The coarse-graining, i.e. maximum number of residues per bead to be tried in the next iteration.")
    
    parser.add_argument("-tf","--topology_file",action="store",dest="topology_file",help="initial topology file") 
    
    parser.add_argument("-bmf","--previous_bead_map_file",action="store",dest="previous_bead_map_file",help="bead map file created in previous iteration")  
                     
    parser.add_argument("-pf","--precisions_file",action="store",dest="precisions_file",help="file containing identity of imprecise beads") 
    
    result = parser.parse_args()
 
    return result

    
def set_next_beadmap():
     # Authors: Shruthi Viswanath
     
    ''' Set the representation for the next iteration of CG.
    '''
  
    # process input
    args=parse_args()
    
  
    if args.usage=="create": 
        
        # first iteration. Get bead map for the highest resolution beads.
        bmb=IMP.optrep.BeadMapBuilder.BeadMapBuilder()
        
        #TODO Note : it assumes that we use the same starting resolution for all proteins which may not be true! 
        bmb.set_bead_map_from_topology_file(topology_file=args.topology_file,resolution=args.next_resolution) 
           
    else:
        
        bmb=IMP.optrep.BeadMapBuilder.BeadMapBuilder()
        
        bmb.set_bead_map_from_beadmap_file(args.previous_bead_map_file)
        
        bmb.set_imprecise_beads_from_file(args.precisions_file)
        
        #TODO Note : it assumes that we use the same resolution for all proteins which may not be true!
        updated = bmb.update_all_bead_maps(int(args.next_resolution))
      
        # make sure the proteins for which we are not changing representation remain exactly same as before
        if not updated:
            print("Process converged.")
            return
            
    #next_iteration_sampling_dir = "r"+str(args.next_resolution)
    
    #os.mkdir(next_iteration_sampling_dir)
    
    ##Irrespective of whether we are creating a new bead map or updating an existing one, write the bead map to a file
    #bmb.write_bead_map_to_file(os.path.join(next_iteration_sampling_dir,'bead_map_'+str(args.next_resolution)+'.txt'))
    
    bmb.write_bead_map_to_file('bead_map_'+str(args.next_resolution)+'.txt')

if __name__ == "__main__" :
    set_next_beadmap()
