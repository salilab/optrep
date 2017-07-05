
import IMP
import IMP.optrep
import IMP.optrep.BeadMapBuilder
import IMP.optrep.GoodScoringModelSelector
import IMP.optrep.SamplingPrecisionEstimator
import os,sys,string,math
import argparse

def parse_args():
    
    parser = argparse.ArgumentParser(description="Get the next iteration's representation (bead map), based on results (precision) of the previous iteration. Example of usage: (1) When creating a new beadmap: set_next_representation.py -r 1 -tf <topology_file> (2) When you have previously created a beadmap and need to update it based on the results of sampling:  set_next_representation.py -r <next_resolution> -bmf <previous_beadmap_file> -pf <precions_file>.Flag -h for more details.")
    
    parser.add_argument("-r","--next_resolution",dest="next_resolution",required=True,type=int,help="The coarse-graining, i.e. maximum number of residues per bead to be tried in the next iteration.")
    
    parser.add_argument("-tf","--topology_file",action="store",dest="topology_file",help="initial topology file") 
    
    parser.add_argument("-bmf","--previous_bead_map_file",action="store",dest="prev_bm_file",help="bead map file created in previous iteration")  
                     
    parser.add_argument("-pf","--precisions_file",dest="precisions_file",help="file containing identity of imprecise beads") 
    
    result = parser.parse_args()
 
    return result

    
def set_next_beadmap():
     # Authors: Shruthi Viswanath
     
    ''' Set the representation for the next iteration of CG.
    '''
  
    # process input
    args=parse_args()
    
  
    if args.next_resolution==1: 
        
        # first iteration. Get bead map for the highest resolution beads.
        bmb=IMP.optrep.BeadMapBuilder.BeadMapBuilder()
        
        bmb.set_bead_map_from_topology_file(topology_file=args.topology_file)
           
    else:
        
        bmb=IMP.optrep.BeadMapBuilder.BeadMapBuilder()
        
        bmb.set_bead_map_from_beadmap_file(arg.previous_bead_map_file)
        
        bmb.set_imprecise_beads_from_file(arg.precisions_file)
        
        updated = bmb.update_all_bead_maps(int(args.next_resolution))
      
        # make sure the proteins for which we are not changing representation remain exactly same as before
        if not updated:
            print "Process converged."
            return
            
    next_iteration_sampling_dir = "r"+str(args.next_resolution)
    
    os.mkdir(next_iteration_sampling_dir)
    
    #Irrespective of whether we are creating a new bead map or updating an existing one, write the bead map to a file
    bmb.write_bead_map_to_file(os.path.join(next_iteration_sampling_dir,'bead_map_'+str(args.next_resolution)+'.txt'))

if __name__ == "__main__" :
    set_next_beadmap()
