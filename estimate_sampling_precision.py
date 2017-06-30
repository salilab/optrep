
import IMP
import IMP.optrep
import IMP.optrep.BeadMapBuilder
import IMP.optrep.GoodScoringModelSelector
import IMP.optrep.SamplingPrecisionEstimator
import os,sys,string,math
import argparse

def parse_args():
    
    parser = argparse.ArgumentParser(description="Get the next iteration's representation (bead map), based on results (precision) of the previous iteration. Example of usage: (1) When creating a new beadmap: set_next_representation.py -r 1 -tf <topology_file> (2) When you have previously created a beadmap and need to update it based on the results of sampling:  set_next_representation.py -r <next_resolution> -pf <previous_beadmap_file> -rd <run_directory_for_sampling> -rp <run_prefix> -gs <grid_size> -xs <xscale>. Flag -h for more details.")
    
    parser.add_argument("-r","--next_resolution",dest="next_resolution",required=True,type=int,help="The coarse-graining, i.e. maximum number of residues per bead to be tried in the next iteration.")
    
    parser.add_argument("-tf","--topology_file",action="store",dest="topology_file",help="initial topology file") 
    
    parser.add_argument("-pf","--previous_bead_map_file",action="store",dest="prev_bm_file",help="bead map file created in previous iteration")  
                     
    parser.add_argument("-rd","--run_directory",dest="run_directory",help="directory in whcih sampling results are stored") 
    
    parser.add_argument("-rp","--run_prefix",dest="run_prefix",help="prefix of runs") 
                        
    parser.add_argument("-gs","--grid_size",dest="grid_size",help="Size of grid for clustering and determining sampling precision") 
    
    parser.add_argument("-xs","--xscale",dest="xscale",help="Sampling precision can be utmost xscale times the representation precision plus 2.0") 
                 
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
        
        gsms=IMP.optrep.GoodScoringModelSelector.GoodScoringModelSelector(arg.run_dir,arg.run_prefix)
        
        ##gsms.get_good_scoring_models(criteria_list=[],keywords_list=[],aggregate_lower_) #TODO specific to application
        
        ### list of (protein, domain) tuples corresponding to topology file
        ### could be changed to PMI selection later
        ##components_to_update=[("A","A_1"),("B","B_1")] #TODO specific to application
        
        #spe=IMP.optrep.SamplingPrecisionEstimator.SamplingPrecisionEstimator(os.path.join(arg.run_dir,"good_scoring_models"),arg.components_to_update) # the components to calculate precision are the 
        ## same as the components to update
        
        #spe.load_coordinates_and_bead_sizes_from_model_files()
        
        #spe.estimate_perbead_sampling_precision(grid_size=arg.grid_size)
        
        #imprecise_beads=spe.get_imprecise_beads(xscale=arg.xscale)
        
        #bmb=IMP.optrep.BeadMapBuilder.BeadMapBuilder()
        
        #bmb.set_bead_map_from_beadmap_file(arg.previous_bead_map_file)
        
        #bmb.set_imprecise_beads(imprecise_beads)
        
        #updated = bmb.update_all_bead_maps(components_to_update,int(args.next_resolution))
      
        ## make sure the proteins for which we are not changing representation remain exactly same as before
        #if not updated:
            #print "Process converged."
            #return
            
    next_iteration_sampling_dir = "r"+str(args.next_resolution)
    
    os.mkdir(next_iteration_sampling_dir)
    
    #Irrespective of whether we are creating a new bead map or updating an existing one, write the bead map to a file
    bmb.write_bead_map_to_file(os.path.join(next_iteration_sampling_dir,'bead_map_'+str(args.next_resolution)+'.txt'))

if __name__ == "__main__" :
    set_next_beadmap()
