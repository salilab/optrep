
# this file contains the setup and tasks for calculating per-bead sampling precision.
#

import IMP
import IMP.optrep

def slave_setup():
    
    
    return()

def master_setup(components_to_update, topology, gsm):
    """ 
    Given a set of protein/domain components, find the total number of beads across proteins.
    """ 
    # outer () is for argument, [] is for vector, inner () is for string pair.

    spe=IMP.optrep.SPE(topology, gsm, components_to_update)
    
    spe.load_coordinates_and_bead_sizes_from_model_files(True) # loads only first model to get the number of beads in the system. 
    
    return(spe.number_of_global_beads_)
    
# Note that setup and tasks are Python callables, i.e. functions (like setup
# above) or classes that implement the __call__ method (like Task below).
# The latter allows for parameters (Python objects) to be passed from the
# master to the slaves.


class SlaveTask(object):

    def __init__(self,grid_size, xscale,start_bead_index,end_bead_index,
                 topology, gsm):
        self.grid_size = grid_size
        self.xscale=xscale
        self.start_bead_index=start_bead_index
        self.end_bead_index=max(end_bead_index, start_bead_index)
        self.topology, self.gsm = topology, gsm
     
        
    def __call__(self):
        """Note that the
        input parameters to this method are those returned by
        the setup function above.
        """
        
        """ 
        Given a set of protein/domain components, load all their coordinates from models into memory.
        Use slaves in a parallel environment to get the sampling precision of each bead. 
        """
        
        components_to_update=IMP.optrep.ProteinDomainList([("B","B_1")])
        
        
        spe=IMP.optrep.SPE(self.topology, self.gsm, components_to_update)
    
        spe.load_coordinates_and_bead_sizes_from_model_files()
        
        bead_precision_output_list = spe.print_precision_for_range_of_beads(self.start_bead_index, self.end_bead_index,
self.grid_size, self.xscale, linear_cutoff=1.0)
                        
        return(bead_precision_output_list)
        #for out_str in bead_precision_output_list:
            #print out_str
            
        # return(1.0)

