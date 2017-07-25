# this file contains the setup and tasks for calculating per-bead sampling precision.
#

import IMP
import IMP.optrep

def slave_setup():
    """ 
    Given a set of protein/domain components, load all their coordinates from models into memory.
    Use slaves in a parallel environment to get the sampling precision of each bead. 
    """
    components_to_update=IMP.optrep.ProteinDomainList([("B","B_1")])

    # outer () is for argument, [] is for vector, inner () is for string pair.

    spe=IMP.optrep.SPE("input/1SYX/1SYX.topology.txt","input/1SYX/good_scoring_models/",components_to_update)

    spe.load_coordinates_and_bead_sizes_from_model_files()
    
    return(spe)

def master_setup():
    """ 
    Given a set of protein/domain components, find the total number of beads across proteins.
    """ 
    components_to_update=IMP.optrep.ProteinDomainList([("B","B_1")])

    # outer () is for argument, [] is for vector, inner () is for string pair.

    spe=IMP.optrep.SPE("input/1SYX/1SYX.topology.txt","input/1SYX/good_scoring_models/",components_to_update)
    
    spe.load_coordinates_and_bead_sizes_from_model_files(True) # loads only first model to get the number of beads in the system. 
    
    return(spe.number_of_global_beads_)
    
# Note that setup and tasks are Python callables, i.e. functions (like setup
# above) or classes that implement the __call__ method (like Task below).
# The latter allows for parameters (Python objects) to be passed from the
# master to the slaves.


class SlaveTask(object):

    def __init__(self,global_bead_index):
        self.global_bead_index = global_bead_index

    def __call__(self, spe):
        """Note that the
        input parameters to this method (spe) are those returned by
        the setup function above."""
        
        bead_precision_output = spe.estimate_and_print_single_bead_precision(self.global_bead_index,
        2.0, 1.0)
        print bead_precision_output

        #return(bead_precision_output)
        return(1.0)

