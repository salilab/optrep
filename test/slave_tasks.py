# this file contains the setup and tasks for calculating per-bead sampling precision.
#

import IMP
import IMP.optrep
def setup():
    """ 
    Given a set of protein/domain components, load all their coordinates from models into memory.
    Use slaves in a parallel environment to get the sampling precision of each bead. 
    """
    
    
    return spe

# Note that setup and tasks are Python callables, i.e. functions (like setup
# above) or classes that implement the __call__ method (like Task below).
# The latter allows for parameters (Python objects) to be passed from the
# master to the slaves.


class SlaveTask(object):

    def __init__(self,global_bead_index):
        self.global_bead_index = global_bead_index

    def __call__(self, spe):
           """Note that the
           input parameters to this method (m, sf, and d) are those returned by
           the setup function above."""
        
        spe.IMP.algebra.Vector3D(0, 0, self.dist))
        return (self.dist, sf.evaluate(False))
