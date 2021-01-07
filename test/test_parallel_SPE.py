from __future__ import print_function
import IMP.test
import IMP
import IMP.atom
import IMP.rmf
import IMP.pmi
import IMP.pmi.tools
import IMP.pmi.topology
import os,sys,string,math
import IMP.optrep
import IMP.parallel
import parallel_tasks

class Tests(IMP.test.TestCase):
    @IMP.test.skip("takes too long to run")
    def test_parallel_spe(self):
        """Test parallel SPE"""
        # Hard-code number of cores for test
        num_cores = 2

        components_to_update=IMP.optrep.ProteinDomainList([("B","B_1")])
        grid_size=2.0
        xscale=1.0

        # Set up a Manager to keep track of slaves and our tasks
        m = IMP.parallel.Manager()

        # Add slaves  
        for i in range(num_cores):
            s = IMP.parallel.LocalWorker()
            m.add_worker(s)
            
        # Generate a context (an environment on each slave in which tasks will be
        # run). Provide a setup function for this context. 
        c = m.get_context(parallel_tasks.slave_setup)

        topology = self.get_input_file_name('1SYX/1SYX.topology.txt')
        gsm = self.get_input_file_name('1SYX/good_scoring_models')
        num_global_beads = parallel_tasks.master_setup(components_to_update,
            topology, gsm)

        # Add tasks with different input parameters
        num_beads_per_core=math.ceil(float(num_global_beads)/float(num_cores))

        start_bead=0
        for i in range(num_cores):
            print(start_bead,min(start_bead+num_beads_per_core-1,num_global_beads-1))
            
            c.add_task(parallel_tasks.SlaveTask(grid_size,xscale,start_bead,min(start_bead+num_beads_per_core-1,num_global_beads-1), topology, gsm))
            start_bead=start_bead+num_beads_per_core
           

        # Run all tasks, distributed between the slaves. Get the results in
        # the order they are returned (not necessarily the order they were created).

        #for bead_precision_list in c.get_results_unordered():
           #for out_string in bead_precision_list:
               #print out_string
               
        for x in c.get_results_unordered():
            print(x)


if __name__ == '__main__':
    IMP.test.main()
