
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

def running_on_cluster():
    import distutils.spawn
    return distutils.spawn.find_executable('qsub') is not None

numCores=int(sys.argv[1])

# Set up a Manager to keep track of slaves and our tasks
m = IMP.parallel.Manager()

# Add slaves  
for i in range(numCores):
    
    if running_on_cluster():
        s = IMP.parallel.SGEQsubSlaveArray()
    else:
        s = IMP.parallel.LocalSlave()
    m.add_slave(s)
    
# Generate a context (an environment on each slave in which tasks will be
# run). Provide a setup function for this context. 
c = m.get_context(parallel_tasks.slave_setup)

num_global_beads = parallel_tasks.master_setup()

# Add tasks with different input parameters
for b in range(num_global_beads):
    c.add_task(parallel_tasks.SlaveTask(b))

# Run all tasks, distributed between the slaves. Get the results in
# the order they are returned (not necessarily the order they were created).
for x in c.get_results_unordered():
    print(x)
