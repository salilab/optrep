import sys
sys.dont_write_bytecode=True

import IMP
import IMP.optrep
import os,string,math
import argparse
import subprocess
import time
import stats_helper

def parse_args():
    
    parser = argparse.ArgumentParser(description="Get the optimal representation by constructing non-uniform representations using incremental coarse-graining. Usage: incremental_coarse_grain.py -s 2IDO -e 1 -r 0. Flag -h for more details.")
    
    parser.add_argument("-s","--system",dest="system",type=str,help="the name of the system")
    
    parser.add_argument("-e","--experiment",dest="experiment",type=str,help="the experiment number")
    
    parser.add_argument("-r","--restart",dest="restart",type=str,help="restart option. e.g. 1-b 1-s 1-g 1-p or 5-b 5-s 5-g 5-p or 0 if start from the beginning. The first number is the coarse-grained bead size you need to start from. The second alphabet stands for what stage we want to restart from: b for bead map generation, s for sampling, g for good-scoring model selection and p for precision estimation.")
  
    result = parser.parse_args()
 
    return result

def create_sampling_qsub_script(bio_system,expt_number,resolution,run_number,imp_dir,sampling_script,cores_per_run,steps_per_run,topo_file,bead_map_file,move_sizes_file,xlinks_file,xlink_avg_distance,ligand_max_trans,ev_weight):

    qsf = open("job_sample.sh","w")
    print >>qsf,"#$ -S /bin/bash"
    print >>qsf,"#$ -N s"+bio_system+"_e"+expt_number+"_r"+resolution+"_rn"+str(run_number)
    print >>qsf,"#$ -o ./"
    print >>qsf,"#$ -e ./"
    print >>qsf,"#$ -r n"
    print >>qsf,"#$ -j n"
    print >>qsf,"#$ -l netappsali=5G"
    print >>qsf,"#$ -l h_rt=5:00:00"
    print >>qsf,"#$ -l arch=linux-x64"
    print >>qsf,"#$ -R yes"
    print >>qsf,"#$ -cwd"
    print >>qsf,"#$ -pe ompi "+cores_per_run
    print >>qsf,"#$ -l hostname='i*'"
    #print >>qsf,"#$ -q lab.q"
    print >>qsf,"module load openmpi-1.6-nodlopen"
    print >>qsf,"module load sali-libraries"
    print >>qsf, "hostname" 
    print >>qsf,"datetime=`date +%x_%T`"
    print >>qsf,"echo Datetime:$datetime"
    print >>qsf,"export PYTHONDONTWRITEBYTECODE=1" 
    print >>qsf,"mpirun -np "+cores_per_run+" "+imp_dir+"/build/setup_environment.sh python -B "+sampling_script+" "+bio_system+" "+topo_file+" "+bead_map_file+" "+move_sizes_file+" "+xlinks_file+" "+xlink_avg_distance+" "+steps_per_run+" "+ligand_max_trans+" "+ev_weight
        
    print >>qsf,"datetime=`date +%x_%T`"
    print >>qsf,"echo Datetime:$datetime"
    
    qsf.close()

def create_precision_qsub_script(bio_system,expt_number,resolution,imp_dir,num_cores,proteins_list,domains_list,topo_file,xscale,linear_cutoff,grid_size,output_prefix):
    
    qsf = open("job_precision.sh","w")
    print >>qsf,"#$ -S /bin/bash"
    print >>qsf,"#$ -N p"+bio_system+"_e"+expt_number+"_r"+resolution
    print >>qsf,"#$ -o ./"
    print >>qsf,"#$ -e ./"
    print >>qsf,"#$ -r n"
    print >>qsf,"#$ -j n"
    print >>qsf,"#$ -l netappsali=5G"
    print >>qsf,"#$ -l h_rt=2:00:00"
    print >>qsf,"#$ -l arch=linux-x64"
    print >>qsf,"#$ -R yes"
    print >>qsf,"#$ -cwd"
    print >>qsf,"#$ -t 1-"+num_cores
    print >>qsf,"#$ -l hostname='i*'"
    #print >>qsf,"#$ -q lab.q"
    print >>qsf,"module load openmpi-1.6-nodlopen"
    print >>qsf,"module load sali-libraries"
    
    # zero-based numbering
    print >>qsf,"j=$(( $SGE_TASK_ID - 1 ))"
   
    print >>qsf,"export PYTHONDONTWRITEBYTECODE=1" 
    print >>qsf,imp_dir+"/build/setup_environment.sh python -B "+imp_dir+"/imp/modules/optrep/pyext/src/estimate_sampling_precision.py -n "+num_cores+" -cn $j -pl "+proteins_list+" -dl "+domains_list+" -rd ./ -tf "+topo_file+" -gs "+grid_size+" -xs "+xscale+" -lc "+linear_cutoff+" -o "+output_prefix
        
    qsf.close()
    
def get_job_id_from_command_output(stdout):
    
    return stdout.strip().split()[2].split('.')[0]

def check_if_jobs_done(job_id_list):
    
    qstat_process = subprocess.Popen(["qstat"],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    qstat_out, qstat_err = qstat_process.communicate()
   
    running_jobs=[]
    for job_line in qstat_out.split('\n') :
        if len(job_line)>1 and not job_line.startswith('job') and not job_line.startswith('-'):
              running_jobs.append(job_line.split()[0])

    done=True
    
    for jid in job_id_list:
        if jid in running_jobs: # both are string! 
            done=False
    
    return done
   
            
def incremental_coarse_grain():
    # Authors: Shruthi Viswanath
     
    # process input
    arg=parse_args()
    
    config_file = os.path.join(os.path.expanduser('~'),"optrep","input",arg.system,arg.system+".config."+arg.experiment) 
    
    config_params = stats_helper.parse_config_file(config_file)
   
    
    sampling_run_prefix = "run."
    
    if arg.restart=="0":
        correctly_restarted = True
    else:
        starting_resolution = arg.restart.split('-')[0]
        starting_stage = arg.restart.split('-')[1]
        correctly_restarted = False # start from the correct previously left off stage
   
    # boolean to say if the optimization process is over
    all_done=False
    
    parent_dir = os.getcwd()
           
    for ires,resolution in enumerate(config_params["RESOLUTIONS_LIST"]): 
        
        curr_resolution_dir = "r"+resolution
        
        if not correctly_restarted and int(starting_resolution) > int(resolution):
            print "skipping resolution",resolution
            continue
                
        # STEP 0. make a directory for the new resolution
        if not os.path.exists(curr_resolution_dir):
            os.mkdir(curr_resolution_dir)
        
        os.chdir(curr_resolution_dir)
        print "Current resolution: ",resolution
            
        if (not correctly_restarted and starting_stage=="b") or correctly_restarted : 
            
            correctly_restarted=True
       
            #Step 1. Get the next level beadmap. 
            if ires==0: 
                # first resolution, create beadmap from topology file
                                    
                ret = subprocess.call([os.path.join(config_params["IMP_DIR"],"build","setup_environment.sh"),"python",os.path.join(config_params["IMP_DIR"],"imp/modules/optrep/pyext/src/set_next_beadmap.py"),"-u","create","-r",resolution,"-tf",os.path.join(config_params["INPUT_DIR"],config_params["TOPOLOGY_FILE"])])
                
            else:
                previous_resolution = config_params["RESOLUTIONS_LIST"][ires-1]
                
                previous_resolution_dir = "../r"+previous_resolution
                
                # subsequent resolutions, create from previous beadmap
                ret = subprocess.call([os.path.join(config_params["IMP_DIR"],"build","setup_environment.sh"),"python",os.path.join(config_params["IMP_DIR"],"imp/modules/optrep/pyext/src/set_next_beadmap.py"),"-u","update","-r",resolution,"-bmf",os.path.join(previous_resolution_dir,"bead_map_"+previous_resolution+".txt"),"-pf",os.path.join(previous_resolution_dir,"bead_precisions_"+previous_resolution+".txt")])
                                    
                if not ret==0:
                    print "Problems creating beadmap in resolution ",resolution," for system ",arg.system
                    exit(1)
          
        # Step 2. Do the sampling
        if (not correctly_restarted and starting_stage=="s") or correctly_restarted : 
            correctly_restarted=True
                        
            sampling_job_ids=[]
            
            # launch a qsub script for each sampling run
            for irun in range(1,int(config_params["NUM_SAMPLING_RUNS"])+1):
        
                # create the new sampling directory
                curr_sampling_dir = sampling_run_prefix+str(irun)
                os.mkdir(curr_sampling_dir)
                os.chdir(curr_sampling_dir)
                    
                # create the new sampling script
                create_sampling_qsub_script(arg.system,arg.experiment,resolution,irun,config_params["IMP_DIR"],config_params["SAMPLING_SCRIPT"],config_params["NUM_CORES_PER_SAMPLING_RUN"],config_params["NUM_STEPS_PER_SAMPLING_RUN"],os.path.join(config_params["INPUT_DIR"],config_params["TOPOLOGY_FILE"]),"../bead_map_"+resolution+".txt",os.path.join(config_params["INPUT_DIR"],config_params["MOVE_SIZES_FILE"]),os.path.join(config_params["INPUT_DIR"],config_params["XLINKS_FILE"]),config_params["XLINK_AVG_DISTANCE"],config_params["LIGAND_MAX_TRANS"],config_params["EV_WEIGHT"])
                
                # run qsub script 
                print "Launching sampling run ",irun
                sample_process = subprocess.Popen(["qsub","job_sample.sh"],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
                sample_out,sample_err = sample_process.communicate()
                            
                if sample_process.returncode!=0: 
                    print sample_err
                    print "Problems launching sampling run ",irun," in  resolution ",resolution," for system ",arg.system
                    exit(1)
            
                submitted_job_id = get_job_id_from_command_output(sample_out)  # get the job ID from the qsub command output
                sampling_job_ids.append(submitted_job_id) 
        
                os.chdir(os.path.join(parent_dir,curr_resolution_dir))
        
            # Step 2.5 Wait for sampling to be done
            sampling_done = check_if_jobs_done(sampling_job_ids)
            
            while not sampling_done:
                
                print "Going back to sleep for 5 mins"
                
                # sleep for 5 minutes
                time.sleep(60)
                
                # try your luck again
                sampling_done = check_if_jobs_done(sampling_job_ids)
             
        # Step 3. Get good-scoring models
        if (not correctly_restarted and starting_stage=="g") or correctly_restarted : 
                   
            correctly_restarted=True
            
            ret = subprocess.call([os.path.join(config_params["IMP_DIR"],"build","setup_environment.sh"),"python",os.path.join(config_params["IMP_DIR"],"imp/modules/optrep/pyext/src/select_good_scoring_models.py"),"-rd","./","-rp",sampling_run_prefix,"-cl",config_params["GOOD_SCORING_MODEL_CRITERIA_LIST"],"-kl",config_params["GOOD_SCORING_MODEL_KEYWORD_LIST"],"-agl",config_params["GOOD_SCORING_MODEL_AGGREGATE_LOWER_THRESHOLDS_LIST"],"-aul",config_params["GOOD_SCORING_MODEL_AGGREGATE_UPPER_THRESHOLDS_LIST"],"-mlt",config_params["GOOD_SCORING_MODEL_MEMBER_LOWER_THRESHOLDS_LIST"]                                                                                                                                                                                     ,"-mut",config_params["GOOD_SCORING_MODEL_MEMBER_UPPER_THRESHOLDS_LIST"]])
            
            if not ret==0:
                print "Problems getting good-scoring models in resolution ",resolution," for system ",arg.system
                exit(1)
            
        # Step 4. Get sampling precision of all beads
        if (not correctly_restarted and starting_stage=="p") or correctly_restarted : 
            
            correctly_restarted=True
                                   
            os.chdir(os.path.join("good_scoring_models"))
        
            create_precision_qsub_script(arg.system,arg.experiment,resolution,config_params["IMP_DIR"],config_params["NUM_CORES_ESTIMATE_PRECISION"],config_params["PROTEINS_TO_OPTIMIZE_LIST"],config_params["DOMAINS_TO_OPTIMIZE_LIST"],os.path.join(config_params["INPUT_DIR"],config_params["TOPOLOGY_FILE"]),config_params["XSCALE"],config_params["LINEAR_CUTOFF"][ires],config_params["GRID_SIZE"],"../bead_precisions_sub")
        
            # run qsub script 
            print "Launching precision run"
            precision_process = subprocess.Popen(["qsub","job_precision.sh"],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
            precision_out,precision_err = precision_process.communicate()
            
            if precision_process.returncode!=0: 
                print precision_err
                print "Problems launching precision run  in  resolution ",resolution," for system ",arg.system
                exit(1)
        
            precision_job_ids=[]
            submitted_job_id = get_job_id_from_command_output(precision_out)  # get the job ID from the qsub command output
            precision_job_ids.append(submitted_job_id) 
        
            # Step 4.5 Wait for the jobs to be done
            precision_done = check_if_jobs_done(precision_job_ids)
        
            while not precision_done:
                # sleep for 5 minutes
                time.sleep(60)
                
                precision_done = check_if_jobs_done(precision_job_ids)
                  
            # Step 5. Collate all precisions in one file
            os.chdir(os.path.join(parent_dir,curr_resolution_dir))
            
            bead_precisions_file = "bead_precisions_"+resolution+".txt"
            
            ret = subprocess.call([os.path.join(config_params["IMP_DIR"],"build","setup_environment.sh"),"python",os.path.join(config_params["IMP_DIR"],"imp/modules/optrep/pyext/src/collate_sampling_precisions.py"),"-n",config_params["NUM_CORES_ESTIMATE_PRECISION"],"-pl",config_params["PROTEINS_TO_OPTIMIZE_LIST"],"-dl",config_params["DOMAINS_TO_OPTIMIZE_LIST"],"-ip","./bead_precisions_sub","-o",bead_precisions_file])
                                    
            if not ret==0:
                print "Problems collating precisions in resolution ",resolution," for system ",arg.system
                exit(1)
        
        # Step 6. Check if done (all beads are precise)
        all_done = stats_helper.no_consecutive_beads_imprecise(bead_precisions_file)            
       
        if all_done:
            break
        
        os.chdir(parent_dir) # do a cd ../
      
 
if __name__ == "__main__" :
    incremental_coarse_grain()
