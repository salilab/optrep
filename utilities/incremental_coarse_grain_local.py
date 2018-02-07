
import IMP
import IMP.optrep
import os,sys,string,math
import argparse
import subprocess
import datetime
import time
import scipy.stats
import stats_helper

def parse_args():
    
    parser = argparse.ArgumentParser(description="Get the optimal representation by constructing non-uniform representations using incremental coarse-graining. Usage: incremental_coarse_grain.py -s 2IDO -e 1 -r 0. Flag -h for more details.")
    
    parser.add_argument("-s","--system",dest="system",type=str,help="the name of the system")
    
    parser.add_argument("-e","--experiment",dest="experiment",type=str,help="the experiment number")
    
    parser.add_argument("-r","--restart",dest="restart",type=str,help="restart option. e.g. 1-b 1-s 1-g 1-p or 5-b 5-s 5-g 5-p or 0 if start from the beginning. The first number is the coarse-grained bead size you need to start from. The second alphabet stands for what stage we want to restart from: b for bead map generation, s for sampling, g for good-scoring model selection and p for precision estimation.")
    
    result = parser.parse_args()
 
    return result
                
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
    
    # always assume we start from the dir of exptNUM/TGT
    
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
            
            sampling_time_file=open("average_sampling_time_res"+resolution+"_expt"+arg.experiment+".txt",'w')
            
            sampling_times_in_seconds=[]
            
            # launch a qsub script for each sampling run
            for irun in range(1,int(config_params["NUM_SAMPLING_RUNS"])+1):
        
                # create the new sampling directory
                curr_sampling_dir = sampling_run_prefix+str(irun)
                os.mkdir(curr_sampling_dir)
                os.chdir(curr_sampling_dir)
            
                sampling_success=False

                while not sampling_success: # do this till a successful sampling run

                   time_start = datetime.datetime.now()
                    
                   ret = subprocess.call(["mpirun","-np",config_params["NUM_CORES_PER_SAMPLING_RUN"],os.path.join(config_params["IMP_DIR"],"build","setup_environment.sh"),"python",os.path.join(config_params["IMP_DIR"],config_params["SAMPLING_SCRIPT"]),arg.system,
                os.path.join(config_params["INPUT_DIR"],config_params["TOPOLOGY_FILE"]),"../bead_map_"+resolution+".txt",os.path.join(config_params["INPUT_DIR"],config_params["MOVE_SIZES_FILE"]),os.path.join(config_params["INPUT_DIR"],config_params["XLINKS_FILE"])
                , config_params["XLINK_AVG_DISTANCE"],config_params["NUM_STEPS_PER_SAMPLING_RUN"],config_params["LIGAND_MAX_TRANS"],config_params["EV_WEIGHT"]])
                    
                   num_lines_sampling_file = int(subprocess.check_output(['wc','-l','output/stat_replica.0.out']).strip().split()[0]) 

                   if num_lines_sampling_file == int(config_params["NUM_STEPS_PER_SAMPLING_RUN"])+1:
                        sampling_success=True
             
                   time_end = datetime.datetime.now()  #.strftime("%m/%d/%Y_%H:%M:%S")
                
                time_diff=time_end-time_start
                
                #time_start = datetime.datetime.now()
                
                #ret = subprocess.call(["mpirun","-np",config_params["NUM_CORES_PER_SAMPLING_RUN"],os.path.join(config_params["IMP_DIR"],"build","setup_environment.sh"),"python",os.path.join(config_params["IMP_DIR"],config_params["SAMPLING_SCRIPT"]),arg.system,os.path.join(config_params["INPUT_DIR"],config_params["TOPOLOGY_FILE"]),"../bead_map_"+resolution+".txt",os.path.join(config_params["INPUT_DIR"],config_params["MOVE_SIZES_FILE"]),os.path.join(config_params["INPUT_DIR"],config_params["XLINKS_FILE"]), config_params["XLINK_AVG_DISTANCE"],config_params["NUM_STEPS_PER_SAMPLING_RUN"]])
                
                #time_end = datetime.datetime.now()  #.strftime("%m/%d/%Y_%H:%M:%S")
                
                time_diff=time_end-time_start
                
                sampling_times_in_seconds.append(time_diff.total_seconds())
                
                os.chdir(os.path.join(parent_dir,curr_resolution_dir))
        
            avg_sampling_time = sum(sampling_times_in_seconds)/float(len(sampling_times_in_seconds))
            std_err_sampling_time = scipy.stats.sem(sampling_times_in_seconds)
            print >>sampling_time_file,"%.2f %.2f" %(avg_sampling_time, std_err_sampling_time)
            sampling_time_file.close()
            
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
        
            bead_precisions_file = "bead_precisions_"+resolution+".txt"
            
            proteins_list = " ".join(config_params["PROTEINS_TO_OPTIMIZE_LIST"]) 

            domains_list = " ".join(config_params["DOMAINS_TO_OPTIMIZE_LIST"])

            ret = subprocess.call([os.path.join(config_params["IMP_DIR"],"build","setup_environment.sh"),"python",os.path.join(config_params["IMP_DIR"],
            "imp/modules/optrep/pyext/src/estimate_sampling_precision_imp_parallel.py"),"-n",config_params["NUM_CORES_ESTIMATE_PRECISION"],"-pl",proteins_list,"-dl",domains_list,
            "-rd","./","-tf",os.path.join(config_params["INPUT_DIR"],config_params["TOPOLOGY_FILE"]),"-gs",config_params["GRID_SIZE"],"-xs",config_params["XSCALE"],"-lc",config_params["LINEAR_CUTOFF"][ires],"-o","../"+bead_precisions_file])
                
            if not ret==0:
                print "Problems getting bead-precisions in resolution ",resolution," for system ",arg.system
                exit(1)

            # Step 5. Collate all precisions in one file
            os.chdir(os.path.join(parent_dir,curr_resolution_dir))
            
            # wait till this is done before going to the next step
            while not os.path.exists(bead_precisions_file):
                time.sleep(60) # sleep for 60 seconds and try again         

        # Step 6. Check if done (all beads are precise)
        all_done = stats_helper.no_consecutive_beads_imprecise(bead_precisions_file)           
       
        if all_done:
            break
        
        os.chdir(parent_dir) # do a cd ../../
      
if __name__ == "__main__" :
    incremental_coarse_grain()
