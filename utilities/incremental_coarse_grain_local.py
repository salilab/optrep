
import IMP
import IMP.optrep
import os,sys,string,math
import argparse
import subprocess
import datetime
import time

def parse_args():
    
    parser = argparse.ArgumentParser(description="Get the optimal representation by constructing non-uniform representations using incremental coarse-graining. Usage: incremental_coarse_grain.py -s 2IDO -e 1. Flag -h for more details.")
    
    parser.add_argument("-s","--system",dest="system",type=str,help="the name of the system")
    
    parser.add_argument("-e","--experiment",dest="experiment",type=str,help="the experiment number")
  
    result = parser.parse_args()
 
    return result


def parse_config_file(config_file):
    
    cf=open(config_file,'r')
    
    configs_dict ={}
    
    for ln in cf.readlines():
    
        if ln.startswith("#") or not ln.strip(): # ignore comment and blank lines
            continue
        
        fields=ln.strip().split("=")
        
        if fields[0]=="RESOLUTIONS_LIST":
            configs_dict[fields[0]]=fields[1].split() # key to list map
            
        else :
            if fields[1].startswith('~'): # location of a file or directory
     
                fields[1]=os.path.expanduser('~')+fields[1].lstrip('~') #os.path.join did not work here for some weird reason!
      
            configs_dict[fields[0]]=fields[1] # just a single key to string map

    cf.close()    
    
    return configs_dict
    
    

def all_beads_precise(bead_precisions_file):
          
    bpf=open(bead_precisions_file,'r')
    
    for ln in bpf.readlines():
        imprecise=int(ln.strip().split()[4])
        if imprecise==1:
            return False
         
    return True
    
            
def incremental_coarse_grain():
    # Authors: Shruthi Viswanath
     
    # process input
    arg=parse_args()
    
    config_file = os.path.join(os.path.expanduser('~'),"optrep","input",arg.system,arg.system+".config."+arg.experiment) 
    
    config_params = parse_config_file(config_file)
    
    # boolean to say if the optimization process is over
    all_done=False
    
    parent_dir = os.getcwd()
           
    for ires,resolution in enumerate(config_params["RESOLUTIONS_LIST"]): 
        
        print "Current resolution: ",resolution
        # assume you are in the directory where you want to store the results
        # STEP 0. make a directory for the new resolution
        
        curr_resolution_dir = "r"+resolution
        os.mkdir(curr_resolution_dir)
        
        os.chdir(curr_resolution_dir)
       
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
        sampling_job_ids=[]
        sampling_run_prefix = "run."
        
        sampling_time_file=open("average_sampling_time_res"+resolution+"_expt"+arg.experiment+".txt",'w')
        
        sampling_times_in_seconds=[]
        
        # launch a qsub script for each sampling run
        for irun in range(1,int(config_params["NUM_SAMPLING_RUNS"])+1):
     
            # create the new sampling directory
            curr_sampling_dir = sampling_run_prefix+str(irun)
            os.mkdir(curr_sampling_dir)
            os.chdir(curr_sampling_dir)
                
            time_start = datetime.datetime.now()
            
            ret = subprocess.call(["mpirun","-np",config_params["NUM_CORES_PER_SAMPLING_RUN"],os.path.join(config_params["IMP_DIR"],"build","setup_environment.sh"),"python",os.path.join(config_params["IMP_DIR"],config_params["SAMPLING_SCRIPT"]),arg.system,
            os.path.join(config_params["INPUT_DIR"],config_params["TOPOLOGY_FILE"]),"../bead_map_"+resolution+".txt",os.path.join(config_params["INPUT_DIR"],config_params["MOVE_SIZES_FILE"]),os.path.join(config_params["INPUT_DIR"],config_params["XLINKS_FILE"])
            , config_params["XLINK_AVG_DISTANCE"],config_params["NUM_STEPS_PER_SAMPLING_RUN"]])
                    
            time_end = datetime.datetime.now()  #.strftime("%m/%d/%Y_%H:%M:%S")
            
            time_diff=time_end-time_start
            
            sampling_times_in_seconds.append(time_diff.total_seconds())
            
            os.chdir(os.path.join(parent_dir,curr_resolution_dir))
    
        avg_sampling_time = sum(sampling_times_in_seconds)/float(len(sampling_times_in_seconds))
        print >>sampling_time_file,avg_sampling_time
        sampling_time_file.close()
        
        # Step 3. Get good-scoring models
       
        ret = subprocess.call([os.path.join(config_params["IMP_DIR"],"build","setup_environment.sh"),"python",os.path.join(config_params["IMP_DIR"],"imp/modules/optrep/pyext/src/select_good_scoring_models.py"),"-rd","./","-rp",sampling_run_prefix,"-cl",config_params["GOOD_SCORING_MODEL_CRITERIA_LIST"],"-kl",config_params["GOOD_SCORING_MODEL_KEYWORD_LIST"],"-agl",config_params["GOOD_SCORING_MODEL_AGGREGATE_LOWER_THRESHOLDS_LIST"],"-aul",config_params["GOOD_SCORING_MODEL_AGGREGATE_UPPER_THRESHOLDS_LIST"],"-mlt",config_params["GOOD_SCORING_MODEL_MEMBER_LOWER_THRESHOLDS_LIST"]                                                                                                                                                                                     ,"-mut",config_params["GOOD_SCORING_MODEL_MEMBER_UPPER_THRESHOLDS_LIST"]])
        
        if not ret==0:
            print "Problems getting good-scoring models in resolution ",resolution," for system ",arg.system
            exit(1)
            
        # Step 4. Get sampling precision of all beads
        os.chdir("good_scoring_models")
     
        bead_precisions_file = "bead_precisions_"+resolution+".txt"
        
        ret = subprocess.call([os.path.join(config_params["IMP_DIR"],"build","setup_environment.sh"),"python",os.path.join(config_params["IMP_DIR"],
        "imp/modules/optrep/pyext/src/estimate_sampling_precision_imp_parallel.py"),"-n",config_params["NUM_CORES_ESTIMATE_PRECISION"],"-pl",config_params["PROTEINS_TO_OPTIMIZE_LIST"],"-dl",config_params["DOMAINS_TO_OPTIMIZE_LIST"],
        "-rd","./","-tf",os.path.join(config_params["INPUT_DIR"],config_params["TOPOLOGY_FILE"]),"-gs",config_params["GRID_SIZE"],"-xs",config_params["XSCALE"],"-o","../"+bead_precisions_file])
                
        # Step 5. Collate all precisions in one file
        os.chdir(os.path.join(parent_dir,curr_resolution_dir))
        
        # wait till this is done before going to the next step
        while not os.path.exists(bead_precisions_file):
            time.sleep(60) # sleep for 60 seconds and try again         

        # Step 6. Check if done (all beads are precise)
        all_done = all_beads_precise(bead_precisions_file)           
       
        if all_done:
            break
        
        os.chdir(parent_dir) # do a cd ../../
      
if __name__ == "__main__" :
    incremental_coarse_grain()
