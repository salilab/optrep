import IMP
import IMP.atom
import IMP.rmf
from subprocess import Popen
import os,sys,string,math
import shutil
 	
class GoodScoringModelSelector(object):
    # Authors: Shruthi Viswanath
    
    ''' Select good-scoring models based on scores and/or data satisfaction.
    Exrtact the corresponding RMFs and put them in a separate directory
    '''
    
    def __init__(self,run_directory,run_prefix):
        """Constructor.
        @param run_directory the directory containing subdirectories of runs
        @param run_prefix the prefix for each run directory. For e.g. if the subdirectories are modeling_run1, modeling_run2, etc. the prefix is modeling_run
        """
        self.run_dir=run_directory
        self.run_prefix=run_prefix
        self.all_good_scoring_models_indices=[]# list with each member as a tuple (run id,replica id,frame id) corresponding to good-scoring models

    def _get_subfields_for_criteria(self,field_headers,keywords_list):
        ''' Given the list of keywords, get all the stat file entries corresponding to each keyword.'''

        fields_for_criteria=[[] for j in range(len(keywords_list))] # list of lists corresponding to field indices for each keyword

        for ki,kw in enumerate(keywords_list):
            for j in range(len(field_headers)):
                if kw in field_headers[j]:
                    fields_for_criteria[ki].append(j)

        return fields_for_criteria

    def get_list_of_good_scoring_models(self,criteria_list=[],keywords_list=[],aggregate_lower_thresholds=[],
                                        aggregate_upper_thresholds=[],member_lower_thresholds=[],member_upper_thresholds=[]): 
        ''' Loops over all stat files in the run directory and populates the 
        list of good-scoring models.
        @param criteria_list is a list of what scores/data are checked to define good-scoring models. For e.g., for checking crosslinks, excluded volume and total score, one would say something like ["crosslinks","EVscore","Total_score"]. Note that all score terms should contain the keyword "score" (case insensitive)
        @param keywords_list is the list of keywords in the PMI stat file that need to be checked for each datatype/score in the criteria list
        @param aggregate_lower_thresholds The list of lower bounds on the values corresponding to fields in the criteria_list. Aggregates are used for terms like % of crosslink satisfaction and thresholds of score terms 
        @param aggregate_upper_thresholds The list of upper bounds on the values corresponding to fields in the criteria_list. Aggregates are used for terms like % of crosslink satisfaction and thresholds of score terms
        @param member_lower_thresholds The list of lower bounds for values of subcomponents of an aggregate term. E.g. for crosslink satisfaction the thresholds are on distances for each individual crosslink. For score terms this can be ignored since thresholds are mentioned in the aggregate fields.
        @param member_upper_thresholds The list of upper bounds for values of subcomponents of an aggregate term. E.g. for crosslink satisfaction the thresholds are on distances for each individual crosslink. For score terms this can be ignored since thresholds are mentioned in the aggregate fields.
        '''

        output_dir=os.path.join(self.run_dir,"good_scoring_models")
        shutil.rmtree(output_dir)
        os.mkdir(output_dir)
        
        outf=open(os.path.join(output_dir,"model_ids_scores.txt"),'w')
      
        for each_run_dir in sorted(os.listdir(self.run_dir)):  
            runid=each_run_dir.split(self.run_prefix)[1]
            
            for each_replica_stat_file in sorted(glob.glob(os.path.join(each_run_dir,"output")+"/stat.*.out")):
                    
                    replicaid=each_replica_stat_file.strip(".out").split(".")[-1]
                    
                    rsf=open(each_replica_stat_file,'r')
                    
                    for linecount,each_model_line in enumerate(sf.readlines()): # for each model in the current replica
                        
                        if linecount==0:
                            field_headers=eval(each_model_line.strip())
                            fields_for_criteria=_get_subfields_for_criteria(field_headers,keywords_list)
                           
                            continue
                        
                        frameid=i-1
                        
                        dat=eval(each_model_line.strip())
                        
                        model_satisfies=False
                        model_criteria_values=[]
                                                
                        for si,score_type in enumerate(criteria_list):
                            if "crosslink" is in score_type.lower():
                                crosslink_distance_values=[dat[j] for j in fields_for_criteria[si]] # TODO : consider ambiguity
                      
                                satisfied_percent,model_satisfies=_get_crosslink_satisfaction(crosslink_distance_values,aggregate_lower_thresholds[si],aggregate_upper_thresholds[si],lower_thresholds[si],upper_thresholds[si])
                                
                                model_criteria_values.append(satisfied_percent)

                            elif "score" is in score_type.lower():
                                score_value=dat[fields_for_criteria[si]]
                                model_satisfies=_get_score_satisfaction(score_value,aggregate_lower_thresholds[si],aggregate_upper_thresholds[si])
                                model_criteria_values.append(score_value)
                      
                            if not model_satisfies:
                                break

                        if model_satisfies:
                            self.all_good_scoring_models_indices.append(runid,replicaid,frameid)
                            print >>outf,len(self.all_good_Scoring_models_indices)-1,runid,replicaid,frameid,
                            for mcv in model_criteria_values:
                                print >>outf,mcv,
                            print >>outf
        
                    rsf.close()

        _extract_models_from_trajectories(output_dir) 

    def _get_crosslink_satisfaction(self,crosslink_distance_values,crosslink_percentage_lower_threshold,crosslink_percentage_upper_threshold,xlink_distance_lower_threshold,xlink_distance_upper_threshold):
        ''' For crosslinks, we want models with atleast x% (e.g. 90%) or more crosslink satisfaction. A crosslink is satisfied if the distance is between the lower and upper distance thresholds 
        @param crosslink_distance_values values of distances in the current model
        @param crosslink_percentage_lower_threshold atleast x% of crosslinks should be within the below distance thresholds
        @param crosslink_percentage_upper_threshold atmost x% of crosslinks should be within the below distance thresholds (usually 100%: dummy parameter)
        @param xlink_distance_lower_threshold a crosslink should be atleast this distance apart (usually 0 Å) to be considered satisfied
        @param xlink_distance_upper_threshold a crosslink should be atmost this distance apart to be considered satisfied
        
        '''
        satisfied_xlinks=0.0
        for d in crosslink_distance_values:
            if d>=xlink_distance_lower_threshold and d<=xlink_distance_upper_threshold:
                satisfied_xlinks+=1.0

        percent_satisfied=satisfied_xlinks/float(len(crosslink_distance_values))

        if percent_satisfied>=crosslink_percentage_lower_threshold and percent_satisfied<=crosslink_percentage_upper_threshold:
            return percent_satisfied,True
        else:
            return percent_satisfied,False

    
    def _get_score_satisfaction(self,score,lower_threshold,upper_threshold):
        ''' Check if the score is within the thresholds
        '''
           
        if score<=upper_threshold and score>=lower_threshold:
            return True
        return False

    def _extract_models_from_trajectories(self,output_dir) 
        ''' Given the list of all good-scoring model indices, extract their frames and store them ordered by the list index.
        '''
        
        for i,gsm in enumerate(self.all_good_scoring_models_indices):
            (runid,replicaid,frameid)=gsm 
            
            trajfile=os.path.join(os.environ['IMP_BIN_DIR'],self.run_dir,self.run_prefix+runid,'output','rmfs',replicaid+'.rmf3')
         
            slice_location=os.path.join(os.environ['IMP_BIN_DIR'],'build','bin','rmf_slice')
                                        
            rmf_slice=Popen(slice_location,trajfile,"-f",frameid,os.path.join(outputdir,str(i)+'.rmf3'))

        
    def _split_good_scoring_models_into_two_subsets(split_type="divide_by_run_ids"):
        ''' Get the listof good scoring models and split them into two samples. Write the model id and sample id onto a file. 
        @param split_type how to split good scoring models into two samples. Current options are:
        (a) divide_by_run_ids : where the list of runids is divided into 2. e.g. if we have runs from 1-50, good scoring models from runs
        1-25 is sample1 and those from runs 26-50 is sample 2. 
        (b) random : split the set of good scoring models into two subsets at random.
        '''
                
        if split_type=="divide_by_run_ids":



        elif split_type=="random":




        # write model and sample IDs to a file

