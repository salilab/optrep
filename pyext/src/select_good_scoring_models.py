
import IMP
import IMP.optrep
import IMP.optrep.BeadMapBuilder
import IMP.optrep.GoodScoringModelSelector
import IMP.optrep.SamplingPrecisionEstimator
import os,sys,string,math
import argparse

def parse_args():
    
    parser = argparse.ArgumentParser(description="List and extract good-scoring models from a set of sampling runs. Example of usage: select_good_scoring_models.py -rd <run_directory_for_sampling> -rp <run_prefix> -gs <grid_size>. Flag -h for more details.")
    
    parser.add_argument("-rd","--run_directory",dest="run_directory",help="directory in whcih sampling results are stored") 
    
    parser.add_argument("-rp","--run_prefix",dest="run_prefix",help="prefix of runs") 
                        
    parser.add_argument("-cl","--criteria_list",nargs='+',type=string,dest="criteria_list",help="list of criteria")
    parser.add_argument("-kl","--keywords_list",nargs='+',type=string,dest="keywords_list",help="list of stat file keywords")
    parser.add_argument("-agl","--aggregate_lower_thresholds",nargs='+',type=float,dest="aggregate_lower_threshold_list",help="aggregate lower threshold values")
    parser.add_argument("-aul","--aggregate_upper_thresholds",nargs='+',type=float,dest="aggregate_upper_threshold_list",help="aggregate upper threshold values")
    parser.add_argument("-mlt","--member_lower_thresholds",nargs='+',type=float,dest="member_lower_thresholds",help="member lower thresholds")
    parser.add_argument("-mut","--member_upper_thresholds",nargs='+',type=float,dest="member_upper_thresholds",help="member upper thresholds")

    result = parser.parse_args()
 
    return result

    
def select_good_scoring_models():
     # Authors: Shruthi Viswanath
     
    # process input
    args=parse_args()
    
    gsms=IMP.optrep.GoodScoringModelSelector.GoodScoringModelSelector(arg.run_dir,arg.run_prefix)
        
    gsms.get_good_scoring_models(criteria_list=arg.criteria_list,keywords_list=arg.keywords_list,aggregate_lower_thresholds=arg.aggregate_lower_thresholds,aggregate_upper_thresholds=arg.aggregate_upper_thresholds,member_lower_thresholds=arg.member_lower_thresholds,member_upper_thresholds=arg.member_upper_thresholds)
        
if __name__ == "__main__" :
    select_good_scoring_models()
