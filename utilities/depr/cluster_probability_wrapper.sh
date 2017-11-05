#!/bin/bash

exptnum=$1

bio=$2

rep=$3 #1, 20, 50,max

for score_type in total_score data_score
do
	for averaging_type in direct_average average_of_average
	do
		for models_to_average in all_gsms all_models retained_clusters
		do
			 ~/imp-clean/build/setup_environment.sh python ~/optrep/code/module/utilities/get_normalized_cluster_probabilities.py $exptnum $bio $rep ~/optrep/expts/ probability $score_type $averaging_type $models_to_average

		done


	done
done	
		


