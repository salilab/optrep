#!/bin/bash

expt=$1
tgt=$2
rep=$3

cd ~/optrep/expts/expt$expt/$tgt/r$rep

numRuns=`ls -d run.* | wc -l`

numReplicas=4

numReplicasMinusOne=3

rm good_scoring_models/all_model_scores.out 

for i in `seq 1 $numRuns`
do 
	for j in `seq 0 $numReplicasMinusOne`
	do

	~/imp-clean/build/setup_environment.sh python ~/imp-clean/imp/modules/pmi/pyext/src/process_output.py -f run.$i/output/stat.$j.out -s rmf_frame_index Total_Score CrossLinkingMassSpectrometryRestraint_Data_Score_XLDSS CrossLinkingMassSpectrometryRestraint_PriorPsi_Score_XLDSS CrossLinkingMassSpectrometryRestraint_PriorSig_Score_XLDSS | grep -v \# | awk -v runid="$i" -v repid="$j" '{print runid,repid,$2,$3,$4+$5+$6,$4}' >> good_scoring_models/all_model_scores.out
	done
done

