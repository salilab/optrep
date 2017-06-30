import IMP
import IMP.atom
import IMP.rmf
import IMP.pmi
import IMP.pmi.tools
import IMP.pmi.topology
import os,sys,string,math
import IMP.optrep.GoodScoringModelSelector

gsms=IMP.optrep.GoodScoringModelSelector.GoodScoringModelSelector("input/1SYX","1SYX_4.")

# test for score terms
#gsms.get_good_scoring_models(criteria_list=["Total_Score"],keywords_list=["Total_Score"],aggregate_lower_thresholds=[1.0],aggregate_upper_thresholds=[11.5],member_lower_thresholds=[0.0],member_upper_thresholds=[0.0])

# test for crosslink terms
gsms.get_good_scoring_models(criteria_list=["Crosslinks"],keywords_list=["CrossLinkingMassSpectrometryRestraint_Distance_"],aggregate_lower_thresholds=[0.9],aggregate_upper_thresholds=[1.0],member_lower_thresholds=[0.0],member_upper_thresholds=[12.0])

