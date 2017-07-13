import IMP
import IMP.atom
import IMP.rmf
import IMP.pmi
import IMP.pmi.tools
import IMP.pmi.topology
import os,sys,string,math
import IMP.optrep.SamplingPrecisionEstimator

components_to_update=[("B","B_1")]

spe=IMP.optrep.SamplingPrecisionEstimator.SamplingPrecisionEstimator("input/1SYX/good_scoring_models",components_to_update)

spe.load_coordinates_and_bead_sizes_from_model_files()

spe.get_all_vs_all_distances(("B","B_1"),0);

#spe.estimate_perbead_sampling_precision(grid_size=2.0)

#spe.get_imprecise_beads(xscale=1.0)

#spe.print_bead_precisions("bead_precisions_python.dat")

