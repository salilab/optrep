
import IMP
import IMP.atom
import IMP.rmf
import IMP.pmi
import IMP.pmi.tools
import IMP.pmi.topology
import os,sys,string,math
import IMP.optrep

components_to_update=IMP.optrep.ProteinDomainList([("B","B_1")]) 
# outer () is for argument, [] is for vector, inner () is for string pair.

spe=IMP.optrep.SPE("input/1SYX/good_scoring_models/",components_to_update)

a=spe._included_protein_domain("Chain_B_full")

print a

b=spe._included_protein_domain("Chain_A_full")

print b

#spe.load_coordinates_and_bead_sizes_from_model_files()

#spe.estimate_perbead_sampling_precision(grid_size=2.0)

#spe.get_imprecise_beads(xscale=1.0)

#spe.print_bead_precisions("bead_precisions_cpp.dat")

