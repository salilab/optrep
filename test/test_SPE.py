from __future__ import print_function, division
import IMP.test
import IMP
import IMP.atom
import IMP.rmf
import IMP.pmi
import IMP.pmi.tools
import IMP.pmi.topology
import os,sys,string,math
import IMP.optrep

class Tests(IMP.test.TestCase):
    def test_spe(self):
        """Test SPE"""
        #components_to_update=IMP.optrep.ProteinDomainList([("B","B_1"),("A","A_1")])

        components_to_update=IMP.optrep.ProteinDomainList([("B","B_1")])

        # outer () is for argument, [] is for vector, inner () is for string pair.
        topology = self.get_input_file_name('1SYX/1SYX.topology.txt')
        gsm = self.get_input_file_name('1SYX/good_scoring_models')

        spe=IMP.optrep.SPE(topology, gsm, components_to_update)

        spe.load_coordinates_and_bead_sizes_from_model_files()

        #for i in range(10):
          #distanceMatrix = spe.get_all_vs_all_distances(i)
          #print i

        # spe.estimate_single_bead_precision(0,grid_size=2.0)

        spe.estimate_all_beads_sampling_precision(grid_size=2.0)

        spe.get_all_imprecise_beads(xscale=1.0, linear_cutoff=1.0)

        spe.print_all_bead_precisions("bead_precisions_cpp.dat")


if __name__ == '__main__':
    IMP.test.main()
