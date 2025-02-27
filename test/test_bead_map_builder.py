import IMP.test
import IMP
import IMP.optrep.BeadMapBuilder

class Tests(IMP.test.TestCase):
    def test_bead_map_builder(self):
        """Test BeadMapBuilder"""
        # test 1 setting from topology file

        bm=IMP.optrep.BeadMapBuilder.BeadMapBuilder()

        topology = self.get_input_file_name("topology.in.txt")
        bm.set_bead_map_from_topology_file(topology)

        bm.write_bead_map_to_file("beadmap_r1.txt")

        #bm.show_bead_map()

        bm2=IMP.optrep.BeadMapBuilder.BeadMapBuilder()

        bm2.set_bead_map_from_beadmap_file("beadmap_r1.txt")

        #bm2.show_bead_map()


        # incrementally CG
        components_to_update_1=[("B","B_1")]

        imprecise_beads_file = self.get_input_file_name("imprecise_beads_B.txt")
        bm.set_imprecise_beads_from_file(imprecise_beads_file)

        updated = bm.update_all_bead_maps(5)

        if updated:
            bm.show_bead_map()

        #another CG
        components_to_update_2=[("A","A_1")]

        imprecise_beads_file = self.get_input_file_name("imprecise_beads_A.txt")

        bm2.set_imprecise_beads_from_file(imprecise_beads_file)

        updated = bm2.update_all_bead_maps(3)


        if updated:
            bm2.show_bead_map()


if __name__ == '__main__':
    IMP.test.main()
