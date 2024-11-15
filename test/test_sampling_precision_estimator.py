import IMP.test
import IMP
import IMP.optrep.SamplingPrecisionEstimator

class Tests(IMP.test.TestCase):
    def test_spe(self):
        """Test sampling precision estimator"""
        components_to_update=[("B","B_1")]

        gsm = self.get_input_file_name('1SYX/good_scoring_models')
        spe=IMP.optrep.SamplingPrecisionEstimator.SamplingPrecisionEstimator(
            gsm, components_to_update)

        spe.load_coordinates_and_bead_sizes_from_model_files()

        #for i in range(10):
            #spe.get_all_vs_all_distances(("B","B_1"),i)
            #print i

        #spe.estimate_single_bead_precision(("B","B_1"),0,grid_size=2.0)
        spe.estimate_perbead_sampling_precision(grid_size=2.0)

        spe.get_imprecise_beads(xscale=1.0)

        spe.print_bead_precisions("bead_precisions_python.dat")


if __name__ == '__main__':
    IMP.test.main()
