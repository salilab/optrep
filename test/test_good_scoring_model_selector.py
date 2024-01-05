import IMP.test
import shutil
import os
import IMP.optrep.GoodScoringModelSelector


class Tests(IMP.test.TestCase):
    def test_gsm(self):
        """Test GoodScoringModelSelector"""
        # Don't trample all over gsmdir used by other tests
        with IMP.test.temporary_directory() as tmpdir:
            gsmdir = os.path.join(tmpdir, "1SYX")
            shutil.copytree(self.get_input_file_name("1SYX"), gsmdir)
            gsms=IMP.optrep.GoodScoringModelSelector.GoodScoringModelSelector(
                gsmdir, "1SYX_4.")

            # test for score terms
            gsms.get_good_scoring_models(criteria_list=["Total_Score"],keywords_list=["Total_Score"],aggregate_lower_thresholds=[1.0],aggregate_upper_thresholds=[11.5],member_lower_thresholds=[0.0],member_upper_thresholds=[0.0])

            # test for crosslink terms
            #gsms.get_good_scoring_models(criteria_list=["Crosslinks"],keywords_list=["CrossLinkingMassSpectrometryRestraint_Distance_"],aggregate_lower_thresholds=[0.9],aggregate_upper_thresholds=[1.0],member_lower_thresholds=[0.0],member_upper_thresholds=[12.0])

if __name__ == '__main__':
    IMP.test.main()
