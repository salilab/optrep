/**
 *  \file SPE.cpp
 *  \brief Estimate particle-wise precision in order to coarse-grain beads
 *
 *  Copyright 2007-2017 IMP Inventors. All rights reserved.
 *
 */
#include <IMP/algebra.h>
#include <IMP/algebra/Vector3D.h>
#include <IMP/atom.h>
#include <IMP/core.h>
#include <iostream>
#include <string>
#include <IMP/rmf.h>
#include <IMP/optrep/SPE.h>

IMPOPTREP_BEGIN_NAMESPACE

SPE::SPE(String gsm_directory,std::vector<std::pair<String, String > > components_calculate_precision) {
        
  models_dir_=gsm_directory;
  components_calculate_precision_=components_calculate_precision;

  // Get the mapping of model index to sample number. 
  get_models_by_sample(gsm_directory+"model_sample_ids.txt");

}

/* Get the mapping of model index to sample number. 
*/
void SPE::get_models_by_sample(String sample_id_file) {

  std::ifstream sifile;
  sifile.open(sample_id_file.c_str());

  int model_index, sample_number;

  while (sifile >> model_index >> sample_number) {
      models_by_sample_[sample_number-1].push_back(model_index);
  } 
                
  total_number_of_models_ = models_by_sample_[0].size() + models_by_sample_[1].size();          
  sifile.close();

}

// String _included_protein_domain(String chain_full_name,protein_domain_list):
//         ''' Check if the current chain should be included for calculating bead precisions, based on the list of protein domains mentioned. 
//         '''
//         for protein_domain in protein_domain_list:
//             
//             protein=protein_domain[0]
//             if protein in chain_full_name:
//                 return protein_domain
// 
//         return False
// 
//     def load_coordinates_and_bead_sizes_from_model_files(self):
//         ''' Load all the coordinates from all good scoring models.
//         Store them by bead type so that it is easy to calculate RMSD/precision.
//         TODO does not do alignment yet
//         TODO does not handle multiscale properly yet? 
//         TODO Take care of ambiguity.
//         TODO Take care of multi-state protein.
//         TODO  There may be memory problem for large assemblies and large number of models. Need to refactor. 
//         TODO Make the code faster by integrating with the GPU based RMSD calculator code.
//         '''
// 
//         for i in range(self.total_num_models): # the models in the good_scoring_model dir are called 0.rmf3, 1.rmf3 and so on.. 
//             m=IMP.Model()
//             mdl_i=os.path.join(self.models_dir,str(i)+'.rmf3')
//     
//             fh_i=RMF.open_rmf_file_read_only(mdl_i)
//     
//             hier_i=IMP.rmf.create_hierarchies(fh_i,m)[0]
// 
//             IMP.rmf.load_frame(fh_i, 0)
// 
//             for state_i in hier_i.get_children():
//                 for chain_i in state_i.get_children():
//                     protein_domain_key = self._included_protein_domain(chain_i.get_name(),self.components_calculate_precision)
//                     bead_index=0
// 
//                     if not protein_domain_key:
//                         continue
//                     
//                     for a in IMP.atom.get_leaves(chain_i):
//                         curr_coords=IMP.core.XYZR(a).get_coordinates()
//                         curr_dia=IMP.core.XYZR(a).get_radius()*2.0
// 
//                         if i==0: #first model, need to initialize the coordinate and radius lists for each bead index
//                             self.bead_coords[protein_domain_key].append([curr_coords])
//                             self.bead_diameter[protein_domain_key].append(curr_dia) # assuming it is same for all models. In future versions this could be different for each model
//                         else:
//                             self.bead_coords[protein_domain_key][bead_index].append(curr_coords)
//                 
//                         # in any case increment bead index
//                         bead_index += 1
// 
//             del m,hier_i

IMPOPTREP_END_NAMESPACE
