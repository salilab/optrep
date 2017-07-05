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

SPE::SPE(std::string gsm_directory,std::vector<std::pair<std::string, std::string > > components_calculate_precision) {
        
  models_dir_=gsm_directory;
  components_calculate_precision_=components_calculate_precision;

  // Get the mapping of model index to sample number. 
  get_models_by_sample(gsm_directory+"model_sample_ids.txt");

}

/* Get the mapping of model index to sample number. 
*/
void SPE::get_models_by_sample(std::string sample_id_file) {

  std::ifstream sifile;
  sifile.open(sample_id_file.c_str());

  int model_index, sample_number;

  while (sifile >> model_index >> sample_number) {
      models_by_sample_[sample_number-1].push_back(model_index);
  } 
                
  total_number_of_models_ = models_by_sample_[0].size() + models_by_sample_[1].size();          
  sifile.close();

}

IMPOPTREP_END_NAMESPACE
