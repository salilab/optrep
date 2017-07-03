/**
 *  \file SamplingPrecisionEstimator.cpp
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

IMPOPTREP_BEGIN_NAMESPACE

SamplingPrecisionEstimator::SamplingPrecisionEstimator(std::string gsm_directory,std::vector<std::pair<std::string, std::string>> components_calculate_precision) {
        
        models_dir_=gsm_directory;
	components_calculate_precision_=components_calculate_precision;

        // Initialize bead precision and diameter dictionaries.
        self.bead_precisions={} #dictionary with key as protein domain and bead precision for each primitive in the protein domain
        self.bead_diameter={} #dictionary with key as protein domain and size of each primitive in the protein domain

        self.bead_imprecise={} # dictionary with key as protein domain and values as booleans that say which primitive is imprecisely sampled

        for protein_domain_key in components_calculate_precision:
            self.bead_precisions[protein_domain_key]=[]
            self.bead_diameter[protein_domain_key]=[]
            self.bead_imprecise[protein_domain_key]=[]
            self.component_coords[protein_domain_key]=[]


        // Get the mapping of model index to sample number. 
        get_models_by_sample(sample_id_file=strcat(gsm_directory,"model_sample_ids.txt"));

}

/* dummy function
*/
void SamplingPrecisionEstimator::show(std::ostream &out) const {
  out << "sampling precision estimator\n";
}

IMPOPTREP_END_NAMESPACE
