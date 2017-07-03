/**
 *  \file SamplingPrecisionEstimator.h
 *  \brief Get per-bead sampling precision
 *
 *  Copyright 2007-2017 IMP Inventors. All rights reserved.
 *
 */

#ifndef IMPOPTREP_SAMPLING_PRECISION_ESTIMATOR_H
#define IMPOPTREP_SAMPLING_PRECISION_ESTIMATOR_H
#include <IMP.h>
#include <IMP/atom.h>
#include <IMP/core.h>

IMPOPTREP_BEGIN_NAMESPACE

//! Obtain bead-wise sampling precisions.
/**
     Get the all vs all distances for every primitive (bead), get bead precisions, and whether a bead is imprecise based on that. 
 */
class IMPOPTREPEXPORT SamplingPrecisionEstimator {
 public:
 /**
	\param[in] gsm_directory the directory containing good-scoring models. Should also contain the sample id file, the file containing the sample identity (1 or 2) for each good-scoring model RMF (model_sample_ids.txt)
     	\param[in] components_calculate_precision is a list of (protein,domain) elements whose precision needs to be calculated (i.e. whose representation needs to be improved).
  */

    SamplingPrecisionEstimator(std::string gsm_directory,std::vector<std::pair<std::string, std::string> > components_calculate_precision);

 protected:
 unsigned int number_of_good_scoring_models;

 // model IDs for models in each sample
 std::vector<int> models_of_sample_1_;
 std::vector<int> models_of_sample_2_;

 /*  Store the coordinates of each model according to the bead index.
  This dictionary is indexed by a protein domain. The elements are lists of lists. 
  Each list corresponds to a bead and is the list of all coordinates of good-scoring models corresponding to that bead.
 */
  std::vector<std::vector<IMP::algebra::Vector3D>>component_coords;  
   
  IMP_OBJECT_METHODS(SamplingPrecisionEstimator);
  IMP_SHOWABLE(SamplingPrecisionEstimator);


}

IMPOPTREP_END_NAMESPACE

#endif /* IMPOPTREP_SAMPLING_PRECISION_ESTIMATOR_H */
