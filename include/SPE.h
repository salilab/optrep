
/**
 *  \file SPE.h
 *  \brief Get per-bead sampling precision
 *
 *  Copyright 2007-2017 IMP Inventors. All rights reserved.
 *
 */

#ifndef IMPOPTREP_SPE_H
#define IMPOPTREP_SPE_H
#include <IMP.h>
#include <IMP/atom.h>
#include <IMP/core.h>
#include "optrep_config.h"

IMPOPTREP_BEGIN_NAMESPACE

//! Obtain bead-wise sampling precisions.
/**
     Get the all vs all distances for every primitive (bead), get bead precisions, and whether a bead is imprecise based on that. 
 */
class IMPOPTREPEXPORT SPE {
 public:
 /**
	\param[in] gsm_directory the directory containing good-scoring models. Should also contain the sample id file, the file containing the sample identity (1 or 2) for each good-scoring model RMF (model_sample_ids.txt)
     	\param[in] components_calculate_precision is a list of (protein,domain) elements whose precision needs to be calculated (i.e. whose representation needs to be improved).
  */

    SPE(String gsm_directory,std::vector<std::pair<String, String> > components_calculate_precision);

 protected:
    /* directory where model RMFs and the model sample IDs are stored. */
    String models_dir_;

    /* list of <protein,domain> names for which representation needs to be optimized*/
    std::vector<std::pair<String, String > > components_calculate_precision_;

 unsigned int total_number_of_models_;

 /* model IDs for models in each sample. array of 2 vectors, one per sample. 
 */
 Ints models_by_sample_[2];

 /*  Store the coordinates of each model according to the bead index.
  Each vector corresponds to a bead and is the list of all coordinates of good-scoring models corresponding to that bead.
 */
  std::vector<std::vector<IMP::algebra::Vector3D> >bead_coords_;  

  /* Note the difference between the rest of the python-implemented classes (e.g. BeadMapBuilder) 
 and the below: the implementation for each data structure here uses a global_bead_index  for the system (includes all proteins and domains) as opposed to the python classes which use a local_bead_index (index specific to a protein and domain). 
  */

  /* vector of precision values for each bead */ 
 Floats bead_precisions_;

 /* vector of diameter values for each bead */
 Floats bead_diameter_;

 /* vector showing which bead is imprecise and needs to be CG'ed */
 std::vector<bool > bead_imprecise_;


 // Methods
 
 void get_models_by_sample(std::string sample_id_file);

  //IMP_OBJECT_METHODS(SPE);


};

IMPOPTREP_END_NAMESPACE

#endif /* IMPOPTREP_SPE_H */
