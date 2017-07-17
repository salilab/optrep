
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

 SPE(String topology_file, String gsm_directory,std::vector<std::pair<String, String> > input_components_calculate_precision);

void load_coordinates_and_bead_sizes_from_model_files(); 

//TODO make protected after testing
protected:
/* directory where model RMFs and the model sample IDs are stored. */
String models_dir_;

/* list of <protein,domain> names for which representation needs to be optimized*/
std::vector<std::pair<String, String > > components_calculate_precision_;

 std::size_t total_number_of_models_;
 
 std::vector<unsigned int> beads_per_protein_domain_;

 /* model IDs for models in each sample. array of 2 vectors, one per sample. 
 */
 Ints models_by_sample_[2];

 /*  Store the coordinates of each model according to the global bead index (protein_domain_index , bead_index for protein domain).
  Each vector corresponds to a bead and is the list of all coordinates of good-scoring model coordinates corresponding to that bead.
 */
  std::vector<std::vector<IMP::algebra::Vector3D> > bead_coords_;  
  
   /* Note the difference between the rest of the python-implemented classes (e.g. BeadMapBuilder) and the below: the implementation for each data structure here uses a global_bead_index  for the system (protein_domain_index , bead_index for protein domain) as opposed to the python classes which use a local_bead_index (index specific to a protein and domain). 
  */

 /* vector of diameter values for each bead */
 Floats bead_diameter_;
 
  /* vector of precision values for each bead */ 
 Floats bead_precisions_;

 /* vector showing which bead is imprecise and needs to be CG'ed */
 std::vector<bool > bead_imprecise_;
 
 // Methods
void order_components_by_topology_file(std::vector<std::pair<String, String > > input_components_calculate_precision, String topology_file);

void get_models_by_sample(String sample_id_file);

int included_protein_domain(String chain_full_name);

IMP::optrep::DistanceMatrix get_all_vs_all_distances(unsigned int global_bead_index);

/*
Float get_sampling_precision(Floats cutoffs,Floats pvals,Floats cramersv,Floats populations);

IMP::optrep::Clusters precision_cluster(Floats distmat,Float rmsd_cutoff);

IMP::algebra::Vector2Ds get_contingency_table(IMP::optrep::Clusters cluster_result);

Float percent_ensemble_explained(IMP::algebra::Vector2Ds ctable);

IMP::optrep::ChiSquareTestResult test_sampling_exhaustiveness(IMP::algebra::Vector2Ds observed_contingency_table);

Float estimate_single_bead_precision(unsigned int global_bead_index,Float grid_size);

bool is_commensurate(Float bead_diameter,Float bead_precision,Float xscale);

void get_imprecise_beads(Float xscale);
	
void print_bead_precisions(std::string out_file_name);

void estimate_perbead_sampling_precision(Float grid_size=1.0);

*/
//IMP_OBJECT_METHODS(SPE);


};

IMPOPTREP_END_NAMESPACE

#endif /* IMPOPTREP_SPE_H */
