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
class DistanceMatrix;

class Cluster;

class ChiSquareTestResult;

class IMPOPTREPEXPORT SPE {
 public:
 /**
	\param[in] gsm_directory the directory containing good-scoring models. Should also contain the sample id file, the file containing the sample identity (1 or 2) for each good-scoring model RMF (model_sample_ids.txt)
    \param[in] components_calculate_precision is a list of (protein,domain) elements whose precision needs to be calculated (i.e. whose representation needs to be improved).
    \param[in] topology_file the same file that was given for sampling (to know the components in the system)
  */

SPE(const String topology_file, const String gsm_directory,const std::vector<std::pair<String, String> > input_components_calculate_precision);

void load_coordinates_and_bead_sizes_from_model_files(bool break_after_first_model=false); 

/* this was for mono (unparallellized) code. Runningall ona a single core */ 
void estimate_all_beads_sampling_precision(const Float grid_size=1.0);

void get_all_imprecise_beads(const Float xscale);:

void print_all_bead_precisions(const std::string out_file_name) const ;

/*end API for mono code*/

/* this is for parallellized code. 1 core per bead. */
String estimate_and_print_single_bead_precision(const unsigned int global_bead_index,const Float grid_size=1.0, const  Float xscale) const;

/* designed return value as a string, as python is capable of parsing it. It could have been an object or an IMP_NAED_TUPLE_3 though */

/*end API for parallel code*/

//TODO make protected after testing
protected:
/* directory where model RMFs and the model sample IDs are stored. */
String models_dir_;

/* list of <protein,domain> names for which representation needs to be optimized*/
std::vector<std::pair<String, String > > components_calculate_precision_;

 std::size_t total_number_of_models_;

 /* required for manager in parallel environment to assign tasks. 
  * Also used by some of the mono methods (non-parallel) */ 
 unsigned int number_of_global_beads_ = 0; 

 std::vector<std::pair<unsigned int,unsigned int > > global_index_to_protein_and_local_bead_index_map_; 

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
void order_components_by_topology_file(const std::vector<std::pair<String, String > > input_components_calculate_precision, const String topology_file);

void get_models_by_sample(const String sample_id_file);

int included_protein_domain(const String chain_full_name) const;

IMP::optrep::DistanceMatrix* get_all_vs_all_distances(const unsigned int  global_bead_index) const;

Float get_sampling_precision(const Floats& cutoffs,const Floats& pvals,const Floats& cramersv,const Floats& populations) const;

IMP::Vector<IMP::optrep::Cluster*> precision_cluster(const Floats& distmat,const Float rmsd_cutoff) const;

IMP::algebra::Vector2Ds get_contingency_table(const IMP::Vector<IMP::optrep::Cluster*> cluster_result) const;

Float percent_ensemble_explained(const IMP::algebra::Vector2Ds& ctable) const;

IMP::optrep::ChiSquareTestResult* test_sampling_exhaustiveness(const IMP::algebra::Vector2Ds& observed_contingency_table) const ;

Float estimate_single_bead_precision(const unsigned int global_bead_index,const Float grid_size) const ;

bool is_commensurate(const Float bead_diameter,const Float bead_precision,const Float xscale) const ;


//IMP_OBJECT_METHODS(SPE);


};

IMPOPTREP_END_NAMESPACE

#endif /* IMPOPTREP_SPE_H */
