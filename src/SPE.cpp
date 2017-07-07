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

number_of_protein_domains_ = components_calculate_precision_.size();

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

/* Check if the current chain should be included for calculating bead precisions, 
based on the list of protein domains mentioned. 
*/
int SPE::included_protein_domain_(String chain_full_name) {
            
    for(unsigned int i=0;i<components_calculate_precision_.size();i++) {
        
        String protein=components_calculate_precision_[i].first;
        
        std::size_t found=chain_full_name.find(protein);
        if (found != String::npos)
            return(i);
    }

    return(-1);
}

void spe::load_coordinates_and_bead_sizes_from_model_files() {
/* Load all the coordinates from all good scoring models.
        Store them by bead type so that it is easy to calculate RMSD/precision.
        TODO does not do alignment yet
        TODO does not handle multiscale properly yet? 
        TODO Take care of ambiguity.
        TODO Take care of multi-state protein.
        TODO  There may be memory problem for large assemblies and large number of models. Need to refactor. 
        TODO Make the code faster by integrating with the GPU based RMSD calculator code.
*/
    /* first create hierarchies from the first frame */
    IMP_NEW(Model, m, ());
    
    String mdl_0 = models_dir + std::to_str(0) + ".rmf3";
    
    RMF::FileConstHandle fh_0=RMF.open_rmf_file_read_only(mdl_0);
    
    IMP::atom::Hierarchies hier_0=IMP.rmf.create_hierarchies(fh_0,m); 
    
    fh_0 = RMF::FileHandle();
    
    
    for (unsigned int i=0;i<total_number_of_models_;i++) { 
    /* the models in the good_scoring_model dir are called 0.rmf3, 1.rmf3 and so on.. */
                
        String mdl_i = models_dir + std::to_str(i) + ".rmf3";
        
        RMF::FileConstHandle fh_i=RMF.open_rmf_file_read_only(mdl_i);
    
        rmf::link_hierarchies(fh_i,hier_0);   

        rmf::load_frame(fh_i, 0);                           //  load frame 0

        // need only from the first hierarchy
        IMP::atom::Hierarchies state_i = hier_0[0].get_children(); 
        for (unsigned int si = 0; si<state_i.size();si++) {
            IMP::atom::Hierarchies chain_i = state_i[si].get_children();
                for (unsigned int ci = 0; ci<chain_i.size();ci++) {
                    
                    protein_domain_index = included_protein_domain_(chain_i[ci]->get_name());
                    
                    if (protein_domain_index == -1) 
                        continue;
                    
                    bead_index=0;
                    
                    atom::Hierarchies beads_i = IMP::atom::get_leaves(chain_i[ci]);
                    for (unsigned int bi = 0; bi<beads_i.size();bi++) { 
                        
                        IMP::algebra::Vector3D curr_coords=IMP::core::XYZR(beads_i[bi]).get_coordinates();
                        Float curr_dia=IMP::core::XYZR(beads_i[bi]).get_radius()*2.0;

                        if (i==0) {                         //first model, need to initialize the coordinate and radius lists for each bead index
                            
                            
                            bead_coords_.push_back([])[protein_domain_index].append([curr_coords])
                            // assuming it is same for all models. In future versions this could be different for each model
                            bead_diameter_.push_back(curr_dia)
                        }
                        else {
                            bead_coords_[protein_domain_index][bead_index].append(curr_coords)
                        }
                
                        // in any case increment bead index
                        bead_index += 1;
                }
        }
        
} 


IMPOPTREP_END_NAMESPACE
