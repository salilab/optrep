
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
#include <algorithm>
#include <cmath>
#include <boost/math/distributions/chi_squared.hpp>
#include <IMP/rmf.h>
#include <IMP/optrep/SPE.h>


IMPOPTREP_BEGIN_NAMESPACE

SPE::SPE(String topology_file, String gsm_directory,std::vector<std::pair<String, String > > input_components_calculate_precision) {
        
models_dir_=gsm_directory;

//components_calculate_precision_=components_calculate_precision;

order_components_by_topology_file(input_components_calculate_precision, topology_file);

// Get the mapping of model index to sample number. 
get_models_by_sample(gsm_directory+"model_sample_ids.txt");

}

void SPE::order_components_by_topology_file(std::vector<std::pair<String, String > > input_components_calculate_precision, String topology_file) {

	std::ifstream tfile;
	tfile.open(topology_file.c_str());
	
	//just need the first two fields
	String line, prot, dom;

	while (std::getline(tfile,line)) {
		std::istringstream ss(line);
		ss >> prot;
		ss >> dom; 

		std::pair<String,String > protein_domain = {prot,dom};

		for(unsigned int i=0;i<input_components_calculate_precision.size();i++){

			if (protein_domain == input_components_calculate_precision[i])	
                components_calculate_precision_.push_back(protein_domain); 
		}
        } 
	
	tfile.close();

	if(input_components_calculate_precision.size()!=components_calculate_precision_.size())
		std::cout<< "some components were not defined in the topology file!" <<std::endl; 

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
        
	// sort the model indices so that you can do binary search on it instead of doing a linear search each time
	for(unsigned int i=0;i<2;i++) {
		std::sort (models_by_sample_[i].begin(), models_by_sample_[i].end());
	}
		        
    total_number_of_models_ = models_by_sample_[0].size() + models_by_sample_[1].size();
    sifile.close();

}

/* Check if the current chain should be included for calculating bead precisions, 
based on the list of protein domains mentioned. 
*/
int SPE::included_protein_domain(String chain_full_name) {
            
    for(unsigned int i=0;i<components_calculate_precision_.size();i++) {
        
        String protein=components_calculate_precision_[i].first;
        
        std::size_t found=chain_full_name.find(protein);
        if (found != String::npos)
            return(i);
    }

    return(-1);
}

void SPE::load_coordinates_and_bead_sizes_from_model_files() {
   /* Load all the coordinates from all good scoring models.
        Store them by bead type so that it is easy to calculate RMSD/precision.
        TODO does not do alignment yet
        TODO does not handle multiscale properly yet? 
        TODO Take care of ambiguity.
        TODO Take care of multi-state protein.
        TODO  There may be memory problem for large assemblies and large number of models. Need to refactor. 
        TODO Make the code faster by integrating with the GPU based RMSD calculator code.
    */
    /* first create hierarchies from the first frame, so that it is faster to link/load frames.
    Assuming same hierarchy for all models of a given sampling run. */
    IMP_NEW(Model, m, ());
    
    String mdl_0 = models_dir_ + std::to_string(0) + ".rmf3";
    
    RMF::FileConstHandle fh_0=RMF::open_rmf_file_read_only(mdl_0);
    
    IMP::atom::Hierarchies hier_0=rmf::create_hierarchies(fh_0,m); 
    
    fh_0 = RMF::FileHandle();
    
    for (unsigned int i=0;i<total_number_of_models_;i++) { 
    /* the models in the good_scoring_model dir are called 0.rmf3, 1.rmf3 and so on.. */
                
        String mdl_i = models_dir_ + std::to_string(i) + ".rmf3";
        
        RMF::FileConstHandle fh_i=RMF::open_rmf_file_read_only(mdl_i);
    
        rmf::link_hierarchies(fh_i,hier_0);   

        rmf::load_frame(fh_i, RMF::FrameID(0));          

	unsigned int  global_bead_index=0;

        /* same for all models: list of (protein/domain,bead) indices for components to calculate precision.
        Can insert beads in this order since the components to calculate precision were sorted by hierarchy/topology file.
        */ 

        // need only from the first hierarchy 
        IMP::atom::Hierarchies states_i = hier_0[0].get_children(); 
        for (unsigned int si = 0; si<states_i.size();si++) {
            IMP::atom::Hierarchies chains_i = states_i[si].get_children();

            for (unsigned int ci = 0; ci<chains_i.size();ci++) {
                    
                    int protein_domain_index = included_protein_domain(chains_i[ci]->get_name());
                
                    if (protein_domain_index == -1) 
                        continue;
                    
                    if (i==0) { // assuming one chain per domain
                        beads_per_protein_domain_.push_back(0);
                    }
		                   
                    IMP::atom::Hierarchies beads_i = IMP::atom::get_leaves(chains_i[ci]);
                    for (unsigned int bi = 0; bi<beads_i.size();bi++) { 
                        
                        IMP::algebra::Vector3D curr_coords = IMP::core::XYZR(beads_i[bi]).get_coordinates();
                        
                        if(i==0) {
                            //create new coords list
                            bead_coords_.push_back( std::vector<IMP::algebra::Vector3D >() );
                            
                            beads_per_protein_domain_.back()+=1;
                                                
                            // assume same bead sizes for all models
                            Float curr_dia=IMP::core::XYZR(beads_i[bi]).get_radius()*2.0;
                            bead_diameter_.push_back(curr_dia);
                        }	
			
                        bead_coords_[global_bead_index].push_back(curr_coords);	

                        // std::cout << global_bead_index<<" " << bead_coords_[global_bead_index][i]<<" "<< curr_coords <<std::endl;

                        global_bead_index++; // we can do this because components to calculate precision are in order of the topology file

                    } //end for beads in protein domain     
			                
              } // end for chain

        } // end for state

 	}// end for each model
	
	
	//std::cout << beads_per_protein_domain_.size() <<" " <<   beads_per_protein_domain_[0] <<  " " <<  beads_per_protein_domain_[1] << std::endl; 
        	
} 

IMP::optrep::DistanceMatrix SPE::get_all_vs_all_distances(unsigned int  global_bead_index) {
        /* Return the distance matrix, minimum and maximum distance per bead.
        */
        IMP::optrep::DistanceMatrix d(total_number_of_models_);

        Float mindist=std::numeric_limits<double>::max();
        Float maxdist=0.0;

        for (unsigned int i = 0;i<total_number_of_models_-1;i++) {
                
            for (unsigned int j = i+1;j<total_number_of_models_;j++) {     
               
                Float dist=IMP::algebra::get_distance(bead_coords_[global_bead_index][i],bead_coords_[global_bead_index][j]);
                d.distmat.push_back(dist);

                if (dist<mindist) {
                    mindist=dist;
                }

                if (dist>maxdist) {
                    maxdist=dist;
                }
            }
        }

         d.mindist = mindist;
         d.maxdist = maxdist;
                    
        return d;

}

Float SPE::get_sampling_precision(Floats cutoffs,Floats pvals,Floats cramersv,Floats populations) {
        /* Given the 3 criteria for each cutoff, returns the sampling precision.
        This is the first cutoff at which all 3 criteria are satisfied.
        */
        Float sampling_precision=std::numeric_limits<double>::max();

        for (unsigned int i = 0;i<cutoffs.size(); i++) {
            if ((populations[i]>80.0) && ((pvals[i]>0.05) || (cramersv[i]<0.10))) {
                if (sampling_precision>cutoffs[i]) { 
                    sampling_precision=cutoffs[i];
                }
            }
            else {                                          // takes care of local minima? We want to choose the first time at which the above condition was  always satisfied
                sampling_precision=std::numeric_limits<double>::max();
            }

        }
        
        return (sampling_precision);

}

IMP::optrep::Clusters SPE::precision_cluster(Floats distmat,Float rmsd_cutoff) {
        /* Perform distance threshold-based clustering given the distance matrix and RMSD cutoff.
         Return the result as a list of clusters.
        */

        IMP::optrep::Clusters cluster_result;
		
        // Populate the neighbors of a given model
        IntsList neighbors;
        neighbors.insert(neighbors.end(),total_number_of_models_, Ints());
        
        for (unsigned int i = 0;i<total_number_of_models_;i++) {
            neighbors[i].push_back(i);            // model is a neighbor of itself 
        }
        
        unsigned int pair_index = 0;
        for (unsigned int i = 0;i<(total_number_of_models_-1);i++) {
            for(unsigned int j = i+1;j<total_number_of_models_;j++) {
                if (distmat[pair_index]<=rmsd_cutoff) {          // accepted to be a neighbor
                        neighbors[i].push_back(j);
                        neighbors[j].push_back(i);
                    }
                pair_index += 1;
            }   
                
        }
                
       // Get the cluster with the most neighbors, and iterate
       // Initially each model is unclustered
       Ints unclustered;
       std::vector<bool> boolUnclustered;
       for (unsigned int i = 0;i<(total_number_of_models_);i++) {
            unclustered.push_back(i);
            boolUnclustered.push_back(true); // 2 redundant pieces of information. Doesn't save time, just doing it for added code clarity
        }

       while(unclustered.size()>0) {
            // get cluster with maximum weight
             unsigned int max_neighbors=0;
             int curr_center=-1;
             for(Ints::iterator eachu=unclustered.begin();eachu!=unclustered.end();eachu++) {  
                 if (neighbors[*eachu].size()>max_neighbors) { // if multiple clusters have same maxweight this tie is broken arbitrarily! 
                     max_neighbors=neighbors[*eachu].size();
                     curr_center=*eachu;
		          }
             }

             //form a new cluster with u and its neighbors
               IMP::optrep::Cluster curr_cluster(curr_center, neighbors[curr_center]);
               cluster_result.push_back(curr_cluster);
                
             //update neighbors 
             for(Ints::iterator ne=neighbors[curr_center].begin();ne!=neighbors[curr_center].end();ne++) {
                 //removes the neighbor from the pool
                 unclustered.erase(std::remove(unclustered.begin(),unclustered.end(),*ne), unclustered.end()); //first occurence of n is removed. 
                 boolUnclustered[*ne]=false; // clustered
              }

             for(Ints::iterator ne=neighbors[curr_center].begin();ne!=neighbors[curr_center].end();ne++) { //for each neighbor
                 for(Ints::iterator unn=neighbors[*ne].begin();unn!=neighbors[*ne].end();unn++) {  //unclustered neighbor of neighbor
                     if (!boolUnclustered[*unn]) {
                         continue;
                        }
                     neighbors[*unn].remove(*ne);
                    }
              }
          }

         return (cluster_result);
 
}
                                                                                                   
 IMP::algebra::Vector2Ds SPE::get_contingency_table(IMP::optrep::Clusters cluster_result) {
	     /* Given the clustering and the identity of models in run1 and run2 creates the contingency table
         with 1 row per cluster and 1 column per sample.
         */
	 
	     IMP::algebra::Vector2Ds full_ctable;
    
         //initially create the table containing all clusters
		 for(unsigned int ic=0;ic<cluster_result.size();ic++)  { //iterate over each cluster
			 full_ctable.push_back(IMP::algebra::Vector2D(0.,0.));
		 }
          
		 for(unsigned int ic=0;ic<cluster_result.size();ic++)  { //iterate over each cluster
			 for(Ints::iterator cm=cluster_result[ic].cluster_members.begin();cm!=cluster_result[ic].cluster_members.end();cm++) {
				 // iterate over each cluster member of the cluster
				
				 if (std::binary_search(models_by_sample_[0].begin(),models_by_sample_[0].end(),*cm)) {
					 full_ctable[ic][0]+=1.0;
				 }
				 else if (std::binary_search(models_by_sample_[1].begin(),models_by_sample_[1].end(),*cm)) { 
				 	full_ctable[ic][1]+=1.0;
				 }
				
			 }
         }
		 
         // reduce the table by eliminating tiny clusters
         IMP::algebra::Vector2Ds reduced_ctable;
         Ints retained_clusters; // unnecessary variable but kept for sake of future use
  
          for(unsigned int ic=0;ic<cluster_result.size();ic++) { 
              if (full_ctable[ic][0]<=10.0 || full_ctable[ic][1]<=10.0)
                  continue;
		 
              reduced_ctable.push_back(IMP::algebra::Vector2D(full_ctable[ic][0],full_ctable[ic][1]));
              retained_clusters.push_back(ic);
  		  }
		
          return(reduced_ctable);
 }

Float SPE::percent_ensemble_explained(IMP::algebra::Vector2Ds ctable) {
          /* Check what proportion of the model space the exhaustiveness was on. 
         */
         if (ctable.size()==0) {
         	return (0.0);
         }
		 
		 Float num_models_in_contingency_table=0.0;
		 // sum of elements of a Vector2D. 
		 for(unsigned int ic=0;ic<ctable.size();ic++) { 
			 num_models_in_contingency_table=num_models_in_contingency_table + ctable[ic][0]+ ctable[ic][1];
		 }
              
          Float percent_clustered=num_models_in_contingency_table*100.0/Float(total_number_of_models_);
  
          return (percent_clustered);
}
 
IMP::optrep::ChiSquareTestResult SPE::test_sampling_exhaustiveness(IMP::algebra::Vector2Ds observed_contingency_table) {
          /* Chi2 test based on the contingency table.
          */
				
    		if (observed_contingency_table.size()==0) {
    	  	  	return(IMP::optrep::ChiSquareTestResult(0.0,1.0));
    		}
			
			else if(observed_contingency_table.size()==1) { //one single cluster
    	  	  	return(IMP::optrep::ChiSquareTestResult(1.0,0.0)); // trivially converged
    		}

			// Step 1. accumulate row and column sums of contingency table
			Floats row_sums;
			Floats column_sums; 
			Float table_total=0.0;
			
			for(unsigned int sc=0;sc<2;sc++) {
				column_sums.push_back(0.0);
			}
			
			for(unsigned int ic=0;ic<observed_contingency_table.size();ic++) { //cluster count
				row_sums.push_back(0.0);
			}
			
			for(unsigned int ic=0;ic<observed_contingency_table.size();ic++) { //cluster count		
				for(unsigned int sc=0;sc<2;sc++) { // sample count (0 and 1)
					row_sums[ic]+=observed_contingency_table[ic][sc];
					column_sums[sc]+=observed_contingency_table[ic][sc];
					table_total+=observed_contingency_table[ic][sc];
				}			
			}	

			// Step 2. calculate expected counts for each cell in the contingency table
			IMP::algebra::Vector2Ds expected_contingency_table;
			
			for(unsigned int ic=0;ic<observed_contingency_table.size();ic++) { //cluster count
				
				expected_contingency_table.push_back(IMP::algebra::Vector2D(0.,0.));
				
				for(unsigned int sc=0;sc<2;sc++) { // sample count (0 and 1)
					
					expected_contingency_table[ic][sc]=row_sums[ic]*column_sums[sc]/table_total;
				}		
			}
			
			Float chisquare=0.0;
			unsigned int degrees_of_freedom=observed_contingency_table.size()-1; // dof=(#rows-1)*(#cols-1)
			
			// Step 3. calculate the chisquare: expected-observed
			for(unsigned int ic=0;ic<observed_contingency_table.size();ic++) { //cluster count		
				for(unsigned int sc=0;sc<2;sc++) { // sample count (0 and 1)
					chisquare+=(observed_contingency_table[ic][sc]-expected_contingency_table[ic][sc])*
						(observed_contingency_table[ic][sc]-expected_contingency_table[ic][sc])/expected_contingency_table[ic][sc] ;
				}			
			}
	
		   boost::math::chi_squared curr_distribution(degrees_of_freedom);
		   Float pvalue = 1.0-boost::math::cdf(curr_distribution,chisquare);  
           Float cramersv=sqrt(chisquare/Float(total_number_of_models_));
      
          return(IMP::optrep::ChiSquareTestResult(pvalue,cramersv));
}

Float SPE::estimate_single_bead_precision(unsigned int global_bead_index,Float grid_size) {
    /* Estimate the sampling precision of the bead_index-th bead of protein domain.
    */
	
    IMP::optrep::DistanceMatrix dm=get_all_vs_all_distances(global_bead_index);

    Floats cutoffs;
	Float curr_cutoff=0.0;
	// Float curr_cutoff=dm.mindist; // the minimum distance is different for different beads, so standardizing it
	while(curr_cutoff<dm.maxdist) {
		cutoffs.append(curr_cutoff);
		curr_cutoff+=grid_size;	
	}
	
    Floats pvals;
	Floats cramersv;
	Floats populations;
	
    for(Floats::iterator c=cutoffs.begin();c!=cutoffs.end();c++) {
        IMP::optrep::Clusters clusters=precision_cluster(dm.distmat,*c);

        IMP::algebra::Vector2Ds ctable=get_contingency_table(clusters);
        IMP::optrep::ChiSquareTestResult ctr=test_sampling_exhaustiveness(ctable);

        Float percent_explained=percent_ensemble_explained(ctable);

        pvals.append(ctr.pvalue);
        cramersv.append(ctr.cramersv);
        populations.append(percent_explained);
	}

    Float sampling_precision = get_sampling_precision(cutoffs,pvals,cramersv,populations);
  
    return(sampling_precision);
	
}
bool SPE::is_commensurate(Float bead_diameter,Float bead_precision,Float xscale) {
    /* Check if the sampling precision of the bead is atmost xscale times the bead_diameter.
    If not the bead is imprecise and needs to be coarse-grained.
    */
    if (bead_precision> xscale*bead_diameter + 2.0) { //not just larger, but significantly larger than the representation precision
        return false;
    } 
    return true;
}
    
void SPE::get_imprecise_beads(Float xscale) {
    /* For each bead check if its size is commensurate with the sampling precision.
    If not, mark it as imprecise. 
    @param xscale used to define imprecise bead. imprecise bead has sampling_precision > xscale*bead_radius.
    */

	unsigned int total_number_of_global_beads=bead_diameter_.size();
	
	for(unsigned int global_bead_index=0;global_bead_index < total_number_of_global_beads;global_bead_index++) {  
        if (!is_commensurate(bead_diameter_[global_bead_index],bead_precisions_[global_bead_index],xscale)) {
        	bead_imprecise_.push_back(true); 
        }
            
        else {
        	bead_imprecise_.push_back(false); 
        }
	
    }
             
}
    

void SPE::print_bead_precisions(std::string out_file_name) {
    /* write bead index, bead precision and whether it is an imprecise bead to a file. 
    */
	FILE* out_file;
	out_file=fopen(out_file_name.c_str(), "w");
   
    unsigned int global_bead_index=0;
    for(unsigned int prot_index=0;prot_index<beads_per_protein_domain_.size();prot_index++) { // for each protein, domain
    	
		for(unsigned int bead_index=0;bead_index<beads_per_protein_domain_[prot];bead_index++) { // for each bead in a domain
			
			fprintf(out_file,"%s %s %u %.3f %d",components_calculate_precision[prot].first,
			components_calculate_precision[prot].second, bead_index, bead_precisions_[global_bead_index],
			int(bead_imprecise_[global_bead_index]));
						
			global_bead_index+=1;
			
		}	
		 
    }
            
    fclose(out_file);
}
    
void SPE::estimate_perbead_sampling_precision(Float grid_size=1.0) {
	/* For each required bead (selection residues mentioned in the class constructor), computes the sampling precision.
    Results are stored in the object's bead_precisions dictionary
    */
	   
	unsigned int total_number_of_global_beads=bead_diameter_.size();
	
	for(unsigned int global_bead_index=0;global_bead_index < total_number_of_global_beads;global_bead_index++) {  
			
			bead_precisions_.push_back(estimate_single_bead_precision(global_bead_index,grid_size));
	
    }
            
}

IMPOPTREP_END_NAMESPACE
