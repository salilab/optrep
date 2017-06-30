import IMP
import IMP.atom
import IMP.rmf
import RMF
from subprocess import Popen
import os,sys,string,math
import shutil
import scipy.stats

class SamplingPrecisionEstimator(object):
    # Authors: Shruthi Viswanath
    
    ''' Given a set of good-scoring model, estimate the bead-wise sampling precision.
    '''
    
    def __init__(self,gsm_directory,components_calculate_precision={}):
        """Constructor.
        @param gsm_directory the directory containing good-scoring models. Should also contain the sample id file,
        the file containing the sample identity (1 or 2) for each good-scoring model RMF (model_sample_ids.txt)
        @param components_calculate_precision is ideally a PMI selection dictionary that selects the coordinates. 
        Right now it is just a list of protein names, e.g. {"A":"A","B":"B"}
        TODO Does not handle different selections of different proteins correctly yet? Easy to integrate this with existing 
        PMI selection tuples Dictionary with keys=groupname, values are tuples for aligning the structures e.g. {"Rpb1": (20,100,"Rpb1"),"Rpb2":"Rpb2"}
        """

        self.models_dir=gsm_directory
        self.models_by_sample={1:[],2:[]}

        self.components_calculate_precision=components_calculate_precision
        self.component_coords={}# store the coordinates of each model acoording to the bead index.
        # This dictionary is indexed by protein. The elements are lists of lists. 
        # Each list corresponds to a bead and is the list of all coordinates of good-scoring models corresponding to that bead.

        '''Initialize bead precision and radii dictionaries.'''
        self.bead_precisions={} #dictionary with key as protein and bead precision for each primitive in the protein
        self.bead_radii={} #dictionary with key as protein and size of each primitive in the protein

        for prot in components_calculate_precision:
            self.bead_precisions[prot]=[]
            self.bead_radii[prot]=[]
            self.component_coords[prot]=[]

        ''' Get the mapping of model index to sample number. '''
        sample_id_file=os.path.join(gsm_directory,"model_sample_ids.txt")

        sif=open(sample_id_file,'r')

        for ln in sif.readlines():
            model_index=int(ln.strip().split()[0])
            sample_number=int(ln.strip().split()[1])
        
            self.models_by_sample[sample_number].append(model_index)
            
        sif.close()

    def _included_protein(self,chain_full_name,protein_list):
        ''' Check if the current chain should be included for calculating bead precisions, based on the list of proteins mentioned. 
        '''
        for prot in protein_list:
            if prot in chain_full_name:
                return True

        return False

    def load_coordinates_and_bead_sizes_from_model_files(self):
        ''' Load all the coordinates from all good scoring models.
        Store them by bead type so that it is easy to calculate RMSD/precision.
        TODO does not do alignment yet
        TODO does not handle multiscale properly yet? 
        TODO Take care of ambiguity.
        TODO Take care of multi-state protein.
        TODO  There may be memory problem for large assemblies and large number of models. Need to refactor. 
        TODO Make the code faster by integrating with the GPU based RMSD calculator code.
        '''

        for ct,i in enumerate(sorted(self.model_sample_ids.keys())):
            m=IMP.Model()
            mdl_i=str(i)+'.rmf3'
    
            fh_i=RMF.open_rmf_file_read_only(mdl_i)
    
            hier_i=IMP.rmf.create_hierarchies(fh_i,m)[0]

            IMP.rmf.load_frame(fh_i, 0)

            bead_index=0

            for state_i in hier_i.get_children():
                for chain_i in state_i.get_children():
                    if _included_protein(chain_i.get_name(),self.components_calculate_precision):
                        for a in IMP.atom.get_leaves(chain_i):
                            curr_coords=IMP.core.XYZ(a).get_coordinates()
                            curr_radius=IMP.core.XYZ(a).get_radius()

                            if ct==1: #first model, need to initialize the coordinate and radius lists for each bead index
                                self.component_coords[prot].append([curr_coords])
                                self.bead_radii[prot].append([curr_radius])
                            else:
                                self.component_coords[prot][bead_index].append(curr_coords)
                                self.bead_radii[prot][bead_index].append(curr_radius)

                            # in any case increment bead index
                            bead_index += 1

            del m,hier_i

    def get_all_vs_all_distances(self,prot,bead_index):
        ''' Return the distance matrix, minimum and maximum distance per bead.
        '''
        num_models=len(self.model_sample_ids)
        
        distmat=numpy.zeros((num_models,num_models))
        #get the all vs all distance matrix, and minimum and maximum distances, of the bead index

        mindist=1000000.0
        maxdist=0.0

        for i in range(num_models-1):
            for j in range(i+1,num_models):                
                dist=IMP.core.XYZ.get_distance(self.component_coords[prot][bead_index][i],self.component_coords[prot][bead_index][j])
                distmat[i][j]=dist

                if dist<mindist:
                    mindist=dist

                if dist>maxdist:
                    maxdist=dist

        return (distmat,mindist,maxdist)


    def get_sampling_precision(self,cutoffs,pvals,cramersv,populations):
        ''' Given the 3 criteria for each cutoff, returns the sampling precision.
        This is the first cutoff at which all 3 criteria are satisfied.
        '''
        sampling_precision=100000.0

        for i in range(len(cutoffs)):
            if populations[i]>80.0 and (pvals[i]>0.05 or cramersv[i]<0.10):
                if sampling_precision>cutoffs[i]:
                    sampling_precision=cutoffs[i]
            else: # takes care of local minima? We want to choose the first time at which the above condition was satisfied
                sampling_precision=100000.0

        return sampling_precision

    def precision_cluster(distmat,rmsd_cutoff):
        ''' Perform distance threshold-based clustering given the distance matrix and RMSD cutoff        '''

        num_models=numpy.shape(distmat)[0]

        #Populate the neighbors ofa given model
        neighbors=[]
        for count in range(num_models):
            neighbors.append([count])  # model is a neighbor of itself

        for i in range(num_models-1):
            for j in range(i+1,num_models):
                if distmat[i][j]<=rmsd_cutoff: # accepted to be a neighbor
                    neighbors[i].append(j)
                    neighbors[j].append(i)

        #Get the cluster with the most neighbors, and iterate
        #Initially each model is unclustered
        unclustered=[]
        boolUnclustered=[]
        for i in range(num_models):
            unclustered.append(i)
            boolUnclustered.append(True)

        cluster_members=[] # list of lists : one list per cluster
        cluster_centers=[]

        while len(unclustered)>0:
            # get cluster with maximum weight
            max_neighbors=0
            currcenter=-1
            for eachu in unclustered:  # if multiple clusters have same maxweight this tie is broken arbitrarily! 
                if len(neighbors[eachu])>max_neighbors:
                    max_neighbors=len(neighbors[eachu])
                    currcenter=eachu

            #form a new cluster with u and its neighbors
            cluster_centers.append(currcenter)
            cluster_members.append([n for n in neighbors[currcenter]])

            #update neighbors 
            for n in neighbors[currcenter]:
                #removes the neighbor from the pool
                unclustered.remove(n) #first occurence of n is removed. 
                boolUnclustered[n]=False # clustered

            for n in neighbors[currcenter]:
                for unn in neighbors[n]: #unclustered neighbor
                    if not boolUnclustered[unn]:
                        continue
                    neighbors[unn].remove(n)

        return cluster_centers,cluster_members
                                                                                                 
    def get_contingency_table(cluster_members,sample1_models,sample2_models):
        ''' Given the clustering and the identity of models in run1 and run2 creates the contingency table
        with 1 row per cluster and 1 column per sample.
        '''

        num_clusters=len(cluster_members)

        # initially create the table containing all clusters
        full_ctable=numpy.zeros((num_clusters,2))

        for ic,cluster in enumerate(cluster_members):
            for member in cluster:
                if member in sample1_models:
                    full_ctable[ic][0]+=1.0
                elif member in sample2_models:
                    full_ctable[ic][1]+=1.0

        # reduce the table by eliminating tiny clusters
        reduced_ctable=[]
        retained_clusters=[]

        for i in range(num_clusters):
            if full_ctable[i][0]<=10.0 or full_ctable[i][1]<=10.0:
                continue
            reduced_ctable.append([full_ctable[i][0],full_ctable[i][1]])
            retained_clusters.append(i)

        return numpy.array(reduced_ctable),retained_clusters

    def percent_ensemble_explained(ctable,total_num_models):
        ''' Check what proportion of the model space the exhaustiveness was on. 
        '''
        if len(ctable)==0:
            return 0.0
        percent_clustered=float(numpy.sum(ctable,axis=0)[0]+numpy.sum(ctable,axis=0)[1])*100.0/float(total_num_models)

        return percent_clustered

    def test_sampling_exhaustiveness(contingency_table,total_num_models):
        ''' Chi2 test based on the contingency table.
        '''

        if len(contingency_table)==0:
            return 0.0,1.0

        ct = numpy.transpose(contingency_table)

        [chisquare,pvalue,dof,expected]=scipy.stats.chi2_contingency(ct)

        if dof==0.0:
            cramersv=0.0 #converged, one single cluster
        else:
            cramersv=math.sqrt(chisquare/float(total_num_models))

        return(pvalue,cramersv)

    def estimate_single_bead_precision(self,prot,bead_index,grid_size):
        ''' Estimate the sampling precision of the bead_index-th bead of protein prot.
        '''
        num_models=len(self.model_sample_ids)
        
        (distmat_bead,mindist_bead,maxdist_bead)=get_all_vs_all_distances(self,prot,bead_index)

        cutoffs=numpy.arange((mindist,maxdist,grid_size)
        
        pvals=[]
        cramersv=[]
        populations=[]

        for c in cutoffs:
            cluster_centers,cluster_members=precision_cluster(distmat_bead,c)

            ctable,retained_clusters=get_contingency_table(cluster_members,self.models_by_sample[1],self.models_by_sample[2]) 

            (pval,crv)=test_sampling_exhaustiveness(ctable,num_models)

            percent_explained= percent_ensemble_explained(ctable,num_models)

            pvals.append(pval)
            cramersv.append(crv)
            populations.append(percent_explained)

        sampling_precision = self.get_sampling_precision(cutoffs,pvals,cramersv,populations)

        return sampling_precision

    def is_commensurate(bead_radius,bead_precision,xscale):
        ''' Check if the sampling precision of the bead is atmost xscale times the bead_radius.
        If not the bead is imprecise and needs to be coarse-grained.
        '''
        if bead_precision>xscale*bead_radius:
            return False
        return True
              
        
    def get_imprecise_beads(self,xscale):
        ''' For each protein, get the list of bead indices for which the bead sizes were not commensurate with the sampling precision.
        @param xscale used to define imprecise bead. imprecise bead has sampling_precision > xscale*bead_radius.
        '''
        imprecise_beads={}
        
        for prot in self.bead_radii:
            imprecise_beads[prot]=[]
            
            for bead_index in range(len(self.bead_radii[prot])):
                 if not self.is_commensurate(self.bead_radii[prot][bead_index],self.bead_precisions[prot][bead_index],xscale):
                     imprecise_beads[prot].append(bead_index)

        return imprecise_beads

    def estimate_perbead_sampling_precision(self,grid_size=2.5):
        ''' For each required bead (selection residues mentioned in the class constructor), computes the sampling precision.
        Results are stored in the object's bead_precisions dictionary
        '''

        self.load_coordinates_and_bead_sizes_from_model_files()

        for prot in self.components_calculate_precision:
            for bead_index in range(len(self.bead_radii[prot])):
                self.bead_precisions[prot][bead_index]=self.estimate_single_bead_precision(prot,bead_index,grid_size)

  
