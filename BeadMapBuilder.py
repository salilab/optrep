import IMP
import IMP.atom
import IMP.rmf
import IMP.pmi
import IMP.pmi.tools
import IMP.pmi.topology
import os,sys,string,math
 	
class BeadMapBuilder(object):
    # Authors: Shruthi Viswanath
    
    ''' Setup the bead to residue map given a topology file, for a given molecule.
    Alter the bead-residue map when given a set of residues to coarse-grain to the next stage.
    '''
    
    def __init__(self):
        """Constructor.
        """
        
        # All the below dictionaries have keys as a tuple: (protein_name,domain_name) 
        self.bead_maps={} # dictionary of beadmaps  for each protein,domain.  
        # The values are a list of lists with each member list representing a bead. 
        # Each member list representing a bead has two members: [start residue index, end residue index]. 
        self.imprecise_beads={} # list of bead indices that should be coarse-grained in the next round, for each (protein,domain) key.
        self.merged_bead={}  # list of booleans, one per bead, for each (protein,domain) key. The boolean says whether that bead was merged with another in the current iteration.
    
            
    def set_bead_map_from_topology_file(self,topology_file):
        '''Given a topology file, read its contents to the bead map object.
        @param topology_file list of proteins, domains, residue ranges, pdb files and chains
        '''
       
        fl=open(topology_file,'r')
              
        for ln in fl.readlines():
            fields=ln.strip().split()
           
            protein=fields[0]
            domain=fields[1]
            start_residue_index=int(fields[3])
            end_residue_index=int(fields[4])
            
            self.set_bead_map_given_residue_range((protein,domain),start_residue_index,end_residue_index)
    
        fl.close()
        
    def set_bead_map_from_beadmap_file(self,bm_file):
        '''Given a bead map file, read its contents to the bead map object.
        @param bm_file bead map used in the previous iteration
        '''
   
        self.bead_maps={}

        fl=open(bm_file,'r')

        for ln in fl.readlines():
            fields=ln.strip().split()
    
            protein=fields[0]
            domain=fields[1]
            bead=[int(fields[2]),int(fields[3])]
       
            protein_domain_key=(protein,domain)

            if protein_domain_key in self.bead_maps:
                self.bead_maps[protein_domain_key].append(bead) #assume the beads are written in order in the file
            else:
                self.bead_maps[protein_domain_key]=[bead]
        
        fl.close()
        
    def set_bead_map_given_residue_range(self,protein_domain_key,start_residue_index,end_residue_index):
        '''Given a range of contiguous residues, set the bead map of the protein domain to this list of residues.
        @param start_residue_index starting residue of the bead map list
        @param end_residue_index end residue of the bead map list
        '''
        # reset the bead map
        self.bead_maps[protein_domain_key]=[]

        for i in range(start_residue_index,end_residue_index+1):
            self.bead_maps[protein_domain_key].append([i,i]) 
    
    def _get_bead_size(self,protein_domain_key,bead_index):
        ''' Given the index of the bead in the list returns the number of residues in the bead.
        @param protein_domain_key is a tuple (prot,domain)
        '''
        return self.bead_maps[protein_domain_key][bead_index][1]-self.bead_maps[protein_domain_key][bead_index][0]+1 # since we include both starting and ending indices
        
    def _number_of_residues_in_protein_domain(self,protein_domain_key):
        ''' Return the number of residues in the prot_domain according to the prot_domain range.
        TODO does not account for missing residues.
        '''
        size = self.bead_maps[protein_domain_key][len(self.bead_maps[protein_domain_key])-1][1]-self.bead_maps[protein_domain_key][0][0]+1
        return size

    def coarsen_bead(self,protein_domain_key,bead_index,cg_bead_size):
        ''' Coarsen the bead with given bead index to the cg_bead_size.
        Use knowledge of the bead list and size of each bead currently in the list.
        Merge along the backbone, once on the left and once on the right till we reach the desired bead size.
        This method makes changes to the bead_maps and merged_bead lists in the object.
        '''
    
        num_residues_bead=self._get_bead_size(protein_domain_key,bead_index)
        
        counter=0 # number of cycles. 1 cycle = 1 bead to the left and right each merged with current bead
        
        while num_residues_bead<cg_bead_size: 
                    
                if num_residues_bead==self._number_of_residues_in_protein_domain(protein_domain_key):
                    break

                counter+=1
                
                for neighbor_index in [bead_index-counter,bead_index+counter]:
                    #print bead_index, neighbor_index,num_residues_bead

                    # bumping into boundaries, or into an already merged neighbor
                    if neighbor_index>=len(self.bead_maps[protein_domain_key]) or neighbor_index<0 or self.merged_bead[protein_domain_key][neighbor_index]:
                        continue
                    
                    # Update the residues contained in the current bead by incorporating residues in the neighboring bead
                    self.bead_maps[protein_domain_key][bead_index][0]=self.bead_maps[protein_domain_key][min(bead_index,neighbor_index)][0]
                
                    self.bead_maps[protein_domain_key][bead_index][1]=self.bead_maps[protein_domain_key][max(bead_index,neighbor_index)][1]
             
                    self.merged_bead[protein_domain_key][neighbor_index]=True # this is the index of the bead that was merged with the current imprecise bead index
                            
                    num_residues_bead=self._get_bead_size(protein_domain_key,bead_index)
                    
                     # check if we did not already reach the maximum CG size       
                    if num_residues_bead>=cg_bead_size: 
                        break
                     
    def update_single_bead_map(self,protein_domain_key,cg_bead_size):
       '''Get the bead map corresponding to the current level of coarse-graining. If we are coarse-graining for the first time, can do highest resolution beads. Else can use the previous bead map and coarse-grain it incrementally. 
       '''
       
       if protein_domain_key not in self.imprecise_beads:
           return False
              
       # Otherwise do incremental CG of the beadmap.
       # Assumes the imprecise_beads are set for the protein
              
       if len(self.imprecise_beads[protein_domain_key])==0:
           return False
           
       if protein_domain_key not in self.bead_maps:
           raise ValueError("Before attempting to coarse-grain, need to set the bead map for %s",protein_domain_key)
           return False
        
       if len(self.bead_maps[protein_domain_key])==0:
           raise ValueError("Before attempting to coarse-grain, need to set the bead map for %s",protein_domain_key)
           return False

       # Initialize the boolean that says which beads were merged
       self.merged_bead[protein_domain_key]=[False for i in range(len(self.bead_maps[protein_domain_key]))] # Boolean that says if a bead index was merged recently in the current round of coarse-graining
              
       for imprecise_bead_index in self.imprecise_beads[protein_domain_key]:
            # Make this bead bigger by combining adjacent beads, one from left and one from right  
            # along sequence till the desired (next level) bead size is reached.

            if self.merged_bead[protein_domain_key][imprecise_bead_index]: # this bead was already merged in the current step. 
                continue
            
            self.coarsen_bead(protein_domain_key,imprecise_bead_index,cg_bead_size) # changes the bead_maps and the merged_bead lists
      
       map_changed=False
       new_bead_map=[]
       
       # Now we create the new coarsened bead map by removing all beads that were marked merged      
       for i in range(len(self.bead_maps[protein_domain_key])):
           if not self.merged_bead[protein_domain_key][i]:
               new_bead_map.append(self.bead_maps[protein_domain_key][i])
        
           else:
               map_changed=True
               
       self.bead_maps[protein_domain_key]=new_bead_map    
   
       # After getting a new bead map, need to reset the imprecise_beads of protein domain
       self.imprecise_beads[protein_domain_key]=[] 
       self.merged_bead[protein_domain_key]=[]
       
       return map_changed

    def update_all_bead_maps(self,protein_domains_to_update,cg_bead_size):
        ''' Given a list of protein domains to update, update the representations of each of these to make the 
        maximum bead size cg_bead_size.'''
        
        updated=False
        
        for protein_domain_key in protein_domains_to_update:
            current_map_updated = self.update_single_bead_map(protein_domain_key,cg_bead_size)
            
            if current_map_updated:
                updated=True

        return updated
    
    def set_imprecise_beads(self,imprecise_bead_dict):
        ''' Since we traverse the protein domain in order of bead index,the imprecise_bead_list is sorted by construction.        
        '''
    
        for protein_domain_key in imprecise_bead_dict:
            self.imprecise_beads[protein_domain_key]=imprecise_bead_dict[protein_domain_key]    
  
        
    def write_bead_map_to_file(self,output_file):
        ''' Write bead maps to file. Need to access this in the case of the next iteration of coarse-graining.
        '''
        f=open(output_file,'w')
        
        for p in self.bead_maps:
            for bead in self.bead_maps[p]:
                print >>f,p[0],p[1],bead[0],bead[1]
             
        f.close()
               
    def show_bead_map(self):
        ''' Output on the screen.'''
        for p in self.bead_maps:
            for bead in self.bead_maps[p]:
                print p[0],p[1],bead[0],bead[1]

