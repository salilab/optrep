from __future__ import print_function
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
    
            
    def set_bead_map_from_topology_file(self,topology_file,resolution=1):
        '''Given a topology file, read its contents to the bead map object.
        @param topology_file list of proteins, domains, residue ranges, pdb files and chains
        @param resolution the size of bead to be created
        '''
       
        fl=open(topology_file,'r')
              
        for ln in fl.readlines():
            fields=ln.strip().split()
           
            protein=fields[0]
            domain=fields[1]
            start_residue_index=int(fields[4])
            end_residue_index=int(fields[5])
            
            self.set_bead_map_given_residue_range((protein,domain),start_residue_index,end_residue_index,resolution)
    
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
        
    def set_bead_map_given_residue_range(self,protein_domain_key,start_residue_index,end_residue_index,resolution):
        '''Given a range of contiguous residues, set the bead map of the protein domain to this list of residues.
        @param start_residue_index starting residue of the bead map list
        @param end_residue_index end residue of the bead map list
        '''
        # reset the bead map
        self.bead_maps[protein_domain_key]=[]

        for i in range(start_residue_index,end_residue_index+1,resolution):
            self.bead_maps[protein_domain_key].append([i,min(i+resolution-1,end_residue_index)]) 
    
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
       
       # current bead on which we are coarse-graining
       curr_bead_index = -1 # uninitialized value
       
       #i is the current bead index
       for i in range(len(self.bead_maps[protein_domain_key])):
           
           # it is an imprecise bead
           if i in self.imprecise_beads[protein_domain_key]: #TODO this can be faster by going serially through imprecise beads
               
               #if we are not cg'ing any bead right now, make this the bead to cg.
               if curr_bead_index==-1:
                   curr_bead_index=i
           
               #cg the current bead to cg by adding residues from the current bead
               self.bead_maps[protein_domain_key][curr_bead_index][1]=self.bead_maps[protein_domain_key][i][1]
               
               # mark bead as merged if it is not the bead we are cg'ing
               if curr_bead_index!=i:
                   self.merged_bead[protein_domain_key][i]=True
                   
               # if maximum bead size reached, or the next bead is precise, need to stop cg'ing and reset the current bead to cg.
               if self._get_bead_size(protein_domain_key,curr_bead_index)>=cg_bead_size or (i+1) not in self.imprecise_beads[protein_domain_key]:
                   curr_bead_index=-1     
               
       # add all non-merged beads to new_bead_map
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

    def update_all_bead_maps(self,cg_bead_size):
        ''' Given a list of protein domains to update, update the representations of each of these to make the 
        maximum bead size cg_bead_size.'''
        
        updated=False
        
        for protein_domain_key in self.imprecise_beads:
            current_map_updated = self.update_single_bead_map(protein_domain_key,cg_bead_size)
            
            if current_map_updated:
                updated=True

        return updated
    
    def set_imprecise_beads_from_file(self,imprecise_beads_file):
        ''' Will only have the beads for components that we want to update. 
        '''
    
        ibf=open(imprecise_beads_file,'r')
        for ln in ibf.readlines():
            fields=ln.strip().split()
            is_imprecise  = int(fields[4])
            if not is_imprecise:
                continue
            protein_domain_key=(fields[0],fields[1])
            if protein_domain_key in self.imprecise_beads:
                self.imprecise_beads[protein_domain_key].append(int(fields[2]))
            else:
                self.imprecise_beads[protein_domain_key]=[int(fields[2])]
                        
        ibf.close()

    def write_bead_map_to_file(self,output_file):
        ''' Write bead maps to file. Need to access this in the case of the next iteration of coarse-graining.
        '''
        with open(output_file, 'w') as f:
            for p in sorted(self.bead_maps):
                for bead in sorted(self.bead_maps[p]):
                    print(p[0],p[1],bead[0],bead[1], file=f)
               
    def show_bead_map(self):
        ''' Output on the screen.'''
        for p in self.bead_maps:
            for bead in self.bead_maps[p]:
                print(p[0],p[1],bead[0],bead[1])

