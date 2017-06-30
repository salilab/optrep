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
        
        self.bead_maps={} # dictionary of beadmap  for each protein domain.  Format : each protein domain in the system that corresponds to a separate line in the topology file is a key. 
        # The values are a list of lists with each member list representing a bead. 
        # Each member list representing the bead has two members: [start residue index, end residue index]. 
        self.imprecise_beads={} # dictionary of bead indices that should be coarse-grained in the next round. Indexed by protein key.
        self.merged_bead={}  # Indexed by protein key. One boolean entry for each bead in the bead map. Says whether it was ignored in the current iteration.
        
        for p in component_list:
            self.bead_maps[p]=[] # each bead map is a list of tuples: (start residue index,end residue index)
            self.imprecise_beads[p]=[]
            self.merged_bead[p]=[] 
            
    def set_bead_map_from_file(self,bm_file,protein):
        '''Given a bead map file, read its contents to the bead map object.
        @param bm_file bead map used in the previous iteration
        '''
        # reset the bead map
        self.bead_maps[protein]=[]

        f=open(bm_file,'r')

        for ln in f.readlines():
            fields=ln.strip().split()
            prot=fields[0]
            
            if prot!=protein:
                continue
            bead=[int(fields[1]),int(fields[2])]
            self.bead_maps[protein].append(bead) #assume the beads are written in order in the file
        f.close()
        
    def set_bead_map_given_residue_range(self,protein,start_residue_index,end_residue_index):
        '''Given a range of contiguous residues, set the bead map of the protein to this list of residues.
        @param start_residue_index starting residue of the bead map list
        @param end_residue_index end residue of the bead map list
        '''
        # reset the bead map
        self.bead_maps[protein]=[]

        for i in range(start_residue_index,end_residue_index+1):
            self.bead_maps[protein].append([i,i]) 
    
    #DEPRECATED METHODS. These are not necessary since we have the topology file. 
    #def _get_residue_indices_from_PDB(self,chain,pdb_file):
        #''' Gets all the residue indices from a model in the PDB file. 
        #TODO assumes 1 protein = 1 PDB, does not account for missing residues.
        #'''

        #if not chain:
            #raise ValueError("No chain specified along with %s",pdb_file)
            #return []

        #mol = IMP.atom.read_pdb(pdb_file,self.model,
                               #IMP.atom.AndPDBSelector(IMP.atom.ChainPDBSelector(chain),IMP.atom.CAlphaPDBSelector()))
        
        ## can do either Selection or Hierarchy
        #selected = IMP.atom.Selection(mol,chain=chain,atom_type=IMP.atom.AtomType('CA'))
        #residue_indices=[]

        #for p in selected.get_selected_particles():
            #res = IMP.atom.Residue(IMP.atom.Atom(p).get_parent())
            #residue_indices.append(res.get_index())

        #return residue_indices
        
    #def _get_residue_indices_from_FASTA(self,chain,fasta_file):
        #''' TODO not implemented/tested. Caution before using. '''

        #if not chain:
            #raise ValueError("No chain specified along with %s",fasta_file)
            #return []

        #fasta_sequences=IMP.pmi.topology.Sequences(fasta_file)

        #residue_indices=[i for i in range(1,len(fasta_sequences[chain])+1)]
        #return residue_indices

    def _get_bead_size(self,protein,bead_index):
        ''' Given the index of the bead in the list returns the number of residues in the bead.
        '''
        return self.bead_maps[protein][bead_index][1]-self.bead_maps[protein][bead_index][0]+1 # since we include both starting and ending indices
        
    def _number_of_residues_in_protein(self,protein):
        ''' Return the number of residues in the protein according to the protein range.
        TODO does not account for missing residues.
        '''
        size = self.bead_maps[protein][len(self.bead_maps[protein])-1][1]-self.bead_maps[protein][0][0]+1
        return size

    def coarsen_bead(self,protein,bead_index,cg_bead_size):
        ''' Coarsen the bead with given bead index to the cg_bead_size.
        Use knowledge of the bead list and size of each bead currently in the list.
        Merge along the backbone, once on the left and once on the right till we reach the desired bead size.
        This method makes changes to the bead_maps and merged_bead lists in the object.
        '''
    
        num_residues_bead=self._get_bead_size(protein,bead_index)
        
        counter=0 # number of cycles. 1 cycle = 1 bead to the left and right each merged with current bead
        
        while num_residues_bead<cg_bead_size: 
                    
                if num_residues_bead==self._number_of_residues_in_protein(protein):
                    break

                counter+=1
                
                for neighbor_index in [bead_index-counter,bead_index+counter]:
                    #print bead_index, neighbor_index,num_residues_bead

                    # bumping into boundaries, or into an already merged neighbor
                    if neighbor_index>=len(self.bead_maps[protein]) or neighbor_index<0 or self.merged_bead[protein][neighbor_index]:
                        continue
                    
                    # Update the residues contained in the current bead by incorporating residues in the neighboring bead
                    self.bead_maps[protein][bead_index][0]=self.bead_maps[protein][min(bead_index,neighbor_index)][0]
                
                    self.bead_maps[protein][bead_index][1]=self.bead_maps[protein][max(bead_index,neighbor_index)][1]
             
                    self.merged_bead[protein][neighbor_index]=True # this is the index of the bead that was merged with the current imprecise bead index
                            
                    num_residues_bead=self._get_bead_size(protein,bead_index)
                    
                     # check if we did not already reach the maximum CG size       
                    if num_residues_bead>=cg_bead_size: 
                        break
                     
    def create_bead_map(self,simplified_topology_file=None,cg_bead_size=1):
       '''Get the bead map corresponding to the current level of coarse-graining. If we are coarse-graining for the first time, can do highest resolution beads. Else can use the previous bead map and coarse-grain it incrementally. 
       '''
       if not protein:
           raise ValueError("Missing protein name in create_bead_map. Make sure you add protein name to constructor")
           return
       
       # First time we are creating beads. Use the PDB/Fasta file
       if cg_bead_size==1:
         
            #if pdb_file:
                #residues=self._get_residue_indices_from_PDB(chain,pdb_file)
                    
            #elif fasta_file:
                #residues=self._get_residue_indices_from_FASTA(chain,fasta_file) 
                    
            #else:
                #raise ValueError("Need either PDB or FASTA file to get the residue indices of %s",protein)
            
            #for i,r in enumerate(residues):
                    #self.bead_maps[protein].append([r,r]) # including starting and ending residue
                    
            

            return
       
       # Otherwise do incremental CG of the beadmap.
       # Assumes the imprecise_beads are set for the protein
              
       if len(self.imprecise_beads[protein])==0:
           raise ValueError("Before attempting to coarse-grain, need to find out the imprecise beads for %s",protein)
       
       # Initialize the boolean that says which beads were merged
       self.merged_bead[protein]=[False for i in range(len(self.bead_maps[protein]))] # Boolean that says if a bead index was merged recently in the current round of coarse-graining
              
       for imprecise_bead_index in self.imprecise_beads[protein]:
            # Make this bead bigger by combining adjacent beads, one from left and one from right  
            # along sequence till the desired (next level) bead size is reached.
      
            if self.merged_bead[protein][imprecise_bead_index]: # this bead was already merged in the current step. 
                continue
            
            self.coarsen_bead(protein,imprecise_bead_index,cg_bead_size) # changes the bead_maps and the merged_bead lists
      
       # Now we create the new coarsened bead map by removing all beads that were marked merged      
       new_bead_map=[self.bead_maps[protein][i] for i in range(len(self.bead_maps[protein])) if not self.merged_bead[protein][i]]
   
       self.bead_maps[protein]=new_bead_map    
            
       # After getting a new bead map, need to reset the imprecise_beads of protein
       self.imprecise_beads[protein]=[] 
       self.merged_bead[protein]=[]

    def set_imprecise_beads(self,sorted_imprecise_bead_list):
        ''' Since we traverse the protein in order of bead index,the imprecise_bead_list is sorted by construction.        
        '''
        for prot in sorted_imprecise_bead_list:
            self.imprecise_beads[protein]=sorted_imprecise_bead_list    
   
     
    def write_bead_map_to_file(self,output_file):
        ''' Write bead maps to file. Need to access this in the case of the next iteration of coarse-graining.
        '''
        f=open(output_file,'w')
        
        for p in self.bead_maps:
            for bead in self.bead_maps[p]:
                print >>f,p,bead[0],bead[1]
             
        f.close()
               
    def show_bead_map(self):
        ''' Output on the screen.'''
        for p in self.bead_maps:
            for bead in self.bead_maps[p]:
                print p,bead[0],bead[1]

