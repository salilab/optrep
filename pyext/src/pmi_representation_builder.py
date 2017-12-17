import IMP
import RMF
import os,string,sys
import IMP.atom
import IMP.rmf
import IMP.pmi
import IMP.pmi.tools
import IMP.pmi.topology
import IMP.pmi.dof
import IMP.pmi.macros
import IMP.optrep.BeadMapBuilder

def add_representation(state,input_dir,fasta_file,topology_file,bead_map_file):
   
    if fasta_file:
        sequences = IMP.pmi.topology.Sequences(fasta_file)
    
    bmb=IMP.optrep.BeadMapBuilder.BeadMapBuilder()
    bmb.set_bead_map_from_beadmap_file(bead_map_file)
    bead_map=bmb.bead_maps
  
    mols=[]
    
    tf=open(topology_file,'r')
    
    for ln in tf.readlines():
      
        protein,domain,protein_chain,fastakey,startres_domain,endres_domain,pdb,pdb_chain,resolution,clr=ln.strip().split()
            
             
        if protein not in state.get_molecules(): #TODO not considered multiple copies yet

            if fastakey!="-" and protein_chain!="-":
                curr_mol=state.create_molecule(protein,chain_id=protein_chain,sequence=sequences[fastakey])
                
            elif fastakey=="-" and protein_chain!="-":
                curr_mol=state.create_molecule(protein,chain_id=protein_chain)    
                
            elif fastakey!="-" and protein_chain=="-":
                curr_mol=state.create_molecule(protein,sequence=sequences[fastakey]) 
            
            else: 
                curr_mol=state.create_molecule(protein)

            
            mols.append(curr_mol)
            
        else: # multiple domains of same protein
            curr_mol=state.get_molecule(protein)
        
        # getcolor in correct format
        try:
            clr=float(clr)
            
        except ValueError:
            if '_' in clr:
                clr = clr.replace('_',' ')
                
        # if structure exists for this domain, add it
        if pdb!="-":

            if pdb_chain=="-":
                pdb_chain=" "

            atomic = curr_mol.add_structure(os.path.join(input_dir,pdb),chain_id=pdb_chain,res_range=(int(startres_domain),int(endres_domain)+1),soft_check=True)
       
        if resolution!="bm": # structure exists, or need not be optimized. Used fixed res representation
            multiscale_res = [int(rsn) for rsn in resolution.split(',')]
           
            curr_mol.add_representation(atomic,resolutions=multiscale_res,color=clr)
            
        # else add the representation from bead map 
        else:
            
            for bead in bead_map[(protein,domain)]:
                startres_bead=bead[0] 
                endres_bead=bead[1]
                                
                #print curr_mol.residue_range(str(startres_bead),str(endres_bead))
                #print curr_mol[startres_bead-1:endres_bead]
                
                #curr_mol.add_representation(curr_mol.residue_range(str(startres_bead),str(endres_bead)),resolutions=[endres_bead-startres_bead+1],color=clr) 
                if pdb!="-":
                    select_residues = curr_mol[startres_bead-1:endres_bead] & atomic 
                else:
                    select_residues = curr_mol[startres_bead-1:endres_bead]
             
                curr_mol.add_representation(select_residues,resolutions=[endres_bead-startres_bead+1],color=clr)
                #TODO not considered multi-scaling yet
                #TODO not considering missing residues in PDB here
     
    return mols
    
