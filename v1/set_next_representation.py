import IMP
import os,sys,string,math
import argparse

def parse_args():
    
    parser = argparse.ArgumentParser(description="Get the next iteration's representation (bead map), based on results (precision) of the previous iteration. Example of usage: set_next_representation.py -r 11  xx xx. Flag -h for more details.")
    
    parser.add_argument("-r","--next_resolution",dest="next_resolution",required=True,type=int,help="The coarse-graining, i.e. maximum number of residues per bead to be tried in the next iteration.")
    
    parser.add_argument("-tf","--topology_file",action="store",dest="topology_file",required=True,help="")
                       
    result = parser.parse_args()
    
    return result

def 
    #TODO Read topology file. 
    #e.g.
    #Protein startres endres pdb chainname
    #Ab 1 34 1AVX.pdb A
    #Ab 35 100 1AVY.pdb A
    #Ab 101 150 - - 
    # The last line above is an unstructured region of protein Ab.
    
    
    return component
    
    
    
def set_next_representation():
     # Authors: Shruthi Viswanath
     
    ''' Set the representation for the next iteration of CG.
    '''
     
    # process input
    args=parse_args()
   
    if args.next_resolution==1: 
        # first iteration. Get bead map for the highest resolution beads.
        bmb=IMP.optrep.BeadMapBuilder(args.topology_file)
        # need protein list and PDB/fasta list
        
    else:
        
        
        
        
        get_protein_list_and_residues
        return
    
    
    
    
    
    
    
    
    
    
    




if __name__=="__init__":
    set_next_representation()
