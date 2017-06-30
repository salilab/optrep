
import IMP
import RMF
import IMP.atom
import IMP.rmf
import IMP.pmi
import IMP.pmi.topology
import IMP.pmi.dof
import IMP.pmi.macros
import IMP.pmi.restraints
import IMP.pmi.restraints.stereochemistry

def getAA_Alphabet(resnam):

	if resnam=='ALA':
		return 'A'
	elif resnam=='CYS':
		return 'C'
	elif resnam=='ASP':
		return 'D'
	elif resnam=='GLU':
		return 'E'
	elif resnam=='PHE':
		return 'F'
	elif resnam=='GLY':
		return 'G'
	elif resnam=='HIS':
		return 'H'
	elif resnam=='ILE':
		return 'I'
	elif resnam=='LYS':
		return 'K'
	elif resnam=='LEU':
		return 'L'
	elif resnam=='MET':
		return 'M'
	elif resnam=='ASN':
		return 'N'
	elif resnam=='PRO':
		return 'P'
	elif resnam=='GLN':
		return 'Q'
	elif resnam=='ARG':
		return 'R'
	elif resnam=='SER':
		return 'S'
	elif resnam=='THR':
		return 'T'
	elif resnam=='VAL':
		return 'V'
	elif resnam=='TRP':
		return 'W'
	elif resnam=='TYR':
		return 'Y'
	else:
		return 'G'
	# default is GLYcine

def getSequenceFromPDB(pdbfile,chain):
        seq = ""
        ms = IMP.Model()
        h = IMP.atom.read_pdb(pdbfile,ms)
        for ch in h.get_children():
             if chain not in ch.get_name():
                continue
	     for res in ch.get_children():
	         seq = seq + getAA_Alphabet(res.get_name())

        print "going",seq
	return seq

lpdb='input/7CEI_l_u.pdb'

# Preliminaries: read sequences, create a system and a state
mdl = IMP.Model()
s = IMP.pmi.topology.System(mdl)
st = s.create_state()

mol=st.create_molecule('B',chain_id='B',sequence=getSequenceFromPDB(lpdb,"B"))


# Slice a molecule (uses python-style 0-ordering)
myres1 = mol[0:10]

# Use PDB numbering with residue_range (inclusive on both ends)
#myres2 = mol.residue_range('446','456')

#mol.add_representation(myres1,
                       #resolutions=[1])
myres2 = mol[10:11] | mol[21:24]

#mol.add_structure(lpdb,chain_id="B",soft_check=True)
mol.add_representation(myres1,resolutions=[4])

# When you have decided all representations, call build()
#  This returns an IMP hierarchy
hier = s.build()

# View your creation with this function
IMP.atom.show_with_representations(hier)



