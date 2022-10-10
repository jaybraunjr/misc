import MDAnalysis as mda
from MDAnalysis import transformations
import numpy as np



def move_prot(ld_file, prot_file, z_dist, outname):

	u = mda.Universe(ld_file)
	halfz = u.dimensions[2] / 2

	UP = u.select_atoms('name P and prop z > %f' %halfz)
	UPz = UP.center_of_mass()[2]
	LP = u.select_atoms('name P and prop z < %f' %halfz)
	LPz = LP.center_of_mass()[2]

	### READ A PROTEIN-ONLY STRUCTURE
	prot = mda.Universe(prot_file).select_atoms('protein')
	prot.positions -= prot.center_of_mass()
	prot.positions += np.array([u.dimensions[0]/2, u.dimensions[1]/2, LPz - z_dist])


	### MERGE membrane and  protein
	newu = mda.Merge(prot, u.atoms)
	newu.dimensions = u.dimensions

	ag = newu.select_atoms('(protein) or (resname POPC DOPE SAPI TRIO) or (resname TIP3 and not byres (resname TIP3 and name OH2 and around 2.8 protein)) or (resname SOD CLA and not byres (resname SOD CLA and around 2.8 protein))')

	return ag.write(outname)
