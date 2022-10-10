import MDAnalysis as mda
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def index(gro_file, traj_file):

	u = mda.Universe(gro_file, traj_file)
	prot = u.select_atoms('protein')
	system = u.select_atoms('resid 1-234467')
	com = prot.center_of_mass()

	halfz = u.dimensions[2]/2 ## Get the half of the z-dimension
	UP = u.select_atoms('name P and prop z > %f' %halfz)
	LP = u.select_atoms('name P and prop z < %f' %halfz)

	with mda.selections.gromacs.SelectionWriter('leaflet.ndx', mode='w') as ndx:
		ndx.write(LP,name='lower_leaflet')
		ndx.write(UP,name='upper_leaflet')
		ndx.write(prot,name='prot')
		ndx.write(system,name='all')