import MDAnalysis as mda
import numpy as np

u = mda.Universe('bilayer.gro')
halfz = u.dimensions[2]/2 ## Get the half of the z-dimension
UP = u.select_atoms('name P and prop z > %f' %halfz) ## Select the phosphorous atoms of the upper leaflet
LP = u.select_atoms('name P and prop z < %f' %halfz) ## Select the phosphorous atoms of the lower leaflet

ag_move = UP.residues.atoms + u.select_atoms('resname TIP3 SOD CLA POT and prop z > %f' %halfz)
## Create a new atomgroup with
## UP.residues.atoms select all the atoms of the upper leaflet
## the latter selects water/ion atoms above the halfz

dz = 80
ag_move.positions += np.array([0, 0, dz])
u.atoms.write('bilayer_moved.gro')