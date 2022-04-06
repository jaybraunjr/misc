import MDAnalysis as mda
from MDAnalysis.analysis import rdf
from MDAnalysis.analysis.rdf import InterRDF
import numpy as np
import math
import pandas as pd
import matplotlib.pyplot as plt








u = mda.Universe('1.21.tpr', '1.21.xtc')





class BeadGroup(object):
    def __init__(self, groups):
        self._groups = groups

    def __len__(self):
        return len(self._groups)

    @property
    def positions(self):
        return np.array([g.center_of_mass() for g in self._groups], dtype=np.float32)

    @property
    def universe(self):
        return self._groups[0].universe


chyo = u.select_atoms('resname CHYO')
chyo_groups = list(chyo.groupby('resids').values())





c1 = BeadGroup(chyo_groups)
c2 = BeadGroup(chyo_groups)







rdf = rdf.InterRDF(c1, c2,
					density=True,
					range=(1,30),
                    exclusion_block=(1, 1))


rdf.run()


plt.plot(rdf.rdf)
plt.xlabel(r"$r$ (Å)")
plt.xlim(0,40)
plt.show()
plt.savefig('rdf_iso.png')

print('finished')
