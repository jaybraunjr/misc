import MDAnalysis as mda
import numpy as np
from MDAnalysis.analysis.distances import distance_array, self_distance_array
from MDAnalysis.core.groups import AtomGroup

class RDF:
    """
    Compute the RDF
    r.bins, r.shell_vol, r.shell_area available
    >>> r = smda.RDF(nbins=100, limits=(0.0, 15.0)[A])

    1) For static atomic groups
    >>> rdf = r.run(ag1, ag2, D=2/3,
    ...             b=0[ns], e=1e10[ns], nblocks=1)
    >>> rdf[:,0] = bins [A]
    >>> rdf[:,1] = RDF (avg) [unitless]
    >>> rdf[:,2] = RDF (std) [unitless]

    2) For dynamic atomic groups
    >>> data = []
    >>> for ts in u.trajectory:
    ...    g1_pos = g1.positions[SELECT]
    ...    g2_pos = g2.positions[SELECT]
    ...    data.append(r.run2d_frame(g1_pos, g2_pos, u.dimensions))
    ...    data.append(r.run3d_frame(g1_pos, g2_pos, u.dimensions))
    >>> rdf = np.transpose([r.bins, np.average(data, axis=0)])
    >>> rdf[:,0] = bins [A]
    >>> rdf[:,1] = RDF  [unitless]
    """

    def __init__(self, nbins=100, limits=(0.0, 15.0)):
        """
        Set up a RDF calculation.
        Define the number of bins and limits.

        Parameter
        ---------
        nbins  = 100  [int]
        limits = (0.0, 15.0) [A]
        """

        self.rdf_settings = {'bins': nbins, 'range': limits}
        self.rmax = limits[1]

        _, edges = np.histogram([-1], **self.rdf_settings)
        self.bins  = 0.5 * (edges[1:] + edges[:-1])
        self.shell_vol  = 4.0/3.0 * np.pi * (np.power(edges[1:], 3) - np.power(edges[:-1], 3))
        self.shell_area = np.pi * (np.power(edges[1:], 2) - np.power(edges[:-1], 2))


    def run(self, ag1, ag2, D = None, bframe=0, eframe=-1, skip=1):
        """
        Run a RDF calculation for static atomic groups

        Parameter
        ---------
        ag1, ag2:   atomic groups
        D  = None   [int]
        bframe = 0  [int] starting frame
        eframe = -1 [int] ending frame
        skip = 1    [int] skip frame

        Output
        ------
        [:,0] = bins [A]
        [:,1] = RDF
        """

        assert isinstance(ag1, AtomGroup)
        assert isinstance(ag2, AtomGroup)
        u = ag1.universe

        rdfs = []
        for ts in u.trajectory[bframe:eframe:skip]:
            if D == 3:
                rdf = self.run3d_frame(ag1.positions, ag2.positions, u.dimensions)
            elif D == 2:
                rdf = self.run2d_frame(ag1.positions, ag2.positions, u.dimensions)

            elif D ==1:
            	rdf = self.run1d_frame(ag1.positions, ag2.positions, u.dimensions)




            else:
                assert 1==0, 'Specify D=2 or D=3'
            rdfs.append(rdf)

        return np.transpose([self.bins, np.average(rdfs, axis=0)])


    def run3d_frame(self, g1_pos, g2_pos, dimensions):
        """
        Run a RDF calculation for one frame.
        Great flexibility as it takes atomic positions

        Parameters
        ----------
        g1_pos = g1.positions [A]
        g2_pos = g2.positions [A]
        dimensions = u.dimensions

        Outputs
        -------
        RDF array [A]
        """

        N = len(g1_pos) * len(g2_pos)
        if N == 0:
            return np.zeros(len(self.bins))

        vol = dimensions[0] * dimensions[1] * dimensions[2]
        density = N / vol

        d = distance_array(g1_pos, g2_pos, box=dimensions)

        count = np.histogram(d[d!=0], **self.rdf_settings)[0]
        count = count.astype(np.float64)
        rdf = count / density / self.shell_vol
        return rdf


    def run2d_frame(self, g1_pos, g2_pos, dimensions):
        """
        Run a RDF calculation for one frame.
        Great flexibility as it takes atomic positions

        Parameters
        ----------
        g1_pos = g1.positions [A]
        g2_pos = g2.positions [A]
        dimensions = u.dimensions

        Outputs
        -------
        RDF array [A]
        """

        N = len(g1_pos) * len(g2_pos)
        if N == 0:
            return np.zeros(len(self.bins))

        g1_pos[:,2] = 0.0
        g2_pos[:,2] = 0.0

        area = dimensions[0] * dimensions[1]
        density = N / area

        d = distance_array(g1_pos, g2_pos, box=dimensions)

        count = np.histogram(d[d!=0], **self.rdf_settings)[0]
        count = count.astype(np.float64)
        rdf = count / density / self.shell_area
        return rdf

    def run1d_frame(self, g1_pos, g2_pos, dimensions):
        """
        Run a RDF calculation for one frame.
        Great flexibility as it takes atomic positions

        Parameters
        ----------
        g1_pos = g1.positions [A]
        g2_pos = g2.positions [A]
        dimensions = u.dimensions

        Outputs
        -------
        RDF array [A]
        """

        N = len(g1_pos) * len(g2_pos)
        if N == 0:
            return np.zeros(len(self.bins))

        area = dimensions[0] * dimensions[2]
        density = N / area

        g1_pos[:,0] = 0.0
        g2_pos[:,0] = 0.0
        g1_pos[:,1] = 0.0
        g2_pos[:,1] = 0.0

        d = distance_array(g1_pos, g2_pos, box=dimensions)

        count = np.histogram(d[d!=0], **self.rdf_settings)[0]
        count = count.astype(np.float64)
        rdf = count / density / self.shell_area
        return rdf


r  = RDF(nbins=100, limits=(0, 25))
u  = mda.Universe('300nsCHYO.tpr', '300nsCHYO.xtc')



#####################################
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

##################################3




### Atomic group for P atoms
C  = u.select_atoms('name C10')

### Atomic group for P atoms in the upper leaflet
#uP = u.select_atoms('name P and prop z > %f' %(u.dimensions[2]/2))

### Atomic group for P atoms in the lower leaflets
#lP = u.select_atoms('name P and prop z < %f' %(u.dimensions[2]/2))


### P-P 3D RDF
rdf = r.run(C, C, D=3)
np.savetxt('C10_3D.dat', rdf)

print('finished 3D')

### P-P 2D RDF
rdf = r.run(C, C, D=2)
np.savetxt('C10_2D.dat', rdf)

print('finished 2D')

### 1D RDF
rdf = r.run(c1, c2, D=1)
np.savetxt('ZX.dat', rdf)

print('finished 1D')

### uP-uP 2D RDF




#rdf = r.run(uP, uP, D=2)
#np.savetxt('test.dat', rdf)

### C-C 2D RDF
### Computing 2D RDF frame-by-frame
### This example doesn't require this, but
### there're some examples that dynamic selections
### are required. Useful for those cases.
data = []
for ts in u.trajectory:
    g1_pos = C.positions
    data.append(r.run2d_frame(g1_pos, g1_pos, u.dimensions))
rdf = np.transpose([r.bins, np.average(data, axis=0)])
np.savetxt('test_ts_2D.dat', rdf)

print('ts_2d310 finished')