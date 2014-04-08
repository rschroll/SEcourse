import os
import time
import fcntl
from subprocess import Popen, PIPE
from StringIO import StringIO
from numpy import loadtxt, zeros, r_, sqrt, cross, dot, fft, argsort, matrix, identity, outer
from numpy.linalg import norm, LinAlgError
from pylab import find
import plot3D #plotMix
try:
    import cPickle as pickle
except ImportError:
    import pickle

plotting = plot3D # By default, but change this on the fly
__version__ = 2

class BaseSurf(object):
    """A base class for abstracting surfaces.  This contains much of what's
    needed for plotting and calculating energy densities.  Descendent classes
    are responsible for providing a load() method which, given a filename,
    sets the vertices, facets, and energies.  The load method is expected to
    create the following attributes:

        verts   N x 3 array containing the Cartesian coordinates of the vertices

        facets  M x 3 array containing the vertices making up each facet

        stretch M x 1 array of the elastic energy for each facet

        f_bend  M x 1 array of the bending energy for each facet

        f_z     M x 1 array of the average out-of-plane height of each facet

    Optionally, it can also create the following attributes:

        bend    N x 1 array containing the bending energy per vertex

        ref_coords  N x 3 array containing the strain-free configuration of the
                    vertices

    Also, f_area should be a M x 1 array of the facet areas, but if things are
    normal, this is automatically calculated.

    Other attributes can also be created, but BaseSurf() won't do anything with
    them.
    """

    def __new__(cls, filename=None, force_load=False, *args, **kw):
        if filename and not force_load:
            pfn = cls.pickle_name(filename)
            try:
                if os.path.getmtime(pfn) < os.path.getmtime(filename):
                    raise OSError
            except OSError:
                pass
            else:
                try:
                    obj = pickle.load(file(pfn, 'r'))
                except (IOError, pickle.UnpicklingError) as e:
                    pass
                else:
                    if obj.__version__ == __version__:
                        return obj
        #return super(BaseSurf, klass).__new__(cls, filename, force_load, *args, **kw)
        return object.__new__(cls) # Above was recommended, but this works the same

    def __init__(self, filename=None, force_load=False, mirror=False, flip_y=False, flip_z=False):
        """Args:    filename    file containing the surface definition

                    force_load  Set to True to force reloading from the
                                original file, instead of the pickled copy.

                    mirror      Whether the surface should be mirrored about
                                x=0.  This becomes an attribute which can be
                                adjusted at anytime.  Defaults to False.

                    flip_y      Whether the y-coordinate should be inverted
                                for display.  Defaults to False.

                    flip_z      Whether the z-coordinate should be inverted
                                for display.  Defaults to False.
        """

        self.__version__ = __version__
        if filename:
            if not hasattr(self, 'verts'):
                self.load(filename)
                self.pickle_save(filename)
        self.mirror = mirror
        if flip_y:
            self.flip_y()
        if flip_z:
            self.flip_z()

    def load(self, filename):
        raise NotImplementedError, "Descendants must implement a load method!"

    def pickle_save(self, filename):
        try:
            pickle.dump(self, file(self.pickle_name(filename), 'w'))
        except (IOError, pickle.PicklingError) as e:
            pass

    @classmethod
    def pickle_name(cls, filename):
        dir, name = os.path.split(filename)
        return os.path.join(dir, '.'+name+'.pkd')

    def calc_area(self):
        """This function is deprecated.  Just used f_area directly."""
        pass

    @property
    def f_area(self):
        if not hasattr(self, '_f_area'):
            self._f_area = r_[[norm(cross(self.verts[f[1],:]-self.verts[f[0],:],
                                          self.verts[f[2],:]-self.verts[f[0],:]))
                               for f in self.facets]] * 0.5
        return self._f_area

    @property
    def bend_density(self):
        return self.f_bend / self.f_area

    @property
    def stretch_density(self):
        return self.stretch / self.f_area

    def power_spectrum(self, slicex=0.5, gridsize=0.005, rot=False):
        """rot True for "rotated" triangulation, where first and last
        rows have twice as many vertices."""
        if self.gridsize:
            gridsize = self.gridsize
        nslices = len(find(self.ref_coords[:,0] == slicex))
        nfreq = len(find(self.ref_coords[:,1] == 0))
        if rot:
            nfreq = (nfreq + 1)/2
        freq = zeros((nslices,nfreq))
        for i in range(nslices):
            ief = find(abs(self.ref_coords[:,1] - self.length*i/(nslices-1)) < gridsize/3)
            a = argsort(self.verts[ief,0])
            if len(a)%2:
                z = r_[self.verts[ief[a],2][:-1], -self.verts[ief[a],2][:-1]]
            else:
                z = r_[self.verts[ief[a],2][1:-1], -self.verts[ief[a],2][1:-1]]
            if rot and i in (0, nslices-1):
                freq[i,:] = abs(fft.rfft(z))[:nfreq]/2 # b/c normalization
            else:
                freq[i,:] = abs(fft.rfft(z))
        return r_[1:nfreq:2]**4 * freq[:,1::2]**2

    def orient_facets(self):
        """Reorder facets such that they have consistent orientation wrt
        the z axis."""
        for i,f in enumerate(self.facets):
            v1 = self.verts[f[1]] - self.verts[f[0]]
            v2 = self.verts[f[2]] - self.verts[f[0]]
            if v1[0]*v2[1] - v1[1]*v2[0] > 0:
                self.facets[i] = [f[0], f[2], f[1]]


    def plot_stretch(self, **kw):
        self.plot_facets(self.stretch_density, "Stretching", **kw)

    def plot_bend(self, **kw):
        self.plot_facets(self.bend_density, "Bending", **kw)

    def plot_height(self, **kw):
        self.plot_facets(self.f_z, "Height", **kw)

    def plot_all(self, title=None, **kw):
        if not hasattr(plotting, "plt"):
            print "Error: plot_all() requires surface.plotting to have matplotlib available."
            print "Try setting surface.plotting = plot3D"
            return
        p = plotting.plt
        p.figure()
        p.subplot(131)
        self.plot_height(colorbar=True, **kw)
        p.subplot(132)
        self.plot_stretch(colorbar=True, **kw)
        p.subplot(133)
        self.plot_bend(colorbar=True, **kw)
        p.subplot(132)
        if title:
            p.title(title)

    def plot_profile(self, slices, axis, fig=None, gridsize=0.005, slicesize=3., xscale=1, yscale=1):
        """Plot a profile of the height z at given slices.

        Inputs: slices      An iterable containing the position at which to
                            plot the profile.

                axis        The axis along which to take the slices:
                                0 - across ridgeline
                                1 - along ridgeline

                fig         The figure on which the plot should be plotted.  If
                            None (default), it is plotted on the current figure.
                            If fig is False in the Boolean sense, a new figure
                            is created.  Otherwise, the figure fig is used.

                gridsize    The spacing of gridpoints.  If self.gridsize is
                            available, this quantity is ignored.

                slicesize   The inverse of the fraction of the gridsize to use
                            as a cutoff for determining of a point is in a slice.

                xscale      Rescale the x-axis while plotting.

                yscale      Rescale the y-axis while plotting.
        """

        if not hasattr(plotting, "plt"):
            print "Error: plot_profile() requires surface.plotting to have matplotlib available."
            print "Try setting surface.plotting = plot3D"
            return
        p = plotting.plt

        if fig is not None:
            if not fig:
                p.figure()
            else:
                p.figure(fig)
        if self.gridsize:
            gridsize = self.gridsize
        for slise in slices:
            inds = find(abs(self.verts[:,abs(axis-1)] - slise) < gridsize/slicesize)
            p.plot(self.verts[inds,axis]/xscale,self.verts[inds,2]/yscale,'.')

    def plot_facets(self, value, label='', ref_coords=False, threeD=False, colorbar=False, **kw):
        """Plot the facets of the surface, colored by value.  Must have
        len(value) == len(facets); otherwise fail."""

        if ref_coords:
            verts = self.ref_coords
        else:
            verts = self.verts

        if len(value) == len(self.verts) and hasattr(plotting, 'plot_vertmesh'):
            plotting.plot_vertmesh(verts, self.facets, value, **kw)
            return

        if len(value) != len(self.facets):
            raise NotImplementedError, "value must be length of facets."

        kw['mirror'] = self.mirror
        if threeD:
            pltfunc = plotting.plot_trimesh
        elif colorbar:
            pltfunc = plotting.plot_trimesh2D_cb
            kw["label"] = label
        else:
            pltfunc = plotting.plot_trimesh2D
        pltfunc(verts, self.facets, value, **kw)

    def plot_points(self, value, ref_coords=False, threeD=False, **kw):
        """Scatter plot value on the surface.  If len(value) == len(verts), plot
        the points at the vertices; if len(value) == len(facets), plot them at
        the centroids of the faces; otherwise fail."""

        inds = range(threeD and 3 or 2)
        if len(value) == len(self.verts):
            if ref_coords:
                pts = tuple(self.ref_coords[:,i] for i in inds)
            else:
                pts = tuple(self.verts[:,i] for i in inds)
        elif len(value) == len(self.facets):
            if ref_coords:
                pts = tuple(self.ref_centroids[:,i] for i in inds)
            else:
                pts = tuple(self.centroids[:,i] for i in inds)
        else:
            raise NotImplementedError, "value must be length of vertices or facets."

        skw = {'edgecolor': 'none', 's': 5, 'c': value}
        skw.update(kw)
        plotting.scatter(*pts, **skw)

    def _get_centroids(self):
        """This function is deprecated.  Just use centroids."""
        pass

    @property
    def centroids(self):
        """The center of each facet."""
        if not hasattr(self, '_centroids'):
            self._centroids = r_[[sum([self.verts[v,:] for v in face],0)/3 for face in self.facets]]
        return self._centroids

    @property
    def ref_centroids(self):
        if not hasattr(self, '_ref_centroids'):
            self._ref_centroids = r_[[sum([self.ref_coords[v,:] for v in face],0)/3 for face in self.facets]]
        return self._ref_centroids

    def flip_y(self):
        """Invert the vertical coordinate: y <-> -y."""
        self.verts *= r_[1,-1,1]
        if hasattr(self, "ref_coords"):
            self.ref_coords *= r_[1,-1,1]

    def flip_z(self):
        """Invert the out-of-plane coordinate: z <-> -z."""
        self.verts *= r_[1,1,-1]
        self.f_z *= -1
        if hasattr(self, "ref_coords"):
            self.ref_coords *= r_[1,1,-1]

class EvolveSurf(BaseSurf):
    """Abstracting a surface from Surface Evolver."""

    def load(self, filename):
        d = loadtxt(filename)
        self.nv, self.nf = map(int, d[0,:2])
        self.thickness, self.gridsize, self.length, self.delta, self.apptension = d[0,2:7]
        self.verts = d[1:self.nv+1,:3]
        self.tension = None
        if d.shape[1] in (8,9):
            self.gcurv = True
            self.hbend = d[1:self.nv+1, 3]
            self.kbend = d[1:self.nv+1, 4]
            self.bend = self.hbend + self.kbend
            if d.shape[1] == 9:
                self.tension = d[1:self.nv+1, 5]
                i = 6
            else:
                i = 5
        elif d.shape[1] == 7:
            self.gcurv = False
            self.bend = d[1:self.nv+1, 3]
            i = 4
        else:
            raise NotImplementedError, "Data file must have 7, 8, or 9 columns"
        self.ref_coords = d[1:self.nv+1, i:]
        self.facets = d[self.nv+1:,:3] - 1 # SE is 1-indexed
        self.stretch = d[self.nv+1:, 3]
        self.form_factors = d[self.nv+1:, 4:7]
        self.f_z = r_[[sum([self.verts[ind,2] for ind in face])/3
                        for face in self.facets]]

        self.v_area = zeros(self.nv)
        for i in range(self.nf):
            for j in range(3):
                self.v_area[self.facets[i][j]] += self.f_area[i]/3
        self.f_bend_density = r_[[sum([self.bend[ind]/self.v_area[ind]
                                        for ind in face])/3
                                   for face in self.facets]]
        self.f_bend = self.f_bend_density * self.f_area
        if self.gcurv:
            self.f_hbend_density = r_[[sum([self.hbend[ind]/self.v_area[ind]
                                             for ind in face])/3
                                        for face in self.facets]]
            self.f_kbend_density = r_[[sum([self.kbend[ind]/self.v_area[ind]
                                             for ind in face])/3
                                        for face in self.facets]]

        if any(self.ref_coords[..., 2]):
            print "Cannot calculate strain when ref_coords do not lie in XY plane. (Sorry.)"
            self.strain = None
        else:
            self.strain = zeros((self.nf, 2, 2))
            for fn in range(self.nf):
                vs = self.facets[fn]
                S = matrix([self.form_factors[fn][:2], self.form_factors[fn][1:]])
                v1 = self.verts[vs[1]] - self.verts[vs[0]]
                v2 = self.verts[vs[2]] - self.verts[vs[0]]
                F = matrix([[dot(v1,v1), dot(v1,v2)], [dot(v1,v2), dot(v2,v2)]])
                # A is Lambda.T in my notation
                A = matrix([[self.ref_coords[vs[1+i]][j] - self.ref_coords[vs[0]][j]
                            for j in range(2)] for i in range(2)])
                try:
                    self.strain[fn, ...] = A.I * 0.5*(F * S.I - identity(2)) * A
                except LinAlgError, e:
                    print "Problem calculating strain:", e

    @property
    def bend_density(self):
        return self.f_bend_density

    def stress(self, pr):
        """Calculate the stress from the strain, given the Poisson ratio."""
        tre = self.strain[:,0,0] + self.strain[:,1,1]
        return self.strain / (1 + pr) + pr / (1 - pr**2) * outer(tre, r_[1,0,0,1]).reshape(-1,2,2)

def read_when_ready(p):
    p.stdin.write('print "DONE READING"\n')
    res = ''
    while not res.count('DONE READING\nEnter command: '):
        try:
            out = p.stdout.read()
        except IOError:
            time.sleep(0.1)
        else:
            if not out:
                if p.poll() is not None:
                    raise IOError, "Evolver hung up unexpectedly:\n" + p.stderr.read()
            res += out
    return res.replace('DONE READING\nEnter command: ', '')

OUTPUTSURF = r"""
printf "%g %g %g %g %g %g %g 0 0\n",
    count(vertices,1), count(facets,1), thicknezz, gridsize, lngth, delta, exttension;
foreach vertices do
    printf "%0.15g %0.15g %0.15g %0.15g %0.15g %0.15g %0.15g %0.15g %0.15g\n",
            x, y, z, bend+freeedgebend, gbend, apptension, ref_coord[1],  ref_coord[2],  ref_coord[3];
foreach facets ff do
    printf "%g %g %g %0.15g %0.15g %0.15g %0.15g 0 0\n",
        ff.vertex[1].id, ff.vertex[2].id, ff.vertex[3].id, stretch,
        form_factors[1],  form_factors[2],  form_factors[3];
"""

class EvolveSurfD(EvolveSurf):
    """Loading from an Evolver dump file."""

    def load(self, filename, outputcode="outputsurf\n"): #OUTPUTSURF):
        p = Popen('evolver -w', shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE, close_fds=True)
        p.stdin.flush()
        p.stdin.write('\npid ::= 0\nload "%s"\n' % filename)
        fcntl.fcntl(p.stdout, fcntl.F_SETFL, os.O_NONBLOCK)
        fcntl.fcntl(p.stderr, fcntl.F_SETFL, os.O_NONBLOCK)
        output = read_when_ready(p)

        p.stdin.write(outputcode)
        output = read_when_ready(p).replace('Enter command: ', '').replace('more> ', '')
        p.stdin.close()

        EvolveSurf.load(self, StringIO(output))

class MinSurf(BaseSurf):
    """Abstracting a surface from the membrane minimizer code."""

    def load(self, filename):
        self.verts = loadtxt(filename)
        self.nv = len(self.verts)
        self.facets = loadtxt(os.path.join(os.path.split(filename)[0],
                                            "XFacetM.txt"), dtype=int) - 1
        self.nf = len(self.facets)
        self.calc_stretch()
        self.f_z = r_[[sum([self.verts[ind,2] for ind in face])/3
                        for face in self.facets]]
        self.calc_bend(filename)

    def calc_stretch(self):
        """Calculate the stretching energy per facet, using the algorithm in
        MEMstretchplot.m, from Eleni.  I hope it's correct."""

        nb12 = nb23 = nb13 = 0.1    # Default length of springs?

        ba12 = self.verts[self.facets[:,1],:] - self.verts[self.facets[:,0],:]
        ba23 = self.verts[self.facets[:,2],:] - self.verts[self.facets[:,1],:]
        ba13 = self.verts[self.facets[:,0],:] - self.verts[self.facets[:,2],:]
        nba12 = sqrt(ba12[:,0]**2 + ba12[:,1]**2 + ba12[:,2]**2)
        nba23 = sqrt(ba23[:,0]**2 + ba23[:,1]**2 + ba23[:,2]**2)
        nba13 = sqrt(ba13[:,0]**2 + ba13[:,1]**2 + ba13[:,2]**2)

        self.stretch = ((nba12-nb12)**2 + (nba23-nb23)**2 + (nba13-nb12)**2)/3
        # What is 1e8 doing at beginning?
        # MEMstretchplot contains some fiddling to the lowest energy facets, but
        # I'm guessing that's just to fix the display.

    def calc_bend(self, filename):
        """Calculate the bending energy per facet, as per MEMstretchplot.m"""

        #list12 = r_[self.facets[:,0], self.facets[:,1]]
        list3 = self.facets[:,2]
        #list23 = r_[self.facets[:,1], self.facets[:,2]]
        list1 = self.facets[:,0]
        #list31 = r_[self.facets[:,2], self.facets[:,0]]
        list2 = self.facets[:,1]

        AdjacentriM = loadtxt(os.path.join(os.path.split(filename)[0],
                                            "XAdjacentriM.txt"), dtype=int) - 1

        facetlist = zeros((self.nf, 3), dtype=int) + len(AdjacentriM)

        for jf in range(self.nf):
            T1 = ((list1[jf] == AdjacentriM[:,0]) * (list2[jf] == AdjacentriM[:,1])) \
                    * (list3[jf] == AdjacentriM[:,2]) + \
                 ((list2[jf] == AdjacentriM[:,0]) * (list1[jf] == AdjacentriM[:,1])) \
                    * (list3[jf] == AdjacentriM[:,3])
            T2 = ((list2[jf] == AdjacentriM[:,0]) * (list3[jf] == AdjacentriM[:,1])) \
                    * (list1[jf] == AdjacentriM[:,2]) + \
                 ((list3[jf] == AdjacentriM[:,0]) * (list2[jf] == AdjacentriM[:,1])) \
                    * (list1[jf] == AdjacentriM[:,3])
            T3 = ((list3[jf] == AdjacentriM[:,0]) * (list1[jf] == AdjacentriM[:,1])) \
                    * (list2[jf] == AdjacentriM[:,2]) + \
                 ((list1[jf] == AdjacentriM[:,0]) * (list3[jf] == AdjacentriM[:,1])) \
                    * (list2[jf] == AdjacentriM[:,3])

            t1 = find(T1)
            t2 = find(T2)
            t3 = find(T3)
            if len(t1) > 1 or len(t2) > 1 or len(t3) > 1:
                print "Error near t1", t1, t2, t3

            if len(t1):
                facetlist[jf,0] = t1

            if len(t2):
                facetlist[jf,1] = t2

            if len(t3):
                facetlist[jf,2] = t3

        Costh0 = zeros(len(AdjacentriM))
        Sinth0 = zeros(len(AdjacentriM))
        Costh = zeros(len(AdjacentriM))
        Sinth = zeros(len(AdjacentriM))

        for ji in range(len(AdjacentriM)):
            n1,n2,n3,n4 = AdjacentriM[ji,:]
            av = self.verts[n1,:]
            bv = self.verts[n2,:]
            cv = self.verts[n3,:]
            cpv = self.verts[n4,:]

            ab = bv - av
            normab = norm(ab)
            ac = cv - av
            acp = cpv - av
            ccp = cpv - cv
            vec1 = cross(ab, ac)
            normv1 = norm(vec1)
            vec2 = cross(acp, ab)
            normv2 = norm(vec2)
            dot12 = dot(ab,ac)
            dot12 = dot(ab,acp)
            dot23 = dot(ac,acp)
            dot11 = dot(ab,ab)
            dot22 = dot(ac,ac)
            dot33 = dot(acp,acp)
            crossba = cross(bv, av)
            crosscpc = cross(cpv, cv)
            crossv1_v2_ab = dot(cross(vec1, vec2), ab)
            costh = dot(vec1, vec2) / (normv1*normv2)
            sinth = crossv1_v2_ab / (normv1*normv2*normab)

            Costh[ji] = costh
            Sinth[ji] = sinth

            Costh0[ji] = 1
            Sinth0[ji] = 0

        CosthA = r_[Costh, 0]
        Costh0A = r_[Costh0, 0]
        SinthA = r_[Sinth, 0]
        Sinth0A = r_[Sinth0, 0]

        self.f_bend = 1 - sum([CosthA[facetlist[:,i]] * Costh0A[facetlist[:,i]]
                               + SinthA[facetlist[:,i]] * Sinth0A[facetlist[:,i]]
                               for i in range(3)], 0)/3
