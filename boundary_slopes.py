import SnapPea, t3m, sys


# helper functions

def dot_product(x,y):
    assert len(x) == len(y)
    dot = 0
    for i in range(len(x)):
        dot += x[i]*y[i]
    return dot

def unique(list):
    sorted_list = list[:]
    sorted_list.sort()
    
    ans = []
    for i in range(len(list) - 1):
        if sorted_list[i] != sorted_list[i + 1]:
            ans.append(sorted_list[i])

    ans += sorted_list[-1:]
    return ans

def is_contained_in(X, Y):
    for x in X:
        if Y.count(x) == 0:
            return 0
    return 1


def weak_normalize_slope( slope ):
    a, b = slope
    if a == b == 0:
        return slope
    g = SnapPea.gcd(a,b)
    a, b = a/g, b/g
    return (a,b)
    
def normalize_slope( slope ):
    a, b = weak_normalize_slope(slope)
    if a == b == 0:
        return slope
    if a < 0:
        a, b = -a, -b
    elif a == 0 and b < 0:
        b = -b
    return (a,b)

def find_distinct_triangulations(M, restarts=5, moves=30):
    tris = [M.clone()]
    for i in range(restarts):
        N = M.clone()
        for j in range(moves):
            N.randomize()
            if not 1 in [ N.same_triangulation(P) for P in tris]:
                tris.append(N)

    return tris
    
# renaming

QuadShifts = (t3m.Shift[t3m.E01], t3m.Shift[t3m.E02], t3m.Shift[t3m.E03])

#---------------------------------------------------------------------
#
# Class: OneCuspedManifold
#
#--------------------------------------------------------------------

# This class is an amalgam of t3m.mcomplex and SnapPea.Triangulation
# It takes the latter to create the class.  Class is setup so that
# calls to many SnapPea.Triangulation methods are mapped through.
# So, e.g. OneCuspedManifold.volume() works

class OneCuspedManifold(t3m.Mcomplex):
    def __init__(self, SnapPea_triangulation):
        if not SnapPea_triangulation.get_triangulation_is_orientable():
            raise ValueError, "Manifold must be orientable"
        if SnapPea_triangulation.get_num_cusps() != 1:
            raise ValueError, "Manifold has more than one cusp"

        tets = t3m.build_tets_from_SnapPea(SnapPea_triangulation, 0)
        t3m.Mcomplex.__init__(self, tets)
        self.SnapPeaTriangulation = SnapPea_triangulation
    
    # should really disable all t3m.Mcomplex methods which
    # change the triangulation as need to keep the two reps
    # in sync.

    # Takes a Surface and computes the shifts along the edges
    # of each tetrahedra.  In terms of the gluing equation variety
    # picture, this is the order of zero of the corresponding shape
    # parameter.

    def find_normal_surfaces(self, modp=0):
        t3m.Mcomplex.find_normal_surfaces(self, modp)
        cusp_eqns = self.get_cusp_equations()[0]
        for S in self.NormalSurfaces:
            S.add_shifts()
            S.add_boundary_slope(cusp_eqns)

    # takes the boundary slopes of the given surfaces,
    # normalizes them, and removes duplicates.

    def boundary_slopes(self):
        surfaces = self.NormalSurfaces
        slopes = unique(map(normalize_slope,[surface.BoundarySlope for surface in surfaces]))
        if (0,0) in slopes:
            slopes.remove( (0,0))
        return  slopes

    # The trick used is that a real incompressible surface can be spun
    # _both_ ways.

    def boundary_slopes_with_trick(self):
        surfaces = self.NormalSurfaces
        directed_slopes = unique(map(weak_normalize_slope,[surface.BoundarySlope for surface in surfaces]))
        normalized_directed_slopes = map(normalize_slope, directed_slopes)
        slopes = self.boundary_slopes()
        return [slope for slope in slopes if normalized_directed_slopes.count(slope) == 2]

# Code for adding many SnapPea.Triangulation calls to OneCuspedManifold

meths = ['Dirichlet', 'core_geodesic',
             'current_fillings_become_meridians', 'cusp_is_fillable',
             'fill_cusp', 'fundamental_group',
             'get_cusp_equation', 'get_cusp_equations', 'get_cusp_is_complete',
             'get_cusp_is_orientable', 'get_cusp_l', 'get_cusp_m',
             'get_cusp_moduli', 'get_cusp_modulus', 'get_cusp_shape',
             'get_cusp_shapes', 'get_drillable_curves', 'get_gluing_data',
             'get_gluing_equations', 'get_name', 'get_num_cusps',
             'get_num_tetrahedra', 'get_solution_type',
             'get_triangulation_is_orientable', 'homology',
             'is_canonical_triangulation',
             'remove_Dehn_fillings', 'set_cusp',
             'set_name', 'shortest_curves_become_meridians', 'symmetry_group',
             'tet_shapes', 'volume']

# there has _got_ to be a better way to do this...

for meth in meths:
    setattr(OneCuspedManifold, meth,
            eval( "lambda self, *rest : apply( getattr(self.SnapPeaTriangulation, '" + meth + "'), rest)"))


# used in next function

def boundary_slopes_from_SnapPea(M):
    N = OneCuspedManifold(M)
    N.find_normal_surfaces()
    return N.boundary_slopes_with_trick()

# tries several different triangulations to find smallest possible number
# of boundary slopes.

def super_boundary_slopes_from_SnapPea(M):
    tris = find_distinct_triangulations(M)
    many_slopes = [ boundary_slopes_from_SnapPea(T) for T in tris]
    final_slopes = []
    for s in many_slopes[0]:
        if not 0 in [ (s in slopes)  for slopes in many_slopes]:
            final_slopes.append(s)
    print len(many_slopes[0]) - len(final_slopes)
    return final_slopes
    

#-------------------------------------------------------------------------------
#
#  Code for testing to make sure conventions match
#
#-------------------------------------------------------------------------------

def quad_equation_from_gluing(eqn):
    nor_eqn = []
    for i in range(len(eqn)/3):
        for k in range(3):
            nor_eqn.append( dot_product( eqn[3*i: 3*(i+1)] , QuadShifts[k]))

    return nor_eqn
            
def normal_surface_equations_from_gluing(M):
    glueeqns =  M.get_gluing_equations()
    return [quad_equation_from_gluing(eqn) for eqn in glueeqns]
    
    
def check_equations(M):
    snapeqn = normal_surface_equations_from_gluing(M)
    snapeqn.sort()
    
    N = t3m.Mcomplex_from_SnapPea(M)
    N.build_matrix()
    t3meqn = N.QuadMatrix.to_list()
    t3meqn.sort()
    
    return snapeqn == t3meqn

#-------------------------------
#
# Code for comparing with old Mathematica version
#
#---------------------------
    
def save_slopes(census, file_name):
    f = open(file_name, "a")

    for M in census:
        if M.get_num_cusps() == 1:
            N = spun_surfaces_from_SnapPea(M)
            bs = tuple(boundary_slopes(N.NormalSurfaces))
            f.write(M.get_name() + "\t" + repr(bs) + "\n" )
            f.flush()
    
def check_against_old_mathematica_prog2(file_name="/Users/dunfield/work/work/old_projects/haken/manifold_data/old_normal_slopes"):
                                      
    f = open(file_name, "r")

    for line in f.xreadlines():
        name, slopes= line.split("\t")
        os = list(eval(slopes))
        os.sort()
        M = SnapPea.get_manifold(name)
        N = spun_surfaces_from_SnapPea(M)
        bs = boundary_slopes(N.NormalSurfaces)
        print name
        if not is_contained_in(bs, os):
            print "WUGGA!!!!!"
            print os
            print bs
        sys.stdout.flush()


    

    
    

