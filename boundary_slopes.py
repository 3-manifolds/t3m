import SnapPea, t3m, sys
from FXrays import find_Xrays
from types import *
import random


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

def is_subset(X, Y):
    for x in X:
        if Y.count(x) == 0:
            return 0
    return 1

def union(X, Y):
    return unique(X + Y)

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

def find_distinct_triangulations(M, restarts=5, moves=40):
    tris = [M.clone()]
    for i in range(restarts):
        N = M.clone()
        for j in range(moves):
            N.randomize()
            if not 1 in [ N.same_triangulation(P) for P in tris]:
                tris.append(N)

    return tris
    
# Used for converting normal surface into tetrahedron edge-shift data.
# QuadShift[k] is the shifts induced by quad Qk3 along edges (E03, E13, E23).
# Note that this follows the convention that the order of the edges
# is the same as the order of the quads, and _not_ in order E01, E02, E03.
#
# Note: this also exists in surface.py and had better be identical there.

QuadShifts = ((0, 1, -1), (-1, 0, 1), (1, -1, 0))


def quad_equation_from_gluing(eqn):
    nor_eqn = []
    for i in range(len(eqn)/3):
        for k in range(3):
            nor_eqn.append( dot_product( eqn[3*i: 3*(i+1)] , QuadShifts[k]))

    return nor_eqn

#  SnapPea uses a diffent convention for the order of its quads than
#  t3m does.  In particular SnapPea is (Q23, Q13, Q03), the exact
#  reverse of t3m.

def convert_quads_SnapPea_to_t3m(eqn):
    assert len(eqn) % 3 == 0
    new_eqn = []
    for i in range(len(eqn)/ 3):
        new_eqn += [eqn[3*i + 2], eqn[3*i + 1], eqn[3*i]]
    return new_eqn
    
        

#---------------------------------------------------------------------
#
# Class: OneCuspedManifold
#
#--------------------------------------------------------------------

# This class is an amalgam of t3m.mcomplex and SnapPea.Triangulation
# It takes the latter to create the class, or alternatively a string
# which specifies the name of a manifold SnapPea knows about.  Class
# is setup so that calls to many SnapPea.Triangulation methods are
# mapped through.  So, e.g. OneCuspedManifold.volume() works

class OneCuspedManifold(t3m.Mcomplex):
    def __init__(self, triangulation):
        if triangulation.__class__ == SnapPea.Triangulation:
            manifold = triangulation
        else:
            manifold = SnapPea.get_manifold(triangulation)
        
        if not manifold.get_triangulation_is_orientable():
            raise ValueError, "Manifold must be orientable"
        if manifold.get_num_cusps() != 1:
            raise ValueError, "Manifold does not have one cusp"

        tets = t3m.build_tets_from_SnapPea(manifold, 0)
        t3m.Mcomplex.__init__(self, tets)
        self.SnapPeaTriangulation = manifold
        self.ClosedSurfaces = []
        self.AngleStructures = []

        # finally set cusp equations

        self.CuspEquations = map(convert_quads_SnapPea_to_t3m, manifold.get_cusp_equations()[0])
        
    # should really disable all t3m.Mcomplex methods which
    # change the triangulation as need to keep the two reps
    # in sync.

    # Takes a Surface and computes the shifts along the edges
    # of each tetrahedra.  In terms of the gluing equation variety
    # picture, this is the order of zero of the corresponding shape
    # parameter.

    def find_normal_surfaces(self, modp=0):
        t3m.Mcomplex.find_normal_surfaces(self, modp)
        for S in self.NormalSurfaces:
            S.add_shifts()
            S.add_boundary_slope(self.CuspEquations)

    # takes the boundary slopes of the given surfaces,
    # normalizes them, and removes duplicates.

    def boundary_slopes(self, surfaces=None):
        if surfaces == None:
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

    def incompressible_slopes(self):
        return self.boundary_slopes([S for S in self.NormalSurfaces if S.Incompressible == 1])
    
    def find_closed_normal_surfaces(self, modp=0):
        for surface in self.ClosedSurfaces:
            surface.erase()
            self.ClosedSurfaces.remove(surface)
        self.build_matrix()

        cusp_eqns = map(quad_equation_from_gluing, self.CuspEquations)
        new_eqns = self.QuadMatrix.matrix.tolist() + cusp_eqns[0] + cusp_eqns[1]
        coeff_list = find_Xrays(self.QuadMatrix.rows + 2,
                                        self.QuadMatrix.columns,
                                        new_eqns, modp)
        for coeff_vector in coeff_list:
            self.ClosedSurfaces.append(t3m.ClosedSurfaceInCusped(self, coeff_vector))

    def  all_closed_surfaces_edge_linking(self):
        for surface in self.ClosedSurfaces:
            if not surface.is_edge_linking_torus()[0]:
                return 0

        return 1

    def has_incompressible_closed_surface(self):
        return 1 in [S.Incompressible for S in self.ClosedSurfaces]

    def closed_normal_surface_info(self):
        out = sys.stdout
        for surface in self.ClosedSurfaces:
            out.write("-------------------------------------\n\n")
            surface.info(out)
            out.write('\n')

    # See class below for more on an angle structure.  The equations
    # are nominally inhomongenious, with (sum angles around edge) = 2
    # pi and (sum 3 angles in triangle) = pi.  So we add an extra
    # dummy variable (essentially "pi") to make them homogeneous.

    def build_angle_matrix(self):
        n = len(self)

        # should really compute these directly from
        # the t3m representation.
        
        gluing_equations = map(convert_quads_SnapPea_to_t3m,
                               self.SnapPeaTriangulation.get_gluing_equations())

        angle_eqns = []
        # these say that the angles 
        for gluing_eqn in gluing_equations:
            angle_eqns += (gluing_eqn + [-2])

        for i in range(n):
            tri_eqn = [0,]*(3*i) + [1,1,1] + [0,]*(3*(n - i - 1)) + [-1]
            angle_eqns += tri_eqn

        self.AngleMatrix = t3m.Matrix(len(gluing_equations) + n, 3*n + 1)
        self.AngleMatrix.matrix = t3m.Numeric.array( angle_eqns, 'i')
        
    def find_angle_structures(self):
        self.build_angle_matrix()
        coeff_list = find_Xrays(self.AngleMatrix.rows,
                                self.AngleMatrix.columns,
                                self.AngleMatrix.matrix, modp=0, filtering=0)
        self.AngleStructures = [AngleStructure(coeff_vector) for coeff_vector in coeff_list]

    def angle_structure_info(self):
        for structure in self.AngleStructures:
            print structure

    def mark_incompressible(self):
        for surfaces in (self.NormalSurfaces, self.ClosedSurfaces):
            for S in surfaces:
                for A in self.AngleStructures:
                    if A.supports_surface(S):
                        assert S.Incompressible != 0
                        S.Incompressible = 1
                        break

        for S in self.ClosedSurfaces:
            if S.is_edge_linking_torus()[0]:
                assert S.Incompressible == None
                S.Incompressible = 0


    def find_all(self):
        self.find_normal_surfaces()
        self.find_closed_normal_surfaces()
        self.find_angle_structures()
        self.mark_incompressible()
        
    
# Code for adding many SnapPea.Triangulation calls to OneCuspedManifold

meths = ['Dirichlet', 'core_geodesic',
             'current_fillings_become_meridians', 'cusp_is_fillable',
             'fill_cusp', 'fundamental_group',  'get_cusp_is_complete',
             'get_cusp_is_orientable', 'get_cusp_l', 'get_cusp_m',
             'get_cusp_moduli', 'get_cusp_modulus', 'get_cusp_shape',
             'get_cusp_shapes', 'get_drillable_curves', 'get_gluing_data',
             'get_name', 'get_num_cusps',
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


class AngleStructure:
    def __init__(self, solution_vector):
        self.Angles = solution_vector[:-1]
        self.Pi = solution_vector[-1]
        self.n = len(self.Angles)/3 

    def format_fraction(self, top):
        g = SnapPea.gcd(top, self.Pi)
        if g != self.Pi:
            return "%d/%d" % (top/g, self.Pi/g)
        else:
            return "%d" % (top/g)
    
    def __repr__(self):
        s =  "<AngleStr: "
        for i in range(self.n):
            s = s + "%s %s %s; " % tuple(map(self.format_fraction, self.Angles[3*i:3*(i + 1)]))
        return s[:-2]+" >"

    # returns true if the edges opposite the quads of the given surface
    # are 0
    
    def supports_surface(self, surface):
        for i in range(self.n):
            if self.Angles[3*i + surface.Quadtypes[i]] != 0:
                return 0

        return 1



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
    

#------------------------------------------------------------------------
#
#  Code for trying to find out as much about the manifold as possible
#
#------------------------------------------------------------------------

class OneCuspedData:
    def __init__(self, name):
        self.Name = name
        self.IncompSlopesSubset = None
        self.IncompSlopesSuperset = None
        self.SmallOrLarge = "?"   # should be one of "small", "?", or "large"

    def update_slopes_subset(self, slopes):
        if self.IncompSlopesSubset == None:
            self.IncompSlopesSubset = slopes
        else:
            self.IncompSlopesSubset = union(self.IncompSlopesSubset, slopes)

    def update_slopes_superset(self, slopes):
        if self.IncompSlopesSuperset == None:
            self.IncompSlopesSuperset = slopes
        else:
            self.IncompSlopesSuperset = [ s for s in self.IncompSlopesSuperset if s in slopes]

    def update_small_or_large(self, change):
        assert change in ("small", "?", "large")
        if change == "small":
            assert self.SmallOrLarge != "large"
        elif change == "large":
            assert self.SmallOrLarge != "small"

        self.SmallOrLarge = change

    def __repr__(self):
        return  self.Name + "\t" + repr(self.IncompSlopesSubset) + "\t" + repr(self.IncompSlopesSuperset) + "\t" + self.SmallOrLarge

def super_surfaces_in_SnapPea(M):
    data = OneCuspedData(M.get_name())
    tris = find_distinct_triangulations(M)

    for T in tris:
        N = OneCuspedManifold(T)
        N.find_all()

        data.update_slopes_subset(N.incompressible_slopes())
        data.update_slopes_superset(N.boundary_slopes_with_trick())
        if N.all_closed_surfaces_edge_linking():
            data.update_small_or_large("small")
        if N.has_incompressible_closed_surface():
            data.update_small_or_large("large")

    return data

            

    

#-------------------------------------------------------------------------
#
#  Code for looking for small one-cusped manifolds
#
#-------------------------------------------------------------------------

def is_small(M):
    N = OneCuspedManifold(M)
    N.find_closed_normal_surfaces()
    return not N.has_non_trivial_closed_surface()

def super_is_small(M):
    tris = find_distinct_triangulations(M)
    for T in tris:
        if is_small(T):
            return 1

    return 0
    
def small_go(file, census):
    f = open(file, "w")
    for M in census:
        if M.get_num_cusps() == 1:
            if super_is_small(M):
                f.write(M.get_name() + "\tsmall\n")
            else:
                f.write(M.get_name() + "\t?\n")
            f.flush()

#small_go("small2", SnapPea.Cusped5tetOrientable)

#-------------------------------------------------------------------------------
#
#  Code for testing to make sure conventions match
#
#-------------------------------------------------------------------------------

def normal_surface_equations_from_gluing(M):
    glueeqns =  M.get_gluing_equations()
    return [quad_equation_from_gluing(convert_quads_SnapPea_to_t3m(eqn))
            for eqn in glueeqns]
    
    
def check_equations(M):
    snapeqn = normal_surface_equations_from_gluing(M)
    snapeqn.sort()
    
    N = t3m.Mcomplex_from_SnapPea(M)
    N.build_matrix()
    t3meqn = N.QuadMatrix.to_list()
    t3meqn.sort()

    #print snapeqn
    #print
    #print t3meqn
    
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
        N = OneCuspedManifold(name)
        N.find_normal_surfaces()
        bs = N.boundary_slopes()
        print name
        if not is_subset(bs, os):
            print "WUGGA!!!!!"
            print os
            print bs
        sys.stdout.flush()

#check_against_old_mathematica_prog2()

# check SnapPea gluing conventions

def evaluate_equation(eqn, edge_params):
    z = 1
    for i in range(len(eqn)):
        z = z*(edge_params[i]**eqn[i])
    return z

                      
def test_gluing_equations(M):
    shapes = M.tet_shapes(1)
    edge_param = []
    for z in [shape["shape rect"] for shape in shapes]:
        edge_param += [z, 1/(1-z), (z-1)/z]

    for eqn in M.get_gluing_equations():
        print evaluate_equation(eqn, edge_param)

    print "cusp"

    for eqn in M.get_cusp_equations()[0]:
        print evaluate_equation(eqn, edge_param)
        
            
    
    

    
    

