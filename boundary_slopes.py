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

# Takes a Surface and computes the shifts along the edges
# of each tetrahedra.  In terms of the gluing equation variety
# picture, this is the order of zero of the corresponding shape
# parameter.

def add_shifts(S):
    shifts = []
    for i in range(len(S.Manifold)):
        shifts += [ S.Coefficients[i] * w for w in QuadShifts[S.Quadtypes[i]]]
    S.Shifts = shifts
    
def add_boundary_slope(surface, cusp_equations):
    assert len(cusp_equations) == 1

    surface.BoundarySlope = (-dot_product(surface.Shifts, cusp_equations[0][1]),
                              dot_product(surface.Shifts, cusp_equations[0][0]) )
    
def spun_surfaces_from_SnapPea(M):
    assert M.get_num_cusps() == 1
    
    N = t3m.Mcomplex_from_SnapPea(M)
    N.find_normal_surfaces()

    cusp_eqn = M.get_cusp_equations()
    for S in N.NormalSurfaces:
        add_shifts(S)
        add_boundary_slope(S, cusp_eqn)

    return N

# takes the boundary slopes of the given surfaces,
# normalizes them, and removes duplicates.

def boundary_slopes(surfaces):
    slopes = unique(map(normalize_slope,[surface.BoundarySlope for surface in surfaces]))
    if (0,0) in slopes:
        slopes.remove( (0,0))
    return  slopes

# The trick used is that a real incompressible surface can be spun
# both ways.

def boundary_slopes_with_trick(surfaces):
    directed_slopes = unique(map(weak_normalize_slope,[surface.BoundarySlope for surface in surfaces]))
    normalized_directed_slopes = map(normalize_slope, directed_slopes)
    slopes = boundary_slopes(surfaces)
    return [slope for slope in slopes if normalized_directed_slopes.count(slope) == 2]


def boundary_slopes_from_SnapPea(M):
    N = spun_surfaces_from_SnapPea(M)
    bs = boundary_slopes_with_trick(N.NormalSurfaces)
    return (N, bs)

def super_boundary_slopes_from_SnapPea(M):
    tris = find_distinct_triangulations(M)
    many_slopes = [ boundary_slopes_from_SnapPea(T)[1] for T in tris]
    final_slopes = []
    for s in many_slopes[0]:
        if not 0 in [ (s in slopes)  for slopes in many_slopes]:
            final_slopes.append(s)
    print len(many_slopes[0]) - len(final_slopes)
    return final_slopes
    
        

    
def (census, file_name):
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


    

    
    

