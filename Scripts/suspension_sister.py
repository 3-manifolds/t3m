from snappea_conversion import *

def  make_test_manifold(shift):
    # replace the five valent edge e3 by all possible
    # suspensions of polygons.  Save the five resulting manifolds
    
    M = SnapPea_to_Mcomplex("m004(1,2)")
    edge = M.Edges[3]
    a = edge.get_arrow()
    t, b = M.suspension_of_polygon(edge.valence())
    t = t[shift : ] + t[ : shift ]
    b = b[shift : ] + b[ : shift ]
    M.replace_star(a, t, b)

    # Might as well do 3->2 moves as these never hurt...

    M.easy_simplify()
    M.rebuild()
    print  len(M), " tets"
    Mcomplex_to_SnapPea(M, "test_mfld%d" % shift)
    return(M)

for i in range(5):
    print i
    print
    M = make_test_manifold(i)
    M.erase()
    print

