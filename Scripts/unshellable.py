from mcomplex import *

# Defines an unshellable triangulation of the 3-ball with all
# vertices on the bdry.    Taken from:
#
#  Ziegler, G. M.  Shelling polyhedral 3-balls and 4-polytopes.
#  Discrete Comput. Geom.  19 (1998) 159 -- 174

#  The following function takes a list of 4-tuples of numbers each of
#  which represents a tetrahedron with labeled vertices.  Assuming
#  that the vertices uniquely determine a tetrahedra, it constructs
#  the corresponding manifold.

def manifold_from_vertices( data ):
    n = len(data)
    T = map( lambda x : Tetrahedron() , range(n) )
    for i in range(n):
        for j in range(i, n):
            match = do_share_face(data[i], data[j])
            if match:
                T[i].attach(TwoSubsimplices[match[0]], T[j], match[1])
    return Mcomplex(T)

# Checks to see if   the two length 4  tuples  have  an overlap  of  3
# elements.  If they do, returns the [face, gluing permutation] for t0
# -> t1.  Otherwise returns None

def do_share_face( t0, t1):
    perm = [-1, ]*4
    for i in range(4):
        for j in range(4):
            if t0[i] == t1[j]:
                perm[i] = j

    if perm.count(-1) != 1:
        return None

    i = perm.index(-1)

    for j in range(4):
        if perm.count(j) == 0:
            perm[i] = j
            break

    return [i, perm]

# Test

def simple_test():
    data = [ [1, 2, 3, 4], [1, 2, 3, 6], [1, 2, 4, 7], [2, 3, 4, 8], [1, 3, 4, 9] ]
    M = manifold_from_vertices(data)
    M.info()
    print len(M)
    M.randomize()
    print len(M)

# The Unshellable Complex.  Initially it has 21 tetrahedra, but the
# program can simplify it down to 16 tets.  

def unshellable():
    data = [ [1, 2, 3, 4], [1, 2, 5, 6], [2, 3, 6, 7], [3, 4, 7, 8],
             [4, 1, 8, 5], [1, 5, 6, 9], [1, 6, 2, 9], [1, 2, 4, 9],
             [1, 4, 8, 9], [1, 8, 5, 9], [2, 5, 6, 0], [2, 6, 7, 0],
             [2, 7, 3, 0], [2, 3, 1, 0], [2, 1, 5, 0], [3, 6, 7, 8],
             [3, 2, 4, 8], [3, 2, 6, 8], [4, 5, 7, 8], [4, 1, 3, 7],
             [4, 1, 5, 7]   ]

    M = manifold_from_vertices(data)
    M.orient()
    #M.info()
    print "Initial num tet: ", len(M)
    #print M.EdgeValences
    #print len(M.Edges), len(M.Vertices)
#    for i in range(10):
#      M.randomize()
#      print "num tet: ",  len(M)
    return M

M = unshellable()
M.info()