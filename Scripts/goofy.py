#$Id: goofy.py,v 1.1.1.1 2002/08/19 16:24:45 t3m Exp $
# Nathan's code for importing and exporting snappea files.
from mcomplex import *
import re

# Converts a SnapPea file to MComplex.  Doesn't really use all the
# structure of the SnapPea file as it relies only on the fact that the
# gluing data for the ith pair of tetrahedra is given by the ith pair
# of lines like:
#
#      2    5    1   34 
#   3120 0321 0132 0132

def SnapPea_to_Mcomplex(file_name):
    data = open(file_name).read()
    count = 0

    neighbors_match = "^\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s*$"
    perm_match = "\s*([0123]{4,4})\s+([0123]{4,4})\s+([0123]{4,4})\s+([0123]{4,4})\s*$"
    snappea_re = re.compile(neighbors_match + perm_match, re.MULTILINE)
    
    fake_tets =[]
    
    curr_poss = 0
    while 1:
        m = snappea_re.search(data, curr_poss)
        if not m:
            break
        else:
            neighbors = map(int, m.group(1,2,3,4))
            perms = []
            for perm in m.group(5,6,7,8):
                perm = map(int, [perm[3], perm[2], perm[1], perm[0]] )
                perms.append(perm)
            fake_tets.append( (neighbors, perms) )
            curr_poss = m.end(8)
        

    return Mcomplex_from_data(fake_tets)

#------------End function SnapPea to Mcomplex--------------------

# Takes a list where the ith element represents the glueing data
# for the ith tetraherda:
#
#  ( [Neighbors], [Glueings] )
#
# and creates the corresponding Mcomplex


def Mcomplex_from_data(fake_tets):
    num_tets = len(fake_tets)
    tets = map(lambda x: Tetrahedron(), range(num_tets))
    for i in range(num_tets):
        neighbors, perms = fake_tets[i]
        for k in range(4):
            tets[i].attach(TwoSubsimplices[k], tets[neighbors[k]], perms[k])

    return Mcomplex(tets)

#-----------End function Mcomplex_from_data--------------------

# Exports an MComplex in SnapPea 2.0 format.
# ASSUMES THAT THE MANIFOLD IS ORIENTABLE, CLOSED, AND THAT THE LINK OF
# ANY VERTEX HAS GENUS AT MOST ONE.

def Mcomplex_to_SnapPea(manifold, file_name ):
    out = open(file_name, "w").write
    out("% Triangulation\n\n" + file_name + "\nnot_attempted\noriented_manifold\nCS_unknown\n\n")
    # Make sure everything is in order
    manifold.rebuild()

    torus_cusps = 0
    for vertex in manifold.Vertices:
        g = vertex.link_genus()
        if g > 1:
            raise ValueError, "Link of vertex has genus more than 1."
        if g == 1:
            torus_cusps = torus_cusps + 1

    # All torus cusps are unfilled
    
    out("%d 0" % torus_cusps)
    for i in range(torus_cusps):
        out( "   torus   0.000000000000   0.000000000000\n" )

    out("\n")

    # The num of tetrahedra

    out("%d\n" % len(manifold))
    
    # Output the tetraheda themselves.
    
    for tet in manifold.Tetrahedra:
        for face in TwoSubsimplices:
            out("    %d" % manifold.Tetrahedra.index( tet.Neighbor[face]))
        out("\n")
        for face in TwoSubsimplices:
            out(" %d%d%d%d" % tet.Gluing[face].tuple())

        out("\n")
        for vert in ZeroSubsimplices:
            if tet.Class[vert].link_genus() == 1:
                out("0 ")
            else:
                out("-1 ")
        out("\n")
        for i in range(4):
            out("0 0 0 0  0 0 0 0   0 0 0 0   0 0 0 0\n")
        out("0.0 0.0\n\n")


