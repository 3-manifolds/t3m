#$Id$
#   t3m - software for studying triangulated 3-manifolds
#   Copyright (C) 2002 Marc Culler, Nathan Dunfield and others
#
#   This program is distributed under the terms of the 
#   GNU General Public License, version 2 or later, as published by
#   the Free Software Foundation.  See the file GPL.txt for details.

import SnapPea
from mcomplex import *
import types


def Mcomplex_from_SnapPea(SnapPeaTriangulation):
    """
    Takes a list where the ith element represents the gluing data
    for the ith tetrahedron:

    ( [Neighbors], [Glueings] )

    and creates the corresponding Mcomplex.
    """
    M = SnapPeaTriangulation
    fill = 1 in [M.cusp_is_fillable(i) for i in range(M.get_num_cusps())]
    gluing_data = SnapPeaTriangulation.get_gluing_data(fill)
    num_tets = len(gluing_data)
    tets = map(lambda x: Tetrahedron(), range(num_tets))
    for i in range(num_tets):
        neighbors, perms = gluing_data[i]
        for k in range(4):
            tets[i].attach(TwoSubsimplices[k], tets[neighbors[k]], perms[k])
    return Mcomplex(tets)


class Manifold_list:
    """
    A Manifold_list is a subscriptable object containing Mcomplexes
    which represent manifolds obtained from SnapPea triangulations.
    """
    def __init__(self, census):
        self.census = census

    def __repr__(self):
        return self.census.__repr__()

    def __len__(self):
        return len(self.census)

    def __getitem__(self, i):
        if type(i) == types.SliceType:
            return Manifold_list(self.census[i])
        else:
            manifold = self.census[i]
            return Mcomplex_from_SnapPea(manifold)    

closed_orientable = Manifold_list( SnapPea.ClosedOrientable)
closed_nonorientable = Manifold_list( SnapPea.ClosedNonorientable)
five_tet_cusped = Manifold_list( SnapPea.Cusped5tet)
six_tet_cusped_orientable = Manifold_list( SnapPea.Cusped6tetOrientable) 
six_tet_cusped_nonorientable = Manifold_list( SnapPea.Cusped6tetNonorientable) 
seven_tet_cusped_orientable = Manifold_list( SnapPea.Cusped7tetOrientable) 
seven_tet_cusped_nonorientable = Manifold_list( SnapPea.Cusped7tetOrientable) 
alternating_knot_ext = Manifold_list(SnapPea.AlternatingKnotExteriors)
nonalternating_knot_ext = Manifold_list(SnapPea.NonalternatingKnotExteriors)

def get_mcomplex(name):
    return Mcomplex_from_SnapPea(SnapPea.get_manifold(name))


__all__ = ('closed_orientable',
           'closed_nonorientable',
           'five_tet_cusped',
           'six_tet_cusped_orientable',
           'six_tet_cusped_nonorientable',
           'seven_tet_cusped_orientable',
           'seven_tet_cusped_nonorientable',
           'alternating_knot_ext',
           'nonalternating_knot_ext',
           'get_mcomplex',
           'Mcomplex_from_SnapPea')

