import SnapPea
from mcomplex import *


def Mcomplex_from_data(SnapPeaTriangulation):
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


class manifold_list:
    """
    A manifold_list is a subscriptable object containing Mcomplexes
    which represent manifolds obtained from SnapPea triangulations.
    """
    def __init__(self, census):
        self.census = census

    def __repr__(self):
        return self.census.__repr__()

    def __len__(self):
        return len(self.census)

    def __getitem__(self, i):
        if not (0 <= i < len(self.census)):
            raise IndexError, "Requested manifold is out of range."
        manifold = self.census[i]
        return Mcomplex_from_data(manifold)    

closed_orientable = manifold_list( SnapPea.ClosedCensus())
five_tet_cusped = manifold_list( SnapPea.CuspedCensus())
six_tet_cusped_orientable = manifold_list( SnapPea.CuspedCensus(6)) 
six_tet_cusped_nonorientable = manifold_list( SnapPea.CuspedCensus(6,'n')) 
seven_tet_cusped_orientable = manifold_list( SnapPea.CuspedCensus(7)) 
seven_tet_cusped_nonorientable = manifold_list( SnapPea.CuspedCensus(7,'n')) 
alternating_knot_ext = manifold_list(SnapPea.AlternatingKnotExteriors)
nonalternating_knot_ext = manifold_list(SnapPea.NonalternatingKnotExteriors)

def get_mcomplex(name):
    return Mcomplex_from_data(SnapPea.get_manifold(name))


__all__ = ('closed_orientable',
           'five_tet_cusped',
           'six_tet_cusped_orientable',
           'six_tet_cusped_nonorientable',
           'seven_tet_cusped_orientable',
           'seven_tet_cusped_nonorientable',
           'alternating_knot_ext',
           'nonalternating_knot_ext',
           'get_mcomplex')

