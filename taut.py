import t3m, SnapPea
from FXrays import find_Xrays
from boundary_slopes import find_distinct_triangulations

EdgeToQuad = {t3m.E01:2, t3m.E02:1, t3m.E03:0, t3m.E12:0, t3m.E13:1, t3m.E23:2}

def get_face_index(self):
    return self.Tetrahedron.Class[self.Face].Index

t3m.Arrow.get_face_index = get_face_index


class Mcomplex_with_taut(t3m.Mcomplex):
    def __init__(self, triangulation):
        if triangulation.__class__ == SnapPea.Triangulation:
            manifold = triangulation
        else:
            manifold = SnapPea.get_manifold(triangulation)
        
        tets = t3m.build_tets_from_SnapPea(manifold, 0)
        t3m.Mcomplex.__init__(self, tets)
        self.Name = manifold.get_name()
        self.SnapPeaTriangulation = manifold
        self.TautStructures = []

    # See class below for more on an angle structure.  The equations
    # are nominally inhomongenious, with (sum angles around edge) = 2
    # pi and (sum 3 angles in triangle) = pi.  So we add an extra
    # dummy variable (essentially "pi") to make them homogeneous.

    def build_angle_matrix(self):
        n = len(self)

        gluing_equations = []
        for edge in self.Edges:
            eqn = [0,]*(3 * n)
            for corner in edge.Corners:
                i = corner.Tetrahedron.Index
                eqn[3*i + EdgeToQuad[corner.Subsimplex]] += 1
            gluing_equations.append(eqn)
                
        angle_eqns = []
        # these say that the angles 
        for gluing_eqn in gluing_equations:
            angle_eqns += (gluing_eqn + [-2])

        for i in range(n):
            tri_eqn = [0,]*(3*i) + [1,1,1] + [0,]*(3*(n - i - 1)) + [-1]
            angle_eqns += tri_eqn

        self.AngleMatrix = t3m.Matrix(len(gluing_equations) + n, 3*n + 1)
        self.AngleMatrix.matrix = t3m.Numeric.array( angle_eqns, 'i')

    def find_taut_structures(self):
        self.build_angle_matrix()
        coeff_list = find_Xrays(self.AngleMatrix.rows,
                                self.AngleMatrix.columns,
                                self.AngleMatrix.matrix, modp=0, filtering=1)
        self.TautStructures = [TautStructure(self, coeff_vector) for coeff_vector in coeff_list]

    def known_to_semifiber(self):
        return 1 in [T.CarriesSemifiber for T in self.TautStructures]


class TautStructure:
    def __init__(self, manifold, angle_vector):
        self.NumTet = len(manifold)
        self.NumFaces = 2 * self.NumTet
        self.PiQuads = [ list(angle_vector[3*i:3*(i+1)]).index(1) for i in range(len(manifold))]
        self.setup_weight_equations(manifold)
        self.Surfaces = []
        self.Carries = None
        self.CarriesSemifiber = None
        self.carries_surface()

    def angle_is_pi(self, tet_index, edge):
        return self.PiQuads[tet_index] == EdgeToQuad[edge]
    
    def setup_weight_equations(self, manifold):
        eqns = t3m.Matrix( len(manifold.Edges), len(manifold.Faces) )
        for i in range(len(manifold.Edges)):
            edge = manifold.Edges[i]
            a = edge.get_arrow()
            while not self.angle_is_pi(a.Tetrahedron.Index, a.Edge):
                a.next()

            branches = [ [], [] ]
            for branch in branches:
                while 1:
                    branch.append(a.get_face_index())
                    a.next()
                    if  self.angle_is_pi(a.Tetrahedron.Index, a.Edge):
                        break

            for t in range(2):
                for f in branches[t]:
                    eqns[i, f] += (-1)**t

        self.WeightEquations = eqns

    def carries_surface(self):
        sols = find_Xrays(self.WeightEquations.rows,
                          self.WeightEquations.columns,
                          self.WeightEquations.matrix, modp=0, filtering=0)
        self.Surfaces = sols

        self.Carries = len(sols) > 0
        if len(sols) == 0:
            self.CarriesSemifiber = 0
        else:
            def max_entry(i):
                return max( [ s[i] for s in sols])
            self.CarriesSemifiber = min( map( max_entry, range(self.NumFaces) ) )  > 0

        
    def __repr__(self):
        s = "<Taut: " + repr(self.PiQuads) + " "
        if self.CarriesSemifiber:
            return s + "s-fiber>"
        elif self.Carries:
            return s + "carries>"
        return s + " >"
    

#---------code for testing SnapPea------------------

def semifiber_SnapPea(M):
    N = Mcomplex_with_taut(M)
    N.find_taut_structures()
    return N.known_to_semifiber()
    
def super_semifiber_SnapPea(M, restarts=15, moves=40):
    if semifiber_SnapPea(M):
        return 1
    
    tris = [M.clone()]
    for i in range(restarts):
        N = M.clone()
        for j in range(moves):
            N.randomize()
            N.simplify()
            if not 1 in [ N.same_triangulation(P) for P in tris]:
                tris.append(N)
                if semifiber_SnapPea(N):
                    return 1

    return 0


def save_semifiber(census):
    f = open("taut_out", "a")
    for M in census:
        if super_semifiber_SnapPea(M):
            f.write(M.get_name() + "\tsemifibers\n")
        else:
            f.write(M.get_name() + "\tunknown\n")
        f.flush()
                            
