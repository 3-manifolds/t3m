#    SnapPea.py
#
#    The Python interface to the SnapPea kernel consists of two parts.
#
#        SnapPea.py (this file) defines a set of objects (Triangulation,
#            AbelianGroup, etc.) purely in Python.
#
#        SnapPeaC.c (in a different file) implements SnapPea.py's methods
#            as a set of wrappers to the standard SnapPea kernel functions,
#            which are written in C.

import SnapPeaC
import os
import string

False = 0
True  = 1


#    Keep track of when memory is released.
def VerifyMyMallocUsage():
    return SnapPeaC.verify_my_malloc_usage()
    

class Triangulation:
    
    def __init__(self, spec, orientable = 1, index = 0):
        #    Triangulation('file_name') reads a file from disk.
        #    Triangulation(num_tet (=5,6,7), orientable (=0(false),1(true)), index)
        #        reads a cusped census manifold.
        #    Triangulation((long) (Trianguation *)) accepts a pointer
        #        to an existing SnapPea Triangulation.
        if spec in range(5,8):
            self.triangulation = SnapPeaC.get_cusped_census_manifold(spec, orientable, index)
        elif type(spec) == type(''):
            self.triangulation = SnapPeaC.get_triangulation(spec)
        else:
            self.triangulation = spec
    
    def __del__(self):
        #    Let the SnapPea kernel free its private representation
        #    of a Triangulation.
        SnapPeaC.free_triangulation(self.triangulation)
    
    def __repr__(self):

        s = '\n'
        s = s + 'name:      %s' % SnapPeaC.get_triangulation_name(self.triangulation) + '\n'
        s = s + 'volume:    %f' % SnapPeaC.volume(self.triangulation)                 + '\n'
        s = s + 'homology:  %s' % self.homology()                                     + '\n'
        s = s + 'Dehn fillings:\n'
        for i in range(SnapPeaC.get_num_cusps(self.triangulation)):
            s = s + '%3i  ' % i
            if SnapPeaC.get_cusp_is_complete(self.triangulation, i):
                s = s + '    complete'
            else:
                s = s + '%f ' % SnapPeaC.get_cusp_m(self.triangulation, i)
                if SnapPeaC.get_cusp_is_orientable(self.triangulation, i):
                    s = s + '%f ' % SnapPeaC.get_cusp_l(self.triangulation, i)
            s = s + '\n'
        
        return s
    
    def save(self, file_name):
        if len(file_name) > 0:
            SnapPeaC.set_triangulation_name(self.triangulation, file_name)
        SnapPeaC.save_triangulation(self.triangulation, file_name)
    
    def clone(self):
        return Triangulation(SnapPeaC.copy_triangulation(self.triangulation))
    
    def get_name(self):
        return SnapPeaC.get_triangulation_name(self.triangulation)
    
    def set_name(self, name):
        SnapPeaC.set_triangulation_name(self.triangulation, name)
        
    def get_triangulation_is_orientable(self):
        return SnapPeaC.get_triangulation_is_orientable(self.triangulation)
    
    def get_solution_type(self):
        return SnapPeaC.get_solution_type(self.triangulation)
    
    def get_num_tetrahedra(self):
        return SnapPeaC.get_num_tetrahedra(self.triangulation)
    
    def get_num_cusps(self):
        return SnapPeaC.get_num_cusps(self.triangulation)
    
    def get_cusp_is_orientable(self, i):
        return SnapPeaC.get_cusp_is_orientable(self.triangulation, i)
    
    def get_cusp_is_complete(self, i):
        return SnapPeaC.get_cusp_is_complete(self.triangulation, i)
    
    def get_cusp_m(self, i):
        return SnapPeaC.get_cusp_m(self.triangulation, i)
    
    def get_cusp_l(self, i):
        return SnapPeaC.get_cusp_l(self.triangulation, i)
    
    def set_cusp(self, index, m = 0, l = 0, recompute = 1):
    
        #    set_cusp(i)
        #        makes cusp i complete.
        #
        #    set_cusp(i, m, l)
        #        does (m,l) Dehn filling on cusp i.
        #        (As a special case, set_cusp(i, 0, 0) is equivalent
        #        to set_cusp(i), and makes the cusp complete.)
        #
        #    set_cusp(i, m, l, 0)
        #        sets the cusp coefficients, but doesn't attempt
        #        to find the hyperbolic structure.
        #        (set_cusp(i, m, l, 1) is equivalent to set_cusp(i, m, l).)
        
        if index in range(SnapPeaC.get_num_cusps(self.triangulation)):
            SnapPeaC.set_cusp_info(self.triangulation, index, m, l, recompute)
        else:
            return 'cusps are numbered 0 through %i\n' % SnapPeaC.get_num_cusps(self.triangulation)
    
    def cusp_is_fillable(self, index):
        return SnapPeaC.cusp_is_fillable(self.triangulation, index)
    
    def remove_Dehn_fillings(self):
        return SnapPeaC.remove_Dehn_fillings(self.triangulation)
    
    def tet_shapes(self, use_fixed_coordinates):
        return SnapPeaC.tet_shapes(self.triangulation, use_fixed_coordinates)
    
    def volume(self):
        return SnapPeaC.volume(self.triangulation)
    
    def homology(self):
        #    If the homology couldn't be computed
        #    (perhaps because some Dehn filling coefficients
        #    aren't integers) return None.
        #    Otherwise return the corresponding AbelianGroup.
        #    Note:  An empty coefficient list is not None,
        #    and will generate an AbelianGroup with
        #    no torsion coefficients.
        theTorsionCoefficients = SnapPeaC.homology(self.triangulation)
        if theTorsionCoefficients == None:
            return None
        else:
            return AbelianGroup(theTorsionCoefficients)
    
    def fundamental_group(    self,
                            simplify_flag = 1,
                            fillings_affect_generators_flag = 1,
                            minimize_num_generators_flag = 0):
        return FundamentalGroup(SnapPeaC.fundamental_group(
                                    self.triangulation,
                                    simplify_flag,
                                    fillings_affect_generators_flag,
                                    minimize_num_generators_flag))

    def fill_cusp(self, index):
        if index not in range(self.get_num_cusps()):
            return 'cusps are numbered 0 through %i\n' % self.get_num_cusps()
        elif self.get_num_cusps() == 1:
            print 'The triangulation must retain at least one cusp.'
        elif not self.cusp_is_fillable(index):
            print 'To permanently fill a cusp, the Dehn filling coefficients must be relatively prime integers.'
        else:
            theFilledTriangulation = SnapPeaC.fill_cusp(self.triangulation, index)
            if (theFilledTriangulation != self.triangulation):
                SnapPeaC.free_triangulation(self.triangulation)
                self.triangulation = theFilledTriangulation
    
    def get_drillable_curves(self, max_segments = 6):
        return SnapPeaC.get_drillable_curves(self.triangulation, max_segments)
    
    def drill_curve(self, index, max_segments = 6):
        theDrilledTriangulation = SnapPeaC.drill_curve(self.triangulation, index, max_segments)
        if (theDrilledTriangulation != 0):
            SnapPeaC.free_triangulation(self.triangulation)
            self.triangulation = theDrilledTriangulation
    
    def get_normal_surfaces(self):
        return SnapPeaC.get_normal_surfaces(self.triangulation)
    
    def split_on_normal_surface(self, index):
        #    SnapPeaC produces bare SnapPea Triangulations.
        #    Convert them to Python Triangulation objects.
        theBareTriangulations = SnapPeaC.split_along_normal_surface(self.triangulation, index)
        theTriangulationObjects = []
        for theBareTriangulation in theBareTriangulations:
            theTriangulationObjects.append(Triangulation(theBareTriangulation))
        return theTriangulationObjects
    
    def core_geodesic(self, index):
        return SnapPeaC.core_geodesic(self.triangulation, index)
    
    def shortest_curves_become_meridians(self):
        SnapPeaC.shortest_curves_become_meridians(self.triangulation)
    
    def current_fillings_become_meridians(self):
        SnapPeaC.current_fillings_become_meridians(self.triangulation)
    
    def simplify(self):
        SnapPeaC.basic_simplification(self.triangulation)
    
    def randomize(self):
        SnapPeaC.randomize_triangulation(self.triangulation)
    
    def reflect(self):
        SnapPeaC.reorient(self.triangulation)
    
    def canonize(self):
        SnapPeaC.proto_canonize(self.triangulation)
    
    def is_canonical_triangulation(self):
        return SnapPeaC.is_canonical_triangulation(self.triangulation)
    
    def symmetry_group(self):

        theRawResults = SnapPeaC.symmetry_group(self.triangulation)

        thePackagedResults = {}
        
        if theRawResults[0] == 0:
            thePackagedResults['manifold'] = None
        else:
            thePackagedResults['manifold'] = SymmetryGroup(theRawResults[0], theRawResults[3])
        
        if theRawResults[1] == 0:
            thePackagedResults['link'] = None
        else:
            thePackagedResults['link'] = SymmetryGroup(theRawResults[1], theRawResults[3])
        
        if theRawResults[2] == 0:
            thePackagedResults['symmetric triangulation'] = None
        else:
            thePackagedResults['symmetric triangulation'] = Triangulation(theRawResults[2])

        return thePackagedResults
    
    def Dirichlet(self, centroid_at_origin=1, maximize_inj_radius=0,
                displacement=(0,0,0)):
        return DirichletDomain(SnapPeaC.Dirichlet(    self.triangulation,
                                                    centroid_at_origin,
                                                    maximize_inj_radius,
                                                    displacement[0],
                                                    displacement[1],
                                                    displacement[2]))


class CuspedCensus:

    def __init__(self, num_tetrahedra = 5, orientable = 'o'):
    
        self.num_tetrahedra    = num_tetrahedra
        self.orientable        = (orientable == 'o' or orientable == 'O')

        if   num_tetrahedra == 5:
            self.num_manifolds = 415
        elif num_tetrahedra == 6:
            if self.orientable:
                self.num_manifolds = 962
            else:
                self.num_manifolds = 259
        elif num_tetrahedra == 7:
            if self.orientable:
                self.num_manifolds = 3552
            else:
                self.num_manifolds = 887
        else:
            self.num_manifolds = 0
    
    def __repr__(self):
        if   self.num_tetrahedra == 5:
            return 'census of all cusped hyperbolic 3-manifolds triangulable with 5 or fewer ideal tetrahedra'
        elif self.num_tetrahedra == 6:
            if self.orientable:
                return 'census of all orientable cusped hyperbolic 3-manifolds triangulable with 6 (but not fewer) ideal tetrahedra'
            else:
                return 'census of all nonorientable cusped hyperbolic 3-manifolds triangulable with 6 (but not fewer) ideal tetrahedra'
        elif self.num_tetrahedra == 7:
            if self.orientable:
                return 'census of all orientable cusped hyperbolic 3-manifolds triangulable with 7 (but not fewer) ideal tetrahedra'
            else:
                return 'census of all nonorientable cusped hyperbolic 3-manifolds triangulable with 7 (but not fewer) ideal tetrahedra'
        else:
            return 'empty census'
    
    def __len__(self):
        return self.num_manifolds
    
    def __getitem__(self, i):
        if not 0 <= i < self.num_manifolds:
            raise IndexError
        else:
            return Triangulation(self.num_tetrahedra, self.orientable, i)

class ClosedCensus:

    def __init__(self, orientable = 'o'):
        self.orientable       = (orientable == 'o' or orientable == 'O')
        if self.orientable:
            filename = 'ClosedOrientableDistinct.txt'
        else:
            filename = 'ClosedNonorientableDistinct.txt'
        self.census_file      = open('Tables/ClosedCensusData' + os.sep + filename)
        self.linelength       = len(self.census_file.readline())
        self.census_file.seek(0,2)
        self.length           = int( (self.census_file.tell() + 1) / self.linelength )

    def __del__(self):
        self.census_file.close()

    def __repr__(self):
        if self.orientable:
            return 'Census of low-volume closed orientable hyperbolic manifolds.' 
        else:
            return 'Census of low-volume closed non-orientable hyperbolic manifolds.' 

    def __getitem__(self, i):
        if not (0 <= i < self.length):
            raise IndexError, "Requested manifold is out of range."
        self.census_file.seek(i*self.linelength)
        volume, number_of_tets, index, m, l = string.split(self.census_file.readline())
        manifold = Triangulation(int(number_of_tets), 1, int(index))
        manifold.set_cusp(0,int(m),int(l),0)
        return manifold    

    def __len__(self):
        return self.length

class AbelianGroup:
    
    def __init__(self, coefficients):
        self.coefficients = coefficients
    
    def __repr__(self):

        coef = self.coefficients
        n = len(coef)

        if (n == 0):
            return 'trivial'
        s = ''
        for i in range(n):
            s = s + 'Z'
            if coef[i] != 0:
                s = s + ('/%i' % coef[i])
            if i != n - 1:
                s = s + ' + '
        return s
    
    def __len__(self):
        return len(self.coefficients)
    
    def __getitem__(self, i):
        return self.coefficients[i]
    
    def Betti_number(self):
        return self.coefficients.count(0)


class FundamentalGroup:

    def __init__(self, fundamental_group):
        #    fundamental_group is an integer that is really
        #    a pointer to a SnapPea kernel FundamentalGroup.
        self.fundamental_group = fundamental_group

    def __del__(self):
        #    Let the SnapPea kernel free its private representation
        #    of a FundamentalGroup.
        SnapPeaC.free_group_presentation(self.fundamental_group)

    def __repr__(self):
        return (    'generators:  ' + self.generators_string() + '\n'
                  + 'relations:\n'
                  + self.relations_string() )
    
    def num_generators(self):
        return SnapPeaC.fg_get_num_generators(self.fundamental_group)
    
    def relations(self):
        return SnapPeaC.fg_get_relations(self.fundamental_group)

    def representation(self, word):
        return SnapPeaC.fg_representation(self.fundamental_group, word)

    def peripheral_curves(self):
        return SnapPeaC.fg_peripheral_curves(self.fundamental_group)
        
    def generators_string(self):
        theAlphabet = 'a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z'
        if self.num_generators() > 0:
            return theAlphabet[0 : 2*self.num_generators() - 1]
        else:
            return 'none'
        
    def relations_string(self):
        theString = ''
        for theWord in self.relations():
            theString = theString + theWord + '\n'
        return theString


class SymmetryGroup:

    def __init__(self, symmetry_group, is_full_group, owns_symmetry_group=True):
        #    symmetry_group is an integer that is really
        #    a pointer to a SnapPea kernel SymmetryGroup.
        self.symmetry_group            = symmetry_group
        self.is_full_group            = is_full_group
        self.owns_symmetry_group    = owns_symmetry_group

    def __del__(self):
        #    Let the SnapPea kernel free its private representation
        #    of a SymmetryGroup.
        if self.owns_symmetry_group:
            SnapPeaC.free_symmetry_group(self.symmetry_group)

    def __repr__(self):

        if self.is_full_group:
            thePretext = ''
        else:
            thePretext = 'at least '

        if self.is_abelian():
            theText = self.abelian_description().__repr__()
        elif self.is_dihedral():
            theText = 'D%d'%(self.order()/2)
        elif self.is_polyhedral():
            theText = self.polyhedral_description()['name']
        elif self.is_S5():
            theText = 'S5'
        elif self.is_direct_product():
            theText =     self.get_factor(0).__repr__()    \
                      + ' X '                            \
                      + self.get_factor(1).__repr__()
        else:
            theText = 'nonabelian group of order %d'%self.order()
        
        return thePretext + theText
    
    def order(self):
        return SnapPeaC.symmetry_group_order(self.symmetry_group)
    
    def is_abelian(self):
        return SnapPeaC.symmetry_group_is_abelian(self.symmetry_group)

    def abelian_description(self):
        return AbelianGroup(SnapPeaC.symmetry_group_abelian_description(self.symmetry_group))
    
    def is_dihedral(self):
        return SnapPeaC.symmetry_group_is_dihedral(self.symmetry_group)
    
    def is_polyhedral(self):
        return SnapPeaC.symmetry_group_is_polyhedral(self.symmetry_group)
    
    def polyhedral_description(self):
        return SnapPeaC.symmetry_group_polyhedral_description(self.symmetry_group)
    
    def is_S5(self):
        return SnapPeaC.symmetry_group_is_S5(self.symmetry_group)
    
    def is_direct_product(self):
        return SnapPeaC.symmetry_group_is_direct_product(self.symmetry_group)
    
    def get_factor(self, i):
        return SymmetryGroup(    SnapPeaC.symmetry_group_factor(self.symmetry_group, i),
                                True,
                                False)
    
    def is_amphicheiral(self):
        return SnapPeaC.symmetry_group_is_amphicheiral(self.symmetry_group)
    
    def is_invertible_knot(self):
        return SnapPeaC.symmetry_group_invertible_knot(self.symmetry_group)
        
    
    #    Please don't confuse abelian_description() (above)
    #    with abelianization() (below).  The former describes
    #    a group which already happens to be abelian, while
    #    the latter gives the quotient of the group by its
    #    commutator subgroup.

    def commutator_subgroup(self):
        return SymmetryGroup(    SnapPeaC.symmetry_group_commutator_subgroup(self.symmetry_group),
                                self.is_full_group,
                                True)

    def abelianization(self):
    
        if self.is_full_group:
        
            #    Compute the abelianization as a pointer to a SnapPea kernel
            #    internal SymmetryGroup data structure.
            theSnapPeaGroup = SnapPeaC.symmetry_group_abelianization(self.symmetry_group)
            
            #    Create the corresponding Python AbelianGroup object.
            theAbelianization = AbelianGroup(SnapPeaC.symmetry_group_abelian_description(theSnapPeaGroup))

            #    Free the SnapPea kernel's data structure.
            SnapPeaC.free_symmetry_group(theSnapPeaGroup)
            
            #    Return the Python object.
            return theAbelianization
            
        else:
            return None

    def center(self):
        if self.is_full_group:
            return SymmetryGroup(    SnapPeaC.symmetry_group_center(self.symmetry_group),
                                    True,
                                    True)
        else:
            return None
    
    def presentation(self):
        if self.is_full_group:
            return SnapPeaC.symmetry_group_presentation(self.symmetry_group)
        else:
            return None
    
    def presentation_text(self):
        if self.is_full_group:
            theAlphabet = 'a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z'
            thePresentation = self.presentation()
            theText = '{'
            theText = theText + theAlphabet[0 : 2*thePresentation['number of generators'] - 1]
            theText = theText + ' |'
            theLowerAlphabet = 'abcdefghijklmnopqrstuvwxyz'
            theUpperAlphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
            theFirstRelationFlag = 1
            for theRelation in thePresentation['relations']:
                if theFirstRelationFlag == 0:
                    theText = theText + ',' 
                else:
                    theFirstRelationFlag = 0
                for theFactor in theRelation:
                    if theFactor[1] < 0:
                        theLetter = theUpperAlphabet[theFactor[0]]
                        thePower  = - theFactor[1]
                    else:
                        theLetter = theLowerAlphabet[theFactor[0]]
                        thePower =   theFactor[1]
                    if thePower > 1:
                        theText = theText + ' %c^%i'%(theLetter, thePower)
                    elif thePower == 1:
                        theText = theText + ' %c'%(theLetter)
                    else:
                        raise RuntimeError, 'zero exponent in symmetry group presentation'
            theText = theText + ' }'
            return theText
        else:
            return None


class DirichletDomain:

    def __init__(self, Dirichlet_domain):
        #    Dirichlet_domain is an integer that is really
        #    a pointer to a SnapPea kernel WEPolyhedron.
        self.Dirichlet_domain = Dirichlet_domain

    def __del__(self):
        #    Let the SnapPea kernel free its private representation
        #    of a WEPolyhedron.
        SnapPeaC.free_Dirichlet_domain(self.Dirichlet_domain)

    def __repr__(self):
        return 'Dirichlet domain v=? e=? f=?'
    
    def off(self):
    
        theText = 'OFF %d %d 0\n\n'%(self.v(), self.f())
        
        theVertices = self.vertices()
        for i in range(self.v()):
            theText = theText + '%11.8f %11.8f %11.8f\n'%theVertices[i]
        theText = theText + '\n'
        
        theFaces = self.faces()
        theFaceColors = self.face_colors()
        for i in range(self.f()):
            theText = theText + '%2d'%len(theFaces[i])
            for j in range(len(theFaces[i])):
                theText = theText + ' %2d'%theFaces[i][j]
            theText = theText + '  %5.3f %5.3f %5.3f 1.000\n'%theFaceColors[i]
        
        return theText
    
    def v(self):
        return SnapPeaC.Dirichlet_num_vertices(self.Dirichlet_domain)
    
    def e(self):
        return SnapPeaC.Dirichlet_num_edges(self.Dirichlet_domain)
    
    def f(self):
        return SnapPeaC.Dirichlet_num_faces(self.Dirichlet_domain)

    def vertices(self):
        return SnapPeaC.Dirichlet_vertices(self.Dirichlet_domain)

    def faces(self):
        return SnapPeaC.Dirichlet_faces(self.Dirichlet_domain)

    def face_colors(self):
        return SnapPeaC.Dirichlet_face_colors(self.Dirichlet_domain)
    
    def face_pairings(self):
        return SnapPeaC.Dirichlet_face_pairings(self.Dirichlet_domain)
