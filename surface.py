#$Id$
#   t3m - software for studying triangulated 3-manifolds
#   Copyright (C) 2002 Marc Culler, Nathan Dunfield and others
#
#   This program is distributed under the terms of the 
#   GNU General Public License, version 2 or later, as published by
#   the Free Software Foundation.  See the file GPL.txt for details.

from simplex import *
from Numeric import *
from LinearAlgebra import *
import sys

# NOTE (1) The functions in this module only make sense for closed
# manifolds.  It will need to be rewritten to accomodate spun normal
# surfaces.  In particular, build_weights tries to compute the
# triangle weights from the quad weights.  We could set them to
# infinity, I suppose, near a torus cusp.
#
# For spun surfaces we will want to compute the boundary slope.  We
# should perhaps decide what the degenerate form of the boundary slope
# is for a spherical vertex link. (0/0?) For a higher genus vertex
# link the object that corresponds to a boundary slope is an element
# of H^1 of the link.  This is the class represented by the shift
# cocycle.  The boundary class generates the kernel of the shift
# class, in the genus 1 case.  Now that I mention it, I think that
# the shift class is the natural object to compute in all cases.
#
# NOTE (2) Many of the functions in this module also assume that the
# triangulation has only one vertex.  If there are more vertices then one
# has to be more careful about saying things like "bounds a thick subcomplex"
# or "bounds a thin subcomplex" - it might bound both.  And an edge linkig
# surface might not be a torus.  It might be good to distinguish surfaces
# that bound regular neighborhoods of graphs from other surfaces that bound
# thin subcomplexes.
#
# NOTE (3) Our plan is to create (at least) three subclasses of the
# Surface class: Closed_Surface, Spun_Surface, Bounded_Surface.

#Incidence dictionaries for quads, triangles and octagons

MeetsQuad = {E01:array((1,1,0)), E02:array((1,0,1)), E21:array((0,1,1)),
             E32:array((1,1,0)), E31:array((1,0,1)), E03:array((0,1,1))}

MeetsTri = {E01:array((1,1,0,0)), E02:array((1,0,1,0)), E21:array((0,1,1,0)),
            E32:array((0,0,1,1)), E31:array((0,1,0,1)), E03:array((1,0,0,1))}

MeetsOct =  {E01:array((1,1,2)), E02:array((1,2,1)), E21:array((2,1,1)),
             E32:array((1,1,2)), E31:array((1,2,1)), E03:array((2,1,1))}

QuadWeights = (array((1,0,0)), array((0,1,0)), array((0,0,1)) )

WeightVector = array([1,1,1])

TypeVector = array([0,1,2]) 

# Used for converting normal surface into tetrahedron edge-shift data.
# QuadShift[k] is the shifts induced by quad Qk3 along edges (E03, E13, E23).
# Note that this follows the convention that the order of the edges
# is the same as the order of the quads, and _not_ in order E01, E02, E03.

QuadShifts = ((0, 1, -1), (-1, 0, 1), (1, -1, 0))

# The format for a coefficient vector is [T0, T1, T2, T3, Q0, Q1, Q2, ...}

NonInteger = 'Error'

# NOTE: The convention is that the order of the quads is (Q03, Q13, Q23)

class Surface:
  Count = 0

  def __init__(self, manifold, quadvector):
    Q = not_equal(quadvector, 0).resize((len(manifold),3))
    A = array(quadvector).resize((len(manifold),3))
    Surface.Count = Surface.Count + 1
    self.Manifold = manifold
    self.Coefficients = matrixmultiply(A,WeightVector)
    self.Quadtypes = matrixmultiply(Q,TypeVector)
    
  def __del__(self):
    Surface.count = Surface.Count - 1

  def erase(self):
    self.Manifold = None

  def type(self):
    if min(self.Coefficients) < 0:
      return "almost-normal"
    else:
      return "normal"

  # computes and records hexagon shift of surface along
  # the edges of each tet.  Order convention of edges in
  # each tet is (E03, E13, E23) the same as the standard
  # order of the quads, and _not_ (E01, E02, E03).

  def add_shifts(self):
    shifts = []
    for i in range(len(self.Manifold)):
        shifts += [ self.Coefficients[i] * w for w in QuadShifts[self.Quadtypes[i]]]
    self.Shifts = shifts

  def info(self, out = sys.stdout):
    M = self.Manifold
    if self.type() == "normal":
      out.write("Normal surface\n")
    for i in range(len(M)):
      quad_weight = self.Coefficients[i]
      if quad_weight == -1:
        weight = "  Quad Type Q%d3, weight: octagon" % self.Quadtypes[i]
      elif quad_weight > 0:
        weight = "  Quad Type  Q%d3, weight %d" % (self.Quadtypes[i], quad_weight)
      else:
        weight = "No quads"
      out.write(weight  + "\n")

class ClosedSurface(Surface):
  Count = 0

  def __init__(self, manifold, quadvector):
    Surface.__init__(self, manifold, quadvector)
    self.Weights = zeros( 7*len(manifold) )
    self.build_weights()
    
  def build_weights(self):

    """Use self.QuadWeights self.QuadTypes vector to construct
    self.Weights and self.EdgeWeights.  The vector self.Weights has size
    7T and gives the weights of triangles and quads in each 3-simplex.
    In each bank of 7 weights, the first 4 are triangle weights and the
    last 3 are quad weights."""

    #   Here is how it works.  We use the system of normal surface
    # equations defined as follows.  For each edge e, and each face f
    # containing e, let t1 and t2 be the two tetrahedra containing f.
    # The weight of e can be expressed as the sum of two quad weights
    # and two triangle weights in each of t1 and t2.  Set these
    # expressions equal to get one equation in the system.  The
    # resulting system can be written as Aw = Bq where w is the
    # (unknown) vector of triangle weights and q is the (known) vector
    # of quad weights.  Note that this system does not have positive
    # coefficients, and the matrix A is singular with null space
    # spanned by the triangle weights of the vertex links, which of
    # course have all quad weights equal to 0.  We can make the system
    # non-singular by restricting to the orthogonal complement of the
    # null space, i.e. by adding an equation for each vertex v which
    # sets the sum of the coefficients of triangles in the support of
    # Link(v) equal to 0.  The resulting system A'w = b is now
    # non-singular, but not square.  Since the linear equation solvers
    # like square systems, we multiply both sides by the transpose of
    # A before solving.  We then adjust the solution vector by adding the
    # unique linear combination of vertex links which produces a non-negative
    # vector w' with a coefficient of 0 on some triangle in each vertex links.
    # Now we can construct self.Weights by interleaving
    # the vector w' and the vector q of quadweights.
    #   The vector self.EdgeWeights is constructed by simply multiplying
    # the incidence matrix for edges meeting quads and triangles times the
    # vector self.Weights.
    
    eqns = []
    constants = []
    edge_matrix = []
    for edge in self.Manifold.Edges:

      # Build the row of the incidence matrix that records which quads
      # and triangles meet this edge.  We can use any tetrahedron that
      # meets this edge to compute the row.
      
      edge_row = zeros( 7*len(self.Manifold) )
      corner = edge.Corners[0]
      j = corner.Tetrahedron.Index
      edge_row[7*j:7*j+4] = MeetsTri[corner.Subsimplex]
      if not self.Coefficients[j] == -1:
        edge_row[7*j+4:7*j+7] = MeetsQuad[corner.Subsimplex]
      else:
        edge_row[7*j+4:7*j+7] = MeetsOct[corner.Subsimplex]
      edge_matrix.append(edge_row)
 
      # Build a row of A and a coefficient of the right hand side for
      # each pair of adjacent 3-simplices that meet this edge.
      for i in range(len(edge.Corners) - 1):
        j = edge.Corners[i].Tetrahedron.Index
        k = edge.Corners[i+1].Tetrahedron.Index
        row = zeros(4*len(self.Manifold))
        row[4*j:4*j+4] = MeetsTri[edge.Corners[i].Subsimplex]
        row[4*k:4*k+4] -= MeetsTri[edge.Corners[i+1].Subsimplex]
        eqns.append(row)
        c = 0
        if self.Coefficients[k] == -1:
          c = MeetsOct[edge.Corners[i+1].Subsimplex][self.Quadtypes[k]]
        else:
          if MeetsQuad[edge.Corners[i+1].Subsimplex][self.Quadtypes[k]]:
            c = self.Coefficients[k]
        if self.Coefficients[j] == -1:
          c -= MeetsOct[edge.Corners[i].Subsimplex][self.Quadtypes[j]]
        else:
          if MeetsQuad[edge.Corners[i].Subsimplex][self.Quadtypes[j]]:
            c -= self.Coefficients[j]
        constants.append(c)
      # Now add the extra equations to kill off the vertex links.
      for vertex in self.Manifold.Vertices:
        eqns.append(vertex.IncidenceVector)
        constants.append(0)

    A =  array(eqns)
    b =  array(constants)
    tA = transpose(A)
    U = matrixmultiply(tA,A)
    v = matrixmultiply(tA,b)

    x = solve_linear_equations(U,v)

    # Subtract off as many vertex links as possible.
    for vertex in self.Manifold.Vertices:
      m = min(compress(vertex.IncidenceVector, x))
      x -= m*vertex.IncidenceVector 

    for i in range(len(self.Manifold)):
      for j in range(4):
        # A hack, since we are not doing the linear algebra over Q.
        if round ( x[4*i+j] ) - x[4*i+j] > .0000001:
          print x
          print self.Coefficients
          print b
          print A
          raise NonInteger, 'Weight is not an integer!'
        self.Weights[7*i + j ] = round( x[4*i + j] )
      if not self.Coefficients[i] == -1:
        self.Weights[7*i + 4: 7*i + 7] = (
          self.Coefficients[i]*QuadWeights[self.Quadtypes[i]] )
      else:
        self.Weights[7*i + 4: 7*i + 7] = QuadWeights[self.Quadtypes[i]]

    self.EdgeWeights = matrixmultiply(array(edge_matrix),self.Weights)


  def euler_characteristic(self):
    # An EdgeValence is the number of tetrahedra that meet the edge.
    # The number of 2-simplices that meet the edge is larger by 1 in
    # the case of a boundary edge.
    valences = array(self.Manifold.EdgeValences)
    for edge in self.Manifold.Edges:
      if edge.IntOrBdry == 'bdry':
        valences[edge.Index] += 1
    V = sum(self.EdgeWeights)
    E = dot(self.EdgeWeights, valences)/2
    F = sum(abs(self.Weights))
    return V - E + F

  # takes either a triangle given as the corresponding vertex
  # or a quad given as an edge disjoint from it.
  
  def get_weight(self, tet_number, subsimplex):
    D = {V0: 0, V1:1, V2:2, V3:3, E03:4, E12:4, E13:5, E02:5, E23:6, E01:6}
    return self.Weights[7*tet_number + D[subsimplex] ]

  def get_edge_weight(self, edge):
    j = edge.Index
    return self.EdgeWeights[j]
  
  # The next function decides if a normal surface bounds a subcomplex.
  # The thing is to note is that given any surface, then there is a unique
  # maximal subcomplex disjoint from it -- consisting of all simplices
  # of any dimension disjoint from it.  (For a normal surface, the boundary
  # of a regular nbhd of this subcomplex is always normal.)   It's not
  # hard to see that a normal surface bounds a subcomplex iff all edge weights
  # are 0 or 2.  The function bounds_subcomplex returns the tuple
  #
  #  (bounds subcomplex, double bounds subcomplex, thick_or_thin)
  #
  # where thick_or_thin describes whether there is a tetrahedron contained
  # in the subcomplex bounded by the surface or its double.

  def bounds_subcomplex(self):
    if self.type() != "normal":
      return (0, 0, None)
    
    bounds_subcomplex = 1
    double_bounds_subcomplex = 1
    for w in self.EdgeWeights:
      if w != 0 and w != 2:
        bounds_subcomplex = 0
      if w != 0 and w != 1:
        double_bounds_subcomplex = 0
      if not (bounds_subcomplex or double_bounds_subcomplex):
        break


    if bounds_subcomplex or double_bounds_subcomplex:
      thick_or_thin = "thin"
      M = self.Manifold
      for i in range(len(M)):
        tet = M.Tetrahedra[i]
        inside = 1
        for e in OneSubsimplices:
          w = self.get_edge_weight(tet.Class[e])
          if w != 0:
            inside = 0
            break

        if inside:
          thick_or_thin = "thick"
          break
        
    else:
      thick_or_thin = None

    return (bounds_subcomplex, double_bounds_subcomplex, thick_or_thin)

###### It is not a torus unless the edge is a loop!
  # A surface is an edge linking torus iff all edge weights are 2 except one which
  # is zero.  Returns pair (is linking torus, edge it links around).
  
  def is_edge_linking_torus(self):
    zeroes = 0
    zero_index = None
    for i in range(len(self.EdgeWeights)):
      w = self.EdgeWeights[i]
      if w == 0:
        if zeroes > 0:
          return (0, None)
        zeroes = 1
        zero_index = i
      elif w != 2:
        return (0, None)
      
    return (1,  self.Manifold.Edges[zero_index])
        
  def info(self, out = sys.stdout):
    M = self.Manifold
    if self.type() == "normal":
      # check if really boring:
      q, e = self.is_edge_linking_torus()
      if q:
        out.write("Normal surface #%d is thin linking torus of edge %s\n"
                   %(M.NormalSurfaces.index(self), e))
        return
      out.write("Normal surface #%d of Euler characteristic %d\n"
                %(M.NormalSurfaces.index(self), self.euler_characteristic()))
      # addional message about bounding subcomplex
      b, d, t = self.bounds_subcomplex()
      if b == 1:
        out.write("  Bounds %s subcomplex\n"  % t)
      elif d == 1:
        out.write("  Double bounds %s subcomplex\n" %t)
      else:
        out.write("  doesn't bound subcomplex\n")
    else:
      out.write("Almost-normal surface #%d of Euler characteristic %d\n"
               % (M.AlmostNormalSurfaces.index(self), 
                  self.euler_characteristic()))
    out.write('\n') 
    for i in range(len(self.Manifold)):
      quad_weight = self.Coefficients[i]
      if quad_weight == -1:
        weight = "  Quad Type Q%d3, weight: octagon" % self.Quadtypes[i]
      elif quad_weight > 0:
        weight = "  Quad Type  Q%d3, weight %d" % (self.Quadtypes[i], quad_weight)
      else:
        weight = "No quads"
      out.write("  In tetrahedron %s :  %s\n" %
                      (self.Manifold.Tetrahedra[i], weight))
      out.write("\tTri weights V0: %d V1: %d V2 : %d V3 : %d\n" 
                % (self.get_weight(i, V0), 
                   self.get_weight(i, V1), 
                   self.get_weight(i, V2),
                   self.get_weight(i, V3)))
      out.write('\n') 

    for i in range(len(self.EdgeWeights)):
      out.write("  Edge %s has weight %d\n" 
                  % (self.Manifold.Edges[i], self.EdgeWeights[i]))

#-----------------end class ClosedSurface---------------------------------------


#-----------------begin class SpunSurface--------------------------------------

def dot_product(x,y):
    assert len(x) == len(y)
    dot = 0
    for i in range(len(x)):
        dot += x[i]*y[i]
    return dot

class SpunSurface(Surface):

  def __init__(self, manifold, quadvector):
    Surface.__init__(self, manifold, quadvector)
    self.Incompressible = None
    self.BoundarySlope = None
  
  def add_boundary_slope(surface, cusp_equations):
    surface.BoundarySlope = (-dot_product(surface.Shifts, cusp_equations[1]),
                             dot_product(surface.Shifts, cusp_equations[0]) )

  def info(self, out = sys.stdout):
    M = self.Manifold
    out.write("SpunSurface: Slope %s;  Incompressible %s\n" % (self.BoundarySlope, self.Incompressible))
    for i in range(len(M)):
      quad_weight = self.Coefficients[i]
      if quad_weight > 0:
        weight = "  Tet %d: Quad Type  Q%d3, weight %d" % (i, self.Quadtypes[i], quad_weight)
      else:
        weight = "  Tet %d: no quads" % i
      out.write(weight  + "\n")


#-------------begin class ClosedSurfaceInCusped------------------------

class ClosedSurfaceInCusped(ClosedSurface):
  def __init__(self, manifold, quadvector):
    ClosedSurface.__init__(self, manifold, quadvector)
    self.Incompressible = None
    self.BoundarySlope = None
  

  def info(self, out = sys.stdout):
    M = self.Manifold
    out.write("ClosedSurfaceInCusped #%d:  Euler %d;  Incompressible %s\n" % (M.ClosedSurfaces.index(self), self.euler_characteristic(), self.Incompressible))
    # check if really boring:
    q, e = self.is_edge_linking_torus()
    if q:
      out.write("    is thin linking surface of edge %s\n" % e)
      return

    # addional message about bounding subcomplex
    b, d, t = self.bounds_subcomplex()
    if b == 1:
      out.write("  Bounds %s subcomplex\n"  % t)
    elif d == 1:
      out.write("  Double bounds %s subcomplex\n" %t)
    else:
      out.write("  Doesn't bound subcomplex\n")

    for i in range(len(self.Manifold)):
      quad_weight = self.Coefficients[i]
      if quad_weight > 0:
        weight = "  Quad Type  Q%d3, weight %d" % (self.Quadtypes[i], quad_weight)
      else:
        weight = "No quads"

      out.write("  In tet %s :  %s\n" %
                      (self.Manifold.Tetrahedra[i], weight))
      out.write("\tTri weights V0: %d V1: %d V2 : %d V3 : %d\n" 
                % (self.get_weight(i, V0), 
                   self.get_weight(i, V1), 
                   self.get_weight(i, V2),
                   self.get_weight(i, V3)))
      out.write('\n') 

    for i in range(len(self.EdgeWeights)):
      out.write("  Edge %s has weight %d\n" 
                  % (self.Manifold.Edges[i], self.EdgeWeights[i]))

