#$Id$
# Loading this module extends the Mcomplex class by adding attributes and
# methods which deal with normal surfaces.

# We work with normal surface equations defined entirely in terms of
# quad types.  The equations specify that the sum of the height shifts
# around an edge is equal to 0. 
# The quad types are 0, 1 or 2 where quad type i separates vertices
# V3 and Vi from the other two.
# The triangle types are 0, 1, 2, 3 where triangle type i cuts off
# vertex i.

from string import *
from mcomplex import *
from surface import *
#from Numeric import *
import os
from array import array
from Solver import solver

# The height shift dictionaries for the three quad types.
Shift = {E01:(-1,1,0), E02:(1,0,-1), E21:(0,-1,1),
         E32:(-1,1,0), E31:(1,0,-1), E03:(0,-1,1)}

# The solution vector templates for the four vertex types.
VertexVector = {V0:(1,0,0,0), V1:(0,1,0,0),
                V2:(0,0,1,0), V3:(0,0,0,1)}


class Matrix:

     def __init__(self, rows, columns):
          self.rows = rows
          self.columns = columns
          self.matrix = array('i', rows*columns*[0]) 

     def __setitem__(self, ij, value):
          i, j = ij
          self.matrix[i*self.columns + j] = value

     def __getitem__(self, ij):
          i, j = ij
          return self.matrix[i*self.columns + j]

     def __repr__(self):
          result = ''
          for i in range(self.rows):
            result += repr(self.matrix[i*self.columns:(i+1)*self.columns].tolist())+'\n'
          return result
     
# Extensions to the Mcomplex class
          
def normal_init(self, tetrahedron_list):
     self.Tetrahedra = tetrahedron_list
     self.Edges                = []
     self.Vertices             = []
     self.NormalSurfaces       = []
     self.AlmostNormalSurfaces = []
     self.build()
Mcomplex.__init__ = normal_init

def build_matrix(self):
    int_edges = [edge for edge in self.Edges if edge.IntOrBdry == 'int']
    self.QuadMatrix = Matrix(len(int_edges), 3*len(self))
    for edge in int_edges:
      for corner in edge.Corners:
        i = int_edges.index(edge)
        j = corner.Tetrahedron.Index
        for k in range(3):
             self.QuadMatrix[i,3*j+k] += Shift[corner.Subsimplex][k]
    self.build_vertex_incidences()
Mcomplex.build_matrix = build_matrix

def build_vertex_incidences(self):
    for vertex in self.Vertices:
       vertex.IncidenceVector = zeros( 4*len(self) )
       for corner in vertex.Corners:
         j = corner.Tetrahedron.Index
         vertex.IncidenceVector[4*j:4*j+4] += VertexVector[corner.Subsimplex]
Mcomplex.build_vertex_incidences = build_vertex_incidences

def find_normal_surfaces(self):
     for surface in self.NormalSurfaces:
          surface.erase()
          self.NormalSurfaces.remove(surface)
     self.build_matrix()
     coeff_list = solver.find_vertices(self.QuadMatrix.rows,
                                  self.QuadMatrix.columns,
                                  self.QuadMatrix.matrix)
     for coeff_vector in coeff_list:
          self.NormalSurfaces.append(Surface(self, coeff_vector))
Mcomplex.find_normal_surfaces = find_normal_surfaces

# We need find_almost_normal_surfaces()


# BROKEN!  We need read/write functions for surfaces
# Load surfaces from file

#def load_surfaces(self, file_name="surf.fifo"):
#     f = open(file_name)
#     exec(f.read())
#     f.close()
#     for datum in normal_surfaces:
#          self.NormalSurfaces.append(Surface(self, datum[0], datum[1]))
#     for datum in almost_normal_surfaces:
#          self.AlmostNormalSurfaces.append(Surface(self, datum[0], datum[1]))
#Mcomplex.load_surfaces = load_surfaces

# Info for printing surfaces

def normal_surface_info(self):
   try:
     out = os.popen('less', 'w')
     for surface in self.NormalSurfaces:
          out.write("-------------------------------------\n\n")
          surface.info(out)
          out.write('\n')
   except IOError:
     pass
Mcomplex.normal_surface_info = normal_surface_info

def almost_normal_surface_info(self):
   try:
     out = os.popen('less','w')
     for surface in self.AlmostNormalSurfaces:
          out.write("-------------------------------------\n\n")
          surface.info(out)
          out.write('\n')
   except IOError:
     pass
Mcomplex.almost_normal_surface_info = almost_normal_surface_info


          
          
