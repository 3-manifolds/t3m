#$Id$
from simplex import *
from tetrahedron import *
from corner import *
from edge import *

class Vertex:

   Count = 0

   def __init__(self):
     Vertex.Count = Vertex.Count + 1
     self.Index = -1
     self.IntOrBdry = ''
     self.Corners = []      # Corners of type "0-simplex in Tetrahedron"
     self.Edges = []        # incident Edges
                             # An Edge will appear twice if both its endpoints
                             # are equal to this Vertex
   def __repr__(self):
     if self.Index > -1:
       return ('v' + str(self.Index) 
            + ' (' + self.IntOrBdry + ')')
     else:
       return '< floating vertex' + str(id(self)) + ' >'

   def __del__(self):
     Vertex.Count = Vertex.Count - 1
#     print 'Deleting vertex at ', id(self)

   def erase(self):
     for corner in self.Corners:
       corner.Tetrahedron.Class[corner.Subsimplex] = None
     for edge in self.Edges:
       try:
         edge.Vertices.remove(self)
       except:
         pass
     self.Index = -1

 
# The link of a vertex in an Mcomplex is a surface
# of arbitrary genus, possibly with non-empty boundary.
# For now I am pretending that links are closed and orientable
   def link_genus(self):
     sum = 12
     for edge in self.Edges:
       sum = sum - 6 + edge.valence()
     return sum/12
