#$Id$
from simplex import *
from tetrahedron import *

# A Corner is a "subsimplex in a tetrahedron".

class Corner:

   Count = 0

   def __init__(self, tetrahedron, subsimplex):
     Corner.Count = Corner.Count + 1
     self.Tetrahedron = tetrahedron
     self.Subsimplex = subsimplex

   def __repr__(self):
     return ('<'+SubsimplexName[self.Subsimplex]+ ' of '+ 
              str(self.Tetrahedron) + '>')

   def __del__(self):
     Corner.Count = Corner.Count - 1
#     print 'Deleting Corner at ', id(self)
