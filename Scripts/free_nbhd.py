#$Id$
from mcomplex import *
import sys

# Creates a free neighborhood with target edge valence N 
# and k layers.

def free_nbhd(N,k):
    Tbase = Tetrahedron()
    Tbase.Name = "base"
    M = Mcomplex([Tbase])

    i = 0
    for one_subsimplex in OneSubsimplices:
      edge = Tbase.Class[one_subsimplex]
      M.add_fan(edge, N[i] - edge.valence())
      i = i+1
     
    return M

def add_layer(M, target_valence):
    for edge in M.Edges:
      if edge.IntOrBdry is "Bdry" 

#-----end of code ------------------------


if __name__ == "__main__":
    M = embedded_nbhd([5, 5, 5, 5, 5, 5])
    print len(M)
    for i in range(1000):
        M.randomize()
        print len(M)
        sys.stdout.flush()


#M.info()
#Mcomplex_to_geo(M, "s3")
        

