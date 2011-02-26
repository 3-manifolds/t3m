#$Id: embedded_nbhds.py,v 1.1.1.1 2002/08/19 16:24:45 t3m Exp $
from mcomplex import *
import sys

# Creates an embedded star of a tetrahedron, where
# the ith edge has valence N[i]

def embedded_nbhd(N):
    Tbase = Tetrahedron()
    Tbase.Name = "base"
    M = Mcomplex([Tbase])

    i = 0
    for one_subsimplex in OneSubsimplices:
      edge = Tbase.Class[one_subsimplex]
      M.add_fan(edge, N[i] - edge.valence())
      i = i+1
     
    return M

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
        

