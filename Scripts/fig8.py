#$Id: fig8.py,v 1.1.1.1 2002/08/19 16:24:45 t3m Exp $
from mcomplex import *

Tets = [Tetrahedron(),Tetrahedron()]

Tets[0].attach(F0,Tets[1],(1,0,2,3))
Tets[0].attach(F1,Tets[1],(0,2,1,3))
Tets[0].attach(F2,Tets[1],(0,1,3,2))
Tets[0].attach(F3,Tets[1],(3,1,2,0))
M = Mcomplex(Tets)
del(Tets)




