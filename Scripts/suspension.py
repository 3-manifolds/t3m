from snappea_conversion import *

def  random_suspension(M):
    # replace the star of a randomly selected edge by a
    # a suspension with a randomly selected cone vertex.
    edge_choice = whrandom.randint(0,len(M.Edges)-1)
    edge = M.Edges[edge_choice]
    valence = edge.valence()
    if not edge.distinct() or valence < 3:
       print "Bad edge"
       return 0
    print "valence = %d" % valence 
    shift = whrandom.randint(0, valence - 1)
    a = edge.get_arrow()
    t, b = M.suspension_of_polygon(valence)
    t = t[shift : ] + t[ : shift ]
    b = b[shift : ] + b[ : shift ]
    if not M.replace_star(a, t, b):
       print 'replace_star failed'

def suspension_test(M,count):
    for i in range(count):
      random_suspension(M)
      M.easy_simplify()
      print "%d of %d" % (M.EdgeValences.count(4), len(M.EdgeValences)), M.EdgeValences
    


