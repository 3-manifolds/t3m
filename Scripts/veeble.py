from generators import *
mflds = orientable_closed_manifolds()

#
#  Returns the first manifold after m which has no self-adjacent edge.
#
def next_after(m):
   while 1:
     m = m+1
     M = mflds[m]
     if M.has_self_adjacent():
       print M.EdgeValences
       M.erase()
     else:
       return M

#
#  A new Mcomplex method to compute the complexity.
#
def complexity(self):
  coefficients = (0,1,2,2)
  cx = 0
  for valence in self.EdgeValences:
     if valence > 3:
       try:
         cx = cx + coefficients[valence-4]
       except IndexError:
         cx = cx + 2
  return cx

Mcomplex.complexity = complexity

#
# A new Mcomplex method to report whether M has a self-adjacent edge.
#
def has_self_adjacent(self):
  for edge in self.Edges:
    if edge.self_adjacent():
      return 1
  return 0

Mcomplex.has_self_adjacent = has_self_adjacent

#
# A new Mcomplex method to list the self-opposite edges.
#
def self_opposite_edges(self):
    so_list = []
    for edge in self.Edges:
      so_list.append(edge.self_opposite())
    return so_list

Mcomplex.self_opposite_edges = self_opposite_edges

def veeble_move(self, edge = None):
    if edge is None:
      edge = self.Edges[ whrandom.randint(0, len(self.Edges) - 1) ]
    print 'move on edge of valence', edge.valence()
    a = edge.get_arrow()
    N,b = self.copy(base_arrow = a)
    c = N.split_star(b.Tetrahedron.Class[b.Edge])
    if not c:
      print "split failed"
      N.erase()
      return 0
    blowdowns = []
    for edge in c.radii():
      if edge.Vertices[0] != edge.Vertices[1]:
        blowdowns.append(N.Edges.index(edge))
    shift = whrandom.randint(0,len(blowdowns)-1)
    blowdowns = blowdowns[shift:] + blowdowns[:shift]
    for bd in blowdowns:
      P = N.copy()
      if not P.smash_star(P.Edges[bd]):
        P.erase()
        continue
      P.easy_simplify()
      if P.complexity() <= self.complexity() :
        self.split_star(a.Tetrahedron.Class[a.Edge])
        self.smash_star(self.Edges[bd])
        self.easy_simplify()
        N.erase()
        P.erase()
        return 1
      P.erase()
    print "All blowdowns increased complexity."
    N.erase()

Mcomplex.veeble_move = veeble_move

def smart_veeble_moves(self, repeats = 1, verbose=1):
   while repeats > 1:
     if self.has_self_adjacent():
       print 'Has self-adjacent edge :^)\n'
       return 1
     if verbose:
       print 'complexity:',  self.complexity()
       print 'valences: ', self.EdgeValences
     self.veeble_move()
     if verbose:
       print 'new complexity: ', self.complexity(), '\n'
     repeats = repeats - 1
   print 'NO self-adjacent edge yet. :^(\n'   
   return 0

Mcomplex.smart_veeble_moves = smart_veeble_moves

# Just do any veeble whatsoever.  Ignore the complexity,
def dumb_veeble_move(self, edge = None):
   if edge is None:
     edge = self.Edges[ whrandom.randint(0, len(self.Edges) -1) ]
   print 'move on edge of valence', edge.valence()
   a = edge.get_arrow()
   N,b = self.copy(base_arrow = a)
   c = N.split_star(b.Tetrahedron.Class[b.Edge])
   if not c:
     print "split failed - edge is self-adjacent"
     N.erase()
     return 0
   blowdowns = []
   for edge in c.radii():
     if edge.Vertices[0] != edge.Vertices[1]:
       blowdowns.append(N.Edges.index(edge))
   if len(blowdowns) == 0:
     print 'smash failed - nothing to smash.'
     N.erase()
     return 0
   shift = whrandom.randint(0,len(blowdowns)-1)
   bd = blowdowns[shift]
   self.split_star(a.Tetrahedron.Class[a.Edge])
   self.smash_star(self.Edges[bd])
   self.easy_simplify()
   N.erase()
   return 1

Mcomplex.dumb_veeble_move = dumb_veeble_move

def dumb_veeble_moves(self, repeats = 1, verbose=1):
   while repeats > 1:
     if self.has_self_adjacent():
       print 'Has self-adjacent edge :^)\n'
       return 1
     if verbose:
       print "complexity:",  self.complexity()
       print 'valences:', self.EdgeValences
     self.dumb_veeble_move()
     if verbose:
       print "new complexity: ", self.complexity(), "\n"
     repeats = repeats - 1
   print 'NO self-adjacent edge yet. :^(\n'   
   return 0

Mcomplex.dumb_veeble_moves = dumb_veeble_moves

# Do a veeble move that will be sure to turn an edge of valence 4
# into an edge of valence 3 so it will bo away when we do the
# 3->2 moves.
def new_veeble_move(self, edge = None):
   if edge is None:
     edge = self.Edges[ whrandom.randint(0, len(self.Edges) -1) ]
   print 'move on edge of valence', edge.valence()
   a = edge.get_arrow()
   N,b = self.copy(base_arrow = a)
   c = N.split_star(b.Tetrahedron.Class[b.Edge])
   if not c:
     print "split failed -- has self_adjacent edge."
     N.erase()
     return 0
   blowdowns = []
   d = c.copy()
   e = c.copy().next().reverse()
   while d.next() != c:
     e.next()
     if ( (d.south_tail().Vertices[0] != d.south_tail().Vertices[1]) and
          (d.north_head().valence() == 4 or  e.south_head().valence() == 4) ):
       print d.north_head().valence(), e.south_head().valence()
       blowdowns.append(N.Edges.index(d.south_tail()))
   if len(blowdowns) == 0:
     print "bad edge"
     N.erase()
     return 0
   shift = whrandom.randint(0,len(blowdowns)-1)
   blowdowns = blowdowns[shift:] + blowdowns[:shift]
   for bd in blowdowns:
     P = N.copy()
     if not P.smash_star(P.Edges[bd]):
       P.erase()
     else:
       break
   self.split_star(a.Tetrahedron.Class[a.Edge])
   self.smash_star(self.Edges[bd])
   self.easy_simplify()
   N.erase()
   P.erase()
   return 1

Mcomplex.new_veeble_move = new_veeble_move

def test():
   testcases = range(len(mflds))
   for i in testcases:
     print 'Manifold', i
     m = mflds[i]
     if m.dumb_veeble_moves(repeats=20,verbose=0) == 0:
       testcases.append(i)
     m.erase()
