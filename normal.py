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
from Numeric import *
import os

# The height shift dictionaries for the three quad types.
Shift = {E01:array((-1,1,0)), E02:array((1,0,-1)), E21:array((0,-1,1)),
         E32:array((-1,1,0)), E31:array((1,0,-1)), E03:array((0,-1,1))}

# The solution vector templates for the four vertex types.
VertexVector = {V0:array((1,0,0,0)), V1:array((0,1,0,0)),
                V2:array((0,0,1,0)), V3:array((0,0,0,1))}

# Extensions to the Mcomplex class

def normal_init(self, tetrahedron_list):
     self.Tetrahedra = tetrahedron_list
     self.Edges                = []
     self.Vertices             = []
     self.NormalSurfaces       = []
     self.AlmostNormalSurfaces = []
     self.build()
Mcomplex.__init__ = normal_init

def build_template(self):
    int_edges = [edge for edge in self.Edges if edge.IntOrBdry == 'int']
    self.Template = zeros( (len(int_edges), 3*len(self) ) )
    for edge in int_edges:
      for corner in edge.Corners:
        i = int_edges.index(edge)
        j = corner.Tetrahedron.Index
        self.Template[i,3*j:3*j+3] += Shift[corner.Subsimplex]
    self.build_vertex_incidences()
Mcomplex.build_template = build_template

def build_vertex_incidences(self):
    for vertex in self.Vertices:
       vertex.IncidenceVector = zeros( 4*len(self) )
       for corner in vertex.Corners:
         j = corner.Tetrahedron.Index
         vertex.IncidenceVector[4*j:4*j+4] += VertexVector[corner.Subsimplex]
Mcomplex.build_vertex_incidences = build_vertex_incidences

def write_pari_script(self, file_name):
    self.build_template()
    if file_name == None:
      gp = os.popen("gp -q", 'w')
      out = gp.write
    else:
      out = open(file_name, "w").write
    m,n = shape(self.Template)
    out('template = Mat([\\\n')
    for i in range(m):
      row = ''
      for j in range(n-1):
        row = row+str(self.Template[i][j])+','
      row = row+str(self.Template[i][j+1])
      if i < m-1:
        row = row + ';\\\n'
      out(row)
    out('])\nread("surf.gp")\ngo()\n')
    if file_name == None:
      gp.close()
Mcomplex.write_pari_script = write_pari_script

# Write a PARI script to find surfaces and then run it.

def run_pari_on_fifo(self):
   
   self.erase_surfaces()
   fd = os.open("surf.fifo", os.O_NONBLOCK)
   fifo = os.fdopen(fd, 'r')
   self.write_pari_script(None)
   exec(fifo)
   fifo.close()
   for datum in pari_normal:
     self.NormalSurfaces.append(Surface(self, datum[0], datum[1]))
   for datum in pari_almost_normal:
     self.AlmostNormalSurfaces.append(Surface(self, datum[0], datum[1]))

# The below preserves the name "surf.fifo" but its no longer a fifo
#  The optional arguement is the name of a file copy the resulting
#  surf.gp so that you don't have to regenerate the surfaces all the
#  time.

def run_pari_via_tempfile(self, save_name=None):
   
   self.erase_surfaces()
   try:
        os.remove("surf.fifo")
   except:
        None
   self.write_pari_script(None)
   fifo = open("surf.fifo")
   exec(fifo.read())
   fifo.close()
   if save_name:
        os.rename("surf.fifo", save_name)
   for datum in pari_normal:
     self.NormalSurfaces.append(Surface(self, datum[0], datum[1]))
   for datum in pari_almost_normal:
     self.AlmostNormalSurfaces.append(Surface(self, datum[0], datum[1]))


#Mcomplex.find_surfaces = run_pari_on_fifo

Mcomplex.find_surfaces = run_pari_via_tempfile

# Load surfaces from file

def load_surfaces(self, file_name="surf.fifo"):
     f = open(file_name)
     exec(f.read())
     f.close()
     for datum in pari_normal:
          self.NormalSurfaces.append(Surface(self, datum[0], datum[1]))
     for datum in pari_almost_normal:
          self.AlmostNormalSurfaces.append(Surface(self, datum[0], datum[1]))

Mcomplex.load_surfaces = load_surfaces

def erase_surfaces(self):
      for surface in self.NormalSurfaces:
        surface.erase()
        self.NormalSurfaces.remove(surface)
      for surface in self.AlmostNormalSurfaces:
        surface.erase()
        self.AlmostNormalSurfaces.remove(surface)
Mcomplex.erase_surfaces = erase_surfaces

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


          
          
