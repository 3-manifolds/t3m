from normal import *
import LinearAlgebra

def qtlist(surface): 
   result = []   
   for i in range(len(surface.Quadtypes)):
     if surface.Coefficients[i] != 0:
        result.append(surface.Quadtypes[i] + 1)
     else:
        result.append(0)
   return result

def compatible(x,y):
   result = 1
   for i in range(len(x)):
     if x[i] != 0 and y[i] != 0 and x[i] != y[i]:
        result = 0
   return result

def build_cells(self):
    qtypes = []
    for surf in self.NormalSurfaces:
      qtypes.append(qtlist(surf))
    self.SolutionCells = [[]]
    for i in range(len(qtypes)):
      self.SolutionCells[0].append((i,))
    while 1:
      newcells = []
      for cell in self.SolutionCells[-1]:
        for i in range( cell[-1]+1, len(qtypes)):
          save = 1
          for n in cell:
             if not compatible(qtypes[n], qtypes[i]):
                save = 0
                break
          if save == 1:
             newcells.append(cell + (i,))
      if len(newcells) == 0:
         return
      else:
         for cell in newcells:
           for i in range(len(cell)):
             try:
                self.SolutionCells[-1].remove(cell[:i]+cell[i+1:])
             except:
                pass
         self.SolutionCells.append(newcells)

Mcomplex.build_cells = build_cells
  
def cell_dimension(self,i,j):
   rows = len(self.SolutionCells[i][j])
   A = zeros((rows,len(self)))
   for k in range(rows):
     A[k] = self.NormalSurfaces[self.SolutionCells[i][j][k]].Coefficients
   return LinearAlgebra.linear_least_squares( A, zeros( (rows,1) ) )[2]
Mcomplex.cell_dimension = cell_dimension  
