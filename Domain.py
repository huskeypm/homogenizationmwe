"""
----------------------------------
Smolfin - a numerical solver of the Smoluchowski equation for interesting geometries
Copyright (C) 2012 Peter Kekenes-Huskey, Ph.D., huskeypm@gmail.com

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA


----------------------------------------------------------------------------
"""
from dolfin import *
from homogutil import *
import numpy as np
class empty:pass

# For marking the unit cell outer boundaries. 
# Assumes the domain is square. 
# Note that the criteria use problem.EPS - if there are problems marking a side correctly (e.g. noisy boundary), try adjusting the 
# EPS value. 

class Domain(object):
  #
  # CLASSES
  #
  class LeftRightBoundary(SubDomain):
          def inside(self, x, on_boundary):
              problem = self.problem 
              scale = problem.scale 
              return ( on_boundary and \
                       # ON EDGE 
                       # (np.abs(x[0]-problem.boundsMin[0]) < EPS or np.abs(x[0]-problem.boundsMax[0]) < EPS)\
                       # BEYOND EDGE
                       (x[0] < (scale[0]*problem.boundsMin[0] + problem.EPS) or x[0] > (scale[0]*problem.boundsMax[0]- problem.EPS))\
                     )

      
  class BackFrontBoundary(SubDomain):
          def inside(self, x, on_boundary):
              problem = self.problem 
              scale = problem.scale 
              return ( on_boundary and \
                       # ON EDGE 
                       # (np.abs(x[1]-problem.boundsMin[1]) < EPS or np.abs(x[1]-problem.boundsMax[1]) < EPS)\
                       # BEYOND EDGE
                       (x[1] < (scale[1]*problem.boundsMin[1] + problem.EPS) or x[1] > (scale[1]*problem.boundsMax[1]- problem.EPS))\
                       )
      
  class TopBottomBoundary(SubDomain):
          def inside(self, x, on_boundary):
              problem = self.problem 
              scale = problem.scale 
              result = ( on_boundary and \
                       # ON EDGE 
                       # (np.abs(x[2]-problem.boundsMin[2]) < EPS or np.abs(x[2]-problem.boundsMax[2]) < EPS)\
                       # BEYOND EDGE
                       (x[2] < (scale[2]*problem.boundsMin[2] + problem.EPS) or x[2] > (scale[2]*problem.boundsMax[2]- problem.EPS))\
                     )
              #print "x[2]:%f %f:%f/%f" % (x[2],problem.boundsMin[2],problem.boundsMin[2],problem.boundsMax[2])
              #print result
              return result 
      
      # Sub domain for Dirichlet boundary condition
  class CenterDomain(SubDomain):
          def inside(self, x, in_boundary):
              problem = self.problem 
              return all(near(x[i], problem.center_coord[i], problem.EPS) for i in range(problem.nDims))



  #
  # FUNCTIONS
  # 
  def __init__(self,type,EPS=1e-1):
    problem = empty()
    problem.scale = np.array([1.0,1.0,1.0])  # can use values less than 1.0 to assign more vertices to boundary (see BCs) 
    problem.gamma = 1.
    problem.volUnitCell= 1.
    problem.init = 1
    problem.EPS = EPS 

    self.utilObj = homogutil(problem)
    self.type = type
    self.problem = problem

  def AssignBC(self,uBoundary=0):
    problem = self.problem


  # Create Dirichlet boundary condition
    bcs = []
    #PKHfixed_center = DirichletBC(problem.V, Constant((0,0,0)), CenterDomain(), "pointwise")
    centerDomain = self.CenterDomain()
    centerDomain.problem = self.problem
    fixed_center = DirichletBC(problem.V, Constant((0,0,0)), centerDomain, "pointwise")
    bcs.append(fixed_center)

    #PKHbc1 = PeriodicBC(problem.V.sub(0), PeriodicLeftRightBoundary())
    leftRightBoundary=self.PeriodicLeftRightBoundary()
    leftRightBoundary.problem = self.problem
    #bc1 = PeriodicBC(problem.V.sub(0), leftRightBoundary)
    #bcs.append(bc1)
    #PKHbc2 = PeriodicBC(problem.V.sub(1), PeriodicBackFrontBoundary())
    backFrontBoundary=self.PeriodicBackFrontBoundary()
    backFrontBoundary.problem = self.problem
    #bc2 = PeriodicBC(problem.V.sub(1), backFrontBoundary)
    #bcs.append(bc2)
    #PKHbc3 = PeriodicBC(problem.V.sub(2), PeriodicTopBottomBoundary())
    topBottomBoundary=self.PeriodicTopBottomBoundary()
    topBottomBoundary.problem = self.problem
    #bc3 = PeriodicBC(problem.V.sub(2), topBottomBoundary)
    #bcs.append(bc3)

    problem.bcs = bcs



  def CalcGeom(self,problem,boundsMin=None,boundsMax=None):
    # SA
    mesh = problem.mesh
    areaExpr = Constant(1.) * ds(domain=mesh)
    # 
    area = assemble(areaExpr)#, mesh=problem.mesh)

    problem.surfaceArea = area
    # this is the 'beta' term in the Goel papers 

    # VOL 
    vol = assemble(Constant(1.) * dx(domain=mesh))
    problem.volume = vol
    #print "SA: %e [um^2]" % area
    #print "Volume: %e [um^3]" % vol  

    (boundsMin,boundIdxMin,boundsMax,boundIdxMax) = self.utilObj.CalcBounds(problem.mesh, boundsMin=boundsMin,boundsMax=boundsMax) 
    diff = boundsMax-boundsMin
    totVol = np.product(diff)
    if MPI.rank(mpi_comm_world())==0:
      print("Total volume (assuming rectangular): %e [um^3]" % totVol)
      print("volume fraction (assuming rectangular): %e [um^3]" % (vol/totVol))    
    problem.volUnitCell = totVol
    

    

