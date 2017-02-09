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
# Defines boundary conditions, etestBC for microdomain
from dolfin import *
from params import *
from Domain import *


#
class MolecularUnitDomain(Domain):
  def __init__(self,fileMesh,\
    boundaryTol=1e-1,\
    outpath="./",\
    name = "Molecular",\
    boxMin=None,
    boxMax=None,
    reflectiveBoundary="none"  # use this to force outer cell boundary as reflective, instead of dirichlet
                               # this is equivalent to having an inclusion boundary coincident with the unit cell boundary
    ):
    super(MolecularUnitDomain,self).__init__(type,EPS=boundaryTol)
    problem = self.problem
    problem.fileMesh = fileMesh
    self.reflectiveBoundary = reflectiveBoundary
    self.boxMin = boxMin
    self.boxMax = boxMax
    problem.name = name
    problem.outpath = outpath

  

  def Setup(self):
    # mesh
    problem = self.problem
    if MPI.rank(mpi_comm_world())==0:
      print "Attempting to load ", problem.fileMesh
    mesh = Mesh(problem.fileMesh)
    problem.mesh = mesh

    utilObj=self.utilObj
    utilObj.GeometryInitializations()
    # PKH PBC Not used anymore 
    #utilObj.DefinePBCMappings()


    problem.V = VectorFunctionSpace(mesh,"CG",1)

    # geom
    self.CalcGeom(problem, 
      boundsMin=self.boxMin,
      boundsMax=self.boxMax) 
   
    


  # bcs
  def AssignBC(self,uBoundary=0):
    problem = self.problem
    nDim = np.shape(problem.mesh.coordinates())[1]

    # Create Dirichlet boundary conditions
    # In order to have a unique solution, we arbitrarily set a point to 0. This 
    # doesn't influence the determination of the effective coefficients, since they
    # are based on the gradient of the correction function solution. 
    bcs = []
    centerDomain = self.CenterDomain()
    centerDomain.problem = self.problem
    fixed_center = DirichletBC(problem.V, Constant(np.zeros(nDim)), centerDomain, "pointwise")
    bcs.append(fixed_center)

    # for solution in the x direction [sub(0)], we set the boundaries along one direction to 0
    # I don't recall why I'm setting one side to 1. though, since I was using these dirichlets to replace a complicated
    # periodic bc. It doesn't appear to influence the solution though  
    leftRightBoundary=self.LeftRightBoundary()
    leftRightBoundary.problem = self.problem
    #PKH 120901 bc1 = PeriodicBC(problem.V.sub(0), leftRightBoundary)
    testBC1 = DirichletBC(problem.V.sub(0), Constant(1.),leftRightBoundary) 
    bc1 = DirichletBC(problem.V.sub(0), Constant(0.),leftRightBoundary)
    if(self.reflectiveBoundary!="leftright"):
      bcs.append(bc1)
    else:
      print "Labeling %s as reflective" % self.reflectiveBoundary

    # same thing for top/bottom along y direction [sub(1)]
    backFrontBoundary=self.BackFrontBoundary()
    backFrontBoundary.problem = self.problem
    testBC2 = DirichletBC(problem.V.sub(1), Constant(1.),backFrontBoundary)
    bc2 = DirichletBC(problem.V.sub(1), Constant(0.),backFrontBoundary)
    if(self.reflectiveBoundary!="backfront"):
      bcs.append(bc2)
    else:
      print "Labeling %s as reflective" % self.reflectiveBoundary
    #print self.reflectiveBoundary
    #quit()

    # same thing for top/bottom along z direction [sub(2)]
    #PKH 120901 topBottomBoundary=self.PeriodicTopBottomBoundary()
    if(nDim>2):
      topBottomBoundary=self.TopBottomBoundary()
      topBottomBoundary.problem = self.problem
      #PKH 120901 bc3 = PeriodicBC(problem.V.sub(2), topBottomBoundary)
      testBC3 = DirichletBC(problem.V.sub(2), Constant(1.),topBottomBoundary)
      bc3 = DirichletBC(problem.V.sub(2), Constant(0.),topBottomBoundary)
      if(self.reflectiveBoundary!="topbottom"):
        bcs.append(bc3)
      else:
        print "Labeling %s as reflective" % self.reflectiveBoundary

    # for visualizing the applied bcs above 
    testBC=True
    if(testBC):
      z = Function(problem.V)
      testBC1.apply(z.vector())
      testBC2.apply(z.vector())
      if(nDim>2):
        testBC3.apply(z.vector())
      File("appliedBCs.pvd") << z
    
    problem.bcs = bcs

