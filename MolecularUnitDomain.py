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
# Defines boundary conditions, etc for microdomain
from dolfin import *
from params import *
from Domain import *


#
class MolecularUnitDomain(Domain):
  def __init__(self,fileMesh,fileSubdomains,\
    # scalar (wrong) or field
    type="field",\
    boundaryTol=1e-1,\
    outpath="./",\
    name = "Molecular",\
    reflectiveBoundary="none"  # use this to force outer cell boundary as reflective, instead of dirichlet 
    ):
    super(MolecularUnitDomain,self).__init__(type,EPS=boundaryTol)
    problem = self.problem
    problem.fileMesh = fileMesh
    problem.fileSubdomains = fileSubdomains
    problem.init = 1
    self.reflectiveBoundary = reflectiveBoundary
    problem.name = name
    problem.outpath = outpath

  

  def Setup(self):
    # mesh
    problem = self.problem
    if MPI.rank(mpi_comm_world())==0:
      print "Attempting to load ", problem.fileMesh
    mesh = Mesh(problem.fileMesh)
    problem.mesh = mesh
    # mesh is in A, so converting to um
    # DISABLED problem.mesh.coordinates()[:] *= parms.Ang_to_um
    # Something is flawed here, since surface area==0 if conversion is used

    utilObj=self.utilObj
    utilObj.GeometryInitializations()
    # PKH PBC Not used anymore 
    #utilObj.DefinePBCMappings()

    problem.V = VectorFunctionSpace(mesh,"CG",1)

    # geom
    self.CalcGeom(problem)


  # bcs
  def AssignBC(self,uBoundary=0):
    problem = self.problem
    nDim = np.shape(problem.mesh.coordinates())[1]

    #print "Probably don't have the right BCs yet"
    # TODO: might need to handle BC at active site as a mixed boundary
    # condition to take into account saturation
    u0 = Constant(np.zeros(nDim))
    u1 = Constant(np.ones(nDim))

    # use user-provided BC instead  
    if(uBoundary != 0):
      u1 = uBoundary

  # Create Dirichlet boundary condition

    bcs = []
    centerDomain = self.CenterDomain()
    centerDomain.problem = self.problem
    #PKHfixed_center = DirichletBC(problem.V, Constant((0,0,0)), centerDomain, "pointwise")
    fixed_center = DirichletBC(problem.V, Constant(np.zeros(nDim)), centerDomain, "pointwise")
    bcs.append(fixed_center)

    #PKH 120901 leftRightBoundary=self.PeriodicLeftRightBoundary()
    leftRightBoundary=self.LeftRightBoundary()
    leftRightBoundary.problem = self.problem
    #PKH 120901 bc1 = PeriodicBC(problem.V.sub(0), leftRightBoundary)
    tc1 = DirichletBC(problem.V.sub(0), Constant(1.),leftRightBoundary)
    bc1 = DirichletBC(problem.V.sub(0), Constant(0.),leftRightBoundary)
    if(self.reflectiveBoundary!="leftright"):
      bcs.append(bc1)
    else:
      print "Labeling %s as reflective" % self.reflectiveBoundary

    #PKH 120901 backFrontBoundary=self.PeriodicBackFrontBoundary()
    backFrontBoundary=self.BackFrontBoundary()
    backFrontBoundary.problem = self.problem
    #PKH 120901 bc2 = PeriodicBC(problem.V.sub(1), backFrontBoundary)
    tc2 = DirichletBC(problem.V.sub(1), Constant(1.),backFrontBoundary)
    bc2 = DirichletBC(problem.V.sub(1), Constant(0.),backFrontBoundary)
    if(self.reflectiveBoundary!="backfront"):
      bcs.append(bc2)
    else:
      print "Labeling %s as reflective" % self.reflectiveBoundary
    #print self.reflectiveBoundary
    #quit()

    #PKH 120901 topBottomBoundary=self.PeriodicTopBottomBoundary()
    if(nDim>2):
      topBottomBoundary=self.TopBottomBoundary()
      topBottomBoundary.problem = self.problem
      #PKH 120901 bc3 = PeriodicBC(problem.V.sub(2), topBottomBoundary)
      tc3 = DirichletBC(problem.V.sub(2), Constant(1.),topBottomBoundary)
      bc3 = DirichletBC(problem.V.sub(2), Constant(0.),topBottomBoundary)
      if(self.reflectiveBoundary!="topbottom"):
        bcs.append(bc3)
      else:
        print "Labeling %s as reflective" % self.reflectiveBoundary

    testBC=1
    if(testBC==1):
      z = Function(problem.V)
      tc1.apply(z.vector())
      tc2.apply(z.vector())
      if(nDim>2):
        tc3.apply(z.vector())
      File("appliedBCs.pvd") << z
    
    problem.bcs = bcs

