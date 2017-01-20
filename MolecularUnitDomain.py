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
import smol

markerActiveSite = 1
markerMolecularBoundary =4
markerOuterBoundary=5


#class TestRL(SubDomain):
#  def inside(self, x, on_boundary):
#    cond = (np.abs(x[0]- -1) < EPS or np.abs(x[0]-1) < EPS)
#    #if(on_boundary):
#    #  print "x",x[0],cond
#    return (on_boundary and cond)
#
#class TestTB(SubDomain):
#  def inside(self, x, on_boundary):
#    cond = (np.abs(x[1]- -1) < EPS or np.abs(x[1]-1) < EPS)
#    #if(on_boundary):
#    #  print "y",x[1],cond
#    return (on_boundary and cond)
#
#class TestBF(SubDomain):
#  def inside(self, x, on_boundary):
#    cond = (np.abs(x[2]- -1) < EPS or np.abs(x[2]-1) < EPS)
#    #if(on_boundary):
#    #  print "z",x[2],cond
#    return (on_boundary and cond)


#
# filePotential - electrostatic potential from APBS, interpolated to FE mesh
class MolecularUnitDomain(Domain):
  def __init__(self,fileMesh,fileSubdomains,filePotential="none",\
    # scalar (wrong) or field
    type="field",\
    boundaryTol=1e-1,\
    # doe mesh come from gamer?
    gamer=1,\
    outpath="./",\
    name = "Molecular",\
    reflectiveBoundary="none",  # use this to force outer cell boundary as reflective, instead of dirichlet 
    q = 2.0,  # charge Ca2+ 
    psi = "none" # usually none, unless psi potential is passed in 
    ):
    super(MolecularUnitDomain,self).__init__(type,EPS=boundaryTol)
    problem = self.problem
    problem.fileMesh = fileMesh
    problem.fileSubdomains = fileSubdomains
    problem.filePotential= filePotential 
    problem.init = 1
    self.reflectiveBoundary = reflectiveBoundary
    #print "Enforcing gamer==1"
    self.gamer = gamer
    problem.name = name
    problem.outpath = outpath
    problem.q = q 
    problem.psi = psi 
    self.markerOuterBoundary = markerOuterBoundary
    problem.smolMode = False

  

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

    if(self.type=="scalar"):
        problem.V = FunctionSpace(mesh,"CG",1)
    elif(self.type=="field"):
        problem.V = VectorFunctionSpace(mesh,"CG",1)

    #print "Not loading subdomains, since don't think they're needed"
    #problem.subdomains = MeshFunction(
    #  "uint", mesh, problem.fileSubdomains)
    #self.markers = [1,4,5]

    # load ESP 
    if(problem.filePotential!="none" or problem.psi!="none"):
      problem.smolMode=True
      self.InitializeElectrostatics()

    # geom
    self.CalcGeom(problem)


  # bcs
  def AssignBC(self,uBoundary=0):
    problem = self.problem
    nDim = np.shape(problem.mesh.coordinates())[1]

    #print "Probably don't have the right BCs yet"
    # TODO: might need to handle BC at active site as a mixed boundary
    # condition to take into account saturation
    if(self.type=="scalar"):
        raise RuntimeError("REPLACE THIS AS DEBUG OPTION. scalar approach is nonsensical") 
        u0 = Constant(0.)
        u1 = Constant(1.)
    elif(self.type=="field"):
        #u0 = Constant((0.,0,0.))
        #u1 = Constant((1.,1.,1.))
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


  def InitializeElectrostatics(self):
      problem = self.problem
      # load in psi from file 
      Vtemp = FunctionSpace(problem.mesh,"CG",1)
      if(problem.psi=="none"): 
        print "Loading electrostatic potential from file"        
        psi = Function(Vtemp,problem.filePotential)
   
      # load in psi from argument 
      else: 
        print "Loading electrostatic potential from argument"
        psi = problem.psi

      smol.ElectrostaticPMF(problem,psi,V=Vtemp,q=problem.q) # pmf stored internally             
      #File("pmftest.pvd") << problem.psi
      #File("pmftest.pvd") << psi
      #quit()
