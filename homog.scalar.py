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
#
# source config.bash
# from homog import *
# (cell,mol) = domost()

from dolfin import *
import numpy as np
from params import *

# classes
from CellularDomain import *
from CellularUnitDomain import *
from MolecularUnitDomain import *
import scalar
import field

class empty:pass


def u0_boundary(x, on_boundary):
  return on_boundary

## solv. homog cell
def solve_homogeneous_unit(problem,type="scalar"):

  ## debug mode
  debug=1
  if(debug==1):
    print "In debug mode - using stored values"
    d_eff = np.array([2.,2.,2.])
    problem.d_eff = d_eff
    V = FunctionSpace(problem.mesh,"CG",1)
    u = Function(V)
    u.vector()[:] = 1
    problem.x = u
    return

  ## solve homog
  if(type=="scalar"):
    scalar.solveHomog(problem)
    D_eff = scalar.compute_eff_diff(problem)
    d_eff = np.array([D_eff,D_eff,D_eff])

  elif (type=="field"):
    field.solveHomog(problem)
    d_eff = field.compute_eff_diff(problem)

  else:
    print "Not supported"
    quit()

  problem.d_eff = d_eff

## Coupled problem (need to check on neumann confs)
# Here I'm just trying to get a general time-dep soln to work. NOthing is
# relevant to my problem yet
def solve_homogenized_whole(wholecell,mol,type="scalar"):

  out  = File(wholecell.name+"_homogenized.pvd")


  print "WARNING: overwriting anistropic D const. Needs to be fixed"
  print "WARNING: must adjust both unit wholecell and whole to enforce VecFunc basis"
  wholecell.d_eff = wholecell.d_eff[0]


  # set up time dep form
  wholecell.u = TrialFunction(wholecell.V)
  wholecell.v = TestFunction(wholecell.V)
  # Assembly of the K, M and A matrices
  K = assemble(wholecell.d_eff * inner(grad(wholecell.u), grad(wholecell.v))*dx,mesh=wholecell.mesh)
  M = assemble(wholecell.u*wholecell.v*dx,mesh=wholecell.mesh)
  #E = assemble(-u*inner(a, grad(v))*dx,mesh=mesh)


# if(1):
#   print "Test, erase me"
#   V = FunctionSpace(wholecell.mesh, "CG", 1)
#   u = TrialFunction(V)
#   v = TestFunction(V)
#   K = assemble(wholecell.d_eff * inner(grad(u), grad(v))*dx,mesh=wholecell.mesh)
#   A = K.copy()
#   A.assign(K)
#   quit()


  u_n = Function(wholecell.V)
  A = K.copy()
  b = Vector(A.size(1))
  b[:]=0.0
  E = assemble(wholecell.dudn*wholecell.v*ds,mesh=wholecell.mesh)
  b += E

  x = u_n.vector()
  #x[:] = wholecell.x.vector()[:] # pass 'x' values from initial solution of homogeneous domains
  x[:] = 0.1 # for init cond
  wholecell.x = Function(wholecell.V)
  wholecell.x.vector()[:] = x[:]
  #wholecell.u.vector()[:] = u0
  #mol.u.vector()[:] = 0.5*u0
  #mol.x = x

  # for multiple bcs
  # for bc in bcs:
  #   bc.apply(A,b)

  dt =0.5
  t = 0.
  tMax = 1
  while (t < tMax):
    print "t %f" % t

    t  += float(dt)

    ## TODO check that flux is correct
    # assume simple flux between compartments
    if(type=="scalar"):
      scalar.CalcConc(wholecell)
      scalar.CalcConc(mol)
    elif(type=="field"):
      field.CalcConc(wholecell)
      field.CalcConc(mol)

    k = 1
    # TODO - need to understand how to get non-zero fluxes and why this matters
    hack = 0.5
    jcell = k*(wholecell.conc - hack*mol.conc)
    print "wholecell: %f" % wholecell.conc
    print "mol: %f" % mol.conc
    jmol  = -jcell

    # TODO add scale factors


    ## TODO add corret time stepping, equations
    # solve cell
    # JOHAN
    A.assign(K)
    A *= parms.d*dt
    A += M

    #  TODOGoel: how/where is the surface defined here? It seems like the boundaries within
    # the unit cell are 'embedded' into the large cell description, so are no longer
    # boundaries
    print" TODO - note: ds needs to be limited to boundary between cell and microdomains [eq use ds(2)]"
    b[:] = 0.0

    # cell + molec coupling
    #print"FOR NOW IGNORING MOLEC CONTRIB"
    #E = cell.d_eff * assemble( jcell * cell.v * ds, mesh=cell.mesh)
    #b += E

    # outer cell boundary
    # TODO verify
    E = assemble(wholecell.dudn*wholecell.v*ds,mesh=wholecell.mesh)
    b += E


    b += M*x
    if(hasattr(wholecell,'bc')):
      wholecell.bc.apply(A,b)

    solve(A,x,b,"gmres","default")

    # store solution 
    wholecell.x.vector()[:] = x[:]


    # store results 
    out << wholecell.x


    # solv mol
    #F = mol.diff_eff * grad(mol.u) * grad(mol.v)
    #F += mol.diff_eff( jmol * mol.v)
    #M = 0
    #solve(F==M,mol.u)
    #write


  quit()

def domost():
  parms.d = 1.
  cell = empty()
  mol  = empty()
  cell.name = "cell"
  mol.name = "mol"

  # TODO need to rescale size
  cell.mesh = UnitCube(8,8,8)
  mol.mesh  = UnitCube(8,8,8)

  solve_homogeneous_unit(cell)
  solve_homogeneous_unit(mol)

  compute_eff_diff(cell)
  compute_eff_diff(mol)

  return (cell,mol)

# test 2
def Test():
  root = "/home/huskeypm/scratch/homog/mol/"

  # celular domain
  prefix = "cell"
  meshFileOuter = root+prefix+"_mesh.xml.gz"
  subdomainFileOuter = root+prefix+"_subdomains.xml.gz"
  cellDomUnit = CellularUnitDomain(meshFileOuter,subdomainFileOuter)
  cellDomUnit.Setup(type="field")
  cellDomUnit.AssignBC()
  solve_homogeneous_unit(cellDomUnit.problem,type="field")


  # molecular domain
  prefix = "mol"
  meshFileInner = root+prefix+"_mesh.xml.gz"
  subdomainFileInner = root+prefix+"_subdomains.xml.gz"
  molDomUnit = MolecularUnitDomain(meshFileInner,subdomainFileInner)
  molDomUnit.Setup(type="field")
  molDomUnit.AssignBC()
  solve_homogeneous_unit(molDomUnit.problem,type="field")



  ## whole cell solutions
  # solve on actual cell domain
  prefix = "wholecell"
  meshFileOuter = root+prefix+"_mesh.xml.gz"
  subdomainFileOuter = root+prefix+"_subdomains.xml.gz"
  cellDomWhole = CellularDomain(meshFileOuter,subdomainFileOuter)
  print "WARNING: using scalar for right now"
  #cellDomWhole.Setup(type="field")
  cellDomWhole.Setup(type="scalar")
  # TODO: Need to replace, since using time-dep Neumann
  cellDomWhole.AssignBC()
  cellDomWhole.problem.d_eff = cellDomUnit.problem.d_eff
  solve_homogenized_whole(cellDomWhole.problem,molDomUnit.problem)


  quit()

# example 2.7 in Goel paper
def GoelEx2p7():
  ## micro
  unitLength = np.array([0.1,0.1,1.0]) # length of microdomain in unit cell
  diff = Constant(2.5) # isotropic

  # our geometry needs to look something like the 'cross' in Fig 2
  # creating this in Blender
  # mesh is crap, doesn't work
  meshFile = "goel.xml.gz"
  problem = empty()
  problem.mesh = Mesh(meshFile)

  #mesh = UnitCube(6,6,6)
  #problem.name ="test"
  #problem.mesh = mesh

  solve_homogeneous_unit(problem)
  quit()


  # note: authors take u0(x=[0.5,0.5,0,5]) = 0 to fix arb. const
  # i think I can treat this like an initial condition and set some
  # location in u0 to be zero.
  u0 = Constant(0)  # this doesn't do anything, but could use as templat ealter
  u_1 = interpolate(u0,V) # applies Expression u0 to FunctionSpace
  # solve stuff
  # u_1.assign(u)


##
## MAIN
##

if __name__ == "__main__":
  msg="script.py <arg>"
  remap = "none"

  #GoelEx2p7()
  Test()
  quit()

  import sys
  if len(sys.argv) < 1:
      raise RuntimeError(msg)


  (cell,mol) = domost()
  solve_homogenized_whole(cell,mol)

