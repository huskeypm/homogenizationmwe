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
# Revisions
#       10.08.10 inception
#
#  source ~/sources/dolfin_smol/config.bash

from dolfin import *
from params import *
parms = params()
parms.d = 1.0 ; #print "WARNING: overriding parms (conflicting w smol on vm)"
# PKH not sure if this will work 
from homogutil import *
import numpy as np



EPS = 1.e-10
################/////////////////


# calculate concentration
def CalcConc(domain):
  problem = domain.problem
  mesh = problem.mesh 
  problem.conc = assemble( problem.x[0] * dx(domain=mesh))

  #problem.conc /= assemble( Constant(1)*dx,mesh=problem.mesh)
  problem.conc /= problem.volume
  

#
# Solve homogenized diff. eqn based on vector field 
#
## solv. homog cell
# See notetaker notes from 2012-04-19 (pg 11) 
# Using first-order representation from Higgins paper ()
# type - scalar is valid only for certain cases
# solverMode - gmres/mumps
def solveHomog(domain,solver="gmres"):     
  solverMode = solver 
  # mesh
  problem = domain.problem
  mesh = problem.mesh
  V = problem.V

  # Define variational problem
  u = TrialFunction(V)
  v = TestFunction(V)
  
  ## LHS terms 
  # Diffusion constant
  Dbulk = parms.d
  nDims = problem.nDims
  Dii  = Constant(Dbulk*np.ones(nDims))
  Aij = diag(Dii)  # for now, but could be anisotropic
  Atilde = Aij 

  
  # Identity matrix 
  Delta = Identity( mesh.ufl_cell().geometric_dimension()) #
 
  # LHS 
  form = inner(Atilde*(grad(u) + Delta), grad(v))*dx
  
  # note: we are mixing linear and bilinear forms, so we need to split
  # these up for assembler 
  LHS = lhs(form)
  RHS = rhs(form)
  a = LHS
  
  ## RHS terms 
  #print "WATNING: still confused about the BC here. Johan says ==0, Gary says we have complicated Jacobian term"
  #n = FacetNormal(mesh)
  #L = inner( n,v ) * ds 
  # Based on Johan's May 14 email, because (d X + I) .n = 0, we do not apply the BC. 
  
  # add in RHS from earlier 
  L = RHS

  
  ## Compute solution
  #solverType="original"
  solverType="krylov"
  #solve(a == L, x, problem.bcs)
  if solverType=="original":
    lvproblem = LinearVariationalProblem(a,L, x, bcs=problem.bcs)
    solver = LinearVariationalSolver(lvproblem)
    solver.parameters["linear_solver"] = "gmres"
    #solver.parameters["preconditioner"] = "ilu"
    #print "Using amg preconditioner instead of ilu"
    solver.parameters["preconditioner"] = "amg"
    #solver.parameters["newton_solver"]["linear_solver"] = "mumps"
    x = Function(V)
    solver.solve(x.vector())

  elif solverType=="krylov":
    #print "NOW USING KRYLOV SOLVER (from stokes-iterative example)" 
    #print "WARNING: made a random qguess as to what to use for the preconditioner [used basic diff eqn, no delta term]"
    # Form for use in constructing preconditioner matrix
    #b = inner(Atilde*(grad(u) + Delta), grad(v))*dx
    b = inner(Atilde*(grad(u)), grad(v))*dx


    # Assemble system
    A, bb = assemble_system(a, L, problem.bcs)

    # Assemble preconditioner system
    P, btmp = assemble_system(b, L, problem.bcs)

    # Create Krylov solver and AMG preconditioner
    solver = KrylovSolver("minres", "amg")
 
    # Associate operator (A) and preconditioner matrix (P)
    solver.set_operators(A, P)

    x = Function(V)
    solver.solve(x.vector(),bb)

  #print solver.parameters


  problem.x = x
  problem.up = Function(problem.V)   

  Vs = FunctionSpace(mesh,"CG",1)
  up = project(x[0],V=Vs)    
  #ar=np.asarray(up.vector())


  problem.up.vector()[:] = problem.x.vector()[:]

  
  # save soln
  #File(problem.name+"_unit.pvd") << problem.up

  # save unprojected soln instead 
  fileName = problem.outpath + problem.name+"_unit.pvd"
  if MPI.rank(mpi_comm_world())==0:
    print "Writing ",fileName
  #File(fileName) <<  problem.x   
  File(fileName) <<  problem.up


  return problem

def compute_eff_diff(domain):
  problem = domain.problem
  mesh = problem.mesh
  dim = mesh.ufl_cell().geometric_dimension()


  dx_int = dx(domain=mesh)

  Vscalar = FunctionSpace(mesh,"CG",1)
  us = TrialFunction(Vscalar)
  vs = TestFunction(Vscalar)

  ## get omega
  # treating each component independtly, following Goel's example in sec 2.7 
  import numpy as np
  omegas = np.zeros(dim)
  x = problem.up
  # Believe it is correct now print "WARNING: verify accuracy of this calculation"
  # I iterate only over the diagonal, since my original diffusion constant is diagonal 
  for i in range(dim):

    # JOHAN 
    D_eff_project = Function(Vscalar)

    ## apply leading exp(-PMF) term to diff. integral for homogenized smol eqn  
    grad_Xi_component = x[i].dx(i)+Constant(1)
    outname = "diff%d.pvd" % i

   
    #print "Solve for estimating Omegas" 
    solve(us*vs*dx_int==grad_Xi_component*vs*dx_int, D_eff_project)
    #File(outname)<<D_eff_project
  
    form = grad_Xi_component * dx_int
  
    omegas[i] = assemble(form)
  
  vol = assemble( Constant(1)*dx_int)#, mesh=mesh )
  

  if MPI.rank(mpi_comm_world())==0:
    print "omegasO ",omegas
  #print omegas/vol
  #problem.omv = omegas/vol
  #problem.domv = parms.d*omegas/vol
  #problem.d_eff = problem.domv 
  #print "WARNING: circumventing use of volUnitCell etc for determining deff"
  #return problem.d_eff 
  #print "vol", vol
  
  
  #print "WARNING: what is gamma?" 
  #omegas /= problem.gamma
  d_eff = parms.d*omegas
  d_eff /= problem.volUnitCell
  if MPI.rank(mpi_comm_world())==0:
   if(dim==3):
    print "d_eff= [%4.2f,%4.2f,%4.2f] for d=%4.2f"%\
      (d_eff[0],d_eff[1],d_eff[2],parms.d)
   else: 
    print "d_eff= [%4.2f,%4.2f] for d=%4.2f"%\
      (d_eff[0],d_eff[1],parms.d)
   print "problem.volUnitCell", problem.volUnitCell

  # thought I was already storing this somewhere
  problem.volFrac = vol / problem.volUnitCell

  # normalize
  #nd_eff= d_eff/ np.linalg.norm(d_eff)
  #print "Deff (normalized) ", nd_eff

  # store 
  problem.d_eff = d_eff
  
  return d_eff

class empty:pass

def doit(fileIn):
  # Create mesh and define function space
  defaultDom = DefaultUnitDomain()
  mesh = UnitCube(8,8,8)
  problem = empty()
  problem.mesh = mesh 
  solveHomog(defaultDom,type="field")
 
import sys

if __name__ == "__main__":
  msg="Purpose: Runs diff. eqn with field as test function "
  msg=msg+"Usage: "
  msg=msg+"script.py <arg>"
  msg=msg+"Notes:"
  remap = "none"


  import sys
  if len(sys.argv) < 2:
      raise RuntimeError(msg)

  fileIn= sys.argv[1]
  if(len(sys.argv)==3):
    1



  print "DOES NOT YIELD CORRECT ANSWERS FOR MPI" 
  doit(fileIn)


