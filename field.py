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

from dolfin import *
from homogutil import *
import numpy as np



EPS = 1.e-10

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
  Dbulk = problem.d  # set diffusion coefficient 
  nDims = problem.nDims
  Dii  = Constant(Dbulk*np.ones(nDims))
  Aij = diag(Dii)  # for now, but could be anisotropic
  Atilde = Aij 

  
  # Identity matrix 
  Delta = Identity( mesh.ufl_cell().geometric_dimension()) #
 
  # LHS. 
  #   dot (A grad u + I )) =0 
  form = inner(Atilde*(grad(u) + Delta), grad(v))*dx  
  # note: we are mixing linear and bilinear forms, so we need to split
  # these up for assembler 
  LHS = lhs(form)
  RHS = rhs(form)
  a = LHS
  
  ## RHS terms 
  # add in RHS from earlier 
  L = RHS

  
  ## Compute solution
  # We can use different solver types, depending on the convergence properies of the system 
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


  ## Store solution 
  problem.x = x

  
  # save unprojected soln instead 
  fileName = problem.outpath + problem.name+"_unit.pvd"
  if MPI.rank(mpi_comm_world())==0:
    print("Writing ",fileName)
  #File(fileName) <<  problem.x   
  File(fileName) <<  problem.x   


  return problem


# Computes effective diffusion constant based on solution from solveHomog

def compute_eff_diff(domain):
  problem = domain.problem
  mesh = problem.mesh
  dim = mesh.ufl_cell().geometric_dimension()


  dx_int = dx(domain=mesh)

  Vscalar = FunctionSpace(mesh,"CG",1)
  us = TrialFunction(Vscalar)
  vs = TestFunction(Vscalar)

  ## get omega (e.g. integration term below) 
  # deff = 1/volUnitCell int( grad u + 1 dx ) 
  # treating each component independtly, following Goel's example in sec 2.7 
  import numpy as np
  omegas = np.zeros(dim)
  x = problem.x  
  # I iterate only over the diagonal, since my original diffusion constant is diagonal 
  for i in range(dim):
    ## This is used for evaluating grad(x)+1
    # x[i].dx(i) is the gradient of the ith component of our solution
    # along the ith direction 
    grad_Xi_component = x[i].dx(i)+Constant(1)
    form = grad_Xi_component * dx_int
    omegas[i] = assemble(form)

    # for visualizing intermediate results 
    #D_eff_project = Function(Vscalar)
    #outname = "diff%d.pvd" % i
    #solve(us*vs*dx_int==grad_Xi_component*vs*dx_int, D_eff_project)
    #File(outname)<<D_eff_project  

  
  
  

  if MPI.rank(mpi_comm_world())==0:
    print("omegasO ",omegas)

  # normalize by volumeUnitCell to get deff
  d_eff = problem.d*omegas
  d_eff /= problem.volUnitCell

  if MPI.rank(mpi_comm_world())==0:
   if(dim==3):
    print("d_eff= [%4.2f,%4.2f,%4.2f] for d=%4.2f"%\
      (d_eff[0],d_eff[1],d_eff[2],problem.d))
   else: 
    print("d_eff= [%4.2f,%4.2f] for d=%4.2f"%\
      (d_eff[0],d_eff[1],problem.d))
   print("problem.volUnitCell", problem.volUnitCell)

  # report volume frac for validation purposes 
  vol = assemble( Constant(1)*dx_int)#, mesh=mesh )  
  problem.volFrac = vol / problem.volUnitCell


  # store 
  problem.d_eff = d_eff
  
  return d_eff
