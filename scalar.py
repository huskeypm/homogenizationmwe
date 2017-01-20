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
# 

from dolfin import *
import numpy as np
from params import *
from dolfin import nabla_grad as grad


class empty:pass

# TODO - add boundary 
# might not be correct
def  CalcArea(mesh, boundary=-1):
  f = Constant(1)
  if(boundary==-1):
    A = f*ds
    a = assemble(A, mesh=mesh)
  else:
    A = f*ds(i)
    a = assemble(A, exterior_facet_domains=self.sub_domains, mesh=self.mesh)

  return a

# calculate concentration
def CalcConc(domain):
  problem = domain.problem 
  if(domain.gamer==0):
    problem.conc = assemble( problem.x * dx,mesh=problem.mesh)
  else:
    problem.conc /= assemble( Constant(1)*dx(1),mesh=problem.mesh)


## solv. homog cell
# See notetaker notes from 2012-04-19 (pg 11) 
# Using first-order representation from Higgins paper ()
# type - scalar is valid only for certain cases
def solveHomog(domain):
  problem = domain.problem 

  # use class definition for BCs, etc
  if(hasattr(problem,'init')): 
    print"Using class definition"

  # load everything w defauit values 
  else:
    print "Overriding - WARNING: Should be in its own class"

    # Create mesh and define function space
    mesh = problem.mesh     
    V = FunctionSpace(mesh, "CG", 1)
    problem.V = V

    # Define boundary conditions
    # dirichlet 
    # TODO verify  
    u0 = Constant(1.)
    bc = DirichletBC(V, u0, u0_boundary)
    problem.bcs = [bc]

  # neumann cond
  # TODO verify 
  # JOHAN 
  dcdn = Constant(-1.0)

  # Define variational problem
  V = problem.V
  u = TrialFunction(V) 
  problem.u = u
  v = TestFunction(V) 
  problem.v = v
  if(problem.gamer==0):
    a = parms.d * inner(grad(u), grad(v))*dx 
  else:
    a = parms.d * inner(grad(u), grad(v))*dx(1) 
    print"Need to double check integration for gamer-derived meshes"
    quit()

  f = Constant(0.0) 
  L = f*v*dx - dcdn*v*ds
  # Compute solution
  x = Function(V) 
  solve(a == L, x, problem.bcs)
  problem.x = x

  # save soln
  File(problem.name+"_unit.pvd") << problem.x
  

## compute effective diff const 
def compute_eff_diff(domain):
  from dolfin import assemble
  problem = domain.problem
  gamma = CalcArea(problem.mesh,boundary=-1)

  # TODO need to figure out how to represent ident

  # TODO FIXME     
  #I = delta
  # JOHAN 
  #integral = assemble( (grad(problem.u) + I ) * dx ) 
  print "not correct"
  integral = 1
  omegac = (1/gamma) * integral

  diffeff = omegac * parms.d
  problem.d_eff = diffeff 

  return diffeff



