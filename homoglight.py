# very simple module for homogenzation 

from dolfin import *
import numpy as np

## Import helper classes
# These classes contain functions the prepare the mesh for boundary condition assignments
# and other initializations. 
from DefaultUnitDomain import *
from MolecularUnitDomain import *

# Contains routines for computing homogenization corrector functions and using those solutions
# to estimate effective coefficients
import field 



## solv. homog cell
def solve_homogeneous_unit(domain,debug=False,solver="gmres"):
  #print "WARNING: assuming D=1."
  problem = domain.problem
  problem.d = 1.0  # diffusion constant

  ## solve 
  field.solveHomog(domain,solver=solver)
  problem.d_eff = field.compute_eff_diff(domain)

def runHomog(fileXML="test.xml",
             verbose=False,\
             boxMin=None,   
             boxMax=None,   
             solver = "gmres",  # default solver, but krylov is available 
             reflectiveBoundary=None):
  # domain set up              
  molDomUnit = MolecularUnitDomain(fileXML,\
                 boxMin=boxMin,
                 boxMax=boxMax,
                 reflectiveBoundary=reflectiveBoundary)          
  molDomUnit.Setup()
  molDomUnit.AssignBC()
  
  # homog solver 
  solve_homogeneous_unit(molDomUnit,solver=solver) 

  # report results
  problem = molDomUnit.problem
  if(verbose and MPI.rank(mpi_comm_world())==0):
    print "From master node:" 
    print "vol domain:", problem.volume
    print "vol unit cell:", problem.volUnitCell
    print "Deff:", problem.d_eff
  return problem 



#!/usr/bin/env python
import sys
#
# Revisions
#       10.08.10 inception
#
def test():
  print "Good morniing!"

if __name__ == "__main__":
  import sys
  scriptName= sys.argv[0]
  msg="""
Purpose: 
  Light-weight wrapper for running homogeniation on xml files 
 
Usage:
"""
  msg+="  %s -file <filename>" % (scriptName)
  msg+="""
  
 
Notes:

"""
  remap = "none"

  if len(sys.argv) < 2:
      raise RuntimeError(msg)

  fileIn= sys.argv[1]
  solver="gmres"
  fileXML = "none" 
  boxMin=None
  boxMax=None
  for i,arg in enumerate(sys.argv):
    if(arg=="-file"):
      fileXML=sys.argv[i+1] 
      
    # I don't think mumps is supported currently, but krylov is fair game   
    if(arg=="-mumps"):
      solver="mumps"



  runHomog(fileXML,verbose=True,solver=solver,boxMin=boxMin,boxMax=boxMax)



