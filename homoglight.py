# very simple module for homogenzation 

#from homog  import *
import MolecularUnitDomain
from dolfin import *



class empty:pass


from dolfin import *
import numpy as np
from params import *
parms = params() # conflicts with params from smol

# classes
from DefaultUnitDomain import *
from MolecularUnitDomain import *
import field 



## solv. homog cell
def solve_homogeneous_unit(domain,debug=False,solver="gmres"):
  #print "WARNING: assuming D=1."
  parms.D = 1.0
  problem = domain.problem

  ## debug mode
  field.solveHomog(domain,solver=solver)
  d_eff = field.compute_eff_diff(domain)

  problem.d_eff = d_eff

def runHomog(fileXML="test.xml",verbose=False,\
             reflectiveBoundary=None):
  fileSubdomains = None   
  molDomUnit = MolecularUnitDomain(fileXML,fileSubdomains,\
                 reflectiveBoundary=reflectiveBoundary)          
  molDomUnit.Setup()
  molDomUnit.AssignBC()
  solve_homogeneous_unit(molDomUnit,solver="gmres") 

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
  for i,arg in enumerate(sys.argv):
    if(arg=="-file"):
      fileXML=sys.argv[i+1] 
    if(arg=="-mumps"):
      solver="mumps"



  runHomog(fileXML,verbose=True,solver=solver)



