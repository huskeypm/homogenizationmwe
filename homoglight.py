# very simple module for homogenzation 

#from homog  import *
import homog
import MolecularUnitDomain
from dolfin import *
def runHomog(fileXML="test.xml",verbose=False,\
             reflectiveBoundary=None):
  fileSubdomains = None   
  molDomUnit = MolecularUnitDomain.MolecularUnitDomain(fileXML,fileSubdomains,\
                 reflectiveBoundary=reflectiveBoundary)          
  molDomUnit.Setup()
  molDomUnit.AssignBC()
  homog.solve_homogeneous_unit(molDomUnit,solver="gmres") 

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



