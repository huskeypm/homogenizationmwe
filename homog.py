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

## TODO
# DONE add in PBC for cell
# use unit cell geometry that is based on the molec geom (right now its sphereically symm)  
# validate against COMSOL
# validate against Goel? (if mesh isn't too difficult)
# verify 'LCC' flux is being applied appopriately 
# Looks like I need to homogenize the fluxes as well (see Higgins) 



from dolfin import *
import numpy as np
from params import *
parms = params() # conflicts with params from smol

# classes
from CellularDomain import *
from DefaultUnitDomain import *
from CellularUnitDomain import *
from MolecularUnitDomain import *

from CellularUnitDomain_TnC import *

import scalar
import field

##
## VAR
## 
debug=0
smolMode = False # tells program that we want to solve the smol equation for the molec domain
smolPsi = "none"   # can pass in Function() containing electrostatic potential values for the molec domain
                   # 'none' means that it will be read from file 
smolq   = 2	   # Ca2+ charge used for smol eqn
case ="none"
molGamer = 0 
validationMode = 0
#outbase= "/tmp/outs/"
#outbase="/home/huskeypm/scratch/homogouts/"
outbase="./"
boundaryTolerance = 1e-1

# rocce
cellPrefix="none"
wholeCellPrefix="none"
molPrefix="none"
root = "/home/huskeypm/scratch/homog/"

# vm
#root = "/home/huskeypm/localTemp/"
 

class empty:pass

## solv. homog cell
def solve_homogeneous_unit(domain,type="field",debug=0,smolMode=False,solver="gmres"):
  #print "WARNING: assuming D=1."
  parms.D = 1.0
  problem = domain.problem

  ## debug mode
  if(debug==1):
  #if(1):
    print "In debug mode - using stored values"
    d_eff = np.array([2.,2.,2.])
    problem.d_eff = d_eff
    V = FunctionSpace(problem.mesh,"CG",1)
    u = Function(V)
    u.vector()[:] = 1
    problem.x = u
    return

  ## solve homog
  field.solveHomog(domain,solver=solver)
  d_eff = field.compute_eff_diff(domain)

  problem.d_eff = d_eff


##
## MAIN
##



if __name__ == "__main__":
  msg="""
\nPurpose: 
  run homogenized problem 
 
Usage:
  homog.py -case myofilament/globular/validationSphere/custom/run <-smol> <-molGamer>
           <-molPrefix molPrefix> <-validation all/etc> <-boundaryTol float>
           or 
           <-file filename.xml>

  where 
    -smol - run molecular domain with electrostatics
    -molGamer - molecule was prepared with Gamer

Notes:
"""

  #GoelEx2p7()

  import sys
  if len(sys.argv) < 2:
      raise RuntimeError(msg)
  
  if MPI.rank(mpi_comm_world())==0:
    print "Writing outputs to %s" % outname 
  #print "WARNING: files arent writing!"

  for i,arg in enumerate(sys.argv):
    if(arg=="-smol"):
      print "Ensabling electrostatics contribution" 
      smolMode = True

    if(arg=="-case"):
      case = sys.argv[i+1]

    if(arg=="-molPrefix"):
      molPrefix = sys.argv[i+1]

    if(arg=="-molGamer"):
      molGamer = 1

    if(arg=="-validation"):
      validationMode=sys.argv[i+1]

    if(arg=="-file"):
      fileName = sys.argv[i+1]

    if(arg=="-boundaryTol"):
      boundaryTolerance=float(sys.argv[i+1])
      print "Using boundary tolerance of %f" % boundaryTolerance


  #
  # Validation 
  #  
  if(validationMode!=0):          
    ValidationPaper(mode=validationMode)
    print "Run paraview.py to view results (on mac/home computer)"
    quit()


  #
  # Normal 
  #  


  #Debug2()
  #quit()

  # ovoerride
  #debug =1
  #cellPrefix = ""
  #molPrefix = "molecular_volmesh"
  #root = "/home/huskeypm/bin/grids/"
  
  #case = "globular"
  #case = "myofilament"
  #case = "validation"


  # globular case
  if(case=="globular"):
    #cellPrefix="mol/cell"
    #wholeCellPrefix="mol/multi_clustered"
    #molPrefix="120529_homog/1CID"
    molPrefix = "mol/1CID"
    SolveHomogSystem(debug=debug,\
      root=root,\
      cellPrefix=cellPrefix, molPrefix=molPrefix,wholeCellPrefix=wholeCellPrefix,\
      molGamer=1)
  
  # TnC/cylindrical case
  elif(case=="myofilament"):
    cellPrefix="mol/cell"
    wholeCellPrefix="mol/multi_clustered"
    molPrefix = "mol/troponin"
    SolveHomogSystem(debug=debug,\
      root=root,\
      cellPrefix=cellPrefix, molPrefix=molPrefix,wholeCellPrefix=wholeCellPrefix,\
      molGamer=molGamer,option="troponin")


  elif(case=="validationSphere"):
    ValidationSphere()

  elif(case=="custom"):
    root = "./"
    results = SolveHomogSystem(debug=debug,\
      root=root,\
      cellPrefix=cellPrefix, molPrefix=molPrefix,wholeCellPrefix=wholeCellPrefix,\
      smolMode = smolMode,
      molGamer=molGamer)
    r=results.molDomUnit.problem.d_eff    
    #r = r/np.linalg.norm(r)
    #print r
  elif(case=="run"):
    root = "./"
    molPrefix = fileName
    results = SolveHomogSystem(debug=debug,\
      root=root,\
      addSuffix=False,
      cellPrefix=cellPrefix, molPrefix=molPrefix,\
      wholeCellPrefix=wholeCellPrefix,\
      smolMode = smolMode,
      molGamer=molGamer)
    r=results.molDomUnit.problem.d_eff


  else:
    msg = "Case " + case + " not understood"   
    raise RuntimeError(msg)


