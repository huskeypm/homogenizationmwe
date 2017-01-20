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


# bcs
def boundary(x,on_boundary):
  return on_boundary

class DefaultUnitDomain(Domain):
  def __init__(self,type="field"):
    super(DefaultUnitDomain,self).__init__(type)

    problem = self.problem
    problem.name = "Default"

  def Setup(self):
    # mesh
    problem = self.problem
    #problem.mesh = UnitCube(8,8,8)
    problem.mesh = UnitSphere(8)
    if(self.type=="scalar"):
        problem.V = FunctionSpace(problem.mesh,"CG",1)
    elif(self.type=="field"):
        problem.V = VectorFunctionSpace(problem.mesh,"CG",1)

    utilObj=self.utilObj
    utilObj.GeometryInitializations()
    utilObj.DefinePBCMappings()

    # geom
    self.CalcGeom(problem)


