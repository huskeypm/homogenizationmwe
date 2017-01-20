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


class params:
  def __init__(self,d=1.): # Be sure to use float 
    self.d = d # [um^2/s] ??  
    self.concInitial = 0.1 # [uM]
    #self.ANG_TO_UM = 1e-4
    #print "WARNING: removed ANG TO UM" 
    self.Ang_to_um = 1     
    self.dt = 0.001
    self.tStep = 20
    self.beta = 1/0.593 # [kcal/mol]

#from params import *
parms = params()

