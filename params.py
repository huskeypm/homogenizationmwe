"""
Overkill, but this allows for sharing params between routines 
"""


class params:
  def __init__(self,d=1.): # Be sure to use float 
    self.d = d # [um^2/s] ??  
#    self.concInitial = 0.1 # [uM]
#    #self.ANG_TO_UM = 1e-4
#    #print "WARNING: removed ANG TO UM" 
#    self.Ang_to_um = 1     
#    self.dt = 0.001
#    self.tStep = 20
#    self.beta = 1/0.593 # [kcal/mol]
#
#from params import *
parms = params()

