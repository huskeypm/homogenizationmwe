
import numpy as np
import homoglight as hl 
import numpy as np
from dolfin import *

def run(fileXML="tests/volFrac_0.50.xml.gz",
  boxMin=None,boxMax=None):
 
  # Run code  

  problem = hl.runHomog(fileXML,verbose=True,
    boxMin=boxMin,boxMax=boxMax)

# Run series of calculations using larger boxes around 
# pore (x,y increasing, z fixed) 
def runSeries(fileXML=None): 
  mesh = Mesh(fileXML)
  V = FunctionSpace(mesh,"CG",1)
  c = mesh.coordinates()
 
  print "Assuming grid is centered at 000"
  daMins = np.min(c,axis=0)
  daMaxs = np.max(c,axis=0)

  incr = 0.005

  volFracs = []
  zDeffs = []
  nCases = 20
  for i in range(nCases):  
    boxMin = daMins - (incr*i)*np.array([1,1,0])
    boxMax = daMaxs + (incr*i)*np.array([1,1,0])
    problem = hl.runHomog(fileXML,verbose=True,
      boxMin=boxMin,boxMax=boxMax)
    #print boxMin, boxMax
    #print problem.volFrac
    #print problem.d_eff   
    volFracs.append(problem.volFrac)
    zDeffs.append(problem.d_eff[2]) 

  # plot 
  import matplotlib.pylab as plt
  plt.plot(zDeffs,volFracs)
  plt.scatter(zDeffs,volFracs)
  plt.ylabel("zDeff")
  plt.xlabel("volFrac") 
  plt.gcf().savefig("cylinder.png",dpi=150) 

 

#
# Message printed when program run without arguments 
#
def helpmsg():
  scriptName= sys.argv[0]
  msg="""
Purpose: 
  Very simple consistency check 
 
Usage:
"""
  msg+="  %s -validation" % (scriptName)
  msg+="""
  
 
Notes:

"""
  return msg

#
# MAIN routine executed when launching this script from command line 
#
if __name__ == "__main__":
  import sys
  msg = helpmsg()

  if len(sys.argv) < 2:
      raise RuntimeError(msg)

  #fileIn= sys.argv[1]
  #if(len(sys.argv)==3):
  #  1
  #  #print "arg"

  # Loops over each argument in the command line 
  for i,arg in enumerate(sys.argv):
    if(arg=="-runUnitCell"):
      fileXML = "tests/cylinder.xml"             
      run(fileXML=fileXML)
      quit()
    elif(arg=="-runSeries"):   
      fileXML = "tests/cylinder.xml"             
      runSeries(fileXML=fileXML)
      quit()
  





  raise RuntimeError("Arguments not understood")




