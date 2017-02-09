
import numpy as np
import homoglight as hl 
import numpy as np
from dolfin import *

def test(): 
  fileXML="tests/squarehole_in_square.xml"  
  problem = hl.runHomog(fileXML,verbose=True)
  print problem.d_eff     # should have smallest volume fraction, hence smallest deff  

  # these two should ahve comparable deffs (occupy about the same vol frac of nearly 1.0)
  fileXML="tests/square_in_square.xml"  
  problem = hl.runHomog(fileXML,verbose=True)
  print problem.d_eff 
  fileXML="tests/squares_in_square.xml"
  problem = hl.runHomog(fileXML,verbose=True)
  print problem.d_eff 

def run(fileXML="tests/volFrac_0.50.xml.gz",
  boxMin=None,boxMax=None):
 
  # Run code  
  mesh = Mesh(fileXML)
  V = FunctionSpace(mesh,"CG",1)
  print "vol",assemble(Constant(1.)*dx(domain=mesh))


  problem = hl.runHomog(fileXML,verbose=True,
    boxMin=boxMin,boxMax=boxMax)
 

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
    # calls 'run' with the next argument following the argument '-validation'
    if(arg=="-validation"):
      run()      
      quit()
    elif(arg=="-run"): 
      run(fileXML=sys.argv[i+1])
      quit()
    #expects
# python validate.py -runUnitCell 0,0,0,50,50,50
    elif(arg=="-runUnitCell"):
      mystr    = sys.argv[i+1]
      #mystr = "0,0,0,50,50,50"
      myspl = mystr.split(',')
      myspl = [np.float(i) for i in myspl]
      boxMin= np.asarray(myspl[0:3])
      boxMax= np.asarray(myspl[3:6])
      print boxMin,boxMax
      #boxMin=np.asarray([0,0,0])
      #boxMax=np.asarray([50.,50.,50.])
      run(boxMin=boxMin,boxMax=boxMax)
      quit()
    elif(arg=="-test"): 
      test()
      quit()
  





  raise RuntimeError("Arguments not understood")




