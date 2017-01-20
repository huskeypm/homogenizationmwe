
import homoglight as hl 
import numpy as np


def doit():                 
  vf = 0.50
  # A 3D mesh with an inclusion that occupies a volume fraction of 0.5 
  fileXML="tests/volFrac_0.50.xml.gz"
 
  # Run code  
  problem = hl.runHomog(fileXML,verbose=True)

  # compare against hashin-shtrikman bound for sphere w 0.5 volume fraction 
  deffHSBound = 2*vf/(3-vf)
  assert(np.abs(problem.d_eff[0]-deffHSBound)<0.01), "Don't commit! somthing changed"
  print "All is ok!"
  quit()

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
    # calls 'doit' with the next argument following the argument '-validation'
    if(arg=="-validation"):
      doit()      
      quit()
  





  raise RuntimeError("Arguments not understood")




