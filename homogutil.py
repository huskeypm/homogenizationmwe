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
# Revisions
#   121009 generalized to 2d
#
from dolfin import *
import numpy as np

###
### Class for dealing with various geometric issues (bounds, pairing between vertices) 
### 

class homogutil: 
  def __init__(self,problem=0):
    self.prob = problem 

  def CalcBounds(self,mesh):
    prob = self.prob
  
    V = FunctionSpace(mesh, "Lagrange", 1)
    prob.nDims = np.shape(mesh.coordinates())[1]
    boundsMin=np.zeros(prob.nDims)
    boundsMax=np.zeros(prob.nDims)

    # need to use 'gather' to poll all the CPUs for their parts of the mesh
    for i in range(prob.nDims):
      tag = "x[%d]"%i
      u = interpolate(Expression(tag), V)
      x = Vector()
      u.vector().gather(x, np.array(range(V.dim()), "intc"))
      boundsMin[i] =np.min(x.array())
      boundsMax[i] =np.max(x.array())

    prob.boundsMin = boundsMin
    prob.boundsMax = boundsMax

    #return (boundsMin,boundIdxMin,boundsMax,boundIdxMax)
    return (boundsMin,-1,boundsMax,-1)
  
  def CalcMidpoint(self,mesh):
    (boundsMin,boundIdxMin,boundsMax,boundIdxMax) = self.CalcBounds(mesh)
    return (boundsMin + boundsMax)/2.
  
  def CalcRanges(self,mesh):
    (boundsMin,boundIdxMin,boundsMax,boundIdxMax) = self.CalcBounds(mesh)
    prob = self.prob
  
   
    ranges = np.zeros(prob.nDims)
    ranges[0] = boundsMax[0]-boundsMin[0]
    ranges[1] = boundsMax[1]-boundsMin[1]
    ranges[2] = boundsMax[2]-boundsMin[2]
  
    return ranges
  
  def CenterMesh(self,mesh):
  
    (boundsMin,boundIdxMin,boundsMax,boundIdxMax) = self.CalcBounds(mesh)
    mp = self.CalcMidpoint(mesh)
    #print CalcMidpoint(mesh)

    center = 1
    if(center):
      for i,c in enumerate( mesh.coordinates() ):
        c -= mp
        mesh.coordinates()[i] = c
  
  def DefinePBCMappings(self):
    prob = self.prob 
    mesh = prob.mesh
    (boundsMin,boundIdxMin,boundsMax,boundIdxMax) = self.CalcBounds(mesh)

    # x component 
    #v0x= np.array(([-1,0,0]))
    #v1x= np.array(([1,0,0]))
    v0x = mesh.coordinates()[boundIdxMin[0],:]
    v1x = mesh.coordinates()[boundIdxMax[0],:]
    prob.targetsx = [tuple(v1x)]
    prob.vert_mapx = {}
    prob.vert_mapx[tuple(v0x)] = tuple(v1x)
    
    # y component 
    v0y = mesh.coordinates()[boundIdxMin[1],:]
    v1y = mesh.coordinates()[boundIdxMax[1],:]
    #v0y= np.array(([0,-1,0]))
    #v1y= np.array(([0,1,0]))
    prob.targetsy = [tuple(v1y)]
    prob.vert_mapy = {}
    prob.vert_mapy[tuple(v0y)] = tuple(v1y)

    print v0x, " maps to " , prob.vert_mapx[ tuple(v0x)  ]
    print v0y, " maps to " , prob.vert_mapy[ tuple(v0y)  ]
    
    # z component 
    #v0z= np.array(([0,0,-1]))
    #v1z= np.array(([0,0,1]))
    if(prob.nDims==3):
      v0z = mesh.coordinates()[boundIdxMin[2],:]
      v1z = mesh.coordinates()[boundIdxMax[2],:]
      prob.targetsz= [tuple(v1z)]
      prob.vert_mapz = {}
      prob.vert_mapz[tuple(v0z)] = tuple(v1z)
      print v0z, " maps to " , prob.vert_mapz[ tuple(v0z)  ]
  

  def GeometryInitializations(self):
    prob = self.prob 
    mesh = prob.mesh
  
    # center
    self.CenterMesh(mesh)
  
    prob.center_coord = self.CalcMidpoint(mesh)
    #coords = mesh.coordinates()
  #
  #
    coords = mesh.coordinates()
    prob.dims = np.array([coords[:,i].max()-coords[:,i].min() for i in range(prob.nDims)])
  #  print "Center:", center_coord
  #  print "Dim:", dims
