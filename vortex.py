import numpy as np
import matplotlib.pyplot as plt
from math import *

N=50  #number of points in each direction
xstart, xend = -2.0,2.0 #x direction boundaries 
ystart, yend = -1.0,1.0
x = np.linspace(xstart,xend,N) # x 1D array 
y = np.linspace(ystart,yend,N) # y 1D array 
X,Y = np.meshgrid(x,y) #generation of the mesh grid 

gamma = 5.0 #vortex strength 
xvortex,yvortex = 0.0, 0.0 #location of vortex 

#function to compute the velocity components of the vortex
def getVelocityVortex(strength,xv,yv,X,Y):
    u = + strength/(2*pi)*(Y-yv)/((X-xv)**2+(Y-yv)**2)
    v = - strength/(2*pi)*(X-xv)/((X-xv)**2+(Y-yv)**2)
    return u,v

#function to compute the stream function of a vortex
def getStreamFunctionVortex(strength,xv,yv,X,Y):
    psi=strength/(4*pi)*np.log((X-xv)**2+(Y-yv)**2)
    return psi

#computing the velocity components on the mesh grid 
uvortex, vvortex = getVelocityVortex(gamma,xvortex,yvortex,X,Y)

#computing the stream function on the mesh grid 
psivortex = getStreamFunctionVortex(gamma,xvortex,yvortex,X,Y)

#plotting 
size = 10
plt.figure(figsize=(size,(yend-ystart)/(xend-xstart)*size))
plt.xlabel('x',fontsize=16)
plt.ylabel ('y', fontsize=16)
plt.xlim(xstart,xend)
plt.ylim(ystart, yend)
plt.streamplot(X,Y,uvortex,vvortex, \
            density = 2.0, linewidth=1, arrowsize = 1, arrowstyle='->')
plt.scatter(xvortex,yvortex,c='r', s=80, marker = 'o');
  

strengthsink= -1.0 #strength of the sink 
xsink,ysink= 0.0, 0.0 #location of the sink 

#function to compute the velocity components of a sink
def getVelocitySink(strength, xs, ys, X,Y):
    u=strength/(2*pi)*(X-xs)/((X-xs)**2+(Y-ys)**2)
    v=strength/(2*pi)*(Y-ys)/((X-xs)**2+(Y-ys)**2)
    return u,v

#function to compute the stream function on the mesh grid 
def getStreamFunctionSink(strength, xs, ys, X, Y):
    psi = strengthsink/(2*pi)*np.arctan2((Y-ys),(X-xs))
    return psi

#computing the velocity components of the mesh grid
usink,vsink=getVelocitySink(strengthsink,xsink,ysink,X,Y)

#computing the stream function
psisink = getStreamFunctionSink(strengthsink, xsink, ysink, X,Y)

#superposition of the sink and the vortex
u=uvortex+usink
v=vvortex+vsink
psi=psivortex+psisink

#plotting
size=10
plt.figure(figsize=(size, (yend-ystart)/(xend-xstart)*size))
plt.xlabel('x', fontsize = 16)
plt.ylabel('y', fontsize = 16)
plt.xlim(xstart, xend)
plt.ylim(ystart, yend)
plt.streamplot(X,Y,u,v,\
        density = 2.0, linewidth=1, arrowsize =1, arrowstyle = '->')
plt.scatter(xvortex,yvortex,c='r',s=80,marker='o');
plt.show()