#lesson 7 Aerodymanics Interferences 
#and method of images

import numpy as np
import matplotlib.pyplot as plt
from math import *
from IPython.core.display import clear_output

N = 50     #number of point in each direction 
xstart,xend = -2.0, 2.0  #x-direction boundaries
ystart,yend = -1.0, 1.0  #y-direction boundaries
x = np.linspace(xstart, xend, N)  #X 1D array 
y = np.linspace(ystart, yend, N)  #y 1D array
X,Y = np.meshgrid(x,y)            #generate mesh

#--------------------source near a wall -----------------------------
class Source: 
    def __init__(self,strength,x,y):
        self.strength = strength
        self.x,self.y = x,y
    #get the velocity field 
    def velocity(self,X,Y):
        self.u = self.strength/(2*pi)*(X-self.x)/((X-self.x)**2+(Y-self.y)**2)
        self.v = self.strength/(2*pi)*(Y-self.y)/((X-self.x)**2+(Y-self.y)**2)
    #get the stream function
    def streamfunction(self,X,Y):
        self.psi = self.strength/(2*pi)*np.arctan2((Y-self.y),(X-self.x))

strengthsource = 1.0
xsource,ysource = 0.0,0.5
source = Source(strengthsource,xsource,ysource)
source.velocity(X,Y)
source.streamfunction(X,Y)

sourceimage = Source(strengthsource,xsource,-ysource)
sourceimage.velocity(X,Y)
sourceimage.streamfunction(X,Y)

#apply superposition of the source and its image
u = source.u + sourceimage.u
v = source.v + sourceimage.v
psi = source.psi + sourceimage.psi

#plotting 
size = 10
plt.figure(num=0,figsize=(size,(yend-ystart)/(xend-xstart)*size))
plt.xlabel('x', fontsize = 16)
plt.ylabel('y', fontsize = 16)
plt.xlim(xstart,xend)
plt.ylim(ystart,yend)
plt.streamplot(X,Y,u,v,\
               density=2.0,linewidth=1,arrowsize=1,arrowstyle='->')
plt.scatter(source.x,source.y,c='r',s=80,marker='o')
plt.scatter(sourceimage.x,sourceimage.y,c='r',s=80,marker='D')
plt.axhline(0.0,color='k',linestyle='--',linewidth=4)
#plt.show()

#--------------------------vortex near a wall------------------

class Vortex:
    def __init__(self,strength,x,y):
        self.strength = strength
        self.x,self.y = x,y
    #get velocity field
    def velocity(self,X,Y):
        self.u = +self.strength/(2*pi)*(Y-self.y)/((X-self.x)**2+(Y-self.y)**2)
        self.v = -self.strength/(2*pi)*(X-self.x)/((X-self.x)**2+(Y-self.y)**2)
    #get stream function
    def streamfunction(self,X,Y):
        self.psi = -self.strength/(4*pi)*np.log((X-self.x)**2+(Y-self.y)**2)

strengthvortex = 1.0
xvortex,yvortex = 0.0,0.5

vortex = Vortex(strengthvortex,xvortex,yvortex)
vortex.velocity(X,Y)
vortex.streamfunction(X,Y)

vorteximage=Vortex(-strengthvortex,xvortex,-yvortex)
vorteximage.velocity(X,Y)
vorteximage.streamfunction(X,Y)

#apply superposition
u = vortex.u + vorteximage.u
v = vortex.v + vorteximage.v
psi = vortex.psi + vorteximage.psi

#plotting
size = 10
plt.figure(num=1,figsize=(size,(yend-ystart)/(xend-xstart)*size))
plt.xlabel('x', fontsize =16)
plt.ylabel('y', fontsize = 16)
plt.xlim(xstart,xend)
plt.ylim(ystart,yend)
plt.streamplot(X,Y,u,v,\
        density=2.0, linewidth=1,arrowsize=1,arrowstyle='->')
plt.scatter(vortex.x,vortex.y, c='m',s=80,marker='o')
plt.scatter(vorteximage.x,vorteximage.y,c='m',s=80,marker='D')
plt.axhline(0.0,color='k',linestyle='--',linewidth=4)
plt.show()
#----------------motion of vortex pair near the ground------------
strengthvortex = 1.0
xvortex1,yvortex1 = -0.1,0.5
xvortex2,yvortex2 = 0.1,0.5

vortex1 = Vortex(+strengthvortex,xvortex1,yvortex1)
vortex2 = Vortex(-strengthvortex,xvortex2,yvortex2)

vortex1.velocity(X,Y)
vortex1.streamfunction(X,Y)
vortex2.velocity(X,Y)
vortex2.streamfunction(X,Y)

vorteximage1 = Vortex(-strengthvortex,xvortex1,-yvortex1)
vorteximage2 = Vortex(+strengthvortex,xvortex2,-yvortex2)

vorteximage1.velocity(X,Y)
vorteximage2.velocity(X,Y)
vorteximage1.streamfunction(X,Y)
vorteximage2.streamfunction(X,Y)

#superposition ofvortex pair and image 
u = vortex1.u + vorteximage1.u + vortex2.u + vorteximage2.u
v = vortex1.v + vorteximage1.v + vortex2.v + vorteximage2.v
psi = vortex1.psi + vorteximage1.psi + vortex2.psi + vorteximage2.psi

#plotting
size = 10
plt.figure(num=3,figsize=(size, (yend-ystart)/(xend-xstart)*size))
plt.xlabel('x', fontsize = 16)
plt.ylabel('y', fontsize = 16)
plt.xlim(xstart,xend)
plt.ylim(ystart,yend)
plt.streamplot(X,Y,u,v,\
        density = 2.0, linewidth=1,arrowsize=1,arrowstyle='->')
plt.scatter(vortex1.x,vortex1.y,c='r',s=80,marker='o')
plt.scatter(vortex2.x,vortex2.y,c='g',s=80,marker='o')
plt.scatter(vorteximage1.x,vorteximage1.y, c='r',s=80,marker='D')
plt.scatter(vorteximage2.x,vorteximage2.y, c='g',s=80,marker='D')
plt.axhline(0.0,color='k',linestyle='--',linewidth=4)
plt.show()



#------------doublet near a wall --------------------------
uinf = 0.0
ufreestream = uinf*np.ones((N,N),dtype=float)
vfreestream = np.zeros((N,N), dtype =float) 
psifreestream = uinf*Y

class Doublet:
    def __init__(self,strength,x,y):
        self.strength = strength
        self.x,self.y = x,y
    #Get velocity field 
    def velocity(self, X,Y):
        self.u = -self.strength/(2*pi)*((X-self.x)**2-(Y-self.y)**2)/((X-self.x)**2+(Y-self.y)**2)**2
        self.v = -self.strength/(2*pi)*2*(X-self.x)*(Y-self.y)/((X-self.x)**2+(Y-self.y)**2)**2
    #get stream function 
    def streamfunction(self,X,Y):
        self.psi = -self.strength/(2*pi)*(Y-self.y)/((X-self.x)**2+(Y-self.y)**2)

strengthdoublet = 1.0
xdoublet,ydoublet = 0.0,0.3
doublet = Doublet(strengthdoublet,xdoublet,ydoublet)
doublet.velocity(X,Y)
doublet.streamfunction(X,Y)

doubletimage = Doublet(strengthdoublet,xdoublet,-ydoublet)
doubletimage.velocity(X,Y)
doubletimage.streamfunction(X,Y)

#superposition of the doublet and its image to the uniform flow 
u = ufreestream + doublet.u + doubletimage.u
v = vfreestream + doublet.v + doubletimage.v
psi = psifreestream + doublet.psi + doubletimage.psi 

#plotting
size = 10
plt.figure(num=2, figsize = (size,(yend-ystart)/(xend-xstart)*size))
plt.xlabel('x', fontsize=16)
plt.ylabel('y', fontsize=16)
plt.xlim(xstart,xend)
plt.ylim(ystart,yend)
plt.streamplot(X,Y,u,v,\
            density=2.0,linewidth=1,arrowsize=1,arrowstyle='->')
plt.scatter(doublet.x,doublet.y,c='r',s=80, marker='o')
plt.scatter(doubletimage.x,doubletimage.y, c='r', s=80, marker='D')
plt.axhline(0.0,color='k',linestyle='--',linewidth=4)
plt.show()




        


