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
plt.show()

        


