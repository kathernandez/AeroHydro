# Lesson 8 Source sheet 
import numpy as np
import matplotlib.pyplot as plt
from math import * 

N = 100                          #number of points in each direction 
xstart,xend = -1.0,1.0           #x-dir boundaries
ystart, yend = -1.5,1.5          #y-dir boundaries 
x = np.linspace(xstart,xend,N)   #x 1D array 
y = np.linspace(ystart,yend,N)  #y 1D array 
X,Y = np.meshgrid(x,y)           #generate meshgrid

uinf = 1.0         #free stream velocity
ufreestream = uinf*np.ones((N,N),dtype=float) 
vfreestream = np.zeros((N,N),dtype=float)


#---------------------finite number of sources along a line---------------------
class Source: 
    def __init__(self, strength, x, y): 
        self.strength = strength
        self.x,self.y = x,y
    #get the velocity field
    def velocity(self,X,Y):
        self.u = self.strength/(2*pi)*(X-self.x)/((X-self.x)**2+(Y-self.y)**2)
        self.v = self.strength/(2*pi)*(Y-self.y)/((X-self.x)**2+(Y-self.y)**2)
    #get the stream function
    def streamfunction(self,X,Y):
        self.psi = self.strength/(2*pi)*np.arctan2((Y-self.y),(X-self.x))
    
#definition of sources
Nsources = 11
strength = 5.0
strengthsource = strength/Nsources
xsource = 0.0
ysource = np.linspace(-1,1,Nsources)

#creation of the source line
sources = np.empty(Nsources,dtype=object)

for i in range(Nsources):
    sources[i] = Source(strengthsource,xsource,ysource[i])
    sources[i].velocity(X,Y)

#superposition
u = ufreestream.copy()
v = vfreestream.copy()

for s in sources:
    u = np.add(u,s.u)
    v = np.add(v,s.v)

#plotting 

size = 6
plt.figure(figsize=(size,(yend-ystart)/(xend-xstart)*size))
plt.grid(True)
plt.xlabel('x', fontsize = 16)
plt.ylabel('y', fontsize = 16) 
plt.xlim(xstart,xend)
plt.ylim(ystart,yend)

plt.streamplot(X,Y,u,v,\
              density=2,linewidth=1,arrowsize=1,arrowstyle='->')
plt.scatter(xsource*np.ones(Nsources,dtype=float),ysource,c='r',s=80,marker='o')
velocity = plt.contourf(X,Y,np.sqrt(u**2+v**2), levels = np.linspace(0.0,0.1,10))
cbar = plt.colorbar(velocity,ticks=[0.0,0.05,0.1],orientation = 'horizontal')
cbar.set_label('Velocity magnitude',fontsize=16);
plt.show()

#-----------infinite line of source-------------------------------

from scipy import integrate

sigma = 2.5 #strength of source sheet

upanel = np.empty((N,N),dtype=float)
vpanel = np.empty((N,N),dtype=float)

#boundaries of the sheet 
ymin,ymax = -1.0,1.0

#Computing the velocity field 
for i in range(N):
    for j in range(N):
        
        integrand = lambda s: X[i,j]/(X[i,j]**2+(Y[i,j]-s)**2)
        upanel[i,j] = sigma/(2*pi)*integrate.quad(integrand,ymin,ymax)[0]
        
        integrand = lambda s: (Y[i,j]-s)/(X[i,j]**2+(Y[i,j]-s)**2)
        vpanel[i,j] = sigma/(2*pi)*integrate.quad(integrand,ymin,ymax)[0]

#superposition to the uniform flow 
u2 = np.add(ufreestream,upanel)
v2 = np.add(vfreestream,vpanel)

# plotting

size = 8
plt.figure(figsize=(size,(yend-ystart)/(xend-xstart)*size))
plt.grid(True)
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xstart,xend)
plt.ylim(ystart,yend)
plt.streamplot(X,Y,u2,v2,\
              density=2, linewidth=1,arrowsize=1,arrowstyle='->')
plt.axvline( 0.0,(ymin-ystart)/(yend-ystart),(ymax-ystart)/(yend-ystart),\
            color='r',linewidth=4)
velocity = plt.contourf(X,Y,np.sqrt(u2**2+v2**2),levels=np.linspace(0.0,0.1,10))
cbar=plt.colorbar(velocity,ticks=[0.0,0.05,0.1],orientation='horizontal')
cbar.set_label('Velocity Magnitude',fontsize=16)
plt.show()



