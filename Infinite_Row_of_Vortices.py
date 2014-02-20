#First, we import python libraries 
import numpy as np
import matplotlib.pyplot as plt
from math import *

#Create a grid mesh 

N=50                          # number of grid points 
xstart,xend = -2.0,2.0        # x boundaries
ystart, yend = -1.0,1.0       # y boudaries 
x=np.linspace(xstart,xend,N)  # x direnction 1D array 
y=np.linspace(ystart,yend,N)  # y direction 1D array 
X,Y=np.meshgrid(x,y)          # mesh grid created 

gamma = 5.0                   #strength of vortex (same for all) 

#define funtion to compute the velocity components
def getVelocityVortex(strength,xv,yv,X,Y):
    u = + strength/(2*pi)*(Y-yv)/((X-xv)**2+(Y-yv)**2)
    v = - strength/(2*pi)*(X-xv)/((X-xv)**2+(Y-yv)**2)
    return u,v

#define function to compute the stream function     
def getStreamFunctionVortex(strength,xv,yv,X,Y):
    psi = strength/(4*pi)*np.log((X-xv)**2+(Y-yv)**2)
    return psi 


#--------Finite row of vortices---------------------------------


#we need to initialize u,v, and psi by assigning 0 arrays 
#Note, np arrays will be used instead of lists because arrays are faster for..
#..computation and require less memory allocation 
uvortex=np.zeros((N,N),dtype=float)
vvortex=np.zeros((N,N),dtype=float) 
psivortex=np.zeros((N,N),dtype=float)

#then we should assign a fixed y-coordinate since its a row only..
#..the xcoordinate will vary but will be spaced evenly 
numbervortices = 8                               #number of vortices to be plotted 
xvortex=np.linspace(xstart,xend,numbervortices)  #location of vortices in the x-direction
yvortex = np.zeros(numbervortices)               #location of vortices in the y-direction  

#apply a for loop to calculate u,v, and psi for every ith vortex 
for i in range(0,numbervortices):
    uvortex_init,vvortex_init = getVelocityVortex(gamma,xvortex[i],yvortex[i],X,Y)
    uvortex= uvortex + uvortex_init
    vvortex= vvortex + vvortex_init
    
    psivortex_init= getStreamFunctionVortex(gamma,xvortex[i],yvortex[i],X,Y)
    psivortex= psivortex + psivortex_init

#plot row of vortices 
size = 10
plt.figure(figsize=(size,(yend-ystart)/(xend-xstart)*size))
plt.xlim(xstart, xend)
plt.ylim(ystart, yend)
plt.xlabel('x', fontsize=18)
plt.ylabel('y', fontsize=18)
plt.streamplot(X,Y,uvortex,vvortex,\
              density=2.0,linewidth=1,arrowsize=1,arrowstyle='->')
plt.scatter(xvortex,yvortex,c='r',s=80,marker = 'o')

plt.show()
        






