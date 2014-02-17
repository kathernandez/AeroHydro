import numpy as np
import matplotlib.pyplot as plt 
from math import * 

N = 200  #Number of points in each direction 
xstart, xend = -4.0, 4.0    #x-dir boundaries
ystart, yend = -2.0, 2.0    #y-dir boundaries
x = np.linspace(xstart, xend, N) #x 1-D array 
y = np.linspace(ystart, xend, N) #y 1-D array 
X,Y = np.meshgrid(x,y)  #generation of grid 
np.shape(X)

uinf = 1.0 #free stream speed
a_degrees = 0.0 #AOA (in degrees)
a = a_degrees*pi/180 

#computing the velocity components on the mesh grid 
ufree=uinf*cos(a)*np.ones((N,N), dtype = float)
vfree=uinf*sin(a)*np.ones((N,N), dtype = float)

#computing the stream function on the mesh grid
psifree = + uinf*cos(a)*Y - uinf*sin(a)*X

#function to compute the velocity field of a source/sink 
def getvelocity (strength,xs,ys,X,Y):
    u=strength/(2*pi)*(X-xs)/((X-xs)**2+(Y-ys)**2)
    v=strength/(2*pi)*(Y-ys)/((X-xs)**2+(Y-ys)**2)
    return u,v

#function to compute the stream-function of a source/sink 
def getstreamfunction(strength,xs,ys,X,Y):
    psi=strength/(2*pi)*np.arctan2((Y-ys),(X-xs))
    return psi

strengthsource = 5.0  #strength of the source 
xsource, ysource = -1.0, 0.0

#computing the velocity components 
usource, vsource = getvelocity(strengthsource, xsource, ysource, X, Y)

#computing the stream function 
psisource = getstreamfunction(strengthsource, xsource, ysource, X, Y)

#print(usource)
#print(vsource)
#print(psisource)

#superposition of the source on the free stream 
u = ufree + usource
v = vfree + vsource
psi= psifree + psisource

#plotting
size = 10
plt.figure(figsize=(size, (yend-ystart)/(xend-xstart) * size))
plt.grid(True)
plt.xlabel('x', fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xstart,xend)
plt.ylim(ystart,yend)
plt.streamplot(X,Y,u,v,density=2.0,linewidth=1,arrowsize=1,arrowstyle='->')
plt.scatter(xsource, ysource, c='#CD2305', s=80, marker = 'o')

#computing the stagnation point
xstag = xsource - strengthsource*cos(a)/(2*pi*uinf)
ystag = ysource - strengthsource*sin(a)/(2*pi*uinf)

#adding the stagnation point to the figure 
plt.scatter(xstag, ystag, c='b', s=80, marker='o')

#adding the dividing line to the figure 
if (a==0.0):
        plt.contour(X,Y,psi,\
        levels = [-strengthsource/2, +strengthsource/2],\
        colors='#CD2305',linewidths=2,linestyles='solid')
        
#plt.show()

#superimpose source and sink in freestream
strengthsink=-5.0 #strength of the sink 
xsink, ysink=1.0,0.0

#computing the velocity feild on the mesh grid 
usink, vsink = getvelocity(strengthsink, xsink, ysink, X,Y)

#computing the stream-function on the grid mesh 
psisink = getstreamfunction(strengthsink, xsink, ysink, X,Y)

u=ufree + usource +usink
v=vfree + vsource +vsink

psi=psifree+psisource+psisink

#plotting 
size = 10
plt.figure(figsize =(size,(yend-ystart)/(xend-xstart)*size))
plt.xlabel('x', fontsize = 16)
plt.ylabel('y', fontsize = 16)
plt.xlim(xstart,xend)
plt.ylim(ystart, yend) 
plt.streamplot(X,Y,u,v,density = 2.0, linewidth=1,arrowsize=1,arrowstyle='->')
plt.scatter([xsource, xsink],[ysource,ysink],c='m',s=80,marker='o')
if (a==0.0):
    plt.contour(X,Y,psi,\
    levels=[0.0],colors='m',linewidths=2,linestyles='solid')
    
#plt.show()

#computing the pressure coefficient 

Cp=1.0-(u**2+v**2)/uinf**2

#Plotting 
size=10
plt.figure(figsize=(1.1*size,(yend-ystart)/(xend-xstart)*size))
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xstart,xend)
plt.ylim(ystart,yend)

contf=plt.contourf(X,Y,Cp,levels=np.linspace(-2.0,1.0,100),extend='both')
cbar = plt.colorbar(contf)
cbar.set_label(r'$C_p$',fontsize=16)
cbar.set_ticks([-2.0,-1.0,0.0,1.0])
plt.scatter([xsource,xsink],[ysource,ysink],c='r',s=80,marker='o')
plt.contour(X,Y,psi,levels=[0.0],colors='r',linewidths=2,linestyles='solid')

plt.show()
    