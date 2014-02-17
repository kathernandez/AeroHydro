import numpy as np
import matplotlib.pyplot as plt
from math import *


N=50 #number of point in each direction
xstart, xend=-2.0,2.0
ystart, yend=-1.0,1.0 
x = np.linspace(xstart,xend,N)
y=np.linspace(ystart, yend, N)
X,Y=np.meshgrid(x,y)

k=1.0 #strength of the doublet
xdoublet, ydoublet = 0.0, 0.0 #location of the doublet

#function to compute the velocity components of a doublet
def getVelocityDoublet(strength, xd, yd, X, Y):
    u = -strength/(2*pi)*((X-xd)**2-(Y-yd)**2)/((X-xd)**2+(Y-yd)**2)**2
    v = -strength/(2*pi)*2*(X-xd)*(Y-yd)/((X-xd)**2+(Y-yd)**2)**2
    return u,v
    
    #function to compute the stream function of doublet
def getStreamFunctionDoublet(strength,xd,yd,X,Y):
    psi = -strength/(2*pi)*(Y-yd)/((X-xd)**2+(Y-yd)**2)
    return psi
        
#computing the velocity components on the mesh grid
udoublet, vdoublet = getVelocityDoublet(k,xdoublet,ydoublet,X,Y)

#computing the stream-function on the mesh grid
psidoublet=getStreamFunctionDoublet(k, xdoublet, ydoublet, X,Y)

#plotting 
size = 10
plt.figure(figsize=(size, (yend-ystart)/(xend-yend)*size))
plt.xlabel('x', fontsize=16)
plt.ylabel('y', fontsize=16)
plt.xlim(xstart,xend)
plt.ylim(ystart,yend)
plt.streamplot(X,Y,udoublet,vdoublet,density=2.0,linewidth=1,arrowsize=1,arrowstyle='->')
plt.scatter(xdoublet, ydoublet, c='r',s=80, marker ='o')
plt.show()
        
#Uniform flow past a doublet
uinf= 1.0 #freestream speed
alpha= 0.0 #angle of attack (indegrees)

ufreestream=uinf*cos(alpha*pi/180)*np.ones((N,N),dtype=float)
vfreestream=uinf*sin(alpha*pi/180)*np.ones((N,N), dtype=float)

psifreestream=uinf*(cos(alpha*pi/180)*Y-sin(alpha*pi/180)*X)

#superposition of the doublet on the freestream flow 
u=ufreestream+udoublet
v=vfreestream+vdoublet
psi=psifreestream+psidoublet

#plotting
size=10
plt.figure(figsize=(size,(yend-ystart)/(xend-xstart)*size))
plt.xlabel('x', fontsize=16)
plt.ylabel('y', fontsize=16)
plt.xlim(xstart,xend)
plt.ylim(ystart, yend)
plt.streamplot(X,Y,u,v, \
            density =2.0,linewidth=1,arrowsize=1,arrowstyle='->')
plt.contour(X,Y,psi,levels=[0.0],colors='r',linewidths=2,linestyles='solid')
plt.scatter(xdoublet,ydoublet,c='r', s=80, marker ='o')
plt.show()

#compute the pressure coefficient cp
Cp=1.0-(u**2+v**2)/uinf**2

#plotting
size=10
plt.figure(num=0,figsize=(size,(yend-ystart)/(xend-xstart)*size))
plt.xlabel('x', fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xstart, xend)
plt.ylim(ystart, yend)
contf=plt.contourf(X,Y,Cp,levels=np.linspace(-2.0,1.0,100),extend='both')
cbar=plt.colorbar(contf)
cbar.set_label(r'$C_p$',fontsize=16)
cbar.set_ticks([-2.0,-1.0,0.01,1.0])
plt.scatter(xdoublet,ydoublet, c='r',s=80,marker='o')
plt.contour(X,Y,psi,\
    levels=[0.0],colors='r',linewidths =2, linestyles ='solid')
    
plt.show()