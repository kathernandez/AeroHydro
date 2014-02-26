import numpy as np
import matplotlib.pyplot as plt
from math import * 

#set up mesh grid 
N=50                           #number of points 
xstart, xend = (-2.0,2.0)     # x-dir boundaries
ystart, yend = (-1.0,1.0)     #y-dir boundaries 
x = np.linspace(xstart,xend,N)  #1D array in x-dir
y = np.linspace(ystart,yend,N)  #1D array in y-dir
X,Y = np.meshgrid(x,y)

k = 1   #strength of doublet 
xdoublet, ydoublet = 0.0, 0.0 #location of the doublet
uinf = 1                      #freestream velocity 

#define functions to compute u,v, and psi

def getVelocityDoublet(strength, xd, yd, X, Y):
    u = -strength/(2*pi) * ((X-xd)**2-(Y-yd)**2)/((X-xd)**2+(Y-yd)**2)**2
    v = -strength/(2*pi) * 2*(X-xd)*(Y-yd)/((X-xd)**2+(Y-yd)**2)**2
    return u,v

def getStreamFunctionDoublet (strength, xd, yd, X, Y):
    psi = -strength/(2*pi)* (Y-yd)/((X-xd)**2+(Y-yd)**2)
    return psi

#compute flow around a cylinder by adding free stream 

udoublet, vdoublet = getVelocityDoublet(k, xdoublet, ydoublet, X, Y)
psidoublet = getStreamFunctionDoublet(k, xdoublet, ydoublet, X, Y)

#compute free stream velocity components and psi

ufreestream = np.ones((N,N), dtype = float)
vfreestream = np.zeros((N,N), dtype = float)
psifreestream = uinf*Y

#apply superposition to create flow around a cylinder 

u = ufreestream + udoublet
v = vfreestream + vdoublet
psi = psidoublet + psifreestream


#Plotting 
size = 10
plt.figure(figsize=(size, (yend-ystart)/(xend-xstart)*size))
plt.xlim(xstart,xend)
plt.ylim(ystart,yend)
plt.xlabel('x')
plt.ylabel('y')
plt.streamplot(X,Y,u,v,\
        density = 2.0, linewidth = 1, arrowsize =1, arrowstyle = '->')
plt.scatter(xdoublet, ydoublet, c='m', s=80, marker = 'o')

#plot the radius of cylinder 
R=sqrt(k/(2*pi*uinf))
circle = plt.Circle((0,0),radius = R, color='m', alpha = 0.5)
plt.gca().add_patch(circle)

#plot stagnation points 
xstag1, ystag1 = + R, 0
xstag2, ystag2 = -R, 0
plt.scatter([xstag1,xstag2],[ystag1,ystag2], c ='g', s=80, marker ='o')
plt.show()


