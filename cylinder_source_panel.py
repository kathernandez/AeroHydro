#lesson 9 flow over a cylinder with source panels 
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
from math import *


uinf = 1

#definition of the cylinder

R = 1.0
theta = np.linspace(0,2*pi,100)
xcylinder,ycylinder = R*np.cos(theta),R*np.sin(theta)

#plot

size = 4
plt.figure(num=0,figsize=(size,size))
plt.grid(True)
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(-1.1,1.1)
plt.ylim(-1.1,1.1)
plt.plot(xcylinder,ycylinder,c='b',ls='-',lw=2)

class Panel:
    def __init__(self,xa,ya,xb,yb):
        self.xa,self.ya = xa,ya   #1st end point
        self.xb,self.yb = xb,yb   #2nd end point
        self.xc,self.yc = (xa+xb)/2,(ya+yb)/2   #center point
        self.length = np.sqrt((xb-xa)**2+(yb-ya)**2)  #length of panel
        
        #orientation of the panel
        if(xb-xa<=0.):self.beta = acos((yb-ya)/self.length)
        elif(xb-xa>0.):self.beta = pi + acos(-(yb-ya)/self.length)
        
        self.sigma = 0.    #source strength
        self.vt = 0.       #tangential velocity
        self.Cp = 0.       #pressure coefficient


Np = 10  #number of panels desired

#defining the end points of the panel
xb = R*np.cos(np.linspace(0,2*pi,Np+1))
yb = R*np.sin(np.linspace(0,2*pi,Np+1))

#defining the panels 
panel = np.empty(Np,dtype=object)
for i in range(Np):
    panel[i] = Panel(xb[i],yb[i],xb[i+1],yb[i+1])

#plotting the discretization 
size = 6
plt.figure(num=1,figsize=(size,size))
plt.grid(True)
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.plot(xcylinder,ycylinder,c='b',ls='-',lw=1)
plt.plot(xb,yb,c='r',ls='-',lw=2)
plt.scatter([p.xa for p in panel],[p.ya for p in panel],c='r',s=40)
plt.scatter([p.xc for p in panel],[p.yc for p in panel],c='k',s=40,zorder=3)
plt.legend(['cylinder','panels','end points','center points'],\
            loc='best',prop={'size':14})
plt.xlim(-1.1,1.1)
plt.ylim(-1.1,1.1);
plt.show()

#function to evaluate the integral 
def I(pi,pj):
    def func(s):
        return(+(pi.xc-(pj.xa-sin(pj.beta)*s))*cos(pi.beta)\
               +(pi.yc-(pj.ya+cos(pj.beta)*s))*sin(pi.beta))\
               /((pi.xc-(pj.xa-sin(pj.beta)*s))**2\
               +(pi.yc-(pj.ya+cos(pj.beta)*s))**2)
    return integrate.quad(lambda s:func(s),0.,pj.length)[0]

A = np.empty((Np,Np),dtype=float)
for i in range(Np):
    for j in range(Np):
        if (i!=j):
            A[i,j] = 0.5/pi*I(panel[i],panel[j])
        else:
            A[i,j] = 0.5

b = -uinf*np.cos([p.beta for p in panel])


#solve linear system
var = np.linalg.solve(A,b)
for i in range(len(panel)):
    panel[i].sigma = var[i]
    
#function to evaluate the integral Iij(zi)
def J(pi,pj):
    def func(s):
        return (-(pi.xc-(pj.xa-sin(pj.beta)*s))*sin(pi.beta)\
              +(pi.yc-(pj.ya+cos(pj.beta)*s))*cos(pi.beta))\
              /((pi.xc-(pj.xa-sin(pj.beta)*s))**2\
              +(pi.yc-(pj.ya+cos(pi.beta)*s))**2)
    return integrate.quad(lambda s:func(s),0.,pj.length)[0]
    
A=np.zeros((Np,Np),dtype=float)
for i in range (Np):
    for j in range (Np):
        if(i!=j):
            A[i,j] = 0.5/pi*J(panel[i],panel[j])

B = -uinf*np.sin([p.beta for p in panel])
sigma = np.array([p.sigma for p in panel])
Vt = np.dot(A,sigma) +B
for i in range(Np):
    panel[i].Vt = Vt[i]

#Calculate pressure coefficient 
for i in range(Np):
    panel[i].Cp = 1-(panel[i].Vt/uinf)**2
    

#plotting Cp
plt.figure(figsize=(10,6))
plt.grid(True)
plt.xlabel('x',fontsize=16)
plt.ylabel('$C_p$',fontsize=16)
plt.plot(xcylinder, 1-4*(ycylinder/R)**2,c='b',ls='-',lw=1,zorder=1)
plt.scatter([p.xc for p in panel],[p.Cp for p in panel],c='r',s=40,zorder=2)
plt.title('Number of panels: %d'%len(panel),fontsize=16)
plt.legend(['analytical','source panel method'],loc = 'best',prop={'size':16})
plt.xlim(-1.0,1.0)
plt.ylim(-4.0,2.0)
plt.show()



#----------Challenge task----------

#generate mesh grid 
N = 100                          #number of points in each direction 
xstart,xend = -1.0,1.0           #x-dir boundaries
ystart, yend = -1.5,1.5          #y-dir boundaries 
x = np.linspace(xstart,xend,N)   #x 1D array 
y = np.linspace(ystart,yend,N)  #y 1D array 
X,Y = np.meshgrid(x,y)           #generate meshgrid

uinf = 1.0         #free stream velocity
ufreestream = uinf*np.ones((N,N),dtype=float) 
vfreestream = np.zeros((N,N),dtype=float)




            
    
        

            




