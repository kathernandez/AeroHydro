#Final Project 
#First attempt to model flow using doublet panels with constant strength
#then try to add the wake behind the airfoil 

#import libraries 
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
from math import *

#read the geometry from a data file x
coords = np.loadtxt('C:/Users/Kat/Documents/Aero_Hydro/AeroPython-master/resources/naca0012.dat')
xp,yp = coords[:,0],coords[:,1]

#plotting the geometry 
valX,valY = 0.1,0.2
xmin,xmax = min(xp),max(xp)
ymin,ymax = min(yp),max(yp)
xstart,xend = xmin-valX*(xmax-xmin), xmax+valX*(xmax-xmin)
ystart,yend = ymin-valY*(ymax-ymin), ymax+valY*(ymax-ymin)
size = 10
plt.figure(figsize=(size,(yend-ystart)/(xend-xstart)*size))
plt.grid(True)
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xstart,xend)
plt.ylim(ystart,yend)
plt.plot(xp,yp,'k-',linewidth=2);
plt.show()

#-----Discretization into Panels-----------
class Panel: 
    def __init__(self,xa,ya,xb,yb):
        self.xa,self.ya = xa,ya                    #1st end-point 
        self.xb,self.yb = xb,yb                    #2nd end-point 
        self.xc,self.yc = (xa+xb)/2,(ya+yb)/2    #control point 
        self.length = sqrt((xb-xa)**2+(yb-ya)**2)  #length of the panel 
        
        #orientation of the panel
        #if (xb-xa<=0.): self.beta = acos((yb-ya)/self.length)
        #elif(xb-xa>0.): self.beta = pi + acos(-(yb-ya)/self.length)
        
        #new orientation of the panel 
        if (xb>xa): self.beta = acos((ya-yb)/self.length)
        else: self.beta = pi+acos((yb-ya)/self.length)
        
        #locations of the panel 
        if (self.beta<=pi): self.loc = 'extrados'
        else: self.loc = 'intrados'
        
        self.K = 0.       #doublet strength 
        self.vt = 0.     #tangential velocity 
        self.Cp = 0.    #pressure coefiicient 

#---Function to discretize the geometry into panels---
def definePanels(N,xp,yp):
    R = (max(xp)-min(xp))/2                                 #radius of circle 
    xcenter = (max(xp)+min(xp))/2                           #x-coord of center    
    xcircle = xcenter + R*np.cos(np.linspace(0,2*pi,N+1))   #x-coord of the circle points
    
    x = np.copy(xcircle) #projection of the xcoord on the surface 
    y = np.empty_like(x) #initialization of the y-coord Numpy Array 
    
    #xp,yp = np.append(xp,xp[0]),np.append(yp,yp[0])         #extend arrays using np.append 
    
    #CLOCKWISE
    xp,yp = np.append(xp[0],xp[::-1]),np.append(yp[0],yp[::-1])
    
    I = 0
    for i in range (N):
        while (I<len(xp)-1):
            if(xp[I]<=x[i]<=xp[I+1] or xp[I+1]<=x[i]<=xp[I]): break
            else: I += 1
        a = (yp[I+1]-yp[I])/(xp[I+1]-xp[I])
        b = yp[I+1]-a*xp[I+1]
        y[i] = a*x[i]+b           #interpolation being used to get y coordinates 
    y[N] = y[0]
    
    panel = np.empty(N,dtype=object)
    for i in range(N):
        panel[i] = Panel(x[i],y[i],x[i+1],y[i+1])
    return panel
    
N = 20  #number of panels
panel = definePanels(N,xp,yp) #discretization of the geometry into panels 

#would I include a panel here where the N = 1 and the start point would be the end
#point of of the airfoil and the end point would be the whate ever the x end is...
#and y would remain constant. Call it the wake panel. Use the Class Panel 

#plotting the geometry with panels 
valX,valY = 0.1,0.2
xmin,xmax = min([p.xa for p in panel]), max([p.xa for p in panel])
ymin,ymax = min([p.ya for p in panel]),max([p.ya for p in panel])
xstart,xend = xmin-valX*(xmax-xmin),xmax+valX*(xmax-xmin)
ystart,yend = ymin-valY*(ymax-ymin),ymax+valY*(ymax-ymin)
size = 10
plt.figure(figsize=(size,(yend-ystart)/(xend-xstart)*size))
size = 10
plt.grid(True) 
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xstart,xend)
plt.ylim(ystart,yend)
plt.plot(xp,yp,'k-',linewidth=2)
plt.plot(np.append([p.xa for p in panel],panel[0].xa),\
        np.append([p.ya for p in panel],panel[0].ya),\
        linestyle='-',linewidth=1,marker='o', markersize=6,color='r'); 
plt.show()

#---Class Freestream containing the freestream conditions ------
class Freestream:
    def __init__(self,uinf,alpha):
        self.uinf = uinf                #velocity magnitude 
        self.alpha = alpha*pi/180       #angle of attack (in degrees) 
        
#definition of the object freestream
uinf = 1.0          #freestream speed 
alpha = 0.0       #angle of attack (in degrees) 
freestream = Freestream(uinf,alpha)     #instantation of the object's freedom 

#---Function to evaluate the integral Iij(zi)---
def I(xci,yci,pj,dxdz,dydz):
    def func(s):
        return(+((yci-(pj.yb+cos(pj.beta)*s))**2\
                -(xci-(pj.xb-sin(pj.beta)*s))**2)*dxdz\
                -(2*(xci-(pj.xb-sin(pj.beta)*s))*(yci-(pj.yb+cos(pj.beta)*s)))*dydz)\
                /(((xci-(pj.xb-sin(pj.beta)*s))**2\
                +(yci-(pj.yb+cos(pj.beta)*s))**2)**2)
    return integrate.quad(lambda s:func(s),0.,pj.length)[0]

#---Function to build doublet matrix---
def doubletMatrix(p):
    N = len(p)
    A = np.empty((N+1,N+1),dtype=float)
    np.fill_diagonal(A,1)
    for i in range (N): 
        for j in range (N): 
            if(i!=j):
                A[i,j] = 0.5/pi * I(p[i].xc, p[i].yc, p[j],+cos(p[i].beta),+sin(p[i].beta))
            if(i==N):
                A[i,j] = 0 
                A[i,0] = 1
                A[i,N] = -1 
    return A 

#----Function to build RHS ------
def buildRHS(p,fs):
    N = len(p)
    B = np.zeros(N+1,dtype = float)
    for i in range(N):
        B[i] = -fs.uinf*cos(fs.alpha-p[i].beta)
        B[N] = 0
    return B

#----Build System---------
A = doubletMatrix(panel)   #here is the problem..panel does not take into account the wake panel 
B = buildRHS(panel,freestream) #need to modify the panel--> change the fuction which defines the geometry 


#---Solve system of linear equations ----
var = np.linalg.solve(A,B)
for i in range (len(panel)):
    panel[i].K = var[i]




        
    

                
    
    
 

