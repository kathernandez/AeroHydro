#Lesson 11 Source-Vortex Panel Methos 
import numpy as np
from scipy import integrate
from math import *
import matplotlib.pyplot as plt 

#-----------Discretization into panels------------

#read of the geometry from the a data file 
coords = np.loadtxt('C:/Users/Kat/Documents/Aero_Hydro/AeroPython-master/resources/naca0012.dat')
xp,yp = coords[:,0], coords[:,1]

# plotting the geometry 
valX,valY = 0.1,0.2
xmin,xmax = min(xp),max(xp)
ymin,ymax = min(yp),max(yp)
xstart,xend = xmin-valX*(xmax-xmin),xmax+valX*(xmax-xmin)
ystart,yend = ymin-valY*(ymax-ymin),ymax+valY*(ymax-ymin)
size = 10
plt.figure(figsize=(size,(yend-ystart)/(xend-xstart)*size))
plt.grid(True)
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xstart,xend)
plt.ylim(ystart,yend)
plt.plot(xp,yp,'k-',linewidth=2);
plt.show()

class Panel:
    def __init__(self,xa,ya,xb,yb):
        self.xa,self.ya = xa,ya  #1st end point
        self.xb,self.xb = xb,yb #second end point 
        self.xc,self.yc = (xa+xb)/2, (ya+yb)/2  #control point 
        self.length = sqrt((xb-xa)**2+(yb-ya)**2)  #length of the panel 
        
        #orientation of the panel 
        if(xb-xa<=0.):self.beta = acos((yb-ya)/self.length)
        elif(xb-xa>0): self.beta = pi+acos(-(yb-ya)/self.length)
        
        #location of panel
        if(self.beta<=pi):self.loc='extrados'
        else:self.loc = 'intrados'
        
        self.sigma = 0.      #source strength
        self.vt = 0.         #tangential velocity 
        self.Cp = 0.         #pressure coefficient 

#function to discretize the geometry into panels 
def definePanels(N,xp,yp):
    R = (max(xp)-min(xp))/2                                 #radius of circle 
    xcenter = (max(xp)+min(xp))/2                           #x-coord of center    
    xcircle = xcenter + R*np.cos(np.linspace(0,2*pi,N+1))   #x-coord of the circle points
    
    x = np.copy(xcircle) #projection of the xcoord on the surface 
    y = np.empty_like(x) #initialization of the y-coord Numpy Array 
    
    xp,yp = np.append(xp,xp[0]),np.append(yp,yp[0])         #extend arrays using np.append 
    
    I = 0
    for i in range (N):
        while (I<len(xp)-1):
            if(xp[I]<=x[i]<=xp[I+1] or xp[I+1]<=x[i]<=xp[I]): break
            else: I += 1
        a = (yp[I+1]-yp[I])/(xp[I+1]-xp[I])
        b = yp[I+1]-a*xp[I+1]
        y[i] = a*x[i]+b
    y[N] = y[0]
    
    panel = np.empty(N,dtype=object)
    for i in range(N):
        panel[i] = Panel(x[i],y[i],x[i+1],y[i+1])
    return panel
    
#----------Define the geometry for the air foil panels 
N = 20                              #number of panels 
panel = definePanels(N,xp,yp)       #discretizaton of the geometry into panels 

#plotting the geometry with the panels 
valX,valY = 0.1, 0.2 
xmin,xmax = min([p.xa for p in panel]), max([p.xa for p in panel]) 
ymin, ymax = min([p.ya for p in panel]), max([p.xa for p in panel]) 
xstart,xend = xmin-valX*(xmax-xmin),xmax+valX*(xmax-xmin)
ystart,yend = ymin-valY*(ymax-ymin),ymax+valY*(ymax-ymin)
size = 10
plt.figure(figsize=(size,(yend-ystart)/(xend-xstart)*size))
plt.grid(True) 
plt.xlabel('x', fontsize = 16)
plt.ylabel('y',fontsize =16)
plt.xlim(xstart,xend)
plt.ylim(ystart,yend)
plt.plot(xp,yp,'k-',linewidth=1)
plt.plot(np.append([p.xa for p in panel],panel[0].xa),\
        np.append([p.ya for p in panel], panel[0].ya),\
        linestyle = '-', linewidth=1,marker='o', markersize = 6, color ='r')
plt.show()

#-------add free stream conditions-----------------------------------
class Freestream:
    def __init__(self,uinf,alpha):
        self.uinf=uinf             #velocity magnitude 
        self.alpha=alpha*pi/180    #angle of attack (in degrees--> radians)

#definition of the object freestream 
uinf = 1    #free stream speed
alpha = 1   #angle of attack (in degrees) 
freestream = Freestream(uinf,alpha)   #instantiation of the object freestream 

#function to evaluate the integral Tij(zi)
def I(xci,yci,pj,dxdz,dydz):
    def func(s):
        return(+(xci-(pj.xa-sin(pj.beta)*s))*dxdz\
        +(yci-(pj.ya+cos(pj.beta)*s))*dydz)\
        /((xci-(pj.xa-sin(pj.beta)*s))**2\
        +(yci-(pj.ya+cos(pj.beta)*s))**2)
    return integrate.quad(lambda s:func(s),0.,pj.length)[0]
    

#function to build the source matrix
def sourceMatrix(p):
    N = len(p)
    A = np.empty((N,N),dtype=float)
    np.fill_diagonal(A,0.5)
    for i in range(N):
        for j in range(N):
            if(i!=j):
                A[i,j] = 0.5/pi*I(p[i].xc, p[i].yc, p[j],+cos(p[i].beta),+sin(p[i].beta))
    return A 



#function to build the vortex array 
def vortexArray(p): 
    N=len(p)
    B=np.zeros(N,dtype=float)
    for i in range(N):
        for j in range(N):
            if(j!=i):
                B[i] -= 0.5/pi*I(p[i].xc,p[i].yc,p[j],+sin(p[i].beta),-cos(p[i].beta))
    return B 


#function to build the kutta array 
def kuttaArray(p):
    N = len(p)
    B = np.zeros(N+1,dtype=float)
    for j in range (N): 
        if (j==0):
            B[j] = 0.5/pi*I(p[N-1].xc, p[N-1].yc, p[j], -sin(p[N-1].beta), +cos(p[N-1].beta))
        elif(j==N-1):
            B[j] = 0.5/pi*I(p[0].xc,p[0].yc,p[j],-sin(p[0].beta),+cos(p[0].beta))
        else: 
            B[j] = 0.5/pi*I(p[0].xc,p[0].yc,p[j],-sin(p[0].beta),+cos(p[0].beta))\
            + 0.5/pi*I(p[N-1].xc, p[N-1].yc,p[j],-sin(p[N-1].beta),+cos(p[N-1].beta))
            
            B[N] -= 0.5/pi*I(p[0].xc,p[0].yc,p[j],+cos(p[0].beta),+sin(p[0].beta))\
                + 0.5/pi*I(p[N-1].xc,p[N-1].yc,p[j],+cos(p[N-1].beta),+sin(p[N-1].beta))
    return B



#function to assemble the global matrix 
def buildMatrix(panel):
    N = len(panel)
    A = np.empty((N+1,N+1),dtype = float)
    AS = sourceMatrix(panel)
    BV = vortexArray(panel)
    BK = kuttaArray(panel)
    A[0:N,0:N],A[0:N,N],A[N,:] = AS[:,:],BV[:],BK[:]
    return A 
    

#funtion to build the right hand side of the linear system 
def buildRHS(p,fs):
    N = len(p)
    B = np.zeros(N+1,dtype = float)
    for i in range(N):
        B[i] = -fs.uinf*cos(fs.alpha-p[i].beta)
    B[N] = -fs.uinf*(sin(fs.alpha-p[0].beta)+sin(fs.alpha-p[N-1].beta))
    return B 


A = buildMatrix(panel)  #calculate the singularity matrix 
B = buildRHS(panel,freestream)  #calculate the freestream RHS


#solve the linear system 
var = np.linalg.solve(A,B)
for i in range (len(panel)):
    panel[i].sigma = var[i]
gamma=var[-1]


#function to calculate the tangential velocity at each control point 


    

    
            
            
    




        