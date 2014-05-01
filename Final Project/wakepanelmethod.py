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
xstart,xend = xmin-valX*(xmax-xmin), xmax+valX*4 #changed the x-end value 
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
    
    panel = np.empty(N+1,dtype=object) #changed the panel size from N to N+1 to include wake panel 
   
    for i in range(N+1):
        if (i==N):
            panel[i] = Panel(x[i],y[i],100000,y[i])#start of the wake panel is the endpoint of the last panel 
        else:                                      #endpoint is 100,000 chord lengths 
           panel[i] = Panel(x[i],y[i],x[i+1],y[i+1]) #including the body panels 
    return panel                                      
    
N = 20  #number of panels
panel = definePanels(N,xp,yp) #discretization of the geometry into panels 

#plotting the geometry with panels 
valX,valY = 0.1,0.2
xmin,xmax = min([p.xa for p in panel]), max([p.xa for p in panel])
ymin,ymax = min([p.ya for p in panel]),max([p.ya for p in panel])
xstart,xend = xmin-valX*(xmax-xmin),xmax+valX*4
ystart,yend = ymin-valY*(ymax-ymin),ymax+valY*(ymax-ymin)
size = 10
plt.figure(figsize=(size,(yend-ystart)/(xend-xstart)*size))
size = 10
plt.grid(True) 
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xstart,xend)
plt.ylim(ystart,yend)
#plt.plot(xp,yp,'k-',linewidth=2)
#plt.plot(np.append([p.xa for p in panel],panel[0].xa),\
       # np.append([p.ya for p in panel],panel[0].ya),\
        #linestyle='-',linewidth=1,marker='o', markersize=6,color='r'); 
plt.plot([p.xa for p in panel],[p.ya for p in panel],\
          linestyle='-',linewidth=1,marker='o',markersize=6,color='r');
plt.plot([p.xb for p in panel],[p.yb for p in panel],\
          linestyle='-',linewidth=1,marker='o',markersize=6,color='r');
plt.xlim(xstart,2) 
plt.show()


#---Class Freestream containing the freestream conditions ------
class Freestream:
    def __init__(self,uinf,alpha):
        self.uinf = uinf                #velocity magnitude 
        self.alpha = alpha*pi/180       #angle of attack (in degrees) 
        
#definition of the object freestream
uinf = 1.0          #freestream speed 
alpha = 0     #angle of attack (in degrees) 
freestream = Freestream(uinf,alpha)     #instantation of the object's freedom 

#---Function to evaluate the integral Iij(zi)---
def I(xci,yci,pj,dxdz,dydz):
    def func(s):
        return(((yci-(pj.yb+cos(pj.beta)*s))**2\
                -(xci-(pj.xb-sin(pj.beta)*s))**2)*dxdz\
                -(2*(xci-(pj.xb-sin(pj.beta)*s))*(yci-(pj.yb+cos(pj.beta)*s)))*dydz)\
                /(((xci-(pj.xb-sin(pj.beta)*s))**2\
                +(yci-(pj.yb+cos(pj.beta)*s))**2)**2)
    return integrate.quad(lambda s:func(s),0.,pj.length)[0]

#---Function to build doublet matrix---
def doubletMatrix(p):
    N = len(p)
    A = np.empty((N+1,N+1),dtype=float)
    np.fill_diagonal(A,0.5)
    for i in range (N): 
        for j in range (N): 
            if(i!=j):
                A[i,j] = 0.5/pi * I(p[i].xc, p[i].yc, p[j],+cos(p[i].beta),+sin(p[i].beta))
            if(i==N):
                A[i,:] = 0 
                A[i,0] = -1            #Apply kutta condition 
                A[i,N-1] = 1
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
A = doubletMatrix(panel) 
B = buildRHS(panel,freestream) 


#---Solve system of linear equations ----
var = np.linalg.solve(A,B)
for i in range (len(panel)):
    panel[i].K = var[i]

#sum of all source sink strengths
print '--> sum of doublet strengths:', sum([p.K*p.length for p in panel])

#calclation of the lift
Cl = p.K*sum([p.length for p in panel])/(0.5*freestream.uinf*(xmax-xmin))
print '--> Lift Coefficient: Cl=',Cl

#-----Generate stream lines-------------

def getVelocityField(panel,freestream,K,X,Y):
    Nx,Ny = X.shape
    u,v = np.empty((Nx,Ny),dtype=float),np.empty((Nx,Ny),dtype=float)
    for i in range(Nx):
        for j in range(Ny):
            u[i,j] = freestream.uinf*cos(freestream.alpha)\
                    +0.5/pi*sum([p.K*I(X[i,j],Y[i,j],p,1,0) for p in panel])
            
            v[i,j] = freestream.uinf*sin(freestream.alpha)\
                    +0.5/pi*sum([p.K*I(X[i,j],Y[i,j],p,0,1) for p in panel])
    return u,v

#definition of mesh grid
Nx,Ny = 20,20
valX,valY = 1.0,2.0

xmin,xmax = min([p.xa for p in panel]),max([p.xa for p in panel])
ymin,ymax = min([p.ya for p in panel]),max([p.ya for p in panel])

xstart,xend = xmin-valX*(xmax-xmin), 2
ystart,yend = ymin-valY*(ymax-ymin),ymax+valY*(ymax-ymin)

X,Y = np.meshgrid(np.linspace(xstart,xend,Nx),np.linspace(ystart,yend,Ny))

# get the velicity field on the mesh grid
u,v = getVelocityField(panel,freestream,gamma,X,Y)

# plotting the velocity field
size=10
plt.figure(figsize=(size,(yend-ystart)/(xend-xstart)*size))
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)

plt.streamplot(X,Y,u,v,density=2,linewidth=1,arrowsize=1,arrowstyle='->')
plt.fill([p.xa for p in panel],[p.ya for p in panel],'ko-',linewidth=2,zorder=2)
plt.xlim(xstart,xend)
plt.ylim(ystart,yend)
plt.title('Streamlines around a NACA 0012 airfoil, '+r'$\alpha=$'+str(alpha));
plt.show()
        
    

                
    
    
 

