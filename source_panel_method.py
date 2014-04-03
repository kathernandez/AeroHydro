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
        if (xb-xa<=0.): self.beta = acos((yb-ya)/self.length)
        elif(xb-xa>0.): self.beta = pi + acos(-(yb-ya)/self.length)
        
        #locations of the panel 
        if (self.beta<=pi): self.loc = 'extrados'
        else: self.loc = 'intados'
        
        self.sigma = 0.  #source strength 
        self.vt = 0.     #tangential vlocity 
        self.Cp = 0.    #pressure coefiicient 

#---Function to discretize the geometry into panels---
def definePanels(N,xp,yp):
    R = (max(xp)-min(xp))/2
    xc,yc = (max(xp)+min(xp))/2,(max(yp)+min(yp))/2
    xcircle = xc + R*np.cos(np.linspace(0,2*pi,N+1))
    ycircle = yc + R*np.sin(np.linspace(0,2*pi,N+1))
    
    x = np.copy(xcircle[0:-1])
    y = np.empty_like(x)
    
    I = 0 
    for i in range (N):
        while(I<len(xp)-1):
            if(xp[I]<=x[i]<=xp[I+1] or xp[I+1]<=x[i]<=xp[I]): break
            else: I+=1
        a = (yp[(I+1)%len(yp)]-yp[I])/(xp[(I+1)%len(yp)]-xp[I])  #percent thing is to obtain the remainder 
        b = yp[(I+1)%len(yp)]-a*xp[(I+1)%len(xp)]
        y[i] = a*x[i]+b
    panel = np.empty(N,dtype=object)
    for i in range (N):
        panel[i] = Panel(x[i],y[i],x[(i+1)%N],y[(i+1)%N])
    
    return panel
    
N = 20  #number of panels
panel = definePanels(N,xp,yp) #discretization of the geometry into panels 

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

#--Class Freestream containg the freestream conditions ---
class Freestream: 
    def __init__(self,uinf,alpha):
        self.uinf = uinf                #velocity magnitude 
        self.alpha = alpha*pi/180       #angle of attack (in degrees) 
        
#definition of the object freestream
uinf = 1.0          #freestream speed 
alpha = 5.0         #angle of attack (in degrees) 
freestream = Freestream(uinf,alpha)     #instantation of the object's freedom 

#function to evaluate the integral Iij(zi)
def I(xci,yci,pj,dxdz,dydz):
    def func(s):
        return (+(xci-(pj.xa-sin(pj.beta)*s))*dxdz\
        +(yci-(pj.ya+cos(pj.beta)*s))*dydz)\
        /((xci-(pj.xa-sin(pj.beta)*s))**2\
        +(yci-(pj.ya+cos(pj.beta)*s))**2)
    return integrate.quad(lambda s:func(s),0.,pj.length)[0]

#------Building the linear system------------------------

#function to build the source matrix
def buildMatrix(p):
    N = len(p)
    A = np.empty((N,N),dtype=float)
    np.fill_diagonal(A,0.5)
    for i in range(N):
        for j in range(N):
            if(i!=j):
                A[i,j] = 0.5/pi*I(p[i].xc,p[i].yc,p[j],cos(p[i].beta),sin(p[i].beta))
    return A

#function to build the right hand side of the linear system 
def buildRHS(p,fs):
    N = len(p)
    B = np.zeros(N,dtype=float)
    for i in range(N):
        B[i] = -fs.uinf*cos(fs.alpha-p[i].beta)
    return B 
    
A = buildMatrix(panel)          #calculate the singularity matrix
B = buildRHS(panel,freestream)  #calculate the freestream RHS

#solve for the linear system 
var = np.linalg.solve(A,B)
for i in range(len(panel)):
    panel[i].sigma = var[i]

#function to calculate the tangential velocity at each control point 
def getTangentVelocity(p,fs,gamma):
    N = len(p)
    A = np.zeros((N,N),dtype = float)
    for i in range(N):
        for j in range(N):
            if(i!=j):
                A[i,j] = 0.5/pi*I(p[i].xc,p[i].yc,p[j],-sin(p[i].beta),cos(p[i].beta))
    B = fs.uinf*np.sin([fs.alpha-pp.beta for pp in p])
    var = np.array([pp.sigma for pp in p])
    vt = np.dot(A,var)+B
    for i in range(N):
        p[i].vt = vt[i]

getTangentVelocity(panel,freestream,gamma)      #get tangential velocity 

#function to get Cp at each control point
def getPressureCoef(p,fs):
    for i in range(len(p)):
        p[i].Cp = 1-(p[i].vt/fs.uinf)**2

getPressureCoef(panel,freestream)

#plotting the coefficient of pressure 
valX,valY = 0.1,0.2
xmin,xmax = min([p.xa for p in panel]),max([p.xa for p in panel])
Cpmin,Cpmax = min([p.Cp for p in panel]), max([p.Cp for p in panel])
xstart,xend = xmin - valX*(xmax-xmin), xmax+valX*(xmax-xmin)
ystart,yend = Cpmin-valY*(Cpmax-Cpmin),Cpmax+valY*(Cpmax-Cpmin)
plt.figure(figsize=(10,6))
plt.grid(True)
plt.xlabel('x',fontsize=16)
plt.ylabel('$C_p$, fontsize=16')
plt.plot([p.xc for p in panel if p.loc=='extrados'],\
        [p.Cp for p in panel if p.loc=='extrados'],'ro-',linewidth=2)  
plt.plot([p.xc for p in panel if p.loc=='intrados'],\
          [p.Cp for p in panel if p.loc=='intrados'],'bo-',linewidth=1)
plt.legend(['extrados','intrados'],'best',prop={'size':14})                
plt.xlim(xstart,xend)
plt.ylim(ystart,yend)
plt.gca().invert_yaxis()
plt.title('Number of panels: %d'%len(panel));
plt.show()

#function to calculate the velocity field of a mesh grid
def getVelocityField(panel,freestream,gamma,X,Y):
    Nx,Ny = X.shape
    u,v = np.empty((Nx,Ny),dtype=float),np.empty((Nx,Ny),dtype=float)
    for i in range(Nx):
        for j in range(Ny):
            u[i,j] = freestream.uinf*cos(freestream.alpha)\
              +0.5/pi*sum([p.sigma*I(X[i,j],Y[i,j],p,1,0)for p in panel])
            v[i,j] = freestream.uinf*sin(freestream.alpha)\
              +0.5/pi*sum([p.sigma*I(X[i,j],Y[i,j],p,0,1)for p in panel])
    return u,v

#definition of the mesh grid 
Nx,Ny = 20,20
valX,valY = 1.0,2.0
xmin,xmax = min([p.xa for p in panel]),max([p.xa for p in panel])
ymin,ymax = min([p.ya for p in panel]),max([p.ya for p in panel])
xstart,xend = xmin-valX*(xmax-xmin),xmax+valX*(xmax-xmin)
ystart,yend = ymin-valY*(ymax-ymin),ymax+valY*(ymax-ymin)
X,Y = np.meshgrid(np.linspace(xstart,xend,Nx),np.linspace(ystart,yend,Ny))

#get the velocity field on the mesh grid 
u,v = getVelocityField(panel,freestream,gamma,X,Y)  

#plotting the velocity field
size=10
plt.figure(figsize=(size,(yend-ystart)/(xend-xstart)*size))
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.streamplot(X,Y,u,v,density=1,linewidth=1,arrowsize=1,arrowstyle='->')
plt.fill([p.xa for p in panel],[p.ya for p in panel],'ko-',linewidth=2,zorder=1)
plt.xlim(xstart,xend)
plt.ylim(ystart,yend)
plt.title('Contour of velocity field');
plt.show()  

#computing the pressure field 
Cp = 1.0-(u**2+v**2)/freestream.uinf**2

#plotting the pressure field 
size = 12
plt.figure(figsize=(1.1*size,(yend-ystart)/(xend-xstart)*size))
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
contf=plt.contourf(X,Y,Cp,levels=np.linspace(-2.0,1.0,100),extend='both')
cbar=plt.colorbar(contf)
cbar.set_label('$C_p$',fontsize=16)
cbar.set-ticks([-2.0,-1.0,0.0,1.0])
plt.fill([p.xc for p in panel],[p.yc for p in panel],'ko-',linewidth=2,zorder=2)
plt.xlim(xstart,xend)
plt.ylim(ystart,yend)
plt.title('Contour of pressure field')
plt.show()


        

        