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
        self.xc, self.yc = (xa+xb)/2, (ya+yb)/2    #control point 
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

    



        

        