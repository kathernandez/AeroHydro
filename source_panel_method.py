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

#-----------------Discretization into Panels--------------
