import numpy as np
import matplotlib.pyplot as plt
from math import * 

#set up mesh grid 
N=50                           #number of points 
xstart, xend = (-2.0,2.0)     # x-dir boundaries
ystart, yend = (-1.0,1.0)     #y-dir boundaries 
x = linspace.np(xstart,xend,N)  #1D array in x-dir
y = linspace.np(ystart,yend,N)  #1D array in y-dir
X,Y = meshgrid.np(x,y)

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


