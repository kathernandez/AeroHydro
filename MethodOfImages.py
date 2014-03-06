#lesson 7 Aerodymanics Interferences 
#and method of images

import numpy as np
import matplotlib.pyplot as plt
from math import *
from Ipython.core.display import clear_output

N = 50     #number of point in each direction 
xstart,xend = -2.0, 2.0  #x-direction boundaries
ystart,yend = -1.0, 1.0  #y-direction boundaries
x = np.linspace(xstart, xend, N)  #X 1D array 
y = np.linspace(ystart, yend, N)  #y 1D array
X,Y = np.meshgrid(x,y)            #generate mesh

class Source: 
    def _init_(self, strength, x, y):
        self.strength = strength
        self.x,self.y = x,y
    #get the velocity field 
    def velocity(self, X,Y):
        self.u = self.strength/(2*pi)*(X-self.x)/((X-self.x)**2+(Y-self.y)**2)
        self.v = self.strength/(2*pi)*(Y-self.y)/((X-self.x)**2+(Y-self.Y)**2)
    #get the stream function
    def streamfunction(self,X,Y):
        self.psi = self.strength/(2*pi)*np.arctan2((Y-self.y),(X-self.x))
        