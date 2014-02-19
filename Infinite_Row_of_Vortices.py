
#Finite row of vortices

#First, we import python libraries 
import numpy as np
import matplotlib.pyplot as plt
from math import *

#Create a grid mesh 

N=50                         # number of grid points 
xstart,xend = -2.0,2.0       # x boundaries
ystart, yend = -1.0,1.0      # y boudaries 
x=np.linpace(xstart,xend,N)  # x direnction 1D array 
y=np.linspace(ystart,yend,N) # y direction 1D array 
X,Y=np.meshgrid(x,y)         # mesh grid created 


