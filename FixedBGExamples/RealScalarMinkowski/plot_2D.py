import numpy as np
import matplotlib.pyplot as plt

# the function that I'm going to plot
t = 0.0
m = 1.0
def z_func(x,y):
 return (4.0*m*np.cos(m*t) * (x*np.sin(m*t)+y*np.cos(m*t)) )
 
x = np.arange(-40.0,40.0,1.0)
y = np.arange(-40.0,40.0,1.0)
X,Y = np.meshgrid(x, y) # grid of point
Z = z_func(X, Y) # evaluation of the function on the grid

contours = plt.contour(X, Y, Z, 3, colors='black')
plt.clabel(contours, inline=True, fontsize=8)

plt.imshow(Z, extent=[0, 5, 0, 5], origin='lower',
           cmap='RdGy', alpha=0.5)
plt.colorbar()
plt.savefig("S_x_theo_x+iy.png")
