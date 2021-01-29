import numpy as np
import math
from matplotlib import pyplot as plt

# plot a graph of actual darkness over expected surface depending on the radius

def r_theta(im, xc, yc):
    
    # returns the radius rr and the angle phi for point (xc,yc)
    ny, nx = im.shape
    yp, xp = np.mgrid[0:ny,0:nx]
    yp = yp - yc
    xp = xp - xc
    rr = np.sqrt(np.power(yp,2.) + np.power(xp,2.))
    phi = np.arctan2(yp, xp)
    
    return(rr, phi)

def planet_array(radius):
    
    # make sure matrix size is odd to center the planet in the array and have padding of "True"
    matrix_size = (radius * 2) + 3

    planet_array = np.ones((matrix_size, matrix_size))
    r_planet,t_planet = r_theta(planet_array, matrix_size/2, matrix_size/2)
    
    return r_planet

analytical = []
numerical = []
radius = []

for planet_radius in range (5,200):
    r_planet = planet_array(planet_radius) > planet_radius
    surface_of_circle = math.pi * (planet_radius ** 2)
    darkness = len(r_planet) ** 2 - np.sum(r_planet)
    numerical.append(darkness)
    analytical.append(surface_of_circle)
comparison = np.array(numerical) / np.array(analytical)
#plt.plot(numerical, label = 'numerical')
#plt.plot(analytical, label = 'analytical')
plt.plot(comparison, label = 'regular / ring')
plt.title("Ratio of numerical surface over analytical surface")
plt.xlabel("Radius (pixels)")
plt.ylabel("Num / Analytical circle surface")

plt.show()
