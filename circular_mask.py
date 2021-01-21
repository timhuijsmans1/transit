import numpy as np
import copy
import math
from matplotlib import pyplot as plt

from kepler_test import coordinates as kepler_output


def r_theta(im, xc, yc):
    
    # returns the radius rr and the angle phi for point (xc,yc)
    ny, nx = im.shape
    yp, xp = np.mgrid[0:ny,0:nx]
    yp = yp - yc
    xp = xp - xc
    rr = np.sqrt(np.power(yp,2.) + np.power(xp,2.))
    phi = np.arctan2(yp, xp)
    
    return(rr, phi)


def star_array_limb(radius, matrix_size):

    star = np.zeros((matrix_size, matrix_size))
    cx = matrix_size / 2
    cy = matrix_size / 2
    cr = radius

    coefficient_array = [0.5651,0.1233,0.2108,-0.0817]
    #coefficient_array = [0.6995,-0.7701,1.5147,-0.6341]
    #coefficient_array = [0.5409,-0.0366,0.5688,-0.3213]
    #coefficient_array = [0.3682,0.7769,-0.6243,0.1868]

    a1 = coefficient_array[0]
    a2 = coefficient_array[1]
    a3 = coefficient_array[2]
    a4 = coefficient_array[3]

    for x in range(matrix_size):
        #print("row number {}" .format(x))
        for y in range(matrix_size): 
            
            r = math.sqrt((x - cx)**2 + (y - cy)**2)
            
            if r < cr:
                mu = math.cos(math.asin(r / cr)) # constant used here is distance to CoRoT-1 in matrix elements
                I = 1 - a1*(1-mu**0.5) - a2*(1-mu) - a3*(1-mu**1.5) - a4*(1-mu**2)
                star[y][x] = I # set the value of this point according to the value of the intensity determined by the non linear formula

    return star

def star_array_no_limb(radius, matrix_size):

    star = np.zeros((matrix_size, matrix_size))
    cx = matrix_size / 2
    cy = matrix_size / 2
    cr = radius

    for x in range(matrix_size):
        #print("row number {}" .format(x))
        for y in range(matrix_size): 
            
            r = math.sqrt((x - cx)**2 + (y - cy)**2)
            
            if r < cr:
                star[y][x] = 1 # set the value of this point according to the value of the intensity determined by the non linear formula
    return star

def planet_array(radius):
    
    # make sure matrix size is odd to center the planet in the array and have padding of "True"
    matrix_size = (radius * 2) + 3

    planet_array = np.ones((matrix_size, matrix_size))
    r_planet,t_planet = r_theta(planet_array, matrix_size/2, matrix_size/2)
    
    return r_planet

def star_mask_by_planet(star_array, cx, cy, r_planet, star_brightness):
    # flip the y-coordinate to correct for the upside down orientation
    cy_flip = len(star_array) - 1 - cy
    matrix_size = len(r_planet)

    # set the part of the matrix to be replaced
    ymin = cy_flip - matrix_size / 2
    ymax = cy_flip + matrix_size / 2
    xmin = cx - matrix_size / 2
    xmax = cx + matrix_size / 2

    # check brightness difference for small part of the matrix
    old_matrix = np.array(star_array[ymin:ymax + 1, xmin:xmax + 1])
    
    # star_array[ymin:ymax + 1, xmin:xmax + 1] = old_matrix * r_planet
    # brightness = np.sum(star_array)
    
    brightness = star_brightness - (np.sum(old_matrix) - np.sum(old_matrix * r_planet))
    star_array[ymin:ymax + 1, xmin:xmax + 1] = old_matrix

    return brightness

if __name__ == "__main__":
    # parameters of star and planet (probably don't want to go much beyond 300 planet and 4000 star)
    planet_radius = 2
    star_radius = 3000
    
    # translated input into matrix sizes
    star_array_size = int(2.8 * star_radius)
    y = star_array_size / 2

    # initialize output
    brightness_profile = []

    # create the star (remove no in case you wish to use limb darkening)
    star = star_array_limb(star_radius, star_array_size)
    star_brightness = np.sum(star)

    # create matrix of the planet
    r_planet = planet_array(planet_radius) > planet_radius

    # because we import the kepler function's variable coordinate, we can call it here the format is a list of coordinate lists: [[x1,y1], [x2,y2], ....]
    kepler_coordinates = np.array(kepler_output)
    scaled_kepler = kepler_coordinates * (star_radius / (0.03))
    translated_kepler = [x+(star_array_size / 2) for x in scaled_kepler]

    # This for loop is now used to move the planet through the matrix. If we have useful Kepler output, we can instead loop over the xy-coordinates in the output
    count = 0
    for coord in translated_kepler:
        brightness = star_mask_by_planet(star, int(coord[0]), int(coord[1]), r_planet, star_brightness)
        brightness_profile.append(brightness)
        count+=1
        print("planet step: {}" .format(count))

    plt.plot(brightness_profile)
    plt.title("Lightcurve of a transit with limb-darkening for r-star: %s pixels and r-planet: %s pixels (w = 180 degrees)" % (star_radius, planet_radius))
    plt.xlabel("Timestep")
    plt.ylabel("Brightness")
    plt.show()



