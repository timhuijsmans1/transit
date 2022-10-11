import numpy as np
import copy
import math
from matplotlib import pyplot as plt

#from kepler_test import coordinates as kepler_output
from kepler_test import kepler_cropped

def r_theta(im, xc, yc):
    
    # returns the radius rr and the angle phi for point (xc,yc)
    ny, nx = im.shape
    yp, xp = np.mgrid[0:ny,0:nx]
    yp = yp - yc
    xp = xp - xc
    rr = np.sqrt(np.power(yp,2.) + np.power(xp,2.))
    phi = np.arctan2(yp, xp)
    
    return(rr, phi)


def ellipse_r_theta(im, xc, yc, axis_ratio, theta):
   
    # axis ratio is b/a
    # theta is the orientation of the major axis input in deg

    ny, nx = im.shape
    yp, xp = np.mgrid[0:ny,0:nx]
     
    theta = math.radians(theta)    
   
    yp = yp - yc
    xp = xp - xc
   
    y_rot = xp * math.cos(theta) - yp* math.sin(theta)    
    x_rot = xp * math.sin(theta) + yp* math.cos(theta)    
   
    x_rot = x_rot * axis_ratio    
   
    rr = np.sqrt(np.power(y_rot,2.) + np.power(x_rot,2.))
    phi = np.arctan2(y_rot, x_rot)
   
    return(rr, phi)


def star_array_limb(radius, matrix_size, temp):

    star = np.zeros((matrix_size, matrix_size))
    cx = matrix_size / 2
    cy = matrix_size / 2
    cr = radius

    if temp == 4000:
        coefficient_array = [0.5651,0.1233,0.2108,-0.0817]
    if temp == 5000:
        coefficient_array = [0.6995,-0.7701,1.5147,-0.6341]

    a1 = coefficient_array[0]
    a2 = coefficient_array[1]
    a3 = coefficient_array[2]
    a4 = coefficient_array[3]

    for x in range(matrix_size):
        print("row number {}" .format(x))
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

def ellipse_array(a, theta, ratio):

    matrix_size = (a * 20) + 3

    planet_array = np.ones((matrix_size, matrix_size))
    r_ellipse, t_ellipse = ellipse_r_theta(planet_array, matrix_size/2, matrix_size/2, ratio, theta)

    return r_ellipse

def planet_array(radius):
    
    # make sure matrix size is odd to center the planet in the array and have padding of "True"
    matrix_size = (radius * 2) + 3

    planet_array = np.ones((matrix_size, matrix_size))
    r_planet,t_planet = r_theta(planet_array, matrix_size/2, matrix_size/2)
    
    return r_planet

def ring_array(radius):
    
    # make sure matrix size is odd to center the planet in the array and have padding of "True"
    matrix_size = (radius * 20) + 3

    planet_array = np.ones((matrix_size, matrix_size))
    r_planet,t_planet = r_theta(planet_array, matrix_size/2, matrix_size/2)
    
    return r_planet

def star_mask_by_planet(star_array, cx, cy, r_planet, star_brightness, planet_step, heatmap):
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

    if heatmap:
        star_array[ymin:ymax + 1, xmin:xmax + 1] = old_matrix * r_planet
        if planet_step == 290:
            plt.imshow(star_array, cmap = 'afmhot')
            plt.colorbar()
            plt.title('Relative brightness of the star with limb darkening')
            plt.show()

    brightness = star_brightness - (np.sum(old_matrix) - np.sum(old_matrix * r_planet))
    star_array[ymin:ymax + 1, xmin:xmax + 1] = old_matrix

    return brightness

# This function loops over the received kepler output
def planet_progression(translated_kepler, star, r_planet, star_brightness, heatmap):    
    count = 0
    brightness_profile = []
    for coord in translated_kepler:
        brightness = star_mask_by_planet(star, int(coord[0]), int(coord[1]), r_planet, star_brightness, count, heatmap)
        brightness_profile.append(brightness)
        count+=1
        print("planet step: {}" .format(count))
    return(brightness_profile)

def retrieve_coordinates(a, i, w_angle, omega, P, T, e, timestep, crop_range, star_radius, star_array_size, kepler_star_radius):

    # set kepler coordinates
    kepler_coordinates = np.array(kepler_cropped(a, i, w_angle, omega, P, T, e, timestep, crop_range)[2])
    scaled_kepler = kepler_coordinates * (star_radius / (kepler_star_radius))
    translated_kepler = [x+(star_array_size / 2) for x in scaled_kepler]
    
    return translated_kepler