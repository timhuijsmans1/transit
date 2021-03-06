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

if __name__ == "__main__":
    # parameters of star and planet
    star_radius = 2000

    """import kepler crop function and set initial values"""
    a = 0.3 #AU
    i = math.radians(90)
    w_angle = math.radians(180) 
    omega = math.radians(0)
    P = 0.5 # period of half a year
    T = -0.01 #time passage through periastron (seconden?) #where the planet starts
    e = 0
    timestep = (0.000003802651) # minutes when P is set to 0.5 yrs

    r_star = 0.01 * a #(Approx 7 solar radii)
    crop_range = r_star * 1.4
    
    # translate input into matrix sizes
    star_array_size = int(3.5 * star_radius)
    y = star_array_size / 2

    # create the star and set the initial brightness (remove no in case you wish to use limb darkening)
    star = star_array_limb(star_radius, star_array_size)
    star_brightness = np.sum(star)

def inclination_plot(normalizer):
    """
    For the inclination plot, we take a larger star radius (Approx. 70 solar radii), since the planet doesn't transit anymore for a smaller star
    """
    
    planet_radius = 20
    r_planet = planet_array(planet_radius) > planet_radius

    # retrieve and translate the Kepler coordinates into x,y coordinates of the matrix
    translated_kepler_1 = retrieve_coordinates(a, math.radians(90), w_angle, omega, P, T, e, timestep, crop_range, star_radius, star_array_size, r_star)
    translated_kepler_2 = retrieve_coordinates(a, math.radians(89.8), w_angle, omega, P, T, e, timestep, crop_range, star_radius, star_array_size, r_star)
    translated_kepler_3 = retrieve_coordinates(a, math.radians(89.6), w_angle, omega, P, T, e, timestep, crop_range, star_radius, star_array_size, r_star)
    translated_kepler_4 = retrieve_coordinates(a, math.radians(89.5), w_angle, omega, P, T, e, timestep, crop_range, star_radius, star_array_size, r_star)

    if normalizer == True:
        max_brightness = star_brightness
    else:
        max_brightness = 1

    # these lines are for creating inclination comparison lightcurves
    brightness_profile_1 = planet_progression(translated_kepler_1, star, r_planet, star_brightness, False) / max_brightness
    brightness_profile_2 = planet_progression(translated_kepler_2, star, r_planet, star_brightness, False) / max_brightness
    brightness_profile_3 = planet_progression(translated_kepler_3, star, r_planet, star_brightness, False) / max_brightness
    brightness_profile_4 = planet_progression(translated_kepler_4, star, r_planet, star_brightness, False) / max_brightness


    plt.plot(brightness_profile_1, label = '90')
    plt.plot(brightness_profile_2, label = '89.8')
    plt.plot(brightness_profile_3, label = '89.6')
    plt.plot(brightness_profile_4, label = '89.5')
    plt.title("Lightcurve of a transit with limb-darkening (from periastron)")
    plt.legend(title= 'Inclination angle (deg)', loc= 'upper center')
    plt.xlim(0, len(translated_kepler_1))
    plt.xlabel("Time (minutes)")
    plt.ylabel("Total image brightness")
    plt.show()

def radii_plot(normalizer):
    planet_radius = [10, 20, 40]

    # create matrix of the planet, 
    planet_matrices = []
    for radius in planet_radius:    
        r_planet = planet_array(radius) > radius
        planet_matrices.append(r_planet)
    
    # retrieve and translate the Kepler coordinates into x,y coordinates of the matrix
    translated_kepler = retrieve_coordinates(a, i, w_angle, omega, P, T, e, timestep, crop_range, star_radius, star_array_size, r_star)
    
    if normalizer == True:
        max_brightness = star_brightness
    else:
        max_brightness = 1

    # these lines are for creating radius comparison lightcurves
    brightness_profile_1 = planet_progression(translated_kepler, star, planet_matrices[0], star_brightness, False) / max_brightness
    brightness_profile_2 = planet_progression(translated_kepler, star, planet_matrices[1], star_brightness, False) / max_brightness
    brightness_profile_3 = planet_progression(translated_kepler, star, planet_matrices[2], star_brightness, False) / max_brightness

    plt.plot(brightness_profile_1, label = '0.35')
    plt.plot(brightness_profile_2, label = '0.7')
    plt.plot(brightness_profile_3, label = '1.4')
    plt.title("Lightcurve of a transit with limb-darkening (from periastron)")
    plt.legend(title= 'radius of planet in earth radii', loc= 'center')
    plt.xlim(0, len(translated_kepler))
    plt.ylim(0.9993, 1)
    plt.xlabel("Time (minutes)")
    plt.ylabel("Relative image brightness")
    plt.show()

def circular_ring_plot(normalizer):
    
    planet_radius = 20

    r_planet_regular = (planet_array(planet_radius) > planet_radius)
    r_planet_ring = (1 - ((ring_array(planet_radius) > planet_radius + 30) * (ring_array(planet_radius) < planet_radius + 40))) * (ring_array(planet_radius) > planet_radius)
    
    translated_kepler = retrieve_coordinates(a, math.radians(90), w_angle, omega, P, T, e, timestep, crop_range, star_radius, star_array_size, r_star)
    
    if normalizer == True:
        max_brightness = star_brightness
    else:
        max_brightness = 1

    brightness_profile_regular = planet_progression(translated_kepler, star, r_planet_regular, star_brightness, False) / max_brightness
    brightness_profile_ring = planet_progression(translated_kepler, star, r_planet_ring, star_brightness, True) / max_brightness

    plt.plot(brightness_profile_regular, label = 'regular')
    plt.plot(brightness_profile_ring, label = 'ring')
    plt.title("Lightcurve of a transit with limb-darkening (from periastron)")
    plt.legend(title= 'type of system', loc= 'center')
    plt.xlim(0, len(translated_kepler))
    plt.xlabel("Time (minutes)")
    plt.ylabel("Total image brightness")
    plt.show()

def elliptical_ring_plot(normalizer):
    
    planet_a = 20

    # set up matrix for regular circular orbit
    r_output = ring_array(planet_a)
    r_planet_regular = (r_output > planet_a)
    
    # set up matrix for elliptical ring and planet
    r_theta_output_1 = ellipse_array(planet_a,45,0.6)
    r_theta_output_2 = ellipse_array(planet_a,90,0.6)
    r_ellipse_1 = (1 - (r_theta_output_1>planet_a+30) * (r_theta_output_1<planet_a+50)) * (r_output>planet_a)
    r_ellipse_2 = (1 - (r_theta_output_2>planet_a+30) * (r_theta_output_2<planet_a+50)) * (r_output>planet_a)
    r_planet_ring = (1 - ((ring_array(planet_a) > planet_a + 30) * (ring_array(planet_a) < planet_a + 50))) * (ring_array(planet_a) > planet_a)
    
    translated_kepler = retrieve_coordinates(a, math.radians(90), w_angle, omega, P, T, e, timestep, crop_range, star_radius, star_array_size, r_star)
    
    if normalizer == True:
        max_brightness = star_brightness
    else:
        max_brightness = 1

    brightness_profile_regular = np.array(planet_progression(translated_kepler, star, r_planet_regular, star_brightness, False)) / max_brightness
    brightness_profile_ellipse_1 = np.array(planet_progression(translated_kepler, star, r_ellipse_1, star_brightness, False)) / max_brightness
    brightness_profile_ellipse_2 = np.array(planet_progression(translated_kepler, star, r_ellipse_2, star_brightness, False)) / max_brightness
    brightness_profile_ring = np.array(planet_progression(translated_kepler, star, r_planet_ring, star_brightness, False)) / max_brightness

    plt.plot(brightness_profile_regular, label = 'regular')
    plt.plot(brightness_profile_ring, label = 'ring')
    plt.plot(brightness_profile_ellipse_1, label = 'ellipse 45 deg')
    plt.plot(brightness_profile_ellipse_2, label = 'ellipse 90 deg')
    plt.title("Lightcurve of a transit with limb-darkening (from periastron)")
    plt.legend(title= 'type of system', loc= 'center')
    plt.xlim(0, len(translated_kepler))
    plt.xlabel("Time (minutes)")
    plt.ylabel("Total image brightness")
    plt.show()

def comparison_plot_ellipses(normalizer):

    planet_a = 20

    # set up matrix for regular circular orbit
    r_output = ring_array(planet_a)
    r_planet_regular = (r_output > planet_a)
    
    # set up matrix for elliptical ring and planet
    r_theta_output_1 = ellipse_array(planet_a,45,0.6)
    r_theta_output_2 = ellipse_array(planet_a,90,0.6)
    r_ellipse_1 = (1 - (r_theta_output_1>planet_a+60) * (r_theta_output_1<planet_a+80)) * (r_output>planet_a)
    r_ellipse_2 = (1 - (r_theta_output_2>planet_a+60) * (r_theta_output_2<planet_a+80)) * (r_output>planet_a)
    r_planet_ring = (1 - ((ring_array(planet_a) > planet_a + 60) * (ring_array(planet_a) < planet_a + 80))) * (ring_array(planet_a) > planet_a)
    
    translated_kepler = retrieve_coordinates(a, math.radians(90), w_angle, omega, P, T, e, timestep, crop_range, star_radius, star_array_size, r_star)
    
    if normalizer == True:
        max_brightness = star_brightness
    else:
        max_brightness = 1

    brightness_profile_ellipse_1 = np.array(planet_progression(translated_kepler, star, r_ellipse_1, star_brightness, True)) / max_brightness
    brightness_profile_ellipse_2 = np.array(planet_progression(translated_kepler, star, r_ellipse_2, star_brightness, True)) / max_brightness

    brightness_profile_comparison = brightness_profile_ellipse_1 / brightness_profile_ellipse_2

    #plt.plot(brightness_profile_ellipse_1, label = '45')
    #plt.plot(brightness_profile_ellipse_2, label = '90')
    plt.plot(brightness_profile_comparison, label = '45 deg / 90 deg')
    plt.title("Ratio 45 deg / 90 degree oriented elliptical transits with limb-darkening (from periastron)")
    #plt.legend(title= 'orientation of ellipse (degrees)', loc= 'center')
    plt.xlim(0, len(translated_kepler))
    plt.xlabel("Time (minutes)")
    plt.ylabel("45 degrees / 90 degrees brightness")
    plt.show()

def comparison_plot_circle(normalizer):

    planet_a = 20
    
    # set up matrix for elliptical ring and planet
    r_planet_ring = (1 - ((ring_array(planet_a) > planet_a + 60) * (ring_array(planet_a) < planet_a + 80))) * (ring_array(planet_a) > planet_a)

    total_masking_ring = (len(r_planet_ring) ** 2) - np.sum(r_planet_ring)
    print(len(r_planet_ring))
    print(total_masking_ring)
    
    translated_kepler = retrieve_coordinates(a, math.radians(90), w_angle, omega, P, T, e, timestep, crop_range, star_radius, star_array_size, r_star)
    
    if normalizer == True:
        max_brightness = star_brightness
    else:
        max_brightness = 1

    # determine the equal brightness regular matrix
    mask_radius = int(math.sqrt(total_masking_ring/math.pi))
    r_planet_regular = planet_array(mask_radius) > mask_radius
    print("Cicular")
    print(mask_radius)

    # circular and ring brightness profile
    #brightness_profile_regular = np.array(planet_progression(translated_kepler, star, r_planet_regular, star_brightness, False)) / max_brightness
    #brightness_profile_ring = np.array(planet_progression(translated_kepler, star, r_planet_ring, star_brightness, False)) / max_brightness

    #brightness_profile_comparison = brightness_profile_regular / brightness_profile_ring

    #plt.plot(brightness_profile_regular, label = 'regular')
    #plt.plot(brightness_profile_ring, label = 'ring')
    #plt.plot(brightness_profile_comparison, label = 'regular / ring')
    plt.title("Ratio of regular / ring transit with limb-darkening (from periastron)")
    #plt.legend(title= 'Type of system', loc= 'center')
    plt.xlim(0, len(translated_kepler))
    plt.xlabel("Time (minutes)")
    plt.ylabel("Regular brightness / ring brightness")
    plt.show()

def comparison_plot_multi(normalizer):

    planet_a = 20
    
    # set up matrix for elliptical ring and planet
    r_planet_ring = (1 - ((ring_array(planet_a) > planet_a + 50) * (ring_array(planet_a) < planet_a + 80))) * (1 - ((ring_array(planet_a) > planet_a + 130) * (ring_array(planet_a) < planet_a + 150))) * (1 - ((ring_array(planet_a) > planet_a + 170) * (ring_array(planet_a) < planet_a + 190))) *  (ring_array(planet_a) > planet_a)

    total_masking_ring = (len(r_planet_ring) ** 2) - np.sum(r_planet_ring)
    print(len(r_planet_ring))
    print(total_masking_ring)
    
    translated_kepler = retrieve_coordinates(a, math.radians(90), w_angle, omega, P, T, e, timestep, crop_range, star_radius, star_array_size, r_star)
    
    if normalizer == True:
        max_brightness = star_brightness
    else:
        max_brightness = 1

    # determine the equal brightness regular matrix
    mask_radius = int(math.sqrt(total_masking_ring/math.pi)) +1
    r_planet_regular = planet_array(mask_radius) > mask_radius
    print(mask_radius)
    print(len(r_planet_regular))
    print((len(r_planet_regular) ** 2) - np.sum(r_planet_regular))

    # circular and ring brightness profile
    brightness_profile_regular = np.array(planet_progression(translated_kepler, star, r_planet_regular, star_brightness, False)) / max_brightness
    brightness_profile_ring = np.array(planet_progression(translated_kepler, star, r_planet_ring, star_brightness, True)) / max_brightness

    brightness_profile_comparison = brightness_profile_regular / brightness_profile_ring

    #plt.plot(brightness_profile_regular, label = 'regular')
    #plt.plot(brightness_profile_ring, label = 'multi ring')
    plt.plot(brightness_profile_comparison, label = 'regular / ring')
    plt.title("Ratio of regular / multi ring transit with limb-darkening (from periastron)")
    #plt.legend(title= 'type of system', loc= 'center')
    plt.xlim(0, len(translated_kepler))
    plt.xlabel("Time (minutes)")
    plt.ylabel("regular / multi ring brighntess")
    plt.show()

def no_planet_system(normalizer):
    """planet moves outside the edge of the star, only the ring transits"""
    planet_a = 20

    # set up matrix for regular circular orbit
    r_output = ring_array(planet_a)
    
    # set up matrix for elliptical ring and planet
    r_theta_output_1 = ellipse_array(planet_a,45,0.6)
    r_theta_output_2 = ellipse_array(planet_a,90,0.6)
    r_ellipse_1 = (1 - (r_theta_output_1>planet_a+100) * (r_theta_output_1<planet_a+140)) * (r_output>planet_a)
    r_ellipse_2 = (1 - (r_theta_output_2>planet_a+100) * (r_theta_output_2<planet_a+140)) * (r_output>planet_a)
    
    translated_kepler = retrieve_coordinates(a, math.radians(90.58), w_angle, omega, P, T, e, timestep, crop_range, star_radius, star_array_size, r_star)
    
    if normalizer == True:
        max_brightness = star_brightness
    else:
        max_brightness = 1

    brightness_profile_ellipse_1 = np.array(planet_progression(translated_kepler, star, r_ellipse_1, star_brightness, True)) / max_brightness
    brightness_profile_ellipse_2 = np.array(planet_progression(translated_kepler, star, r_ellipse_2, star_brightness, True)) / max_brightness

    relative_comparison = brightness_profile_ellipse_1 / brightness_profile_ellipse_2

    #plt.plot(brightness_profile_ellipse_1, label = 'ellipse 45 deg')
    #plt.plot(brightness_profile_ellipse_2, label = 'ellipse 90 deg')
    plt.plot(relative_comparison, label = 'ellipse 90 deg')
    plt.title("Ratio of 45 degree / 90 degree no companion system with limb-darkening (from periastron)")
    #plt.legend(title= 'type of system', loc= 'right center')
    plt.xlim(0, len(translated_kepler))
    plt.xlabel("Time (minutes)")
    plt.ylabel("45 degrees / 90 degrees brightness")
    plt.show()


"""
Choose the type of plot you want to generate here
"""
#radii_plot(True)
#inclination_plot(True)
#elliptical_ring_plot(True)
#comparison_plot_ellipses(True)
comparison_plot_circle(True)
#comparison_plot_multi(True)
#no_planet_system(True)
