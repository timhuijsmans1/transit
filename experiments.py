import numpy as np
import copy
import math
from matplotlib import pyplot as plt

#from kepler_test import coordinates as kepler_output
from kepler_test import kepler_cropped
from matrix_generation import *

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
    
    # set up matrix for ring and planet
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
    brightness_profile_regular = np.array(planet_progression(translated_kepler, star, r_planet_regular, star_brightness, True)) / max_brightness
    brightness_profile_ring = np.array(planet_progression(translated_kepler, star, r_planet_ring, star_brightness, True)) / max_brightness

    #brightness_profile_comparison = brightness_profile_regular / brightness_profile_ring

    plt.plot(brightness_profile_regular, label = 'regular')
    plt.plot(brightness_profile_ring, label = 'ring')
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

def circle_temperatures():
    """planet moves outside the edge of the star, only the ring transits"""
    star_4000 = star_array_limb(star_radius, star_array_size, 4000)
    star_brightness_4000 = np.sum(star_4000)

    star_5000 = star_array_limb(star_radius, star_array_size, 5000)
    star_brightness_5000 = np.sum(star_5000)
    
    planet_a = 20

    # set up matrix for regular circular orbit
    r_output = ring_array(planet_a)
    
    # set up matrix for ring and planet
    r_planet_ring = (1 - ((ring_array(planet_a) > planet_a + 60) * (ring_array(planet_a) < planet_a + 80))) * (ring_array(planet_a) > planet_a)
    
    translated_kepler = retrieve_coordinates(a, math.radians(90), w_angle, omega, P, T, e, timestep, crop_range, star_radius, star_array_size, r_star)

    brightness_profile_ellipse_1 = np.array(planet_progression(translated_kepler, star_4000, r_planet_ring, star_brightness_4000, False))
    brightness_profile_ellipse_2 = np.array(planet_progression(translated_kepler, star_5000, r_planet_ring, star_brightness_5000, False))

    relative_comparison = brightness_profile_ellipse_1 / brightness_profile_ellipse_2

    #plt.plot(brightness_profile_ellipse_1, label = 'ellipse 45 deg')
    #plt.plot(brightness_profile_ellipse_2, label = 'ellipse 90 deg')
    plt.plot(relative_comparison, label = 'ellipse 90 deg')
    plt.title("Ratio of 4000k / 5000K for a circular ring (from periastron)")
    #plt.legend(title= 'type of system', loc= 'right center')
    plt.xlim(0, len(translated_kepler))
    plt.xlabel("Time (minutes)")
    plt.ylabel("Brightness profile ratio for 4000K / 5000K star")
    plt.show()

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
    star = star_array_limb(star_radius, star_array_size, 4000)
    star_brightness = np.sum(star)

    #radii_plot(True)
    #inclination_plot(True)
    #elliptical_ring_plot(True)
    #comparison_plot_ellipses(True)
    comparison_plot_circle(True)
    #comparison_plot_multi(True)
    #no_planet_system(True)
    #circle_temperatures()