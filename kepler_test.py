import numpy as np
from matplotlib import pyplot as plt
import math
import copy

 

#Step 1: Calculate the Thiele Innes parameters, independent of eccentricity
def Thiele_Innes(a, i, w, omega):

    A = a*(math.cos(w) * math.cos(omega) - math.sin(w) * math.sin(omega) * math.cos(i))
    B = a*(math.cos(w) * math.sin(omega) + math.sin(w) * math.cos(omega) * math.cos(i))

    F = a*(-math.sin(w) * math.cos(omega) - math.cos(w) * math.sin(omega) * math.cos(i))
    G = a*(-math.sin(w) * math.sin(omega) + math.cos(w) * math.cos(omega) * math.cos(i))

    return A, B, F, G

 
#Step 2: Calculate the mean anomaly M, before you can calculate eccentric anomaly E
def mean_anomaly(P, t, T):

    M = (2*math.pi/P) * (t - T)

    return M

 
#Step 3: calculate the eccentric anomaly E with the mean anomaly 
def eccentric_anomaly(e, P, M):    

    E_0 = M + e * math.sin(M) + e ** 2 / M * math.sin(2*M)
    E_value = [E_0]

    i = 0
    dif = 1

    while dif > 0.01:
        M_index = E_value[i] - e * math.sin(E_value[i])
        E_index = E_value[i] + ((M - M_index) / (1 - e * math.cos(E_value[i])))
        E_value.append(E_index)
        dif = abs(E_value[i + 1] - E_value[i])
        i += 1

    return E_value[-1]

#Step 4: calculate the rectengular coordinates with results form step 1, 2 and 3
def rect_coord(TI, e, E):

    X = math.cos(E) - e
    Y = math.sqrt(1 - e ** 2) * math.sin(E)

    x_pos = TI[0]*X + TI[2]*Y
    y_pos = TI[1]*X + TI[3]*Y
    return x_pos,y_pos

def kepler_cropped(a, i, w, omega, P, T, e, timestep):# loop over timesteps until one transit has passed

    TI = Thiele_Innes(a, i, w, omega)

    coordinates = []
    x_coordinates = []
    y_coordinates = []
        
    for t in np.arange (0, 1, timestep):
        M = mean_anomaly(P, t, T)
        E = eccentric_anomaly(e, P, M)
        
        # This is a little hacky because the first coordinate will never be added alhough it might be moving positive x,
        # but it works for now as long as we choose T accordingly
        if t == 0:    
            old_x = rect_coord(TI, e, E)[0]

        # cut out only the small transit part
        if rect_coord(TI, e, E)[0] > old_x and t != 0 and rect_coord(TI, e, E)[0] > (1.2 * (-2) * (10 ** -7)) and rect_coord(TI, e, E)[0] < (1.2 * (2) * (10 ** -7)):
            
            x_coordinates.append(rect_coord(TI, e, E)[0])
            y_coordinates.append(rect_coord(TI, e, E)[1])
            
            coordinates.append([rect_coord(TI, e, E)[0], rect_coord(TI, e, E)[1]]) # list of x,y coordinate points

            old_x = rect_coord(TI, e, E)[0]

    return x_coordinates, y_coordinates, coordinates

def kepler(a, i, w, omega, P, T, e, timestep):# loop over timesteps until one transit has passed
   
    TI = Thiele_Innes(a, i, w, omega)

    coordinates = []
    x_coordinates = []
    y_coordinates = []
    E_list = []
    M_list = []
    t_list = []

    for t in np.arange (0, 1, timestep):
        M = mean_anomaly(P, t, T)
        E = eccentric_anomaly(e, P, M)
        E_list.append(E)
        M_list.append(M)
        t_list.append(t)


        x_coordinates.append(rect_coord(TI, e, E)[0])
        y_coordinates.append(rect_coord(TI, e, E)[1])
            
        coordinates.append([rect_coord(TI, e, E)[0], rect_coord(TI, e, E)[1]])
        
    plt.plot(t_list, E_list)
    plt.show()
    return x_coordinates, y_coordinates, coordinates

# initial conditions
a = 0.3 #boogseconden
i = math.radians(0)
w = math.radians(0) 
omega = math.radians(0)
P = 1 
T = 0.15 #time passage through periastron (seconden?) #where the planet starts
e = 0.28
timestep = 1./24

# you can choose here if you would like to run the cropped version or the regular orbit
#coord = kepler_cropped(a, i, w, omega, P, T, e, timestep)
coord = kepler(a, i, w, omega, P, T, e, timestep)
x_coordinates = coord[0]
y_coordinates = coord[1]
coordinates = coord[2]

# if __name__ == "__main__":
#     for i in range(0,len(coordinates)):
#         plt.plot(0,0,'yo', markersize=10)
#         plt.plot(x_coordinates[0], y_coordinates[0], 'ro', markersize=5,label='start')
#         plt.plot(x_coordinates[:i],y_coordinates[:i],'bo', markersize=1)
#         plt.plot(x_coordinates[i], y_coordinates[i], "go", markersize = 10)
#         plt.draw()
#         plt.pause(0.001)
#         plt.clf()
#         plt.xlim(-0.6,0.4)
#         plt.ylim(-0.4,0.4)
#     plt.show()