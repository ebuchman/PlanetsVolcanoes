#Script for simulating temperature and density of a planet to evaluate volcanic activity
import sys, os
import pickle, json

from scipy import *
from numpy import *
from pylab import *

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt


sys.setrecursionlimit(5000)

LOCK_BOUNDARIES = True

#Spherical Shell class represents a spherically symmetrical segment of the planet
class Shell:

        #Initialize a spherical shell for the planet
        def __init__(self, temp, r, conductivity, rho, thickness):
                self.temperature = temp
                self.radius = r
                self.conductivity = conductivity
                self.density = rho
                self.thickness = thickness

#Planet composed of spherical shells
class Planet:
        
        shells = []
        num_shells = 0
        radius = 0

        #Set boundaries for planetary regions (upper bounds of regions in terms of Earth radii)
        inner_core = 0.191
        outer_core = 0.544
        inner_mantle = 0.882
        outer_mantle = .995
        crust = 1.0
        atmosphere =1.002
        
        #Initialize with empty (parameterless) shells
        def __init__(self, num_shells, radius):
                #Compute thickness of shells in the planet
                self.shell_thickness = radius / float(num_shells)
                
                #Planet needs at least two shells
                if (num_shells < 2):
                        num_shells = 2
                
                #Initialize an array of shells to fill the planet
                #Shells initialized with T=0, R=i*dr, k=0, rho=0 
                for i in range(0, num_shells):
                        new_shell = Shell(temp = 0, r = self.shell_thickness * i, conductivity = 0, rho = 1, thickness = self.shell_thickness)
                        self.shells.append(new_shell)
                
                #Store the data
                self.num_shells = num_shells
                self.radius = radius
        
        #Recursively compute the new temperature of target shell using the two neighbour shells (inner and outer shells)
        #Target_shell should be the index number of the shell you want (in the shells[] array)
        #Squared_change accumulates how much the temperature overall changed in a given round of recursions
        #Meant to be called with initial target_shell at 0 (i.e. goes from core outwards)
        def compute_temperature(self, target_shell, T_list = None, squared_change = 0):
                if T_list is None: 
                    T_list = []
                #Set inner and outer shell indices
                inner_shell = target_shell - 1
                outer_shell = target_shell + 1
                
                #Check outer boundary
                if (outer_shell >= self.num_shells and LOCK_BOUNDARIES):
                        #Don't change the temperature at outer edge (it's a boundary condition)
                        new_temperature = self.shells[target_shell].temperature
                        #T_list = []
                #Check inner boundary
                elif (target_shell == 0 and LOCK_BOUNDARIES):
                        new_temperature = self.shells[target_shell].temperature
                        #T_list = []
                else: 
                        #Compute new temperature as a weighted average of the neighbours
                        inner_weighting = self.shells[inner_shell].conductivity*4*pi*(self.shells[inner_shell].thickness + self.shells[inner_shell].radius)**2
                        outer_weighting = self.shells[outer_shell].conductivity*4*pi*(self.shells[outer_shell].radius)**2
                        norm_factor = (1.0/(inner_weighting + outer_weighting))
                        new_temperature = norm_factor*(inner_weighting*self.shells[inner_shell].temperature + outer_weighting*self.shells[outer_shell].temperature)
                
                #Accumulate the temperature change
                squared_change += (self.shells[target_shell].temperature - new_temperature)**2
                
                #Set the new temperature
                self.shells[target_shell].temperature = new_temperature


                #Recur the function outward
                if (outer_shell < self.num_shells):
                        change, _ = self.compute_temperature(outer_shell, T_list)
                        squared_change += change
                        
                T_list.insert(0, new_temperature)

                #This is used in the main loop to determine when the temperature has equilibrated
                return squared_change, T_list
        
        #Compute the densities and conductivities of the shells
        def compute_density_and_conductivity(self):
              
                
                #Iterate through each shell
                for shell in self.shells:
                        #Find which region of the planet the current shell is in
                        region = shell.radius / self.radius
                        
                        #Apply a different function depending which region of the planet you're in
                        if (region < self.inner_core):
                               shell.density = 13.1 - 0.00035*shell.radius
                               shell.conductivity = 200
                        elif (region < self.outer_core):
                               shell.density =12.2 - 0.001*shell.radius
                               shell.conductivity = 10 - 0.006*(shell.radius - 0.191)
                        elif (region < self.inner_mantle):
                               shell.density = 5.6 - 0.00055*shell.radius
                               shell.conductivity = 5
                        elif (region < self.outer_mantle):
                                shell.density = 4.4 - 0.0013*shell.radius
                                shell.conductivity = 5
                        elif (region < self.crust):
                                shell.density = 2.9 - 0.023*shell.radius
                                shell.conductivity = 6.5
                        else:
                                shell.density = 0.00121 * exp((shell.radius - 1)/0.0012543)
                                shell.conductivity = 0.025
        

def experiment(mass = 1.0, base_num_shells = 250, max_time = 10000, init_core_temp = 7000):
    #******Simulate the Planet******#

    # make directory for experimental results
    experiment_name = 'm%.2f_t%d.png'%(mass, max_time)
    dir_name = os.path.normpath(os.path.join('results', experiment_name) )
    os.mkdir(dir_name)
        
    #Set radius and mass the planet (in Earth units)
    radius = mass**(1.0/3.0)

    num_shells = int(base_num_shells* radius)

    #Create a planet
    planet = Planet(num_shells, radius) 

    #Set density of planet
    planet.compute_density_and_conductivity()

    #Set the initial temp in Kelvins
    for i in range (0, num_shells - 2):
            planet.shells[i].temperature = init_core_temp

    #Set boundary conditions for the temperature initialization
    planet.shells[0].temperature = init_core_temp
    planet.shells[num_shells - 1].temperature =  0

    #Track time-evolved temperature of certain shells
    inner_core_shell = int(num_shells * planet.inner_core)
    outer_core_shell = int(num_shells * planet.outer_core)
    inner_mantle_shell = int(num_shells * planet.inner_mantle)
    outer_mantle_shell = int(num_shells * planet.outer_mantle)
    crust_shell = int(num_shells - 2)

    # initialize 2D matrix to store temperature at each radial location over time
    T_data = zeros((num_shells, max_time))

    #Iteratively compute the temperature
    for i in range (0, max_time):
            change, temp_list = planet.compute_temperature(0)            
            T_data[:, i] = asarray(temp_list)
                            
    #Build arrays for parameters and simulated variables
    # density and conductivity are constant
    density = zeros(num_shells)
    conductivity = zeros(num_shells)

    density_array = zeros(num_shells)
    conductivity_array = zeros(num_shells)
    
    for i in range(0, num_shells):
            density_array[i], conductivity_array[i] = planet.shells[i].density, planet.shells[i].conductivity
    radii = asarray([i*planet.shell_thickness for i in xrange(num_shells)])
    time = asarray(range(max_time))
    final_temp = T_data[:, -1]

    # make plots for final values
    plt.figure(0)
    ax1 = plt.subplot(5, 1, 1)
    plt.plot(radii, final_temp)
    ax1.set_xlabel("Radius (Earth radii)")
    ax1.set_ylabel("Temp (arb.)")
    ax1.set_title("Radial Final Temperature Profile")

    ax2 = plt.subplot(5,1,3)
    plt.plot(radii , density)
    ax2.set_xlabel("Radius (Earth radii)")
    ax2.set_ylabel("Density (arb.)")
    ax2.set_title("Radial Density Profile")

    ax3 = plt.subplot(5,1,5)
    plt.plot(radii , conductivity)
    ax3.set_xlabel("Radius (Earth radii)")
    ax3.set_ylabel("Conductivity (arb.)")
    ax3.set_title("Radial Conductivity Profile")

    plt.savefig(os.path.normpath(os.path.join(dir_name, 'final_data.png')))
    #plt.show()

    plt.clf()

    # plot temperate data over time as heat map
    plt.figure(1)
    plt.imshow(T_data, aspect = 'auto')
    plt.axhline(inner_core_shell, color='black')
    plt.axhline(outer_core_shell, color='black')
    plt.axhline(inner_mantle_shell, color='black')
    plt.axhline(outer_mantle_shell, color='black')
    plt.axhline(crust_shell, color='black')
    plt.colorbar()

    plt.show()
    plt.savefig(os.path.join(dir_name, 'heatmap.png'))


    # save parameters and results
    init_params = {'mass': mass , 'base_num_shells' : base_num_shells, 'max_time': max_time, 'init_core_temp': init_core_temp}
    params = [T_data, radii, time, density, conductivity]

    f = open(os.path.join(dir_name, 'results.pkl'), 'wb')
    cPickle.dump(params, f)
    f.close()
    
    f = open(os.path.join(dir_name, 'params'), 'wb')
    json.dump(f, init_params)
    f.close()


experiment()
quit()

radii, time = meshgrid(radii, time)

fig = plt.figure()
ax = fig.gca(projection='3d')

surf = ax.plot_surface(radii, time, T_data)
plt.show()


'''
#Plot temperatures of tracked shells
temps = loadtxt("temp_track.txt")

plot(temps[:,0] , temps[:,1])
xlabel("Time (arb.)")
ylabel("Temp (arb.)")
title("Time Evolved Temperature (Inner Core)")
show()

plot(temps[:,0] , temps[:,2])
xlabel("Time (arb.)")
ylabel("Temp (arb.)")
title("Time Evolved Temperature (Outer Core)")
show()

plot(temps[:,0] , temps[:,3])
xlabel("Time (arb.)")
ylabel("Temp (arb.)")
title("Time Evolved Temperature (Inner Mantle)")
show()

plot(temps[:,0] , temps[:,4])
xlabel("Time (arb.)")
ylabel("Temp (arb.)")
title("Time Evolved Temperature (Crust)")
show()
'''
