#Script for simulating temperature and density of a planet to evaluate volcanic activity
import sys, os
import cPickle, json

from numpy import *

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm

sys.setrecursionlimit(5000)

LOCK_BOUNDARIES = True

SAVE = False
SAVE_PARAMS = False

#Spherical Shell class represents a spherically symmetrical segment of the planet
class Shell:

        #Initialize a spherical shell for the planet
        def __init__(self, temp, r, conductivity, rho, thickness):
                self.temperature = temp
                self.radius = r
                self.conductivity = conductivity
                self.density = rho
                self.diffusivity = 0
                self.thickness = thickness

#Planet composed of spherical shells
class Planet:
        
        shells = []
        num_shells = 0
        radius = 0

        #Set boundaries for planetary regions (upper bounds of regions in terms of Earth radii)
        inner_core = 0.192
        outer_core = 0.546
        inner_mantle = 0.895
        trans_mantle = 0.937
        outer_mantle = .965
        lvz = .996
        crust = 1.0
        
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
                        inner_weighting = self.shells[inner_shell].diffusivity*4*pi*(self.shells[inner_shell].thickness + self.shells[inner_shell].radius)**2
                        outer_weighting = self.shells[outer_shell].diffusivity*4*pi*(self.shells[outer_shell].radius)**2
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
                                shell.density = 13.0885 - 8.8381*shell.radius**2
                                shell.conductivity = 80
                        elif (region < self.outer_core):
                                shell.density = 12.5815 - 1.2638*shell.radius - 3.6426*shell.radius**2 - 5.5281*shell.radius**3
                                shell.conductivity = 10  #-0.01257193448*shell.radius+48.35033200
                        elif (region < self.inner_mantle):
                                shell.density = 7.9565 - 6.4761*shell.radius + 5.5283*shell.radius**2 -3.0807*shell.radius**3
                                shell.conductivity =5
                        elif (region < self.trans_mantle):
                                shell.density = 11.2494- 8.0298*shell.radius
                                shell.conductivity =  5
                        elif (region < self.outer_mantle):
                                shell.density = 7.1089 - 3.8045*shell.radius
                                shell.conductivity = 5
                        elif (region < self.lvz):
                                shell.density = 2.691 + 0.6924*shell.radius
                                shell.conductivity = 5 
                        elif (region <= self.crust):
                                shell.density = 2.75
                                shell.conductivity = 4.5
                        else:
                                shell.density = 2.75
                                shell.conductivity = 4.5
                        
                        shell.diffusivity = shell.conductivity/shell.density

        #Compute mass of planet
        def compute_mass(self):
                mass = 0
                for i in arange(self.num_shells):
                        mass += self.shells[i].density*((self.shells[i].radius*637810000)**2)*4*pi*(self.shells[i].thickness*637810000)
                return mass / 5.97E27
        
        

def experiment(mass = 1.0, base_num_shells = 250, max_time = 100, init_core_temp = 7000, epsilon = 0.01, name = 'Planet'):
    #******Simulate the Planet******#

    # make directory for experimental results
    experiment_name = '%s_m%.2f_t%d'%(name, mass, max_time)
    dir_name = os.path.join('results', 'experiment_'+experiment_name) 
    if not os.path.exists('results'):
        os.mkdir('results')
    if not os.path.exists(dir_name):
        os.mkdir(dir_name)
        
    #Set radius and mass the planet (in Earth units)
    radius = mass**(1.0/3.0)

    num_shells = int(base_num_shells* radius)

    #Create a planet
    planet = Planet(num_shells, radius) 

    #Set density and conductivity of planet
    planet.compute_density_and_conductivity()

    #Set the initial temp in Kelvins
    for i in range (0, num_shells - 2):
            planet.shells[i].temperature = init_core_temp

    #Set boundary conditions for the temperature initialization
    planet.shells[0].temperature = init_core_temp
    planet.shells[num_shells - 1].temperature =  255

    # initialize 2D matrix to store temperature at each radial location over time
    T_data = zeros((num_shells, max_time))

    equilibrium_time = max_time
    
    #Iteratively compute the temperature
    for i in range (0, max_time):
            change, temp_list = planet.compute_temperature(0)            
            T_data[:, i] = asarray(temp_list)
            
            if equilibrium_time == max_time and (change**0.5)/num_shells < epsilon:
                equilibrium_time = i
    
    print "Equilibrium Time: " + str(equilibrium_time)
    print "Total Planet Mass: " + str(planet.compute_mass())
    
    #Build arrays for parameters and simulated variables
    # density and conductivity are constant
    density = zeros(num_shells)
    conductivity = zeros(num_shells)
    for i in range(0, num_shells):
            density[i], conductivity[i] = planet.shells[i].density, planet.shells[i].conductivity
    radii = asarray([i*planet.shell_thickness for i in xrange(0, num_shells, 1)])
    time = asarray(range(0,max_time, 1))
    final_temp = T_data[:, -1]
    
    # plot data and results
    plotData(radii, final_temp, density, conductivity)
    if SAVE: 
        plt.savefig(os.path.join(dir_name, name+': final_data.png'))
        plt.clf()
    else: plt.show()
    
    plotHeatMap(num_shells, planet, T_data, equilibrium_time)
    if SAVE: 
        plt.savefig(os.path.join(dir_name, name+': heatmap.png'))
        plt.clf()
    else: plt.show()

    plotTemperatureSurface(radii, time, T_data, r_step = 1, t_step = int(time.shape[0]/10.0))
    if SAVE:
        plt.savefig(os.path.join(dir_name, name+': temp_surface.png'))
        plt.clf()    
    else: plt.show()


    if SAVE_PARAMS:
        # save parameters and results
        init_params = {'mass': mass, 'base_num_shells' : base_num_shells, 'max_time': max_time, 'init_core_temp': init_core_temp}
        params = [T_data, radii, time, density, conductivity]

        f = open(os.path.join(dir_name, 'results.pkl'), 'wb')
        cPickle.dump(params, f)
        f.close()
        
        f = open(os.path.join(dir_name, 'params'), 'wb')
        json.dump(init_params, f)
        f.close()


#######################
## Plotting Functions##
#######################

def plotData(radii, final_temp, density, conductivity):
    # make plots for final values
    plt.figure(0)
    ax1 = plt.subplot(3, 1, 1)
    plt.plot(radii, final_temp)
    ax1.set_xlabel("Radius (Earth radii)")
    ax1.set_ylabel("Temp (arb.)")
    #ax1.set_title("Radial Final Temperature Profile")

    # a better implementation of the density/conductivity joint plot will be required....
    ax2 = plt.subplot(3,1,3)
    plt.plot(radii, density, 'blue')
    plt.plot(radii, conductivity, 'green')
    ax2.set_xlabel("Radius (Earth radii)")
    ax2.set_ylabel("Density and Conductivity (arb.)")
    #ax2.set_title("Radial Density and Conductivity Profile")

    '''
    ax3 = plt.subplot(5,1,5)
    plt.plot(radii, conductivity)
    ax3.set_xlabel("Radius (Earth radii)")
    ax3.set_ylabel("Conductivity (arb.)")
    ax3.set_title("Radial Conductivity Profile")
    '''

def plotHeatMap(num_shells, planet, T_data, eq_time):

    #Track time-evolved temperature of certain shells
    inner_core_shell = int(num_shells * planet.inner_core)
    outer_core_shell = int(num_shells * planet.outer_core)
    inner_mantle_shell = int(num_shells * planet.inner_mantle)
    trans_mantle_shell = int(num_shells * planet.trans_mantle)
    outer_mantle_shell = int(num_shells * planet.outer_mantle)
    lvz_shell = int(num_shells * planet.lvz)
    crust_shell = int(num_shells - 2)
    
    # plot temperate data over time as heat map
    plt.figure(1)
    plot = plt.imshow(T_data, aspect = 'auto')
    plt.xlabel("Time (Iterations)")
    plt.ylabel("Radius (Spherical Shells)")
    plt.axhline(inner_core_shell, color='black')
    plt.axhline(outer_core_shell, color='black')
    plt.axhline(inner_mantle_shell, color='black')
    plt.axhline(trans_mantle_shell, color='black')
    plt.axhline(outer_mantle_shell, color='black')
    plt.axhline(lvz_shell, color='black')
    plt.axhline(crust_shell, color='black')
    plt.axvline(eq_time, color='black')
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('Temperature (Kelvin)')

    return plot

def plotTemperatureSurface(radii, time, T_data, r_step = 1, t_step = 10):
    T_data=asarray([T_data[:, i] for i in xrange(0, T_data.shape[1], t_step)])
    T_data=asarray([T_data[j, :] for j in xrange(0, T_data.shape[0], r_step)])

    #print time.shape

    radii = asarray([radii[i] for i in xrange(0, radii.shape[0], r_step)])
    time = asarray(range(0,time.shape[0], t_step))
    radii, time = meshgrid(radii, time)


    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.set_xlabel("Radius (Earth Radii)")
    ax.set_ylabel("Time (Iterations)")
    ax.set_zlabel("Temperature (Kelvin)")
    surf = ax.plot_surface(radii, time, T_data, rstride = 1, cstride = 15, cmap =cm.coolwarm)
    
    return surf


###############
##Experiments##
###############

#experiment(mass = 1.0, base_num_shells = 250, max_time = 25000, init_core_temp = 7000,  name = 'Earth_diffusivity')
#experiment(mass = 0.055, base_num_shells = 250, max_time = 25000, init_core_temp = 7000, name = 'Mercury')
#experiment(mass = 0.81, base_num_shells = 250, max_time = 25000, init_core_temp = 7000, name = 'Venus')
#experiment(mass = 0.0123, base_num_shells = 250, max_time = 25000, init_core_temp = 7000, name = 'Moon')
#experiment(mass = 0.107, base_num_shells = 250, max_time = 25000, init_core_temp = 7000, name = 'Mars')

#experiment(max_time = 10000)


experiment(mass = 0.5, base_num_shells = 250, max_time = 25000, init_core_temp = 7000)
#experiment(mass = 1.0, base_num_shells = 250, max_time = 25000, init_core_temp = 7000)
#experiment(mass = 1.0, base_num_shells = 250, max_time = 25000, init_core_temp = 7000, name = 'Planet')
#experiment(mass = 1.5, base_num_shells = 250, max_time = 25000, init_core_temp = 7000, name = 'Planet')
#experiment(mass = 2.0, base_num_shells = 250, max_time = 25000, init_core_temp = 7000, name = 'Planet')
#experiment(mass = 2.5, base_num_shells = 250, max_time = 25000, init_core_temp = 7000, name = 'Planet')
#experiment(mass = 3.0, base_num_shells = 250, max_time = 25000, init_core_temp = 7000, name = 'Planet')
#experiment(mass = 3.5, base_num_shells = 250, max_time = 40000, init_core_temp = 7000, name = 'Planet')

