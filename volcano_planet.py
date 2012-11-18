#Script for simulating temperature and density of a planet to evaluate volcanic activity
import sys

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
        
                
#******Simulate the Planet******#

#Set radius and mass the planet (in Earth units)
mass = 1.0
radius = mass**(1.0/3.0)

#Set how many shells we want
num_shells = int(250 * radius)

#Set max simulation time
max_time = 10000



#Create a planet
planet = Planet(num_shells, radius) 

#Set density of planet
planet.compute_density_and_conductivity()

#Set the initial temp in Kelvins
init_core_temp = 7000
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

#Open file for temperature tracking
file_temp_track = open("temp_track.txt", "w")


#Open file for time tracking of all parameters
time_track = open('time_track.txt', 'w')


T_data = zeros((num_shells, max_time))

#Iteratively compute the temperature
for i in range (0, max_time):
        #Track temperature of specific shells
        change, temp_list = planet.compute_temperature(0)            
        T_data[:, i] = asarray(temp_list)
        
file_temp_track.close()
        

#Check initialized parameters and save to file
file_init_params = open("final_params.txt", "w")
for i in range(0, num_shells):
        s = str(planet.shells[i].radius) + " " + str(planet.shells[i].temperature) + " " + str(planet.shells[i].density) + " " + str(planet.shells[i].conductivity) +  "\n"
        file_init_params.write(s)
file_init_params.close()

#Plot initialized parameters
params = loadtxt("final_params.txt")

plt.figure(0)

ax1 = plt.subplot(5, 1, 1)
plt.plot(params[:,0] , params[:,1])
ax1.set_xlabel("Radius (Earth radii)")
ax1.set_ylabel("Temp (arb.)")
ax1.set_title("Radial Initial Temperature Profile")

ax2 = plt.subplot(5,1,3)
plt.plot(params[:,0] , params[:,2])
ax2.set_xlabel("Radius (Earth radii)")
ax2.set_ylabel("Density (arb.)")
ax2.set_title("Radial Initial Density Profile")

ax3 = plt.subplot(5,1,5)
plt.plot(params[:,0] , params[:,3])
ax3.set_xlabel("Radius (Earth radii)")
ax3.set_ylabel("Conductivity (arb.)")
ax3.set_title("Radial Initial Conductivity Profile")

plt.savefig('data_r%f_t%d.png'%(radius, max_time))
#plt.show()

plt.clf()

plt.figure(1)

radii = asarray([i*planet.shell_thickness for i in xrange(num_shells)])
time = asarray(range(max_time))

plt.imshow(T_data, aspect = 'auto')
plt.axhline(inner_core_shell, color='black')
plt.axhline(outer_core_shell, color='black')
plt.axhline(inner_mantle_shell, color='black')
plt.axhline(outer_mantle_shell, color='black')
plt.axhline(crust_shell, color='black')

#plt.show()
plt.savefig('heatmap_r%f_t%d.png'%(radius, max_time))
'''
#Track time-evolved temperature of certain shells
inner_core_shell = int(num_shells * planet.inner_core)
outer_core_shell = int(num_shells * planet.outer_core)
inner_mantle_shell = int(num_shells * planet.inner_mantle)
outer_mantle_shell = int(num_shells * planet.outer_mantle)
crust_shell = int(num_shells - 2)

'''



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
