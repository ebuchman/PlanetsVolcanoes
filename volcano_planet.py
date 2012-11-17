#Script for simulating temperature and density of a planet to evaluate volcanic activity

from scipy import *
from numpy import *
from pylab import *

LOCK_BOUNDARIES = True

#Spherical Shell class represents a spherically symmetrical segment of the planet
class Shell:
        
        temperature = 0
        radius = 0
        conductivity = 0
        density = 0
        thickness = 0

        #Initialize a spherical shell for the planet
        def __init__(self, temp = temperature, r = radius, conductivity = conductivity, rho = density, thickness = thickness):
                self.temperature = temperature
                self.radius = radius
                self.conductivity = conductivity
                self.density = density
                self.thickness = thickness
        
        #Sets the temperature of a shell
        def set_temperature(self, temperature):
                self.temperature = temperature
        
        #Sets the thickness of a shell
        def set_thickness(self, thickness):
                self.thickness = thickness
        
        #Sets the density of a shell
        def set_density(self, density):
                self.density = density

        #Sets the conductivity of a shell
        def set_conductivity(self, conductivity):
                self.conductivity = conductivity
        
        #Gets the temperature of a shell
        def get_temperature(self):
                return self.temperature
                
        #Gets the radius of a shell
        def get_radius(self):
                return self.radius
                
        #Gets the thickness of a shell
        def get_thickness(self):
                return self.thickness

        #Gets the density of a shell
        def get_density(self):
                return self.density

        #Gets the conductivity of a shell
        def get_conductivity(self):
                return self.conductivity
        
        

#Planet composed of spherical shells
class Planet:
        
        shells = []
        num_shells = 0
        radius = 0


        
        #Initialize with empty (parameterless) shells
        def __init__(self, num_shells, radius):
                #Compute thickness of shells in the planet
                shell_thickness = radius / float(num_shells)
                
                #Planet needs at least two shells
                if (num_shells < 2):
                        num_shells = 2
                
                #Initialize an array of shells to fill the planet
                #Shells initialized with T=0, R=i*dr, k=0, rho=0 
                for i in range(0, num_shells):
                        new_shell = Shell(temp = 0, r = shell_thickness * i, conduictivity = 0, rho = 1, thickness = shell_thickness)
                        self.shells.append(new_shell)
                
                #Store the data
                self.num_shells = num_shells
                self.radius = radius
        
        #Recursively compute the new temperature of target shell using the two neighbour shells (inner and outer shells)
        #Target_shell should be the index number of the shell you want (in the shells[] array)
        #Squared_change accumulates how much the temperature overall changed in a given round of recursions
        #Meant to be called with initial target_shell at 0 (i.e. goes from core outwards)
        def compute_init_temperature(self, target_shell, squared_change = 0):
                
                #Set inner and outer shell indices
                inner_shell = target_shell - 1
                outer_shell = target_shell + 1
                
                #Check outer boundary
                if (outer_shell >= self.num_shells and LOCK_BOUNDARIES):
                        #Don't change the temperature at outer edge (it's a boundary condition)
                        new_temperature = self.shells[target_shell].get_temperature()
                #Check inner boundary
                elif (target_shell == 0 and LOCK_BOUNDARIES):
                        new_temperature = self.shells[target_shell].get_temperature()
                else: 
                        #Compute new temperature as a weighted average of the neighbours
                        inner_weighting = 4*pi*(self.shells[inner_shell].get_thickness() + self.shells[inner_shell].get_radius())**2
                        outer_weighting = 4*pi*(self.shells[outer_shell].get_radius())**2
                        norm_factor = (1.0/(inner_weighting + outer_weighting))
                        new_temperature = norm_factor*(inner_weighting*self.shells[inner_shell].get_temperature() + outer_weighting*self.shells[outer_shell].get_temperature())
                
                #Accumulate the temperature change
                squared_change += (self.shells[target_shell].get_temperature() - new_temperature)**2
                
                #Set the new temperature
                self.shells[target_shell].set_temperature(new_temperature)
                
                #Recur the function outward
                if (outer_shell < self.num_shells):
                        squared_change += self.compute_init_temperature(outer_shell)
                
                #This is used in the main loop to determine when the temperature has equilibrated
                return squared_change
        
        #Compute the densities and conductivities of the shells
        def compute_density_and_conductivity(self):
                #Set boundaries for planetary regions (upper bounds of regions in terms of Earth radii)
                inner_core = 0.191
                outer_core = 0.544
                inner_mantle = 0.882
                outer_mantle = .995
                crust = 1.0
                atmosphere =1.002
                
                #Iterate through each shell
                for shell in self.shells:
                        #Find which region of the planet the current shell is in
                        region = shell.get_radius() / self.radius
                        
                        #Apply a different function depending which region of the planet you're in
                        if (region < inner_core):
                               shell.set_density(13.1 - 0.00035*shell.get_radius())
                               shell.set_conductivity(200)
                        elif (region < outer_core):
                               shell.set_density(12.2 - 0.001*shell.get_radius())
                               shell.set_conductivity(5)
                        elif (region < inner_mantle):
                               shell.set_density(5.6 - 0.00055*shell.get_radius())
                               shell.set_conductivity(5)
                        elif (region < outer_mantle):
                                shell.set_density(4.4 - 0.0013*shell.get_radius())
                                shell.set_conductivity(5)
                        elif (region < crust):
                                shell.set_density(2.9 - 0.023*shell.get_radius())
                                shell.set_conductivity(6.5)
                        else:
                                shell.set_density(0.00121 * exp((shell.get_radius() - 1)/0.0012543))
                                shell.set_conductivity(0.025)
        
        #Set temperature of a particular shell
        def set_shell_temperature(self, index, temperature):
                self.shells[index].set_temperature(temperature)
        
        #Get a shell temperature
        def get_shell_temperature(self, index):
                return self.shells[index].get_temperature()

        def get_shell_density(self, index):
                return self.shells[index].get_density()
        
        #Get a shell conductivity
        def get_shell_conductivity(self, index):
                return self.shells[index].get_conductivity()
        
        

                
                
#******Simulate the Planet******#

#Set how many shells we want
num_shells = 100

#Set radius of the planet (in Earth radii)
radius = 1

#Create a planet
planet = Planet(num_shells, radius) 

#Set density of planet
planet.compute_density_and_conductivity()

#Set the initial temp in Kelvins
init_core_temp = 7000
for i in range (0, num_shells - 2):
        planet.set_shell_temperature(i, init_core_temp)

#Set boundary conditions for the temperature initialization
planet.set_shell_temperature(0, init_core_temp)
planet.set_shell_temperature(num_shells - 1, 0)

#Iteratively compute the temperature
change = 10000
change_cond = 1000
while(change > change_cond):
        change = planet.compute_init_temperature(1)

#Check initialized parameters and save to file
file_init_params = open("init_params.txt", "w")
for i in range(0, num_shells):
        s = str(i) + " " + str(planet.get_shell_temperature(i)) + " " + str(planet.get_shell_density(i)) + " " + str(planet.get_shell_conductivity(i)) +  "\n"
        file_init_params.write(s)
file_init_params.close()

#Plot initialized parameters
params = loadtxt("init_params.txt")

plot(params[:,0] , params[:,1])
xlabel("shell (#)")
ylabel("Temp (arb.)")
title("Radial Initial Temperature Profile")
show()

plot(params[:,0] , params[:,2])
xlabel("shell (#)")
ylabel("Density (arb.)")
title("Radial Initial Density Profile")
show()

plot(params[:,0] , params[:,3])
xlabel("shell (#)")
ylabel("Conductivity (arb.)")
title("Radial Initial Conductivity Profile")
show()


