import numpy as np
import matplotlib.pyplot as plt
import random
import copy
import seaborn as sns

L=10
num_particles = 864
temperature = 1
density=0.85
epsilon=1
sigma=1
h=0.032
mass = 1



#defining the shape and volume of the box
def cube(L):
    volume = L ** 3
    return volume



#applying the periodic boundary condition
def PBC(position, L):
    for i in range(3):  #considering a 3D box (i=x,y,z)
        if position[i] < 0: #agar position zare ma az samte chape mehvare morede nazar(masalan i=x) kharej shod b andaze tule zel e an mehvar b an ezafe kon , dar in sorat engar az rast vared mishavad
            position[i] += L
        elif position[i] >= L:#agar position zare ma az samte raste mehvare morede nazar(masalan i=x) kharej shod b andaze tule zel e an mehvar az an kam kon , dar in sorat engar az chap vared mishavad
            position[i] -= L
    return position


#initializing the positions based on PBC
def initialize_system_with_pbc(num_particles, L):
    # Generate random particle positions within the box
    particle_positions = [
        [L * random.random(), L * random.random(), L * random.random()]
        for _ in range(num_particles)
    ]

    # Apply PBC to ensure particles stay within the box
    for i in range(len(particle_positions)):
        particle_positions[i] = PBC(particle_positions[i], L)

    return np.array(particle_positions)
#initializing the velocities based on the reduced units
def initialize_velocities(num_particles, temperature):
    # Reduced mass and energy units
    mass = 1.0
    epsilon = 1.0

    # Boltzmann constant in reduced units
    k_boltzmann = 1.0

    # Calculate the standard deviation of the Maxwell-Boltzmann distribution
    std_dev = np.sqrt(temperature / k_boltzmann)

    # Generate random velocities from a normal distribution
    velocity = np.random.normal(loc=0, scale=std_dev, size=(num_particles, 3))

    return velocity

#now defining the Lennard Jones potential  &  Force
def Lennard_Jones_Potential(r,epsilon, sigma):
    potential=4*epsilon*((sigma/r)**12-(sigma/r)**6)
    return potential

def lj_force(positions,epsilon,sigma):
    num_particles = len(positions)
    forces = np.zeros_like(positions)
    for i in range(num_particles):
        for j in range(i + 1, num_particles):
            r = np.linalg.norm(positions[i] - positions[j])
            direction = (positions[j] - positions[i]) / r
            force_magnitude = 48 * epsilon * (((sigma)**12/r**13)-1/2*((sigma)**6)/r**7)
            forces[i] += force_magnitude * direction
            forces[j] -= force_magnitude * direction

    return forces

def velocity_verlet(positions, velocities, forces, h, mass):
    half_h = 0.5 * h
    positions += velocities * h + 0.5 * forces * (h ** 2) / mass # r(t+h)
    new_forces = lj_force(positions, epsilon, sigma) # F(t+h)
    velocities += 0.5 * (forces + new_forces) * half_h / mass # v(t+h)
    return positions, velocities, new_forces




particle_positions = initialize_system_with_pbc(num_particles, L)
velocities = initialize_velocities(num_particles, temperature)
forces = lj_force(particle_positions, epsilon, sigma)
positions_list = [copy.deepcopy(particle_positions)]
velocity_list=[copy.deepcopy(velocities)]
forces_list=[copy.deepcopy(forces)]
# Perform a few integration steps
for _ in range(3):
    particle_positions, velocities, forces = velocity_verlet(particle_positions, velocities, forces, h, mass)
    positions_list.append(copy.deepcopy(particle_positions))
    velocity_list.append(copy.deepcopy(velocities))
    forces_list.append(copy.deepcopy(forces))

#Display the first few velocities
print("Initialized Positions :", particle_positions)
# print("Initialized Velocities:",velocities)
# print("Box Volume:",cube(L))
# print(len(particle_positions))
# print(len(velocities))
print(positions_list)

