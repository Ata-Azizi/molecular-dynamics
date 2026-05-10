import numpy as np
import matplotlib.pyplot as plt

# Parameters
n_particles = 864
density = 0.85
temperature = 2.0
time_steps = 15.0
time_step_size = 0.032
cutoff = 2.5
sigma = 1.0

# Initialization
box_size = (n_particles / density) ** (1 / 3)
positions = np.zeros((n_particles, 3))
velocities = np.zeros((n_particles, 3))

# Initialize positions with FCC lattice
lattice_spacing = box_size / (4 ** (1 / 3))
particle_count = 0
for x in range(4):
    for y in range(4):
        for z in range(4):
            if particle_count < n_particles:
                positions[particle_count] = np.array([x, y, z]) * lattice_spacing
                particle_count += 1

# Initialize velocities with Maxwell distribution
velocities = np.random.normal(0, np.sqrt(temperature), size=(n_particles, 3))


# Velocity Verlet integration
def lj_force(r, sigma, cutoff):
    r_magnitude = np.linalg.norm(r, axis=-1)
    r_mask = r_magnitude < cutoff
    r_cut = r[r_mask]
    r_magnitude_cut = r_magnitude[r_mask]

    force_magnitude = 48 * ((sigma / r_magnitude_cut) ** 12 - 0.5 * (sigma / r_magnitude_cut) ** 6) / r_magnitude_cut
    force_magnitude[np.isnan(force_magnitude)] = 0  # Handle division by zero

    forces = np.zeros_like(r)
    for i in range(3):
        forces[r_mask, i] = force_magnitude * r_cut[:, i] / r_magnitude_cut

    return forces


# Simulation
num_steps = int(time_steps / time_step_size)
trajectory = np.zeros((num_steps, n_particles, 3))

for step in range(num_steps):
    # Velocity Verlet algorithm
    positions += time_step_size * (velocities + 0.5 * time_step_size * lj_force(positions, sigma, cutoff))
    velocities += 0.5 * time_step_size * lj_force(positions, sigma, cutoff)
    velocities += 0.5 * time_step_size * lj_force(positions, sigma, cutoff)

    # Apply periodic boundary conditions
    positions = np.mod(positions, box_size)

    # Store the trajectory for analysis or visualization
    trajectory[step] = positions.copy()

# Visualization (plot the trajectory)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_box_aspect([1, 1, 1])

for particle in range(n_particles):
    ax.plot(trajectory[:, particle, 0], trajectory[:, particle, 1], trajectory[:, particle, 2])

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.show()



