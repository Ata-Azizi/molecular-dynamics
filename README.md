# Molecular Dynamics Simulation — Lennard-Jones System

A three-dimensional molecular dynamics simulation of interacting particles using the **Lennard-Jones potential**, periodic boundary conditions, and the Velocity Verlet integration algorithm.

This project demonstrates classical many-body simulation, numerical integration, statistical mechanics, and energy-conservation analysis using Python.

## Physical Model

Particles interact through the Lennard-Jones potential:

$$
V(r)=4\epsilon\left[\left(\frac{\sigma}{r}\right)^{12}-\left(\frac{\sigma}{r}\right)^6\right]
$$

where:

- \(r\) is the distance between two particles.
- \(\epsilon\) determines the interaction energy scale.
- \(\sigma\) determines the characteristic particle size.

The corresponding forces are calculated from the negative gradient of the potential.

## Numerical Methods

### FCC Lattice Initialization

Particles are initially arranged on a face-centered cubic (FCC) lattice to avoid unphysical overlaps.

Initial velocities are sampled from a Maxwell–Boltzmann distribution, with center-of-mass motion removed.

### Periodic Boundary Conditions

Periodic boundary conditions approximate an extended bulk system by wrapping particles across the simulation boundaries.

The minimum-image convention is used to calculate pairwise separations.

### Velocity Verlet Integration

Newton's equations of motion are integrated using the Velocity Verlet algorithm:

$$
\mathbf r(t+\Delta t)=\mathbf r(t)+\mathbf v(t)\Delta t+\frac{\mathbf F(t)}{2m}\Delta t^2
$$

$$
\mathbf v(t+\Delta t)=\mathbf v(t)+\frac{\mathbf F(t)+\mathbf F(t+\Delta t)}{2m}\Delta t
$$

This method is widely used in molecular dynamics because of its numerical stability and favorable long-term energy behavior.

## Notebook

**`Molecular_Dynamics_Lennard_Jones.ipynb`**

The notebook includes:

- Three-dimensional FCC lattice generation
- Lennard-Jones force and potential calculations
- Periodic boundary conditions
- Velocity Verlet integration
- Kinetic, potential, and total energy calculations
- Temperature evolution
- Energy-conservation diagnostics
- Three-dimensional particle visualizations

The default configuration uses 108 particles for convenient execution. A larger 864-particle system can be selected by changing the FCC lattice configuration.

## Requirements

- Python 3
- NumPy
- Matplotlib
- Jupyter Notebook or VS Code with Jupyter support

Install dependencies:

```bash
pip install numpy matplotlib jupyter
```

## Running the Simulation

Clone the repository:

```bash
git clone https://github.com/Ata-Azizi/molecular-dynamics.git
cd molecular-dynamics
```

Open `Molecular_Dynamics_Lennard_Jones.ipynb` in Jupyter Notebook or VS Code and execute the cells sequentially.

## Scientific Relevance

This project demonstrates numerical techniques used in computational physics, including many-particle interactions, differential-equation integration, statistical mechanics, and simulation validation.

These techniques provide a foundation for more advanced computational modeling and scientific computing.

## Author

**Ata Azizi Seraji**

M.Sc. Theoretical Physics — University of Göttingen

[GitHub Profile](https://github.com/Ata-Azizi)