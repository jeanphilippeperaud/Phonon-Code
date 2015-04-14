# Phonon-Code
A code to simulate linearized phonon transport

This code is a tool for solving the linearized phonon Boltzmann equation in the relaxation time approximation, in 2D.
In the present version, the user specifies the geometry, the source terms, the "detectors" and various inputs (number of particles, temperature of linearization). The code returns the Monte Carlo estimates obtained with the desired number of particles.

# Input files
The input files are all .txt format. They can be ranked in three categories: geometry, sources and detectors.
<h3>- Geometry files</h3>
As to now, there are three geometry files, each corresponding to a type of boundary condition: "prescribed.txt" (handling prescribed temperature walls, sometimes otherwise called isothermal walls), "periodic.txt" (handling periodic boundary conditions), "reflective.txt" (handling reflective boundaries).

<h3>- Source files </h3>

<h3>- Detectors </h3>
