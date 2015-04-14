# Phonon-Code
A code to simulate linearized phonon transport

This code is a tool for solving the linearized phonon Boltzmann equation in the relaxation time approximation, in 2D.
In the present version, the user specifies the geometry, the source terms, the "detectors" and various inputs (number of particles, temperature of linearization). The code returns the Monte Carlo estimates obtained with the desired number of particles.

# Input files
The input files are all .txt format. They can be ranked in four categories: base, geometry, sources and detectors. In the current version, the names of the input files should not be changed (only the contents should be adapted by the user).

<h3> - Base input file </h3>
The input file named "input.txt" is just used for specifying:
- The number of particles to be used (first entry)
- The maximum number of collision events that a particle can undergo. By "collision", we mean here any "randomizing event", namely, a scattering event or a diffuse reflection on a boundary. This number is necessary only for <b> steady state </b> calculations <b>with no absorption term</b>, where particles need a termination criterion. Typically, this is the case for the study of 2D periodic systems. In other cases, particle trajectories get "naturally" terminated when they exit the system through absorbing boundaries or because the maximum estimate time has been reached.
- Reference temperature for linearizing the system

<h3>- Geometry files</h3>
As to now, there are three geometry files, each corresponding to a type of boundary condition: "prescribed.txt" (handling prescribed temperature walls, sometimes otherwise called isothermal walls), "periodic.txt" (handling periodic boundary conditions), "reflective.txt" (handling reflective boundaries). Each boundary is implemented as an <b> oriented </b> segment. The orientation of the segment (let's call it vector s) is such that, if vector n is its normal pointing inward the system, s and n must form a direct basis, in the sense that their determinant is positive. Please refer to the provided example files.
- "prescribed.txt" is organized as follows: the first entry is the total number. After this first entry, the segments are specified at each line as "x0 y0 x1 y1 T_0 T_eq"

<h3>- Source files </h3>
There are four source files:
- "volumetric.txt" This corresponds to volumetric heating by Joule effect.
<h3>- Detectors </h3>

# Output format

# Troubleshoot, frequent mistakes and comments
