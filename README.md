# Phonon-Code
A code to simulate linearized phonon transport

This code is a tool for solving the linearized phonon Boltzmann equation in the relaxation time approximation, in 2D.
In the present version, the user specifies the geometry, the source terms, the "detectors" and various inputs (number of particles, temperature of linearization). The code returns the Monte Carlo estimates obtained with the desired number of particles.

This code is Open Source and provided freely. We shall appreciate if scientific work done using this code contains an explicit acknowledgement and the following references which served as a basis for this code:

<p><img src="../../Publications/whiteball.gif" alt="" align="bottom">
J-P. Péraud and N. G. Hadjiconstantinou,
<b>"An alternative approach to efficient simulation of micro/nanoscale phonon transport"</b>,
<em>Applied Physics Letters,</em> <b>101</b>, 153114, 2012.
</p>

<p><img src="../../Publications/whiteball.gif" alt="" align="bottom">
J-P. Péraud and N. G. Hadjiconstantinou,
<b>"Efficient simulation of multidimensional phonon transport using energy-based variance-reduced Monte Carlo formulations"</b>,
<em>Physical Review B,</em> <b>84</b>, 205331, 2011.

# Compiling and running
Typing "make" in the command line should most of the time compile the code and output an executable named "MC.out".
Note that the compiler requires the use of the boost library. For convenience, all the files of the boost library that are needed for compiling are provided (in directory boost_dir/ ).
The complete boost library is available at www.boost.org . I am in no way affiliated to boost.org .
Once compiled, you may run Example 1 by placing the executable in the Example 1 folder and running it from there. The calculation will produce two files (results_T.txt and results_H.txt) as outputs.


# Input files
The input files are all .txt format. They can be ranked in four categories: base, geometry, sources and detectors. In the current version, the names of the input files should not be changed (only the contents should be adapted by the user).

<h3> - Base input file </h3>
The input file named "input.txt" is just used for specifying:
- The type of simulation (FORWARD or ADJOINT)
- The number of particles to be used (first entry)
- The maximum number of collision events that a particle can undergo. By "collision", we mean here any "randomizing event", namely, a scattering event or a diffuse reflection on a boundary. This number is necessary only for <b> steady state </b> calculations <b>with no absorption term</b>, where particles need a termination criterion. Typically, this is the case for the study of 2D periodic systems. In other cases, particle trajectories get "naturally" terminated when they exit the system through absorbing boundaries or because the maximum estimate time has been reached.
- Reference temperature for linearizing the system

The second base input file specifies the material properties. In this version, we provide the properties (assumed isotropic) for silicon. The file is named "dataSi.txt". The first line defines the number of frequency cells. The following lines describe the properties for each of these cells, in the following order: radial frequency - density of states - group velocity - "Delta" of frequency (width of the frequency cell) - relaxation time - polarization . The polarization is coded by 1 for LA, 2 for TA.

<h3>- Geometry files</h3>
As to now, there are three geometry files, each corresponding to a type of boundary condition: "prescribed.txt" (handling prescribed temperature walls, sometimes otherwise called isothermal walls), "periodic.txt" (handling periodic boundary conditions), "reflective.txt" (handling reflective boundaries). Each boundary is implemented as an <b> oriented </b> segment. The orientation of the segment (let's call it vector s) is such that, if vector n is its normal pointing inward the system, s and n must form a direct basis, in the sense that their determinant is positive. Please refer to the provided example files.
- "prescribed.txt" is organized as follows: the first entry is the total number. After this first entry, the segments are specified at each line as "x0 y0 x1 y1 T_0 T_eq", where the first two coordinates define the starting point and x1 and y1 define the ending point of the segment. T_0 refers to the actual temperature of the boundary (for instance, 301 K) and T_eq is the reference temperature (typically, 300 K)
- "periodic.txt". The first line only includes the number of periodic boundaries. The next lines are of the form "x0 y0 x1 y1 xa ya". The first four coordinates define the segment as in "prescribed.txt", while the two last coordinates refers to the vector translation that should applied to the particle when it hits the boundary.
- "reflective.txt": the first line indicates the number of these boundaries. Thean each line is of the form "x0 y0 x1 y1 p". "p" refers to the degree of specularity of the boundary (p=1 indicates a perfectly specular boundary).

<h3>- Source files </h3>
There are three source files. Each of them define sources that are spreads through 2D regions. Each 2D region is defined as a quadrilateral, which are defined in the source files by the coordinates of four points.
- "volumetric.txt": this corresponds to volumetric heating by Joule effect. As in the other files, the first entry defines the number of such sources. Then, at each line, there are 10 numbers. The first 8 numbers describe the coordinates of the points defining the quadrilateral. The quadrilateral must be <b> convex</b>  for the code to work stably. The points must be ordered in the direct direction. In the current version, the two last number of each line are a temperature and the reference temperature. They indicate the rate of temperature increase that should be engendered by the source in a closed system (K per second). For instance, if the last two numbers of a line are the values 310 and 300, then the system undergoes a temperature increase of 10K per second.
- "body_force.txt": this corresponds to the source term that emerges when a spatial variable control is applied to the Boltzmann equation. Seen differently, it corresponds to an "imposed" temperature gradient. The first line indicates the number of such sources, while the following lines include 10 numbers. The first 8 numbers define the quadrilateral, while the two last number define the vector representing the "imposed" temperature gradient. This type of source term is used essentially to simulate periodic systems.
- "initial.txt": this corresponds to the initial condition. It is useful only for transient cases. Again, the first line indicates the number of elements, while the next lines are made of 10 numbers. The 8 first numbers define the quadrilateral, while the last two numbers refer to a temperature T and the reference temperature T_eq. T simply defines the initial temperature in this region.
 
<b> Note: </b> while the regions are primarily coded as quadrilaterals, triangles might be used as well. To do so, one just has to define a quadrialeral with two identical points.

<h3>- Detectors </h3>
Detectors are of two types: temperature detectors and heat flux detectors. They are respectively defined by the files "T_detectors.txt" and "H_detectors.txt". Both files follow the same structure as previously, in the sense that the first number defines the number of such detectors and the following lines define quadrilaterals that define 2D regions. Each estimate is therefore a spatial average over the defined region. In the case of "H_detectors.txt", the basic 8 quadrilater coordinates are completed with two further coordinates. These two coordinates define the vector refering to the components over which the heat flux is projected. For instance, in order to obtain the y-component of the heat flux, these two numbers must be "0 1". In order to obtain the x-component, it must be "1 0".

There is a third detector file called "times.txt". This file refers to the measurement times of the detectors, for transient cases. <b> In order to run a transient calculation, this file must be included in the same folder as the executable. In order to run a steady state calculation, this file must be removed.</b> The first entry of this file must be the number of measurement times. The following entries are the measurement times.

# Output format
In the FORWARD version, the code produces two files:
- "results_T.txt" returns the temperature estimates.
- "results_H.txt" returns the heat flux estimates.
In the BACKWARD version, the code produces:
- "adj_T_results.txt"
- "adj_H_results.txt"
In both cases, the estimates appear as a single column if the calculation is steady state. If the calculation is transient, then each column correspond to a measurement time.

# Examples
There currently are five examples in the repository, respectively in directory Example1, Example1adjoint, Example2, Example2adjoint and thin-film. 

# Visualization
In folder "Tools", the python script "draw.py" should (normally) show you your detectors along with the associated temperature.
After you have run the simulation. Place the python file in your simulation folder, then type "python draw.py".
Limited to steady-state cases for now.
Note: this requires the matplotlib library (www.matplotlib.org) and numpy (www.numpy.org)

# Troubleshoot, frequent mistakes and comments

It is advised to double check if the geometry is well defined and if all boundaries coherently assemble to define the compuational domain. The latest update of the source file Forward_MC.cpp includes a small function which outputs a "geometry.m" file which the user can run in matlab and which will plot the geometry that was entered in the simulation. With this tool, one may check the correct spatial definition of the boundaries, sources and detectors.

# Running an "adjoint" simulation
You need to compile "adjoint_MC.cpp" with "g++ -I boost_dir adjoint_MC.cpp -o adj.out" .
Move the executable adj.out to the folder in which your input files are located.
At line 1 of "input.txt", indicate "ADJOINT". Line 2 is then the total number of particles to be used, line 3 is the max number of collisions and line 4 is the linearization temperature.
IMPORTANT NOTE: this is still an early version, far from being thoroughly tested.


