# Ray-Tracing-Many-Spheres

This code performs Monte Carlo ray tracing to calculate the View Factors or Radiation Distribution Factors (RDF) in large groups of uniform-sized spheres, such as particle beds. 

Operating systems: Built using Ubuntu 16.04, using C++ with MPI. For running on Windows, you can enable “Ubuntu in Windows”. 

The View Factor is more commonly known, being the fraction of rays emitted from surface 1 which strike surface 2; no reflections are accounted for. On the other hand, the RDF is the fraction of rays emitted from surface 1 which are _absorbed_ by surface 2, after any reflections off nearby spheres. Thus, the View Factor is the same as the RDF for the special case of absorptivity = 1, so this code can be used for either. Reflections are modeled as diffuse, though specular reflections may be added in the future. Spheres are opaque (transmissivity = 0), and absorptivity is specified by the user.


How it works:

The code takes the xyz positions of the particle centers (from a csv file) as the main input for ray tracing, which makes it work seamlessly with Discrete Element Method (DEM) particle simulations. DEM is a modeling technique where many spheres are inserted into space and simulated as they collide off each other and any walls. It is a good way to generate realistic domains of particles for studying packed beds or flowing groups of particles. DEM is used in many applications, including powder beds for sintering, chemical processing, nuclear pebble bed reactors, and concentrating solar power. The open source software LIGGGHTS was used in the development of this code, but any other DEM code could be used instead. To study a non-random arrangement of spheres (e.g. simple cubic packing), the sphere center positions could be calculated and written as a text file (with a program such as Matlab or Octave). 

The ray tracing process is as follows: Rays are traced from certain specified "home" (or "emitting") particles. The photon is released from a random position on the emitting sphere's surface, in a random direction in the outward-facing hemisphere from the emission location. The photon is traced until it hits another sphere, at which point the photon either absorbs or reflects from the surface. If the photon is reflected, a new random angle from the reflection location is chosen, and the path of the photon continues. When the photon is eventually absorbed, the photon's path is stopped and counted as an absorption by the absorbing sphere. Then a new photon is released from the emitting sphere, until many photons (typically 10^5 or more) have been emitted. This is repeated for any spheres specified as "home" spheres. At the end, the RDF between two spheres is calculated as the number of absorptions by the receiving sphere divided by the total number emitted from the "home" sphere. 

Optionally, a single wall can also be specified, so the RDF or View Factor from each particle to the wall can be found as well. 

This code uses geometric optics to simulate rays, which is valid for particles which are large in comparison to the wavelength of light. This is true for most macro sized particles (e.g. from marbles down to a fine sand), but it may not be valid for particles only a few microns across. The code is built under the assumption of gray, diffuse emissions and reflections.

This was developed during my PhD work at the Middle East Technical University in Ankara, Turkey. The relevant chapter from the thesis is uploaded to GitHub, which has more details on how to use the code. The code was validated using a simulation in Fluent. Further details can be found in the thesis and in the references below. If you use the code or the equations from it, please cite one of the following papers:

[1] E.F. Johnson, İ. Tarı, D. Baker, Radiative heat transfer in the discrete element method using distance based approximations, Powder Technol. 380 (2021) 164–182. doi:10.1016/j.powtec.2020.11.050.

[2] E. Johnson, İ. Tarı, D. Baker, A Monte Carlo method to solve for radiative effective thermal conductivity for particle beds of various solid fractions and emissivities, J. Quant. Spectrosc. Radiat. Transf. 250 (2020). doi:https://doi.org/10.1016/j.jqsrt.2020.107014. 

The motivation for publishing this code is to provide a tool for researchers and students studying particle systems and powder beds. Others in literature have clearly developed similar codes, but they are generally not publically available. My hope is that others can use this code to develop or validate radiation models for particle beds. 

A note about the code writing: This was my first significant effort in writing in C or C++, which is complicated by using MPI for parallel processing. So it is certainly full of syntax that is not “best practice”, and there are sure to be better ways that this code could have been written. It works, it’s accurate, and it’s fast, but the code itself isn’t pretty. It is, however, very well annotated.

How to use it? 
More details are given in the thesis, but the basic steps are as follows. The only prerequisite in Ubuntu is to install openmpi (with, for example: sudo apt-get install openmpi-bin).
1)	Download the source code and directory structure
2)	Place a text file in the “pos” directory with the xyz positions of the particle centers (often found with DEM, but an example pos file is given)
3)	Place a “home_id” text file in the “pos” directory, specifying which particles will be the home (emitting) particles
4)	Open a terminal in the main RayTracingManySpheres folder
5)	Issue a command to compile the code: mpicxx -std=c++11 PW_MCRT_1.6.cpp -o PW_MCRT_1.6
6)	Run the executable (here, 4 is the number of processors specified): mpirun -np 4 PW_MCRT_1.6
7)	Answer the questions for the simulation parameters: number of photons to send per emitting particle, the total number of particles, the absorptivity, and the radius

The outputs are text files giving: 
1) A square matrix with the RDF values from each sphere to every other sphere
2) A 4-column matrix with columns of Emitting Particle ID, Absorbing Particle ID, Center-to-Center Distance, and RDF.
3) For simulations with a wall specified, a 3-column matrix with columns of Emitting Particle ID, Particle-center to Wall Distance, and Particle-Wall RDF

Recommended additional software: Open the "pos" file in ParaView to visualize the sphere positions. 
