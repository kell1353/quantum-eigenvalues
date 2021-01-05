# Overview
This program provides numerical and analytical solutions for the energies at a few different potentials. The infinite well, harmonic 
oscillator and the Woods-Saxon potential. It writes these solutions into files for each potential for a few different states.
 
This program calculates numerically eigenenergies and eigen-wavefunctions for the Woods Saxon potential and provides output for analysis.


# Compilation Instructions
Compilation is all done using the makefile in the repository. Type `make` into your command line to compile the files.

Initially sets woods_saxon as the executable name.
- Then creates the types object file.
- Then creates the eigen_solver object file using the types object file.
- Then creates the hamiltonian object file using the types object file.
- Then creates the qm_solver object file using the types, eigen_solver, and hamiltonian object files.
- Then creates the read_write object file using the types, and qm_solver object files.
- Then creates the main object file using the types, read_write, and qm_solver object files.


# Usage Instructions 
Once you have compiled everything execute the program (./woods_saxon)
The program will prompt you to type in values for the size of the box, number of lattice points and radius of Woods-Saxon potential. Make sure the file is located in the source files.
It will then write the normalized probabilty densities for each scenario into three seperate files.
It will also write a file containing the Woods-Saxon energies as a function over radius.
Once created you will be able to run the Jupyter notebook on these data files. 


# Expected Behavior
Once you have typed the file name, it then will perform calculation, print the analytical and numerical energies comparison and write the results of the program into 
three probability density files and a file with the Woods-Saxon energy results.

The probability density files should have 4 columns. 'x', 'ground state', '1st excited', '2nd excited'.

- x: This column will contain the sample points in the box.
- ground state: This column will contain the ground state probability densities.
- 1st excited: This column will contain the 1st excited probability densities.
- 2nd excited: This column will contain the 2nd excited probability densities.

The other file should have 4 columns. 'radius', 'ground state', '1st excited', '2nd excited'.

- radius: radius of the Woods Saxon potential
- ground state: This column will contain the ground state probability densities.
- 1st excited: This column will contain the 1st excited probability densities.
- 2nd excited: This column will contain the 2nd excited probability densities.
