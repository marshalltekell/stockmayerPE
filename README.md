# About this project

Author: Marshall Tekell
Contact: mct2180@columbia.edu

This repository is meant to support a recent publication in Macromolecules (in revision) and to provide examples of particular algorithms for other researchers in the field.

There are really two parts to this project.
-  (1) Generating starting configurations of polymer electrolytes with the appropriate scaling and force-field parameters for molecular dynamics simulations in LAMMPS.  
-  (2) Analyzing the data from the associated simulations in a computationally efficient manner.

The code associated with (1) is contained in ./run. This has four directories — SPE, neat, one_cat, one_an – which correspond to concentrated electrolyte, neat polymer, single cation, and single anion simulations, respectively.

Each of these starting configurations has an example_series folder where one variable (here it is the magnitude of the dipole moment) is systematically varied.

The master.py file allows the user to change this "series variable" and to systematically vary any revelent parameter. For each unique simulation, master.py creates a new subdirectory and populates it with the files need from the templates folder. The files in the templates folder are then altered based on the series variable (and all the other variables contained in master.py).

Then, for each particular force-field implementatation, master.py also generates sbatch scripts for submission to a shared high-performance computing resource. This assumes that LAMMPS is downloaded and available as a binary with the appropriate packages at the location indicated by the pathname. The only required user input is to simply type "sh sub.sh" to to submit all of the production run scripts for the force-field implementations that were systematically varied with a single variable.

Once you have obtained all of the data for a particular force-field, you would import the analysis scripts available in src. 

This scripts make use of SWIG (Simplified Wrapper and Interface Generator) to compile analysis packages written in C and to make them available as shared libraries in Python.

This allows you to, say, calculate the mean-squared displacement with Pythonic syntax but with the speed-up (usually orders of magnitude) afforded by a low-level language like C. For calculations like that partial static structure factor, where you loop over N1 particles, over N2 particles, over 300 q vectors, for let's say 50 different q magnitudes, over T timesteps, these speed-ups are absolutely necessary. 

There are some particulars in getting these things to work together. Specifically, these C functions all operate on numpy arrays in-place and do not return values (i.e. void functions). The syntax for the numpy arrays has to be specified explicitly in the cfunctions.i file. This could be improved in future implementations. 
 
You will need to make the C shared library on your device as the binary is hardware specific (to the best of my knowledge). To do this, simply type "sh swigit.sh" which makes a change to an environment variable (numpy library) and then runs the setup.py script.


If you have any questions, please email me.

 
