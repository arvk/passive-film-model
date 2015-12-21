passive-film-model
==================
Note: This is a work in progress

The passive film model is a coupled kinetic Monte Carlo (kMC) and phase field (PF) method to simulate the multiscale behavior of multi-phase passive films in electrochemical environments.


Dependencies
============
* [PETSc](http://www.mcs.anl.gov/petsc/) 3.2 for solving phase-field equations
* [SPPARKS](http://spparks.sandia.gov/), compiled as a library, for performing 2D kMC modeling at interfaces between phases
* [VisIt](https://wci.llnl.gov/simulation/computer-codes/visit/) for visualizing output .vts files.


Compiling and Running
=====================
After choosing an appropriate Fortran compiler in the Makefile and verifying the values of environment variables `PETSC_DIR` and `SPPARKS_LIB_DIR`, compile with GNU Make using `make -f Makefile` and run using `./model.exe`

Documentation
=============
The code is partly documented using [FORD](https://github.com/cmacmackin/ford). To compile the documentation, from the main directory, run

`ford sample_parameters/ford_project_file`

More extensive documentation is on the way.

Contributing
============
Please report any errors and issues [here](issues)

To do
=====
1. Update code for PETSc 3.6
2. Document literature sources for simulation parameters

License
=======
The passive film model source code and documentation are distributed 'as is' under the [GNU General Public License v2](doc/LICENSE).
