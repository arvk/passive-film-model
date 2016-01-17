passive-film-model
==================
Note: This is a work in progress

The passive film model is a coupled kinetic Monte Carlo (kMC) and phase field (PF) method to simulate the multiscale behavior of multi-phase passive films in electrochemical environments.


Dependencies
============
* [autotools](http://www.gnu.org/software/autoconf/autoconf.html) for building the program
* [PETSc](http://www.mcs.anl.gov/petsc/) (v3.2-3.5) for solving phase-field equations
* Custom [fork](https://github.com/arvk/spparks-pfm) of [SPPARKS](http://spparks.sandia.gov/), compiled as a library, for performing 2D kMC modeling at interfaces between phases
* [VisIt](https://wci.llnl.gov/simulation/computer-codes/visit/) for visualizing output .vts files.


Compiling and Running
=====================
1. Run autoreconf using `autoreconf -i`
2. Run the configure script and provide the path the the `SPPARKS` library using the `--with-spparks` flag
   * If the PETSc installation is successful, the environmental variables `PETSC_DIR` and `PETSC_ARCH` are available. In this case, run `./configure --with-spparks=<PATH_TO_SPPARKS>`
   * If these environment variables are not currently defined, you can define them in the configure command using the `./configure --with-spparks=<PATH_TO_SPPARKS> --with-petscdir=<PATH_TO_PETSC> --with-petscarch=<PETSC_ARCH>`
3. Run `make`

The executable, `pfm` is in `src/`

Documentation
=============
The code is partly documented using [FORD](https://github.com/cmacmackin/ford). To compile the documentation, from the main directory, run

`ford sample_parameters/ford_project_file`

More extensive documentation is on the way.

Contributing
============
Please report any errors and issues using the issues tab. 

To do
=====
1. Update code for PETSc 3.6
2. Document literature sources for simulation parameters

License
=======
The passive film model source code and documentation are distributed 'as is' under the [GNU General Public License v2](doc/LICENSE).
