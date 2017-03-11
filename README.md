# WCHR-regent

## Overview
This code implements the WCHR scheme (see `docs/` for documentation about the scheme) to solve the Navier-Stokes equations for compressible fluid flow.

## Pre-requisites
- The code is written in the [Regent](http://regent-lang.org/ "Regent programming language") programming language based on the [legion](https://github.com/StanfordLegion/legion "Legion runtime system") runtime.

- The code also uses [SuperLU](http://crd-legacy.lbl.gov/~xiaoye/SuperLU/ "SuperLU") to solve block-tridagonal matrices. In `superlu_util.t`, change the path to the SuperLU library.

- HDF5 for file I/O. Make sure Regent is also installed with HDF5 support. In `IO.rg`, change the location to the HDF5 library.

## Running the code
To run the code, first navigate to the `/src/` directory. Ensure that the `problem.rg` symlink is linked to the right problem file. In the problem file, you can set the number of grid points in each direction and set the initial conditions. Also set the `HDF_ROOT` environment variable to the path of your HDF5 library installation to use the file I/O features.

Run the code using `<path/to/regent.py> main.rg`.

## Changing the problem
To change the problem you want to run, simply link the `problem.rg` symlink to a new problem file using `ln -s <path/to/new/problem.rg> problem.rg`.

## TODO
- Add non-periodic boundary conditions
- Parallelization
- Explicit-Compact interpolation
- Implement CFL time step condition
- Different Runge-Kutta time stepping schemes

### Authors
- Akshay Subramaniam (akshays@stanford.edu)
- Man Long Wong (wongml@stanford.edu)
