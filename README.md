# WCHR-regent

## Overview
This code implements the WCHR scheme (see `docs/` for documentation about the scheme) to solve the Navier-Stokes equations for compressible fluid flow.

## Pre-requisites
- The code is written in the [Regent](http://regent-lang.org/ "Regent programming language") programming language based on the [legion](https://github.com/StanfordLegion/legion "Legion runtime system") runtime.

- [SuperLU](http://crd-legacy.lbl.gov/~xiaoye/SuperLU/ "SuperLU") to solve block-tridagonal matrices.
  - If you're unsure about how to compile SuperLU, make a `build` directory in the SuperLU source directory and from it run `cmake .. -DBUILD_SHARED_LIBS=ON -DUSE_XSDK_DEFAULTS=ON -DCMAKE_INSTALL_PREFIX=<path-to-superlu-installation> -DCMAKE_BUILD_TYPE=Release` and then do `make; make install`.
- HDF5 for file I/O. Make sure Regent is also installed with HDF5 support.

## Compiling the code
To compile the code, first navigate to the `/src/` directory. Ensure that the `problem.rg` symlink is linked to the right problem file in the `src/problems/` directory. In the problem file, you can set the number of grid points in each direction, the time stepping settings and set the initial conditions. Also set the `SUPERLU_PATH` environment variable to the location of the SuperLU installation. To use the HDF5 I/O features, set `USE_IO=1` and set the `HDF_ROOT` environment variable to the path of your HDF5 library installation to use the file I/O features.

Compile the code to an executable using `SAVEOBJ=1 <path/to/regent.py> main.rg`. This will generate the executable called `wchr`. Note that the problem size is a compile time constant and changing the problem size requires recompilation.

## Running the code
If the code is compiled to an executable, run the code using `./wchr -ll:cpu <number of CPUs> [OPTIONS]`. Else, run the code directly using `<path/to/regent.py> main.rg -ll:cpu <number of CPUs> [OPTIONS]`.
```
OPTIONS
  -h               : Print the usage and exit.
  -prefix {prefix} : Use {prefix} as prefix for file I/O.
  -stats {nstats}  : Print stats only every {nstats} steps.
  -p {value}       : Set the number of parallel tasks to {value} (default = 1).
```

## Changing the problem
To change the problem you want to run, simply link the `problem.rg` symlink to a new problem file using `ln -s <path/to/new/problem.rg> problem.rg`.

## TODO
- Add non-periodic boundary conditions (possibly with ghost cells)
- Custom block-tridiagonal solver to remove SuperLU dependency and reduce round-off errors
- Get the SPMD transformation working
- Make the code GPU compatible (at least the interpolation and solve routines)

### Authors
- Akshay Subramaniam (akshays@stanford.edu)
- Man Long Wong (wongml@stanford.edu)
