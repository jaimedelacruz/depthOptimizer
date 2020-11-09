# depthOptimizer
A depth optimization tool for 3D rMHD simulations.

## Description
This code precalculates a log tau_500 scale, optimizes the depth-scale for radiative transfer calculations and optionally fills the atom population densities of a 6-level H atom.

The optimization algorithm is based on those found in the MULTI code (Carlsson 1986), but with the addition of a search for gradients in velocity. The code operates in parallel in a share memory machine.

The input/ouput is done via FITS files. The file has dimensions [Ny, Nx, 16, Nz], where optimization is done along the z-axis. The variable ordering along the 3rd axis is: [z, ltau (output), Tg, Vz, Vturb, Bx, By, Bz, Ne, rho, nHI_0, nHI_1, nHI_2, nH1_3, nH1_4, NHII].

All variables should be in CGS units. If nH is blanck, the code will fill it up assuming LTE.

## Usage
```
./depthOpt.x filein.fits fileout.fits nthreads 0 smooth_window temperature_max ltau_max
```

If we assume a machine with 63 cores, and a smooth window of 15 grid cells for the gradients (not for the variables), a temperature cut of 50000 and a maximum continuum log tau of 2:
```
./depthOpt.x filein.fits fileout.fits nthreads 0 15 50000 2
```

## Compilations and dependencies
The code requires a C++-17 compiler: GCC from version 7, clang++ from version 6 and intel compilers from 2019. The code requires the cfitsio library, which is included in most linux distributions. If you have a Mac, it can be installed via, e.g., MacPorts or Homebrew.

Just type make in the src/ folder.

