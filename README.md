# POLRAD 2.0 Radiative Corrections

FORTRAN code for radiative correction calculations used in CLAS12 analysis, specifically for inclusive DIS electron scattering on nuclear targets.

## Attribution and Contact
**Original Developer:** Dr. Alexander Ilyichev <ily@hep.by>
- This code was provided by Alexander Ilyichev for radiative corrections analysis of CLAS12 RG-D data
- Includes setup for Carbon target with example input file

## Purpose
Calculate internal radiative corrections for:
- Electron-proton scattering (DIS continuum + elastic radiative tail)
- Electron-nucleus scattering (DIS continuum + elastic tail + quasi-elastic radiative tail)
- Supports analysis of LD2 and nuclear targets (C, Cu, Sn) for CLAS12 RG-D experiment

## Files
- `polrad20.f` - Main Fortran 77 source code
- `Makefile` - Build configuration
- `input.dat` - Input parameters and configuration (example for Carbon)

## Building and Running
```bash
make
./polrad20 < input.dat
