# GISTEMP Uncertainty Ensemble Source Code
Nathan Lenssen, Gavin Schmidt, Michael Hendrickson, Reto Ruedy, NASA GISTEMP Team

## Quick start guide

## Outline of Codebase

`Code/` Contains the `R` scripts used to run the analysis

`Figures/` Output of final and intermediate Figures

`Namelists/` Contains parameters for the analysis. Most changes to the behavior of the analysis is done by changing values in the namelist. In particular, the data directory paths are set here (see next section for more information)

`Scripts/` (Work in progress) Contains master scripts to run the full analysis with one command.

## Data Directories



## GHCN-ERSST-GISTEMP Ensemble
* Currently being run on Discover
* 100 member ensemble of operational python GISTEMP run with input from the GHCN and ERSST uncertainty ensembles
* Contains all relevant sources of uncertainty other than the LSAT sampling uncertainty which will be added using the following R codebase

## Steps

### Step 1: Process ERA Data

### Step 2: Process GHCN Data

### Step 3: Interpolate ERA Data

* Create 1200 km interpolated ERA5 fields on merraGrid NetCDF and convert from native grid to 2x2 grid

### Step 4: Spatial Analysis of Sampling Errors

*  Calculation of sampling reconstruction error fields as in `diffMat` and `dataInds` from `spatialUncertaintyCoarse.R` at 2x2 resolution

### Step 5: Generate Full Ensemble

* Pull data from discover, generate full ensemble (ERRST, GHCN + Sampling) NetCDF file at 2x2 grid

### Step 6: Process Ensemble

### Step 7: Run the LSAT decomponsition analyses

### Step 8: Make Figures
