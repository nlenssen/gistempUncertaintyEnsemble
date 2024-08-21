# GISTEMP Uncertainty Ensemble Source Code
Nathan Lenssen, Gavin Schmidt, Michael Hendrickson, Reto Ruedy, NASA GISTEMP Team

## Quick start guide

## Outline of Codebase

`Code/` Contains the `R` scripts used to run the analysis

`Figures/` Output of final and intermediate Figures

`Namelists/` Contains parameters for the analysis. Most changes to the behavior of the analysis is done by changing values in the namelist. In particular, the data directory paths are set here (see next section for more information)

`Scripts/` Contains a master script to run the full analysis with one command.

## Data Directories

`scratchDir/` Contains the majority of the analysis files. This will be 50+ GB so will often not fit in a home directory on a remote cluster
`ensembleOutDir/` The directory where the final ensemble and related analyses are saved
`plotdir/` The directory figures (in PDF format) are written to


## GHCN-ERSST-GISTEMP Ensemble
* Currently being run on NASA Discover
* 100 member ensemble of operational python GISTEMP run with input from the GHCN and ERSST uncertainty ensembles
* Contains all relevant sources of uncertainty other than the LSAT sampling uncertainty which will be added using the following R codebase

## Steps

### Step 1: Process ERA Data

* Load and process all of the ERA5 data to be used in the sampling uncertainty estimation (Step 4)

### Step 2: Process GHCN Data

* Load and process the operational GHCNv4 data to determine LSAT historical coverage

### Step 3: Interpolate ERA Data

* Create 1200 km interpolated ERA5 fields on merraGrid NetCDF and convert from native grid to 2x2 grid

### Step 4: Spatial Analysis of Sampling Errors

*  Empirical estimation of sampling error using the reconstruction error fields as in `diffMat` and `dataInds` from `spatialUncertaintyCoarse.R` at 2x2 resolution

### Step 5: Generate Full Ensemble

* Generate full ensemble (ERRST, GHCN + Sampling) NetCDF file at 2x2 grid

### Step 6: Process Ensemble

* Step 6a calculates gridded statistics from the full 200-member ensemble
* Step 6b calculates time-series statistics from the full 200-member ensemble

### Step 7: Run the LSAT decomponsition analyses

This is the analysis used to decompose the LSAT uncertainty into sampling and homogenization components as in 

* Step 6a calculates gridded statistics from the full 200-member ensemble
* Step 6b calculates time-series statistics from the full 200-member ensemble

### Step 8: Create Figures

* Create all figures in Lenssen et al. 2024