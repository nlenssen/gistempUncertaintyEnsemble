# GISTEMP Uncertainty Ensemble Source Code
Nathan Lenssen, Gavin Schmidt, Michael Hendrickson, Peter Jacobs, Matthew Menne, Reto Ruedy, NASA GISTEMP Team

Current version: v1.0.0

## Outline of Codebase

`Code/` Contains the `R` scripts used to run the analysis

`Figures/` Output of final and intermediate Figures

`Namelists/` Contains parameters for the analysis. Most changes to the behavior of the analysis is done by changing values in the namelist. In particular, the data directory paths are set here (see next section for more information)

`Scripts/` Contains a master script to run the full analysis with one command.

## Directories

`scratchDir/` Contains the majority of the analysis files. This will be 50+ GB so will often not fit in a home directory on a remote cluster. This directory should have `Raw/`, `Intermediate/`, and `Output/` subdirectories.

`ensembleOutDir/` The directory where the final ensemble and related analyses are saved
`plotdir/` The directory figures (in PDF format) are written to

## Data Archive
The full `scratchDir/` has been archived for Lenssen et al. (2024). Two repositories were used due to the 50GB limit:
* [The results and raw data directories](https://zenodo.org/doi/10.5281/zenodo.13343335)
* [The intermediate data products, including the GHCN-ERSST-GISTEMP ensemble](https://zenodo.org/doi/10.5281/zenodo.13344579)

## Recommended Architecture

This analysis was run on an AWS hpc6a.48xlarge (96 core, 384 GB RAM) instance with 300GB of gp2 (SSD) mounted to the instance. As this analysis is often I/O limited, a SSD can dramatically improve the wall-time for running the full analysis. 

The analysis can be run on smaller clusters, but the number of cores used for various steps will need to be adjusted in `Namelists/awsNamelist.Rnl`. The limiting factor will generally be the available memory, not the number of avaiable cores.

A huge thank you to Hoot Thompson, Dorian Crockrell, Dan'l Pierce, Garrison Vaughan, and Dan Duffy at the NASA Center for Climate Simulation (NCCS) for supercomputing and cloud computing support. This research has made use of the NASA Goddard Science Managed Cloud Environment (SMCE), which is a service of the Computational & Information Sciences and Technology Office at the NASA Goddard Space Flight Center.

## GHCN-ERSST-GISTEMP Ensemble
* Currently being run on NASA Discover
* A 100 member ensemble using operational python GISTEMP with input from the GHCN and ERSST uncertainty ensembles
* Contains all relevant sources of uncertainty other than the LSAT sampling uncertainty which will be added using the following R codebase
* Downloadable as part of the [Intermediate repository on Zenodo](https://doi.org/10.5281/zenodo.13344579) for Lenssen et al. (2024)

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

* Create all figures in Lenssen et al. (2024)

## Versions

v1.0.0 Public release of codebase alongside Lenssen et al. (2024)
