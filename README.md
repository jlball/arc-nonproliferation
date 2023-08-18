# arc-nonproliferation
An OpenMC-based repo for studying proliferation issues in ARC-class fusion reactors.

## Installation and Setup
This repo requires the latest development version of OpenMC (as of 12/24/22) and its depenencies. I am fairly confident that all of the OpenMC dependencies cover the required python packages needed for the rest of the code. To learn how to install OpenMC from source, visit the official website.

The repo is built around a custom python package. After installing OpenMC from source, navigate to the repo directory and run:

`pip install -e .`

to install the custom API localy.

## File Structure and Organization

This repository is built around a custom python api called `arc-nonproliferation` which includes classes, functions, and data which is helpful in conducting analyses of ARC-class reactors. This API was adpoted from one originally developed by Ethan Peterson for the Fall 2022 edition of MIT's course 22.63 Engineering Principles for Fusion Reactors, taught by Dennis Whyte.

### arc_nonproliferation

#### device.py
The device class inherits from the OpenMC model class, and is designed to allow for the rapid parametric generation of ARC OpenMC models with FLiBe doped with a specified amount of fertile material. This is achieved using the `generate_device()` function. This allows for many models with varying fertile inventories, dopant types, and lithium enrichments to be generated and simulated in a single script.

### data

This directory contains supporting data files, like outputs from stochastic volume calculations, lists of points which specify vacuum vessel and blanket geometries, and OpenMC chain files.

### openmc_scripts
The openmc_scripts directory contains all scripts which run OpenMC simulations and analyze the resulting data. The different types of simulations are each in separate directories. They are described below:

| Directory   | Description |
| ----------- | ----------- |
| scan        | scans over fertile inventory with simple transport calculations, no depletion or time effects are resolved.  |
| depletion_scan   | scans over fertile inventory with coupled depletion calculations, evaluating at least one transport calculation at each timestep.         |
| independent_depletion  | scans over fertile inventory with independent depletion calculations  |
| geometry_plots | scripts for using OpenMC's built in geometry plotting tools to visualize model geometry |
| depletion | runs a single coupled depletion calculation for a single `device` instance. |

Each directory contains a script which runs a calculation and generates an output directory structure. These scripts take as an argument when run from the command line the name of the directory in which to store the output. There is an accompianing "post-processing" script which has the same name as the script which ran the simulation, but with "pp-" placed before the name. To run the post processing script, pass in the name of the directory you wish to analyse as a command line argument. The figures generated by the post processing script will be stored in a directory called "figures" inside of the main directory for the run.

By default, all output directories are placed in the same directory as the scripts. Some of the scripts, like `independent-depletion-arc-1.py` and `depletion-scan-arc-1.py` take a second command line argument, the enrichment of Li-6 in percent.

Included in each directory is an output folder which contains the data and figures used in the publication.


