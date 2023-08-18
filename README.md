# arc-nonproliferation
An OpenMC-based repo for studying proliferation issues in ARC-class fusion reactors.

## Installation and Setup
This repo requires the latest development version of OpenMC (as of 12/24/22) and its depenencies. I am fairly confident that all of the OpenMC dependencies cover the required python packages needed for the rest of the code. To learn how to install OpenMC from source, visit the official website.

The repo is built around a custom python package. After installing OpenMC from source, navigate to the repo directory and run:

`pip install -e .`

to install the custom API localy.

## File Structure and Organization

This repository is built around a custom python api called `arc-nonproliferation` which includes classes, functions, and data which is helpful in conducting analyses of ARC-class reactors. This API was adpoted from one originally developed by Ethan Peterson for the Fall 2022 edition of MIT's course 22.63 Engineering Principles for Fusion Reactors, taught by Dennis Whyte.

### The Device Class
The device class inherits from the OpenMC model class, and is designed to allow for the rapid parametric generation of ARC OpenMC models with FLiBe doped with a specified amount of fertile material. This is achieved using the `generate_device()` function. This allows for many models with varying fertile inventories, dopant types, and lithium enrichments to be generated and simulated in a single script.

### openmc_scripts
The openmc_scripts directory contains all scripts which run OpenMC simulations and analyze the resulting data. The different types of simulations are each in separate directories. They are described below:

| Directory   | Description |
| ----------- | ----------- |
| scan        | scans over fertile inventory with simple transport calculations, no depletion or time effects are resolved.  |
| depletion_scan   | scans over fertile inventory with coupled depletion calculations, evaluating at least one transport calculation at each timestep.         |