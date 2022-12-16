#!/usr/bin/env python
import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="arc-nonproliferation",
    version="0.1.0-dev",
    description="OpenMC-based API for studying fissile breeding in fusion liquid immersion blankets",
    long_description=long_description,
    packages=setuptools.find_packages()
)