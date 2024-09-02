---
layout: default
nav_exclude: false
title: Installation
parent: Instructions
nav_order: 10
---

# Installation

### Pip install

Ashlar can be installed in most Python environments using `pip`:
``` bash
pip install ashlar
```

### Using a conda environment

If you don't already have [miniconda](https://docs.conda.io/en/latest/miniconda.html)
or [Anaconda](https://www.anaconda.com/products/individual), download Anaconda and
install. Then, run the following commands from a terminal (Linux/Mac) or Anaconda
command prompt (Windows):

Create a named conda environment with python 3.12:
```bash
conda create -y -n ashlar python=3.12
```

Activate the conda environment:
```bash
conda activate ashlar
```

In the activated environment, install dependencies and ashlar itself:
```bash
conda install -y -c conda-forge numpy scipy matplotlib networkx scikit-image scikit-learn tifffile zarr pyjnius blessed
pip install ashlar
```

### Docker image
The docker image of ashlar is on DockerHub at [labsyspharm/ashlar](https://hub.docker.com/r/labsyspharm/ashlar) and should be suitable for many use cases.

**Return to the [quick start guide](./) to learn more about how to use ASHLAR.**
