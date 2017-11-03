# ASHLAR: Alignment by Simultaneous Harmonization of Layer/Adjacency Registration

## Installation (Linux)

ASHLAR is written in Python, but has python dependencies which make use of C.

The recommended way to install ASHLAR is to make use of a fairly comprehensive
distribution of python binaries called miniconda.

Note: Users with an existing python environment may wish to skip the final step
adding the conda bin directory to the `PATH` as this will override the current
for all shells. Instead, the `PATH` can be added as and when needed.

```bash
# Define an installation directory, e.g.
export CONDA_ROOT=/opt/ashlarconda

# Download Minconda (Python 2) installer
wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh

# Run miniconda installer
sh Miniconda2-latest-Linux-x86_64.sh -b -u -p ${CONDA_ROOT}

# Install conda packages required
${CONDA_ROOT}/bin/conda install -q -y numpy matplotlib networkx pathlib2 \
    scipy scikit-image scikit-learn

# Install conda packages from other channels required
${CONDA_ROOT}/bin/conda install -q -y -c conda-forge pyfftw

# Install pip packages required (not available or up-to-date in conda)
${CONDA_ROOT}/bin/pip install -q -U ModestImage javabridge python-bioformats

# Install ASHLAR
${CONDA_ROOT}/bin/pip install -q -U \
    git+https://github.com/dpwrussell/ashlar@module

# Make this conda install available on your PATH
echo ${CONDA_ROOT}/bin':$PATH' >> ~/.bashrc
```

## Installation (Windows)

TBD

## Usage

```bash
mosaic <scan_path> [--ffc <flat_field_path>]
```
