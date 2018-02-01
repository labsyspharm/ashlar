# ASHLAR: Alignment by Simultaneous Harmonization of Layer/Adjacency Registration

## Usage

```
ashlar [-h] [-q] [--ffp FFP_FILE] FILE [FILE ...]

positional arguments:
  FILE            an image file to be processed

optional arguments:
  -h, --help      show this help message and exit
  -q, --quiet     suppress progress display
  --ffp FFP_FILE  path to flat field profile image
```

## Installation

### Linux

On Linux, installation is fairly straightforward, but you must first install
your distribution's development package for libfftw3 and a JDK (any one of
versions 1.6, 1.7 or 1.8). You will also have to manually run `pip install
numpy` before `pip install ashlar` due to a requirement of one of the required
packages.

### Using Anaconda Python for MacOS (or Linux)

On MacOS, obtaining and configuring the necessary native libraries can be a bit
challenging. For users on those platforms, or for Linux users having trouble
with the instructions above, the Anaconda Python distribution can simplify the
process.

If you don't already have Anaconda, download it from
https://www.anaconda.com/download/ and install. Then, run the following
commands from a terminal:

```bash
conda install -q -y -c conda-forge pyfftw
pip install -q -U ashlar
```

### Windows

The pyfftw dependency is not currently supported on Windows. We are currently
investigating a workaround.
