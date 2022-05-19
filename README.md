# ASHLAR: Alignment by Simultaneous Harmonization of Layer/Adjacency Registration

## Whole-slide microscopy image stitching and registration in Python

**Ashlar** performs fast, high-quality stitching of microscopy images. It also
co-registers multiple rounds of cyclic imaging for methods such as CyCIF and
CODEX. Ashlar can read image data directly from BioFormats-supported microscope
vendor file formats as well as a directory of plain TIFF files. Output is saved
as pyramidal, tiled OME-TIFF.

Note that Ashlar requires unstitched individual "tile" images as input, so it is
not suitable for microscopes or slide scanners that only provide pre-stitched
images.

**Visit [labsyspharm.github.io/ashlar/](https://labsyspharm.github.io/ashlar/) for the most up-to-date information on ASHLAR.**

## Usage

```
ashlar [-h] [-o PATH] [-c CHANNEL] [--flip-x] [--flip-y]
       [--flip-mosaic-x] [--flip-mosaic-y]
       [--output-channels CHANNEL [CHANNEL ...]] [-m SHIFT]
       [--filter-sigma SIGMA] [--tile-size PIXELS]
       [--ffp FILE [FILE ...]] [--dfp FILE [FILE ...]] [--plates] [-q]
       [--version]
       FILE [FILE ...]

Stitch and align multi-tile cyclic microscope images

positional arguments:
  FILE                  Image file(s) to be processed, one per cycle

optional arguments:
  -h, --help            Show this help message and exit
  -o PATH, --output PATH
                        Output file. If PATH ends in .ome.tif a pyramidal OME-
                        TIFF will be written. If PATH ends in just .tif and
                        includes {cycle} and {channel} placeholders, a series
                        of single-channel plain TIFF files will be written. If
                        PATH starts with a relative or absolute path to
                        another directory, that directory must already exist.
                        (default: ashlar_output.ome.tif)
  -c CHANNEL, --align-channel CHANNEL
                        Reference channel number for image alignment.
                        Numbering starts at 0. (default: 0)
  --flip-x              Flip tile positions left-to-right
  --flip-y              Flip tile positions top-to-bottom
  --flip-mosaic-x       Flip output image left-to-right
  --flip-mosaic-y       Flip output image top-to-bottom
  --output-channels CHANNEL [CHANNEL ...]
                        Output only specified channels for each cycle.
                        Numbering starts at 0. (default: all channels)
  -m SHIFT, --maximum-shift SHIFT
                        Maximum allowed per-tile corrective shift in microns
                        (default: 15)
  --filter-sigma SIGMA  Filter images before alignment using a Gaussian kernel
                        with s.d. of SIGMA pixels (default: no filtering)
  --tile-size PIXELS    Pyramid tile size for OME-TIFF output (default: 1024)
  --ffp FILE [FILE ...]
                        Perform flat field illumination correction using the
                        given profile image. Specify one common file for all
                        cycles or one file for every cycle. Channel counts
                        must match input files. (default: no flat field
                        correction)
  --dfp FILE [FILE ...]
                        Perform dark field illumination correction using the
                        given profile image. Specify one common file for all
                        cycles or one file for every cycle. Channel counts
                        must match input files. (default: no dark field
                        correction)
  --plates              Enable plate mode for HTS data
  -q, --quiet           Suppress progress display
  --version             Show program's version number and exit
```

## Installation

### Pip install

Ashlar can be installed in most Python environments using `pip`:
``` bash
pip install ashlar
```

### Using a conda environment

If you don't already have [miniconda](https://docs.conda.io/en/latest/miniconda.html)
or [Anaconda](https://www.anaconda.com/products/individual), download the python
3.x version and install. Then, run the following commands from a terminal (Linux/Mac)
or command prompt (Windows):

Create a named conda environment with python 3.10:
```bash
conda create -y -n ashlar python=3.10
```

Activate the conda environment:
```bash
conda activate ashlar
```

In the activated environment, install dependencies and ashlar itself:
```bash
conda install -y -c conda-forge numpy scipy matplotlib networkx scikit-image=0.19 scikit-learn "tifffile>=2022.4.8" zarr pyjnius blessed
pip install ashlar
```

### Docker image

The docker image of ashlar is on DockerHub at `labsyspharm/ashlar` which should be
suitable for many use cases.
