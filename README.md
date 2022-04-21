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
       [--output-channels [CHANNEL [CHANNEL ...]]] [-m SHIFT]
       [--filter-sigma SIGMA] [-f FORMAT] [--pyramid]
       [--tile-size PIXELS] [--ffp [FILE [FILE ...]]]
       [--dfp [FILE [FILE ...]]] [--plates] [-q] [--version]
       FILE [FILE ...]

Stitch and align one or more multi-series images

positional arguments:
  FILE                  an image file to be processed (one file per cycle)

optional arguments:
  -h, --help            show this help message and exit
  -o PATH, --output PATH
                        write output to PATH; default is
                        ashlar_output.ome.tif. If value ends in .ome.tif an
                        OME-TIFF with tiled image pyramid will be written. If
                        value ends in just .tif and includes {cycle} and
                        {channel} placeholders a series of single-channel TIFF
                        files will be written. Otherwise value will be
                        interpreted as a directory and the '-f' and '--
                        pyramid' arguments will control the file names and
                        format.
  -c CHANNEL, --align-channel CHANNEL
                        align images using channel number CHANNEL; numbering
                        starts at 0
  --flip-x              flip tile positions left-to-right to account for
                        unusual microscope configurations
  --flip-y              flip tile positions top-to-bottom to account for
                        unusual microscope configurations
  --flip-mosaic-x       flip output image horizontally
  --flip-mosaic-y       flip output image vertically
  --output-channels [CHANNEL [CHANNEL ...]]
                        output only channels listed in CHANNELS; numbering
                        starts at 0
  -m SHIFT, --maximum-shift SHIFT
                        maximum allowed per-tile corrective shift in microns
  --filter-sigma SIGMA  width in pixels of Gaussian filter to apply to images
                        before alignment; default is 0 which disables
                        filtering
  -f FORMAT, --filename-format FORMAT
                        use FORMAT to generate output filenames, with {cycle}
                        and {channel} as required placeholders for the cycle
                        and channel numbers; default is
                        cycle_{cycle}_channel_{channel}.tif (DEPRECATED: Use
                        the '-o' argument to specify the output filename
                        format.)
  --pyramid             write output as a single pyramidal OME-TIFF
  --tile-size PIXELS    set tile width and height to PIXELS (pyramid output
                        only); default is 1024
  --ffp [FILE [FILE ...]]
                        read flat field profile image from FILES; if specified
                        must be one common file for all cycles or one file for
                        each cycle
  --dfp [FILE [FILE ...]]
                        read dark field profile image from FILES; if specified
                        must be one common file for all cycles or one file for
                        each cycle
  --plates              enable plate mode for HTS data
  -q, --quiet           suppress progress display
  --version             print version
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

Create a named conda environment with python 3.7:
```bash
conda create -y -n ashlar python=3.7
```

Activate the conda environment:
```bash
conda activate ashlar
```

In the activated environment, install dependencies and ashlar itself:
```bash
conda install -y -c conda-forge numpy scipy matplotlib networkx scikit-image=0.16.2 scikit-learn tifffile zarr pyjnius=1.2.1
pip install ashlar
```

### Docker image

The docker image of ashlar is on DockerHub at `labsyspharm/ashlar` which should be
suitable for many use cases.
