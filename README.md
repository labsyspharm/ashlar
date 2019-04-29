# ASHLAR: Alignment by Simultaneous Harmonization of Layer/Adjacency Registration

## Usage

```
ashlar [-h] [-o DIR] [-c [CHANNEL]]
       [--output-channels [CHANNEL [CHANNEL ...]]] [-m SHIFT]
       [--filter-sigma SIGMA] [-f FORMAT] [--pyramid]
       [--tile-size PIXELS] [--ffp [FILE [FILE ...]]]
       [--dfp [FILE [FILE ...]]] [--plates] [-q] [--version]
       [FILE [FILE ...]]

Stitch and align one or more multi-series images

positional arguments:
  FILE                  an image file to be processed (one file per cycle)

optional arguments:
  -h, --help            show this help message and exit
  -o DIR, --output DIR  write output image files to DIR; default is the
                        current directory
  -c [CHANNEL], --align-channel [CHANNEL]
                        align images using channel number CHANNEL; numbering
                        starts at 0
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
                        cycle_{cycle}_channel_{channel}.tif
  --pyramid             write output as a single pyramidal TIFF
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
investigating a workaround. There is an experimental Docker image on DockerHub
at `sorgerlab/ashlar` which should be suitable for many use cases.
