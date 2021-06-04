# ASHLAR: Alignment by Simultaneous Harmonization of Layer/Adjacency Registration

## Usage

```
ashlar [-h] [-o DIR] [-c [CHANNEL]] [--flip-x] [--flip-y]
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
conda install -y -c conda-forge numpy scipy matplotlib networkx scikit-image=0.16.2 scikit-learn pyjnius
pip install ashlar
```

### Docker image

The docker image of ashlar is on DockerHub at `sorgerlab/ashlar` which should be
suitable for many use cases.
