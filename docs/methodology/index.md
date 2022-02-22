---
layout: page
title: Tutorial
nav_order: 8
has_children: true
---

# ASHLAR: Alignment by Simultaneous Harmonization of Layer/Adjacency Registration

For detailed instructions on how to install ASHLAR, see [installation](https://hxu-hms.github.io/ashlar/methodology/installation.html).

## Stitch & register whole-slide microscopy images in Python

**Ashlar** performs fast, high-quality stitching of microscopy images. It also
co-registers multiple rounds of cyclic imaging for methods such as CyCIF and
CODEX. Ashlar can read image data directly from BioFormats-supported microscope
vendor file formats as well as a directory of plain TIFF files. Output is saved
as pyramidal, tiled OME-TIFF.

*Note:* Ashlar requires unstitched individual "tile" images as input, so it is
not suitable for microscopes or slide scanners that only provide pre-stitched
images.

## Parameters

For detailed information on tuning these parameters, see [Parameter Tuning](https://hxu-hms.github.io/ashlar/methodology/parameter-tuning.html).

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


