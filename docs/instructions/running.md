---
layout: default
nav_exclude: false
title: Running ASHLAR
parent: Instructions
nav_order: 12
---

# Running ASHLAR

This page serves as a reference guide for ASHLAR syntax. Refer to the [quick start guide](./) for step-by-step instructions.

## Description
**Ashlar** performs fast, high-quality stitching of microscopy images. It also co-registers multiple rounds of cyclic imaging for methods such as CyCIF and CODEX. 

## Input
An ```.ome.tiff``` file of **unstitched**\* images

{: .fs-3 }
{: .fw-300 }
> Many BioFormats-supported microscope vendor file formats are supported. Further information on how to determine if your file format is compatible will be added soon. 

{: .fs-3 }
{: .fw-300 }
> \*Ashlar requires unstitched individual "tile" images as input - it is not suitable for microscopes or slide scanners that only provide pre-stitched images. A list of microscopes that we know are compatible with ASHLAR will be coming soon. 

## Output
A pyramidal, tiled ```.ome.tif```

## Usage:
Stitch and align one or more multi-series images
```
ashlar FILE [FILE ...] [OPTIONS] 
```

### Required arguments

| Name | Description |
|---|---|
| ```FILE``` | an image file to be processed (one file per cycle) |

### Optional arguments

|  Name; Shorthand | Description | Default|
|---|---|---|
|```--help; -h```| Show this help message and exit| |
|```--output DIR; -o DIR```|Write image files to DIR|Current directory|
|```--align-channel CHANNEL; -c CHANNEL```| Align images using channel number CHANNEL | Numbering starts at 0|
|```--flip-x```|Flip tile positions left-to-right to account for unusual microscope configurations | |
|```--flip-y```|Flip tile positions top-to-bottom to account for unusual microscope configurations | |
|```--flip-mosaic-x```|Flip mosaic image horizontally||
|```--flip-mosaic-y```|Flip mosaic image vertically||
|```--output-channels CHANNEL [CHANNEL...]```|Output only channels listed in CHANNELS|Numbering starts at 0|
|```--maximum-shift SHIFT; -m SHIFT```|Maximum allowed per-tile corrective shift in microns||
|```--filter-sigma SIGMA```|Width in pixels of Gaussian filter to apply to images before alignment| Default is 0 (which disables filtering)|
|```--filename-format FORMAT; -f FORMAT```|Use FORMAT to generate output filenames, with {cycle} and {channel} as required placeholders for the cycle and channel numbers | default is cycle\_{cycle}\_channel\_{channel}.tif|
|```--pyramid```|Write output as a single pyramidal TIFF||
|```--tile-size PIXELS```|Set tile width and height to PIXELS (pyramid output only)|Default is 1024|
|```--ffp FILE [FILE ...]```|Read flat field profile image from FILES|If specified must be one common file for all cycles or one file for each cycle|
|```--dfp FILE [FILE ...]```|Read dark field profile image from FILES|If specified must be one common file for all cycles or one file for each cycle|
|```--plates```|Enable mode for multi-well plates (for high-throughput screening assays)||
|```--quiet; -q```|Suppress progress display||
|```--version```|Print version||

  
{: .fs-3 }
{: .fw-300 }
**Note:** Detailed information on how to tune these parameters will be added soon.

