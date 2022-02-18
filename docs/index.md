---
layout: default
title: Home
nav_order: 1
description: ""
hero_heading: "ASHLAR"
hero_body: "A software package that combines multi-tile microscopy images into a high dimensional mosaic image."
hero_ctas:
  - label: "OVERVIEW"
    link: "overview-land.html"
  - label: "DATASET"
    link: "dataset.html"
last_modified_date: 2022-02-17
---

# Key Features

## Precise stitching for single-cell level accuracy
ASHLAR uses a novel alignment method that results in improved accuracy across imaging cycles.  
The first cycle of images are combined by aligning and stitching overlapping regions, then images between cycles are aligned by finding the relative position of images from each subsequent cycle to the first cycle image. This minimizes alignment errors between cycles. 

## Designed for multiplexed images
ASHLAR has been optimized for computationally efficient generation of high dimensional mosaic images.

## Supports irregular sample geometries
ASHLAR allows reconstruction of irregular and discontinuous sample edges, as is often seen in biological slide samples.

## Seamlessly compatibile with standard microscopy workflows
ASHLAR directly reads image formats from most commercial microscopes and slide scanners and outputs standard OME-TIFF images. 

## Open source
ASHLAR is written in Python and is available under the MIT License. It can be downloaded at: [https://github.com/labsyspharm/ashlar](https://github.com/labsyspharm/ashlar).

## Integrated into the Multiple Choice Microscopy Pipeline
For more information about the MCMICRO pipeline, see [mcmicro.org](mcmicro.org).







