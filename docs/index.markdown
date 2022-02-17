---
layout: default
title: Home
nav_order: 1
description: ""
hero_heading: "Alignment by Simultaneous Harmonization of Layer/Adjacency Registration (ASHLAR)"
hero_body: "ASHLAR is a software package that combines multi-tile microscopy images into a high dimensional mosaic image."
hero_ctas:
  - label: "SUMMARY"
    link: "summary.html"
  - label: "DATASET"
    link: "dataset.html"
last_modified_date: 2022-02-17
---

# Key Features

## Precise Stitching for single-cell level accuracy
ASHLAR uses a novel alignment method that results in improved accuracy across imaging cycles. In short, the first cycle of images are combined by aligning the overlapping regions. Next, ASHLAR aligns images between cycles by finding the relative of images from each subsequent cycle to the first cycle image. This results in minimal alignment errors between cycles. 

## Designed for multiplexed images
ASHLAR has been optimized for computationally efficient generation of high dimensional mosaic images.

## Able to handle samples with irregular or discontinuous edges
ASHLAR allows reconstruction of irregular and discontinuous sample edges, as is often seen in biological slide samples.

## Integrated into the Multiple Choice Microscopy Pipeline
For more information about the MCMICRO pipeline, see [mcmicro.org](mcmicro.org).

## Seamless Compatibility
ASHLAR directly reads image formats from most commercial microscopes and slide scanners and outputs standard OME-TIFF images. 

## Open source
ASHLAR is written in Python and is available under the MIT License. It can be downloaded at: [https://github.com/labsyspharm/ashlar](https://github.com/labsyspharm/ashlar).







