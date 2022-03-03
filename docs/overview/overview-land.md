---
layout: default
title: Overview
nav_order: 2
has_children: true
last_modified_date: 2022-02-18
---

# Overview 

{:.no_toc }

<details open markdown="block">
  <summary>
    Table of contents
  </summary>
  {: .text-delta }
  - TOC
{:toc}
</details>

# What's ASHLAR?
**ASHLAR: Alignment by Simultaneous Harmonization of Layer / Adjacency Registration** is a software package that combines multi-tile microscopy images into a high-dimensional mosaic image.

# Why ASHLAR?

Recent imaging advances have enabled scientists to cyclically image pathology samples for 20-60 unique proteins, resulting in a tremendous amount of data per pathology slide. These sample areas are generally quite large, so high-magnification images must ‘tile’ to cover the whole slide. Together, this results in highly multiplexed, biologically-rich, image sets that encompass many sample positions and many proteins. To visualize and analyze this data, the images must be combined to seamlessly align adjacent images within cycles and corresponding locations between cycles. To date, this has been a challenge.  

ASHLAR is a new software package that fills this need and enables multiplexed image alignment with sub-pixel resolution accuracy. 

# How does it work?
ASHLAR combines multiplex images through a 3-step process that:  
>1.  Stitches together adjacent images within the first imaging cycle,  
>2. registers the relative location of images from subsequent cycles to the first cycle, and  
>3. uses these relative positions to generate a final multidimensional mosaic image. 

ASHLAR reads image formats from most commercial microscopes and outputs standard OME-TIFF images. 
ASHLAR is written in Python and is available under the MIT License. Access it at: [https://github.com/labsyspharm/ashlar](https://github.com/labsyspharm/ashlar).

## Step 0: Collect multidimensional image data   
**Collect overlapping image tiles from sample.**

![Grid of overlapping image tiles, with each tile labeled with a consecutive tile identifier that indicates the microscope path. Each tile corresponds to a multi-channel image, which may encompass 4-6 channels of data.]({{ site.baseurl }}/assets/images/Step0.png)

Each image tile is assigned an identifier based on sample location. Location data can be extracted directly from BioFormats image metadata ([Li et al., 2016](https://doi.org/10.1016/j.ymeth.2015.10.006)) that is produced by many commercial microscopes. ASHLAR can also process images from microscopes that do not support BioFormats if they follow consistent naming convention, acquisition order, and tile overlap.

> **Note:** *There are a number of methods to collect multidimensional images. ASHLAR is compatible with many of these methods. For more information on CyCIF, the primary multiplexed image acquisition method used by the Laboratory of Systems Pharmacology at HMS, read the [CyCIF manuscript](https://doi.org/10.7554/eLife.31657) or [protocol](https://dx.doi.org/10.17504/protocols.io.bjiukkew).*

## Step 1: Stitching
**Align adjacent images across a single cycle by using overlapping regions.**

![Grid of overlapping image tiles, with an inset that zooms into the region of overlap between 4 tiles. The first image in the inset shows that the original positions of the overlapped tiles are misaligned. The second image shows the corrected alignment, and final stitched image.]({{ site.baseurl }}/assets/images/Step1.png)


*Use one representative channel and propagate locations to other channels per tile. (DNA marker is recommended.)* 


## Step 2: Registration
**Align corresponding tiles from each subsequent cycle.**

![Top image shows several image tiles from different cycles that correspond to the same position within the sample. An inset shows the overlapping region between the first and second cycle of imaging in one region. The first image of the inset shows misalignment. The second image shows the corrected alignment with altered positions. This altered position is then saved and propagated to other channels within the image tile. This process is repeated for each subsequent image in the cycle.]({{ site.baseurl }}/assets/images/Step2.png)


*Use one representative channel and propagate locations to other channels per tile. (DNA marker is recommended.)* 

## Step 3: Mosaic Image Generation
**All corrected tile positions from Step 1 (stitching) and Step 2 (registration) are combined.**

![Representation of the mosaic image data. Each cycle contains stitched image data for multiple channels, and these aligned cycles are stacked into a single image file.]({{ site.baseurl }}/assets/images/Step3.png)

# Learn More
**View the [detailed computational methods](./DetCompMethods.html) for more information on how each step is performed.**

*For more details, read the preprint manuscript here: [https://doi.org/10.1101/2021.04.20.440625](https://doi.org/10.1101/2021.04.20.440625).*

## Integration with MCMICRO Pipeline
**The Multiple-choice microscopy pipeline (MCMICRO)** is an end-to-end processing pipeline to transform large, multi-channel whole slide images into single-cell data. See [https://mcmicro.org/](https://mcmicro.org/) for more information on the MCMICRO pipeline documentation, implementation, troubleshooting, and more.

![Visual overview of the MC MICRO pipeline components: Basic for illumination correction, Ashlar for alignment and stitching, Coreograph for TMA Core detection, UnMicst or S3 segmenter for segmentation, MC Quant for image quantification.]({{ site.baseurl }}/assets/images/mcmicro-pipeline-two-rows-v2.png)
