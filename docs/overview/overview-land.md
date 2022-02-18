---
layout: default
title: Overview
nav_order: 2
has_children: true
last_modified_date: 2022-02-18
---

# Overview

### Why ASHLAR?
{: .no_toc }

Recent imaging advances have enabled scientists to cyclically image pathology samples for 20-60 unique proteins, resulting in a tremendous amount of data per pathology slide. These sample areas are generally quite large, so high-magnification images must ‘tile’ to cover the whole slide. Together, this results in highly multiplexed, biologically-rich, image sets that encompass many sample positions and many proteins. To visualize and analyze this data, the images must be combined to seamlessly align adjacent images within cycles and corresponding locations between cycles. To date, this has been a challenge.  
ASHLAR is a new software package that fills this need and enables multiplexed image alignment with sub-pixel resolution accuracy. 

### How does it work?
ASHLAR combines multiplex images through a 3-step process that:  
i) stitches together adjacent images within the first imaging cycle,  
ii) registers the relative location of images from subsequent cycles to the first cycle, and  
iii) uses these relative positions to generate a final multidimensional mosaic image. 

ASHLAR reads image formats from most commercial microscopes and outputs standard OME-TIFF images.  
ASHLAR is written in Python and is available under the MIT License. Download ASHLAR at: [https://github.com/labsyspharm/ashlar](https://github.com/labsyspharm/ashlar).

## Step 0: Collect multidimensional image data   
**Collect overlapping image tiles from sample.**

![Grid of overlapping image tiles, with each tile labeled with a consecutive tile identifier that indicates the microsope path. Each tile corresponds to a multi-channel image, which may encompass 4-6 channels of data.]({{ site.baseurl }}/assets/images/Step0.png)

Each image tile is assigned an identifier based on sample location. Location data can be extracted directly from BioFormats image metadata ([Li et al., 2016](https://doi.org/10.1016/j.ymeth.2015.10.006)) that is produced by many commercial microscopes. ASHLAR can also process images from microscopes that do not support BioFormats if they follow consistent naming convention, acquisition order, and tile overlap.

*For more information on how multidimensional images are collected experimentally, see [https://dx.doi.org/10.17504/protocols.io.bjiukkew](https://dx.doi.org/10.17504/protocols.io.bjiukkew).*

## Step 1: Stitching
**Align adjacent images across a single cycle by using overlapping regions.**

![Grid of overlapping image tiles, with an inset that zooms into the region of overlap between 4 tiles. The first image in the inset shows that the original postitions of the overlapped tiles are misaligned. The second image shows the corrected alignment, and final stitched image.]({{ site.baseurl }}/assets/images/Step1.png)


*Use one representative channel and propagate locations to other channels per tile. (DNA marker is recommended.)* 


## Step 2: Registration
**Align corresponding tiles from each subsequent cycle.**

![Top image shows several image tiles from different cycles that correspond to the same position within the sample. An inset shows the overlapping region between the first and second cycle of imaging in one region. The first image of the inset shows misalignment. The second image shows the corrected alignment with altered positions. This altered position is then saved and propogated to other channels within the image tile. This process is repeated for each subsequent image in the cycle.]({{ site.baseurl }}/assets/images/Step2.png)


*Use one representative channel and propagate locations to other channels per tile. (DNA marker is recommended.)* 

## Step 3: Mosaic Image Generation
**All corrected tile positions from Step 1 (stitching) and Step 2 (registration) are combined.**

![Representation of the mosaic image data. Each cycle contains stiched image data for multiple channels, and these aligned cycles are stacked into a single image file.]({{ site.baseurl }}/assets/images/Step3.png)

## Learn More
**View the [detailed computational methods](./overview/DetCompMethods.html) for more information on how each step is performed.**

*For more details, read the pre-print manuscript here: [https://doi.org/10.1101/2021.04.20.440625](https://doi.org/10.1101/2021.04.20.440625).*

<details open markdown="block">
  <summary>
    Table of contents
  </summary>
  {: .text-delta }
1. TOC
{:toc}
</details>
