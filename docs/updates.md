---
layout: default
title: Updates
nav_exclude: false
nav_order: 95

---

<details open markdown="block">
  <summary>
    Table of contents
  </summary>
  {: .text-delta }
  - TOC
{:toc}
</details>

# Software Updates

ASHLAR is periodically being updated to improve functionality. To see up-to-date changes, view the [GitHub release notes](https://github.com/labsyspharm/ashlar/releases).

## Past Updates:

### v1.14.0
#### New features
* Added `--flip-mosaic-x` and `--flip-mosaic-y` arguments to allow flipping output mosaic image. (#90) 

### v1.13.1
#### Internal enhancements
* Update README to include `--flip-x` and `--flip-y` arguments in the usage summary and refresh the installation instructions. (#83, #89)
* Allow `ashlar.reg.plot_*` functions that can display the stitched image to use a downsampled image. This makes these functions practical for large images which would cause `matplotlib.imshow` to struggle. (#86)
* The `preview_slide` script now downsamples images by a factor of 10 by default for improved performance. The factor can be adjusted with the `-d` argument. (#86)

### v1.13.0
#### New features
* FileSeriesReader adds new parameters `layout` and `direction`. Layout can be set to `raster` (default) or `snake`. Direction can be set to `horizontal` (default) or `vertical`. (#76)
#### Internal enhancements
* Add layer quality plot `ashlar.reg.plot_layer_quality`. (#81)
#### Bug fixes
* Pin random state by default for reproducible stitching results. (#77)
* Multi-tile .CZI files and .CZI files with extra images (such as slide labels) are now read correctly. (#80)
* Catch rare error caused when assembling mosaic. (#79)

# Website Updates

This website was last updated at the time of publication (February 2022).
