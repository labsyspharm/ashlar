# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

---

## [afr-2024-9-1] - 2024-09-14

### Added

- `rcashlar combine` command: : Added option to flip the combined image
  vertically and horizontally. The default behavior is to flip horizontally.

### Changed

- `rcashlar subtract` command:
  - Matching of channels is now based on excitation and emission wavelengths.
  - The command now supports scans with a different number of channels.

### Fixed

- Microscope model: Updated the hard-coded microscope model
- `rcashlar register` command: Fixed an issue where the command would not
  function properly when run on the same scan.
- Dependencies: Updated versions of `palom` to resolve memory leak issue and add
  `lxml` because it's not a dependency in `ome-types` >= 0.5.
  