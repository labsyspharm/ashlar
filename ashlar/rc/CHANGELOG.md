# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

---

## [afr-2025.5.1] - 2025-05-05

### Changed

- Applied pre-smoothing for angle registration to improve accuracy.
- Switched to less aggressive window masking (Kaiser with `beta=2`, was Hann) to
  enhance phase correlation accuracy for low signal-to-noise (S/N) images.

### Improved

- Print refined rotation angle during alignment for better debugging and
  analysis.
- Enhanced thumbnail contrast in QC plots for better visualization.

## [afr-2025.4.1] - 2025-04-03

### Fixed

- Restore the `--add_camera_bias_back` operation: this functionality was
  unintentionally broken in commit 8e4a0e.
- Skip camera bias correction for the fiducial channel: For non-fiducial
  channels, camera bias correction will always be applied, regardless of the
  presence of a BG channel.

## [afr-2025.2.1] - 2025-02-03

### Changed

- Update Bioformats JAR to v8.0.1
- Use GitHub hosted bioformats instead of from the UoD server

## [afr-2024.12.1] - 2024-12-27

### Added

- `rcashlar version` command: Print version then exit.

### Changed

- Update versioning to be compatible with semantic versioning
- `max_shift` in `stitch` and `register` command default to 30 µm (was 15)
- `rcashlar combine` command:
  - Replace `--dna_filename` with `--dna_file`
  - The DNA channel in the `--dna_file` will be the first channel of the
    combined image; other channels in that image will not be written to the
    combined image
  - DNA channel in other images will not be written to the combined image

## [afr-2024-9-2] - 2024-09-30

### Changed

- Microscope serial number metadata: use "serialNumber" field in rcjob file,
  fall back to "jobHostname" when not available for backward compatibility

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
  