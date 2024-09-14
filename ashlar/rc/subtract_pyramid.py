import itertools
import pathlib

import numpy as np
import tifffile
import zarr
import joblib

from .. import _version, reg


class SubtractPyramid(reg.PyramidWriter):
    def __init__(
        self,
        bg_mosaic,
        ab_mosaic,
        path,
        as_float=False,
        fiducial_channel=0,
        bg_intensity_scaling_factor=None,
        camera_bias=0,
        add_camera_bias_back=False,
        peak_size=1024,
        verbose=False,
        subtraction_config=None,
    ):
        super().__init__(
            [bg_mosaic, ab_mosaic], path, peak_size=peak_size, verbose=verbose
        )
        self.bg_mosaic = bg_mosaic
        self.ab_mosaic = ab_mosaic
        assert self.bg_mosaic.dtype == self.ab_mosaic.dtype

        self.as_float = as_float
        self.fiducial_channel = fiducial_channel
        self.bg_intensity_scaling_factor = bg_intensity_scaling_factor
        self.camera_bias = camera_bias
        self.add_camera_bias_back = add_camera_bias_back

        if bg_intensity_scaling_factor is None:
            self.bg_intensity_scaling_factor = np.ones(self.num_channels)
        assert len(self.bg_intensity_scaling_factor) == self.num_channels

        if subtraction_config is not None:
            # FIXME re-assign bg_mosaic.channels and handle None
            self.bg_mosaic.channels = [
                cc["bg_channel"].get("channel_index", None) for cc in subtraction_config
            ]
            bg_int_factor = [
                np.divide(
                    cc["exposure_time"], cc["bg_channel"].get("exposure_time", np.nan)
                )
                for cc in subtraction_config
            ]
            self.bg_intensity_scaling_factor = np.nan_to_num(bg_int_factor, nan=0)

        if self.fiducial_channel is not None:
            if self.fiducial_channel in self.ab_mosaic.channels:
                self.bg_mosaic.channels[
                    self.ab_mosaic.channels.index(self.fiducial_channel)
                ] = None

    def assemble_all_channels(self):
        self.cache_path = f"{pathlib.Path(self.path)}.zarr"
        root = zarr.open(self.cache_path, mode="w")
        root.create_groups("ab_mosaic", "bg_mosaic")
        tasks = []
        for idx, (ch_ab, ch_bg) in enumerate(
            zip(self.ab_mosaic.channels, self.bg_mosaic.channels)
        ):
            root["ab_mosaic"].zeros(
                idx,
                shape=self.ab_mosaic.shape,
                dtype=self.dtype,
                chunks=(1024, 1024),
            )
            root["bg_mosaic"].zeros(
                idx,
                shape=self.bg_mosaic.shape,
                dtype=self.bg_mosaic.dtype,
                chunks=(1024, 1024),
            )
            tasks.append((self.ab_mosaic.assemble_channel, ch_ab, root["ab_mosaic"][idx]))
            if ch_bg is None:
                continue
            tasks.append((self.bg_mosaic.assemble_channel, ch_bg, root["bg_mosaic"][idx]))
        n_jobs = min(len(tasks), joblib.cpu_count())
        verboses = [True] + [False] * (len(tasks) - 1)
        print("Generating mosaics in parallel")
        _ = joblib.Parallel(n_jobs=n_jobs, verbose=0)(
            joblib.delayed(m_func)(channel, out=out_zarr, verbose=v)
            for (m_func, channel, out_zarr), v in zip(tasks, verboses)
        )
        self.mosaics_zarr = root

    def base_tiles(self):
        h, w = self.base_shape
        th, tw = self.tile_shapes[0]
        is_int_dtype = np.issubdtype(self.dtype, np.integer)
        if is_int_dtype:
            iinfo = np.iinfo(self.dtype)
            dmin, dmax = iinfo.min, iinfo.max
        for idx, (ch_ab, ch_bg, int_scaling_factor) in enumerate(
            zip(
                self.ab_mosaic.channels,
                self.bg_mosaic.channels,
                self.bg_intensity_scaling_factor,
            )
        ):
            print(f"Channel {ch_ab}:")
            if ch_bg is not None:
                print("    Ab & Bg Image")
            else:
                print("    Ab Image")

            bg_img = self.mosaics_zarr["bg_mosaic"][idx]
            ab_img = self.mosaics_zarr["ab_mosaic"][idx]
            for y, x in itertools.product(range(0, h, th), range(0, w, tw)):
                corrected_ab_tile = (
                    ab_img[y : y + th, x : x + tw].astype(np.float32) - self.camera_bias
                )
                if ch_bg is None:
                    corrected_bg_tile = np.zeros_like(corrected_ab_tile)
                else:
                    corrected_bg_tile = (
                        bg_img[y : y + th, x : x + tw].astype(np.float32)
                        - self.camera_bias
                    )
                subtracted = corrected_ab_tile - corrected_bg_tile * int_scaling_factor
                if self.as_float or (not is_int_dtype):
                    yield subtracted
                else:
                    yield np.clip(np.round(subtracted), dmin, dmax).astype(self.dtype)
            # Allow img to be freed immediately to avoid keeping it in
            # memory while the next loop iteration calls assemble_channel.
            ab_img = None
            bg_img = None

    @property
    def num_channels(self):
        return len(self.ab_mosaic.channels)

    @property
    def dtype(self):
        return self.ab_mosaic.dtype

    @property
    def channels(self):
        return self.ab_mosaic.channels

    def run(self):
        if not hasattr(self, "mosaics_zarr"):
            self.assemble_all_channels()
        dtype = np.float32 if self.as_float else self.dtype
        compression = dict(compression="adobe_deflate", predictor=True)
        if dtype == np.float32:
            compression = dict(compression=None)
        pixel_size = self.ref_mosaic.aligner.metadata.pixel_size
        resolution_cm = 10000 / pixel_size
        software = f"Ashlar v{_version.get_versions()['version']}"
        metadata = {
            "Creator": software,
            "Pixels": {
                "PhysicalSizeX": pixel_size,
                "PhysicalSizeXUnit": "\u00b5m",
                "PhysicalSizeY": pixel_size,
                "PhysicalSizeYUnit": "\u00b5m",
            },
        }
        with tifffile.TiffWriter(self.path, ome=True, bigtiff=True) as tiff:
            tiff.write(
                data=self.base_tiles(),
                metadata=metadata,
                software=software.encode("utf-8"),
                shape=self.level_full_shapes[0],
                subifds=int(self.num_levels - 1),
                dtype=dtype,
                tile=self.tile_shapes[0],
                resolution=(resolution_cm, resolution_cm, "centimeter"),
                # FIXME Propagate this from input files (especially RGB).
                photometric="minisblack",
                **compression,
            )
            if self.verbose:
                print("Generating pyramid")
            for level, (shape, tile_shape) in enumerate(
                zip(self.level_full_shapes[1:], self.tile_shapes[1:]), 1
            ):
                if self.verbose:
                    print(f"    Level {level} ({shape[2]} x {shape[1]})")
                tiff.write(
                    data=self.subres_tiles(level),
                    shape=shape,
                    subfiletype=1,
                    dtype=dtype,
                    tile=tile_shape,
                    **compression,
                )
                if self.verbose:
                    print()

    def cleanup(self):
        self.mosaics_zarr.store.rmdir()
