import pathlib

import numpy as np
import tifffile
import zarr
from joblib import Parallel, delayed

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
        peak_size=1024,
        verbose=False
    ):
        super().__init__(
            [bg_mosaic, ab_mosaic], path, peak_size=peak_size, verbose=verbose
        )
        self.bg_mosaic = bg_mosaic
        self.ab_mosaic = ab_mosaic
        assert self.bg_mosaic.channels == self.ab_mosaic.channels
        assert self.bg_mosaic.dtype == self.ab_mosaic.dtype

        self.as_float = as_float
        self.fiducial_channel = fiducial_channel
        self.bg_intensity_scaling_factor = bg_intensity_scaling_factor
        
        if bg_intensity_scaling_factor is None:
            self.bg_intensity_scaling_factor = np.ones(self.num_channels)
        assert len(self.bg_intensity_scaling_factor) == len(self.num_channels)
    
    def assemble_all_channels(self):
        self.cache_path = f"{pathlib.Path(self.path)}.zarr"
        root = zarr.open(self.cache_path, mode='w')
        root.create_groups('ab_mosaic', 'bg_mosaic')
        to_assemble = []
        for channel in self.channels:
            root['ab_mosaic'].zeros(
                channel, 
                shape=self.ab_mosaic.shape,
                dtype=self.dtype,
                chunks=(1024, 1024)
            )
            root['bg_mosaic'].zeros(
                channel,
                shape=self.bg_mosaic.shape,
                dtype=self.bg_mosaic.dtype,
                chunks=(1024, 1024)
            )
            to_assemble.append((self.ab_mosaic.assemble_channel, channel, root['ab_mosaic'][channel]))
            if channel != self.fiducial_channel:
                to_assemble.append((self.bg_mosaic.assemble_channel, channel, root['bg_mosaic'][channel]))

        Parallel(n_jobs=len(to_assemble), verbose=1)(
            delayed(m_func)(channel, out=out_zarr)
            for m_func, channel, out_zarr in to_assemble
        )
        self.mosaics_zarr = root

    def base_tiles(self):
        h, w = self.base_shape
        th, tw = self.tile_shapes[0]
        is_int_dtype = np.issubdtype(self.dtype, np.integer)
        if is_int_dtype:
            iinfo = np.iinfo(self.dtype)
            dmin, dmax = iinfo.min, iinfo.max
        for channel, int_scaling_factor in zip(
            self.channels, self.bg_intensity_scaling_factor
        ):
            print(f"Channel {channel}:")

            if channel != self.fiducial_channel:
                print("    Ab & Bg Image")
                bg_img = self.mosaics_zarr['bg_mosaic'][channel]
            else:
                print("    Ab Image")
                bg_img = np.zeros(self.ab_mosaic.shape, dtype=self.dtype)
            ab_img = self.mosaics_zarr['ab_mosaic'][channel]
            for y in range(0, h, th):
                for x in range(0, w, tw):
                    # Returning a copy makes the array contiguous, avoiding
                    # a severely unoptimized code path in ndarray.tofile.
                    subtracted = (
                        ab_img[y:y+th, x:x+tw].astype(np.float32) - 
                        bg_img[y:y+th, x:x+tw].astype(np.float32) * int_scaling_factor
                    )
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
        return len(self.bg_mosaic.channels)
    
    @property
    def dtype(self):
        return self.bg_mosaic.dtype

    @property
    def channels(self):
        return self.bg_mosaic.channels
    
    def run(self):
        if not hasattr(self, 'mosaics_zarr'):
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
                "PhysicalSizeX": pixel_size, "PhysicalSizeXUnit": "\u00b5m",
                "PhysicalSizeY": pixel_size, "PhysicalSizeYUnit": "\u00b5m"
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
                **compression
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
                    **compression
                )
                if self.verbose:
                    print()
    
    def cleanup(self):
        self.mosaics_zarr.store.rmdir()