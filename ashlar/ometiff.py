import threading
import logging
import tifffile
from . import reg


class _SubstringFilter(logging.Filter):
    def __init__(self, *substrings):
        self.substrings = substrings

    def filter(self, record):
        return not any(s in record.getMessage() for s in self.substrings)

logging.getLogger('tifffile').addFilter(
    _SubstringFilter('OME series contains invalid TiffData index')
)


class OmeTiffReader(reg.BioformatsReader):
    """Uses BioFormats for metadata and tifffile for pixel data.

    Compared to BioformatsReader, this can be significantly faster for
    compressed OME-TIFF files because tifffile uses native decompression
    libraries. Thread-local TiffFile handles allow parallel reads without
    locking.
    """

    def __init__(self, path, plate=None, well=None):
        super().__init__(path, plate, well)
        self._thread_local = threading.local()

    def __getstate__(self):
        state = super().__getstate__()
        del state['_thread_local']
        return state

    def __setstate__(self, state):
        super().__setstate__(state)
        self._thread_local = threading.local()

    @property
    def _tif(self):
        try:
            return self._thread_local.tif
        except AttributeError:
            self._thread_local.tif = tifffile.TiffFile(self.path)
            return self._thread_local.tif

    def read(self, series, c):
        actual_series = self.metadata.active_series[series]
        tif_series = self._tif.series
        num_images = self.metadata.num_images

        if len(tif_series) >= num_images:
            # One tifffile series per BioFormats image; index by channel only.
            series_obj = tif_series[actual_series]
            axes = series_obj.axes
            page_axes = [a for a in axes if a not in ('Y', 'X')]
            coords = {'T': 0, 'C': c, 'Z': 0}
            page_idx = 0
            stride = 1
            for ax in reversed(page_axes):
                page_idx += coords.get(ax, 0) * stride
                stride *= series_obj.shape[axes.index(ax)]
        else:
            # All tiles and channels are flattened into a single tifffile
            # series (e.g. axes='IYX'). Assumes tile-major, channel-minor
            # ordering (OME DimensionOrder "XYZCT"), which is the most common.
            series_obj = tif_series[0]
            n_pages_per_series = len(series_obj.pages) // num_images
            page_idx = actual_series * n_pages_per_series + c

        return series_obj.pages[page_idx].asarray()
