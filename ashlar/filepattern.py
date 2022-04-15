import re
import pathlib
import numpy as np
import skimage.io
from . import reg


# Classes for reading datasets consisting of TIFF files with a naming pattern.
# The pattern must include row and column as integers, and optionally a channel
# name or number.
#
# This code is experimental and probably still has a lot of rough edges.


class FilePatternMetadata(reg.Metadata):

    def __init__(self, path, pattern, overlap, pixel_size):
        # The pattern argument uses the Python Format String syntax with
        # required "row" and "col" and optionally "channel" fields. A width
        # specification with leading zeros must be used for any fields that are
        # zero-padded. The pattern is used both to parse the filenames upon
        # initialization as well as to synthesize filenames when reading images.
        # Example pattern: 'img_{channel}_r{row:03}_c{col:03}.tif'
        self.path = pathlib.Path(path)
        self.pattern = pattern
        self.overlap = overlap
        self._pixel_size = pixel_size
        self._enumerate_tiles()

    def _enumerate_tiles(self):
        # Translate a restricted subset of the "format" pattern language to
        # a matching regex with named capture.
        regex = re.sub(r'{([^:}]+)(?:[^}]*)}', r'(?P<\1>.*?)',
                       self.pattern.replace('.', '\.'))
        rows = set()
        cols = set()
        channels = set()
        n = 0
        for p in self.path.iterdir():
            match = re.match(regex, p.name)
            if match:
                gd = match.groupdict()
                rows.add(int(gd['row']))
                cols.add(int(gd['col']))
                channels.add(gd.get('channel'))
                n += 1
        if n != len(rows) * len(cols) * len(channels):
            raise Exception("Tiles do not form a full rectangular grid")
        self._actual_num_images = len(rows) * len(cols)
        self.channel_map = dict(enumerate(sorted(channels)))
        self.height = len(rows)
        self.width = len(cols)
        self.row_offset = min(rows)
        self.col_offset = min(cols)
        path = self.path / self.pattern.format(
            row=self.row_offset, col=self.col_offset,
            channel=self.channel_map[0]
        )
        img = skimage.io.imread(str(path))
        self._tile_size = np.array(img.shape[:2])
        self.multi_channel_tiles = False
        # Handle multi-channel tiles (pattern must not include channel).
        if len(self.channel_map) == 1 and img.ndim == 3:
            self.channel_map = {c: None for c in range(img.shape[2])}
            self.multi_channel_tiles = True
        self._num_channels = len(self.channel_map)

    @property
    def _num_images(self):
        return self._actual_num_images

    @property
    def num_channels(self):
        return self._num_channels

    @property
    def pixel_size(self):
        return self._pixel_size

    @property
    def pixel_dtype(self):
        return np.uint16

    def tile_position(self, i):
        row, col = self.tile_rc(i)
        return [row, col] * self.tile_size(i) * (1 - self.overlap)

    def tile_size(self, i):
        return self._tile_size

    def tile_rc(self, i):
        row = i // self.width + self.row_offset
        col = i % self.width + self.col_offset
        return row, col


class FilePatternReader(reg.Reader):

    def __init__(self, path, pattern, overlap, pixel_size=1.0):
        # See FilePatternMetadata for an explanation of the pattern syntax.
        self.path = pathlib.Path(path)
        self.pattern = pattern
        self.overlap = overlap
        self.metadata = FilePatternMetadata(
            self.path, self.pattern, overlap, pixel_size
        )

    def read(self, series, c):
        path = str(self.path / self.filename(series, c))
        kwargs = {}
        if self.metadata.multi_channel_tiles:
            kwargs['key'] = c
        else:
            # In case of multi-plane images, only take the first plane. The
            # processing code only handles 2D image arrays!
            kwargs['key'] = 0
        return skimage.io.imread(path, **kwargs)

    def filename(self, series, c):
        row, col = self.metadata.tile_rc(series)
        c = self.metadata.channel_map[c]
        return self.pattern.format(row=row, col=col, channel=c)
