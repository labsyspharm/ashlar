import re
import itertools
import pathlib
import numpy as np
import skimage.io
from . import reg


# Classes for reading datasets consisting of TIFF files with a naming pattern.
# The pattern must include an integer series number, and optionally a channel
# name or number, and well name.
#
# This code is experimental and probably still has a lot of rough edges.


def format_to_regex(s):
    # Translate a restricted subset of the "format" pattern language to
    # a matching regex with named capture.
    s = s.replace('.', '\.')
    regex = re.sub(r'{([^:}]+):?([^}]*)}', f2r_repl, s)
    return regex

def f2r_repl(m):
    r = '(?P<' + m.group(1) + '>.'
    if re.match(r'^\d+$', m.group(2)):
        r += '{' + m.group(2) + '}'
    else:
        r += '*?'
    r += ')'
    return r


class FileSeriesMetadata(reg.PlateMetadata):

    def __init__(
        self, path, pattern, overlap, width, height, layout, direction,
        pixel_size,
    ):
        # The pattern argument uses the Python Format String syntax with
        # required "series" and optionally "channel" fields. A width
        # specification with leading zeros must be used for any fields that are
        # zero-padded. The pattern is used both to parse the filenames upon
        # initialization as well as to synthesize filenames when reading images.
        # Example pattern: 'img_s{series}_w{channel}.tif'
        super(FileSeriesMetadata, self).__init__()
        self.path = pathlib.Path(path)
        self.pattern = pattern
        self.overlap = overlap
        self.width = width
        self.height = height
        self.layout = layout
        self.direction = direction
        self._pixel_size = pixel_size
        self._enumerate_tiles()

    def _enumerate_tiles(self):
        regex = format_to_regex(self.pattern)
        wells = set()
        series = set()
        channels = set()
        n = 0
        self.filename_components = {}
        for p in self.path.iterdir():
            match = re.match(regex, p.name)
            if match:
                gd = match.groupdict()
                w = gd.get('well')
                s = int(gd['series'])
                c = gd.get('channel')
                wells.add(w)
                series.add(s)
                channels.add(c)
                self.filename_components[w, s, c] = gd
                n += 1
        if len(self.filename_components) != len(wells) * len(series) * len(channels):
            raise Exception("Missing images detected")
        # Build sorted list of (well, series) tuples for all wells.
        self.all_series = sorted(set(
            k[:2] for k in self.filename_components.keys()
        ))
        self.well_map = dict(enumerate(sorted(wells)))
        self._actual_num_images = len(series) * len(wells)
        self.channel_map = dict(enumerate(sorted(channels)))
        path = self.path / self.filename(0, 0)
        img = skimage.io.imread(str(path))
        self._tile_size = np.array(img.shape[:2])
        self._dtype = img.dtype
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

    # This only supports a single plate at the moment.
    @property
    def num_plates(self):
        return 1

    @property
    def num_wells(self):
        return [len(self.well_map)]

    @property
    def plate_well_series(self):
        return [
            list([a[0] for a in v] for k, v in itertools.groupby(enumerate(self.all_series), key=lambda x: x[1][0]))
        ]

    def plate_name(self, i):
        assert i == 0, "Plate index out of range"
        return "Plate_1"

    def well_name(self, plate, i):
        assert plate == 0, "Plate index out of range"
        return self.well_map[i]

    @property
    def pixel_size(self):
        return self._pixel_size

    @property
    def pixel_dtype(self):
        return self._dtype

    def tile_position(self, i):
        row, col = self.tile_rc(i)
        return [row, col] * self.tile_size(i) * (1 - self.overlap)

    def tile_size(self, i):
        return self._tile_size

    def tile_rc(self, i):
        if self.direction == "horizontal":
            row = i // self.width
            col = i % self.width
            if self.layout == "snake" and row % 2 == 1:
                col = self.width - 1 - col
        else:
            row = i % self.height
            col = i // self.height
            if self.layout == "snake" and col % 2 == 1:
                row = self.height - 1 - row
        return row, col

    def filename(self, series, c):
        well, series = self.all_series[self.active_series[series]]
        c = self.channel_map[c]
        components = self.filename_components[well, series, c]
        return self.pattern.format(**components)


class FileSeriesReader(reg.PlateReader):

    def __init__(
        self, path, pattern, overlap, width, height, layout="raster",
        direction="horizontal", pixel_size=1.0, plate=None, well=None
    ):
        # See FileSeriesMetadata for an explanation of the pattern syntax.
        if layout not in ("raster", "snake"):
            raise ValueError("layout must be 'raster' or 'snake'")
        if direction not in ("horizontal", "vertical"):
            raise ValueError("direction must be 'horizontal' or 'vertical'")
        self.path = pathlib.Path(path)
        self.pattern = pattern
        self.metadata = FileSeriesMetadata(
            self.path, self.pattern, overlap, width, height, layout, direction,
            pixel_size
        )
        self.metadata.set_active_plate_well(plate, well)

    def read(self, series, c):
        # TODO: Address tension between non-plate and plate-aware modes
        # here and in Metadata class.
        path = str(self.path / self.metadata.filename(series, c))
        kwargs = {}
        if self.metadata.multi_channel_tiles:
            kwargs['key'] = c
        else:
            # In case of multi-plane images, only take the first plane. The
            # processing code only handles 2D image arrays!
            kwargs['key'] = 0
        return skimage.io.imread(path, **kwargs)
