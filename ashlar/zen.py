import pathlib
from urllib.parse import unquote as urllib_unquote
import xml.etree.ElementTree
import numpy as np
import skimage.io
from . import reg


class ZenMetadata(reg.Metadata):

    def __init__(self, path):
        self.path = pathlib.Path(path).resolve()
        self.base_path = self.path.parent
        self._init_metadata()

    def _init_metadata(self):
        tree = xml.etree.ElementTree.parse(str(self.path))
        positions_map = {}
        tile_size = None
        self._num_channels = 0
        self.image_paths = {}
        for image in tree.findall('Image'):
            path = urllib_unquote(image.findtext('Filename'))
            bounds = image.find('Bounds')
            series = int(bounds.attrib['StartM'])
            channel = int(bounds.attrib['StartC'])
            start_x = int(bounds.attrib['StartX'])
            start_y = int(bounds.attrib['StartY'])
            size_x = int(bounds.attrib['SizeX'])
            size_y = int(bounds.attrib['SizeY'])
            position = start_y, start_x
            size = size_y, size_x
            if channel == 0:
                positions_map[series] = position
            else:
                if positions_map[series] != position:
                    print(
                        "Position for series %d, channel %d doesn't match"
                        " channel 0." % (series, channel)
                    )
            if tile_size is None:
                tile_size = size
            else:
                if tile_size != size:
                    raise Exception(
                        "Size for series %d, channel %d doesn't match"
                        " first image." % (series, channel)
                    )
            self.image_paths[series, channel] = path
            self._num_channels = max(self._num_channels, channel + 1)
        positions = [pos for series, pos in sorted(positions_map.items())]
        self._positions = np.array(positions, np.float64)
        self._tile_size = np.array(tile_size, np.int64)

    @property
    def _num_images(self):
        return len(self._positions)

    @property
    def num_channels(self):
        return self._num_channels

    @property
    def pixel_size(self):
        return 1.0

    @property
    def pixel_dtype(self):
        return np.dtype(np.uint16)

    def tile_size(self, i):
        return self._tile_size

    def image_path(self, series, c):
        return self.base_path / self.image_paths[series, c]


class ZenReader(reg.Reader):

    def __init__(self, path):
        self.metadata = ZenMetadata(path)
        self.path = pathlib.Path(path)

    def read(self, series, c):
        path = self.metadata.image_path(series, c)
        img = skimage.io.imread(str(path), key=0)
        return img
