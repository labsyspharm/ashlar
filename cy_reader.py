import reg
import re
import collections
import pathlib2 as pathlib
import numpy as np
import skimage.io
import skimage.external.tifffile


ImageMetadata = collections.namedtuple('ImageMetadata',
                                       'path position size')


class Metadata(reg.Metadata):

    def __init__(self, path):
        path = pathlib.Path(path).resolve()
        f = open(str(path.joinpath('TileConfiguration.txt')))
        self.images = []
        for line in f:
            if line[0] == '#':
                continue
            data = re.split(' *; *', line.strip())
            if len(data) != 3:
                continue
            filename, dummy, position = data
            image_path = str(path.joinpath(filename))
            position_str = re.split(r' *, *', position.strip('()'))
            position = np.array([float(p) for p in position_str])[::-1]
            with skimage.external.tifffile.TiffFile(image_path) as tf:
                size = np.array(tf.pages[0].shape)
            im = ImageMetadata(path=image_path, position=position, size=size)
            self.images.append(im)

    @property
    def num_images(self):
        return len(self.images)

    @property
    def num_channels(self):
        return 1

    @property
    def pixel_size(self):
        return np.array([1.0, 1.0])

    def tile_position(self, i):
        return self.images[i].position

    def tile_size(self, i):
        return self.images[i].size


class Reader(reg.Reader):

    def __init__(self, path):
        self.path = path
        self.metadata = Metadata(self.path)

    def read(self, series, c):
        assert c == 0
        return skimage.io.imread(self.metadata.images[series].path)


