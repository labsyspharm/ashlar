from __future__ import division
import atexit
import bioformats
import javabridge
import numpy as np
import scipy.ndimage
import skimage.util
import skimage.feature
import skimage.filters


def _init_bioformats():
    if javabridge._javabridge.get_vm().is_active():
        return
    javabridge.start_vm(class_path=bioformats.JARS)
    bioformats.init_logger()
    # Hack module to fix py3 assumptions which break XML parsing.
    bioformats.omexml.str = unicode

def _deinit_bioformats():
    javabridge.kill_vm()


class Metadata(object):

    def __init__(self, path):
        _init_bioformats()
        ome_xml = bioformats.get_omexml_metadata(path)
        self._metadata = bioformats.OMEXML(ome_xml)

    @property
    def num_images(self):
        return self._metadata.image_count

    @property
    def pixel_size(self):
        px_node = self._metadata.image(0).Pixels.node
        return np.array([
            float(px_node.get('PhysicalSize%s' % d)) for d in 'Y', 'X'
        ])

    def tile_position(self, i):
        plane = self._metadata.image(i).Pixels.Plane(0)
        position_microns = np.array([plane.PositionY, plane.PositionX])
        position_pixels = position_microns / self.pixel_size
        return position_pixels

    def tile_size(self, i):
        pixels = self._metadata.image(i).Pixels
        return np.array([pixels.SizeY, pixels.SizeX])

    @property
    def positions(self):
        return np.vstack([
            self.tile_position(i) for i in range(self.num_images)
        ])

    @property
    def centers(self):
        return self.positions + self.sizes / 2

    @property
    def origin(self):
        return self.positions.min(axis=0)

    @property
    def sizes(self):
        return np.vstack([
            self.tile_size(i) for i in range(self.num_images)
        ])


class Reader(object):

    def __init__(self, path):
        _init_bioformats()
        self.path = path
        self.metadata = Metadata(self.path)
        self.ir = bioformats.ImageReader(self.path)

    def read(self, series, c):
        return self.ir.read(c=c, series=series)


class Aligner(object):

    def __init__(self, reader):
        self.reader = reader
        self._cache = {}

    def register(self, t1, t2):
        print '  %d -> %d' % (t1, t2),
        key = tuple(sorted((t1, t2)))
        try:
            shift, error = self._cache[key]
            print '<cached>',
        except KeyError:
            im1 = skimage.filters.laplace(self.reader.read(series=t1, c=0))
            im2 = skimage.filters.laplace(self.reader.read(series=t2, c=0))
            shift, error, _ = skimage.feature.register_translation(im1, im2, 10)
            self._cache[key] = (shift, error)
        print
        return shift, error


def paste(target, img, pos, debug=False):
    """Composite img into target using maximum intensity projection.

    target: uint
    img: float

    """
    pos_f, pos_i = np.modf(pos)
    yi, xi = pos_i.astype('i8')
    # Clip img to the edges of the mosaic.
    if yi < 0:
        img = img[-yi:]
        yi = 0
    if xi < 0:
        img = img[:, -xi:]
        xi = 0
    # This is a bit wrong on the edges in the subpixel shift direction. The
    # fractional pixels that would be shifted off the edges of the image are
    # actually discarded. However since the images being tiled in this
    # application have far more overlap than a single pixel, it's irrelevant.
    target_slice = target[yi:yi+img.shape[0], xi:xi+img.shape[1]]
    img = crop_like(img, target_slice)
    img = scipy.ndimage.shift(img, pos_f)
    np.clip(img, 0, 1, img)
    img = skimage.util.dtype.convert(img, target.dtype)
    target_slice[:, :] = np.maximum(target_slice, img)
    if debug:
        # Render a faint outline of the pasted image.
        # TODO 6000 is arbitrary and should be calculated from the data.
        # Also these lines can cause registration problems, so ideally
        # this step should be performed on the final images by using the
        # accumulated list of per-tile offsets.
        target_slice[0, :] += 6000
        target_slice[-1, :] += 6000
        target_slice[1:-1, 0] += 6000
        target_slice[1:-1, -1] += 6000
        np.clip(target_slice[:, :], 0, np.iinfo(target.dtype).max)


def crop_like(img, target):
    if (img.shape[0] > target.shape[0]):
        img = img[:target.shape[0], :]
    if (img.shape[1] > target.shape[1]):
        img = img[:, :target.shape[1]]
    return img
