from __future__ import division
import atexit
import bioformats
import javabridge
import numpy as np
import scipy.ndimage
import skimage.util
import skimage.feature
import skimage.filters
import matplotlib.pyplot as plt


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
        s0 = self.tile_size(0)
        if any(any(self.tile_size(i) != s0) for i in range(1, self.num_images)):
            raise ValueError("Image series must all have the same dimensions")
        self.size = s0

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
    def grid_dimensions(self):
        pos = self.positions
        shape = np.array([len(set(pos[:, d])) for d in range(2)])
        if np.prod(shape) != self.num_images:
            raise ValueError("Series positions do not form a grid")
        return shape


    @property
    def positions(self):
        return np.vstack([
            self.tile_position(i) for i in range(self.num_images)
        ])

    @property
    def centers(self):
        return self.positions + self.size / 2

    @property
    def origin(self):
        return self.positions.min(axis=0)


class Reader(object):

    def __init__(self, path):
        _init_bioformats()
        self.path = path
        self.metadata = Metadata(self.path)
        self.ir = bioformats.ImageReader(self.path)

    def read(self, series, c):
        return np.flipud(self.ir.read(c=c, series=series))


class EdgeAligner(object):

    def __init__(self, reader):
        self.reader = reader
        self.max_shift = 0.05
        self._cache = {}

    def register(self, t1, t2):
        #print '  %d -> %d' % (t1, t2),
        key = tuple(sorted((t1, t2)))
        try:
            shift, error = self._cache[key]
            #print '<cached>',
        except KeyError:
            # Register nearest-pixel image overlaps.
            img1, img2 = self.overlap(t1, t2)
            img1 = whiten(img1)
            img2 = whiten(img2)
            shift, error, _ = skimage.feature.register_translation(img1, img2,
                                                                   10)
            # Add fractional part of offset back in.
            offset1, offset2, _ = self.intersection(t1, t2)
            shift += np.modf(offset1 - offset2)[0]
            # Constrain shift.
            if any(np.abs(shift) > self.max_shift * self.reader.metadata.size):
                shift[:] = 0
                error = 1
            self._cache[key] = (shift, error)
        #print
        if t1 > t2:
            shift = -shift
        # Return copy of shift to prevent corruption of cached values.
        return shift.copy(), error

    def intersection(self, t1, t2):
        corners1 = self.reader.metadata.positions[[t1, t2]]
        corners2 = corners1 + self.reader.metadata.size
        return intersection(corners1, corners2)

    def crop(self, tile, offset, shape):
        img = self.reader.read(series=tile, c=0)
        return crop(img, offset, shape)

    def overlap(self, t1, t2):
        offset1, offset2, shape = self.intersection(t1, t2)
        img1 = self.crop(t1, offset1, shape)
        img2 = self.crop(t2, offset2, shape)
        return img1, img2

    def debug(self, t1, t2):
        shift, _ = self.register(t1, t2)
        o1, o2 = self.overlap(t1, t2)
        w1 = whiten(o1)
        w2 = whiten(o2)
        corr = np.fft.fftshift(np.abs(np.fft.ifft2(
            np.fft.fft2(w1) * np.fft.fft2(w2).conj()
        )))
        stack = np.vstack
        rows, cols = 3, 1
        if corr.shape[0] > corr.shape[1]:
            stack = np.hstack
            rows, cols = cols, rows
        plt.figure()
        plt.subplot(rows, cols, 1)
        plt.imshow(stack([o1, o2]))
        ax = plt.subplot(rows, cols, 2)
        ax.set_xticks([])
        ax.set_yticks([])
        plt.imshow(stack([w1, w2]))
        ax = plt.subplot(rows, cols, 3)
        ax.set_xticks([])
        ax.set_yticks([])
        plt.imshow(corr)
        origin = np.array(corr.shape) / 2
        plt.plot(origin[1], origin[0], 'r+')
        shift += origin
        plt.plot(shift[1], shift[0], 'rx')
        plt.tight_layout(0, 0, 0)


class LayerAligner(object):

    def __init__(self, reader, reference_reader, reference_positions,
                 reference_shifts):
        self.reader = reader
        self.reference_reader = reference_reader
        self.reference_positions = reference_positions
        self.reference_shifts = reference_shifts
        self.max_shift = 0.05
        self.positions = reader.metadata.positions - reader.metadata.origin
        dist = scipy.spatial.distance.cdist(self.reference_positions,
                                            self.positions)
        self.reference_idx = np.argmin(dist, 0)

    def register(self, t):
        ref_img, img = self.overlap(t)
        ref_img = whiten(ref_img)
        img = whiten(img)
        shift, error, _ = skimage.feature.register_translation(ref_img, img, 10)
        # Add fractional part of offset back in.
        offset1, offset2, _ = self.intersection(t)
        shift += np.modf(offset1 - offset2)[0]
        new_position = self.positions[t] + shift
        # Constrain shift.
        rel_shift = shift - self.reference_shifts[t]
        if any(np.abs(rel_shift) > self.max_shift * self.reader.metadata.size):
            print "\n%s > %s" % (np.abs(shift), self.max_shift * self.reader.metadata.size)
            new_position = self.positions[t]
            error = 1
        return new_position, error

    def ref_pos(self, t):
        return self.reference_positions[self.reference_idx[t]]

    def intersection(self, t):
        corners1 = np.vstack([self.ref_pos(t), self.positions[t]])
        corners2 = corners1 + self.reader.metadata.size
        offset1, offset2, shape = intersection(corners1, corners2)
        shape = shape // 32 * 32
        return offset1, offset2, shape

    def overlap(self, t):
        offset1, offset2, shape = self.intersection(t)
        ref_t = self.reference_idx[t]
        img1 = self.reference_reader.read(series=ref_t, c=0)
        img2 = self.reader.read(series=t, c=0)
        ov1 = crop(img1, offset1, shape)
        ov2 = crop(img2, offset2, shape)
        return ov1, ov2

    def debug(self, t):
        shift, _ = self.register(t)
        o1, o2 = self.overlap(t)
        w1 = whiten(o1)
        w2 = whiten(o2)
        corr = np.fft.fftshift(np.abs(np.fft.ifft2(
            np.fft.fft2(w1) * np.fft.fft2(w2).conj()
        )))
        plt.figure()
        plt.subplot(1, 3, 1)
        plt.imshow(np.vstack([o1, o2]))
        ax = plt.subplot(1, 3, 2)
        ax.set_xticks([])
        ax.set_yticks([])
        plt.imshow(np.vstack([w1, w2]))
        ax = plt.subplot(1, 3, 3)
        ax.set_xticks([])
        ax.set_yticks([])
        plt.imshow(corr)
        shift -= self.ref_pos(t)
        origin = np.array(corr.shape) / 2
        plt.plot(origin[1], origin[0], 'r+')
        shift += origin
        plt.plot(shift[1], shift[0], 'rx')
        plt.tight_layout(0, 0, 0)


def whiten(img):
    img = skimage.filters.laplace(img)
    # Other possible whitening functions:
    #img = skimage.filters.roberts(img)
    #img = skimage.filters.scharr(img)
    #img = skimage.filters.sobel(img)
    #img = np.log(img)
    #img = img - scipy.ndimage.filters.gaussian_filter(img, 2) + 0.5
    return img


def intersection(corners1, corners2):
    position = corners1.max(axis=0)
    shape = np.ceil(corners2.min(axis=0) - position).astype(int)
    if any(shape <= 0):
        raise ValueError("Tiles do not intersect")
    offset1, offset2 = corners1 - position
    return offset1, offset2, shape


def crop(img, offset, shape):
    # Note that this only crops to the nearest whole-pixel offset.
    start = -offset.astype(int)
    end = start + shape
    img = img[start[0]:end[0], start[1]:end[1]]
    return img


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
