from __future__ import division
import sys
import warnings
import bioformats
import javabridge
import numpy as np
import scipy.ndimage
import skimage.util
import skimage.feature
import skimage.filters
import skimage.restoration.uft
import skimage.io
import sklearn.linear_model
import pyfftw
import networkx as nx
import queue
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import modest_image

# Patch np.fft to use pyfftw so skimage utilities can benefit.
np.fft = pyfftw.interfaces.numpy_fft


def _init_bioformats():
    if javabridge._javabridge.get_vm().is_active():
        return
    javabridge.start_vm(class_path=bioformats.JARS)
    bioformats.init_logger()
    # Hack module to fix py3 assumptions which break XML parsing.
    bioformats.omexml.str = unicode

def _deinit_bioformats():
    javabridge.kill_vm()


# TODO:
# - Write tables with summary information about alignments.


class Metadata(object):

    @property
    def num_images(self):
        raise NotImplementedError

    @property
    def num_channels(self):
        raise NotImplementedError

    @property
    def pixel_size(self):
        raise NotImplementedError

    def tile_position(self, i):
        raise NotImplementedError

    def tile_size(self, i):
        raise NotImplementedError

    @property
    def grid_dimensions(self):
        pos = self.positions
        shape = np.array([len(set(pos[:, d])) for d in range(2)])
        if np.prod(shape) != self.num_images:
            raise ValueError("Series positions do not form a grid")
        return shape

    @property
    def positions(self):
        if not hasattr(self, '_positions'):
            self._positions = np.vstack([
                self.tile_position(i) for i in range(self.num_images)
            ])
        return self._positions

    @property
    def size(self):
        if not hasattr(self, '_size'):
            s0 = self.tile_size(0)
            image_ids = range(1, self.num_images)
            if any(any(self.tile_size(i) != s0) for i in image_ids):
                raise ValueError("Image series must all have the same dimensions")
            self._size = s0
        return self._size

    @property
    def centers(self):
        return self.positions + self.size / 2

    @property
    def origin(self):
        return self.positions.min(axis=0)


class Reader(object):

    def read(self, series, c):
        raise NotImplementedError


class BioformatsMetadata(Metadata):

    def __init__(self, path):
        _init_bioformats()
        ome_xml = bioformats.get_omexml_metadata(path)
        self._metadata = bioformats.OMEXML(ome_xml)

    @property
    def num_images(self):
        return self._metadata.image_count

    @property
    def num_channels(self):
        return self._metadata.image(0).Pixels.channel_count

    @property
    def pixel_size(self):
        px_node = self._metadata.image(0).Pixels.node
        return np.array([
            float(px_node.get('PhysicalSize%s' % d)) for d in 'Y', 'X'
        ])

    def tile_position(self, i):
        plane = self._metadata.image(i).Pixels.Plane(0)
        # Invert Y so that stage position coordinates and image pixel
        # coordinates are aligned.
        position_microns = np.array([-plane.PositionY, plane.PositionX])
        position_pixels = position_microns / self.pixel_size
        return position_pixels

    def tile_size(self, i):
        pixels = self._metadata.image(i).Pixels
        return np.array([pixels.SizeY, pixels.SizeX])


class BioformatsReader(Reader):

    def __init__(self, path):
        _init_bioformats()
        self.path = path
        self.metadata = BioformatsMetadata(self.path)
        self.ir = bioformats.ImageReader(self.path)

    def read(self, series, c):
        return self.ir.read(c=c, series=series, rescale=False)


# TileStatistics = collections.namedtuple(
#     'TileStatistics',
#     'scan tile x_original y_original x y shift_x shift_y error'
# )


@property
def neighbors_graph(aligner):
    """Return graph of neighboring (overlapping) tiles.

    Tiles are considered neighbors if the 'city block' distance between them
    is less than the largest tile dimension.

    """
    # FIXME: This should properly test for overlap, possibly via
    # intersection of bounding rectangles.
    if not hasattr(aligner, '_neighbors_graph'):
        pdist = scipy.spatial.distance.pdist(aligner.metadata.positions,
                                             metric='cityblock')
        sp = scipy.spatial.distance.squareform(pdist)
        max_distance = aligner.metadata.size.max()
        edges = zip(*np.nonzero((sp > 0) & (sp < max_distance)))
        graph = nx.from_edgelist(edges)
        aligner._neighbors_graph = graph
    return aligner._neighbors_graph


class EdgeAligner(object):

    def __init__(self, reader, verbose=False):
        self.reader = reader
        self.verbose = verbose
        self.max_shift = 0.05
        self._cache = {}

    neighbors_graph = neighbors_graph

    def run(self):
        self.register_all()
        self.build_spanning_tree()
        self.calculate_positions()
        self.fit_model()

    def register_all(self):
        n = self.neighbors_graph.size()
        for i, (t1, t2) in enumerate(self.neighbors_graph.edges_iter(), 1):
            if self.verbose:
                sys.stdout.write('\r    aligning edge %d/%d' % (i, n))
                sys.stdout.flush()
            self.register_pair(t1, t2)
        if self.verbose:
            print

    def build_spanning_tree(self):
        line_graph = nx.line_graph(self.neighbors_graph)
        spanning_tree = nx.Graph()
        fringe = queue.PriorityQueue()
        start_edge = self.best_edge
        fringe.put((self.register_pair(*start_edge)[1], start_edge))
        while not fringe.empty():
            _, edge = fringe.get()
            if edge[0] in spanning_tree and edge[1] in spanning_tree:
                continue
            spanning_tree.add_edge(*edge)
            for next_edge in set(line_graph.neighbors(edge)):
                fringe.put((self.register_pair(*next_edge)[1], next_edge))
        self.spanning_tree = spanning_tree

    def calculate_positions(self):
        # Use the source node of the edge with the best alignment quality as the
        # reference tile against which all others will be aligned.
        reference_node = self.best_edge[0]
        shifts = {reference_node: np.array([0, 0])}
        for edge in nx.traversal.dfs_edges(self.spanning_tree, reference_node):
            source, dest = edge
            if source not in shifts:
                source, dest = dest, source
            shifts[dest] = shifts[source] + self.register_pair(source, dest)[0]
        self.shifts = np.array(zip(*sorted(shifts.items()))[1])
        self.positions = self.metadata.positions + self.shifts

    def fit_model(self):
        # Build list of connected components from spanning tree with
        # infinite-error edges discarded.
        forest = self.spanning_tree.copy()
        forest.remove_edges_from(
            e for e, v in self._cache.items() if v[1] == np.inf
        )
        components = sorted((list(c) for c in nx.connected_components(forest)),
                            key=len, reverse=True)
        # Fit LR model on positions of largest connected component.
        cc0 = components[0]
        self.lr = sklearn.linear_model.LinearRegression()
        self.lr.fit(self.metadata.positions[cc0], self.positions[cc0])
        # Adjust position of remaining components so their centroids match
        # the predictions of the model.
        for nodes in components[1:]:
            centroid_m = np.mean(self.metadata.positions[nodes], axis=0)
            centroid_f = np.mean(self.positions[nodes], axis=0)
            shift = self.lr.predict([centroid_m])[0] - centroid_f
            self.positions[nodes] += shift
        # Adjust positions and model intercept to put origin at 0,0.
        self.origin = self.positions.min(axis=0)
        self.positions -= self.origin
        self.lr.intercept_ -= self.origin
        self.centers = self.positions + self.metadata.size / 2


    def register_pair(self, t1, t2):
        """Return relative shift between images and the alignment error."""
        key = tuple(sorted((t1, t2)))
        try:
            shift, error = self._cache[key]
        except KeyError:
            # Register nearest-pixel image overlaps.
            img1, img2 = self.overlap(t1, t2)
            img1_f = fft2(whiten(img1))
            img2_f = fft2(whiten(img2))
            shift, error, _ = skimage.feature.register_translation(
                img1_f, img2_f, 10, 'fourier'
            )
            # Add fractional part of offset back in.
            offset1, offset2, _ = self.intersection(t1, t2)
            shift += np.modf(offset1 - offset2)[0]
            # Constrain shift.
            if any(np.abs(shift) > self.max_shift * self.metadata.size):
                shift[:] = 0
                error = np.inf
            self._cache[key] = (shift, error)
        if t1 > t2:
            shift = -shift
        # Return copy of shift to prevent corruption of cached values.
        return shift.copy(), error

    def intersection(self, t1, t2):
        corners1 = self.metadata.positions[[t1, t2]]
        corners2 = corners1 + self.metadata.size
        return intersection(corners1, corners2)

    def crop(self, tile, offset, shape):
        img = self.reader.read(series=tile, c=0)
        return crop(img, offset, shape)

    def overlap(self, t1, t2):
        offset1, offset2, shape = self.intersection(t1, t2)
        img1 = self.crop(t1, offset1, shape)
        img2 = self.crop(t2, offset2, shape)
        return img1, img2

    @property
    def best_edge(self):
        ordered_keys = sorted(self._cache, key=lambda k: self._cache[k][1])
        return ordered_keys[0]

    @property
    def metadata(self):
        return self.reader.metadata

    @property
    def mosaic_shape(self):
        upper_corners = self.positions + self.metadata.size
        max_dimensions = upper_corners.max(axis=0)
        return np.ceil(max_dimensions).astype(int)

    def debug(self, t1, t2):
        shift, _ = self.register_pair(t1, t2)
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
        plt.imshow(stack([w1, w2]).real)
        ax = plt.subplot(rows, cols, 3)
        ax.set_xticks([])
        ax.set_yticks([])
        plt.imshow(corr)
        origin = np.array(corr.shape) / 2
        plt.plot(origin[1], origin[0], 'r+')
        shift += origin
        plt.plot(shift[1], shift[0], 'rx')
        plt.colorbar()
        plt.tight_layout(0, 0, 0)


class LayerAligner(object):

    def __init__(self, reader, reference_aligner, verbose=False):
        self.reader = reader
        self.reference_aligner = reference_aligner
        self.verbose = verbose
        self.max_shift = 0.05
        self.tile_positions = self.metadata.positions - reference_aligner.origin
        reference_positions = (reference_aligner.metadata.positions
                               - reference_aligner.origin)
        dist = scipy.spatial.distance.cdist(reference_positions,
                                            self.tile_positions)
        self.reference_idx = np.argmin(dist, 0)
        self.reference_positions = reference_positions[self.reference_idx]

    neighbors_graph = neighbors_graph

    def run(self):
        self.register_all()
        self.calculate_positions()

    def register_all(self):
        n = self.metadata.num_images
        self.shifts = np.empty((n, 2))
        self.errors = np.empty(n)
        for i in range(n):
            if self.verbose:
                sys.stdout.write("\r    aligning tile %d/%d" % (i + 1, n))
                sys.stdout.flush()
            shift, error = self.register(i)
            self.shifts[i] = shift
            self.errors[i] = error
        if self.verbose:
            print

    def calculate_positions(self):
        self.positions = self.reference_aligner.positions + self.shifts
        self.constrain_positions()
        self.centers = self.positions + self.metadata.size / 2

    def constrain_positions(self):
        # Computed shifts of exactly 0,0 seem to result from failed
        # registration. We need to throw those out for this purpose.
        discard = (self.shifts == 0).all(axis=1)
        # Take the mean of registered shifts to determine the offset
        # (translation) from the reference image to this one.
        offset = np.nan_to_num(np.mean(self.shifts[~discard], axis=0))
        # Here we assume the fitted linear model from the reference image is
        # still appropriate, apart from the extra offset we just computed.
        predictions = self.reference_aligner.lr.predict(self.metadata.positions)
        # Discard any tile registration that's too far from the linear model,
        # replacing it with the relevant model prediction.
        distance = np.linalg.norm(self.positions - predictions - offset, axis=1)
        max_dist = self.max_shift * max(self.metadata.size)
        extremes = distance > max_dist
        # Recalculate the mean shift, also ignoring the extreme values.
        discard |= extremes
        self.offset = np.nan_to_num(np.mean(self.shifts[~discard], axis=0))
        # Fill in discarded shifts from the predictions.
        self.positions[discard] = predictions[discard] + self.offset

    def register(self, t):
        """Return relative shift between images and the alignment error."""
        ref_img, img = self.overlap(t)
        ref_img_f = fft2(whiten(ref_img))
        img_f = fft2(whiten(img))
        shift, error, _ = skimage.feature.register_translation(
            ref_img_f, img_f, 10, 'fourier'
        )
        # Add reported difference in stage positions.
        offset1, _, _ = self.intersection(t)
        shift -= offset1
        return shift, error

    def intersection(self, t):
        corners1 = np.vstack([self.reference_positions[t],
                              self.tile_positions[t]])
        corners2 = corners1 + self.reader.metadata.size
        offset1, offset2, shape = intersection(corners1, corners2)
        shape = shape // 32 * 32
        return offset1, offset2, shape

    def overlap(self, t):
        offset1, offset2, shape = self.intersection(t)
        ref_t = self.reference_idx[t]
        img1 = self.reference_aligner.reader.read(series=ref_t, c=0)
        img2 = self.reader.read(series=t, c=0)
        ov1 = crop(img1, offset1, shape)
        ov2 = crop(img2, offset2, shape)
        return ov1, ov2

    @property
    def metadata(self):
        return self.reader.metadata

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
        plt.imshow(np.vstack([w1, w2]).real)
        ax = plt.subplot(1, 3, 3)
        ax.set_xticks([])
        ax.set_yticks([])
        plt.imshow(corr)
        origin = np.array(corr.shape) / 2
        plt.plot(origin[1], origin[0], 'r+')
        shift += origin
        plt.plot(shift[1], shift[0], 'rx')
        plt.tight_layout(0, 0, 0)


class Mosaic(object):

    def __init__(self, aligner, shape, filename_format, channels=None,
                 verbose=False):
        self.aligner = aligner
        self.shape = shape
        self.filename_format = filename_format
        self.channels = self._sanitize_channels(channels)
        self.verbose = verbose
        self.filenames = []

    def _sanitize_channels(self, channels):
        all_channels = range(self.aligner.metadata.num_channels)
        if channels is None:
            channels = all_channels
        invalid_channels = sorted(set(channels) - set(all_channels))
        if invalid_channels:
            raise ValueError("invalid channels: %s" % invalid_channels)
        return channels

    def run(self):
        num_tiles = len(self.aligner.positions)
        for channel in self.channels:
            if self.verbose:
                print '    Channel %d:' % channel
            mosaic_image = np.zeros(self.shape, dtype=np.uint16)
            for tile, position in enumerate(self.aligner.positions):
                if self.verbose:
                    sys.stdout.write('\r        merging tile %d/%d'
                                     % (tile + 1, num_tiles))
                    sys.stdout.flush()
                tile_image = self.aligner.reader.read(c=channel, series=tile)
                paste(mosaic_image, tile_image, position)
            if self.verbose:
                print
            filename = self.filename_format % {'channel': channel}
            self.filenames.append(filename)
            if self.verbose:
                print "        writing %s" % filename
            with warnings.catch_warnings():
                warnings.filterwarnings(
                    'ignore', r'.* is a low contrast image', UserWarning,
                    '^skimage\.io'
                )
                skimage.io.imsave(filename, mosaic_image)


def fft2(img):
    return pyfftw.builders.fft2(img, planner_effort='FFTW_ESTIMATE',
                                avoid_copy=True, auto_align_input=True,
                                auto_contiguous=True)()


# Pre-calculate the Laplacian operator kernel. We'll always be using 2D images.
_laplace_kernel = skimage.restoration.uft.laplacian(2, (3, 3))[1]

def whiten(img):
    # Copied from skimage.filters.edges, with explicit aligned output from
    # convolve. Also the mask option was dropped.
    img = skimage.img_as_float(img)
    output = pyfftw.empty_aligned(img.shape, 'complex64')
    output.imag[:] = 0
    scipy.ndimage.convolve(img, _laplace_kernel, output.real)
    return output

    # Other possible whitening functions:
    #img = skimage.filters.roberts(img)
    #img = skimage.filters.scharr(img)
    #img = skimage.filters.sobel(img)
    #img = np.log(img)
    #img = img - scipy.ndimage.filters.gaussian_filter(img, 2) + 0.5


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


# TODO:
# - Deal with ringing from high-frequency elements. The wrapped edges of the
#   image are especially bad, where the wrapping introduces sharp
#   discontinuities. The edge artifacts could be dealt with in several ways
#   (extend the trailing image edge via mirroring, throw away some of the
#   trailing edge of the shifted result) but edges in the "true" image content
#   would require proper pre-filtering. What filter to use, and how to apply it
#   quickly?
# - Can we use real FFT for a ~50% overall speedup? Fourier-space matrices will
#   all be half-size in the last dimension, so FFT is around 50% faster and our
#   fshift calculations will be too.
# - Trailing edge pixels should be zeroed to match the behavior of
#   scipy.ndimage.shift, which we rely on in our maximum-intensity projection.
def fourier_shift(img, shift):
    # Ensure properly aligned complex64 data (fft requires complex to avoid
    # reallocation and copying).
    img = skimage.util.dtype.convert(img, dtype=np.float32)
    img = pyfftw.byte_align(img, dtype=np.complex64)
    # Compute per-axis frequency values according to the Fourier shift theorem.
    # (Read "w" here as "omega".) We pre-multiply as many scalar values as
    # possible on these vectors to avoid operations on the full w matrix below.
    v = np.fft.fftfreq(img.shape[0])
    wy = (2 * np.pi * v * shift[0]).astype(np.float32).reshape(-1, 1)
    u = np.fft.fftfreq(img.shape[1])
    wx = (2 * np.pi * u * shift[1]).astype(np.float32)
    # Add column and row vector to get full expanded matrix of frequencies.
    w = wy + wx
    # We perform an explicit application of Euler's formula with careful
    # management of output arrays to avoid extra memory allocations and copies,
    # squeezing out some speed over the obvious np.exp(-1j*w).
    fshift = np.empty_like(img, dtype=np.complex64)
    np.cos(w, out=fshift.real)
    np.sin(w, out=fshift.imag)
    np.negative(fshift.imag, out=fshift.imag)
    # Perform the FFT, multiply in-place by the shift matrix, then IFFT.
    freq = pyfftw.builders.fft2(img, planner_effort='FFTW_ESTIMATE',
                                avoid_copy=True, auto_align_input=True,
                                auto_contiguous=True)()
    freq *= fshift
    img_s = pyfftw.builders.ifft2(freq, planner_effort='FFTW_ESTIMATE',
                                  avoid_copy=True, auto_align_input=True,
                                  auto_contiguous=True)()
    # Any non-zero imaginary component of the resulting array is due to
    # numerical error, so we can just return the real part.
    # FIXME need to zero out row(s) and column(s) we shifted away from,
    # since at this point we have a cyclic rotation rather than a shift.
    return img_s.real


def paste(target, img, pos):
    """Composite img into target using maximum intensity projection."""
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
    if np.issubdtype(img.dtype, float):
        np.clip(img, 0, 1, img)
    img = skimage.util.dtype.convert(img, target.dtype)
    target_slice[:, :] = np.maximum(target_slice, img)


def crop_like(img, target):
    if (img.shape[0] > target.shape[0]):
        img = img[:target.shape[0], :]
    if (img.shape[1] > target.shape[1]):
        img = img[:, :target.shape[1]]
    return img


def plot_edge_shifts(aligner, mosaic):
    plt.figure()
    ax = plt.gca()
    modest_image.imshow(ax, mosaic)
    h, w = aligner.reader.metadata.size
    # Bounding boxes denoting new tile positions.
    for xy in np.fliplr(aligner.positions):
        rect = mpatches.Rectangle(xy, w, h, color='black', fill=False, lw=0.5)
        ax.add_patch(rect)
    # Compute per-edge relative shifts from tile positions.
    edges = np.array(aligner.spanning_tree.edges())
    dist = aligner.metadata.positions - aligner.positions
    shifts = dist[edges[:, 0]] - dist[edges[:, 1]]
    shift_distances = np.linalg.norm(shifts, axis=1)
    # Spanning tree with nodes at new tile positions, edges colored by shift
    # distance (brighter = farther).
    nx.draw(
        aligner.spanning_tree, ax=ax, with_labels=True,
        pos=np.fliplr(aligner.centers), edge_color=shift_distances,
        edge_cmap=plt.get_cmap('Blues_r'), width=2, node_size=100, font_size=6
    )

def plot_edge_quality(aligner, mosaic):
    centers = aligner.reader.metadata.centers - aligner.reader.metadata.origin
    nrows, ncols = 1, 2
    if mosaic.shape[1] * 2 / mosaic.shape[0] < 4 / 3:
        nrows, ncols = ncols, nrows
    plt.figure()
    ax = plt.subplot(nrows, ncols,1)
    modest_image.imshow(ax, mosaic)
    error = np.array([aligner._cache[tuple(sorted(e))][1]
                      for e in aligner.neighbors_graph.edges()])
    # Manually center and scale data to 0-1, except infinity which is set to -1.
    # This lets us use the purple-green diverging color map to color the graph
    # edges and cause the "infinity" edges to disappear into the background
    # (which is itself purple).
    infs = error == np.inf
    error[infs] = -1
    error_f = error[~infs]
    emin = np.min(error_f)
    emax = np.max(error_f)
    error[~infs] = (error_f - emin) / (emax - emin)
    # Neighbor graph colored by edge alignment quality (brighter = better).
    nx.draw(
        aligner.neighbors_graph, ax=ax, with_labels=True,
        pos=np.fliplr(centers), edge_color=error, edge_vmin=-1, edge_vmax=1,
        edge_cmap=plt.get_cmap('PRGn'), width=2, node_size=100, font_size=6
    )
    ax = plt.subplot(nrows, ncols, 2)
    modest_image.imshow(ax, mosaic)
    # Spanning tree with nodes at original tile positions.
    nx.draw(
        aligner.spanning_tree, ax=ax, with_labels=True,
        pos=np.fliplr(centers), edge_color='royalblue',
        width=2, node_size=100, font_size=6
    )

def plot_layer_shifts(aligner, mosaic):
    plt.figure()
    ax = plt.gca()
    modest_image.imshow(ax, mosaic)
    h, w = aligner.metadata.size
    # Bounding boxes denoting new tile positions.
    for xy in np.fliplr(aligner.positions):
        rect = mpatches.Rectangle(xy, w, h, color='black', fill=False, lw=0.5)
        ax.add_patch(rect)
    # Neighbor graph with edges hidden, i.e. just show nodes.
    nx.draw(
        aligner.neighbors_graph, ax=ax, with_labels=True,
        pos=np.fliplr(aligner.centers), edge_color='none',
        node_size=100, font_size=6
    )
