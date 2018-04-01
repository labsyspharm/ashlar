from __future__ import division, print_function
import sys
import warnings
import xml.etree.ElementTree
try:
    import pathlib
except ImportError:
    import pathlib2 as pathlib
import jnius_config
import numpy as np
import scipy.ndimage
import skimage.util
import skimage.feature
import skimage.filters
import skimage.restoration.uft
import skimage.io
import skimage.exposure
import sklearn.linear_model
import pyfftw
import networkx as nx
import queue
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
try:
    import modest_image
except ImportError:
    modest_image = None

# Patch np.fft to use pyfftw so skimage utilities can benefit.
np.fft = pyfftw.interfaces.numpy_fft


if not jnius_config.vm_running:
    pkg_root = pathlib.Path(__file__).parent.resolve()
    bf_jar_path = pkg_root / 'jars' / 'loci_tools.jar'
    if not bf_jar_path.exists():
        raise RuntimeError("loci_tools.jar missing from distribution"
                           " (expected it at %s)" % bf_jar_path)
    jnius_config.add_classpath(str(bf_jar_path))

import jnius

JString = jnius.autoclass('java.lang.String')
DebugTools = jnius.autoclass('loci.common.DebugTools')
IFormatReader = jnius.autoclass('loci.formats.IFormatReader')
MetadataStore = jnius.autoclass('loci.formats.meta.MetadataStore')
ServiceFactory = jnius.autoclass('loci.common.services.ServiceFactory')
OMEXMLService = jnius.autoclass('loci.formats.services.OMEXMLService')
ImageReader = jnius.autoclass('loci.formats.ImageReader')

DebugTools.setRootLevel("ERROR")


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

    @property
    def pixel_dtype(self):
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

    _pixel_dtypes = {
        'uint16': np.uint16,
    }

    def __init__(self, path):
        self.path = path
        self._init_metadata()

    def __getstate__(self):
        state = self.__dict__.copy()
        del state['_metadata']
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)
        self._init_metadata()

    def _init_metadata(self):

        factory = ServiceFactory()
        service = jnius.cast(OMEXMLService, factory.getInstance(OMEXMLService))
        metadata = service.createOMEXMLMetadata()
        reader = ImageReader()
        reader.setMetadataStore(metadata)
        # FIXME Workaround for pyjnius #300 until there is a new release.
        # Passing a python string directly here corrupts the value under Python
        # 3, but explicitly converting it into a Java string works.
        path_jstring = JString(self.path)
        reader.setId(path_jstring)

        xml_content = metadata.dumpXML()
        self._metadata = metadata
        self._omexml_root = xml.etree.ElementTree.fromstring(xml_content)
        self.format_name = reader.getFormat()

    @property
    def num_images(self):
        count = self._metadata.imageCount
        # Skip final overview slide in Metamorph Slide Scan data if present.
        if (self.format_name == 'Metamorph STK'
            and 'overview' in self._metadata.getImageName(count - 1).lower()):
            count -= 1
        return count

    @property
    def num_channels(self):
        return self._metadata.getChannelCount(0)

    @property
    def pixel_size(self):
        values = []
        for dim in ('Y', 'X'):
            method = getattr(self._metadata, 'getPixelsPhysicalSize%s' % dim)
            v = method(0).value().doubleValue()
            values.append(v)
        if values[0] != values[1]:
            raise Exception("Can't handle non-square pixels (%f, %f)"
                            % tuple(values))
        return values[0]

    @property
    def pixel_dtype(self):
        return self._pixel_dtypes[self._metadata.getPixelsType(0).value]

    def tile_position(self, i):
        values = []
        for dim in ('Y', 'X'):
            method = getattr(self._metadata, 'getPlanePosition%s' % dim)
            # FIXME verify all planes have the same X,Y position.
            v = method(i, 0).value().doubleValue()
            values.append(v)
        position_microns = np.array(values, dtype=float)
        if self.format_name != 'Metamorph STK':
            # Invert Y so that stage position coordinates and image pixel
            # coordinates are aligned.
            position_microns *= [-1, 1]
        position_pixels = position_microns / self.pixel_size
        return position_pixels

    def tile_size(self, i):
        values = []
        for dim in ('Y', 'X'):
            method = getattr(self._metadata, 'getPixelsSize%s' % dim)
            v = method(i).value
            values.append(v)
        return np.array(values, dtype=int)


class BioformatsReader(Reader):

    def __init__(self, path):
        self.path = path
        self.metadata = BioformatsMetadata(self.path)
        self._init_ir()

    def __getstate__(self):
        state = self.__dict__.copy()
        del state['ir']
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)
        self._init_ir()

    def _init_ir(self):
        ImageReader = jnius.autoclass('loci.formats.ImageReader')
        self._reader = ImageReader()
        # FIXME Workaround for pyjnius #300 (see above for details).
        path_jstring = JString(self.path)
        self._reader.setId(path_jstring)

    def read(self, series, c):
        self._reader.setSeries(series)
        index = self._reader.getIndex(0, c, 0)
        byte_array = self._reader.openBytes(index)
        dtype = self.metadata.pixel_dtype
        shape = self.metadata.tile_size(series)
        img = np.frombuffer(byte_array.tostring(), dtype=dtype).reshape(shape)
        return img


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

    def __init__(self, reader, channel=0, max_shift=15, verbose=False):
        self.reader = reader
        self.channel = channel
        self.verbose = verbose
        # Unit is micrometers.
        self.max_shift = max_shift
        self.max_shift_pixels = self.max_shift / self.metadata.pixel_size
        self._cache = {}

    neighbors_graph = neighbors_graph

    def run(self):
        self.register_all()
        self.build_spanning_tree()
        self.calculate_positions()
        self.fit_model()

    def register_all(self):
        n = self.neighbors_graph.size()
        for i, (t1, t2) in enumerate(self.neighbors_graph.edges(), 1):
            if self.verbose:
                sys.stdout.write('\r    aligning edge %d/%d' % (i, n))
                sys.stdout.flush()
            self.register_pair(t1, t2)
        if self.verbose:
            print()

    def build_spanning_tree(self):
        # Note that this may be disconnected, so it's technically a forest.
        g = nx.Graph()
        g.add_nodes_from(self.neighbors_graph)
        g.add_weighted_edges_from(
            (t1, t2, error)
            for (t1, t2), (_, error) in self._cache.items()
            if np.isfinite(error)
        )
        spanning_tree = nx.Graph()
        spanning_tree.add_nodes_from(g)
        for cc in nx.connected_component_subgraphs(g):
            center = nx.center(cc)[0]
            paths = nx.single_source_dijkstra_path(cc, center).values()
            for path in paths:
                spanning_tree.add_path(path)
        self.spanning_tree = spanning_tree

    def calculate_positions(self):
        shifts = {}
        for cc in nx.connected_component_subgraphs(self.spanning_tree):
            center = nx.center(cc)[0]
            shifts[center] = np.array([0, 0])
            for edge in nx.traversal.bfs_edges(cc, center):
                source, dest = edge
                if source not in shifts:
                    source, dest = dest, source
                shift = self.register_pair(source, dest)[0]
                shifts[dest] = shifts[source] + shift
        self.shifts = np.array([s for _, s in sorted(shifts.items())])
        self.positions = self.metadata.positions + self.shifts

    def fit_model(self):
        components = sorted(
            nx.connected_component_subgraphs(self.spanning_tree),
            key=len, reverse=True
        )
        # Fit LR model on positions of largest connected component.
        cc0 = list(components[0])
        self.lr = sklearn.linear_model.LinearRegression()
        self.lr.fit(self.metadata.positions[cc0], self.positions[cc0])
        # Adjust position of remaining components so their centroids match
        # the predictions of the model.
        for cc in components[1:]:
            nodes = list(cc)
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
            # Phase correlation can only report a shift of up to half of the
            # image size in any given direction. Here we calculate the tile
            # overlap image size that will support observing the maximum allowed
            # shift.
            max_its_size = self.max_shift_pixels * 2
            results = []
            # Perform registration with increasingly larger minimum overlap
            # sizes, starting from the "native" overlap and repeatedly doubling
            # it until we reach or exceed the maximum required overlap
            # calculated above. We then choose the result with the lowest
            # error. This is necessary because the phase correlation can't
            # detect a shift greater than half of the overlap size.
            native_its = self.intersection(t1, t2, None)
            min_size = native_its.shape
            while True:
                shift, error = self._register(t1, t2, min_size)
                # Constrain shift.
                if any(np.abs(shift) > self.max_shift_pixels):
                    shift[:] = 0
                    error = np.inf
                results.append([shift, error])
                if all(min_size > max_its_size):
                    break
                min_size *= 2
            # Take result with lowest error.
            shift, error = min(results, key=lambda x: x[1])
            self._cache[key] = (shift, error)
        if t1 > t2:
            shift = -shift
        # Return copy of shift to prevent corruption of cached values.
        return shift.copy(), error

    def _register(self, t1, t2, min_size):
        # Perform the actual pixel-based alignment, given a minimum size for the
        # overlap region.
        its, img1, img2 = self.overlap(t1, t2, min_size)
        img1_f = fft2(whiten(img1))
        img2_f = fft2(whiten(img2))
        shift, error, _ = skimage.feature.register_translation(
            img1_f, img2_f, 10, 'fourier'
        )
        # Account for padding, flipping the sign depending on the direction
        # between the tiles.
        p1, p2 = self.metadata.positions[[t1, t2]]
        sx = 1 if p1[1] >= p2[1] else -1
        sy = 1 if p1[0] >= p2[0] else -1
        shift += its.padding * [sy, sx]
        return shift, error

    def intersection(self, t1, t2, min_size):
        corners1 = self.metadata.positions[[t1, t2]]
        corners2 = corners1 + self.metadata.size
        return Intersection(corners1, corners2, min_size)

    def crop(self, tile, offset, shape):
        img = self.reader.read(series=tile, c=self.channel)
        return crop(img, offset, shape)

    def overlap(self, t1, t2, min_size):
        its = self.intersection(t1, t2, min_size)
        img1 = self.crop(t1, its.offsets[0], its.shape)
        img2 = self.crop(t2, its.offsets[1], its.shape)
        return its, img1, img2

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

    def debug(self, t1, t2, min_size=None):
        shift, _ = self._register(t1, t2, min_size)
        its, o1, o2 = self.overlap(t1, t2, min_size)
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
        origin = np.array(corr.shape) // 2
        plt.plot(origin[1], origin[0], 'r+')
        # FIXME This is wrong when t1 > t2.
        shift += origin + its.padding
        plt.plot(shift[1], shift[0], 'rx')
        plt.colorbar()
        plt.tight_layout(0, 0, 0)


class LayerAligner(object):

    def __init__(self, reader, reference_aligner, channel=None, max_shift=15,
                 verbose=False):
        self.reader = reader
        self.reference_aligner = reference_aligner
        if channel is None:
            channel = reference_aligner.channel
        self.channel = channel
        # Unit is micrometers.
        self.max_shift = max_shift
        self.max_shift_pixels = self.max_shift / self.metadata.pixel_size
        self.verbose = verbose
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
            print()

    def calculate_positions(self):
        self.positions = self.reference_aligner.positions + self.shifts
        self.constrain_positions()
        self.centers = self.positions + self.metadata.size / 2

    def constrain_positions(self):
        # Computed shifts of exactly 0,0 seem to result from failed
        # registration. We need to throw those out for this purpose.
        discard = (self.shifts == 0).all(axis=1)
        # Take the median of registered shifts to determine the offset
        # (translation) from the reference image to this one.
        offset = np.nan_to_num(np.median(self.shifts[~discard], axis=0))
        # Here we assume the fitted linear model from the reference image is
        # still appropriate, apart from the extra offset we just computed.
        predictions = self.reference_aligner.lr.predict(self.metadata.positions)
        # Discard any tile registration that's too far from the linear model,
        # replacing it with the relevant model prediction.
        distance = np.linalg.norm(self.positions - predictions - offset, axis=1)
        max_dist = self.max_shift_pixels
        extremes = distance > max_dist
        # Recalculate the mean shift, also ignoring the extreme values.
        discard |= extremes
        self.offset = np.nan_to_num(np.mean(self.shifts[~discard], axis=0))
        # Fill in discarded shifts from the predictions.
        self.positions[discard] = predictions[discard] + self.offset

    def register(self, t):
        """Return relative shift between images and the alignment error."""
        its, ref_img, img = self.overlap(t)
        ref_img_f = fft2(whiten(ref_img))
        img_f = fft2(whiten(img))
        shift, error, _ = skimage.feature.register_translation(
            ref_img_f, img_f, 10, 'fourier'
        )
        # Add reported difference in stage positions.
        shift += its.offsets[0]
        return shift, error

    def intersection(self, t):
        corners1 = np.vstack([self.reference_positions[t],
                              self.tile_positions[t]])
        corners2 = corners1 + self.reader.metadata.size
        its = Intersection(corners1, corners2)
        its.shape = its.shape // 32 * 32
        return its

    def overlap(self, t):
        its = self.intersection(t)
        ref_t = self.reference_idx[t]
        img1 = self.reference_aligner.reader.read(series=ref_t, c=0)
        img2 = self.reader.read(series=t, c=self.channel)
        ov1 = crop(img1, its.offsets[0], its.shape)
        ov2 = crop(img2, its.offsets[1], its.shape)
        return its, ov1, ov2

    @property
    def metadata(self):
        return self.reader.metadata

    def debug(self, t):
        shift, _ = self.register(t)
        its, o1, o2 = self.overlap(t)
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
        origin = np.array(corr.shape) // 2
        plt.plot(origin[1], origin[0], 'r+')
        shift += origin - its.offsets[0]
        plt.plot(shift[1], shift[0], 'rx')
        plt.tight_layout(0, 0, 0)


class Intersection(object):

    def __init__(self, corners1, corners2, min_size=None):
        if min_size is None:
            min_size = np.zeros(2)
        self._calculate(corners1, corners2, min_size)

    def _calculate(self, corners1, corners2, min_size):
        max_shape = (corners2 - corners1).max(axis=0)
        min_size = min_size.clip(0, max_shape)
        position = corners1.max(axis=0)
        initial_shape = np.ceil(corners2.min(axis=0) - position).astype(int)
        if any(initial_shape <= 0):
            raise ValueError("Tiles do not intersect")
        clipped_shape = np.maximum(initial_shape, min_size)
        self.shape = np.ceil(clipped_shape).astype(int)
        self.padding = self.shape - initial_shape
        self.offsets = np.maximum(position - corners1 - self.padding, 0)

    def __repr__(self):
        s = 'shape: {0.shape}\npadding: {0.padding}\noffsets:\n{0.offsets}'
        return s.format(self)


class Mosaic(object):

    def __init__(self, aligner, shape, filename_format, channels=None,
                 ffp_path=None, dfp_path=None, verbose=False):
        self.aligner = aligner
        self.shape = tuple(shape)
        self.filename_format = filename_format
        self.channels = self._sanitize_channels(channels)
        self.dtype = np.uint16
        self._load_correction_profiles(dfp_path, ffp_path)
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

    def _load_correction_profiles(self, dfp_path, ffp_path):
        if dfp_path or ffp_path:
            c = max(self.channels) + 1
            self.dfp = skimage.io.imread(dfp_path) if dfp_path else np.zeros(c)
            self.ffp = skimage.io.imread(ffp_path) if ffp_path else np.ones(c)
            self.dfp /= np.iinfo(self.dtype).max
            self.do_correction = True
        else:
            self.do_correction = False

    def run(self, mode='write', debug=False):
        if mode not in ('write', 'return'):
            raise ValueError('Invalid mode')
        num_tiles = len(self.aligner.positions)
        all_images = []
        if debug:
            node_colors = nx.greedy_color(self.aligner.neighbors_graph)
            num_colors = max(node_colors.values()) + 1
            if num_colors > 3:
                raise ValueError("neighbor graph requires more than 3 colors")
        for channel in self.channels:
            if self.verbose:
                print('    Channel %d:' % channel)
            if not debug:
                mosaic_image = np.zeros(self.shape, self.dtype)
            else:
                mosaic_image = np.zeros(self.shape + (3,), np.float32)
            for tile, position in enumerate(self.aligner.positions):
                if self.verbose:
                    sys.stdout.write('\r        merging tile %d/%d'
                                     % (tile + 1, num_tiles))
                    sys.stdout.flush()
                tile_image = self.aligner.reader.read(c=channel, series=tile)
                tile_image = self.correct_illumination(tile_image, channel)
                if debug:
                    color_channel = node_colors[tile]
                    rgb_image = np.zeros(tile_image.shape + (3,),
                                         tile_image.dtype)
                    rgb_image[:,:,color_channel] = tile_image
                    tile_image = rgb_image
                func = pastefunc_blend if not debug else np.add
                paste(mosaic_image, tile_image, position, func=func)
            if debug:
                np.clip(mosaic_image, 0, 1, out=mosaic_image)
                w = int(1e6)
                mi_flat = mosaic_image.reshape(-1, 3)
                for p in np.arange(0, mi_flat.shape[0], w, dtype=int):
                    mi_flat[p:p+w] = skimage.exposure.adjust_gamma(
                        mi_flat[p:p+w], 1/2.2
                    )
            if self.verbose:
                print()
            if mode == 'write':
                filename = self.filename_format.format(channel=channel)
                self.filenames.append(filename)
                if self.verbose:
                    print("        writing %s" % filename)
                with warnings.catch_warnings():
                    warnings.filterwarnings(
                        'ignore', r'.* is a low contrast image', UserWarning,
                        '^skimage\.io'
                    )
                    skimage.io.imsave(filename, mosaic_image)
            elif mode == 'return':
                all_images.append(mosaic_image)
        if mode == 'return':
            return all_images

    def correct_illumination(self, img, channel):
        if self.do_correction:
            img = skimage.img_as_float(img, force_copy=True)
            img -= self.dfp[..., channel]
            img /= self.ffp[..., channel]
            img.clip(0, 1, out=img)
        return img


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


def crop(img, offset, shape):
    # Note that this only crops to the nearest whole-pixel offset.
    start = offset.astype(int)
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


def paste(target, img, pos, func=None):
    """Composite img into target."""
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
    if img.ndim == 2:
        img = scipy.ndimage.shift(img, pos_f)
    else:
        for c in range(img.shape[2]):
            img[...,c] = scipy.ndimage.shift(img[...,c], pos_f)
    # Drop one pixel on all edges to avoid artifacts from the subpixel shift.
    # FIXME We should only drop the edges opposite the shift direction.
    img = img[1:-1,1:-1]
    target_slice = target_slice[1:-1,1:-1]
    if np.issubdtype(img.dtype, np.floating):
        np.clip(img, 0, 1, img)
    img = skimage.util.dtype.convert(img, target.dtype)
    if func is None:
        target_slice[:] = img
    elif isinstance(func, np.ufunc):
        func(target_slice, img, out=target_slice)
    else:
        target_slice[:] = func(target_slice, img)


def pastefunc_blend(target, img):
    # Linear blend based on distance to unfilled space in target.
    dist = scipy.ndimage.distance_transform_cdt(target)
    dmax = dist.max()
    if dmax == 0:
        alpha = 0
    else:
        alpha = dist / dist.max()
    return target * alpha + img * (1 - alpha)


def crop_like(img, target):
    if (img.shape[0] > target.shape[0]):
        img = img[:target.shape[0], :]
    if (img.shape[1] > target.shape[1]):
        img = img[:, :target.shape[1]]
    return img


def plot_edge_shifts(aligner, img=None, bounds=True):
    fig = plt.figure()
    ax = plt.gca()
    draw_mosaic_image(ax, aligner, img)
    h, w = aligner.reader.metadata.size
    if bounds:
        # Bounding boxes denoting new tile positions.
        for xy in np.fliplr(aligner.positions):
            rect = mpatches.Rectangle(xy, w, h, color='black', fill=False,
                                      lw=0.5)
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
    fig.set_facecolor('black')


def plot_edge_quality(aligner, img=None):
    centers = aligner.reader.metadata.centers - aligner.reader.metadata.origin
    nrows, ncols = 1, 2
    if aligner.mosaic_shape[1] * 2 / aligner.mosaic_shape[0] > 4 / 3:
        nrows, ncols = ncols, nrows
    fig = plt.figure()
    ax = plt.subplot(nrows, ncols,1)
    draw_mosaic_image(ax, aligner, img)
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
    draw_mosaic_image(ax, aligner, img)
    # Spanning tree with nodes at original tile positions.
    nx.draw(
        aligner.spanning_tree, ax=ax, with_labels=True,
        pos=np.fliplr(centers), edge_color='royalblue',
        width=2, node_size=100, font_size=6
    )
    fig.set_facecolor('black')


def plot_layer_shifts(aligner, img=None):
    fig = plt.figure()
    ax = plt.gca()
    draw_mosaic_image(ax, aligner, img)
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
    fig.set_facecolor('black')


def draw_mosaic_image(ax, aligner, img):
    if img is not None:
        if sys.version_info[0] != 2:
            warnings.warn('ModestImage module (required for image display)'
                          ' is only compatible with Python 2')
            img = None
        elif modest_image is None:
            warnings.warn('Please install ModestImage for image display')
            img = None
    if img is not None:
        modest_image.imshow(ax, img)
    else:
        h, w = aligner.mosaic_shape
        # Draw a single-pixel image in the lowest color in the colormap,
        # stretched across the same extent that the real image would be.
        # This makes the graph edge colors visible even if there's no image.
        ax.imshow([[0]], extent=(-0.5, w-0.5, h-0.5, -0.5))
