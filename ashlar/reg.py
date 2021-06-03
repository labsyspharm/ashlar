import sys
import warnings
import re
import xml.etree.ElementTree
import io
import uuid
import struct
import pathlib
import jnius_config
import numpy as np
import scipy.spatial.distance
import scipy.fft
import skimage.util
import skimage.util.dtype
import skimage.io
import skimage.exposure
import skimage.transform
import sklearn.linear_model
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.cm as mcm
import matplotlib.patches as mpatches
import matplotlib.patheffects as mpatheffects
from . import utils
from . import thumbnail
from . import __version__ as _version


if not jnius_config.vm_running:
    pkg_root = pathlib.Path(__file__).parent.resolve()
    bf_jar_path = pkg_root / 'jars' / 'loci_tools.jar'
    if not bf_jar_path.exists():
        raise RuntimeError("loci_tools.jar missing from distribution"
                           " (expected it at %s)" % bf_jar_path)
    jnius_config.add_classpath(str(bf_jar_path))

import jnius

DebugTools = jnius.autoclass('loci.common.DebugTools')
IFormatReader = jnius.autoclass('loci.formats.IFormatReader')
MetadataRetrieve = jnius.autoclass('ome.xml.meta.MetadataRetrieve')
ServiceFactory = jnius.autoclass('loci.common.services.ServiceFactory')
OMEXMLService = jnius.autoclass('loci.formats.services.OMEXMLService')
ChannelSeparator = jnius.autoclass('loci.formats.ChannelSeparator')
DynamicMetadataOptions = jnius.autoclass('loci.formats.in.DynamicMetadataOptions')
UNITS = jnius.autoclass('ome.units.UNITS')
DebugTools.enableLogging("ERROR")


# TODO:
# - Write tables with summary information about alignments.


class Metadata(object):

    @property
    def _num_images(self):
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
    def num_images(self):
        return self._num_images

    @property
    def positions(self):
        if not hasattr(self, '_positions'):
            self._positions = np.vstack([
                self.tile_position(i) for i in range(self._num_images)
            ])
        return self._positions

    @property
    def size(self):
        if not hasattr(self, '_size'):
            s0 = self.tile_size(0)
            image_ids = range(1, self._num_images)
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


class PlateMetadata(Metadata):

    def __init__(self):
        super(PlateMetadata, self).__init__()
        self.set_active_plate_well(None, None)

    @property
    def num_plates(self):
        raise NotImplementedError

    @property
    def num_wells(self):
        raise NotImplementedError

    @property
    def plate_well_series(self):
        raise NotImplementedError

    def plate_name(self, i):
        raise NotImplementedError

    def well_name(self, plate, i):
        raise NotImplementedError

    def set_active_plate_well(self, plate, well):
        if (plate is None) ^ (well is None):
            raise ValueError("plate and well must be both set or both None")
        self.active_plate = plate
        self.active_well = well

    @property
    def active_series(self):
        if self.active_plate is None:
            return range(self._num_images)
        else:
            return self.plate_well_series[self.active_plate][self.active_well]

    @property
    def plate_names(self):
        if not hasattr(self, '_plate_names'):
            self._plate_names = [
                self.plate_name(i) for i in range(self.num_plates)
            ]
        return self._plate_names

    @property
    def well_names(self):
        if not hasattr(self, '_well_names'):
            self._well_names = [
                [self.well_name(p, i) for i in range(num_plate_wells)]
                for p, num_plate_wells in enumerate(self.num_wells)
            ]
        return self._well_names

    @Metadata.num_images.getter
    def num_images(self):
        return len(self.active_series)

    @Metadata.positions.getter
    def positions(self):
        return Metadata.positions.fget(self)[self.active_series]

    # FIXME Metadata.grid_dimensions should be overriden here or removed.


class Reader(object):

    def read(self, series, c):
        raise NotImplementedError


class PlateReader(Reader):
    # No API here, just a way to signal that a subclass's metadata class
    # inherits from PlateMetadata. This is probably a sign that the
    # architectural split between Metadata and Reader should be reconsidered.
    pass


class BioformatsMetadata(PlateMetadata):

    _pixel_dtypes = {
        'uint8': np.dtype(np.uint8),
        'uint16': np.dtype(np.uint16),
    }

    _ome_dtypes = {v: k for k, v in _pixel_dtypes.items()}

    def __init__(self, path):
        super(BioformatsMetadata, self).__init__()
        self.path = path
        self._init_metadata()

    def __getstate__(self):
        state = self.__dict__.copy()
        del state['_reader'], state['_metadata'], state['_omexml_root']
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)
        self._init_metadata()

    def _init_metadata(self):

        factory = ServiceFactory()
        service = jnius.cast(OMEXMLService, factory.getInstance(OMEXMLService))
        metadata = service.createOMEXMLMetadata()
        self._reader = ChannelSeparator()
        self._reader.setMetadataStore(metadata)
        # For multi-scene .CZI files, we need raw tiles instead of the
        # auto-stitched mosaic and we don't want labels or overview images
        options = DynamicMetadataOptions()
        options.setBoolean('zeissczi.autostitch', False)
        options.setBoolean('zeissczi.attachments', False)
        self._reader.setMetadataOptions(options)
        self._reader.setId(self.path)

        xml_content = service.getOMEXML(metadata)
        self._metadata = jnius.cast(MetadataRetrieve, metadata)
        self._omexml_root = xml.etree.ElementTree.fromstring(xml_content)
        self.format_name = self._reader.getFormat()

    @property
    def _num_images(self):
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
    def num_plates(self):
        return self._metadata.getPlateCount()

    @property
    def num_wells(self):
        return [self._metadata.getWellCount(i) for i in range(self.num_plates)]

    @property
    def plate_well_series(self):
        if hasattr(self, '_plate_well_series'):
            return self._plate_well_series
        # FIXME Store slice objects to save resources where possible.
        series = [
            [
                np.array([
                    self._metadata.getWellSampleIndex(p, w, s).value
                    for s in range(self._metadata.getWellSampleCount(p, w))
                ], dtype=int)
                for w in range(num_wells)
            ]
            for p, num_wells in enumerate(self.num_wells)
        ]
        return series

    @property
    def pixel_size(self):
        values = []
        for dim in ('Y', 'X'):
            method = getattr(self._metadata, 'getPixelsPhysicalSize%s' % dim)
            v_units = method(0)
            if v_units is None:
                warn_data(
                    "Pixel size undefined; falling back to 1.0 \u03BCm."
                )
                value = 1.0
            else:
                value = v_units.value(UNITS.MICROMETER).doubleValue()
            values.append(value)
        if values[0] != values[1]:
            raise Exception("Can't handle non-square pixels (%f, %f)"
                            % tuple(values))
        return values[0]

    @property
    def pixel_dtype(self):
        return self._pixel_dtypes[self._metadata.getPixelsType(0).value]

    def plate_name(self, i):
        return self._metadata.getPlateName(i)

    @property
    def well_naming(self):
        if not hasattr(self, '_well_naming'):
            _well_naming = []
            for p in range(self.num_plates):
                row_nc = self._metadata.getPlateRowNamingConvention(p)
                column_nc = self._metadata.getPlateColumnNamingConvention(p)
                if row_nc is not None:
                    row_nc = row_nc.value
                else:
                    row_nc = 'letter'
                if column_nc is not None:
                    column_nc = column_nc.value
                else:
                    column_nc = 'number'
                if row_nc not in ('letter', 'number') or column_nc != 'number':
                    raise RuntimeError(
                        "Can't handle well naming convention row={} column={}"
                        .format(row_nc, column_nc)
                    )
                _well_naming.append([row_nc, column_nc])
            self._well_naming = _well_naming
        return self._well_naming

    def well_name(self, plate, i):
        row = self._metadata.getWellRow(plate, i).value
        column = self._metadata.getWellColumn(plate, i).value
        row_nc, column_nc = self.well_naming[plate]
        # FIXME Support formatting with 384/1536-well plates.
        assert row_nc in ('letter', 'number')
        assert column_nc == 'number'
        if row_nc == 'number':
            row_fmt = '{:02}'.format(row + 1)
        else:
            row_fmt = chr(ord('A') + row)
        column_fmt = '{:02}'.format(column + 1)
        return row_fmt + column_fmt

    def tile_position(self, i):
        planeCount = self._metadata.getPlaneCount(i)
        values = []
        for dim in ('Y', 'X'):
            method = getattr(self._metadata, 'getPlanePosition%s' % dim)
            # FIXME verify all planes have the same X,Y position.
            if planeCount > 0:
                # Returns None if planePositionX/Y not defined.
                v_units = method(i, 0)
            else:
                # Simple file formats don't have planes at all.
                v_units = None
            if v_units is None:
                warn_data(
                    "Stage coordinates undefined; falling back to (0, 0)."
                )
                values = [0.0, 0.0]
                break
            else:
                v = v_units.value(UNITS.MICROMETER)
                if v is None:
                    # Conversion failed, which usually happens when the unit is
                    # "reference frame". Proceed as if it's actually microns but
                    # emit a warning.
                    warn_data(
                        "Stage coordinates' measurement unit is undefined;"
                        " assuming \u03BCm."
                    )
                    v = v_units.value()
                value = v.doubleValue()
            values.append(value)
        position_microns = np.array(values, dtype=float)
        # Invert Y so that stage position coordinates and image pixel
        # coordinates are aligned (most formats seem to work this way).
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


class BioformatsReader(PlateReader):

    def __init__(self, path, plate=None, well=None):
        self.path = path
        self.metadata = BioformatsMetadata(self.path)
        self.metadata.set_active_plate_well(plate, well)

    def read(self, series, c):
        self.metadata._reader.setSeries(self.metadata.active_series[series])
        index = self.metadata._reader.getIndex(0, c, 0)
        byte_array = self.metadata._reader.openBytes(index)
        dtype = self.metadata.pixel_dtype
        shape = self.metadata.tile_size(series)
        img = np.frombuffer(byte_array.tostring(), dtype=dtype).reshape(shape)
        return img


class CachingReader(Reader):
    """Wraps a reader to provide tile image caching."""

    def __init__(self, reader, channel):
        self.reader = reader
        self.channel = channel
        self._cache = {}

    @property
    def metadata(self):
        return self.reader.metadata

    def read(self, series, c):
        if c == self.channel and series in self._cache:
            img = self._cache[series]
        else:
            img = self.reader.read(series, c)
        if c == self.channel and series not in self._cache:
            self._cache[series] = img
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
        max_distance = aligner.metadata.size.max() + 1
        edges = zip(*np.nonzero((sp > 0) & (sp < max_distance)))
        graph = nx.from_edgelist(edges)
        graph.add_nodes_from(range(aligner.metadata.num_images))
        aligner._neighbors_graph = graph
    return aligner._neighbors_graph


class EdgeAligner(object):

    def __init__(
        self, reader, channel=0, max_shift=15, false_positive_ratio=0.01,
        randomize=False, filter_sigma=0.0, do_make_thumbnail=True, verbose=False
    ):
        self.channel = channel
        self.reader = CachingReader(reader, self.channel)
        self.verbose = verbose
        # Unit is micrometers.
        self.max_shift = max_shift
        self.max_shift_pixels = self.max_shift / self.metadata.pixel_size
        self.false_positive_ratio = false_positive_ratio
        self.randomize = randomize
        self.filter_sigma = filter_sigma
        self.do_make_thumbnail = do_make_thumbnail
        self._cache = {}

    neighbors_graph = neighbors_graph

    def run(self):
        self.make_thumbnail()
        self.check_overlaps()
        self.compute_threshold()
        self.register_all()
        self.build_spanning_tree()
        self.calculate_positions()
        self.fit_model()

    def make_thumbnail(self):
        if not self.do_make_thumbnail:
            return
        self.reader.thumbnail = thumbnail.make_thumbnail(
            self.reader, channel=self.channel
        )

    def check_overlaps(self):
        # This might be better addressed by removing the +1 from the
        # neighbors_graph max_distance calculation and ensuring the graph is
        # fully connected.
        pos = self.metadata.positions
        overlaps = np.array([
            self.metadata.size - abs(pos[t1] - pos[t2])
            for t1, t2 in self.neighbors_graph.edges
        ])
        failures = np.any(overlaps < 1, axis=1) if len(overlaps) else []
        if len(failures) and all(failures):
            warn_data("No tiles overlap, attempting alignment anyway.")
        elif any(failures):
            warn_data("Some neighboring tiles have zero overlap.")

    def compute_threshold(self):
        # Compute error threshold for rejecting aligments. We generate a
        # distribution of error scores for many known non-overlapping image
        # regions and take a certain percentile as the maximum allowable error.
        # The percentile becomes our accepted false-positive ratio.
        edges = self.neighbors_graph.edges
        num_tiles = self.metadata.num_images
        # If not enough tiles overlap to matter, skip this whole thing.
        if len(edges) <= 1:
            self.errors_negative_sampled = np.empty(0)
            self.max_error = np.inf
            return
        widths = np.array([
            self.intersection(t1, t2).shape.min()
            for t1, t2 in edges
        ])
        w = widths.max()
        max_offset = self.metadata.size[0] - w
        # Number of possible pairs minus number of actual neighbor pairs.
        num_distant_pairs = num_tiles * (num_tiles - 1) // 2 - len(edges)
        # Reduce permutation count for small datasets -- there are fewer
        # possible truly distinct strips with fewer tiles. The calculation here
        # is just a heuristic, not rigorously derived.
        n = 1000 if num_distant_pairs > 8 else (num_distant_pairs + 1) * 10
        pairs = np.empty((n, 2), dtype=int)
        offsets = np.empty((n, 2), dtype=int)
        # Generate n random non-overlapping image strips. Strips are always
        # horizontal, across the entire image width.
        max_tries = 100
        if self.randomize is False:
            random_state = np.random.RandomState(0)
        else:
            random_state = np.random.RandomState()
        for i in range(n):
            # Limit tries to avoid infinite loop in pathological cases.
            for current_try in range(max_tries):
                t1, t2 = random_state.randint(self.metadata.num_images, size=2)
                o1, o2 = random_state.randint(max_offset, size=2)
                # Check for non-overlapping strips and abort the retry loop.
                if t1 != t2 and (t1, t2) not in edges:
                    # Different, non-neighboring tiles -- always OK.
                    break
                elif t1 == t2 and abs(o1 - o2) > w:
                    # Same tile OK if strips don't overlap within the image.
                    break
                elif (t1, t2) in edges:
                    # Neighbors OK if either strip is entirely outside the
                    # expected overlap region (based on nominal positions).
                    its = self.intersection(t1, t2, np.repeat(w, 2))
                    ioff1, ioff2 = its.offsets[:, 0]
                    if (
                        its.shape[0] > its.shape[1]
                        or o1 < ioff1 - w or o1 > ioff1 + w
                        or o2 < ioff2 - w or o2 > ioff2 + w
                    ):
                        break
            else:
                # Retries exhausted. This should be very rare.
                warn_data(
                    "Could not find non-overlapping strips in {max_tries} tries"
                )
            pairs[i] = t1, t2
            offsets[i] = o1, o2
        errors = np.empty(n)
        for i, ((t1, t2), (offset1, offset2)) in enumerate(zip(pairs, offsets)):
            if self.verbose and (i % 10 == 9 or i == n - 1):
                sys.stdout.write(
                    '\r    quantifying alignment error %d/%d' % (i + 1, n)
                )
                sys.stdout.flush()
            img1 = self.reader.read(t1, self.channel)[offset1:offset1+w, :]
            img2 = self.reader.read(t2, self.channel)[offset2:offset2+w, :]
            _, errors[i] = utils.register(img1, img2, self.filter_sigma, upsample=1)
        if self.verbose:
            print()
        self.errors_negative_sampled = errors
        self.max_error = np.percentile(errors, self.false_positive_ratio * 100)

    def register_all(self):
        n = self.neighbors_graph.size()
        for i, (t1, t2) in enumerate(self.neighbors_graph.edges, 1):
            if self.verbose:
                sys.stdout.write('\r    aligning edge %d/%d' % (i, n))
                sys.stdout.flush()
            self.register_pair(t1, t2)
        if self.verbose:
            print()
        self.all_errors = np.array([x[1] for x in self._cache.values()])
        # Set error values above the threshold to infinity.
        for k, v in self._cache.items():
            if v[1] > self.max_error or any(np.abs(v[0]) > self.max_shift_pixels):
                self._cache[k] = (v[0], np.inf)

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
        for c in nx.connected_components(g):
            cc = g.subgraph(c)
            center = nx.center(cc)[0]
            paths = nx.single_source_dijkstra_path(cc, center).values()
            for path in paths:
                nx.add_path(spanning_tree, path)
        self.spanning_tree = spanning_tree

    def calculate_positions(self):
        shifts = {}
        for c in nx.connected_components(self.spanning_tree):
            cc = self.spanning_tree.subgraph(c)
            center = nx.center(cc)[0]
            shifts[center] = np.array([0, 0])
            for edge in nx.traversal.bfs_edges(cc, center):
                source, dest = edge
                if source not in shifts:
                    source, dest = dest, source
                shift = self.register_pair(source, dest)[0]
                shifts[dest] = shifts[source] + shift
        if shifts:
            self.shifts = np.array([s for _, s in sorted(shifts.items())])
            self.positions = self.metadata.positions + self.shifts
        else:
            # TODO: fill in shifts and positions with 0x2 arrays
            raise NotImplementedError("No images")

    def fit_model(self):
        components = sorted(
            nx.connected_components(self.spanning_tree),
            key=len, reverse=True
        )
        # Fit LR model on positions of largest connected component.
        cc0 = list(components[0])
        self.lr = sklearn.linear_model.LinearRegression()
        self.lr.fit(self.metadata.positions[cc0], self.positions[cc0])
        # Fix up degenerate transform matrix (e.g. when we have only one tile).
        if (self.lr.coef_ == 0).all():
            self.lr.coef_ = np.diag(np.ones(2))
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
            # We test a series of increasing overlap window sizes to help avoid
            # missing alignments when the stage position error is large relative
            # to the tile overlap. Simply using a large overlap in all cases
            # limits the maximum achievable correlation thus increasing the
            # error metric, leading to worse overall results. The window size
            # starts at the nominal size and doubles until it's at least 10% of
            # the tile size. If the nominal overlap is already 10% or greater,
            # we only use that one size.
            smin = self.intersection(key[0], key[1]).shape
            smax = np.round(self.metadata.size * 0.1)
            sizes = [smin]
            while any(sizes[-1] < smax):
                sizes.append(sizes[-1] * 2)
            results = [self._register(key[0], key[1], s) for s in sizes]
            # Use the shift from the window size that gave the lowest error.
            shift, _ = min(results, key=lambda r: r[1])
            # Extract the images from the nominal overlap window but with the
            # shift applied to the second tile's position, and compute the error
            # metric on these images. This should be even lower than the error
            # computed above.
            _, o1, o2 = self.overlap(key[0], key[1], shift=shift)
            error = utils.nccw(o1, o2, self.filter_sigma)
            self._cache[key] = (shift, error)
        if t1 > t2:
            shift = -shift
        # Return copy of shift to prevent corruption of cached values.
        return shift.copy(), error

    def _register(self, t1, t2, min_size=0):
        its, img1, img2 = self.overlap(t1, t2, min_size)
        # Account for padding, flipping the sign depending on the direction
        # between the tiles.
        p1, p2 = self.metadata.positions[[t1, t2]]
        sx = 1 if p1[1] >= p2[1] else -1
        sy = 1 if p1[0] >= p2[0] else -1
        padding = its.padding * [sy, sx]
        shift, error = utils.register(img1, img2, self.filter_sigma)
        shift += padding
        return shift, error

    def intersection(self, t1, t2, min_size=0, shift=None):
        corners1 = self.metadata.positions[[t1, t2]]
        if shift is not None:
            corners1[1] += shift
        corners2 = corners1 + self.metadata.size
        return Intersection(corners1, corners2, min_size)

    def crop(self, tile, offset, shape):
        img = self.reader.read(series=tile, c=self.channel)
        return utils.crop(img, offset, shape)

    def overlap(self, t1, t2, min_size=0, shift=None):
        its = self.intersection(t1, t2, min_size, shift)
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

    def debug(self, t1, t2, min_size=0):
        shift, _ = self._register(t1, t2, min_size)
        its, o1, o2 = self.overlap(t1, t2, min_size)
        w1 = utils.whiten(o1, self.filter_sigma)
        w2 = utils.whiten(o2, self.filter_sigma)
        corr = scipy.fft.fftshift(np.abs(scipy.fft.ifft2(
            scipy.fft.fft2(w1) * scipy.fft.fft2(w2).conj()
        )))
        corr /= (np.linalg.norm(w1) * np.linalg.norm(w2))
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
        plt.imshow(corr, vmin=np.exp(-10))
        cbar = plt.colorbar()
        cbar.ax.yaxis.set_major_locator(
            plt.FixedLocator(cbar.mappable.get_clim())
        )
        cbar.ax.yaxis.set_major_formatter(
            plt.FuncFormatter(lambda x, pos: "{:.2f}".format(-np.log(x)))
        )
        origin = np.array(corr.shape) // 2
        plt.plot(origin[1], origin[0], 'r+')
        # FIXME This is wrong when t1 > t2.
        shift += origin + its.padding
        plt.plot(shift[1], shift[0], 'rx')
        plt.tight_layout()


class LayerAligner(object):

    def __init__(self, reader, reference_aligner, channel=None, max_shift=15,
                 filter_sigma=0.0, verbose=False):
        self.reader = reader
        self.reference_aligner = reference_aligner
        if channel is None:
            channel = reference_aligner.channel
        self.channel = channel
        # Unit is micrometers.
        self.max_shift = max_shift
        self.max_shift_pixels = self.max_shift / self.metadata.pixel_size
        self.filter_sigma = filter_sigma
        self.verbose = verbose
        # FIXME Still a bit muddled here on the use of metadata positions vs.
        # corrected positions from the reference aligner. We probably want to
        # use metadata positions to find the cycle-to-cycle tile
        # correspondences, but the corrected positions for computing our
        # corrected positions.

    neighbors_graph = neighbors_graph

    def run(self):
        self.make_thumbnail()
        self.coarse_align()
        self.register_all()
        self.calculate_positions()

    def make_thumbnail(self):
        self.reader.thumbnail = thumbnail.make_thumbnail(
            self.reader, channel=self.channel
        )

    def coarse_align(self):
        self.cycle_offset = thumbnail.calculate_cycle_offset(
            self.reference_aligner.reader, self.reader
        )
        self.corrected_nominal_positions = self.metadata.positions + self.cycle_offset
        reference_positions = self.reference_aligner.metadata.positions
        dist = scipy.spatial.distance.cdist(reference_positions,
                                            self.corrected_nominal_positions)
        self.reference_idx = np.argmin(dist, 0)
        self.reference_positions = reference_positions[self.reference_idx]
        self.reference_aligner_positions = self.reference_aligner.positions[self.reference_idx]

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
        self.positions = (
            self.corrected_nominal_positions
            + self.shifts
            + self.reference_aligner_positions
            - self.reference_positions
        )
        self.constrain_positions()
        self.centers = self.positions + self.metadata.size / 2

    def constrain_positions(self):
        # Discard camera background registration which will shift target
        # positions to reference aligner positions, due to strong
        # self-correlation of the sensor dark current pattern which dominates in
        # low-signal images.
        position_diffs = np.absolute(
            self.positions - self.reference_aligner_positions
        )
        # Round the diffs to one decimal point because the subpixel shifts are
        # calculated by 10x upsampling.
        position_diffs = np.rint(position_diffs * 10) / 10
        discard = (position_diffs == 0).all(axis=1)
        # Discard any tile registration that error is infinite
        discard |= np.isinf(self.errors)
        # Take the median of registered shifts to determine the offset
        # (translation) from the reference image to this one.
        if discard.all():
            offset = 0
        else:
            offset = np.nan_to_num(np.median(self.shifts[~discard], axis=0))
        # Here we assume the fitted linear model from the reference image is
        # still appropriate, apart from the extra offset we just computed.
        predictions = self.reference_aligner.lr.predict(self.corrected_nominal_positions)
        # Discard any tile registration that's too far from the linear model,
        # replacing it with the relevant model prediction.
        distance = np.linalg.norm(self.positions - predictions - offset, axis=1)
        max_dist = self.max_shift_pixels
        extremes = distance > max_dist
        # Recalculate the mean shift, also ignoring the extreme values.
        discard |= extremes
        self.discard = discard
        if discard.all():
            self.offset = 0
        else:
            self.offset = np.nan_to_num(np.mean(self.shifts[~discard], axis=0))
        # Fill in discarded shifts from the predictions.
        self.positions[discard] = predictions[discard] + self.offset

    def register(self, t):
        """Return relative shift between images and the alignment error."""
        its, ref_img, img = self.overlap(t)
        if np.any(np.array(its.shape) == 0):
            return (0, 0), np.inf
        shift, error = utils.register(ref_img, img, self.filter_sigma)
        # We don't use padding and thus can skip the math to account for it.
        assert (its.padding == 0).all(), "Unexpected non-zero padding"
        return shift, error

    def intersection(self, t):
        corners1 = np.vstack([self.reference_positions[t],
                              self.corrected_nominal_positions[t]])
        corners2 = corners1 + self.reader.metadata.size
        its = Intersection(corners1, corners2)
        its.shape = its.shape // 32 * 32
        return its

    def overlap(self, t):
        its = self.intersection(t)
        ref_t = self.reference_idx[t]
        img1 = self.reference_aligner.reader.read(
            series=ref_t, c=self.reference_aligner.channel
        )
        img2 = self.reader.read(series=t, c=self.channel)
        ov1 = utils.crop(img1, its.offsets[0], its.shape)
        ov2 = utils.crop(img2, its.offsets[1], its.shape)
        return its, ov1, ov2

    @property
    def metadata(self):
        return self.reader.metadata

    def debug(self, t):
        shift, _ = self.register(t)
        its, o1, o2 = self.overlap(t)
        w1 = utils.whiten(o1, self.filter_sigma)
        w2 = utils.whiten(o2, self.filter_sigma)
        corr = scipy.fft.fftshift(np.abs(scipy.fft.ifft2(
            scipy.fft.fft2(w1) * scipy.fft.fft2(w2).conj()
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
        shift += origin
        plt.plot(shift[1], shift[0], 'rx')
        plt.tight_layout(0, 0, 0)


class Intersection(object):

    def __init__(self, corners1, corners2, min_size=0):
        if np.isscalar(min_size):
            min_size = np.repeat(min_size, 2)
        self._calculate(corners1, corners2, min_size)

    def _calculate(self, corners1, corners2, min_size):
        max_shape = (corners2 - corners1).max(axis=0)
        min_size = min_size.clip(1, max_shape)
        position = corners1.max(axis=0)
        initial_shape = np.floor(corners2.min(axis=0) - position).astype(int)
        clipped_shape = np.maximum(initial_shape, min_size)
        self.shape = np.ceil(clipped_shape).astype(int)
        self.padding = self.shape - initial_shape
        self.offsets = np.maximum(position - corners1 - self.padding, 0)

    def __repr__(self):
        s = 'shape: {0.shape}\npadding: {0.padding}\noffsets:\n{0.offsets}'
        return s.format(self)


class Mosaic(object):

    def __init__(
            self, aligner, shape, filename_format, channels=None,
            ffp_path=None, dfp_path=None, flip_mosaic_x=False, flip_mosaic_y=False,
            combined=False, tile_size=None, first=False, verbose=False
    ):
        self.aligner = aligner
        self.shape = tuple(shape)
        self.filename_format = filename_format
        self.channels = self._sanitize_channels(channels)
        self.flip_mosaic_x = flip_mosaic_x
        self.flip_mosaic_y = flip_mosaic_y
        self.combined = combined
        self.tile_size = tile_size
        self.first = first
        self.dtype = aligner.metadata.pixel_dtype
        self._load_correction_profiles(dfp_path, ffp_path)
        self.verbose = verbose

    def _sanitize_channels(self, channels):
        all_channels = range(self.aligner.metadata.num_channels)
        if channels is None:
            channels = all_channels
        invalid_channels = sorted(set(channels) - set(all_channels))
        if invalid_channels:
            raise ValueError("invalid channels: %s" % invalid_channels)
        return channels

    def _load_single_profile(self, path, num_channels, img_size, profile_type):
        """Load, normalize, and validate illumination profile.

        Parameters
        ----------
        path : str
            Path to the image being loaded.
        num_channels : int
            Expected number of channels in the profile image.
        img_size : tuple
            Shape of a 2D image in (row, column).
        profile_type : str
            Type of profile, only accepts 'dark' and 'flat'.

        Returns
        ----------
        ndarray
            Image as numpy array in the (channel, row, column) arrangement.
            If ``path`` is ``None``, return an array in (channel, 1, 1) shape.
            The values in the array are 0 and 1 for dark- and flat-field profile, respectively.
        """
        assert profile_type in ('dark', 'flat'), "profile_type must be either 'dark' or 'flat'."
        if path is None:
            profile_shape = (num_channels, 1, 1)
            return (
                np.zeros(profile_shape)
                    if profile_type == 'dark'
                    else np.ones(profile_shape)
            )

        expected_ndim = 2 if num_channels == 1 else 3
        profile = skimage.io.imread(path)
        if profile.ndim != expected_ndim:
            raise ValueError(
                'Expect dimensionality is {} for {}-field profile but {} has {} dimensions.'.format(
                    expected_ndim, profile_type, path, profile.ndim
                )
            )

        profile = np.atleast_3d(profile)
        # skimage.io.imread convert images with 3 and 4 channels into (Y, X, C) shape,
        # but as (C, Y, X) for images with other channel numbers. We normalize
        # image-shape to (C, Y, X) regardless of the number of channels in the image.
        if num_channels in (1, 3, 4):
            profile = np.moveaxis(profile, 2, 0)
        if profile.shape != (num_channels,) + img_size:
            raise ValueError(
                '{}-field profile shape {} does not match target image shape {}.'.format(
                    profile_type.capitalize(), profile.shape, img_size
                )
            )

        return profile

    def _load_correction_profiles(self, dfp_path, ffp_path):
        if dfp_path is None and ffp_path is None:
            self.do_correction = False
        else:
            num_channels = self.aligner.metadata.num_channels
            img_size = tuple(self.aligner.metadata.size)
            self.dfp = self._load_single_profile(dfp_path, num_channels, img_size, 'dark')
            self.ffp = self._load_single_profile(ffp_path, num_channels, img_size, 'flat')

            # FIXME This assumes integer dtypes. Do we need to support floats?
            self.dfp /= np.iinfo(self.dtype).max
            self.do_correction = True

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
        for ci, channel in enumerate(self.channels):
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
                func = utils.pastefunc_blend if not debug else np.add
                utils.paste(mosaic_image, tile_image, position, func=func)
            if debug:
                np.clip(mosaic_image, 0, 1, out=mosaic_image)
                w = int(1e6)
                mi_flat = mosaic_image.reshape(-1, 3)
                for p in np.arange(0, mi_flat.shape[0], w, dtype=int):
                    mi_flat[p:p+w] = skimage.exposure.adjust_gamma(
                        mi_flat[p:p+w], 1/2.2
                    )
            if self.flip_mosaic_x:
                mosaic_image = np.fliplr(mosaic_image)
            if self.flip_mosaic_y:
                mosaic_image = np.flipud(mosaic_image)
            if self.verbose:
                print()
            if mode == 'write':
                filename = self.filename_format.format(channel=channel)
                kwargs = {}
                if self.combined:
                    kwargs['bigtiff'] = True
                    # FIXME Propagate this from input files (esp. RGB).
                    kwargs['photometric'] = 'minisblack'
                    resolution = np.round(10000 / self.aligner.reader.metadata.pixel_size)
                    # FIXME Switch to "CENTIMETER" once we use tifffile directly.
                    kwargs['resolution'] = (resolution, resolution, 'cm')
                    kwargs['metadata'] = None
                    if self.first and ci == 0:
                        # Set description to a short placeholder that will fit
                        # within the IFD. We'll check for this string later.
                        kwargs['description'] = '!!xml!!'
                        kwargs['software'] = (
                            'Ashlar v{} (Glencoe/Faas pyramid output)'
                            .format(_version)
                        )
                    else:
                        # Overwite if first channel of first cycle.
                        kwargs['append'] = True
                if self.tile_size:
                    kwargs['tile'] = (self.tile_size, self.tile_size)
                if self.verbose:
                    print("        writing to %s" % filename)
                utils.imsave(filename, mosaic_image, **kwargs)
            elif mode == 'return':
                all_images.append(mosaic_image)
        if mode == 'return':
            return all_images

    def correct_illumination(self, img, channel):
        if self.do_correction:
            img = skimage.util.img_as_float(img, force_copy=True)
            img -= self.dfp[channel, ...]
            img /= self.ffp[channel, ...]
            img.clip(0, 1, out=img)
        return img


def build_pyramid(
        path, num_channels, shape, dtype, pixel_size, tile_size, verbose=False
):
    max_level = 0
    shapes = [shape]
    while any(s > tile_size for s in shape):
        prev_level = max_level
        max_level += 1
        if verbose:
            print("    Level %d:" % max_level)
        for i in range(num_channels):
            if verbose:
                sys.stdout.write('\r        processing channel %d/%d'
                                 % (i + 1, num_channels))
                sys.stdout.flush()
            img = skimage.io.imread(path, series=prev_level, key=i)
            img = skimage.transform.pyramid_reduce(img, multichannel=False)
            img = skimage.util.dtype.convert(img, dtype)
            utils.imsave(
                path, img, bigtiff=True, metadata=None, append=True,
                tile=(tile_size, tile_size), photometric='minisblack'
            )
        shapes.append(img.shape)
        if verbose:
            print()
        shape = img.shape
    # Now that we have the number and dimensions of all levels, we can generate
    # the corresponding OME-XML and patch it into the Image Description tag of
    # the first IFD.
    filename = pathlib.Path(path).name
    img_uuid = uuid.uuid4().urn
    ome_dtype = BioformatsMetadata._ome_dtypes[dtype]
    ifd = 0
    xml = io.StringIO()
    xml.write(u'<?xml version="1.0" encoding="UTF-8"?>')
    xml.write(
        (u'<OME xmlns="http://www.openmicroscopy.org/Schemas/OME/2016-06"'
         ' xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"'
         ' UUID="{uuid}"'
         ' xsi:schemaLocation="http://www.openmicroscopy.org/Schemas/OME/2016-06'
         ' http://www.openmicroscopy.org/Schemas/OME/2016-06/ome.xsd">')
        .format(uuid=img_uuid)
    )
    for level in range(max_level + 1):
        shape = shapes[level]
        if level == 0:
            psize_xml = (
                u'PhysicalSizeX="{0}" PhysicalSizeXUnit="\u00b5m"'
                u' PhysicalSizeY="{0}" PhysicalSizeYUnit="\u00b5m"'
                .format(pixel_size)
            )
        else:
            psize_xml = u''
        xml.write(u'<Image ID="Image:{}">'.format(level))
        xml.write(
            (u'<Pixels BigEndian="false" DimensionOrder="XYZCT"'
             ' ID="Pixels:{level}" {psize_xml} SizeC="{num_channels}" SizeT="1"'
             ' SizeX="{sizex}" SizeY="{sizey}" SizeZ="1" Type="{ome_dtype}">')
            .format(
                level=level, psize_xml=psize_xml, num_channels=num_channels,
                sizex=shape[1], sizey=shape[0], ome_dtype=ome_dtype
            )
        )
        for channel in range(num_channels):
            xml.write(
                (u'<Channel ID="Channel:{level}:{channel}"'
                 + (u' Name="Channel {channel}"' if level == 0 else u'')
                 + u' SamplesPerPixel="1"><LightPath/></Channel>')
                .format(level=level, channel=channel)
            )
        for channel in range(num_channels):
            xml.write(
                (u'<TiffData FirstC="{channel}" FirstT="0" FirstZ="0"'
                 ' IFD="{ifd}" PlaneCount="1">'
                 '<UUID FileName="{filename}">{uuid}</UUID>'
                 '</TiffData>')
                .format(
                    channel=channel, ifd=ifd, filename=filename, uuid=img_uuid
                )
            )
            ifd += 1
        if level == 0:
            for channel in range(num_channels):
                xml.write(
                    u'<Plane TheC="{channel}" TheT="0" TheZ="0"/>'
                    .format(channel=channel)
                )
        xml.write(u'</Pixels>')
        xml.write(u'</Image>')
    xml.write(u'</OME>')
    xml_bytes = xml.getvalue().encode('utf-8') + b'\x00'
    # Append the XML and patch up the Image Description tag in the first IFD.
    with open(path, 'rb+') as f:
        f.seek(0, io.SEEK_END)
        xml_offset = f.tell()
        f.write(xml_bytes)
        f.seek(0)
        ifd_block = f.read(500)
        match = re.search(b'!!xml!!\x00', ifd_block)
        if match is None:
            raise RuntimeError("Did not find placeholder string in IFD")
        f.seek(match.start() - 8)
        f.write(struct.pack('<Q', len(xml_bytes)))
        f.write(struct.pack('<Q', xml_offset))


class DataWarning(UserWarning):
    """Warnings about the content of user-provided image data."""
    pass


def warn_data(message):
    warnings.warn(message, DataWarning)


def plot_edge_shifts(aligner, img=None, bounds=True, im_kwargs=None):
    if im_kwargs is None:
        im_kwargs = {}
    fig = plt.figure()
    ax = plt.gca()
    draw_mosaic_image(ax, aligner, img, **im_kwargs)
    h, w = aligner.reader.metadata.size
    if bounds:
        # Bounding boxes denoting new tile positions.
        for xy in np.fliplr(aligner.positions):
            rect = mpatches.Rectangle(xy, w, h, color='black', fill=False,
                                      lw=0.5)
            ax.add_patch(rect)
    # Compute per-edge relative shifts from tile positions.
    edges = np.array(list(aligner.spanning_tree.edges))
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


def plot_edge_quality(
    aligner, img=None, show_tree=True, pos='metadata', im_kwargs=None, nx_kwargs=None
):
    if pos == 'metadata':
        centers = aligner.metadata.centers - aligner.metadata.origin
    elif pos == 'aligner':
        centers = aligner.centers
    else:
        raise ValueError("pos must be either 'metadata' or 'aligner'")
    if im_kwargs is None:
        im_kwargs = {}
    if nx_kwargs is None:
        nx_kwargs = {}
    final_nx_kwargs = dict(width=2, node_size=100, font_size=6)
    final_nx_kwargs.update(nx_kwargs)
    if show_tree:
        nrows, ncols = 1, 2
        if aligner.mosaic_shape[1] * 2 / aligner.mosaic_shape[0] > 2 * 4 / 3:
            nrows, ncols = ncols, nrows
    else:
        nrows, ncols = 1, 1
    fig = plt.figure()
    ax = plt.subplot(nrows, ncols, 1)
    draw_mosaic_image(ax, aligner, img, **im_kwargs)
    error = np.array([aligner._cache[tuple(sorted(e))][1]
                      for e in aligner.neighbors_graph.edges])
    # Manually center and scale data to 0-1, except infinity which is set to -1.
    # This lets us use the purple-green diverging color map to color the graph
    # edges and cause the "infinity" edges to disappear into the background
    # (which is itself purple).
    infs = error == np.inf
    error[infs] = -1
    if not infs.all():
        error_f = error[~infs]
        emin = np.min(error_f)
        emax = np.max(error_f)
        if emin == emax:
            # Always true when there's only one edge. Otherwise it's unlikely
            # but theoretically possible.
            erange = 1
        else:
            erange = emax - emin
        error[~infs] = (error_f - emin) / erange
    # Neighbor graph colored by edge alignment quality (brighter = better).
    nx.draw(
        aligner.neighbors_graph, ax=ax, with_labels=True,
        pos=np.fliplr(centers), edge_color=error, edge_vmin=-1, edge_vmax=1,
        edge_cmap=plt.get_cmap('PRGn'), **final_nx_kwargs
    )
    if show_tree:
        ax = plt.subplot(nrows, ncols, 2)
        draw_mosaic_image(ax, aligner, img, **im_kwargs)
        # Spanning tree with nodes at original tile positions.
        nx.draw(
            aligner.spanning_tree, ax=ax, with_labels=True,
            pos=np.fliplr(centers), edge_color='royalblue',
            **final_nx_kwargs
        )
    fig.set_facecolor('black')


def plot_edge_scatter(aligner, annotate=True):
    import seaborn as sns
    xdata = aligner.all_errors
    ydata = np.clip(
        [np.linalg.norm(v[0]) for v in aligner._cache.values()], 0.01, np.inf
    )
    pdata = np.clip(aligner.errors_negative_sampled, 0, 10)
    g = sns.JointGrid(xdata, ydata)
    g.plot_joint(sns.scatterplot, alpha=0.5)
    _, xbins = np.histogram(np.hstack([xdata, pdata]), bins=40)
    sns.distplot(
        xdata, ax=g.ax_marg_x, kde=False, bins=xbins, norm_hist=True
    )
    sns.distplot(
        pdata, ax=g.ax_marg_x, kde=False, bins=xbins, norm_hist=True,
        hist_kws=dict(histtype='step')
    )
    g.ax_joint.axvline(aligner.max_error, c='k', ls=':')
    g.ax_joint.axhline(aligner.max_shift_pixels, c='k', ls=':')
    g.ax_joint.set_yscale('log')
    g.set_axis_labels('error', 'shift')
    if annotate:
        for pair, x, y in zip(aligner.neighbors_graph.edges, xdata, ydata):
            plt.annotate(str(pair), (x, y), alpha=0.1)
    plt.tight_layout()


def plot_layer_shifts(aligner, img=None, im_kwargs=None):
    if im_kwargs is None:
        im_kwargs = {}
    fig = plt.figure()
    ax = plt.gca()
    draw_mosaic_image(ax, aligner, img, **im_kwargs)
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


def plot_layer_quality(
    aligner, img=None, scale=1.0, artist='patches', annotate=True, im_kwargs=None
):
    if im_kwargs is None:
        im_kwargs = {}
    fig = plt.figure()
    ax = plt.gca()
    draw_mosaic_image(ax, aligner, img, **im_kwargs)

    h, w = aligner.metadata.size
    positions, centers, shifts = aligner.positions, aligner.centers, aligner.shifts

    if scale != 1.0:
        h, w, positions, centers, shifts = [
            scale * i for i in [h, w, positions, centers, shifts]
        ]

    # Bounding boxes denoting new tile positions.
    color_index = skimage.exposure.rescale_intensity(
        aligner.errors, out_range=np.uint8
    ).astype(np.uint8)
    color_map = mcm.magma_r
    for xy, c_idx in zip(np.fliplr(positions), color_index):
        rect = mpatches.Rectangle(
            xy, w, h, color=color_map(c_idx), fill=False, lw=0.5
        )
        ax.add_patch(rect)
    
    # Annotate tile numbering.
    if annotate:
        for idx, (x, y) in enumerate(np.fliplr(positions)):
            text = plt.annotate(str(idx), (x+0.1*w, y+0.9*h), alpha=0.7)
            # Add outline to text for better contrast in different background color.
            text_outline = mpatheffects.Stroke(linewidth=1, foreground='#AAA')
            text.set_path_effects(
                [text_outline, mpatheffects.Normal()]
            )

    if artist == 'quiver':
        ax.quiver(
            *centers.T[::-1], *shifts.T[::-1], aligner.discard,
            units='dots', width=2, scale=1, scale_units='xy', angles='xy',
            cmap='Greys'
        )
    if artist == 'patches':
        for xy, dxy, is_discarded in zip(
            np.fliplr(centers), np.fliplr(shifts), aligner.discard
        ):
            arrow = mpatches.FancyArrowPatch(
                xy, np.array(xy) + np.array(dxy), 
                arrowstyle='->', color='0' if is_discarded else '1',
                mutation_scale=8,
                )
            ax.add_patch(arrow)
    ax.axis('off')


def draw_mosaic_image(ax, aligner, img, **kwargs):
    if img is None:
        img = [[0]]
    h, w = aligner.mosaic_shape
    ax.imshow(img, extent=(-0.5, w-0.5, h-0.5, -0.5), **kwargs)
