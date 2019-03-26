from __future__ import division, print_function
import sys
import warnings
import re
import itertools
import xml.etree.ElementTree
import io
import uuid
import struct
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
import skimage.transform
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
from . import __version__ as _version

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
ChannelSeparator = jnius.autoclass('loci.formats.ChannelSeparator')
UNITS = jnius.autoclass('ome.units.UNITS')

# Another workaround for pyjnius #300 (see below).
DebugTools.enableLogging(JString("ERROR"))


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


class BioformatsMetadata(PlateMetadata):

    _pixel_dtypes = {
        'uint8': np.uint8,
        'uint16': np.uint16,
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
        # FIXME Workaround for pyjnius #300 until there is a new release.
        # Passing a python string directly here corrupts the value under Python
        # 3, but explicitly converting it into a Java string works.
        path_jstring = JString(self.path)
        self._reader.setId(path_jstring)

        xml_content = metadata.dumpXML()
        self._metadata = metadata
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
            v = method(0).value(UNITS.MICROMETER).doubleValue()
            values.append(v)
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
        values = []
        for dim in ('Y', 'X'):
            method = getattr(self._metadata, 'getPlanePosition%s' % dim)
            # FIXME verify all planes have the same X,Y position.
            v_units = method(i, 0)
            v = v_units.value(UNITS.MICROMETER)
            if v is None:
                # Conversion failed, which usually happens when the unit is
                # "reference frame". Proceed as if it's actually microns but
                # emit a warning.
                warnings.warn(
                    "Stage coordinates' measurement unit is undefined;"
                    " assuming micrometers."
                )
                v = v_units.value()
            values.append(v.doubleValue())
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
        aligner._neighbors_graph = graph
    return aligner._neighbors_graph


class EdgeAligner(object):

    def __init__(
        self, reader, channel=0, max_shift=15, false_positive_ratio=0.01,
        filter_sigma=0.0, verbose=False
    ):
        self.reader = reader
        self.channel = channel
        self.verbose = verbose
        # Unit is micrometers.
        self.max_shift = max_shift
        self.max_shift_pixels = self.max_shift / self.metadata.pixel_size
        self.false_positive_ratio = false_positive_ratio
        self.filter_sigma = filter_sigma
        self._cache = {}

    neighbors_graph = neighbors_graph

    def run(self):
        self.check_overlaps()
        self.compute_threshold()
        self.register_all()
        self.build_spanning_tree()
        self.calculate_positions()
        self.fit_model()

    def check_overlaps(self):
        # This might be better addressed by removing the +1 from the
        # neighbors_graph max_distance calculation and ensuring the graph is
        # fully connected.
        pos = self.metadata.positions
        overlaps = np.array([
            self.metadata.size - abs(pos[t1] - pos[t2])
            for t1, t2 in self.neighbors_graph.edges
        ])
        failures = np.any(overlaps < 1, axis=1)
        if all(failures):
            warnings.warn("No tiles overlap, attempting alignment anyway.")
        elif any(failures):
            warnings.warn("Some neighboring tiles have zero overlap.")

    def compute_threshold(self):
        # Compute error threshold for rejecting aligments. We generate a
        # distribution of error scores for many known non-overlapping image
        # regions and take a certain percentile as the maximum allowable error.
        # The percentile becomes our accepted false-positive ratio.
        edges = self.neighbors_graph.edges()
        min_size = np.repeat(self.max_shift_pixels * 2, 2)
        widths = np.array([
            self.intersection(t1, t2, min_size).shape.min()
            for t1, t2 in edges
        ])
        w = widths.max()
        max_offset = self.metadata.size[0] - w
        n = 1000
        pairs = np.empty((n, 2), dtype=int)
        offsets = np.empty((n, 2), dtype=int)
        # Generate n random image pairs and strip offsets for alignment. Even
        # with a small number of tiles, we can easily get 1000 non-repeating
        # strips due to also varying the offset.
        i = 0
        for i in range(n):
            # Ensure pair is two different tiles and tiles are not neighbors.
            # This is more conservative than necessary -- we could admit
            # neighbors as long as the chosen offsets don't correspond to the
            # actual overlap region.
            while True:
                t1, t2 = np.random.randint(self.metadata.num_images, size=2)
                if t1 != t2 and (t1, t2) not in edges:
                    break
            pairs[i] = t1, t2
            offsets[i] = np.random.randint(max_offset, size=2)
            i += 1
        errors = np.empty(n)
        for i, ((t1, t2), (offset1, offset2)) in enumerate(zip(pairs, offsets)):
            if self.verbose and (i % 10 == 9 or i == n - 1):
                sys.stdout.write(
                    '\r    quantifying alignment error %d/%d' % (i + 1, n)
                )
                sys.stdout.flush()
            img1 = self.reader.read(t1, self.channel)[offset1:offset1+w, :]
            img2 = self.reader.read(t2, self.channel)[offset2:offset2+w, :]
            _, errors[i] = register(img1, img2, self.filter_sigma)
        if self.verbose:
            print()
        self.errors_negative_sampled = errors
        self.max_error = np.percentile(errors, self.false_positive_ratio * 100)

    def register_all(self):
        n = self.neighbors_graph.size()
        for i, (t1, t2) in enumerate(self.neighbors_graph.edges(), 1):
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
        if shifts:
            self.shifts = np.array([s for _, s in sorted(shifts.items())])
            self.positions = self.metadata.positions + self.shifts
        else:
            # TODO: fill in shifts and positions with 0x2 arrays
            raise NotImplementedError("No images")

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
            max_its_size = np.repeat(self.max_shift_pixels * 2, 2)
            shift, error = self._register(t1, t2, max_its_size)
            self._cache[key] = (shift, error)
        if t1 > t2:
            shift = -shift
        # Return copy of shift to prevent corruption of cached values.
        return shift.copy(), error

    def _register(self, t1, t2, min_size):
        its, img1, img2 = self.overlap(t1, t2, min_size)
        # Account for padding, flipping the sign depending on the direction
        # between the tiles.
        p1, p2 = self.metadata.positions[[t1, t2]]
        sx = 1 if p1[1] >= p2[1] else -1
        sy = 1 if p1[0] >= p2[0] else -1
        padding = its.padding * [sy, sx]
        shift, error = register(img1, img2, self.filter_sigma)
        shift += padding
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
        w1 = whiten(o1, self.filter_sigma)
        w2 = whiten(o2, self.filter_sigma)
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
        self.tile_positions = self.metadata.positions - reference_aligner.origin
        reference_positions = reference_aligner.positions
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
        self.positions = self.tile_positions + self.shifts
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
        shift, error = register(ref_img, img, self.filter_sigma)
        # We don't use padding and thus can skip the math to account for it.
        assert (its.padding == 0).all(), "Unexpected non-zero padding"
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
        img1 = self.reference_aligner.reader.read(
            series=ref_t, c=self.reference_aligner.channel
        )
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
        w1 = whiten(o1, self.filter_sigma)
        w2 = whiten(o2, self.filter_sigma)
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
        shift += origin
        plt.plot(shift[1], shift[0], 'rx')
        plt.tight_layout(0, 0, 0)


class Intersection(object):

    def __init__(self, corners1, corners2, min_size=None):
        if min_size is None:
            min_size = np.zeros(2)
        self._calculate(corners1, corners2, min_size)

    def _calculate(self, corners1, corners2, min_size):
        max_shape = (corners2 - corners1).max(axis=0)
        min_size = min_size.clip(1, max_shape)
        position = corners1.max(axis=0)
        initial_shape = np.ceil(corners2.min(axis=0) - position).astype(int)
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
            ffp_path=None, dfp_path=None, combined=False, tile_size=None,
            first=False, verbose=False
    ):
        self.aligner = aligner
        self.shape = tuple(shape)
        self.filename_format = filename_format
        self.channels = self._sanitize_channels(channels)
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

    def _load_correction_profiles(self, dfp_path, ffp_path):
        if dfp_path or ffp_path:
            c = max(self.channels) + 1
            self.dfp = np.atleast_3d(
                skimage.io.imread(dfp_path) if dfp_path else np.zeros(c)
            )
            self.ffp = np.atleast_3d(
                skimage.io.imread(ffp_path) if ffp_path else np.ones(c)
            )
            # FIXME This assumes integer dtypes. Do we need to support floats?
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
                kwargs = {}
                if self.combined:
                    kwargs['bigtiff'] = True
                    # FIXME Propagate this from input files (esp. RGB).
                    kwargs['photometric'] = 'minisblack'
                    resolution = np.round(10000 / self.aligner.reader.metadata.pixel_size)
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
                with warnings.catch_warnings():
                    warnings.filterwarnings(
                        'ignore', r'.* is a low contrast image', UserWarning,
                        '^skimage\.io'
                    )
                    skimage.io.imsave(filename, mosaic_image, **kwargs)
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
            with warnings.catch_warnings():
                warnings.filterwarnings(
                    'ignore', 'Possible precision loss', UserWarning,
                    '^skimage\.util\.dtype'
                )
                warnings.filterwarnings(
                    'ignore', '.* is a low contrast image', UserWarning,
                    '^skimage\.io'
                )
                img = skimage.util.dtype.convert(img, dtype)
                skimage.io.imsave(
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


def fft2(img):
    return pyfftw.builders.fft2(img, planner_effort='FFTW_ESTIMATE',
                                avoid_copy=True, auto_align_input=True,
                                auto_contiguous=True)()


# Pre-calculate the Laplacian operator kernel. We'll always be using 2D images.
_laplace_kernel = skimage.restoration.uft.laplacian(2, (3, 3))[1]

def whiten(img, sigma):
    # Copied from skimage.filters.edges, with explicit aligned output from
    # convolve. Also the mask option was dropped.
    img = skimage.img_as_float(img)
    output = pyfftw.empty_aligned(img.shape, 'complex64')
    output.imag[:] = 0
    if sigma == 0:
        scipy.ndimage.convolve(img, _laplace_kernel, output.real)
    else:
        scipy.ndimage.gaussian_laplace(img, sigma, output=output.real)
    return output

    # Other possible whitening functions:
    #img = skimage.filters.roberts(img)
    #img = skimage.filters.scharr(img)
    #img = skimage.filters.sobel(img)
    #img = np.log(img)
    #img = img - scipy.ndimage.filters.gaussian_filter(img, 2) + 0.5


def register(img1, img2, sigma):
    img1w = whiten(img1, sigma)
    img2w = whiten(img2, sigma)
    img1_f = fft2(img1w)
    img2_f = fft2(img2w)
    img1w = img1w.real
    img2w = img2w.real
    shift, _, _ = skimage.feature.register_translation(
        img1_f, img2_f, 10, 'fourier'
    )
    # At this point we may have a shift in the wrong quadrant since the FFT
    # assumes the signal is periodic. We test all four possibilities and return
    # the shift that gives the highest direct correlation (sum of products).
    shape = np.array(img1.shape)
    shift_pos = (shift + shape) % shape
    shift_neg = shift_pos - shape
    shifts = list(itertools.product(*zip(shift_pos, shift_neg)))
    correlations = [
        np.sum(img1w * scipy.ndimage.shift(img2w, s))
        for s in shifts
    ]
    idx = np.argmax(correlations)
    shift = shifts[idx]
    correlation = correlations[idx]
    total_amplitude = np.linalg.norm(img1w) * np.linalg.norm(img2w)
    if correlation > 0 and total_amplitude > 0:
        error = -np.log(correlation / total_amplitude)
    else:
        error = np.inf
    return shift, error


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
    target_slice = target[yi:yi+img.shape[0], xi:xi+img.shape[1]]
    img = crop_like(img, target_slice)
    if img.ndim == 2:
        img = scipy.ndimage.shift(img, pos_f)
    else:
        for c in range(img.shape[2]):
            img[...,c] = scipy.ndimage.shift(img[...,c], pos_f)
    # For any axis where there is a non-zero subpixel shift, crop out the last
    # row or column of pixels on the "losing" side. These pixels will be darker
    # than normal and will introduce artifacts in most blending modes.
    y1 = None if pos_f[0] <= 0 else 1
    y2 = None if pos_f[0] >= 0 else -1
    x1 = None if pos_f[1] <= 0 else 1
    x2 = None if pos_f[1] >= 0 else -1
    img = img[y1:y2, x1:x2]
    target_slice = target_slice[y1:y2, x1:x2]
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


def plot_edge_quality(
    aligner, img=None, show_tree=True, pos='metadata', use_mi=True,
    nx_kwargs=None
):
    if pos == 'metadata':
        centers = aligner.metadata.centers - aligner.metadata.origin
    elif pos == 'aligner':
        centers = aligner.centers
    else:
        raise ValueError("pos must be either 'metadata' or 'aligner'")
    if nx_kwargs is None:
        nx_kwargs = {}
    final_nx_kwargs = dict(width=2, node_size=100, font_size=6)
    final_nx_kwargs.update(nx_kwargs)
    if show_tree:
        nrows, ncols = 1, 2
        if aligner.mosaic_shape[1] * 2 / aligner.mosaic_shape[0] > 4 / 3:
            nrows, ncols = ncols, nrows
    else:
        nrows, ncols = 1, 1
    fig = plt.figure()
    ax = plt.subplot(nrows, ncols, 1)
    draw_mosaic_image(ax, aligner, img, use_mi)
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
        edge_cmap=plt.get_cmap('PRGn'), **final_nx_kwargs
    )
    if show_tree:
        ax = plt.subplot(nrows, ncols, 2)
        draw_mosaic_image(ax, aligner, img, use_mi)
        # Spanning tree with nodes at original tile positions.
        nx.draw(
            aligner.spanning_tree, ax=ax, with_labels=True,
            pos=np.fliplr(centers), edge_color='royalblue',
            **final_nx_kwargs
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


def draw_mosaic_image(ax, aligner, img, use_mi=True):
    if use_mi and img is not None:
        if sys.version_info[0] != 2:
            warnings.warn('ModestImage module (required for image display)'
                          ' is only compatible with Python 2')
            img = None
        elif modest_image is None:
            warnings.warn('Please install ModestImage for image display')
            img = None
    if img is not None:
        if use_mi:
            modest_image.imshow(ax, img)
        else:
            ax.imshow(img)
    else:
        h, w = aligner.mosaic_shape
        # Draw a single-pixel image in the lowest color in the colormap,
        # stretched across the same extent that the real image would be.
        # This makes the graph edge colors visible even if there's no image.
        ax.imshow([[0]], extent=(-0.5, w-0.5, h-0.5, -0.5))
