from __future__ import division
import sys
import time
import collections
import queue
import bioformats
import javabridge
import numpy as np
import pandas as pd
try:
    import pyfftw
    np.fft = pyfftw.interfaces.numpy_fft
except ImportError:
    print "pyfftw not found, falling back to numpy.fft"
import scipy.ndimage
import networkx as nx
import skimage.feature
import skimage.io
from skimage.restoration.uft import laplacian
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import modest_image

import reg


def graph_from_positions(positions, max_distance):
    pdist = scipy.spatial.distance.pdist(positions, metric='cityblock')
    sp = scipy.spatial.distance.squareform(pdist)
    edges = zip(*np.nonzero((sp > 0) & (sp < max_distance)))
    graph = nx.from_edgelist(edges)
    return graph


TileStatistics = collections.namedtuple(
    'TileStatistics',
    'scan tile x_original y_original x y shift_x shift_y error'
)


filepaths = sys.argv[1:]
assert len(filepaths) > 0
assert all(p.endswith('.rcpnl') for p in filepaths)

reader0 = reg.Reader(filepaths[0])
metadata = reader0.metadata
aligner = reg.EdgeAligner(reader0)

positions = metadata.positions - metadata.origin
centers = metadata.centers - metadata.origin
neighbor_max_distance = metadata.size.max()
graph = graph_from_positions(positions, neighbor_max_distance)

num_edges = graph.size()
for i, (t1, t2) in enumerate(graph.edges_iter(), 1):
    sys.stdout.write('\raligning: %d/%d' % (i, num_edges))
    sys.stdout.flush()
    aligner.register(t1, t2)
print

lg = nx.line_graph(graph)
spanning_tree = nx.Graph()
fringe = queue.PriorityQueue()
start_edge = sorted(aligner._cache.keys(), key=lambda k: aligner._cache[k][1])[0]
shifts = {start_edge[0]: np.array([0, 0])}
fringe.put((aligner.register(*start_edge)[1], start_edge))
while not fringe.empty():
    _, edge = fringe.get()
    if edge[0] in spanning_tree and edge[1] in spanning_tree:
        continue
    spanning_tree.add_edge(*edge)
    source, dest = edge
    if source not in shifts:
        source, dest = dest, source
    shifts[dest] = shifts[source] + aligner.register(source, dest)[0]
    for next_edge in set(lg.neighbors(edge)):
        fringe.put((aligner.register(*next_edge)[1], next_edge))

shifts = np.array(zip(*sorted(shifts.items()))[1])
new_positions = positions + shifts
new_positions -= new_positions.min(axis=0)
new_centers = new_positions + metadata.size / 2
mshape = np.ceil((new_positions + metadata.size).max(axis=0)).astype(int)
mosaic0 = np.zeros(mshape, dtype=np.uint16)
for i, npos in enumerate(new_positions):
    sys.stdout.write("\rScan 0: merging %d/%d" % (i + 1, metadata.num_images))
    sys.stdout.flush()
    reg.paste(mosaic0, reader0.read(c=0, series=i), npos)
print
skimage.io.imsave('scan_0_0.tif', mosaic0)

for s, filepath in enumerate(filepaths[1:], 1):
    reader = reg.Reader(filepath)
    metadata = reader.metadata
    mosaic = np.zeros(mshape, dtype=np.uint16)
    aligner_l = reg.LayerAligner(reader, reader0, new_positions, shifts)
    for i in range(metadata.num_images):
        sys.stdout.write("\rScan %d: merging %d/%d"
                         % (s, i + 1, metadata.num_images))
        sys.stdout.flush()
        npos, error = aligner_l.register(i)
        reg.paste(mosaic, reader.read(c=0, series=i), npos)
    print
    skimage.io.imsave('scan_%d_0.tif' % s, mosaic)


try:
    __IPYTHON__
except:
    reg._deinit_bioformats()


"""

plt.figure()
ax = plt.gca()
modest_image.imshow(ax, mosaic)
h, w = metadata.size
for xy in np.fliplr(new_positions):
    ax.add_patch(mpatches.Rectangle(xy, w, h, color='black', fill=False, lw=0.5))
nx.draw(spanning_tree, ax=ax, pos=np.fliplr(new_centers), with_labels=True,
    edge_color=np.sum(np.array([aligner._cache[tuple(sorted(e))][0] for e in spanning_tree.edges()]) ** 2, axis=1) ** 0.5,
    edge_cmap=plt.get_cmap('Blues_r'), width=2, node_size=100, font_size=6)


nrows, ncols = 1, 2
if mosaic.shape[1] / mosaic.shape[0] / 2 < 4 / 3:
    nrows, ncols = ncols, nrows
plt.figure()
ax = plt.subplot(nrows, ncols,1)
modest_image.imshow(ax, mosaic)
nx.draw(
    graph, ax=ax, pos=np.fliplr(centers), with_labels=True,
    edge_color=[aligner._cache[tuple(sorted(e))][1] for e in graph.edges()],
    edge_cmap=plt.get_cmap('hot_r'), width=2, node_size=100, font_size=6
)
ax = plt.subplot(nrows, ncols, 2)
modest_image.imshow(ax, mosaic)
nx.draw(
    spanning_tree, ax=ax, pos=np.fliplr(centers), with_labels=True,
    edge_color='royalblue', width=2, node_size=100, font_size=6
)


"""
