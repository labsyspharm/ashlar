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


filepath = sys.argv[1]
assert filepath.endswith('.rcpnl')

reader = reg.Reader(filepath)
metadata = reader.metadata
aligner = reg.Aligner(reader)

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


new_positions = positions + np.array(zip(*sorted(shifts.items()))[1])
new_positions -= new_positions.min(axis=0)
new_centers = new_positions + metadata.size / 2
mshape = np.ceil((new_positions + metadata.size).max(axis=0)).astype(int)
mosaic = np.zeros(mshape, dtype=np.uint16)
for i, npos in enumerate(new_positions):
    sys.stdout.write("\rLoading %d/%d" % (i + 1, metadata.num_images))
    sys.stdout.flush()
    reg.paste(mosaic, reader.read(c=0, series=i), npos)
print


try:
    __IPYTHON__
except:
    reg._deinit_bioformats()


"""

plt.figure()
ax = plt.gca()
modest_image.imshow(ax, mosaic)
nx.draw(spanning_tree, ax=ax, pos=np.fliplr(new_centers), with_labels=True,
    edge_color='royalblue', width=2, node_size=100, font_size=6)


plt.figure()
ax = plt.subplot(121)
modest_image.imshow(ax, mosaic)
nx.draw(
    graph, ax=ax, pos=np.fliplr(centers), with_labels=True,
    edge_color=[aligner._cache[tuple(sorted(e))][1] for e in graph.edges()],
    edge_cmap=plt.get_cmap('hot_r'), width=2, node_size=100, font_size=6
)
ax = plt.subplot(122)
modest_image.imshow(ax, mosaic)
nx.draw(
    spanning_tree, ax=ax, pos=np.fliplr(centers), with_labels=True,
    edge_color='royalblue', width=2, node_size=100, font_size=6
)


"""
