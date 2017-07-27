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

neighbor_max_distance = metadata.sizes.max(axis=1).min()
graph = graph_from_positions(metadata.positions, neighbor_max_distance)

fringe = queue.PriorityQueue()
fringe.put((0, 0))
costs = {0: 0}
came_from = {}
shifts = {}

while not fringe.empty():
    print 'scored: %d  fringe: %d' % (len(costs), fringe.qsize())
    _, cur_tile = fringe.get()
    for next_tile in graph.neighbors(cur_tile):
        shift, error = aligner.register(cur_tile, next_tile)
        new_cost = costs[cur_tile] + error
        if next_tile not in costs or new_cost < costs[next_tile]:
            costs[next_tile] = new_cost
            fringe.put((new_cost, next_tile))
            came_from[next_tile] = cur_tile
            shifts[next_tile] = shift
spanning_tree = nx.from_edgelist(came_from.items())


try:
    __IPYTHON__
except:
    reg._deinit_bioformats()
