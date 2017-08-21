from __future__ import division
import sys
import collections
import queue
import numpy as np
import pandas as pd
import skimage.io
import reg


TileStatistics = collections.namedtuple(
    'TileStatistics',
    'scan tile x_original y_original x y shift_x shift_y error'
)


filepaths = sys.argv[1:]
assert len(filepaths) > 0
assert all(p.endswith('.rcpnl') for p in filepaths)

reader0 = reg.Reader(filepaths[0])
metadata = reader0.metadata
aligner = reg.EdgeAligner(reader0, verbose=True)
aligner.run()

mshape = np.ceil((aligner.positions + metadata.size).max(axis=0)).astype(int)
mosaic0 = np.zeros(mshape, dtype=np.uint16)
for i, npos in enumerate(aligner.positions):
    sys.stdout.write("\rScan 0: merging %d/%d" % (i + 1, metadata.num_images))
    sys.stdout.flush()
    reg.paste(mosaic0, reader0.read(c=0, series=i), npos)
print
skimage.io.imsave('scan_0_0.tif', mosaic0)

for s, filepath in enumerate(filepaths[1:], 1):
    reader = reg.Reader(filepath)
    metadata = reader.metadata
    mosaic = np.zeros(mshape, dtype=np.uint16)
    aligner_l = reg.LayerAligner(reader, reader0, aligner)
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
