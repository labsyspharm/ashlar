import sys
try:
    import pathlib2 as pathlib
except:
    pass
import numpy as np
import skimage.io
import reg


filepaths = sys.argv[1:]
assert len(filepaths) > 0
if len(filepaths) == 1:
    path = pathlib.Path(filepaths[0])
    if path.is_dir():
        filepaths = sorted(str(p) for p in path.glob('*rcpnl'))
assert all(p.endswith('.rcpnl') for p in filepaths)

print 'Scan 0:'
print '    reading %s' % filepaths[0]
reader0 = reg.Reader(filepaths[0])
metadata = reader0.metadata
aligner0 = reg.EdgeAligner(reader0, verbose=True)
aligner0.run()
print '    merging...'
mshape = aligner0.mosaic_shape
mosaic = np.zeros(mshape, dtype=np.uint16)
for i, npos in enumerate(aligner0.positions):
    reg.paste(mosaic, reader0.read(c=0, series=i), npos)
skimage.io.imsave('scan_0_0.tif', mosaic)

aligners = []
for scan, filepath in enumerate(filepaths[1:], 1):
    print 'Scan %d:' % scan
    print '    reading %s' % filepath
    reader = reg.Reader(filepath)
    aligner = reg.LayerAligner(reader, aligner0, verbose=True)
    aligner.run()
    aligners.append(aligner)
    print '    merging...'
    mosaic = np.zeros(mshape, dtype=np.uint16)
    for i, npos in enumerate(aligner.positions):
        reg.paste(mosaic, reader.read(c=0, series=i), npos)
    skimage.io.imsave('scan_%d_0.tif' % scan, mosaic)


try:
    __IPYTHON__
except:
    reg._deinit_bioformats()
