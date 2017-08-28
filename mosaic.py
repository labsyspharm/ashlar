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

aligners = []
mosaics = []

print 'Scan 0:'
print '    reading %s' % filepaths[0]
reader = reg.Reader(filepaths[0])
aligner = reg.EdgeAligner(reader, verbose=True)
aligner.run()
mshape = aligner.mosaic_shape
mosaic = reg.Mosaic(aligner, mshape, 'scan_0_%(channel)d.tif', verbose=True)
mosaic.run()
aligners.append(aligner)
mosaics.append(mosaic)

for scan, filepath in enumerate(filepaths[1:], 1):
    print 'Scan %d:' % scan
    print '    reading %s' % filepath
    reader = reg.Reader(filepath)
    aligner = reg.LayerAligner(reader, aligners[0], verbose=True)
    aligner.run()
    filename_format = 'scan_%d_%%(channel)d.tif' % scan
    mosaic = reg.Mosaic(aligner, mshape, filename_format, verbose=True)
    mosaic.run()
    aligners.append(aligner)
    mosaics.append(mosaic)


try:
    __IPYTHON__
except:
    reg._deinit_bioformats()
