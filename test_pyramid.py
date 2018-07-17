import os
import numpy as np
import skimage.io
import skimage.transform
from skimage.util.dtype import convert as skimage_convert

orig_dir = os.getcwd()
os.chdir(os.path.dirname(os.path.abspath(__file__)) + '/output/BP40')

out_filename = 'test_pyramid.tif'
TILE_SIZE = 1024

print "level 0"
for i in range(4):
    print "  channel", i
    kwargs = {}
    if i == 0:
        with open('ome.xml') as f:
            kwargs['description'] = f.read()
    else:
        kwargs['append'] = True
    img = skimage.io.imread('cycle_0_channel_{}.tif'.format(i))
    skimage.io.imsave(
        out_filename, img, bigtiff=True, metadata=None,
        tile=(TILE_SIZE, TILE_SIZE), photometric='minisblack', **kwargs
    )

dtype = img.dtype
level_shape = img.shape
level = 0

while any(s > TILE_SIZE for s in level_shape):
    prev_level = level
    level += 1
    print "level", level
    for i in range(4):
        print "  channel", i
        img = skimage.io.imread(
            out_filename, series=prev_level, key=i
        )
        img = skimage.transform.pyramid_reduce(img, multichannel=False)
        img = skimage_convert(img, dtype)
        skimage.io.imsave(
            out_filename, img, bigtiff=True, metadata=None, append=True,
            tile=(TILE_SIZE, TILE_SIZE), photometric='minisblack'
        )
    level_shape = img.shape

os.chdir(orig_dir)
