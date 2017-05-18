from __future__ import division
import bioformats
import javabridge
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import shift
from skimage.feature import register_translation
from skimage.filters import laplace
from skimage.io import imsave


def lowpass(im, sigma):
    dtype = im.dtype
    im = im.astype('f8')
    im -= scipy.ndimage.gaussian_filter(im, sigma)
    im = np.maximum(im, 0)
    im = im.astype(dtype)
    return im


def get_position_range(metadata, dimension):
    attr = 'Position' + dimension.upper()
    values = set(getattr(metadata.image(i).Pixels.Plane(0), attr)
                 for i in range(metadata.image_count))
    return max(values) - min(values)


def get_position(metadata, i):
    plane = metadata.image(i).Pixels.Plane(0)
    return np.array([plane.PositionY, plane.PositionX])


def paste(target, im, pos):
    pos_f, pos_i = np.modf(pos)
    yi, xi = pos_i.astype('i8')
    # Clip against the top and left of the mosaic.
    if yi < 0:
        im = im[-yi:]
        yi = 0
    if xi < 0:
        im = im[:, -xi:]
        xi = 0
    # This is a bit wrong on the edges in the subpixel shift direction. The
    # fractional pixels that would be shifted off the edges of the image are
    # actually discarded. However since the images being tiled in this
    # application have far more overlap than a single pixel, it's irrelevant.
    target_slice = target[yi:yi+im.shape[0], xi:xi+im.shape[1]]
    target_slice[:, :] = np.maximum(target_slice, shift(im, pos_f))

    target_slice[0, :] += 6000
    target_slice[-1, :] += 6000
    target_slice[1:-1, 0] += 6000
    target_slice[1:-1, -1] += 6000


MARGIN_FACTOR = 1.1


javabridge.start_vm(class_path=bioformats.JARS)
# Hack module to fix py3 assumptions which break XML parsing.
bioformats.omexml.str = unicode

filename = '1_40X_BACKGROUND_ONLY_Scan_20170425_191309_01x4x00176.rcpnl'
#filename = '2_40X_LY6C_CD8_CD68_Scan_20170427_134107_01x4x00176.rcpnl'
filepath = '../../dloads/' + filename
ir = bioformats.ImageReader(filepath)
metadata = bioformats.OMEXML(bioformats.get_omexml_metadata(filepath))

img0 = np.flipud(ir.read(c=0, series=0, rescale=False))
x_range = get_position_range(metadata, 'x')
y_range = get_position_range(metadata, 'y')
px_node = metadata.image(0).Pixels.node
pixel_size_x = float(px_node.get('PhysicalSizeX')) * 1.02
pixel_size_y = float(px_node.get('PhysicalSizeY')) * 1.02
pixel_sizes = np.array([pixel_size_y, pixel_size_x])
mosaic_width = (x_range / pixel_size_x + img0.shape[1]) * MARGIN_FACTOR
mosaic_height = (y_range / pixel_size_y + img0.shape[0]) * MARGIN_FACTOR
mosaic_shape = np.trunc([mosaic_height, mosaic_width]).astype(int)

mosaic = np.zeros(mosaic_shape, dtype=img0.dtype)
# FIXME: Assumes we always start in the corner.
paste(mosaic, img0, (0, 0))
plane0_meta = metadata.image(0).Pixels.Plane(0)
pos0 = np.array([plane0_meta.PositionY, plane0_meta.PositionX])

for i in range(1, metadata.image_count):
    print "Registering %d/%d" % (i, metadata.image_count)
    img = np.flipud(ir.read(c=0, series=i, rescale=False))
    h, w = img.shape
    sy, sx = ((get_position(metadata, i) - pos0) / pixel_sizes).astype(int)
    mtile = mosaic[sy:sy+h, sx:sx+w]
    offsets , _, _ = register_translation(laplace(mtile), laplace(img), 100)
    #print "  offset: %f,%f" % (offsets[1], offsets[0])
    pos = offsets + [sy, sx]
    paste(mosaic, img, pos)

gamma_corrected = (mosaic/mosaic.max())**(1/2.2)
imsave('m2.jpg', gamma_corrected, qualty=90)
