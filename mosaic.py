from __future__ import division
import sys
import gc
from collections import defaultdict
import bioformats
import javabridge
import matplotlib.pyplot as plt
import numpy as np
import scipy.ndimage
import skimage.feature
import skimage.filters
import skimage.io


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


def paste(target, img, pos, debug=False):
    """Composite img into target using maximum intensity projection."""
    pos_f, pos_i = np.modf(pos)
    yi, xi = pos_i.astype('i8')
    # Clip img to the edges of the mosaic.
    # TODO Also need to handle the bottom and right edges. We have padding
    # there for now so it's not important yet.
    if yi < 0:
        img = img[-yi:]
        yi = 0
    if xi < 0:
        img = img[:, -xi:]
        xi = 0
    # This is a bit wrong on the edges in the subpixel shift direction. The
    # fractional pixels that would be shifted off the edges of the image are
    # actually discarded. However since the images being tiled in this
    # application have far more overlap than a single pixel, it's irrelevant.
    target_slice = target[yi:yi+img.shape[0], xi:xi+img.shape[1]]
    img_subpixel_shifted = scipy.ndimage.shift(img, pos_f)
    target_slice[:, :] = np.maximum(target_slice, img_subpixel_shifted)

    if debug:
        # Render a faint outline of the pasted image.
        # TODO 6000 is arbitrary and should be calculated from the data.
        target_slice[0, :] += 6000
        target_slice[-1, :] += 6000
        target_slice[1:-1, 0] += 6000
        target_slice[1:-1, -1] += 6000


def difference(a, b):
    return np.where(a >= b, a - b, 0)


def gamma_correct(a, gamma):
    max = np.iinfo(a.dtype).max
    return (a / max)**(1 / gamma)




MARGIN_FACTOR = 1.1
PX_SCALE_X = 1.02 #1.010
PX_SCALE_Y = 1.02 #1.011
GAMMA = 2.2
#WX, WY = 4000, 4000
WX, WY = 999999, 999999

javabridge.start_vm(class_path=bioformats.JARS)
# Hack module to fix py3 assumptions which break XML parsing.
bioformats.omexml.str = unicode

filenames = (
    '1_40X_BACKGROUND_ONLY_Scan_20170425_191309_01x4x00176.rcpnl',
    '2_40X_LY6C_CD8_CD68_Scan_20170427_134107_01x4x00176.rcpnl',
    '3_40X_BACKGROUND_ONLY_Scan_20170428_121003_01x4x00176.rcpnl',
    '4_40X_B220_CD4_CD49B_Scan_20170501_120526_01x4x00176.rcpnl',
    '5_40X_BACKGROUND_ONLY_can_20170502_213630_01x4x00176.rcpnl',
    '6_40X_CD11B_FOXP3_VIMENTIN_Scan_20170505_113103_01x4x00176.rcpnl',
)

positions = defaultdict(dict)

for scan, filename in enumerate(filenames, 1):

    print "Scan %d\n==========" % scan
    filepath = '../../dloads/' + filename
    ir = bioformats.ImageReader(filepath)
    metadata = bioformats.OMEXML(bioformats.get_omexml_metadata(filepath))

    # bg = np.flipud(ir.read(c=0, series=175, rescale=False))
    # bg = scipy.ndimage.gaussian_filter(bg, 100)
    #bg = 0

    img0 = np.flipud(ir.read(c=0, series=0, rescale=False))
    #img0 = difference(np.flipud(ir.read(c=0, series=0, rescale=False)), bg)
    x_range = get_position_range(metadata, 'x')
    y_range = get_position_range(metadata, 'y')
    px_node = metadata.image(0).Pixels.node
    pixel_size_x = float(px_node.get('PhysicalSizeX')) * PX_SCALE_X
    pixel_size_y = float(px_node.get('PhysicalSizeY')) * PX_SCALE_Y
    pixel_sizes = np.array([pixel_size_y, pixel_size_x])
    mosaic_width = (x_range / pixel_size_x + img0.shape[1]) * MARGIN_FACTOR
    mosaic_height = (y_range / pixel_size_y + img0.shape[0]) * MARGIN_FACTOR
    mosaic_shape = np.trunc([mosaic_height, mosaic_width]).astype(int)

    mosaic = np.zeros(mosaic_shape, dtype=img0.dtype)
    if scan == 1:
        reference = mosaic
        # FIXME: Assumes we always start in the corner.
        pos = (0, 0)
        positions[scan][0] = pos
        paste(mosaic, img0, pos)
        plane0_meta = metadata.image(0).Pixels.Plane(0)
        pos0 = np.array([plane0_meta.PositionY, plane0_meta.PositionX])
        first_image = 1
    else:
        first_image = 0

    print
    for i in range(first_image, metadata.image_count):
        print "\rRegistering tile %d/%d" % (i + 1, metadata.image_count),
        sys.stdout.flush()
        sy, sx = ((get_position(metadata, i) - pos0) / pixel_sizes).astype(int)
        if sx > WX or sy > WY:
            continue
        img = np.flipud(ir.read(c=0, series=i, rescale=False))
        #img = difference(np.flipud(ir.read(c=0, series=i, rescale=False)), bg)
        h, w = img.shape
        reftile_f = skimage.filters.laplace(reference[sy:sy+h, sx:sx+w])
        img_f = skimage.filters.laplace(img)
        shift, _, _ = skimage.feature.register_translation(reftile_f, img_f, 10)
        #print "  shift: %f,%f" % (shift[1], shift[0])
        pos = shift + [sy, sx]
        positions[scan][i] = pos
        paste(mosaic, img, pos)
        gc.collect()
    print

    gamma_corrected = gamma_correct(mosaic[:WY,:WX], GAMMA)
    skimage.io.imsave('scan_%d.jpg' % scan, gamma_corrected, quality=95)

    for c in range(1, metadata.image(0).Pixels.channel_count):
        print "Aligning channel %d" % c
        mosaic = np.zeros_like(reference)
        for i in range(0, metadata.image_count):
            print "\r  tile %d/%d" % (i + 1, metadata.image_count),
            sys.stdout.flush()
            try:
                pos = positions[scan][i]
            except:
                # Temporary while working with WX,WY
                continue
            img = np.flipud(ir.read(c=c, series=i, rescale=False))
            paste(mosaic, img, pos)
            gc.collect()
        print

        gamma_corrected = gamma_correct(mosaic[:WY,:WX], GAMMA)
        skimage.io.imsave('scan_%d_%d.jpg' % (scan, c), gamma_corrected,
                          quality=95)
