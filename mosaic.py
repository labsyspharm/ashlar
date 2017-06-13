from __future__ import division
import sys
import gc
import time
from collections import defaultdict
import bioformats
import javabridge
import matplotlib.pyplot as plt
import numpy as np
try:
    import pyfftw
    np.fft = pyfftw.interfaces.numpy_fft
except ImportError:
    print "pyfftw not found, falling back to numpy.fft"
import scipy.ndimage
import skimage.feature
import skimage.io
from skimage.restoration.uft import laplacian


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
    img = crop_like(img, target_slice)
    img_subpixel_shifted = scipy.ndimage.shift(img, pos_f)
    result = np.maximum(target_slice, img_subpixel_shifted)
    result = np.minimum(result, 1.0)
    target_slice[:, :] = result

    if debug:
        # Render a faint outline of the pasted image.
        # TODO .1 is arbitrary and should be calculated from the data.
        target_slice[0, :] += .1
        target_slice[-1, :] += .1
        target_slice[1:-1, 0] += .1
        target_slice[1:-1, -1] += .1
        target_slice[:, :] = np.minimum(target_slice[:, :], 1)


def empty_aligned(shape, dtype, align=32):
    """Get an empty array with the specified byte alignment."""
    dtype = np.dtype(dtype)
    n = dtype.itemsize * np.prod(shape)
    a = np.empty(n + (align - 1), dtype=np.uint8)
    data_align = a.ctypes.data % align
    offset = 0 if data_align == 0 else (align - data_align)
    b = a[offset : offset + n]
    return b.view(dtype).reshape(shape)


def laplace(image, ksize=3):
    """Find the edges of an image using the Laplace operator.


    Copied from skimage.filters.edges, with explicit aligned output from
    convolve. Also the mask option was dropped.

    """

    image = skimage.img_as_float(image)
    # Create the discrete Laplacian operator - We keep only the real part of the
    # filter.
    _, laplace_op = laplacian(image.ndim, (ksize, ) * image.ndim)
    output = empty_aligned(image.shape, 'complex64')
    output.imag[:] = 0
    scipy.ndimage.convolve(image, laplace_op, output.real)
    return output


def subtract(a, b):
    return np.where(a >= b, a - b, 0)


def crop_like(img, target):
    if (img.shape[0] > target.shape[0]):
        img = img[:target.shape[0], :]
    if (img.shape[1] > target.shape[1]):
        img = img[:, :target.shape[1]]
    return img


def correct_illumination(img, ff):
    output = img / ff
    np.minimum(output, 1, output)
    return output


def read_ff(channel):
    ff = skimage.io.imread('flat_field_ch%d.tif' % channel).astype('u2').astype('f8')
    ff = scipy.ndimage.gaussian_filter(ff, 200)
    ff /= ff.mean()
    return ff


def gamma_correct(a, gamma):
    #max = np.iinfo(a.dtype).max
    #return (a / max)**(1 / gamma)
    return a ** (1 / gamma)


def read_image(reader, **kwargs):
    img = reader.read(**kwargs)
    img = np.flipud(img)
    return img


MARGIN_FACTOR = 1.1
PX_SCALE_X = 1.02 #1.010
PX_SCALE_Y = 1.02 #1.011
GAMMA = 2.2
#WX, WY = 4000, 4000
WX, WY = 2000, 2000

javabridge.start_vm(class_path=bioformats.JARS)
#signal.signal(signal.SIGINT, sigint_handler)
# Hack module to fix py3 assumptions which break XML parsing.
bioformats.omexml.str = unicode

filenames = (
    '1_40X_BACKGROUND_ONLY_Scan_20170425_191309_01x4x00176.rcpnl',
    '2_40X_LY6C_CD8_CD68_Scan_20170427_134107_01x4x00176.rcpnl',
    # '3_40X_BACKGROUND_ONLY_Scan_20170428_121003_01x4x00176.rcpnl',
    # '4_40X_B220_CD4_CD49B_Scan_20170501_120526_01x4x00176.rcpnl',
    # '5_40X_BACKGROUND_ONLY_can_20170502_213630_01x4x00176.rcpnl',
    # '6_40X_CD11B_FOXP3_VIMENTIN_Scan_20170505_113103_01x4x00176.rcpnl',
)

positions = defaultdict(dict)

for scan, filename in enumerate(filenames, 1):

    print "Scan %d\n==========" % scan
    filepath = sys.argv[1] + '/' + filename
    ir = bioformats.ImageReader(filepath)
    metadata = bioformats.OMEXML(bioformats.get_omexml_metadata(filepath))

    ff = read_ff(0)
    img0 = read_image(ir, c=0, series=0)
    # Warm up fftw.
    skimage.feature.register_translation(laplace(img0), laplace(img0), 10, 'fourier')
    img0 = correct_illumination(img0, ff)
    x_range = get_position_range(metadata, 'x')
    y_range = get_position_range(metadata, 'y')
    px_node = metadata.image(0).Pixels.node
    pixel_size_x = float(px_node.get('PhysicalSizeX')) * PX_SCALE_X
    pixel_size_y = float(px_node.get('PhysicalSizeY')) * PX_SCALE_Y
    pixel_sizes = np.array([pixel_size_y, pixel_size_x])
    mosaic_width = (x_range / pixel_size_x + img0.shape[1]) * MARGIN_FACTOR
    mosaic_height = (y_range / pixel_size_y + img0.shape[0]) * MARGIN_FACTOR
    # FIXME Temporary clipping to subregion [0,0]-[WX-1,WY-1].
    mosaic_width = min(mosaic_width, WX)
    mosaic_height = min(mosaic_height, WY)
    mosaic_shape = np.trunc([mosaic_height, mosaic_width]).astype(int)

    t_mos_start = time.time()
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
    reg_time = 0
    n_regs = 0
    for i in range(first_image, metadata.image_count):
        print "\rRegistering tile %d/%d" % (i + 1, metadata.image_count),
        sys.stdout.flush()
        sy, sx = ((get_position(metadata, i) - pos0) / pixel_sizes).astype(int)
        # FIXME Temporary optimization for working with subregions when testing.
        if sx > mosaic_width or sy > mosaic_height:
            continue
        img = read_image(ir, c=0, series=i)
        img = correct_illumination(img, ff)
        h, w = img.shape
        reftile_f = pyfftw.builders.fft2(laplace(reference[sy:sy+h, sx:sx+w]),
                                         avoid_copy=True)()
        img = crop_like(img, reftile_f)
        img_f = pyfftw.builders.fft2(laplace(img), avoid_copy=True)()
        t_start = time.time()
        shift, _, _ = skimage.feature.register_translation(reftile_f, img_f,
                                                           10, 'fourier')
        reg_time += time.time() - t_start
        n_regs += 1
        #print "  shift: %f,%f" % (shift[1], shift[0])
        pos = shift + [sy, sx]
        positions[scan][i] = pos
        paste(mosaic, img, pos)
        gc.collect()
    print
    print ("Registration: %g ms/frame (%g seconds, %d frames)"
           % (reg_time / n_regs * 1000, reg_time, n_regs))
    print "Total mosaicing time: %g s" % (time.time() - t_mos_start)

    gamma_corrected = gamma_correct(mosaic, GAMMA)
    skimage.io.imsave('scan_%d.jpg' % scan, gamma_corrected, quality=95)
    del gamma_corrected
    gc.collect()

    if 'background' in filename.lower():

        ir_bg = ir

    else:

        for c in range(1, metadata.image(0).Pixels.channel_count):

            print "Aligning channel %d" % c

            mosaic = np.zeros_like(reference)
            mosaic_bg = np.zeros_like(reference)
            for mode, s, r, m in (('bg', scan-1, ir_bg, mosaic_bg),
                                  ('fg', scan, ir, mosaic)):
                ff = read_ff(c)
                for i in range(0, metadata.image_count):
                    print "\r  %s tile %d/%d" % (mode, i + 1,
                                                 metadata.image_count),
                    sys.stdout.flush()
                    try:
                        pos = positions[s][i]
                    except KeyError:
                        # Temporary while working with subregions.
                        continue
                    img = read_image(r, c=c, series=i)
                    img = correct_illumination(img, ff)
                    paste(m, img, pos)
                    gc.collect()
                print
                skimage.io.imsave('scan_%d_%d.jpg' % (s, c), m, quality=95)

            mosaic = subtract(mosaic, mosaic_bg)
            del mosaic_bg
            gc.collect()
            mosaic = gamma_correct(mosaic, GAMMA)
            skimage.io.imsave('cycle_%d_%d.jpg' % (scan//2, c), mosaic,
                              quality=95)
            del mosaic
            gc.collect()

javabridge.kill_vm()
