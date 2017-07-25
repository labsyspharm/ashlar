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


def get_position_range(metadata, dimension):
    attr = 'Position' + dimension.upper()
    values = set(getattr(metadata.image(i).Pixels.Plane(0), attr)
                 for i in range(metadata.image_count))
    return max(values) - min(values)


def get_position(metadata, i):
    plane = metadata.image(i).Pixels.Plane(0)
    return np.array([plane.PositionY, plane.PositionX])


def paste(target, img, pos, debug=False):
    """Composite img into target using maximum intensity projection.

    target: uint
    img: float

    """
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
    img = scipy.ndimage.shift(img, pos_f)
    np.clip(img, 0, 1, img)
    img = skimage.util.dtype.convert(img, target.dtype)
    target_slice[:, :] = np.maximum(target_slice, img)
    if debug:
        # Render a faint outline of the pasted image.
        # TODO 6000 is arbitrary and should be calculated from the data.
        # Also these lines can cause registration problems, so ideally
        # this step should be performed on the final images by using the
        # accumulated list of per-tile offsets.
        target_slice[0, :] += 6000
        target_slice[-1, :] += 6000
        target_slice[1:-1, 0] += 6000
        target_slice[1:-1, -1] += 6000
        np.clip(target_slice[:, :], 0, np.iinfo(target.dtype).max)


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

    image: any dtype
    returns: float

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
    """
    img: float
    ff: float (mean-normalized)
    returns: float
    """
    if ff is None:
        return img
    output = img / ff
    np.clip(output, 0, 1, output)
    return output


def read_ff(channel):
    """
    returns: float
    """
    try:
        ff = skimage.io.imread('flat_field_ch%d.tif' % channel)
    except IOError:
        print ("WARNING: No flat-field image for channel %d; not correcting"
               % channel)
        return None
    ff = ff / ff.mean()
    return ff


def gamma_correct(a, gamma):
    """
    a: any dtype
    returns: uint16
    """
    a = skimage.img_as_float(a)
    a **= 1 / gamma
    return skimage.img_as_uint(a)


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
bioformats.init_logger()
#signal.signal(signal.SIGINT, sigint_handler)
# Hack module to fix py3 assumptions which break XML parsing.
bioformats.omexml.str = unicode

filepaths = sorted(sys.argv[1:])
assert all(p.endswith('.rcpnl') for p in filepaths)

positions = defaultdict(dict)
ir_bg = None

for scan, filepath in enumerate(filepaths, 1):

    print "Scan %d\n==========" % scan
    ir = bioformats.ImageReader(filepath)
    metadata = bioformats.OMEXML(bioformats.get_omexml_metadata(filepath))

    ff = read_ff(0)
    img0 = read_image(ir, c=0, series=0)
    # Warm up fftw to make future timing accurate.
    #skimage.feature.register_translation(laplace(img0), laplace(img0), 10, 'fourier')
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

    #t_mos_start = time.time()
    mosaic = np.zeros(mosaic_shape, dtype=np.uint16)
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

    #print
    #reg_time = 0
    #n_regs = 0
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
        #t_start = time.time()
        shift, _, _ = skimage.feature.register_translation(reftile_f, img_f,
                                                           10, 'fourier')
        #reg_time += time.time() - t_start
        #n_regs += 1
        #print "  shift: %f,%f" % (shift[1], shift[0])
        pos = shift + [sy, sx]
        positions[scan][i] = pos
        paste(mosaic, img, pos)
        gc.collect()
    print
    #print ("Registration: %g ms/frame (%g seconds, %d frames)"
    #       % (reg_time / n_regs * 1000, reg_time, n_regs))
    #print "Total mosaicing time: %g s" % (time.time() - t_mos_start)

    skimage.io.imsave('scan_%d_0.tif' % scan, mosaic)

    if 'background' in filepath.lower():

        ir_bg = ir

    else:

        for c in range(1, metadata.image(0).Pixels.channel_count):

            print "Aligning channel %d" % c

            mosaic = np.zeros_like(reference)
            mosaic_bg = np.zeros_like(reference) if ir_bg else None
            for mode, s, r, m in (('bg', scan-1, ir_bg, mosaic_bg),
                                  ('fg', scan, ir, mosaic)):
                if r is None:
                    continue
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
                skimage.io.imsave('scan_%d_%d.tif' % (s, c), m)

            if ir_bg:
                mosaic = subtract(mosaic, mosaic_bg)
                del mosaic_bg
                gc.collect()
                skimage.io.imsave('cycle_%d_%d.tif' % (scan//2, c), mosaic)
                del mosaic
            gc.collect()

        if ir_bg:
            ir_bg.close()
        ir.close()

try:
    __IPYTHON__
except:
    javabridge.kill_vm()
