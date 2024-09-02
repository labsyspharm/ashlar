import sys
import pathlib
import numpy as np
from . import utils
from skimage.transform import rescale
from skimage.registration import phase_cross_correlation
import tifffile


def calculate_scale(reader, default_scale=0.05, min_size=1000):
    """Return scaling factor for a thumbnail with a minimum size constraint."""
    positions = reader.metadata.positions - reader.metadata.origin
    full_shape = (positions + reader.metadata.size).max(axis=0)
    thumbnail_shape = full_shape * default_scale
    if any(thumbnail_shape >= min_size):
        scale = default_scale
    else:
        # Compute the scale needed to make the thumbnail min_size on the longest
        # side, but don't scale smaller images up.
        scale = min(min_size / max(full_shape), 1)
    return scale


def make_thumbnail(reader, channel=0, scale=0.05):
    positions = reader.metadata.positions - reader.metadata.origin
    full_shape = (positions + reader.metadata.size).max(axis=0)
    mshape = np.ceil(full_shape * scale).astype(int)
    mosaic = np.zeros(mshape, dtype=np.uint16)
    total = reader.metadata.num_images
    for i, pos_s in enumerate(positions * scale):
        sys.stdout.write("\r    assembling thumbnail %d/%d" % (i + 1, total))
        sys.stdout.flush()
        img = reader.read(c=channel, series=i)
        # We don't need anti-aliasing as long as the coarse features in the
        # images are bigger than the scale factor. This speeds up the rescaling
        # dramatically.
        img_s = rescale(img, scale, anti_aliasing=False)
        utils.paste(mosaic, img_s, pos_s, utils.pastefunc_blend)
    print()
    return mosaic


def calculate_image_offset(img1, img2, upsample_factor=1):
    ref = utils.window(utils.whiten(img1, 1))
    test = utils.window(utils.whiten(img2, 1))
    shift = phase_cross_correlation(
        ref,
        test,
        upsample_factor=upsample_factor,
        normalization=None
    )[0]
    return shift


def calculate_cycle_offset(reader1, reader2, scale=0.05):
    if not hasattr(reader1, 'thumbnail'):
        raise ValueError('reader1 does not have a thumbnail')
    if not hasattr(reader2, 'thumbnail'):
        raise ValueError('reader2 does not have a thumbnail')
    img1 = reader1.thumbnail
    img2 = reader2.thumbnail
    if img1.shape != img2.shape:
        padded_shape = np.array((img1.shape, img2.shape)).max(axis=0)
        padded_img1, padded_img2 = np.zeros(padded_shape), np.zeros(padded_shape)
        utils.paste(padded_img1, img1, [0, 0])
        utils.paste(padded_img2, img2, [0, 0])
        img1 = padded_img1
        img2 = padded_img2
    img_offset = calculate_image_offset(img1, img2, int(1/scale)) / scale
    img_offset -= (reader2.metadata.origin - reader1.metadata.origin)
    print(
        '\r    estimated cycle offset [y x] =',
        img_offset
    )
    return img_offset


def _save_as_tif(img, file_path, post_fix=''):
    input_path = pathlib.Path(file_path)
    if input_path.is_dir():
        raise RuntimeError("file_path must point to a file not a directory")
    filename = input_path.name.replace('.', post_fix + '.', 1)
    out_path = input_path.parent / filename
    out_path = out_path.with_suffix('.tif')
    tifffile.imwrite(out_path, img)
