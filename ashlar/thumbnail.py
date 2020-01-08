from __future__ import print_function, division
import sys
import pathlib
import numpy as np
from . import utils
from skimage.transform import rescale
from skimage.feature import register_translation
from skimage.io import imsave
import warnings


def make_thumbnail(reader, channel=0, scale=0.05):
    
    metadata = reader.metadata
    positions = metadata.positions - metadata.origin
    coordinate_max = (positions + metadata.size).max(axis=0) 
    mshape = ((coordinate_max+ 1) * scale).astype(int)
    mosaic = np.zeros(mshape, dtype=np.uint16)

    total = reader.metadata.num_images
    for i in range(total):
        sys.stdout.write("\r        Loading %d/%d" % (i + 1, total))
        sys.stdout.flush()
        utils.paste(
            mosaic, rescale(reader.read(c=channel, series=i), scale, multichannel=False), 
            positions[i] * scale, np.maximum
        )
    print()

    return mosaic


def calculate_image_offset(img1, img2, upsample_factor=1):

    ref = utils.fft2(utils.whiten(img1, 0))
    test = utils.fft2(utils.whiten(img2, 0))

    shift, error, _ = register_translation(ref, test, upsample_factor, 'fourier')

    return shift


def calculate_cycle_offset(reader1, reader2, scale=0.05):
      
    if not (hasattr(reader1, 'thumbnail') and hasattr(reader2, 'thumbnail')):
        raise ValueError(
            'Thumbnail image needs to be assigned to the thumbnail attribute of both'
            ' readers before endering this function.'
        )

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
        '\r        Estimated offset [y x] =',
        img_offset
    )

    return img_offset
  

def _save_as_tif(img, file_path, post_fix=''):
    
    input_path = pathlib.Path(file_path)

    if input_path.is_dir():
        raise RuntimeError("file_path musts point to a file not a directory")

    filename = input_path.name.replace('.', post_fix + '.', 1)
    out_path = input_path.parent / filename
    out_path = out_path.with_suffix('.tif')

    utils.imsave(out_path, img)
