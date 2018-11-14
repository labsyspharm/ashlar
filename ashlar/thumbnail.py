from __future__ import print_function, division
import sys
try:
    import pathlib
except ImportError:
    import pathlib2 as pathlib
import numpy as np
from . import reg
from skimage.transform import rescale
from skimage.feature import register_translation
from skimage.io import imsave


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
        reg.paste(
            mosaic, rescale(reader.read(c=channel, series=i), scale), positions[i] * scale, np.maximum
        )
    print()

    return mosaic


def calculate_image_offset(img1, img2, upsample_factor=1):

    ref = reg.fft2(reg.whiten(img1))
    test = reg.fft2(reg.whiten(img2))

    shift, error, _ = register_translation(ref, test, upsample_factor, 'fourier')

    return shift


def calculate_cycle_offset(reader1, reader2, channel=0, scale=0.05, save=(False, False)):
    
    img1 = reader1.thumbnail_img \
        if hasattr (reader1, 'thumbnail_img') \
        else make_thumbnail(reader1, channel=channel, scale=scale)
    img2 = make_thumbnail(reader2, channel=channel, scale=scale)
    reader1.thumbnail_img = img1

    for file_path, img in zip(
        [reader1.path, reader2.path],
        [img1*save[0], img2*save[1]]
    ):
        if img.any(): _save_as_tif(img, file_path, post_fix='-thumbnail')

    img_offset = calculate_image_offset(img1, img2, int(1/scale)) * 1/scale
    ori_diff = reader2.metadata.origin - reader1.metadata.origin

    print(
        '\r        Estimated offset [y x] =',
        img_offset - ori_diff
    )

    return img_offset - ori_diff


def thumbnail(reader, post_fix='-thumbnail', channel=0, scale=0.05):

    img = make_thumbnail(reader, channel, scale)
    _save_as_tif(img, reader.path, post_fix)
    
    return img
    

def _save_as_tif(img, file_path, post_fix=''):
    
    input_path = pathlib.Path(file_path)

    if input_path.is_dir():
        raise RuntimeError("file_path musts point to a file not a directory")

    filename = input_path.name.replace('.', post_fix + '.', 1)
    out_path = input_path.parent / filename
    out_path = out_path.with_suffix('').with_suffix('.tif')

    imsave(out_path, img)