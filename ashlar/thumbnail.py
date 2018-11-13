from __future__ import print_function, division
import sys
import numpy as np
from . import reg
from skimage.transform import rescale
from skimage.feature import register_translation


def make_thumbnail(reader, channel=0, scale=0.05):
    
    metadata = reader.metadata
    positions = metadata.positions - metadata.origin
    coordinate_max = (positions + metadata.size).max(axis=0) 
    mshape = ((coordinate_max+ 1) * scale).astype(int)
    mosaic = np.zeros(mshape, dtype=np.uint16)

    total = reader.metadata.num_images
    for i in range(total):
        sys.stdout.write("\rLoading %d/%d" % (i + 1, total))
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

def calculate_cycle_offset(reader1, reader2, channel=0, scale=0.05):
    
    img1 = make_thumbnail(reader1, channel=channel, scale=scale)
    img2 = make_thumbnail(reader2, channel=channel, scale=scale)

    img_offset = calculate_image_offset(img1, img2, int(1/scale)) * 1/scale
    ori_diff = reader2.metadata.origin - reader1.metadata.origin

    print(
        'Estimated offset [ y, x] =',
        img_offset - ori_diff
    )

    return img_offset - ori_diff

