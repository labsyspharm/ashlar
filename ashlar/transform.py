import numpy as np
from skimage.transform import warp


def _barrel_mapping(xy, center, k):
    x, y = xy.T.astype(float)
    y0, x0 = center
    x -= x0
    y -= y0
    f = 1 + k  * (x ** 2 + y ** 2)
    xy[..., 0] = x * f + x0
    xy[..., 1] = y * f + y0
    return xy


def barrel_correction(
    image,
    k,
    center=None,
    output_shape=None,
    order=1,
    mode=None,
    cval=0,
    clip=True,
    preserve_range=False,
):

    if mode is None:
        mode = "constant"

    if center is None:
        center = np.array(image.shape)[:2] / 2

    warp_args = {"center": center, "k": k}

    return warp(image, _barrel_mapping, map_args=warp_args,
                output_shape=output_shape, order=order, mode=mode, cval=cval,
                clip=clip, preserve_range=preserve_range)
