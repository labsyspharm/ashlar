import numpy as np
from skimage.transform import warp


def _cart2polar_mapping(xy, center):
    x, y = xy.T
    y0, x0 = center

    rho = np.sqrt((x - x0) ** 2 + (y - y0) ** 2)
    theta = np.arctan2(y - y0, x - x0)

    # Normalize
    rho *= x.max() / rho.max()
    theta += np.pi
    theta *= y.max() / (2*np.pi)

    xy[..., 0] = rho
    xy[..., 1] = theta
    return xy


def cart2polar(image, center=None, output_shape=None, order=1, mode=None,
               cval=0, clip=True, preserve_range=False):

    """
    Perform a cartesian to polar coordinates transform.

    Parameters
    ----------
    image : ndarray
        Input image.
    center : (row, column) tuple or (2,) ndarray, optional
        Center coordinate of transformation.

    Returns
    -------
    polar : ndarray
        Polar version of the input.

    Other parameters
    ----------------
    output_shape : tuple (rows, cols), optional
        Shape of the output image generated. By default the shape of the input
        image is preserved.
    order : int, optional
        The order of the spline interpolation, default is 1. The order has to
        be in the range 0-5. See `skimage.transform.warp` for detail.
    mode : {'constant', 'edge', 'symmetric', 'reflect', 'wrap'}, optional
        Points outside the boundaries of the input are filled according
        to the given mode, with 'constant' used as the default. Modes match
        the behaviour of `numpy.pad`.
    cval : float, optional
        Used in conjunction with mode 'constant', the value outside
        the image boundaries.
    clip : bool, optional
        Whether to clip the output to the range of values of the input image.
        This is enabled by default, since higher order interpolation may
        produce values outside the given input range.
    preserve_range : bool, optional
        Whether to keep the original range of values. Otherwise, the input
        image is converted according to the conventions of `img_as_float`.
    """
    if mode is None:
        mode = "constant"

    if center is None:
        center = np.array(image.shape)[:2] / 2

    warp_args = {"center": center}

    return warp(image, _cart2polar_mapping, map_args=warp_args,
                output_shape=output_shape, order=order, mode=mode, cval=cval,
                clip=clip, preserve_range=preserve_range)


def _polar2cart_mapping(xy, center):
    x, y = xy.T
    x0, y0 = center

    # Rescale y value to match 0:2pi
    y = y / y.max() * 2 * np.pi
    xx = x * np.cos(y) + x0
    yy = x * np.sin(y) + y0

    xy[..., 0] = yy
    xy[..., 1] = xx
    return xy


def polar2cart(image, center=None, output_shape=None, order=1, mode=None,
               cval=0, clip=True, preserve_range=False):
    """
    Perform a polar to cartesian coordinates transform.

    Parameters
    ----------
    image : ndarray
        Input image.
    center : (row, column) tuple or (2,) ndarray, optional
        Center coordinate of transformation.

    Returns
    -------
    cartesian : ndarray
        Cartesian version of the input.
        Rows correspond to radius and columns to angle values.

    Other parameters
    ----------------
    output_shape : tuple (rows, cols), optional
        Shape of the output image generated. By default the shape of the input
        image is preserved.
    order : int, optional
        The order of the spline interpolation, default is 1. The order has to
        be in the range 0-5. See `skimage.transform.warp` for detail.
    mode : {'constant', 'edge', 'symmetric', 'reflect', 'wrap'}, optional
        Points outside the boundaries of the input are filled according
        to the given mode, with 'constant' used as the default. Modes match
        the behaviour of `numpy.pad`.
    cval : float, optional
        Used in conjunction with mode 'constant', the value outside
        the image boundaries.
    clip : bool, optional
        Whether to clip the output to the range of values of the input image.
        This is enabled by default, since higher order interpolation may
        produce values outside the given input range.
    preserve_range : bool, optional
        Whether to keep the original range of values. Otherwise, the input
        image is converted according to the conventions of `img_as_float`.
    """
    if mode is None:
        mode = "constant"

    if center is None:
        center = np.array(image.shape)[:2] / 2

    warp_args = {"center": center}

    return warp(image, _polar2cart_mapping, map_args=warp_args,
                output_shape=output_shape, order=order, mode=mode, cval=cval,
                clip=clip, preserve_range=preserve_range)


def _barrel_mapping(xy, center, k):
    x, y = xy.T.astype(float)
    y0, x0 = center
    x -= x0
    y -= y0
    # The warp mapping function defines the INVERSE map to apply, which in the
    # case of a distortion correction is the FORWARD map of the distortion
    # itself. See Fitzgibbon 2001 eq 1 for details.
    r2 = x ** 2 + y ** 2
    f = 1 + k * r2
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
    """Perform a transform to correct for barrel or pincushion distortion.

    Parameters
    ----------
    image : ndarray
        Input image.
    k : float
        Distortion parameter
    center : (row, column) tuple or (2,) ndarray, optional
        Center coordinate of transformation.

    Returns
    -------
    cartesian : ndarray
        Cartesian version of the input.
        Rows correspond to radius and columns to angle values.

    Other parameters
    ----------------
    output_shape : tuple (rows, cols), optional
        Shape of the output image generated. By default the shape of the input
        image is preserved.
    order : int, optional
        The order of the spline interpolation, default is 1. The order has to
        be in the range 0-5. See `skimage.transform.warp` for detail.
    mode : {'constant', 'edge', 'symmetric', 'reflect', 'wrap'}, optional
        Points outside the boundaries of the input are filled according
        to the given mode, with 'constant' used as the default. Modes match
        the behaviour of `numpy.pad`.
    cval : float, optional
        Used in conjunction with mode 'constant', the value outside
        the image boundaries.
    clip : bool, optional
        Whether to clip the output to the range of values of the input image.
        This is enabled by default, since higher order interpolation may
        produce values outside the given input range.
    preserve_range : bool, optional
        Whether to keep the original range of values. Otherwise, the input
        image is converted according to the conventions of `img_as_float`.

    Notes
    -----
    Radial distortion is modeled here using a one-parameter model described in
    [1]_ (Equation 1).

    References
    ----------
    .. [1] A. W. Fitzgibbon, "Simultaneous linear estimation of multiple view
       geometry and lens distortion," Proceedings of the 2001 IEEE Computer
       Society Conference on Computer Vision and Pattern Recognition. CVPR 2001,
       Kauai, HI, USA, 2001, pp. I-I. :DOI:`10.1109/CVPR.2001.990465`.

    """

    if mode is None:
        mode = "constant"

    if center is None:
        center = np.array(image.shape)[:2] / 2

    warp_args = {"center": center, "k": k}

    return warp(image, _barrel_mapping, map_args=warp_args,
                output_shape=output_shape, order=order, mode=mode, cval=cval,
                clip=clip, preserve_range=preserve_range)
