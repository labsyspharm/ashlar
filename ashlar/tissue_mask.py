import numpy as np
import skimage.filters
import skimage.transform
import tqdm
import zarr

from . import thumbnail
from compress_bg.main import entropy_img_to_masks, local_entropy, plot_tissue_mask


def make_mask(
    thumbnail: np.ndarray,
    entropy_kernel_size: int = 5,
    dilation_radius: int = 5,
    plot: bool = False,
    level_center: float = 0.5,
    level_adjust: int = 0,
    figure_title: str = "",
):
    _tt = np.array(thumbnail)
    np.clip(_tt, np.percentile(_tt[_tt > 0], 0).astype(_tt.dtype), None, out=_tt)
    entropy_thumbnail = local_entropy(np.log1p(_tt), kernel_size=entropy_kernel_size)
    thumbnail = np.log1p(thumbnail)

    erange = np.ptp(entropy_thumbnail)
    _threshold = skimage.filters.threshold_otsu(entropy_thumbnail)
    _max = entropy_thumbnail.max() - 0.1 * erange
    _min = entropy_thumbnail.min() + 0.1 * erange
    _threshold = np.clip(_threshold, _min, _max)

    level_center = np.clip(
        level_center, (_min - _threshold) / erange, (_max - _threshold) / erange
    )

    # forcing threshold to be within 10th and 90th percent of the range
    _threshold += level_center * erange

    thresholds = np.concatenate(
        [
            np.linspace(entropy_thumbnail.min(), _threshold, 4)[1:],
            np.linspace(_threshold, entropy_thumbnail.max(), 4)[1:-1],
        ]
    )
    level_adjusts = np.arange(-2, 3, 1)
    masks = entropy_img_to_masks(
        thumbnail, entropy_thumbnail, thresholds, dilation_radius
    )
    mask = masks[list(level_adjusts).index(level_adjust)]

    if plot:
        fig = plot_tissue_mask(
            thumbnail,
            entropy_thumbnail,
            masks,
            thresholds,
            list(level_adjusts).index(level_adjust),
        )
        fig.suptitle(figure_title)

    return mask


def make_reader_mask(reader, thumbnail_px_size=20, **kwargs):
    Affine = skimage.transform.AffineTransform

    metadata = reader.metadata
    px_size = metadata.pixel_size
    factor = px_size / thumbnail_px_size
    timg = thumbnail.make_thumbnail(reader, scale=factor)
    mask = make_mask(timg, **kwargs)

    positions = factor * (metadata.positions - metadata.positions.min(axis=0))

    tile_mask = zarr.group()
    for ii, pp in enumerate(tqdm.tqdm(positions, desc="Making tissue mask")):
        h, w = metadata.tile_size(ii)
        out_shape = np.ceil(factor * metadata.tile_size(ii)).astype("int")
        tform = Affine(translation=pp[::-1])
        tmask = skimage.transform.warp(mask, tform, output_shape=out_shape, order=0)
        tile_mask[str(ii)] = skimage.transform.rescale(tmask, 1 / factor, order=0)[
            :h, :w
        ]
    return tile_mask
