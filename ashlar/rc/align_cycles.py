import numpy as np
import scipy.spatial
import skimage.transform
import tqdm

from .. import thumbnail, utils
from . import rotation_utils


def correct_position(layer_aligner, angle):
    cycle_tform = thumbnail.align_cycles(
        layer_aligner.reference_aligner.reader,
        layer_aligner.reader,
        scale=layer_aligner.reference_aligner.reader.thumbnail_scale,
        angle=angle
    )
    layer_aligner.cycle_tform = cycle_tform
    layer_aligner.corrected_nominal_positions = (
        np.fliplr(cycle_tform(np.fliplr(layer_aligner.reader.metadata.positions)))
    )


def set_pairs(layer_aligner):
    reference_positions = layer_aligner.reference_aligner.metadata.positions
    dist = scipy.spatial.distance.cdist(reference_positions,
                                        layer_aligner.corrected_nominal_positions)
    layer_aligner.reference_idx = np.argmin(dist, 0)
    layer_aligner.reference_positions = reference_positions[layer_aligner.reference_idx]
    layer_aligner.reference_aligner_positions = layer_aligner.reference_aligner.positions[layer_aligner.reference_idx]


def register(layer_aligner, t):
    its, ref_img, img = layer_aligner.overlap(t)
    if np.any(np.array(its.shape) == 0):
        return np.nan
    return utils.register_angle(ref_img, img, layer_aligner.filter_sigma)


def tile_edge_score(reader, channel):
    score = np.zeros(reader.metadata.num_images)
    for i in tqdm.trange(reader.metadata.num_images):
        score[i] = rotation_utils.var_of_laplacian(reader.read(i, channel))
    return score


def refine_angle(layer_aligner, rank=None, top_k=None):
    tiles = np.arange(layer_aligner.reader.metadata.num_images)
    if rank is not None:
        tiles = tiles[rank]
    if top_k is None:
        top_k = len(tiles)
    top_k = min(top_k, len(tiles))
    tiles = tiles[:top_k]
    from joblib import Parallel, cpu_count, delayed
    img_pairs = filter(lambda x: x[0].size > 0, [
        layer_aligner.overlap(t)[1:3]
        for t in tiles
    ])
    angles = Parallel(verbose=1, n_jobs=cpu_count())(
        delayed(utils.register_angle)(img1, img2, layer_aligner.filter_sigma)
        for img1, img2 in img_pairs
    )
    return np.nanmedian(angles)


import skimage.transform

from .. import reg
from . import preproc_reader, rotation_utils


def process_rotated_reader(
    reader, edge_aligner,
    channel=None, max_shift=15,
    filter_sigma=0.0
):
    c2r = reader
    c1e = edge_aligner
    c21l = reg.LayerAligner(
        c2r, c1e, verbose=True,
        channel=channel, max_shift=max_shift,
        filter_sigma=filter_sigma
    )
    c21l.make_thumbnail()

    SCALE = c1e.reader.thumbnail_scale

    correct_position(c21l, angle=None)
    set_pairs(c21l)

    if c21l.cycle_tform.rotation == 0:
        c21l.register_all()
        c21l.calculate_positions()
        c21l.mosaic_shape = c1e.mosaic_shape
        return c21l

    edgy_scores = tile_edge_score(c2r, c21l.channel)
    angle = refine_angle(c21l, rank=np.argsort(edgy_scores)[::-1], top_k=30)
    print(f'\r    refined cycle rotation = {angle:.4f} degrees')

    ori_shape = c2r.metadata.size
    rotation_slice = rotation_utils.compute_slice(ori_shape, angle)
    crop_shape = np.zeros(ori_shape)[rotation_slice].shape

    cycle_tform = thumbnail.align_cycles(
        c21l.reference_aligner.reader,
        c21l.reader,
        scale=SCALE,
        angle=angle,
    )
    corrected_positions = (
        np.fliplr(cycle_tform(np.fliplr(c2r.metadata.centers)))
    )
    corrected_positions += np.multiply(-0.5, crop_shape)

    c2rr = preproc_reader.PreprocBioformatsReader(
        c2r.path,
        angle=angle,
        center_crop_shape=crop_shape
    )
    c2rr.metadata._positions = corrected_positions
    c2rr.metadata.extent = (
        c2rr.metadata.positions.max(axis=0) + c2rr.metadata.size - c2rr.metadata.origin
    )

    rthumbnail = skimage.transform.rotate(c21l.reader.thumbnail, angle=angle, resize=True)
    t_offset = .5*np.subtract(rthumbnail.shape, SCALE * c2rr.metadata.extent)
    ro, co = np.around(t_offset).astype(int)
    slice_r = slice(None) if ro == 0 else slice(ro, -ro)
    slice_c = slice(None) if co == 0 else slice(co, -co)
    c2rr.thumbnail = rthumbnail[slice_r, slice_c]


    c21lr = reg.LayerAligner(
        c2rr, c1e, verbose=True,
        channel=channel, max_shift=max_shift,
        filter_sigma=filter_sigma
    )
    c21lr.corrected_nominal_positions = corrected_positions
    set_pairs(c21lr)
    c21lr.register_all()
    c21lr.calculate_positions()
    c21lr.mosaic_shape = c1e.mosaic_shape

    return c21lr
