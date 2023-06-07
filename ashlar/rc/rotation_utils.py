import skimage.morphology
import skimage.measure
import skimage.transform
import numpy as np
import cv2


def compute_slice(shape, angle):
    assert(len(shape)) == 2
    mask = skimage.transform.rotate(np.ones(shape), angle) == 1
    return slice_rotated_internal_rect(mask)


def slice_rotated_internal_rect(mask):
    labeled = skimage.morphology.label(~mask)
    if labeled.max() == 0:
        return slice(None), slice(None)
    bboxes = skimage.measure.regionprops_table(labeled, properties=('bbox',))
    assert len(bboxes) == 4
    rows = np.concatenate([bboxes['bbox-0'], bboxes['bbox-2']])
    cols = np.concatenate([bboxes['bbox-1'], bboxes['bbox-3']])
    urows = np.unique(rows)
    ucols = np.unique(cols)
    if len(urows) < 4:
        urows = np.repeat(urows, 2)
    if len(ucols) < 4:
        ucols = np.repeat(ucols, 2)
    return slice(*urows[[1, -2]]), slice(*ucols[[1, -2]])


def var_of_laplacian(img, sigma=0):
    assert img.ndim == 2
    if sigma != 0:
        img = cv2.GaussianBlur(img.astype(np.float32), ksize=(0, 0), sigmaX=sigma)
    return cv2.Laplacian(img, cv2.CV_32F, ksize=1).var()