import numpy as np
from ashlar import utils
import skimage.data
import skimage.filters
import skimage.transform

import matplotlib.pyplot as plt
import tqdm

def blob_noise(fraction, noise_sd=0, blob_img_seed=None):
    blob_base = skimage.data.binary_blobs(
        100,
        blob_size_fraction=5/100,
        volume_fraction=fraction/(100*100),
        seed=blob_img_seed
    ).astype(float)
    rgn = np.random.default_rng()
    noise = rgn.normal(0, noise_sd, 100*100).reshape(100, 100)
    return blob_base + noise

# radial distortion
def radial_distort(xy, warp_center, k1=0.01, k2=0.002):
    assert warp_center in ['left', 'center', 'right']
    half = xy.mean(axis=0)
    if warp_center == 'right':
        center = [0, half[1]]
    elif warp_center == 'center':
        center = half
    elif warp_center == 'left':
        center = [2*half[0], half[1]]
    xy -= center
    xy /= half
    r = np.linalg.norm(xy, axis=1)
    m_r = 1 + k1*r + k2*r**2
    xy /= m_r.reshape(-1, 1)
    return xy * half + center


# test_range = np.linspace(0, 100*100, 1000)

# simulates two modes, one mode contains very few objects in overlapping blocks
# and accounts for 30% of total overlaps, rest of the overlaps has good quality
# for phase correlation
test_range = np.sort([
    *np.random.default_rng().normal(5, 10, 300),
    *np.random.default_rng().normal(100, 20, 700)
])
test_range = test_range[test_range > 0]


# testing and visualizatino function
# left panels shows results from original approach (w/o adding noise to
# laplacian filtered image) while gaussian noise is added to ALL the laplacian
# filtered image in the right panels

# the 1-percentile error cutoff tends to fail when the image has little noise
# (when `NOISE_SD` is low) or/and the distortion is significant (`K2` is high)
def plot_tests(
    test_range,
    SIGMA=1,
    NOISE_SD=0.005,
    K1=0.01,
    K2=0.002
):

    permutation_errors = np.empty(test_range.shape)
    edge_amplitudes = np.empty((*test_range.shape, 2))

    for idx, i in enumerate(tqdm.tqdm(test_range, desc='Computing edge amp', ascii=True)):
        # find edge amplitude threshold using "non-overlapping" blocks
        img1 = blob_noise(i, noise_sd=NOISE_SD)
        img2 = blob_noise(i, noise_sd=NOISE_SD)
        permutation_errors[idx] = utils.register(img1, img2, SIGMA, upsample=1)[1]
        edge_amplitudes[idx] = [
            utils.edge_amplitude(img1, SIGMA),
            utils.edge_amplitude(img2, SIGMA)
        ]
        
    noise_factor = skimage.filters.threshold_triangle(edge_amplitudes)

    permutation_errors_noise = np.empty(test_range.shape)
    for idx, i in enumerate(tqdm.tqdm(test_range, desc='Computing errors', ascii=True)):
        # calculate permutation errors w/ and w/o added noise
        img1 = blob_noise(i, noise_sd=NOISE_SD)
        img2 = blob_noise(i, noise_sd=NOISE_SD)
        permutation_errors_noise[idx] = utils.register(
            img1, img2, SIGMA,
            upsample=1, noise_factor=noise_factor
        )[1]


    errors = np.empty(test_range.shape)
    shifts = np.empty((*test_range.shape, 2))

    errors_noise = np.empty(test_range.shape)
    shifts_noise = np.empty((*test_range.shape, 2))

    for idx, i in enumerate(tqdm.tqdm(test_range, desc='Registering imgs', ascii=True)):
        # synthesize "overlapping blocks" and add gaussian noise and barrel
        # distortion
        img1 = blob_noise(i, noise_sd=NOISE_SD, blob_img_seed=1001)
        img2 = blob_noise(i, noise_sd=NOISE_SD, blob_img_seed=1001)
        img1 = skimage.transform.warp(
            img1, radial_distort, map_args=dict(warp_center='left', k1=K1, k2=K2)
        )
        img2 = skimage.transform.warp(
            img2, radial_distort, map_args=dict(warp_center='right', k1=K1, k2=K2)
        )

        shifts[idx], errors[idx] = utils.register(img1, img2, SIGMA, upsample=1)
        shifts_noise[idx], errors_noise[idx] = utils.register(
            img1, img2, SIGMA,
            upsample=1, noise_factor=noise_factor
        )

    passed = errors < np.percentile(permutation_errors, 1)
    passed_noise = errors_noise < np.percentile(permutation_errors_noise, 1)

    fig, axs = plt.subplots(2, 2, sharex=True, sharey=False)
    fig.suptitle(f'SIGMA={SIGMA}, NOISE_SD={NOISE_SD}, K1={K1}, K2={K2}')

    kwargs = dict(linewidths=0, s=8, alpha=0.5)
    axs[0][0].set_title('error w/o noise')
    axs[0][0].scatter(test_range, permutation_errors, c='#666666', **kwargs)
    axs[0][0].scatter(test_range, errors, c=passed, cmap='PiYG', **kwargs)
    axs[0][0].axhline(np.percentile(permutation_errors, 1), c='k', lw=1)
    axs[0][1].set_title('error w/ noise')
    axs[0][1].scatter(test_range, permutation_errors_noise, c='#666666', **kwargs)
    axs[0][1].scatter(test_range, errors_noise, c=passed_noise, cmap='PiYG', **kwargs)
    axs[0][1].axhline(np.percentile(permutation_errors_noise, 1), c='k', lw=1)

    axs[1][0].set_title('shift distance w/o noise')
    axs[1][0].scatter(test_range, np.linalg.norm(shifts, axis=1), c=passed, cmap='PiYG', **kwargs)
    axs[1][1].set_title('shift distance w/ noise')
    axs[1][1].scatter(test_range, np.linalg.norm(shifts_noise, axis=1), c=passed_noise, cmap='PiYG', **kwargs)
