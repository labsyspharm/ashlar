import warnings
import sys
import glob
import re
import subprocess
import gc
import numpy as np
import skimage.transform
import skimage.exposure
from .. import utils


# Target width.
TW = 1920

FORMAT = "frame_%04d.jpg"
MOVIE_FILENAME = 'alignment.mp4'

CMD = ("ffmpeg -v error -r 5 -y -i " + FORMAT + " -an"
       " -vcodec libx264 -profile:v main -level 3 -pix_fmt yuv420p -crf 18"
       " " + MOVIE_FILENAME)


def main(argv=sys.argv):
    # Sort paths by cycle number, numerically.
    paths = sorted(glob.glob('cycle_*_channel_0.tif'),
                   key=lambda s: int(s.split('_')[1]))
    for i, in_path in enumerate(paths):
        out_path = FORMAT % i
        print('%s -> %s' % (in_path, out_path))
        img = skimage.io.imread(in_path)
        if i == 0:
            h, w = img.shape
            # Tall images will be rotated by 90 degrees.
            rotate = h > w
            if rotate:
                w, h = h, w
            scale = w / TW
            # Ensure final dimensions are multiples of 2, per video codec
            # restrictions.
            w = int(w / scale // 2 * 2)
            h = int(h / scale // 2 * 2)
        if rotate:
            img = np.rot90(img)
        img_new = skimage.transform.resize(img, (h, w), mode='reflect')
        # Free this memory as soon as possible.
        del img
        gc.collect()
        bins = np.linspace(0, img_new.max(), 3000)
        counts, _ = np.histogram(img_new, bins=bins)
        # Find peak, skipping first and last bin.
        vmin = bins[np.argmax(counts[1:-1]) + 1]
        vmax = np.percentile(img_new, 99.5)
        img_new = skimage.exposure.rescale_intensity(img_new, (vmin, vmax))
        img_new = skimage.exposure.adjust_gamma(img_new, 1/2.2)
        utils.imsave(out_path, img_new)
    print('rendering frames to %s' % (MOVIE_FILENAME))
    subprocess.call(CMD.split(' '))
    return 0


if __name__ == '__main__':
    main()
