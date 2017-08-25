from __future__ import division, print_function
import warnings
import sys
import glob
import re
import subprocess
import numpy as np
import skimage.io
import skimage.transform


# Target width.
TW = 1920

FORMAT = "frame_%04d.jpg"

CMD = ("ffmpeg -r 5 -i " + FORMAT + " -y -vcodec libx264 -profile:v main "
       "-level 3 -pix_fmt yuv420p -crf 18 -an alignment.mp4")

def main(args):
    # Sort paths by scan number, numerically.
    paths = sorted(glob.glob('scan_*_0.tif'),
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
        # TODO: contrast stretching
        with warnings.catch_warnings():
            warnings.filterwarnings(
                'ignore', r'Possible precision loss', UserWarning,
                '^skimage\.util\.dtype'
            )
            skimage.io.imsave(out_path, img_new)
    subprocess.call(CMD.split(' '))
    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
