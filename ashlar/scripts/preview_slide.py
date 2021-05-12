import warnings
import sys
import argparse
import numpy as np
import skimage.transform
import skimage.util
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.text as mtext
from .. import reg, utils


def main(argv=sys.argv):

    parser = argparse.ArgumentParser(
        description=(
            "Stitch a multi-tile image using the nominal stage positions and"
            " display the resulting image."
        )
    )
    parser.add_argument(
        "input", help="Path to image (BioFormats supported formats only)"
    )
    parser.add_argument(
        "-n", "--numbers", action="store_true", help="Display tile numbers"
    )
    parser.add_argument(
        "-b", "--bounds", action="store_true", help="Display tile bounds"
    )
    parser.add_argument(
        "-c", "--channel", type=int, default=0,
        help="Channel number to display; default: 0",
    )
    parser.add_argument(
        "-d", "--downsample", metavar="FACTOR", type=float, default=10,
        help="Downsample the image resolution by this factor for display (reduces"
        " memory requirements and improves performance); default: 10",
    )
    parser.add_argument(
        "-l", "--log", action="store_true", help="Log-transform pixel intensities"
        " (helps visualize dim images)",
    )
    args = parser.parse_args()

    reader = reg.BioformatsReader(args.input)
    metadata = reader.metadata

    resolution_scale = 1 / args.downsample
    positions = (metadata.positions - metadata.origin) * resolution_scale
    pmax = (positions + metadata.size * resolution_scale).max(axis=0)
    mshape = (pmax + 0.5).astype(int)
    mosaic = np.zeros(mshape, dtype=np.uint16)

    total = reader.metadata.num_images
    for i in range(total):
        sys.stdout.write("\rLoading %d/%d" % (i + 1, total))
        sys.stdout.flush()
        img = reader.read(c=args.channel, series=i)
        img = skimage.transform.rescale(img, resolution_scale, anti_aliasing=False)
        img = skimage.img_as_uint(img)
        if args.log:
            intensity_scale = 65535 / np.log(65535)
            img = (np.log(np.maximum(img, 1)) * intensity_scale).astype(np.uint16)
        # Round position so paste will skip the expensive subpixel shift.
        pos = np.round(positions[i])
        utils.paste(mosaic, img, pos, np.maximum)
    print()

    ax = plt.gca()

    plt.imshow(X=mosaic, axes=ax, extent=(0, pmax[1], pmax[0], 0))

    h, w = metadata.size * resolution_scale
    for i, (x, y) in enumerate(np.fliplr(positions)):
        if args.bounds:
            rect = mpatches.Rectangle((x, y), w, h, color='black', fill=False)
            ax.add_patch(rect)
        if args.numbers:
            xc = x + w / 2
            yc = y + h / 2
            circle = mpatches.Circle((xc, yc), w / 5, color='salmon', alpha=0.5)
            text = mtext.Text(
                xc, yc, str(i), color='k', size=10, ha='center', va='center'
            )
            ax.add_patch(circle)
            ax.add_artist(text)

    plt.show()


if __name__ == '__main__':
    main()
