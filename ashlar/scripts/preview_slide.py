from __future__ import print_function, division
import warnings
import sys
import argparse
import numpy as np
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
    args = parser.parse_args()

    try:
        from modest_image import imshow
    except ImportError:
        warnings.warn("Please install ModestImage to speed up image rendering")
        imshow = plt.imshow

    reader = reg.BioformatsReader(args.input)
    metadata = reader.metadata

    positions = metadata.positions - metadata.origin
    mshape = ((positions + metadata.size).max(axis=0) + 1).astype(int)
    mosaic = np.zeros(mshape, dtype=np.uint16)

    total = reader.metadata.num_images
    for i in range(total):
        sys.stdout.write("\rLoading %d/%d" % (i + 1, total))
        sys.stdout.flush()
        img = reader.read(c=args.channel, series=i)
        utils.paste(mosaic, img, positions[i], np.maximum)
    print()

    ax = plt.gca()

    imshow(X=mosaic, axes=ax)

    h, w = metadata.size
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
