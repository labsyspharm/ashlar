from __future__ import print_function
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from .. import reg


def main(argv=sys.argv):

    try:
        import modest_image
    except ImportError:
        print("Please install ModestImage (Python 2.7 only) to use this script.")
        exit()

    filepath = argv[1]

    channel = 0
    if len(argv) >= 3:
        channel = int(argv[2])

    reader = reg.BioformatsReader(filepath)
    metadata = reader.metadata

    positions = metadata.positions - metadata.origin
    mshape = ((metadata.positions + metadata.size - metadata.origin).max(axis=0) + 1).astype(int)
    mosaic = np.zeros(mshape, dtype=np.uint16)

    total = reader.metadata.num_images
    for i in range(total):
        sys.stdout.write("\rLoading %d/%d" % (i + 1, total))
        sys.stdout.flush()
        reg.paste(mosaic, reader.read(c=channel, series=i), positions[i])
    print()

    ax = plt.gca()

    modest_image.imshow(ax, mosaic)

    h, w = metadata.size
    for xy in np.fliplr(positions):
        ax.add_patch(mpatches.Rectangle(xy, w, h, color='black', fill=False))

    plt.show()


    try:
        __IPYTHON__
    except:
        reg._deinit_bioformats()


if __name__ == '__main__':
    main()
