import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import modest_image
from .. import reg


def main(argv=sys.argv):
    filepath = argv[1]
    assert filepath.endswith('.rcpnl')

    reader = reg.Reader(filepath)
    metadata = reader.metadata

    positions = metadata.positions - metadata.origin
    mshape = ((metadata.positions + metadata.size - metadata.origin).max(axis=0) + 1).astype(int)
    mosaic = np.zeros(mshape, dtype=np.uint16)

    total = reader.metadata.num_images
    for i in range(total):
        sys.stdout.write("\rLoading %d/%d" % (i + 1, total))
        sys.stdout.flush()
        reg.paste(mosaic, reader.read(c=0, series=i), positions[i])
    print

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
