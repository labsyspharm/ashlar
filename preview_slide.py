import sys
import numpy as np
import matplotlib.pyplot as plt
import modest_image
import reg

filepath = sys.argv[1]
assert filepath.endswith('.rcpnl')

reader = reg.Reader(filepath)
metadata = reader.metadata

positions = metadata.positions - metadata.origin
mshape = ((metadata.positions + metadata.sizes - metadata.origin).max(axis=0) + 1).astype(int)
mosaic = np.zeros(mshape, dtype=np.uint16)
modest_image.imshow(plt.gca(), mosaic, cmap='gray', vmin=0, vmax=65535)

total = reader.metadata.num_images
for i in range(total):
    sys.stdout.write("\rLoading %d/%d" % (i, total))
    sys.stdout.flush()
    reg.paste(mosaic, np.flipud(reader.read(c=0, series=i)), positions[i])
    plt.pause(0.01)
print
plt.show()

try:
    __IPYTHON__
except:
    reg._deinit_bioformats()
