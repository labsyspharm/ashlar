from __future__ import division
import sys
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import skimage.io
import skimage.color


MAX_L = 60.0
MAX_C = 100.0


paths = sorted(sys.argv[1:])
img0 = skimage.io.imread(paths[0])
rgb = np.zeros(img0.shape + (3,))

lch = np.empty_like(rgb, np.float64)
hues = np.linspace(0, 2*np.pi, len(paths)+1)[:-1]
for path, hue in zip(paths, hues):
    print path
    img = skimage.io.imread(path).astype(np.float64)
    img /= img.max()
    lch[:, :, 0] = img * MAX_L
    lch[:, :, 1] = img * MAX_C
    lch[:, :, 2] = hue
    rgb += skimage.color.lab2rgb(skimage.color.lch2lab(lch))
np.clip(rgb, 0, 1, out=rgb)

cmap_lch = np.empty((len(hues), 1, 3))
cmap_lch[:, 0, 0] = MAX_L
cmap_lch[:, 0, 1] = MAX_C
cmap_lch[:, 0, 2] = hues
cmap_rgb = skimage.color.lab2rgb(skimage.color.lch2lab(cmap_lch))
cmap = mpl.colors.ListedColormap(cmap_rgb.reshape((len(hues), 3)))

skimage.io.imsave('colorized.jpg', rgb, quality=95)

plt.imshow(rgb, cmap=cmap)
plt.colorbar()
plt.show()
