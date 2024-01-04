import skimage.data
import skimage.util
import numpy as np
import tifffile
import skimage.filters

img = skimage.data.astronaut()[..., 0]
imgs = skimage.util.view_as_windows(img, window_shape=(100, 100), step=90)

h, w, *_ = imgs.shape
pixel_size = 0.3

positions = np.mgrid[:h, :w].reshape(2, -1).T * 90
random_gaussian_sigma = np.random.random(h * w) * 3

with tifffile.TiffWriter("01-ref.ome.tif", bigtiff=True) as tif:
    for ii, p, ss in zip(
        imgs.reshape(-1, *(100, 100)), positions, random_gaussian_sigma
    ):
        ii = skimage.filters.gaussian(ii, ss)
        ii = skimage.util.img_as_ubyte(ii)
        ii = np.array([ii] * 5, dtype=np.uint16)
        metadata = {
            "Pixels": {
                "PhysicalSizeX": pixel_size,
                "PhysicalSizeXUnit": "µm",
                "PhysicalSizeY": pixel_size,
                "PhysicalSizeYUnit": "µm",
            },
            "Plane": {
                "PositionX": [p[1] * pixel_size] * ii.shape[0],
                "PositionY": [p[0] * pixel_size] * ii.shape[0],
            },
        }
        tif.write(ii, metadata=metadata)

with tifffile.TiffWriter("02-bg.ome.tif", bigtiff=True) as tif:
    for ii, p, ss in zip(
        imgs.reshape(-1, *(100, 100)), positions, random_gaussian_sigma
    ):
        ii = skimage.filters.gaussian(ii, ss)
        ii = skimage.util.img_as_ubyte(ii)
        ii = ii.astype(np.uint16)
        ii = np.array([ii + 105] * 5, dtype=np.uint16)
        metadata = {
            "Pixels": {
                "PhysicalSizeX": pixel_size,
                "PhysicalSizeXUnit": "µm",
                "PhysicalSizeY": pixel_size,
                "PhysicalSizeYUnit": "µm",
            },
            "Plane": {
                "PositionX": [p[1] * pixel_size] * ii.shape[0],
                "PositionY": [p[0] * pixel_size] * ii.shape[0],
            },
        }
        tif.write(ii, metadata=metadata)

with tifffile.TiffWriter("03-ab.ome.tif", bigtiff=True) as tif:
    for ii, p, ss in zip(
        imgs.reshape(-1, *(100, 100)), positions, random_gaussian_sigma
    ):
        ii = skimage.filters.gaussian(ii, ss)
        ii = skimage.util.img_as_ubyte(ii)
        ii = ii.astype(np.float32)
        ii = np.array([ii, ii * 5, ii * 10, ii / 5, ii / 10])
        ii += 105
        ii = np.round(ii).astype(np.uint16)
        metadata = {
            "Pixels": {
                "PhysicalSizeX": pixel_size,
                "PhysicalSizeXUnit": "µm",
                "PhysicalSizeY": pixel_size,
                "PhysicalSizeYUnit": "µm",
            },
            "Plane": {
                "PositionX": [p[1] * pixel_size] * ii.shape[0],
                "PositionY": [p[0] * pixel_size] * ii.shape[0],
            },
        }
        tif.write(ii, metadata=metadata)
