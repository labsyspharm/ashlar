from __future__ import division
import bioformats
import javabridge
import numpy as np
import scipy.ndimage
import skimage.feature
import matplotlib.pyplot as plt
from PIL import Image


def lowpass(im, sigma):
    dtype = im.dtype
    im = im.astype('f8')
    im -= scipy.ndimage.gaussian_filter(im, sigma)
    im = np.maximum(im, 0)
    im = im.astype(dtype)
    return im


javabridge.start_vm(class_path=bioformats.JARS)
# Hack module to fix py3 assumptions which break XML parsing.
bioformats.omexml.str = unicode

filename = '2_40X_LY6C_CD8_CD68_Scan_20170427_134107_01x4x00176.rcpnl'
filepath = '../../dloads/' + filename
ir = bioformats.ImageReader(filepath)
metadata = bioformats.OMEXML(bioformats.get_omexml_metadata(filepath))

#metadata.image(0).Pixels.Plane(0).PositionX

im1 = ir.read(c=0, series=93, rescale=False) # 28054.0, 14227.267
im2 = ir.read(c=0, series=94, rescale=False) # 28441.7, 14227.267
im3 = ir.read(c=0, series=82, rescale=False)
im4 = ir.read(c=0, series=83, rescale=False)
sigma = 100
im1a = lowpass(ir.read(c=1, series=93, rescale=False), sigma)
im2a = lowpass(ir.read(c=1, series=94, rescale=False), sigma)
im3a = lowpass(ir.read(c=1, series=82, rescale=False), sigma)
im4a = lowpass(ir.read(c=1, series=83, rescale=False), sigma)
im1b = lowpass(ir.read(c=2, series=93, rescale=False), sigma)
im2b = lowpass(ir.read(c=2, series=94, rescale=False), sigma)
im3b = lowpass(ir.read(c=2, series=82, rescale=False), sigma)
im4b = lowpass(ir.read(c=2, series=83, rescale=False), sigma)
im1c = lowpass(ir.read(c=3, series=93, rescale=False), sigma)
im2c = lowpass(ir.read(c=3, series=94, rescale=False), sigma)
im3c = lowpass(ir.read(c=3, series=82, rescale=False), sigma)
im4c = lowpass(ir.read(c=3, series=83, rescale=False), sigma)

im1f = skimage.filters.laplace(im1)
im2f = skimage.filters.laplace(im2)
im3f = skimage.filters.laplace(im3)
im4f = skimage.filters.laplace(im4)
shift12, _, _ = skimage.feature.register_translation(im1f, im2f, 100)
shift13, _, _ = skimage.feature.register_translation(im1f, im3f, 100)
shift24, _, _ = skimage.feature.register_translation(im2f, im4f, 100)

shift = np.vstack((shift12, shift13, shift24))
shift[2] += shift[0]
shift_f, shift_i = np.modf(shift)
shift_i = shift_i[:, ::-1].astype('i8')

m = Image.new('I;16', (im1.shape[1]*2, im1.shape[0]*2))
m.paste(Image.fromarray(im1), (0, 0))
m.paste(Image.fromarray(scipy.ndimage.shift(im2, shift_f[0])), tuple((1280, 0) + shift_i[0]))
m.paste(Image.fromarray(scipy.ndimage.shift(im3, shift_f[1])), tuple((0, 1080) + shift_i[1]))
m.paste(Image.fromarray(scipy.ndimage.shift(im4, shift_f[2])), tuple((1280, 1080) + shift_i[2]))
plt.imshow(m)
plt.show()

ma = Image.new('I;16', m.size)
ma.paste(Image.fromarray(im1a), (0, 0))
ma.paste(Image.fromarray(scipy.ndimage.shift(im2a, shift_f[0])), tuple((1280, 0) + shift_i[0]))
ma.paste(Image.fromarray(scipy.ndimage.shift(im3a, shift_f[1])), tuple((0, 1080) + shift_i[1]))
ma.paste(Image.fromarray(scipy.ndimage.shift(im4a, shift_f[2])), tuple((1280, 1080) + shift_i[2]))

mb = Image.new('I;16', m.size)
mb.paste(Image.fromarray(im1b), (0, 0))
mb.paste(Image.fromarray(scipy.ndimage.shift(im2b, shift_f[0])), tuple((1280, 0) + shift_i[0]))
mb.paste(Image.fromarray(scipy.ndimage.shift(im3b, shift_f[1])), tuple((0, 1080) + shift_i[1]))
mb.paste(Image.fromarray(scipy.ndimage.shift(im4b, shift_f[2])), tuple((1280, 1080) + shift_i[2]))

mc = Image.new('I;16', m.size)
mc.paste(Image.fromarray(im1c), (0, 0))
mc.paste(Image.fromarray(scipy.ndimage.shift(im2c, shift_f[0])), tuple((1280, 0) + shift_i[0]))
mc.paste(Image.fromarray(scipy.ndimage.shift(im3c, shift_f[1])), tuple((0, 1080) + shift_i[1]))
mc.paste(Image.fromarray(scipy.ndimage.shift(im4c, shift_f[2])), tuple((1280, 1080) + shift_i[2]))

mm = np.dstack((ma, mb, mc)).astype('f8')
mm_max = mm.max(0).max(0)
mm_min = mm.min(0).min(0)
mm_norm = (((mm - mm_min) / (mm_max - mm_min)) ** (1/2.2) * 255).astype('u1')
mrgb = Image.fromarray(mm_norm)
plt.figure()
plt.imshow(mrgb)
