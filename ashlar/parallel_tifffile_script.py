from ashlar import reg
import pathlib
import joblib


ometiffs = sorted(pathlib.Path('.').glob('*.ome.tiff'))

c1r = reg.BioformatsReader(str(ometiffs.pop(0)))
c1e = reg.EdgeAligner(c1r, max_shift=30, verbose=True)
c1e.run()


def wrap(reader_path, reference_aligner):
    reader = reg.BioformatsReader(str(reader_path))
    cx1l = reg.LayerAligner(reader, reference_aligner, verbose=True)
    cx1l.run()
    return cx1l

aligned = joblib.Parallel(n_jobs=-1, verbose=1)(joblib.delayed(wrap)(o, c1e) for o in ometiffs)

aligned.insert(0, c1e)

mosaics = [
    reg.Mosaic(
        aligner, c1e.mosaic_shape, 'zzzzz.ome.tif', channels=[0],
        combined=True, tile_size=(1024, 1024), verbose=True
    ) 
    for aligner in aligned
]

reg.write_pyramid(
    mosaics[:2]
)
