import json
import pathlib
import pickle

import matplotlib.pyplot as plt
import numpy as np

from .. import reg
from . import align_cycles, subtract_pyramid


def stitch(
    path, raw_endwith='pysed.ome.tif',
    channel=0, max_shift=15, alpha=0.01, max_error=None,
    filter_sigma=1.0, is_cli=True
):
    path = pathlib.Path(path)
    raws = sorted(path.glob(f"*{raw_endwith}"))
    assert len(raws) == 1
    raw = raws[0]

    out_path = path / f"{raw.stem}.ashlar.pkl"

    c1r = reg.BioformatsReader(str(raw))
    if 'rcpnl' not in raw_endwith:
        _ = c1r.metadata.positions
        c1r.metadata._positions *= [-1, 1]
    c1e = reg.EdgeAligner(
        c1r, verbose=True,
        channel=channel, max_shift=max_shift,
        alpha=alpha, max_error=max_error,
        filter_sigma=filter_sigma
    )
    c1e.run()

    reg.plot_edge_quality(c1e, img=c1e.reader.thumbnail)
    fig = plt.gcf()
    fig.suptitle(path.name, color='white')
    fig.set_size_inches(fig.get_size_inches()*2)    
    fig.tight_layout()
    fig.savefig(path / f"{raw.stem}.ashlarqc.pdf", bbox_inches='tight')
    plt.close(fig)

    reg.plot_edge_scatter(c1e)
    fig = plt.gcf()
    fig.suptitle(path.name)
    fig.set_size_inches(fig.get_size_inches()*2)    
    fig.tight_layout()
    fig.savefig(path / f"{raw.stem}.ashlarqcsctr.pdf", bbox_inches='tight')
    plt.close(fig)

    c1e.reader._cache = {}
    with open(out_path, 'wb') as f:
        pickle.dump(c1e, f)

    if is_cli:
        return 0
    return c1e


def register(
    ref_path, moving_path, raw_endwith='pysed.ome.tif',
    channel1=None, channel2=None,
    max_shift=15, filter_sigma=0.0, is_cli=True
):
    ref_path, moving_path = pathlib.Path(ref_path), pathlib.Path(moving_path)
    c1e = _load_ashlar_pkl(ref_path)
    
    if channel1 is not None:
        c1e.channel = channel1

    raws = sorted(moving_path.glob(f"*{raw_endwith}"))
    assert len(raws) == 1
    raw = raws[0]

    out_path = moving_path / f"{raw.stem}.ashlar.pkl"

    c2r = reg.BioformatsReader(str(raw))
    if 'rcpnl' not in raw_endwith:
        _ = c2r.metadata.positions
        c2r.metadata._positions *= [-1, 1]
    c21l = align_cycles.process_rotated_reader(
        c2r, c1e,
        channel=channel2, max_shift=max_shift,
        filter_sigma=filter_sigma
    )

    reg.plot_layer_quality(c21l, img=c21l.reader.thumbnail)
    fig = plt.gcf()
    fig.suptitle(raw.name)
    fig.set_size_inches(fig.get_size_inches()*2)    
    fig.tight_layout()
    fig.savefig(moving_path / f"{raw.stem}.ashlarqc.pdf", bbox_inches='tight')
    plt.close('all')

    c1e.reader._cache = {}
    with open(out_path, 'wb') as f:
        pickle.dump(c21l, f)
    
    if is_cli:
        return 0
    return c21l


def assemble(path, channels=None, is_cli=True):
    path = pathlib.Path(path)
    aligner = _load_ashlar_pkl(path)
    out_path = path / f"{aligner.from_pickle.stem}.ome.tif"
    
    mosaic = reg.Mosaic(aligner, shape=aligner.mosaic_shape, verbose=False, channels=channels)
    writer = reg.PyramidWriter([mosaic], out_path, parallel_assemble=True, verbose=True)
    writer.run()
    if is_cli:
        return 0
    return out_path


def subtract(
        bg_path, ab_path,
        fiducial_channel=0, as_float=False, is_cli=True
    ):
    bg_path, ab_path = pathlib.Path(bg_path), pathlib.Path(ab_path)
    bg_aligner = _load_ashlar_pkl(bg_path)
    ab_aligner = _load_ashlar_pkl(ab_path)
    
    mosaic_pair = [
        reg.Mosaic(aligner, bg_aligner.mosaic_shape, verbose=False)
        for aligner in (bg_aligner, ab_aligner)
    ]

    out_path = ''
    out_path = ab_path / f"{ab_aligner.from_pickle.stem}-subtracted.ome.tif"
    scaling_factors = _exposure_time_factor(bg_path, ab_path)

    sp = subtract_pyramid.SubtractPyramid(
        mosaic_pair,
        str(out_path),
        verbose=True,
        scaling_factors=scaling_factors,
        fiducial_channel=fiducial_channel,
        as_float=as_float
    )
    sp.assemble_all_channels()
    sp.run()
    sp.cleanup()
    if is_cli:
        return 0
    return out_path


def _load_ashlar_pkl(path):
    ashlar_pkls = sorted(path.glob('*.ashlar.pkl'))
    if len(ashlar_pkls) == 0:
        raise FileNotFoundError(
            f"Scan {path} has not been stitched or registered. "
            f"Run stitch/register command first."
        )
    with open(ashlar_pkls[-1], 'rb') as f:
        aligner = pickle.load(f)
    aligner.from_pickle = ashlar_pkls[-1]
    raw_path = None
    reader = aligner.reader
    if issubclass(type(aligner.reader), reg.CachingReader):
        reader = aligner.reader.reader
    if not issubclass(type(reader), reg.BioformatsReader):
        raise NotImplementedError
    raw_path = pathlib.Path(reader.path)
    alt_raw_path = path / raw_path.name
    if not raw_path.exists():
        print(f"{raw_path} does not exist.")
        if not alt_raw_path.exists():
            print(f"{alt_raw_path} does not exist.")
            raise FileNotFoundError
        reader.path = str(alt_raw_path)
    return aligner


def _exposure_time_factor(bg_path, ab_path):
    bg_jobs = sorted(bg_path.glob('*.rcjob'))
    ab_jobs = sorted(ab_path.glob('*.rcjob'))
    assert len(bg_jobs) == 1
    assert len(ab_jobs) == 1
    bg = _exposure_time_rcjob(bg_jobs[0])
    ab = _exposure_time_rcjob(ab_jobs[0])
    return np.divide(ab, bg)
        

def _exposure_time_rcjob(path):
    with open(path) as f:
        cfg = json.load(f)
    return cfg['scanner']['assay']['exposures']


def main():
    import fire
    fire.Fire({
        'stitch': stitch,
        'register': register,
        'assemble': assemble,
        'subtract': subtract
    })


