import json
import pathlib
import pickle
import sys
from typing import Iterable

import matplotlib.pyplot as plt
import numpy as np

from .. import reg
from . import align_cycles, subtract_pyramid


def stitch(
    path: str | pathlib.Path,
    raw_endwith: str = 'pysed.ome.tif',
    channel: int = 0,
    max_shift: float = 15,
    alpha: float = 0.01,
    max_error: float | None = None,
    filter_sigma: float = 1.0,
    is_cli: bool = True
):
    path = pathlib.Path(path).absolute()
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
    ref_path: str | pathlib.Path,
    moving_path: str | pathlib.Path,
    raw_endwith: str = 'pysed.ome.tif',
    channel_ref: int | None = None,
    channel_moving: int | None = None,
    max_shift: float = 15,
    filter_sigma: float = 1.0,
    is_cli: bool = True
):
    ref_path = pathlib.Path(ref_path).absolute()
    moving_path = pathlib.Path(moving_path).absolute()
    c1e = _load_ashlar_pkl(ref_path)
    
    if issubclass(type(c1e), reg.LayerAligner):
        ref_edgealigner_path = c1e.reference_aligner.reader.reader.path
        ref_edgealigner_path = pathlib.Path(ref_edgealigner_path).parent
        reg.warn_data(
            "\n\n"
            f"Image in {ref_path} was aligned to {ref_edgealigner_path}.\n"
            f"Using {ref_path} as reference propogates errors.\n"
            f"Consider using {ref_edgealigner_path} as reference instead.\n"
        )
        c1e.lr = c1e.reference_aligner.lr
    if channel_ref is not None:
        c1e.channel = channel_ref

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
        channel=channel_moving, max_shift=max_shift,
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


def assemble(
    path: str | pathlib.Path,
    channels: list[int] | None = None,
    is_cli: bool = True
):
    path = pathlib.Path(path).absolute()
    aligner = _load_ashlar_pkl(path)
    out_path = path / f"{aligner.from_pickle.stem}.ome.tif"
    
    mosaic = reg.Mosaic(aligner, shape=aligner.mosaic_shape, verbose=False, channels=channels)
    writer = reg.PyramidWriter([mosaic], out_path, parallel_assemble=True, verbose=True)
    writer.run()
    if is_cli:
        return 0
    return out_path


def subtract(
    bg_path: str | pathlib.Path,
    ab_path: str | pathlib.Path,
    fiducial_channel: int = 0,
    bg_intensity_scaling_factor: str | Iterable[int] | None = 'rcjob',
    as_float: bool = False,
    is_cli: bool = True
):
    bg_path = pathlib.Path(bg_path).absolute()
    ab_path = pathlib.Path(ab_path).absolute()
    bg_aligner = _load_ashlar_pkl(bg_path)
    ab_aligner = _load_ashlar_pkl(ab_path)
    
    bg_mosaic, ab_mosaic = [
        reg.Mosaic(aligner, bg_aligner.mosaic_shape, verbose=False)
        for aligner in (bg_aligner, ab_aligner)
    ]

    out_path = ab_path / f"{ab_aligner.from_pickle.stem}-subtracted.ome.tif"
    
    if type(bg_intensity_scaling_factor) == str:
        assert bg_intensity_scaling_factor == 'rcjob'
        bg_intensity_scaling_factor = _exposure_time_factor(bg_path, ab_path)

    sp = subtract_pyramid.SubtractPyramid(
        bg_mosaic, ab_mosaic,
        str(out_path),
        verbose=True,
        bg_intensity_scaling_factor=bg_intensity_scaling_factor,
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
    combine_cycles = None
    if sys.argv[1] == 'combine':
        from .combine import combine_cycles

    import fire
    fire.Fire({
        'stitch': stitch,
        'register': register,
        'assemble': assemble,
        'subtract': subtract,
        'combine': combine_cycles
    })


if __name__ == '__main__':
    sys.exit(main())