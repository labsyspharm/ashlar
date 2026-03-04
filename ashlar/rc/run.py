import pathlib
import pickle
import sys

import matplotlib.pyplot as plt
import numpy as np
import tifffile

from .. import reg, utils
from .. import __version__ as VERSION
from . import align_cycles, ome_metadata


def stitch(
    path: str | pathlib.Path,
    raw_endwith: str = "pysed.ome.tif",
    channel: int = 0,
    max_shift: float = 30,
    alpha: float = 0.01,
    max_error: float | None = None,
    filter_sigma: float = 1.0,
    is_cli: bool = True,
):
    path = pathlib.Path(path).absolute()
    raws = sorted(path.glob(f"*{raw_endwith}"))
    assert len(raws) == 1
    raw = raws[0]

    output_path = path / f"{raw.stem}.ashlar.pkl"

    c1r = reg.BioformatsReader(str(raw))
    if "rcpnl" not in raw_endwith:
        _ = c1r.metadata.positions
        c1r.metadata._positions *= [-1, 1]
    c1e = reg.EdgeAligner(
        c1r,
        verbose=True,
        channel=channel,
        max_shift=max_shift,
        alpha=alpha,
        max_error=max_error,
        filter_sigma=filter_sigma,
    )
    c1e.run()

    reg.plot_edge_quality(
        c1e,
        img=utils.visualize_long_tail_image(c1e.reader.thumbnail),
        im_kwargs=dict(
            interpolation="nearest",
            vmin=np.percentile(
                utils.visualize_long_tail_image(c1e.reader.thumbnail), 25
            ),
        ),
    )
    fig = plt.gcf()
    fig.suptitle(path.name, color="white")
    fig.set_size_inches(fig.get_size_inches() * 4)
    fig.tight_layout()
    fig.savefig(path / f"{raw.stem}.ashlarqc.pdf", bbox_inches="tight")
    plt.close(fig)

    reg.plot_edge_scatter(c1e)
    fig = plt.gcf()
    fig.suptitle(path.name)
    fig.set_size_inches(fig.get_size_inches() * 2)
    fig.tight_layout()
    fig.savefig(path / f"{raw.stem}.ashlarqcsctr.pdf", bbox_inches="tight")
    plt.close(fig)

    c1e.reader._cache = {}
    with open(output_path, "wb") as f:
        pickle.dump(c1e, f)

    if is_cli:
        return 0
    return c1e


def register(
    ref_path: str | pathlib.Path,
    moving_path: str | pathlib.Path,
    raw_endwith: str = "pysed.ome.tif",
    channel_ref: int | None = None,
    channel_moving: int | None = None,
    max_shift: float = 30,
    filter_sigma: float = 1.0,
    is_cli: bool = True,
):
    ref_path = pathlib.Path(ref_path).absolute()
    moving_path = pathlib.Path(moving_path).absolute()
    c1e = _load_ashlar_pkl(ref_path)

    if isinstance(c1e, reg.LayerAligner):
        ref_edgealigner_path = _get_ref_path(c1e)
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

    output_path = moving_path / f"{raw.stem}.ashlar.pkl"

    c2r = reg.BioformatsReader(str(raw))
    if "rcpnl" not in raw_endwith:
        _ = c2r.metadata.positions
        c2r.metadata._positions *= [-1, 1]
    c21l = align_cycles.process_rotated_reader(
        c2r, c1e, channel=channel_moving, max_shift=max_shift, filter_sigma=filter_sigma
    )

    reg.plot_layer_quality(
        c21l,
        img=utils.visualize_long_tail_image(c21l.reader.thumbnail),
        im_kwargs=dict(
            interpolation="nearest",
            vmin=np.percentile(
                utils.visualize_long_tail_image(c21l.reader.thumbnail), 25
            ),
        ),
    )
    fig = plt.gcf()
    fig.suptitle(raw.name)
    fig.set_size_inches(fig.get_size_inches() * 2)
    fig.tight_layout()
    fig.savefig(moving_path / f"{raw.stem}.ashlarqc.pdf", bbox_inches="tight")
    plt.close("all")

    c1e.reader._cache = {}
    with open(output_path, "wb") as f:
        pickle.dump(c21l, f)

    if is_cli:
        return 0
    return c21l


def assemble(
    path: str | pathlib.Path,
    output_path: str | pathlib.Path = None,
    channels: list[int] | None = None,
    is_cli: bool = True,
):
    path = pathlib.Path(path).absolute()
    aligner = _load_ashlar_pkl(path)
    if output_path is None:
        output_path = path / f"{aligner.from_pickle.stem}.ome.tif"
    else:
        output_path = _custom_output_path(output_path)
    mosaic = reg.Mosaic(
        aligner, shape=aligner.mosaic_shape, verbose=False, channels=channels
    )
    writer = reg.PyramidWriter(
        [mosaic], output_path, parallel_assemble=True, verbose=True
    )
    writer.run()

    ref_path = _get_ref_path(aligner)
    ome = ome_metadata._assemble_metadata(
        ref_path, _get_aligner_reader_path(aligner), output_path, channels
    )
    tifffile.tiffcomment(output_path, ome.to_xml().encode())
    if is_cli:
        return 0
    return output_path


def _load_ashlar_pkl(path):
    ashlar_pkls = sorted(path.glob("*.ashlar.pkl"))
    if len(ashlar_pkls) == 0:
        raise FileNotFoundError(
            f"Scan {path} has not been stitched or registered. "
            f"Run stitch/register command first."
        )
    with open(ashlar_pkls[-1], "rb") as f:
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


def _custom_output_path(output_path):
    output_path = pathlib.Path(output_path)
    assert output_path.name.endswith(
        ".ome.tif"
    ), f"`output_path` must be a file path ends with .ome.tif; not {output_path}"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    return output_path


def _get_ref_path(aligner):
    if not isinstance(aligner, reg.LayerAligner):
        return None
    ref = aligner.reference_aligner
    if isinstance(ref.reader, reg.CachingReader):
        return ref.reader.reader.path
    return ref.reader.path


def _get_aligner_reader_path(aligner):
    if not isinstance(aligner, (reg.EdgeAligner, reg.LayerAligner)):
        return
    reader = aligner.reader
    if isinstance(reader, reg.CachingReader):
        return reader.reader.path
    return reader.path


def print_version():
    print(VERSION)
    return 0


def main():
    def combine_cycles():
        return

    if len(sys.argv) > 1:
        if sys.argv[1] == "combine":
            from .combine import combine_cycles

    import fire

    fire.Fire(
        {
            "stitch": stitch,
            "register": register,
            "assemble": assemble,
            "combine": combine_cycles,
            "version": print_version,
        }
    )


if __name__ == "__main__":
    sys.exit(main())
