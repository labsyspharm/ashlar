import json
import pathlib
import pickle
import sys
from typing import Iterable

import matplotlib.pyplot as plt
import numpy as np
import tifffile

from .. import reg
from . import align_cycles, subtract_pyramid, ome_metadata


def stitch(
    path: str | pathlib.Path,
    raw_endwith: str = "pysed.ome.tif",
    channel: int = 0,
    max_shift: float = 15,
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

    reg.plot_edge_quality(c1e, img=c1e.reader.thumbnail)
    fig = plt.gcf()
    fig.suptitle(path.name, color="white")
    fig.set_size_inches(fig.get_size_inches() * 2)
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
    max_shift: float = 15,
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

    reg.plot_layer_quality(c21l, img=c21l.reader.thumbnail)
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


def subtract(
    bg_path: str | pathlib.Path,
    ab_path: str | pathlib.Path,
    output_path: str | pathlib.Path = None,
    fiducial_channel: int = 0,
    bg_intensity_scaling_factor: str | Iterable[float] | None = "rcjob",
    channel_matching_by: str = "rcjob",
    camera_bias: float = 105.0,
    add_camera_bias_back: bool = False,
    as_float: bool = False,
    is_cli: bool = True,
):
    assert channel_matching_by in ["rcjob", "index"]

    bg_path = pathlib.Path(bg_path).absolute()
    ab_path = pathlib.Path(ab_path).absolute()
    bg_aligner = _load_ashlar_pkl(bg_path)
    ab_aligner = _load_ashlar_pkl(ab_path)

    bg_mosaic, ab_mosaic = [
        reg.Mosaic(aligner, bg_aligner.mosaic_shape, verbose=False)
        for aligner in (bg_aligner, ab_aligner)
    ]

    if output_path is None:
        output_path = ab_path / f"{ab_aligner.from_pickle.stem}-subtracted.ome.tif"
    else:
        output_path = _custom_output_path(output_path)

    if isinstance(bg_intensity_scaling_factor, str):
        assert bg_intensity_scaling_factor == "rcjob"
        bg_intensity_scaling_factor = _exposure_time_factor(bg_path, ab_path)

    if channel_matching_by == "index":
        # match bg and ab channels by their indecies; number of channels must match
        assert bg_mosaic.channels == ab_mosaic.channels
        assert len(bg_intensity_scaling_factor) == len(ab_mosaic.channels)
        subtraction_config = _subtraction_config_by_index(
            ab_mosaic.channels, bg_intensity_scaling_factor
        )
    else:
        # use excitation and emission info in rcjob file to map channel between
        # bg and ab cycles; do not subtract an ab channel if no matching channel
        # is present in the bg cycle
        subtraction_config = _subtraction_config_by_rcjob(bg_path, ab_path)
        bg_intensity_scaling_factor = None

    sp = subtract_pyramid.SubtractPyramid(
        bg_mosaic,
        ab_mosaic,
        str(output_path),
        verbose=True,
        bg_intensity_scaling_factor=bg_intensity_scaling_factor,
        camera_bias=camera_bias,
        add_camera_bias_back=add_camera_bias_back,
        fiducial_channel=fiducial_channel,
        as_float=as_float,
        subtraction_config=subtraction_config,
    )
    sp.assemble_all_channels()
    sp.run()
    sp.cleanup()

    ome = ome_metadata._subtract_metadata(
        bg_path=_get_aligner_reader_path(bg_aligner),
        ab_path=_get_aligner_reader_path(ab_aligner),
        bg_ref_path=_get_ref_path(bg_aligner),
        ab_ref_path=_get_ref_path(ab_aligner),
        output_path=output_path,
        fiducial_channel=fiducial_channel,
        subtraction_config=subtraction_config,
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


def _exposure_time_factor(bg_path, ab_path):
    bg_jobs = sorted(bg_path.glob("*.rcjob"))
    ab_jobs = sorted(ab_path.glob("*.rcjob"))
    assert len(bg_jobs) == 1
    assert len(ab_jobs) == 1
    bg = _exposure_time_rcjob(bg_jobs[0])
    ab = _exposure_time_rcjob(ab_jobs[0])
    if len(ab) != len(bg):
        print(
            "\n"
            f"Number of channels does not match\n"
            f"\t{len(bg)} channels in {bg_path}\n"
            f"\t{len(ab)} channels in {ab_path}\n"
        )
        return None
    return np.divide(ab, bg)


def _exposure_time_rcjob(path):
    with open(path) as f:
        cfg = json.load(f)
    return cfg["scanner"]["assay"]["exposures"]


def _subtraction_config_by_index(channels, bg_intensity_scaling_factor):
    subtraction_config = []
    for cc, ff in zip(channels, bg_intensity_scaling_factor):
        subtraction_config.append(
            {
                "channel_index": cc,
                "exposure_time": ff,
                # artificially set bg channel exposure time to 1 - the intensity
                # scaling factor is already calculated
                "bg_channel": {"channel_index": cc, "exposure_time": 1},
            }
        )
    return subtraction_config


def _subtraction_config_by_rcjob(bg_path, ab_path):
    bg_jobs = sorted(bg_path.glob("*.rcjob"))
    ab_jobs = sorted(ab_path.glob("*.rcjob"))
    assert len(bg_jobs) == 1
    assert len(ab_jobs) == 1
    job_bg = bg_jobs[0]
    job_ab = ab_jobs[0]
    settings_bg = _channel_settings_rcjob(job_bg)
    settings_ab = _channel_settings_rcjob(job_ab)

    for ab in settings_ab:
        ab["bg_channel"] = {}
        for bg in settings_bg:
            if ab["ex_em"] == bg["ex_em"]:
                ab["bg_channel"] = bg
                # break at first match
                break
    return settings_ab


def _channel_settings_rcjob(path):
    with open(path) as f:
        cfg = json.load(f)
    assay = cfg["scanner"]["assay"]

    settings = []
    for idx, (ex, em, et) in enumerate(
        zip(assay["excitations"], assay["emissions"], assay["exposures"])
    ):
        settings.append({"channel_index": idx, "ex_em": (ex, em), "exposure_time": et})
    return settings


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
            "subtract": subtract,
            "combine": combine_cycles,
        }
    )


if __name__ == "__main__":
    sys.exit(main())
