import argparse
import pathlib
import pickle
import sys

import lxml.etree
import matplotlib.pyplot as plt
import numpy as np
import ome_types
import tifffile
from ome_types import _conversion

from .. import __version__ as VERSION
from .. import reg, utils
from . import align_cycles, run


# ---------------------------------------------------------------------------- #
#                              single-command func                             #
# ---------------------------------------------------------------------------- #
def run_orion(
    paths: list[str | pathlib.Path],
    output_path: str | pathlib.Path = None,
    channel: int = 0,
    max_shift: float = 30,
    filter_sigma: float = 1.0,
    output_channels: list[int] | None = None,
    alpha: float = 0.01,
    max_error: float | None = None,
    n_jobs: int = 10,
    no_mask_background: bool = False,
):

    paths = [pathlib.Path(pp) for pp in paths]
    for pp in paths:
        assert pp.exists(), f"{pp} does not exist"
    output_path = pathlib.Path(output_path)

    assert output_path.name.endswith(".ome.tif")
    output_path.parent.mkdir(exist_ok=True, parents=True)

    ref, *movings = paths

    c1e = run.stitch(
        path=ref.parent,
        raw_endwith=ref.name,
        channel=channel,
        max_shift=max_shift,
        alpha=alpha,
        max_error=max_error,
        filter_sigma=filter_sigma,
        is_cli=False,
    )

    aligners = [c1e]

    for mm in movings:
        raw = mm.absolute()

        pickle_path = raw.parent / f"{raw.stem}.ashlar.pkl"

        c2r = reg.BioformatsReader(str(raw))
        if "rcpnl" not in raw.name:
            _ = c2r.metadata.positions
            c2r.metadata._positions *= [-1, 1]
        c21l = align_cycles.process_rotated_reader(
            c2r, c1e, channel=channel, max_shift=max_shift, filter_sigma=filter_sigma
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
        fig.savefig(raw.parent / f"{raw.stem}.ashlarqc.pdf", bbox_inches="tight")
        plt.close("all")

        c1e.reader._cache = {}
        with open(pickle_path, "wb") as f:
            pickle.dump(c21l, f)

        aligners.append(c21l)

    mosaic_shape = c1e.mosaic_shape
    mosaics = [
        reg.Mosaic(aa, shape=mosaic_shape, verbose=False, channels=output_channels)
        for aa in aligners
    ]

    do_mask_tissue = not no_mask_background
    writer = reg.PyramidWriter(
        mosaics,
        output_path,
        parallel_assemble=True,
        do_mask_tissue=do_mask_tissue,
        n_jobs=n_jobs,
        verbose=True,
    )
    writer.run()

    # ----------------------------- add channel name ----------------------------- #
    channel_names = [_channel_name_from_tif(pp) for pp in paths]

    names = []
    for nn in channel_names:
        if output_channels is None:
            names.extend(nn)
            continue
        names.extend(list(np.asarray(nn)[output_channels]))

    _add_channel_name(output_path, names)
    return


def _channel_name_from_tif(img_path):
    img_path = pathlib.Path(img_path)
    assert img_path.name.endswith(".pysed.ome.tif")

    xml = _conversion.tiff2xml(img_path)
    root = lxml.etree.fromstring(xml)
    image_el = root.find("{*}Image")

    if image_el is None:
        return

    ome_image = ome_types.from_xml(lxml.etree.tostring(image_el))
    return [cc.name for cc in ome_image.pixels.channels]


def _add_channel_name(tiff_path, channel_names):
    ome = ome_types.from_tiff(tiff_path)
    n_channels = len(ome.images[0].pixels.channels)
    n_names = len(channel_names)
    assert n_channels == n_names, (
        f"Number of channels ({n_channels}) in '{tiff_path}' does not match number of channel names ({n_names})."
    )

    for channel, name in zip(ome.images[0].pixels.channels, channel_names):
        channel.name = name

    tifffile.tiffcomment(tiff_path, ome.to_xml().encode())
    return


def main(argv=sys.argv):
    parser = argparse.ArgumentParser(
        description="Align and assemble Orion images",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {VERSION}",
    )
    parser.add_argument(
        "paths",
        nargs="+",
        type=pathlib.Path,
        help="Input file paths (one or more).",
    )
    parser.add_argument(
        "--output-channels",
        type=int,
        nargs="+",
        default=None,
        help="Optional list of channels to include in the output file.",
    )
    parser.add_argument(
        "--channel",
        type=int,
        default=0,
        help="Channel index to process.",
    )
    parser.add_argument(
        "--max-shift",
        type=float,
        default=30.0,
        help="Maximum allowed shift.",
    )
    parser.add_argument(
        "--filter-sigma",
        type=float,
        default=1.0,
        help="Gaussian filter sigma.",
    )
    parser.add_argument(
        "--output-path",
        type=pathlib.Path,
        default=None,
        help="Path to save output results. Required if multiple inputs.",
    )
    parser.add_argument(
        "--alpha",
        type=float,
        default=0.01,
        help="Alpha parameter.",
    )
    parser.add_argument(
        "--max-error",
        type=float,
        default=None,
        help="Maximum error tolerance.",
    )
    parser.add_argument(
        "--n-jobs",
        type=int,
        default=10,
        help="Number of parallel jobs.",
    )
    parser.add_argument(
        "--no-mask-background",
        action="store_true",
        help="Do not automatically mask out background region",
    )

    args = parser.parse_args(argv[1:])

    # Handle output path logic
    if len(args.paths) > 1 and args.output_path is None:
        parser.error(
            "When providing multiple input files, you must specify --output-path."
        )
    elif len(args.paths) == 1 and args.output_path is None:
        # Auto-generate output path: <input>-orion.ome.tif
        input_file = args.paths[0]
        args.output_path = input_file.with_name(input_file.stem + "-orion.ome.tif")

    run_orion(
        paths=args.paths,
        output_path=args.output_path,
        channel=args.channel,
        max_shift=args.max_shift,
        filter_sigma=args.filter_sigma,
        output_channels=args.output_channels,
        alpha=args.alpha,
        max_error=args.max_error,
        n_jobs=args.n_jobs,
        no_mask_background=args.no_mask_background,
    )


if __name__ == "__main__":
    sys.exit(main())
