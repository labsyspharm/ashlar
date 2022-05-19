import warnings
import sys
import re
import argparse
import pathlib
import blessed
from .. import __version__ as VERSION
from .. import reg
from ..reg import PlateReader, BioformatsReader
from ..filepattern import FilePatternReader
from ..fileseries import FileSeriesReader
from ..zen import ZenReader


def main(argv=sys.argv):

    parser = argparse.ArgumentParser(
        description='Stitch and align multi-tile cyclic microscope images',
        formatter_class=HelpFormatter,
    )
    parser.add_argument(
        'filepaths', metavar='FILE', nargs='+',
        help='Image file(s) to be processed, one per cycle',
    )
    parser.add_argument(
        '-o', '--output', dest='output', default='ashlar_output.ome.tif',
        metavar='PATH',
        help=("Output file. If PATH ends in .ome.tif a pyramidal OME-TIFF will"
              " be written. If PATH ends in just .tif and includes {cycle} and"
              " {channel} placeholders, a series of single-channel plain TIFF"
              " files will be written. If PATH starts with a relative or"
              " absolute path to another directory, that directory must already"
              " exist."),
    )
    parser.add_argument(
        '-c', '--align-channel', dest='align_channel', type=int,
        default='0', metavar='CHANNEL',
        help=('Reference channel number for image alignment. Numbering starts'
              ' at 0.'),
    )
    parser.add_argument(
        '--flip-x', default=False, action='store_true',
        help='Flip tile positions left-to-right',
    )
    parser.add_argument(
        '--flip-y', default=False, action='store_true',
        help='Flip tile positions top-to-bottom',
    )
    parser.add_argument(
        '--flip-mosaic-x', default=False, action='store_true',
        help='Flip output image left-to-right',
    )
    parser.add_argument(
        '--flip-mosaic-y', default=False, action='store_true',
        help='Flip output image top-to-bottom',
    )
    parser.add_argument(
        '--output-channels', nargs='+', type=int, metavar='CHANNEL',
        help=('Output only specified channels for each cycle. Numbering starts'
              ' at 0. (default: all channels)'),
    )
    parser.add_argument(
        '-m', '--maximum-shift', type=float, default=15, metavar='SHIFT',
        help='Maximum allowed per-tile corrective shift in microns',
    )
    parser.add_argument(
        '--filter-sigma', type=float, default=0, metavar='SIGMA',
        help=('Filter images before alignment using a Gaussian kernel with s.d.'
              ' of SIGMA pixels (default: no filtering)'),
    )
    parser.add_argument(
        '-f', '--filename-format', dest='filename_format',
        default='cycle_{cycle}_channel_{channel}.tif', help=argparse.SUPPRESS,
    )
    parser.add_argument(
        '--pyramid', default=False, action='store_true', help=argparse.SUPPRESS
    )
    parser.add_argument(
        '--tile-size', type=int, default=1024, metavar='PIXELS',
        help='Pyramid tile size for OME-TIFF output',
    )
    parser.add_argument(
        '--ffp', metavar='FILE', nargs='+',
        help=("Perform flat field illumination correction using the given"
              " profile image. Specify one common file for all cycles or one"
              " file for every cycle. Channel counts must match input files."
              " (default: no flat field correction)"),
    )
    parser.add_argument(
        '--dfp', metavar='FILE', nargs='+',
        help=("Perform dark field illumination correction using the given"
              " profile image. Specify one common file for all cycles or one"
              " file for every cycle. Channel counts must match input files."
              " (default: no dark field correction)"),
    )
    parser.add_argument(
        '--plates', default=False, action='store_true',
        help='Enable plate mode for HTS data',
    )
    parser.add_argument(
        '-q', '--quiet', dest='quiet', default=False,
        action='store_true', help='Suppress progress display',
    )
    parser.add_argument(
        '--version', action='version', version=f"ashlar {VERSION}"
    )
    args = parser.parse_args(argv[1:])

    configure_terminal()
    configure_warning_format()

    filepaths = args.filepaths

    output_path = pathlib.Path(args.output)
    op_tiff = bool(re.search(r"\.tiff?$", output_path.name, re.IGNORECASE))
    ff_default = args.filename_format == parser.get_default("filename_format")
    if op_tiff and ff_default:
        # Standard usage: -o includes a .tif filename, -f not included.
        args.filename_format = output_path.name
        output_path = output_path.parent
    else:
        # Old, deprecated usage: -o is a directory and/or -f was specified.
        if ff_default:
            warnings.warn(
                "The output path must include a filename with a .tif or .tiff"
                " suffix. Specifying only a directory path with -o/--output has"
                " been deprecated and will be disabled in a future version. See"
                " the -o documentation for details.",
                reg.Warning,
            )
        else:
            warnings.warn(
                "The -f/--filename-format argument has been deprecated and its"
                " functionality merged into the -o argument. See the -o"
                " documentation for details.",
                reg.Warning,
            )
        if op_tiff and not output_path.is_dir():
            # Checking is_dir() avoids erroring out in the strange but legal
            # situation where output_path is a DIRECTORY that ends in .tif !
            print_error(
                "Filename may be appended to the output path specified by"
                " -o/--output, or specified separately with"
                " -f/--filename-format, but not both at the same time"
            )
            return 1
        if not re.search(r"\.tiff?$", args.filename_format, re.IGNORECASE):
            print_error(
                f"Filename format does not end in .tif: {args.filename_format}"
            )
            return 1
    if not output_path.is_dir():
        print_error(
            "Output location does not exist or is not a directory:"
            f" {output_path}"
        )
        return 1
    if re.search(r"\.ome\.tiff?$", args.filename_format, re.IGNORECASE):
        args.pyramid = True

    if args.tile_size != parser.get_default("tile_size") and not args.pyramid:
        print_error("--tile-size can only be used with OME-TIFF output")
        return 1

    ffp_paths = args.ffp
    if ffp_paths:
        if len(ffp_paths) not in (0, 1, len(filepaths)):
            print_error(
                "Wrong number of flat-field profiles. Must be 1, or {}"
                " (number of input files)".format(len(filepaths))
            )
            return 1
        if len(ffp_paths) == 1:
            ffp_paths = ffp_paths * len(filepaths)

    dfp_paths = args.dfp
    if dfp_paths:
        if len(dfp_paths) not in (0, 1, len(filepaths)):
            print_error(
                "Wrong number of dark-field profiles. Must be 1, or {}"
                " (number of input files)".format(len(filepaths))
            )
            return 1
        if len(dfp_paths) == 1:
            dfp_paths = dfp_paths * len(filepaths)

    aligner_args = {}
    aligner_args['channel'] = args.align_channel
    aligner_args['verbose'] = not args.quiet
    aligner_args['max_shift'] = args.maximum_shift
    aligner_args['filter_sigma'] = args.filter_sigma

    mosaic_args = {}
    if args.output_channels:
        mosaic_args['channels'] = args.output_channels
    if args.pyramid:
        mosaic_args['tile_size'] = args.tile_size
    if args.quiet is False:
        mosaic_args['verbose'] = True
    mosaic_args['flip_mosaic_x'] = args.flip_mosaic_x
    mosaic_args['flip_mosaic_y'] = args.flip_mosaic_y

    try:
        if args.plates:
            return process_plates(
                filepaths, output_path, args.filename_format, args.flip_x,
                args.flip_y, ffp_paths, dfp_paths, aligner_args, mosaic_args,
                args.pyramid, args.quiet
            )
        else:
            mosaic_path_format = str(output_path / args.filename_format)
            return process_single(
                filepaths, mosaic_path_format, args.flip_x, args.flip_y,
                ffp_paths, dfp_paths, aligner_args, mosaic_args, args.pyramid,
                args.quiet
            )
    except ProcessingError as e:
        print_error(str(e))
        return 1


def process_single(
    filepaths, output_path_format, flip_x, flip_y, ffp_paths, dfp_paths,
    aligner_args, mosaic_args, pyramid, quiet, plate_well=None
):

    mosaic_args = mosaic_args.copy()
    writer_args = {}
    if pyramid:
        writer_args["tile_size"] = mosaic_args.pop("tile_size", None)
    mosaics = []

    if not quiet:
        print("Stitching and registering input images")
        print('Cycle 0:')
        print('    reading %s' % filepaths[0])
    reader = build_reader(filepaths[0], plate_well=plate_well)
    process_axis_flip(reader, flip_x, flip_y)
    ea_args = aligner_args.copy()
    if len(filepaths) == 1:
        ea_args['do_make_thumbnail'] = False
    edge_aligner = reg.EdgeAligner(reader, **ea_args)
    edge_aligner.run()
    mshape = edge_aligner.mosaic_shape
    mosaic_args_final = mosaic_args.copy()
    if ffp_paths:
        mosaic_args_final['ffp_path'] = ffp_paths[0]
    if dfp_paths:
        mosaic_args_final['dfp_path'] = dfp_paths[0]
    mosaics.append(reg.Mosaic(edge_aligner, mshape, **mosaic_args_final))

    for cycle, filepath in enumerate(filepaths[1:], 1):
        if not quiet:
            print('Cycle %d:' % cycle)
            print('    reading %s' % filepath)
        reader = build_reader(filepath, plate_well=plate_well)
        process_axis_flip(reader, flip_x, flip_y)
        layer_aligner = reg.LayerAligner(reader, edge_aligner, **aligner_args)
        layer_aligner.run()
        mosaic_args_final = mosaic_args.copy()
        if ffp_paths:
            mosaic_args_final['ffp_path'] = ffp_paths[cycle]
        if dfp_paths:
            mosaic_args_final['dfp_path'] = dfp_paths[cycle]
        mosaics.append(reg.Mosaic(layer_aligner, mshape, **mosaic_args_final))

    # Disable reader caching to save memory during mosaicing and writing.
    edge_aligner.reader = edge_aligner.reader.reader

    if not quiet:
        print()
        print(f"Merging tiles and writing to {output_path_format}")
    writer_class = reg.PyramidWriter if pyramid else reg.TiffListWriter
    writer = writer_class(
        mosaics, output_path_format, verbose=not quiet, **writer_args
    )
    writer.run()

    return 0


def process_plates(
    filepaths, output_path, filename_format, flip_x, flip_y, ffp_paths,
    dfp_paths, aligner_args, mosaic_args, pyramid, quiet
):

    temp_reader = build_reader(filepaths[0])
    metadata = temp_reader.metadata
    if metadata.num_plates == 0:
        # FIXME raise ProcessingError here instead?
        print("Dataset does not contain plate information.")
        return 1

    for p, plate_name in enumerate(metadata.plate_names):
        print("Plate {} ({})\n==========\n".format(p, plate_name))
        for w, well_name in enumerate(metadata.well_names[p]):
            print("Well {}\n-----".format(well_name))
            if len(metadata.plate_well_series[p][w]) > 0:
                well_path = output_path / plate_name / well_name
                well_path.mkdir(parents=True, exist_ok=True)
                mosaic_path_format = str(well_path / filename_format)
                process_single(
                    filepaths, mosaic_path_format, flip_x, flip_y,
                    ffp_paths, dfp_paths, aligner_args, mosaic_args, pyramid,
                    quiet, plate_well=(p, w)
                )
            else:
                print("Skipping -- No images found.")
            print()
        print()

    return 0


def process_axis_flip(reader, flip_x, flip_y):
    metadata = reader.metadata
    # Trigger lazy initialization.
    _ = metadata.positions
    sx = -1 if flip_x else 1
    sy = -1 if flip_y else 1
    metadata._positions *= [sy, sx]


readers = {
    'filepattern': FilePatternReader,
    'fileseries': FileSeriesReader,
    'bioformats': BioformatsReader,
    'zen': ZenReader,
}

# This is a short-term hack to provide a way to specify alternate reader
# classes and pass specific args to them.
def build_reader(path, plate_well=None):
    # Default to BioformatsReader if name not specified.
    reader_class = BioformatsReader
    kwargs = {}
    match = re.match(
        r'(?P<reader>\w+)\|(?P<path>.*?)(\|(?P<kwargs>.*))?$', path
    )
    if match:
        path = match.group('path')
        reader_name = match.group('reader')
        reader_class = readers.get(reader_name)
        if reader_class is None:
            raise ProcessingError("Unknown reader: {}".format(reader_name))
        kwargs.update(parse_kwargs_string(match.group('kwargs')))
    if plate_well is not None:
        if not issubclass(reader_class, PlateReader):
            raise ProcessingError(
                "The %s reader does not support plate/well processing"
                % reader_class.__name__
            )
        kwargs.update(plate=plate_well[0], well=plate_well[1])
    reader = reader_class(path, **kwargs)
    return reader


def parse_kwargs_string(string):
    kwargs = {}
    if string is not None:
        for piece in string.split('|'):
            name, value = piece.split('=')
            # Optimistically parse as float.
            try:
                value = float(value)
            except ValueError:
                pass
            kwargs[name] = value
    return kwargs


def configure_terminal():
    global terminal
    terminal = blessed.Terminal()


def print_error(message):
    print(terminal.bright_red("ERROR:"), message)


def warning_formatter(message, category, filename, lineno, line=None):
    if issubclass(category, reg.Warning):
        return terminal.bright_yellow("WARNING:") + f" {message}\n"
    else:
        return _old_formatwarning(message, category, filename, lineno, line)


def configure_warning_format():
    global _old_formatwarning
    _old_formatwarning = warnings.formatwarning
    warnings.formatwarning = warning_formatter


class HelpFormatter(argparse.HelpFormatter):
    """Help message formatter which adds default values to argument help.

    Based on argparse.ArgumentDefaultsHelpFormatter.
    """

    def _get_help_string(self, action):
        help = action.help
        if isinstance(action, (argparse._HelpAction, argparse._VersionAction)):
            help = help.capitalize()
        elif (
            not isinstance(action, argparse._StoreTrueAction)
            and "%(default)" not in help
            and "(default:" not in help
            and action.default is not argparse.SUPPRESS
        ):
            defaulting_nargs = [argparse.OPTIONAL, argparse.ZERO_OR_MORE]
            if action.option_strings or action.nargs in defaulting_nargs:
                help += " (default: %(default)s)"
        return help


class ProcessingError(RuntimeError):
    pass


if __name__ == '__main__':
    sys.exit(main())
