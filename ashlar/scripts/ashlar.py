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
        description='Stitch and align one or more multi-series images'
    )
    parser.add_argument(
        'filepaths', metavar='FILE', nargs='+',
        help='an image file to be processed (one file per cycle)'
    )
    parser.add_argument(
        '-o', '--output', dest='output', default='ashlar-output.ome.tif',
        metavar='PATH',
        help="write output to PATH; default is ashlar_output.ome.tif. If value"
        " ends in .ome.tif an OME-TIFF with tiled image pyramid will be written."
        " If value ends in just .tif and includes {cycle} and {channel}"
        " placeholders a series of single-channel TIFF files will be written."
        " Otherwise value will be interpreted as a directory and the '-f' and"
        " '--pyramid' arguments will control the file names and format."
    )
    parser.add_argument(
        '-c', '--align-channel', dest='align_channel', type=int,
        default='0', metavar='CHANNEL',
        help=('align images using channel number CHANNEL; numbering starts'
              ' at 0')
    )
    parser.add_argument(
        '--flip-x', default=False, action='store_true',
        help=('flip tile positions left-to-right to account for unusual'
              ' microscope configurations')
    )
    parser.add_argument(
        '--flip-y', default=False, action='store_true',
        help=('flip tile positions top-to-bottom to account for unusual'
              ' microscope configurations')
    )
    parser.add_argument(
        '--flip-mosaic-x', default=False, action='store_true',
        help=('flip output image horizontally')
    )
    parser.add_argument(
        '--flip-mosaic-y', default=False, action='store_true',
        help=('flip output image vertically')
    )
    parser.add_argument(
        '--output-channels', nargs='*', type=int, metavar='CHANNEL',
        help=('output only channels listed in CHANNELS; numbering starts at 0')
    )
    parser.add_argument(
        '-m', '--maximum-shift', type=float, default=15, metavar='SHIFT',
        help='maximum allowed per-tile corrective shift in microns'
    )
    parser.add_argument(
        '--filter-sigma', type=float, default=0.0, metavar='SIGMA',
        help=('width in pixels of Gaussian filter to apply to images before'
              ' alignment; default is 0 which disables filtering')
    )
    arg_f_default = 'cycle_{cycle}_channel_{channel}.tif'
    parser.add_argument(
        '-f', '--filename-format', dest='filename_format',
        default=arg_f_default, metavar='FORMAT',
        help="use FORMAT to generate output filenames, with {cycle} and"
        " {channel} as required placeholders for the cycle and channel"
        f" numbers; default is {arg_f_default} (DEPRECATED: Use the '-o'"
        " argument to specify the output filename format.)"
    )
    parser.add_argument(
        '--pyramid', default=False, action='store_true',
        help='write output as a single pyramidal OME-TIFF'
    )
    # Implement default-value logic ourselves so we can detect when the user
    # has explicitly set a value.
    tile_size_default = 1024
    parser.add_argument(
        '--tile-size', type=int, default=None, metavar='PIXELS',
        help=('set tile width and height to PIXELS (pyramid output only);'
              ' default is {default}'.format(default=tile_size_default))
    )
    parser.add_argument(
        '--ffp', metavar='FILE', nargs='*',
        help=('read flat field profile image from FILES; if specified must'
              ' be one common file for all cycles or one file for each cycle')
    )
    parser.add_argument(
        '--dfp', metavar='FILE', nargs='*',
        help=('read dark field profile image from FILES; if specified must'
              ' be one common file for all cycles or one file for each cycle')
    )
    parser.add_argument(
        '--plates', default=False, action='store_true',
        help='enable plate mode for HTS data'
    )
    parser.add_argument(
        '-q', '--quiet', dest='quiet', default=False, action='store_true',
        help='suppress progress display'
    )
    parser.add_argument(
        '--version', dest='version', default=False, action='store_true',
        help='print version'
    )
    args = parser.parse_args(argv[1:])

    configure_terminal()
    configure_warning_format()

    if args.version:
        print('ashlar {}'.format(VERSION))
        return 0

    if len(args.filepaths) == 0:
        parser.print_usage()
        return 1

    filepaths = args.filepaths

    output_path = pathlib.Path(args.output)
    if re.search(r"\.tiff?$", output_path.name):
        if args.filename_format != arg_f_default:
            print_error(
                "Filename may be appended to the output path specified by"
                " -o/--output, or specified separately with"
                " -f/--filename-format, but not both at the same time"
            )
            return 1
        if re.search(r"\.ome\.tiff?$", output_path.name):
            args.pyramid = True
        args.filename_format = output_path.name
        output_path = output_path.parent
    if output_path.is_dir() and not output_path.exists():
        print_error("Output directory '{}' does not exist".format(output_path))
        return 1

    if args.tile_size and not args.pyramid:
        print_error("--tile-size can only be used with --pyramid")
        return 1
    if args.tile_size is None:
        # Implement default value logic as mentioned in argparser setup above.
        args.tile_size = tile_size_default

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
    if issubclass(category, reg.DataWarning):
        return terminal.bright_yellow("WARNING:") + f" {message}\n"
    else:
        return _old_formatwarning(message, category, filename, lineno, line)


def configure_warning_format():
    global _old_formatwarning
    _old_formatwarning = warnings.formatwarning
    warnings.formatwarning = warning_formatter


class ProcessingError(RuntimeError):
    pass


if __name__ == '__main__':
    sys.exit(main())
