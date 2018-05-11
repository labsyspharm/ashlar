import sys
import argparse
try:
    import pathlib
except ImportError:
    import pathlib2 as pathlib
from .. import __version__ as VERSION
from .. import reg


def main(argv=sys.argv):
    try:
        return(main_inner(argv))
    finally:
        # FIXME This whole main wrapper is unneeded with the new jnius
        # implementation, but I'm leaving this block intact for now to avoid
        # potential rebase problems.
        pass


def main_inner(argv):

    parser = argparse.ArgumentParser(
        description='Stitch and align one or more multi-series images'
    )
    parser.add_argument(
        'filepaths', metavar='FILE', nargs='*',
        help='an image file to be processed (one file per cycle)'
    )
    parser.add_argument(
        '-o', '--output', dest='output', default='.', metavar='DIR',
        help='write output image files to DIR; default is the current directory'
    )
    parser.add_argument(
        '-c', '--align-channel', dest='align_channel', nargs='?', type=int,
        default='0', metavar='CHANNEL',
        help=('align images using channel number CHANNEL; numbering starts'
              ' at 0')
    )
    parser.add_argument(
        '--output-channels', nargs='*', type=int, metavar='CHANNEL',
        help=('output only channels listed in CHANNELS; numbering starts at 0')
    )
    parser.add_argument(
        '-m', '--maximum-shift', type=float, default=15, metavar='SHIFT',
        help='maximum allowed per-tile corrective shift in microns'
    )
    arg_f_default = 'cycle_{cycle}_channel_{channel}.tif'
    parser.add_argument(
        '-f', '--filename-format', dest='filename_format',
        default=arg_f_default, metavar='FORMAT',
        help=('use FORMAT to generate output filenames, with {{cycle}} and'
              ' {{channel}} as required placeholders for the cycle and channel'
              ' numbers; default is {default}'.format(default=arg_f_default))
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
        '-q', '--quiet', dest='quiet', default=False, action='store_true',
        help='suppress progress display'
    )
    parser.add_argument(
        '--version', dest='version', default=False, action='store_true',
        help='print version'
    )
    args = parser.parse_args(argv[1:])

    if args.version:
        print('ashlar {}'.format(VERSION))
        return 0

    if len(args.filepaths) == 0:
        parser.print_usage()
        return 1

    filepaths = args.filepaths
    if len(filepaths) == 1:
        path = pathlib.Path(filepaths[0])
        if path.is_dir():
            filepaths = sorted(str(p) for p in path.glob('*rcpnl'))

    output_path = pathlib.Path(args.output)
    if not output_path.exists():
        print("Output directory '{}' does not exist".format(output_path))
        return 1

    ffp_paths = args.ffp
    if ffp_paths:
        if len(ffp_paths) not in (0, 1, len(filepaths)):
            print("Wrong number of flat-field profiles. Must be 1, or {}"
                  " (number of input files)".format(len(filepaths)))
            return 1
        if len(ffp_paths) == 1:
            ffp_paths = ffp_paths * len(filepaths)

    dfp_paths = args.dfp
    if dfp_paths:
        if len(dfp_paths) not in (0, 1, len(filepaths)):
            print("Wrong number of dark-field profiles. Must be 1, or {}"
                  " (number of input files)".format(len(filepaths)))
            return 1
        if len(dfp_paths) == 1:
            dfp_paths = dfp_paths * len(filepaths)

    aligners = []
    mosaics = []

    if not args.quiet:
        print('Cycle 0:')
        print('    reading %s' % filepaths[0])

    aargs = {}
    aargs['channel'] = args.align_channel
    aargs['verbose'] = not args.quiet
    aargs['max_shift'] = args.maximum_shift
    reader = reg.BioformatsReader(filepaths[0])
    aligner = reg.EdgeAligner(reader, **aargs)
    aligner.run()
    mshape = aligner.mosaic_shape

    margs = {}
    if args.output_channels:
        margs['channels'] = args.output_channels
    if ffp_paths:
        margs['ffp_path'] = ffp_paths[0]
    if dfp_paths:
        margs['dfp_path'] = dfp_paths[0]
    if args.quiet is False:
        margs['verbose'] = True
    m_format = str(output_path / args.filename_format)
    mosaic = reg.Mosaic(aligner, mshape, format_cycle(m_format, 0), **margs)
    mosaic.run()

    aligners.append(aligner)
    mosaics.append(mosaic)

    for cycle, filepath in enumerate(filepaths[1:], 1):
        if not args.quiet:
            print('Cycle %d:' % cycle)
            print('    reading %s' % filepath)
        reader = reg.BioformatsReader(filepath)
        aligner = reg.LayerAligner(reader, aligners[0], **aargs)
        aligner.run()
        if ffp_paths:
            margs['ffp_path'] = ffp_paths[cycle]
        if dfp_paths:
            margs['dfp_path'] = dfp_paths[cycle]
        mosaic = reg.Mosaic(aligner, mshape, format_cycle(m_format, cycle),
                            **margs)
        mosaic.run()
        aligners.append(aligner)
        mosaics.append(mosaic)

    return 0


def format_cycle(f, cycle):
    return f.format(cycle=cycle, channel='{channel}')


if __name__ == '__main__':
    sys.exit(main())
