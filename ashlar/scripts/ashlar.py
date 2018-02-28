import sys
import argparse
try:
    import pathlib
except ImportError:
    import pathlib2 as pathlib
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
        'filepaths', metavar='FILE', nargs='+',
        help='an image file to be processed (one file per cycle)'
    )
    parser.add_argument(
        '-o', '--output', dest='output', default='.', metavar='DIR',
        help='write output image files to DIR; default is the current directory'
    )
    arg_f_default = 'cycle_{cycle}_channel_{channel}.tif'
    parser.add_argument(
        '-f', '--output-format', dest='output_format',
        default=arg_f_default, metavar='FORMAT',
        help=('use FORMAT to generate output filenames, with {{cycle}} and'
              ' {{channel}} as required placeholders for the cycle and channel'
              ' numbers; default is {default}'.format(default=arg_f_default)),
    )
    parser.add_argument(
        '--ffp', metavar='FILE',
        help='read flat field profile image from FILE'
    )
    parser.add_argument(
        '-q', '--quiet', dest='quiet', default=False, action='store_true',
        help='suppress progress display'
    )
    args = parser.parse_args(argv[1:])

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

    aligners = []
    mosaics = []

    if not args.quiet:
        print('Cycle 0:')
        print('    reading %s' % filepaths[0])
    reader = reg.BioformatsReader(filepaths[0])
    aligner = reg.EdgeAligner(reader, verbose=not args.quiet)
    aligner.run()
    mshape = aligner.mosaic_shape
    margs = {}
    if args.ffp:
        margs['ffp_path'] = args.ffp
    if args.quiet is False:
        margs['verbose'] = True

    m_format = str(output_path / args.output_format)
    mosaic = reg.Mosaic(aligner, mshape, format_cycle(m_format, 0), **margs)
    mosaic.run()
    aligners.append(aligner)
    mosaics.append(mosaic)

    for cycle, filepath in enumerate(filepaths[1:], 1):
        if not args.quiet:
            print('Cycle %d:' % cycle)
            print('    reading %s' % filepath)
        reader = reg.BioformatsReader(filepath)
        aligner = reg.LayerAligner(reader, aligners[0],
                                   verbose=not args.quiet)
        aligner.run()
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
