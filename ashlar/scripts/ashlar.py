import sys
import argparse
try:
    import pathlib
except ImportError:
    import pathlib2 as pathlib
from .. import reg


def main(argv=sys.argv):

    parser = argparse.ArgumentParser(
        description='Stitch and align one or more multi-series images'
    )
    parser.add_argument(
        'filepaths', metavar='FILE', nargs='+',
        help='an image file to be processed'
    )
    parser.add_argument(
        '-q', '--quiet', dest='quiet', default=False, action='store_true',
        help='suppress progress display'
    )
    parser.add_argument(
        '--ffp', metavar='FFP_FILE',  help='path to flat field profile image'
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

    aligners = []
    mosaics = []

    if not args.quiet:
        print('Scan 0:')
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
    mosaic = reg.Mosaic(aligner, mshape, 'scan_0_%(channel)d.tif', **margs)
    mosaic.run()
    aligners.append(aligner)
    mosaics.append(mosaic)

    for scan, filepath in enumerate(filepaths[1:], 1):
        if not args.quiet:
            print('Scan %d:' % scan)
            print('    reading %s' % filepath)
        reader = reg.BioformatsReader(filepath)
        aligner = reg.LayerAligner(reader, aligners[0],
                                   verbose=not args.quiet)
        aligner.run()
        filename_format = 'scan_%d_%%(channel)d.tif' % scan
        mosaic = reg.Mosaic(aligner, mshape, filename_format, **margs)
        mosaic.run()
        aligners.append(aligner)
        mosaics.append(mosaic)


    try:
        __IPYTHON__
    except:
        reg._deinit_bioformats()


if __name__ == '__main__':
    main()
