import pathlib
import sys

import palom
import tifffile

from . import ome_metadata


def combine_cycles(
    input_dir: str | pathlib.Path = ".",
    glob_file_pattern: str = "*/*.ashlar*.ome.tif",
    output_path: str | pathlib.Path = None,
    # if `input_files` is specified, `input_dir` and `glob_file_pattern` are
    # ignored for input file list generation
    input_files: list[str | pathlib.Path] = None,
    dna_file: str = None,
    # set `dna_channel_number` to None to keep all the channels in the input
    # files
    dna_channel_number: int = 0,
    flip_image_x: bool = True,
    flip_image_y: bool = False,
    overwrite: bool = False,
    _level: int = 0,
):
    assert dna_channel_number >= 0, "`dna_channel_number` must be >= 0"
    assert dna_file is not None, "`dna_file` is required"
    dna_file = pathlib.Path(dna_file)
    assert dna_file.exists(), f"{dna_file} does not exist"
    input_files = _get_input_files(input_dir, glob_file_pattern, input_files)
    output_path = _get_output_path(input_dir, output_path)
    if not overwrite:
        assert (
            not output_path.exists()
        ), f"output path ({output_path}) exists and `overwrite` is False"
    # prepend input files with dna file
    input_files.insert(0, dna_file)
    _input_files = "\n\t".join(
        [
            f"{p} (for DNA channel)" if idx == 0 else f"{p}"
            for idx, p in enumerate(input_files)
        ]
    )

    print_str = f"""
Processing:
\t{_input_files}

Writing to:
\t{output_path}

DNA channel number:
\t{dna_channel_number}

Flip image:
\tX: {flip_image_x}
\tY: {flip_image_y}
"""
    print(print_str)

    readers = [palom.reader.OmePyramidReader(p) for p in input_files]
    _check_xy_dimension(readers, _level)
    _check_dna_channel(readers, dna_channel_number)

    mosaics = []

    for idx, reader in enumerate(readers):
        _mosaic = reader.pyramid[_level]
        if flip_image_x:
            _mosaic = _mosaic[:, :, ::-1]
        if flip_image_y:
            _mosaic = _mosaic[:, ::-1, :]
        channels = list(range(_mosaic.shape[0]))
        channels.pop(dna_channel_number)
        if idx == 0:
            channels = [dna_channel_number]
        assert len(channels) > 0, (
            f"no channel left in {reader.path} when setting "
            f"`dna_channel_number={dna_channel_number}`"
        )
        mosaics.append(_mosaic[channels])

    try:
        tif_tags = _src_tif_tags(input_files[1])
    except Exception:
        tif_tags = {}

    palom.pyramid.write_pyramid(
        mosaics=mosaics,
        output_path=output_path,
        pixel_size=readers[1].pixel_size,
        downscale_factor=2,
        compression="zlib",
        tile_size=1024,
        save_RAM=True,
        kwargs_tifffile=tif_tags,
    )

    ome = ome_metadata._combine_metadata(
        input_files=input_files,
        output_path=output_path,
        dna_file_index=0,
        dna_channel_number=dna_channel_number,
    )
    tifffile.tiffcomment(output_path, ome.to_xml().encode())

    return 0


def _check_xy_dimension(readers, level):
    ref_img = readers[0].pyramid[level]
    ref_shape = ref_img.shape[1:]
    for rr in readers:
        _curr_shape = rr.pyramid[level].shape[1:]
        assert ref_shape == _curr_shape, (
            f"image x-y dimension does not match: "
            f"{ref_shape} - {readers[0].path} and "
            f"{_curr_shape} - {rr.path}"
        )


def _check_dna_channel(readers, dna_channel_number):
    for rr in readers:
        num_channels = rr.pyramid[0].shape[0]
        assert dna_channel_number < num_channels, (
            f"`dna_channel_number={dna_channel_number}` is out of range; "
            f"{rr.path} has {num_channels} channels"
        )


def _get_output_path(input_dir, output_path):
    if output_path is not None:
        output_path = pathlib.Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
    else:
        input_dir = pathlib.Path(input_dir).absolute()
        output_path = input_dir / f"{input_dir.name}-combined.ome.tif"
    return output_path


def _get_input_files(input_dir, glob_file_pattern, input_files):
    if input_files is not None:
        input_files = [pathlib.Path(p) for p in input_files]
        for file in input_files:
            assert file.exists(), f"input file {file} does not exist"
    else:
        input_dir = pathlib.Path(input_dir)
        assert input_dir.exists(), f"`input_dir` ({input_dir}) does not exist"
        input_files = sorted(pathlib.Path(input_dir).glob(glob_file_pattern))
        assert len(input_files) > 1, (
            f"{len(input_files)} file found, nothing to combine. "
            f"`input_dir`: {input_dir}; "
            f"`glob_file_pattern`: {glob_file_pattern}"
        )
    return input_files


def _src_tif_tags(img_path):
    kwargs_tifffile = {}
    with tifffile.TiffFile(img_path) as tif:
        kwargs_tifffile.update(
            dict(
                photometric=tif.pages[0].photometric.value,
                resolution=tif.pages[0].resolution,
                resolutionunit=tif.pages[0].resolutionunit.value,
                software=tif.pages[0].software,
            )
        )
    return kwargs_tifffile


def main():
    import fire

    fire.Fire(combine_cycles)


if __name__ == "__main__":
    sys.exit(main())

"""
NAME
    combine.py

SYNOPSIS
    combine.py <flags>

FLAGS
    --input_dir=INPUT_DIR
        Type: str | pathlib.Path
        Default: '.'
    -g, --glob_file_pattern=GLOB_FILE_PATTERN
        Type: str
        Default: '*/*.ashlar*.ome.tif'
    --output_path=OUTPUT_PATH
        Type: Optional[str | pathlib.Path]
        Default: None
    --input_files=INPUT_FILES
        Type: Optional[list]
        Default: None
    --dna_filename=DNA_FILENAME
        Type: Optional[str]
        Default: None
    --dna_channel_number=DNA_CHANNEL_NUMBER
        Type: int
        Default: 0
    --overwrite=OVERWRITE
        Type: bool
        Default: False
    -_, --_level=_LEVEL
        Type: int
        Default: 0
"""
