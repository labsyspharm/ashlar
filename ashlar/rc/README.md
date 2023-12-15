# Running ASHLAR reactively

## Install

If you don't already have [miniconda](https://docs.conda.io/en/latest/miniconda.html)
or [Anaconda](https://www.anaconda.com/products/individual), download the python
3.x version and install. Then, run the following commands from a terminal (Linux/Mac)
or command prompt (Windows):

Create a named conda environment with python 3.10:

```bash
conda create -y -n rcashlar python=3.10 openjdk pyqt
```

Activate the conda environment and pip install aslar from GitHub:

```bash
conda activate rcashlar
python -m pip install "ashlar @ git+https://github.com/yu-anchen/ashlar@afr-2023-12-1"
```

## Workflow

1. 1-target-folder mode
    - pickle file (`.ashlar.pkl`) not present (❌)
        - `stitch` raw image, will generate:
            - 1 pickle file (`.ashlar.pkl`)
            - 2 QC plot PDFs
    - pickle file (`.ashlar.pkl`) present (✅)
        - `assemble` mosaic, will generate:
            - 1 pyramidal image (`.ashlar.ome.tif`)

2. 2-target-folder mode
    - pickle files (`.ashlar.pkl`): ✅ ❌
        - `register` the second to the first, will generate
            - 1 pickle file (`.ashlar.pkl`) for the second folder
            - 1 QC plot PDF
    - pickle files (`.ashlar.pkl`): ✅ ✅
        - `subtract` the first mosaic (background/autofluorescence image) from
          the second mosaic (antibody image), will generate
            - 1 pyramidal image in the second folder (`.ashlar-subtracted.ome.tif`)

3. (Optional) Combine multiple OME-TIFFs into one single OME-TIFF

    ```text
    NAME
        rcashlar combine

    SYNOPSIS
        rcashlar combine <flags>

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
    ```

    __Parameter details:__
    - Use `--input_dir` and `--glob_file_pattern` to search and sort
      automatically in lexicographical order.
    - If `--output_path` is not provided, the default output path is
      `{input-dir}/{name-of-input-dir}-combined.ome.tif`.
    - Use `--input_files` to customize files and file ordering to combine. In
      the format of:

      ```bash
      --input_files "[r'path/to/file1.ome.tif', r'path/to/file2.ome.tif']"
      ```

    - Use `--dna_filename` to select which input file get to keep its DNA
      channel in the output file. __Provide just the "file name" not "file
      path".__ If not specified, the first file in the file list will be used.
    - Use `--dna_channel_number` to specify DNA channel location; default is 0
      (first channel). Only one DNA channel image will be written to the
      combined OME-TIFF.

## Example scenario and commands

NOTE: The following commands are executed in the `rcashlar` conda environment,
which can be activated by running

```bash
conda activate rcashlar
```

1. Slide being imaged for the first time, this scan serves as a reference, and all the
   successive scans will be registered to it.
   - Run `pysed`
   - Run `rcashlar stitch`

        ```bash
        rcashlar stitch /project/demo/21-2D474_001-2D01@20230523_173947_224759
        ```

   - (Optional) run `rcashlar assemble`

        ```bash
        rcashlar assemble /project/demo/21-2D474_001-2D01@20230523_173947_224759
        ```

2. Slide being imaged for the second time. This scan does not have antibody
   stains. It serves as background and autofluorescence.
   - Run `pysed`
   - Run `rcashlar register` to the first scan

        ```bash
        rcashlar register /project/demo/21-2D474_001-2D01@20230523_173947_224759 /project/demo/21-2D474_001-2D01@20230523_180743_729614
        ```

3. Slide stained with antibodies and being imaged for the third time.
    - Run `pysed`
    - Run `rcashlar register` to the first scan

        ```bash
        rcashlar register /project/demo/21-2D474_001-2D01@20230523_173947_224759 /project/demo/21-2D474_001-2D01@20230523_183243_352969
        ```

    - Run `rcashlar subtract` against the second scan

        ```bash
        rcashlar subtract /project/demo/230523_3scans/21-2D474_001-2D01@20230523_180743_729614 /project/demo/230523_3scans/21-2D474_001-2D01@20230523_183243_352969
        ```

4. Combine the assembled images from the first and the subtracted image into a
   single OME-TIFF.
    - Run `rcashlar combine` on the `demo/` directory

        ```bash
        rcashlar combine --input_dir /project/demo
        ```

## Command reference

1. Stitch

    ```text
    NAME
        rcashlar stitch

    SYNOPSIS
        rcashlar stitch PATH <flags>

    POSITIONAL ARGUMENTS
        PATH
            Type: str | pathlib.Path

    FLAGS
        -r, --raw_endwith=RAW_ENDWITH
            Type: str
            Default: 'pysed.ome.tif'
        -c, --channel=CHANNEL
            Type: int
            Default: 0
        --max_shift=MAX_SHIFT
            Type: float
            Default: 15
        -a, --alpha=ALPHA
            Type: float
            Default: 0.01
        --max_error=MAX_ERROR
            Type: Optional[float | None]
            Default: None
        -f, --filter_sigma=FILTER_SIGMA
            Type: float
            Default: 1.0
    ```

1. Register

    ```text
    NAME
        rcashlar register

    SYNOPSIS
        rcashlar register REF_PATH MOVING_PATH <flags>

    POSITIONAL ARGUMENTS
        REF_PATH
            Type: str | pathlib.Path
        MOVING_PATH
            Type: str | pathlib.Path

    FLAGS
        -r, --raw_endwith=RAW_ENDWITH
            Type: str
            Default: 'pysed.ome.tif'
        --channel_ref=CHANNEL_REF
            Type: Optional[int | None]
            Default: None
        --channel_moving=CHANNEL_MOVING
            Type: Optional[int | None]
            Default: None
        -m, --max_shift=MAX_SHIFT
            Type: float
            Default: 15
        -f, --filter_sigma=FILTER_SIGMA
            Type: float
            Default: 1.0
    ```

1. Assemble

    ```text
    NAME
        rcashlar assemble

    SYNOPSIS
        rcashlar assemble PATH <flags>

    POSITIONAL ARGUMENTS
        PATH
            Type: str | pathlib.Path

    FLAGS
        -o, --output_path=OUTPUT_PATH
            Type: Optional[str | pathlib.Path]
            Default: None
        -c, --channels=CHANNELS
            Type: Optional[list[int] | None]
            Default: None
    ```

1. Subtract

    ```text
    NAME
        rcashlar subtract

    SYNOPSIS
        rcashlar subtract BG_PATH AB_PATH <flags>

    POSITIONAL ARGUMENTS
        BG_PATH
            Type: str | pathlib.Path
        AB_PATH
            Type: str | pathlib.Path

    FLAGS
        -o, --output_path=OUTPUT_PATH
            Type: Optional[str | pathlib.Path]
            Default: None
        -f, --fiducial_channel=FIDUCIAL_CHANNEL
            Type: int
            Default: 0
        -b, --bg_intensity_scaling_factor=BG_INTENSITY_SCALING_FACTOR
            Type: Union
            Default: 'rcjob'
        -c, --camera_bias=CAMERA_BIAS
            Type: float
            Default: 105.0
        --add_camera_bias_back=ADD_CAMERA_BIAS_BACK
            Type: bool
            Default: False
        --as_float=AS_FLOAT
            Type: bool
            Default: False
    ```
