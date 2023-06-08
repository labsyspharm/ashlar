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
python -m pip install "ashlar @ git+https://github.com/yu-anchen/ashlar@dd08fadaf19bf7300a6b5d162e0fb32c4212a2dc"
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

## Example scenario and commands

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
        rcashlar rcashlar register /project/demo/21-2D474_001-2D01@20230523_173947_224759 /project/demo/21-2D474_001-2D01@20230523_183243_352969
        ```

    - Run `rcashlar subtract` against the second scan

        ```bash
        rcashlar subtract /project/demo/230523_3scans/21-2D474_001-2D01@20230523_180743_729614 /project/demo/230523_3scans/21-2D474_001-2D01@20230523_183243_352969
        ```
