1. Generate conda env spec

```bash
docker run -it --rm --platform linux/amd64 -v "$(pwd)":/data mambaorg/micromamba:1.5.1 bash


micromamba create -y -n rcashlar -c conda-forge openjdk=11 python=3.10 "scikit-image<0.20" scikit-learn "zarr<2.15" tifffile imagecodecs matplotlib networkx seaborn pyjnius blessed tqdm scipy dask numpy loguru=0.5.3 "ome-types>0.3" "pydantic<2" pint napari-lazy-openslide yamale fire termcolor git wget unzip procps-ng

micromamba activate rcashlar

python -m pip install --dry-run "ashlar @ git+https://github.com/yu-anchen/ashlar@afr-2023-10-2"
# Would install ashlar-1.17.0+48.gdb3e2e3 opencv-python-4.8.1.78 palom-2023.9.2

micromamba env export --explicit > /data/docker-env.lock

# optional
python -m pip install --no-deps "ashlar @ git+https://github.com/yu-anchen/ashlar@afr-2023-10-2"
python -m pip install --no-deps palom==2023.9.2 opencv-python-headless==4.8.0.76
```

1. Test build

```bash
docker build --platform linux/amd64 --tag afr-2023-10-2 .
```

1. Test run

```bash
docker run -it --rm --platform linux/amd64 afr-2023-10-2 bash
```
