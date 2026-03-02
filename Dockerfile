FROM mambaorg/micromamba:1.5.1

COPY --chown=$MAMBA_USER:$MAMBA_USER docker-env.lock /tmp/docker-env.lock
RUN micromamba install --name base --yes --file /tmp/docker-env.lock \
    && micromamba clean --trash -aflp --yes

# add conda path to PATH to allow entrypoint overwrite
ENV PATH="${PATH}:/opt/conda/bin"

# pip install packages that are not available/problematic on conda-forge
RUN python -m pip install \
    --no-deps \
    opencv-python-headless==4.8.0.76 \
    palom==2023.9.2 \
    && python -m pip cache purge

RUN python -m pip install \
    --no-deps \
    "ashlar @ git+https://github.com/yu-anchen/ashlar@afr-2023-10-2" \
    && python -m pip cache purge
