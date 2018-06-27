FROM ubuntu

RUN apt-get update \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y \
        build-essential \
        ffmpeg \
        git \
        libfftw3-dev \
        libfftw3-double3 \
        libfftw3-long3 \
        libfftw3-single3 \
        openjdk-8-jdk-headless \
        python-dev \
        python-pip \
        python-tk \
    && rm -rf /var/lib/apt/lists/*

RUN pip install -q -U \
    cython \
    numpy \
    pip \
    scipy

COPY / /app/ashlar/
RUN pip install /app/ashlar

ENV JAVA_HOME /usr/lib/jvm/java-8-openjdk-amd64
ENV OMP_NUM_THREADS 1

VOLUME /data
VOLUME /data2
VOLUME /data3

WORKDIR /data
