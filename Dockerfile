FROM ubuntu

RUN apt-get update
RUN apt-get install -y libfftw3-dev libfftw3-single3 libfftw3-long3 \
  libfftw3-double3
RUN apt-get install -y python-pip python-dev build-essential \
  openjdk-8-jdk-headless python-tk
RUN pip install -q -U pip
RUN pip install -q numpy
RUN pip install -q scipy
RUN pip install ashlar

VOLUME /input
VOLUME /output
VOLUME /ffc

WORKDIR /output
ENTRYPOINT ["mosaic", "/input"]
