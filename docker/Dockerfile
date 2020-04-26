FROM nvidia/cuda:10.0-devel-ubuntu18.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-key adv --keyserver keys.openpgp.org --recv-key 612DEFB798507F25 \
    && echo "deb [ arch=amd64 ] https://downloads.skewed.de/apt bionic main" | tee -a /etc/apt/sources.list \
    && apt-get update \
    && apt-get -y install sudo libcanberra-gtk3-module python3-graph-tool python3-cairo python3-matplotlib python3-gi python3-gi-cairo python3-pip \
    && pip3 install pycuda pyopenms Cython

RUN groupadd -g 1000 user \
    && useradd -m -g 1000 -u 1000 user \
    && passwd -d user \
    && usermod -a -G sudo user
