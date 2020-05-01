FROM nvidia/cuda:10.0-devel-ubuntu18.04


ENV DEBIAN_FRONTEND=noninteractive


# Install graph-tool, ignore recommends as it installs packages for Python 2, not Python 3
RUN apt-key adv --keyserver keys.openpgp.org --recv-key 612DEFB798507F25 \
    && echo "deb [ arch=amd64 ] https://downloads.skewed.de/apt bionic main" | tee -a /etc/apt/sources.list \
    && apt-get update \
    && apt-get -y install python3-gv python3-tk python3-cairo python3-gi-cairo python3-gi gir1.2-gtk-3.0 \
    && apt-get -y install --no-install-recommends python3-graph-tool


# Install pip and dependencies, upgrade scipy, numpy and matplotlib
RUN apt-get -y install curl python3 python3-dev python3-distutils \
    && curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py \
    && python3 get-pip.py \
    && pip install pycuda pyopenms Cython scipy numpy matplotlib --upgrade


# Install sudo and add non-root user
RUN apt-get -y install sudo \
    && groupadd -g 1000 user \
    && useradd -m -g 1000 -u 1000 user \
    && passwd -d user \
    && usermod -a -G sudo user
