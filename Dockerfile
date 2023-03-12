FROM ubuntu:18.04

MAINTAINER Dmitrii Smirnov "Dmitrii.Smirnov@skoltech.edu"

SHELL ["/bin/bash", "-c"]

RUN apt-get update && apt-get -y --no-install-recommends install \
    build-essential \
    clang \
    cmake \
    mpich \
    libboost-all-dev \
    wget








