FROM hub.opensciencegrid.org/htc/ubuntu:24.04

RUN apt-get update && apt-get install -y \
    build-essential \
    cmake


COPY src ../../RCC_cpp/oneStep
WORKDIR /src

RUN cmake -B build
RUN cmake --build build --config Release

