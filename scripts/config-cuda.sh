#!/usr/bin/env bash



../forestclaw/configure \
    --enable-cudaclaw \
    --without-blas \
    --enable-debug \
    --with-p4est=/home/donnacalhoun/projects/ForestClaw/code/p4est-build-serial/local \
    --with-sc=/home/donnacalhoun/projects/ForestClaw/code/p4est-build-serial/local \
    CPP="gcc -E" \
    CC=gcc \
    CFLAGS="-O0 -g -std=c99" \
    CXX=g++ \
    CXXFLAGS="-O0 -g" \
    FC=gfortran \
    F77=gfortran \
    FFLAGS="-O0 -g" \
    NVCC=nvcc \
    CUDA_CFLAGS=\
    CUDA_LDFLAGS=\
    LIBS="-L/cm/shared/apps/cuda90/toolkit/lib64/ -lcuda -lcudart" \
    --disable-shared

#   --disable-shared
#     FLIBS="-L/cm/shared/apps/gcc/4.8.1/lib64 -lgfortran" \
    