#!/bin/bash

rm -rf assoc-mac
mkdir assoc-mac

LIB=/usr/local/lib
LIBGFORTRAN=${LIB}/libgfortran.a
LIBQUADMATH=${LIB}/libquadmath.a
LIBGCC=${LIB}/gcc/x86_64-apple-darwin17.5.0/8.1.0/libgcc.a

g++ *.cpp -o assoc-mac/assoc -s -O2 -std=c++11 -lopenblas $LIBGFORTRAN $LIBQUADMATH $LIBGCC
