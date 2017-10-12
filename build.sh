#!/bin/bash

rm -rf assoc-$1
mkdir assoc-$1

if [ $1 == "lnx64" ]; then

    g++ *.cpp -o assoc-$1/assoc -s -O2 -std=c++11 -static -lopenblasp -lgfortran -lquadmath -lpthread

elif [ $1 == "win32" ]; then

    i686-w64-mingw32-g++ *.cpp -o assoc-$1/assoc.exe -s -O2 -std=c++11 -static -L./mingw32/lib -lopenblas -lgfortran -lquadmath

elif [ $1 == "win64" ]; then

    x86_64-w64-mingw32-g++ *.cpp -o assoc-$1/assoc.exe -s -O2 -std=c++11 -static -L./mingw64/lib -lopenblas -lgfortran -lquadmath

fi
