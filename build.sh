#!/bin/bash

# -lopenblasp / -lopenblaso
g++ *.cpp -o assoc -s -O2 -std=c++11 -static -lopenblas -lgfortran -lquadmath
