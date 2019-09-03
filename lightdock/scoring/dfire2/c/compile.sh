#!/bin/bash

rm -rf *.so
#OPT="-g3" CFLAGS="-O0 -g3" LDFLAGS="-g3" python3 setup.py build_ext --inplace
#OPT="-g -O0" python3 setup.py build_ext --inplace
OPT="-g -O3" python3 setup.py build_ext --inplace
#OPT="-g -m64 -Ofast -flto -march=corei7-avx -funroll-loops -ftree-vectorize" python3 setup.py build_ext --inplace
#python3 setup.py build_ext --inplace
python3 setup.py clean
