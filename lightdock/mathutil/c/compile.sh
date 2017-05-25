#!/bin/bash

rm -rf *.so
# Compilation flags
OPT="-O3 -g" python setup.py build_ext --inplace

# Other possible compilation flags:
#OPT="-g -pg" CFLAGS="-O0 -g -pg" LDFLAGS="-g -pg" python setup.py build_ext --inplace
#OPT="-g -m64 -Ofast -flto -march=corei7-avx -funroll-loops -ftree-vectorize" python setup.py build_ext --inplace

python setup.py clean
