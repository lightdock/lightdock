#!/bin/bash

rm -rf *.so
OPT="-g -O3" python3 setup.py build_ext --inplace
python3 setup.py clean
