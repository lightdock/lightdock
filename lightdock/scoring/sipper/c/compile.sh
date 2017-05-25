#!/bin/bash

rm -rf *.so
OPT="-g -O3" python setup.py build_ext --inplace
python setup.py clean
