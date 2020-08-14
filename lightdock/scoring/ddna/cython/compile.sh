#!/bin/bash

rm -rf *.so *.c
python3 setup.py build_ext --inplace
python3 setup.py clean
