#!/bin/bash

rm -rf *.so *.c
python setup.py build_ext --inplace
python setup.py clean
