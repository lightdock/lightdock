from distutils.core import setup, Extension
import numpy as np

setup(
    ext_modules=[Extension("cdfire", ["cdfire.c"])],
    include_dirs = [np.get_include()]
)
