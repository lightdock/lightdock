import numpy as np
from setuptools import setup, Extension
from Cython.Build import cythonize

setup(
    ext_modules=cythonize([
        Extension(
            "cdfire2",
            ["cdfire2.pyx"],
            include_dirs=[np.get_include()],
            define_macros=[("CYTHON_TRACE", "1")],
        )
    ])
)
