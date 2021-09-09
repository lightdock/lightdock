import numpy as np
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

setup(
    cmdclass={"build_ext": build_ext},
    ext_modules=[
        Extension(
            "cdfire2",
            ["cdfire2.pyx"],
            include_dirs=[np.get_include()],
            define_macros=[("CYTHON_TRACE", "1")],
        )
    ],
)
