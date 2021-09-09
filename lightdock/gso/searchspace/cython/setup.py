from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

setup(
    cmdclass={"build_ext": build_ext},
    ext_modules=[
        Extension("j1", ["j1.pyx"]),
        Extension("j2", ["j2.pyx"]),
        Extension("j3", ["j3.pyx"]),
        Extension("j4", ["j4.pyx"]),
        Extension("j5", ["j5.pyx"]),
    ],
)
