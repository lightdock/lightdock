from distutils.core import setup, Extension

setup(
    ext_modules=[Extension("quaternion", ["quaternion.c"])]
)
