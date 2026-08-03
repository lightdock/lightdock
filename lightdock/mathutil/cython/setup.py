from setuptools import setup, Extension
from Cython.Build import cythonize

# Define the extensions
extensions = [
    Extension("cutil", ["cutil.pyx"]),
    Extension("quaternion", ["quaternion.pyx"]),
]

# Pass the extensions to cythonize() inside setup()
setup(
    ext_modules=cythonize(extensions),
)
