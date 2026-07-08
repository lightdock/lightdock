from setuptools import setup, Extension
from Cython.Build import cythonize

extensions = [
    Extension(f"j{i}", [f"j{i}.pyx"]) 
    for i in range(1, 6)
]

setup(
    ext_modules=cythonize(extensions),
)
