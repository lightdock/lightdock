import setuptools
from Cython.Build import cythonize
from distutils.core import setup, Extension
import numpy

with open("README.md", "r") as fh:
    long_description = fh.read()

exts = [Extension(name='cutil',
                  sources=["lightdock/mathutil/cython/cutil.pyx"],
                  include_dirs=[numpy.get_include()]),
        Extension(name='quaternion',
                  sources=["lightdock/mathutil/cython/quaternion.pyx"],
                  include_dirs=[numpy.get_include()]),
        Extension(name='j1',
                  sources=["lightdock/gso/searchspace/cython/j1.pyx"],
                  include_dirs=[numpy.get_include()]),
        Extension(name='j2',
                  sources=["lightdock/gso/searchspace/cython/j2.pyx"],
                  include_dirs=[numpy.get_include()]),
        Extension(name='j3',
                  sources=["lightdock/gso/searchspace/cython/j3.pyx"],
                  include_dirs=[numpy.get_include()]),
        Extension(name='j4',
                  sources=["lightdock/gso/searchspace/cython/j4.pyx"],
                  include_dirs=[numpy.get_include()]),
        Extension(name='j5',
                  sources=["lightdock/gso/searchspace/cython/j5.pyx"],
                  include_dirs=[numpy.get_include()]),
        Extension(name='cpisa',
                  sources=["lightdock/scoring/pisa/cython/cpisa.pyx"],
                  include_dirs=[numpy.get_include()]),
        Extension(name='cdfire2',
                  sources=["lightdock/scoring/dfire2/c/cdfire2.c"],
                  include_dirs=[numpy.get_include()]),
        Extension(name='sd',
                  sources=["lightdock/scoring/sd/energy/c/sd.c"],
                  include_dirs=[numpy.get_include()]),
        Extension(name='fastdfire',
                  sources=["lightdock/scoring/fastdfire/c/cdfire.c"],
                  include_dirs=[numpy.get_include()]),
        Extension(name='cpydock',
                  sources=["lightdock/scoring/cpydock/energy/c/cpydock.c"],
                  include_dirs=[numpy.get_include()]),
        Extension(name='cvdw',
                  sources=["lightdock/scoring/vdw/energy/c/cvdw.c"],
                  include_dirs=[numpy.get_include()]),
        Extension(name='cdna',
                  sources=["lightdock/scoring/dna/energy/c/cdna.c"],
                  include_dirs=[numpy.get_include()]),
        Extension(name='sipper',
                  sources=["lightdock/scoring/sipper/c/sipper.c"],
                  include_dirs=[numpy.get_include()]),
]

setuptools.setup(
    name='lightdock',
    version='0.0.3',
    description="A macromolecular docking framework",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://lightdock.org/",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires=['numpy', 'scipy', 'cython', 'biopython', 'freesasa', 'prody'],
    scripts=['bin/ant_thony.py','bin/lgd_add_chain.py','bin/lgd_calculate_diameter.py','bin/lgd_calculate_reference_points.py','bin/lgd_calculate_scoring.py','bin/lgd_calculate_surface_density.py','bin/lgd_cluster_bsas.py','bin/lgd_filter_membrane.py','bin/lgd_filter_restraints.py','bin/lgd_generate_conformations.py','bin/lgd_generate_glowworm_positions.py','bin/lgd_generate_trajectory.py','bin/lgd_gso_to_csv.py','bin/lgd_move_anm.py','bin/lgd_prepare_new_simulation.py','bin/lgd_quaternion_to_euler.py','bin/lgd_rank.py','bin/lgd_stats.py','bin/lgd_success_rate.py','bin/lgd_top.py','bin/lightdock3.py','bin/lightdock3_setup.py'],
    ext_modules=cythonize(exts)
)
