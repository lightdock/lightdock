import setuptools
from distutils.core import setup, Extension
from setuptools.command.test import test as TestCommand
import numpy


# Inspired by the example at https://pytest.org/latest/goodpractises.html
class NoseTestCommand(TestCommand):
    def finalize_options(self):
        TestCommand.finalize_options(self)
        self.test_args = []
        self.test_suite = True

    def run_tests(self):
        import nose
        nose.run_exit(argv=['nosetests'])


with open("README.md", "r") as fh:
    long_description = fh.read()

exts = [Extension(name='lightdock.mathutil.cython.cutil',
                  sources=["lightdock/mathutil/cython/cutil.c"],
                  include_dirs=[numpy.get_include()]),
        Extension(name='lightdock.mathutil.cython.quaternion',
                  sources=["lightdock/mathutil/cython/quaternion.c"],
                  include_dirs=[numpy.get_include()]),
        Extension(name='lightdock.gso.searchspace.cython.j1',
                  sources=["lightdock/gso/searchspace/cython/j1.c"],
                  include_dirs=[numpy.get_include()]),
        Extension(name='lightdock.gso.searchspace.cython.j2',
                  sources=["lightdock/gso/searchspace/cython/j2.c"],
                  include_dirs=[numpy.get_include()]),
        Extension(name='lightdock.gso.searchspace.cython.j3',
                  sources=["lightdock/gso/searchspace/cython/j3.c"],
                  include_dirs=[numpy.get_include()]),
        Extension(name='lightdock.gso.searchspace.cython.j4',
                  sources=["lightdock/gso/searchspace/cython/j4.c"],
                  include_dirs=[numpy.get_include()]),
        Extension(name='lightdock.gso.searchspace.cython.j5',
                  sources=["lightdock/gso/searchspace/cython/j5.c"],
                  include_dirs=[numpy.get_include()]),
        Extension(name='lightdock.scoring.pisa.cython.cpisa',
                  sources=["lightdock/scoring/pisa/cython/cpisa.c"],
                  include_dirs=[numpy.get_include()]),
        Extension(name='lightdock.scoring.dfire.cython.cdfire',
                  sources=["lightdock/scoring/dfire/cython/cdfire.c"],
                  include_dirs=[numpy.get_include()]),
        Extension(name='lightdock.scoring.dfire2.c.cdfire2',
                  sources=["lightdock/scoring/dfire2/c/cdfire2.c"],
                  include_dirs=[numpy.get_include()]),
        Extension(name='lightdock.scoring.sd.energy.c.sd',
                  sources=["lightdock/scoring/sd/energy/c/sd.c"],
                  include_dirs=[numpy.get_include()]),
        Extension(name='lightdock.scoring.fastdfire.c.cdfire',
                  sources=["lightdock/scoring/fastdfire/c/cdfire.c"],
                  include_dirs=[numpy.get_include()]),
        Extension(name='lightdock.scoring.cpydock.energy.c.cpydock',
                  sources=["lightdock/scoring/cpydock/energy/c/cpydock.c"],
                  include_dirs=[numpy.get_include()]),
        Extension(name='lightdock.scoring.vdw.energy.c.cvdw',
                  sources=["lightdock/scoring/vdw/energy/c/cvdw.c"],
                  include_dirs=[numpy.get_include()]),
        Extension(name='lightdock.scoring.dna.energy.c.cdna',
                  sources=["lightdock/scoring/dna/energy/c/cdna.c"],
                  include_dirs=[numpy.get_include()]),
        Extension(name='lightdock.scoring.sipper.c.sipper',
                  sources=["lightdock/scoring/sipper/c/sipper.c"],
                  include_dirs=[numpy.get_include()]),
]

setuptools.setup(
    name='lightdock',
    version='0.8.0b2',
    description="A macromolecular docking framework",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://lightdock.org/",
    packages=setuptools.find_packages(),
    include_package_data=True,
    license='GPLv3 License',
    classifiers=[
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Intended Audience :: Education",
        "Intended Audience :: Science/Research",
        "License :: Free For Educational Use",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: POSIX",
        "Programming Language :: Python :: Implementation :: CPython",
        "Topic :: Scientific/Engineering :: Artificial Intelligence",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Chemistry"
    ],
    python_requires='>=3.5',
    install_requires=[
        'numpy>=1.17.1', 
        'scipy>=1.3.1', 
        'cython>=0.29.13', 
        'biopython>=1.74',
        'pyparsing>=2.4.2',
        'prody>=1.10.11',
        'freesasa>=2.0.3',
    ],
    scripts=[
        'bin/ant_thony.py','bin/lgd_add_chain.py','bin/lgd_calculate_diameter.py',
        'bin/lgd_calculate_reference_points.py','bin/lgd_calculate_scoring.py',
        'bin/lgd_calculate_surface_density.py','bin/lgd_cluster_bsas.py',
        'bin/lgd_filter_membrane.py','bin/lgd_filter_restraints.py',
        'bin/lgd_generate_conformations.py','bin/lgd_generate_glowworm_positions.py',
        'bin/lgd_generate_trajectory.py','bin/lgd_gso_to_csv.py','bin/lgd_move_anm.py',
        'bin/lgd_prepare_new_simulation.py','bin/lgd_quaternion_to_euler.py',
        'bin/lgd_rank.py','bin/lgd_stats.py','bin/lgd_success_rate.py','bin/lgd_top.py',
        'bin/lightdock3.py','bin/lightdock3_setup.py'
    ],
    setup_requires=[
        'nose'
    ],
    cmdclass={'test': NoseTestCommand},
    ext_modules=exts,
    zip_safe=False
)
