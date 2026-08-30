import setuptools
from setuptools import Extension
import numpy as np

numpy_include = np.get_include()

extension_names = [
    "lightdock.mathutil.cython.cutil",
    "lightdock.mathutil.cython.quaternion",
    "lightdock.gso.searchspace.cython.j1",
    "lightdock.gso.searchspace.cython.j2",
    "lightdock.gso.searchspace.cython.j3",
    "lightdock.gso.searchspace.cython.j4",
    "lightdock.gso.searchspace.cython.j5",
    "lightdock.scoring.pisa.cython.cpisa",
    "lightdock.scoring.dfire.cython.cdfire",
    "lightdock.scoring.ddna.cython.cddna",
    "lightdock.scoring.dfire2.c.cdfire2",
    "lightdock.scoring.sd.energy.c.sd",
    "lightdock.scoring.fastdfire.c.cdfire",
    "lightdock.scoring.cpydock.energy.c.cpydock",
    "lightdock.scoring.vdw.energy.c.cvdw",
    "lightdock.scoring.dna.energy.c.cdna",
    "lightdock.scoring.sipper.c.sipper",
]

exts = [
    Extension(
        name=name,
        sources=[f"{name.replace('.', '/')}.c"],
        include_dirs=[numpy_include]
    )
    for name in extension_names
]

setuptools.setup(
    include_package_data=True,
    ext_modules=exts,
    zip_safe=False,
)