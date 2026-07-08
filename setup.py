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
    packages=setuptools.find_namespace_packages(),
    include_package_data=True,
    ext_modules=exts,
    zip_safe=False,
    scripts=[
        "bin/ant_thony.py",
        "bin/lgd_calculate_diameter.py",
        "bin/lgd_calculate_reference_points.py",
        "bin/lgd_calculate_scoring.py",
        "bin/lgd_cluster_bsas.py",
        "bin/lgd_copy_structures.py",
        "bin/lgd_create_membrane.py",
        "bin/lgd_dummify.py",
        "bin/lgd_filter_membrane.py",
        "bin/lgd_filter_restraints.py",
        "bin/lgd_flatten.py",
        "bin/lgd_generate_conformations.py",
        "bin/lgd_generate_glowworm_positions.py",
        "bin/lgd_generate_trajectory.py",
        "bin/lgd_gso_to_csv.py",
        "bin/lgd_map_contacts.py",
        "bin/lgd_move_anm.py",
        "bin/lgd_rank.py",
        "bin/lgd_rank_swarm.py",
        "bin/lgd_run.py",
        "bin/lgd_setup.py",
        "bin/lgd_top.py",
    ],
)