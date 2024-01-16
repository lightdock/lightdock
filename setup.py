import os
import setuptools
from distutils.core import Extension


# MDAnalysis NumPy delay on setup.py
def abspath(file):
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), file)


class LDExtension(Extension, object):
    """Derived class for handling setup-time (numpy) dependencies."""

    def __init__(self, name, sources, *args, **kwargs):
        self._ld_include_dirs = []
        super(LDExtension, self).__init__(name, sources, *args, **kwargs)

    @property
    def include_dirs(self):
        if not self._ld_include_dirs:
            for item in self._ld_include_dir_args:
                try:
                    self._ld_include_dirs.append(item())
                except TypeError:
                    item = abspath(item)
                    self._ld_include_dirs.append((item))
        return self._ld_include_dirs

    @include_dirs.setter
    def include_dirs(self, val):
        self._ld_include_dir_args = val


def get_numpy_include():
    import builtins

    builtins.__NUMPY_SETUP__ = False
    try:
        import numpy as np
    except ImportError:
        raise SystemExit("LightDock requires NumPy for setup")
    return np.get_include()


exts = [
    LDExtension(
        name="lightdock.mathutil.cython.cutil",
        sources=["lightdock/mathutil/cython/cutil.c"],
        include_dirs=[get_numpy_include],
    ),
    LDExtension(
        name="lightdock.mathutil.cython.quaternion",
        sources=["lightdock/mathutil/cython/quaternion.c"],
        include_dirs=[get_numpy_include],
    ),
    LDExtension(
        name="lightdock.gso.searchspace.cython.j1",
        sources=["lightdock/gso/searchspace/cython/j1.c"],
        include_dirs=[get_numpy_include],
    ),
    LDExtension(
        name="lightdock.gso.searchspace.cython.j2",
        sources=["lightdock/gso/searchspace/cython/j2.c"],
        include_dirs=[get_numpy_include],
    ),
    LDExtension(
        name="lightdock.gso.searchspace.cython.j3",
        sources=["lightdock/gso/searchspace/cython/j3.c"],
        include_dirs=[get_numpy_include],
    ),
    LDExtension(
        name="lightdock.gso.searchspace.cython.j4",
        sources=["lightdock/gso/searchspace/cython/j4.c"],
        include_dirs=[get_numpy_include],
    ),
    LDExtension(
        name="lightdock.gso.searchspace.cython.j5",
        sources=["lightdock/gso/searchspace/cython/j5.c"],
        include_dirs=[get_numpy_include],
    ),
    LDExtension(
        name="lightdock.scoring.pisa.cython.cpisa",
        sources=["lightdock/scoring/pisa/cython/cpisa.c"],
        include_dirs=[get_numpy_include],
    ),
    LDExtension(
        name="lightdock.scoring.dfire.cython.cdfire",
        sources=["lightdock/scoring/dfire/cython/cdfire.c"],
        include_dirs=[get_numpy_include],
    ),
    LDExtension(
        name="lightdock.scoring.ddna.cython.cddna",
        sources=["lightdock/scoring/ddna/cython/cddna.c"],
        include_dirs=[get_numpy_include],
    ),
    LDExtension(
        name="lightdock.scoring.dfire2.c.cdfire2",
        sources=["lightdock/scoring/dfire2/c/cdfire2.c"],
        include_dirs=[get_numpy_include],
    ),
    LDExtension(
        name="lightdock.scoring.sd.energy.c.sd",
        sources=["lightdock/scoring/sd/energy/c/sd.c"],
        include_dirs=[get_numpy_include],
    ),
    LDExtension(
        name="lightdock.scoring.fastdfire.c.cdfire",
        sources=["lightdock/scoring/fastdfire/c/cdfire.c"],
        include_dirs=[get_numpy_include],
    ),
    LDExtension(
        name="lightdock.scoring.cpydock.energy.c.cpydock",
        sources=["lightdock/scoring/cpydock/energy/c/cpydock.c"],
        include_dirs=[get_numpy_include],
    ),
    LDExtension(
        name="lightdock.scoring.vdw.energy.c.cvdw",
        sources=["lightdock/scoring/vdw/energy/c/cvdw.c"],
        include_dirs=[get_numpy_include],
    ),
    LDExtension(
        name="lightdock.scoring.dna.energy.c.cdna",
        sources=["lightdock/scoring/dna/energy/c/cdna.c"],
        include_dirs=[get_numpy_include],
    ),
    LDExtension(
        name="lightdock.scoring.sipper.c.sipper",
        sources=["lightdock/scoring/sipper/c/sipper.c"],
        include_dirs=[get_numpy_include],
    ),
]

setuptools.setup(
    packages=setuptools.find_namespace_packages(),
    include_package_data=True,
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
    ext_modules=exts,
    zip_safe=False,
)
