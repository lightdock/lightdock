![LightDock](media/lightdock_banner.png "LightDock")

<hr>

# Official Documentation

#### Table of Contents

- [Introduction](#1-introduction)
- [Setup a simulation](#2-setup-a-simulation)
- [Run a simulation](#3-run-a-simulation)
- [Generate models](#4-generate-models)
- [Clustering](#5-clustering)
- [Custom Scoring Functions](#6-custom-scoring-functions)


## 1. Introduction

### 1.1. What is LightDock?
LightDock is a protein-protein and protein-DNA docking protocol based on the [Glowworm Swarm Optimization](https://link.springer.com/article/10.1007/s11721-008-0021-5) (GSO) algorithm. The LightDock protocol has been published in Oxford Bioinformatics:

**LightDock: a new multi-scale approach to protein–protein docking**<br>
[Brian Jiménez-García](https://scholar.google.es/citations?user=eVN1WVYAAAAJ&hl=en), Jorge Roel-Touris, Miguel Romero-Durana, Miquel Vidal, Daniel Jiménez-González and Juan Fernández-Recio<br>
*Bioinformatics*, Volume 34, Issue 1, 1 January 2018, Pages 49–55, [https://doi.org/10.1093/bioinformatics/btx555](https://doi.org/10.1093/bioinformatics/btx555)

Specific details of the protocol can be found in the publication mentioned above, but in summary, LightDock is:

- *Ab initio* protocol, which means that only requires of the 3D coordinates of the protein partners for predicting the protein-protein or protein-DNA complex.

- Capable of modeling protein-protein and protein-DNA complexes in rigid-body fashion or modeling backbone flexibility using [Anisotropic Network Model](https://en.wikipedia.org/wiki/Anisotropic_Network_Model) (ANM). If ANM mode is activated, LightDock calculates the Ca-Ca ANM model using the awesome [ProDy](http://prody.csb.pitt.edu/) Python library. By default, the first 10 non-trivial normal modes are calculated for both receptor and ligand (in every residue backbone, extended to side-chains). See [Prody ANM documentation](http://prody.csb.pitt.edu/tutorials/enm_analysis/anm.html) for an example.

- Customizable by the user. LightDock is not only a protocol, but a framework for testing and developing custom scoring functions. The GSO optimization algorithm is agnostic of the force-field used, so in theory LightDock is capable of minimizing the docking energies in any force-field given by the user. See *Custom Scoring Functions* section for more details.

- Prepared to scale in [HPC](https://en.wikipedia.org/wiki/Supercomputer) architectures. LightDock nature is *embarrassingly parallel* as each swarm is treated as an independent simulation. This property makes LightDock to scale up to a number of CPU cores equal to the number of swarms simulated. Two implementations are given: 1) [multiprocessing](https://docs.python.org/2/library/multiprocessing.html) (by default) and 2) MPI (using [mpi4py](http://mpi4py.scipy.org/docs/) library).

- Capable of using multiple scoring functions during the minimization. Instead of specifiying a single scoring function, a file containing the weight and the name of the scoring function can be given as an input. LightDock will combine the different scoring functions as a linear combination of their value multiplied by the weight specified in the file.


### 1.2. Swarms
In LightDock, the receptor molecule is keep fixed (despite atoms positions can move if ANM mode is enabled). Over its surface, a set of points is calculated. Each of these points is a swarm center which represent an independent simulation. For example, for complex [1VFB](https://www.rcsb.org/structure/1VFB), 400 swarms are calculated:

![Swarms centers](media/1vfb_cswarms_centers.png "Swarm centers")

For each of these swarm centers, a number *N* of glowworms, the algorithm agents, are disposed in a random way. Every glowworm represents a possible ligand conformation. In the following figure a set of 300 glowworms are displayed in a single swarm:

![Swarms and glowworms](media/swarm.png "Swarm and glowworms")

More in detail, each glowworm is represented as a 3D-axis object in its center of mass and oriented as the actual 3D-axis orientation:

![Swarms and glowworms](media/1e6e_swarm.png "Swarm and glowworms")


## 2. Setup a simulation

In LightDock versions prior to 0.5.0, this step was optional. From version 0.5.0, a setup step is mandatory for any simulation. This new Lightdock setup reduces complexity of the main program and it will allow the user to apply distance restraints.

If you execute `lightdock_setup` several options will appear:

```bash
usage: lightdock_setup [-h] [--seed_points STARTING_POINTS_SEED]
                       [-ft ftdock_file] [--noxt] [-anm] [--seed_anm ANM_SEED]
                       [-anm_rec ANM_REC] [-anm_lig ANM_LIG]
                       receptor_pdb_file ligand_pdb_file swarms glowworms
lightdock_setup: error: too few arguments

```

The mandatory arguments are the PDB file containing the receptor structure, the PDB file containing the ligand structure, the number of swarms of the simulations and the number or glowworms of each swarm in this order.

In the original publication, for each complex of the [Protein-Protein Benchmark v5](https://zlab.umassmed.edu/benchmark/) analyzed, the number of swarms and glowworms were:

```
Number of swarms: 400
Number of glowworms per swarm: 300
```

We found reasonable those values for number of swarms and glowworms per swarm in order to ensure a certain density and exhaustiveness of sampling, despite for smaller molecules, a smaller number of swarms could be completely valid.

Below there is a description of the rest of accepted paramenters of `lightdock_setup`:

- **seed_points** *STARTING_POINTS_SEED*: An integer can be specified as the seed used in the random number generator of the initial random poses of the ligand.
- **ft** *ftdock_file*: LightDock can read the output of the venerable [FTDock](http://www.sbg.bio.ic.ac.uk/docking/ftdock.html) software in order to use the FTDock rigid-body predictions as the starting poses of the LightDock simulation. In order to do so, `lightdock_setup` classifies the different FTDock predictions according to its translation into the corresponding swarm over the surface of the receptor.
- **noxt**: If this option is enabled, LightDock ignores OXT atoms. This is useful for several scoring functions which don't understand this special type of atom.
- **anm**: If this option is enabled, the ANM mode is activated and backbone flexibility is modeled using ANM (via ProDy).
- **seed_anm** *ANM_SEED*: An integer can be specified as the seed used in the random number generator of ANM normal modes extents. 
- **anm_rec** *ANM_REC*: The number of non-trivial normal modes calculated for the recepetor in the ANM mode.
- **anm_lig** *ANM_LIG*: The number of non-trivial normal modes calculated for the ligand in the ANM mode.

### 2.1. Results of the setup

After the execution of `lightdock_setup` script, several files and directories will be generated in the root of your docking project:

- `init`: A directory containing initialization data, see below:
    - `cluster_centers.pdb`: A file in PDB format containing dummy atoms which coordinates correspond to a swarm center.
    - `initial_positions_X.dat`: One file like this for each swarm (where X indicates the id of the swarm), containing the translation, quaternion and normal modes extents (if ANM is activated) for each of the glowworms in this swarm.
    - `starting_positions_X.pdb`: Same as before, but a PDB file with the glowworm translation vector expressed as a dummy atom.
    - `starting_poses_X.bild`: Same as before, but this time every glowworm pose in rotational space is expressed as a BILD 3D object which can be readed by UCSF Chimera. More information about the format can be found [here](https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/bild.html).
- `swarm_0, ..., swarm_(N-1)`: For each of the N swarms specified by the `swarms` argument, a directory is created. Note that the first directory is called `swarm_0` and not `swarm_1`.
- `setup.json`: A file with a summary of the setup step in JSON format. This file will be necessary for running the simulation in the next step.
- `lightdock_*.pdb`: LightDock parses the input PDB structures and creates two new PDB files, one for the receptor and one for the ligand.
- `*.xyz.npy`: Two files with this extension, one for the receptor and one for the ligand, which contains information about the minimum ellipsoid containing each of the structures in NumPy format.
- `lightdock_rec.nm.npy`and `lightdock_lig.nm.npy`: If ANM is enabled, two files are created containing the normal modes information calculated by the ProDy library.

### 2.2. Tips and tricks

- As a general rule of thumb, the receptor structure is the bigger molecule and the ligand the smaller. With bigger and smaller, the metric used is the longest diameter of the minimum ellipsoid containing the molecule. A script called `lgd_calculate_diameter.py` can be found in `$LIGHTDOCK_HOME/bin/support` path in order to calculate an approximation of that metric.

- If the `init` directory exists, `lightdock_setup` makes use of the information contained in that directory.

- The file `setup.json` is intended for reproducibility of the results, but not to be modified by the user.

## 3. Run a simulation

TBC


## 4. Generate models

TBC


## 5. Clustering

TBC


## 6. Custom Scoring Functions

TBC