# LightDock

## 1. Introduction
LightDock is a protein-protein, protein-peptide and protein-DNA docking framework based on the [Glowworm Swarm Optimization](https://link.springer.com/article/10.1007/s11721-008-0021-5) (GSO) algorithm.

The LightDock framework is highly versatile, with many options that can be further developed and optimized by the users: it can accept any user-defined scoring function, can use local gradient-free minimization, the simulation can be restrained from the beginning to focus on user-assigned interacting regions, **it supports residue restraints in both receptor and ligand partners** and it has support for the use of pre-calculated conformers for both receptor and ligand.

## 2. Reference
The first version of the LightDock protocol was published in [Oxford Bioinformatics](https://academic.oup.com/bioinformatics) journal. Please cite this reference if you use LightDock in your research:

**LightDock: a new multi-scale approach to protein–protein docking**<br>
[Brian Jiménez-García](http://bjimenezgarcia.com), Jorge Roel-Touris, Miguel Romero-Durana, Miquel Vidal, Daniel Jiménez-González and Juan Fernández-Recio<br>
*Bioinformatics*, Volume 34, Issue 1, 1 January 2018, Pages 49–55, [https://doi.org/10.1093/bioinformatics/btx555](https://doi.org/10.1093/bioinformatics/btx555)

A second article about the implementation details and performance of the new protocol for including residue restraints is avaiable:

**LightDock goes information-driven**<br>
Jorge Roel-Touris, Alexandre M.J.J. Bonvin, Brian Jiménez-García<br>
*Bioinformatics*, , btz642; doi: [https://doi.org/10.1093/bioinformatics/btz642](https://doi.org/10.1093/bioinformatics/btz642)


## 3. Installation
### 3.1. Dependencies
LightDock has the following dependencies:

* NumPy (<http://www.numpy.org/>)
* Scipy (<http://www.scipy.org/>)
* Cython (<http://cython.org/>)
* BioPython (<http://biopython.org>)
* ProDy (<http://prody.csb.pitt.edu/>)

Optional dependencies:

* MPI4py (if you plan to use MPI support which is experimental at the moment, <http://pythonhosted.org/mpi4py/>)
* Freesasa (only if `cpydock` scoring function is used, <http://freesasa.github.io/>)

#### 3.1.1. Installing NumPy, Scipy, Cython and Biopython

```bash
pip3 install numpy, scipy, cython, biopython, pyparsing, prody
```

**Make sure all libraries are from the same Python 3.5+ series.**


### 3.2. Install LightDock
The fastest way to install LightDock is to use `pip`:

```bash
pip3 install lightdock
```

### 3.3. Alternative installations
Please visit the [repository documentation on GitHub](https://github.com/brianjimenez/lightdock/tree/python3) for alternative installation.


## 4. Documentation

The complete documentation about how to run the LightDock protocol can be found at:

<p align="center">
	<b><a href="https://brianjimenez.github.io/lightdock">https://brianjimenez.github.io/lightdock/</a>
	</b>
</p>


## 5. Get Help

LightDock is being actively developed and some issues may arise or you may get some extra help to run LightDock. In those cases, there are two main ways to get help:

1. Read the [FAQ](https://github.com/brianjimenez/lightdock/blob/python3/docs/FAQ.md) in case your problem is known
2. Open a [new issue in this repository](https://github.com/brianjimenez/lightdock/issues/new)
3. Or write an email to <b.jimenezgarcia@uu.nl>
