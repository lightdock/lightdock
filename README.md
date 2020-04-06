[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![PyPi version](https://img.shields.io/pypi/v/lightdock.svg)](https://pypi.org/project/lightdock/)
[![PyPi Downloads](https://img.shields.io/pypi/dm/lightdock?label=PyPI%20Downloads)](https://pypistats.org/packages/lightdock)
[![Supported versions](https://img.shields.io/pypi/pyversions/lightdock.svg)](https://pypi.org/project/lightdock)
[![Build Status](https://travis-ci.com/lightdock/lightdock.svg?branch=master)](https://travis-ci.com/lightdock/lightdock)
[![Code Coverage](https://codecov.io/gh/lightdock/lightdock/branch/master/graph/badge.svg)](https://codecov.io/gh/lightdock/lightdock)
[![Commit activity](https://img.shields.io/github/commit-activity/m/lightdock/lightdock.svg)](https://github.com/lightdock/lightdock/commits/master)
 
<p align="center">
    <img src="https://lightdock.org/assets/images/lightdock_logo.png">
</p>

## 1. Synopsis
**LightDock** is a protein-protein, protein-peptide and protein-DNA docking framework based on the [Glowworm Swarm Optimization](https://link.springer.com/article/10.1007/s11721-008-0021-5) (GSO) algorithm.

The LightDock framework is highly versatile, with many options that can be further developed and optimized by the users: it can accept any user-defined scoring function, can use local gradient-free minimization, the simulation can be restrained from the beginning to focus on user-assigned interacting regions, **it supports residue restraints in both receptor and ligand partners**.

## 2. Reference
LightDock protocol and the updates to make use of residue restraints have been published in [Oxford Bioinformatics](https://academic.oup.com/bioinformatics) journal. Please cite these references if you use LightDock in your research:

**LightDock: a new multi-scale approach to protein–protein docking**<br>
[Brian Jiménez-García](http://bjimenezgarcia.com), Jorge Roel-Touris, Miguel Romero-Durana, Miquel Vidal, Daniel Jiménez-González and Juan Fernández-Recio<br>
*Bioinformatics*, Volume 34, Issue 1, 1 January 2018, Pages 49–55, [https://doi.org/10.1093/bioinformatics/btx555](https://doi.org/10.1093/bioinformatics/btx555)

**LightDock goes information-driven**<br>
Jorge Roel-Touris, Alexandre M.J.J. Bonvin, [Brian Jiménez-García](http://bjimenezgarcia.com)<br>
*Bioinformatics*, btz642; doi: [https://doi.org/10.1093/bioinformatics/btz642](https://doi.org/10.1093/bioinformatics/btz642)


## 3. Installation
### 3.1. Dependencies
LightDock has the following dependencies:

* NumPy (<http://www.numpy.org/>)
* Scipy (<http://www.scipy.org/>)
* Cython (<http://cython.org/>)
* BioPython (<http://biopython.org>)
* ProDy (<http://prody.csb.pitt.edu/>)
* Freesasa (<http://freesasa.github.io/>)

Optional dependencies are:

* MPI4py (if you plan to use MPI support which is experimental at the moment, <http://pythonhosted.org/mpi4py/>)


#### 3.1.1. Installing NumPy, Scipy, Cython, Biopython and ProDy

```bash
pip3 install numpy, scipy, cython, biopython, pyparsing, prody
```

| WARNING: Issues installing ProDy in macOS have been reported. In case you are not able to install ProDy using `pip`, you may try to install it from the latest version in GitHub. Please see instructions [here](https://github.com/prody/ProDy/issues/864). |
| --- |


### 3.2. Install LightDock
The fastest way to install LightDock is to use `pip`:

```bash
pip3 install lightdock
```

### 3.3. Alternative installations
Please visit the old Python-2.7 [repository documentation on GitHub](https://github.com/lightdock/lightdock-python2.7/tree/master) for alternative installation.


## 4. Documentation

The complete documentation about how to run the LightDock protocol can be found at [https://lightdock.org](https://lightdock.org).


## 5. Get Help

LightDock is being actively developed and some issues may arise or you may get some extra help to run LightDock. In those cases, there are two main ways to get help:

1. Read the [FAQ](https://lightdock.org/tutorials/faq) in case your problem is known
2. Open a [new issue in this repository](https://github.com/lightdock/lightdock/issues/new)
3. Or write an email to <b.jimenezgarcia@uu.nl>

## 6. LICENSE

LightDock is available under GPLv3 License. See LICENSE document for more details.
