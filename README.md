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

Lightdock software is compatible and it has been tested with the followings OS:

* **macOS**: El Capitan, Sierra, High Sierra, Mojave, Catalina.
* **GNU/Linux**: Ubuntu 16+, Debian Stretch+, Scientific Linux 6+, CentOS 6+.

Microsoft Windows is not officially supported, despite many parts of the protocol might be able to run. Please use it at your own risk. If you wish to contribute testing and developing LightDock for Windows, please contact us.

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


### 3.2. Install LightDock
The fastest way to install LightDock is to use `pip`:

```bash
pip3 install lightdock
```

## 4. Development
For development and extension of the LightDock code, please follow these instructions:

### 4.1. Clone
Clone this repository:

```bash
git clone https://github.com/lightdock/lightdock.git
```

### 4.2. Compile Python C and Cython extensions

Please make sure dependencies are already installed (via pip, package manager, etc.):

* numpy>=1.17.1
* scipy>=1.3.1
* cython>=0.29.13
* biopython>=1.74
* pyparsing>=2.4.2
* prody>=1.10.11
* freesasa>=2.0.3

There is as bash script to compile all the extensions:

```bash
cd lightdock
./setup.sh
```

### 4.3. Add Lightdock to your path

Add the following lines to your `~/.bashrc` file, don't forget to change `/path/to/lightdock`:

```bash
# LightDock
export LIGHTDOCK_HOME="/path/to/lightdock"
export PATH=$PATH:$LIGHTDOCK_HOME/bin
export PYTHONPATH=$PYTHONPATH:$LIGHTDOCK_HOME
```

Don't forget to apply the changes:

```bash
source ~/.bashrc
```

### 4.4. Testing

You can run LightDock tests:

```bash
cd lightdock
nosetests-3.8 
```

## 5. Documentation

The complete documentation about how to run the LightDock protocol and several tutorials and use cases can be found at [https://lightdock.org/tutorials](https://lightdock.org/tutorials).


## 6. Get Help

LightDock is being actively developed and some issues may arise or you may need extra help to run LightDock. In those cases, there are two main ways to get help:

1. Read the [FAQ](https://lightdock.org/tutorials/faq) in case your problem is known
2. Open a [new issue in this repository](https://github.com/lightdock/lightdock/issues/new)
3. Or write an email to <b.jimenezgarcia@uu.nl>

## 7. LICENSE

LightDock is available under GPLv3 License. See [LICENSE](LICENSE) document for more details.
