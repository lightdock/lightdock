![LightDock](docs/media/lightdock_banner.png "LightDock")

#### Table of Contents

- [1. Introduction](#1-introduction)
- [2. Reference](#2-reference)
- [3. Installation](#3-installation)
- [4. Documentation](#4-documentation)
- [5. Get Help](#5-get-help)

## 1. Introduction
LightDock is a protein-protein, protein-peptide and protein-DNA docking framework based on the [Glowworm Swarm Optimization](https://link.springer.com/article/10.1007/s11721-008-0021-5) (GSO) algorithm.

**The framework is written in the Python programming language (version 2.7) and allows the users to incorporate their own scoring function.**

The LightDock framework is highly versatile, with many options that can be further developed and optimized by the users: it can accept any user-defined scoring function, can use local gradient-free minimization, the simulation can be restrained from the beginning to focus on user-assigned interacting regions, **it supports residue restraints in both receptor and ligand partners** and it has support for the use of pre-calculated conformers for both receptor and ligand.

## 2. Reference
The first version of the LightDock protocol was published in [Oxford Bioinformatics](https://academic.oup.com/bioinformatics) journal. Please cite this reference if you use LightDock in your research:

**LightDock: a new multi-scale approach to protein–protein docking**<br>
[Brian Jiménez-García](http://bjimenezgarcia.com), Jorge Roel-Touris, Miguel Romero-Durana, Miquel Vidal, Daniel Jiménez-González and Juan Fernández-Recio<br>
*Bioinformatics*, Volume 34, Issue 1, 1 January 2018, Pages 49–55, [https://doi.org/10.1093/bioinformatics/btx555](https://doi.org/10.1093/bioinformatics/btx555)

A preprint about the implementation details and performance of the new protocol for including residue restraints is avaiable at [biorxiv](https://www.biorxiv.org/content/10.1101/595983v1):

**LightDock goes information-driven**<br>
Jorge Roel-Touris, Alexandre M.J.J. Bonvin, Brian Jiménez-García<br>
*bioRxiv* 595983; doi: [https://doi.org/10.1101/595983](https://doi.org/10.1101/595983)

## 3. Installation
### 3.1. Dependencies
LightDock has the following dependencies:

* **Python 2.7.x**
* Nose (<http://nose.readthedocs.io/en/latest/>)
* NumPy (<http://www.numpy.org/>)
* Scipy (<http://www.scipy.org/>)
* Cython (<http://cython.org/>)
* BioPython (<http://biopython.org>)
* MPI4py (<http://pythonhosted.org/mpi4py/>)
* ProDy (<http://prody.csb.pitt.edu/>)
* Freesasa (only if `cpydock` scoring function is used and to run the complete test set, <http://freesasa.github.io/>)

#### 3.1.1. Installing NumPy, Scipy, Cython, Biopython, Nose and MPI4py
NumPy, Scipy, Cython, Biopython, Nose and MPI4py libraries are usually available as packages in most of GNU/Linux distributions. For example, to install them in Ubuntu execute:

```bash
sudo apt-get update && apt-get install python-numpy python-scipy cython python-biopython python-nose python-nose2 python-mpi4py
```

**Make sure all libraries are from the Python 2.7.x series.**

#### 3.1.2. Installing ProDy
To install ProDy library, the simplest way is to use pip (you can use sudo to install it system-wide):

```bash
pip install -U ProDy
```

You may also need to install `pyparsing` dependency:

```bash
pip install pyparsing
```

More instructions on how to install it can be found in the official documentation (<http://prody.csb.pitt.edu/downloads/>).


#### 3.1.3. Installing FreeSASA (optional)
In case of using `cpydock` scoring function or to execute the tests, **Freesasa** library has to be installed and compiled with the python-binding options. Tested version in 
LightDock is 1.1 (<https://github.com/mittinatten/freesasa/tree/1.1>). To install freesasa 1.1, please follow these instructions (change `path/to/install`):

```bash
git clone https://github.com/mittinatten/freesasa.git
cd freesasa
git checkout tags/1.1
autoreconf -i
./configure --enable-python-bindings --prefix=path/to/install
make
make install
```

For more recent versions of freesasa, please check the instructions for installing it on its Github (<https://github.com/mittinatten/freesasa>). 

### 3.2. Download LightDock
The fastest way to install LightDock is to use `git` to clone the repository from GitHub:

```bash
git clone https://github.com/brianjimenez/lightdock.git
```

A directory called `lightdock` is now available. This is the path necessary for setting the enviroment variable `LIGHTDOCK_HOME` in your bash. Change your `~/.bashrc` accordingly (add the following lines to your `~/.bashrc`):

```bash
export LIGHTDOCK_HOME=/path/to/lightdock/folder
export PATH=$PATH:$LIGHTDOCK_HOME/bin:$LIGHTDOCK_HOME/bin/post:$LIGHTDOCK_HOME/bin/support
export PYTHONPATH=$PYTHONPATH:$LIGHTDOCK_HOME
```

### 3.3. Compilation of high-intensive calculation source code
Once the library dependencies are installed, a compilation of some high-intensive calculation parts is required. To do so, a script is provided:

```bash
cd ${LIGHTDOCK_HOME}/bin/setup
./setup.sh
```

### 3.4. Testing the framework (optional)
LightDock makes use of nosetests library for testing the different parts of the framework. There are two levels of testing: unitary and regression. 

Library unit tests. To run them execute:

```bash
cd $LIGHTDOCK_HOME
./run_tests.sh lib
```

Regression short tests. To run them execute:

```bash
cd $LIGHTDOCK_HOME
./run_tests.sh reg
```

Regression long tests. To run them execute (this may take several minutes):

```bash
cd $LIGHTDOCK_HOME
export LIGHTDOCK_LONG_TEST=true
./run_tests.sh reg
```

**NOTE**: some tests may fail depending on float accuracy in your installation. This is probably not relevant unless you plan to run lightdock using that piece of code. Please, open an issue in this GitHub repository to get help.


## 4. Documentation

This GitHub [README.md](README.md) file is intended as an installation guide.

The complete documentation about how to run the LightDock protocol can be found at:

<p align="center">
	<b><a href="https://brianjimenez.github.io/lightdock">https://brianjimenez.github.io/lightdock/</a>
	</b>
</p>


## 5. Get Help

LightDock is being actively developed and some issues may arise or you may get some extra help to run LightDock. In those cases, there are two main ways to get help:

1. Read the [FAQ](docs/FAQ.md) in case your problem is known
2. Open a [new issue in this repository](https://github.com/brianjimenez/lightdock/issues/new)
3. Or write an email to <b.jimenezgarcia@uu.nl>



