![LightDock](docs/media/lightdock_banner.png "LightDock")

#### Table of Contents

- [Introduction](#1-introduction)
- [Installation](#2-installation)
- [LightDock for the impatient](#3-lightdock-for-the-impatient)
- [Documentation](#4-documentation)

## 1. Introduction
LightDock is a protein-protein docking framework based on the [Glowworm Swarm Optimization](https://link.springer.com/article/10.1007/s11721-008-0021-5) (GSO) algorithm.

**The framework is written in the Python programming language (version 2.7) and allows the users to incorporate their own scoring function.**

The LightDock framework is highly versatile, with many options that can be further developed and optimized by the users: it can accept any user-defined scoring function, can use local gradient-free minimization, the simulation can be restrained from the beginning to focus on user-assigned interacting regions, and it has support for the use of pre-calculated conformers for both receptor and ligand.

## 2. Installation
### 2.1. Dependencies
LightDock has the following dependencies:

* Python 2.7.x
* Nose (<http://nose.readthedocs.io/en/latest/>)
* NumPy (<http://www.numpy.org/>)
* Scipy (<http://www.scipy.org/>)
* Cython (<http://cython.org/>)
* BioPython (<http://biopython.org>)
* MPI4py (<http://pythonhosted.org/mpi4py/>)
* ProDy (<http://prody.csb.pitt.edu/>)
* Freesasa (only if `cpydock` scoring function is used and to execute the tests, <http://freesasa.github.io/>)

NumPy, Scipy, Cython, Biopython, Nose and MPI4py libraries are usually available as packages in most of GNU/Linux distributions. To install them in Ubuntu execute:

```bash
sudo apt-get update && apt-get install python-numpy python-scipy cython python-biopython python-nose2 python-mpi4py
```

**Make sure all libraries are from the Python 2.7.x series.**

To install ProDy library, the simplest way is to use pip (you can use sudo to install it system-wide):

```bash
pip install -U ProDy
```

More instructions on how to install it can be found in the official documentation (<http://prody.csb.pitt.edu/downloads/>).

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

### 2.2. Download LightDock
The fastest way to install LightDock is to use git to clone the repository in GitHub:

```bash
git clone https://github.com/brianjimenez/lightdock.git
```
A directory called `lightdock` is now available. This is the path for `LIGHTDOCK_HOME`. Change your `~/.bashrc` accordingly (add the following lines to your `~/.bashrc`):

```bash
export LIGHTDOCK_HOME=/path/to/lightdock/folder
export PATH=$PATH:$LIGHTDOCK_HOME/bin:$LIGHTDOCK_HOME/lightdock/bin/post:$LIGHTDOCK_HOME/lightdock/bin/support
export PYTHONPATH=$PYTHONPATH:$LIGHTDOCK_HOME/lightdock/
```

### 2.3. Compilation of high-intensive calculation source code
Once the library dependencies are installed, a compilation of some high-intensive calculation parts is required. To do so, a script is provided:

```bash
cd ${LIGHTDOCK_HOME}/bin/setup
./setup.sh
```

### 2.4. Testing the framework
LightDock makes use of nosetests library for testing the different parts of the framework. There are two levels of testing: unitary and regression. 

Library unit tests:

```bash
cd $LIGHTDOCK_HOME
./run_tests.sh lib
```

Regression short tests:

```bash
cd $LIGHTDOCK_HOME
./run_tests.sh reg
```

Regression long tests (this may take several minutes):

```bash
cd $LIGHTDOCK_HOME
export LIGHTDOCK_LONG_TEST=true
./run_tests.sh reg
```

## 3. LightDock for the impatient
The simplest way to perform a protein-protein docking in LightDock is to use default parameters and to only provide two [PDB](http://www.rcsb.org/pdb/static.do?p=file_formats/pdb/index.html) files for both receptor and ligand proteins.

### 3.1. Simplest example
You fill find in the [examples](examples/) a folder [2UUY](examples/2UUY/) with a PDB file for the receptor [2UUY_rec.pdb](examples/2UUY/2UUY_rec.pdb) and the ligand [2UUY_lig.pdb](examples/2UUY/2UUY_lig.pdb).

In previous versions of LightDock, a setup step was not required. This has changed from **version 0.5.0** and now a simulation setup is required. Follow the next steps in order to perform your first protein-protein docking with LightDock:

#### 3.1.1. Copying data
Create a directory and copy the sample data provided:

```bash
cd $LIGHTDOCK_HOME
cd examples
mkdir test
cd test
cp ../2UUY/2UUY*.pdb .
```

#### 3.1.2. LightDock setup
Execute `lightdock_setup` script in order to prepare your LightDock simulation:

```bash
lightdock_setup 2UUY_rec.pdb 2UUY_lig.pdb 1 10
```

The previous command executes LightDock with the default parameters and only calculating 1 swarm containing 10 glowworms.

Here it is the output:

```
@> ProDy is configured: verbosity='info'
[lightdock_setup] INFO: Reading structure from 2UUY_rec.pdb PDB file...
[lightdock_setup] INFO: 1628 atoms, 223 residues read.
[lightdock_setup] INFO: Reading structure from 2UUY_lig.pdb PDB file...
[lightdock_setup] INFO: 415 atoms, 55 residues read.
[lightdock_setup] INFO: Calculating reference points for receptor 2UUY_rec.pdb...
[lightdock_setup] INFO: Done.
[lightdock_setup] INFO: Calculating reference points for ligand 2UUY_lig.pdb...
[lightdock_setup] INFO: Done.
[lightdock_setup] INFO: Saving processed structure to PDB file...
[lightdock_setup] INFO: Done.
[lightdock_setup] INFO: Saving processed structure to PDB file...
[lightdock_setup] INFO: Done.
[lightdock_setup] INFO: Calculating starting positions...
[lightdock_setup] INFO: Generated 1 positions files
[lightdock_setup] INFO: Done.
[lightdock_setup] INFO: Preparing environment
[lightdock_setup] INFO: Done.
[lightdock_setup] INFO: LightDock setup OK 
```

**A new file called** `setup.json` **has been generated with the simulation information** .

If you execute `lightdock_setup` without arguments a little help is displayed:

```bash
lightdock_setup
@> ProDy is configured: verbosity='info'
usage: lightdock_setup [-h] [--seed_points STARTING_POINTS_SEED]
                       [-ft ftdock_file] [--noxt] [-anm] [--seed_anm ANM_SEED]
                       [-anm_rec ANM_REC] [-anm_lig ANM_LIG]
                       receptor_pdb_file ligand_pdb_file swarms glowworms
lightdock_setup: error: too few arguments
```

#### 3.1.3. LightDock run
Execute `lightdock` script in order to run your first LightDock simulation:

```bash
lightdock setup.json 10
```

You will see an output like this:

```
@> ProDy is configured: verbosity='info'
[lightdock] INFO: simulation parameters saved to ./lightdock.info
[lightdock_setup] INFO: Reading structure from 2UUY_rec.pdb PDB file...
[lightdock_setup] INFO: 1628 atoms, 223 residues read.
[lightdock_setup] INFO: Reading structure from 2UUY_lig.pdb PDB file...
[lightdock_setup] INFO: 415 atoms, 55 residues read.
[lightdock] INFO: Loading scoring function...
[lightdock] INFO: Using DFIRE scoring function
[lightdock] INFO: Done.
[kraken] WARNING: Number of cores has not been specified or is incorrect. Using available cores.
[kraken] INFO: Kraken has 4 tentacles (cpu cores)
[kraken] INFO: Tentacle ready with 1 tasks
[kraken] INFO: Tentacle ready with 0 tasks
[kraken] INFO: Tentacle ready with 0 tasks
[kraken] INFO: Tentacle ready with 0 tasks
[kraken] INFO: 1 ships ready to be smashed
[lightdock] INFO: Monster spotted
[kraken] INFO: Release the Kraken!
[kraken] INFO: folding tentacle Tentacle-2
[0] step 1
[kraken] INFO: folding tentacle Tentacle-4
[kraken] INFO: folding tentacle Tentacle-3
[0] step 2
[0] step 3
[0] step 4
[0] step 5
[0] step 6
[0] step 7
[0] step 8
[0] step 9
[0] step 10
[kraken] INFO: folding tentacle Tentacle-1
[kraken] INFO: 1 ships destroyed
[lightdock] INFO: Finished.
```

By default, LightDock makes use of the DFIRE scoring function. There is a warning on the number of CPU cores used. By default, LightDock will look for the total number of cores. If you want to specify a different number, use the flag <code>-c NUMBER_CORES</code>. **Note that MPI is also supported using the -mpi flag**.

For each of the swarms, there is a folder called <code>swarm_X</code>. In our example, we use only one swarm so there is only a folder <code>swarm_0</code>. Inside, we can find the file containing the result of the simulation `gso_10.out`. See an example of this file from the *2UUY* example: [gso_5.out](examples/2UUY/cluster_0/gso_5.out). 

In this file, every line corresponds to a glowworm agent in the algorithm:

```
#Coordinates  RecID  LigID  Luciferin  Neighbor's number  Vision Range  Scoring
(31.4171143,  1.8570079, -6.3956223, -0.1058407, -0.4849369,  0.5997430, -0.6276482)    0    0  10.79432165  0 2.200   7.52191884
```

Finally, to generate the final docked PDB structures, we will use the script **lgd_generate_conformations.py**:

```bash
cd swarm_0
lgd_generate_conformations.py ../2UUY_rec.pdb ../2UUY_lig.pdb gso_10.out 10
```

Output is:

```
@> ProDy is configured: verbosity='info'
[generate_conformations] INFO: Reading ../lightdock_2UUY_rec.pdb receptor PDB file...
[generate_conformations] INFO: 1628 atoms, 223 residues read.
[generate_conformations] INFO: Reading ../lightdock_2UUY_lig.pdb ligand PDB file...
[generate_conformations] INFO: 415 atoms, 55 residues read.
[generate_conformations] INFO: Read 10 coordinate lines
[generate_conformations] INFO: Generated 10 conformations
```

Inside the <code>swarm_0</code> folder 10 new PDB structures corresponding to the 10 glowworm agents used in the example have been generated.


## 4. Documentation

This GitHub [README.md](README.md) file is intended as an installation guide and to show the simplest protein-protein docking example.

The complete information can be found at [http://brianjimenez.github.io/lightdock](http://brianjimenez.github.io/lightdock)

