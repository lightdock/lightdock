# LightDock
## Protein-protein docking framework
Protein Interactions and Docking Group <http://life.bsc.es/pid/pidweb/>

Life Sciences Department - Barcelona Supercomputing Center <http://www.bsc.es/>

## 1. Introduction
LightDock is a protein-protein docking framework based on the Glowworm Swarm Optimization (GSO) algorithm.

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
* Freesasa (only if PyDock scoring function is used, <http://freesasa.github.io/>)

NumPy, Scipy, Cython, Biopython, Nose and MPI4py libraries are usually available as packages in most of GNU/Linux distributions. To install them in Ubuntu execute:

`sudo apt-get update && apt-get install python-numpy python-scipy cython python-biopython python-nose2 python-mpi4py`

Make sure all libraries are from Python 2.7.x series.

To install ProDy library, the simplest way is to use pip:

`pip install -U ProDy`

More instructions on how to install it can be found in the official documentation (<http://prody.csb.pitt.edu/downloads/>).

In case of using PyDock scoring function, **Freesasa** library has to be installed and compiled with the python-binding options.
Please, check the instructions for installing it on its Github (<https://github.com/mittinatten/freesasa>). Tested version in 
LightDock is 1.1 (<https://github.com/mittinatten/freesasa/tree/1.1>).

### 2.2. Compilation of high-intensive calculation source code
Once the library dependencies are installed, a compilation of some high-intensive calculation parts is required. To do so, a script is provided:

```
cd bin/setup
./setup.sh
```

### 2.3. Testing the framework
LightDock makes use of nosetests library for testing the different parts of the framework. There are two levels of testing: unitary and regression. Before running the tests, please add LightDock folder to your path environment variable. 

In bash (**complete the LIGHTDOCK_HOME variable according to the path where you have downloaded and unpacked LightDock**):

```
export LIGHTDOCK_HOME=/path/to/lightdock/folder
export PATH=$PATH:$LIGHTDOCK_HOME/bin:$LIGHTDOCK_HOME/lightdock/bin/post
export PYTHONPATH=$PYTHONPATH:$LIGHTDOCK_HOME/lightdock/:$LIGHTDOCK_HOME/lightdock/lightdock
```

Then you will be able to run the tests. To do so:

Library unit tests:

```
./run_tests.sh lib
```

Regression short tests:

```
./run_tests.sh reg
```

Regression long tests (this may take several minutes):

```
export LIGHTDOCK_LONG_TEST=true
./run_tests.sh reg
```

## 3. Executing LightDock
The simplest way to perform a protein-protein docking in LightDock is to use default parameters and to only provide two [PDB](http://www.rcsb.org/pdb/static.do?p=file_formats/pdb/index.html) files for both receptor and ligand proteins.

### 3.1. Simplest example
You fill find in the [examples](examples/) a folder [2UUY](examples/2UUY/) with a PDB file for the receptor [2UUY_rec.pdb](examples/2UUY/2UUY_rec.pdb) and the ligand [2UUY_lig.pdb](examples/2UUY/2UUY_lig.pdb).

Execute LightDock with the default parameters and only calculating 1 initial glowworm group with 10 glowworms and for 5 steps of the algorithm:

<code>
lightdock 2UUY_rec.pdb 2UUY_lig.pdb 1 10 5
</code>

By default, LightDock makes use of the DFIRE scoring function. The output of the execution should be something similar to this:

```
[lightdock] INFO: Parameters are:
[lightdock] INFO: Receptor: 2UUY_rec.pdb
[lightdock] INFO: Ligand: 2UUY_lig.pdb
[lightdock] INFO: Number of clusters: 1
[lightdock] INFO: Number of glowworms per cluster: 10
[lightdock] INFO: Simulation steps: 5
[lightdock] INFO: GSO seed: 324324
[lightdock] INFO: Starting points seed: 324324
[lightdock] INFO: Translation step: 0.5
[lightdock] INFO: Rotation step: 0.5
[lightdock] INFO: lightdock parameters saved to ./lightdock.info
[lightdock] INFO: Reading 2UUY_rec.pdb receptor PDB file...
[lightdock] INFO: 1628 atoms, 223 residues read.
[lightdock] INFO: Reading 2UUY_lig.pdb ligand PDB file...
[lightdock] INFO: 415 atoms, 55 residues read.
[lightdock] INFO: Saving processed structures to PDB files...
[lightdock] INFO: Done.
[lightdock] INFO: Calculating starting positions...
[lightdock] INFO: Generated 1 positions files
[lightdock] INFO: Done.
[lightdock] INFO: Loading scoring function...
[lightdock] INFO: Using DFIRE scoring function
[lightdock] INFO: Done.
[lightdock] INFO: Preparing environment
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
[kraken] INFO: folding tentacle Tentacle-3
[0] step 1
[kraken] INFO: folding tentacle Tentacle-4
[0] step 2
[0] step 3
[0] step 4
[0] step 5
[kraken] INFO: folding tentacle Tentacle-1
[kraken] INFO: 1 ships destroyed
[lightdock] INFO: Finished.
```

There is a warning on the number of CPU cores used. By default, LightDock will look for the total number of cores. If you want to specify a different number, use the flag <code>-c NUMBER_CORES</code>.

For each of the initial glowworm groups, there is a folder called <code>cluster_X</code>. In our example, we use only one initial glowworm group so there is only a folder <code>cluster_0</code>. Inside, we can find the file containing the result of the simulation: [gso_5.out](examples/2UUY/cluster_0/gso_5.out). 

In this file, every line corresponds to a glowworm agent in the algorithm:

```
#Coordinates  RecID  LigID  Luciferin  Neighbor's number  Vision Range  Scoring
(31.4171143,  1.8570079, -6.3956223, -0.1058407, -0.4849369,  0.5997430, -0.6276482)    0    0  10.79432165  0 2.200   7.52191884
```

Finally, to generate the final docked PDB structures, we will use the script **generate_conformations.py**:

```
generate_conformations.py 2UUY_rec.pdb 2UUY_lig.pdb cluster_0/gso_5.out 10 
@> ProDy is configured: verbosity='info'
[generate_conformations] INFO: Reading lightdock_2UUY_rec.pdb receptor PDB file...
[generate_conformations] INFO: 1628 atoms, 223 residues read.
[generate_conformations] INFO: Reading lightdock_2UUY_lig.pdb ligand PDB file...
[generate_conformations] INFO: 415 atoms, 55 residues read.
[generate_conformations] INFO: Read 10 coordinate lines
[generate_conformations] INFO: Generated 10 conformations
```

In the <code>cluster_0</code> folder we will find the 10 structures corresponding to the 10 glowworm agents used in the example.
