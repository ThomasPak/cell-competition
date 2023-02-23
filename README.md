# Modelling cell competition

This repository contains the code developed and used in (i) the thesis
**Modelling cell competition**, submitted by Thomas Pak for the degree of
*Doctor of Philosophy* at the University of Oxford (Trinity 2021), and (ii) the
paper **A mathematical framework for the emergence of winners and losers in
cell competition** (review pending).

* The `well_mixed` and `vertex` directories contain the code for the well-mixed
  and vertex models, respectively.
* The `sims` directory contain the scripts used to run the simulation suites
  described in the thesis and in the paper.
* The `utils` directory contains scripts to collate the data generated by the
  well-mixed and vertex simulations.
* The `data` directory contains the data generated by simulations and used for
  generating the thesis and paper figures.
* The `figures` directory contains the scripts used to generate the thesis and
  paper figures.
* The `Makefile` file contains commands for generating the thesis and paper
  figures.  To generate the thesis and paper figures, simply run the command
  `make`.  The generated figures will be placed in the new directory `images`.

## Well-mixed model

The file `well_mixed/well_mixed_death_clock.py` implements the class
`WellMixedSimulator` for running simulations of the well-mixed death clock
model, and the class `WellMixedSimulationData` for analysing the raw data
generated by well-mixed simulations.

### Prerequisites

* [Python 3](https://www.python.org/)
* [NumPy](https://numpy.org/)
* [SciPy](https://scipy.org/)
* [pandas](https://pandas.pydata.org/)

### Testing

Run `well_mixed/well_mixed_death_clock.py` as a script to run a test suite,
comprising of unit tests and simulation tests, that verifies whether the
implemented classes function as expected:

```
$ python3 well_mixed/well_mixed_death_clock.py
```

If this command completes without raising errors, then the tests were
successful.

## Vertex model

The directory `vertex` is a Chaste user project (see
[https://chaste.cs.ox.ac.uk/trac/wiki/ChasteGuides/UserProjects](https://chaste.cs.ox.ac.uk/trac/wiki/ChasteGuides/UserProjects)
for information on Chaste user projects).  This means that the directory
`vertex` needs to be linked into an existing Chaste installation, and that the
vertex model is built as part of a Chaste build.

In the instructions below, we assume that the environment variables
`$CHASTE_SRC` and `$CHASTE_BUILD` are set to the Chaste source and build
directories, respectively.

### Prerequisites

* [Chaste](https://chaste.cs.ox.ac.uk/trac/wiki/GettingStarted)

Please note that this code has only been tested for the **Chaste 2021.1
release**.  Ensure that you are using the Chaste 2021.1 release with the
following commands.

```
$ cd $CHASTE_SRC
$ git checkout 2021.1
```

### Patches

We implemented an ad hoc fix for handling the case where one of two
neighbouring triangular cells on the boundary performs a T2 swap.

We also implemented an ad hoc fix for handling the case where a T1 swap results
in a concave triangular element.

These ad hoc fixes were included as a patch to the Chaste source code, which
can be applied as follows:

```
$ cd $CHASTE_SRC
$ git apply /path/to/cell-competition/vertex/patch/chaste-2021.1-ad-hoc-fixes.patch
```

### Linking the user project

Link the user project directory to the `projects` directory in Chaste so that
CMake can find the user project in the build step.

```
$ ln -s /path/to/cell-competition/vertex $CHASTE_SRC/projects/cell-competition
```

### Building

After running CMake, the target `project_cell-competition` builds all the
executables and tests related to the user project.

```
$ cd $CHASTE_BUILD
$ cmake $CHASTE_SRC
$ make project_cell-competition
```

### Testing

After building Chaste with the `cell-competition` user project, run CTest with
the label `Continuous_project-cell-competition` to test the vertex model build.

```
$ cd $CHASTE_BUILD
$ ctest -L Continuous_project_cell-competition
```
