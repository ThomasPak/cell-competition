# Simulation scripts

This directory contains scripts used to run the simulation suites used in the
thesis.

## Well-mixed simulations

The well-mixed simulations are executed by the `run-well-mixed-xxx.py` scripts.
The working directory should be the same as the directory containing the
simulation script.  For example, to run the simulations in
`ch4-effective-g1-duration`:

```
$ cd ch4-effective-g1-duration
$ python3 run-well-mixed-effective-g1-duration.py
```

## Vertex simulations

The vertex simulations are executed by the `run-vertex-xxx.sh` scripts.  The
working directory should be the same as the directory containing the simulation
script.  The script accepts a simulation number as the argument.  For example,
to run simulation 69 in `ch2-mechanical-study`:

```
$ cd ch2-mechanical-study
$ bash run-vertex-mechanical-study.sh 69
```
