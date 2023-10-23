# Simulation scripts

This directory contains scripts used to run the simulation suites used in the
thesis.

## Well-mixed simulations

The well-mixed simulations are executed by the `run-well-mixed-xxx.py` scripts.
The working directory should be the same as the directory containing the
simulation script.  For example, to run the simulations in
`ch4-effective-g1-duration`:

```
$ cd sims/thesis/ch4-effective-g1-duration
$ python3 run-well-mixed-effective-g1-duration.py
```

To specify a particular seed, provide a non-negative integer as the first
command-line argument.  For instance, to use seed zero:

```
$ python3 run-well-mixed-effective-g1-duration.py 0
```

Furthermore, we can invoke the script with a range to simulate only a subset
of simulations.

```
$ python3 run-well-mixed-effective-g1-duration.py 0 0 10
```

The second and third arguments specify the start and the end of the range,
respectively.  Note that we use zero-based indexing and specify the range as a
half-open interval.  Hence, the command above runs simulations 0 though 9.

Splitting up the simulation suite in this way allows for efficient
parallellisation.

## Vertex simulations

The vertex simulations are executed by the `run-vertex-xxx.sh` scripts.  The
working directory should be the same as the directory containing the simulation
script.  The script accepts a simulation number as the argument.  For example,
to run simulation 69 in `ch2-mechanical-study`:

```
$ cd sims/thesis/ch2-mechanical-study
$ bash run-vertex-mechanical-study.sh 69
```
