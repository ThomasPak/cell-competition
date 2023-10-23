# Paper simulation scripts

This directory contains scripts used to run the simulation suites used in the
paper.  The scripts are grouped by the corresponding paper section, e.g.
the scripts for the simulation suite described in Section 4.1.2 are found
in `section-4.1.2`.

## Well-mixed simulations

The well-mixed simulations are executed by the `run-well-mixed-xxx.py` scripts.
The working directory should be the same as the directory containing the
simulation script.  For example, to run the simulations in
`section-4.1.2`:

```
$ cd section-4.1.2
$ python3 run-well-mixed-exponential-proliferation-regimes.py
```

## Vertex simulations

The vertex simulations are executed by the `run-vertex-xxx.sh` scripts.  The
working directory should be the same as the directory containing the simulation
script.  The script accepts a simulation number as the argument.  For example,
to run simulation 69 in `section-2.2`:

```
$ cd section-2.2
$ bash run-vertex-mechanical-latin-study.sh 69
```

## Note for Windows users

Most of the files in this directory are Linux-type symbolic links to
directories in `/sims/thesis`.  These do not work in a Windows environment.
Instead, please navigate to the relevant directory using the following table.

| Symbolic link   | Target directory                                     |
|-----------------|------------------------------------------------------|
| `section-s4`    | `/sims/thesis/ch5-exponential-survival-probability`  |
| `section-4.1.2` | `/sims/thesis/ch5-exponential-proliferation-regimes` |
| `section-s6`    | `/sims/thesis/ch6-heterotypic-survival-difference`   |
| `section-4.2.3` | `/sims/thesis/ch6-heterotypic-proliferation-regimes` |
