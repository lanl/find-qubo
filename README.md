find-qubo
=========

Description
-----------

The goal of a [*quadratic unconstrained binary optimization*](https://en.wikipedia.org/wiki/Quadratic_unconstrained_binary_optimization) (QUBO) problem is to find the values of the <i>x<sub>i</sub></i> ∈ [0, 1] that minimize the quadratic pseudo-Boolean expression

![qubo](https://bit.ly/3rULTP0)

where <i>a<sub>i</sub></i> ∈ ℝ, and <i>b<sub>i,j</sub></i> ∈ ℝ are given.  There is always one minimum, but there can be more.  For example, <i>f</i>(<b><i>x</i></b>) = 3<i>x</i><sub>3</sub> + <i>x</i><sub>1</sub><i>x</i><sub>2</sub> − 2<i>x</i><sub>1</sub><i>x</i><sub>3</sub> − 2<i>x</i><sub>2</sub><i>x</i><sub>3</sub> is minimized at {<i>x</i><sub>1</sub>, <i>x</i><sub>2</sub>, <i>x</i><sub>3</sub>} = {0, 0, 0}, {0, 1, 0}, {1, 0, 0}, and {1, 1, 1}:

| <i>x</i><sub>1</sub> | <i>x</i><sub>2</sub> | <i>x</i><sub>3</sub> | <i>f</i>(<b><i>x</i></b>) |
| :------------------: | :------------------: | :------------------: | :-----------------------: |
|                    0 |                    0 |                    0 |                         0 |
|                    0 |                    0 |                    1 |                         3 |
|                    0 |                    1 |                    0 |                         0 |
|                    0 |                    1 |                    1 |                         1 |
|                    1 |                    0 |                    0 |                         0 |
|                    1 |                    0 |                    1 |                         1 |
|                    1 |                    1 |                    0 |                         1 |
|                    1 |                    1 |                    1 |                         0 |

A surprising number of problems can be expressed in this form, and it is the native input format for specialized computers based on classical or quantum annealing (such as the those developed by [D-Wave Systems](https://www.dwavesys.com/)).  One nice property of QUBOs is that they are additive in the sense that <i>f</i>(<b><i>x</i></b>) + <i>g</i>(<b><i>x</i></b>) is minimized at the intersection of the set of  <b><i>x</i></b> where <i>f</i>(<b><i>x</i></b>) is minimized and the set of  <b><i>x</i></b> where <i>g</i>(<b><i>x</i></b>) is minimized, assuming this intersection is non-empty.

The challenge that **find-qubo** addresses is the inverse problem: Given a set of <b><i>x</i></b> that should be minimized, return a corresponding <i>f</i>(<b><i>x</i></b>).  What makes this difficult is that many problems are not directly solvable.  In this case, additional ("ancillary") <i>x<sub>i</sub></i> variables must be introduced to produce a valid <i>f</i>(<b><i>x</i></b>).

The importance of a tool like **find-qubo** is that it can be used to construct primitive QUBOs that can then be combined to formulate more complicated optimization problems.

Installation
------------

Multiple implementations of **find-qubo** are provided in this repository:

| Tool            | Language | Features | Non-standard dependencies |
| :-------------- | :------- | :------- | :------------------------ |
| `find-qubo`     | Go       | Choice of fast heuristic or slow brute force, [MPI](https://en.wikipedia.org/wiki/Message_Passing_Interface) parallelism | [CLP](https://www.coin-or.org/Clp/) |
| `find-qubo-pm`  | Python   | Faster than brute force | [Ocean](https://ocean.dwavesys.com/) |
| `find-qubo-smt` | Python   | Very fast on modest-sized problems | [Z3](https://github.com/Z3Prover/z3) |

(In the above, `pm` stands for [penalty model](https://docs.ocean.dwavesys.com/projects/penaltymodel/), and `smt` stands for [satisfiability modulo theories](https://en.wikipedia.org/wiki/Satisfiability_modulo_theories).)

All three implementations iterate over the number of ancillary variables, first trying zero, then one, then two, and so forth.  Each variable added makes the search for a solution exponentially slower.  The first implementation, `find-qubo`, is the most experimental.  It investigates a heuristic that often finds a solution extremely fast but may overlook a solution at a given number of ancillae and require more ancillae than necessary, becoming  very slow.  For comparison purposes, it also provides a brute-force mode, which can find either a single solution or a tally of all solutions (to gauge problem difficulty) and is extremely slow.  `find-qubo-smt` is the most recent implementation and can be very fast on problems that do not contain too many variables.

**find-qubo** provides a `Makefile` for building `find-qubo`.  Remove `-tags mpi` from the `GOFLAGS` line to disable support for MPI-based parallelism.  `find-qubo-pm` and `find-qubo-smt` do not need to be compiled, but their dependencies must be met before they can run.  Use the [`pip`](https://docs.python.org/3/installing/index.html) command to install [`dwave-ocean-sdk`](https://pypi.org/project/dwave-ocean-sdk/) (for `find-qubo-pm`) or [`z3-solver`](https://pypi.org/project/z3-solver/) (for `find-qubo-smt`).

Usage
-----

All three implementations require the name of a file containing a truth table.  (See [the `testcases` directory](testcases/) for examples.)  `find-qubo` accepts a large number of options, and `find-qubo-smt` accepts a few.  Run with `--help` for information.

Citation
--------

The following peer-reviewed publication introduces **find-qubo**:

> Scott Pakin, "A Simple Heuristic for Expressing a Truth Table as a Quadratic Pseudo-Boolean Function".  Proceedings of the 2021 IEEE International Conference on Quantum Computing and Engineering (QCE), Broomfield, Colorado, USA, 17–22 October 2021, pp. 218–224.  ISBN: 978-1-6654-1691-7, DOI: [10.1109/QCE52317.2021.00039](https://doi.org/10.1109/QCE52317.2021.00039).

Legal statement
---------------

Copyright © 2021 Triad National Security, LLC.
All rights reserved.

This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S.  Department of Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.

This program is open source under the [BSD-3 License](LICENSE.md).  Its LANL-internal identifier is C21038.

Author
------

Scott Pakin, *pakin@lanl.gov*
