# Code for "On the emergence of single versus multi-state allostery". 
## Eric Rouviere, Rama Ranganathan, Olivier Rivoire.



### Installation instructions (Mac or Linux).

All code is written in julia1.7.2.  Download and install julia1.7.2 at from

https://julialang.org/downloads/oldreleases/

Next, open up a terminal type,

`julia`

A julia REPL should be running. To install the needed Julia packages, type,

`] add PyPlot, LaTeXStrings, Revise`

kill that REPL session.

### Evolving elastic networks

The julia notebook `evolve_single_net.ipynb` contains the basic functions to build a network, evolve it under a selection for allostery, and analyze the result. All code is sourced from the directory `src/`.


