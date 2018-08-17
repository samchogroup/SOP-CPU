# SOP-CPU - modified on 04/12/2017

Cell array and Hybrid method implementation by Sajant Anand under the direction of Prof. Sam Cho. The Cell Array method breaks the simulation space into cubic cells and places each bead into one of these. We then only calculate the interaction of beads in neighboring cells each simulation iteration. For the Hybrid method, we use the same cell array structure, but we place interactions that are in neighboring cells into the neighbor list. We then use this neighbor list for X simulation interations to calculate interactions and energy updates (as is done in the neighbor list method).

The Benchmark folder contains both the input files necessary to run timings tests and timing files with the results used in the Anand 2018 CS Thesis. To run these tests, move the Benchmark folder to the same parent directory of the SOP folder and then run 'runme.sh'.
