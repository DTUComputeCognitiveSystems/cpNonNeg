cpNonNeg
=======

This is an implementation of the CANDECOMP/PARAFAC model for tensor factorization of non-negative data. Missing values are handled by marginalization, i.e. ignored during optimization. The code was used in a project at Technical University of Denmark as part of one of the authors master degree, that ultimately lead to the publication [Non-negative Tensor Factorization with missing data for the modeling of gene expressions in the Human Brain](http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=6958919&tag=1 "NTF for missing data"). 
All implementation was done in MATLAB.

Inlcudes:
-------

* cpNonNeg.m - main MATLAB function
* cpNonNeg_sub.m - NMF solver for CP-subproblem
* krprod.m - Kathri-Rao product for tensors
* matricizing.m - matricizing operation
* tmult.m - tensor multiplication (mode specific)
* unmatricizing.m - tensor reconstruction from matrix

Written by: Søren Føns Vind Nielsen and Morten Mørup
CogSys, Technical University of Denmark, May 2014
