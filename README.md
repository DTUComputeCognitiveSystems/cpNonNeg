cpNonNeg
=======

This is an implementation of the CANDECOMP/PARAFAC model for tensor factorization of non-negative data. Missing values are handled by marginalization, i.e. ignored during optimization. 
All implementation was done in MATLAB.

Inlcudes:
-------

* cpNonNeg.m - main MATLAB function
* cpNonNeg_sub.m - NMF solver for CP-subproblem
* krprod.m - Kathri-Rao product for tensors
* matricizing.m - matricizing operation
* tmult.m - tensor multiplication (mode specific)
* unmatricizing.m - tensor recontruction from matrix
* faces.mat - test data set (2-way face database - 19x19 pictures)
* faces_tensor.mat - same test set as above but 3-way

Written by: Søren Føns Vind Nielsen and Morten Mørup
CogSys, Technical University of Denmark, May 2014
