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
 
Example function call:
--------
The following script (available in the repository) shows a basic usage of the code. 
NB! The code will terminate quickly and not neccessarily give meaningful results (due to the random data).
```matlab
% example script
X = randn(100,100,10); % 100 x 100 x 10 normally-random tensor
ND = ndims(X);
N = size(X);
D = 5; % Latent factors
Finit = cell(ND,1); % initialization of factors (default)
scale = std(X(:)); % scale of data

for i = 1:ND
   Finit{i}=(scale.^(1/ND))*rand(N(i),D); 
end

% options
options.maxiter = 250; % number of iterations
options.mu = 0; % no multiplicative update steps are taken
options.hals = 1; % hierarchical alternating least sqaures steps are taken

% run
[FACT,SSEv,CPUt]=cpNonNeg(X,D,Finit,options);

% FACT gives back factors in a cell array just as Finit was initialized
% SSEv is the sum of squared (reconstruction) errors in each iteration
% (vector)
% CPUt is the CPU time used in each iteration
```

Written by: Søren Føns Vind Nielsen and Morten Mørup
CogSys, Technical University of Denmark, May 2014
