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
%% Generate synthetic data
N = [1000 50 25]; % Tensor dimensions
Nx = length(N);
F = cell(Nx,1);
for i = 1:Nx
        F{i} = rand(N(i),D);
end

% Diagonal identity tensor
I=zeros(D*ones(1,Nx));
for j=1:D
        I(j,j,j)=1;
end

Y=tmult(I,F{1},1);
% Data tensor
for ip = 2:Nx
        Y=tmult(Y,F{ip},ip);
end

sig2 = 0.1; % noise level
X = Y + sqrt(sig2)*randn(N);

%% Holdout missing data
p = 0.05; % holdout percentage (missing data)
NE = prod(size(X)); % number of elements in tensor
R = rand(NE,1)>(1-p); % holdout logical indices
X(R) = nan; % missing values are treated as NaN

%% Model specification
D = 5; % number of latent componenents in the model
Finit = cell(ND,1); % initialization of factors (default)
scale = std(X(:)); % scale of data

for i = 1:ND
   Finit{i}=(scale.^(1/ND))*rand(N(i),D); 
end

% options
options.maxiter = 250; % number of iterations
options.mu = 0; % no multiplicative update steps are taken
options.hals = 1; % hierarchical alternating least sqaures steps are taken

%% Run
[FACT,SSEv,CPUt]=cpNonNeg(X,D,Finit,options);

% FACT gives back factors in a cell array just as Finit was initialized
% SSEv is the sum of squared (reconstruction) errors in each iteration
% (vector)
% CPUt is the CPU time used in each iteration
```

Written by: Søren Føns Vind Nielsen and Morten Mørup
CogSys, Technical University of Denmark, May 2014
