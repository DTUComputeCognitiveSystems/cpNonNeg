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