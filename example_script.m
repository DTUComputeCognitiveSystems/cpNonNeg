% Example script for using the CP-model with non-negativity constraints 
% Written by Søren Føns Vind Nielsen (sfvn at dtu dot dk)
% - Don't worry that the output looks weird...the input is weird as well

%% Generate synthetic data
D_true = 5;
N = [1000 50 25]; % Tensor dimensions
Nx = length(N);
F = cell(Nx,1);
for i = 1:Nx
        F{i} = rand(N(i),D_true);
end

% Diagonal identity tensor
I=zeros(D_true*ones(1,Nx));
for j=1:D_true
        I(j,j,j)=1;
end

Y=tmult(I,F{1},1);
% Data tensor
for ip = 2:Nx
        Y=tmult(Y,F{ip},ip);
end

sig2 = 0.5; % noise level
C = 5; % affine transformation to ensure non-negatitivty
X = Y + sqrt(sig2)*randn(N) + C*ones(N);

assert(min(X(:))>0);

%% Holdout missing data
p = 0.20; % holdout fraction (missing data)
NE = prod(size(X)); % number of elements in tensor
R = rand(NE,1)>(1-p); % holdout logical indices
X(R) = nan; % missing values are treated as NaN

%% Model specification
D = 5; % number of latent componenents in the model
Finit = cell(Nx,1); % initialization of factors (default)
scale = nanstd(X(:)); % scale of data

for i = 1:Nx
   Finit{i}=(scale.^(1/Nx))*rand(N(i),D); 
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