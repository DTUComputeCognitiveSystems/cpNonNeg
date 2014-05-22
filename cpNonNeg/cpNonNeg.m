function [FACT,SSEv,CPUt]=cpNonNeg(X,D,Finit,options)
% The CandeComp/PARAFAC with Non-Negativity constraints
% 
% The function fits a CP-model to the tensor data X with number of
% components equal to D, with inital factors Finit. The subproblem in this
% ALS approach is solved by NMF, with either Multiplicative Updates (MU) or
% Hierarchical Alternating Least Squares (HALS) (or a combination). Script
% is optimized for missing values. 
%
% Usage:
% [FACT,SSEv,CPUt] = cpNonNeg(X,D[,Finit,options])
%
% Input:
% X             n-way array to decompose. If data includes missing values
%               these should be NaN.
% D             integer number of components
% Finit         cell-array of length n with initial cp-factor matrices
%               (optional). Default is uniformly distributed random
%               matrices scaled accordin to n'th root of standard deviation
%               of data.
% options       struct with the following fields (optional)
%                   - options.maxiter: Number of outer iterations by als
%                   (default is 250)
%                   - options.mu: binary - if 1 MU steps
%                   are taken in the NMF-subproblem (default is 0)
%                   - options.hals: binary - if 1 HALS steps
%                   are taken in the NMF-subproblem (default is 1)
% Output:
% FACT          cell array: FACT{i} is the factors found for the i'th
%               modality
% SSEv          vector with SSE in each iteration
% CPUt          vector with CPU time after start in each iteration
%
% Written by: Søren Føns Vind Nielsen and Morten Mørup
% CogSys, Technical University of Denmark, May 2014

%% Initialization
global X2 scale

if nargin<4
    options.maxiter = 250;
    options.hals = 1;
    options.mu = 0;
end
hals = options.hals;
mu = options.mu;
maxiter = options.maxiter;

% Missing data
R = ~isnan(X); % index of present data
X(~R) = 0;

Nx=ndims(X);
N=size(X);
scale = std(X(:)); % Proper scaled initialization is essential for convergence
Xm = cell(Nx,1);

% Initial factors
for i=1:Nx
    if nargin < 3
        FACT{i}=(scale.^(1/Nx))*rand(N(i),D);
    else
        FACT{i} = Finit{i};
    end
    Xm{i}=matricizing(X,i);
    Rm{i}=matricizing(R,i);
end


SST=sum(X(:).^2); X2 = SST;
SSE=inf;
dSSE=inf;
SSEv=[];
CPUt=[];

disp([' '])
disp(['Non Negativity CP optimization'])
disp(['A ' num2str(D) ' component model will be fitted']);
disp([' '])
disp(['To stop algorithm press control C'])
disp([' ']);
dheader = sprintf('%12s | %12s | %12s | %12s ','Iteration','Expl. var.','dSSE','Time');
dline = sprintf('-------------+--------------+--------------+--------------+');

start = tic;
tic;
iter=0;
while dSSE>=1e-6*SSE && iter<maxiter
    if mod(iter,100)==0
        disp(dline); disp(dheader); disp(dline);
    end
    
    iter=iter+1;
    SSE_old=SSE;
    for i=1:Nx
        ind=1:Nx;
        ind(i)=[];
        kr=FACT{ind(1)};
        for z=ind(2:end)
            kr=krprod(FACT{z}, kr);
        end
        % Minimize cost function wrt. to mode using NMF
        if D>=50 && N(i)>30000 % Memory problems - split up in two subproblems
            nn = N(i)/2;
            [FACT{i}(1:nn,:), SSE1]=cpNonNeg_sub(Xm{i}(1:nn,:),FACT{i}(1:nn,:),kr,Rm{i}(1:nn,:),iter==1,hals,mu);
            [FACT{i}(nn+1:end,:), SSE2]=cpNonNeg_sub(Xm{i}(nn+1:end,:),FACT{i}(nn+1:end,:),kr,Rm{i}(nn+1:end,:),iter==1,hals,mu);
            SSE = SSE1 + SSE2;
        else
            [FACT{i}, SSE]=cpNonNeg_sub(Xm{i},FACT{i},kr,Rm{i},iter==1,hals,mu);
        end
    end   
    dSSE=SSE_old-SSE;
    CPUt = [CPUt toc(start)];
    SSEv = [SSEv SSE];
    if mod(iter,5)==0
        disp(sprintf('%12.0f | %12.4f | %12.4e | %12.4e |',iter, (SST-SSE)/SST,dSSE,toc));
        tic;
    end
end
% Display final iteration
disp(sprintf('%12.0f | %12.4f | %12.4e | %12.4e |',iter, (SST-SSE)/SST,dSSE,toc));
disp(['Total time [s]: ' num2str(toc(start))])


%eof
end