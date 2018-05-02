function [W, SSE] = cpNonNeg_sub(Xn,W,Ht,R,initialize,hals,mu,maxiter)
% Solves the NMF subproblem for the Non-Negative constrained CP-model
% See cpNonNeg.m for documentation
global X2 scale
%% Initialization
[~,D] = size(W);
[N,M] = size(Xn);

if nargin < 8 || isempty(maxiter)
	maxiter=M;
end

tol = 1e-9;
eps = 1e-12; % Regularization

k = 0; % while loop counter

%% main algorithm
% indexing for use in updates
h2ind = [1];
for i = 2:D
    h2ind(i) = h2ind(i-1)+ D-i+2;
end

IND = nan(D,D-1);
IND(1,:) = 2:D;
for i = 2:D
    if i == D
        IND(D,:) = IND(1:(D-1),D-1)';
    else
        IND(i,1:(i-1)) = IND(1:(i-1),i-1)';
        IND(i,i:end) = h2ind(i)+1:h2ind(i+1)-1;
    end
end

idxGp = nan(D);
idx_not_d = nan(D,D-1);
for d=1:D
    idx_t = 1:D; idx_t(d) = [];      
    idx_not_d(d,:) = idx_t; 
    idxGp(d,:) = sort([h2ind(d) IND(d,:)],'ascend');
end

% Premultiplication
XHt = Xn*Ht;
RHH = zeros(N,D*(D+1)/2);
for i = 1:D
    if i == D
        RHH(:,end) = R*bsxfun(@times,Ht(:,D),Ht(:,D));
    else
        RHH(:,h2ind(i):h2ind(i+1)-1) = R*bsxfun(@times,Ht(:,i:end),Ht(:,i));
    end
end

Gp = zeros(N,D);
for d = 1:D
     Gp(:,d) = sum(RHH(:,idxGp(d,:)).*W,2);
end
cost_eval =X2+sum(sum(W.*Gp))-2*sum(sum(XHt.*W));
cost_diff = inf;

while (k < maxiter && cost_diff >= tol*cost_eval) % Updating W
    cost_eval_old = cost_eval;
    % Multiplicative Update
    if initialize || mu
        w_up = (XHt+scale*eps)./(Gp+scale*eps);
        W = W.*w_up;
        W(W==0 & Gp<XHt)=eps*scale; % inadmissible zeros and 'stuck'-fix
    end
    
    % HALS update
    if ~initialize && hals
        for d = randperm(D) % columnwise update            
            W(:,d) = max(0,(XHt(:,d)-(sum(RHH(:,IND(d,:)).*W(:,idx_not_d(d,:)),2)))./RHH(:,h2ind(d)));
        end
    end
    
    for d = 1:D
        Gp(:,d) = sum(RHH(:,idxGp(d,:)).*W,2); % Gp recalculated for cost
    end    
    costpart = sum(sum(W.*Gp))-2*sum(sum(XHt.*W));
    cost_eval = X2+costpart;
    cost_diff = cost_eval_old - cost_eval; % convergence check
    
    k = k+1;
end
SSE=cost_eval;

%eof
end