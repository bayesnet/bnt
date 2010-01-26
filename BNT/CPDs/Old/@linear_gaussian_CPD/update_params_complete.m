function CPD = update_params_complete(CPD, self_ev, pev)
% UPDATE_PARAMS_COMPLETE Bayesian parameter updating given completely observed data (linear_gaussian)
% CPD = update_params_complete(CPD, self_ev, pev)
%
% self_ev{m} is the evidence on this node in case m.
% pev{i,m} is the evidence on the i'th parent in case m
%
% We update the hyperparams and set the params to the mean of the posterior.

y = cat(1, self_ev{:});
X = cell2num(pev)';
[N k] = size(X); % each row is a case

n0 = CPD.prior.n;
th0 = CPD.prior.theta;
CPD.prior.theta = inv(n0 + X'*X)*(n0*th0 + X'*y);
thn = CPD.prior.theta;
CPD.prior.beta = CPD.prior.beta + 0.5*(y-X*thn)'*y + 0.5*(th0-thn)'*n0*th0;
CPD.prior.alpha = CPD.prior.alpha + 0.5*N;
CPD.prior.n = CPD.prior.n + X'*X;

  
% set params to their mean
CPD.theta = CPD.prior.theta;
%CPD.sigma = CPD.prior.beta/CPD.prior.alpha; % mean of Gamma is E[lambda] = alpha/beta 
