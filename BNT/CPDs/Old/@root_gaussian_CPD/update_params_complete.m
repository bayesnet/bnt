function CPD = update_params_complete(CPD, self_ev, pev)
% UPDATE_PARAMS_COMPLETE Bayesian parameter updating given completely observed data (root_gaussian)
% CPD = update_params_complete(CPD, self_ev, pev)
%
% self_ev{m} is the evidence on this node in case m.
% pev{i,m} is the evidence on the i'th parent in case m (ignored)
%
% We update the hyperparams and set the params to the mean of the posterior.

X = cell2num(self_ev);
[k N] = size(X); % each column is a case

one = ones(N,1);
xbar = X*one / N; % = mean(X')'
S = X*(eye(N) - one*one'/N)*X';

n0 = CPD.prior.n;
nn = 1/(n0 + N);
mu0 = CPD.prior.mu;
CPD.prior.mu = nn*(n0*mu0 + N*xbar);
CPD.prior.alpha = CPD.prior.alpha + 0.5*N;
CPD.prior.beta = CPD.prior.beta + 0.5*S + 0.5*nn*N*n0*(mu0-xbar)*(mu0-xbar)';
CPD.prior.n = CPD.prior.n + N;

% set params to their mean
CPD.mu = CPD.prior.mu;
% E[Cov] = E inv(n lambda) = 1/(n (alpha-(k+1)/2)) beta
CPD.Sigma = CPD.prior.beta /(CPD.prior.n * (CPD.prior.alpha - (k+1)/2));
  
