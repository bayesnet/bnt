function [net, gamma, logev] = evidence_weighted(net, x, t, eso_w, num)
%EVIDENCE Re-estimate hyperparameters using evidence approximation.
%
%	Description
%	[NET] = EVIDENCE(NET, X, T) re-estimates the hyperparameters ALPHA
%	and BETA by applying Bayesian re-estimation formulae for NUM
%	iterations. The hyperparameter ALPHA can be a simple scalar
%	associated with an isotropic prior on the weights, or can be a vector
%	in which each component is associated with a group of weights as
%	defined by the INDEX matrix in the NET data structure. These more
%	complex priors can be set up for an MLP using MLPPRIOR. Initial
%	values for the iterative re-estimation are taken from the network
%	data structure NET passed as an input argument, while the return
%	argument NET contains the re-estimated values.
%
%	[NET, GAMMA, LOGEV] = EVIDENCE(NET, X, T, NUM) allows the re-
%	estimation  formula to be applied for NUM cycles in which the re-
%	estimated values for the hyperparameters from each cycle are used to
%	re-evaluate the Hessian matrix for the next cycle.  The return value
%	GAMMA is the number of well-determined parameters and LOGEV is the
%	log of the evidence.
%
%	See also
%	MLPPRIOR, NETGRAD, NETHESS, DEMEV1, DEMARD
%

%	Copyright (c) Ian T Nabney (1996-9)

errstring = consist(net, '', x, t);
if ~isempty(errstring)
  error(errstring);
end

ndata = size(x, 1);
if nargin == 4
  num = 1;
end

if isfield(net,'beta')
    beta = net.beta;
else 
    beta = 1;
end;

% Extract weights from network
pakstr = [net.type, 'pak'];
w = feval(pakstr, net);

% Evaluate data-dependent contribution to the Hessian matrix.
[h, dh] = nethess_weighted(w, net, x, t, eso_w); 

% Now set the negative eigenvalues to zero.
[evec, evl] = eig(dh);
evl = evl.*(evl > 0);
% safe_evl is used to avoid taking log of zero
safe_evl = evl + eps.*(evl <= 0);

% Do the re-estimation. 
for k = 1 : num
  [e, edata, eprior] = neterr_weighted(w, net, x, t, eso_w);
  h = nethess_weighted(w, net, x, t, eso_w, dh);
  % Re-estimate alpha.
  if size(net.alpha) == [1 1]
    % Evaluate number of well-determined parameters.
    if k == 1
      % Form vector of eigenvalues
      evl = diag(evl);
      safe_evl = diag(safe_evl);
    end
    B = beta*evl;
    gamma = sum(B./(B + net.alpha));       
    net.alpha = 0.5*gamma/eprior;
       
    % Partially evaluate log evidence
    logev = e - 0.5*sum(log(safe_evl)) + 0.5*net.nwts*log(net.alpha) - ...
      0.5*ndata*log(2*pi);
  else
    ngroups = size(net.alpha, 1);
    gams = zeros(1, ngroups);
    logas = zeros(1, ngroups);
    traces = zeros(1, ngroups);
    % Reconstruct data hessian with negative eigenvalues set to zero.
    dh = evec*evl*evec';
    hinv = inv(nethess_weighted(w, net, x, t, eso_w, dh));
    for m = 1 : ngroups
      group_nweights = sum(net.index(:, m));
      gams(m) = group_nweights - ...
	        net.alpha(m)*sum(diag(hinv).*net.index(:,m));
      net.alpha(m) = real(gams(m)/(2*eprior(m)));
      % Weight alphas by number of weights in group
      logas(m) = 0.5*group_nweights*log(net.alpha(m));
      % Compute sum of evalues corresponding to group
      traces(m) = sum(log(safe_evl*net.index(:,m)));
    end 
    gamma = sum(gams, 2);
    logev = e - 0.5*sum(traces) + sum(logas) - 0.5*ndata*log(2*pi);
  end
  % Re-estimate beta.
  if isfield(net, 'beta')
      net.beta = 0.5*(net.nout*ndata - gamma)/edata;
  end
  logev = logev + 0.5*ndata*log(beta);
end

