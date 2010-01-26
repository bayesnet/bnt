function CPD = linear_gaussian_CPD(bnet, self, theta, sigma, theta0, n0, alpha0, beta0)
% LINEAR_GAUSSIAN_CPD Make a linear Gaussian distrib.
%
% CPD = linear_gaussian_CPD(bnet, self, theta, lambda)
% This defines the distribution P(Y|X) =  N(y | theta'*x, sigma),
% where y (self) is a scalar, theta is a regression vector, and sigma is the variance.
% Pass in [] to generate a default random value for a parameter.
%
% CPD = linear_gaussian_CPD(bnet, self, [], [], theta0, n0, alpha0, beta0)
% defines a Normal-Gamma prior over the parameters:
%   P(theta | lambda) = N(theta | theta0, n0*lambda)
%   P(lambda) = Gamma(lambda | alpha0, beta0)
% where lambda = 1/sigma is the precision for y.
% n0 is a precision matrix, beta0 is a scale factor.
% Pass in [] to generate a default value for a hyperparameter.
% theta and sigma will be set to their prior expected values.
% See "Bayesian Theory", Bernardo and Smith (2000), p442.


if nargin==0
  % This occurs if we are trying to load an object from a file.
  CPD = init_fields;
  CPD = class(CPD, 'linear_gaussian_CPD', generic_CPD(0));
  return;
elseif isa(bnet, 'linear_gaussian_CPD')
  % This might occur if we are copying an object.
  CPD = bnet;
  return;
end
CPD = init_fields;


ns = bnet.node_sizes;
ps = parents(bnet.dag, self);
d = sum(ns(ps));
assert(ns(self)==1);


if nargin < 5,
  prior = [];
  if isempty(theta), theta = randn(d, 1); end
  if isempty(sigma), sigma = 1; end
else
  
  %if isempty(theta0), theta0 = zeros(d, 1); end
  %if isempty(n0), n0 = 0.1*eye(d); end
  %if isempty(alpha0), alpha0 = 0.1; end
  %if isempty(beta0), beta0 = 0.1; end
   
  % use non-informative priors
  if isempty(theta0), theta0 = zeros(d, 1); end
  if isempty(n0), n0 = 0.001*ones(d); end
  if isempty(alpha0), alpha0 = -d/2 + 0.001; end
  if isempty(beta0), beta0 = 0.001; end

  prior.theta = theta0;
  prior.n = n0;
  prior.alpha = alpha0;
  prior.beta = beta0;
  
  % set params to their mean
  theta = prior.theta;
  %sigma = prior.beta/prior.alpha; % mean of Gamma is E[lambda] = alpha/beta 
end


CPD.self = self;
CPD.theta = theta;
CPD.sigma = sigma;
CPD.prior = prior;


clamped = 0;
CPD = class(CPD, 'linear_gaussian_CPD', generic_CPD(clamped));


%%%%%%%%%%%

function CPD = init_fields()
% This ensures we define the fields in the same order 
% no matter whether we load an object from a file,
% or create it from scratch. (Matlab requires this.)

CPD.self = [];
CPD.theta = [];
CPD.sigma = [];
CPD.prior = [];
