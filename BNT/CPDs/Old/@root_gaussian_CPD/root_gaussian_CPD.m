function CPD = root_gaussian_CPD(bnet, self, mu, Sigma, mu0, n0, alpha0, beta0)
% ROOT_GAUSSIAN_CPD Make an unconditional Gaussian distrib.
%
% CPD = root_gaussian_CPD(bnet, self, mu, Sigma)
% This defines the distribution Y ~ N(mu, Sigma),
% Pass in [] to generate a default random value for a parameter.
%
% CPD = root_gaussian_CPD(bnet, self, [], [], mu0, n0, alpha0, beta0)
% defines a Normal-Wishart prior over the parameters:
%   P(mu | lambda) = N(mu | mu0, n0*lambda)
%   P(lambda) = Wishart(lambda | alpha0, beta0)
% where lambda = inv(Sigma) is the precision matrix of mu.
% n0 is a scale factor, beta0 is a precision matrix.
% Pass in [] to generate a default value for a hyperparameter.
% mu and Sigma will be set to their prior expected values.
% See "Bayesian Theory", Bernardo and Smith (2000), p441.


if nargin==0
  % This occurs if we are trying to load an object from a file.
  CPD = init_fields;
  CPD = class(CPD, 'root_gaussian_CPD', generic_CPD(0));
  return;
elseif isa(bnet, 'root_gaussian_CPD')
  % This might occur if we are copying an object.
  CPD = bnet;
  return;
end
CPD = init_fields;


ns = bnet.node_sizes;
d = ns(self);

if nargin < 5,
  prior = [];
  if isempty(mu), mu = randn(d, 1); end
  if isempty(Sigma), Sigma = eye(d); end
else
  if isempty(mu0), mu0 = zeros(d, 1); end
  if isempty(n0), n0 = 0.1; end
  if isempty(alpha0), alpha0 = (d-1)/2 + 1; end % Wishart requires 2 alpha > d-1
  if isempty(beta0), beta0 = eye(d); end
  
  prior.mu = mu0;
  prior.n = n0;
  prior.alpha = alpha0;
  prior.beta = beta0;
  
  % set params to their mean
  mu = prior.mu;
  Sigma = prior.beta/prior.alpha; % mean of Wishart is E[lambda] = alpha*inv(beta)
end

CPD.self = self;
CPD.mu = mu;
CPD.Sigma = Sigma;
CPD.prior = prior;

clamped = 0;
CPD = class(CPD, 'root_gaussian_CPD', generic_CPD(clamped));


%%%%%%%%%%%

function CPD = init_fields()
% This ensures we define the fields in the same order 
% no matter whether we load an object from a file,
% or create it from scratch. (Matlab requires this.)

CPD.self = [];
CPD.mu = [];
CPD.Sigma = [];
CPD.prior = [];
