function CPD = gaussian_CPD(varargin)
% GAUSSIAN_CPD Make a conditional linear Gaussian distrib.
%
% To define this CPD precisely, call the continuous (cts) parents (if any) X,
% the discrete parents (if any) Q, and this node Y. Then the distribution on Y is:
% - no parents: Y ~ N(mu, Sigma)
% - cts parents : Y|X=x ~ N(mu + W x, Sigma)
% - discrete parents: Y|Q=i ~ N(mu(i), Sigma(i))
% - cts and discrete parents: Y|X=x,Q=i ~ N(mu(i) + W(i) x, Sigma(i))
%
% CPD = gaussian_CPD(bnet, node, ...) will create a CPD with random parameters,
% where node is the number of a node in this equivalence class.
%
% The list below gives optional arguments [default value in brackets].
% (Let ns(i) be the size of node i, X = ns(X), Y = ns(Y) and Q = prod(ns(Q)).)
%
% mean       - mu(:,i) is the mean given Q=i [ randn(Y,Q) ]
% cov        - Sigma(:,:,i) is the covariance given Q=i [ repmat(eye(Y,Y), [1 1 Q]) ]
% weights    - W(:,:,i) is the regression matrix given Q=i [ randn(Y,X,Q) ]
% cov_type   - if 'diag', Sigma(:,:,i) is diagonal [ 'full' ]
% tied_cov   - if 1, we constrain Sigma(:,:,i) to be the same for all i [0]
% clamp_mean - if 1, we do not adjust mu(:,i) during learning [0]
% clamp_cov  - if 1, we do not adjust Sigma(:,:,i) during learning [0]
% clamp_weights - if 1, we do not adjust W(:,:,i) during learning [0]
% cov_prior_weight - weight given to I prior for estimating Sigma [0.01]
%
% e.g., CPD = gaussian_CPD(bnet, i, 'mean', [0; 0], 'clamp_mean', 'yes')
%
% For backwards compatibility with BNT2, you can also specify the parameters in the following order
%   CPD = gaussian_CPD(bnet, self, mu, Sigma, W, cov_type, tied_cov, clamp_mean, clamp_cov, clamp_weight)
%
% Sometimes it is useful to create an "isolated" CPD, without needing to pass in a bnet.
% In this case, you must specify the discrete and cts parents (dps, cps) and the family sizes, followed
% by the optional arguments above:
%   CPD = gaussian_CPD('self', i, 'dps', dps, 'cps', cps, 'sz', fam_size, ...)


if nargin==0
  % This occurs if we are trying to load an object from a file.
  CPD = init_fields;
  clamp = 0;
  CPD = class(CPD, 'gaussian_CPD', generic_CPD(clamp));
  return;
elseif isa(varargin{1}, 'gaussian_CPD')
  % This might occur if we are copying an object.
  CPD = varargin{1};
  return;
end
CPD = init_fields;
 
CPD = class(CPD, 'gaussian_CPD', generic_CPD(0));


% parse mandatory arguments
if ~isstr(varargin{1}) % pass in bnet
  bnet = varargin{1};
  self = varargin{2};
  args = varargin(3:end);
  ns = bnet.node_sizes;
  ps = parents(bnet.dag, self);
  dps = myintersect(ps, bnet.dnodes);
  cps = myintersect(ps, bnet.cnodes);
  fam_sz = ns([ps self]);
else
  disp('parsing new style')
  for i=1:2:length(varargin)
    switch varargin{i},
     case 'self', self = varargin{i+1}; 
     case 'dps',  dps = varargin{i+1};
     case 'cps',  cps = varargin{i+1};
     case 'sz',   fam_sz = varargin{i+1};
    end
  end
  ps = myunion(dps, cps);
  args = varargin;
end

CPD.self = self;
CPD.sizes = fam_sz;

% Figure out which (if any) of the parents are discrete, and which cts, and how big they are
% dps = discrete parents, cps = cts parents
CPD.cps = find_equiv_posns(cps, ps); % cts parent index
CPD.dps = find_equiv_posns(dps, ps);
ss = fam_sz(end);
psz = fam_sz(1:end-1);
dpsz = prod(psz(CPD.dps));
cpsz = sum(psz(CPD.cps));

% set default params
CPD.mean = randn(ss, dpsz);
CPD.cov = 100*repmat(eye(ss), [1 1 dpsz]);    
CPD.weights = randn(ss, cpsz, dpsz);
CPD.cov_type = 'full';
CPD.tied_cov = 0;
CPD.clamped_mean = 0;
CPD.clamped_cov = 0;
CPD.clamped_weights = 0;
CPD.cov_prior_weight = 0.01;

nargs = length(args);
if nargs > 0
  if ~isstr(args{1})
    % gaussian_CPD(bnet, self, mu, Sigma, W, cov_type, tied_cov, clamp_mean, clamp_cov, clamp_weights)
    if nargs >= 1 & ~isempty(args{1}), CPD.mean = args{1}; end
    if nargs >= 2 & ~isempty(args{2}), CPD.cov = args{2}; end
    if nargs >= 3 & ~isempty(args{3}), CPD.weights = args{3}; end
    if nargs >= 4 & ~isempty(args{4}), CPD.cov_type = args{4}; end
    if nargs >= 5 & ~isempty(args{5}) & strcmp(args{5}, 'tied'), CPD.tied_cov = 1; end
    if nargs >= 6 & ~isempty(args{6}), CPD.clamped_mean = 1; end
    if nargs >= 7 & ~isempty(args{7}), CPD.clamped_cov = 1; end
    if nargs >= 8 & ~isempty(args{8}), CPD.clamped_weights = 1; end
  else
    CPD = set_fields(CPD, args{:});
  end
end

% Make sure the matrices have 1 dimension per discrete parent.
% Bug fix due to Xuejing Sun 3/6/01
CPD.mean = myreshape(CPD.mean, [ss ns(dps)]);
CPD.cov = myreshape(CPD.cov, [ss ss ns(dps)]);
CPD.weights = myreshape(CPD.weights, [ss cpsz ns(dps)]);
  
CPD.init_cov = CPD.cov;  % we reset to this if things go wrong during learning

% expected sufficient statistics 
CPD.Wsum = zeros(dpsz,1);
CPD.WYsum = zeros(ss, dpsz);
CPD.WXsum = zeros(cpsz, dpsz);
CPD.WYYsum = zeros(ss, ss, dpsz);
CPD.WXXsum = zeros(cpsz, cpsz, dpsz);
CPD.WXYsum = zeros(cpsz, ss, dpsz);

% For BIC
CPD.nsamples = 0;
switch CPD.cov_type
  case 'full',
    ncov_params = ss*(ss-1)/2; % since symmetric (and positive definite)
  case 'diag',
    ncov_params = ss;
  otherwise
    error(['unrecognized cov_type ' cov_type]);
end
% params = weights + mean + cov
if CPD.tied_cov
  CPD.nparams = ss*cpsz*dpsz + ss*dpsz + ncov_params;
else
  CPD.nparams = ss*cpsz*dpsz + ss*dpsz + dpsz*ncov_params;
end



clamped = CPD.clamped_mean & CPD.clamped_cov & CPD.clamped_weights;
CPD = set_clamped(CPD, clamped);

%%%%%%%%%%%

function CPD = init_fields()
% This ensures we define the fields in the same order 
% no matter whether we load an object from a file,
% or create it from scratch. (Matlab requires this.)

CPD.self = [];
CPD.sizes = [];
CPD.cps = [];
CPD.dps = [];
CPD.mean = [];
CPD.cov = [];
CPD.weights = [];
CPD.clamped_mean = [];
CPD.clamped_cov = [];
CPD.clamped_weights = [];
CPD.init_cov = [];
CPD.cov_type = [];
CPD.tied_cov = [];
CPD.Wsum = [];
CPD.WYsum = [];
CPD.WXsum = [];
CPD.WYYsum = [];
CPD.WXXsum = [];
CPD.WXYsum = [];
CPD.nsamples = [];
CPD.nparams = [];            
CPD.cov_prior_weight = [];
