function CPD = mlp_CPD(bnet, self, nhidden, w1, b1, w2, b2, clamped, max_iter, verbose, wthresh,  llthresh)
% MLP_CPD Make a CPD from a Multi Layer Perceptron (i.e., feedforward neural network)
%
% We use a different MLP for each discrete parent combination (if there are any discrete parents).
% We currently assume this node (the child) is discrete.
%
% CPD = mlp_CPD(bnet, self, nhidden)
% will create a CPD with random parameters, where self is the number of this node and nhidden the number of the hidden nodes.
% The params are drawn from N(0, s*I), where s = 1/sqrt(n+1), n = length(X).
%
% CPD = mlp_CPD(bnet, self, nhidden, w1, b1, w2, b2) allows you to specify the params, where
%  w1 = first-layer weight matrix
%  b1 = first-layer bias vector
%  w2 = second-layer weight matrix
%  b2 = second-layer bias vector
% These are assumed to be the same for each discrete parent combination.
% If any of these are [], random values will be created.
%
% CPD = mlp_CPD(bnet, self, nhidden, w1, b1, w2, b2, clamped) allows you to prevent the params from being
% updated during learning (if clamped = 1). Default: clamped = 0.
%
% CPD = mlp_CPD(bnet, self, nhidden, w1, b1, w2, b2, clamped, max_iter, verbose, wthresh,  llthresh)
% alllows you to specify params that control the M step:
%  max_iter - the maximum number of steps to take (default: 10)
%  verbose - controls whether to print (default: 0 means silent).
%  wthresh - a measure of the precision required for the value of
%     the weights W at the solution. Default: 1e-2.
%  llthresh - a measure of the precision required of the objective
%     function (log-likelihood) at the solution.  Both this and the previous condition must
%     be satisfied for termination. Default: 1e-2.
%
% For learning, we use a weighted version of scaled conjugated gradient in the M step.

if nargin==0
  % This occurs if we are trying to load an object from a file.
  CPD = init_fields;
  CPD = class(CPD, 'mlp_CPD', discrete_CPD(0,[]));
  return;
elseif isa(bnet, 'mlp_CPD')
  % This might occur if we are copying an object.
  CPD = bnet;
  return;
end
CPD = init_fields;

assert(myismember(self, bnet.dnodes));
ns = bnet.node_sizes;

ps = parents(bnet.dag, self);
dnodes = mysetdiff(1:length(bnet.dag), bnet.cnodes);
dps = myintersect(ps, dnodes);
cps = myintersect(ps, bnet.cnodes);
dpsz = prod(ns(dps));
cpsz = sum(ns(cps));
self_size = ns(self);

% discrete/cts parent index - which ones of my parents are discrete/cts?
CPD.dpndx = find_equiv_posns(dps, ps); 
CPD.cpndx = find_equiv_posns(cps, ps);

CPD.mlp = cell(1,dpsz);
for i=1:dpsz
    CPD.mlp{i} = mlp(cpsz, nhidden, self_size, 'softmax');
    if nargin >=4 & ~isempty(w1)
        CPD.mlp{i}.w1 = w1;
    end
    if nargin >=5 & ~isempty(b1)
        CPD.mlp{i}.b1 = b1; 
    end
    if nargin >=6 & ~isempty(w2)
        CPD.mlp{i}.w2 = w2; 
    end
    if nargin >=7 & ~isempty(b2)
        CPD.mlp{i}.b2 = b2; 
    end
    W1app(:,:,i)=CPD.mlp{i}.w1;
    W2app(:,:,i)=CPD.mlp{i}.w2;
    b1app(i,:)=CPD.mlp{i}.b1;
    b2app(i,:)=CPD.mlp{i}.b2;
end
if nargin < 8, clamped = 0; end
if nargin < 9, max_iter = 10; end
if nargin < 10, verbose = 0; end
if nargin < 11, wthresh = 1e-2; end
if nargin < 12, llthresh = 1e-2; end

CPD.self = self;
CPD.max_iter = max_iter;
CPD.verbose = verbose;
CPD.wthresh = wthresh;
CPD.llthresh = llthresh;

% sufficient statistics 
% Since MLP is not in the exponential family, we must store all the raw data.
%
CPD.W1=W1app;                     % Extract all the parameters of the node for handling discrete obs parents
CPD.W2=W2app;                     %
nparaW=[size(W1app) size(W2app)]; %
CPD.b1=b1app;                     %
CPD.b2=b2app;                     %
nparab=[size(b1app) size(b2app)]; %

CPD.sizes=bnet.node_sizes(:);   % used in CPD_to_table to pump up the node sizes

CPD.parent_vals = [];        % X(l,:) = value of cts parents in l'th example

CPD.eso_weights=[];          % weights used by the SCG algorithm 

CPD.self_vals = [];          % Y(l,:) = value of self in l'th example

% For BIC
CPD.nsamples = 0;   
CPD.nparams=prod(nparaW)+prod(nparab);
CPD = class(CPD, 'mlp_CPD', discrete_CPD(clamped, ns([ps self])));

%%%%%%%%%%%

function CPD = init_fields()
% This ensures we define the fields in the same order 
% no matter whether we load an object from a file,
% or create it from scratch. (Matlab requires this.)

CPD.mlp = {};
CPD.self = [];
CPD.max_iter = [];
CPD.verbose = [];
CPD.wthresh = [];
CPD.llthresh = [];
CPD.approx_hess = [];
CPD.W1 = [];
CPD.W2 = [];
CPD.b1 = [];
CPD.b2 = [];
CPD.sizes = [];
CPD.parent_vals = [];
CPD.eso_weights=[];
CPD.self_vals = [];
CPD.nsamples = [];
CPD.nparams = [];
