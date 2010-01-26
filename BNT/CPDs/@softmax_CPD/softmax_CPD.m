function CPD = softmax_CPD(bnet, self, varargin)
% SOFTMAX_CPD Make a softmax (multinomial logit) CPD
%
% To define this CPD precisely, let W be an (m x n) matrix with W(i,:) = {i-th row of B} 
% => we can define the following vectorial function:
%    
%                                   softmax: R^n |--> R^m  
%                  softmax(z,i-th)=exp(W(i,:)*z)/sum_k(exp(W(k,:)*z))      
%
% (this constructor augments z with a one at the beginning to introduce an offset term (=bias, intercept))                                   
% Now call the continuous (cts) and always observed (obs) parents X,
% the discrete parents (if any) Q, and this node Y then we use the discrete parent(s) just  to index
% the parameter vectors (c.f., conditional Gaussian nodes); that is:
%                 prob(Y=i | X=x, Q=j) = softmax(x,i-th|j)
% where '|j' means that we are using the j-th (m x n) parameters matrix W(:,:,j).
% If there are no discrete parents, this is a regular softmax node.
% If Y is binary, this is a logistic (sigmoid) function.
%
% CPD = softmax_CPD(bnet, node_num, ...) will create a softmax CPD with random parameters,
% where node is the number of a node in this equivalence class.
%
% The following optional arguments can be specified in the form of name/value pairs:
% [default value in brackets]
% (Let ns(i) be the size of node i, X = ns(X), Y = ns(Y), Q1=ns(dps(1)), Q2=ns(dps(2)), ...
% where dps are the discrete parents; if there are no discrete parents, we set Q1=1.)
%
% discrete - the discrete parents that we want to treat like the cts ones [ [] ]. 
%            This can be used to define sigmoid belief network - see below the reference.             
%            For example suppose that Y has one cts parents X and two discrete ones: Q, C1 where:
%            -> Q is binary (1/2) and used just to index the parameters of 'self'
%            -> C1 is ternary (1/2/3) and treated as a cts node <=> its values appear into the linear 
%               part of the softmax function
%            then:
%                     prob(Y|X=x, Q=q, C1=c1)= softmax(W(:,:,q)' * y)
%            where y = [1 | delta(C1,1) delta(C1,2) delta(C1,3) | x(:)']' and delta(Y,a)=indicator(Y=a).
% weights - (w(:,j,a,b,...) - w(:,j',a,b,...)) is ppn to dec. boundary
%           between j,j' given Q1=a,Q2=b,... [ randn(X,Y,Q1,Q2,...) ]
% offset  - (b(j,a,b,...) - b(j',a,b,...)) is the offset to dec. boundary
%           between j,j' given Q1=a,Q2=b,... [ randn(Y,Q1,Q2,...) ]
%
% e.g., CPD = softmax_CPD(bnet, i, 'offset', zeros(ns(i),1));
%
% The following fields control the behavior of the M step, which uses 
% a weighted version of the Iteratively Reweighted Least Squares (WIRLS) if dps_as_cps=[]; or
% a weighted SCG otherwise, as implemented in Netlab, and modified by Pierpaolo Brutti.
%
% clamped     - 'yes' means don't adjust params during learning ['no']
% max_iter    - the maximum number of steps to take [10]
% verbose     - 'yes' means print the LL at each step of IRLS ['no']
% wthresh     - convergence threshold for weights [1e-2]
% llthresh    - convergence threshold for log likelihood [1e-2]
% approx_hess - 'yes' means approximate the Hessian for speed ['no']
%
% For backwards compatibility with BNT2, you can also specify the parameters in the following order
%   softmax_CPD(bnet, self, w, b, clamped, max_iter, verbose, wthresh, llthresh, approx_hess)
%
% REFERENCE
% For details on the sigmoid belief nets, see:
% - Neal (1992). Connectionist learning of belief networks, Artificial Intelligence, 56, 71-113.
% - Saul, Jakkola, Jordan (1996). Mean field theory for sigmoid belief networks, Journal of Artificial Intelligence Reseach (4), pagg. 61-76.
%
% For details on the M step, see:
% - K. Chen, L. Xu, H. Chi (1999). Improved learning algorithms for mixtures of experts in multiclass 
%       classification. Neural Networks 12, pp. 1229-1252.
% - M.I. Jordan, R.A. Jacobs (1994). Hierarchical Mixtures of Experts and the EM algorithm. 
%       Neural Computation 6, pp. 181-214.
% - S.R. Waterhouse, A.J. Robinson (1994). Classification Using Hierarchical Mixtures of Experts. In Proc. IEEE
%       Workshop on Neural Network for Signal Processing IV, pp. 177-186

if nargin==0
  % This occurs if we are trying to load an object from a file.
  CPD = init_fields;
  CPD = class(CPD, 'softmax_CPD', discrete_CPD(0, []));
  return;
elseif isa(bnet, 'softmax_CPD')
  % This might occur if we are copying an object.
  CPD = bnet;
  return;
end
CPD = init_fields;

assert(myismember(self, bnet.dnodes));
ns = bnet.node_sizes;
ps = parents(bnet.dag, self);
dps = myintersect(ps, bnet.dnodes);
cps = myintersect(ps, bnet.cnodes);

clamped = 0;
CPD = class(CPD, 'softmax_CPD', discrete_CPD(clamped, ns([ps self])));

dps_as_cpssz = 0;
dps_as_cps = [];
% determine if any discrete parents are to be treated as cts
if nargin >= 3 & isstr(varargin{1}) % might have passed in 'discrete'
  for i=1:2:length(varargin)
    if strcmp(varargin{i}, 'discrete')
      dps_as_cps = varargin{i+1};
      assert(myismember(dps_as_cps, dps));
      dps = mysetdiff(dps, dps_as_cps);         % put out the dps treated as cts
      CPD.dps_as_cps.ndx = find_equiv_posns(dps_as_cps, ps);
      CPD.dps_as_cps.separator = [0 cumsum(ns(dps_as_cps(1:end-1)))]; % concatenated dps_as_cps dims separators
      dps_as_cpssz = sum(ns(dps_as_cps));
      break;
    end
  end
end
assert(~isempty(union(cps, dps_as_cps)));   % It have to be at least a cts or a dps_as_cps parents
self_size = ns(self); 
cpsz = sum(ns(cps));  
glimsz = prod(ns(dps));
CPD.dpndx = find_equiv_posns(dps, ps);  % it contains only the indeces of the 'pure' dps
CPD.cpndx = find_equiv_posns(cps, ps);

CPD.self  = self;
CPD.solo  = (length(ns)<=2);
CPD.sizes = bnet.node_sizes([ps self]);

% set default params
CPD.max_iter = 10;
CPD.verbose = 0;
CPD.wthresh = 1e-2;
CPD.llthresh = 1e-2;
CPD.approx_hess = 0;
CPD.glim = cell(1,glimsz);
for i=1:glimsz
  CPD.glim{i} = glm(dps_as_cpssz + cpsz, self_size, 'softmax');
end

if nargin >= 3
  args = varargin;
  nargs = length(args);
  if ~isstr(args{1})
    %   softmax_CPD(bnet, self, w, b, clamped, max_iter, verbose, wthresh, llthresh, approx_hess)
    if nargs >= 1 & ~isempty(args{1}), CPD = set_fields(CPD, 'weights', args{1}); end
    if nargs >= 2 & ~isempty(args{2}), CPD = set_fields(CPD, 'offset', args{2});  end
    if nargs >= 3 & ~isempty(args{3}), CPD = set_clamped(CPD, args{3});           end
    if nargs >= 4 & ~isempty(args{4}), CPD.max_iter    = args{4}; end
    if nargs >= 5 & ~isempty(args{5}), CPD.verbose     = args{5}; end
    if nargs >= 6 & ~isempty(args{6}), CPD.wthresh     = args{6}; end
    if nargs >= 7 & ~isempty(args{7}), CPD.llthresh   = args{7}; end
    if nargs >= 8 & ~isempty(args{8}), CPD.approx_hess = args{8}; end
  else
    CPD = set_fields(CPD, args{:});
  end
end

% sufficient statistics 
% Since dsoftmax is not in the exponential family, we must store all the raw data.
CPD.parent_vals = [];         % X(l,:) = value of cts parents in l'th example
CPD.self_vals = [];           % Y(l,:) = value of self in l'th example

CPD.eso_weights=[];           % weights used by the WIRLS algorithm

% For BIC
CPD.nsamples = 0;   
if ~adjustable_CPD(CPD),
   CPD.nparams=0;
else
   [W, b] = extract_params(CPD);
   CPD.nparams= prod(size(W)) + prod(size(b));
end

%%%%%%%%%%%

function CPD = init_fields()
% This ensures we define the fields in the same order 
% no matter whether we load an object from a file,
% or create it from scratch. (Matlab requires this.)

CPD.glim = {};
CPD.self = [];
CPD.solo = [];
CPD.max_iter = [];
CPD.verbose = [];
CPD.wthresh = [];
CPD.llthresh = [];
CPD.approx_hess = [];
CPD.sizes = [];
CPD.parent_vals = [];
CPD.eso_weights=[];
CPD.self_vals = [];
CPD.nsamples = [];
CPD.nparams = [];
CPD.dpndx = [];
CPD.cpndx = [];
CPD.dps_as_cps.ndx = [];
CPD.dps_as_cps.separator = [];
