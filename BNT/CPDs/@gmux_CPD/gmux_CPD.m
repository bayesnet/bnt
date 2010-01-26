function CPD = gmux_CPD(bnet, self, varargin)
% GMUX_CPD Make a Gaussian multiplexer node
%
% CPD = gmux_CPD(bnet, node, ...) is used similarly to gaussian_CPD,
% except we assume there is exactly one discrete parent (call it M)
% which is used to select which cts parent to pass through to the output.
% i.e., we define P(Y=y|M=m, X1, ..., XK) = N(y | W(m)*x(m) + mu(m), Sigma(m))
% where Y represents this node, and the Xi's are the cts parents.
% All the Xi must have the same size, and the num values for M must be K.
%
% Currently the params for this kind of CPD cannot be learned.
%
% Optional arguments [ default in brackets ]
%
% mean       - mu(:,i) is the mean given M=i [ zeros(Y,K) ]
% cov        - Sigma(:,:,i) is the covariance given M=i [ repmat(1*eye(Y,Y), [1 1 K]) ]
% weights    - W(:,:,i) is the regression matrix given M=i [ randn(Y,X,K) ]

if nargin==0
  % This occurs if we are trying to load an object from a file.
  CPD = init_fields;
  clamp = 0;
  CPD = class(CPD, 'gmux_CPD', generic_CPD(clamp));
  return;
elseif isa(bnet, 'gmux_CPD')
  % This might occur if we are copying an object.
  CPD = bnet;
  return;
end
CPD = init_fields;
 
CPD = class(CPD, 'gmux_CPD', generic_CPD(1));

ns = bnet.node_sizes;
ps = parents(bnet.dag, self);
dps = myintersect(ps, bnet.dnodes);
cps = myintersect(ps, bnet.cnodes);
fam_sz = ns([ps self]);

CPD.self = self;
CPD.sizes = fam_sz;

% Figure out which (if any) of the parents are discrete, and which cts, and how big they are
% dps = discrete parents, cps = cts parents
CPD.cps = find_equiv_posns(cps, ps); % cts parent index
CPD.dps = find_equiv_posns(dps, ps);
if length(CPD.dps) ~= 1
  error('gmux must have exactly 1 discrete parent')
end
ss = fam_sz(end);
cpsz = fam_sz(CPD.cps(1)); % in gaussian_CPD, cpsz = sum(fam_sz(CPD.cps))
if ~all(fam_sz(CPD.cps) == cpsz)
  error('all cts parents must have same size')
end
dpsz = fam_sz(CPD.dps);
if dpsz ~= length(cps)
  error(['the arity of the mux node is ' num2str(dpsz) ...
	 ' but there are ' num2str(length(cps)) ' cts parents']);
end

% set default params
%CPD.mean = zeros(ss, 1);
%CPD.cov = eye(ss);
%CPD.weights = randn(ss, cpsz);
CPD.mean = zeros(ss, dpsz);
CPD.cov = 1*repmat(eye(ss), [1 1 dpsz]);    
CPD.weights = randn(ss, cpsz, dpsz);

args = varargin;
nargs = length(args);
for i=1:2:nargs
  switch args{i},
   case 'mean',        CPD.mean = args{i+1}; 
   case 'cov',         CPD.cov = args{i+1}; 
   case 'weights',    CPD.weights = args{i+1}; 
   otherwise,  
    error(['invalid argument name ' args{i}]);
  end
end

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

