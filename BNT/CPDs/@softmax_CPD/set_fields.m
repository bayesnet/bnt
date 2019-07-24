function CPD = set_params(CPD, varargin)
% SET_PARAMS Set the parameters (fields) for a softmax_CPD object
% CPD = set_params(CPD, name/value pairs)
%
% The following optional arguments can be specified in the form of name/value pairs:
% (Let ns(i) be the size of node i, X = ns(X), Y = ns(Y), Q1=ns(dps(1)), Q2=ns(dps(2)), ...
%   where dps are the discrete parents; if there are no discrete parents, we set Q1=1.)
%
% weights - (W(:,j,a,b,...) - W(:,j',a,b,...)) is ppn to dec. boundary
%           between j,j' given Q1=a,Q2=b,... [ randn(X,Y,Q1,Q2,...) ]
% offset  - (offset(j,a,b,...) - offset(j',a,b,...)) is the offset to dec. boundary
%           between j,j' given Q1=a,Q2=b,... [ randn(Y,Q1,Q2,...) ]
% clamped     - 'yes' means don't adjust params during learning ['no']
% max_iter    - the maximum number of steps to take [10]
% verbose     - 'yes' means print the LL at each step of IRLS ['no']
% wthresh     - convergence threshold for weights [1e-2]
% llthresh    - convergence threshold for log likelihood [1e-2]
% approx_hess - 'yes' means approximate the Hessian for speed ['no']
%
% e.g., CPD = set_params(CPD,'offset', zeros(ns(i),1));

args = varargin;
nargs = length(args);
glimsz = prod(CPD.sizes(CPD.dpndx));
for i=1:2:nargs
  switch args{i},
   case 'discrete',     str='nothing to do';   
   case 'clamped',      CPD = set_clamped(CPD, strcmp(args{i+1}, 'yes'));
   case 'max_iter',     CPD.max_iter = args{i+1};
   case 'verbose',      CPD.verbose = strcmp(args{i+1}, 'yes');
   case 'max_iter',     CPD.max_iter = args{i+1};
   case 'wthresh',      CPD.wthresh = args{i+1};
   case 'llthresh',     CPD.llthresh = args{i+1};
   case 'approx_hess',  CPD.approx_hess = strcmp(args{i+1}, 'yes');
   case 'weights',      for q=1:glimsz, CPD.glim{q}.w1 = args{i+1}(:,:,q); end; 
   case 'offset',
    if glimsz == 1
      CPD.glim{1}.b1 = args{i+1};
    else
      for q=1:glimsz, CPD.glim{q}.b1 = args{i+1}(:,q); end; 
    end
   otherwise,  
    error(['invalid argument name ' args{i}]);       
  end
end
