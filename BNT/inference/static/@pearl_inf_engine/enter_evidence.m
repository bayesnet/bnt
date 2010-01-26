function [engine, loglik, niter] = enter_evidence(engine, evidence, varargin)
% ENTER_EVIDENCE Add the specified evidence to the network (pearl)
% [engine, loglik, num_iter] = enter_evidence(engine, evidence, ...)
% evidence{i} = [] if if X(i) is hidden, and otherwise contains its observed value (scalar or column vector)
%
% The following optional arguments can be specified in the form of name/value pa irs:
% [default value in brackets]
%
% maximize - if 1, does max-product instead of sum-product [0]
% 'filename' -  msgs will be printed to this file, so you can assess convergence while it runs [engine.filename]
%
% e.g., engine = enter_evidence(engine, ev, 'maximize', 1)
%     
% For discrete nodes, loglik is the negative Bethe free energy evaluated at the final beliefs.
% For Gaussian nodes, loglik is currently always 0.
%
% 'num_iter' returns the number of iterations used.

maximize = 0;
filename = engine.filename;

% parse optional params
args = varargin;
nargs = length(args);
if nargs > 0
  for i=1:2:nargs
    switch args{i},
     case 'maximize', maximize = args{i+1};
     case 'filename', filename = args{i+1};
     otherwise,
      error(['invalid argument name ' args{i}]);
    end
  end
end
    

if maximize
  error('can''t handle max-prop yet')
end

engine.maximize = maximize;
engine.filename = filename;
engine.bel = []; % reset if necessary

bnet = bnet_from_engine(engine);
N = length(bnet.dag);
ns = bnet.node_sizes(:);

observed_bitv = ~isemptycell(evidence);
disconnected = find(engine.disconnected_nodes_bitv);
if ~all(observed_bitv(disconnected))
  error(['The following discrete nodes must be observed: ' num2str(disconnected)])
end
msg = init_pearl_msgs(engine.msg_type, engine.msg_dag, ns, evidence);

niter = 1;
switch engine.protocol
 case 'parallel', [msg, niter] = parallel_protocol(engine, evidence, msg);
 case 'tree', msg = tree_protocol(engine, evidence, msg);
 otherwise,
  error(['unrecognized protocol ' engine.protocol])
end
engine.niter = niter;

engine.marginal = cell(1,N);
nodes = find(~engine.disconnected_nodes_bitv);
for n=nodes(:)'
  engine.marginal{n} = compute_bel(engine.msg_type, msg{n}.pi, msg{n}.lambda);
end

engine.evidence = evidence; % needed by marginal_nodes and marginal_family
engine.msg = msg;  % needed by marginal_family

if (nargout >= 2)
  if (engine.msg_type == 'd')
    loglik = bethe_free_energy(engine, evidence);
  else
    loglik = 0;
  end
end



%%%%%%%%%%%

function msg =  init_pearl_msgs(msg_type, dag, ns, evidence)
% INIT_MSGS Initialize the lambda/pi message and state vectors
% msg =  init_msgs(dag, ns, evidence)
%

N = length(dag);
msg = cell(1,N);
observed = ~isemptycell(evidence);
lam_msg = 1;

for n=1:N
  ps = parents(dag, n);
  msg{n}.pi_from_parent = cell(1, length(ps));
  for i=1:length(ps)
    p = ps(i);
    msg{n}.pi_from_parent{i} = mk_msg(msg_type, ns(p));
  end
  
  cs = children(dag, n);
  msg{n}.lambda_from_child = cell(1, length(cs));
  for i=1:length(cs)
    c = cs(i);
    msg{n}.lambda_from_child{i} = mk_msg(msg_type, ns(n), lam_msg);
  end

  msg{n}.lambda = mk_msg(msg_type, ns(n), lam_msg);
  msg{n}.pi = mk_msg(msg_type, ns(n));
  
  if observed(n)
    msg{n}.lambda_from_self = mk_msg_with_evidence(msg_type, ns(n), evidence{n});
  else
    msg{n}.lambda_from_self = mk_msg(msg_type, ns(n), lam_msg);
  end
end



%%%%%%%%%

function msg =  mk_msg(msg_type, sz, is_lambda_msg)

if nargin < 3, is_lambda_msg = 0; end

switch msg_type
 case 'd', msg = ones(sz, 1);
 case 'g', 
  if is_lambda_msg
    msg.precision = zeros(sz, sz);
    msg.info_state = zeros(sz, 1);
  else
    msg.Sigma = zeros(sz, sz);
    msg.mu = zeros(sz,1);
  end
end

%%%%%%%%%%%%

function msg = mk_msg_with_evidence(msg_type, sz, val)

switch msg_type
 case 'd',
  msg = zeros(sz, 1);
  msg(val) = 1;
 case 'g',
  %msg.observed_val = val(:);
  msg.precision = inf;
  msg.mu = val(:);
end
