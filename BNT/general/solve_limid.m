function [strategy, MEU, niter] = solve_limid(engine, varargin)
% SOLVE_LIMID Find the (locally) optimal strategy for a LIMID
% [strategy, MEU, niter] = solve_limid(inf_engine, ...)
%
% strategy{d} = stochastic policy for node d (a decision node)
% MEU = maximum expected utility
% niter = num iterations used
%
% The following optional arguments can be specified in the form of name/value pairs:
% [default in brackets]
%
% max_iter - max. num. iterations [ 1 ]
% tol - tolerance required of consecutive MEU values, used to assess convergence [1e-3]
% order - order in which decision nodes are optimized [ reverse numerical order ]
%
% e.g., solve_limid(engine, 'tol', 1e-2, 'max_iter', 10)

bnet = bnet_from_engine(engine);

% default values
max_iter = 1;
tol = 1e-3;
D = bnet.decision_nodes;
order = D(end:-1:1);

args = varargin;
nargs = length(args);
for i=1:2:nargs
  switch args{i},
   case 'max_iter', max_iter  = args{i+1}; 
   case 'tol',      tol = args{i+1}; 
   case 'order',    order = args{i+1}; 
   otherwise,  
    error(['invalid argument name ' args{i}]);       
  end
end

CPDs = bnet.CPD;
ns = bnet.node_sizes;
N = length(ns);
evidence = cell(1,N);
strategy = cell(1, N);

iter = 1;
converged = 0;
oldMEU = 0;
while ~converged & (iter <= max_iter)
  for d=order(:)'
    engine = enter_evidence(engine, evidence, 'exclude', d);
    [m, pot] = marginal_family(engine, d);
    %pot = marginal_family_pot(engine, d);
    [policy, score] = upot_to_opt_policy(pot);    
    e = bnet.equiv_class(d);
    CPDs{e} = set_fields(CPDs{e}, 'policy', policy);
    engine = update_engine(engine, CPDs);
    strategy{d} = policy;
  end  
  engine = enter_evidence(engine, evidence);
  [m, pot] = marginal_nodes(engine, []);
  %pot = marginal_family_pot(engine, []);
  [dummy, MEU] = upot_to_opt_policy(pot);    
  if approxeq(MEU, oldMEU, tol)
    converged = 1;
  end
  oldMEU = MEU;
  iter = iter + 1;
end
niter = iter - 1;
