function [engine, loglik] = enter_evidence(engine, evidence, varargin)
% ENTER_EVIDENCE Add the specified evidence to the network (global_joint)
% [engine, loglik] = enter_evidence(engine, evidence, ...)
%
% evidence{i} = [] if if X(i) is hidden, and otherwise contains its observed value.
%
% Warning: Computing the log likelihood requires marginalizing all the nodes and can be slow.
%
% The list below gives optional arguments [default value in brackets].      
%
% exclude - list of nodes whose potential will not be included in the joint [ [] ]
%
% e.g., engine = enter_evidence(engine, ev, 'exclude', 3)

exclude = [];
maximize = 0;

if nargin >= 3
  args = varargin;
  nargs = length(args);
  for i=1:2:nargs
    switch args{i},
     case 'exclude', exclude = args{i+1};
     case 'maximize', maximize = args{i+1};
     otherwise,
      error(['invalid argument name ' args{i}]);
    end
  end
end
  
assert(~maximize)
bnet = bnet_from_engine(engine);
N = length(bnet.node_sizes);
%[engine.jpot, loglik] = compute_joint_pot(bnet, mysetdiff(1:N, exclude), evidence, 1:N);
[engine.jpot] = compute_joint_pot(bnet, mysetdiff(1:N, exclude), evidence, 1:N);
% jpot should not be normalized, otherwise it gives wrong resutls for limids like asia_dt1
if nargout == 2
  [m] = marginal_nodes(engine, []);
  [T, lik] = normalize(m.T);
  loglik = log(lik);
end     
%[engine.jpot loglik] = normalize_pot(engine.jpot);
