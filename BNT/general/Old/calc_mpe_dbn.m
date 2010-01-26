function [mpe, ll] = calc_mpe_dbn(engine, evidence, break_ties)
% CALC_MPE Computes the most probable explanation of the evidence
% [mpe, ll] = calc_mpe_dbn(engine, evidence, break_ties)
%
% INPUT
% engine must support max-propagation
% evidence{i,t} is the observed value of node i in slice t, or [] if hidden
%
% OUTPUT
% mpe{i,t} is the most likely value of node i (cell array!)
% ll is the log-likelihood of the globally best assignment
%
% This currently only works when all hidden nodes are discrete

if nargin < 3, break_ties = 0; end

if break_ties
  disp('warning: break ties is ignored')
end

[engine, ll] = enter_evidence(engine, evidence, 'maximize', 1);

observed = ~isemptycell(evidence);
[ss T] = size(evidence);
scalar = 1;
N = length(evidence);
mpe = cell(ss,T);
bnet = bnet_from_engine(engine);
ns = bnet.node_sizes;
for t=1:T
  for i=1:ss
    m = marginal_nodes(engine, i, t);
    % observed nodes are all set to 1 inside the inference engine, so we must undo this
    if observed(i,t)
      mpe{i,t} = evidence{i,t};
    else
      assert(length(m.T) == ns(i));
      mpe{i,t} = argmax(m.T);
    end
  end
end
