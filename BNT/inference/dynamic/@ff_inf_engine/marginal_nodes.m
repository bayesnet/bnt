function marginal = marginal_nodes(engine, nodes, t)
% MARGINAL_NODES Compute the marginal on the specified query nodes (ff)
% marginal = marginal_nodes(engine, i, t)
% returns Pr(X(i,t) | Y(1:T)), where X(i,t) is the i'th node in the t'th slice.
% If enter_evidence used filtering instead of smoothing, this will return  Pr(X(i,t) | Y(1:t)).

if nargin < 3, t = 1; end
assert(length(nodes)==1);
i = nodes(end);
if myismember(i, engine.hnodes)
  marginal = pot_to_marginal(engine.marginals{i,t});
else
  marginal = pot_to_marginal(dpot(i, 1, 1)); % observed
end

bnet = bnet_from_engine(engine);
ss = length(bnet.intra);
% we convert the domain to the unrolled numbering system
% so that update_ess extracts the right evidence.
marginal.domain = nodes+(t-1)*ss;   
