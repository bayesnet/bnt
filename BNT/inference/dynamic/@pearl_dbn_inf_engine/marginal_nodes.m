function marginal = marginal_nodes(engine, nodes, t)
% MARGINAL_NODES Compute the marginal on the specified query nodes (pearl_dbn)
% marginal = marginal_nodes(engine, i, t)
% returns Pr(X(i,t) | Y(1:T)), where X(i,t) is the i'th node in the t'th slice.
% If enter_evidence used filtering instead of smoothing, this will return  Pr(X(i,t) | Y(1:t)).

if nargin < 3, t = 1; end
assert(length(nodes)==1);
i = nodes(end);
if ~myismember(i, engine.onodes)
  marginal.T = engine.marginal{i,t};
else
  marginal.T = 1; % observed
end

% we convert the domain to the unrolled numbering system
% so that update_ess extracts the right evidence.
marginal.domain = nodes+(t-1)*engine.ss;       
