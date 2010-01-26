function marginal = marginal_nodes(engine, nodes, t, add_ev)
% MARGINAL_NODES Compute the marginal on the specified query nodes (hmm)
% marginal = marginal_nodes(engine, nodes, t, add_ev)
%
% 'nodes' must be a single node.
% t is the time slice.

if nargin < 3, t = 1; end
if nargin < 4, add_ev = 0; end

assert(length(nodes)==1)
ss = engine.slice_size;

i = nodes(1);
bigT = engine.one_slice_marginal(:,t);
dom = i + (t-1)*ss;

ns = engine.eff_node_sizes(:);
bigdom = 1:ss;
marginal.T = marg_table(bigT, bigdom + (t-1)*ss, ns(bigdom), dom, engine.maximize);

marginal.domain = dom;
marginal.mu = [];
marginal.Sigma = [];

if add_ev
  marginal = add_ev_to_dmarginal(marginal, engine.evidence, engine.node_sizes);
end    
 
