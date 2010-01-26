function marginal = marginal_nodes(engine, nodes, t)
% MARGINAL_NODES Compute the marginal on the specified query nodes (frontier)
% marginal = marginal_nodes(engine, nodes, t)
%
% 't' specifies the time slice of the earliest node in 'nodes'.
% 'nodes' cannot span more than 2 time slices.
%
% Example:
% Consider a DBN with 2 nodes per slice.
% Then t=2, nodes=[1 3] refers to node 1 in slice 2 and node 1 in slice 3,
% i.e., nodes 3 and 5 in the unrolled network,

if nargin < 3, t = 1; end
assert(length(nodes)==1);
i = nodes(1);
bigpot = engine.fwdback{i,t};
bnet = bnet_from_engine(engine);
ss  = length(bnet.intra);
nodes = nodes + (t-1)*ss;
%if t > 1, nodes = nodes + ss; end
marginal = pot_to_marginal(marginalize_pot(bigpot, nodes));
