function marginal = marginal_nodes(engine, nodes, t)
% MARGINAL_NODES Compute the marginal on the specified query nodes (jtree_unrolled_dbn)
% marginal = marginal_nodes(engine, nodes, t)
%
% 't' specifies the time slice of the earliest node in 'nodes'.
% 'nodes' must occur in some clique.
%
% Example:
% Consider a DBN with 2 nodes per slice.
% Then t=2, nodes=[1 3] refers to node 1 in slice 2 and node 1 in slice 3,
% i.e., nodes 3 and 5 in the unrolled network,

if nargin < 3, t = 1; end

bnet = bnet_from_engine(engine);
ss = length(bnet.intra);
query = nodes + (t-1)*ss;
marginal = marginal_nodes(engine.sub_engine, query);    
