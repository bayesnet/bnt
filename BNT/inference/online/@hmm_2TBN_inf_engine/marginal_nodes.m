function marginal = marginal_nodes(engine, b, nodes, t, add_ev)
% MARGINAL_NODES Compute the marginal on the specified nodes (hmm_2TBN)
% marginal = marginal_nodes(engine, b, nodes, t, add_ev)
%
% nodes must be a singleton set 

assert(length(nodes)==1)
ss = engine.slice_size;

i = nodes(1);
bigT = b.gamma;
dom = i + (t-1)*ss;

%id = engine.marg_singleton_ndx_id(i);
%global SD_NDX
%ndx = SD_NDX{id};
%marginal.T = marg_table_ndxSD(bigT, engine.maximize, ndx);

ns = engine.eff_node_sizes(:);
bigdom = 1:ss;
marginal.T = marg_table(bigT, bigdom + (t-1)*ss, ns(bigdom), dom, engine.maximize);

marginal.domain = dom;
assert(~add_ev);
%if add_ev
%  marginal = add_ev_to_dmarginal(marginal, engine.evidence, engine.node_sizes);
%end    
