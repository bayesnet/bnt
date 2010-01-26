function marginal = marginal_nodes(engine, nodes, t)
% MARGINAL_NODES Compute the marginal on the specified query nodes (bk_ff_hmm)
% marginal = marginal_nodes(engine, i, t)

assert(length(nodes)==1);
i = nodes(end);
%assert(myismember(i, engine.hnodes));
marginal = engine.marginals{i,t};
bnet = bnet_from_engine(engine);
ss = length(bnet.intra);
marginal.domain = i + (t-1)*ss;
