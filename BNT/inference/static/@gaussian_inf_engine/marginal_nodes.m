function marginal = marginal_nodes(engine, query)
% MARGINAL_NODES Compute the marginal on the specified query nodes (gaussian)
% marginal = marginal_nodes(engine, query)

% Compute sum_{Hsum} Pr(Hkeep, Hsum | o)
H = engine.hnodes;
bnet = bnet_from_engine(engine);
ns = bnet.node_sizes;
Hkeep = myintersect(H, query);
Hsum = mysetdiff(H, Hkeep);

[marginal.mu, marginal.Sigma] = marginalize_gaussian(engine.Hmu, engine.HSigma, Hkeep, Hsum, ns);
marginal.domain = query;
marginal.T = 1;

