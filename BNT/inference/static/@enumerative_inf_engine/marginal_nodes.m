function [marginal, loglik] = marginal_nodes(engine, query)
% MARGINAL_NODES Compute the marginal on the specified query nodes (enumerative_inf)
% [marginal, loglik] = marginal_nodes(engine, query)


if isempty(query) & nargout < 2
  marginal.T = 1;
  marginal.domain = [];
  return;
end

evidence = engine.evidence;
bnet = bnet_from_engine(engine);
assert(isempty(bnet.cnodes));
n = length(bnet.dag);
observed = ~isemptycell(evidence);
vals = cat(1,evidence{observed});
vals = vals(:)';
ns = bnet.node_sizes;

sz = ns(query);
T = 0*myones(sz);
p = 0;
for i=1:prod(ns)
  inst = ind2subv(ns, i); % i'th instantiation
  if isempty(vals) | inst(observed) == vals % agrees with evidence
    prob = exp(log_lik_complete(bnet, num2cell(inst(:))));
    p = p + prob;
    v = inst(query);
    j = subv2ind(sz, v);
    T(j) = T(j) + prob;
  end
end

[T, lik] = normalise(T);
lik = p;
loglik = log(lik);

Tsmall = shrink_obs_dims_in_table(T, query, evidence);
marginal.domain = query;
marginal.T = Tsmall;
