function marginal = marginal_nodes(engine, nodes)
% MARGINAL_NODES Compute the marginal on the specified query nodes (likelihood_weighting)
% marginal = marginal_nodes(engine, nodes)

bnet = bnet_from_engine(engine);
ddom = myintersect(nodes, bnet.dnodes);
cdom = myintersect(nodes, bnet.cnodes);
nsamples = size(engine.samples, 1);
ns = bnet.node_sizes;

%w = normalise(engine.weights);
w = engine.weights;
if mysubset(nodes, ddom)
  T = 0*myones(ns(nodes));
  P = prod(ns(nodes));
  indices = ind2subv(ns(nodes), 1:P);
  samples = reshape(cat(1, engine.samples{:,nodes}), nsamples, length(nodes));
  for j = 1:P
    rows = find_rows(samples, indices(j,:));
    T(j) = sum(w(rows));
  end
  T = normalise(T);
  marginal.T = T;
elseif subset(nodes, cdom)
  samples = reshape(cat(1, engine.samples{:,nodes}), nsamples*sum(ns(nodes)), length(nodes));
  [marginal.mu, marginal.Sigma] =  wstats(samples', normalise(w));
else
  error('can''t handle mixed marginals yet');
end

marginal.domain = nodes;

%%%%%%%%%

function rows = find_rows(M, v)
% FINDROWS Find rows which are equal to a specified vector
% rows = findrows(M, v)
% Each row of M is a sample

temp = abs(M - repmat(v, size(M, 1), 1));
rows = find(sum(temp,2) == 0);      

%%%%%%%%

function [mu, Sigma] = wstats(X, w)

% Computes the weighted mean and weighted covariance matrix for a given
% set of observations X(:,i), and a set of normalised weights w(i).
% Each column of X is a sample.

d = X - repmat(X * w', 1, size(X, 2));
mu = sum(X .* repmat(w, size(X, 1), 1), 2);
Sigma = d * diag(w) * d';          
