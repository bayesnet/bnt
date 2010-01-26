function [engine, loglik] = enter_evidence(engine, evidence, varargin)
% ENTER_EVIDENCE Add the specified evidence to the network (gaussian_inf_engine)
% [engine, loglik] = enter_evidence(engine, evidence, ...)
%
% evidence{i} = [] if if X(i) is hidden, and otherwise contains its observed value (scalar or column vector)

bnet = bnet_from_engine(engine);
ns = bnet.node_sizes;
O = find(~isemptycell(evidence));
H = find(isemptycell(evidence));
vals = cat(1, evidence{O});

% Compute Pr(H|o)
[Hmu, HSigma, loglik] = condition_gaussian(engine.mu, engine.Sigma, H, O, vals(:), ns);

engine.Hmu = Hmu;
engine.HSigma = HSigma;
engine.hnodes = H;

%%%%%%%%

function [mu2, Sigma2, loglik] = condition_gaussian(mu, Sigma, X, Y, y, ns)
% CONDITION_GAUSSIAN Compute Pr(X|Y=y) where X and Y are jointly Gaussian.
% [mu2, Sigma2, ll] = condition_gaussian(mu, Sigma, X, Y, y, ns)

if isempty(y)
  mu2 = mu;
  Sigma2 = Sigma;
  loglik = 0;
  return;
end

use_log = 1;

if length(Y)==length(mu) % instantiating every variable
  mu2 = y;
  Sigma2 = zeros(length(y));
  loglik = gaussian_prob(y, mu, Sigma, use_log);
  return;
end

[muX, muY, SXX, SXY, SYX, SYY] = partition_matrix_vec(mu, Sigma, X, Y, ns);
K = SXY*inv(SYY);
mu2 = muX + K*(y-muY);
Sigma2 = SXX - K*SYX;
loglik = gaussian_prob(y, muY, SYY, use_log);
