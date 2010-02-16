function m = pot_to_marginal(pot)
% POT_TO_MARGINAL Convert a mpot to a marginal structure.
% m = pot_to_marginal(pot)

m.domain = pot.domain;
m.T = exp(pot.logp);
m.mu = pot.mu;
m.Sigma = pot.Sigma;

if isvectorBNT(m.T)
  m.T = m.T(:)';
end

