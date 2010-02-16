function m = pot_to_marginal(pot)
% POT_TO_MARGINAL Convert a cgpot to a marginal structure.
% m = pot_to_marginal(pot)

pot = cg_can_to_mom(pot);
m.domain = pot.domain;
n = pot.csize;
d = length(pot.mom);
if n==0
  m.mu = [];
  m.Sigma = [];
else
  m.mu = zeros(n, d);
  m.Sigma = zeros(n, n, d);
end
m.T = 0*myones(pot.dsizes);
for i=1:pot.dsize
  s = struct(pot.mom{i}); % violate privacy of object
  if n > 0
    m.mu(:,i) = s.mu;
    m.Sigma(:,:,i) = s.Sigma;
  end
  m.T(i) = exp(s.logp);
end     
if isvectorBNT(m.T)
  m.T = m.T(:)';
end
