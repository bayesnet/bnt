function m = pot_to_marginal(pot)
% POT_TO_MARGINAL Convert a scgpot to a marginal structure.
% m = pot_to_marginal(pot)

assert(isempty(pot.ctaildom))
m.domain = pot.domain;
n = pot.cheadsize;
d = pot.dsize;

if n==0
  m.mu = [];
  m.Sigma = [];
else
  m.mu = zeros(n, d);
  m.Sigma = zeros(n, n, d);
end
%m.T = 0*myones(pot.dsizes);
m.T = 0*myones(pot.dsize);
for i=1:pot.dsize
  potc = struct(pot.scgpotc{i}); % violate privacy of object
  if n > 0
    m.mu(:,i) = potc.A;
    m.Sigma(:,:,i) = potc.C;
  end
  m.T(i) = potc.p;
end     
if isvectorBNT(m.T)
  m.T = m.T(:)';
end
