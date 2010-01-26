function m = pot_to_marginal(pot)
% POT_TO_MARGINAL Convert a dpot to a marginal structure.
% m = pot_to_marginal(pot)

m.domain = pot.domain;
m.T = pot.T;
m.mu = [];
m.Sigma = [];  

%if isvector(m.T)
%  m.T = m.T(:);
%end
