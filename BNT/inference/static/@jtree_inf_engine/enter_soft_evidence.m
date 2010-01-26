function [clpot, loglik] = enter_soft_evidence(engine, clique, potential, onodes, pot_type)
% ENTER_SOFT_EVIDENCE Add the specified potentials to the network (jtree)
% [clpot, loglik] = enter_soft_evidence(engine, clique, potential, onodes, pot_type, maximize)
%
% We multiply potential{i} onto clique(i) before propagating.
% We return all the modified clique potentials.

% only used by BK!

[clpot, seppot] = init_pot(engine, clique, potential, pot_type, onodes);
[clpot, seppot] = collect_evidence(engine, clpot, seppot);
[clpot, seppot] = distribute_evidence(engine, clpot, seppot);

C = length(clpot);
ll = zeros(1, C);
for i=1:C
  [clpot{i}, ll(i)] = normalize_pot(clpot{i});
end
loglik = ll(1); % we can extract the likelihood from any clique


