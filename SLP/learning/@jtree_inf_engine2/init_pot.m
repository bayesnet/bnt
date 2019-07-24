function [clpot, seppot] = init_pot(engine, clqs, pots, pot_type, onodes, ndx)
% INIT_POT Initialise potentials with evidence (jtree_inf)
% function [clpot, seppot] = init_pot(engine, clqs, pots, pot_type, onodes)

cliques = engine.cliques;
bnet = bnet_from_engine(engine);
% Set the clique potentials to all 1s
C = length(cliques);
clpot = cell(1,C);
for i=1:C
  clpot{i} = mk_initial_pot(pot_type, cliques{i}, bnet.node_sizes(:), bnet.cnodes(:), onodes);
end

% Multiply on specified potentials
for i=1:length(clqs)
  c = clqs(i);
  clpot{c} = multiply_by_pot(clpot{c}, pots{i});
end

seppot = cell(C,C); % implicitely initialized to 1
