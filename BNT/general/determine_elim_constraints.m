function partial_order = determine_elim_constraints(bnet, onodes)
% DETERMINE_ELIM_CONSTRAINTS Determine what the constraints are (if any) on the elimination ordering.
% partial_order = determine_elim_constraints(bnet, onodes)
%
% A graph with different kinds of nodes (e.g., discrete and cts, or decision and rnd) is called marked. 
% A strong root is guaranteed to exist if the marked graph is triangulated and does not have any paths of
% the form discrete -> cts -> discrete. In general we need to add extra edges to
% the moral graph to ensure this (see example in Lauritzen (1992) fig 3b).
% However, a simpler sufficient condition is to eliminate all the cts nodes before the discrete ones,
% because then, as we move from the leaves to the root, the cts nodes get marginalized away
% and we are left with purely discrete cliques.
%
% partial_order(i,j)=1 if we must marginalize j *before* i
% (so i will be nearer the strong root).
% If the hidden nodes are either all discrete or all cts, we set partial_order = [].
%
% For details, see
% - Jensen, Jensen and Dittmer, "From influence diagrams to junction trees", UAI 94.
% - Lauritzen, "Propgation of probabilities, means, and variances in mixed graphical
%   association models", JASA 87(420):1098--1108, 1992.
% - K. Olesen, "Causal probabilistic networks with both discrete and continuous variables",
%      IEEE Pami 15(3), 1993


n = length(bnet.dag);
pot_type = determine_pot_type(bnet, onodes);
if (pot_type == 'd') | (pot_type == 'g')
  partial_order = [];
  return;
end


partial_order = sparse(n,n);
partial_order(bnet.dnodes, bnet.cnodes) = 1;

% Integrate out cts nodes before their discrete parents - see Olesen (1993) p9
% This method gives the wrong results on cg1.m!
if 0
for i=bnet.cnodes(:)'
  dps = myintersect(parents(bnet.dag, i), bnet.dnodes);
  partial_order(dps, i)=1;
end
end
