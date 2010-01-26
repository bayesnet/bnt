function [jpot, loglik] = compute_joint_pot(bnet, nodes, evidence, domain)
% COMPUTE_JOINT_POT Compute the global joint potential of a Bayes net
% function jpot = compute_joint_pot(bnet, nodes, evidence, domain)

if nargin < 4, domain = nodes; end

onodes = find(~isemptycell(evidence));
pot_type = determine_pot_type(bnet, onodes, domain);

jpot = mk_initial_pot(pot_type, domain, bnet.node_sizes, bnet.cnodes, onodes);
for i=nodes(:)'
  e = bnet.equiv_class(i);
  fam = family(bnet.dag, i);
  pot = convert_to_pot(bnet.CPD{e}, pot_type, fam(:), evidence);
  jpot = multiply_by_pot(jpot, pot);
end                                                  
%[jpot, loglik] = normalize_pot(jpot); % causes errors in asia_dt1 etc
