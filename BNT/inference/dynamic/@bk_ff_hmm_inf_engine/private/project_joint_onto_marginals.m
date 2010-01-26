function [marginals, marginalsT] = project_joint_onto_marginals(joint, hnodes, ns)

ss = length(ns);
jointT = dpot(hnodes, ns(hnodes), joint);
marginalsT = cell(1, ss);
marginals = cell(1,ss);
for i=hnodes(:)'
  marginalsT{i} = marginalize_pot(jointT, i);
  m = pot_to_marginal(marginalsT{i});
  marginals{i} = m.T(:);
end
