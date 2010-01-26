function joint = combine_marginals_into_joint(marginalsT, hnodes, ns)

jointT = dpot(hnodes, ns(hnodes));
for i=hnodes(:)'
  jointT = multiply_by_pot(jointT, marginalsT{i});
end
m = pot_to_marginal(jointT);
joint = m.T(:);           
