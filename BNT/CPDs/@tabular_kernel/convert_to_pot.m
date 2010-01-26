function pot = convert_to_pot(CPD, pot_type, domain, evidence)
% CONVERT_TO_POT Convert a tabular CPD to one or more potentials
% pot = convert_to_pot(CPD, pot_type, domain, evidence)

% This is the same as discrete_CPD/convert_to_pot,
% except we didn't want to the kernel to inherit methods like sample_node etc.

sz = CPD.sz;
ns = zeros(1, max(domain));
ns(domain) = sz;

odom = domain(~isemptycell(evidence(domain)));
T = convert_to_table(CPD, domain, evidence);

switch pot_type
 case 'u',
  pot = upot(domain, sz, T, 0*myones(sz));  
 case 'd',
  ns(odom) = 1;
  pot = dpot(domain, ns(domain), T);          
 case 'c',
  % Since we want the output to be a Gaussian, the whole family must be observed.
  % In other words, the potential is really just a constant.
  p = T.p;
  %p = prob_node(CPD, evidence(domain(end)), evidence(domain(1:end-1)));
  ns(domain) = 0;
  pot = cpot(domain, ns(domain), log(p));       
 case 'cg',
  T = T(:);
  ns(odom) = 1;
  can = cell(1, length(T));
  for i=1:length(T)
    can{i} = cpot([], [], log(T(i)));
  end
  pot = cgpot(domain, [], ns, can);   
end

