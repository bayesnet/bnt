function pot = convert_to_pot(CPD, pot_type, domain, evidence)
% CONVERT_TO_POT Convert a tabular utility node to one or more potentials
% pot = convert_to_pot(CPD, pot_type, domain, evidence)

switch pot_type
 case 'u',
  sz = [CPD.sizes 1]; % the utility node itself has size 1
  pot = upot(domain, sz, 1*myones(sz), myreshape(CPD.T, sz));   
 otherwise,
  error(['can''t convert a utility node to a ' pot_type ' potential']);
end
