function pot = CPD_to_pot(pot_type, CPD, domain, ns, cnodes, evidence)
% CPD_TO_POT Convert a CPD to a potential of the specified form, incorporating any evidence
% pot = CPD_to_pot(pot_type, CPD, domain, node_sizes, cnodes, evidence)
%
% pot_type is one of 'd', 'g', or 'cg'.
% domain is the domain of CPD.
% node_sizes(i) is the size of node i.
% cnodes = the cts nodes
% evidence{i} is the evidence on the i'th node.

switch pot_type
 case 'd',
  pot = CPD_to_dpot(CPD, domain, ns, cnodes, evidence);
 case 'g',
  pot = CPD_to_cpot(CPD, domain, ns, cnodes, evidence);
 case 'cg',
  pot = CPD_to_cgpot(CPD, domain, ns, cnodes, evidence);
 otherwise,
  error(['can''t handle pot_type ' pot_type]);
end
              
