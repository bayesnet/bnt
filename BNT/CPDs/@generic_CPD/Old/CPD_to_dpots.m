function pots = CPD_to_dpots(CPD, domain, ns, cnodes, evidence)
% CPD_TO_DPOTS Convert the CPD to several discrete potentials, for different instantiations (generic)
% pots = CPD_to_dpots(CPD, domain, ns, cnodes, evidence)
%
% domain(:,i) is the domain of the i'th instantiation of CPD.
% node_sizes(i) is the size of node i.
% cnodes = all the cts nodes
% evidence{i} is the evidence on the i'th node.
%
% This just calls CPD_to_dpot for each domain.
    
nCPDs = size(domain,2);
pots = cell(1,nCPDs);
for i=1:nCPDs
  pots{i} = CPD_to_dpot(CPD, domain(:,i), ns, cnodes, evidence);
end
