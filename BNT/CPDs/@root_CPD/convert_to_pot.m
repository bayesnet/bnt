function pot = convert_to_pot(CPD, pot_type, domain, evidence)
% CONVERT_TO_POT Convert a root CPD to one or more potentials
% pots = convert_to_pot(CPD, pot_type, domain, evidence)

assert(length(domain)==1);
assert(~isempty(evidence(domain)));
T = 1;   

sz = CPD.sizes;
ns = zeros(1, max(domain));
ns(domain) = sz;

switch pot_type
 case 'u',
  pot = upot(domain, 1, T, 0);
 case 'd',
  ns(domain) = 1;
  pot = dpot(domain, ns(domain), T);          
 case {'c','g'},
  ns(domain) = 0;
  pot = cpot(domain, ns(domain), 0);
 case 'cg',
  ddom = [];
  cdom = domain; % we assume the root node is cts
  %pot = cgpot(ddom, cdom, ns, {cpot([],[],0)});
  pot = cgpot(ddom, cdom, ns);
end

