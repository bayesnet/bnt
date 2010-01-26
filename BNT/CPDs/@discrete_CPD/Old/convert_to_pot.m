function pot = convert_to_pot(CPD, pot_type, domain, evidence)
% CONVERT_TO_POT Convert a tabular CPD to one or more potentials
% pots = convert_to_pot(CPD, pot_type, domain, evidence)
%
% pots{i} = CPD evaluated using evidence(domain(:,i))
% If 'domains' is a single row vector, pots will be an object, not a cell array.

ncases = size(domain,2);
assert(ncases==1); % not yet vectorized

sz = dom_sizes(CPD);
ns = zeros(1, max(domain));
ns(domain) = sz;

local_ev = evidence(domain);
obs_bitv = ~isemptycell(local_ev);
odom = domain(obs_bitv);
T = convert_to_table(CPD, domain, local_ev, obs_bitv);

switch pot_type
 case 'u',
  pot = upot(domain, sz, T, 0*myones(sz));  
 case 'd',
  ns(odom) = 1;
  pot = dpot(domain, ns(domain), T);          
 case {'c','g'},
  % Since we want the output to be a Gaussian, the whole family must be observed.
  % In other words, the potential is really just a constant.
  p = T;
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
 otherwise,
  error(['unrecognized pot type ' pot_type])
end

