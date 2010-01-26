function pot = convert_to_pot(CPD, pot_type, domain, evidence)
% CONVERT_TO_POT Convert a softmax CPD to a potential
% pots = convert_to_pot(CPD, pot_type, domain, evidence)
%
% pots = CPD evaluated using evidence(domain)

ncases = size(domain,2);
assert(ncases==1); % not yet vectorized

sz = dom_sizes(CPD);
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
  ps = domain(1:end-1);
  dps = ps(CPD.dpndx);
  cps = ps(CPD.cpndx);
  ddom = [dps CPD.self];
  cdom = cps;
  pot = cgpot(ddom, cdom, ns, can);   
  
 case 'scg'
  T = T(:);
  ns(odom) = 1;
  pot_array = cell(1, length(T));
  for i=1:length(T)
    pot_array{i} = scgcpot([], [], T(i));
  end
  pot = scgpot(domain, [], [], ns, pot_array);   

 otherwise,
  error(['unrecognized pot type ' pot_type])
end

