function pot = convert_to_pot(CPD, pot_type, domain, evidence)
% CONVERT_TO_POT Convert a discrete CPD to a potential
% pot = convert_to_pot(CPD, pot_type, domain, evidence)
%
% pots = CPD evaluated using evidence(domain)

ncases = size(domain,2);
assert(ncases==1); % not yet vectorized

sz = dom_sizes(CPD);
ns = zeros(1, max(domain));
ns(domain) = sz;

CPT1 = CPD_to_CPT(CPD);
spar = issparse(CPT1);
odom = domain(~isemptycell(evidence(domain)));
if spar
   T = convert_to_sparse_table(CPD, domain, evidence);
else 
   T = convert_to_table(CPD, domain, evidence);
end

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
    if T(i) == 0 
      can{i} = cpot([], [], -Inf); % bug fix by Bob Welch 20/2/04
    else
      can{i} = cpot([], [], log(T(i)));
    end;
  end
  pot = cgpot(domain, [], ns, can); 
  
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

