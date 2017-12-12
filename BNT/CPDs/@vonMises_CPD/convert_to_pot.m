function pot = convert_to_pot( CPD, pot_type, domain, evidence )
%CONVERT_TO_POT Summary of this function goes here
%   Detailed explanation goes here


sz = CPD.sizes;
ns = zeros(1, max(domain));
ns(domain) = sz;

odom = domain(~isemptycell(evidence(domain)));
ps = domain(1:end-1);
cps = ps(CPD.cps);
dps = ps(CPD.dps);
self = domain(end);
cdom = [cps(:)' self];
ddom = dps;
cnodes = cdom;
  
switch pot_type
 case 'u'
  error('vonMises utility potentials not yet supported');
 
 case 'd'
  T = convert_to_table(CPD, domain, evidence);
  ns(odom) = 1;
  pot = dpot(domain, ns(domain), T);          

 case {'c','g'}
 [m,k,w] = vonMises_CPD_params_given_dps(CPD, domain, evidence);
  pot = linear_vonMises_to_cpot(m, k, w, domain, ns, cnodes, evidence);
 case 'cg'
  error('pot type not yet supported');

 case 'scg'
  error('pot type not yet supported');

 otherwise
  error(['unrecognized pot_type' pot_type])
end

end

