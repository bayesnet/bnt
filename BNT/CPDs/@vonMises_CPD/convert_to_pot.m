function pot = convert_to_pot( CPD, pot_type, domain, evidence )
%CONVERT_TO_POT Convert a Von Mises CPD to one or more potentials depending
%on the pot_type.
%   Currently the case 'd' and case 'c','g' are functioning correctly and
%   have been tested. However, case 'cg' is in the midst of being tested.

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
   error('vonMises cg section still testing...');
  [m, k, W] = vonMises_CPD_params_given_dps(CPD, domain, evidence);
  % Convert each conditional Von Mises to a canonical potential
  cobs = myintersect(cdom, odom);
  dobs = myintersect(ddom, odom);
  ens = ns; % effective node size
  ens(cobs) = 0;
  ens(dobs) = 1;
  dpsize = prod(ens(dps));
  can = cell(1, dpsize);
  for i=1:dpsize
    if isempty(W)
      can{i} = linear_vonMises_to_cpot(m(:,i), k(:,:,i), [], cdom, ns, cnodes, evidence);
    else
      can{i} = linear_vonMises_to_cpot(m(:,i), k(:,:,i), W(:,:,i), cdom, ns, cnodes, evidence);
    end
  end
  pot = cgpot(ddom, cdom, ens, can);

 case 'scg'
  error('pot type not yet supported');

 otherwise
  error(['unrecognized pot_type' pot_type])
end

end

