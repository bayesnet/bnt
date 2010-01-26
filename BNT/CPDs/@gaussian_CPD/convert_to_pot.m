function pot = convert_to_pot(CPD, pot_type, domain, evidence)
% CONVERT_TO_POT Convert a Gaussian CPD to one or more potentials
% pot = convert_to_pot(CPD, pot_type, domain, evidence)

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
 case 'u',
  error('gaussian utility potentials not yet supported');
 
 case 'd',
  T = convert_to_table(CPD, domain, evidence);
  ns(odom) = 1;
  pot = dpot(domain, ns(domain), T);          

 case {'c','g'},
  [m, C, W] = gaussian_CPD_params_given_dps(CPD, domain, evidence);
  pot = linear_gaussian_to_cpot(m, C, W, domain, ns, cnodes, evidence);

 case 'cg',
  [m, C, W] = gaussian_CPD_params_given_dps(CPD, domain, evidence);
  % Convert each conditional Gaussian to a canonical potential
  cobs = myintersect(cdom, odom);
  dobs = myintersect(ddom, odom);
  ens = ns; % effective node size
  ens(cobs) = 0;
  ens(dobs) = 1;
  dpsize = prod(ens(dps));
  can = cell(1, dpsize);
  for i=1:dpsize
    if isempty(W)
      can{i} = linear_gaussian_to_cpot(m(:,i), C(:,:,i), [], cdom, ns, cnodes, evidence);
    else
      can{i} = linear_gaussian_to_cpot(m(:,i), C(:,:,i), W(:,:,i), cdom, ns, cnodes, evidence);
    end
  end
  pot = cgpot(ddom, cdom, ens, can);

 case 'scg',
  [m, C, W] = gaussian_CPD_params_given_dps(CPD, domain, evidence);
  cobs = myintersect(cdom, odom);
  dobs = myintersect(ddom, odom);
  ens = ns; % effective node size
  ens(cobs) = 0;
  ens(dobs) = 1;
  dpsize = prod(ens(dps));
  cpsize = size(W, 2); % cts parents size
  ss = size(m, 1); % self size
  cheaddom = self;
  ctaildom = cps(:)';
  pot_array = cell(1, dpsize);
  for i=1:dpsize
    pot_array{i} = scgcpot(ss, cpsize, 1, m(:,i), W(:,:,i), C(:,:,i));
  end
  pot = scgpot(ddom, cheaddom, ctaildom, ens, pot_array);

 otherwise,
  error(['unrecognized pot_type' pot_type])
end

