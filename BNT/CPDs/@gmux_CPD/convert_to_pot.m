function pot = convert_to_pot(CPD, pot_type, domain, evidence)
% CONVERT_TO_POT Convert a gmux CPD to a Gaussian potential
% pot = convert_to_pot(CPD, pot_type, domain, evidence)
  
switch pot_type
 case {'d', 'u', 'cg', 'scg'},
  error(['can''t convert gmux to potential of type ' pot_type])

 case {'c','g'},
  % We create a large weight matrix with zeros in all blocks corresponding
  % to the non-chosen parents, since they are effectively disconnected.
  % The chosen parent is determined by the value, m,  of the discrete parent.
  % Thus the potential is as large as the whole family.
  ps = domain(1:end-1);
  dps = ps(CPD.dps); % CPD.dps is an index, not a node number (because of param tying)
  cps = ps(CPD.cps);
  m = evidence{dps};
  if isempty(m)
    error('gmux node must have observed discrete parent')
  end
  bs = CPD.sizes(CPD.cps);
  b = block(m, bs);
  sum_cpsz = sum(CPD.sizes(CPD.cps));
  selfsz = CPD.sizes(end);
  W = zeros(selfsz, sum_cpsz);
  W(:,b) = CPD.weights(:,:,m);

  ns = zeros(1, max(domain));
  ns(domain) = CPD.sizes;
  self = domain(end);
  cdom = [cps(:)' self];
  pot = linear_gaussian_to_cpot(CPD.mean(:,m), CPD.cov(:,:,m), W, domain, ns, cdom, evidence);
  
 otherwise,
  error(['unrecognized pot_type' pot_type])
end

