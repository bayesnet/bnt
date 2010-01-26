function L = log_prior(CPD)
% LOG_PRIOR Return log P(theta) for a tabular CPD 
% L = log_prior(CPD)

switch CPD.prior_type
 case 'none',
  L = 0;
 case 'dirichlet',
  D = CPD.dirichlet(:);
  L = sum(log(D + (D==0)));
 case 'entropic',
  % log-prior = log exp(-H(theta)) = sum_i theta_i log (theta_i)
  fam_sz = CPD.sizes;
  psz = prod(fam_sz(1:end-1));
  ssz = fam_sz(end);
  C = reshape(CPD.CPT, psz, ssz);
  L = sum(sum(C .* log(C + (C==0))));
end
