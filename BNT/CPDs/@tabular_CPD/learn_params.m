function CPD = learn_params(CPD, fam, data, ns, cnodes)
%function CPD = learn_params(CPD, local_data)
% LEARN_PARAMS Compute the ML/MAP estimate of the params of a tabular CPD given complete data
% CPD = learn_params(CPD, local_data)
%
% local_data(i,m) is the value of i'th family member in case m (can be cell array).

local_data = data(fam, :); 
if iscell(local_data)
  local_data = cell2num(local_data);
end
counts = compute_counts(local_data, CPD.sizes);
switch CPD.prior_type
 case 'none', CPD.CPT = mk_stochastic(counts); 
 case 'dirichlet', CPD.CPT = mk_stochastic(counts + CPD.dirichlet); 
 otherwise, error(['unrecognized prior ' CPD.prior_type])
end
