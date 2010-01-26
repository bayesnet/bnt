function CPD = update_params(CPD, ev, counts)
% UPDATE_PARAMS Update the Dirichlet pseudo counts and compute the new MAP param estimates (tabular)
%
% CPD = update_params(CPD, ev) uses the evidence on the family from a single case.
%
% CPD = update_params(CPD, [], counts) does a batch update using the specified suff. stats.

if nargin < 3
  n = length(ev);
  data = cat(1, ev{:}); % convert to a vector of scalars
  counts = compute_counts(data(:)', 1:n, mysize(CPD.CPT));
end
  
CPD.prior = CPD.prior + counts;
CPD.CPT = mk_stochastic(CPD.prior);
