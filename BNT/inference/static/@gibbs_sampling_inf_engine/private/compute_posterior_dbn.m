function post = compute_posterior_dbn(bnet, state, i, n, strides, families, ...
				  CPT)
% COMPUTE_POSTERIOR
%
% post = compute_posterior(bnet, state, i, n, strides, families,
% cpts)
%
% Compute the posterior distribution on node X_i^n of a DBN,
% conditional on evidence in the cell array state
%
% strides is the cached result of compute_strides(bnet)
% families is the cached result of compute_families(bnet)
% cpt is the cached result of get_cpts(bnet)
%
% post is a one-dimensional table



% First multiply in the cpt of the node itself
post = get_slice_dbn(bnet, state, i, n, i, n, strides, families, CPT);
post = post(:);

% Then multiply in CPTs of children that are in this slice
for j = children(bnet.intra, i)
  slice = get_slice_dbn(bnet, state, j, n, i, n, strides, families, CPT);
  post = post.*slice(:);
end

% Finally, if necessary, multiply in CPTs of children in the next
% slice 
if (n < size(state,2))
  for j = children(bnet.inter, i)
    slice = get_slice_dbn(bnet, state, j, n+1, i, n, strides, families, ...
			    CPT);
    post = post.*slice(:);
  end
end

post = normalise(post);




















