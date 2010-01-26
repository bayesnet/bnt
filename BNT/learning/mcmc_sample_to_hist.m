function mcmc_post = mcmc_sample_to_hist(sampled_graphs, dags)
% MCMC_SAMPLE_TO_HIST Convert a set of sampled dags into a histogram over dags
% hist = mcmc_sample_to_hist(sampled_graphs, dags)
%
% sampled_graphs{m} is the m'th sampled dag
% dags{i} is the i'th dag in the hypothesis space
% hist(i) = Pr(model i | data)

ndags = length(dags);
nsamples = length(sampled_graphs);
nnodes = length(dags{1});
% sampled_bitv(m, :) is the m'th sampled graph represented as a vector of n^2 bits, computed
% by stacking the columns of the adjacency matrix vertically.
sampled_bitvs = zeros(nsamples, nnodes*nnodes);
for m=1:nsamples
  sampled_bitvs(m, :) = sampled_graphs{m}(:)';
end
  
[ugraphs, I, J] = unique(sampled_bitvs, 'rows');  % each row of ugraphs is a unique bit vector
sampled_indices  = subv2ind(2*ones(1,nnodes*nnodes), ugraphs+1);
counts = hist(J, 1:size(ugraphs,1)); % counts(i) = number of times graphs(i,:) occurs in the sample

mcmc_post = zeros(1, ndags);
for i=1:ndags
  bitv = dags{i}(:)';
  % Find the samples that corresponds to this graph by converting the graphs to bitvectors and
  % then to integers.
  ndx = subv2ind(2*ones(1,nnodes*nnodes), bitv+1);
  locn = find(ndx == sampled_indices);
  if ~isempty(locn)
    mcmc_post(i) = counts(locn);
  end
end
mcmc_post = normalise(mcmc_post);


