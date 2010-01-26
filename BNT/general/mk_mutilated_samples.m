function [data, clamped] = mk_mutilated_samples(bnet, ncases, max_clamp, usecell)
% GEN_MUTILATED_SAMPLES Do random interventions and then draw random samples
% [data, clamped] = gen_mutilated_samples(bnet, ncases, max_clamp, usecell)
%
% At each step, we pick a random subset of size 0 .. max_clamp, and 
% clamp these nodes to random values.
%
% data(i,m) is the value of node i in case m.
% clamped(i,m) = 1 if node i in case m was set by intervention.

if nargin < 4, usecell = 1; end

ns = bnet.node_sizes;
n = length(bnet.dag);
if usecell
  data = cell(n, ncases);
else
  data = zeros(n, ncases);
end
clamped = zeros(n, ncases);

csubsets = subsets(1:n, max_clamp, 0); % includes the empty set
distrib_cset = normalise(ones(1, length(csubsets)));

for m=1:ncases
  cset = csubsets{sample_discrete(distrib_cset)};
  nvals = prod(ns(cset));
  distrib_cvals = normalise(ones(1, nvals));
  cvals = ind2subv(ns(cset), sample_discrete(distrib_cvals));
  mutilated_bnet = do_intervention(bnet, cset, cvals);
  ev = sample_bnet(mutilated_bnet);
  if usecell
    data(:,m) = ev;
  else
    data(:,m) = cell2num(ev);
  end
  clamped(cset,m) = 1;
end
