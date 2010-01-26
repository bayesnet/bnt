function [dag, order] = mk_rnd_dag(N, max_fan_in)
% MY_MK_RND_DAG  Create a random directed acyclic graph
%
% [dag, order] = my_mk_rnd_dag(N, max_fan_in)
%  max_fan_in defaults to N.
%  order is the random topological order that was chosen

% Modified by Sonia Leach 2/25/02

if nargin < 2, max_fan_in = N; end

order = randperm(N);
dag = zeros(N,N);
for i=2:N
  j = order(i);
  %k = sample_discrete(normalise(ones(1, min(i-1, max_fan_in))));
  k = sample_discrete(normalise(ones(1, min(i-1, max_fan_in)+1))) - 1; % min = 0 (bug fix due to
                                                                       % Pedrito, 7/28/04)
  SS = order(1:i-1);          % get Set of possible parentS
  p  = randperm(length(SS));  % permute order of set
  dag(SS(p(1:k)),j) = 1;      % take first k in permuted order

  % Kevin had:
  %SS = subsets(order(1:i-1), k, k);
  %p = sample_discrete(normalise(ones(1, length(SS))));
  %dag(SS{p}, j) = 1;
end
