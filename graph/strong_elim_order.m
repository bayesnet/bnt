function order = strong_elim_order(G, node_sizes, partial_order)
% STRONG_ELIM_ORDER Find an elimination order to produce a strongly triangulated graph.
% order = strong_elim_order(moral_graph, node_sizes, partial_order)
% 
% partial_order(i,j)=1 if we must marginalize i *after* j
% (so i will be nearer the strong root).
% e.g., if j is a decision node and i is its information set:
%   we cannot maximize j if we have marginalized out some of i
% e.g., if j is a continuous child and i is its discrete parent:
%   we want to integrate out the cts nodes before the discrete ones,
%   so that the marginal is strong.
%
% For details, see
% - Jensen, Jensen and Dittmer, "From influence diagrams to junction trees", UAI 94.
% - Lauritzen, "Propgation of probabilities, means, and variances in mixed graphical
%   association models", JASA 87(420):1098--1108, 1992.
%
% On p369 of the Jensen paper, they state "the reverse of the elimination order must be some
% extension of [the partial order] to a total order".
% We make no attempt to find the best such total ordering, in the sense of minimizing the weight
% of the resulting cliques.

% Example from the Jensen paper:
% Let us number the nodes in Fig 1 from top to bottom, left to right,
% so a=1,b=2,D1=3,c=4,...,l=14,j=15,k=16.
% The elimination ordering they propose on p370 is [14 15 16 11 12 1 4 5 10 8 13 9 7 6 3 2];

if 0
  total_order = topological_sort(partial_order);
  order = total_order(end:-1:1); % no attempt to find an optimal constrained ordering!
  return;
end

% The following implementation is due to Ilya Shpitser and seems to give wrong
% results on cg1

n = length(G);
MG = G; % copy the original graph
uneliminated = ones(1,n);
order = zeros(1,n);

for i=1:n
  roots = [];
  k = 1;
  for j=1:n
    if sum(partial_order(j,:)) == 0
      roots(k) = j;
      k = k + 1;
    end
  end
  U = find(uneliminated);
  valid = myintersect(U, roots);
  % Choose the best node from the set of valid candidates
  score1 = zeros(1,length(valid));
  score2 = zeros(1,length(valid));
  for j=1:length(valid)
    k = valid(j);
    ns = myintersect(neighbors(G, k), U);
    l = length(ns);
    M = MG(ns,ns);
    score1(j) = l^2 - sum(M(:)); % num. added edges
    score2(j) = prod(node_sizes([k ns])); % weight of clique
  end
  j1s = find(score1==min(score1));
  j = j1s(argmin(score2(j1s)));
  k = valid(j);
  uneliminated(k) = 0;
  order(i) = k;
  ns = myintersect(neighbors(G, k), U);
  if ~isempty(ns)
    G(ns,ns) = 1;
    G = setdiag(G,0);
  end
  partial_order(:,k) = 0;
end
