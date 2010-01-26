function inter = learn_struct_dbn_reveal(seqs, ns, max_fan_in, penalty)
% LEARN_STRUCT_DBN_REVEAL Learn inter-slice adjacency matrix given fully observable discrete time series
% inter = learn_struct_dbn_reveal(seqs, node_sizes, max_fan_in, penalty)
% 
% seqs{l}{i,t} = value of node i in slice t of time-series l.
%   If you have a single time series in an N*T array D, use
%      seqs = { num2cell(D) }.
%   If you have L time series, each of length T, in an N*T*L array D, use
%      seqs= cell(1,L); for l=1:L, seqs{l} = num2cell(D(:,:,l)); end
%   or, in vectorized form,
%      seqs = squeeze(num2cell(num2cell(D),[1 2]));
% Currently the data is assumed to be discrete (1,2,...)
%
% node_sizes(i) is the number of possible values for node i
% max_fan_in is the largest number of parents we allow per node (default: N)
% penalty is weight given to the complexity penalty (default: 0.5)
%  A penalty of 0.5 gives the BIC score.
%  A penalty of 0 gives the ML score.
%  Maximizing likelihood is equivalent to maximizing mutual information between parents and child.
%
% inter(i,j) = 1 iff node in slice t connects to node j in slice t+1
%
% The parent set for each node in slice 2 is computed by evaluating all subsets of nodes in slice 1,
% and picking the largest scoring one. This takes O(n^k) time per node, where n is the num. nodes
% per slice, and k <= n is the max fan in.
% Since all the nodes are observed, we do not need to use an inference engine.
% And since we are only learning the inter-slice matrix, we do not need to check for cycles.
%
% This algorithm is described in
% - "REVEAL: A general reverse engineering algorithm for inference of genetic network
%      architectures", Liang et al. PSB 1998
% - "Extended dependency analysis of large systems",
%       Roger Conant, Intl. J. General Systems, 1988, vol 14, pp 97-141
% - "Learning the structure of DBNs", Friedman, Murphy and Russell, UAI 1998.

n = length(ns);

if nargin < 3, max_fan_in = n; end
if nargin < 4, penalty = 0.5; end

inter = zeros(n,n);

if ~iscell(seqs)
  data{1} = seqs;
end

nseq = length(seqs);
nslices = 0;
data = cell(1, nseq);
for l=1:nseq
  nslices = nslices + size(seqs{l}, 2);
  data{l} = cell2num(seqs{l})'; % each row is a case
end
ndata = nslices - nseq; % subtract off the initial slice of each sequence

% We concatenate the sequences as in the following example.
% Let there be 2 sequences of lengths 4 and 5, with n nodes per slice,
% and let i be the target node.
% Then we construct following matrix D 
%
% s{1}{1,1} ... s{1}{1,3}     s{2}{1,1} ... s{2}{1,4}
% ....
% s{1}{n,1} ... s{1}{n,3}     s{2}{n,1} ... s{2}{n,4}
% s{1}{i,2} ... s{1}{i,4}     s{2}{i,2} ... s{2}{i,5}
%
% D(1:n, i) is the i'th input and D(n+1, i) is the i'th output.
% 
% We concatenate each sequence separately to avoid treating the transition
% from the end of one sequence to the beginning of another as a "normal" transition.


for i=1:n
  D = [];
  for l=1:nseq
    T = size(seqs{l}, 2);
    A = cell2num(seqs{l}(:, 1:T-1));
    B = cell2num(seqs{l}(i, 2:T));
    C = [A;B];
    D = [D C];
  end
  SS = subsets(1:n, max_fan_in, 1); % skip the empty set 
  nSS = length(SS);
  bic_score = zeros(1, nSS);
  ll_score = zeros(1, nSS);
  target = n+1;
  ns2 = [ns ns(i)];
  for h=1:nSS
    ps = SS{h};
    dom = [ps target];
    counts = compute_counts(D(dom, :), ns2(dom));
    CPT = mk_stochastic(counts);
    [bic_score(h), ll_score(h)] = bic_score_family(counts, CPT, ndata);
  end
  if penalty == 0
    h = argmax(ll_score);
  else
    h = argmax(bic_score);
  end
  ps = SS{h};
  inter(ps, i) = 1;
end
