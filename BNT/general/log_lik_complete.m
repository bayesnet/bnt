function L = log_lik_complete(bnet, cases, clamped)
% LOG_LIK_COMPLETE Compute sum_m sum_i log P(x(i,m)| x(pi_i,m), theta_i) for a completely observed data set
% L = log_lik_complete(bnet, cases, clamped)
%
% If there is a missing data, you must use an inference engine.
% cases(i,m) is the value assigned to node i in case m.
% (If there are vector-valued nodes, cases should be a cell array.)
% clamped(i,m) = 1 if node i was set by intervention in case m (default: clamped = zeros)
% Clamped nodes contribute a factor of 1.0 to the likelihood.

if iscell(cases), usecell = 1; else usecell = 0; end

n = length(bnet.dag);
ncases = size(cases, 2);
if n ~= size(cases, 1)
  error('data should be of size nnodes * ncases');
end

if nargin < 3, clamped = zeros(n,ncases); end

L = 0;
for i=1:n
  ps = parents(bnet.dag, i);
  e = bnet.equiv_class(i);
  u = find(clamped(i,:)==0);
  ll = log_prob_node(bnet.CPD{e}, cases(i,u), cases(ps,u));
  if approxeq(exp(ll), 0),  fprintf('node %d has very low likelihood\n'); end
  L = L + ll;
end

