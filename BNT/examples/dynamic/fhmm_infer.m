function [loglik, gamma] = fhmm_infer(inter, CPTs_slice1, CPTs, obsmat, node_sizes)
% FHMM_INFER Exact inference for a factorial HMM.
% [loglik, gamma] = fhmm_infer(inter, CPTs_slice1, CPTs, obsmat, node_sizes)
%
% Inputs:
% inter - the inter-slice adjacency matrix
% CPTs_slice1{s}(j) = Pr(Q(s,1) = j) where Q(s,t) = hidden node s in slice t
% CPT{s}(i1, i2, ..., j) = Pr(Q(s,t) = j | Pa(s,t-1) = i1, i2, ...),
% obsmat(i,t) = Pr(y(t) | Q(t)=i)
% node_sizes is a vector with the cardinality of the hidden nodes
%
% Outputs:
% gamma(i,t) = Pr(X(t)=i | O(1:T)) as in an HMM,
% except that i is interpreted as an M digit, base-K number (if there are M chains each of cardinality K).
%
%
% For M chains each of cardinality K, the frontiers  (i.e., cliques)
% contain M+1 nodes, and it takes M steps to advance the frontier by one time step,
% so the run time is O(T M K^(M+1)).
% An HMM takes O(T S^2) where S is the size of the state space.
% Collapsing the FHMM to an HMM results in S = K^M.
% For details, see
%   "The Factored Frontier Algorithm for Approximate Inference in DBNs",
%    Kevin Murphy and Yair Weiss, submitted to NIPS 2000.
%
% The frontier algorithm makes the following topological assumptions:
% 
%  - All nodes are persistent (connect to the next slice)
%  - No connections within a timeslice
%  - There is a single observation variable, which depends on all the hidden nodes
%  - Each node can have several parents in the previous time slice (generalizes a FHMM slightly)
%

% The forwards pass of the frontier algorithm can be explained with the following example.
% Suppose we have 3 hidden nodes per slice, A, B, C.
% The goal is to compute alpha(j, t) = Pr( (A_t,B_t,C_t)=j | Y(1:t))
% We move alpha from t to t+1 one node at a time, as follows.
% We define the following quantities:
% s([a1 b1 c1], 1) = Prob(A(t)=a1, B(t)=b1, C(t)=c1 | Y(1:t)) = alpha(j, t)
% s([a2 b1 c1], 2) = Prob(A(t+1)=a2, B(t)=b1, C(t)=c1 | Y(1:t))
% s([a2 b2 c1], 3) = Prob(A(t+1)=a2, B(t+1)=b2, C(t)=c1 | Y(1:t))
% s([a2 b2 c2], 4) = Prob(A(t+1)=a2, B(t+1)=b2, C(t+1)=c2 | Y(1:t))
% s([a2 b2 c2], 5) = Prob(A(t+1)=a2, B(t+1)=b2, C(t+1)=c2 | Y(1:t+1)) = alpha(j, t+1)
%
% These can be computed recursively as follows:
%
% s([a2 b1 c1], 2) = sum_{a1} P(a2|a1) s([a1 b1 c1], 1)
% s([a2 b2 c1], 3) = sum_{b1} P(b2|b1) s([a2 b1 c1], 2)
% s([a2 b2 c2], 4) = sum_{c1} P(c2|c1) s([a2 b2 c1], 1)
% s([a2 b2 c2], 5) = normalise( s([a2 b2 c2], 4) .* P(Y(t+1)|a2,b2,c2)


[kk,ll,mm] = make_frontier_indices(inter, node_sizes); % can pass in as args

scaled = 1;

M = length(node_sizes);
S = prod(node_sizes);
T = size(obsmat, 2);

alpha = zeros(S, T);
beta = zeros(S, T);
gamma = zeros(S, T);
scale = zeros(1,T);
tiny = exp(-700);


alpha(:,1) = make_prior_from_CPTs(CPTs_slice1, node_sizes);
alpha(:,1) = alpha(:,1) .* obsmat(:, 1);

if scaled
  s = sum(alpha(:,1));
  if s==0, s = s + tiny; end
  scale(1) = 1/s;
else
  scale(1) = 1;
end
alpha(:,1) = alpha(:,1) * scale(1);

%a = zeros(S, M+1);
%b = zeros(S, M+1);
anew = zeros(S,1);
aold = zeros(S,1);
bnew = zeros(S,1);
bold = zeros(S,1);

for t=2:T
  %a(:,1) = alpha(:,t-1);
  aold =  alpha(:,t-1);
  
  c = 1;
  for i=1:M
    ns = node_sizes(i);
    cpt = CPTs{i};
    for j=1:S
      s = 0;
      for xx=1:ns
	%k = kk(xx,j,i);
	%l = ll(xx,j,i);
	k = kk(c);
	l = ll(c);
	c = c + 1;
	% s = s + a(k,i) * CPTs{i}(l);
	s = s + aold(k) * cpt(l);
      end
      %a(j,i+1) = s;
      anew(j) = s;
    end
    aold = anew;
  end
  
  %alpha(:,t) = a(:,M+1) .* obsmat(:, obs(t));
  alpha(:,t) = anew .* obsmat(:, t);

  if scaled
    s = sum(alpha(:,t));
    if s==0, s = s + tiny; end
    scale(t) = 1/s;
  else
    scale(t) = 1;
  end
  alpha(:,t) = alpha(:,t) * scale(t);

end


beta(:,T) = ones(S,1) * scale(T);
for t=T-1:-1:1
  %b(:,1) = beta(:,t+1) .* obsmat(:, obs(t+1));
  bold = beta(:,t+1) .* obsmat(:, t+1);

  c = 1;
  for i=1:M
    ns = node_sizes(i);
    cpt = CPTs{i};
    for j=1:S
      s = 0;
      for xx=1:ns
	%k = kk(xx,j,i);
	%m = mm(xx,j,i);
	k = kk(c);
	m = mm(c);
	c = c + 1;
	% s = s + b(k,i) * CPTs{i}(m);
	s = s + bold(k) * cpt(m);
      end
      %b(j,i+1) = s;
      bnew(j) = s;
    end
    bold = bnew;
  end
  % beta(:,t) = b(:,M+1) * scale(t);
  beta(:,t) = bnew * scale(t);
end


if scaled
  loglik = -sum(log(scale)); % scale(i) is finite
else
  lik = alpha(:,1)' * beta(:,1);
  loglik = log(lik+tiny);
end

for t=1:T
  gamma(:,t) = normalise(alpha(:,t) .* beta(:,t));
end

%%%%%%%%%%%

function [kk,ll,mm] = make_frontier_indices(inter, node_sizes)
%
% Precompute indices for use in the frontier algorithm.
% These only depend on the topology, not the parameters or data.
% Hence we can compute them outside of fhmm_infer.
% This saves a lot of run-time computation.

M = length(node_sizes);
S = prod(node_sizes);

mns = max(node_sizes);
kk = zeros(mns, S, M);
ll = zeros(mns, S, M);
mm = zeros(mns, S, M);

for i=1:M
  for j=1:S
    u = ind2subv(node_sizes, j);
    x = u(i);
    for xx=1:node_sizes(i)
      uu = u;
      uu(i) = xx;
      k = subv2ind(node_sizes, uu);
      kk(xx,j,i) = k;
      ps = find(inter(:,i)==1);
      ps = ps(:)';
      l = subv2ind(node_sizes([ps i]), [uu(ps) x]); % sum over parent
      ll(xx,j,i) = l;
      m = subv2ind(node_sizes([ps i]), [u(ps) xx]); % sum over child
      mm(xx,j,i) = m;
    end
  end
end

%%%%%%%%%

function prior=make_prior_from_CPTs(indiv_priors, node_sizes)
%
% composite_prior=make_prior(individual_priors, node_sizes)
% Make the prior for the first node in a Markov chain
% from the priors on each node in the equivalent DBN.
% prior{i}(j) = Pr(X_i=j), where X_i is the i'th node in slice 1.
% composite_prior(i) = Pr(slice1 = i).

n = length(indiv_priors);
S = prod(node_sizes);
prior = zeros(S,1);
for i=1:S
  vi = ind2subv(node_sizes, i);
  p = 1;
  for k=1:n
    p = p * indiv_priors{k}(vi(k));
  end
  prior(i) = p;
end



%%%%%%%%%%%

function [loglik, alpha, beta] = FHMM_slow(inter, CPTs_slice1, CPTs, obsmat, node_sizes, data)
% 
% Same as the above, except we don't use the optimization of computing the indices outside the loop.


scaled = 1;

M = length(node_sizes);
S = prod(node_sizes);
[numex T] = size(data);

obs = data;

alpha = zeros(S, T);
beta = zeros(S, T);
a = zeros(S, M+1);
b = zeros(S, M+1);
scale = zeros(1,T);

alpha(:,1) = make_prior_from_CPTs(CPTs_slice1, node_sizes);
alpha(:,1) = alpha(:,1) .* obsmat(:, obs(1));
if scaled
  s = sum(alpha(:,1));
  if s==0, s = s + tiny; end
  scale(1) = 1/s;
else
  scale(1) = 1;
end
alpha(:,1) = alpha(:,1) * scale(1);

for t=2:T
  fprintf(1, 't %d\n', t);
  a(:,1) = alpha(:,t-1);
  for i=1:M
    for j=1:S
      u = ind2subv(node_sizes, j);
      xnew = u(i);
      s = 0;
      for xold=1:node_sizes(i)
	uold = u;
	uold(i) = xold;
	k = subv2ind(node_sizes, uold);
	ps = find(inter(:,i)==1);
	ps = ps(:)';
	l = subv2ind(node_sizes([ps i]), [uold(ps) xnew]);
	s = s + a(k,i) * CPTs{i}(l);
      end
      a(j,i+1) = s;
    end
  end
  alpha(:,t) = a(:,M+1) .* obsmat(:, obs(t));

  if scaled
    s = sum(alpha(:,t));
    if s==0, s = s + tiny; end
    scale(t) = 1/s;
  else
    scale(t) = 1;
  end
  alpha(:,t) = alpha(:,t) * scale(t);

end


beta(:,T) = ones(S,1) * scale(T);
for t=T-1:-1:1
  fprintf(1, 't %d\n', t);
  b(:,1) = beta(:,t+1) .* obsmat(:, obs(t+1));
  for i=1:M
    for j=1:S
      u = ind2subv(node_sizes, j);
      xold = u(i);
      s = 0;
      for xnew=1:node_sizes(i)
	unew = u;
	unew(i) = xnew;
	k = subv2ind(node_sizes, unew);
	ps = find(inter(:,i)==1);
	ps = ps(:)';
	l = subv2ind(node_sizes([ps i]), [u(ps) xnew]);
	s = s + b(k,i) * CPTs{i}(l);
      end
      b(j,i+1) = s;
    end
  end
  beta(:,t) = b(:,M+1) * scale(t);
end


if scaled
  loglik = -sum(log(scale)); % scale(i) is finite
else
  lik = alpha(:,1)' * beta(:,1);
  loglik = log(lik+tiny);
end
