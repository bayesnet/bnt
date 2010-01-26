function [clpot, loglik] = enter_soft_evidence(engine, CPDpot, observed, pot_type, filter)
% ENTER_SOFT_EVIDENCE Add the specified soft evidence to the network (bk)
% [clpot, loglik] = enter_soft_evidence(engine, CPDpot, observed, pot_type, filter)

[ss T] = size(CPDpot);
C = length(engine.clusters);
S = length(engine.separators);
Q = length(cliques_from_engine(engine.sub_engine));
Q1 = length(cliques_from_engine(engine.sub_engine1));
clpot = cell(Q,T);
alpha = cell(C,T);

% Forwards
% The method is a generalization of the following HMM equation:
% alpha(j,t) = normalise( (sum_i alpha(i,t-1) * transmat(i,j)) * obsmat(j,t) )
% where alpha(j,t) = Pr(Q(t)=j | y(1:t))
t = 1;
[clpot(1:Q1,t), logscale(t)] = enter_soft_evidence(engine.sub_engine1, engine.clq_ass_to_node1(:), ...
					   CPDpot(:,1), find(observed(:,1)), pot_type);
for c=1:C
  k = engine.clq_ass_to_cluster1(c);
  alpha{c,t} = marginalize_pot(clpot{k,t}, engine.clusters{c});
end

%=== FH: For each separator s, divide some cluster potential by s's potential
alpha_orig = alpha(:,t);
for s=1:S
  c = engine.cluster_ass_to_separator(s);
  alpha{c,t} = divide_by_pot(alpha{c,t}, marginalize_pot(alpha_orig{c}, engine.separators{s}));
end

% For filtering, clpot{1} contains evidence on slice 1 only

%fprintf('alphas t=%d\n', t);
%for c=1:8
%  temp = pot_to_marginal(alpha{c,t});
%  temp.T
%end

% clpot{t} contains evidence from slices t-1, t for t > 1
clqs = [engine.clq_ass_to_cluster(:,1); engine.clq_ass_to_node(:,2)];
for t=2:T
  pots = [alpha(:,t-1); CPDpot(:,t)];
  [clpot(:,t), logscale(t)] = enter_soft_evidence(engine.sub_engine, clqs, pots, find(observed(:,t-1:t)),  pot_type);
  for c=1:C
    k = engine.clq_ass_to_cluster(c,2);
    cl = engine.clusters{c};
    alpha{c,t} = marginalize_pot(clpot{k,t}, cl+ss); % extract slice 2 posterior
    alpha{c,t} = set_domain_pot(alpha{c,t}, cl); % shift back to slice 1 for re-use as prior
  end
  %=== FH: For each separator s, divide some cluster potential by s's potential
  alpha_orig = alpha(:,t);
  for s=1:S
    c = engine.cluster_ass_to_separator(s);
    alpha{c,t} = divide_by_pot(alpha{c,t}, marginalize_pot(alpha_orig{c}, engine.separators{s}));
  end
end

loglik = sum(logscale); 

if filter
  return;
end

% Backwards
% The method is a generalization of the following HMM equation:
% beta(i,t) = (sum_j transmat(i,j) * obsmat(j,t+1) * beta(j,t+1))
% where beta(i,t) = Pr(y(t+1:T) | Q(t)=i)
t = T;
bnet = bnet_from_engine(engine);
beta = cell(C,T);
for c=1:C
  beta{c,t} = mk_initial_pot(pot_type, engine.clusters{c} + ss, bnet.node_sizes(:), bnet.cnodes(:), ...
			     find(observed(:,t-1:t)));
end
%=== FH: For each separator s, divide some cluster potential by s's potential
beta_orig = beta(:,t);
for s=1:S
  c = engine.cluster_ass_to_separator(s);
  beta{c,t} = divide_by_pot(beta{c,t}, marginalize_pot(beta_orig{c}, engine.separators{s}+ss));
end

for t=T-1:-1:1
  clqs = [engine.clq_ass_to_cluster(:,2); engine.clq_ass_to_node(:,2)];
  pots = [beta(:,t+1); CPDpot(:,t+1)];
  temp = enter_soft_evidence(engine.sub_engine, clqs, pots, find(observed(:,t:t+1)),  pot_type);
  for c=1:C
    k = engine.clq_ass_to_cluster(c,1);
    cl = engine.clusters{c};
    beta{c,t} = marginalize_pot(temp{k}, cl); % extract slice 1
    beta{c,t} = set_domain_pot(beta{c,t}, cl + ss); % shift fwd to slice 2
  end
  %=== FH: For each separator s, divide some cluster potential by s's potential
  beta_orig = beta(:,t);
  for s=1:S
    c = engine.cluster_ass_to_separator(s);
    beta{c,t} = divide_by_pot(beta{c,t}, marginalize_pot(beta_orig{c}, engine.separators{s}+ss));
  end
end

% Combine
% The method is a generalization of the following HMM equation:
% xi(i,j,t) = normalise( alpha(i,t) * transmat(i,j) * obsmat(j,t+1) * beta(j,t+1) )
% where xi(i,j,t) = Pr(Q(t)=i, Q(t+1)=j | y(1:T))
for t=1:T-1
  clqs = [engine.clq_ass_to_cluster(:); engine.clq_ass_to_node(:,2)];
  pots = [alpha(:,t); beta(:,t+1); CPDpot(:,t+1)];
  clpot(:,t+1) = enter_soft_evidence(engine.sub_engine, clqs, pots, find(observed(:,t:t+1)),  pot_type);
end
% for smoothing, clpot{1} is undefined
for k=1:Q1
  clpot{k,1} = []; 
end


