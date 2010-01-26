function [gamma, loglik, marginals, marginalsT] = bk_ff_fb(prior, transmat, obslik, filter_only, hnodes, ns)
% BK_FF_FB Fully factored Boyen-Koller version of forwards-backwards
% [gamma, loglik, marginals, marginalsT] = bk_ff_hmm(prior, transmat, obslik, filter_only, hnodes, ns)

ss  = length(ns);
S = length(prior);
T = size(obslik, 2);
marginals = cell(ss,T);
marginalsT = cell(ss,T);
scale = zeros(1,T);
alpha = zeros(S, T);

transmat2 = transmat';
for t=1:T
  if t==1
    [alpha(:,t), scale(t)] = normalise(prior(:) .* obslik(:,t));
  else
    [alpha(:,t), scale(t)] = normalise((transmat2 * alpha(:,t-1)) .* obslik(:,t));
  end
  [marginals(:,t), marginalsT(:,t)] = project_joint_onto_marginals(alpha(:,t), hnodes, ns);
  alpha(:,t) = combine_marginals_into_joint(marginalsT(:,t), hnodes, ns);
  %fprintf('alpha t=%d\n', t);
  %celldisp(marginals(1:8,t))
end
loglik = sum(log(scale));

if filter_only
  gamma = alpha;
  return;
end

beta = zeros(S,T);
gamma = zeros(S,T);
t = T;
beta(:,t) = ones(S,1);
gamma(:,t) = normalise(alpha(:,t) .* beta(:,t));
[marginals(:,t), marginalsT(:,t)] = project_joint_onto_marginals(gamma(:,t), hnodes, ns);

for t=T-1:-1:1
  b = beta(:,t+1) .* obslik(:,t+1); 
  beta(:,t) = normalise((transmat * b));
  [junk, tempT] = project_joint_onto_marginals(beta(:,t), hnodes, ns);
  beta(:,t) = combine_marginals_into_joint(tempT, hnodes, ns);
  %gamma(:,t) = normalise(alpha(:,t) .* beta(:,t));
  %[marginals(:,t), marginalsT(:,t)] = project_joint_onto_marginals(gamma(:,t), hnodes, ns);
end

gamma2 = zeros(S,T);
for t=T-1:-1:1
  b = beta(:,t+1) .* obslik(:,t+1); 
  xi(:,:,t) = normalise((transmat .* (alpha(:,t) * b')));      
  if t==T-1
    gamma2(:,T) = sum(xi(:,:,T-1), 1)';
  end
  gamma2(:,t) = sum(xi(:,:,t), 2);
  [marginals(:,t), marginalsT(:,t)] = project_joint_onto_marginals(gamma2(:,t), hnodes, ns);
end


